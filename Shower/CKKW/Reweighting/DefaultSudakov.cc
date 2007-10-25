// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DefaultSudakov class.
//

#include "DefaultSudakov.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Repository/Repository.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Config/Constants.h"

#include "DefaultReweighter.h"
#include "DefaultJetMeasure.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "DefaultSudakov.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

  namespace gsl_interface {

    tDefaultJetMeasurePtr _resolution;
    tSplittingFnPtr _splittingFunction;
    tShowerAlphaPtr _showerAlpha;
    IdList _ids;
    bool _useMassive;
    bool _initial;

    double _accuracy;

    double _integrand (double z, void * vq2) {

      Energy2 q2 = (*(double*)(vq2)) * GeV2;

      pair<Energy2,double> shower = _resolution->invertClustering(_ids,q2,z,_initial);

      double value = 1.;

      if (!_initial) {
	value = _showerAlpha->value(sqr(shower.second*(1.-shower.second))*shower.first)/(2.*Constants::pi);
	value *= _splittingFunction->P(shower.second,shower.second*(1.-shower.second)*shower.first,_ids,_useMassive);
	value *= _resolution->showerJacobian(_ids,q2,z,_initial)*GeV2/shower.first;
      } else {
	value = _showerAlpha->value(sqr(1.-shower.second)*shower.first)/(2.*Constants::pi);
	value *= _splittingFunction->P(shower.second,(1.-shower.second)*shower.first/shower.second,_ids,_useMassive);
	value *= _resolution->showerJacobian(_ids,q2,z,_initial)*GeV2/shower.first;
      }

      return value;

    }

    double _z_integral (double q2, void * vQ2) {

      Energy2 Q2 = (*(double*)(vQ2)) * GeV2;

      pair<double,double> zlims = _resolution->zLimits(q2*GeV2,Q2,_ids,_initial);

      gsl_function integrand;
      integrand.function = &_integrand;
      integrand.params = &q2;

      double value, err;

      gsl_integration_workspace * workspace = gsl_integration_workspace_alloc(1000);

      gsl_integration_qag (&integrand,zlims.first,zlims.second,
			   0,_accuracy,1000,GSL_INTEG_GAUSS61,workspace,
			   &value,&err);

      gsl_integration_workspace_free (workspace);

      return value;


    }

  }

DefaultSudakov::~DefaultSudakov() {
  if (_data_allocated) {
    delete _qvalues;
    delete _ivalues;
  }
  if (_interpolation_allocated) gsl_spline_free (_spline);
  if (_interpolation_allocated) gsl_interp_accel_free (_acc);
}

void DefaultSudakov::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
  os << _reweighter << _splittingFunction << _ids << _initial << _useMassive;
}

void DefaultSudakov::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
  is >> _reweighter >> _splittingFunction >> _ids >> _initial >> _useMassive;
}

ClassDescription<DefaultSudakov> DefaultSudakov::initDefaultSudakov;
// Definition of the static class description member.

void DefaultSudakov::Init() {

  static ClassDocumentation<DefaultSudakov> documentation
    ("DefaultSudakov is the class which performs all the "
     "numerics for Sudakov reweighting in a standard CKKW "
     "approach.");

}

string DefaultSudakov::dataFile () const {
  ostringstream name;

  name << _reweighter->sudakovDataPath() << "/"
       << _splittingFunction->name() << "_";

  for (IdList::const_iterator i = _ids.begin(); i != _ids.end(); ++i)
    name << *i << "_";

  name << _initial << "_"
       << _reweighter->resolution()->name() << "_"
       << _reweighter->resolution()->resolution()/GeV2 << "_"
       << _reweighter->showerAlpha()->value(_reweighter->resolution()->resolution())
       << ".dat";

  return name.str();
}

void DefaultSudakov::initialize (bool run) {

#ifdef HERWIG_DEBUG_CKKW_EXTREME
  if (run) {
    Repository::clog() << "== DefaultSudakov::initialize (run)" << endl;
    Repository::clog() << flush;
  }
  else
    Repository::clog() << "== DefaultSudakov::initialize (init)" << endl;
#endif

  bool dointegration = false;

  ifstream test (dataFile().c_str());

  // if not in run mode, look for data file.
  // if not existing, compute the interpolation points.

  if (!run) {
    if (!test) {
      Repository::clog() 
	<< "CKKW : DefaultSudakov: "
	<< "Interpolation data file\n" + dataFile () + "\nmissing, computing interpolation points. This may take some time."
	<< endl;
      test.close ();
      dointegration = true;
    } else {
      test.close();
    }
  }

  // if in run mode and no data file, abort

  if (run && !test) {
    test.close();
    throw Exception() << "CKKW : DefaultSudakov::initialize () : "
		      << "data file\n" + dataFile() + "\nmissing." << Exception::runerror;
  }

  // if the jet resolution can't handle the branching
  // abort the run

  if (!dynamic_ptr_cast<DefaultJetMeasurePtr>(_reweighter->resolution())->canHandle(_ids,_initial))
    Throw<InitException>() << "CKKW : DefaultSudakov::initialize () : "
			   << "jet resolution " << _reweighter->resolution()->name()
			   << " cannot handle the branching for computing "
			   << dataFile();

  // allocate double arrays for gsl

  unsigned long interpolationPoints =
    (unsigned long)((_reweighter->sudakovMaxScale()
		     - _reweighter->resolution()->minResolvableScale(_ids[0],_initial))
		  /_reweighter->interpolationSpacing());

  _qvalues = new double [interpolationPoints];
  _ivalues = new double [interpolationPoints];

  _data_allocated = true;

  // if we are to do the integrations, calculate and write out interpolation points

  if (dointegration) {

#ifdef HERWIG_DEBUG_CKKW_EXTREME
    Repository::clog() << "starting integration" << endl;
    Repository::clog() << flush;
#endif

    gsl_integration_workspace * workspace = gsl_integration_workspace_alloc(1000);

    // do we need a mutex here?

    gsl_interface::_resolution = dynamic_ptr_cast<DefaultJetMeasurePtr>(_reweighter->resolution());
    gsl_interface::_splittingFunction = _splittingFunction;
    gsl_interface::_showerAlpha = _reweighter->showerAlpha();
    gsl_interface::_ids = _ids;
    gsl_interface::_initial = _initial;
    gsl_interface::_useMassive = _useMassive;
    gsl_interface::_accuracy = _reweighter->integrationAccuracy();

    gsl_function sintegrand;

    sintegrand.function = &gsl_interface::_z_integral;

    double stepLength = _reweighter->interpolationSpacing()/GeV2;

    ofstream odata (dataFile().c_str());

    for (unsigned int i = 0; i<interpolationPoints; i++) {
      _qvalues[i] =
	_reweighter->resolution()->minResolvableScale (_ids[0],_initial)/GeV2 + i*stepLength;

      double value, err;

      if (i%(interpolationPoints/80)==0) Repository::clog() << ".";

      if (i==0) { value = 0.; err = 0.;  }

      if (i !=0) {

	sintegrand.params = &_qvalues[i];

	gsl_integration_qag (&sintegrand,
			     _reweighter->resolution()->minResolvableScale (_ids[0],_initial)/GeV2,
			     _qvalues[i],
			     0,_reweighter->integrationAccuracy(),1000,GSL_INTEG_GAUSS61,workspace,
			     &value,&err);

      }

      _ivalues[i] = value;

      odata << _qvalues[i] << "\t" << _ivalues[i] << "\n";

    }

    Repository::clog() << endl;

    gsl_integration_workspace_free (workspace);

  }

  // in run mode, we load the data file
  // some exception handling missing here

  if (run) {

#ifdef HERWIG_DEBUG_CKKW_EXTREME
    Repository::clog() << "loading interpolation data" << endl;
#endif

    ifstream indata (dataFile().c_str());

    for (unsigned int i = 0; i<interpolationPoints; i++) {
      indata >> _qvalues[i] >> _ivalues[i];
    }

    // continue to initialize the interpolation
    
    _acc = gsl_interp_accel_alloc ();
    _spline = gsl_spline_alloc (gsl_interp_cspline, interpolationPoints);
    _interpolation_allocated = true;

    gsl_spline_init (_spline, _qvalues, _ivalues, interpolationPoints);

  }

}

double DefaultSudakov::interpolate (Energy2 q) {
  if (q> _reweighter->sudakovMaxScale()) throw Exception() << "CKKW : DefaultSudakov::interpolate(Energy2) : "
							   << "scale exceeds maximum scale" << Exception::eventerror;

  if (q< _reweighter->resolution()->minResolvableScale (_ids[0],_initial)) return 1.;

  double exponent = gsl_spline_eval (_spline, q/GeV2, _acc);

#ifdef HERWIG_DEBUG_CKKW_CHECK_SUDAKOVS
  _sudakov_calls.push_back(make_pair(q,exponent));
#endif

  return exp(-exponent);
}

#ifdef HERWIG_DEBUG_CKKW_CHECK_SUDAKOVS

void DefaultSudakov::dumpSudakovCalls() {
  string checkData = dataFile() + ".check";
  ofstream out (checkData.c_str());
  for (vector<pair<Energy2,double> >::iterator c = _sudakov_calls.begin();
       c != _sudakov_calls.end(); ++c) {
    out << c->first/GeV2 << " " << c->second << endl;
  }
}

#endif
