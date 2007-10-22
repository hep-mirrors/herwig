// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DefaultReweighter class.
//

#include "DefaultReweighter.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Repository/EventGenerator.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "DefaultReweighter.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#ifdef HERWIG_DEBUG_CKKW
#include "ThePEG/Repository/Repository.h"
#endif

using namespace Herwig;

DefaultReweighter::~DefaultReweighter() {}

void DefaultReweighter::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
  os << _splittingMap << _sudakovMap << ounit(_interpolationSpacing,GeV2) 
     << _integrationAccuracy << ounit(_sudakovMaxScale,GeV2)
     << _sudakovDataPath << _useMassive << _sudakovUnweight;
}

void DefaultReweighter::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
  is >> _splittingMap >> _sudakovMap >> iunit(_interpolationSpacing,GeV2)
     >> _integrationAccuracy >> iunit(_sudakovMaxScale,GeV2)
     >> _sudakovDataPath >> _useMassive >> _sudakovUnweight;
}

AbstractClassDescription<DefaultReweighter> DefaultReweighter::initDefaultReweighter;
// Definition of the static class description member.

void DefaultReweighter::Init() {

  static ClassDocumentation<DefaultReweighter> documentation
    ("DefaultReweighter performs the Sudakov reweighting for standard CKKW approaches.");




  static Parameter<DefaultReweighter,Energy2> interfaceInterpolationSpacing
    ("InterpolationSpacing",
     "The spacing between interpolation points.",
     &DefaultReweighter::_interpolationSpacing, GeV2, 1.0*GeV2, 0.0*GeV2, 0*GeV2,
     false, false, Interface::lowerlim);


  static Parameter<DefaultReweighter,double> interfaceIntegrationAccuracy
    ("IntegrationAccuracy",
     "The relative accuracy to which Sudakov exponents are to be computed.",
     &DefaultReweighter::_integrationAccuracy, 1e-3, 0, 1,
     false, false, Interface::limited);


  static Parameter<DefaultReweighter,Energy2> interfaceSudakovMaxScale
    ("SudakovMaxScale",
     "The maximum scale up to which interpolation poiints for Sudakov form factors are to be computed.",
     &DefaultReweighter::_sudakovMaxScale, GeV2, 1000000.0*GeV2, 0.0*GeV2, 0*GeV2,
     false, false, Interface::lowerlim);


  static Parameter<DefaultReweighter,string> interfaceSudakovDataPath
    ("SudakovDataPath",
     "The path where Sudakov interpolation data is stored.",
     &DefaultReweighter::_sudakovDataPath, "./",
     false, false);


  static Switch<DefaultReweighter,bool> interfaceUseMassive
    ("UseMassive",
     "Wether or not to use massive splitting functions",
     &DefaultReweighter::_useMassive, true, false, false);
  static SwitchOption interfaceUseMassiveUseMassiveOn
    (interfaceUseMassive,
     "On",
     "Use massive splitting functions",
     true);
  static SwitchOption interfaceUseMassiveUseMassiveOff
    (interfaceUseMassive,
     "Off",
     "Do not use massive splitting functions",
     false);


  static Switch<DefaultReweighter,bool> interfaceSudakovUnweight
    ("SudakovUnweight",
     "Wether or not to use the Sudakov weight divided by the maximum Sudakov weight.",
     &DefaultReweighter::_sudakovUnweight, true, false, false);
  static SwitchOption interfaceSudakovUnweightOn
    (interfaceSudakovUnweight,
     "On",
     "Divide Sudakov weights by maximum weight",
     true);
  static SwitchOption interfaceSudakovUnweightOff
    (interfaceSudakovUnweight,
     "Off",
     "Do not divide Sudakov weight by the maximum Sudakov weight",
     false);

}

double DefaultReweighter::sudakovReweight (CascadeHistory history, unsigned int mult, unsigned int minmult) {

#ifdef HERWIG_DEBUG_CKKW_EXTREME
  generator()->log() << "== DefaultReweighter::sudakovReweight" << endl;
#endif

  double weight = 1.;

  if (mult == minmult && _sudakovUnweight) return weight;

  for (vector<ClusteringParticlePtr>::iterator p = history.reconstructed.begin();
       p != history.reconstructed.end(); ++p) {
    if (!(**p).noReweight() && (**p).pData().colour != (**p).pData().antiColour) {
#ifdef HERWIG_DEBUG_CKKW_EXTREME
      generator()->log() << "reweighting for clustering particle "
			 << *p << endl
			 << "from Q/GeV = " << sqrt((**p).productionScale())/GeV
			 << " to q/GeV = " << sqrt((**p).splittingScale())/GeV << endl;
#endif
      pair<long,bool> sudakovKey;
      sudakovKey.first = (**p).pData().partonId.PDGId;
      // for antiquarks use quark sudakov
      if (abs(sudakovKey.first)<7 && sudakovKey.first<0) sudakovKey.first = -sudakovKey.first;
      sudakovKey.second = ((**p).pData().partonId.state == ClusteringParticleState::initial);
      bool gotone=false;
      for(multimap<pair<long,bool>,DefaultSudakovPtr>::iterator s = _sudakovMap.lower_bound(sudakovKey);
	  s != _sudakovMap.upper_bound(sudakovKey); ++s) {
	if (!gotone) gotone = true;
	double sweight = (*(s->second))((**p).productionScale(),(**p).splittingScale());
#ifdef HERWIG_DEBUG_CKKW_EXTREME
	generator()->log() << "using sudakov " << s->second << " weight = " << sweight << endl;
#endif
	weight *= sweight;
      }
      if (!gotone && ((**p).pData().colour != (**p).pData().antiColour))
	generator()->log() << "CKKW : DefaultReweighter::sudakovReweight() : No Sudakov"
			   << " form factor found for " << (**p).pData().partonId.PDGId
			   << " " << (**p).pData().partonId.state << endl;
    }
  }

  // now divide by the maximum weight, i.e. the weight associated with
  // the hard process

  if (_sudakovUnweight) {

#ifdef HERWIG_DEBUG_CKKW_EXTREME
    generator()->log() << "dividing by hard process weight" << endl;
#endif

    for (vector<tClusteringParticlePtr>::iterator p = history.hardProcess.begin();
	 p != history.hardProcess.end(); ++p) {
      if ((**p).pData().colour != (**p).pData().antiColour) {
#ifdef HERWIG_DEBUG_CKKW_EXTREME
	generator()->log() << "reweighting for clustering particle "
			   << *p << endl
			   << "from Q/GeV = " << sqrt((**p).productionScale())/GeV << endl;
#endif
	pair<long,bool> sudakovKey;
	sudakovKey.first = (**p).pData().partonId.PDGId;
	// for antiquarks use quark sudakov
	if (abs(sudakovKey.first)<7 && sudakovKey.first<0) sudakovKey.first = -sudakovKey.first;
	sudakovKey.second = ((**p).pData().partonId.state == ClusteringParticleState::initial);
	bool gotone=false;
	for(multimap<pair<long,bool>,DefaultSudakovPtr>::iterator s = _sudakovMap.lower_bound(sudakovKey);
	    s != _sudakovMap.upper_bound(sudakovKey); ++s) {
	  if (!gotone) gotone = true;
	  double sweight = (*(s->second))((**p).productionScale(),
					  resolution()
					  ->minResolvableScale(sudakovKey.first,sudakovKey.second));
#ifdef HERWIG_DEBUG_CKKW_EXTREME
	  generator()->log() << "using sudakov " << s->second << " weight = " << sweight << endl;
#endif
	  weight /= sweight;
	}
	if (!gotone && ((**p).pData().colour != (**p).pData().antiColour))
	  generator()->log() << "CKKW : DefaultReweighter::sudakovReweight() : No Sudakov"
			     << " form factor found for " << (**p).pData().partonId.PDGId
			     << " " << (**p).pData().partonId.state << endl;
      }
    }

  }

  if (weight > 1.) throw Exception() << "CKKW : DefaultReweighter::sudakovReweight() : weight larger than one"
				     << Exception::eventerror;
  return weight;

}

void DefaultReweighter::setup (CascadeReconstructorPtr reconstructor) {

#ifdef HERWIG_DEBUG_CKKW
  Repository::clog() << "== DefaultReweighter::setup" << endl
		     << "using " << _splittingMap.size() << " splittings" << endl;
#endif

  Reweighter::setup (reconstructor);

  // make sure alpha_s is initialized, befor we start to calculate Sudakovs
  showerAlpha()->initialize();

  // now build Sudakovs from splitting functions.

  DefaultSudakovPtr sudakov;

  for (map<splittingKey,SplittingFnPtr>::iterator s = _splittingMap.begin();
       s != _splittingMap.end(); ++s) {

    if (!s->second)
      Throw<InitException> () << "CKKW : DefaultReweighter::setup() : null pointer to splitting function";

#ifdef HERWIG_DEBUG_CKKW
    Repository::clog() << "Building Sudakov for "
		       << (s->first.initial ? "initial " : "final ")
		       << "state splitting ";
    for(IdList::const_iterator i = s->first.ids.begin(); i != s->first.ids.end(); ++i)
      Repository::clog() << *i << " ";
    Repository::clog() << "\n"
		       << "using splitting function "
		       << s->second->name() << "(" << s->second << ")" << endl;
#endif
    sudakov = new_ptr(DefaultSudakov(this,s->second,s->first.ids,_useMassive,s->first.initial));
    sudakov->initialize(false);
    _sudakovMap.insert(make_pair(make_pair(s->first.ids[0],s->first.initial),sudakov));
  }

}

void DefaultReweighter::initialize () {

#ifdef HERWIG_DEBUG_CKKW
  Repository::clog() << "== DefaultReweighter::initialize" << endl;
#endif

  // initialize sudakovs

  for (multimap<pair<long,bool>, DefaultSudakovPtr>::iterator s = _sudakovMap.begin();
       s != _sudakovMap.end(); ++s) {
#ifdef HERWIG_DEBUG_CKKW
    Repository::clog() << "initializing sudakov " << s->second
		       << " for reweighting " << (s->first.second ? "initial" : "final")
		       << " PDGId " << s->first.first
		       << " using splitting function " << s->second->splittingFunction()->name() << endl;
#endif
    s->second->initialize(true);
  }

}
