// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SudakovFormFactor class.
//

#include "SudakovFormFactor.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include <cassert>

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "SudakovFormFactor.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

SudakovFormFactor::~SudakovFormFactor() {}

void SudakovFormFactor::persistentOutput(PersistentOStream & os) const {
  os << _splittingFn << _alpha << _variables << _pdfmax;
}

void SudakovFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> _splittingFn >> _alpha >> _variables >> _pdfmax;
}

ClassDescription<SudakovFormFactor> SudakovFormFactor::initSudakovFormFactor;
// Definition of the static class description member.

void SudakovFormFactor::Init() {

  static ClassDocumentation<SudakovFormFactor> documentation
    ("Class that, given a splitting function, returns values for Sudakov.");

  static Reference<SudakovFormFactor,SplittingFunction>
    interfaceSplittingFunction("SplittingFunction",
			       "A reference to the SplittingFunction object",
			       &Herwig::SudakovFormFactor::_splittingFn,
			       false, false, true, false);

  static Reference<SudakovFormFactor,ShowerAlpha>
    interfaceAlpha("Alpha",
		   "A reference to the Alpha object",
		   &Herwig::SudakovFormFactor::_alpha,
		   false, false, true, false);

  static Parameter<SudakovFormFactor,double> interfacePDFmax
    ("PDFmax",
     "Maximum value of PDF weight. ",
     &SudakovFormFactor::_pdfmax, 35.0, 1.0, 1000.0,
     false, false, Interface::limited);
}

void SudakovFormFactor::setupLookupTables() {}

void SudakovFormFactor::
guessTimeLike(double &z, double &z0, double &z1, Energy2 &t, const Energy2 t0,
	      const Energy kinCutoff, const bool glueEmits) const {
  Energy2 told = t;
  // the larger PS-boundary in z (could be part of ShowerKinematics, really!)
  if (glueEmits) 
    {
      // no emission possible
      if(t<16.*t0)
	{
	  t=-1.;
	  return;
	}
      z0 = 0.5*(1.-sqrt(1.-4.*sqrt(t0/t)));
      z1 = 1.-z0;
    } 
  else 
    {
      z0 = 0.5*sqrt(t0/t);
      z1 = 1.-0.5*kinCutoff/sqrt(t);
    }
  // guess values of t and z
  t = guesst(told, z0, z1); 
  z = guessz(z0, z1); 
  // actual values for z-limits
  if (glueEmits) 
    {
      z0 = 0.5*(1.-sqrt(1.-4.*sqrt(t0/t)));
      z1 = 1.-z0;
    } 
  else 
    {
      z0 = 0.5*sqrt(t0/t)/2.;
      z1 = 1.-0.5*kinCutoff/sqrt(t);
    }
} 

void SudakovFormFactor::
guessSpaceLike(double &z, double &z0, double &z1, Energy2 &t, 
	       const Energy2 t0, const Energy kinCutoff, const double x) 
  const {
  
  Energy2 told = t;
  // the overestimated PS-boundary in z
  z0 = x;
  double yy = 1.+sqr(kinCutoff)/t/2.;
  z1 = yy - sqrt(sqr(yy)-1.); 
  if (z1 < z0) {
    t = -1.0*GeV;
    // we can return here, 
    // if t=-1 the calling function will return q=-1 anyway
    // no matter what the other variables are
    // dgrell: look at logic of this calling stack again
    return;
  }
  // guess values of t and z
  t = guesst(told, z0, z1, true); 
  z = guessz(z0, z1); 
  // actual values for z-limits
  yy = 1.+sqr(kinCutoff)/t/2.;
  z1 = yy - sqrt(sqr(yy)-1.); 
  // if new upper limit less than lower
  if (z1 < z0) t = -1.0*GeV;
} 

bool SudakovFormFactor::PSVeto(const double z, const double z0, 
			       const double z1, const Energy2 t, 
			       const Energy2 tmin, const Energy2 t0,
			       const Energy kinCutoff, bool glueEmits) {
  // still inside PS?
  if((z < z0 || z > z1) && t > tmin) return true;
  else {
    // still REALLY inside PS? (overestimated allowed PS)
    if (glueEmits) {if (sqr(z*(1.-z))*t < t0) return true;}
    else if (sqr(1.-z)*(t*sqr(z) - t0) < z*sqr(kinCutoff)) return true;
  }
  return false;
}


 
bool SudakovFormFactor::PDFVeto(const double z, const Energy2 t, 
				const double x,
				const tcPDPtr parton0, 
				const tcPDPtr parton1) const {
  assert(_variables->currentPDF());
  // remember: pdf's q is cut in pdf class.  shoudl probably be done here! 
  // this would correspond to QSPAC in F-HERWIG. 
  double ratio = 
    _variables->currentPDF()->xfx(_variables->beamParticle(),parton0,t,x/z)/
    _variables->currentPDF()->xfx(_variables->beamParticle(),parton1,t,x);
  // ratio / PDFMax must be a probability <= 1.0
  if (ratio > _pdfmax) {
    generator()->log() << "PDFVeto warning: Ratio (" << ratio 
			    << ") > " << name() << ":PDFmax ("
			    <<_pdfmax <<")\n";
  }
  return ratio < UseRandom::rnd()*_pdfmax;
}
 
Energy SudakovFormFactor::generateNextTimeBranching(const Energy startingScale,
						    const IdList &ids,
						    const bool reverseAO)
{
  // First reset the internal kinematics variables that can
  // have been eventually set in the previous call to thie method.
  _q = Energy();
  _z = 0.0;
  _phi = 0.0; 
  // perform initialisation
  initialize(ids);
  // reverse angular ordering isn't implemented
  if (reverseAO) 
    throw Exception() << "Reverse Ordering not implemented in " 
		      << "SudakovFormFactor::generateNextTimeBranching"
		      << Exception::runerror;
  // normal ordering
  double z0,z1;
  Energy2 t,t0,tmax,tmin;
  Energy kinCutoff;


  bool glueEmits = (ids[0]==ParticleID::g);
  Energy m0 = getParticleData(ids[0])->mass();
  Energy m1 = getParticleData(ids[1])->mass();




  tmax = sqr(startingScale);
  t = tmax;
  if(glueEmits) initialize(t0,tmin,tmax,kinCutoff,m1);
  else initialize(t0,tmin,tmax,kinCutoff,m0);
  if(_q < 0.) return _q;
  do 
    { 
      guessTimeLike(_z,z0,z1,t,t0,kinCutoff,glueEmits);
      // Our t is too low now, terminate guesses
      if(tVeto(t,tmin)) break;
    } 
  while(PSVeto(_z,z0,z1,t,tmin,t0,kinCutoff,glueEmits) ||
	SplittingFnVeto(_z,t,ids) || 
	alphaSVeto(sqr(_z*(1.-_z))*t));
  if(t > 0) _q = sqrt(t);
  else _q = -1.;
  _phi = 2.*pi*UseRandom::rnd();
  return _q;
}

Energy SudakovFormFactor::generateNextSpaceBranching(const Energy startingQ,
			                             const IdList &ids,
						     double x,
						     const bool revOrd) {
  _q = Energy();
  _z = 0.0;
  _phi = 0.0;
  double z0,z1;
  // reverse angular ordering isn't implemented
  if(revOrd) 
    throw Exception() << "Reverse Ordering not implemented in " 
		      << "SudakovFormFactor::generateNextSpaceBranching"
		      << Exception::runerror;
  // normal ordering
  // All the variables needed
  Energy2 t,t0,tmax,tmin;
  Energy kinCutoff;
  bool glueEmits = (ids[0]==ParticleID::g);
  
  // Different order, incoming parton is id =  1, outgoing are id=0,2
  tcPDPtr parton0 = getParticleData(ids[0]);
  tcPDPtr parton1 = getParticleData(ids[1]);
  Energy m0 = parton0->mass();
  Energy m1 = parton1->mass();
  tmax = sqr(startingQ);
  t = tmax;
  // Initialize the variables
  if(glueEmits) initialize(t0,tmin,tmax,kinCutoff,m1);
  else initialize(t0,tmin,tmax,kinCutoff,m0);
  if(_q < 0.) return _q;
  
  // Now do the veto algorithm
  do 
    { 
      guessSpaceLike(_z,z0,z1,t,t0,kinCutoff,x);
      if(z0 > z1 || tVeto(t,tmin)) break;
    } 
  while(_z > z1 || 
	SplittingFnVeto(_z,t/sqr(_z),ids) || 
	alphaSVeto(sqr(1.-_z)*t) || 
	PDFVeto(_z,t,x,parton0,parton1));
  if(t > 0 && z0 < z1) _q = sqrt(t);
  else _q = -1.;

  _phi = 2.*pi*UseRandom::rnd();
  
  return _q;
}

void SudakovFormFactor::initialize(const IdList & ids)
{
  _ids=ids;
  _masses.clear();
  _masssquared.clear();
  unsigned int ix;
  for(ix=0;ix<_ids.size();++ix)
    _masses.push_back(getParticleData(_ids[ix])->constituentMass());
  Energy kinCutoff=
    _variables->kinematicCutOff(kinScale(),
				*std::max_element(_masses.begin(),_masses.end()));
  for(ix=0;ix<_masses.size();++ix)
    {
      _masses[ix]=max(kinCutoff,_masses[ix]);
      _masssquared.push_back(_masses[ix]);
    }
}
