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

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "SudakovFormFactor.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

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

bool SudakovFormFactor::guessTimeLike(Energy2 &t,Energy2 tmin) const
{
  Energy2 told = t;
  // calculate limits on z and if lower>upper return
  if(!computeTimeLikeLimits(t)) return false;
  // guess values of t and z
  t = guesst(told,0); 
  _z = guessz(); 
  // actual values for z-limits
  if(!computeTimeLikeLimits(t)) return false;
  if(t<tmin)
    {
      t=-1.0*GeV;
      return false;
    }
  else
    return true; 
} 

bool SudakovFormFactor::guessSpaceLike(Energy2 &t, Energy2 tmin, const double x) const
{
  Energy2 told = t;
  // calculate limits on z if lower>upper return
  if(!computeSpaceLikeLimits(t,x)) return false;
  // guess values of t and z
  t = guesst(told,1); 
  _z = guessz(); 
  // actual values for z-limits
  if(!computeSpaceLikeLimits(t,x)) return false;
  if(t<tmin)
    {
      t=-1.0*GeV;
      return false;
    }
  else
    return true; 
} 

bool SudakovFormFactor::PSVeto(const Energy2 t) {
  // still inside PS, return true if outside
  // check vs overestimated limits
  if(_z < _zlimits.first || _z > _zlimits.second) return true;
  // compute the pt
  Energy2 pt2=sqr(_z*(1.-_z))*t-_masssquared[1]*(1.-_z)-_masssquared[2]*_z;
  if(_ids[0]!=ParticleID::g) pt2+=_z*(1.-_z)*_masssquared[0];
  // if pt2<0 veto
  if(pt2<0.) return true;
  // otherwise calculate pt and return
  _pt=sqrt(pt2);
  return false;
}
 
bool SudakovFormFactor::PDFVeto(const double z, const Energy2 t, 
				const double x,
				const tcPDPtr parton0, 
				const tcPDPtr parton1) const {
  assert(_variables->currentPDF());
  // remember: pdf's q is cut in pdf class.  should probably be done here! 
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
						    const IdList &ids)
{
  // First reset the internal kinematics variables that can
  // have been eventually set in the previous call to the method.
  _q = Energy();
  _z = 0.0;
  _phi = 0.0; 
  // perform initialization
  Energy2 tmax(sqr(startingScale)),tmin;
  initialize(ids,tmin);
  // check max > min
  if(tmax<=tmin)
    {
      _q=-1.;
      return _q;
    }
  // calculate next value of t using veto algorithm
  Energy2 t(tmax);
  do  
    if(!guessTimeLike(t,tmin)) break;
  while(PSVeto(t) || SplittingFnVeto(_z,t,ids) || 
	alphaSVeto(sqr(_z*(1.-_z))*t));
  if(t > 0) _q = sqrt(t);
  else _q = -1.;
  _phi = 2.*pi*UseRandom::rnd();
  return _q;
}

Energy SudakovFormFactor::generateNextSpaceBranching(const Energy startingQ,
			                             const IdList &ids,
						     double x) {
  // First reset the internal kinematics variables that can
  // have been eventually set in the previous call to the method.
  _q = Energy();
  _z = 0.0;
  _phi = 0.0;
  // perform the initialization
  Energy2 tmax(sqr(startingQ)),tmin;
  initialize(ids,tmin);
  // check max > min
  if(tmax<=tmin)
    {
      _q=-1.;
      return _q;
    }
  // extract the partons which are needed for the PDF veto
  // Different order, incoming parton is id =  1, outgoing are id=0,2
  tcPDPtr parton0 = getParticleData(ids[0]);
  tcPDPtr parton1 = getParticleData(ids[1]);
  // calculate next value of t using veto algorithm
  Energy2 t(tmax),pt2(0.);
  do
    {
      if(!guessSpaceLike(t,tmin,x)) break;
      pt2=sqr(1.-_z)*t-_z*sqr(_kinCutoff);
    }
  while(_z > _zlimits.second || 
	SplittingFnVeto(_z,t/sqr(_z),ids) || 
	alphaSVeto(sqr(1.-_z)*t) || 
	PDFVeto(_z,t,x,parton0,parton1)||pt2<0);
  if(t > 0 && _zlimits.first < _zlimits.second) 
    _q = sqrt(t);
  else
    {
      _q = -1.;
      return _q;
    }
  _phi = 2.*pi*UseRandom::rnd();
  _pt=sqrt(pt2);
  return _q;
}

void SudakovFormFactor::initialize(const IdList & ids, Energy2 & tmin)
{
  _ids=ids;
  _masses.clear();
  _masssquared.clear();
  tmin=0.;
  unsigned int ix;
  for(ix=0;ix<_ids.size();++ix)
    _masses.push_back(getParticleData(_ids[ix])->mass());
  _kinCutoff=
    _variables->kinematicCutOff(kinScale(),
				*std::max_element(_masses.begin(),_masses.end()));
  for(ix=0;ix<_masses.size();++ix)
    {
      _masses[ix]=max(_kinCutoff,_masses[ix]);
      _masssquared.push_back(sqr(_masses[ix]));
      if(ix>0) tmin=max(_masssquared[ix],tmin);
    }
}

Energy SudakovFormFactor::generateNextDecayBranching(const Energy startingScale,
						     const Energy stoppingScale,
						     const Energy minmass,
						     const IdList &ids)
{
  // First reset the internal kinematics variables that can
  // have been eventually set in the previous call to this method.
  _q = ShowerVariables::HUGEMASS;
  _z = 0.0;
  _phi = 0.0; 
  // perform initialisation
  Energy2 tmax(sqr(stoppingScale)),tmin;
  initialize(ids,tmin);
  tmin=sqr(startingScale);
  // check some branching possible
  if(tmax<=tmin)
    {
      _q=-1.;
      return _q;
    }
  // perform the evolution
  Energy2 t(tmin);
  do 
    if(!guessDecay(t,tmax,minmass)) break;
  while(SplittingFnVeto(_z,t/sqr(_z),ids)|| 
	alphaSVeto(sqr(1.-_z)*t)||
	sqr(1.-_z)*(t-_masssquared[0])<_z*sqr(_kinCutoff)||
	t*(1.-_z)>_masssquared[0]-sqr(minmass));
  if(t > 0)
    {
      _q = sqrt(t);
      _pt = sqrt(sqr(1.-_z)*(t-_masssquared[0])-_z*sqr(_kinCutoff));
    }
  else _q = -1.;
  _phi = 2.*pi*UseRandom::rnd();
  return _q;
}

bool SudakovFormFactor::guessDecay(Energy2 &t,Energy2 tmax, Energy minmass) const
{
  // previous scale
  Energy2 told = t;
  // overestimated limits on z
  _zlimits.first  = sqr(minmass/_masses[0]);
  _zlimits.second = 1.-_kinCutoff/sqrt(tmax-_masssquared[0])
    +0.5*sqr(_kinCutoff)/(tmax-_masssquared[0]);
  // guess values of t and z
  t = guesst(told,2); 
  _z = guessz(); 
  // actual values for z-limits
  _zlimits.first  = 0.;
  _zlimits.second = 1.-_kinCutoff/sqrt(t-_masssquared[0])
    +0.5*sqr(_kinCutoff)/(t-_masssquared[0]);
  if(t>tmax)
    {
      t=-1.0*GeV;
      return false;
    }
  else
    return true; 
} 

bool SudakovFormFactor::computeTimeLikeLimits(Energy2 & t) const
{
  // special case for gluon radiating
  if(_ids[0]==ParticleID::g)
    {
      // no emission possible
      if(t<16.*_masssquared[1])
	{
	  t=-1.;
	  return false;
	}
      // overestimate of the limits
      _zlimits.first  = 0.5*(1.-sqrt(1.-4.*sqrt(_masssquared[1]/t)));
      _zlimits.second = 1.-_zlimits.first;
    }
  // special case for radiated particle is gluon 
  else if(_ids[2]==ParticleID::g)
    {
      _zlimits.first  = 0.5*sqrt(_masssquared[1]/t);
      _zlimits.second = 1.-0.5*sqrt(_masssquared[2]/t);
    }
  else if(_ids[1]==ParticleID::g)
    {
      _zlimits.second  = 0.5*sqrt(_masssquared[2]/t);
      _zlimits.first   = 1.-0.5*sqrt(_masssquared[1]/t);
    }
  else
    {throw Exception() << "SudakovFormFactor::computeTimeLikeLimits() " 
			<< "general case not implemented " << Exception::runerror;}
  return true;
}

bool SudakovFormFactor::computeSpaceLikeLimits(Energy2 & t, double x) const
{
  // compute the limits
  _zlimits.first = x;
  double yy = 1.+0.5*sqr(_kinCutoff)/t;
  _zlimits.second = yy - sqrt(sqr(yy)-1.); 
  // return false if lower>upper
  if(_zlimits.second<_zlimits.first)
    {
      t=-1.*GeV;
      return false;
    }
  else
    return true;
}
