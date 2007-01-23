// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ForcedSplitting class.
//

#include "ForcedSplitting.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include <ThePEG/Interface/Reference.h>
#include <ThePEG/Interface/Parameter.h>
#include <ThePEG/PDT/EnumParticles.h>
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Remnant.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "ThePEG/Utilities/Timer.h"
#include <cassert>

using namespace Herwig;

void ForcedSplitting::persistentOutput(PersistentOStream & os) const {
  os << _kinCutoff << _range << _qspac << _zbin << _ybin << _nbinmax << _alpha;
}

void ForcedSplitting::persistentInput(PersistentIStream & is, int) {
  is >> _kinCutoff >> _range >> _qspac >> _zbin >> _ybin >> _nbinmax >> _alpha;
}

ClassDescription<ForcedSplitting> ForcedSplitting::initForcedSplitting;
// Definition of the static class description member.

void ForcedSplitting::Init() {

  static ClassDocumentation<ForcedSplitting> documentation
    ("This class is responsible for correctly tying the parton shower to "
     "the remaining flavours in the hadron and producing the correct remnant");

  static Parameter<ForcedSplitting,Energy> interfaceKinCutoff
    ("KinCutoff",
     "Parameter kinCutoff used to constrain qtilde",
     &ForcedSplitting::_kinCutoff, GeV, 0.75*GeV, 0.5*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<ForcedSplitting,double> interfaceEmissionRange
    ("EmissionRange",
     "Factor above the minimum possible value in which the forced splitting is allowed.",
     &ForcedSplitting::_range, 1.1, 1.0, 10.0,
     false, false, Interface::limited);

  static Parameter<ForcedSplitting,double> interfaceZBinSize
    ("ZBinSize",
     "The size of the vbins in z for the interpolation of the splitting function.",
     &ForcedSplitting::_zbin, 0.05, 0.001, 0.1,
     false, false, Interface::limited);

  static Parameter<ForcedSplitting,int> interfaceMaxBin
    ("MaxBin",
     "Maximum number of z bins",
     &ForcedSplitting::_nbinmax, 100, 10, 1000,
     false, false, Interface::limited);

  static Reference<ForcedSplitting,ShowerAlpha> interfaceAlphaS
    ("AlphaS",
     "Pointer to object to calculate the strong coupling",
     &ForcedSplitting::_alpha, false, false, true, false, false);


  static Parameter<ForcedSplitting,Energy> interfaceQSpac
    ("QSpac",
     "The starting scale for the evolution in the forced splitting",
     &ForcedSplitting::_qspac, GeV, 2.5*GeV, 0.5*GeV, 10.0*GeV,
     false, false, Interface::limited);

}

tPVector ForcedSplitting::split(const tPVector & tagged, tStepPtr pstep) {
  Timer<5000> timer("ForcedSplitting::split");
  // extract the remnants
  tPVector rems;
  tPVector output;
  for(unsigned int ix=0;ix<tagged.size();++ix) {
    if(tagged[ix]->id()!=ExtraParticleID::Remnant) output.push_back(tagged[ix]);
    else rems.push_back(tagged[ix]);
  }
  // if none return
  if(rems.size()==0) return output;
  // must be two
  if(rems.size()!=2) throw Exception() << "Must be either no Remnants or two in "
				       << "ForcedSplitting::split()"
				       << Exception::eventerror;
  try {
    // Find beam particles
    PPair beam = generator()->currentEventHandler ()->currentCollision()->incoming();
    tPPtr rem[2];
    tPPtr part1,part2;
    PVector::const_iterator it= beam.first->children().begin();
    // remnant from first beam
    for(; it != beam.first->children().end();it++) {
      if((*it)->children().size()==0) rem[0] = *it;
      else part1 = *it;
    }
    // remnant from second beam
    for(it = beam.second->children().begin(); it !=beam.second->children().end();
	it++) {
      if((*it)->children().size()==0) rem[1] = *it;
      else part2 = *it;
    }
    // momentum fractions of parton going into the shower
    double x1(part1->momentum().rho()/beam.first->momentum().rho());
    double x2(part2->momentum().rho()/beam.second->momentum().rho());
    // split the remnants
    if(rem[0]) {
      _beam=dynamic_ptr_cast<Ptr<BeamParticleData>::const_pointer>
	(beam.first->dataPtr());
      split(rem[0],part1,pstep,x1);
    }
    if(rem[1]) {
      _beam=dynamic_ptr_cast<Ptr<BeamParticleData>::const_pointer>
	(beam.second->dataPtr());
      split(rem[1],part2,pstep,x2);
    }
    // get the masses of the remnants
    Energy mrem[2];
    Lorentz5Momentum ptotal,pnew[2];
    for(unsigned int ix=0;ix<2;++ix) {
      pnew[ix]=Lorentz5Momentum();
      for(unsigned int iy=0;iy<rem[ix]->children().size();++iy) {
	pnew[ix]+=rem[ix]->children()[iy]->momentum();
      }
      mrem[ix]=sqrt(pnew[ix].m2());
    }
    // now find the remnant remnant cmf frame
    Lorentz5Momentum prem[2]={rem[0]->momentum(),rem[1]->momentum()};
    ptotal=prem[0]+prem[1];
    ptotal.rescaleMass();
    // boost momenta to this frame
    if(ptotal.m()<0) throw Exception() << "Space-Like Remnant in " 
				       << "ForcedSplitting::handle() " 
				       << Exception::eventerror;
    Hep3Vector boostv(-ptotal.boostVector());
    ptotal.boost(boostv);
    for(unsigned int ix=0;ix<2;++ix) {
      prem[ix].boost(boostv);
      // set the masses and energies,
      prem[ix].setMass(mrem[ix]);
      prem[ix].setE(0.5/ptotal.m()*(sqr(ptotal.m())+sqr(mrem[ix])-sqr(mrem[1-ix])));
      prem[ix].rescaleRho();
      // boost back to the lab
      prem[ix].boost(-boostv);
      // set the momenta of the remnants
      rem[ix]->set5Momentum(prem[ix]);
    }
    // boost the decay products
    for(unsigned int ix=0;ix<2;++ix) {
      // factors for Lorentz transform
      Energy ea(pnew[ix].e()),eb(prem[ix].e()),mr(prem[ix].mass());
      Energy2 mr2(sqr(mr));
      long double beta1 = -pnew[ix].pz()/pnew[ix].e();
      long double beta2 =  prem[ix].pz()/prem[ix].e();
      long double gamma = (ea*eb)/mr2;
      long double sum   = beta1+beta2;
      long double prod  = 1.+beta1*beta2;
      // small approx for accuracy
      if(abs(beta1)<1e-5||abs(beta2)<1e-5) {
	prod = 0.5*mr2*(sqr(1./ea)+sqr(1./eb))
	  +0.125*sqr(mr2)*sqr(1./ea+1./eb)*sqr(1./ea-1./eb);
	sum = 0.5*mr2*(1./ea+1./eb)*(1./ea-1./eb)*
	  (1.+0.25*mr2*(sqr(1./ea)+sqr(1./eb)));
	if(beta1>0) sum*=-1.;
      }
      // boost the children
      for(unsigned int iy=0;iy<rem[ix]->children().size();++iy) {
	Energy pz(rem[ix]->children()[iy]->momentum().pz());
	Energy ee(rem[ix]->children()[iy]->momentum().e());
	Energy mm(rem[ix]->children()[iy]->momentum().mass());
	Lorentz5Momentum pold=rem[ix]->children()[iy]->momentum();
	Lorentz5Momentum ptemp(pold.px(),pold.py(),gamma*(ee*sum+pz*prod),
			       gamma*(pz*sum+ee*prod),mm);
	rem[ix]->children()[iy]->set5Momentum(ptemp);
      }
    }
  }
  catch(std::exception & e) 
    {throw Exception() << "Caught exception\n"
		       << e.what() 
		       <<  "\nin ForcedSplitting::handle() "
		       << Exception::eventerror;}
  // insert the products of the splitting in the vector
  for(unsigned int ix=0;ix<rems.size();++ix) {
    for(unsigned int iy=0;iy<rems[ix]->children().size();++iy) {
      output.push_back(rems[ix]->children()[iy]);
    }
  }
  return output;
}

/****
 * Now we need to force the final (up to two) splittings so that we are left
 * with only the valence quarks that make up the incoming hadron. If we have
 * terminated the shower on a sea quark, then we need two splittings, one
 * to a gluon and one to a valence quark. If we are on a gluon we just go
 * to a valence quark. We must also choose qtilda and z for each splitting so
 * that the kinematics can be reconstructed properly. This is no longer 
 * sampled according to the splitting functions, as they no longer have space
 * in the virtuality (since the shower has terminated). Instead we use a new
 * distribution.
 * NOTE: temporarily chosen linearly in z and logarithmically in qtilda, this
 * may be changed later.
 ****/
void ForcedSplitting::split(const tPPtr rem,const tPPtr part, 
			    const tStepPtr step,const double xin) {
  long hadronId = rem->parents()[0]->id();
  // return if not from a hadron
  if(abs(hadronId)<99) return;
  int maxIdx = 3;
  int idx = -1;
  long lg(ParticleID::g),currentPart(part->id());
  Lorentz5Momentum usedMomentum=Lorentz5Momentum(),lastp = part->momentum();
  PPtr lastColour(part),newPart;
  ColinePtr x1, x2;
  long quarks[3];
  quarks[0] = hadronId % 10;
  quarks[1] = (hadronId/10)%10;
  quarks[2] = (hadronId/100)%10;
  
  // NOTE TODO: Must make sure that the sign is correct for the meson quarks
  if(quarks[2] == 0) maxIdx = 2; // we have a meson
  
  // initial x and scale
  Energy oldx = xin; 
  Energy oldQ=_qspac;

  // Look first at sea quarks, these must go to a gluon, we then handle
  // the gluons in the next step
  if(currentPart != quarks[0] && currentPart != quarks[1] && 
     currentPart != quarks[2] && currentPart != ParticleID::g) {
    // Create the new parton with its momentum and parent/child relationship set
    newPart = forceSplit(rem, -currentPart, oldQ, oldx, lastp, 
			 usedMomentum,1, step);
    // Set the proper colour connections
    x1 = new_ptr(ColourLine());
    if(lastColour->colourLine()) x1->addAntiColoured(newPart);
    else if(lastColour->antiColourLine()) x1->addColoured(newPart);
    currentPart = lg;
  }
  // We now handle the gluons, either it is where the shower terminated or
  // it has been created by splitting a sea quark
  // gluon
  if(currentPart == ParticleID::g) { 
    // Create new particles, splitting is q->g q
    // First choose which q
    idx = UseRandom::irnd(maxIdx);
    Lorentz5Momentum s = rem->momentum()-usedMomentum;
    // Generate the next parton, with s momentum remaining in the remnant.
    oldQ=_qspac;
    newPart = forceSplit(rem, quarks[idx], oldQ, oldx, lastp, usedMomentum,2,step);
    // Several colour connection cases...
    if(x1) {
      bool npC = newPart->hasColour();
      if(npC) {
	x2 = lastColour->antiColourLine();
	if(!x1->antiColoured().empty())   x1->addColoured(newPart);
	else if(!x1->coloured().empty())  x2->addColoured(newPart);
      } 
      else {
	ColinePtr x3 = lastColour->colourLine();
	if(!x1->coloured().empty()) x1->addAntiColoured(newPart);
	else if(!x1->antiColoured().empty()) x3->addAntiColoured(newPart);
	x2=x1;
	x1=x3;
      }
    } 
    else {
      x1 = lastColour->colourLine();
      x2 = lastColour->antiColourLine();
      if(newPart->hasColour()) x2->addColoured(newPart);
      else if(newPart->hasAntiColour()) x1->addAntiColoured(newPart);
    }
    currentPart = quarks[idx];
  } 
  // find the extracted quark if not known
  if(idx<0) {
    unsigned int ix=0;
    do {
      if(quarks[ix]==currentPart){idx=ix;}
      ++ix;
    }
    while(idx<0&&ix<3);
  }
  // Lastly, do the final split into the (di)quark and a parton
  newPart = finalSplit(rem,maxIdx,quarks,idx,usedMomentum, step);
  // Set colour connections, case 1, no other forced splittings
  if(!x1 || !x2) {
    if(rem->colourLine()) rem->colourLine()->addColoured(newPart);
    else if(rem->antiColourLine()) 
      rem->antiColourLine()->addAntiColoured(newPart);
  } 
  else {
    if(getParticleData(currentPart)->hasColour()) x1->addAntiColoured(newPart);
    else x2->addColoured(newPart);
  }
}

// This creates the parton to split and sets it momentum and parent/child
// relationships
PPtr ForcedSplitting::forceSplit(const tPPtr rem, long child, Energy &oldQ, 
				 double &oldx, Lorentz5Momentum &pf, 
				 Lorentz5Momentum &p,
				 const unsigned int iopt,
				 const tStepPtr step) {
  Lorentz5Momentum beam = rem->parents()[0]->momentum();
  PPtr parton = new_ptr(Particle(getParticleData(child)));
  Lorentz5Momentum partonp = emit(beam,oldQ,oldx,parton,pf,iopt);
  p += partonp;
  parton->set5Momentum(partonp);
  step->addDecayProduct(rem,parton);
  return parton;
}

// This forces the final output of the remnant ((di)quark) and sets the
// momentum and parent/child relationships
PPtr ForcedSplitting::finalSplit(const tPPtr rem, int maxIdx, 
				 long quarks[3], int idx, 
				 Lorentz5Momentum usedMomentum, 
				 const tStepPtr step) {
  tPPtr hadron;
  if(rem->parents().size() == 1) hadron = rem->parents()[0];
  // First decide what the remnant is
  long remId(0);
  int sign, spin;// Meson hadronic state
  if(maxIdx == 2) { 
    remId = quarks[(idx+1)%2];
  }
  // Baryonic hadron 
  else { 
    // Get the other 2 elements of the array
    // use modulus to simplify things. idx is the array entry for the
    // parton which eventually leads to the hard process, the other two
    // elements of the array constitute the remnant.
    long id1 = quarks[(idx+1)%3];
    long id2 = quarks[(idx+2)%3];
    if (abs(id1) > abs(id2)) swap(id1, id2);
    sign = (id1 < 0) ? -1 : 1; // Needed for the spin 0/1 part
    remId = id2*1000+id1*100;
    
    // Now decide if we have spin 0 diquark or spin 1 diquark
    if(id1 == id2) spin = 3; // spin 1
    else spin = 1; // otherwise spin 0
    remId += sign*spin;
  }
   
  // Create the remnant and set its momentum, also reset all of the decay 
  // products from the hadron
  PPtr remnant = new_ptr(Particle(getParticleData(remId)));
  Lorentz5Momentum prem(rem->momentum()-usedMomentum);
  prem.setMass(getParticleData(remId)->constituentMass());
  prem.rescaleEnergy();
  remnant->set5Momentum(prem);
  // Add the remnant to the step, this will be changed again if the
  // shower is vetoed. Set the colour connections as well
  step->addDecayProduct(rem,remnant);
  return remnant;
}

// This defines the momentum for an emitted parton, currently no pt is
// given to the produced partons, z is generated uniformly.
Lorentz5Momentum ForcedSplitting::emit(const Lorentz5Momentum &par,
				       Energy &lastQ, double &lastx, 
				       PPtr parton,
				       Lorentz5Momentum &pf,
				       const unsigned int iopt) {
  assert(iopt==1||iopt==2);
  // the last scale is minimum of last value and upper limit
  Energy minQ=_range*_kinCutoff*sqrt(lastx)/(1-lastx);
  if(minQ>lastQ) lastQ=minQ;
  // generate the new value of qtilde
  // weighted towards the lower value: dP/dQ = 1/Q -> Q(R) =
  // Q0 (Qmax/Q0)^R
  Energy q;
  double zmin,zmax,yy;
  do {
    q = minQ*pow(lastQ/minQ,UseRandom::rnd());
    zmin = lastx;
    yy   = 1.+0.5*sqr(_kinCutoff/q);
    zmax = yy - sqrt(sqr(yy)-1.);
  }
  while(zmax<zmin);
  // now generate z as in FORTRAN HERWIG
  // use y = ln(z/(1-z)) as integration variable
  double ymin=log(zmin/(1.-zmin));
  double ymax=log(zmax/(1.-zmax));
  double dely=ymax-ymin;
  unsigned int nz=std::min(int(_ybin*dely+1),_nbinmax);
  dely/=nz;
  yy=ymin+0.5*dely;
  double psum(0.);
  tcPDPtr gluon=getParticleData(ParticleID::g);
  vector<double> prob;
  for(unsigned int iz=0;iz<nz;++iz) {
    double ez=exp(yy);
    double wr=1.+ez;
    double zr=wr/ez;
    double wz=1./wr;
    double zz=wz*ez;
    double az=wz*zz*_alpha->value(max(wz*q,_kinCutoff));
    // g -> q qbar
    if(iopt==1) {
      // calculate splitting function
      double pdf=_beam->pdf()->xfx(_beam,gluon,sqr(q),lastx*zr);
      psum+=pdf*az*0.5*(sqr(zz)+sqr(wz));
      prob.push_back(psum);
    }
    // q -> q g
    else {
      // calculate splitting function
      double pdf=_beam->pdf()->xfx(_beam,parton->dataPtr(),sqr(q),lastx*zr);
      psum+=pdf*az*4./3.*(1.+sqr(wz))*zr;
      prob.push_back(psum);
    }
    yy+=dely;
  }
  // choose z
  double pval=psum*UseRandom::rnd();
  unsigned int iz=0;
  for(;iz<prob.size();++iz) {
    if(prob[iz]>pval) break;
  }
  if(iz==prob.size()) --iz;
  double ey=exp(ymin+dely*(float(iz+1)-UseRandom::rnd()));
  double z=ey/(1.+ey);
  Energy pt2=sqr((1.-z)*q)- z*sqr(_kinCutoff);
  Energy2 emittedm2 = sqr(parton->dataPtr()->constituentMass());
  // Now boost pcm and pf to z only frame
  Lorentz5Momentum p       = Lorentz5Momentum(0.0,  par.vect());
  Lorentz5Momentum n       = Lorentz5Momentum(0.0, -par.vect());
  // generate phi and compute pt of branching
  double phi = 2.*pi*UseRandom::rnd();
  Energy pt=sqrt(pt2);
  Lorentz5Momentum qt   = LorentzMomentum(pt*cos(phi), pt*sin(phi), 0.0, 0.0);
  // compute alpha for previous particle
  Energy2 p_dot_n  = p*n;
  double lastalpha = pf*n/p_dot_n;
  Lorentz5Momentum qtout=qt;
  Energy2 qtout2=-qt*qt;
  double alphaout=(1.-z)/z*lastalpha;
  double betaout=0.5*(emittedm2+qtout2)/alphaout/p_dot_n;
  Lorentz5Momentum k=alphaout*p+betaout*n+qtout;
  k.rescaleMass();
  pf+=k;
  lastQ=q;
  lastx/=z;
  return k;
}

