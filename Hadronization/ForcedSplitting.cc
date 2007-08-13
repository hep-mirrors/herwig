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
#include "CheckId.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "ThePEG/Utilities/Timer.h"
#include <cassert>

using namespace Herwig;

void ForcedSplitting::persistentOutput(PersistentOStream & os) const {
  os << ounit(_kinCutoff,GeV) << _range << ounit(_qspac,GeV) 
     << _zbin << _ybin << _nbinmax << _alpha;
}

void ForcedSplitting::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_kinCutoff,GeV) >> _range >> iunit(_qspac,GeV) 
     >> _zbin >> _ybin >> _nbinmax >> _alpha;
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

void ForcedSplitting::split(PVector & currentlist) {
  Timer<5000> timer("ForcedSplitting::split");
  // extract the remnants
  tPVector rems;
  PVector output;
  for(unsigned int ix=0;ix<currentlist.size();++ix) {
    if(currentlist[ix]->id()!=ExtraParticleID::Remnant) 
      output.push_back(currentlist[ix]);
    else 
      rems.push_back(currentlist[ix]);
  }
  // if none return
  if(rems.size()==0) return;
  // must be two 
  /// \todo what about ep collisions?
  if(rems.size()!=2) 
    throw Exception() << "Must be either no Remnants or two in "
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

    assert(rem[0] && rem[1]);

    // momentum fractions of parton going into the shower
    double x1(part1->momentum().rho()/beam.first->momentum().rho());
    double x2(part2->momentum().rho()/beam.second->momentum().rho());
    // split the remnants
    _beam = dynamic_ptr_cast<Ptr<BeamParticleData>::const_pointer>
      (beam.first->dataPtr());
    split(rem[0],part1,x1);
    
    
    _beam = dynamic_ptr_cast<Ptr<BeamParticleData>::const_pointer>
      (beam.second->dataPtr());
    split(rem[1],part2,x2);
      
    // get the masses of the remnants
    Energy mrem[2];
    Lorentz5Momentum pnew[2];
    for(unsigned int ix=0; ix < 2; ++ix) {
      pnew[ix]=Lorentz5Momentum();
      for(unsigned int iy=0; iy < rem[ix]->children().size(); ++iy) {
	pnew[ix] += rem[ix]->children()[iy]->momentum();
      }
      mrem[ix] = sqrt(pnew[ix].m2());
    }
    // now find the remnant remnant cmf frame
    Lorentz5Momentum prem[2]={rem[0]->momentum(),
			      rem[1]->momentum()};
    Lorentz5Momentum ptotal = prem[0] + prem[1];
    ptotal.rescaleMass();
    // boost momenta to this frame
    if(ptotal.m() < Energy()) throw Exception() << "Space-Like Remnant in " 
						<< "ForcedSplitting::handle() " 
						<< Exception::eventerror;
    Boost boostv(-ptotal.boostVector());
    ptotal.boost(boostv);
    for(unsigned int ix=0; ix<2; ++ix) {
      prem[ix].boost(boostv);
      // set the masses and energies,
      prem[ix].setMass(mrem[ix]);
      prem[ix].setE(0.5/ptotal.m() * (sqr(ptotal.m()) 
				      + sqr(mrem[ix])
				      - sqr(mrem[1-ix]))
		    );
      prem[ix].rescaleRho();
      // boost back to the lab
      prem[ix].boost(-boostv);
      // set the momenta of the remnants
      rem[ix]->set5Momentum(prem[ix]);
    }
    // boost the decay products
    for(unsigned int ix=0;ix<2;++ix) { 
      Boost btorest(-pnew[ix].boostVector()); 
      Boost bfmrest( prem[ix].boostVector()); 
      for(unsigned int iy=0;iy<rem[ix]->children().size();++iy) { 
	Lorentz5Momentum ptemp = rem[ix]->children()[iy]->momentum(); 
	ptemp.boost(btorest); 
	ptemp.boost(bfmrest); 
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
  swap(output,currentlist);
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
			    const double xin) {
  long hadronId = rem->parents()[0]->id();
  // return if not from a hadron
  /// \todo where else would it come from? throw here?
  if(abs(hadronId)<99) return;

  tcPDPtr currentParton = part->dataPtr();
  Lorentz5Momentum usedMomentum;
  Lorentz5Momentum lastp = part->momentum();
  PPtr lastColouredParton(part);

  long quarks[3];
  quarks[0] = (hadronId/10)   % 10;
  quarks[1] = (hadronId/100)  % 10;
  quarks[2] = (hadronId/1000) % 10;
  
  /// \todo Must make sure that the sign is correct for the meson quarks

  const int maxIdx = (quarks[2] == 0) ? 2 : 3; // we have a meson / baryon
  
  // initial x and scale
  double oldx = xin; 
  Energy oldQ=_qspac;

  // Look first at sea quarks, these must go to a gluon, we then handle
  // the gluons in the next step
  ColinePtr x1,x2;

  if(currentParton->id() != quarks[0] && currentParton->id() != quarks[1] && 
     currentParton->id() != quarks[2] && currentParton->id() != ParticleID::g) {
    // Create the new parton with its momentum and parent/child relationship set
    PPtr newParton = forceSplit(rem, currentParton->CC(), oldQ, oldx, lastp, 
				usedMomentum,1);
    // Set the proper colour connections
    x1 = new_ptr(ColourLine());
    if(lastColouredParton->colourLine()) x1->addAntiColoured(newParton);
    else if(lastColouredParton->antiColourLine()) x1->addColoured(newParton);
    currentParton = getParticleData(ParticleID::g);
  }
  // We now handle the gluons, either it is where the shower terminated or
  // it has been created by splitting a sea quark
  int idx = -1;
  if(currentParton->id() == ParticleID::g) { 
    // Create new particles, splitting is q->g q
    // First choose which q
    idx = UseRandom::irnd(maxIdx);
    Lorentz5Momentum s = rem->momentum() - usedMomentum;
    // Generate the next parton, with s momentum remaining in the remnant.
    oldQ = _qspac;
    PPtr newParton = forceSplit(rem,getParticleData(quarks[idx]), oldQ, oldx, 
				lastp, usedMomentum,2);
    // Several colour connection cases...
    if(x1) {
      bool npC = newParton->hasColour();
      if(npC) {
	x2 = lastColouredParton->antiColourLine();
	if(!x1->antiColoured().empty())   x1->addColoured(newParton);
	else if(!x1->coloured().empty())  x2->addColoured(newParton);
      } 
      else {
	ColinePtr x3 = lastColouredParton->colourLine();
	if(!x1->coloured().empty()) x1->addAntiColoured(newParton);
	else if(!x1->antiColoured().empty()) x3->addAntiColoured(newParton);
	x2=x1;
	x1=x3;
      }
    } 
    else {
      x1 = lastColouredParton->colourLine();
      x2 = lastColouredParton->antiColourLine();
      if(newParton->hasColour()) x2->addColoured(newParton);
      else if(newParton->hasAntiColour()) x1->addAntiColoured(newParton);
    }
    currentParton = getParticleData(quarks[idx]);
  } 
  // find the extracted quark if not known
  if(idx == -1) {
    // must have extracted valence quark
    if      (currentParton->id() ==  quarks[0]) idx = 0;
    else if (currentParton->id() ==  quarks[1]) idx = 1;
    else if (currentParton->id() ==  quarks[2]) idx = 2;
    else assert(false && "Should have got valence quark in ForcedSplitting::split()");
  }
  // Lastly, do the final split into the (di)quark and a parton
  PPtr newParton = finalSplit(rem,maxIdx,quarks,idx,usedMomentum);
  // Set colour connections, case 1, no other forced splittings
  if(!x1 || !x2) {
    if(rem->colourLine()) rem->colourLine()->addColoured(newParton);
    else if(rem->antiColourLine()) 
      rem->antiColourLine()->addAntiColoured(newParton);
  } 
  else {
    if(currentParton->hasColour()) x1->addAntiColoured(newParton);
    else x2->addColoured(newParton);
  }
}

// This creates the parton to split and sets it momentum and parent/child
// relationships
PPtr ForcedSplitting::forceSplit(const tPPtr rem, tcPDPtr child, Energy &oldQ, 
				 double &oldx, Lorentz5Momentum &pf, 
				 Lorentz5Momentum &p,
				 const unsigned int iopt) {
  Lorentz5Momentum beam = rem->parents()[0]->momentum();
  PPtr parton = child->produceParticle();
  Lorentz5Momentum partonp = emit(beam,oldQ,oldx,parton,pf,iopt);
  parton->set5Momentum(partonp);
  p += partonp;
  //  step->addDecayProduct(rem,parton);
  rem->addChild(parton);
  return parton;
}

// This forces the final output of the remnant ((di)quark) and sets the
// momentum and parent/child relationships
PPtr ForcedSplitting::finalSplit(const tPPtr rem, int maxIdx, 
				 long quarks[3], int idx, 
				 Lorentz5Momentum usedMomentum) {
  // First decide what the remnant is
  PPtr remnant = PPtr ();

  if(maxIdx == 2) { 
    // Meson
    assert(idx == 0 || idx == 1);
    long remId = quarks[(idx+1)%2];
    remnant = getParticle(remId);
  }
  // Baryon
  else { 
    // Get the other 2 elements of the array
    // use modulus to simplify things. idx is the array entry for the
    // parton which eventually leads to the hard process, the other two
    // elements of the array constitute the remnant.
    long id1 = quarks[(idx+1)%3];
    long id2 = quarks[(idx+2)%3];

    tcPDPtr remData = getParticleData(CheckId::makeDiquarkID(id1,id2));
    remnant = remData->produceParticle();
  }
   
  // Create the remnant and set its momentum, also reset all of the decay 
  // products from the hadron

  Lorentz5Momentum prem(rem->momentum() - usedMomentum);
  prem.setMass(remnant->dataPtr()->constituentMass());
  prem.rescaleEnergy();

  remnant->set5Momentum(prem);
  // Add the remnant to the step, this will be changed again if the
  // shower is vetoed 
  /// \todo Is this comment still true?. 
  // Set the colour connections as well
  //  step->addDecayProduct(rem,remnant);
  rem->addChild(remnant);
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
    double az=wz*zz*_alpha->value(sqr(max(wz*q,_kinCutoff)));
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
  Energy2 pt2=sqr((1.-z)*q)- z*sqr(_kinCutoff);
  Energy2 emittedm2 = sqr(parton->dataPtr()->constituentMass());
  // Now boost pcm and pf to z only frame
  Lorentz5Momentum p       = Lorentz5Momentum(0.0*MeV,  par.vect());
  Lorentz5Momentum n       = Lorentz5Momentum(0.0*MeV, -par.vect());
  // generate phi and compute pt of branching
  double phi = Constants::twopi*UseRandom::rnd();
  Energy pt=sqrt(pt2);
  Lorentz5Momentum qt   = LorentzMomentum(pt*cos(phi), pt*sin(phi), 0.0*MeV, 0.0*MeV);
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

