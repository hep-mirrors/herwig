// -*- C++ -*-
//
// ForcedSplitting.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
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
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "Herwig++/Utilities/EnumParticles.h"
#include "Herwig++/Shower/ShowerHandler.h"
#include <cassert>

using namespace Herwig;

void ForcedSplitting::persistentOutput(PersistentOStream & os) const {
  os << ounit(_kinCutoff, GeV) << _range
     << _zbin << _ybin << _nbinmax << _alpha;
}

void ForcedSplitting::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_kinCutoff, GeV) >> _range
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

}


// This creates the parton to split and sets it momentum and parent/child
// relationships
PPtr ForcedSplitting::forceSplit(const tRemPPtr rem, long child, Energy &oldQ, 
				 double &oldx, Lorentz5Momentum &pf, 
				 Lorentz5Momentum &p,
				 const unsigned int iopt,
				 const tStepPtr step) const {
  Lorentz5Momentum beam = _beam->momentum();
  PPtr parton = new_ptr(Particle(getParticleData(child)));
  Lorentz5Momentum partonp = emit(beam,oldQ,oldx,parton,pf,iopt);
  p += partonp;
  parton->set5Momentum(partonp);
  step->addDecayProduct(rem,parton,false);
  return parton;
}

// This forces the final output of the remnant ((di)quark) and sets the
// momentum and parent/child relationships
PPtr ForcedSplitting::finalSplit(const tRemPPtr rem, long remID, 
				 Lorentz5Momentum usedMomentum, 
				 const tStepPtr step) const {
  // Create the remnant and set its momentum, also reset all of the decay 
  // products from the hadron
  PPtr remnant = new_ptr(Particle(getParticleData(remID)));
  Lorentz5Momentum prem(rem->momentum()-usedMomentum);
  prem.setMass(getParticleData(remID)->constituentMass());
  prem.rescaleEnergy();
  remnant->set5Momentum(prem);
  // Add the remnant to the step, but don't do colour connections
  step->addDecayProduct(rem,remnant,false);
  return remnant;
}

// This defines the momentum for an emitted parton, currently no pt is
// given to the produced partons, z is generated uniformly.
Lorentz5Momentum ForcedSplitting::emit(const Lorentz5Momentum &par,
				       Energy &lastQ, double &lastx, 
				       PPtr parton,
				       Lorentz5Momentum &pf,
				       const unsigned int iopt) const {
  assert(iopt==1||iopt==2);
  Ptr<BeamParticleData>::const_pointer beam = 
    dynamic_ptr_cast<Ptr<BeamParticleData>::const_pointer>(_beam->dataPtr());

  if(!_pdf)
    throw Exception() << "No PDF object present in "
		      << "ForcedSplitting::emit(...)"
		      << Exception::runerror; 

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
  //  unsigned int nz=std::min(int(_ybin*dely+1),_nbinmax);
  unsigned int nz=_nbinmax;
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
      double pdfval(0.0);
      // SG modified this, should be x/z rather than x/(1-z)! 
      // SP as q is always less than forcedSplitScale, the pdf scale is fixed
      // pdfval=_pdf->xfx(beam,gluon,sqr(q),lastx*zr);
      pdfval=_pdf->xfx(beam,gluon,sqr(_forcedSplitScale),lastx*zr);
      // SG: this is symmetric in z <-> 1-z
      psum+=pdfval*az*0.5*(sqr(zz)+sqr(wz));
    }
    // q -> q g
    else {
      // calculate splitting function
      double pdfval(0.0);
      // SG modified this, should be x/z rather than x/(1-z)! 
      // SP as q is always less than forcedSplitScale, the pdf scale is fixed
      // pdfval=_pdf->xfx(beam,parton->dataPtr(),sqr(q),lastx*zr);
      pdfval=_pdf->xfx(beam,parton->dataPtr(),sqr(_forcedSplitScale),lastx*zr);
      // SG this splitting function has to have a 1/z pole! 
      psum+=pdfval*az*4./3.*(1.+sqr(wz))*zr;
    }
    prob.push_back(psum);
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
  Lorentz5Momentum p = Lorentz5Momentum(0.0*MeV,  par.vect());
  Lorentz5Momentum n = Lorentz5Momentum(0.0*MeV, -par.vect());
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

