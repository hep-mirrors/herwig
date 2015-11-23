// -*- C++ -*-
//
// FRVDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FRVDecayer class.
//

#include "FRVDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig/Utilities/Kinematics.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

IBPtr FRVDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr FRVDecayer::fullclone() const {
  return new_ptr(*this);
}

void FRVDecayer::doinit() {
  perturbativeVertex_ = dynamic_ptr_cast<RFVVertexPtr>        (getVertex());
  abstractVertex_     = dynamic_ptr_cast<AbstractRFVVertexPtr>(getVertex());
  GeneralTwoBodyDecayer::doinit();
}

void FRVDecayer::persistentOutput(PersistentOStream & os) const {
  os << abstractVertex_ << perturbativeVertex_;
}

void FRVDecayer::persistentInput(PersistentIStream & is, int) {
  is >> abstractVertex_ >> perturbativeVertex_;
}

double FRVDecayer::me2(const int , const Particle & inpart,
		       const ParticleVector & decay, 
		       MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1Half,PDT::Spin3Half,PDT::Spin1)));
  // decaying fermion or antifermion
  bool ferm = inpart.id() > 0;
  // initialize
  if(meopt==Initialize) {
    // spinors and rho
    if(ferm) {
      SpinorWaveFunction   ::calculateWaveFunctions(wave_,rho_,
						    const_ptr_cast<tPPtr>(&inpart),
						    incoming);
      if(wave_[0].wave().Type() != u_spinortype)
	for(unsigned int ix = 0; ix < 2; ++ix) wave_   [ix].conjugate();
    }
    else {
      SpinorBarWaveFunction::calculateWaveFunctions(wavebar_,rho_,
						    const_ptr_cast<tPPtr>(&inpart),
						    incoming);
      if(wavebar_[0].wave().Type() != v_spinortype)
	for(unsigned int ix = 0; ix < 2; ++ix) wavebar_[ix].conjugate();
    }
  }
  // setup spin info when needed
  if(meopt==Terminate) {
    // for the decaying particle
    if(ferm) {
      SpinorWaveFunction::
	constructSpinInfo(wave_,const_ptr_cast<tPPtr>(&inpart),incoming,true);
      RSSpinorBarWaveFunction::constructSpinInfo(RSwavebar_,decay[0],outgoing,true);
    }
    else {
      SpinorBarWaveFunction::
	constructSpinInfo(wavebar_,const_ptr_cast<tPPtr>(&inpart),incoming,true);
      RSSpinorWaveFunction::constructSpinInfo(RSwave_,decay[0],outgoing,true);
    }
    VectorWaveFunction::
      constructSpinInfo(vector_,decay[1],outgoing,true,false);
  }
  Energy2 scale(sqr(inpart.mass()));
  if(ferm)
    RSSpinorBarWaveFunction::
      calculateWaveFunctions(RSwavebar_,decay[0],outgoing);
  else
    RSSpinorWaveFunction::
      calculateWaveFunctions(RSwave_   ,decay[0],outgoing);
  bool massless = decay[1]->dataPtr()->mass()==ZERO;
  VectorWaveFunction::
    calculateWaveFunctions(vector_,decay[1],outgoing,massless);
  // loop over helicities
  for(unsigned int if1 = 0; if1 < 2; ++if1) {
    for(unsigned int if2 = 0; if2 < 4; ++if2) {
      for(unsigned int vhel = 0; vhel < 3; ++vhel) {
	if(massless && vhel == 1) ++vhel;
	if(ferm)
	  (*ME())(if1, if2,vhel) = 
	    abstractVertex_->evaluate(scale,wave_[if1],
				      RSwavebar_[if2],vector_[vhel]);
	else
	  (*ME())(if1, if2, vhel) = 
	    abstractVertex_->evaluate(scale,RSwave_[if2],
				      wavebar_[if1],vector_[vhel]);
	
      }
    }
  }
  double output=(ME()->contract(rho_)).real()/scale*UnitRemoval::E2;
  // test
//   Energy m1(inpart.mass()),m2(decay[0]->mass()),m3(decay[1]->mass());
//   Energy2 m12(m1*m1),m22(m2*m2),m32(m3*m3);
//   Energy Qp(sqrt(sqr(m1+m2)-sqr(m3))),Qm(sqrt(sqr(m1-m2)-sqr(m3)));
//   double r2(sqrt(2.)),r3(sqrt(3.));
//   Energy pcm(Kinematics::pstarTwoBodyDecay(m1,m2,m3));
//   vector<Complex> left  = perturbativeVertex_-> left();
//   vector<Complex> right = perturbativeVertex_->right();
//   Complex A1 = 0.5*(left [0]+right[0])*perturbativeVertex_-> norm();
//   Complex B1 = 0.5*(right[0]- left[0])*perturbativeVertex_-> norm();
//   complex<InvEnergy> A2 = 0.5*(left [1]+right[1])*perturbativeVertex_-> norm()*UnitRemoval::InvE;
//   complex<InvEnergy> B2 = 0.5*(right[1]- left[1])*perturbativeVertex_-> norm()*UnitRemoval::InvE;
//   complex<InvEnergy2> A3 = 0.5*(left [2]+right[2])*perturbativeVertex_-> norm()*UnitRemoval::InvE2;
//   complex<InvEnergy2> B3 = 0.5*(right[2]- left[2])*perturbativeVertex_-> norm()*UnitRemoval::InvE2;
//   complex<Energy> h1(-2.*Qp*A1),h2(2.*Qm*B1);
//   complex<Energy> h3(-2./r3*Qp*(A1-Qm*Qm/m2*A2));
//   complex<Energy> h4( 2./r3*Qm*(B1-Qp*Qp/m2*B2));
//   complex<Energy> h5(ZERO),h6(ZERO);
//   if(decay[1]->mass()>ZERO) {
//     h5 = -2.*r2/r3/m2/m3*Qp*(0.5*(m12-m22-m32)*A1+0.5*Qm*Qm*(m1+m2)*A2
//     			   +m12*pcm*pcm*A3);
//     h6 =  2.*r2/r3/m2/m3*Qm*(0.5*(m12-m22-m32)*B1-0.5*Qp*Qp*(m1-m2)*B2
// 					   +m12*pcm*pcm*B3);
//   }
//   cout << "testing 1/2->3/2 1 " << inpart.id() << " "
//        << output << "   " 
//        << 0.25*(h1*conj(h1)+h2*conj(h2)+h3*conj(h3)+
// 		h4*conj(h4)+h5*conj(h5)+h6*conj(h6))/sqr(inpart.mass()) << "   " 
//        << 0.25*(h1*conj(h1)+h2*conj(h2)+h3*conj(h3)+
// 		h4*conj(h4)+h5*conj(h5)+h6*conj(h6))/sqr(inpart.mass())/output << endl;
  // colour and identical particle factors
  output *= colourFactor(inpart.dataPtr(),decay[0]->dataPtr(),decay[1]->dataPtr());
  // return the answer
  return output;
}

Energy FRVDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return ZERO;
  if(perturbativeVertex_) {
    Energy m1(inpart.second),m2(outa.second),m3(outb.second);
    Energy2 m12(m1*m1),m22(m2*m2),m32(m3*m3);
    Energy Qp(sqrt(sqr(m1+m2)-sqr(m3))),Qm(sqrt(sqr(m1-m2)-sqr(m3)));
    double r2(sqrt(2.)),r3(sqrt(3.));
    Energy pcm(Kinematics::pstarTwoBodyDecay(m1,m2,m3));
    // couplings
    tcPDPtr in = inpart.first->CC() ? tcPDPtr(inpart.first->CC()) : inpart.first;
    perturbativeVertex_->setCoupling(sqr(inpart.second), outa.first, 
				     in,  outb.first);
    vector<Complex> left  = perturbativeVertex_-> left();
    vector<Complex> right = perturbativeVertex_->right();
    Complex A1 = 0.5*(left [0]+right[0])*perturbativeVertex_-> norm();
    Complex B1 = 0.5*(right[0]- left[0])*perturbativeVertex_-> norm();
    complex<InvEnergy> A2 = 0.5*(left [1]+right[1])*perturbativeVertex_-> norm()*UnitRemoval::InvE;
    complex<InvEnergy> B2 = 0.5*(right[1]- left[1])*perturbativeVertex_-> norm()*UnitRemoval::InvE;
    complex<InvEnergy2> A3 = 0.5*(left [2]+right[2])*perturbativeVertex_-> norm()*UnitRemoval::InvE2;
    complex<InvEnergy2> B3 = 0.5*(right[2]- left[2])*perturbativeVertex_-> norm()*UnitRemoval::InvE2;
    complex<Energy> h1(-2.*Qp*A1),h2(2.*Qm*B1);
    complex<Energy> h3(-2./r3*Qp*(A1-Qm*Qm/m2*A2));
    complex<Energy> h4( 2./r3*Qm*(B1-Qp*Qp/m2*B2));
    complex<Energy> h5(ZERO),h6(ZERO);
    if(outb.second>ZERO) {
      h5 = -2.*r2/r3/m2/m3*Qp*(0.5*(m12-m22-m32)*A1+0.5*Qm*Qm*(m1+m2)*A2
			       +m12*pcm*pcm*A3);
      h6 =  2.*r2/r3/m2/m3*Qm*(0.5*(m12-m22-m32)*B1-0.5*Qp*Qp*(m1-m2)*B2
			       +m12*pcm*pcm*B3);
    }
    double me2 = 0.25*real(h1*conj(h1)+h2*conj(h2)+h3*conj(h3)+
			   h4*conj(h4)+h5*conj(h5)+h6*conj(h6))/sqr(inpart.second);
    Energy output = me2*pcm/8./Constants::pi;
    // colour factor
    output *= colourFactor(inpart.first,outa.first,outb.first);
    return output;
  }
  else {
    return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
  }
}

ClassDescription<FRVDecayer> FRVDecayer::initFRVDecayer;
// Definition of the static class description member.

void FRVDecayer::Init() {

  static ClassDocumentation<FRVDecayer> documentation
    ("The FRVDecayer class handles the decay of a fermion to "
     "a spin-3/2 particle and a vector boson.");

}
