// -*- C++ -*-
//
// SRFDecayer.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SRFDecayer class.
//

#include "SRFDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig++/Utilities/Kinematics.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

IBPtr SRFDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr SRFDecayer::fullclone() const {
  return new_ptr(*this);
}

void SRFDecayer::doinit() {
  perturbativeVertex_ = dynamic_ptr_cast<RFSVertexPtr>        (getVertex());
  abstractVertex_     = dynamic_ptr_cast<AbstractRFSVertexPtr>(getVertex());
  GeneralTwoBodyDecayer::doinit();
}

void SRFDecayer::persistentOutput(PersistentOStream & os) const {
  os << abstractVertex_ << perturbativeVertex_;
}

void SRFDecayer::persistentInput(PersistentIStream & is, int) {
  is >> abstractVertex_ >> perturbativeVertex_;
}

ClassDescription<SRFDecayer> SRFDecayer::initSRFDecayer;
// Definition of the static class description member.

void SRFDecayer::Init() {

  static ClassDocumentation<SRFDecayer> documentation
    ("This class implements to decay of a scalar to a spin-3/2 and"
     " spin-1/2 fermion");

}

double SRFDecayer::me2(const int , const Particle & inpart,
		       const ParticleVector & decay,MEOption meopt) const {
  unsigned int irs=0,ifm=1;
  if(decay[0]->dataPtr()->iSpin()==PDT::Spin1Half) swap(irs,ifm);
  bool ferm = decay[ifm]->id()<0;
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(rho_,const_ptr_cast<tPPtr>(&inpart),incoming);
    swave_ = ScalarWaveFunction(inpart.momentum(),inpart.dataPtr(),incoming);
    if(irs==0)
      ME(DecayMatrixElement(PDT::Spin0,PDT::Spin3Half,PDT::Spin1Half));
    else
      ME(DecayMatrixElement(PDT::Spin0,PDT::Spin1Half,PDT::Spin3Half));
  }
  if(meopt==Terminate) {
    ScalarWaveFunction::
      constructSpinInfo(const_ptr_cast<tPPtr>(&inpart),incoming,true);
    if(ferm) {
      RSSpinorBarWaveFunction::
	constructSpinInfo(RSwavebar_,decay[irs],outgoing,true);
      SpinorWaveFunction::
	constructSpinInfo(wave_     ,decay[ifm],outgoing,true);
    }
    else {
      RSSpinorWaveFunction::
	constructSpinInfo(RSwave_ ,decay[irs],outgoing,true);
      SpinorBarWaveFunction::
	constructSpinInfo(wavebar_,decay[ifm],outgoing,true);
    }
    return 0.;
  }
  if(ferm) {
    RSSpinorBarWaveFunction::
      calculateWaveFunctions(RSwavebar_,decay[irs],outgoing);
    SpinorWaveFunction::
      calculateWaveFunctions(wave_     ,decay[ifm],outgoing);
  }
  else {
    RSSpinorWaveFunction::
      calculateWaveFunctions(RSwave_ ,decay[irs],outgoing);
    SpinorBarWaveFunction::
      calculateWaveFunctions(wavebar_,decay[ifm],outgoing);
  }
  Energy2 scale(sqr(inpart.mass()));
  for(unsigned int ifm = 0; ifm < 4; ++ifm){
    for(unsigned int ia = 0; ia < 2; ++ia) {
      if(irs==0) {
	if(ferm)
	  ME()(0, ifm, ia) = abstractVertex_->evaluate(scale,wave_[ia],
						       RSwavebar_[ifm],swave_);
	else
	  ME()(0, ifm, ia) = abstractVertex_->evaluate(scale,RSwave_[ifm],
						       wavebar_[ia],swave_);
      }
      else {
	if(ferm)
	  ME()(0, ia, ifm) = abstractVertex_->evaluate(scale,wave_[ia],
						       RSwavebar_[ifm],swave_);
	else
	  ME()(0, ia, ifm) = abstractVertex_->evaluate(scale,RSwave_[ifm],
						       wavebar_[ia],swave_);
      }
    }
  }
  double output = (ME().contract(rho_)).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(inpart.dataPtr(),decay[irs]->dataPtr(),
			 decay[ifm]->dataPtr());
  // return the answer
  return output;
}

Energy SRFDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return ZERO;
  if(perturbativeVertex_) {
    Energy q = inpart.second;
    Energy m1 = outa.second, m2 = outb.second;
    // couplings
    tcPDPtr in = inpart.first->CC() ? tcPDPtr(inpart.first->CC()) : inpart.first;
    if(outa.first->iSpin()==PDT::Spin1Half) {
      swap(m1,m2);
      perturbativeVertex_->setCoupling(sqr(inpart.second),outb.first,
				       outa.first, in);
    }
    else {
      perturbativeVertex_->setCoupling(sqr(inpart.second),outa.first,
				       outb.first, in);
    }
    Complex left  = perturbativeVertex_-> left()*perturbativeVertex_-> norm();
    Complex right = perturbativeVertex_->right()*perturbativeVertex_-> norm();
    complex<InvEnergy> A1 = 0.5*(left+right)*UnitRemoval::InvE;
    complex<InvEnergy> B1 = 0.5*(right-left)*UnitRemoval::InvE;
    Energy2 q2(q*q),m12(m1*m1),m22(m2*m2);
    Energy2 pcm2(0.25*(q2*(q2-2.*m12-2.*m22)+(m12-m22)*(m12-m22))/q2);
    Energy pcm(sqrt(pcm2));
    Energy Qp(sqrt(-sqr(m2+m1)+q2)),Qm(sqrt(-sqr(m2-m1)+q2));
    double r23(sqrt(2./3.));
    complex<Energy> h1(-2.*r23*pcm*q/m1*Qm*B1);
    complex<Energy> h2( 2.*r23*pcm*q/m1*Qp*A1);
    double me2 = real(h1*conj(h1)+h2*conj(h2))/2./sqr(inpart.second);
    Energy output = me2*pcm/8./Constants::pi;
    // colour factor
    output *= colourFactor(inpart.first,outa.first,outb.first);
    // return the answer
    return output;
  }
  else {
    return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
  }
}

