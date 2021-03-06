// -*- C++ -*-
//
// FRSDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FRSDecayer class.
//

#include "FRSDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig/Utilities/Kinematics.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

IBPtr FRSDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr FRSDecayer::fullclone() const {
  return new_ptr(*this);
}

void FRSDecayer::setDecayInfo(PDPtr incoming, PDPair outgoing,
			      vector<VertexBasePtr> vertex,
			      map<ShowerInteraction,VertexBasePtr> &,
			      const vector<map<ShowerInteraction,VertexBasePtr> > &,
			      map<ShowerInteraction,VertexBasePtr>) {
  decayInfo(incoming,outgoing);
  for(auto vert : vertex) {
    vertex_            .push_back(dynamic_ptr_cast<AbstractRFSVertexPtr>(vert));
    perturbativeVertex_.push_back(dynamic_ptr_cast<RFSVertexPtr>        (vert));
  }
}

void FRSDecayer::persistentOutput(PersistentOStream & os) const {
  os << perturbativeVertex_ << vertex_;
}

void FRSDecayer::persistentInput(PersistentIStream & is, int) {
  is >> perturbativeVertex_ >> vertex_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<FRSDecayer,GeneralTwoBodyDecayer>
describeHerwigFRSDecayer("Herwig::FRSDecayer", "Herwig.so");

void FRSDecayer::Init() {

  static ClassDocumentation<FRSDecayer> documentation
    ("The FRSDecayer class implements the decay of a fermion to "
     "a spin-3/2 fermion and a scalar.");

}

void FRSDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  bool ferm = part.id() > 0;
  // for the decaying particle
  if(ferm) {
    SpinorWaveFunction::
      constructSpinInfo(wave_,const_ptr_cast<tPPtr>(&part),incoming,true);
    RSSpinorBarWaveFunction::constructSpinInfo(RSwavebar_,decay[0],outgoing,true);
  }
  else {
    SpinorBarWaveFunction::
      constructSpinInfo(wavebar_,const_ptr_cast<tPPtr>(&part),incoming,true);
    RSSpinorWaveFunction::constructSpinInfo(RSwave_,decay[0],outgoing,true);
  }
  ScalarWaveFunction::constructSpinInfo(decay[1],outgoing,true);
}

double FRSDecayer::me2(const int,const Particle & part,
		       const tPDVector & outgoing,
		       const vector<Lorentz5Momentum> & momenta,
		       MEOption meopt) const {
  bool ferm = part.id() > 0;
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1Half,PDT::Spin3Half,PDT::Spin0)));
  if(meopt==Initialize) {
    // spinors and rho
    if(ferm) {
      SpinorWaveFunction   ::calculateWaveFunctions(wave_,rho_,
						    const_ptr_cast<tPPtr>(&part),
						    incoming);
      if(wave_[0].wave().Type() != SpinorType::u)
	for(unsigned int ix = 0; ix < 2; ++ix) wave_   [ix].conjugate();
    }
    else {
      SpinorBarWaveFunction::calculateWaveFunctions(wavebar_,rho_,
						    const_ptr_cast<tPPtr>(&part),
						    incoming);
      if(wavebar_[0].wave().Type() != SpinorType::v)
	for(unsigned int ix = 0; ix < 2; ++ix) wavebar_[ix].conjugate();
    }
    // fix rho if no correlations
    fixRho(rho_);
  }
  if(ferm)
    RSSpinorBarWaveFunction::
      calculateWaveFunctions(RSwavebar_,momenta[0],outgoing[1],Helicity::outgoing);
  else
    RSSpinorWaveFunction::
      calculateWaveFunctions(RSwave_   ,momenta[0],outgoing[1],Helicity::outgoing);
  ScalarWaveFunction scal(momenta[1],outgoing[1],Helicity::outgoing);
  Energy2 scale(sqr(part.mass()));
  for(unsigned int if1 = 0; if1 < 2; ++if1) {
    for(unsigned int if2 = 0; if2 < 4; ++if2) {
      (*ME())(if1, if2, 0) = 0.;
      for(auto vert : vertex_) {
	if(ferm) (*ME())(if1, if2, 0) +=
		   vert->evaluate(scale,wave_[if1],RSwavebar_[if2],scal);
	else     (*ME())(if1, if2, 0) += 
		   vert->evaluate(scale,RSwave_[if2],wavebar_[if1],scal);
      }
    }
  }
  double output = (ME()->contract(rho_)).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(part.dataPtr(),outgoing[0],outgoing[1]);
  // test code
//   Energy q = part.mass();
//   Energy m1 = decay[0]->mass();
//   Energy m2 = decay[1]->mass();
//   Energy2 q2(q*q),m12(m1*m1),m22(m2*m2);
//   Energy2 pcm2(0.25*(q2*(q2-2.*m12-2.*m22)+(m12-m22)*(m12-m22))/q2);
//   Energy pcm(sqrt(pcm2));
//   Energy Qp(sqrt((q+m1)*(q+m1)-m22)),Qm(sqrt((q-m1)*(q-m1)-m22));
//   double r23(sqrt(2./3.));
//   // couplings
//   Complex left  = perturbativeVertex_-> left()*perturbativeVertex_-> norm();
//   Complex right = perturbativeVertex_->right()*perturbativeVertex_-> norm();
//   complex<InvEnergy> A1 = 0.5*(left+right)*UnitRemoval::InvE;
//   complex<InvEnergy> B1 = 0.5*(right-left)*UnitRemoval::InvE;
//   complex<Energy> h1(-2.*r23*pcm*q/m1*Qm*B1);
//   complex<Energy> h2( 2.*r23*pcm*q/m1*Qp*A1);
//   cout << "testing 1/2->3/2 0 "
//        << output*scale/GeV2 << "   " 
//        << real(h1*conj(h1)+h2*conj(h2))/4./GeV2     << "   " 
//        << real(h1*conj(h1)+h2*conj(h2))/4./(output*scale) << endl;
  // return the answer
  return output;
}

Energy FRSDecayer::partialWidth(PMPair inpart, PMPair outa,
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return ZERO;
  if(perturbativeVertex_.size()==1 &&
     perturbativeVertex_[0]) {
    Energy q = inpart.second;
    Energy m1 = outa.second;
    Energy m2 = outb.second;
    Energy2 q2(q*q),m12(m1*m1),m22(m2*m2);
    Energy2 pcm2(0.25*(q2*(q2-2.*m12-2.*m22)+(m12-m22)*(m12-m22))/q2);
    Energy pcm(sqrt(pcm2));
    Energy Qp(sqrt((q+m1)*(q+m1)-m22)),Qm(sqrt((q-m1)*(q-m1)-m22));
    double r23(sqrt(2./3.));
    // couplings
    tcPDPtr in = inpart.first->CC() ? tcPDPtr(inpart.first->CC()) : inpart.first;
    perturbativeVertex_[0]->setCoupling(sqr(inpart.second), outa.first, 
				     in,  outb.first);
    Complex left  = perturbativeVertex_[0]-> left()*perturbativeVertex_[0]-> norm();
    Complex right = perturbativeVertex_[0]->right()*perturbativeVertex_[0]-> norm();
    complex<InvEnergy> A1 = 0.5*(left+right)*UnitRemoval::InvE;
    complex<InvEnergy> B1 = 0.5*(right-left)*UnitRemoval::InvE;
    complex<Energy> h1(-2.*r23*pcm*q/m1*Qm*B1);
    complex<Energy> h2( 2.*r23*pcm*q/m1*Qp*A1);
    double me2 = real(h1*conj(h1)+h2*conj(h2))/4./sqr(inpart.second);
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

