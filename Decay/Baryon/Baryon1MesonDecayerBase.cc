// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Baryon1MesonDecayerBase class.
//

#include "Baryon1MesonDecayerBase.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig/Utilities/Kinematics.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/RSSpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/RSSpinorBarWaveFunction.h"
#include "ThePEG/Helicity/LorentzPolarizationVector.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

AbstractNoPIOClassDescription<Baryon1MesonDecayerBase> 
Baryon1MesonDecayerBase::initBaryon1MesonDecayerBase;
// Definition of the static class description member.

void Baryon1MesonDecayerBase::Init() {

  static ClassDocumentation<Baryon1MesonDecayerBase> documentation
    ("The Baryon1MesonDecayerBase class is the base class for"
     " the decays of the baryons to a baryon and a pseudoscalar or vector meson.");

}

// return the matrix element squared for a given mode and phase-space channel
// (inherited from DecayIntegrator and implemented here)
double Baryon1MesonDecayerBase::me2(const int ichan,
				    const Particle & inpart,
				    const ParticleVector & decay,
				    MEOption meopt) const {
  double me(0.);
  // decide which matrix element we are doing
  // incoming spin-1/2 particle
  if(inpart.dataPtr()->iSpin()==2) {
    // decay to spin-1/2 particle
    if(decay[0]->dataPtr()->iSpin()==2) {
      // scalar meson
      if(decay[1]->dataPtr()->iSpin()==1)
	me=halfHalfScalar(ichan,inpart,decay,meopt);
      // vector meson
      else if(decay[1]->dataPtr()->iSpin()==3)
	me=halfHalfVector(ichan,inpart,decay,meopt);
      else
	throw DecayIntegratorError() << "Unknown outgoing meson spin in "
				     << "Baryon1MesonDecayerBase::me2()" 
				     << Exception::abortnow;
    }
    // decay to spin-3/2 particle
    else if(decay[0]->dataPtr()->iSpin()==4) {
      // scalar meson
      if(decay[1]->dataPtr()->iSpin()==1)
	me=halfThreeHalfScalar(ichan,inpart,decay,meopt);
      // vector meson
      else if(decay[1]->dataPtr()->iSpin()==3)
	me=halfThreeHalfVector(ichan,inpart,decay,meopt);
      else
	throw DecayIntegratorError() << "Unknown outgoing meson spin in "
				     << "Baryon1MesonDecayerBase::me2()" 
				     << Exception::abortnow;
    }
    // unknown
    else
      throw DecayIntegratorError() << "Unknown outgoing baryon spin in "
				   << "Baryon1MesonDecayerBase::me2()" 
				   << Exception::abortnow;
  }
  // incoming spin-3/2 particle
  else if(inpart.dataPtr()->iSpin()==4) {
    // decay to spin-1/2 particle
    if(decay[0]->dataPtr()->iSpin()==2) {
      // scalar meson
      if(decay[1]->dataPtr()->iSpin()==1)
	me=threeHalfHalfScalar(ichan,inpart,decay,meopt);
      // vector meson
      else if(decay[1]->dataPtr()->iSpin()==3)
	me=threeHalfHalfVector(ichan,inpart,decay,meopt);
      else
	throw DecayIntegratorError() << "Unknown outgoing meson spin in "
				     << "Baryon1MesonDecayerBase::me2()" 
				     << Exception::abortnow;
    }
    // decay to spin-3/2 particle
    else if(decay[0]->dataPtr()->iSpin()==4) {
      // scalar meson
      if(decay[1]->dataPtr()->iSpin()==1)
	me=threeHalfThreeHalfScalar(ichan,inpart,decay,meopt);
      else
	throw DecayIntegratorError() << "Unknown outgoing meson spin in "
				     << "Baryon1MesonDecayerBase::me2()" 
				     << Exception::abortnow;
    }
    // unknown
    else
      throw DecayIntegratorError() << "Unknown outgoing baryon spin in "
				   << "Baryon1MesonDecayerBase::me2()" 
				   << Exception::abortnow;
  }
  // unknown
  else
    throw DecayIntegratorError() << "Unknown incoming spin in "
				 << "Baryon1MesonDecayerBase::me2()" 
				 << Exception::abortnow;
  return me;
}

// matrix element for the decay of a spin-1/2 fermion to a spin-1/2 fermion and
// a pseudoscalar meson
double Baryon1MesonDecayerBase::
halfHalfScalar(const int,const Particle & inpart,
	       const ParticleVector & decay,MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin0)));
  // spinors etc for the decaying particle
  if(meopt==Initialize) {
    // spinors and rho
    if(inpart.id()>0)
      SpinorWaveFunction   ::calculateWaveFunctions(_inHalf,_rho,
						    const_ptr_cast<tPPtr>(&inpart),
						    incoming);
    else
      SpinorBarWaveFunction::calculateWaveFunctions(_inHalfBar,_rho,
						    const_ptr_cast<tPPtr>(&inpart),
						    incoming);
    // matrix element
  }
  // setup spin info when needed
  if(meopt==Terminate) {
    // for the decaying particle
    if(inpart.id()>0) {
      SpinorWaveFunction::
	constructSpinInfo(_inHalf,const_ptr_cast<tPPtr>(&inpart),incoming,true);
      SpinorBarWaveFunction::constructSpinInfo(_inHalfBar,decay[0],outgoing,true);
    }
    else {
      SpinorBarWaveFunction::
	constructSpinInfo(_inHalfBar,const_ptr_cast<tPPtr>(&inpart),incoming,true);
      SpinorWaveFunction::constructSpinInfo(_inHalf,decay[0],outgoing,true);
    }
    ScalarWaveFunction::constructSpinInfo(decay[1],outgoing,true);
    return 0.;
  }
  // spinors for the decay product
  if(inpart.id()>0) {
    SpinorBarWaveFunction::calculateWaveFunctions(_inHalfBar,decay[0],outgoing);
  }
  else {
    SpinorWaveFunction   ::calculateWaveFunctions(_inHalf,decay[0],outgoing);
  }
  // get the couplings
  Complex A,B;
  halfHalfScalarCoupling(imode(),inpart.mass(),decay[0]->mass(),decay[1]->mass(),A,B);
  Complex left,right,meout;
  // coupling for an incoming particle
  if(inpart.id()>0) {
    left  = (A-B);
    right = (A+B);
  }
  // coupling for an incoming antiparticle
  else {
    left  = conj(A+B);
    right = conj(A-B);
  }
  // calculate the matrix element
  vector<unsigned int> ispin(3,0);
  unsigned int ix,iy;
//   Complex output(0.);
  for(ix=0;ix<2;++ix) {
    for(iy=0;iy<2;++iy) {
      if(decay[0]->id()>0){ispin[0]=iy;ispin[1]=ix;}
      else{ispin[0]=ix;ispin[1]=iy;}
      (*ME())(ispin)=_inHalf[iy].generalScalar(_inHalfBar[ix],left,right)/inpart.mass();
//       output += norm(ME()(ispin));
    }
  }
  // test of the matrix elemen// t
//   Energy m1(inpart.mass()),m2(decay[0]->mass()),m3(decay[1]->mass());
//   Energy Qp(sqrt(sqr(m1+m2)-sqr(m3))),Qm(sqrt(sqr(m1-m2)-sqr(m3)));
//   Complex h1(2.*Qp*A/inpart.mass()),h2(-2.*Qm*B/inpart.mass());
//   generator()->log() << "testing 1/2->1/2 0 " 
//  		     << 0.5*output << "   " 
//  		     << 0.25*(h1*conj(h1)+h2*conj(h2)) << "   " 
//  		     << 0.5*(h1*conj(h1)+h2*conj(h2))/output << endl;
//   generator()->log() << "testing alpha " << 
//     (norm(0.5*(h1+h2))-norm(0.5*(h1-h2)))/
//     (norm(0.5*(h1+h2))+norm(0.5*(h1-h2))) << "\n";
//   Energy pcm(Kinematics::pstarTwoBodyDecay(m1,m2,m3));
//   generator()->log() << "testing masses " << m1/GeV << " " << m2/GeV << " " << m3/GeV
// 		     << "\n";
//   generator()->log() << "testing partial " << pcm*0.5*output/8./Constants::pi/MeV
// 		     << "\n";
  // store the matrix element
  return (ME()->contract(_rho)).real();
}

// matrix element for the decay of a spin-1/2 fermion to a spin-1/2 fermion and
// a vector meson
double Baryon1MesonDecayerBase::
halfHalfVector(const int,const Particle & inpart,
	       const ParticleVector & decay,MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin1)));
  // check if the outgoing meson is really a photon
  bool photon=decay[1]->id()==ParticleID::gamma;
  // spinors etc for the decaying particle
  if(meopt==Initialize) {
    // spinors and rho
    if(inpart.id()>0)
      SpinorWaveFunction   ::calculateWaveFunctions(_inHalf,_rho,
						    const_ptr_cast<tPPtr>(&inpart),
						    incoming);
    else
      SpinorBarWaveFunction::calculateWaveFunctions(_inHalfBar,_rho,
						    const_ptr_cast<tPPtr>(&inpart),
						    incoming);
    // matrix element
  }
  // setup spin info when needed
  if(meopt==Terminate) {
    // for the decaying particle
    if(inpart.id()>0) {
      SpinorWaveFunction::
	constructSpinInfo(_inHalf,const_ptr_cast<tPPtr>(&inpart),incoming,true);
      SpinorBarWaveFunction::constructSpinInfo(_inHalfBar,decay[0],outgoing,true);

    }
    else {
      SpinorBarWaveFunction::
	constructSpinInfo(_inHalfBar,const_ptr_cast<tPPtr>(&inpart),incoming,true);
      SpinorWaveFunction::constructSpinInfo(_inHalf,decay[0],outgoing,true);
    }
    VectorWaveFunction::constructSpinInfo(_inVec,decay[1],outgoing,true,photon);
    return 0.;
  }
  // spinors for the decay product
  if(inpart.id()>0) {
    SpinorBarWaveFunction::calculateWaveFunctions(_inHalfBar,decay[0],outgoing);
  }
  else {
    SpinorWaveFunction   ::calculateWaveFunctions(_inHalf,decay[0],outgoing);
  }
  VectorWaveFunction::calculateWaveFunctions(_inVec,decay[1],outgoing,photon);
  // get the couplings
  Complex A1,A2,B1,B2;
  halfHalfVectorCoupling(imode(),inpart.mass(),decay[0]->mass(),decay[1]->mass(),
			 A1,A2,B1,B2);
  Complex lS,rS,lV,rV;
  complex<Energy> scalar;
  // couplings for an incoming particle
  if(inpart.id()>0) {
    lS = (A2-B2);
    rS = (A2+B2);
    lV = (A1-B1);
    rV = (A1+B1);
  }
  else {
    lS = -conj(A2+B2);
    rS = -conj(A2-B2);
    lV =  conj(A1-B1);
    rV =  conj(A1+B1);
  }
  // calculate the matrix element
  // decide which type of mode to do
  Energy msum(inpart.mass()+decay[0]->mass());
  vector<unsigned int> ispin(3);
  LorentzVector<complex<Energy> > svec;
  Complex prod;
//   Complex output(0.);
  unsigned int ix,iy;
  for(ix=0;ix<2;++ix) {
    for(iy=0;iy<2;++iy) {
      // scalar like piece
      scalar = _inHalf[iy].generalScalar(_inHalfBar[ix],lS,rS);
      // vector like piece
      svec   = _inHalf[iy].generalCurrent(_inHalfBar[ix],lV,rV);
      if(decay[0]->id()>0) {
	ispin[0] = iy;
	ispin[1] = ix;
      }
      else {
	ispin[0] = ix;
	ispin[1] = iy;
      }
      for(ispin[2]=0;ispin[2]<3;++ispin[2]) {
	ispin[2]=ispin[2];
	prod=_inVec[ispin[2]].dot(inpart.momentum())/msum;
	(*ME())(ispin)=(svec.dot(_inVec[ispin[2]])+prod*scalar)/inpart.mass();
// 	output += norm(ME()(ispin));
      }
    }
  }
  // test of the matrix element
//   Energy m1(inpart.mass()),m2(decay[0]->mass()),m3(decay[1]->mass());
//   Energy Qp(sqrt(sqr(m1+m2)-sqr(m3))),Qm(sqrt(sqr(m1-m2)-sqr(m3)));
//   double r2(sqrt(2.));
//   Energy pcm(Kinematics::pstarTwoBodyDecay(m1,m2,m3));
//   Complex h1(2.*r2*Qp*B1/inpart.mass()),h2(-2.*r2*Qm*A1/inpart.mass()),
//     h3(2./m3*(Qp*(m1-m2)*B1-Qm*m1*B2*pcm/(m1+m2))/inpart.mass()),
//     h4(2./m3*(Qm*(m1+m2)*A1+Qp*m1*A2*pcm/(m1+m2))/inpart.mass());
//   generator()->log() << "testing 1/2->1/2 1 " 
// 		     << 0.5*output << "   " 
// 		     << 0.25*(h1*conj(h1)+h2*conj(h2)+h3*conj(h3)+h4*conj(h4)) << "   " 
// 		     << 0.50*(h1*conj(h1)+h2*conj(h2)+h3*conj(h3)+h4*conj(h4))/output 
// 		     << "\n";
//   generator()->log() << "alpha = " << 2.*(norm(h3)+norm(h4))/(norm(h1)+norm(h2))-1.
// 		     << "\n";
  // return the answer
  return (ME()->contract(_rho)).real();
}

// matrix element for the decay of a spin-1/2 fermion to a spin-3/2 fermion and
// a scalar meson
double Baryon1MesonDecayerBase::halfThreeHalfScalar(const int,
						    const Particle & inpart,
						    const ParticleVector & decay,
						    MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half,PDT::Spin3Half,PDT::Spin0)));
  // spinors etc for the decaying particle
  if(meopt==Initialize) {
    // spinors and rho
    if(inpart.id()>0)
      SpinorWaveFunction   ::calculateWaveFunctions(_inHalf,_rho,
						    const_ptr_cast<tPPtr>(&inpart),
						    incoming);
    else
      SpinorBarWaveFunction::calculateWaveFunctions(_inHalfBar,_rho,
						    const_ptr_cast<tPPtr>(&inpart),
						    incoming);
    // matrix element
  }
  // setup spin info when needed
  if(meopt==Terminate) {
    // for the decaying particle
    if(inpart.id()>0) {
      SpinorWaveFunction::
	constructSpinInfo(_inHalf,const_ptr_cast<tPPtr>(&inpart),incoming,true);
      RSSpinorBarWaveFunction::constructSpinInfo(_inThreeHalfBar,
						 decay[0],outgoing,true);

    }
    else {
      SpinorBarWaveFunction::
	constructSpinInfo(_inHalfBar,const_ptr_cast<tPPtr>(&inpart),incoming,true);
      RSSpinorWaveFunction::constructSpinInfo(_inThreeHalf,
					      decay[0],outgoing,true);
    }
    ScalarWaveFunction::constructSpinInfo(decay[1],outgoing,true);
    return 0.;
  }
  // spinors for the decay product
  LorentzPolarizationVector in=UnitRemoval::InvE*inpart.momentum();
  if(inpart.id()>0) {
    RSSpinorBarWaveFunction::
      calculateWaveFunctions(_inThreeHalfBar,decay[0],outgoing);
    _inHalfBar.resize(_inThreeHalfBar.size());
    for(unsigned int ix=0;ix<_inThreeHalfBar.size();++ix)
      _inHalfBar[ix] = _inThreeHalfBar[ix].dot(in);
  }
  else {
    RSSpinorWaveFunction::
      calculateWaveFunctions(_inThreeHalf,decay[0],outgoing);
    _inHalf.resize(_inThreeHalf.size());
    for(unsigned int ix=0;ix<_inThreeHalf.size();++ix)
      _inHalf[ix] = _inThreeHalf[ix].dot(in);
  }
  // get the couplings
  Complex A,B,left,right;
  Energy msum(inpart.mass()+decay[0]->mass());
  halfThreeHalfScalarCoupling(imode(),inpart.mass(),decay[0]->mass(),decay[1]->mass(),
			      A,B);
  // incoming particle
  if(inpart.id()>0) {
    left=(A-B);
    right=(A+B);
  }
  // incoming anti-particle
  else {
    left=conj(A+B);
    right=conj(A-B);
  }
  vector<unsigned int> ispin(3,0);
  //Complex output(0.);
  for(unsigned ixa=0;ixa<2;++ixa) {
    for(unsigned int iya=0;iya<4;++iya) {
      unsigned int ix(iya),iy(ixa);
      if(decay[0]->id()<0) swap(ix,iy);
      ispin[0]=ixa;
      ispin[1]=iya;
      complex<double> value = _inHalf[iy].generalScalar(_inHalfBar[ix],left,right)
	*UnitRemoval::E/inpart.mass()/msum;
      (*ME())(ispin) = value;
      //output+= norm(ME()(ispin));
    }
  }
  double output = (ME()->contract(_rho)).real();
  // test of the matrix element
//   Energy m1(inpart.mass()),m2(decay[0]->mass()),m3(decay[1]->mass());
//   Energy Qp(sqrt(sqr(m1+m2)-sqr(m3))),Qm(sqrt(sqr(m1-m2)-sqr(m3)));
//   double r23(sqrt(2./3.));
//   Energy pcm(Kinematics::pstarTwoBodyDecay(m1,m2,m3));
//   complex<Energy> h1(-2.*r23*pcm*m1/m2*Qm*B/(m1+m2)),h2( 2.*r23*pcm*m1/m2*Qp*A/(m1+m2));
//   cout << "testing 1/2->3/2 0 " << inpart.id() << " "
//        << output << "   " 
//        << 0.25*(h1*conj(h1)+h2*conj(h2))/sqr(inpart.mass()) << "   " 
//        << 0.25*(h1*conj(h1)+h2*conj(h2))/sqr(inpart.mass())/output << endl;
  // return the answer
  return output;
}

// matrix element for the decay of a spin-1/2 fermion to a spin-3/2 fermion and
// a vector meson
double Baryon1MesonDecayerBase::
halfThreeHalfVector(const int,const Particle & inpart,
		    const ParticleVector & decay, MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half,PDT::Spin3Half,PDT::Spin1)));
  // check if the outgoing meson is really a photon
  bool photon=decay[1]->id()==ParticleID::gamma;
  // spinors etc for the decaying particle
  if(meopt==Initialize) {
    // spinors and rho
    if(inpart.id()>0)
      SpinorWaveFunction   ::calculateWaveFunctions(_inHalf,_rho,
						    const_ptr_cast<tPPtr>(&inpart),
						    incoming);
    else
      SpinorBarWaveFunction::calculateWaveFunctions(_inHalfBar,_rho,
						    const_ptr_cast<tPPtr>(&inpart),
						    incoming);
    // matrix element
  }
  // setup spin info when needed
  if(meopt==Terminate) {
    // for the decaying particle
    if(inpart.id()>0) {
      SpinorWaveFunction::
	constructSpinInfo(_inHalf,const_ptr_cast<tPPtr>(&inpart),incoming,true);
      RSSpinorBarWaveFunction::constructSpinInfo(_inThreeHalfBar,
						 decay[0],outgoing,true);

    }
    else {
      SpinorBarWaveFunction::
	constructSpinInfo(_inHalfBar,const_ptr_cast<tPPtr>(&inpart),incoming,true);
      RSSpinorWaveFunction::constructSpinInfo(_inThreeHalf,
					      decay[0],outgoing,true);
    }
    VectorWaveFunction::constructSpinInfo(_inVec,decay[1],outgoing,true,photon);
    return 0.;
  }
  LorentzPolarizationVector in=UnitRemoval::InvE*inpart.momentum();
  if(inpart.id()>0) {
    RSSpinorBarWaveFunction::
      calculateWaveFunctions(_inThreeHalfBar,decay[0],outgoing);
    _inHalfBar.resize(_inThreeHalfBar.size());
    for(unsigned int ix=0;ix<_inThreeHalfBar.size();++ix)
      _inHalfBar[ix] = _inThreeHalfBar[ix].dot(in);
    }
  else {
    RSSpinorWaveFunction::
      calculateWaveFunctions(_inThreeHalf,decay[0],outgoing);
    _inHalf.resize(_inThreeHalf.size());
    for(unsigned int ix=0;ix<_inThreeHalf.size();++ix)
      _inHalf[ix] = _inThreeHalf[ix].dot(in);
  }
  ME()->zero();
  VectorWaveFunction::calculateWaveFunctions(_inVec,decay[1],outgoing,photon);
  // get the couplings
  Complex A1,A2,A3,B1,B2,B3;
  halfThreeHalfVectorCoupling(imode(),inpart.mass(),decay[0]->mass(),decay[1]->mass(),
			      A1,A2,A3,B1,B2,B3);
  Energy msum(inpart.mass()+decay[0]->mass());
  Complex lS,rS,lV,rV,left,right;
  // incoming particle
  if(inpart.id()>0) {
    lS=(A3-B3);rS=(A3+B3);
    lV=(A2-B2);rV=(A2+B2);
    left=(A1-B1);right=(A1+B1);
  }
  // incoming anti-particle
  else {
    lS=conj(A3+B3);rS=conj(A3-B3);
    lV=-conj(A2-B2);rV=-conj(A2+B2);
    left=conj(A1+B1);right=conj(A1-B1);
  }
  // compute the matrix element
  vector<unsigned int> ispin(3);
  LorentzVector<complex<Energy> > svec;
  Complex prod;
  complex<Energy> scalar;
  LorentzSpinor<SqrtEnergy> stemp;
  LorentzSpinorBar<SqrtEnergy> sbtemp;
  for(unsigned iya=0;iya<4;++iya) {
    ispin[1]=iya;
    // piece where the vector-spinor is dotted with the momentum of the
    // incoming fermion
    for(unsigned ixa=0;ixa<2;++ixa) {
      unsigned int ix(iya),iy(ixa);
      if(decay[0]->id()<0) swap(ix,iy);
      scalar = _inHalf[iy].generalScalar (_inHalfBar[ix],lS,rS);
      svec   = _inHalf[iy].generalCurrent(_inHalfBar[ix],lV,rV);
      ispin[0]=ixa;
      for(unsigned int iz=0;iz<3;++iz) {
	ispin[2]=iz;
	prod=_inVec[iz].dot(inpart.momentum())/msum;
	(*ME())(ispin) += (svec.dot(_inVec[iz])+prod*scalar)*
	  UnitRemoval::E/msum/inpart.mass();
      }
    }
    // the piece where the vector spinor is dotted with the polarization vector
    for(unsigned int iz=0;iz<3;++iz) {
      ispin[2]=iz;
      if(decay[0]->id()>0) sbtemp = _inThreeHalfBar[iya].dot(_inVec[iz]);
      else                 stemp  = _inThreeHalf[iya].dot(_inVec[iz]);
      for(unsigned int ixa=0;ixa<2;++ixa) {
	ispin[0]=ixa;
	if(decay[0]->id()>0) stemp  = _inHalf[ixa];
	else                 sbtemp = _inHalfBar[ixa];
	(*ME())(ispin) += stemp.generalScalar(sbtemp,left,right)/inpart.mass();
      }
    }
  }
  double output = (ME()->contract(_rho)).real();
  // test of the matrix element
//   Energy m1(inpart.mass()),m2(decay[0]->mass()),m3(decay[1]->mass());
//   Energy2 m12(m1*m1),m22(m2*m2),m32(m3*m3);
//   Energy Qp(sqrt(sqr(m1+m2)-sqr(m3))),Qm(sqrt(sqr(m1-m2)-sqr(m3)));
//   double r2(sqrt(2.)),r3(sqrt(3.));
//   Energy pcm(Kinematics::pstarTwoBodyDecay(m1,m2,m3));
//   complex<Energy> h1(-2.*Qp*A1),h2(2.*Qm*B1);
//   complex<Energy> h3(-2./r3*Qp*(A1-Qm*Qm/m2*A2/msum));
//   complex<Energy> h4( 2./r3*Qm*(B1-Qp*Qp/m2*B2/msum));
//   complex<Energy> h5(-2.*r2/r3/m2/m3*Qp*(0.5*(m12-m22-m32)*A1+0.5*Qm*Qm*(m1+m2)*A2/msum
// 					 +m12*pcm*pcm*A3/msum/msum));
//   complex<Energy> h6( 2.*r2/r3/m2/m3*Qm*(0.5*(m12-m22-m32)*B1-0.5*Qp*Qp*(m1-m2)*B2/msum
// 					 +m12*pcm*pcm*B3/msum/msum));
//   cout << "testing 1/2->3/2 1 " << inpart.id() << " "
//        << output << "   " 
//        << 0.25*(h1*conj(h1)+h2*conj(h2)+h3*conj(h3)+
// 		h4*conj(h4)+h5*conj(h5)+h6*conj(h6))/sqr(inpart.mass()) << "   " 
//        << 0.25*(h1*conj(h1)+h2*conj(h2)+h3*conj(h3)+
// 		h4*conj(h4)+h5*conj(h5)+h6*conj(h6))/sqr(inpart.mass())/output << endl;
  // return the answer
  return output;
}


// matrix element for the decay of a spin-3/2 fermion to a spin-1/2 fermion and
// a scalar meson
double Baryon1MesonDecayerBase::
threeHalfHalfScalar(const int,const Particle & inpart,
		    const ParticleVector & decay, MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin3Half,PDT::Spin1Half,PDT::Spin0)));
  // spinors etc for the decaying particle
  if(meopt==Initialize) {
    // spinors and rho
    if(inpart.id()>0) {
      RSSpinorWaveFunction   ::calculateWaveFunctions(_inThreeHalf,_rho,
						      const_ptr_cast<tPPtr>(&inpart),
						      incoming);
    }
    else {
      RSSpinorBarWaveFunction::calculateWaveFunctions(_inThreeHalfBar,_rho,
						      const_ptr_cast<tPPtr>(&inpart),
						      incoming);
    }
    // matrix element
  }
  // setup spin info when needed
  if(meopt==Terminate) {
    // for the decaying particle
    if(inpart.id()>0) {
      RSSpinorWaveFunction::
	constructSpinInfo(_inThreeHalf,const_ptr_cast<tPPtr>(&inpart),incoming,true);
      SpinorBarWaveFunction::constructSpinInfo(_inHalfBar,decay[0],outgoing,true);

    }
    else {
      RSSpinorBarWaveFunction::
	constructSpinInfo(_inThreeHalfBar,const_ptr_cast<tPPtr>(&inpart),incoming,true);
      SpinorWaveFunction::constructSpinInfo(_inHalf,decay[0],outgoing,true);
    }
    ScalarWaveFunction::constructSpinInfo(decay[1],outgoing,true);
    return 0.;
  }
  // spinors for the decay product
  if(inpart.id()>0) {
    SpinorBarWaveFunction::calculateWaveFunctions(_inHalfBar,decay[0],outgoing);
  }
  else {
    SpinorWaveFunction   ::calculateWaveFunctions(_inHalf,decay[0],outgoing);
  }
  LorentzPolarizationVector out=UnitRemoval::InvE*decay[0]->momentum();
  if(inpart.id()>0) {
    _inHalf.resize(_inThreeHalf.size());
    for(unsigned int ix=0;ix<_inThreeHalf.size();++ix)
      _inHalf[ix] = _inThreeHalf[ix].dot(out);
  }
  else {
    _inHalfBar.resize(_inThreeHalfBar.size());
    for(unsigned int ix=0;ix<_inThreeHalfBar.size();++ix)
      _inHalfBar[ix] = _inThreeHalfBar[ix].dot(out);
  }
  // get the couplings
  Complex A,B;
  Energy msum=inpart.mass()+decay[0]->mass();
  threeHalfHalfScalarCoupling(imode(),inpart.mass(),decay[0]->mass(),decay[1]->mass(),
			      A,B);
  Complex left,right;
  // incoming particle
  if(inpart.id()>0) {
    left=(A-B);
    right=(A+B);
  }
  // incoming anti-particle
  else {
    left=conj(A+B);
    right=conj(A-B);
  }
  // compute the matrix element
  vector<unsigned int> ispin(3,0);
  for(unsigned ixa=0;ixa<2;++ixa) {
    for(unsigned iya=0;iya<4;++iya) {
      unsigned int iy=iya,ix=ixa;
      if(decay[0]->id()<0) swap(ix,iy);
      ispin[0]=iya;
      ispin[1]=ixa;
      (*ME())(ispin) = _inHalf[iy].generalScalar(_inHalfBar[ix],left,right)*
	UnitRemoval::E/msum/inpart.mass();
    }
  }
  double output = (ME()->contract(_rho)).real();
  // test of the matrix element
//   Energy m1(inpart.mass()),m2(decay[0]->mass()),m3(decay[1]->mass());
//   Energy Qp(sqrt(sqr(m1+m2)-sqr(m3))),Qm(sqrt(sqr(m1-m2)-sqr(m3)));
//   double r23(sqrt(2./3.));
//   Energy pcm(Kinematics::pstarTwoBodyDecay(m1,m2,m3));
//   complex<Energy> h1(-2.*r23*pcm*Qm*B/(m1+m2)),h2( 2.*r23*pcm*Qp*A/(m1+m2));
//   cout << "testing 3/2->1/2 0 " << inpart.id() << " "
//        << output << "   " 
//        << 0.125*(h1*conj(h1)+h2*conj(h2))/sqr(inpart.mass()) << "   " 
//        << 0.125*(h1*conj(h1)+h2*conj(h2))/sqr(inpart.mass())/output << endl;
  // return the answer
  return output;
}

// matrix element for the decay of a spin-3/2 fermion to a spin-3/2 fermion and
// a scalar meson
double Baryon1MesonDecayerBase::threeHalfThreeHalfScalar(const int,
							 const Particle & inpart,
							 const ParticleVector & decay,
							 MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin3Half,PDT::Spin3Half,PDT::Spin0)));
  // spinors etc for the decaying particle
  if(meopt==Initialize) {
    // spinors and rho
    if(inpart.id()>0) {
      RSSpinorWaveFunction   ::calculateWaveFunctions(_inThreeHalf,_rho,
						      const_ptr_cast<tPPtr>(&inpart),
						      incoming);
    }
    else {
      RSSpinorBarWaveFunction::calculateWaveFunctions(_inThreeHalfBar,_rho,
						      const_ptr_cast<tPPtr>(&inpart),
						      incoming);
    }
    // matrix element
  }
  // setup spin info when needed
  if(meopt==Terminate) {
    // for the decaying particle
    if(inpart.id()>0) {
      RSSpinorWaveFunction::
	constructSpinInfo(_inThreeHalf,const_ptr_cast<tPPtr>(&inpart),incoming,true);
      RSSpinorBarWaveFunction::constructSpinInfo(_inThreeHalfBar,
						 decay[0],outgoing,true);
      
    }
    else {
      RSSpinorBarWaveFunction::
	constructSpinInfo(_inThreeHalfBar,const_ptr_cast<tPPtr>(&inpart),incoming,true);
      RSSpinorWaveFunction::constructSpinInfo(_inThreeHalf,
					      decay[0],outgoing,true);
    }
    ScalarWaveFunction::constructSpinInfo(decay[1],outgoing,true);
    return 0.;
  }
  // spinors for the decay product
  LorentzPolarizationVector in=UnitRemoval::InvE*inpart.momentum();
  if(inpart.id()>0) {
    RSSpinorBarWaveFunction::
      calculateWaveFunctions(_inThreeHalfBar,decay[0],outgoing);
  }
  else {
    RSSpinorWaveFunction::
      calculateWaveFunctions(_inThreeHalf,decay[0],outgoing);
  }
  _inHalf.resize(_inThreeHalf.size());
  _inHalfBar.resize(_inThreeHalfBar.size());
  for(unsigned int ix=0;ix<_inThreeHalf.size();++ix) {
    _inHalf[ix] = _inThreeHalf[ix].dot(in);
    _inHalfBar[ix] = _inThreeHalfBar[ix].dot(in);
  }
  // get the couplings
  Complex A1,B1,A2,B2;
  Energy msum(inpart.mass()+decay[0]->mass());
  threeHalfThreeHalfScalarCoupling(imode(),inpart.mass(),decay[0]->mass(),
				   decay[1]->mass(),A1,A2,B1,B2);
  Complex left1,right1,left2,right2;
  // incoming particle
  if(inpart.id()>0) {
    left1=(A1-B1); right1=(A1+B1);
    left2=(A2-B2); right2=(A2+B2);
  }
  // incoming anti-particle
  else {
    left1=(A1+B1); right1=(A1-B1);
    left2=(A2+B2); right2=(A2-B2);
  }
  // compute the matrix element
  vector<unsigned int> ispin(3,0);
  for(unsigned ixa=0;ixa<4;++ixa) {
    for(unsigned iya=0;iya<4;++iya) {
      unsigned int iy=iya,ix=ixa;
      if(decay[0]->id()<0) swap(ix,iy);
      ispin[0]=iya;
      ispin[1]=ixa;
      (*ME())(ispin)=(_inThreeHalf[iy].generalScalar(_inThreeHalfBar[ix],left1,right1)
		   +_inHalf[iy].generalScalar( _inHalfBar[ix],left2,right2)
		   *UnitRemoval::E2/sqr(msum))/inpart.mass();
    }
  }
  // return the answer
  return (ME()->contract(_rho)).real();
}

// matrix element for the decay of a spin-3/2 fermion to a spin-1/2 fermion and
// a vector meson

double Baryon1MesonDecayerBase::
threeHalfHalfVector(const int,const Particle & inpart,
		    const ParticleVector & decay,MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin3Half,PDT::Spin1Half,PDT::Spin1)));
  // check if the outgoing meson is really a photon
  bool photon=decay[1]->id()==ParticleID::gamma;
  // spinors etc for the decaying particle
  if(meopt==Initialize) {
    // spinors and rho
    if(inpart.id()>0) {
      RSSpinorWaveFunction   ::calculateWaveFunctions(_inThreeHalf,_rho,
						      const_ptr_cast<tPPtr>(&inpart),
						      incoming);
    }
    else {
      RSSpinorBarWaveFunction::calculateWaveFunctions(_inThreeHalfBar,_rho,
						      const_ptr_cast<tPPtr>(&inpart),
						      incoming);
    }
    // matrix element
  }

  // setup spin info when needed
  if(meopt==Terminate) {
    // for the decaying particle
    if(inpart.id()>0) {
      RSSpinorWaveFunction::
	constructSpinInfo(_inThreeHalf,const_ptr_cast<tPPtr>(&inpart),incoming,true);
      SpinorBarWaveFunction::constructSpinInfo(_inHalfBar,decay[0],outgoing,true);

    }
    else {
      RSSpinorBarWaveFunction::
	constructSpinInfo(_inThreeHalfBar,const_ptr_cast<tPPtr>(&inpart),incoming,true);
      SpinorWaveFunction::constructSpinInfo(_inHalf,decay[0],outgoing,true);
    }
    VectorWaveFunction::constructSpinInfo(_inVec,decay[1],outgoing,true,photon);
    return 0.;
  }
  // spinors for the decay product
  if(inpart.id()>0) {
    SpinorBarWaveFunction::calculateWaveFunctions(_inHalfBar,decay[0],outgoing);
  }
  else {
    SpinorWaveFunction   ::calculateWaveFunctions(_inHalf,decay[0],outgoing);
  }
  LorentzPolarizationVector out=UnitRemoval::InvE*decay[0]->momentum();
  if(inpart.id()>0) {
    _inHalf.resize(_inThreeHalf.size());
    for(unsigned int ix=0;ix<_inThreeHalf.size();++ix)
      _inHalf[ix] = _inThreeHalf[ix].dot(out);
  }
  else {
    _inHalfBar.resize(_inThreeHalfBar.size());
    for(unsigned int ix=0;ix<_inThreeHalfBar.size();++ix)
      _inHalfBar[ix] = _inThreeHalfBar[ix].dot(out);
  }
  ME()->zero();
  VectorWaveFunction::calculateWaveFunctions(_inVec,decay[1],outgoing,photon);
  // get the couplings
  Complex A1,A2,A3,B1,B2,B3,prod,meout;
  threeHalfHalfVectorCoupling(imode(),inpart.mass(),decay[0]->mass(),decay[1]->mass(),
			      A1,A2,A3,B1,B2,B3);
  Energy msum(inpart.mass()+decay[0]->mass());
  Complex lS,rS,lV,rV,left,right;
  // incoming particle
  if(inpart.id()>0) {
      lS=(A3-B3);rS=(A3+B3);
      lV=(A2-B2);rV=(A2+B2);
      left=(A1-B1);right=(A1+B1);
  }
  // incoming anti-particle
  else {
    lS=conj(A3+B3);rS=conj(A3-B3);
    lV=-conj(A2-B2);rV=-conj(A2+B2);
    left=conj(A1+B1);right=conj(A1-B1);
  }
  // compute the matrix element
  vector<unsigned int> ispin(3);
  LorentzVector<complex<Energy> > svec;
  LorentzSpinor<SqrtEnergy> stemp;
  LorentzSpinorBar<SqrtEnergy> sbtemp;
  complex<Energy> scalar;
  for(unsigned iya=0;iya<4;++iya) {
    ispin[0]=iya;
    for(unsigned ixa=0;ixa<2;++ixa) {
      unsigned int iy=iya,ix=ixa;
      if(decay[0]->id()<0) swap(ix,iy);
      scalar = _inHalf[iy].generalScalar( _inHalfBar[ix],lS,rS);
      svec   = _inHalf[iy].generalCurrent(_inHalfBar[ix],lV,rV);
      ispin[1]=ixa;
      for(unsigned int iz=0;iz<3;++iz) {
	ispin[2]=iz;
	prod=_inVec[iz].dot(decay[0]->momentum())/msum;
	(*ME())(ispin) += (svec.dot(_inVec[iz])+prod*scalar)*
	  UnitRemoval::E/msum/inpart.mass();
      }
    }
    // the piece where the vector spinor is dotted with the polarization vector
    for(unsigned iz=0;iz<3;++iz) {
      ispin[2]=iz;
      if(decay[0]->id()>0) stemp  = _inThreeHalf[iya].dot(_inVec[iz]);
      else                 sbtemp = _inThreeHalfBar[iya].dot(_inVec[iz]);
      for(unsigned int ixa=0;ixa<2;++ixa) {
	ispin[1]=ixa;
	if(decay[0]->id()>0) sbtemp = _inHalfBar[ixa];
	else                 stemp  = _inHalf[ixa];
	(*ME())(ispin) += stemp.generalScalar(sbtemp,left,right)/inpart.mass();
      }
    }
  }
  double output = (ME()->contract(_rho)).real();
  // testing code
//   Energy m1(inpart.mass()),m2(decay[0]->mass()),m3(decay[1]->mass());
//   Energy2 m12(m1*m1),m22(m2*m2),m32(m3*m3);
//   Energy Qp(sqrt(sqr(m1+m2)-sqr(m3))),Qm(sqrt(sqr(m1-m2)-sqr(m3)));
//   double r2(sqrt(2.)),r3(sqrt(3.));
//   Energy pcm(Kinematics::pstarTwoBodyDecay(m1,m2,m3));
//   complex<Energy> h1(-2.*Qp*A1),h2(2.*Qm*B1);
//   complex<Energy> h3(-2./r3*Qp*(A1-Qm*Qm/m1*A2/msum));
//   complex<Energy> h4( 2./r3*Qm*(B1-Qp*Qp/m1*B2/msum));
//   complex<Energy> h5(-2.*r2/r3/m1/m3*Qp*(0.5*(m22-m12-m32)*A1+0.5*Qm*Qm*(m1+m2)*A2/msum
//  					 +m12*pcm*pcm*A3/msum/msum));
//   complex<Energy> h6( 2.*r2/r3/m1/m3*Qm*(0.5*(m22-m12-m32)*B1-0.5*Qp*Qp*(m2-m1)*B2/msum
// 					 +m22*pcm*pcm*B3/msum/msum));
//   cout << "testing 3/2->1/2 1 " << inpart.id() << " "
//        << output << "   " 
//        << 0.25*(h1*conj(h1)+h2*conj(h2)+h3*conj(h3)+
// 		h4*conj(h4)+h5*conj(h5)+h6*conj(h6))/sqr(inpart.mass()) << "   " 
//        << 0.25*(h1*conj(h1)+h2*conj(h2)+h3*conj(h3)+
// 		h4*conj(h4)+h5*conj(h5)+h6*conj(h6))/sqr(inpart.mass())/output << endl;
  // return the answer
  return output;
}


void Baryon1MesonDecayerBase::halfHalfScalarCoupling(int,Energy,Energy,Energy,
						     Complex&,Complex&) const {
  throw DecayIntegratorError() << "Baryon1MesonDecayerBase::halfHalfScalarCoupling()"
			      << " called from base class this must be implemented"
			      << " in the inheriting class" << Exception::abortnow;
}

void Baryon1MesonDecayerBase::halfHalfVectorCoupling(int,Energy,Energy,Energy,
						     Complex&,Complex&,
						     Complex&,Complex&) const {
  throw DecayIntegratorError() << "Baryon1MesonDecayerBase::halfHalfVectorCoupling()" 
			       << " called from base class this must be implemented " 
			       << "in the inheriting class" << Exception::abortnow;
}

void Baryon1MesonDecayerBase::halfThreeHalfScalarCoupling(int,Energy,Energy,Energy,
							  Complex&,Complex&) const {
  throw DecayIntegratorError() << "Baryon1MesonDecayerBase::halfThreeHalfScalarCoupling"
			       << "() called from base class this must be implemented"
			       << " in the inheriting class" << Exception::abortnow;
}

void Baryon1MesonDecayerBase::halfThreeHalfVectorCoupling(int,Energy,Energy,Energy,
							  Complex&,Complex&,
							  Complex&,Complex&,
							  Complex&,Complex&) const {
  throw DecayIntegratorError() << "Baryon1MesonDecayerBase::halfThreeHalfVectorCoupling"
			       << "() called from base class this must be implemented "
			       << "in the inheriting class" << Exception::abortnow;
}

void Baryon1MesonDecayerBase::threeHalfHalfScalarCoupling(int,Energy,Energy,Energy,
							  Complex&,Complex&) const {
  throw DecayIntegratorError() << "Baryon1MesonDecayerBase::threeHalfHalfScalarCoupling"
			       << "() called from base class this must be implemented"
			       << " in the inheriting class" << Exception::abortnow;
}

void Baryon1MesonDecayerBase::threeHalfHalfVectorCoupling(int,Energy,Energy,Energy,
							  Complex&,Complex&,
							  Complex&,Complex&,
							  Complex&,Complex&) const {
  throw DecayIntegratorError() << "Baryon1MesonDecayerBase::threeHalfHalfVectorCoupling"
			       << "() called from base class this must be implemented "
			       << "in the inheriting class" << Exception::abortnow;
}

void Baryon1MesonDecayerBase::threeHalfThreeHalfScalarCoupling(int,Energy,Energy,
							       Energy,Complex&,
							       Complex&,Complex&,
							       Complex&) const {
  throw DecayIntegratorError() << "Baryon1MesonDecayerBase::threeHalfThreeHalfScalar"
			       << "Coupling() called from base class this must be "
			       << "implemented in the inheriting class" 
			       << Exception::abortnow;
}

bool Baryon1MesonDecayerBase::twoBodyMEcode(const DecayMode & dm,int & mecode,
					    double & coupling) const {
  coupling=1.;
  unsigned int inspin(dm.parent()->iSpin()),outspin,outmes;
  ParticleMSet::const_iterator pit(dm.products().begin());
  bool order; 
  if((**pit).iSpin()%2==0) {
    order=true;
    outspin=(**pit).iSpin();
    ++pit;outmes=(**pit).iSpin();
  }
  else {
    order=false;
    outmes=(**pit).iSpin();++pit;
    outspin=(**pit).iSpin();
  }
  mecode=-1;
  if(inspin==2) {
    if(outspin==2){if(outmes==1){mecode=101;}else{mecode=102;}}
    else if(outspin==4){if(outmes==1){mecode=103;}else{mecode=104;}}
  }
  else if(inspin==4) {
    if(outspin==2){if(outmes==1){mecode=105;}else{mecode=106;}}
    else if(outspin==4){if(outmes==1){mecode=107;}else{mecode=108;}}
  }
  return order;
}

void Baryon1MesonDecayerBase::dataBaseOutput(ofstream & os,bool header) const {
  DecayIntegrator::dataBaseOutput(os,header);
}
