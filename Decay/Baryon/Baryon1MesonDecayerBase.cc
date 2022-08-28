// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Baryon1MesonDecayerBase class.
//

#include "Baryon1MesonDecayerBase.h"
#include "ThePEG/Utilities/DescribeClass.h"
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
#include "ThePEG/Helicity/HelicityFunctions.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractNoPIOClass<Baryon1MesonDecayerBase,DecayIntegrator>
describeHerwigBaryon1MesonDecayerBase("Herwig::Baryon1MesonDecayerBase", "HwBaryonDecay.so");

void Baryon1MesonDecayerBase::Init() {

  static ClassDocumentation<Baryon1MesonDecayerBase> documentation
    ("The Baryon1MesonDecayerBase class is the base class for"
     " the decays of the baryons to a baryon and a pseudoscalar or vector meson.");

}

void Baryon1MesonDecayerBase::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  // for the decaying particle
  if(part.id()>0) {
    // incoming particle
    if(part.dataPtr()->iSpin()==PDT::Spin1Half) {
      SpinorWaveFunction::
	constructSpinInfo(_inHalf,const_ptr_cast<tPPtr>(&part),incoming,true);
    }
    else if(part.dataPtr()->iSpin()==PDT::Spin3Half) {
      RSSpinorWaveFunction::
	constructSpinInfo(_inThreeHalf,const_ptr_cast<tPPtr>(&part),incoming,true);
    }
    else
      assert(false);
    // outgoing fermion
    if(decay[0]->dataPtr()->iSpin()==PDT::Spin1Half) {
      SpinorBarWaveFunction::constructSpinInfo(_inHalfBar,decay[0],outgoing,true);
    }
    else if(decay[0]->dataPtr()->iSpin()==PDT::Spin3Half) {
      RSSpinorBarWaveFunction::constructSpinInfo(_inThreeHalfBar,
						 decay[0],outgoing,true);
    }
    else
      assert(false);
  }
  else {
    // incoming particle
    if(part.dataPtr()->iSpin()==PDT::Spin1Half) {
      SpinorBarWaveFunction::
	constructSpinInfo(_inHalfBar,const_ptr_cast<tPPtr>(&part),incoming,true);
    }
    else if(part.dataPtr()->iSpin()==PDT::Spin3Half) {
      RSSpinorBarWaveFunction::
	constructSpinInfo(_inThreeHalfBar,const_ptr_cast<tPPtr>(&part),incoming,true);
    }
    else
      assert(false);
    // outgoing fermion
    if(decay[0]->dataPtr()->iSpin()==PDT::Spin1Half) {
      SpinorWaveFunction::constructSpinInfo(_inHalf,decay[0],outgoing,true);
    }
    else if(decay[0]->dataPtr()->iSpin()==PDT::Spin3Half) {
      RSSpinorWaveFunction::constructSpinInfo(_inThreeHalf,
					      decay[0],outgoing,true);
    }
    else
      assert(false);
  }
  // outgoing meson
  if(decay[1]->dataPtr()->iSpin()==PDT::Spin0) {
    ScalarWaveFunction::constructSpinInfo(decay[1],outgoing,true);
  }
  else if(decay[1]->dataPtr()->iSpin()==PDT::Spin1) {
    VectorWaveFunction::constructSpinInfo(_inVec,decay[1],outgoing,true,
					  decay[1]->id()==ParticleID::gamma);
  }
  else
    assert(false);
}

// return the matrix element squared for a given mode and phase-space channel
// (inherited from DecayIntegrator and implemented here)
double Baryon1MesonDecayerBase::me2(const int ichan,const Particle & part,
				    const tPDVector & outgoing,
				    const vector<Lorentz5Momentum> & momenta,
				    MEOption meopt) const {
  double me(0.);
  // decide which matrix element we are doing
  // incoming spin-1/2 particle
  if(part.dataPtr()->iSpin()==2) {
    // decay to spin-1/2 particle
    if(outgoing[0]->iSpin()==2) {
      // scalar meson
      if(outgoing[1]->iSpin()==1)
   	me=halfHalfScalar(ichan,part,outgoing,momenta,meopt);
      // vector meson
      else if(outgoing[1]->iSpin()==3)
   	me=halfHalfVector(ichan,part,outgoing,momenta,meopt);
      else
   	throw DecayIntegratorError() << "Unknown outgoing meson spin in "
				      << "Baryon1MesonDecayerBase::me2()" 
				      << Exception::abortnow;
    }
    // decay to spin-3/2 particle
    else if(outgoing[0]->iSpin()==4) {
      // scalar meson
      if(outgoing[1]->iSpin()==1)
   	me=halfThreeHalfScalar(ichan,part,outgoing,momenta,meopt);
      // vector meson
      else if(outgoing[1]->iSpin()==3)
   	me=halfThreeHalfVector(ichan,part,outgoing,momenta,meopt);
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
  else if(part.dataPtr()->iSpin()==4) {
    // decay to spin-1/2 particle
    if(outgoing[0]->iSpin()==2) {
      // scalar meson
      if(outgoing[1]->iSpin()==1)
   	me=threeHalfHalfScalar(ichan,part,outgoing,momenta,meopt);
      // vector meson
      else if(outgoing[1]->iSpin()==3)
   	me=threeHalfHalfVector(ichan,part,outgoing,momenta,meopt);
      else
   	throw DecayIntegratorError() << "Unknown outgoing meson spin in "
				      << "Baryon1MesonDecayerBase::me2()" 
				      << Exception::abortnow;
    }
    // decay to spin-3/2 particle
    else if(outgoing[0]->iSpin()==4) {
      // scalar meson
      if(outgoing[1]->iSpin()==1)
   	me=threeHalfThreeHalfScalar(ichan,part,outgoing,momenta,meopt);
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
halfHalfScalar(const int,const Particle & part, const tPDVector & outgoing,
	       const vector<Lorentz5Momentum> & momenta,
	       MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin0)));
  // spinors etc for the decaying particle
  if(meopt==Initialize) {
    // spinors and rho
    if(part.id()>0)
      SpinorWaveFunction   ::calculateWaveFunctions(_inHalf,_rho,
						    const_ptr_cast<tPPtr>(&part),
						    incoming);
    else
      SpinorBarWaveFunction::calculateWaveFunctions(_inHalfBar,_rho,
						    const_ptr_cast<tPPtr>(&part),
						    incoming);
    // matrix element
  }
  // spinors for the decay product
  if(part.id()>0) {
    _inHalfBar.resize(2);
    for(unsigned int ix=0;ix<2;++ix)
      _inHalfBar[ix] = HelicityFunctions::dimensionedSpinorBar(-momenta[0],ix,Helicity::outgoing);
  }
  else {
    _inHalf.resize(2);
    for(unsigned int ix=0;ix<2;++ix)
      _inHalf[ix] = HelicityFunctions::dimensionedSpinor   (-momenta[0],ix,Helicity::outgoing);
  }
  // get the couplings
  Complex A,B;
  halfHalfScalarCoupling(imode(),part.mass(),momenta[0].mass(),momenta[1].mass(),A,B);
  Complex left,right;
  // coupling for an incoming particle
  if(part.id()>0) {
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
  // Complex output(0.);
  for(ix=0;ix<2;++ix) {
    for(iy=0;iy<2;++iy) {
      if(outgoing[0]->id()>0){ispin[0]=iy;ispin[1]=ix;}
      else{ispin[0]=ix;ispin[1]=iy;}
      (*ME())(ispin)=Complex(_inHalf[iy].generalScalar(_inHalfBar[ix],left,right)/part.mass());
    }
  }
  // test of the matrix element
  // Energy m1(part.mass()),m2(momenta[0].mass()),m3(momenta[1].mass());
  // Energy Qp(sqrt(sqr(m1+m2)-sqr(m3))),Qm(sqrt(sqr(m1-m2)-sqr(m3)));
  // Complex h1(2.*Qp*A/part.mass()),h2(-2.*Qm*B/part.mass());
  // cout << "testing 1/2->1/2 0 " 
  //      << 0.5*output << "   " 
  //      << 0.25*(h1*conj(h1)+h2*conj(h2)) << "   " 
  //      << 0.5*(h1*conj(h1)+h2*conj(h2))/output << endl;
  // cout << "testing alpha " << 
  //   (norm(0.5*(h1+h2))-norm(0.5*(h1-h2)))/
  //   (norm(0.5*(h1+h2))+norm(0.5*(h1-h2))) << "\n";
  // store the matrix element
  return (ME()->contract(_rho)).real();
}

// matrix element for the decay of a spin-1/2 fermion to a spin-1/2 fermion and
// a vector meson
double Baryon1MesonDecayerBase::
halfHalfVector(const int,const Particle & part, const tPDVector & outgoing,
	       const vector<Lorentz5Momentum> & momenta, MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin1)));
  // check if the outgoing meson is really a photon
  bool photon=outgoing[1]->id()==ParticleID::gamma;
  // spinors etc for the decaying particle
  if(meopt==Initialize) {
    // spinors and rho
    if(part.id()>0)
      SpinorWaveFunction   ::calculateWaveFunctions(_inHalf,_rho,
						    const_ptr_cast<tPPtr>(&part),
						    incoming);
    else
      SpinorBarWaveFunction::calculateWaveFunctions(_inHalfBar,_rho,
						    const_ptr_cast<tPPtr>(&part),
						    incoming);
    // matrix element
  }
  // spinors for the decay product
  if(part.id()>0) {
    _inHalfBar.resize(2);
    for(unsigned int ix=0;ix<2;++ix)
      _inHalfBar[ix] = HelicityFunctions::dimensionedSpinorBar(-momenta[0],ix,Helicity::outgoing);
  }
  else {
    _inHalf.resize(2);
    for(unsigned int ix=0;ix<2;++ix)
      _inHalf[ix] = HelicityFunctions::dimensionedSpinor   (-momenta[0],ix,Helicity::outgoing);
  }
  _inVec.resize(3);
  for(unsigned int ix=0;ix<3;++ix) {
    if(photon && ix==1) continue;
    _inVec[ix] = HelicityFunctions::polarizationVector(-momenta[1],ix,Helicity::outgoing);
  }
  // get the couplings
  Complex A1,A2,B1,B2;
  halfHalfVectorCoupling(imode(),part.mass(),momenta[0].mass(),momenta[1].mass(),
			 A1,A2,B1,B2);
  Complex lS,rS,lV,rV;
  complex<Energy> scalar;
  // couplings for an incoming particle
  if(part.id()>0) {
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
  Energy msum(part.mass()+momenta[0].mass());
  vector<unsigned int> ispin(3);
  LorentzVector<complex<Energy> > svec;
  Complex prod;
  // Complex output(0.);
  unsigned int ix,iy;
  for(ix=0;ix<2;++ix) {
    for(iy=0;iy<2;++iy) {
      // scalar like piece
      scalar = _inHalf[iy].generalScalar(_inHalfBar[ix],lS,rS);
      // vector like piece
      svec   = _inHalf[iy].generalCurrent(_inHalfBar[ix],lV,rV);
      if(outgoing[0]->id()>0) {
	ispin[0] = iy;
	ispin[1] = ix;
      }
      else {
	ispin[0] = ix;
	ispin[1] = iy;
      }
      for(ispin[2]=0;ispin[2]<3;++ispin[2]) {
	ispin[2]=ispin[2];
	prod=_inVec[ispin[2]].dot(part.momentum())/msum;
	(*ME())(ispin)=(svec.dot(_inVec[ispin[2]])+prod*scalar)/part.mass();
 	// output += norm((*ME())(ispin));
      }
    }
  }
  // test of the matrix element
  // Energy m1(part.mass()),m2(momenta[0].mass()),m3(momenta[1].mass());
  // Energy Qp(sqrt(sqr(m1+m2)-sqr(m3))),Qm(sqrt(sqr(m1-m2)-sqr(m3)));
  // double r2(sqrt(2.));
  // Energy pcm(Kinematics::pstarTwoBodyDecay(m1,m2,m3));
  // if(m3!=ZERO) {
  //   Complex h1(2.*r2*Qp*B1/part.mass()),h2(-2.*r2*Qm*A1/part.mass()),
  //     h3(2./m3*(Qp*(m1-m2)*B1-Qm*m1*B2*pcm/(m1+m2))/part.mass()),
  //     h4(2./m3*(Qm*(m1+m2)*A1+Qp*m1*A2*pcm/(m1+m2))/part.mass());
  //   generator()->log() << "testing 1/2->1/2 1 " 
  // 		       << 0.5*output << "   " 
  // 		       << 0.25*(h1*conj(h1)+h2*conj(h2)+h3*conj(h3)+h4*conj(h4)) << "   " 
  // 		       << 0.50*(h1*conj(h1)+h2*conj(h2)+h3*conj(h3)+h4*conj(h4))/output 
  // 		       << "\n";
  //   generator()->log() << "alpha = " << 2.*(norm(h3)+norm(h4))/(norm(h1)+norm(h2))-1.
  // 		       << "\n";
  //   generator()->log() << "masses " << m1/GeV << " " << m2/GeV << " " << m3/GeV << "\n";
  // }
  // return the answer
  return (ME()->contract(_rho)).real();
}

// matrix element for the decay of a spin-1/2 fermion to a spin-3/2 fermion and
// a scalar meson
double Baryon1MesonDecayerBase::
halfThreeHalfScalar(const int,const Particle & part, const tPDVector & outgoing,
		    const vector<Lorentz5Momentum> & momenta, MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half,PDT::Spin3Half,PDT::Spin0)));
  // spinors etc for the decaying particle
  if(meopt==Initialize) {
    // spinors and rho
    if(part.id()>0)
      SpinorWaveFunction   ::calculateWaveFunctions(_inHalf,_rho,
						    const_ptr_cast<tPPtr>(&part),
						    incoming);
    else
      SpinorBarWaveFunction::calculateWaveFunctions(_inHalfBar,_rho,
						    const_ptr_cast<tPPtr>(&part),
						    incoming);
    // matrix element
  }
  // spinors for the decay product
  LorentzPolarizationVector in=UnitRemoval::InvE*part.momentum();
  if(part.id()>0) {
    RSSpinorBarWaveFunction swave(momenta[0],outgoing[0],Helicity::outgoing);
    _inThreeHalfBar.resize(4);
    _inHalfBar.resize(4);
    for(unsigned int ihel=0;ihel<4;++ihel) {
      swave.reset(ihel);
      _inThreeHalfBar[ihel]=swave.dimensionedWf();
      _inHalfBar[ihel] = _inThreeHalfBar[ihel].dot(in);
    }
  }
  else {
    RSSpinorWaveFunction swave(momenta[0],outgoing[0],Helicity::outgoing);
    _inThreeHalf.resize(4);
    _inHalf.resize(4);
    for(unsigned int ihel=0;ihel<4;++ihel) {
      swave.reset(ihel);
      _inThreeHalf[ihel]=swave.dimensionedWf();
      _inHalf[ihel] = _inThreeHalf[ihel].dot(in);
    }
  }
  // get the couplings
  Complex A,B,left,right;
  Energy msum(part.mass()+momenta[0].mass());
  halfThreeHalfScalarCoupling(imode(),part.mass(),momenta[0].mass(),momenta[1].mass(),
			      A,B);
  // incoming particle
  if(part.id()>0) {
    left=(A-B);
    right=(A+B);
  }
  // incoming anti-particle
  else {
    left=conj(A+B);
    right=conj(A-B);
  }
  vector<unsigned int> ispin(3,0);
  for(unsigned ixa=0;ixa<2;++ixa) {
    for(unsigned int iya=0;iya<4;++iya) {
      unsigned int ix(iya),iy(ixa);
      if(outgoing[0]->id()<0) swap(ix,iy);
      ispin[0]=ixa;
      ispin[1]=iya;
      complex<double> value = _inHalf[iy].generalScalar(_inHalfBar[ix],left,right)
	*UnitRemoval::E/part.mass()/msum;
      (*ME())(ispin) = value;
    }
  }
  double output = (ME()->contract(_rho)).real();
  // test of the matrix element
   // Energy m1(part.mass()),m2(momenta[0].mass()),m3(momenta[1].mass());
   // Energy Qp(sqrt(sqr(m1+m2)-sqr(m3))),Qm(sqrt(sqr(m1-m2)-sqr(m3)));
   // double r23(sqrt(2./3.));
   // Energy pcm(Kinematics::pstarTwoBodyDecay(m1,m2,m3));
   // complex<Energy> h1(-2.*r23*pcm*m1/m2*Qm*B/(m1+m2)),h2( 2.*r23*pcm*m1/m2*Qp*A/(m1+m2));
   // cout << "testing 1/2->3/2 0 " << part.id() << " "
   //      << output << "   " 
   //      << 0.25*(h1*conj(h1)+h2*conj(h2))/sqr(part.mass()) << "   " 
   //      << 0.25*(h1*conj(h1)+h2*conj(h2))/sqr(part.mass())/output << endl;
  // return the answer
  return output;
}

// matrix element for the decay of a spin-1/2 fermion to a spin-3/2 fermion and
// a vector meson
double Baryon1MesonDecayerBase::
halfThreeHalfVector(const int,const Particle & part, const tPDVector & outgoing,
		    const vector<Lorentz5Momentum> & momenta, MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half,PDT::Spin3Half,PDT::Spin1)));
  // check if the outgoing meson is really a photon
  bool photon=outgoing[1]->id()==ParticleID::gamma;
  // spinors etc for the decaying particle
  if(meopt==Initialize) {
    // spinors and rho
    if(part.id()>0)
      SpinorWaveFunction   ::calculateWaveFunctions(_inHalf,_rho,
						    const_ptr_cast<tPPtr>(&part),
						    incoming);
    else
      SpinorBarWaveFunction::calculateWaveFunctions(_inHalfBar,_rho,
						    const_ptr_cast<tPPtr>(&part),
						    incoming);
    // matrix element
  }
  LorentzPolarizationVector in=UnitRemoval::InvE*part.momentum();
  // wavefunctions for outgoing fermion
  if(part.id()>0) {
    RSSpinorBarWaveFunction swave(momenta[0],outgoing[0],Helicity::outgoing);
    _inThreeHalfBar.resize(4);
    _inHalfBar.resize(4);
    for(unsigned int ihel=0;ihel<4;++ihel) {
      swave.reset(ihel);
      _inThreeHalfBar[ihel]=swave.dimensionedWf();
      _inHalfBar[ihel] = _inThreeHalfBar[ihel].dot(in);
    }
  }
  else {
    RSSpinorWaveFunction swave(momenta[0],outgoing[0],Helicity::outgoing);
    _inThreeHalf.resize(4);
    _inHalf.resize(4);
    for(unsigned int ihel=0;ihel<4;++ihel) {
      swave.reset(ihel);
      _inThreeHalf[ihel]=swave.dimensionedWf();
      _inHalf[ihel] = _inThreeHalf[ihel].dot(in);
    }
  }
  ME()->zero();
  _inVec.resize(3);
  for(unsigned int ix=0;ix<3;++ix) {
    if(photon && ix==1) continue;
    _inVec[ix] = HelicityFunctions::polarizationVector(-momenta[1],ix,Helicity::outgoing);
  }
  // get the couplings
  Complex A1,A2,A3,B1,B2,B3;
  halfThreeHalfVectorCoupling(imode(),part.mass(),momenta[0].mass(),momenta[1].mass(),
			      A1,A2,A3,B1,B2,B3);
  Energy msum(part.mass()+momenta[0].mass());
  Complex lS,rS,lV,rV,left,right;
  // incoming particle
  if(part.id()>0) {
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
      if(outgoing[0]->id()<0) swap(ix,iy);
      scalar = _inHalf[iy].generalScalar (_inHalfBar[ix],lS,rS);
      svec   = _inHalf[iy].generalCurrent(_inHalfBar[ix],lV,rV);
      ispin[0]=ixa;
      for(unsigned int iz=0;iz<3;++iz) {
	ispin[2]=iz;
	prod=_inVec[iz].dot(part.momentum())/msum;
	(*ME())(ispin) += (svec.dot(_inVec[iz])+prod*scalar)*
	  UnitRemoval::E/msum/part.mass();
      }
    }
    // the piece where the vector spinor is dotted with the polarization vector
    for(unsigned int iz=0;iz<3;++iz) {
      ispin[2]=iz;
      if(outgoing[0]->id()>0) sbtemp = _inThreeHalfBar[iya].dot(_inVec[iz]);
      else                 stemp  = _inThreeHalf[iya].dot(_inVec[iz]);
      for(unsigned int ixa=0;ixa<2;++ixa) {
	ispin[0]=ixa;
	if(outgoing[0]->id()>0) stemp  = _inHalf[ixa];
	else                    sbtemp = _inHalfBar[ixa];
	(*ME())(ispin) += Complex(stemp.generalScalar(sbtemp,left,right)/part.mass());
      }
    }
  }
  double output = (ME()->contract(_rho)).real();
  // test of the matrix element
  // Energy m1(part.mass()),m2(momenta[0].mass()),m3(momenta[1].mass());
  // Energy2 m12(m1*m1),m22(m2*m2),m32(m3*m3);
  // Energy Qp(sqrt(sqr(m1+m2)-sqr(m3))),Qm(sqrt(sqr(m1-m2)-sqr(m3)));
  // double r2(sqrt(2.)),r3(sqrt(3.));
  // Energy pcm(Kinematics::pstarTwoBodyDecay(m1,m2,m3));
  // complex<Energy> h1(-2.*Qp*A1),h2(2.*Qm*B1);
  // complex<Energy> h3(-2./r3*Qp*(A1-Qm*Qm/m2*A2/msum));
  // complex<Energy> h4( 2./r3*Qm*(B1-Qp*Qp/m2*B2/msum));
  // complex<Energy> h5(-2.*r2/r3/m2/m3*Qp*(0.5*(m12-m22-m32)*A1+0.5*Qm*Qm*(m1+m2)*A2/msum
  // 					 +m12*pcm*pcm*A3/msum/msum));
  // complex<Energy> h6( 2.*r2/r3/m2/m3*Qm*(0.5*(m12-m22-m32)*B1-0.5*Qp*Qp*(m1-m2)*B2/msum
  // 					 +m12*pcm*pcm*B3/msum/msum));
  // cout << "testing 1/2->3/2 1 " << part.id() << " "
  //      << output << "   " 
  //      << 0.25*(h1*conj(h1)+h2*conj(h2)+h3*conj(h3)+
  // 		h4*conj(h4)+h5*conj(h5)+h6*conj(h6))/sqr(part.mass()) << "   " 
  //      << 0.25*(h1*conj(h1)+h2*conj(h2)+h3*conj(h3)+
  // 		h4*conj(h4)+h5*conj(h5)+h6*conj(h6))/sqr(part.mass())/output << endl;
  // return the answer
  return output;
}


// matrix element for the decay of a spin-3/2 fermion to a spin-1/2 fermion and
// a scalar meson
double Baryon1MesonDecayerBase::
threeHalfHalfScalar(const int,const Particle & part, const tPDVector & outgoing,
		    const vector<Lorentz5Momentum> & momenta, MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin3Half,PDT::Spin1Half,PDT::Spin0)));
  // spinors etc for the decaying particle
  if(meopt==Initialize) {
    // spinors and rho
    if(part.id()>0) {
      RSSpinorWaveFunction   ::calculateWaveFunctions(_inThreeHalf,_rho,
						      const_ptr_cast<tPPtr>(&part),
						      incoming);
    }
    else {
      RSSpinorBarWaveFunction::calculateWaveFunctions(_inThreeHalfBar,_rho,
						      const_ptr_cast<tPPtr>(&part),
						      incoming);
    }
  }
  // spinors for the decay product
  if(part.id()>0) {
    _inHalfBar.resize(2);
    for(unsigned int ix=0;ix<2;++ix)
      _inHalfBar[ix] = HelicityFunctions::dimensionedSpinorBar(-momenta[0],ix,Helicity::outgoing);
  }
  else {
    _inHalf.resize(2);
    for(unsigned int ix=0;ix<2;++ix)
      _inHalf[ix] = HelicityFunctions::dimensionedSpinor   (-momenta[0],ix,Helicity::outgoing);
  }
  LorentzPolarizationVector out=UnitRemoval::InvE*momenta[0];
  if(part.id()>0) {
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
  Energy msum=part.mass()+momenta[0].mass();
  threeHalfHalfScalarCoupling(imode(),part.mass(),momenta[0].mass(),momenta[1].mass(),
			      A,B);
  Complex left,right;
  // incoming particle
  if(part.id()>0) {
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
      if(outgoing[0]->id()<0) swap(ix,iy);
      ispin[0]=iya;
      ispin[1]=ixa;
      (*ME())(ispin) = Complex(_inHalf[iy].generalScalar(_inHalfBar[ix],left,right)*
			       UnitRemoval::E/msum/part.mass());
    }
  }
  // test of the matrix element
  // double test = (ME()->contract(RhoDMatrix(PDT::Spin3Half))).real();
  // Energy m1(part.mass()),m2(momenta[0].mass()),m3(momenta[1].mass());
  // Energy Qp(sqrt(sqr(m1+m2)-sqr(m3))),Qm(sqrt(sqr(m1-m2)-sqr(m3)));
  // double r23(sqrt(2./3.));
  // Energy pcm(Kinematics::pstarTwoBodyDecay(m1,m2,m3));
  // complex<Energy> h1(-2.*r23*pcm*Qm*B/(m1+m2)),h2( 2.*r23*pcm*Qp*A/(m1+m2));
  // generator()->log() << "testing 3/2->1/2 0 " << part.id() << " "
  // 		     << test << "   " 
  // 		     << 0.125*(h1*conj(h1)+h2*conj(h2))/sqr(part.mass()) << "   " 
  // 		     << 0.125*(h1*conj(h1)+h2*conj(h2))/sqr(part.mass())/test << endl;
  // return the answer
  return (ME()->contract(_rho)).real();;
}

// matrix element for the decay of a spin-3/2 fermion to a spin-3/2 fermion and
// a scalar meson
double Baryon1MesonDecayerBase::threeHalfThreeHalfScalar(const int, const Particle & part,
							 const tPDVector & outgoing,
							 const vector<Lorentz5Momentum> & momenta,
							 MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin3Half,PDT::Spin3Half,PDT::Spin0)));
  // spinors etc for the decaying particle
  if(meopt==Initialize) {
    // spinors and rho
    if(part.id()>0) {
      RSSpinorWaveFunction   ::calculateWaveFunctions(_inThreeHalf,_rho,
						      const_ptr_cast<tPPtr>(&part),
						      incoming);
    }
    else {
      RSSpinorBarWaveFunction::calculateWaveFunctions(_inThreeHalfBar,_rho,
						      const_ptr_cast<tPPtr>(&part),
						      incoming);
    }
    // matrix element
  }
  // spinors for the decay product
  // wavefunctions for outgoing fermion
  if(part.id()>0) {
    RSSpinorBarWaveFunction swave(momenta[0],outgoing[0],Helicity::outgoing);
    _inThreeHalfBar.resize(4);
    for(unsigned int ihel=0;ihel<4;++ihel) {
      swave.reset(ihel);
      _inThreeHalfBar[ihel]=swave.dimensionedWf();
    }
  }
  else {
    RSSpinorWaveFunction swave(momenta[0],outgoing[0],Helicity::outgoing);
    _inThreeHalf.resize(4);
    for(unsigned int ihel=0;ihel<4;++ihel) {
      swave.reset(ihel);
      _inThreeHalf[ihel]=swave.dimensionedWf();
    }
  }
  LorentzPolarizationVector in = UnitRemoval::InvE*part.momentum();
  _inHalf.resize(_inThreeHalf.size());
  _inHalfBar.resize(_inThreeHalfBar.size());
  for(unsigned int ix=0;ix<_inThreeHalf.size();++ix) {
    _inHalf[ix] = _inThreeHalf[ix].dot(in);
    _inHalfBar[ix] = _inThreeHalfBar[ix].dot(in);
  }
  // get the couplings
  Complex A1,B1,A2,B2;
  Energy msum(part.mass()+momenta[0].mass());
  threeHalfThreeHalfScalarCoupling(imode(),part.mass(),momenta[0].mass(),
				   momenta[1].mass(),A1,A2,B1,B2);
  Complex left1,right1,left2,right2;
  // incoming particle
  if(part.id()>0) {
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
      if(outgoing[0]->id()<0) swap(ix,iy);
      ispin[0]=iya;
      ispin[1]=ixa;
      (*ME())(ispin)=Complex((_inThreeHalf[iy].generalScalar(_inThreeHalfBar[ix],left1,right1)
			      +_inHalf[iy].generalScalar( _inHalfBar[ix],left2,right2)
			      *UnitRemoval::E2/sqr(msum))/part.mass());
    }
  }
  // return the answer
  return (ME()->contract(_rho)).real();
}

// matrix element for the decay of a spin-3/2 fermion to a spin-1/2 fermion and
// a vector meson

double Baryon1MesonDecayerBase::
threeHalfHalfVector(const int,const Particle & part, const tPDVector & outgoing,
		    const vector<Lorentz5Momentum> & momenta, MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin3Half,PDT::Spin1Half,PDT::Spin1)));
  // check if the outgoing meson is really a photon
  bool photon=outgoing[1]->id()==ParticleID::gamma;
  // spinors etc for the decaying particle
  if(meopt==Initialize) {
    // spinors and rho
    if(part.id()>0) {
      RSSpinorWaveFunction   ::calculateWaveFunctions(_inThreeHalf,_rho,
						      const_ptr_cast<tPPtr>(&part),
						      incoming);
    }
    else {
      RSSpinorBarWaveFunction::calculateWaveFunctions(_inThreeHalfBar,_rho,
						      const_ptr_cast<tPPtr>(&part),
						      incoming);
    }
  }
  if(part.id()>0) {
    _inHalfBar.resize(2);
    for(unsigned int ix=0;ix<2;++ix)
      _inHalfBar[ix] = HelicityFunctions::dimensionedSpinorBar(-momenta[0],ix,Helicity::outgoing);
  }
  else {
    _inHalf.resize(2);
    for(unsigned int ix=0;ix<2;++ix)
      _inHalf[ix] = HelicityFunctions::dimensionedSpinor   (-momenta[0],ix,Helicity::outgoing);
  }
  LorentzPolarizationVector out=UnitRemoval::InvE*momenta[0];
  if(part.id()>0) {
    _inHalf.resize(_inThreeHalf.size());
    for(unsigned int ix=0;ix<_inThreeHalf.size();++ix)
      _inHalf[ix] = _inThreeHalf[ix].dot(out);
  }
  else {
    _inHalfBar.resize(_inThreeHalfBar.size());
    for(unsigned int ix=0;ix<_inThreeHalfBar.size();++ix)
      _inHalfBar[ix] = _inThreeHalfBar[ix].dot(out);
  }
  // wavefunctions for outgoing fermion
  ME()->zero();
  _inVec.resize(3);
  for(unsigned int ix=0;ix<3;++ix) {
    if(photon && ix==1) continue;
    _inVec[ix] = HelicityFunctions::polarizationVector(-momenta[1],ix,Helicity::outgoing);
  }
  // get the couplings
  Complex A1,A2,A3,B1,B2,B3,prod;
  threeHalfHalfVectorCoupling(imode(),part.mass(),momenta[0].mass(),momenta[1].mass(),
			      A1,A2,A3,B1,B2,B3);
  Energy msum(part.mass()+momenta[0].mass());
  Complex lS,rS,lV,rV,left,right;
  // incoming particle
  if(part.id()>0) {
    lS    = (A3-B3);
    rS    = (A3+B3);
    lV    = (A2-B2);
    rV    = (A2+B2);
    left  = (A1-B1);
    right = (A1+B1);
  }
  // incoming anti-particle
  else {
    lS    = conj(A3+B3);
    rS    = conj(A3-B3);
    lV    =-conj(A2-B2);
    rV    =-conj(A2+B2);
    left  = conj(A1+B1);
    right = conj(A1-B1);
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
      if(outgoing[0]->id()<0) swap(ix,iy);
      scalar = _inHalf[iy].generalScalar( _inHalfBar[ix],lS,rS);
      svec   = _inHalf[iy].generalCurrent(_inHalfBar[ix],lV,rV);
      ispin[1]=ixa;
      for(unsigned int iz=0;iz<3;++iz) {
	ispin[2]=iz;
	prod=_inVec[iz].dot(momenta[0])/msum;
	(*ME())(ispin) += (svec.dot(_inVec[iz])+prod*scalar)*
	  UnitRemoval::E/msum/part.mass();
      }
    }
    // the piece where the vector spinor is dotted with the polarization vector
    for(unsigned iz=0;iz<3;++iz) {
      ispin[2]=iz;
      if(outgoing[0]->id()>0) stemp  = _inThreeHalf[iya].dot(_inVec[iz]);
      else                    sbtemp = _inThreeHalfBar[iya].dot(_inVec[iz]);
      for(unsigned int ixa=0;ixa<2;++ixa) {
	ispin[1]=ixa;
	if(outgoing[0]->id()>0) sbtemp = _inHalfBar[ixa];
	else                 stemp  = _inHalf[ixa];
	(*ME())(ispin) += Complex(stemp.generalScalar(sbtemp,left,right)/part.mass());
      }
    }
  }
  double output = (ME()->contract(_rho)).real();
  // testing code
   // Energy m1(part.mass()),m2(momenta[0].mass()),m3(momenta[1].mass());
   // Energy2 m12(m1*m1),m22(m2*m2),m32(m3*m3);
   // Energy Qp(sqrt(sqr(m1+m2)-sqr(m3))),Qm(sqrt(sqr(m1-m2)-sqr(m3)));
   // double r2(sqrt(2.)),r3(sqrt(3.));
   // Energy pcm(Kinematics::pstarTwoBodyDecay(m1,m2,m3));
   // complex<Energy> h1(-2.*Qp*A1),h2(2.*Qm*B1);
   // complex<Energy> h3(-2./r3*Qp*(A1-Qm*Qm/m1*A2/msum));
   // complex<Energy> h4( 2./r3*Qm*(B1-Qp*Qp/m1*B2/msum));
   // complex<Energy> h5(ZERO),h6(ZERO);
   // if(m3!=ZERO) {
   //   h5 = (-2.*r2/r3/m1/m3*Qp*(0.5*(m22-m12-m32)*A1+0.5*Qm*Qm*(m1+m2)*A2/msum
   // 					 +m12*pcm*pcm*A3/msum/msum));
   //   h6 = ( 2.*r2/r3/m1/m3*Qm*(0.5*(m22-m12-m32)*B1-0.5*Qp*Qp*(m2-m1)*B2/msum
   // 			       +m22*pcm*pcm*B3/msum/msum));
   // }
   // cout << "testing 3/2->1/2 1 " << part.id() << " "
   //      << output << "   " 
   //      << 0.25*(h1*conj(h1)+h2*conj(h2)+h3*conj(h3)+
   // 		h4*conj(h4)+h5*conj(h5)+h6*conj(h6))/sqr(part.mass()) << "   " 
   //      << 0.25*(h1*conj(h1)+h2*conj(h2)+h3*conj(h3)+
   // 		h4*conj(h4)+h5*conj(h5)+h6*conj(h6))/sqr(part.mass())/output << endl;
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
