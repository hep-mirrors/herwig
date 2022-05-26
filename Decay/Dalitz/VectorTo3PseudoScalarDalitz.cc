// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorTo3PseudoScalarDalitz class.
//

#include "VectorTo3PseudoScalarDalitz.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/epsilon.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/PDT/DecayMode.h"

using namespace Herwig;

IBPtr VectorTo3PseudoScalarDalitz::clone() const {
  return new_ptr(*this);
}

IBPtr VectorTo3PseudoScalarDalitz::fullclone() const {
  return new_ptr(*this);
}

void VectorTo3PseudoScalarDalitz::persistentOutput(PersistentOStream & os) const {
}

void VectorTo3PseudoScalarDalitz::persistentInput(PersistentIStream & is, int) {
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<VectorTo3PseudoScalarDalitz,DalitzBase>
describeHerwigVectorTo3PseudoScalarDalitz("Herwig::VectorTo3PseudoScalarDalitz", "HwDalitzDecay.so");

void VectorTo3PseudoScalarDalitz::Init() {

  static ClassDocumentation<VectorTo3PseudoScalarDalitz> documentation
    ("The VectorTo3PseudoScalarDalitz class provides a base class "
     "for the decay of vector mesons to 3 pseudoscalar mesons");

}

void VectorTo3PseudoScalarDalitz::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  VectorWaveFunction::constructSpinInfo(vectors_,const_ptr_cast<tPPtr>(&part),
					incoming,true,false);
  // set up the spin information for the decay products
  for(unsigned int ix=0;ix<3;++ix)
    ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
}

double VectorTo3PseudoScalarDalitz::me2(const int ichan, const Particle & part,
				    const tPDVector & ,
				    const vector<Lorentz5Momentum> & momenta,
				    MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1,PDT::Spin0,PDT::Spin0,PDT::Spin0)));
  useMe();
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(vectors_,rho_,
  						const_ptr_cast<tPPtr>(&part),
  						incoming,false);
  }
  // compute the matrix element
  // // work out the prefactor
  // complex<InvEnergy2> pre(ZERO);
  // Complex resfact,ii(0.,1.);
  // if(ichan<0){pre=_ccoupling[imode()][3];}
  // Energy pcm;
  // // work out the direct invariant masses needed
  // Energy mrho0(sqrt(momenta[1].m2(momenta[2])));
  // Energy mrhop(sqrt(momenta[1].m2(momenta[0])));
  // Energy mrhom(sqrt(momenta[2].m2(momenta[0])));
  // // contribution of the resonances
  // int ichannow(-3);
  // for(unsigned int ix=0;ix<3;++ix) {
  //   ichannow+=3;
  //   if((ix==0 && _rho1wgt[imode()]>0.) || (ix==1 && _rho2wgt[imode()]>0.) ||
  //      (ix==2 && _rho3wgt[imode()]>0.)) {
  //     if(ichan<0) {
  // 	// rho0 contribution
  // 	pcm = Kinematics::pstarTwoBodyDecay(mrho0,_mpic,_mpic);
  // 	resfact = _rhomass2[imode()][ix]/
  // 	  (mrho0*mrho0-_rhomass2[imode()][ix]
  // 	   +ii*pcm*pcm*pcm*_rho0const[imode()][ix]/mrho0);
  // 	// rho+ contribution
  // 	pcm = Kinematics::pstarTwoBodyDecay(mrhop,_mpic,_mpi0);
  // 	resfact+= _rhomass2[imode()][ix]/
  // 	  (mrhop*mrhop-_rhomass2[imode()][ix]
  // 	   +ii*pcm*pcm*pcm*_rhocconst[imode()][ix]/mrhop);
  // 	// rho- contribution
  // 	pcm = Kinematics::pstarTwoBodyDecay(mrhom,_mpic,_mpi0);
  // 	resfact+= _rhomass2[imode()][ix]/
  // 	  (mrhom*mrhom-_rhomass2[imode()][ix]
  // 	   +ii*pcm*pcm*pcm*_rhocconst[imode()][ix]/mrhom);
  // 	// add the contribution
  //     }
  //     else if(ichan==ichannow) {
  // 	pcm = Kinematics::pstarTwoBodyDecay(mrho0,_mpic,_mpic);
  // 	resfact = _rhomass2[imode()][ix]/
  // 	  (mrho0*mrho0-_rhomass2[imode()][ix]
  // 	   +ii*pcm*pcm*pcm*_rho0const[imode()][ix]/mrho0);
  //     }
  //     else if(ichan==ichannow+1) {
  // 	pcm = Kinematics::pstarTwoBodyDecay(mrhop,_mpic,_mpi0);
  // 	resfact+= _rhomass2[imode()][ix]/
  // 	  (mrhop*mrhop-_rhomass2[imode()][ix]
  // 	   +ii*pcm*pcm*pcm*_rhocconst[imode()][ix]/mrhop);
  //     }
  //     else if(ichan==ichannow+2) {
  // 	pcm = Kinematics::pstarTwoBodyDecay(mrhom,_mpic,_mpi0);
  // 	resfact+= _rhomass2[imode()][ix]/
  // 	  (mrhom*mrhom-_rhomass2[imode()][ix]
  // 	   +ii*pcm*pcm*pcm*_rhocconst[imode()][ix]/mrhom);
  //     }
  //     pre += resfact * _ccoupling[imode()][ix];
  //     ichannow+=3;
  //   }
  // }
  // // polarization vector piece
  // LorentzPolarizationVector 
  //   scalar=_coupling[imode()]*pre*epsilon(momenta[0],
  // 					  momenta[1],
  // 					  momenta[2]);
  // // compute the matrix element
  // for(unsigned int ix=0;ix<3;++ix) {
  //   (*ME())(ix,0,0,0)=scalar.dot(vectors_[ix]);
  // }
  // // return the answer
  // return ME()->contract(_rho).real();
}

double VectorTo3PseudoScalarDalitz::
threeBodyMatrixElement(const int imode, const Energy2 q2,
		       const  Energy2 s3, const Energy2 s2, const Energy2 s1, const 
		       Energy , const Energy , const Energy ) const {
  // Lorentz5Momentum p1,p2,p3; Energy2 ee1,ee2,ee3;Energy pp1,pp2,pp3;
  // Energy q(sqrt(q2));
  // Energy2 mpi2c(_mpic*_mpic),mpi20(_mpi0*_mpi0);
  // p1.setE(0.5*(q2+mpi20-s1)/q); ee1=p1.e()*p1.e(); pp1=sqrt(ee1-mpi20);
  // p2.setE(0.5*(q2+mpi2c-s2)/q); ee2=p2.e()*p2.e(); pp2=sqrt(ee2-mpi2c);
  // p3.setE(0.5*(q2+mpi2c-s3)/q); ee3=p3.e()*p3.e(); pp3=sqrt(ee3-mpi2c);
  // // take momentum of 1 parallel to z axis
  // p1.setX(ZERO);p1.setY(ZERO);p1.setZ(pp1);
  // // construct 2 
  // double cos2(0.5*(ee1+ee2-ee3-mpi20)/pp1/pp2);
  // p2.setX(pp2*sqrt(1.-cos2*cos2)); p2.setY(ZERO); p2.setZ(-pp2*cos2);
  // // construct 3
  // double cos3(0.5*(ee1-ee2+ee3-mpi20)/pp1/pp3);
  // p3.setX(-pp3*sqrt(1.-cos3*cos3)); p3.setY(ZERO); p3.setZ(-pp3*cos3); 
  // // compute the prefactor
  // complex<InvEnergy2> pre(_ccoupling[imode][3]);
  // Complex resfact,ii(0.,1.);
  // // rho0 contribution
  // Energy pcm,mrho1(sqrt(s1)),mrho2(sqrt(s2)),mrho3(sqrt(s3));
  // for(unsigned int ix=0;ix<3;++ix) {
  //   // rho0 contribution
  //   pcm = Kinematics::pstarTwoBodyDecay(mrho1,_mpic,_mpic);
  //   resfact = _rhomass2[imode][ix]/(mrho1*mrho1-_rhomass2[imode][ix]
  // 				    +ii*pcm*pcm*pcm*_rho0const[imode][ix]/mrho1);
  //   // rho+ contribution
  //   pcm = Kinematics::pstarTwoBodyDecay(mrho2,_mpic,_mpi0);
  //   resfact+= _rhomass2[imode][ix]/(mrho2*mrho3-_rhomass2[imode][ix]
  // 				    +ii*pcm*pcm*pcm*_rhocconst[imode][ix]/mrho2);
  //   // rho- contribution
  //   pcm = Kinematics::pstarTwoBodyDecay(mrho3,_mpic,_mpi0);
  //   resfact+= _rhomass2[imode][ix]/(mrho3*mrho3-_rhomass2[imode][ix]
  // 				    +ii*pcm*pcm*pcm*_rhocconst[imode][ix]/mrho3);
  //   // add the contribution
  //   pre+=resfact *_ccoupling[imode][ix];
  // }
  // LorentzPolarizationVector current =
  //   _coupling[imode]*(pre*epsilon(p1,p2,p3));
  // Complex temp(current.dot(current.conjugate()));
  // return -temp.real()/3.;
} 

WidthCalculatorBasePtr 
VectorTo3PseudoScalarDalitz::threeBodyMEIntegrator(const DecayMode & dm) const {
  int imode=0;
  // construct the integrator
  vector<double> inweights;
  inweights.reserve(resonances().size());
  int iloc=-1;
  for(unsigned int ix=0;ix<resonances().size();++ix) {
    tPDPtr resonance = getParticleData(resonances()[ix].id);
    if(resonance) {
      ++iloc;
  //     mode->addChannel((PhaseSpaceChannel(mode),0,resonance,0,resonances()[ix].spectator+1,1,resonances()[ix].daughter1+1,1,resonances()[ix].daughter2+1));
  //     resetIntermediate(resonance,resonances()[ix].mass,abs(resonances()[ix].width));
  //     ++ix;
    }
    }
    
//   vector<double> inweights(3,1./3.);
//   vector<int> intype;intype.push_back(1);intype.push_back(2);intype.push_back(3);
//   Energy mrho(getParticleData(ParticleID::rhoplus)->mass());
//   Energy wrho(getParticleData(ParticleID::rhoplus)->width());
//   vector<Energy> inmass(3,mrho);
//   vector<Energy> inwidth(3,wrho);
//   vector<double> inpow(2,0.0);
//   //tcDecayIntegratorPtr decayer(this);
//   WidthCalculatorBasePtr output(
//     new_ptr(ThreeBodyAllOnCalculator<VectorTo3PseudoScalarDalitz>
// 	    (inweights,intype,inmass,inwidth,inpow,
// 	     *this,imode,_mpi0,_mpic,_mpic)));
//   return output;
}

void VectorTo3PseudoScalarDalitz::dataBaseOutput(ofstream & output,
					     bool header) const {
  if(header){output << "update decayers set parameters=\"";}
  // parameters for the DalitzBase base class
  DalitzBase::dataBaseOutput(output,false);
  if(header){output << "\n\" where BINARY ThePEGName=\"" 
   		    << fullName() << "\";" << endl;}
}
