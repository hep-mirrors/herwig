// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the OmegaPionsCurrent class.
//

#include "OmegaPionsCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/epsilon.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"

using namespace Herwig;
using namespace ThePEG::Helicity;
using namespace ThePEG::Helicity;

void OmegaPionsCurrent::persistentOutput(PersistentOStream & os) const {
  os << _fpi << _gfact << _mrho << _ma1 << _localparameters << _fa << _cfact << _mpi 
     << _Bcoup << _Dcoup;
}

void OmegaPionsCurrent::persistentInput(PersistentIStream & is, int) {
  is >> _fpi >> _gfact >> _mrho >> _ma1 >> _localparameters >> _fa >> _cfact >> _mpi
     >> _Bcoup >> _Dcoup;
}

ClassDescription<OmegaPionsCurrent> OmegaPionsCurrent::initOmegaPionsCurrent;
// Definition of the static class description member.

void OmegaPionsCurrent::Init() {

  static ClassDocumentation<OmegaPionsCurrent> documentation
    ("The OmegaPionsCurrent class implements the currents for omega and up to"
     "three pions");

  static Parameter<OmegaPionsCurrent,Energy> interfaceFpi
    ("Fpi",
     "The pion decay constant (sqrt(2) times PDG)",
     &OmegaPionsCurrent::_fpi, MeV, 186.*MeV, 0.0*MeV, 200.0*MeV,
     false, false, true);

  static Parameter<OmegaPionsCurrent,double> interfaceg
    ("g",
     "The g coupling",
     &OmegaPionsCurrent::_gfact, 0.35, 0.0, 10.0,
     false, false, true);

  static Parameter<OmegaPionsCurrent,Energy> interfaceRhoMass
    ("RhoMass",
     "The mass of the rho",
     &OmegaPionsCurrent::_mrho, MeV, 769.9*MeV, 0.0*MeV, 1000.0*MeV,
     false, false, true);

  static Parameter<OmegaPionsCurrent,Energy> interfacea1Mass
    ("a1Mass",
     "The mass of the a1",
     &OmegaPionsCurrent::_ma1, MeV, 1230.*MeV, 0.0*MeV, 2000.0*MeV,
     false, false, true);

  static Switch<OmegaPionsCurrent,bool> interfaceLocalParameters
    ("LocalParameters",
     "Whether or not to use local values of the rho and a1 masses",
     &OmegaPionsCurrent::_localparameters, true, false, false);
  static SwitchOption interfaceLocalParametersLocal
    (interfaceLocalParameters,
     "Local",
     "Use local values",
     true);
  static SwitchOption interfaceLocalParametersParticleData
    (interfaceLocalParameters,
     "ParticleData",
     "Use values from the particle data objects",
     false);
}

PDVector OmegaPionsCurrent::particles(int icharge, unsigned int imode, int , int ) {
  PDVector output(1,getParticleData(ParticleID::omega));
  tPDPtr pim,pip,pi0(getParticleData(ParticleID::pi0));
  if(icharge==3) {
    pim=getParticleData(ParticleID::piminus);
    pip=getParticleData(ParticleID::piplus);
  }
  else {
    pim=getParticleData(ParticleID::piplus);
    pip=getParticleData(ParticleID::piminus);
  }
  output.push_back(pip);
  if(imode==1) {
    output.push_back(pi0);
  }
  else if(imode==2){
    output.push_back(pi0);
    output.push_back(pi0);
  }
  else if(imode==3) {
    output.push_back(pip);
    output.push_back(pim);
  }
  return output;
}

bool OmegaPionsCurrent::accept(vector<int> id) {
  unsigned int npi0(0),npim(0),npip(0),nomega(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
      if(id[ix]==ParticleID::omega)        ++nomega;
      else if(id[ix]==ParticleID::piplus)  ++npip;
      else if(id[ix]==ParticleID::piminus) ++npim;
      else if(id[ix]==ParticleID::pi0)     ++npi0;
    }
  if(id.size()-nomega-npip-npim-npi0!=0||nomega!=1) return false;
  int net=npip-npim;net=abs(net);
  return net==1&&id.size()<=4;
}

unsigned int OmegaPionsCurrent::decayMode(vector<int> id) {
  unsigned int npip(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(abs(id[ix])==ParticleID::piplus) ++npip;
  }
  unsigned int imode(id.size()-2);
  if(npip==3) ++imode;
  return imode;
}

bool OmegaPionsCurrent::createMode(int icharge,unsigned int imode,
				   DecayPhaseSpaceModePtr mode,
				   unsigned int iloc,unsigned int ires,
				   DecayPhaseSpaceChannelPtr phase,Energy upp) {
  if(abs(icharge)!=3) return false;
  Energy min(getParticleData(ParticleID::omega)->massMin()+_mpi);
  if(imode==1)     min+=_mpi;
  else if(imode>1) min+=2.*_mpi;
  if(min<upp) return false;
  tPDPtr W(getParticleData(icharge/3*ParticleID::Wplus));
  tPDPtr rhop,rhom,rho0(getParticleData(ParticleID::rho0)),a1,
    a10(getParticleData(ParticleID::a_10));
  if(icharge==3) {
    rhom=getParticleData(ParticleID::rhoplus);
    rhop=getParticleData(ParticleID::rhominus);
    a1  =getParticleData(ParticleID::a_1plus);
  }
  else {
    rhop=getParticleData(ParticleID::rhoplus);
    rhom=getParticleData(ParticleID::rhominus);
    a1  =getParticleData(ParticleID::a_1minus);
  }
  DecayPhaseSpaceChannelPtr newchannel;
  if(imode==0) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
    newchannel->addIntermediate(W   ,0,0.0,iloc,iloc+1);
    mode->addChannel(newchannel);
  }
  else if(imode==1) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
    newchannel->addIntermediate(a1  ,0,0.0,-ires-1,iloc+2);
    newchannel->addIntermediate(rhom,0,0.0,iloc   ,iloc+1);
    mode->addChannel(newchannel);
    newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
    newchannel->addIntermediate(a1  ,0,0.0,-ires-1,iloc+1);
    newchannel->addIntermediate(rho0,0,0.0,iloc   ,iloc+2);
    mode->addChannel(newchannel);
  }
  else if(imode==2) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
    newchannel->addIntermediate(W   ,0,0.0,-ires-1,iloc+2);
    newchannel->addIntermediate(a1  ,0,0.0,-ires-2,iloc+1);
    newchannel->addIntermediate(rho0,0,0.0,iloc   ,iloc+3);
    mode->addChannel(newchannel);
    newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
    newchannel->addIntermediate(W   ,0,0.0,-ires-1,iloc+3);
    newchannel->addIntermediate(a1  ,0,0.0,-ires-2,iloc+1);
    newchannel->addIntermediate(rho0,0,0.0,iloc   ,iloc+2);
    mode->addChannel(newchannel);
    newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
    newchannel->addIntermediate(W   ,0,0.0,-ires-1,iloc+2);
    newchannel->addIntermediate(a1  ,0,0.0,-ires-2,iloc+3);
    newchannel->addIntermediate(rhom,0,0.0,iloc   ,iloc+1);
    mode->addChannel(newchannel);
    newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
    newchannel->addIntermediate(W   ,0,0.0,-ires-1,iloc+3);
    newchannel->addIntermediate(a1  ,0,0.0,-ires-2,iloc+2);
    newchannel->addIntermediate(rhom,0,0.0,iloc   ,iloc+1);
    mode->addChannel(newchannel);
  }
  else if(imode==3) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
    newchannel->addIntermediate(W   ,0,0.0,-ires-1,iloc+1);
    newchannel->addIntermediate(a10 ,0,0.0,-ires-2,iloc+2);
    newchannel->addIntermediate(rhop,0,0.0,iloc   ,iloc+3);
    mode->addChannel(newchannel);
    newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
    newchannel->addIntermediate(W   ,0,0.0,-ires-1,iloc+2);
    newchannel->addIntermediate(a10 ,0,0.0,-ires-2,iloc+1);
    newchannel->addIntermediate(rhop,0,0.0,iloc   ,iloc+3);
    mode->addChannel(newchannel);
    newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
    newchannel->addIntermediate(W   ,0,0.0,-ires-1,iloc+1);
    newchannel->addIntermediate(a10 ,0,0.0,-ires-2,iloc+3);
    newchannel->addIntermediate(rhom,0,0.0,iloc   ,iloc+2);
    mode->addChannel(newchannel);
    newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
    newchannel->addIntermediate(W   ,0,0.0,-ires-1,iloc+2);
    newchannel->addIntermediate(a10 ,0,0.0,-ires-2,iloc+3);
    newchannel->addIntermediate(rhom,0,0.0,iloc   ,iloc+1);
    mode->addChannel(newchannel);
  }
  else {
    throw Exception() << "Unknown mode in OmegaPionsCurrent::createMode() " 
		      << Exception::abortnow;
  }
  // reset resonance parameters for the integration
  mode->resetIntermediate(a1,_ma1,a1Width(_ma1));
  mode->resetIntermediate(rho0,_mrho,rhoWidth(_mrho));
  mode->resetIntermediate(rhop,_mrho,rhoWidth(_mrho));
  mode->resetIntermediate(rhom,_mrho,rhoWidth(_mrho));
  return true;
}

vector<LorentzPolarizationVector> 
OmegaPionsCurrent::current(bool vertex,const int imode,const int ichan,Energy & scale, 
			   const ParticleVector & decay) const {
  vector<LorentzPolarizationVector> temp;
  Lorentz5Momentum ptotal(decay[0]->momentum());
  // create the wavefunctions
  VectorWaveFunction(temp,decay[0],outgoing,true,false,vertex);
  unsigned int ix;
  for(ix=1;ix<decay.size();++ix) {
    ScalarWaveFunction(decay[ix],outgoing,true,vertex);
    ptotal+=decay[ix]->momentum();
  }
  ptotal.rescaleMass();
  scale=ptotal.mass();
  // calculate the current
  Energy2 q2(ptotal.m2());
  // current for omega pi
  if(imode==0) {
    for(ix=0;ix<3;++ix) {
      temp[ix] = 3./sqr(pi)/_gfact/_fpi*rhoBreitWigner(0,q2)*
	epsilon(decay[0]->momentum(),temp[ix],decay[1]->momentum());
    }
  }
  else if(imode==1) {
    // prefactor for direct diagrams
    double pre1(6.*scale/sqr(pi*_fpi)/_gfact*(1.-2.*_cfact/_gfact));
    double pre2(3.*scale/_fpi/_fa/sqr(pi*_gfact));
    // propagator for the a_1
    Complex prop(a1BreitWigner(0,q2));
    // propagators for the rho
    Energy2 
      k2 ((ptotal-decay[2]->momentum()).m2()),
      kp2((ptotal-decay[1]->momentum()).m2());
    Complex rho1(rhoBreitWigner(1, k2)); 
    Complex rho2(rhoBreitWigner(1,kp2));
    if(imode==0){rho2=0.;}
    else if(imode==1){rho1=0.;}
    // vectors and dot products
    LorentzPolarizationVector vec1,vec2;
    Energy2 dot1(ptotal*decay[2]->momentum()),dot2(ptotal*decay[1]->momentum());
    // coupling factors
    Energy A1(Acoupling(q2,k2)),A2(Acoupling(q2,kp2));
    Complex prod1,prod2;
    for(ix=0;ix<3;++ix) {
      vec1=epsilon(temp[ix],decay[0]->momentum(),decay[2]->momentum());
      vec2=epsilon(temp[ix],decay[0]->momentum(),decay[1]->momentum());
      prod1=vec1*ptotal;
      prod2=vec2*ptotal;
      temp[ix] =
	// first diagram 
	pre1*(2.*(2.*_cfact/_gfact+prop)*prod2/q2*ptotal-vec2+vec1)+
	// second diagram
	2.*pre2*prop*(
		      prod2*(A1*rho1+A2*rho2-_Bcoup*(dot2*rho2+dot1*rho1))/q2*ptotal
		      +_Bcoup*prod2*(rho2*decay[1]->momentum()+
				     rho1*decay[2]->momentum())+
		      A2*rho2*vec1-A1*rho1*vec2);
    }
  }
  else if(imode==2) {
    Complex pre(3.*q2/(sqr(pi)*_gfact*_fpi)*rhoBreitWigner(0,q2));
    LorentzPolarizationVector vec1,vec2,vec3;
    for(ix=0;ix<3;++ix) {
      vec1=epsilon(decay[2]->momentum(),decay[0]->momentum(),temp[ix]);
      vec2=epsilon(decay[3]->momentum(),decay[0]->momentum(),temp[ix]);
      vec3=epsilon(decay[1]->momentum(),decay[0]->momentum(),temp[ix]);
      temp[ix]=
	pre*(
	     ffunction(vec3,q2,ptotal,decay[2]->momentum(),decay[3]->momentum())
	     +ffunction(vec3,q2,ptotal,decay[3]->momentum(),decay[2]->momentum())
	     -ffunction(vec1,q2,ptotal,decay[3]->momentum(),decay[1]->momentum())
	     -ffunction(vec2,q2,ptotal,decay[2]->momentum(),decay[1]->momentum())
	     );
    }
  }
  else if(imode==3) {
    Complex pre(3.*q2/(sqr(pi)*_gfact*_fpi)*rhoBreitWigner(0,q2));
    LorentzPolarizationVector vec1,vec2,vec3;
    for(ix=0;ix<3;++ix) {
      vec1=epsilon(decay[1]->momentum(),decay[0]->momentum(),temp[ix]);
      vec2=epsilon(decay[3]->momentum(),decay[0]->momentum(),temp[ix]);
      vec3=epsilon(decay[2]->momentum(),decay[0]->momentum(),temp[ix]);
      temp[ix]=
	pre*(
	     ffunction(vec3,q2,ptotal,decay[1]->momentum(),decay[3]->momentum())
	     +ffunction(vec1,q2,ptotal,decay[2]->momentum(),decay[3]->momentum())
	     -ffunction(vec2,q2,ptotal,decay[1]->momentum(),decay[2]->momentum())
	     -ffunction(vec2,q2,ptotal,decay[2]->momentum(),decay[1]->momentum())
	     );
    }
  }
  // return the answer
  return temp;
}

void OmegaPionsCurrent::dataBaseOutput(ofstream & output,bool header,bool create) const {
  if(header) {
    output << "update decayers set parameters=\"";
  }
  if(create) {
    output << "create Herwig++::OmegaPionsCurrent " << fullName() << " \n";
  }
  output << "set " << fullName() << ":Fpi " << _fpi/MeV << "\n";
  output << "set " << fullName() << ":g " <<_gfact << "\n";
  output << "set " << fullName() << ":RhoMass " << _mrho/MeV << "\n";
  output << "set " << fullName() << ":a1Mass " << _ma1/MeV << "\n";
  output << "set " << fullName() << ":LocalParameters " << _localparameters << "\n";
  WeakDecayCurrent::dataBaseOutput(output,false,false);
  if(header) {
    output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
  }
}
  
LorentzPolarizationVector 
OmegaPionsCurrent::ffunction(LorentzPolarizationVector & alpha,Energy q2,
			     const Lorentz5Momentum & q,const Lorentz5Momentum & p1,
			     const Lorentz5Momentum & p2) const {
  // momentum and scale for the a_1
  Lorentz5Momentum k1(q-p1);
  Energy2 k12(k1.m2());
  Energy2 r2((k1-p2).m2());
  // rho breit wigner
  Complex prop(rhoBreitWigner(1,r2));
  // a_1 propagator and scale dependent couplings
  Complex BWa1(a1BreitWigner(1,k12)),Aq(Acoupling(k12,q2)),Ar(Acoupling(k12,r2));
  //momentum dot products
  Energy2 k1p2(k1*p2),k1p1(k1*p1),p1p2(p1*p2),qp2(q*p2),qp1(q*p1);
  // functions for the vector
  // a_1 contributions
  Complex f(BWa1*Aq*Ar);
  Complex f11(BWa1*Ar*_Bcoup);
  Complex f12(BWa1*(_Dcoup*(Ar+Aq)+(k1p2-k1p1)*_Bcoup*_Dcoup
		    +p1p2*sqr(_Bcoup)-k12*sqr(_Dcoup))
	      -BWa1/sqr(_ma1)*(-Aq+k1p2*_Bcoup+k12*_Dcoup)*(Ar+k1p2*_Bcoup-k12*_Dcoup));
  Complex f21(0.);
  Complex f22(BWa1*Aq*_Bcoup);
  // F1 and F2 functions needed for some diagrams
  Complex F2(sqr(2./_fpi)*(2.*sqr(_cfact/_gfact)
			   -0.75/_gfact*sqr((1.-2.*_cfact/_gfact)/pi)
			   -1./sqr(2.*pi*_gfact)*(1.-2.*_cfact/_gfact)));
  Complex f1pre(2./_gfact);
  complex<InvEnergy2> f1fact((sqr(1.-2.*_cfact/_gfact)-sqr(2.*pi*_cfact))/sqr(pi*_fpi));
  Complex F1k1p2(f1pre*(1.-k1p2*f1fact)),F1k1p1(f1pre*(1.+k1p1*f1fact));
  // first piece
  f12 += 
    0.5*sqr(_fa/(pi*_gfact*sqr(_ma1)))*
    (-Aq+k1p1*_Bcoup+k12*_Dcoup)*(Ar+k1p2*_Bcoup-k12*_Dcoup)
    +0.5/_gfact/_fpi/sqr(pi*_ma1)*_fa*(1.-2.*_cfact/_gfact)*(-Aq+k1p1*_Bcoup+k12*_Dcoup)*
    (2.*F1k1p2-k12*F2)
    -0.5/_gfact/_fpi/sqr(pi*_ma1)*_fa*(1.-2.*_cfact/_gfact)*( Ar+k1p2*_Bcoup-k12*_Dcoup)*
    (2.*F1k1p1-k12*F2)
    -0.5*sqr((1.-2.*_cfact/_gfact)/(pi*_fpi))*(2.*F1k1p2-k12*F2)*(2.*F1k1p2-k12*F2);
  // 1c diagrams
  f12 += -1./k12*(-4.*F1k1p1*F1k1p2+2.*k12*F2*F1k1p2+2.*k12*F2*F1k1p1-sqr(k12*F2));
  // 1a diagrams
  f   +=sqr(2./_fpi)*(sqr(_ma1)*(1.-0.5/sqr(pi*_gfact))-sqr(_mrho)
		      -2.*p1p2*(-2.*sqr(_cfact/_gfact)
				+sqr((1.-2.*_cfact/_gfact)/(2.*pi*_gfact))));
  f11 +=-2./sqr(pi*_gfact*_fpi)*(1.-2.*_cfact/_gfact);
  f12 +=-8./sqr(_fpi)*(-2.*sqr(_cfact/_gfact)
		       +0.75/sqr(pi*_gfact)*sqr(1.-2.*_cfact/_gfact)
		       +1.50/sqr(pi*_gfact)*(1.-2.*_cfact/_gfact));
  f21 += 8./sqr(_fpi)*(-sqr(2.*_cfact/_gfact)+0.75/sqr(pi*_gfact)*(1.-2.*_cfact/_gfact));
  f22 +=-2./sqr(pi*_gfact*_fpi)*(1.-2.*_cfact/_gfact);
  // 1d diagrams
  Complex rho(rhoBreitWigner(1,(p1+p2).m2())),
    rfact((1.*p1p2/sqr(pi*_fpi)*(sqr(1.-2.*_cfact/_gfact)-sqr(2.*pi*_cfact))));
  f   +=  8./sqr(_gfact)*rho*(qp2-qp1)*rfact;
  f12 += 16./sqr(_gfact)*rho*rfact;
  f21 -= 16./sqr(_gfact)*rho*rfact;
  // calculate the current
  Complex dot1(alpha*p1),dot2(alpha*p2);
  return prop*LorentzPolarizationVector(f*alpha+(f11*dot1+f21*dot2)*p1
					+(f12*dot1+f22*dot2)*p2);
}

