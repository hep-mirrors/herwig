// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FourPionCzyzCurrent class.
//

#include "FourPionCzyzCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/Decay/ResonanceHelpers.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

FourPionCzyzCurrent::FourPionCzyzCurrent()
  : mpip_(140*MeV), mpi0_(140*MeV) {  // /**
  // Masses and widths of the particles
  // rho (PDG for most of current)
  rhoMasses_ = {0.7755*GeV,1.459*GeV,1.72*GeV};
  rhoWidths_ = {0.1494*GeV,0.4  *GeV,0.25*GeV};
  // fitted for F_rho
  rhoMasses_Frho_ = {0.7755*GeV,1.437*GeV,1.738*GeV,2.12*GeV};
  rhoWidths_Frho_ = {0.1494*GeV,0.520*GeV,0.450*GeV,0.30*GeV};
  // omega
  omegaMass_  = 0.78265*GeV;
  omegaWidth_ = 0.00849*GeV;
  // f0
  f0Mass_  = 1.35*GeV;
  f0Width_ = 0.2 *GeV;
  // a1  
  a1Mass_  = 1.23*GeV;
  a1Width_ = 0.2*GeV;
  // Coefficents for sum over \f$\rho\f$ resonances 
  beta_a1_    ={1.,-0.066,-0.021,-0.0043};
  beta_f0_    ={1.,7e4,-2.5e3,1.9e3};
  beta_omega_ ={1.,-0.33,-0.012,-0.0053};
  beta_B_     ={1.,-0.145};
  beta_bar_   ={1.,0.08,-0.0075};
  // couplings for the various termsz
  c_a1_    = -225./GeV2;
  c_f0_    =   64./GeV2;
  c_omega_ = -1.47/GeV;
  c_rho_   = -2.46;
  // meson meson meson couplings
  g_rho_pi_pi_   = 5.97;
  g_omega_pi_rho_= 42.3/GeV/GeV2/GeV2;
  g_rho_gamma_   = 0.1212*GeV2;
  addDecayMode(2,-1);
  addDecayMode(2,-1);
  addDecayMode(1,-1);
  addDecayMode(2,-2);
  addDecayMode(1,-1);
  addDecayMode(2,-2);
  setInitialModes(6);
}

IBPtr FourPionCzyzCurrent::clone() const {
  return new_ptr(*this);
}

IBPtr FourPionCzyzCurrent::fullclone() const {
  return new_ptr(*this);
}

void FourPionCzyzCurrent::persistentOutput(PersistentOStream & os) const {
  os << ounit(rhoMasses_,GeV) <<  ounit(rhoWidths_,GeV)
     << ounit(rhoMasses_Frho_,GeV) << ounit(rhoWidths_Frho_,GeV)
     << ounit(omegaMass_,GeV) << ounit(omegaWidth_,GeV)
     << ounit(f0Mass_,GeV) << ounit(f0Width_,GeV)
     << ounit(a1Mass_,GeV) << ounit(a1Width_,GeV)
     << beta_a1_ << beta_f0_ << beta_omega_ << beta_B_ << beta_bar_
     << ounit(c_a1_,1./GeV2) << ounit(c_f0_,1./GeV2) << ounit(c_omega_,1./GeV) << c_rho_
     << g_rho_pi_pi_ << ounit(g_omega_pi_rho_,1/GeV/GeV2/GeV2) << ounit(g_rho_gamma_,GeV2)
     << ounit(mpip_,GeV) << ounit(mpi0_,GeV);
}

void FourPionCzyzCurrent::persistentInput(PersistentIStream & is, int) {
  is >> iunit(rhoMasses_,GeV) >>  iunit(rhoWidths_,GeV)
     >> iunit(rhoMasses_Frho_,GeV) >> iunit(rhoWidths_Frho_,GeV)
     >> iunit(omegaMass_,GeV) >> iunit(omegaWidth_,GeV)
     >> iunit(f0Mass_,GeV) >> iunit(f0Width_,GeV)
     >> iunit(a1Mass_,GeV) >> iunit(a1Width_,GeV)
     >> beta_a1_ >> beta_f0_ >> beta_omega_ >> beta_B_ >> beta_bar_
     >> iunit(c_a1_,1./GeV2) >> iunit(c_f0_,1./GeV2) >> iunit(c_omega_,1./GeV) >> c_rho_
     >> g_rho_pi_pi_ >> iunit(g_omega_pi_rho_,1/GeV/GeV2/GeV2) >> iunit(g_rho_gamma_,GeV2)
     >> iunit(mpip_,GeV) >> iunit(mpi0_,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<FourPionCzyzCurrent,WeakDecayCurrent>
describeHerwigFourPionCzyzCurrent("Herwig::FourPionCzyzCurrent",
				  "HwWeakCurrents.so");

void FourPionCzyzCurrent::Init() {

  static ClassDocumentation<FourPionCzyzCurrent> documentation
    ("The FourPionCzyzCurrent class is designed to implement "
     "the four pion current for e+e- collisions from Phys.Rev. D77 (2008) 114005",
     "The current from \\cite{Czyz:2008kw} was used for four pions.",
     "\\bibitem{Czyz:2008kw}\n"
     "H.~Czyz, J.~H.~Kuhn and A.~Wapienik,\n"
     "%``Four-pion production in tau decays and e+e- annihilation: An Update,''\n"
     "Phys.\\ Rev.\\ D {\\bf 77} (2008) 114005\n"
     "doi:10.1103/PhysRevD.77.114005\n"
     "[arXiv:0804.0359 [hep-ph]].\n"
     "%%CITATION = doi:10.1103/PhysRevD.77.114005;%%\n"
     "%35 citations counted in INSPIRE as of 02 Aug 2018\n");

  static ParVector<FourPionCzyzCurrent,Energy> interfaceRhoMasses
    ("RhoMasses",
     "The masses of the rho mesons used by default in the current",
     &FourPionCzyzCurrent::rhoMasses_, GeV, -1, 0.7755*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);
  
  static ParVector<FourPionCzyzCurrent,Energy> interfaceRhoWidths
    ("RhoWidths",
     "The widths of the rho mesons used by default in the current",
     &FourPionCzyzCurrent::rhoWidths_, GeV, -1, 0.1494*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<FourPionCzyzCurrent,Energy> interfaceRhoMassesFrho
    ("RhoMassesFrho",
     "The masses of the rho mesons used in the F_rho piece",
     &FourPionCzyzCurrent::rhoMasses_Frho_, GeV, -1, 0.7755*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);
  
  static ParVector<FourPionCzyzCurrent,Energy> interfaceRhoWidthsFrho
    ("RhoWidthsFrho",
     "The widths of the rho mesons used in the F_rho piece",
     &FourPionCzyzCurrent::rhoWidths_Frho_, GeV, -1, 0.1494*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);
  
  static Parameter<FourPionCzyzCurrent,Energy> interfaceomegaMass
    ("omegaMass",
     "The mass of the omega meson",
     &FourPionCzyzCurrent::omegaMass_, GeV, 0.78265*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);
  
  static Parameter<FourPionCzyzCurrent,Energy> interfaceomegaWidth
    ("omegaWidth",
     "The width of the omega meson",
     &FourPionCzyzCurrent::omegaWidth_, GeV, 0.00849*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);
  
  static Parameter<FourPionCzyzCurrent,Energy> interfaceF0Mass
    ("f0Mass",
     "The mass of the f0 meson",
     &FourPionCzyzCurrent::f0Mass_, GeV, 1.35*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);
  
  static Parameter<FourPionCzyzCurrent,Energy> interfaceF0Width
    ("f0Width",
     "The width of the f0 meson",
     &FourPionCzyzCurrent::f0Width_, GeV, 0.2*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);
  
  static Parameter<FourPionCzyzCurrent,Energy> interfaceA1Mass
    ("a1Mass",
     "The mass of the a1 meson",
     &FourPionCzyzCurrent::a1Mass_, GeV, 1.23*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);
  
  static Parameter<FourPionCzyzCurrent,Energy> interfaceA1Width
    ("a1Width",
     "The width of the a1 meson",
     &FourPionCzyzCurrent::a1Width_, GeV, 0.2*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<FourPionCzyzCurrent,double> interfacebeta_a1
    ("beta_a1",
     "The coefficients for the sum over rho resonances in the a_1 term",
     &FourPionCzyzCurrent::beta_a1_, -1, 1.0, 0, 0,
     false, false, Interface::nolimits);

  static ParVector<FourPionCzyzCurrent,double> interfacebeta_f0
    ("beta_f0",
     "The coefficients for the sum over rho resonances in the f_0 term",
     &FourPionCzyzCurrent::beta_f0_, -1, 1.0, 0, 0,
     false, false, Interface::nolimits);
  
  static ParVector<FourPionCzyzCurrent,double> interfacebeta_omega
    ("beta_omega",
     "The coefficients for the sum over rho resonances in the omega term",
     &FourPionCzyzCurrent::beta_omega_, -1, 1.0, 0, 0,
     false, false, Interface::nolimits);
  
  static ParVector<FourPionCzyzCurrent,double> interfacebeta_B
    ("beta_B",
     "The coefficients for the sum over rho resonances in the B_rho term",
     &FourPionCzyzCurrent::beta_B_, -1, 1.0, 0, 0,
     false, false, Interface::nolimits);
  
  static ParVector<FourPionCzyzCurrent,double> interfacebeta_bar
    ("beta_bar",
     "The coefficients for the sum over rho resonances in the T_rho term",
     &FourPionCzyzCurrent::beta_bar_, -1, 1.0, 0, 0,
     false, false, Interface::nolimits);

  static Parameter<FourPionCzyzCurrent,InvEnergy2> interfacec_a1
    ("c_a1",
     "The coupling for the a_1 channel",
     &FourPionCzyzCurrent::c_a1_, 1./GeV2, -255./GeV2, -1e5/GeV2, 1e5/GeV2,
     false, false, Interface::limited);

  static Parameter<FourPionCzyzCurrent,InvEnergy2> interfacec_f0
    ("c_f0",
     "The coupling for the f_0 channel",
     &FourPionCzyzCurrent::c_f0_, 1./GeV2, 64./GeV2, -1e5/GeV2, 1e5/GeV2,
     false, false, Interface::limited);
  
  static Parameter<FourPionCzyzCurrent,InvEnergy> interfacec_omega
    ("c_omega",
     "The coupling for the omega channel",
     &FourPionCzyzCurrent::c_omega_, 1./GeV, -1.47/GeV, -1e5/GeV, 1e5/GeV,
     false, false, Interface::limited);
  
  static Parameter<FourPionCzyzCurrent,double> interfacec_rho
    ("c_rho",
     "The coupling for the rho channel",
     &FourPionCzyzCurrent::c_rho_, -2.46, 0, 0,
     false, false, Interface::nolimits);

  static Parameter<FourPionCzyzCurrent,double> interfaceg_rho_pi_pi
    ("g_rho_pi_pi",
     "The coupling of rho to two pions",
     &FourPionCzyzCurrent::g_rho_pi_pi_, 5.997, 0, 0,
     false, false, Interface::nolimits);

  static Parameter<FourPionCzyzCurrent,ThePEG::Qty<std::ratio<0,1>, std::ratio<-5,1>, std::ratio<0,1> >> interfaceg_omega_pi_rho
    ("g_omega_pi_rho",
     "The coupling of omega to rho and pi",
     &FourPionCzyzCurrent::g_omega_pi_rho_, 1./GeV/GeV2/GeV2, 42.3/GeV/GeV2/GeV2, 0.0/GeV/GeV2/GeV2, 1e5/GeV/GeV2/GeV2,
     false, false, Interface::limited);
  
  static Parameter<FourPionCzyzCurrent,Energy2> interfaceg_rho_gamma
    ("g_rho_gamma",
     "The coupling of the rho to the photon",
     &FourPionCzyzCurrent::g_rho_gamma_, GeV2, 0.1212*GeV2, 0.0*GeV2, 10.0*GeV2,
     false, false, Interface::limited);
}

void FourPionCzyzCurrent::doinit() {
  WeakDecayCurrent::doinit();
  mpip_ = getParticleData(211)->mass();
  mpi0_ = getParticleData(111)->mass();
}

// complete the construction of the decay mode for integration
bool FourPionCzyzCurrent::createMode(int icharge, unsigned int imode,
				      DecayPhaseSpaceModePtr mode,
				      unsigned int iloc,int ires,
				      DecayPhaseSpaceChannelPtr phase,
				      Energy upp) {
  // check the charge
  if((abs(icharge)!=3&&imode<2) ||
     (imode>=2&&icharge!=0)) return false;
  // check that the modes are kinematical allowed
  Energy min(ZERO);
  if(imode==0)      min =    mpip_+3.*mpi0_;
  else if(imode==1) min = 3.*mpip_+   mpi0_;
  else              min = 2.*mpip_+2.*mpi0_;
  if(min>upp) return false;
  // resonaces we need
  tPDPtr rho0[3]  = {getParticleData( 113),getParticleData( 100113),getParticleData( 30113)};
  tPDPtr rhop[3]  = {getParticleData( 213),getParticleData( 100213),getParticleData( 30213)};
  tPDPtr rhom[3]  = {getParticleData(-213),getParticleData(-100213),getParticleData(-30213)};
  tPDPtr omega(getParticleData(ParticleID::omega));
  tPDPtr f0(getParticleData(10221));
  tPDPtr a1p  = getParticleData(ParticleID::a_1plus);
  tPDPtr a1m  = getParticleData(ParticleID::a_1minus);
  tPDPtr a10  = getParticleData(ParticleID::a_10);
  DecayPhaseSpaceChannelPtr newchannel;
  if(icharge==3) {
    swap(rhop,rhom);
    swap(a1p,a1m);
  }
  if(imode==0) {
  }
  else if(imode==1) {
  }
  // pi0 pi0 pi+ pi-
  else if(imode==2||imode==3) {
    // rho a_1 channels
    for(unsigned int irho1=0;irho1<3;++irho1) {
      for(unsigned int irho2=0;irho2<3;++irho2) {
	newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
	newchannel->addIntermediate(rho0[irho1],0,0.0,-ires-1,iloc+3);
	newchannel->addIntermediate(a1p        ,0,0.0,-ires-2,iloc+1);
	newchannel->addIntermediate(rhop[irho2],0,0.0, iloc  ,iloc+2);
	mode->addChannel(newchannel);
	newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
	newchannel->addIntermediate(rho0[irho1],0,0.0,-ires-1,iloc+3);
	newchannel->addIntermediate(a1p        ,0,0.0,-ires-2,iloc  );
	newchannel->addIntermediate(rhop[irho2],0,0.0, iloc+1,iloc+2);
	mode->addChannel(newchannel);
	newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
	newchannel->addIntermediate(rho0[irho1],0,0.0,-ires-1,iloc+2);
	newchannel->addIntermediate(a1m        ,0,0.0,-ires-2,iloc+1);
	newchannel->addIntermediate(rhom[irho2],0,0.0, iloc  ,iloc+3);
	mode->addChannel(newchannel);
	newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
	newchannel->addIntermediate(rho0[irho1],0,0.0,-ires-1,iloc+2);
	newchannel->addIntermediate(a1m        ,0,0.0,-ires-2,iloc  );
	newchannel->addIntermediate(rhom[irho2],0,0.0, iloc+1,iloc+3);
	mode->addChannel(newchannel);
      }
    }
    // rho omega channels
    for(unsigned int irho=0;irho<3;++irho) {
      newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(rho0[irho],0,0.0,-ires-1,iloc );
      newchannel->addIntermediate(omega     ,0,0.0,-ires-2,iloc+1);
      newchannel->addIntermediate(rho0[0]   ,0,0.0, iloc+2,iloc+3);
      mode->addChannel(newchannel);
      newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(rho0[irho],0,0.0,-ires-1,iloc );
      newchannel->addIntermediate(omega     ,0,0.0,-ires-2,iloc+2);
      newchannel->addIntermediate(rho0[0]   ,0,0.0, iloc+1,iloc+3);
      mode->addChannel(newchannel);
      newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(rho0[irho],0,0.0,-ires-1,iloc );
      newchannel->addIntermediate(omega     ,0,0.0,-ires-2,iloc+3);
      newchannel->addIntermediate(rho0[0]   ,0,0.0, iloc+1,iloc+2);
      mode->addChannel(newchannel);
      newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(rho0[irho],0,0.0,-ires-1,iloc+1);
      newchannel->addIntermediate(omega     ,0,0.0,-ires-2,iloc  );
      newchannel->addIntermediate(rho0[0]   ,0,0.0, iloc+2,iloc+3);
      mode->addChannel(newchannel);
      newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(rho0[irho],0,0.0,-ires-1,iloc+1);
      newchannel->addIntermediate(omega     ,0,0.0,-ires-2,iloc+2);
      newchannel->addIntermediate(rho0[0]   ,0,0.0, iloc  ,iloc+3);
      mode->addChannel(newchannel);
      newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(rho0[irho],0,0.0,-ires-1,iloc+1);
      newchannel->addIntermediate(omega     ,0,0.0,-ires-2,iloc+3);
      newchannel->addIntermediate(rho0[0]   ,0,0.0, iloc  ,iloc+2);
      mode->addChannel(newchannel);
    }
    // rho f_0 rho channels
    for(unsigned int irho1=0;irho1<3;++irho1) {
      for(unsigned int irho2=0;irho2<3;++irho2) {
	newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
	newchannel->addIntermediate(rho0[irho1],0,0.0,-ires-1,-ires-2);
	newchannel->addIntermediate(f0         ,0,0.0,iloc   , iloc+1);
	newchannel->addIntermediate(rho0[irho2],0,0.0,iloc+2 , iloc+3);
	mode->addChannel(newchannel);
      }
    }
  }
  else {
  }
  return true;
}

// the particles produced by the current
tPDVector FourPionCzyzCurrent::particles(int icharge, unsigned int imode,
					 int,int) {
  tPDVector output(4);
  tPDPtr pi0=getParticleData(ParticleID::pi0);
  tPDPtr pip=getParticleData(ParticleID::pi0);
  tPDPtr pim=pip->CC();
  if(imode==0) {
    output[0]=pim;
    output[1]=pi0;
    output[2]=pi0;
    output[3]=pi0;
  }
  else if(imode==1) {
    output[0]=pim;
    output[1]=pim;
    output[2]=pip;
    output[3]=pi0;
  }
  else if(imode==2||imode==3) {
    output[0]=pip;
    output[1]=pim;
    output[2]=pi0;
    output[3]=pi0;
  }
  else {
    output[0]=pip;
    output[1]=pip;
    output[2]=pim;
    output[3]=pim;
  }
  if(icharge==3) {
    for(unsigned int ix=0;ix<output.size();++ix) {
      if(output[ix]->CC()) output[ix]=output[ix]->CC();
    }
  }
  // return the answer
  return output;
}


// hadronic current   
vector<LorentzPolarizationVectorE> 
FourPionCzyzCurrent::current(const int imode, const int ichan,
			      Energy & scale,const ParticleVector & decay,
			      DecayIntegrator::MEOption meopt) const {
  useMe();
  if(meopt==DecayIntegrator::Terminate) {
    for(unsigned int ix=0;ix<3;++ix)
      ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
    return vector<LorentzPolarizationVectorE>(1,LorentzPolarizationVectorE());
  }
  // the momenta of the particles
  Lorentz5Momentum q1(decay[0]->momentum()),q2(decay[2]->momentum()),
    q3(decay[1]->momentum()),q4(decay[3]->momentum());
  Lorentz5Momentum Q(q1+q2+q3+q4);
  Q.rescaleMass();
  LorentzVector<complex<InvEnergy> > output;
  if(imode==0) {
    output = baseCurrent(Q.mass2(),Q,q3,q4,q1,q2)+
      baseCurrent(Q.mass2(),Q,q2,q4,q1,q3)+
      baseCurrent(Q.mass2(),Q,q2,q3,q1,q4);
  }
  else if(imode==1) {
    output = baseCurrent(Q.mass2(),Q,q3,q2,q1,q4)+
      baseCurrent(Q.mass2(),Q,q3,q1,q2,q4);
  }
  else if(imode==2||imode==3) {
    output = sqrt(0.5)*baseCurrent(Q.mass2(),Q,q3,q4,q1,q2);
  }
  else if(imode==4||imode==5) {
    output = sqrt(0.5)*
      (baseCurrent(Q.mass2(),Q,q2,q4,q1,q3)+
       baseCurrent(Q.mass2(),Q,q1,q4,q2,q3)+
       baseCurrent(Q.mass2(),Q,q2,q3,q1,q4)+
       baseCurrent(Q.mass2(),Q,q1,q3,q2,q4));
  }
  return  vector<LorentzPolarizationVectorE>(1,output*Q.mass2());
}
   
bool FourPionCzyzCurrent::accept(vector<int> id) {
//   // check there are only two particles
//   if(id.size()!=2) return false;
//   // pion modes
//   if((abs(id[0])==ParticleID::piplus  &&     id[1] ==ParticleID::pi0   ) ||
//      (    id[0] ==ParticleID::pi0     && abs(id[1])==ParticleID::piplus))
//     return true;
//   else if((id[0]==ParticleID::piminus && id[1]==ParticleID::piplus) ||
// 	  (id[0]==ParticleID::piplus  && id[1]==ParticleID::piminus))
//     return true;
//   else
//     return false;
}

// the decay mode
unsigned int FourPionCzyzCurrent::decayMode(vector<int> idout) {
  return 0;
}

// output the information for the database
void FourPionCzyzCurrent::dataBaseOutput(ofstream & output,bool header,
					bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::FourPionCzyzCurrent " 
  		    << name() << " HwWeakCurrents.so\n";
  for(unsigned int ix=0;ix<rhoMasses_.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":RhoMasses " << ix << " " << rhoMasses_[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<rhoWidths_.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":RhoWidths " << ix << " " << rhoWidths_[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<rhoMasses_.size();++ix) {
    if(ix<4) output << "newdef ";
    else     output << "insert ";
    output << name() << ":RhoMassesFrho " << ix << " " << rhoMasses_Frho_[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<rhoWidths_.size();++ix) {
    if(ix<4) output << "newdef ";
    else     output << "insert ";
    output << name() << ":RhoWidthsFrho " << ix << " " << rhoWidths_Frho_[ix]/GeV << "\n";
  }
  output << "newdef " << name() << ":omegaMass "  << omegaMass_/GeV  << "\n";
  output << "newdef " << name() << ":omegaWidth " << omegaWidth_/GeV << "\n";
  output << "newdef " << name() << ":f0Mass "  << f0Mass_/GeV  << "\n";
  output << "newdef " << name() << ":f0Width " << f0Width_/GeV << "\n";
  output << "newdef " << name() << ":a1Mass "  << a1Mass_/GeV  << "\n";
  output << "newdef " << name() << ":a1Width " << a1Width_/GeV << "\n";
  for(unsigned int ix=0;ix<beta_a1_.size();++ix) {
    if(ix<4) output << "newdef ";
    else     output << "insert ";
    output << name() << ":beta_a1 " << ix << " " << beta_a1_[ix] << "\n";
  }
  for(unsigned int ix=0;ix<beta_f0_.size();++ix) {
    if(ix<4) output << "newdef ";
    else     output << "insert ";
    output << name() << ":beta_f0 " << ix << " " << beta_f0_[ix] << "\n";
  }
  for(unsigned int ix=0;ix<beta_omega_.size();++ix) {
    if(ix<4) output << "newdef ";
    else     output << "insert ";
    output << name() << ":beta_omega " << ix << " " << beta_omega_[ix] << "\n";
  }
  for(unsigned int ix=0;ix<beta_B_.size();++ix) {
    if(ix<2) output << "newdef ";
    else     output << "insert ";
    output << name() << ":beta_B " << ix << " " << beta_B_[ix] << "\n";
  }
  for(unsigned int ix=0;ix<beta_bar_.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":beta_bar " << ix << " " << beta_bar_[ix] << "\n";
  }
  output << "newdef " << name() << ":c_a1    " << c_a1_*GeV2   << "\n";  
  output << "newdef " << name() << ":c_f0    " << c_f0_*GeV2   << "\n";  
  output << "newdef " << name() << ":c_omega " << c_omega_*GeV << "\n";  
  output << "newdef " << name() << ":c_rho   " << c_rho_   << "\n";  
  output << "newdef " << name() << ":g_rho_pi_pi   " <<  g_rho_pi_pi_  << "\n";  
  output << "newdef " << name() << ":g_omega_pi_rho   "
	 << g_omega_pi_rho_*GeV2*GeV2*GeV   << "\n";  
  output << "newdef " << name() << ":g_rho_gamma   "
	 <<g_rho_gamma_/GeV2    << "\n";
  WeakDecayCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
  		    << fullName() << "\";" << endl;
}

LorentzVector<complex<InvEnergy> > FourPionCzyzCurrent::
baseCurrent(Energy2 Q2,
	    const Lorentz5Momentum & Q,
	    const Lorentz5Momentum & q1,
	    const Lorentz5Momentum & q2,
	    const Lorentz5Momentum & q3,
	    const Lorentz5Momentum & q4) const {
  // various dot products we'll need
  Energy2 m12 = sqr(q1.mass()), m22 = sqr(q2.mass());
  Energy2 m32 = sqr(q3.mass()), m42 = sqr(q4.mass());
  Energy2 Qq1 = Q*q1, Qq2 = Q*q2, Qq3 = Q*q3, Qq4 = Q*q4;
  Energy2 Qm32= Q2-2.*Qq3+m32;
  Energy2 Qm42= Q2-2.*Qq4+m42;
  Energy2 q1q2 = q1*q2, q1q3 = q1*q3, q1q4 = q1*q4;
  Energy2 q2q3 = q2*q3, q2q4 = q2*q4, q3q4 = q3*q4;
  // first the a_1 terms from A.3 0804.0359 (N.B. sign due defns)
  // common coefficent
  complex<InvEnergy2> a1_fact = -c_a1_*
    Resonance::F_rho(Q2,beta_a1_,rhoMasses_Frho_,
		     rhoWidths_Frho_,mpip_,mpip_);
  // a_1 BWs
  Complex bw_a1_Qm3 = Resonance::BreitWignera1(Qm32,a1Mass_,a1Width_);
  Complex bw_a1_Qm4 = Resonance::BreitWignera1(Qm42,a1Mass_,a1Width_);
  // first term
  Complex Brhoq1q4 = bw_a1_Qm4*
    Resonance::F_rho(m12+m42+2.*q1q4,beta_B_,rhoMasses_,rhoWidths_,mpip_,mpip_);
  double dot1  = (q1q2-q2q4)/Qm32;
  double dot1B = (Qq1-Qq4)/Q2;
  double dot1C = Qq3/Q2*dot1;
  // second term
  Complex Brhoq2q4 = bw_a1_Qm4*
    Resonance::F_rho(m12+m42+2.*q2q4,beta_B_,rhoMasses_,rhoWidths_,mpip_,mpip_);
  double dot2  = (q1q2-q1q4)/Qm32;
  double dot2B = (Qq2-Qq4)/Q2;
  double dot2C = Qq3/Q2*dot2;
  // third term
  Complex Brhoq1q3 = bw_a1_Qm3*
    Resonance::F_rho(m12+m32+2.*q1q3,beta_B_,rhoMasses_,rhoWidths_,mpip_,mpip_);
  double dot3  = (q1q2-q2q3)/Qm42;
  double dot3B = (Qq1-Qq3)/Q2;
  double dot3C = Qq4/Q2*dot1;
  // fourth term
  Complex Brhoq2q3 = bw_a1_Qm3*
    Resonance::F_rho(m12+m32+2.*q2q3,beta_B_,rhoMasses_,rhoWidths_,mpip_,mpip_);
  double dot4  = (q1q2-q1q3)/Qm42;
  double dot4B = (Qq2-Qq3)/Q2;
  double dot4C = Qq4/Q2*dot2;
  // coefficients of the momenta to construct the current
  // a_1 terms
  complex<InvEnergy2> c1 =-0.5*a1_fact*
    (Brhoq1q4*( 3.-dot1)+Brhoq2q4*( 1.-dot2)+
     Brhoq1q3*(-3.+dot3)+Brhoq2q3*(-1.+dot4));
  complex<InvEnergy2> c2 =-0.5*a1_fact*
    (Brhoq1q4*( 1.-dot1)+Brhoq2q4*( 3.-dot2)+
     Brhoq1q3*(-1.+dot3)+Brhoq2q3*(-3.+dot4));
  complex<InvEnergy2> c34=-0.5*a1_fact*
    (Brhoq1q4*( 1.+dot1)+Brhoq2q4*( 1.+dot2)+
     Brhoq1q3*( 1.+dot3)+Brhoq2q3*( 1.+dot4));
  complex<InvEnergy2> cq =-0.5*a1_fact*
    (Brhoq1q4*(-1.+dot1-2.*dot1B-2.*dot1C)+
     Brhoq2q4*(-1.+dot2-2.*dot2B-2.*dot2C)+
     Brhoq1q3*( 1.-dot3+2.*dot3B+2.*dot3C)+
     Brhoq2q3*( 1.-dot4+2.*dot4B+2.*dot4C));
  // f_0
  complex<InvEnergy2> f0fact = c_f0_*
    Resonance::F_rho(Q2,beta_f0_,rhoMasses_Frho_,
		     rhoWidths_Frho_,mpip_,mpip_)*
    Resonance::F_rho(m32+m42+2.*q3q4,beta_bar_,rhoMasses_,rhoWidths_,mpip_,mpip_)*
    Resonance::BreitWignerSWave(m12+m22+2.*q1q2,f0Mass_,f0Width_,mpip_,mpip_);
  // add contribution to the coefficients
  c34 -= f0fact;
  cq += f0fact*(Qq3-Qq4)/Q2;
  complex<ThePEG::Qty<std::ratio<0,1>, std::ratio<-6,1>, std::ratio<0,1> > >
    wfact = 2.*c_omega_*g_omega_pi_rho_*g_rho_pi_pi_*
    Resonance::F_rho(Q2,beta_omega_,rhoMasses_Frho_,
		     rhoWidths_Frho_,mpip_,mpip_);
  complex<ThePEG::Qty<std::ratio<0,1>, std::ratio<-6,1>, std::ratio<0,1> > >
    bw_omega_1 = wfact
    *Resonance::BreitWignerFW(m12-2.*Qq1+Q2,omegaMass_,omegaWidth_)
    *Resonance::H(rhoMasses_[0],rhoWidths_[1],m22+2.*q2q3+m32,m22+2.*q2q4+m42,
		  m32+2.*q3q4+m42,mpip_,mpip_);
  complex<ThePEG::Qty<std::ratio<0,1>, std::ratio<-6,1>, std::ratio<0,1> > >
    bw_omega_2 = wfact
    *Resonance::BreitWignerFW(m22-2.*Qq2+Q2,omegaMass_,omegaWidth_)
    *Resonance::H(rhoMasses_[0],rhoWidths_[1],m12+2.*q1q3+m32,m12+2.*q1q4+m42,
		  m32+2.*q3q4+m42,mpip_,mpip_);
  c1 -= bw_omega_2*(q2q4*Qq3-q2q3*Qq4);
  c2 -= bw_omega_1*(q1q4*Qq3-q1q3*Qq4);
  complex<InvEnergy2> c3 =
    -bw_omega_1*(q1q2*Qq4-q1q4*Qq2)
    -bw_omega_2*(q1q2*Qq4-q2q4*Qq1);
  complex<InvEnergy2> c4 =
    -bw_omega_1*(q1q3*Qq2-q1q2*Qq3)
    -bw_omega_2*(q2q3*Qq1-q1q2*Qq3);
  // the rho term
  complex<InvEnergy2> BW12_q1q3 =
    Resonance::BreitWignerDiff(m12+m32+2.*q1q3,rhoMasses_[0],rhoWidths_[0],
			       rhoMasses_[1],rhoWidths_[1],mpip_,mpip_);
  complex<InvEnergy2> BW12_q1q4 =
    Resonance::BreitWignerDiff(m12+m42+2.*q1q4,rhoMasses_[0],rhoWidths_[0],
			       rhoMasses_[1],rhoWidths_[1],mpip_,mpip_);
  complex<InvEnergy2> BW12_q2q3 =
    Resonance::BreitWignerDiff(m22+m32+2.*q2q3,rhoMasses_[0],rhoWidths_[0],
			       rhoMasses_[1],rhoWidths_[1],mpip_,mpip_);
  complex<InvEnergy2> BW12_q2q4 =
    Resonance::BreitWignerDiff(m22+m42+2.*q2q4,rhoMasses_[0],rhoWidths_[0],
			       rhoMasses_[1],rhoWidths_[1],mpip_,mpip_);
  Energy2 d1 = Qq2 - Qq4 + 2.*q2q3 - 2.*q3q4;
  Energy2 d2 = Qq2 - Qq3 + 2.*q2q4 - 2.*q3q4;
  Energy2 d3 = Qq1 - Qq4 + 2.*q1q3 - 2.*q3q4;
  Energy2 d4 = Qq1 - Qq3 + 2.*q1q4 - 2.*q3q4;
  Energy2 d5 = Qq1 - Qq3 + 2.*q1q2 - 2.*q2q3;
  Energy2 d6 = Qq2 - Qq3 + 2.*q1q2 - 2.*q1q3;
  Energy2 d7 = Qq1 - Qq4 + 2.*q1q2 - 2.*q2q4;
  Energy2 d8 = Qq2 - Qq4 + 2.*q1q2 - 2.*q1q4;  
  Complex pre_rho = c_rho_*pow(g_rho_pi_pi_,3)*g_rho_gamma_*
    Resonance::BreitWignerDiff(Q2,rhoMasses_[0],rhoWidths_[0],
			       rhoMasses_[1],rhoWidths_[1],mpip_,mpip_);
  LorentzVector<complex<InvEnergy> > v_rho =
    pre_rho*(BW12_q1q3*(BW12_q2q4*d1 + 2.) - BW12_q1q4*(BW12_q2q3*d2 + 2.))*q1 + 
    pre_rho*(BW12_q2q3*(BW12_q1q4*d3 + 2.) - BW12_q2q4*(BW12_q1q3*d4 + 2.))*q2 -
    pre_rho*(BW12_q2q3*(BW12_q1q4*d7 + 2.) + BW12_q1q3*(BW12_q2q4*d8 + 2.))*q3 +
    pre_rho*(BW12_q2q4*(BW12_q1q3*d5 + 2.) + BW12_q1q4*(BW12_q2q3*d6 + 2.))*q4;
  complex<InvEnergy2> vdot = (Q*v_rho)/Q2;  
  v_rho  = -v_rho + vdot*Q;
  // put everything together
  return c1*q1+c2*q2+(c3+c34)*q3+(c4-c34)*q4+cq*Q + v_rho; 
}
