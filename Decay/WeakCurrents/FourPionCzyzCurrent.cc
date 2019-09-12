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
  : mpip_(140*MeV), mpi0_(140*MeV), channelMap_(6,vector<int>()) {
  // Masses and widths of the particles
  // rho (PDG for most of current)
  rhoMasses_ = {0.7755*GeV,1.459*GeV,1.72*GeV};
  rhoWidths_ = {0.1494*GeV,0.4  *GeV,0.25*GeV};
  // fitted for F_rho
  rhoMasses_Frho_ = {0.7755*GeV,1.437*GeV             ,1.738*GeV             ,2.12*GeV};
  rhoWidths_Frho_ = {0.1494*GeV,0.6784258847826656*GeV,0.8049153117715863*GeV,0.20924569673342333*GeV};
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
  beta_a1_    ={1.,-0.051864694215520334, -0.0416090742847935,-0.0018940213981381284 };
  beta_f0_    ={1.,    73673.55747104406, -26116.259838644528,     332.84652898870786};
  beta_omega_ ={1., -0.3668543468911129,0.03624673983924276 , -0.0047186283455050064 };
  beta_B_     ={1.,-0.145};
  beta_bar_   ={1.,0.08,-0.0075};  
  // couplings for the various terms
  c_a1_    = -201.7929015686851/GeV2;
  c_f0_    =  124.09577796839454/GeV2;
  c_omega_ = -1.5792052826500913/GeV;
  c_rho_   = -2.308949096570198;
  // meson meson meson couplings
  g_rho_pi_pi_    = 5.997;
  g_omega_pi_rho_ = 42.3/GeV/GeV2/GeV2;
  g_rho_gamma_    = 0.1212*GeV2;
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
     << ounit(mpip_,GeV) << ounit(mpi0_,GeV) << channelMap_;
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
     >> iunit(mpip_,GeV) >> iunit(mpi0_,GeV) >> channelMap_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<FourPionCzyzCurrent,WeakCurrent>
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
  WeakCurrent::doinit();
  mpip_ = getParticleData(211)->mass();
  mpi0_ = getParticleData(111)->mass();
  // test of the current for a fixed momentum configuration
  // Lorentz5Momentum q1(0.13061870567796208*GeV,-0.21736300316234394*GeV,
  // 		      0.51725254282500699*GeV,0.59167288008090657*GeV);
  // Lorentz5Momentum q2(-1.1388573867484255 *GeV,      0.37727761929037396 *GeV,      0.31336706796993302 *GeV,    1.2472979400305677*GeV);
  // Lorentz5Momentum q3(0.11806773412672231 *GeV,      0.17873024832600765   *GeV,    0.10345508827447017   *GeV,   0.27580297667647757*GeV );
  // Lorentz5Momentum q4 (7.7487017488620830E-002*GeV, 0.16118198754624435    *GeV,    6.5813182962706746E-002*GeV, 0.23620982448991124*GeV );
  // q1.rescaleMass();
  // q2.rescaleMass();
  // q2.rescaleMass();
  // q3.rescaleMass();
  // cerr << q1/GeV << " " << q1.mass()/GeV << "\n";
  // cerr << q2/GeV << " " << q2.mass()/GeV << "\n";
  // cerr << q3/GeV << " " << q3.mass()/GeV << "\n";
  // cerr << q4/GeV << " " << q4.mass()/GeV << "\n";
  // Lorentz5Momentum Q(q1+q2+q3+q4);
  // Q.rescaleMass();
 
  // LorentzVector<complex<InvEnergy> > base = baseCurrent(Q.mass2(),tcPDPtr(),-1,Q,q1,q2,q3,q4);
  // LorentzVector<complex<InvEnergy> > test( Complex(  376.35697290395467     ,  66.446392015809550     )/GeV,
  // 					   Complex( -135.73591401998152     ,  112.36912660073307     )/GeV,
  // 					   Complex(  83.215302375273723     , -54.986430577097920     )/GeV,
  // 					   Complex( -123.56434266557559     , -22.465096431505703     )/GeV);
  // cerr << "current test X " << (base.x()-test.x())/(base.x()+test.x()) << "\n";
  // cerr << "current test Y " << (base.y()-test.y())/(base.x()+test.y()) << "\n";
  // cerr << "current test Z " << (base.z()-test.z())/(base.x()+test.z()) << "\n";
  // cerr << "current test T " << (base.t()-test.t())/(base.x()+test.t()) << "\n";
}

void FourPionCzyzCurrent::createChannels(unsigned int imode,
					 int icharge,  tcPDPtr resonance,
					 unsigned int iloc,int ires,
					 tPDVector outgoing, PhaseSpaceModePtr mode,
					 PhaseSpaceChannel phase,
					 unsigned int j1, unsigned int j2,
					 unsigned int j3, unsigned int j4,
					 int & nchan) {
  tPDPtr rho0[3]  = {getParticleData( 113),getParticleData( 100113),getParticleData( 30113)};
  tPDPtr rhop[3]  = {getParticleData( 213),getParticleData( 100213),getParticleData( 30213)};
  tPDPtr rhom[3]  = {getParticleData(-213),getParticleData(-100213),getParticleData(-30213)};
  tPDPtr omega(getParticleData(ParticleID::omega));
  tPDPtr f0(getParticleData(10221));
  tPDPtr a1p  = getParticleData(ParticleID::a_1plus);
  tPDPtr a1m  = getParticleData(ParticleID::a_1minus);
  tPDPtr a10  = getParticleData(ParticleID::a_10);
  if(icharge==3) {
    swap(rhop,rhom);
    swap(a1p,a1m);
  }
  int rhoCharge;
  tPDPtr rho;
  tPDPtr rhoin;
  // first the a1 channels
  for(unsigned int irho=0;irho<3;++irho) {
    tPDPtr rhoin = icharge==0 ? rho0[irho] : rhom[irho];
    if(resonance && rhoin!=resonance) {
      nchan+=8;
      continue;
    }
    int a1Charge;
    tPDPtr a1;
    for(unsigned int irho2=0;irho2<2;++irho2) {
      // // find the a1
      a1Charge = outgoing[j1-1]->iCharge()+outgoing[j2-1]->iCharge()+outgoing[j4-1]->iCharge();
      a1 = a1Charge==0 ? a10 : a1p;
      if(a1->iCharge()!=a1Charge) a1=a1m;
      assert(abs(a1Charge)<=3);
      // find the rho
      rhoCharge = outgoing[j1-1]->iCharge()+outgoing[j4-1]->iCharge();
      rho = rhoCharge==0 ? rho0[irho2] : rhop[irho2];
      if(rho->iCharge()!=rhoCharge) rho = rhom[irho2];
      assert(abs(rhoCharge)<=3);
      mode->addChannel((PhaseSpaceChannel(phase),ires,rhoin,ires+1,a1,ires+1,iloc+j3,
      			ires+2,rho,ires+2,iloc+j2,ires+3,iloc+j1,ires+3,iloc+j4));
      ++nchan; channelMap_[imode].push_back(nchan);
      // find the rho
      if(imode!=4) {
	rhoCharge = outgoing[j2-1]->iCharge()+outgoing[j4-1]->iCharge();
	rho = rhoCharge==0 ? rho0[irho2] : rhop[irho2];
	if(rho->iCharge()!=rhoCharge) rho = rhom[irho2];
	assert(abs(rhoCharge)<=3);
	mode->addChannel((PhaseSpaceChannel(phase),ires,rhoin,ires+1,a1,ires+1,iloc+j3,
			  ires+2,rho,ires+2,iloc+j1,ires+3,iloc+j2,ires+3,iloc+j4));
	++nchan; channelMap_[imode].push_back(nchan);
      }
      else ++nchan;
      // find the second a1
      a1Charge = outgoing[j1-1]->iCharge()+outgoing[j2-1]->iCharge()+outgoing[j3-1]->iCharge();
      a1 = a1Charge==0 ? a10 : a1p;
      if(a1->iCharge()!=a1Charge) a1=a1m;
      assert(abs(a1Charge)<=3);
      // find the rho
      if(imode!=4) {
	rhoCharge = outgoing[j1-1]->iCharge()+outgoing[j3-1]->iCharge();
	rho = rhoCharge==0 ? rho0[irho2] : rhop[irho2];
	if(rho->iCharge()!=rhoCharge) rho = rhom[irho2];
	assert(abs(rhoCharge)<=3);
	mode->addChannel((PhaseSpaceChannel(phase),ires,rhoin,ires+1,a1,ires+1,iloc+j4,
			  ires+2,rho,ires+2,iloc+j2,ires+3,iloc+j1,ires+3,iloc+j3));
	++nchan; channelMap_[imode].push_back(nchan);
      }
      else
	++nchan;
      // find the rho
      if(imode!=1) {
	rhoCharge = outgoing[j2-1]->iCharge()+outgoing[j3-1]->iCharge();
	rho = rhoCharge==0 ? rho0[irho2] : rhop[irho2];
	if(rho->iCharge()!=rhoCharge) rho = rhom[irho2];
	assert(abs(rhoCharge)<=3);
	mode->addChannel((PhaseSpaceChannel(phase),ires,rhoin,ires+1,a1,ires+1,iloc+j4,
			  ires+2,rho,ires+2,iloc+j1,ires+3,iloc+j2,ires+3,iloc+j3));
	++nchan; channelMap_[imode].push_back(nchan);
      }
      else
	++nchan;
    }
  }
  // now the f_0 channel
  for(unsigned int irho=0;irho<3;++irho) {
    tPDPtr rhoin = icharge==0 ? rho0[irho] : rhom[irho];
    if(resonance && rhoin!=resonance) continue;
    rhoCharge = outgoing[j1-1]->iCharge()+outgoing[j2-1]->iCharge();
    for(unsigned int irho2=0;irho2<3;++irho2) {
      rho= rhoCharge==0 ? rho0[irho2] : rhop[irho2];
      if(rho->iCharge()!=rhoCharge) rho = rhom[irho2];
      assert(abs(rhoCharge)<=3);
      mode->addChannel((PhaseSpaceChannel(phase),ires,rhoin,ires+1,rho,ires+1,f0,
  			ires+2,iloc+j1,ires+2,iloc+j2,ires+3,iloc+j3,ires+3,iloc+j4));
      ++nchan; channelMap_[imode].push_back(nchan);
    }
  }
  // now the omega channels
  if(imode>=1&&imode<=3) {
    for(unsigned int irho=0;irho<3;++irho) {
      tPDPtr rhoin = icharge==0 ? rho0[irho] : rhom[irho];
      if(resonance && rhoin!=resonance) continue;
      if(imode!=1) {
	rhoCharge = outgoing[j2-1]->iCharge()+outgoing[j3-1]->iCharge();
	rho= rhoCharge==0 ? rho0[0] : rhop[0];
	if(rho->iCharge()!=rhoCharge) rho = rhom[0];
	assert(abs(rhoCharge)<=3);
	mode->addChannel((PhaseSpaceChannel(phase),ires,rhoin,ires+1,omega,ires+1,iloc+j1,
			  ires+2,rho,ires+2,iloc+j4,ires+3,iloc+j2,ires+3,iloc+j3));
	++nchan; channelMap_[imode].push_back(nchan);
	rhoCharge = outgoing[j2-1]->iCharge()+outgoing[j4-1]->iCharge();
	rho= rhoCharge==0 ? rho0[0] : rhop[0];
	if(rho->iCharge()!=rhoCharge) rho = rhom[0];
	assert(abs(rhoCharge)<=3);
	mode->addChannel((PhaseSpaceChannel(phase),ires,rhoin,ires+1,omega,ires+1,iloc+j1,
			  ires+2,rho,ires+2,iloc+j3,ires+3,iloc+j2,ires+3,iloc+j4));
	++nchan; channelMap_[imode].push_back(nchan);
	rhoCharge = outgoing[j3-1]->iCharge()+outgoing[j4-1]->iCharge();
	rho= rhoCharge==0 ? rho0[0] : rhop[0];
	if(rho->iCharge()!=rhoCharge) rho = rhom[0];
	assert(abs(rhoCharge)<=3);
	mode->addChannel((PhaseSpaceChannel(phase),ires,rhoin,ires+1,omega,ires+1,iloc+j1,
			  ires+2,rho,ires+2,iloc+j2,ires+3,iloc+j3,ires+3,iloc+j4));
	++nchan; channelMap_[imode].push_back(nchan);
      }
      else nchan+=3;
      rhoCharge = outgoing[j1-1]->iCharge()+outgoing[j3-1]->iCharge();
      rho= rhoCharge==0 ? rho0[0] : rhop[0];
      if(rho->iCharge()!=rhoCharge) rho = rhom[0];
      assert(abs(rhoCharge)<=3);
      mode->addChannel((PhaseSpaceChannel(phase),ires,rhoin,ires+1,omega,ires+1,iloc+j2,
			ires+2,rho,ires+2,iloc+j4,ires+3,iloc+j1,ires+3,iloc+j3));
      ++nchan; channelMap_[imode].push_back(nchan);
      rhoCharge = outgoing[j1-1]->iCharge()+outgoing[j4-1]->iCharge();
      rho= rhoCharge==0 ? rho0[0] : rhop[0];
      if(rho->iCharge()!=rhoCharge) rho = rhom[0];
      assert(abs(rhoCharge)<=3);
      mode->addChannel((PhaseSpaceChannel(phase),ires,rhoin,ires+1,omega,ires+1,iloc+j2,
			ires+2,rho,ires+2,iloc+j3,ires+3,iloc+j1,ires+3,iloc+j4));
      ++nchan; channelMap_[imode].push_back(nchan);
      rhoCharge = outgoing[j3-1]->iCharge()+outgoing[j4-1]->iCharge();
      rho= rhoCharge==0 ? rho0[0] : rhop[0];
      if(rho->iCharge()!=rhoCharge) rho = rhom[0];
      assert(abs(rhoCharge)<=3);
      mode->addChannel((PhaseSpaceChannel(phase),ires,rhoin,ires+1,omega,ires+1,iloc+j2,
			ires+2,rho,ires+2,iloc+j1,ires+3,iloc+j3,ires+3,iloc+j4));
      ++nchan; channelMap_[imode].push_back(nchan);
    }
  }
  else
    nchan+=18;
  // the rho channels cancel for -000 and ++--
  if(imode==0 || imode>3) {
    nchan+=16;
    return;
  }
  // rho channels
  for(unsigned int irho=0;irho<2;++irho) {
    tPDPtr rhoin = icharge==0 ? rho0[irho] : rhom[irho];
    if(resonance && rhoin!=resonance) continue;
    for(unsigned int irho1=0;irho1<2;++irho1) {
      for(unsigned int irho2=0;irho2<2;++irho2) {
  	int rho1Charge =  outgoing[j1-1]->iCharge()+outgoing[j3-1]->iCharge();
  	tPDPtr rho1= rho1Charge==0 ? rho0[irho1] : rhop[irho1];
  	if(rho1->iCharge()!=rho1Charge) rho1 = rhom[irho1];
  	assert(abs(rho1Charge)<=3);
  	int rho2Charge =  outgoing[j2-1]->iCharge()+outgoing[j4-1]->iCharge();
  	tPDPtr rho2= rho2Charge==0 ? rho0[irho2] : rhop[irho2];
  	if(rho2->iCharge()!=rho2Charge) rho2 = rhom[irho2];
  	assert(abs(rho2Charge)<=3);
  	mode->addChannel((PhaseSpaceChannel(phase),ires,rhoin,ires+1,rho1,ires+1,rho2,
  			  ires+2,iloc+j1,ires+2,iloc+j3,ires+3,iloc+j2,ires+3,iloc+j4));
	++nchan; channelMap_[imode].push_back(nchan);
  	assert(icharge==rho1Charge+rho2Charge);
	if(imode!=1) {
	  rho1Charge =  outgoing[j1-1]->iCharge()+outgoing[j4-1]->iCharge();
	  rho1= rho1Charge==0 ? rho0[irho1] : rhop[irho1];
	  if(rho1->iCharge()!=rho1Charge) rho1 = rhom[irho1];
	  assert(abs(rho1Charge)<=3);
	  rho2Charge =  outgoing[j2-1]->iCharge()+outgoing[j3-1]->iCharge();
	  rho2= rho2Charge==0 ? rho0[irho2] : rhop[irho2];
	  if(rho2->iCharge()!=rho2Charge) rho2 = rhom[irho2];
	  assert(abs(rho2Charge)<=3);
	  mode->addChannel((PhaseSpaceChannel(phase),ires,rhoin,ires+1,rho1,ires+1,rho2,
			    ires+2,iloc+j1,ires+2,iloc+j4,ires+3,iloc+j2,ires+3,iloc+j3));
	  ++nchan; channelMap_[imode].push_back(nchan);
	  assert(icharge==rho1Charge+rho2Charge);
	}
	else
	  ++nchan;
      }
    }
  }
}

// complete the construction of the decay mode for integration
bool FourPionCzyzCurrent::createMode(int icharge, tcPDPtr resonance,
				     FlavourInfo flavour,
				     unsigned int imode,PhaseSpaceModePtr mode,
				     unsigned int iloc,int ires,
				     PhaseSpaceChannel phase, Energy upp ) {
  // check the charge
  if((abs(icharge)!=3&&imode<2) ||
     (imode>=2&&icharge!=0)) return false;
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown) {
    if(flavour.I!=IsoSpin::IOne) return false;
  }
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
    case IsoSpin::I3Zero:
      if(imode==0) return false;
      break;
    case IsoSpin::I3One:
      if(imode==1 || icharge ==-3) return false;
      break;
    case IsoSpin::I3MinusOne:
      if(imode==1 || icharge ==3) return false;
      break;
    default:
      return false;
    }
  }
  // and other flavour
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero) return false;
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero      ) return false;
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero      ) return false;
  // check that the modes are kinematical allowed
  Energy min(ZERO);
  if(imode==0)                min =    mpip_+3.*mpi0_;
  else if(imode==1)           min = 3.*mpip_+   mpi0_;
  else if(imode==2||imode==3) min = 2.*mpip_+2.*mpi0_;
  else                        min = 4.*mpip_;
  if(min>upp) return false;
  // resonances we need
  // get the external particles
  int iq(0),ia(0);
  tPDVector outgoing = particles(icharge,imode,iq,ia);
  int nchan(-1);
  channelMap_[imode].clear();
  if(imode==0) {
    createChannels(imode,icharge,resonance,iloc,ires,outgoing,mode,phase,3,4,1,2,nchan);
    createChannels(imode,icharge,resonance,iloc,ires,outgoing,mode,phase,2,4,1,3,nchan);
    createChannels(imode,icharge,resonance,iloc,ires,outgoing,mode,phase,2,3,1,4,nchan);
  }
  else if(imode==1) {
    createChannels(imode,icharge,resonance,iloc,ires,outgoing,mode,phase,3,2,1,4,nchan);
    createChannels(imode,icharge,resonance,iloc,ires,outgoing,mode,phase,3,1,2,4,nchan);
  }
  // pi0 pi0 pi+ pi-
  else if(imode==2||imode==3) {
    createChannels(imode,icharge,resonance,iloc,ires,outgoing,mode,phase,3,4,1,2,nchan);
  }
  else {
    createChannels(imode,icharge,resonance,iloc,ires,outgoing,mode,phase,2,4,1,3,nchan);
    createChannels(imode,icharge,resonance,iloc,ires,outgoing,mode,phase,1,4,2,3,nchan);
    createChannels(imode,icharge,resonance,iloc,ires,outgoing,mode,phase,2,3,1,4,nchan);
    createChannels(imode,icharge,resonance,iloc,ires,outgoing,mode,phase,1,3,2,4,nchan);
  }
  return true;
}

// the particles produced by the current
tPDVector FourPionCzyzCurrent::particles(int icharge, unsigned int imode,
					 int,int) {
  tPDVector output(4);
  tPDPtr pi0=getParticleData(ParticleID::pi0);
  tPDPtr pip=getParticleData(ParticleID::piplus);
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
FourPionCzyzCurrent::current(tcPDPtr resonance,
			     FlavourInfo flavour,
			     const int imode, const int ichan,Energy & scale,
			     const tPDVector & outgoing,
			     const vector<Lorentz5Momentum> & momenta,
			     DecayIntegrator::MEOption) const {
  int icharge(0);
  for(tPDPtr out : outgoing) icharge+=out->iCharge();
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown) {
    if(flavour.I!=IsoSpin::IOne) return vector<LorentzPolarizationVectorE>();
  }
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
    case IsoSpin::I3Zero:
      if(imode==0) return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3One:
      if(imode==1 || icharge ==-3) return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3MinusOne:
      if(imode==1 || icharge ==3) return vector<LorentzPolarizationVectorE>();
      break;
    default:
      return vector<LorentzPolarizationVectorE>();
    }
  }
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero) return vector<LorentzPolarizationVectorE>();
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero      ) return vector<LorentzPolarizationVectorE>();
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero      ) return vector<LorentzPolarizationVectorE>();
  useMe();
  // the momenta of the particles
  Lorentz5Momentum q1(momenta[0]);
  Lorentz5Momentum q2(momenta[1]);
  Lorentz5Momentum q3(momenta[2]);
  Lorentz5Momentum q4(momenta[3]);
  Lorentz5Momentum Q(q1+q2+q3+q4);
  Q.rescaleMass();
  scale = Q.mass();
  LorentzVector<complex<InvEnergy> > output;
  assert(ichan<int(channelMap_[imode].size()));
  int ichannelB = ichan<0 ? -1 : channelMap_[imode][ichan];
  if(imode==0) {
    if(ichannelB<51)
      output  += baseCurrent(Q.mass2(),resonance,ichannelB   ,Q,q3,q4,q1,q2);
    if(ichannelB<0 || (ichannelB>=67&&ichannelB<133))
       output += baseCurrent(Q.mass2(),resonance,ichannelB-67,Q,q2,q4,q1,q3);
    if(ichannelB<0 ||  ichannelB>=133) 
       output += baseCurrent(Q.mass2(),resonance,ichannelB-133,Q,q2,q3,q1,q4);
    output *= sqrt(1./3.);
  }
  else if(imode==1) {
    if(ichannelB<117)
      output += baseCurrent(Q.mass2(),resonance,ichannelB,Q,q3,q2,q1,q4);
    if(ichannelB<0||ichannelB>=67)
      output += baseCurrent(Q.mass2(),resonance,ichannelB-67,Q,q3,q1,q2,q4);
  }
  else if(imode==2||imode==3) {
    output = baseCurrent(Q.mass2(),resonance,ichannelB,Q,q3,q4,q1,q2);
  }
  else if(imode==4||imode==5) {
    if(ichannelB<67)
      output += baseCurrent(Q.mass2(),resonance,ichannelB    ,Q,q2,q4,q1,q3);
    if(ichannelB<0 || (ichannelB>=67&&ichannelB<133))
      output += baseCurrent(Q.mass2(),resonance,ichannelB-67 ,Q,q1,q4,q2,q3);
    if(ichannelB<0 || (ichannelB>=133&&ichannelB<200)) 
      output += baseCurrent(Q.mass2(),resonance,ichannelB-133,Q,q2,q3,q1,q4);
    if(ichannelB<0 ||  ichannelB>=200)
      output += baseCurrent(Q.mass2(),resonance,ichannelB-200,Q,q1,q3,q2,q4);
  }
  return  vector<LorentzPolarizationVectorE>(1,output*Q.mass2());
}
   
bool FourPionCzyzCurrent::accept(vector<int> id) {
  bool allowed(false);
  // check four products
  if(id.size()!=4){return false;}
  int npiminus=0,npiplus=0,npi0=0;
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(id[ix]==ParticleID:: piplus)      ++npiplus;
    else if(id[ix]==ParticleID::piminus) ++npiminus;
    else if(id[ix]==ParticleID::pi0)     ++npi0;
  }
  if(npiminus==2&&npiplus==1&&npi0==1)      allowed=true;
  else if(npiminus==1&&npi0==3)             allowed=true;
  else if(npiplus==2&&npiminus==1&&npi0==1) allowed=true;
  else if(npiplus==1&&npi0==3)              allowed=true;
  else if(npiplus==2&&npiminus==2)          allowed=true;
  else if(npiplus==1&&npiminus==1&&npi0==2) allowed=true;
  return allowed;
}

// the decay mode
unsigned int FourPionCzyzCurrent::decayMode(vector<int> idout) {
  unsigned int npi(0);
  for(unsigned int ix=0;ix<idout.size();++ix) {
    if(abs(idout[ix])==ParticleID::piplus) ++npi;
  }
  if(npi==1) return 0;
  else if(npi==2) return 2;
  else if(npi==3) return 1;
  else return 4;
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
  WeakCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
  		    << fullName() << "\";" << endl;
}

LorentzVector<complex<InvEnergy> > FourPionCzyzCurrent::
baseCurrent(Energy2 Q2, tcPDPtr resonance,const int ichan,
	    const Lorentz5Momentum & Q , const Lorentz5Momentum & q1,
	    const Lorentz5Momentum & q2, const Lorentz5Momentum & q3,
	    const Lorentz5Momentum & q4) const {
  // check the resonance
  int ires0(-1);
  if(resonance) {
    switch(resonance->id()/1000) {
    case 0:
      ires0=0;
      break;
    case 100:
      ires0=1;
      break;
    case 30 :
      ires0=2;
      break;
    default:
      assert(false);
    }
  }
  // various dot products we'll need
  Energy2 m12 = sqr(q1.mass()), m22 = sqr(q2.mass());
  Energy2 m32 = sqr(q3.mass()), m42 = sqr(q4.mass());
  Energy2 Qq1 = Q*q1, Qq2 = Q*q2, Qq3 = Q*q3, Qq4 = Q*q4;
  Energy2 Qm32= Q2-2.*Qq3+m32;
  Energy2 Qm42= Q2-2.*Qq4+m42;
  Energy2 q1q2 = q1*q2, q1q3 = q1*q3, q1q4 = q1*q4;
  Energy2 q2q3 = q2*q3, q2q4 = q2*q4, q3q4 = q3*q4;
  // coefficients of the momenta
  complex<InvEnergy2> c1(ZERO),c2(ZERO),c34(ZERO),cq(ZERO),c3(ZERO),c4(ZERO);
  // first the a_1 terms from A.3 0804.0359 (N.B. sign due defns)
  // common coefficent
  if(ichan<24) {
    int ires1(-1),ires2(-1),iterm(-1);
    if(ichan>=0) {
      ires1  = ichan/8;
      ires2  = (ichan/4)%2;
      iterm = ichan%4;
    }
    if(ires0>=0 && ires1<0) ires1=ires0;
    complex<InvEnergy2> a1_fact;
    if(ires1<0) {
      a1_fact = -c_a1_*
  	Resonance::F_rho(Q2,beta_a1_,rhoMasses_Frho_,
  			 rhoWidths_Frho_,mpip_,mpip_);
    }
    else {
      if(ires0<0 || ires0==ires1) {
  	a1_fact = -c_a1_*beta_a1_[ires1]*
  	  Resonance::BreitWignerPWave(Q2,rhoMasses_Frho_[ires1],rhoWidths_Frho_[ires1],mpip_,mpip_)
  	  /std::accumulate(beta_a1_.begin(),beta_a1_.end(),0.);
      }
      else
  	a1_fact=ZERO;
    }
    // first 2 terms
    if(iterm<2) {
      Complex bw_a1_Qm3 = Resonance::BreitWignera1(Qm32,a1Mass_,a1Width_);
      // first term
      if(iterm<=0) {
  	Complex Brhoq1q4(0.);
  	if(ires2<0) {
  	  Brhoq1q4 = bw_a1_Qm3*
  	  Resonance::F_rho(m12+m42+2.*q1q4,beta_B_,rhoMasses_,rhoWidths_,mpip_,mpip_);
  	}
  	else {
  	  Brhoq1q4 = bw_a1_Qm3*beta_B_[ires2]*
  	    Resonance::BreitWignerPWave(m12+m42+2.*q1q4,rhoMasses_[ires2],rhoWidths_[ires2],mpip_,mpip_)/
  	    std::accumulate(beta_B_.begin(),beta_B_.end(),0.);
  	}
  	double dot1  = (q1q2-q2q4)/Qm32;
  	double dot1B = (Qq1-Qq4)/Q2;
  	double dot1C = Qq3/Q2*dot1;
  	// coefficients of the momenta to construct the current
  	c1 += 0.5*a1_fact*Brhoq1q4*( 3.-dot1);
  	c2 += 0.5*a1_fact*Brhoq1q4*( 1.-dot1);
  	c34+= 0.5*a1_fact*Brhoq1q4*( 1.+dot1);
  	cq += 0.5*a1_fact*Brhoq1q4*(-1.+dot1-2.*dot1B-2.*dot1C);
      }
      // second term
      if(iterm<0||iterm==1) {
      	Complex Brhoq2q4(0.);
      	if(ires2<0) {
      	  Brhoq2q4= bw_a1_Qm3*
      	    Resonance::F_rho(m12+m42+2.*q2q4,beta_B_,rhoMasses_,rhoWidths_,mpip_,mpip_);
      	}
      	else {
      	  Brhoq2q4= bw_a1_Qm3*beta_B_[ires2]*
      	    Resonance::BreitWignerPWave(m12+m42+2.*q2q4,rhoMasses_[ires2],rhoWidths_[ires2],mpip_,mpip_)/
      	    std::accumulate(beta_B_.begin(),beta_B_.end(),0.);
      	}
      	double dot2  = (q1q2-q1q4)/Qm32;
      	double dot2B = (Qq2-Qq4)/Q2;
      	double dot2C = Qq3/Q2*dot2;
      	// coefficients of the momenta to construct the current
      	// a_1 terms
      	c1 += 0.5*a1_fact*Brhoq2q4*( 1.-dot2);
      	c2 += 0.5*a1_fact*Brhoq2q4*( 3.-dot2);
      	c34+= 0.5*a1_fact*Brhoq2q4*( 1.+dot2);
      	cq += 0.5*a1_fact*Brhoq2q4*(-1.+dot2-2.*dot2B-2.*dot2C);
      }
    }
    // second 2 terms
    if(iterm<0||iterm>=2) {
      Complex bw_a1_Qm4 = Resonance::BreitWignera1(Qm42,a1Mass_,a1Width_);
      // third term
      if(iterm<0||iterm==2) {
      	Complex Brhoq1q3(0.);
      	if(ires2<0) {
      	  Brhoq1q3 = bw_a1_Qm4*
      	    Resonance::F_rho(m12+m32+2.*q1q3,beta_B_,rhoMasses_,rhoWidths_,mpip_,mpip_);
      	}
      	else {
      	  Brhoq1q3 = bw_a1_Qm4*beta_B_[ires2]*
      	    Resonance::BreitWignerPWave(m12+m32+2.*q1q3,rhoMasses_[ires2],rhoWidths_[ires2],mpip_,mpip_)/
      	    std::accumulate(beta_B_.begin(),beta_B_.end(),0.);
      	}
      	double dot3  = (q1q2-q2q3)/Qm42;
      	double dot3B = (Qq1-Qq3)/Q2;
      	double dot3C = Qq4/Q2*dot3;
      	// coefficients of the momenta to construct the current
      	// a_1 terms
      	c1 += 0.5*a1_fact*Brhoq1q3*(-3.+dot3);
      	c2 += 0.5*a1_fact*Brhoq1q3*(-1.+dot3);
      	c34+= 0.5*a1_fact*Brhoq1q3*( 1.+dot3);
      	cq += 0.5*a1_fact*Brhoq1q3*( 1.-dot3+2.*dot3B+2.*dot3C);
      }
      // fourth term
      if(iterm<0||iterm==3) {
  	Complex Brhoq2q3(0.);
  	if(ires2<0) {
  	  Brhoq2q3 = bw_a1_Qm4*
  	  Resonance::F_rho(m12+m32+2.*q2q3,beta_B_,rhoMasses_,rhoWidths_,mpip_,mpip_);
  	}
  	else {
  	  Brhoq2q3 = bw_a1_Qm4*beta_B_[ires2]*
  	  Resonance::BreitWignerPWave(m12+m32+2.*q2q3,rhoMasses_[ires2],rhoWidths_[ires2],mpip_,mpip_)/
  	    std::accumulate(beta_B_.begin(),beta_B_.end(),0.);
  	}
  	double dot4  = (q1q2-q1q3)/Qm42;
  	double dot4B = (Qq2-Qq3)/Q2;
  	double dot4C = Qq4/Q2*dot4;
  	// coefficients of the momenta to construct the current
  	// a_1 terms
  	c1 += 0.5*a1_fact*Brhoq2q3*(-1.+dot4);
  	c2 += 0.5*a1_fact*Brhoq2q3*(-3.+dot4);
  	c34+= 0.5*a1_fact*Brhoq2q3*( 1.+dot4);
  	cq += 0.5*a1_fact*Brhoq2q3*( 1.-dot4+2.*dot4B+2.*dot4C);
      }
    }
  }
  // f_0
  if(ichan<0 || (ichan>=24 && ichan<33) ) {
    int ires1(-1),ires2(-2);
    if(ichan>0) {
      ires1 = (ichan-24)/3;
      ires2 = (ichan-24)%3;
    }
    Complex rho1(0.);
    if(ires0>=0 && ires1<0) ires1=ires0;
    if(ires1<0) {
      rho1 = Resonance::F_rho(Q2,beta_f0_,rhoMasses_Frho_,rhoWidths_Frho_,mpip_,mpip_);
    }
    else {
      if(ires0<0 || ires0==ires1) {
  	rho1 = beta_f0_[ires1]*Resonance::BreitWignerPWave(Q2,rhoMasses_Frho_[ires1],rhoWidths_Frho_[ires1],mpip_,mpip_)/
  	  std::accumulate(beta_f0_.begin(),beta_f0_.end(),0.);
      }
      else
  	rho1=0.;
    }
    Complex rho2(0.);
    if(ires2<0) {
      rho2 = Resonance::F_rho(m32+m42+2.*q3q4,beta_bar_,rhoMasses_,rhoWidths_,mpip_,mpip_);
    }
    else {
      rho2 = beta_bar_[ires2]*Resonance::BreitWignerPWave(m32+m42+2.*q3q4,rhoMasses_[ires2],rhoWidths_[ires2],mpip_,mpip_)/
  	std::accumulate(beta_bar_.begin(),beta_bar_.end(),0.);
    }
    complex<InvEnergy2> f0fact = c_f0_*rho1*rho2*
      Resonance::BreitWignerSWave(m12+m22+2.*q1q2,f0Mass_,f0Width_,mpip_,mpip_);
    // add contribution to the coefficients
    c34 -= f0fact;
    cq  += f0fact*(Qq3-Qq4)/Q2;
  }
  // omega
  if(ichan<0 || (ichan>=33&&ichan<=50) ) {
    int ires1(-1),iterm(-1);
    if(ichan>0) {
      ires1 = (ichan-33)/6;
      iterm = (ichan-33)%6;
    }
    Complex rho1(0.);
    if(ires0>=0 && ires1<0) ires1=ires0;
    if(ires1<0) {
      rho1 = Resonance::F_rho(Q2,beta_omega_,rhoMasses_Frho_,rhoWidths_Frho_,mpip_,mpip_);
    }
    else {
      if(ires0<0 || ires0==ires1) {
	rho1 = beta_omega_[ires1]*Resonance::BreitWignerPWave(Q2,rhoMasses_Frho_[ires1],rhoWidths_Frho_[ires1],mpip_,mpip_)/
	  std::accumulate(beta_omega_.begin(),beta_omega_.end(),0.);
      }
      else
	rho1=0.;
    }
    complex<ThePEG::Qty<std::ratio<0,1>, std::ratio<-6,1>, std::ratio<0,1> > >
      wfact = 2.*c_omega_*g_omega_pi_rho_*g_rho_pi_pi_*rho1;
    if(iterm<3) {
      Complex Hterm(0.);
      if(iterm<0) {
	Hterm = Resonance::H(rhoMasses_[0],rhoWidths_[0],m22+2.*q2q3+m32,m22+2.*q2q4+m42,
			     m32+2.*q3q4+m42,mpip_,mpip_);
      }
      else if(iterm==0)
	Hterm = Resonance::BreitWignerPWave(m22+2.*q2q3+m32,rhoMasses_[0],rhoWidths_[0],mpip_,mpip_);
      else if(iterm==1)
	Hterm = Resonance::BreitWignerPWave(m22+2.*q2q4+m42,rhoMasses_[0],rhoWidths_[0],mpip_,mpip_);
      else if(iterm==2)
	Hterm = Resonance::BreitWignerPWave(m32+2.*q3q4+m42,rhoMasses_[0],rhoWidths_[0],mpip_,mpip_);
      else
	assert(false);
      complex<ThePEG::Qty<std::ratio<0,1>, std::ratio<-6,1>, std::ratio<0,1> > >
	bw_omega_1 = wfact*Resonance::BreitWignerFW(m12-2.*Qq1+Q2,omegaMass_,omegaWidth_)*Hterm;
      c2 -=  bw_omega_1*(q1q4*Qq3-q1q3*Qq4);
      c3 -=  bw_omega_1*(q1q2*Qq4-q1q4*Qq2);
      c4 -=  bw_omega_1*(q1q3*Qq2-q1q2*Qq3);
    }
    if(iterm<0||iterm>=3) {
      Complex Hterm(0.);
      if(iterm<0) {
	Hterm = Resonance::H(rhoMasses_[0],rhoWidths_[0],m12+2.*q1q3+m32,m12+2.*q1q4+m42,
			     m32+2.*q3q4+m42,mpip_,mpip_);
      }
      else if(iterm==3)
	Hterm = Resonance::BreitWignerPWave(m12+2.*q1q3+m32,rhoMasses_[0],rhoWidths_[0],mpip_,mpip_);
      else if(iterm==4)
	Hterm = Resonance::BreitWignerPWave(m12+2.*q1q4+m42,rhoMasses_[0],rhoWidths_[0],mpip_,mpip_);
      else if(iterm==5)
	Hterm = Resonance::BreitWignerPWave(m32+2.*q3q4+m42,rhoMasses_[0],rhoWidths_[0],mpip_,mpip_);
      else
	assert(false);
      complex<ThePEG::Qty<std::ratio<0,1>, std::ratio<-6,1>, std::ratio<0,1> > >
	bw_omega_2 = wfact*Resonance::BreitWignerFW(m22-2.*Qq2+Q2,omegaMass_,omegaWidth_)*Hterm;
      c1 -= bw_omega_2*(q2q4*Qq3-q2q3*Qq4);
      c3 -= bw_omega_2*(q1q2*Qq4-q2q4*Qq1);
      c4 -= bw_omega_2*(q2q3*Qq1-q1q2*Qq3);
    }
  }
  // the rho term
  LorentzVector<complex<InvEnergy> > v_rho;
  if(ichan<0||ichan>=51) {
    int ires1(-1),ires2(-1),ires3(-1),iterm(-1);
    if(ichan>0) {
      ires1 = (ichan-51)/8;
      ires2 = ((ichan-51)/4)%2;
      ires3 = ((ichan-51)/2)%2;
      iterm = (ichan-51)%2;
    }
    // prefactor
    complex<InvEnergy2> rho1;
    if(ires0>=0 && ires1<0) ires1=ires0;
    if(ires1<0) {
      rho1 = Resonance::BreitWignerDiff(Q2,rhoMasses_[0],rhoWidths_[0],
   					rhoMasses_[1],rhoWidths_[1],mpip_,mpip_);
    }
    else {
      if(ires0<0 || ires0==ires1) {
   	if(ires1==0)
   	  rho1 =  Resonance::BreitWignerPWave(Q2,rhoMasses_[0],rhoWidths_[0],mpip_,mpip_)/sqr(rhoMasses_[0]);
   	else if(ires1==1)
   	  rho1 = -Resonance::BreitWignerPWave(Q2,rhoMasses_[1],rhoWidths_[1],mpip_,mpip_)/sqr(rhoMasses_[1]);
   	else
   	  rho1 = ZERO;
      }
      else
   	rho1 = ZERO;
    }
    Complex pre_rho = c_rho_*pow(g_rho_pi_pi_,3)*g_rho_gamma_*rho1;
    complex<InvEnergy2> BW12_q1q3_A(ZERO),BW12_q1q4_A(ZERO),BW12_q2q3_B(ZERO),BW12_q2q4_B(ZERO);
    if(ires2<0) {
      BW12_q1q3_A = Resonance::BreitWignerDiff(m12+m32+2.*q1q3,rhoMasses_[0],rhoWidths_[0],
  					       rhoMasses_[1],rhoWidths_[1],mpip_,mpip_);
      BW12_q1q4_A = Resonance::BreitWignerDiff(m12+m42+2.*q1q4,rhoMasses_[0],rhoWidths_[0],
  					       rhoMasses_[1],rhoWidths_[1],mpip_,mpip_);
      BW12_q2q3_B = Resonance::BreitWignerDiff(m22+m32+2.*q2q3,rhoMasses_[0],rhoWidths_[0],
  					       rhoMasses_[1],rhoWidths_[1],mpip_,mpip_);
      
      BW12_q2q4_B = Resonance::BreitWignerDiff(m22+m42+2.*q2q4,rhoMasses_[0],rhoWidths_[0],
  					       rhoMasses_[1],rhoWidths_[1],mpip_,mpip_);
      
    }
    else {
      if(ires2==0) {
  	BW12_q1q3_A = Resonance::BreitWignerPWave(m12+m32+2.*q1q3,rhoMasses_[0],rhoWidths_[0],mpip_,mpip_)/sqr(rhoMasses_[0]);
  	BW12_q1q4_A = Resonance::BreitWignerPWave(m12+m42+2.*q1q4,rhoMasses_[0],rhoWidths_[0],mpip_,mpip_)/sqr(rhoMasses_[0]);
      }
      else {
  	BW12_q1q3_A = -Resonance::BreitWignerPWave(m12+m32+2.*q1q3,rhoMasses_[1],rhoWidths_[1],mpip_,mpip_)/sqr(rhoMasses_[1]);
  	BW12_q1q4_A = -Resonance::BreitWignerPWave(m12+m42+2.*q1q4,rhoMasses_[1],rhoWidths_[1],mpip_,mpip_)/sqr(rhoMasses_[1]);
      }
    }
    if(ires3>=0) {
      if(ires3==0) {
  	BW12_q2q3_B = Resonance::BreitWignerPWave(m22+m32+2.*q2q3,rhoMasses_[0],rhoWidths_[0],mpip_,mpip_)/sqr(rhoMasses_[0]);
  	BW12_q2q4_B = Resonance::BreitWignerPWave(m22+m42+2.*q2q4,rhoMasses_[0],rhoWidths_[0],mpip_,mpip_)/sqr(rhoMasses_[0]);
      }
      else {
  	BW12_q2q3_B = -Resonance::BreitWignerPWave(m22+m32+2.*q2q3,rhoMasses_[1],rhoWidths_[1],mpip_,mpip_)/sqr(rhoMasses_[1]);
  	BW12_q2q4_B = -Resonance::BreitWignerPWave(m22+m42+2.*q2q4,rhoMasses_[1],rhoWidths_[1],mpip_,mpip_)/sqr(rhoMasses_[1]);
      }
    }
    // now the various terms
    if(iterm<=0) {
      Energy2 d1 = Qq2 - Qq4 + 2.*q2q3 - 2.*q3q4;
      v_rho += pre_rho*BW12_q1q3_A*(BW12_q2q4_B*d1 + 2.)*q1;
    }
    if(iterm<=0) {
      Energy2 d5 = Qq1 - Qq3 + 2.*q1q2 - 2.*q2q3;
      v_rho += pre_rho*BW12_q2q4_B*(BW12_q1q3_A*d5 + 2.)*q4;
    }
    if(iterm<0||iterm==1) {
      Energy2 d2 = Qq2 - Qq3 + 2.*q2q4 - 2.*q3q4;
      v_rho -= pre_rho*BW12_q1q4_A*(BW12_q2q3_B*d2 + 2.)*q1;
    }
    if(iterm<0||iterm==1) {
      Energy2 d7 = Qq1 - Qq4 + 2.*q1q2 - 2.*q2q4;
      v_rho -= pre_rho*BW12_q2q3_B*(BW12_q1q4_A*d7 + 2.)*q3;
    }
    if(iterm<0||iterm==1) {
      Energy2 d3 = Qq1 - Qq4 + 2.*q1q3 - 2.*q3q4;
      v_rho += pre_rho*BW12_q2q3_B*(BW12_q1q4_A*d3 + 2.)*q2;
    }
    if(iterm<0||iterm==1) {
      Energy2 d6 = Qq2 - Qq3 + 2.*q1q2 - 2.*q1q3;
      v_rho += pre_rho*BW12_q1q4_A*(BW12_q2q3_B*d6 + 2.)*q4;
    }
    if(iterm<=0) {
      Energy2 d4 = Qq1 - Qq3 + 2.*q1q4 - 2.*q3q4;
      v_rho -= pre_rho*BW12_q2q4_B*(BW12_q1q3_A*d4 + 2.)*q2;
    }
    if(iterm<=0) {
      Energy2 d8 = Qq2 - Qq4 + 2.*q1q2 - 2.*q1q4;  
      v_rho -= pre_rho*BW12_q1q3_A*(BW12_q2q4_B*d8 + 2.)*q3;
    }
    complex<InvEnergy2> vdot = (Q*v_rho)/Q2;  
    v_rho  = -v_rho + vdot*Q;
  }
  // put everything together
  return c1*q1+c2*q2+(c3+c34)*q3+(c4-c34)*q4+cq*Q+v_rho;
}
