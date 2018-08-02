// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ThreePionCzyzCurrent class.
//

#include "ThreePionCzyzCurrent.h"
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

static const InvEnergy3 InvGeV3 = pow<-3,1>(GeV);

ThreePionCzyzCurrent::ThreePionCzyzCurrent()
  : mpip_(140*MeV), mpi0_(140*MeV) {
  // parameters for I=0
  // masses and widths
  rhoMasses_ = {0.77609*GeV,1.465*GeV,1.7  *GeV};
  rhoWidths_ = {0.14446*GeV,0.31 *GeV,0.235*GeV};
  omegaMasses_ = {782.4*MeV,1375*MeV,1631*MeV};
  omegaWidths_ = {8.69 *MeV, 250*MeV, 245*MeV};
  phiMass_     = 1019.24*MeV;
  phiWidth_    = 4.14*MeV;
  // couplings
  coup_I0_ = {18.20*InvGeV3,-0.87*InvGeV3,-0.77*InvGeV3,
	      -1.12*InvGeV3,-0.72*InvGeV3,-0.59*InvGeV3};
  // parameters for I=1
  rhoMasses_I1_ = {0.77609*GeV,1.7*GeV };
  rhoWidths_I1_ = {0.14446*GeV,0.26*GeV};
  omegaMass_I1_ = 782.59*MeV;
  omegaWidth_I1_= 8.49*MeV;
  // couplings
  sigma_ = -0.1;
  GW_pre_    = 1.55/sqrt(2.)*12.924*0.266/GeV;
  g_omega_pi_pi_ = 0.185;
  addDecayMode(1,-1);
  addDecayMode(2,-2);
  setInitialModes(2);
}

IBPtr ThreePionCzyzCurrent::clone() const {
  return new_ptr(*this);
}

IBPtr ThreePionCzyzCurrent::fullclone() const {
  return new_ptr(*this);
}

void ThreePionCzyzCurrent::persistentOutput(PersistentOStream & os) const {
  os << ounit(rhoMasses_,GeV) << ounit(rhoWidths_,GeV)
     << ounit(mpip_,GeV) << ounit(mpi0_,GeV)
     << ounit(omegaMasses_,GeV) << ounit(omegaWidths_,GeV)
     << ounit(phiMass_,GeV) <<  ounit(phiWidth_,GeV) << ounit(coup_I0_,InvGeV3)
     << ounit(rhoMasses_I1_,GeV) << ounit(rhoWidths_I1_,GeV)
     << ounit(omegaMass_I1_,GeV) << ounit(omegaWidth_I1_,GeV)
     << sigma_ << ounit(GW_pre_,1./GeV) << g_omega_pi_pi_ << ounit(GW_,GeV);
}

void ThreePionCzyzCurrent::persistentInput(PersistentIStream & is, int) {
  is >> iunit(rhoMasses_,GeV) >> iunit(rhoWidths_,GeV)
     >> iunit(mpip_,GeV) >> iunit(mpi0_,GeV)
     >> iunit(omegaMasses_,GeV) >> iunit(omegaWidths_,GeV)
     >> iunit(phiMass_,GeV) >>  iunit(phiWidth_,GeV) >> iunit(coup_I0_,InvGeV3)
     >> iunit(rhoMasses_I1_,GeV) >> iunit(rhoWidths_I1_,GeV)
     >> iunit(omegaMass_I1_,GeV) >> iunit(omegaWidth_I1_,GeV)
     >> sigma_ >> iunit(GW_pre_,1./GeV) >> g_omega_pi_pi_ >> iunit(GW_,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<ThreePionCzyzCurrent,WeakDecayCurrent>
describeHerwigThreePionCzyzCurrent("Herwig::ThreePionCzyzCurrent",
				     "HwWeakCurrents.so");

void ThreePionCzyzCurrent::Init() {

  static ClassDocumentation<ThreePionCzyzCurrent> documentation
    ("The ThreePionCzyzCurrent class is designed to implement "
     "the three pion current for e+e- collisions from Eur.Phys.J. C47 (2006) 617-624",
     "The current from \\cite{Czyz:2005as} was used for $\\pi^+\\pi^-\\pi^0$",
     "\\bibitem{Czyz:2005as}\n"
     "H.~Czyz, A.~Grzelinska, J.~H.~Kuhn and G.~Rodrigo,\n"
     "%``Electron-positron annihilation into three pions and the radiative return,''\n"
     "Eur.\\ Phys.\\ J.\\ C {\\bf 47} (2006) 617\n"
     "doi:10.1140/epjc/s2006-02614-7\n"
     "[hep-ph/0512180].\n"
     "%%CITATION = doi:10.1140/epjc/s2006-02614-7;%%\n"
     "%32 citations counted in INSPIRE as of 01 Aug 2018\n"
     );

  static ParVector<ThreePionCzyzCurrent,Energy> interfaceRhoMassesI0
    ("RhoMassesI0",
     "The rho masses for the I=0 part of the current",
     &ThreePionCzyzCurrent::rhoMasses_, GeV, -1, 0.766*GeV, 0*GeV, 0*GeV,
     false, false, Interface::nolimits);
  
  static ParVector<ThreePionCzyzCurrent,Energy> interfaceRhoWidthsI0
    ("RhoWidthsI0",
     "The rho masses for the I=0 part of the current",
     &ThreePionCzyzCurrent::rhoWidths_, GeV, -1, 0.766*GeV, 0*GeV, 0*GeV,
     false, false, Interface::nolimits);
  
  static ParVector<ThreePionCzyzCurrent,Energy> interfaceOmegaMassesI0
    ("OmegaMassesI0",
     "The omega masses for the I=0 part of the current",
     &ThreePionCzyzCurrent::omegaMasses_, GeV, -1, 0.766*GeV, 0*GeV, 0*GeV,
     false, false, Interface::nolimits);
  
  static ParVector<ThreePionCzyzCurrent,Energy> interfaceOmegaWidthsI0
    ("OmegaWidthsI0",
     "The omega masses for the I=0 part of the current",
     &ThreePionCzyzCurrent::omegaWidths_, GeV, -1, 0.766*GeV, 0*GeV, 0*GeV,
     false, false, Interface::nolimits);
 
  static Parameter<ThreePionCzyzCurrent,Energy> interfacePhiMass
    ("PhiMass",
     "The mass of the phi meson",
     &ThreePionCzyzCurrent::phiMass_, GeV, 1.0*GeV, 0*GeV, 0*GeV,
     false, false, Interface::nolimits);
  
  static Parameter<ThreePionCzyzCurrent,Energy> interfacePhiWidth
    ("PhiWidth",
     "The width of the phi meson",
     &ThreePionCzyzCurrent::phiWidth_, GeV, 1.0*GeV, 0*GeV, 0*GeV,
     false, false, Interface::nolimits);
  
  static ParVector<ThreePionCzyzCurrent,InvEnergy3> interfaceCouplingsI0
    ("CouplingsI0",
     "The couplings for the I=0 component",
     &ThreePionCzyzCurrent::coup_I0_, InvGeV3, -1, 1.0*InvGeV3, 0*InvGeV3, 0*InvGeV3,
     false, false, Interface::nolimits);

  static ParVector<ThreePionCzyzCurrent,Energy> interfaceRhoMassesI1
    ("RhoMassesI1",
     "The rho masses for the I=1 part of the current",
     &ThreePionCzyzCurrent::rhoMasses_I1_, GeV, -1, 0.766*GeV, 0*GeV, 0*GeV,
     false, false, Interface::nolimits);
  
  static ParVector<ThreePionCzyzCurrent,Energy> interfaceRhoWidthsI1
    ("RhoWidthsI1",
     "The rho masses for the I=0 part of the current",
     &ThreePionCzyzCurrent::rhoWidths_I1_, GeV, -1, 0.766*GeV, 0*GeV, 0*GeV,
     false, false, Interface::nolimits);
   
  static Parameter<ThreePionCzyzCurrent,Energy> interfaceOmegaMass
    ("OmegaMass",
     "The mass of the omega meson",
     &ThreePionCzyzCurrent::omegaMass_I1_, GeV, 0.78259*GeV, 0*GeV, 0*GeV,
     false, false, Interface::nolimits);
  
  static Parameter<ThreePionCzyzCurrent,Energy> interfaceOmegaWidth
    ("OmegaWidth",
     "The width of the omega meson",
     &ThreePionCzyzCurrent::omegaWidth_I1_, GeV, 0.00849*GeV, 0*GeV, 0*GeV,
     false, false, Interface::nolimits);
  
  static Parameter<ThreePionCzyzCurrent,double> interfacesigma
    ("sigma",
     "The sigma parameter for the I=1 component",
     &ThreePionCzyzCurrent::sigma_, -0.1, -10., 10.0,
     false, false, Interface::limited);

  static Parameter<ThreePionCzyzCurrent,InvEnergy> interfaceGWPrefactor
    ("GWPrefactor",
     "The prefactor for the G omega coupling",
     &ThreePionCzyzCurrent::GW_pre_, 1./GeV,  1.55/sqrt(2.)*12.924*0.266/GeV, 0./GeV, 1e5/GeV,
     false, false, Interface::limited);

  static Parameter<ThreePionCzyzCurrent,double> interfaceg_omega_pipi
    ("g_omega_pipi",
     "The coupling of the omega meson to two pions",
     &ThreePionCzyzCurrent::g_omega_pi_pi_, 0.185, 0.0, 1.0,
     false, false, Interface::limited);

}

void ThreePionCzyzCurrent::doinit() {
  WeakDecayCurrent::doinit();
  GW_ = GW_pre_*sqr(rhoMasses_I1_[0])*g_omega_pi_pi_;
  mpip_ = getParticleData(211)->mass();
  mpi0_ = getParticleData(111)->mass();
}

// complete the construction of the decay mode for integration
bool ThreePionCzyzCurrent::createMode(int icharge, unsigned int imode,
				      DecayPhaseSpaceModePtr mode,
				      unsigned int iloc,unsigned int ires,
				      DecayPhaseSpaceChannelPtr phase,
				      Energy upp) {
  if(imode>2 || icharge != 0) return false;
  tPDPtr pip = getParticleData(ParticleID::piplus);
  tPDPtr pim = getParticleData(ParticleID::piminus);
  tPDPtr pi0 = getParticleData(ParticleID::pi0);
  if(2*pip->mass()+pi0->mass()>upp) return false;
  // resonaces we need
  tPDPtr omega[4] = {getParticleData( 223),getParticleData( 100223),getParticleData( 30223),
		      getParticleData( 333)};
  tPDPtr rho0[3]  = {getParticleData( 113),getParticleData( 100113),getParticleData( 30113)};
  tPDPtr rhop[3]  = {getParticleData( 213),getParticleData( 100213),getParticleData( 30213)};
  tPDPtr rhom[3]  = {getParticleData(-213),getParticleData(-100213),getParticleData(-30213)};
  DecayPhaseSpaceChannelPtr newchannel;
  // I=0 channels
  // omega/omega -> rho pi
  for(unsigned int ix=0;ix<4;++ix) {
    newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
    newchannel->addIntermediate(omega[ix],0,0.0,-ires-1,iloc);
    newchannel->addIntermediate(rhom[0]  ,0,0.0, iloc+1,iloc+2);
    mode->addChannel(newchannel);
    newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
    newchannel->addIntermediate(omega[ix],0,0.0,-ires-1,iloc+1);
    newchannel->addIntermediate(rhop[0]  ,0,0.0, iloc  ,iloc+2);
    mode->addChannel(newchannel);
    newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
    newchannel->addIntermediate(omega[ix],0,0.0,-ires-1,iloc+2);
    newchannel->addIntermediate(rho0[0]  ,0,0.0, iloc  ,iloc+1);
    mode->addChannel(newchannel);
  }
  // phi rho 1450
  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
  newchannel->addIntermediate(omega[3],0,0.0,-ires-1,iloc);
  newchannel->addIntermediate(rhom[1]  ,0,0.0, iloc+1,iloc+2);
  mode->addChannel(newchannel);
  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
  newchannel->addIntermediate(omega[3],0,0.0,-ires-1,iloc+1);
  newchannel->addIntermediate(rhop[1]  ,0,0.0, iloc  ,iloc+2);
  mode->addChannel(newchannel);
  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
  newchannel->addIntermediate(omega[3],0,0.0,-ires-1,iloc+2);
  newchannel->addIntermediate(rho0[1]  ,0,0.0, iloc  ,iloc+1);
  mode->addChannel(newchannel);
  // omega 1650 rho 1700
  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
  newchannel->addIntermediate(omega[2],0,0.0,-ires-1,iloc);
  newchannel->addIntermediate(rhom[2]  ,0,0.0, iloc+1,iloc+2);
  mode->addChannel(newchannel);
  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
  newchannel->addIntermediate(omega[2],0,0.0,-ires-1,iloc+1);
  newchannel->addIntermediate(rhop[2]  ,0,0.0, iloc  ,iloc+2);
  mode->addChannel(newchannel);
  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
  newchannel->addIntermediate(omega[2],0,0.0,-ires-1,iloc+2);
  newchannel->addIntermediate(rho0[2]  ,0,0.0, iloc  ,iloc+1);
  mode->addChannel(newchannel);
  return true;
}

// the particles produced by the current
tPDVector ThreePionCzyzCurrent::particles(int icharge, unsigned int,
					     int,int) {
  assert(icharge==0);
  tPDVector extpart(3);
  extpart[0]=getParticleData(ParticleID::piplus);
  extpart[1]=getParticleData(ParticleID::piminus);
  extpart[2]=getParticleData(ParticleID::pi0);
  // return the answer
  return extpart;
}


// hadronic current   
vector<LorentzPolarizationVectorE> 
ThreePionCzyzCurrent::current(const int imode, const int ichan,
			      Energy & scale,const ParticleVector & decay,
			      DecayIntegrator::MEOption meopt) const {
  useMe();
  if(meopt==DecayIntegrator::Terminate) {
    for(unsigned int ix=0;ix<3;++ix)
      ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
    return vector<LorentzPolarizationVectorE>(1,LorentzPolarizationVectorE());
  }
  // calculate q2,s1,s2,s3
  Lorentz5Momentum q;
  for(unsigned int ix=0;ix<decay.size();++ix)
    q+=decay[ix]->momentum();
  q.rescaleMass();
  scale=q.mass();
  Energy2 q2=q.mass2();
  Energy2 sm = (decay[1]->momentum()+decay[2]->momentum()).m2();
  Energy2 sp = (decay[0]->momentum()+decay[2]->momentum()).m2();
  Energy2 s0 = (decay[0]->momentum()+decay[1]->momentum()).m2();
  // isospin = 0
  complex<InvEnergy3> F_I1 =
    Resonance::H(rhoMasses_[0],rhoWidths_[0],sp,sm,s0,mpip_,mpi0_)*
    (coup_I0_[0]*Resonance::BreitWignerFW(q2,omegaMasses_[0],omegaWidths_[0])+
     coup_I0_[1]*Resonance::BreitWignerFW(q2,phiMass_       ,phiWidth_      )+
     coup_I0_[2]*Resonance::BreitWignerFW(q2,omegaMasses_[1],omegaWidths_[1])+
     coup_I0_[3]*Resonance::BreitWignerFW(q2,omegaMasses_[2],omegaWidths_[2]))
    +coup_I0_[4]*Resonance::H(rhoMasses_[1],rhoWidths_[1],sp,sm,s0,mpip_,mpi0_)*
    Resonance::BreitWignerFW(q2,phiMass_,phiWidth_)
    +coup_I0_[5]*Resonance::H(rhoMasses_[2],rhoWidths_[2],sp,sm,s0,mpip_,mpi0_)*
    Resonance::BreitWignerFW(q2,omegaMasses_[2],omegaWidths_[3]);
  // isospin = 1
  complex<InvEnergy3> F_I0 = GW_*
    Resonance::BreitWignerFW(q2,omegaMass_I1_,omegaWidth_I1_)/sqr(omegaMass_I1_)*
    (Resonance::BreitWignerPWave(s0,rhoMasses_I1_[0],
				 rhoWidths_I1_[0],mpip_,mpip_)/sqr(rhoMasses_I1_[0])+
     sigma_*Resonance::BreitWignerPWave(s0,rhoMasses_I1_[1],
					rhoWidths_I1_[1],mpip_,mpip_)/sqr(rhoMasses_I1_[0]));
  // the current
  LorentzPolarizationVector vect = (F_I0+F_I1)*
    Helicity::epsilon(decay[0]->momentum(),
		      decay[1]->momentum(),
		      decay[2]->momentum());
  // factor to get dimensions correct
  return vector<LorentzPolarizationVectorE>(1,q.mass()*vect);
}
   
bool ThreePionCzyzCurrent::accept(vector<int> id) {
  // check there are only two particles
  if(id.size()!=2) return false;
  // pion modes
  if((abs(id[0])==ParticleID::piplus  &&     id[1] ==ParticleID::pi0   ) ||
     (    id[0] ==ParticleID::pi0     && abs(id[1])==ParticleID::piplus))
    return true;
  else if((id[0]==ParticleID::piminus && id[1]==ParticleID::piplus) ||
	  (id[0]==ParticleID::piplus  && id[1]==ParticleID::piminus))
    return true;
  else
    return false;
}

// the decay mode
unsigned int ThreePionCzyzCurrent::decayMode(vector<int> idout) {
  return 0;
}

// output the information for the database
void ThreePionCzyzCurrent::dataBaseOutput(ofstream & output,bool header,
					bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::ThreePionCzyzCurrent " 
  		    << name() << " HwWeakCurrents.so\n";
  for(unsigned int ix=0;ix<rhoMasses_.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":RhoMassesI0 " << ix << " " << rhoMasses_[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<rhoWidths_.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":RhoWidthsI0 " << ix << " " << rhoWidths_[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<omegaMasses_.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":OmegaMassesI0 " << ix << " " << omegaMasses_[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<omegaWidths_.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":OmegaWidthsI0 " << ix << " " << omegaWidths_[ix]/GeV << "\n";
  }
  output << "newdef " << name() << ":PhiMass "  << phiMass_/GeV  << "\n";
  output << "newdef " << name() << ":PhiWidth " << phiWidth_/GeV << "\n";
  for(unsigned int ix=0;ix<coup_I0_.size();++ix) {
    if(ix<6) output << "newdef ";
    else     output << "insert ";
    output << name() << ":CouplingsI0 " << ix << " " << coup_I0_[ix]*GeV*GeV2 << "\n";
  }
  
  for(unsigned int ix=0;ix<rhoMasses_I1_.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":RhoMassesI1 " << ix << " " << rhoMasses_I1_[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<rhoWidths_I1_.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":RhoWidthsI1 " << ix << " " << rhoWidths_I1_[ix]/GeV << "\n";
  }
  output << "newdef " << name() << ":OmegaMass "  << omegaMass_I1_/GeV  << "\n";
  output << "newdef " << name() << ":OmegaWidth " << omegaWidth_I1_/GeV << "\n";
  output << "newdef " << name() << ":sigma "      << sigma_     << "\n";  
  output << "newdef " << name() << ":GWPrefactor "      << GW_pre_*GeV     << "\n";  
  output << "newdef " << name() << ":g_omega_pipi "      << g_omega_pi_pi_ << "\n";
  WeakDecayCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
  		    << fullName() << "\";" << endl;
}
