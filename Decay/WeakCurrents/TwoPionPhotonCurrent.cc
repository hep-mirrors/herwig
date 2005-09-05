// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TwoPionPhotonCurrent class.
//
//  Author: Peter Richardson
//

#include "TwoPionPhotonCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "TwoPionPhotonCurrent.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"

namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;
using Helicity::ScalarWaveFunction;
using Helicity::VectorWaveFunction;
using Helicity::outgoing;

TwoPionPhotonCurrent::~TwoPionPhotonCurrent() {}

void TwoPionPhotonCurrent::persistentOutput(PersistentOStream & os) const {
  os << _grho << _grhoomegapi << _resweights << _rhoparameters << _rhomasses 
     << _rhowidths << _omegaparameters << _omegamass << _omegawidth << _intmass 
     << _intwidth ;
}

void TwoPionPhotonCurrent::persistentInput(PersistentIStream & is, int) { 
  is >> _grho >> _grhoomegapi >> _resweights >> _rhoparameters >> _rhomasses 
     >> _rhowidths >> _omegaparameters >> _omegamass >> _omegawidth >> _intmass
     >> _intwidth;
}

ClassDescription<TwoPionPhotonCurrent> TwoPionPhotonCurrent::initTwoPionPhotonCurrent;
// Definition of the static class description member.

void TwoPionPhotonCurrent::Init() {

  static ParVector<TwoPionPhotonCurrent,double> interfacereswgt
    ("Weights",
     "The weights of the different resonances for the decay tau -> nu pi pi gamma",
     &TwoPionPhotonCurrent::_resweights,
     0, 0, 0, -1000, 1000, false, false, true);
    
  static Switch<TwoPionPhotonCurrent,bool> interfaceRhoParameters
    ("RhoParameters",
     "Use local values for the rho meson masses and widths",
     &TwoPionPhotonCurrent::_rhoparameters, true, false, false);

  static SwitchOption interfaceRhoParameterstrue
    (interfaceRhoParameters,
     "Local",
     "Use local values",
     true);
  static SwitchOption interfaceRhoParametersParticleData
    (interfaceRhoParameters,
     "ParticleData",
     "Use the value from the particle data objects",
     false);
  
  static Switch<TwoPionPhotonCurrent,bool> interfaceomegaParameters
    ("omegaParameters",
     "Use local values for the omega meson masses and widths",
     &TwoPionPhotonCurrent::_omegaparameters, true, false, false);
  static SwitchOption interfaceomegaParameterstrue
    (interfaceomegaParameters,
     "Local",
     "Use local values",
     true);
  static SwitchOption interfaceomegaParametersParticleData
    (interfaceomegaParameters,
     "ParticleData",
     "Use the value from the particle data objects",
     false);

  static ParVector<TwoPionPhotonCurrent,Energy> interfaceRhoMasses
    ("RhoMasses",
     "The masses of the different rho resonances for the decay tau ->  pi pi photon",
     &TwoPionPhotonCurrent::_rhomasses, MeV, -1, 773.*MeV, 0.0*MeV, 10000.*MeV,
     false, false, true);

  static ParVector<TwoPionPhotonCurrent,Energy> interfaceRhoWidths
    ("RhoWidths",
     "The widths of the different rho resonances for the decay tau -> nu pi pi photon",
     &TwoPionPhotonCurrent::_rhowidths, MeV, -1, 145.*MeV, 0.0*MeV, 1000.*MeV,
     false, false, true);

  static Parameter<TwoPionPhotonCurrent,Energy> interfaceomegamass
    ("omegamass",
     "The mass of the omega",
     &TwoPionPhotonCurrent::_omegamass, GeV, 0.782*GeV, 0.0*GeV, 1.0*GeV,
     false, false, true);
  
  static Parameter<TwoPionPhotonCurrent,Energy> interfaceomegawidth
    ("omegawidth",
     "The width of the omega for the decay tau- -> pi pi photon",
     &TwoPionPhotonCurrent::_omegawidth, GeV, 0.0085*GeV, 0.*GeV, 1.*GeV,
     false, false, false);
  
  static ClassDocumentation<TwoPionPhotonCurrent> documentation
    ("The \\classname{TwoPionPhotonCurrent} class implements the decay "
     "tau+/- -> pi+/- pi0 gamma via an omega.");
  
  static Parameter<TwoPionPhotonCurrent,Energy2> interfacegrho
    ("grho",
     "The rho meson decay constant.",
     &TwoPionPhotonCurrent::_grho, GeV2, 0.11238947*GeV2, -1.*GeV2, 1.*GeV2,
     false, false, false);

  static Parameter<TwoPionPhotonCurrent,InvEnergy> interfacegrhoomegapi
    ("grhoomegapi",
     "The rho-omega-pi coupling",
     &TwoPionPhotonCurrent::_grhoomegapi, 1./GeV, 12.924/GeV,
     -100./GeV, 100./GeV,
     false, false, false);

  static Parameter<TwoPionPhotonCurrent,Energy> interfaceIntegrationMass
    ("IntegrationMass",
     "Mass of the pseudoresonance used to improve integration effciency",
     &TwoPionPhotonCurrent::_intmass, GeV, 1.4*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<TwoPionPhotonCurrent,Energy> interfaceIntegrationWidth
    ("IntegrationWidth",
     "Width of the pseudoresonance used to improve integration effciency",
     &TwoPionPhotonCurrent::_intwidth, GeV, 0.5*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);
}

// complete the construction of the decay mode for integration
bool TwoPionPhotonCurrent::createMode(int icharge, unsigned int imode,
				      DecayPhaseSpaceModePtr mode,
				      unsigned int iloc,unsigned int ires,
				      DecayPhaseSpaceChannelPtr phase,Energy upp)
{
  if(icharge!=3&&icharge!=-3){return false;}
  bool kineallowed(true);
  // check that the mode is are kinematical allowed
  Energy min(getParticleData(ParticleID::piplus)->mass()+
	     getParticleData(ParticleID::pi0)->mass());
  if(min>upp){kineallowed=false;}
  if(kineallowed==false){return kineallowed;}
  // set up the integration channels;
  tPDPtr omega(getParticleData(ParticleID::omega));
  tPDPtr W(getParticleData(ParticleID::Wplus));
  if(icharge<0){W=W->CC();}
  DecayPhaseSpaceChannelPtr newchannel; 
  newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
  newchannel->addIntermediate(W,0,0.0,-ires-1,iloc);
  newchannel->addIntermediate(omega,0,0.0,iloc+1,iloc+2);
  mode->addChannel(newchannel);
  // reset the masses and widths of the resonances if needed
  mode->resetIntermediate(W,_intmass,_intwidth);
  // set up the omega masses and widths
  if(_omegaparameters){mode->resetIntermediate(omega,_omegamass,_omegawidth);}
  return kineallowed;
}

// the particles produced by the current
PDVector TwoPionPhotonCurrent::particles(int icharge, unsigned int imode,int iq,int ia)
{
  PDVector extpart;
  if(icharge==3)
    {
      extpart.push_back(getParticleData(ParticleID::piplus));
      extpart.push_back(getParticleData(ParticleID::pi0));
      extpart.push_back(getParticleData(ParticleID::gamma));
    }
  else if(icharge==-3)
    {
      extpart.push_back(getParticleData(ParticleID::piminus));
      extpart.push_back(getParticleData(ParticleID::pi0));
      extpart.push_back(getParticleData(ParticleID::gamma));
    }
  return extpart;
}


// the hadronic currents    
vector<LorentzPolarizationVector> 
TwoPionPhotonCurrent::current(bool vertex, const int imode, const int ichan, 
			      Energy & scale,const ParticleVector & decay) const
{
  vector<LorentzPolarizationVector> temp;
  // locate the particles
  Lorentz5Momentum pout(decay[1]->momentum()+decay[2]->momentum()+decay[0]->momentum());
  // overall hadronic mass
  pout.rescaleMass();
  scale=pout.mass();
  Energy2 q2(pout.m2());
  // mass of the omega
  pout = decay[1]->momentum()+decay[2]->momentum();
  pout.rescaleMass();
  Energy s2(pout.m2());
  // compute the prefactor
  Complex prefactor(-FFunction(0.)*FFunction(q2)*scale*
		    sqrt(2.*pi*generator()->standardModel()->alphaEM())*
		    BreitWigner(s2,10));
  // dot products which don't depend on the polarization vector
  Energy2 dot12(decay[2]->momentum()*decay[1]->momentum());
  Energy2 dot13(decay[2]->momentum()*decay[0]->momentum());
  Energy2 dot23(decay[1]->momentum()*decay[0]->momentum());
  Energy2 mpi2 (decay[0]->mass()*decay[0]->mass());

  // workaround for gcc 3.2.3 bug
  // construct the spininfomation objects for the decay products
  //ALB ScalarWaveFunction(decay[0],outgoing,true,vertex);
  //ALB ScalarWaveFunction(decay[1],outgoing,true,vertex);
  PPtr mytemp = decay[0];
  ScalarWaveFunction(mytemp,outgoing,true,vertex);
  PPtr mytemp_bis = decay[1];
  ScalarWaveFunction(mytemp_bis,outgoing,true,vertex);

  VectorWaveFunction(temp,decay[2],outgoing,true,true,vertex);
  Complex dote2,dote3,coeffa,coeffb,coeffc;
  for(unsigned int ix=0;ix<3;++ix)
    {
      if(ix!=1)
	{
	  // obtain the dot products we need
	  dote2 = temp[ix]*decay[1]->momentum();
	  dote3 = temp[ix]*decay[0]->momentum();
	  // now compute the coefficients
	  coeffa = mpi2*dot13-dot12*(dot23-dot13);
	  coeffb = dote2*dot13-dote3*dot12;
	  coeffc = dote2*dot23-dote3*(mpi2+dot12);
	  // finally compute the current
	  temp[ix]= prefactor*(coeffa*temp[ix]-coeffb*decay[1]->momentum()+
			       coeffc*decay[2]->momentum());
	}
      else{temp[ix]=LorentzPolarizationVector();}
    }
  return temp;
}

bool TwoPionPhotonCurrent::accept(vector<int> id)
{
  if(id.size()!=3){return false;}
  unsigned int npiplus(0),npi0(0),ngamma(0);
  for(unsigned int ix=0;ix<id.size();++ix)
    {
      if(abs(id[ix])==ParticleID::piplus){++npiplus;}
      else if(id[ix]==ParticleID::gamma){++ngamma;}
      else if(id[ix]==ParticleID::pi0){++npi0;}
    }
  return npiplus==1&&ngamma==1&&npi0==1;
}

unsigned int TwoPionPhotonCurrent::decayMode(vector<int> idout){return 0;}

// output the information for the database
void TwoPionPhotonCurrent::dataBaseOutput(ofstream & output) const
{
  output << "create /Herwig++/TwoPionPhotonCurrent " << fullName() << " \n";
  output << "set " << fullName() << ":RhoParameters "    << _rhoparameters << "\n";
  output << "set " << fullName() << ":omegaParameters "    << _omegaparameters << "\n";
  output << "set " << fullName() << ":omegamass "    << _omegamass/GeV << "\n";
  output << "set " << fullName() << ":omegawidth "    << _omegawidth/GeV << "\n";
  output << "set " << fullName() << ":grho "    << _grho/GeV2 << "\n";
  output << "set " << fullName() << ":grhoomegapi "    << _grhoomegapi*GeV << "\n";
  output << "set " << fullName() << ":IntegrationMass "  << _intmass/GeV  << "\n";
  output << "set " << fullName() << ":IntegrationWidth " << _intwidth/GeV  << "\n";
  unsigned int ix;
  for(ix=0;ix<_resweights.size();++ix)
    {
      if(ix<3){output << "set " << fullName() << ":Weights " << ix 
		      << " " << _resweights[ix] << "\n";}
      else{output << "insert " << fullName() << ":Weights " << ix 
		  << " " << _resweights[ix] << "\n";}
    }
  for(ix=0;ix<_rhomasses.size();++ix)
    {
      if(ix<2){output << "set " << fullName() << ":RhoMasses " << ix 
		      << " " << _rhomasses[ix]/MeV << "\n";}
      else{output << "insert " << fullName() << ":RhoMasses " << ix 
		  << " " << _rhomasses[ix]/MeV << "\n";}
    }
  for(ix=0;ix<_rhowidths.size();++ix)
    {
      if(ix<2){output << "set " << fullName() << ":RhoWidths " << ix 
		      << " " << _rhowidths[ix]/MeV << "\n";}
      else{output << "insert " << fullName() << ":RhoWidths " << ix 
		  << " " << _rhowidths[ix]/MeV << "\n";}
    }
}
}
