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
#include "ThePEG/Helicity/ScalarSpinInfo.h"
#include "ThePEG/Helicity/VectorSpinInfo.h"

namespace Herwig {
using namespace ThePEG;
using ThePEG::Helicity::ScalarSpinInfo;
using ThePEG::Helicity::VectorSpinInfo;
using ThePEG::Helicity::VectorSpinPtr;
using ThePEG::Helicity::tcVectorSpinPtr;
using Helicity::VectorWaveFunction;
using Helicity::Direction;
using Helicity::incoming;
using Helicity::outgoing;

TwoPionPhotonCurrent::~TwoPionPhotonCurrent() {}

void TwoPionPhotonCurrent::persistentOutput(PersistentOStream & os) const {
  os << _grho << _grhoomegapi << _resweights << _rhoparameters << _rhomasses 
     << _rhowidths << _omegaparameters << _omegamass << _omegawidth ;
}

void TwoPionPhotonCurrent::persistentInput(PersistentIStream & is, int) { 
  is >> _grho >> _grhoomegapi >> _resweights >> _rhoparameters >> _rhomasses 
     >> _rhowidths >> _omegaparameters >> _omegamass >> _omegawidth ;
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

  static ParVector<TwoPionPhotonCurrent,double> interfaceRhoMasses
    ("RhoMasses",
     "The masses of the different rho resonances for the decay tau ->  pi pi photon",
     &TwoPionPhotonCurrent::_rhomasses,
     0, 0, 0, -10000, 10000, false, false, true);
  
  static ParVector<TwoPionPhotonCurrent,double> interfaceRhoWidths
    ("RhoWidths",
     "The widths of the different rho resonances for the decay tau -> nu pi pi photon",
     &TwoPionPhotonCurrent::_rhowidths,
     0, 0, 0, -10000, 10000, false, false, true);

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
     &TwoPionPhotonCurrent::_grho, GeV2, 0.11238947*GeV2, -1.0e12*GeV2, 1.0e12*GeV2,
     false, false, false);

  static Parameter<TwoPionPhotonCurrent,InvEnergy> interfacegrhoomegapi
    ("grhoomegapi",
     "The rho-omega-pi coupling",
     &TwoPionPhotonCurrent::_grhoomegapi, 1./GeV, 12.924/GeV,
     -1.0e12*1./GeV, 1.0e12*1./GeV,
     false, false, false);
}


// complete the construction of the decay mode for integration
bool TwoPionPhotonCurrent::createMode(int icharge, unsigned int imode,
				      DecayPhaseSpaceModePtr mode,
				      unsigned int iloc,unsigned int ires,
				      DecayPhaseSpaceChannelPtr phase,Energy upp)
{
  if(icharge!=3&&icharge!=-3){return false;}
  bool kineallowed=true;
  // check that the mode is are kinematical allowed
  Energy min=getParticleData(ParticleID::piplus)->mass()
    +getParticleData(ParticleID::pi0)->mass();
  if(min>upp){kineallowed=false;}
  if(kineallowed==false){return kineallowed;}
  // set up the integration channels;
  tPDPtr omega=getParticleData(ParticleID::omega);
  tPDPtr temp;
  DecayPhaseSpaceChannelPtr newchannel; 
  for(unsigned int ix=0;ix<3;++ix)
    {
      if(ix==0){temp = getParticleData(-213);}
      else if(ix==1){temp = getParticleData(-100213);}
      else if(ix==2){temp = getParticleData(-30213);}
      if(temp)
	{
	  newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(temp,0,0.0,-ires-1,iloc);
	  newchannel->addIntermediate(omega,0,0.0,iloc+1,iloc+2);
	  newchannel->init();
	  mode->addChannel(newchannel);
	}
    }
  // reset the masses and widths of the resonances if needed
  // set up the rho masses and widths
  for(unsigned int ix=0;ix<3;++ix)
    {
      if(ix==0){temp = getParticleData(-213);}
      else if(ix==1){temp = getParticleData(-100213);}
      else if(ix==2){temp = getParticleData(-30213);}
      // if using local values
      if(_rhoparameters&&ix<_rhomasses.size())
	{mode->resetIntermediate(temp,_rhomasses[ix],_rhowidths[ix]);}
    }
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
			      const Particle & inpart,
			      const ParticleVector & decay) const
{
  vector<LorentzPolarizationVector> temp;
  LorentzPolarizationVector vect;
  // locate the particles
  unsigned int ipi0=decay.size()-2,ipic=decay.size()-3,igamma=decay.size()-1;
  Lorentz5Momentum pout=decay[ipi0]->momentum()+decay[igamma]->momentum()
    +decay[ipic]->momentum();
  // overall hadronic mass
  Energy2 q2 = pout.m2();
  // mass of the omega
  pout = decay[igamma]->momentum()+decay[ipi0]->momentum();
  Energy s2  = pout.m2();
  // compute the prefactor
  complex<double> prefactor=-FFunction(0.,-1)*FFunction(q2,ichan)*inpart.mass()*
    sqrt(2.*pi*generator()->standardModel()->alphaEM())*BreitWigner(s2,10);
  // dot products which don't depend on the polarization vector
  Energy2 dot12 = decay[igamma]->momentum().dot(decay[ipi0]->momentum());
  Energy2 dot13 = decay[igamma]->momentum().dot(decay[ipic]->momentum());
  Energy2 dot23 = decay[ipi0]->momentum().dot(decay[ipic]->momentum());
  Energy2 mpi2  = decay[ipic]->mass()*decay[ipic]->mass();
  // construct the spininfomation objects for the decay products
  VectorSpinPtr hwtemp;
  if(vertex)
    {
      SpinPtr spi0=new_ptr(ScalarSpinInfo(decay[ipi0]->momentum(),true));
      decay[ipi0]->spinInfo(spi0);
      SpinPtr spic=new_ptr(ScalarSpinInfo(decay[ipic]->momentum(),true));
      decay[ipic]->spinInfo(spic);
      SpinPtr sgamma=new_ptr(VectorSpinInfo(decay[igamma]->momentum(),true));
      decay[igamma]->spinInfo(sgamma);
      hwtemp= dynamic_ptr_cast<VectorSpinPtr>(sgamma);
    }
  // loop over the photon polarizations
  VectorWaveFunction photon(decay[igamma]->momentum(),decay[igamma]->dataPtr(),
			    outgoing);
  complex<double> dote2,dote3,coeffa,coeffb,coeffc;
  for(int ihel=-1;ihel<2;++ihel)
    {
      if(ihel!=0)
	{
	  photon.reset(ihel);
	  if(vertex){hwtemp->setBasisState(ihel,photon.Wave());}
	  // obtain the dot products we need
	  dote2 =
	    +photon.t()*decay[ipi0]->momentum().e()
	    -photon.x()*decay[ipi0]->momentum().px()
	    -photon.y()*decay[ipi0]->momentum().py()
	    -photon.z()*decay[ipi0]->momentum().pz();
	  dote3 =
	    +photon.t()*decay[ipic]->momentum().e()
	    -photon.x()*decay[ipic]->momentum().px()
	    -photon.y()*decay[ipic]->momentum().py()
	    -photon.z()*decay[ipic]->momentum().pz();
	  // now compute the coefficients
	  coeffa = mpi2*dot13-dot12*(dot23-dot13);
	  coeffb = dote2*dot13-dote3*dot12;
	  coeffc = dote2*dot23-dote3*(mpi2+dot12);
	  // finally compute the current
	  for(unsigned int ix=0;ix<4;++ix)
	    {vect[ix]=
		 photon(ix)*coeffa
	       -decay[ipi0]->momentum()[ix]*coeffb
	       +decay[igamma]->momentum()[ix]*coeffc;
	      vect[ix]*=prefactor;
	    }
	  temp.push_back(LorentzPolarizationVector(vect[0],vect[1],
						   vect[2],vect[3]));
	}
      else
	{
	  LorentzPolarizationVector veczero=LorentzPolarizationVector();
	  temp.push_back(veczero);
	  if(vertex){hwtemp->setBasisState(ihel,photon.Wave());}
	}
    }
  return temp;
}


bool TwoPionPhotonCurrent::accept(vector<int> id)
{
  if(id.size()!=3){return false;}
  unsigned int npiplus(0),npi0(0),ngamma(0);
  for(unsigned int ix=0;ix<id.size();++ix)
    {
      if(abs(id[ix])==ParticleID:: piplus){++npiplus;}
      else if(id[ix]==ParticleID::gamma){++ngamma;}
      else if(id[ix]==ParticleID::pi0){++npi0;}
    }
  return npiplus==1&&ngamma==1&&npi0==1;
}

unsigned int TwoPionPhotonCurrent::decayMode(vector<int> idout){return 0;}

// output the information for the database
void TwoPionPhotonCurrent::dataBaseOutput(ofstream & output)
{
  output << "create /Herwig++/TwoPionPhotonCurrent " << fullName() << " \n";
  output << "set " << fullName() << ":RhoParameters "    << _rhoparameters << "\n";
  output << "set " << fullName() << ":omegaParameters "    << _omegaparameters << "\n";
  output << "set " << fullName() << ":omegamass "    << _omegamass/GeV << "\n";
  output << "set " << fullName() << ":omegawidth "    << _omegawidth/GeV << "\n";
  output << "set " << fullName() << ":grho "    << _grho/GeV2 << "\n";
  output << "set " << fullName() << ":grhoomegapi "    << _grhoomegapi*GeV << "\n";
  for(unsigned int ix=0;ix<_resweights.size();++ix)
    {output << "insert " << fullName() << ":Weights " << ix 
	    << " " << _resweights[ix] << "\n";}
  for(unsigned int ix=0;ix<_rhomasses.size();++ix)
    {output << "insert " << fullName() << ":RhoMasses " << ix 
	    << " " << _rhomasses[ix] << "\n";}
  for(unsigned int ix=0;ix<_rhowidths.size();++ix)
    {output << "insert " << fullName() << ":RhoWidths " << ix 
	    << " " << _rhowidths[ix] << "\n";}
}
}
