// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TwoMesonRhoKStarCurrent class.
//
//  Author: Peter Richardson
//

#include "TwoMesonRhoKStarCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/ScalarSpinInfo.h"


#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "TwoMesonRhoKStarCurrent.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"


namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;
using ThePEG::Helicity::ScalarSpinInfo;

TwoMesonRhoKStarCurrent::~TwoMesonRhoKStarCurrent() {}

void TwoMesonRhoKStarCurrent::persistentOutput(PersistentOStream & os) const {
  os << _pimodel << _Kmodel << _piwgt << _kwgt << _rhoparameters
     << _Kstarparameters << _rhomasses << _rhowidths << _Kstarmasses << _Kstarwidths
     <<  _mass << _width << _mass2 << _massw << _massa <<_massb << _mom << _dhdq2
     << _hm2 << _dparam;
}

void TwoMesonRhoKStarCurrent::persistentInput(PersistentIStream & is, int) {
  is >> _pimodel >> _Kmodel >> _piwgt >> _kwgt >> _rhoparameters
     >> _Kstarparameters >> _rhomasses >> _rhowidths >> _Kstarmasses >> _Kstarwidths
     >> _mass >> _width >> _mass2 >> _massw >> _massa >> _massb >> _mom >> _dhdq2
     >> _hm2 >> _dparam;
}

ClassDescription<TwoMesonRhoKStarCurrent> TwoMesonRhoKStarCurrent::initTwoMesonRhoKStarCurrent;
// Definition of the static class description member.

void TwoMesonRhoKStarCurrent::Init() {

  static ParVector<TwoMesonRhoKStarCurrent,double> interfaceRhoMasses
    ("RhoMasses",
     "The masses of the different rho resonances for the pi pi channel",
     &TwoMesonRhoKStarCurrent::_rhomasses,
     0, 0, 0, -10000, 10000, false, false, true);
  
  static ParVector<TwoMesonRhoKStarCurrent,double> interfaceRhoWidths
    ("RhoWidths",
     "The widths of the different rho resonances for the pi pi channel",
     &TwoMesonRhoKStarCurrent::_rhowidths,
     0, 0, 0, -10000, 10000, false, false, true);
  
  static ParVector<TwoMesonRhoKStarCurrent,double> interfaceKstarMasses
    ("KstarMasses",
     "The masses of the different rho resonances for the pi pi channel",
     &TwoMesonRhoKStarCurrent::_Kstarmasses,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<TwoMesonRhoKStarCurrent,double> interfaceKstarWidths
    ("KstarWidths",
     "The widths of the different rho resonances for the pi pi channel",
     &TwoMesonRhoKStarCurrent::_Kstarwidths,
     0, 0, 0, -10000, 10000, false, false, true);
  
  static Switch<TwoMesonRhoKStarCurrent,bool> interfaceRhoParameters
    ("RhoParameters",
     "Use local values for the rho meson masses and widths",
     &TwoMesonRhoKStarCurrent::_rhoparameters, true, false, false);
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

  static Switch<TwoMesonRhoKStarCurrent,bool> interfaceKstarParameters
    ("KstarParameters",
     "Use local values for the Kstar meson masses and widths",
     &TwoMesonRhoKStarCurrent::_Kstarparameters, true, false, false);
  static SwitchOption interfaceKstarParameterstrue
    (interfaceKstarParameters,
     "Local",
     "Use local values",
     true);
  static SwitchOption interfaceKstarParametersParticleData
    (interfaceKstarParameters,
     "ParticleData",
     "Use the value from the particle data objects",
     false);

  static ParVector<TwoMesonRhoKStarCurrent,double> interfacepiwgt
    ("PiWeight",
     "The weights of the different resonances for the pi pi channel",
     &TwoMesonRhoKStarCurrent::_piwgt,
     0, 0, 0, -1000, 1000, false, false, true);
  
  static ParVector<TwoMesonRhoKStarCurrent,double> interfacekwgt
    ("KWeight",
     "The weights of the different resonances for the K pi channel",
     &TwoMesonRhoKStarCurrent::_kwgt,
     0, 0, 0, -1000, 1000, false, false, true);
  
  static Switch<TwoMesonRhoKStarCurrent,int> interfacePiModel
    ("PiModel",
     "The model to use for the propagator for the pion modes.",
     &TwoMesonRhoKStarCurrent::_pimodel, 0, false, false);
  static SwitchOption interfacePiModelKuhn
    (interfacePiModel,
     "Kuhn",
     "The model of Kuhn and Santamaria",
     0);
  static SwitchOption interfacePiModelGounaris
    (interfacePiModel,
     "Gounaris",
     "The model of Gounaris and Sakurai.",
     1);
  
  static Switch<TwoMesonRhoKStarCurrent,int> interfaceKModel
    ("KModel",
     "The model to use for the propagator for the kaon modes.",
     &TwoMesonRhoKStarCurrent::_Kmodel, 0, false, false);
  static SwitchOption interfaceKModelKuhn
    (interfaceKModel,
     "Kuhn",
     "The model of Kuhn and Santamaria",
     0);
  static SwitchOption interfaceKModelGounaris
    (interfaceKModel,
     "Gounaris",
     "The model of Gounaris and Sakurai.",
     1);

  static ClassDocumentation<TwoMesonRhoKStarCurrent> documentation
    ("The \\classname{TwoMesonRhoKStarCurrent} class is designed to implement weak"
     "decay to two scalar mesons using the models of either Kuhn and "
     "Santamaria (Z. Phys. C48, 445 (1990)) or Gounaris and Sakurai Phys. Rev. "
     "Lett. 21, 244 (1968).  The mixing parameters are taken from "
     "Phys. Rev. D61:112002,2000 (CLEO), although the PDG values for the "
     "masses and widths are used, for the decay pi+/- pi0."
     " The decay K pi is assumed to  be dominated by the lowest lying K* resonance.");

}

// complete the construction of the decay mode for integration
bool TwoMesonRhoKStarCurrent::createMode(int icharge, unsigned int imode,
					 DecayPhaseSpaceModePtr mode,
					 unsigned int iloc,unsigned int ires,
					 DecayPhaseSpaceChannelPtr phase,Energy upp)
{
  if(abs(icharge)!=3){return false;}
  // make sure that the decays are kinematically allowed
  bool kineallowed(true);
  tPDPtr part[2];
  if(imode==0)
    {
      part[0]=getParticleData(ParticleID::piplus);
      part[1]=getParticleData(ParticleID::pi0);
    }
  else if(imode==1)
    {
      part[0]=getParticleData(ParticleID::Kplus);
      part[1]=getParticleData(ParticleID::pi0);
    }
  else if(imode==2)
    {
      part[0]=getParticleData(ParticleID::K0);
      part[1]=getParticleData(ParticleID::piplus);
    }
  else if(imode==3)
    {
      part[0]=getParticleData(ParticleID::Kplus);
      part[1]=getParticleData(ParticleID::K0);
    }
  else if(imode==4)
    {
      part[0]=getParticleData(ParticleID::eta);
      part[1]=getParticleData(ParticleID::Kplus);
    }
  else
    {kineallowed=false;}
  Energy min=
    +part[0]->mass()-part[0]->widthLoCut()
    +part[1]->mass()-part[1]->widthLoCut();
  if(min>upp){kineallowed=false;}
  if(kineallowed==false){return kineallowed;}
  DecayPhaseSpaceChannelPtr newchannel;
  // set the resonances
  // two pion or  K+ K0 decay
  tPDPtr res[3];
  if(imode==0||imode==3)
    {
      if(icharge==3)
	{
	  res[0]=getParticleData(213);
	  res[1]=getParticleData(100213);
	  res[2]=getParticleData(30213);
	}
      else
	{
	  res[0]=getParticleData(-213);
	  res[1]=getParticleData(-100213);
	  res[2]=getParticleData(-30213);
	}
    }
  // K+ pi0 or K0 pi+ or K eta decay
  else if(imode==1||imode==2||imode==4)
    {
      if(icharge==3)
	{
	  res[0]=getParticleData(323);
	  res[1]=getParticleData(100323);
	  res[2]=getParticleData(30323);
	}
      else
	{
	  res[0]=getParticleData(-323);
	  res[1]=getParticleData(-100323);
	  res[2]=getParticleData(-30323);
	}
    }
  else
    {cerr << "Failure of initialisation in TwoMesonRhoKStarCurrent" << endl;}
  // create the channels
  for(unsigned int ix=0;ix<3;++ix)
    {
      if(res[ix])
	{
	  newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(res[ix],0,0.0,iloc,iloc+1);
	  newchannel->init();
	  mode->addChannel(newchannel);
	}
    }
  // reset the masses in the intergrators if needed
  // for the rho 
  if(_rhoparameters&&(imode==0||imode==3))
    {for(unsigned int ix=0;ix<3;++ix)
	{if(ix<_rhomasses.size())
	    {mode->resetIntermediate(res[ix],_rhomasses[ix],_rhowidths[ix]);}}}
  // for the K*
  else if(_Kstarparameters)
    {for(unsigned int ix=0;ix<3;++ix)
	{if(ix<_Kstarmasses.size())
	    {mode->resetIntermediate(res[ix],_Kstarmasses[ix],_Kstarwidths[ix]);}}}
  // return if successful
  return kineallowed;
}

// the particles produced by the current
  PDVector TwoMesonRhoKStarCurrent::particles(int icharge, unsigned int imode,
					      int iq,int ia)
{
  PDVector output(2);
  if(icharge==3)
    {
      if(imode==0)
	{
	  output[0]=getParticleData(ParticleID::piplus);
	  output[1]=getParticleData(ParticleID::pi0);
	}
      else if(imode==1)
	{
	  output[0]=getParticleData(ParticleID::Kplus);
	  output[1]=getParticleData(ParticleID::pi0);
	}
      else if(imode==2)
	{
	  output[0]=getParticleData(ParticleID::K0);
	  output[1]=getParticleData(ParticleID::piplus);
	}
      else if(imode==3)
	{
	  output[0]=getParticleData(ParticleID::Kplus);
	  output[1]=getParticleData(ParticleID::Kbar0);
	}
      else
	{
	  output[0]=getParticleData(ParticleID::eta);
	  output[1]=getParticleData(ParticleID::Kplus);
	}
    }
  else if(icharge==-3)
    {
      if(imode==0)
	{
	  output[0]=getParticleData(ParticleID::piminus);
	  output[1]=getParticleData(ParticleID::pi0);
	}
      else if(imode==1)
	{
	  output[0]=getParticleData(ParticleID::Kminus);
	  output[1]=getParticleData(ParticleID::pi0);
	}
      else if(imode==2)
	{
	  output[0]=getParticleData(ParticleID::Kbar0);
	  output[1]=getParticleData(ParticleID::piminus);
	}
      else if(imode==3)
	{
	  output[0]=getParticleData(ParticleID::Kminus);
	  output[1]=getParticleData(ParticleID::K0);
	}
      else
	{
	  output[0]=getParticleData(ParticleID::eta);
	  output[1]=getParticleData(ParticleID::Kminus);
	}
    }
  return output;
}

// hadronic current   
vector<LorentzPolarizationVector> 
TwoMesonRhoKStarCurrent::current(bool vertex, const int imode, const int ichan, 
				 const Particle & inpart,
				 const ParticleVector & outpart) const
{
  // find the mesons
  unsigned int imes[2]={outpart.size()-2,outpart.size()-1};
  // momentum difference and sum of the mesons
  Lorentz5Momentum pdiff=outpart[imes[0]]->momentum()-outpart[imes[1]]->momentum();
  Lorentz5Momentum psum =outpart[imes[0]]->momentum()+outpart[imes[1]]->momentum();
  // mass2 of vector intermediate state
  Energy2 q2=psum.m2();
  double dot = psum.dot(pdiff)/q2;
  psum *=dot;
  LorentzPolarizationVector vect;
  // calculate the current
  Complex FPI(0.),denom(0.);
  // pion mode
  if(imode==0)
    {
      if(ichan<0)
	{
	  for(unsigned int ix=0;ix<_piwgt.size()&&ix<3;++ix)
	    {
	      FPI+=_piwgt[ix]*BreitWigner(q2,_pimodel,0,ix);
	      denom+=_piwgt[ix];
	    }
	}
      else if(ichan<int(_piwgt.size())&&ichan<3)
	{FPI=_piwgt[ichan]*BreitWigner(q2,_pimodel,0,ichan);denom=1.;}
      FPI *= sqrt(2.0)/denom;
    }
  // single kaon modes
  else if(imode==1||imode==2)
    {
      if(ichan<0)
	{
	  for(unsigned int ix=0;ix<_kwgt.size()&&ix<3;++ix)
	    {
	      FPI+=_kwgt[ix]*BreitWigner(q2,_Kmodel,1,ix);
	      denom+=_kwgt[ix];
	    }
	}
      else if (ichan<int(_kwgt.size())&&ichan<3)
	{FPI=_kwgt[ichan]*BreitWigner(q2,_Kmodel,1,ichan);denom=1.;}
      if(imode==1){FPI *= sqrt(2./3.)/denom;}
      else{FPI *= sqrt(4./3.)/denom;}
      // compute the current
      pdiff-=psum;
    }
  // two kaon modes
  else if(imode==3)
    {
      if(ichan<0)
	{
	  for(unsigned int ix=0;ix<_piwgt.size()&&ix<3;++ix)
	    {
	      FPI+=_piwgt[ix]*BreitWigner(q2,_pimodel,0,ix);
	      denom+=_piwgt[ix];
	    }
	}
      else if(ichan<int(_piwgt.size())&&ichan<3)
	{FPI=_piwgt[ichan]*BreitWigner(q2,_pimodel,0,ichan);denom=1.;}
      FPI *= sqrt(2.0)/denom;
      pdiff-=psum;
    }
  // the kaon eta mode
  else if(imode==4)
    {
      if(ichan<0)
	{
	  for(unsigned int ix=0;ix<_kwgt.size()&&ix<3;++ix)
	    {
	      FPI+=_kwgt[ix]*BreitWigner(q2,_Kmodel,1,ix);
	      denom+=_kwgt[ix];
	    }
	}
      else if (ichan<int(_kwgt.size())&&ichan<3)
	{FPI=_kwgt[ichan]*BreitWigner(q2,_Kmodel,1,ichan);denom=1.;
      FPI *=2.*denom;}
      pdiff-=psum;
    }
  // storage for the current
  vector<LorentzPolarizationVector> temp;temp.push_back(FPI*pdiff);
  // set up the spininfo for the decay products
  if(vertex)
    {
      for(unsigned int ix=0;ix<2;++ix)
	{
	  SpinPtr stemp = new_ptr(ScalarSpinInfo(outpart[imes[ix]]->momentum(),true));
	  outpart[imes[ix]]->spinInfo(stemp);
	}
    }
  return temp;
}
   
bool TwoMesonRhoKStarCurrent::accept(vector<int> id)
{
  bool allowed(false);
  // check there are only two particles
  if(id.size()!=2){return false;}
  // pion modes
  if((id[0]==ParticleID::piminus && id[1]==ParticleID::pi0)     ||
     (id[0]==ParticleID::pi0     && id[1]==ParticleID::piminus) ||
     (id[0]==ParticleID::piplus  && id[1]==ParticleID::pi0)     ||
     (id[0]==ParticleID::pi0     && id[1]==ParticleID::piplus))
    {allowed=true;}
  // single charged kaon
  else if((id[0]==ParticleID::Kminus && id[1]==ParticleID::pi0)    ||
	  (id[0]==ParticleID::pi0    && id[1]==ParticleID::Kminus) ||
	  (id[0]==ParticleID::Kplus  && id[1]==ParticleID::pi0)    ||
	  (id[0]==ParticleID::pi0    && id[1]==ParticleID::Kplus))
    {allowed=true;}
  // single neutral kaon
  else if((id[0]==ParticleID::piminus && id[1]==ParticleID::Kbar0)   ||
	  (id[0]==ParticleID::Kbar0   && id[1]==ParticleID::piminus) ||
	  (id[0]==ParticleID::piplus  && id[1]==ParticleID::K0)      ||
	  (id[0]==ParticleID::K0      && id[1]==ParticleID::piplus))
    {allowed=true;}
  // two kaons
  else if((id[0]==ParticleID::Kminus && id[1]==ParticleID::K0)     ||
	  (id[0]==ParticleID::K0     && id[1]==ParticleID::Kminus) ||
	  (id[0]==ParticleID::Kplus  && id[1]==ParticleID::Kbar0)  ||
	  (id[0]==ParticleID::Kbar0  && id[1]==ParticleID::Kplus))
    {allowed=true;}
  // charged kaon and eta
  else if((id[0]==ParticleID::Kminus && id[1]==ParticleID::eta)    ||
	  (id[0]==ParticleID::eta    && id[1]==ParticleID::Kminus) ||
	  (id[0]==ParticleID::Kplus  && id[1]==ParticleID::eta)    ||
	  (id[0]==ParticleID::eta    && id[1]==ParticleID::Kplus))
    {allowed=true;}
  return allowed;
}

// the decay mode
unsigned int TwoMesonRhoKStarCurrent::decayMode(vector<int> idout)
{
  unsigned int imode(0),nkaon(0);
  for(unsigned int ix=0;ix<idout.size();++ix)
    {
      if(idout[ix]==ParticleID::K0){imode=1;++nkaon;}
      else if (abs(idout[ix])==ParticleID::Kplus){imode=2;++nkaon;}
      else if (idout[ix]==ParticleID::eta){imode=4;}
    }
  if(nkaon==2){imode=3;}
  return imode;
}

// output the information for the database
void TwoMesonRhoKStarCurrent::dataBaseOutput(ofstream & output)
{
  output << "create /Herwig++/TwoMesonRhoKStarCurrent " << fullName() << " \n";
  for(unsigned int ix=0;ix<_rhomasses.size();++ix)
    {output << "insert " << fullName() << ":RhoMasses " << ix 
	    << " " << _rhomasses[ix] << "\n";}
  for(unsigned int ix=0;ix<_rhowidths.size();++ix)
    {output << "insert " << fullName() << ":RhoWidths " << ix 
	    << " " << _rhowidths[ix] << "\n";}
  for(unsigned int ix=0;ix<_Kstarmasses.size();++ix)
    {output << "insert " << fullName() << ":KstarMasses " << ix 
	    << " " << _Kstarmasses[ix] << "\n";}
  for(unsigned int ix=0;ix<_Kstarwidths.size();++ix)
    {output << "insert " << fullName() << ":KstarWidths " << ix 
	    << " " << _Kstarwidths[ix] << "\n";}
  output << "set " << fullName() << ":RhoParameters " << _rhoparameters << "\n";
  output << "set " << fullName() << ":KstarParameters " << _Kstarparameters << "\n";
  for(unsigned int ix=0;ix<_piwgt.size();++ix)
    {output << "insert " << fullName() << ":PiWeight " << ix 
	    << " " << _piwgt[ix] << "\n";}
  for(unsigned int ix=0;ix<_kwgt.size();++ix)
    {output << "insert " << fullName() << ":KWeight " << ix 
	    << " " << _kwgt[ix] << "\n";}
  output << "set " << fullName() << ":PiModel " << _pimodel << "\n";
  output << "set " << fullName() << ":KModel  " << _Kmodel << "\n";
}

}
