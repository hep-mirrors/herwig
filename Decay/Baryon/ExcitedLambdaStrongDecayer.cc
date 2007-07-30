// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ExcitedLambdaStrongDecayer class.
//

#include "ExcitedLambdaStrongDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/RSSpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/RSSpinorBarWaveFunction.h"
#include "Herwig++/Utilities/Kinematics.h"
#include "Herwig++/PDT/WidthCalculatorBase.h"

using namespace Herwig;
using namespace ThePEG::Helicity;
using namespace ThePEG::Helicity;

ExcitedLambdaStrongDecayer::ExcitedLambdaStrongDecayer() {
  // pion decay constant
  _fpi = 132*MeV;
  // g_2 coupling.
  _g2=0.570;
  // h_2 coupling.
  _h2=0.572;
  // h_8 coupling
  _h8=3.59e-3/MeV;
  // Lambda_c1
  _incoming.push_back(14122);_outgoing.push_back(4122);_charged.push_back(2);
  _incoming.push_back(14122);_outgoing.push_back(4122);_charged.push_back(0);
  // Lambda_c1*
  _incoming.push_back(4124);_outgoing.push_back(4122);_charged.push_back(2);
  _incoming.push_back(4124);_outgoing.push_back(4122);_charged.push_back(0);
  // Lambda_b1
  _incoming.push_back(15122);_outgoing.push_back(5122);_charged.push_back(2);
  _incoming.push_back(15122);_outgoing.push_back(5122);_charged.push_back(0);
  // Lambda_b1*
  _incoming.push_back(5124);_outgoing.push_back(5122);_charged.push_back(2);
  _incoming.push_back(5124);_outgoing.push_back(5122);_charged.push_back(0);
  // Xi_c1+ to Xi_c+ pi+ pi-
  _incoming.push_back(14322);_outgoing.push_back(4232);_charged.push_back(2);
  _incoming.push_back(14324);_outgoing.push_back(4232);_charged.push_back(2);
  // Xi_c1+ to Xi_c0 pi+ pi0
  _incoming.push_back(14322);_outgoing.push_back(4132);_charged.push_back(1);
  _incoming.push_back(14324);_outgoing.push_back(4132);_charged.push_back(1);
  // Xi_c1+ to Xi_c+ pi0 pi0
  _incoming.push_back(14322);_outgoing.push_back(4232);_charged.push_back(0);
  _incoming.push_back(14324);_outgoing.push_back(4232);_charged.push_back(0);
  // Xi_c10 to Xi_c0 pi+ pi-
  _incoming.push_back(14312);_outgoing.push_back(4132);_charged.push_back(2);
  _incoming.push_back(14314);_outgoing.push_back(4132);_charged.push_back(2);
  // Xi_c10 to Xi_c+ pi- pi0
  _incoming.push_back(14312);_outgoing.push_back(4232);_charged.push_back(1);
  _incoming.push_back(14314);_outgoing.push_back(4232);_charged.push_back(1);
  // Xi_c10 to Xi_c0 pi0 pi0
  _incoming.push_back(14312);_outgoing.push_back(4132);_charged.push_back(0);
  _incoming.push_back(14314);_outgoing.push_back(4132);_charged.push_back(0);
  // Xi_b10 to Xi_b0 pi+ pi-
  _incoming.push_back(15322);_outgoing.push_back(5232);_charged.push_back(2);
  _incoming.push_back(15324);_outgoing.push_back(5232);_charged.push_back(2);
  // Xi_b10 to Xi_b- pi+ pi0
  _incoming.push_back(15322);_outgoing.push_back(5132);_charged.push_back(1);
  _incoming.push_back(15324);_outgoing.push_back(5132);_charged.push_back(1);
  // Xi_b10 to Xi_b0 pi0 pi0
  _incoming.push_back(15322);_outgoing.push_back(5232);_charged.push_back(0);
  _incoming.push_back(15324);_outgoing.push_back(5232);_charged.push_back(0);
  // Xi_b1- to Xi_b- pi+ pi-
  _incoming.push_back(15312);_outgoing.push_back(5132);_charged.push_back(2);
  _incoming.push_back(15314);_outgoing.push_back(5132);_charged.push_back(2);
  // Xi_b1- to Xi_b0 pi- pi0
  _incoming.push_back(15312);_outgoing.push_back(5232);_charged.push_back(1);
  _incoming.push_back(15314);_outgoing.push_back(5232);_charged.push_back(1);
  // Xi_b1- to Xi_b- pi0 pi0
  _incoming.push_back(15312);_outgoing.push_back(5132);_charged.push_back(0);
  _incoming.push_back(15314);_outgoing.push_back(5132);_charged.push_back(0);
}

void ExcitedLambdaStrongDecayer::doinit() throw(InitException) {
  DecayIntegrator::doinit();
  // check consistency of the parameters
  unsigned int isize(_incoming.size());
  if(isize!=_outgoing.size()||isize!=_charged.size()) {
    throw InitException() << "Inconsistent parameters in ExcitedLambdaStrongDecayer"
			  << "::doinit()" << Exception::abortnow;
  }
  // set-up for the phase-space
  PDVector extpart(4);
  DecayPhaseSpaceModePtr mode;
  DecayPhaseSpaceChannelPtr channel;
  vector<double>::iterator start,end;
  double maxweight;
  vector<double> channelwgts;
  // particle data for the pions
  tPDPtr pip(getParticleData(ParticleID::piplus)),
    pi0(getParticleData(ParticleID::pi0)),pim(getParticleData(ParticleID::piminus));
  // particle data for the sigma_c's
  tPDPtr sigmac[3] ={getParticleData(4222),getParticleData(4212),getParticleData(4112)};
  tPDPtr sigmacs[3]={getParticleData(4224),getParticleData(4214),getParticleData(4114)};
  // particle data for the sigma_b's
  tPDPtr sigmab[3] ={getParticleData(5222),getParticleData(5212),getParticleData(5112)};
  tPDPtr sigmabs[3]={getParticleData(5224),getParticleData(5214),getParticleData(5114)};
  // particle data for xi_c's
  tPDPtr xic[2]  = {getParticleData(5322),getParticleData(5312)};
  tPDPtr xics[2] = {getParticleData(5324),getParticleData(5314)};
  // particle data for xi_b's
  tPDPtr xib[2]  = {getParticleData(5322),getParticleData(5312)};
  tPDPtr xibs[2] = {getParticleData(5324),getParticleData(5314)};
  // create the phase-space channels
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    // external particles for the mode
    extpart[0]=getParticleData(_incoming[ix]);
    extpart[1]=getParticleData(_outgoing[ix]);
    if(_charged[ix]==2) {
      extpart[2]=pip;
      extpart[3]=pim;
    }
    else if(_charged[ix]==1) {
      if(extpart[0]->iCharge()-extpart[1]->iCharge()==3) {
	extpart[2]=pip;
	extpart[3]=pi0;
      }
      else {
	extpart[2]=pim;
	extpart[3]=pi0;
      }
    }
    else {
      extpart[2]=pi0;
      extpart[3]=pi0;
    }
    // create the mode
    mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
    // channels for excited lambda_c's
    if(_incoming[ix]==14122||_incoming[ix]==4124) {
      // channels for charged pions
      if(_charged[ix]==2) {
	// sigma_c0(*) channels 
	channel=new_ptr(DecayPhaseSpaceChannel(mode));
	channel->addIntermediate(extpart[0],0,0.0,2,-1);
	channel->addIntermediate(sigmac[2],0,0.0,1,3);
	mode->addChannel(channel);
	channel=new_ptr(DecayPhaseSpaceChannel(mode));
	channel->addIntermediate(extpart[0],0,0.0,2,-1);
	channel->addIntermediate(sigmacs[2],0,0.0,1,3);
	mode->addChannel(channel);
	// sigma_c++(*) channels
	channel=new_ptr(DecayPhaseSpaceChannel(mode));
	channel->addIntermediate(extpart[0],0,0.0,3,-1);
	channel->addIntermediate(sigmac[0],0,0.0,1,2);
	mode->addChannel(channel);
	channel=new_ptr(DecayPhaseSpaceChannel(mode));
	channel->addIntermediate(extpart[0],0,0.0,3,-1);
	channel->addIntermediate(sigmacs[0],0,0.0,1,2);
	mode->addChannel(channel);
      }
      // channels for neutral pions
      else {
	// sigma_c+(*) channels 
	channel=new_ptr(DecayPhaseSpaceChannel(mode));
	channel->addIntermediate(extpart[0],0,0.0,2,-1);
	channel->addIntermediate(sigmac[1],0,0.0,1,3);
	mode->addChannel(channel);
	channel=new_ptr(DecayPhaseSpaceChannel(mode));
	channel->addIntermediate(extpart[0],0,0.0,2,-1);
	channel->addIntermediate(sigmacs[1],0,0.0,1,3);
	mode->addChannel(channel);
	// sigma_c+(*) channels (exchange of neutral pions)
	channel=new_ptr(DecayPhaseSpaceChannel(mode));
	channel->addIntermediate(extpart[0],0,0.0,3,-1);
	channel->addIntermediate(sigmac[1],0,0.0,1,2);
	mode->addChannel(channel);
	channel=new_ptr(DecayPhaseSpaceChannel(mode));
	channel->addIntermediate(extpart[0],0,0.0,3,-1);
	channel->addIntermediate(sigmacs[1],0,0.0,1,2);
	mode->addChannel(channel);
      }
    }
    // channels for excited lambda_b's
    else if(_incoming[ix]==15122||_incoming[ix]==5124) {
      // channels for charged pions
      if(_charged[ix]==2) {
	// sigma_b-(*) channels 
	channel=new_ptr(DecayPhaseSpaceChannel(mode));
	channel->addIntermediate(extpart[0],0,0.0,2,-1);
	channel->addIntermediate(sigmab[2],0,0.0,1,3);
	mode->addChannel(channel);
	channel=new_ptr(DecayPhaseSpaceChannel(mode));
	channel->addIntermediate(extpart[0],0,0.0,2,-1);
	channel->addIntermediate(sigmabs[2],0,0.0,1,3);
	mode->addChannel(channel);
	// sigma_b+(*) channels
	channel=new_ptr(DecayPhaseSpaceChannel(mode));
	channel->addIntermediate(extpart[0],0,0.0,3,-1);
	channel->addIntermediate(sigmab[0],0,0.0,1,2);
	mode->addChannel(channel);
	channel=new_ptr(DecayPhaseSpaceChannel(mode));
	channel->addIntermediate(extpart[0],0,0.0,3,-1);
	channel->addIntermediate(sigmabs[0],0,0.0,1,2);
	mode->addChannel(channel);
      }
      // channels for neutral pions
      else {
	// sigma_b0(*) channels 
	channel=new_ptr(DecayPhaseSpaceChannel(mode));
	channel->addIntermediate(extpart[0],0,0.0,2,-1);
	channel->addIntermediate(sigmab[1],0,0.0,1,3);
	mode->addChannel(channel);
	channel=new_ptr(DecayPhaseSpaceChannel(mode));
	channel->addIntermediate(extpart[0],0,0.0,2,-1);
	channel->addIntermediate(sigmabs[1],0,0.0,1,3);
	mode->addChannel(channel);
	// sigma_b0(*) channels (exchange of neutral pions)
	channel=new_ptr(DecayPhaseSpaceChannel(mode));
	channel->addIntermediate(extpart[0],0,0.0,3,-1);
	channel->addIntermediate(sigmab[1],0,0.0,1,2);
	mode->addChannel(channel);
	channel=new_ptr(DecayPhaseSpaceChannel(mode));
	channel->addIntermediate(extpart[0],0,0.0,3,-1);
	channel->addIntermediate(sigmabs[1],0,0.0,1,2);
	mode->addChannel(channel);
      }
    }
    else if(_incoming[ix]==14322||_incoming[ix]==14324) {
      // channels for two charged pions
      if(_charged[ix]==2) {
	channel=new_ptr(DecayPhaseSpaceChannel(mode));
	channel->addIntermediate(extpart[0],0,0.0,2,-1);
	if(_incoming[ix]==14322) {
	  channel->addIntermediate(xic[1],0,0.0,1,3);
	}
	else {
	  channel->addIntermediate(xics[1],0,0.0,1,3);
	}
	mode->addChannel(channel);
      }
      // channels for one charged pion
      else if(_charged[ix]==1) {
	channel=new_ptr(DecayPhaseSpaceChannel(mode));
	channel->addIntermediate(extpart[0],0,0.0,2,-1);
	if(_incoming[ix]==14322) {
	  channel->addIntermediate(xic[1],0,0.0,1,3);
	}
	else {
	  channel->addIntermediate(xics[1],0,0.0,1,3);
	}
	mode->addChannel(channel);
	channel=new_ptr(DecayPhaseSpaceChannel(mode));
	channel->addIntermediate(extpart[0],0,0.0,3,-1);
	if(_incoming[ix]==14322) {
	  channel->addIntermediate(xic[0],0,0.0,1,2);
	}
	else {
	  channel->addIntermediate(xics[0],0,0.0,1,2);
	}
	mode->addChannel(channel);
      }  
      // channels for two neutral pions
      else {
	channel=new_ptr(DecayPhaseSpaceChannel(mode));
	channel->addIntermediate(extpart[0],0,0.0,2,-1);
	if(_incoming[ix]==14322) {
	  channel->addIntermediate(xic[0],0,0.0,1,3);
	}
	else {
	  channel->addIntermediate(xics[0],0,0.0,1,3);
	}
	mode->addChannel(channel);
	channel=new_ptr(DecayPhaseSpaceChannel(mode));
	channel->addIntermediate(extpart[0],0,0.0,3,-1);
	if(_incoming[ix]==14322) {
	  channel->addIntermediate(xic[0],0,0.0,1,2);
	}
	else {
	  channel->addIntermediate(xics[0],0,0.0,1,2);
	}
	mode->addChannel(channel);
      }
    }
    else if(_incoming[ix]==14312||_incoming[ix]==14314) {
      // channels for two charged pions
      if(_charged[ix]==2) {
	channel=new_ptr(DecayPhaseSpaceChannel(mode));
	channel->addIntermediate(extpart[0],0,0.0,3,-1);
	if(_incoming[ix]==14312) {
	  channel->addIntermediate(xic[0],0,0.0,1,2);
	}
	else {
	  channel->addIntermediate(xics[0],0,0.0,1,2);
	}
	mode->addChannel(channel);
      }
      // channels for one charged pion
      else if(_charged[ix]==1) {
	channel=new_ptr(DecayPhaseSpaceChannel(mode));
	channel->addIntermediate(extpart[0],0,0.0,2,-1);
	if(_incoming[ix]==14312) {
	  channel->addIntermediate(xic[0],0,0.0,1,3);
	}
	else {
	  channel->addIntermediate(xics[0],0,0.0,1,3);
	}
	mode->addChannel(channel);
	channel=new_ptr(DecayPhaseSpaceChannel(mode));
	channel->addIntermediate(extpart[0],0,0.0,3,-1);
	if(_incoming[ix]==14312) {
	  channel->addIntermediate(xic[1],0,0.0,1,2);
	}
	else {
	  channel->addIntermediate(xics[1],0,0.0,1,2);
	}
	mode->addChannel(channel);
      }
      // channels for two neutral pions
      else {
	channel=new_ptr(DecayPhaseSpaceChannel(mode));
	channel->addIntermediate(extpart[0],0,0.0,2,-1);
	if(_incoming[ix]==14312) {
	  channel->addIntermediate(xic[1],0,0.0,1,3);
	}
	else {
	  channel->addIntermediate(xics[1],0,0.0,1,3);
	}
	mode->addChannel(channel);
	channel=new_ptr(DecayPhaseSpaceChannel(mode));
	channel->addIntermediate(extpart[0],0,0.0,3,-1);
	if(_incoming[ix]==14312) {
	  channel->addIntermediate(xic[1],0,0.0,1,2);
	}
	else {
	  channel->addIntermediate(xics[1],0,0.0,1,2);
	}
	mode->addChannel(channel);
      }
    }
    else if(_incoming[ix]==15322||_incoming[ix]==15324) {
      // channels for two charged pions
      if(_charged[ix]==2) {
	channel=new_ptr(DecayPhaseSpaceChannel(mode));
	channel->addIntermediate(extpart[0],0,0.0,2,-1);
	if(_incoming[ix]==15322) {
	  channel->addIntermediate(xib[1],0,0.0,1,3);
	}
	else {
	  channel->addIntermediate(xibs[1],0,0.0,1,3);
	}
	mode->addChannel(channel);
      }
      // channels for one charged pion
      else if(_charged[ix]==1) {
	channel=new_ptr(DecayPhaseSpaceChannel(mode));
	channel->addIntermediate(extpart[0],0,0.0,2,-1);
	if(_incoming[ix]==15322) {
	  channel->addIntermediate(xib[1],0,0.0,1,3);
	  }
	else {
	  channel->addIntermediate(xibs[1],0,0.0,1,3);
	}
	mode->addChannel(channel);
	channel=new_ptr(DecayPhaseSpaceChannel(mode));
	channel->addIntermediate(extpart[0],0,0.0,3,-1);
	if(_incoming[ix]==15322) {
	  channel->addIntermediate(xib[0],0,0.0,1,2);
	}
	else {
	  channel->addIntermediate(xibs[0],0,0.0,1,2);
	}
	mode->addChannel(channel);
      }
      // channels for two neutral pions
      else {
	channel=new_ptr(DecayPhaseSpaceChannel(mode));
	channel->addIntermediate(extpart[0],0,0.0,2,-1);
	if(_incoming[ix]==15322) {
	  channel->addIntermediate(xib[0],0,0.0,1,3);
	}
	else {
	  channel->addIntermediate(xibs[0],0,0.0,1,3);
	}
	mode->addChannel(channel);
	channel=new_ptr(DecayPhaseSpaceChannel(mode));
	channel->addIntermediate(extpart[0],0,0.0,3,-1);
	if(_incoming[ix]==15322) {
	  channel->addIntermediate(xib[0],0,0.0,1,2);
	}
	  else {
	    channel->addIntermediate(xibs[0],0,0.0,1,2);
	  }
	mode->addChannel(channel);
      }
    }
    else if(_incoming[ix]==15312||_incoming[ix]==15314) {
      // channels for two charged pions
      if(_charged[ix]==2) {
	channel=new_ptr(DecayPhaseSpaceChannel(mode));
	channel->addIntermediate(extpart[0],0,0.0,3,-1);
	if(_incoming[ix]==15312) {
	  channel->addIntermediate(xib[0],0,0.0,1,2);
	}
	else {
	  channel->addIntermediate(xibs[0],0,0.0,1,2);
	}
	mode->addChannel(channel);
      }
      // channels for one charged pion
      else if(_charged[ix]==1) {
	channel=new_ptr(DecayPhaseSpaceChannel(mode));
	channel->addIntermediate(extpart[0],0,0.0,2,-1);
	if(_incoming[ix]==15312) {
	  channel->addIntermediate(xib[0],0,0.0,1,3);
	  }
	else {
	  channel->addIntermediate(xibs[0],0,0.0,1,3);
	}
	mode->addChannel(channel);
	channel=new_ptr(DecayPhaseSpaceChannel(mode));
	channel->addIntermediate(extpart[0],0,0.0,3,-1);
	if(_incoming[ix]==15312) {
	  channel->addIntermediate(xib[1],0,0.0,1,2);
	}
	else {
	  channel->addIntermediate(xibs[1],0,0.0,1,2);
	}
	mode->addChannel(channel);
      }
      // channels for two neutral pions
      else {
	channel=new_ptr(DecayPhaseSpaceChannel(mode));
	channel->addIntermediate(extpart[0],0,0.0,2,-1);
	if(_incoming[ix]==15312) {
	  channel->addIntermediate(xib[1],0,0.0,1,3);
	}
	else {
	  channel->addIntermediate(xibs[1],0,0.0,1,3);
	}
	mode->addChannel(channel);
	channel=new_ptr(DecayPhaseSpaceChannel(mode));
	channel->addIntermediate(extpart[0],0,0.0,3,-1);
	if(_incoming[ix]==15312) {
	  channel->addIntermediate(xib[1],0,0.0,1,2);
	}
	else {
	  channel->addIntermediate(xibs[1],0,0.0,1,2);
	}
	mode->addChannel(channel);
      }
    }
    else {
      throw InitException() << "Unknown mode in ExcitedLambdaStrongDecayer::doinit()"
			    << Exception::abortnow;
    }
    // add the mode
    if(_wgtmax.size()>numberModes()) maxweight=_wgtmax[numberModes()];
    else                             maxweight=0.;
    // the weights for the channel
    if(_wgtloc.size()>numberModes()&&
       _wgtloc[numberModes()]+mode->numberChannels()<=_weights.size()) {
      start=_weights.begin()+_wgtloc[numberModes()];
      end  = start+mode->numberChannels();
      channelwgts=vector<double>(start,end);
    }
    else {
      channelwgts.resize(mode->numberChannels(),1./(mode->numberChannels()));
    }
    // finally add the mode
    addMode(mode,maxweight,channelwgts);
  }
  generateIntermediates(true);
}

void ExcitedLambdaStrongDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(_fpi,GeV) << _g2 << _h2 << ounit(_h8,1./GeV) 
     << _incoming << _outgoing << _charged << _wgtloc 
     << _wgtmax << _weights;
}

void ExcitedLambdaStrongDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_fpi,GeV) >> _g2 >> _h2 >> iunit(_h8,1./GeV) 
     >> _incoming >> _outgoing >> _charged >> _wgtloc 
     >> _wgtmax >> _weights;
}

ClassDescription<ExcitedLambdaStrongDecayer> ExcitedLambdaStrongDecayer::initExcitedLambdaStrongDecayer;
// Definition of the static class description member.

void ExcitedLambdaStrongDecayer::Init() {
  static ClassDocumentation<ExcitedLambdaStrongDecayer> documentation
    ("The ExcitedLambdaStrongDecayer class is designed for the decay of excited Lambda"
     " baryons containing a heavy quark to two pions and the lightest Lambda of that"
     " species");

  static Parameter<ExcitedLambdaStrongDecayer,Energy> interfaceFpi
    ("Fpi",
     "The pion decay constant.",
     &ExcitedLambdaStrongDecayer::_fpi, MeV, 130.7*MeV, 100.0*MeV, 200.0*MeV,
     false, false, true);

  static Parameter<ExcitedLambdaStrongDecayer,double> interfaceg2
    ("g2",
     "The g2 coupling",
     &ExcitedLambdaStrongDecayer::_g2, 0.53852, -10.0, 10.0,
     false, false, true);

  static Parameter<ExcitedLambdaStrongDecayer,double> interfaceh2
    ("h2",
     "The h2 coupling",
     &ExcitedLambdaStrongDecayer::_h2, 0.4899, -10.0, 10.0,
     false, false, true);

  static Parameter<ExcitedLambdaStrongDecayer,InvEnergy> interfaceh8
    ("h8",
     "The h8 coupling",
     &ExcitedLambdaStrongDecayer::_h8, 1./GeV, 0./GeV, -10./GeV, 10./GeV,
     false, false, true);

  static ParVector<ExcitedLambdaStrongDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code of the incoming baryon",
     &ExcitedLambdaStrongDecayer::_incoming, -1, 0, 0, 1000000,
     false, false, true);

  static ParVector<ExcitedLambdaStrongDecayer,int> interfaceOutgoing
    ("Outgoing",
     "The PDG code of the outgoing baryon",
     &ExcitedLambdaStrongDecayer::_outgoing, -1, 0, 0, 1000000,
     false, false, true);

  static ParVector<ExcitedLambdaStrongDecayer,int> interfaceWeightLocation
    ("WeightLocation",
     "The locations of the weights for a given channel in the vector",
     &ExcitedLambdaStrongDecayer::_wgtloc,
     0, 0, 0, 0, 10000, false, false, true);

  static ParVector<ExcitedLambdaStrongDecayer,double> interfaceWeightMax
    ("MaximumWeight",
     "The maximum weight for a given channel.",
     &ExcitedLambdaStrongDecayer::_wgtmax,
     0, 0, 0, 0., 100., false, false, true);

  static ParVector<ExcitedLambdaStrongDecayer,double> interfaceWeights
    ("Weights",
     "The weights for the integration.",
     &ExcitedLambdaStrongDecayer::_weights,
     0, 0, 0, 0., 1., false, false, true);

  static ParVector<ExcitedLambdaStrongDecayer,int> interfaceCharged
    ("Charged",
     "Number of charged pions produced in the decay",
     &ExcitedLambdaStrongDecayer::_charged, -1, 0, 0, 2,
     false, false, true);

}


// the matrix element
double ExcitedLambdaStrongDecayer::me2(bool vertex, const int ichan,
				       const Particle & part,
				       const ParticleVector & decay) const {
  vector<LorentzSpinor<SqrtEnergy> > sp;
  vector<LorentzRSSpinor<SqrtEnergy> > Rsp;
  vector<LorentzSpinorBar<SqrtEnergy> > sbar;
  vector<LorentzRSSpinorBar<SqrtEnergy> > Rsbar;
  RhoDMatrix temp;
  // workaround for gcc 3.2.3 bug
  // construct or obtain the spin information
  //ALB ScalarWaveFunction(decay[1],outgoing,true,vertex);
  //ALB ScalarWaveFunction(decay[2],outgoing,true,vertex);
  PPtr mytemp1 = decay[1];
  ScalarWaveFunction(mytemp1,outgoing,true,vertex);
  PPtr mytemp2 = decay[2];
  ScalarWaveFunction(mytemp2,outgoing,true,vertex);

  if(part.id()>0) {
    if(part.dataPtr()->iSpin()==2) {
      SpinorWaveFunction(sp,temp,const_ptr_cast<tPPtr>(&part),
			 incoming,true,vertex);
    }
    else {
      RSSpinorWaveFunction(Rsp,temp,const_ptr_cast<tPPtr>(&part),
			   incoming,true,vertex);
    }
    SpinorBarWaveFunction(sbar,decay[0],outgoing,true,vertex);
  }
  else {
    if(part.dataPtr()->iSpin()==2) {
      SpinorBarWaveFunction(sbar,temp,const_ptr_cast<tPPtr>(&part),incoming,true,
			    vertex);
    }
    else {
      RSSpinorBarWaveFunction(Rsbar,temp,const_ptr_cast<tPPtr>(&part),
			      incoming,true,vertex);
    }
    SpinorWaveFunction(sp,decay[0],outgoing,true,vertex);
  }
  // calculate the widths of the different intermediates
  // production of Lambda_c+
  tcPDPtr res[6];
  Energy gam[6];
  if(abs(decay[0]->id())==ParticleID::Lambda_cplus) {
    res[0]=getParticleData(4222);
    res[1]=getParticleData(4212);
    res[2]=getParticleData(4112);
    res[3]=getParticleData(4224);
    res[4]=getParticleData(4214);
    res[5]=getParticleData(4114);
  }
  // production of Lambda_b+
  else if(abs(decay[0]->id())==ParticleID::Lambda_b0) {
    res[0]=getParticleData(5222);
    res[1]=getParticleData(5212);
    res[2]=getParticleData(5112);
    res[3]=getParticleData(5224);
    res[4]=getParticleData(5214);
    res[5]=getParticleData(5114);
  }
  // production of Xi_c
  else if(abs(decay[0]->id())==ParticleID::Xi_cplus||
	  abs(decay[0]->id())==ParticleID::Xi_c0) {
    res[0]=getParticleData(4322);
    res[1]=getParticleData(4312);
    res[2]=getParticleData(4324);
    res[3]=getParticleData(4314);
    res[4]=getParticleData(4324);
    res[5]=getParticleData(4314);
  }
  // production of Xi_b
  else if(abs(decay[0]->id())==ParticleID::Xi_b0||
	  abs(decay[0]->id())==ParticleID::Xi_bminus) {
    res[0]=getParticleData(5322);
    res[1]=getParticleData(5312);
    res[2]=getParticleData(5324);
    res[3]=getParticleData(5314);
    res[4]=getParticleData(5324);
    res[5]=getParticleData(5314);
  }
  // compute the widths
  Energy ppi,mpi(decay[1]->mass());
  for(unsigned int ix=0;ix<6;++ix) {
    ppi = Kinematics::pstarTwoBodyDecay(res[ix]->mass(),mpi,decay[0]->mass());
    gam[ix] = 0.5*_g2*_g2*decay[0]->mass()*ppi*ppi*ppi/Constants::pi/_fpi/_fpi/res[ix]->mass();
  }
  double output;
  complex<InvEnergy> prop[4];
  Complex ii(0.,1.);
  Energy delta(part.mass()-decay[0]->mass()),E1(decay[1]->momentum().e()),
    E2(decay[2]->momentum().e()),mout(decay[0]->mass());
  Energy2  p12(E1*E1-mpi*mpi),p22(E2*E2-mpi*mpi);
  Energy2 dot12(decay[1]->momentum()*decay[2]->momentum());
  if(part.dataPtr()->iSpin()==2) {
    DecayMatrixElement newME(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin0,PDT::Spin0);
    Complex A(0.),B(0.);
    // coefficients for the matrix element
    if(abs(decay[0]->id())==ParticleID::Lambda_cplus||
       abs(decay[0]->id())==ParticleID::Lambda_b0) {
      // propagators
      if(abs(decay[1]->id())==ParticleID::piplus) {
	prop[0] = 1./(delta-res[2]->mass()+mout-E1+ii*0.5*gam[2]);
	prop[1] = 1./(delta-res[0]->mass()+mout-E2+ii*0.5*gam[0]);
	prop[2] = 1./(delta-res[5]->mass()+mout-E1+ii*0.5*gam[5]);
	prop[3] = 1./(delta-res[3]->mass()+mout-E2+ii*0.5*gam[3]);
      }
      else {
	prop[0] = 1./(delta-res[1]->mass()+mout-E1+ii*0.5*gam[1]);
	prop[1] = 1./(delta-res[1]->mass()+mout-E2+ii*0.5*gam[1]);
	prop[2] = 1./(delta-res[4]->mass()+mout-E1+ii*0.5*gam[4]);
	prop[3] = 1./(delta-res[4]->mass()+mout-E2+ii*0.5*gam[4]);
      }
      if(ichan==0) {
	prop[1]=0./MeV;
	prop[2]=0./MeV;
	prop[3]=0./MeV;
      }
      else if(ichan==1) {
	prop[0]=0./MeV;
	prop[1]=0./MeV;
	prop[3]=0./MeV;
      }
      else if(ichan==2) {
	prop[0]=0./MeV;
	prop[2]=0./MeV;
	prop[3]=0./MeV;
      }
      else if(ichan==3) {
	prop[0]=0./MeV;
	prop[1]=0./MeV;
	prop[2]=0./MeV;
      }
      // the coefficients
      A = _h2*E1*prop[0]-2.*_h8*p12/3.*prop[2]+2.*_h8*(E1*E2-dot12)*prop[3];
      B = _h2*E2*prop[1]-2.*_h8*p22/3.*prop[3]+2.*_h8*(E1*E2-dot12)*prop[2];
    }
    // Xi_c1+ to pi+ pi- Xi_c1+
    else if((abs(decay[0]->id())==ParticleID::Xi_cplus||
	     abs(decay[0]->id())==ParticleID::Xi_b0     )&&_charged[imode()]==2) {
      prop[0] = 1./(delta-res[1]->mass()+mout-E1+ii*0.5*gam[1]);
      A = 0.5*_h2*E1*prop[0];
      B = 0.;
    }
    // Xi_c10 to pi+ pi- Xi_c10
    else if((abs(decay[0]->id())==ParticleID::Xi_c0||
	     abs(decay[0]->id())==ParticleID::Xi_bminus)&&_charged[imode()]==2) {
      prop[0] = 1./(delta-res[0]->mass()+mout-E2+ii*0.5*gam[0]);
      A = 0.;
      B = 0.5*_h2*E2*prop[0];
    }
    // Xi_c1+ to pi+ pi0 Xi_c10
    else if((abs(decay[0]->id())==ParticleID::Xi_c0||
	     abs(decay[0]->id())==ParticleID::Xi_bminus)&&_charged[imode()]==1) {
      prop[0] = 1./(delta-res[1]->mass()+mout-E1+ii*0.5*gam[1]);
      prop[1] = 1./(delta-res[0]->mass()+mout-E2+ii*0.5*gam[0]);
      if(ichan==0)      prop[1]=0./MeV;
      else if(ichan==1) prop[0]=0./MeV;
      A = 0.35355339*_h2*E1*prop[0];
      B = 0.35355339*_h2*E2*prop[1];
    }
    // Xi_c10 to pi- pi0 Xi_c+
    else if((abs(decay[0]->id())==ParticleID::Xi_cplus||
	     abs(decay[0]->id())==ParticleID::Xi_b0    )&&_charged[imode()]==1) {
      prop[0] = 1./(delta-res[0]->mass()+mout-E1+ii*0.5*gam[0]);
      prop[1] = 1./(delta-res[1]->mass()+mout-E2+ii*0.5*gam[1]);
      if(ichan==0)      prop[1]=0./MeV;
      else if(ichan==1) prop[0]=0./MeV;
      A = 0.35355339*_h2*E1*prop[0];
      B = 0.35355339*_h2*E2*prop[1];
    }
    // Xi_c1+ to pi0 pi0 Xi_c1+
    else if((abs(decay[0]->id())==ParticleID::Xi_cplus||
	     abs(decay[0]->id())==ParticleID::Xi_b0     )&&_charged[imode()]==0) {
      prop[0] = 1./(delta-res[0]->mass()+mout-E1+ii*0.5*gam[0]);
      prop[1] = 1./(delta-res[0]->mass()+mout-E2+ii*0.5*gam[0]);
      if(ichan==0)      prop[1]=0./MeV;
      else if(ichan==1) prop[0]=0./MeV;
      A = 0.25*_h2*E1*prop[0];
      B = 0.25*_h2*E2*prop[1];
    }
    // Xi_c1+ to pi0 pi0 Xi_c1+
    else if((abs(decay[0]->id())==ParticleID::Xi_c0||
	     abs(decay[0]->id())==ParticleID::Xi_bminus)&&_charged[imode()]==0) {
      prop[0] = 1./(delta-res[1]->mass()+mout-E1+ii*0.5*gam[1]);
      prop[1] = 1./(delta-res[1]->mass()+mout-E2+ii*0.5*gam[1]);
      if(ichan==0)      prop[1]=0./MeV;
      else if(ichan==1) prop[0]=0./MeV;
      A = 0.25*_h2*E1*prop[0];
      B = 0.25*_h2*E2*prop[1];
    }
    // compute the matrix element
    unsigned int ixa,iya;
    vector<unsigned int> ihel(4);ihel[2]=0;ihel[3]=0;
    LorentzPolarizationVectorE current;
    for(ixa=0;ixa<2;++ixa) {
      for(iya=0;iya<2;++iya) {
	if(part.id()>0) {
	  ihel[0]=ixa;
	  ihel[1]=iya;
	}
	else {
	  ihel[0]=iya;
	  ihel[1]=ixa;
	}
	current = sp[iya].generalCurrent(sbar[ixa],-1.,1.);
	newME(ihel) = (A*(current*decay[2]->momentum())
	  +B*(current*decay[1]->momentum()))/sqr(_fpi);
      }
    }
    output=(newME.contract(temp)).real()*_g2*_g2;
    /*
    // full matrix element
    Energy2 dot01(part.momentum()*decay[1]->momentum());
    Energy2 dot02(part.momentum()*decay[2]->momentum());
    Energy2 dot03(part.momentum()*decay[0]->momentum());
    Energy2 dot12(decay[1]->momentum()*decay[2]->momentum());
    Energy2 dot13(decay[1]->momentum()*decay[0]->momentum());
    Energy2 dot23(decay[2]->momentum()*decay[0]->momentum());
    Energy2 m12(decay[1]->mass()*decay[1]->mass()),
    m22(decay[2]->mass()*decay[2]->mass());
    Energy m3(decay[0]->mass()),m0(part.mass());
    double me = 
    4.*(A*conj(A)).real()*(2.*dot23*dot02-m22*(dot03+m0*m3))+
    4.*(B*conj(B)).real()*(2.*dot13*dot01-m12*(dot03+m0*m3))+
    8.*(A*conj(B)).real()*(dot23*dot01+dot13*dot02-dot12*(dot03+m0*m3));
    me *=0.5*_g2*_g2/_fpi/_fpi/_fpi/_fpi;
    //cout << "testing " <<   me/output << endl;
    */
    if(decay[1]->id()==decay[2]->id()){output*=0.5;}
    ME(newME);
  }
  else {
    DecayMatrixElement newME(PDT::Spin3Half,PDT::Spin1Half,PDT::Spin0,PDT::Spin0);
    Complex C(0.),E(0.);
    complex<InvEnergy2> D(0./MeV2),F(0./MeV2);
    if(abs(decay[0]->id())==ParticleID::Lambda_cplus||
       abs(decay[0]->id())==ParticleID::Lambda_b0) {
      // propagators
      if(abs(decay[1]->id())==ParticleID::piplus) {
	prop[0] = 1./(delta-res[2]->mass()+mout-E1+ii*0.5*gam[2]);
	prop[1] = 1./(delta-res[0]->mass()+mout-E2+ii*0.5*gam[0]);
	prop[2] = 1./(delta-res[5]->mass()+mout-E1+ii*0.5*gam[5]);
	prop[3] = 1./(delta-res[3]->mass()+mout-E2+ii*0.5*gam[3]);
      }
      else {
	prop[0] = 1./(delta-res[1]->mass()+mout-E1+ii*0.5*gam[1]);
	prop[1] = 1./(delta-res[1]->mass()+mout-E2+ii*0.5*gam[1]);
	prop[2] = 1./(delta-res[4]->mass()+mout-E1+ii*0.5*gam[4]);
	prop[3] = 1./(delta-res[4]->mass()+mout-E2+ii*0.5*gam[4]);
	    }
      if(ichan==0) {
	prop[1]=0./MeV;
	prop[2]=0./MeV;
	prop[3]=0./MeV;
      }
      else if(ichan==1) {
	prop[0]=0./MeV;
	prop[1]=0./MeV;
	prop[3]=0./MeV;
      }
      else if(ichan==2) {
	prop[0]=0./MeV;
	prop[2]=0./MeV;
	prop[3]=0./MeV;
      }
      else if(ichan==3) {
	prop[0]=0./MeV;
	prop[1]=0./MeV;
	prop[2]=0./MeV;
      }
      // coefficients
      C = (_h2*E2-2./3.*_h8*p22)*prop[3]
	+2./3.*_h8*(E1*E2-dot12)*(prop[0]+2.*prop[2]);
      D = 2./3.*_h8*(-prop[0]+prop[2]);
      E = (_h2*E1-2./3.*_h8*p12)*prop[2]
	+2./3.*_h8*(E1*E2-dot12)*(prop[1]+2.*prop[3]);
      F =-2./3.*_h8*(-prop[1]+prop[3]);	  
    }
    // Xi_c1+ to pi+ pi- Xi_c1+
    else if((abs(decay[0]->id())==ParticleID::Xi_cplus||
	     abs(decay[0]->id())==ParticleID::Xi_b0     )&&_charged[imode()]==2) {
      prop[0] = 1./(delta-res[3]->mass()+mout-E1+ii*0.5*gam[3]);
      E = 0.5*_h2*E1*prop[0];
      C = 0.;
    }
    // Xi_c10 to pi+ pi- Xi_c10
    else if((abs(decay[0]->id())==ParticleID::Xi_c0||
	     abs(decay[0]->id())==ParticleID::Xi_bminus)&&_charged[imode()]==2) {
      prop[0] = 1./(delta-res[2]->mass()+mout-E2+ii*0.5*gam[2]);
      E = 0.;
      C = 0.5*_h2*E2*prop[0];
    }
    // Xi_c1+ to pi+ pi0 Xi_c10
    else if((abs(decay[0]->id())==ParticleID::Xi_c0||
	     abs(decay[0]->id())==ParticleID::Xi_bminus)&&_charged[imode()]==1) {
      prop[0] = 1./(delta-res[3]->mass()+mout-E1+ii*0.5*gam[3]);
      prop[1] = 1./(delta-res[2]->mass()+mout-E2+ii*0.5*gam[2]);
      if(ichan==0)      prop[1]=0./MeV;
      else if(ichan==1) prop[0]=0./MeV;
      E = 0.35355339*_h2*E1*prop[0];
      C = 0.35355339*_h2*E2*prop[1];
    }
    // Xi_c10 to pi- pi0 Xi_c+
    else if((abs(decay[0]->id())==ParticleID::Xi_cplus||
	     abs(decay[0]->id())==ParticleID::Xi_b0    )&&_charged[imode()]==1) {
      prop[0] = 1./(delta-res[2]->mass()+mout-E1+ii*0.5*gam[2]);
      prop[1] = 1./(delta-res[3]->mass()+mout-E2+ii*0.5*gam[3]);
      if(ichan==0)      prop[1]=0./MeV;
      else if(ichan==1) prop[0]=0./MeV;
      E = 0.35355339*_h2*E1*prop[0];
      C = 0.35355339*_h2*E2*prop[1];
    }
    // Xi_c1+ to pi0 pi0 Xi_c1+
    else if((abs(decay[0]->id())==ParticleID::Xi_cplus||
	     abs(decay[0]->id())==ParticleID::Xi_b0     )&&_charged[imode()]==0) {
      prop[0] = 1./(delta-res[2]->mass()+mout-E1+ii*0.5*gam[2]);
      prop[1] = 1./(delta-res[2]->mass()+mout-E2+ii*0.5*gam[2]);
      if(ichan==0)      prop[1]=0./MeV;
      else if(ichan==1) prop[0]=0./MeV;
      E = 0.25*_h2*E1*prop[0];
      C = 0.25*_h2*E2*prop[1];
    }
    // Xi_c1+ to pi0 pi0 Xi_c1+
    else if((abs(decay[0]->id())==ParticleID::Xi_c0||
	     abs(decay[0]->id())==ParticleID::Xi_bminus)&&_charged[imode()]==0) {
      prop[0] = 1./(delta-res[3]->mass()+mout-E1+ii*0.5*gam[3]);
      prop[1] = 1./(delta-res[3]->mass()+mout-E2+ii*0.5*gam[3]);
      if(ichan==0)      prop[1]=0./MeV;
      else if(ichan==1) prop[0]=0./MeV;
      E = 0.25*_h2*E1*prop[0];
      C = 0.25*_h2*E2*prop[1];
    }
    vector<unsigned int> ihel(4);ihel[2]=0;ihel[3]=0;
    complex<Energy2> Cm,Em;
    complex<Energy4> Dm,Fm;
    Energy  msum(part.mass()+decay[0]->mass());
    Energy2 dot(decay[1]->momentum()*(part.momentum()+decay[0]->momentum()));
    if(part.id()>0) {
      LorentzSpinor<SqrtEnergy> sp1,sp2;
      for(ihel[0]=0;ihel[0]<4;++ihel[0]) {
	sp1 = Rsp[ihel[0]].dot(decay[1]->momentum());
	sp2 = Rsp[ihel[0]].dot(decay[2]->momentum());
	for(ihel[1]=0;ihel[1]<2;++ihel[1]) {
	  Cm = sp1.scalar(sbar[ihel[1]])*UnitRemoval::E;
 	  Dm = -complex<Energy4>(msum*(sp1.vectorCurrent(sbar[ihel[1]]))*decay[1]->momentum()*UnitRemoval::E)+complex<Energy4>(dot*Cm);
	  Em = sp2.scalar(sbar[ihel[1]])*UnitRemoval::E;
	  Fm = -complex<Energy4>(msum*(sp2.vectorCurrent(sbar[ihel[1]])).dot(decay[1]->momentum())*UnitRemoval::E)
	    +complex<Energy4>(dot*Em);
	  newME(ihel) = (C*Cm+D*Dm+E*Em+F*Fm)/sqr(_fpi);
	}
      }
    }
    else {
      LorentzSpinorBar<SqrtEnergy> sp1,sp2;
      for(ihel[0]=0;ihel[0]<4;++ihel[0]) {
	sp1 = Rsbar[ihel[0]].dot(decay[1]->momentum());
	sp2 = Rsbar[ihel[0]].dot(decay[2]->momentum());
	for(ihel[1]=0;ihel[1]<2;++ihel[1]) {
	  Cm = sp[ihel[1]].scalar(sp1)*UnitRemoval::E;
	  Dm = -complex<Energy4>(msum*(sp[ihel[1]].vectorCurrent(sp1))*decay[1]->momentum()*UnitRemoval::E)
	    +complex<Energy4>(dot*Cm);
	  Em = sp[ihel[1]].scalar(sp2)*UnitRemoval::E;
	  Fm =  complex<Energy4>(msum*(sp[ihel[1]].vectorCurrent(sp2))*decay[1]->momentum()*UnitRemoval::E)
	    -complex<Energy4>(dot*Em);
	  newME(ihel) = (C*Cm+D*Dm+E*Em+F*Fm)/sqr(_fpi);
	}
      }
    }
    // final calculation
    output=3.*(newME.contract(temp)).real()*_g2*_g2;
    ME(newME);
    /*
      double test = 4.*part.mass()*decay[0]->mass()/_fpi/_fpi/_fpi/_fpi*_g2*_g2*
      (p12*(C*conj(C)).real()+p22*(E*conj(E)).real()
      +2.*(E1*E2-dot12)*(C*conj(E)).real()
      +(p12*p22-(E1*E2-dot12)*(E1*E2-dot12))*(+p12*(D*conj(D)).real()
      +p22*(F*conj(F)).real()
      -(C*conj(F)).real()+(D*conj(E)).real()
      +2.*(E1*E2-dot12)*(D*conj(F)).real()));
      Energy2 dot01(part.momentum()*decay[1]->momentum());
      Energy2 dot02(part.momentum()*decay[2]->momentum());
      Energy2 dot03(part.momentum()*decay[0]->momentum());
      Energy2 dot12(decay[1]->momentum()*decay[2]->momentum());
      Energy2 dot13(decay[1]->momentum()*decay[0]->momentum());
      Energy2 dot23(decay[2]->momentum()*decay[0]->momentum());
      Energy2 m12(decay[1]->mass()*decay[1]->mass()),
      m22(decay[2]->mass()*decay[2]->mass());
      Energy m3(decay[0]->mass()),m0(part.mass());
      double me = 
	8.*p12*(dot03+m0*m3)*(C*conj(C)).real()+
	8.*p22*(dot03+m0*m3)*(E*conj(E)).real()+
	8.*(m0*m3*(m12*m22-dot12*dot12)-dot03*dot12*dot12
	    +2.*dot12*(dot01*dot23+dot02*dot13)-2.*m12*dot02*dot23-2.*m22*dot01*dot13
	    +m12*m22*dot03)*(p12*(D*conj(D)).real()+p22*(F*conj(F)).real())+
	16.*(dot01*dot23-dot13*dot02)*(p12*(C*conj(D)).real()+p22*(E*conj(F)).real())
	+8.*(2.*(m0*m3+dot03)*(E1*E2-dot12)+dot13*dot02-dot23*dot01
	     +m0*(-dot13*E2+dot23*E1))*(C*conj(E)).real()
	+8.*(2.*(E1*E2-dot12)*(dot01*dot23-dot02*dot13)
	     +dot12*dot12*dot03-2.*dot02*dot13*dot12-2.*dot01*dot23*dot12
	     +2.*m22*dot01*dot13+2.*m12*dot02*dot23-m12*m22*dot03
	     +m0*m3*(dot12*dot12-m12*m22)
	     -m3*(E2*(dot01*dot12-m12*dot02)+E1*(dot02*dot12-m22*dot01))
	     +m0*(E2*(dot13*dot12-m12*dot23)+E1*(dot23*dot12-m22*dot13))
	     )*(C*conj(F)).real()
	-8.*(-2.*(E1*E2-dot12)*(dot01*dot23-dot02*dot13)
	     +dot12*dot12*dot03-2.*dot02*dot13*dot12-2.*dot01*dot23*dot12
	     +2.*m22*dot01*dot13+2.*m12*dot02*dot23-m12*m22*dot03
	     +m0*m3*(dot12*dot12-m12*m22)
	     -m3*(E2*(dot01*dot12-m12*dot02)+E1*(dot02*dot12-m22*dot01))
	     +m0*(E2*(dot13*dot12-m12*dot23)+E1*(dot23*dot12-m22*dot13))
	     )*(D*conj(E)).real()
	+8.*( 2.*(E1*E2-dot12)*(m0*m3*(m12*m22-dot12*dot12)-2.*m12*dot02*dot23
				-2.*m22*dot01*dot13+m12*m22*dot03-dot03*dot12*dot12
				+2.*dot12*(dot01*dot23+dot02*dot13))
	      +4.*(dot12*dot12-m12*m22)*(m0*(E1*dot23-E2*dot13)
					 +dot02*dot13-dot01*dot23))*(D*conj(F)).real();
      me *=0.25*_g2*_g2/_fpi/_fpi/_fpi/_fpi;

      cout << "testingB " << me/output << endl;
    */
    if(decay[1]->id()==decay[2]->id()){output*=0.5;}
  }
  return output;
}

WidthCalculatorBasePtr 
ExcitedLambdaStrongDecayer::threeBodyMEIntegrator(const DecayMode & ) const {
  return WidthCalculatorBasePtr();
}

double ExcitedLambdaStrongDecayer::threeBodyMatrixElement(int ,Energy2 ,
							  Energy2 ,Energy2 ,
							  Energy2 ,Energy ,
							  Energy ,Energy ) {
  return -1;
}

void ExcitedLambdaStrongDecayer::dataBaseOutput(ofstream & , bool) const{
}

int ExcitedLambdaStrongDecayer::modeNumber(bool & cc,const DecayMode & dm) const {
  int idout(0),id,imode(-1);
  unsigned int npi0(0),ix(0);
  ParticleMSet::const_iterator pit(dm.products().begin());
  for( ;pit!=dm.products().end();++pit) {
    id=(**pit).id();
    if(id==ParticleID::pi0){++npi0;}
    else if(id!=ParticleID::piplus&&id!=ParticleID::piminus){idout=id;}
  }
  int charged(2-npi0);
  id=dm.parent()->id();
  do {
    if(id==_incoming[ix]&&idout==_outgoing[ix]&&charged==_charged[ix])
      {imode=ix;cc=false;}
    else if(id==-_incoming[ix]&&idout==-_outgoing[ix]&&charged==_charged[ix])
      {imode=ix;cc=true;}
    ++ix;
  }
  while(imode<0&&ix<_incoming.size());
  return imode;
}
