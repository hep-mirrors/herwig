// -*- C++ -*-
//
// TauDecayer2.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TauDecayer2 class.
//
//  Author: Peter Richardson
//

#include "TauDecayer2.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Helicity/VectorSpinInfo.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig/Decay/DecayVertex.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"
#include "ThePEG/Helicity/FermionSpinInfo.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void TauDecayer2::doinit() {
  DecayIntegrator2::doinit();
  assert(false);
  // // make sure the current got initialised
  // current_->init();
  // // set up the phase-space channels
  // tPDPtr tau = getParticleData(ParticleID::tauminus);
  // tPDPtr nu  = getParticleData(ParticleID::nu_tau);
  // Energy mtau(tau->mass());
  // double maxweight;
  // vector<double> channelwgts;
  // modeMap_.clear();
  // vector<double>::iterator start,end;
  // for(unsigned int ix=0;ix<current_->numberOfModes();++ix) {
  //   // get the external particles for this mode
  //   tPDVector extpart = {tau,nu};
  //   int iq(0),ia(0);
  //   tPDVector ptemp  = current_->particles(-3,ix,iq,ia);
  //   extpart.insert(std::end(extpart), std::begin(ptemp), std::end(ptemp));
  //   // create the mode
  //   DecayPhaseSpaceModePtr mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
  //   // create the first piece of the channel
  //   DecayPhaseSpaceChannelPtr channel = new_ptr(DecayPhaseSpaceChannel(mode));
  //   channel->addIntermediate(extpart[0],0,0.0,-1,1);
  //   if(!current_->createMode(-3,ix,mode,2,1,channel,mtau)) continue;
  //   // the maximum weight and the channel weights
  //   // the maximum
  //   maxweight = wgtMax_.size()>numberModes() ? wgtMax_[numberModes()] : 0;
  //   // the weights for the channel
  //   if(wgtLoc_.size()>numberModes()&&
  //      wgtLoc_[numberModes()]+mode->numberChannels()<=weights_.size()) {
  //     start=weights_.begin()+wgtLoc_[numberModes()];
  //     end  = start+mode->numberChannels();
  //     channelwgts=vector<double>(start,end);
  //   }
  //   else {
  //     channelwgts.resize(mode->numberChannels(),1./(mode->numberChannels()));
  //   }
  //   modeMap_.push_back(ix);
  //   // special for the two body modes
  //   if(extpart.size()==3) {
  //     channelwgts.clear();
  //     mode=new_ptr(DecayPhaseSpaceMode(extpart,this));
  //   }
  //   addMode(mode,maxweight,channelwgts);
  // }
  // current_->reset();
  // current_->touch();
  // current_->update();
}

void TauDecayer2::doinitrun() {
  assert(false);
  // current_->initrun();
  // DecayIntegrator2::doinitrun();
  // if(initialize()) {
  //   weights_.clear();wgtLoc_.clear();wgtMax_.clear();
  //   unsigned int ix,iy;
  //   for(ix=0;ix<numberModes();++ix) {
  //     wgtMax_.push_back(mode(ix)->maxWeight());
  //     wgtLoc_.push_back(weights_.size());
  //     for(iy=0;iy<mode(ix)->numberChannels();++iy) {
  // 	weights_.push_back(mode(ix)->channelWeight(iy));
  //     }
  //   }
  // }
}

bool TauDecayer2::accept(tcPDPtr parent, const tPDVector & children) const {
  bool allowed(false);
  // find the neutrino 
  int idnu(0),idtemp,idin(parent->id());
  vector<int> idother;
  tPDVector::const_iterator pit  = children.begin();
  tPDVector::const_iterator pend = children.end();
  for( ; pit!=pend;++pit) {
    idtemp=(**pit).id();
    if(abs(idtemp)==16) idnu=idtemp; 
    else                idother.push_back(idtemp);
  }
  if((idnu==ParticleID::nu_tau    && idin==ParticleID::tauminus)||
     (idnu==ParticleID::nu_taubar && idin==ParticleID::tauplus )) {
    allowed=current_->accept(idother);
  }
  return allowed;
}


int TauDecayer2::modeNumber(bool & cc,tcPDPtr parent, const tPDVector & children) const {
  int imode(-1);
  tPDVector::const_iterator pit  = children.begin();
  tPDVector::const_iterator pend = children.end();
  int idtemp;vector<int> idother;
  for( ; pit!=pend;++pit) {
    idtemp=(**pit).id();
    if(abs(idtemp)!=16) idother.push_back(idtemp);
  }
  unsigned int itemp=current_->decayMode(idother);
  for(unsigned int ix=0;ix<modeMap_.size();++ix) {
    if(modeMap_[ix]==itemp) imode=ix;
  }
  // perform the decay
  cc=parent->id()==ParticleID::tauplus;
  return imode;
}


void TauDecayer2::persistentOutput(PersistentOStream & os) const {
  os << modeMap_ << current_ << wgtLoc_ 
     << wgtMax_ << weights_ << polOpt_ << tauMpol_ << tauPpol_;
}

void TauDecayer2::persistentInput(PersistentIStream & is, int) {
  is >> modeMap_ >> current_ >> wgtLoc_ 
     >> wgtMax_ >> weights_ >> polOpt_ >> tauMpol_ >> tauPpol_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<TauDecayer2,DecayIntegrator2>
describeHerwigTauDecayer2("Herwig::TauDecayer2", "HwTauDecay.so");

void TauDecayer2::Init() {

  static ClassDocumentation<TauDecayer2> documentation
    ("The TauDecayer2 class is designed to use a weak current"
     " to perform the decay of the tau.");

  static Reference<TauDecayer2,WeakCurrent> interfaceWeakCurrent
    ("WeakCurrent",
     "The reference for the decay current to be used.",
     &TauDecayer2::current_, false, false, true, false, false);

  static ParVector<TauDecayer2,int> interfaceWeightLocation
    ("WeightLocation",
     "The locations of the weights for a given channel in the vector",
     &TauDecayer2::wgtLoc_,
     0, 0, 0, 0, 10000, false, false, true);

  static ParVector<TauDecayer2,double> interfaceWeightMax
    ("MaximumWeight",
     "The maximum weight for a given channel.",
     &TauDecayer2::wgtMax_,
     0, 0, 0, 0., 100., false, false, true);

  static ParVector<TauDecayer2,double> interfaceWeights
    ("Weights",
     "The weights for the integration.",
     &TauDecayer2::weights_,
     0, 0, 0, 0., 1., false, false, true);

  static Switch<TauDecayer2,bool> interfacePolarizationOption
    ("PolarizationOption",
     "Option of forcing the polarization of the tau leptons, N.B. you"
     " should only use this option for making distributions for"
     " comparision if you really know what you are doing.",
     &TauDecayer2::polOpt_, false, false, false);
  static SwitchOption interfacePolarizationOptionDefault
    (interfacePolarizationOption,
     "Default",
     "Don't force the polarization use the full spin density matrices"
     " to get the right answer",
     false);
  static SwitchOption interfacePolarizationOptionForce
    (interfacePolarizationOption,
     "Force",
     "Force the polarizations",
     true);

  static Parameter<TauDecayer2,double> interfaceTauMinusPolarization
    ("TauMinusPolarization",
     "The polarization of the tau-, left=-1, right=+1 if this is forced.",
     &TauDecayer2::tauMpol_, 0.0, -1.0, 1.0,
     false, false, Interface::limited);


  static Parameter<TauDecayer2,double> interfaceTauPlusPolarization
    ("TauPlusPolarization",
     "The polarization of the tau+, left=-1, right=+1 if this is forced.",
     &TauDecayer2::tauPpol_, 0.0, -1.0, 1.0,
     false, false, Interface::limited);

}

void TauDecayer2::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  assert(false);
}
// combine the currents to give the matrix element
double TauDecayer2::me2(const int ichan, const Particle & part,
				    const tPDVector & ,
				    const vector<Lorentz5Momentum> & momenta,
				    MEOption meopt) const {
  assert(false);
  // // map the mode to those in the current
  // int mode(modeMap_[imode()]);
  // // get the particles for the hadronic current
  // ParticleVector hadpart(decay.begin()+1,decay.end());
  // Energy q;
  // // extract info on the decaying particle
  // if(meopt==Initialize) {
  //   // spin density matrix for the decaying particle
  //   rho_ = RhoDMatrix(PDT::Spin1Half);
  //   if(inpart.id()==ParticleID::tauminus)
  //     SpinorWaveFunction   ::calculateWaveFunctions(inSpin_,rho_,
  // 						    const_ptr_cast<tPPtr>(&inpart),
  // 						    incoming);
  //   else
  //     SpinorBarWaveFunction::calculateWaveFunctions(inBar_ ,rho_,
  // 						    const_ptr_cast<tPPtr>(&inpart),
  // 						    incoming);
  //   // fix rho if no correlations
  //   fixRho(rho_);
  //   if(polOpt_) {
  //     rho_(0,1) = rho_(1,0) = 0.;
  //     if(inpart.id()==ParticleID::tauminus) {
  // 	rho_(0,0) = 0.5*(1.-tauMpol_);
  // 	rho_(1,1) = 0.5*(1.+tauMpol_);
  //     }
  //     else {
  // 	rho_(0,0) = 0.5*(1.+tauPpol_);
  // 	rho_(1,1) = 0.5*(1.-tauPpol_);
  //     }
  //   }
  //   // work out the mapping for the hadron vector
  //   constants_ = vector<unsigned int>(decay.size()+1);
  //   iSpin_     = vector<PDT::Spin   >(decay.size());
  //   int itemp(1);
  //   unsigned int ix(decay.size());
  //   do {
  //     --ix;
  //     iSpin_[ix]     = decay[ix]->data().iSpin();
  //     itemp         *= iSpin_[ix];
  //     constants_[ix] = itemp;
  //   }
  //   while(ix>0);
  //   constants_[decay.size()] = 1;
  //   constants_[0           ] = constants_[1];
  // }
  // if(!ME())
  //   ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1Half,iSpin_)));  
  // // connect the spininfo up if needed
  // if(meopt==Terminate) {
  //   if(inpart.id()==ParticleID::tauminus) {
  //     SpinorWaveFunction   ::
  // 	constructSpinInfo(inSpin_,const_ptr_cast<tPPtr>(&inpart),incoming,true);
  //     SpinorBarWaveFunction::
  // 	constructSpinInfo(inBar_,decay[0],outgoing,true);
  //   }
  //   else {
  //     SpinorBarWaveFunction::
  // 	constructSpinInfo(inBar_ ,const_ptr_cast<tPPtr>(&inpart),incoming,true);
  //     SpinorWaveFunction::
  // 	constructSpinInfo(inSpin_,decay[0],outgoing,true);
  //   }
  //   current_->current(mode,ichan,q,hadpart,meopt);
  //   return 0.;
  // }
  // // calculate the spinors for the decay products
  // if(inpart.id()==ParticleID::tauminus)
  //   SpinorBarWaveFunction::calculateWaveFunctions(inBar_ ,decay[0],outgoing);
  // else   
  //   SpinorWaveFunction   ::calculateWaveFunctions(inSpin_,decay[0],outgoing);
  // // calculate the hadron current
  // vector<LorentzPolarizationVectorE> 
  //   hadron(current_->current(mode,ichan,q,hadpart,meopt));
  // // prefactor
  // double pre = sqr(pow(inpart.mass()/q,int(hadpart.size()-2)));
  // // calculate the lepton current
  // LorentzPolarizationVectorE lepton[2][2];
  // for(unsigned ix=0;ix<2;++ix) {
  //   for(unsigned iy=0;iy<2;++iy) {
  //     if(inpart.id()==15) 
  // 	lepton[ix][iy]=2.*inSpin_[ix].leftCurrent(inBar_[iy]); 
  //     else                
  // 	lepton[iy][ix]=2.*inSpin_[ix].leftCurrent(inBar_[iy]); 
  //   }
  // }
  // // compute the matrix element
  // vector<unsigned int> ihel(decay.size()+1);
  // for(unsigned int hhel=0;hhel<hadron.size();++hhel) {
  //   // map the index for the hadrons to a helicity state
  //   for(unsigned int ix=decay.size();ix>1;--ix) {
  //     ihel[ix]=(hhel%constants_[ix-1])/constants_[ix];
  //   }
  //   // loop over the helicities of the tau and neutrino and set up the matrix 
  //   // element
  //   for(ihel[1]=0;ihel[1]<2;++ihel[1]){
  //     for(ihel[0]=0;ihel[0]<2;++ihel[0]) {
  // 	(*ME())(ihel)= lepton[ihel[0]][ihel[1]].dot(hadron[hhel])*
  // 	  SM().fermiConstant();
  //     }
  //   }
  // }
  // // multiply by the CKM element
  // int iq,ia;
  // current_->decayModeInfo(mode,iq,ia);
  // double ckm(1.);
  // if(iq<=6) {
  //   if(iq%2==0) ckm = SM().CKM(iq/2-1,(abs(ia)-1)/2);
  //   else        ckm = SM().CKM(abs(ia)/2-1,(iq-1)/2);
  // }
  // return 0.5*pre*ckm*(ME()->contract(rho_)).real();
}
  
// output the setup information for the particle database
void TauDecayer2::dataBaseOutput(ofstream & output,bool header) const {
  unsigned int ix;
  if(header) output << "update decayers set parameters=\"";
  DecayIntegrator2::dataBaseOutput(output,false);
  for(ix=0;ix<wgtLoc_.size();++ix) {
    output << "insert " << name() << ":WeightLocation " << ix << " " 
	   << wgtLoc_[ix] << "\n";
  }
  for(ix=0;ix<wgtMax_.size();++ix) {
    output << "insert " << name() << ":MaximumWeight "  << ix << " " 
	   << wgtMax_[ix] << "\n";
  }
  for(ix=0;ix<weights_.size();++ix) {
    output << "insert " << name() << ":Weights "        << ix << " " 
	   << weights_[ix] << "\n";
  }
  current_->dataBaseOutput(output,false,true);
  output << "newdef " << name() << ":WeakCurrent " << current_->name() << " \n";
  output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";\n";
}
