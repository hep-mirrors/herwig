// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TauDecayer class.
//
//  Author: Peter Richardson
//

#include "TauDecayer.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Helicity/VectorSpinInfo.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig++/Decay/DecayVertex.h"
#include "ThePEG/Helicity/FermionSpinInfo.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

namespace Herwig {
using namespace ThePEG;
using Helicity::SpinorWaveFunction;
using Helicity::SpinorBarWaveFunction;
using ThePEG::Helicity::LorentzPolarizationVector;
using ThePEG::Helicity::LorentzSpinor;
using ThePEG::Helicity::LorentzSpinorBar;
using ThePEG::Helicity::RhoDMatrix;
using Helicity::Direction;
using Helicity::incoming;
using Helicity::outgoing;

void TauDecayer::doinit() throw(InitException) {
  DecayIntegrator::doinit();
  // make sure the current got initialised
  _current->init();
  // set up the phase-space channels
  DecayPhaseSpaceModePtr mode;
  DecayPhaseSpaceChannelPtr channel;
  PDVector extpart,ptemp;
  extpart.push_back(getParticleData(ParticleID::tauminus));
  extpart.push_back(getParticleData(ParticleID::nu_tau));
  Energy mtau(extpart[0]->mass());
  double maxweight;
  vector<double> channelwgts;
  int iq(0),ia(0);
  _modemap.clear();
  unsigned int ix,iy;
  bool done;
  vector<double>::iterator start,end;
  for(ix=0;ix<_current->numberOfModes();++ix) {
    // get the external particles for this mode
    extpart.resize(2);
    ptemp=_current->particles(-3,ix,iq,ia);
    for(iy=0;iy<ptemp.size();++iy) extpart.push_back(ptemp[iy]);
    // create the mode
    mode=new_ptr(DecayPhaseSpaceMode(extpart,this));
    // create the first piece of the channel
    channel = new_ptr(DecayPhaseSpaceChannel(mode));
    channel->addIntermediate(extpart[0],0,0.0,-1,1);
    done=_current->createMode(-3,ix,mode,2,1,channel,mtau);
    if(done) {
      // the maximum weight and the channel weights
      // the maximum
      maxweight = _wgtmax.size()>numberModes() ? _wgtmax[numberModes()] : 0;
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
      _modemap.push_back(ix);
      // special for the two body modes
      if(extpart.size()==3) {
	channelwgts.clear();
	mode=new_ptr(DecayPhaseSpaceMode(extpart,this));
      }
      addMode(mode,maxweight,channelwgts);
    }
  }
  _current->reset();
  _current->touch();
  _current->update();
}
  
bool TauDecayer::accept(const DecayMode & dm) const {
  bool allowed(false);
  // find the neutrino 
  int idnu(0),idtemp,idin(dm.parent()->id());
  vector<int> idother;
  ParticleMSet::const_iterator pit  = dm.products().begin();
  ParticleMSet::const_iterator pend = dm.products().end();
  for( ; pit!=pend;++pit) {
    idtemp=(**pit).id();
    if(abs(idtemp)==16) idnu=idtemp; 
    else                idother.push_back(idtemp);
  }
  if((idnu==ParticleID::nu_tau    && idin==ParticleID::tauminus)||
     (idnu==ParticleID::nu_taubar && idin==ParticleID::tauplus )) {
    allowed=_current->accept(idother);
  }
  return allowed;
}


int TauDecayer::modeNumber(bool & cc,const DecayMode & dm) const {
  int imode(-1);
  ParticleMSet::const_iterator pit = dm.products().begin();
  ParticleMSet::const_iterator pend = dm.products().end();
  int idtemp;vector<int> idother;
  for( ; pit!=pend;++pit) {
    idtemp=(**pit).id();
    if(abs(idtemp)!=16) idother.push_back(idtemp);
  }
  unsigned int itemp=_current->decayMode(idother);
  for(unsigned int ix=0;ix<_modemap.size();++ix) {
    if(_modemap[ix]==itemp) imode=ix;
  }
  // perform the decay
  cc=dm.parent()->id()==ParticleID::tauplus;
  return imode;
}


void TauDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(_gf,1/GeV2) << _modemap << _current << _wgtloc << _wgtmax << _weights;
}

void TauDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_gf,1/GeV2) >> _modemap >> _current >> _wgtloc >> _wgtmax >> _weights;
}

ClassDescription<TauDecayer> TauDecayer::initTauDecayer;
// Definition of the static class description member.

void TauDecayer::Init() {

  static ClassDocumentation<TauDecayer> documentation
    ("The TauDecayer class is designed to use a weak current"
     " to perform the decay of the tau.");

  static Parameter<TauDecayer,InvEnergy2> interfaceGFermi
    ("GFermi",
     "The Fermi coupling constant",
     &TauDecayer::_gf, 1./GeV2, 1.16637E-5/GeV2, 0./GeV2, 1.e-3/GeV2,
     false, false, false);

  static Reference<TauDecayer,WeakDecayCurrent> interfaceWeakCurrent
    ("WeakCurrent",
     "The reference for the decay current to be used.",
     &TauDecayer::_current, false, false, true, false, false);

  static ParVector<TauDecayer,int> interfaceWeightLocation
    ("WeightLocation",
     "The locations of the weights for a given channel in the vector",
     &TauDecayer::_wgtloc,
     0, 0, 0, 0, 10000, false, false, true);

  static ParVector<TauDecayer,double> interfaceWeightMax
    ("MaximumWeight",
     "The maximum weight for a given channel.",
     &TauDecayer::_wgtmax,
     0, 0, 0, 0., 100., false, false, true);

  static ParVector<TauDecayer,double> interfaceWeights
    ("Weights",
     "The weights for the integration.",
     &TauDecayer::_weights,
     0, 0, 0, 0., 1., false, false, true);
}

// combine the currents to give the matrix element
double TauDecayer::me2(bool vertex, const int ichan,const Particle & inpart,
		       const ParticleVector & decay) const {
  // spin density matrix for the decaying particle
  RhoDMatrix temp(PDT::Spin1Half);temp.average();
  // storage for the wavefunctions of the tau and neutrino
  vector<LorentzSpinor<SqrtEnergy> > wave;
  vector<LorentzSpinorBar<SqrtEnergy> > wavebar;
  // calculate or extract the wavefunctions of the neutrino
  if(inpart.id()==ParticleID::tauminus) {
    SpinorWaveFunction   (wave   ,temp,const_ptr_cast<tPPtr>(&inpart),
			  incoming,true,vertex);
    SpinorBarWaveFunction(wavebar,decay[0],outgoing,true,vertex);
  }
  else {
    SpinorBarWaveFunction(wavebar,temp,const_ptr_cast<tPPtr>(&inpart),
			  incoming,true,vertex);
    SpinorWaveFunction   (wave   ,decay[0],outgoing,true,vertex);
  }
  // map the mode to those in the current
  int mode(_modemap[imode()]);
  // get the particles for the hadronic current
  ParticleVector::const_iterator start(decay.begin()+1),end(decay.end());
  ParticleVector hadpart(start,end);
  // calculate the hadron current
  Energy q;
  vector<LorentzPolarizationVectorE> 
    hadron(_current->current(vertex,mode,ichan,q,hadpart));
  // prefactor
  double pre(pow(inpart.mass()/q,int(hadpart.size()-2)));pre*=pre;
  // work out the mapping for the hadron vector
  vector<unsigned int> constants(decay.size()+1),ihel(decay.size()+1);
  vector<PDT::Spin> ispin(decay.size());
  int itemp(1);
  unsigned int hhel,ix(decay.size()),iy;
  do {
    --ix;
    ispin[ix]=decay[ix]->data().iSpin();
    itemp*=ispin[ix];constants[ix]=itemp;
  }
  while(ix>0);
  constants[decay.size()]=1;
  constants[0]=constants[1];
  // calculate the lepton current
  LorentzPolarizationVectorE lepton[2][2];
  for(ix=0;ix<2;++ix) {
    for(iy=0;iy<2;++iy) {
      if(inpart.id()==15) 
	lepton[ix][iy]=2.*wave[ix].leftCurrent(wavebar[iy]); 
      else                
	lepton[iy][ix]=2.*wave[ix].leftCurrent(wavebar[iy]); 
    }
  }
  // compute the matrix element
  DecayMatrixElement newME(PDT::Spin1Half,ispin);
  for(hhel=0;hhel<hadron.size();++hhel) {
    // map the index for the hadrons to a helicity state
    for(ix=decay.size();ix>1;--ix) {
      ihel[ix]=(hhel%constants[ix-1])/constants[ix];
    }
    // loop over the helicities of the tau and neutrino and set up the matrix 
    // element
    for(ihel[1]=0;ihel[1]<2;++ihel[1]){
      for(ihel[0]=0;ihel[0]<2;++ihel[0]) {
	newME(ihel)= lepton[ihel[0]][ihel[1]].dot(hadron[hhel])*UnitRemoval::InvE2;
      }
    }
  }
  // store the matrix element
  ME(newME);
  // multiply by the CKM element
  int iq,ia;
  _current->decayModeInfo(mode,iq,ia);
  double ckm(1.);
  if(iq<=6) {
    if(iq%2==0) ckm = SM().CKM(iq/2-1,(abs(ia)-1)/2);
    else        ckm = SM().CKM(abs(ia)/2-1,(iq-1)/2);
  }
  return 0.5*pre*ckm*(newME.contract(temp)).real()*_gf*_gf*UnitRemoval::E4;
}
  
// output the setup information for the particle database
void TauDecayer::dataBaseOutput(ofstream & output,bool header) const {
  unsigned int ix;
  if(header) output << "update decayers set parameters=\"";
  DecayIntegrator::dataBaseOutput(output,false);
  output << "set " << fullName() << ":GFermi "    << _gf*GeV2 << "\n";
  for(ix=0;ix<_wgtloc.size();++ix) {
    output << "insert " << fullName() << ":WeightLocation " << ix << " " 
	   << _wgtloc[ix] << "\n";
  }
  for(ix=0;ix<_wgtmax.size();++ix) {
    output << "insert " << fullName() << ":MaximumWeight "  << ix << " " 
	   << _wgtmax[ix] << "\n";
  }
  for(ix=0;ix<_weights.size();++ix) {
    output << "insert " << fullName() << ":Weights "        << ix << " " 
	   << _weights[ix] << "\n";
  }
  _current->dataBaseOutput(output,false,true);
  output << "set " << fullName() << ":WeakCurrent " << _current->fullName() << " \n";
  output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";\n";
}
}
