// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMHiggsWWDecayer class.
//

#include "SMHiggsWWDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/PDT/ParticleData.h"

using namespace Herwig;
typedef Selector<tDMPtr> DecaySelector;

ClassDescription<SMHiggsWWDecayer> SMHiggsWWDecayer::initSMHiggsWWDecayer;
// Definition of the static class description member.

void SMHiggsWWDecayer::Init() {

  static ClassDocumentation<SMHiggsWWDecayer> documentation
    ("The SMHiggsWWDecayer class performs the decay of the Standard Model Higgs"
     " boson to W+w- and Z0Z0");

  static Parameter<SMHiggsWWDecayer,double> interfaceWMaximum
    ("WMaximum",
     "The maximum weight for H-> W+W- decays",
     &SMHiggsWWDecayer::_wmax, 7.0, 0.0001, 1000.,
     false, false, Interface::limited);

  static Parameter<SMHiggsWWDecayer,double> interfaceZMaximum
    ("ZMaximum",
     "The maximum weight for H-> Z0Z0 decays",
     &SMHiggsWWDecayer::_zmax, 0.4, 0.0001, 1000.,
     false, false, Interface::limited);

  static Switch<SMHiggsWWDecayer,bool> interfaceBreitWigner
    ("BreitWigner",
     "Whether to generate the boson masses using a Breit-Wigner or"
     " power-law distribution.",
     &SMHiggsWWDecayer::_breit, true, false, false);
  static SwitchOption interfaceBreitWignerBreitWigner
    (interfaceBreitWigner,
     "BreitWigner",
     "Use a Breit-Wigner",
     true);
  static SwitchOption interfaceBreitWignerPowerLaw
    (interfaceBreitWigner,
     "PowerLaw",
     "Use a power law",
     false);

  static Parameter<SMHiggsWWDecayer,double> interfacePower
    ("Power",
     "The power to use for the power law",
     &SMHiggsWWDecayer::_power, 0.0, -5.0, 5.0,
     false, false, Interface::limited);

}

SMHiggsWWDecayer::SMHiggsWWDecayer() : _wmax(7.0), _zmax(0.4), 
				       _breit(true),_power(0.)
{}

void SMHiggsWWDecayer::doinit() throw(InitException) {
  DecayIntegrator::doinit();
  // get the vertices from the Standard Model object
  tcHwSMPtr hwsm=dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  if(!hwsm) 
    throw InitException() << "SMHiggsWWDecayer needs the StandardModel class"
			  << " to be either the Herwig++ one or a class inheriting"
			  << " from it";
  _theFFWVertex = hwsm->vertexFFW();
  _theFFZVertex = hwsm->vertexFFZ();
  _theHVVVertex = hwsm->vertexWWH();
  // phase space options
  unsigned int iopt = _breit ? 0 : 1;
  // the W+W- decays
  tPDPtr h0     = getParticleData(ParticleID::h0);
  tPDPtr wplus  = getParticleData(ParticleID::Wplus);
  tPDPtr wminus = getParticleData(ParticleID::Wminus);
  DecaySelector wpDecay =  wplus->decaySelector();
  DecaySelector wmDecay = wminus->decaySelector();
  PDVector extpart(5);
  extpart[0]=h0;
  DecayPhaseSpaceModePtr mode;
  DecayPhaseSpaceChannelPtr newchannel;
  vector<double> wgt(1,1.);
  unsigned int imode=0;
  for(DecaySelector::const_iterator wp=wpDecay.begin();wp!=wpDecay.end();++wp) {
    // extract the decay products of W+
    PDVector prod=(*wp).second->orderedProducts();
    if(prod[0]->id()<prod[1]->id()) swap(prod[0],prod[1]);
    extpart[1]=prod[0];
    extpart[2]=prod[1];
    for(DecaySelector::const_iterator wm=wmDecay.begin();wm!=wmDecay.end();++wm) {
      // extract the decay products of W-
      PDVector prod=(*wm).second->orderedProducts();
      if(prod[0]->id()<prod[1]->id()) swap(prod[0],prod[1]);
      extpart[3]=prod[0];
      extpart[4]=prod[1];
      // create the new mode
      mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
      // create the phase space channel
      newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
      newchannel->addIntermediate(extpart[0],0,0.0,-1,-2);
      newchannel->addIntermediate(wplus     ,iopt,_power, 1, 2);
      newchannel->addIntermediate(wminus    ,iopt,_power, 3, 4);
      mode->addChannel(newchannel);
      addMode(mode,_wmax,wgt);
      // insert mode into selector
      _ratio.push_back(wp->second->brat()*wm->second->brat());
      _wdecays.insert (_ratio.back(),imode);
      ++imode;
    }
  }
  // the Z0Z0 decays
  tPDPtr Z0=getParticleData(ParticleID::Z0);
  DecaySelector Z0Decay = Z0->decaySelector();
  for(DecaySelector::const_iterator z1=Z0Decay.begin();z1!=Z0Decay.end();++z1) {
    // extract the decay products of W+
    PDVector prod=(*z1).second->orderedProducts();
    if(prod[0]->id()<prod[1]->id()) swap(prod[0],prod[1]);
    extpart[1]=prod[0];
    extpart[2]=prod[1];
    for(DecaySelector::const_iterator z2=Z0Decay.begin();z2!=Z0Decay.end();++z2) {
      // extract the decay products of W-
      PDVector prod=(*z2).second->orderedProducts();
      if(prod[0]->id()<prod[1]->id()) swap(prod[0],prod[1]);
      extpart[3]=prod[0];
      extpart[4]=prod[1];
      // create the new mode
      mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
      // create the phase space channel
      newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
      newchannel->addIntermediate(extpart[0],0,0.0,-1,-2);
      newchannel->addIntermediate(Z0        ,iopt,_power, 1, 2);
      newchannel->addIntermediate(Z0        ,iopt,_power, 3, 4);
      mode->addChannel(newchannel);
      addMode(mode,_zmax,wgt);
      // insert mode into selector
      _ratio.push_back(z1->second->brat()*z2->second->brat());
      _zdecays.insert (_ratio.back(),imode);
      ++imode;
    }
  }
}

bool SMHiggsWWDecayer::accept(tcPDPtr parent, const PDVector & children) const {
  // if not two decay products return false
  if(children.size()!=2) return false;
  // if not decaying higgs return false
  if(parent->id()!=ParticleID::h0) return false;
  PDVector::const_iterator pit = children.begin();
  int id1=(**pit).id();
  ++pit;
  int id2=(**pit).id();
  if((id1==-id2&&abs(id1)==ParticleID::Wplus)||
     (id1== id2&&    id1 ==ParticleID::Z0))
    return true;
  else
    return false;
}

void SMHiggsWWDecayer::persistentOutput(PersistentOStream & os) const {
  os << _theFFWVertex << _theFFZVertex << _theHVVVertex 
     << _wdecays << _zdecays << _ratio << _wmax << _zmax << _breit << _power;
}

void SMHiggsWWDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _theFFWVertex >> _theFFZVertex >> _theHVVVertex 
     >> _wdecays >> _zdecays >> _ratio >> _wmax >> _zmax >> _breit >> _power;
}

ParticleVector SMHiggsWWDecayer::decay(const Particle & parent,
				       const PDVector & children) const {
  // select the decay modes of the gauge bosons
  unsigned int imode;
  if(abs(children[0]->id())==ParticleID::Wplus)
    imode=_wdecays.select(UseRandom::rnd());
  else
    imode=_zdecays.select(UseRandom::rnd());
  return generate(true,false,imode,parent);
}

double SMHiggsWWDecayer::me2(bool vertex, const int, const Particle & inpart,
			     const ParticleVector & decay) const {
  RhoDMatrix rhoin;
  // check if the incoming particle has a spin info and if not create it
  ScalarWaveFunction inwave = ScalarWaveFunction(const_ptr_cast<tPPtr>(&inpart),
						 rhoin,incoming,true,vertex);
  // check if Z or W decay
  bool Z0=decay[0]->id()==-decay[1]->id();
  // get the intermediates and vertex
  tcPDPtr inter[2];
  FFVVertexPtr vert;
  if(Z0) {
    inter[0]=getParticleData(ParticleID::Z0);
    inter[1]=inter[0];
    vert=_theFFZVertex;
  }
  else {
    inter[0]=getParticleData(ParticleID::Wplus);
    inter[1]=getParticleData(ParticleID::Wminus);
    vert=_theFFWVertex;
  }
  // construct the spinors for the outgoing particles
  vector<SpinorWaveFunction   > awave1,awave2;
  vector<SpinorBarWaveFunction> fwave1,fwave2;
  SpinorBarWaveFunction(fwave1,decay[0],outgoing,true,vertex);
  SpinorWaveFunction   (awave1,decay[1],outgoing,true,vertex);
  SpinorBarWaveFunction(fwave2,decay[2],outgoing,true,vertex);
  SpinorWaveFunction   (awave2,decay[3],outgoing,true,vertex);
  // compute the matrix element
  DecayMatrixElement newme(PDT::Spin0,PDT::Spin1Half,PDT::Spin1Half,
			   PDT::Spin1Half,PDT::Spin1Half);
  Energy2 scale0(sqr(inpart.mass()));
  Energy2 scale1((decay[0]->momentum()+decay[1]->momentum()).m2());
  Energy2 scale2((decay[2]->momentum()+decay[3]->momentum()).m2());
  // compute the boson currents
  VectorWaveFunction curr1[2][2],curr2[2][2];
  unsigned int ohel1,ohel2,ohel3,ohel4;
  for(ohel1=0;ohel1<2;++ohel1) {
    for(ohel2=0;ohel2<2;++ohel2) {
      curr1[ohel1][ohel2]=vert->evaluate(scale1,1,inter[0],
					 awave1[ohel2],fwave1[ohel1]);
      curr2[ohel1][ohel2]=vert->evaluate(scale2,1,inter[1],
					 awave2[ohel2],fwave2[ohel1]);
    }
  }
  // compute the matrix element
  for(ohel1=0;ohel1<2;++ohel1) {
    for(ohel2=0;ohel2<2;++ohel2) {
      for(ohel3=0;ohel3<2;++ohel3) {
	for(ohel4=0;ohel4<2;++ohel4) {
	  newme(0,ohel1,ohel2,ohel3,ohel4)=
	    _theHVVVertex->evaluate(scale0,curr1[ohel1][ohel2],
				    curr2[ohel3][ohel4],inwave);
	}
      }
    }
  }
  ME(newme);
  double output=(newme.contract(rhoin)).real()*scale0*UnitRemoval::InvE2;
  // set up the colour flows
  if(decay[0]->coloured()) {
    output*=3.;
    decay[0]->antiColourNeighbour(decay[1]);
  }
  if(decay[2]->coloured()) {
    output*=3.;
    decay[2]->antiColourNeighbour(decay[3]);
  }
  // divide out the gauge boson branching ratios
  output/=_ratio[imode()];
  // if Z0 decays identical particle factor
  if(Z0) output*=0.5;
  // return the answer
  return output;
}








