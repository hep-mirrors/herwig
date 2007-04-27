// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GoityRobertsDecayer class.
//

#include "GoityRobertsDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig++/Helicity/EpsFunction.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

using namespace Herwig;
using namespace Herwig::Helicity;
using namespace ThePEG::Helicity;

void GoityRobertsDecayer::persistentOutput(PersistentOStream & os) const {
  os << _current << _GF << _includeDstar << _includeDstarstar << _beta1S << _beta2S 
     << _beta1P << _beta1D << _fpi << _deltaM2S << _deltaM1P << _deltaM1D << _gamma2S 
     << _gamma1P << _gamma1D << _Lambdabar << _g << _alpha1 << _alpha2 << _alpha3 
     << _wgtloc << _wgtmax << _weights;
}

void GoityRobertsDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _current >> _GF >> _includeDstar >> _includeDstarstar >> _beta1S >> _beta2S 
     >> _beta1P >> _beta1D >> _fpi >> _deltaM2S >> _deltaM1P >> _deltaM1D >> _gamma2S 
     >> _gamma1P >> _gamma1D >> _Lambdabar >> _g >> _alpha1 >> _alpha2 >> _alpha3 
     >> _wgtloc >> _wgtmax >> _weights;
}

ClassDescription<GoityRobertsDecayer> GoityRobertsDecayer::initGoityRobertsDecayer;
// Definition of the static class description member.

void GoityRobertsDecayer::Init() {

  static ClassDocumentation<GoityRobertsDecayer> documentation
    ("The GoityRobertsDecayer class implements the model of PRD51 3459 for semi-leptonic"
     " B decays to D(*) pi l nu ");

  static Reference<GoityRobertsDecayer,LeptonNeutrinoCurrent> interfaceCurrent
    ("Current",
     "The current for the leptons produced in the decay.",
     &GoityRobertsDecayer::_current, true, true, true, false, false);

  static Parameter<GoityRobertsDecayer,InvEnergy2> interfaceGFermi
    ("GFermi",
     "The Fermi coupling constant",
     &GoityRobertsDecayer::_GF, 1./GeV2, 1.16639E-5/GeV2,
     0./GeV2, 1.e-4/GeV2, false, false, false);

  static Switch<GoityRobertsDecayer,bool> interfaceIncludeDstar
    ("IncludeDstar",
     "Include the D* resonance in the decay B -> D pi",
     &GoityRobertsDecayer::_includeDstar, false, false, false);
  static SwitchOption interfaceIncludeDstarInclude
    (interfaceIncludeDstar,
     "Include",
     "Include the D*",
     true);
  static SwitchOption interfaceIncludeDstarOmit
    (interfaceIncludeDstar,
     "Omit",
     "Don't include the D*",
     false);

  static Switch<GoityRobertsDecayer,bool> interfaceIncludeDstarstar
    ("IncludeD**",
     "Include the D** resonances in the decays",
     &GoityRobertsDecayer::_includeDstarstar, false, false, false);
  static SwitchOption interfaceIncludeDstarstarInclude
    (interfaceIncludeDstarstar,
     "Include",
     "Include the D**",
     true);
  static SwitchOption interfaceIncludeDstarstarOmit
    (interfaceIncludeDstarstar,
     "Omit",
     "Don't include the D**",
     false);

  static Parameter<GoityRobertsDecayer,Energy> interfaceBeta1S
    ("Beta1S",
     "The beta parameter for the 1S states",
     &GoityRobertsDecayer::_beta1S, GeV, 0.29*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<GoityRobertsDecayer,Energy> interfaceBeta2S
    ("Beta2S",
     "The beta parameter for the 2S states",
     &GoityRobertsDecayer::_beta2S, GeV, 0.29*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<GoityRobertsDecayer,Energy> interfaceBeta1P
    ("Beta1P",
     "The beta parameter for the 1P states",
     &GoityRobertsDecayer::_beta1P, GeV, 0.28*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<GoityRobertsDecayer,Energy> interfaceBeta1D
    ("Beta1D",
     "The beta parameter for the 1D states",
     &GoityRobertsDecayer::_beta1D, GeV, 0.26*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<GoityRobertsDecayer,Energy> interfacefpi
    ("fpi",
     "The pion decay constant",
     &GoityRobertsDecayer::_fpi, MeV, 92.4*MeV, 0.*MeV, 200.*MeV,
     false, false, false); 

  static Parameter<GoityRobertsDecayer,Energy> interfaceDeltaM2S
    ("DeltaM2S",
     "The mass difference for the 2S states.",
     &GoityRobertsDecayer::_deltaM2S, GeV, 0.56*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<GoityRobertsDecayer,Energy> interfaceDeltaM1P
    ("DeltaM1P",
     "The mass difference for the 1P states.",
     &GoityRobertsDecayer::_deltaM1P, GeV, 0.39*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<GoityRobertsDecayer,Energy> interfaceDeltaM1D
    ("DeltaM1D",
     "The mass difference for the 1D states.",
     &GoityRobertsDecayer::_deltaM1D, GeV, 0.71*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<GoityRobertsDecayer,Energy> interfaceGamma2S
    ("Gamma2S",
     "The width for the 2S states",
     &GoityRobertsDecayer::_gamma2S, GeV, 0.191*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<GoityRobertsDecayer,Energy> interfaceGamma1P
    ("Gamma1P",
     "The width for the 1P states",
     &GoityRobertsDecayer::_gamma1P, GeV, 1.040*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<GoityRobertsDecayer,Energy> interfaceGamma1D
    ("Gamma1D",
     "The width for the 1D states",
     &GoityRobertsDecayer::_gamma1D, GeV, 0.405*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<GoityRobertsDecayer,Energy> interfaceLambdabar
    ("Lambdabar",
     "The Lambdabar mass difference parameter",
     &GoityRobertsDecayer::_Lambdabar, GeV, 0.75*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<GoityRobertsDecayer,double> interfaceg
    ("g",
     "The coupling for the 1S states",
     &GoityRobertsDecayer::_g, 0.50, -10.0, 10.0,
     false, false, true);

  static Parameter<GoityRobertsDecayer,double> interfaceAlpha1
    ("Alpha1",
     "The coupling for the excited 1P states",
     &GoityRobertsDecayer::_alpha1, -1.43, -10.0, 10.0,
     false, false, true);

  static Parameter<GoityRobertsDecayer,double> interfaceAlpha2
    ("Alpha2",
     "The coupling for the excited 1D states",
     &GoityRobertsDecayer::_alpha2, -0.14, -10.0, 10.0,
     false, false, true);

  static Parameter<GoityRobertsDecayer,double> interfaceAlpha3
    ("Alpha3",
     "The coupling for the excited 2S states",
     &GoityRobertsDecayer::_alpha3, 0.69, -10.0, 10.0,
     false, false, true);

  static ParVector<GoityRobertsDecayer,int> interfaceWeightLocation
    ("WeightLocation",
     "The locations of the weights for a given channel in the vector",
     &GoityRobertsDecayer::_wgtloc,
     0, 0, 0, 0, 10000, false, false, true);

  static ParVector<GoityRobertsDecayer,double> interfaceWeightMax
    ("MaximumWeight",
     "The maximum weight for a given channel.",
     &GoityRobertsDecayer::_wgtmax,
     0, 0, 0, 0., 100., false, false, true);

  static ParVector<GoityRobertsDecayer,double> interfaceWeights
    ("Weights",
     "The weights for the integration.",
     &GoityRobertsDecayer::_weights,
     0, 0, 0, 0., 1., false, false, true);
}


void GoityRobertsDecayer::doinit() throw(InitException) {
  DecayIntegrator::doinit();
  // make sure the current got initialized
  _current->init();
  // intermediate particles for the various modes
  tPDPtr B01mbar=getParticleData(-513);
  tPDPtr B00pbar=getParticleData(-10511);
  tPDPtr B01pbar=getParticleData(-10513);
  tPDPtr D01m   =getParticleData(423);
  tPDPtr D00p   =getParticleData(10421);
  tPDPtr D01p   =getParticleData(10423);
  tPDPtr Bm1m   =getParticleData(-523);
  tPDPtr Bm0p   =getParticleData(-10521);
  tPDPtr Bm1p   =getParticleData(-10523);
  tPDPtr Dp1m   =getParticleData(413);
  tPDPtr Dp0p   =getParticleData(10411);
  tPDPtr Dp1p   =getParticleData(10413);
  unsigned int ix,iy,iz;
  int Wcharge,iq,ia;
  bool done;
  Energy min;
  double maxweight;
  vector<double> channelwgts;
  vector<double>::iterator start,end;
  // the channels
  vector<DecayPhaseSpaceChannelPtr> channel;
  DecayPhaseSpaceModePtr mode;
  PDVector extpart(3),ptemp;
  // first B- to D+ pi-
  extpart[0]=getParticleData(ParticleID::Bminus);
  extpart[1]=getParticleData(ParticleID::Dplus);
  extpart[2]=getParticleData(ParticleID::piminus);
  min = extpart[0]->massMax()-extpart[1]->massMin()-extpart[2]->massMin();
  Wcharge = extpart[0]->iCharge()-extpart[1]->iCharge()-extpart[2]->iCharge();
  for(iy=0;iy<_current->numberOfModes();++iy) {
    channel.resize(0);extpart.resize(3);
    // get the particles for the current
    _current->decayModeInfo(iy,iq,ia);
    ptemp=_current->particles(Wcharge,iy,iq,ia);
    for(iz=0;iz<ptemp.size();++iz) extpart.push_back(ptemp[iz]);
    // create the mode
    mode=new_ptr(DecayPhaseSpaceMode(extpart,this));
    // create the different possible phase space channels
    channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
    channel.back()->addIntermediate(extpart[0],0,0.0,2,-1);
    channel.back()->addIntermediate(B01mbar,0,0.0,1,-2);
    channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
    channel.back()->addIntermediate(extpart[0],0,0.0,2,-1);
    channel.back()->addIntermediate(B00pbar,0,0.0,1,-2);
    channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
    channel.back()->addIntermediate(extpart[0],0,0.0,2,-1);
    channel.back()->addIntermediate(B01pbar,0,0.0,1,-2);
    // add the D* diagrams if needed
    if(_includeDstar) {
      channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
      channel.back()->addIntermediate(extpart[0],0,0.0,-1,-2);
      channel.back()->addIntermediate(D01m,0,0.0,1,2);
    }
    // add the D** diagrams if needed
    if(_includeDstarstar) {
      channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
      channel.back()->addIntermediate(extpart[0],0,0.0,-1,-2);
      channel.back()->addIntermediate(D00p,0,0.0,1,2);
      channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
      channel.back()->addIntermediate(extpart[0],0,0.0,-1,-2);
      channel.back()->addIntermediate(D01p,0,0.0,1,2);
    }
    for(ix=0;ix<channel.size();++ix) {
      done=_current->createMode(Wcharge,iy,mode,3,1,channel[ix],min);
    }
    if(_wgtmax.size()>numberModes()) maxweight=_wgtmax[numberModes()];
    else                             maxweight=0.;
    if(_wgtloc.size()>numberModes()&&
       _wgtloc[numberModes()]+mode->numberChannels()<=_weights.size()) {
      start=_weights.begin()+_wgtloc[numberModes()];
      end  = start+mode->numberChannels();
      channelwgts=vector<double>(start,end);
    }
    else {
      channelwgts.resize(mode->numberChannels(),1./(mode->numberChannels()));
    }
    addMode(mode,maxweight,channelwgts);
  }
  // then B- to D0 pi0
  extpart[0]=getParticleData(ParticleID::Bminus);
  extpart[1]=getParticleData(ParticleID::D0);
  extpart[2]=getParticleData(ParticleID::pi0);
  min = extpart[0]->massMax()-extpart[1]->massMin()-extpart[2]->massMin();
  Wcharge = extpart[0]->iCharge()-extpart[1]->iCharge()-extpart[2]->iCharge();
  for(iy=0;iy<_current->numberOfModes();++iy) {
    channel.resize(0);extpart.resize(3);
    // get the particles for the current
    _current->decayModeInfo(iy,iq,ia);
    ptemp=_current->particles(Wcharge,iy,iq,ia);
    for(iz=0;iz<ptemp.size();++iz) extpart.push_back(ptemp[iz]);
    // create the mode
    mode=new_ptr(DecayPhaseSpaceMode(extpart,this));
    // create the different possible phase space channels
    channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
    channel.back()->addIntermediate(extpart[0],0,0.0,2,-1);
    channel.back()->addIntermediate(Bm1m,0,0.0,1,-2);
    channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
    channel.back()->addIntermediate(extpart[0],0,0.0,2,-1);
    channel.back()->addIntermediate(Bm0p,0,0.0,1,-2);
    channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
    channel.back()->addIntermediate(extpart[0],0,0.0,2,-1);
    channel.back()->addIntermediate(Bm1p,0,0.0,1,-2);
    // add the D* diagrams if needed
    if(_includeDstar) {
      channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
      channel.back()->addIntermediate(extpart[0],0,0.0,-1,-2);
      channel.back()->addIntermediate(D01m,0,0.0,1,2);
    }
    // add the D** diagrams if needed
    if(_includeDstarstar) {
      channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
      channel.back()->addIntermediate(extpart[0],0,0.0,-1,-2);
      channel.back()->addIntermediate(D00p,0,0.0,1,2);
      channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
      channel.back()->addIntermediate(extpart[0],0,0.0,-1,-2);
      channel.back()->addIntermediate(D01p,0,0.0,1,2);
    }
    for(ix=0;ix<channel.size();++ix) {
      done=_current->createMode(Wcharge,iy,mode,3,1,channel[ix],min);
    }
    if(_wgtmax.size()>numberModes()) maxweight=_wgtmax[numberModes()];
    else                             maxweight=0.;
    if(_wgtloc.size()>numberModes()&&
       _wgtloc[numberModes()]+mode->numberChannels()<=_weights.size()) {
      start=_weights.begin()+_wgtloc[numberModes()];
      end  = start+mode->numberChannels();
      channelwgts=vector<double>(start,end);
    }
    else {
      channelwgts.resize(mode->numberChannels(),1./(mode->numberChannels()));
    }
    addMode(mode,maxweight,channelwgts);
  }
  // then B0bar  to D0 pi+
  extpart[0]=getParticleData(ParticleID::Bbar0);
  extpart[1]=getParticleData(ParticleID::D0);
  extpart[2]=getParticleData(ParticleID::piplus);
  min = extpart[0]->massMax()-extpart[1]->massMin()-extpart[2]->massMin();
  Wcharge = extpart[0]->iCharge()-extpart[1]->iCharge()-extpart[2]->iCharge();
  for(iy=0;iy<_current->numberOfModes();++iy) {
    channel.resize(0);extpart.resize(3);
    // get the particles for the current
    _current->decayModeInfo(iy,iq,ia);
    ptemp=_current->particles(Wcharge,iy,iq,ia);
    for(iz=0;iz<ptemp.size();++iz) extpart.push_back(ptemp[iz]);
    // create the mode
    mode=new_ptr(DecayPhaseSpaceMode(extpart,this));
    // create the different possible phase space channels
    channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
    channel.back()->addIntermediate(extpart[0],0,0.0,2,-1);
    channel.back()->addIntermediate(Bm1m,0,0.0,1,-2);
    channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
    channel.back()->addIntermediate(extpart[0],0,0.0,2,-1);
    channel.back()->addIntermediate(Bm0p,0,0.0,1,-2);
    channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
    channel.back()->addIntermediate(extpart[0],0,0.0,2,-1);
    channel.back()->addIntermediate(Bm1p,0,0.0,1,-2);
    // add the D* diagrams if needed
    if(_includeDstar) {
      channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
      channel.back()->addIntermediate(extpart[0],0,0.0,-1,-2);
      channel.back()->addIntermediate(Dp1m,0,0.0,1,2);
    }
    // add the D** diagrams if needed
    if(_includeDstarstar) {
      channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
      channel.back()->addIntermediate(extpart[0],0,0.0,-1,-2);
      channel.back()->addIntermediate(Dp0p,0,0.0,1,2);
      channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
      channel.back()->addIntermediate(extpart[0],0,0.0,-1,-2);
	  channel.back()->addIntermediate(Dp1p,0,0.0,1,2);
    }
    for(ix=0;ix<channel.size();++ix) {
      done=_current->createMode(Wcharge,iy,mode,3,1,channel[ix],min);
    }
    if(_wgtmax.size()>numberModes()) maxweight=_wgtmax[numberModes()];
    else                             maxweight=0.;
    if(_wgtloc.size()>numberModes()&&
       _wgtloc[numberModes()]+mode->numberChannels()<=_weights.size()) {
      start=_weights.begin()+_wgtloc[numberModes()];
      end  = start+mode->numberChannels();
      channelwgts=vector<double>(start,end);
    }
    else {
      channelwgts.resize(mode->numberChannels(),1./(mode->numberChannels()));
    }
    addMode(mode,maxweight,channelwgts);
  }
  // first B- to D+ pi0
  extpart[0]=getParticleData(ParticleID::Bbar0);
  extpart[1]=getParticleData(ParticleID::Dplus);
  extpart[2]=getParticleData(ParticleID::pi0);
  min = extpart[0]->massMax()-extpart[1]->massMin()-extpart[2]->massMin();
  Wcharge = extpart[0]->iCharge()-extpart[1]->iCharge()-extpart[2]->iCharge();
  for(iy=0;iy<_current->numberOfModes();++iy) {
    channel.resize(0);extpart.resize(3);
    // get the particles for the current
    _current->decayModeInfo(iy,iq,ia);
    ptemp=_current->particles(Wcharge,iy,iq,ia);
    for(iz=0;iz<ptemp.size();++iz) extpart.push_back(ptemp[iz]);
    // create the mode
    mode=new_ptr(DecayPhaseSpaceMode(extpart,this));
    // create the different possible phase space channels
    channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
    channel.back()->addIntermediate(extpart[0],0,0.0,2,-1);
    channel.back()->addIntermediate(B01mbar,0,0.0,1,-2);
    channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
    channel.back()->addIntermediate(extpart[0],0,0.0,2,-1);
    channel.back()->addIntermediate(B00pbar,0,0.0,1,-2);
    channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
    channel.back()->addIntermediate(extpart[0],0,0.0,2,-1);
    channel.back()->addIntermediate(B01pbar,0,0.0,1,-2);
    // add the D* diagrams if needed
    if(_includeDstar) {
      channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
      channel.back()->addIntermediate(extpart[0],0,0.0,-1,-2);
      channel.back()->addIntermediate(Dp1m,0,0.0,1,2);
    }
    // add the D** diagrams if needed
    if(_includeDstarstar) {
      channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
      channel.back()->addIntermediate(extpart[0],0,0.0,-1,-2);
      channel.back()->addIntermediate(Dp0p,0,0.0,1,2);
      channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
      channel.back()->addIntermediate(extpart[0],0,0.0,-1,-2);
      channel.back()->addIntermediate(Dp1p,0,0.0,1,2);
    }
    for(ix=0;ix<channel.size();++ix) {
      done=_current->createMode(Wcharge,iy,mode,3,1,channel[ix],min);
    }
    if(_wgtmax.size()>numberModes()) maxweight=_wgtmax[numberModes()];
    else                             maxweight=0.;
    if(_wgtloc.size()>numberModes()&&
       _wgtloc[numberModes()]+mode->numberChannels()<=_weights.size()) {
      start=_weights.begin()+_wgtloc[numberModes()];
      end  = start+mode->numberChannels();
      channelwgts=vector<double>(start,end);
    }
    else {
      channelwgts.resize(mode->numberChannels(),1./(mode->numberChannels()));
    }
    addMode(mode,maxweight,channelwgts);
  }
  // B- to D*+ pi-
  extpart[0]=getParticleData(ParticleID::Bminus);
  extpart[1]=getParticleData(ParticleID::Dstarplus);
  extpart[2]=getParticleData(ParticleID::piminus);
  min = extpart[0]->massMax()-extpart[1]->massMin()-extpart[2]->massMin();
  Wcharge = extpart[0]->iCharge()-extpart[1]->iCharge()-extpart[2]->iCharge();
  for(iy=0;iy<_current->numberOfModes();++iy) {
    channel.resize(0);extpart.resize(3);
    // get the particles for the current
    _current->decayModeInfo(iy,iq,ia);
    ptemp=_current->particles(Wcharge,iy,iq,ia);
    for(iz=0;iz<ptemp.size();++iz) extpart.push_back(ptemp[iz]);
    // create the mode
    mode=new_ptr(DecayPhaseSpaceMode(extpart,this));
    // create the different possible phase space channels
    channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
    channel.back()->addIntermediate(extpart[0],0,0.0,2,-1);
    channel.back()->addIntermediate(B01mbar,0,0.0,1,-2);
    channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
    channel.back()->addIntermediate(extpart[0],0,0.0,2,-1);
    channel.back()->addIntermediate(B00pbar,0,0.0,1,-2);
    channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
    channel.back()->addIntermediate(extpart[0],0,0.0,2,-1);
    channel.back()->addIntermediate(B01pbar,0,0.0,1,-2);
    // add the D* diagrams if needed
    if(_includeDstar) {
      channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
      channel.back()->addIntermediate(extpart[0],0,0.0,-1,-2);
      channel.back()->addIntermediate(D01m,0,0.0,1,2);
    }
    // add the D** diagrams if needed
    if(_includeDstarstar) {
      channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
      channel.back()->addIntermediate(extpart[0],0,0.0,-1,-2);
      channel.back()->addIntermediate(D00p,0,0.0,1,2);
    }
    for(ix=0;ix<channel.size();++ix) {
      done=_current->createMode(Wcharge,iy,mode,3,1,channel[ix],min);
    }
    if(_wgtmax.size()>numberModes()) maxweight=_wgtmax[numberModes()];
    else                             maxweight=0.;
    if(_wgtloc.size()>numberModes()&&
       _wgtloc[numberModes()]+mode->numberChannels()<=_weights.size()) {
      start=_weights.begin()+_wgtloc[numberModes()];
      end  = start+mode->numberChannels();
      channelwgts=vector<double>(start,end);
    }
    else {
      channelwgts.resize(mode->numberChannels(),1./(mode->numberChannels()));
    }
    addMode(mode,maxweight,channelwgts);
  }
  // then B- to D*0 pi0
  extpart[0]=getParticleData(ParticleID::Bminus);
  extpart[1]=getParticleData(ParticleID::Dstar0);
  extpart[2]=getParticleData(ParticleID::pi0);
  min = extpart[0]->massMax()-extpart[1]->massMin()-extpart[2]->massMin();
  Wcharge = extpart[0]->iCharge()-extpart[1]->iCharge()-extpart[2]->iCharge();
  for(iy=0;iy<_current->numberOfModes();++iy)  {
    channel.resize(0);extpart.resize(3);
    // get the particles for the current
    _current->decayModeInfo(iy,iq,ia);
    ptemp=_current->particles(Wcharge,iy,iq,ia);
    for(iz=0;iz<ptemp.size();++iz) extpart.push_back(ptemp[iz]);
    // create the mode
    mode=new_ptr(DecayPhaseSpaceMode(extpart,this));
    // create the different possible phase space channels
    channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
    channel.back()->addIntermediate(extpart[0],0,0.0,2,-1);
    channel.back()->addIntermediate(Bm1m,0,0.0,1,-2);
    channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
    channel.back()->addIntermediate(extpart[0],0,0.0,2,-1);
    channel.back()->addIntermediate(Bm0p,0,0.0,1,-2);
    channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
    channel.back()->addIntermediate(extpart[0],0,0.0,2,-1);
    channel.back()->addIntermediate(Bm1p,0,0.0,1,-2);
    // add the D* diagrams if needed
    if(_includeDstar) {
      channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
      channel.back()->addIntermediate(extpart[0],0,0.0,-1,-2);
      channel.back()->addIntermediate(D01m,0,0.0,1,2);
    }
    // add the D** diagrams if needed
    if(_includeDstarstar) {
      channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
      channel.back()->addIntermediate(extpart[0],0,0.0,-1,-2);
      channel.back()->addIntermediate(D00p,0,0.0,1,2);
    }
    for(ix=0;ix<channel.size();++ix) {
      done=_current->createMode(Wcharge,iy,mode,3,1,channel[ix],min);
    }
    if(_wgtmax.size()>numberModes()) maxweight=_wgtmax[numberModes()];
    else                             maxweight=0.;
    if(_wgtloc.size()>numberModes()&&
       _wgtloc[numberModes()]+mode->numberChannels()<=_weights.size()) {
      start=_weights.begin()+_wgtloc[numberModes()];
      end  = start+mode->numberChannels();
      channelwgts=vector<double>(start,end);
    }
    else {
      channelwgts.resize(mode->numberChannels(),1./(mode->numberChannels()));
    }
    addMode(mode,maxweight,channelwgts);
  }
  // first B0bar  to D*0 pi+
  extpart[0]=getParticleData(ParticleID::Bbar0);
  extpart[1]=getParticleData(ParticleID::Dstar0);
  extpart[2]=getParticleData(ParticleID::piplus);
  min = extpart[0]->massMax()-extpart[1]->massMin()-extpart[2]->massMin();
  Wcharge = extpart[0]->iCharge()-extpart[1]->iCharge()-extpart[2]->iCharge();
  for(iy=0;iy<_current->numberOfModes();++iy) {
    channel.resize(0);extpart.resize(3);
    // get the particles for the current
    _current->decayModeInfo(iy,iq,ia);
    ptemp=_current->particles(Wcharge,iy,iq,ia);
    for(iz=0;iz<ptemp.size();++iz) extpart.push_back(ptemp[iz]);
    // create the mode
    mode=new_ptr(DecayPhaseSpaceMode(extpart,this));
    // create the different possible phase space channels
    channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
    channel.back()->addIntermediate(extpart[0],0,0.0,2,-1);
    channel.back()->addIntermediate(Bm1m,0,0.0,1,-2);
    channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
    channel.back()->addIntermediate(extpart[0],0,0.0,2,-1);
    channel.back()->addIntermediate(Bm0p,0,0.0,1,-2);
    channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
    channel.back()->addIntermediate(extpart[0],0,0.0,2,-1);
    channel.back()->addIntermediate(Bm1p,0,0.0,1,-2);
    // add the D* diagrams if needed
    if(_includeDstar) {
      channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
      channel.back()->addIntermediate(extpart[0],0,0.0,-1,-2);
      channel.back()->addIntermediate(Dp1m,0,0.0,1,2);
    }
    // add the D** diagrams if needed
    if(_includeDstarstar) {
      channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
      channel.back()->addIntermediate(extpart[0],0,0.0,-1,-2);
      channel.back()->addIntermediate(Dp0p,0,0.0,1,2);
    }
    for(ix=0;ix<channel.size();++ix) {
      done=_current->createMode(Wcharge,iy,mode,3,1,channel[ix],min);
    }
    if(_wgtmax.size()>numberModes()) maxweight=_wgtmax[numberModes()];
    else                             maxweight=0.;
    if(_wgtloc.size()>numberModes()&&
       _wgtloc[numberModes()]+mode->numberChannels()<=_weights.size()) {
      start=_weights.begin()+_wgtloc[numberModes()];
      end  = start+mode->numberChannels();
      channelwgts=vector<double>(start,end);
    }
    else {
      channelwgts.resize(mode->numberChannels(),1./(mode->numberChannels()));
    }
    addMode(mode,maxweight,channelwgts);
  }
  // first B- to D*+ pi0
  extpart[0]=getParticleData(ParticleID::Bbar0);
  extpart[1]=getParticleData(ParticleID::Dstarplus);
  extpart[2]=getParticleData(ParticleID::pi0);
  min = extpart[0]->massMax()-extpart[1]->massMin()-extpart[2]->massMin();
  Wcharge = extpart[0]->iCharge()-extpart[1]->iCharge()-extpart[2]->iCharge();
  for(iy=0;iy<_current->numberOfModes();++iy) {
    channel.resize(0);extpart.resize(3);
    // get the particles for the current
    _current->decayModeInfo(iy,iq,ia);
    ptemp=_current->particles(Wcharge,iy,iq,ia);
    for(iz=0;iz<ptemp.size();++iz) extpart.push_back(ptemp[iz]);
    // create the mode
    mode=new_ptr(DecayPhaseSpaceMode(extpart,this));
    // create the different possible phase space channels
    channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
    channel.back()->addIntermediate(extpart[0],0,0.0,2,-1);
    channel.back()->addIntermediate(B01mbar,0,0.0,1,-2);
    channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
    channel.back()->addIntermediate(extpart[0],0,0.0,2,-1);
    channel.back()->addIntermediate(B00pbar,0,0.0,1,-2);
    channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
    channel.back()->addIntermediate(extpart[0],0,0.0,2,-1);
    channel.back()->addIntermediate(B01pbar,0,0.0,1,-2);
    // add the D* diagrams if needed
    if(_includeDstar) {
      channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
      channel.back()->addIntermediate(extpart[0],0,0.0,-1,-2);
      channel.back()->addIntermediate(Dp1m,0,0.0,1,2);
    }
    // add the D** diagrams if needed
    if(_includeDstarstar) {
      channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
      channel.back()->addIntermediate(extpart[0],0,0.0,-1,-2);
      channel.back()->addIntermediate(Dp0p,0,0.0,1,2);
    }
    for(ix=0;ix<channel.size();++ix) {
      done=_current->createMode(Wcharge,iy,mode,3,1,channel[ix],min);
    }
    if(_wgtmax.size()>numberModes()) maxweight=_wgtmax[numberModes()];
    else                             maxweight=0.;
    if(_wgtloc.size()>numberModes()&&
       _wgtloc[numberModes()]+mode->numberChannels()<=_weights.size()) {
      start=_weights.begin()+_wgtloc[numberModes()];
      end  = start+mode->numberChannels();
      channelwgts=vector<double>(start,end);
    }
    else {
      channelwgts.resize(mode->numberChannels(),1./(mode->numberChannels()));
    }
    addMode(mode,maxweight,channelwgts);
  }
  generateIntermediates(_includeDstar||_includeDstarstar);
}

void GoityRobertsDecayer::dataBaseOutput(ofstream & output, bool header) const {
  unsigned int ix;
  if(header){output << "update decayers set parameters=\"";}
  DecayIntegrator::dataBaseOutput(output,false);
  _current->dataBaseOutput(output,false,true);
  output << "set " << fullName() << ":Current " << _current->fullName() << " \n";
  output << "set " << fullName() << ":IncludeDstar " << _includeDstar << " \n";
  output << "set " << fullName() << ":IncludeD** " << _includeDstarstar << " \n";
  output << "set " << fullName() << ":Beta1S " << _beta1S/GeV << " \n";
  output << "set " << fullName() << ":Beta2S " << _beta2S/GeV << " \n";
  output << "set " << fullName() << ":Beta1P " << _beta1P/GeV << " \n";
  output << "set " << fullName() << ":Beta1D " << _beta1D/GeV << " \n";
  output << "set " << fullName() << ":fpi " << _fpi/MeV << " \n";
  output << "set " << fullName() << ":DeltaM2S " << _deltaM2S/GeV << " \n";
  output << "set " << fullName() << ":DeltaM1P " << _deltaM1P/GeV << " \n";
  output << "set " << fullName() << ":DeltaM1D " << _deltaM1D/GeV << " \n";
  output << "set " << fullName() << ":Gamma2S " << _gamma2S/GeV << " \n";
  output << "set " << fullName() << ":Gamma1P " << _gamma1P/GeV << " \n";
  output << "set " << fullName() << ":Gamma1D " << _gamma1D/GeV << " \n";
  output << "set " << fullName() << ":Lambdabar " << _Lambdabar/GeV << " \n";
  output << "set " << fullName() << ":g " << _g << " \n";
  output << "set " << fullName() << ":Alpha1 " << _alpha1 << " \n";
  output << "set " << fullName() << ":Alpha2 " << _alpha2 << " \n";
  output << "set " << fullName() << ":Alpha3 " << _alpha3 << " \n";
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
  output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}

bool GoityRobertsDecayer::accept(const DecayMode & dm) const {
  bool allowed(false);
  int id0(dm.parent()->id()),idtemp,idD(0),idp(0);
  vector<int> idother;
  // check number of decay products
  if(dm.products().size()!=4){return false;}
  ParticleMSet::const_iterator pit(dm.products().begin()),pend(dm.products().end());
  for( ; pit!=pend;++pit) {
    idtemp=(**pit).id();
    if(abs(idtemp)<=16){idother.push_back(idtemp);}
    else if(abs(idtemp/100)==4){idD=idtemp;}
    else {idp=idtemp;}
  }
  if(idother.size()!=2||idD==0||idp==0){return false;}
  // check the mesons
  if(id0==ParticleID::Bbar0) {
    if((idD==ParticleID::D0        &&idp==ParticleID::piplus )||
       (idD==ParticleID::Dplus     &&idp==ParticleID::pi0    )||
       (idD==ParticleID::Dstar0    &&idp==ParticleID::piplus )||
       (idD==ParticleID::Dstarplus &&idp==ParticleID::pi0    )) allowed=true;
  }
  else if(id0==ParticleID::B0) {
    if((idD==ParticleID::Dbar0     &&idp==ParticleID::piminus)||
       (idD==ParticleID::Dminus    &&idp==ParticleID::pi0   )||
       (idD==ParticleID::Dstarbar0 &&idp==ParticleID::piminus)||
       (idD==ParticleID::Dstarminus&&idp==ParticleID::pi0   )) allowed=true;
  }
  else if(id0==ParticleID::Bminus) {
    if((idD==ParticleID::D0        &&idp==ParticleID::pi0    )||
       (idD==ParticleID::Dplus     &&idp==ParticleID::piminus)||
       (idD==ParticleID::Dstar0    &&idp==ParticleID::pi0    )||
       (idD==ParticleID::Dstarplus &&idp==ParticleID::piminus)) allowed=true;
  }
  else if(id0==ParticleID::Bplus) {
    if((idD==ParticleID::Dbar0     &&idp==ParticleID::pi0    )||
       (idD==ParticleID::Dminus    &&idp==ParticleID::piplus )||
       (idD==ParticleID::Dstarbar0 &&idp==ParticleID::pi0    )||
       (idD==ParticleID::Dstarminus&&idp==ParticleID::piplus )) allowed=true;
  }
  if(!allowed){return false;}
  // and the leptons
  return _current->accept(idother);
}

int  GoityRobertsDecayer::modeNumber(bool & cc,const DecayMode & dm) const {
  // find the ids of the particles for the decay current
  ParticleMSet::const_iterator pit = dm.products().begin();
  ParticleMSet::const_iterator pend = dm.products().end();
  int id0(dm.parent()->id()),idD(0),idtemp;
  vector<int> idother;
  // find the pdg code for the D meson
  for( ; pit!=pend;++pit) {
    idtemp=(**pit).id();
    if(abs(idtemp/100)==4)   idD=idtemp;
    else if(abs(idtemp)<=16) idother.push_back(idtemp);
  }
  if(idD==0) return -1;
  int imode(-1);
  if(     abs(id0)==ParticleID::Bplus&&abs(idD)==ParticleID::Dplus    ) imode=0;
  else if(abs(id0)==ParticleID::Bplus&&abs(idD)==ParticleID::D0       ) imode=1;
  else if(abs(id0)==ParticleID::B0   &&abs(idD)==ParticleID::D0       ) imode=2;
  else if(abs(id0)==ParticleID::B0   &&abs(idD)==ParticleID::Dplus    ) imode=3;
  else if(abs(id0)==ParticleID::Bplus&&abs(idD)==ParticleID::Dstarplus) imode=4;
  else if(abs(id0)==ParticleID::Bplus&&abs(idD)==ParticleID::Dstar0   ) imode=5;
  else if(abs(id0)==ParticleID::B0   &&abs(idD)==ParticleID::Dstar0   ) imode=6;
  else if(abs(id0)==ParticleID::B0   &&abs(idD)==ParticleID::Dstarplus) imode=7;
  if(imode<0) return imode;
  imode = 3*imode+_current->decayMode(idother);
  cc = id0>0;
  return imode;
}

double GoityRobertsDecayer::me2(bool vertex, const int ichan,const Particle & inpart,
				const ParticleVector & decay) const {
  // spin info for the decaying particle
  ScalarWaveFunction(const_ptr_cast<tPPtr>(&inpart),incoming,true,vertex);
  // spin info for the outgoing pion
  ScalarWaveFunction(decay[1],outgoing,true,vertex);
  // calculate some common variables
  double omega(inpart.momentum()*decay[0]->momentum()/inpart.mass()/decay[0]->mass());
  Energy mb(inpart.mass()),md(decay[0]->mass());
  Complex ii(0.,1.);
  // dot products we will need
  Energy 
     dotv(inpart.momentum()*decay[1]->momentum()/inpart.mass()),
    dotvp(decay[0]->momentum()*decay[1]->momentum()/decay[0]->mass());
  // delta M parameters
  complex<Energy> dmt1(_deltaM1P-0.5*ii*_gamma1P);
  complex<Energy> dmt2(_deltaM1D-0.5*ii*_gamma1D);
  complex<Energy> dmt3(_deltaM2S-0.5*ii*_gamma2S);
  // calculate the mass splittings for the lightest multiplet
  tcPDPtr Rstar,R;
  // calculate the mass differences
  if(abs(inpart.id())==ParticleID::B0||abs(inpart.id())==ParticleID::Bstar0) {
    R     = getParticleData(ParticleID::B0);
    Rstar = getParticleData(ParticleID::Bstar0);
  }
  else {
    R     = getParticleData(ParticleID::Bplus);
    Rstar = getParticleData(ParticleID::Bstarplus);
  }
  complex<Energy> dmb=Rstar->mass()-R->mass()-ii*0.5*Rstar->width();
  if(abs(decay[0]->id())==ParticleID::Dplus&&
     decay[1]->id()==ParticleID::pi0) {
    R     = getParticleData(ParticleID::Dplus);
    Rstar = getParticleData(ParticleID::Dstarplus);
  }
  else if(abs(decay[0]->id())==ParticleID::D0&&
	  abs(decay[1]->id())==ParticleID::piplus) {
    R     = getParticleData(ParticleID::D0);
    Rstar = getParticleData(ParticleID::Dstarplus);
  }
  else if(abs(decay[0]->id())==ParticleID::Dplus&&
	   abs(decay[1]->id())==ParticleID::piplus) {
    R     = getParticleData(ParticleID::Dplus);
    Rstar = getParticleData(ParticleID::Dstar0);
  }
  else if(abs(decay[0]->id())==ParticleID::D0&&
	  abs(decay[1]->id())==ParticleID::pi0) {
    R     = getParticleData(ParticleID::D0);
    Rstar = getParticleData(ParticleID::Dstar0);
  }
  complex<Energy> dmd=Rstar->mass()-R->mass()-ii*0.5*Rstar->width();
  // get the IW form factors
  double xi,xi1,rho1,rho2;
  calculateFormFactors(omega,xi,xi1,rho1,rho2);
  // hadronic current
  vector<LorentzPolarizationVector> hadron;
  // calculate the matrix element
  // for D pi
  if(decay[0]->dataPtr()->iSpin()==PDT::Spin0) {
    // spnin info for the outgoing scalar
    ScalarWaveFunction(decay[0],outgoing,true,vertex);
    // non-resonant form factors ( without D*)
    Complex hnr  = 0.5*_g/_fpi*xi/mb/md/(dotv+dmb);
    Complex A1nr =-0.5*_g/_fpi*xi*(1.+omega)/(dotv+dmb);
    Complex A2nr = 0.5*_g/_fpi*xi/mb*(dotv+dotvp)/(dotv+dmb);
    Complex A3nr =0.;
    // add D* if needed
    if(_includeDstar) {
      hnr  -= 0.5*_g/_fpi*xi/mb/md/(dotvp-dmd);
      A1nr += 0.5*_g/_fpi*xi*(1.+omega)/(dotvp-dmd);
      A3nr +=-0.5*_g/_fpi*xi/md*(dotv+dotvp)/(dotvp-dmd);
    }
    // resonant pieces
    Complex hr  = _alpha2*rho2/6./_fpi/mb/md*(omega-1.)/(dotv+dmt2)
      +0.5*_alpha3/_fpi/md/mb*xi1/(dotv+dmt3);
    Complex A1r =-_alpha2*rho2/6./_fpi*(omega*omega-1.)/(dotv+dmt2)
      -0.5*_alpha3*xi1/_fpi*(1.+omega)/(dotv+dmt3);
    Complex A2r = 0.5*_alpha1*rho1/_fpi/mb*dotv/(dotv+dmt1)
      +_alpha2*rho2/_fpi/mb/(dotv+dmt2)/3.*(0.5*(omega*dotvp-dotv)+dotvp-omega*dotv)
      +0.5*_alpha3*xi1/_fpi/mb*(dotv+dotvp)/(dotv+dmt3);
    Complex A3r =-0.5*_alpha1*rho1/_fpi/md*dotv/(dotv+dmt1)
      -0.5*_alpha2*rho2/_fpi/md/(dotv+dmt2)*(dotvp-omega*dotv);
    // add D** resonances
    if(_includeDstarstar) {
      hr  -=_alpha2*rho2/6./_fpi/mb/md*(omega-1.)/(dotvp-dmt2)
	+0.5*_alpha3/_fpi/md/mb*xi1/(dotvp-dmt3);
      A1r +=_alpha2*rho2/6./_fpi*(omega*omega-1.)/(dotvp-dmt2)
	+0.5*_alpha3*xi1/_fpi*(1.+omega)/(dotvp-dmt3);
      A2r += 0.5*_alpha1*rho1/_fpi/mb*dotvp/(dotvp-dmt1)
	+0.5*_alpha2*rho2/_fpi/mb*(dotv-omega*dotvp)/(dotvp-dmt2);
      A3r +=-0.5*_alpha1*rho1/_fpi/md*dotvp/(dotvp-dmt1)
	-_alpha2*rho2/_fpi/md/3./(dotvp-dmt2)*
	(0.5*(omega*dotv-dotvp)+dotv-omega*dotvp)
	-0.5*_alpha3*xi1/_fpi/md*(dotv+dotvp)/(dotvp-dmt3);
    }
    // construct the hadron vector
    hadron.push_back(-ii*(hnr+hr)*EpsFunction::product(inpart.momentum(),
						       decay[0]->momentum(),
						       decay[1]->momentum())
		     +(A1nr+A1r)*decay[1]->momentum()
		     +(A2nr+A2r)*inpart.momentum()
		     +(A3nr+A3r)*decay[0]->momentum());
  }
  // for D* pi
  else {
    // construct the spin info and get the polarization vectors
    vector<LorentzPolarizationVector> eps;
    VectorWaveFunction(eps,decay[0],outgoing,true,false,vertex);
    // non-resonant pieces (D and D* left in)
    Complex f3nr(0.),f4nr(0.),g1nr(0.),g2nr(0.),g3nr(0.),g5nr(0.);
    Complex h1nr = -_g*xi/_fpi/mb/md*dotv/(dotv+dmb);
    Complex h2nr = -_g*xi/_fpi/mb/(dotv+dmb);
    Complex h3nr = -_g*xi/_fpi/md*(1./(dotv+dmb)-(1.+omega)/dotvp);
    Complex f1nr = -0.5*_g*xi/_fpi/mb*(1./(dotv+dmb)-1./(dotvp+dmd));
    Complex f2nr = mb/md*f1nr;
    Complex f5nr = 0.5*_g*xi/_fpi/mb/md*(1.+dotv/(dotv+dmb));
    Complex f6nr = 0.5*_g*xi/_fpi/mb*(1./(dotv+dmb)-1./dotvp);
    Complex knr  = 0.5*_g*xi/_fpi*((dotvp-omega*dotv)/(dotv+dmb)+
				   (dotv-omega*dotvp)/dotvp);
    Complex g4nr = _g*xi/_fpi/md/dotvp;
    // resonant part (D** removed)
    Complex f3r(0.),f4r(0.),g1r(0.),g2r(0.),g5r(0.);
    Complex h1r =-_alpha1*rho1/_fpi/md/mb*dotv/(dotv+dmt1)
      +_alpha2*rho2/3./_fpi/mb/md*(dotv*(1.+2.*omega)-dotvp)/(dotv+dmt2)
      -_alpha3*xi1/_fpi/mb/md*dotv/(dotv+dmt3);
    Complex h2r = -_alpha2*rho2*(1.+omega)/3./_fpi/mb/(dotv+dmt2)
      -_alpha3*xi1/_fpi/mb/(dotv+dmt3);
    Complex h3r = _alpha2*rho2/3./_fpi/md*(1.+omega)/(dotv+dmt2)
      -_alpha3*xi1/_fpi/md/(dotv+dmt3);
    Complex f1r = -_alpha2*rho2*(omega-1.)/6./_fpi/mb/(dotv+dmt2)
      -0.5*_alpha3*xi1/_fpi/mb/(dotv+dmt3);
    Complex f5r = 
      0.5*_alpha1*rho1/_fpi/mb/md*dotv/(dotv+dmt1)+
      0.5*_alpha2*rho2/_fpi/mb/md*(dotvp-dotv*(1.+2.*omega)/3.)/(dotv+dmt2)+
      0.5*_alpha3*xi1 /_fpi/mb/md*dotv/(dotv+dmt3);
    Complex f6r = 
      _alpha2*rho2*(omega-1.)/6./_fpi/mb/(dotv+dmt2)+
      _alpha3*xi1            /2./_fpi/mb/(dotv+dmt3);
    Complex kr = 
      -_alpha1/2.*rho1*(omega-1.)/_fpi*dotv              /(dotv+dmt1)
      -_alpha2/3.*rho2*(omega-1.)/_fpi*(dotvp-omega*dotv)/(dotv+dmt2)
      +_alpha3/2.*xi1            /_fpi*(dotvp-omega*dotv)/(dotv+dmt3);
    Complex g4r = 2.*_alpha2*rho2/3./_fpi/md/(dotv+dmt2);
    // D** part of the resonant part
    if(_includeDstarstar) {
      h1r += -_alpha1*rho1/_fpi/md/mb*dotvp/(dotvp-dmt1)	
	+0.5*_alpha2*rho2/_fpi/mb/md*(omega*dotvp-dotv)/(dotvp-dmt2);
      h3r += -_alpha2*rho2/6./_fpi/md*(omega*omega-1.)/(dotvp-dmt2)
	+_alpha3*xi1/_fpi/md*(1.+omega)/(dotvp-dmt3);
      f1r +=  _alpha2*rho2*(omega-1.)/6./_fpi/mb/(dotvp-dmt2)
	+0.5*_alpha3*xi1/_fpi/mb/(dotvp-dmt3);
      f5r += 
	0.5*_alpha1*rho1/_fpi/mb/md*dotvp/(dotvp-dmt1)+
	0.5*_alpha2*rho2/_fpi/mb/md*(dotv-dotvp*(1.+2.*omega)/3.)/(dotvp-dmt2)+
	0.5*_alpha3*xi1 /_fpi/mb/md*dotvp/(dotvp-dmt3);
      f6r += 
	-_alpha2*rho2*(omega-1.)/6./_fpi/mb/(dotvp-dmt2)
	-_alpha3*xi1            /2./_fpi/mb/(dotvp-dmt3);
      kr  += 
	-_alpha1/2.*rho1*(omega-1.)/_fpi*dotvp             /(dotvp-dmt1)
	-_alpha2/3.*rho2*(omega-1.)/_fpi*(dotv-omega*dotvp)/(dotvp-dmt2)
	+_alpha3/2.*xi1            /_fpi*(dotv-omega*dotvp)/(dotvp-dmt3);
      g4r +=-_alpha2*rho2/3./_fpi/md*(1.+0.5*omega)/(dotvp-dmt2)+
	_alpha3*xi1/_fpi/md/(dotvp-dmt3);
    }
    Complex f2r=mb/md*f1r;
    Complex g3r=-g2r;
    complex<Energy> pieps,Beps; 
    for(unsigned int ix=0;ix<3;++ix) {
      pieps = eps[ix]*decay[1]->momentum();
      Beps  = eps[ix]*inpart.momentum();
      hadron.push_back((knr+kr)*eps[ix]+
		       (f1nr+f1r)*pieps*inpart.momentum()+
		       (f2nr+f2r)*pieps*decay[0]->momentum()+
		       (f3nr+f3r)*pieps*decay[1]->momentum()+
		       (f4nr+f4r)*Beps *inpart.momentum()+
		       (f5nr+f5r)*Beps *decay[0]->momentum()+
		       (f6nr+f6r)*Beps *decay[1]->momentum()
		       +ii*0.5*((g1r+g1nr)*pieps+(g2r+g2nr)*Beps)/mb/md*
		       EpsFunction::product(inpart.momentum(),decay[0]->momentum(),
					    decay[1]->momentum())
		       +ii*0.5*(eps[ix]*EpsFunction::product(inpart.momentum(),
							     decay[0]->momentum(),
							     decay[1]->momentum()))
		       /mb/md*((g3r+g3nr)*inpart.momentum()+
			       (g4r+g4nr)*decay[0]->momentum()+
			       (g5r+g5nr)*decay[1]->momentum())
		       +ii*0.5*
		       ((h1r+h1nr)*EpsFunction::product(eps[ix],
							inpart.momentum(),
							decay[0]->momentum())+
			(h2r+h2nr)*EpsFunction::product(eps[ix],
							inpart.momentum(),
						       decay[1]->momentum())+
			(h3r+h3nr)*EpsFunction::product(eps[ix],
							decay[0]->momentum(),
							decay[1]->momentum()))
		       );
    }
  }
  // construct the lepton current
  Energy scale;
  ParticleVector leptons;
  leptons.push_back(decay[2]);leptons.push_back(decay[3]);
  int mode=(abs(decay[1]->id())-11)/2;
  vector<LorentzPolarizationVector> lepton(_current->current(vertex,mode,
							     ichan,scale,leptons));
  // calculate the matrix elemet
  DecayMatrixElement newME(PDT::Spin0,decay[0]->dataPtr()->iSpin(),PDT::Spin0,
			   PDT::Spin1Half,PDT::Spin1Half);
  vector<unsigned int> ihel(decay.size()+1);
  ihel[0]=0;
  ihel[2]=0;
  for(ihel[1]=0;ihel[1]<hadron.size();++ihel[1]) {
    for(ihel[3]=0;ihel[3]<2;++ihel[3]) {
      for(ihel[4]=0;ihel[4]<2;++ihel[4]) {
	newME(ihel)= lepton[2*ihel[3]+ihel[4]]*hadron[ihel[1]];
      }
    }
  }
  RhoDMatrix temp(PDT::Spin0); temp.average();
  // store the matrix element
  ME(newME);
  double me=0.5*(newME.contract(temp)).real()*_GF*_GF*md*mb*SM().CKM(1,2);
  if(abs(decay[1]->id())==ParticleID::piplus) me*=2.;
  return me*inpart.mass()*inpart.mass();
}
