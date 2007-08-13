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
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/epsilon.h"
#include "ThePEG/Helicity/LorentzTensor.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

using namespace Herwig;
using namespace ThePEG::Helicity;
using namespace ThePEG::Helicity;

GoityRobertsDecayer::GoityRobertsDecayer() {
  // fermi constant
  _GF = 1.16637E-5/GeV2;
  //  Include the D* in the B to D pi decay
  _includeDstar=false;
  //  Include the D** in the B to D(*) pi decay
  _includeDstarstar=false;
  // wavefunction beta parameters
  _beta1S=0.285*GeV;
  _beta2S=0.285*GeV;
  _beta1P=0.28*GeV;
  _beta1D=0.26*GeV;
  // The pion decay constant
  _fpi = 92.4*MeV;
  // The mass difference for the mesons
  _deltaM2S=0.563*GeV;
  _deltaM1P=0.392*GeV;
  _deltaM1D=0.709*GeV;
  // the widths for the mesons
  _gamma2S= 191*MeV;
  _gamma1P=1040*MeV;
  _gamma1D= 405*MeV;
  // Lambdabar parameter for the form factors.
  _lambdabar=0.75*GeV;
  // The couplings for the decays
  _g      = 0.50;
  _alpha1 =-1.43;
  _alpha2 =-0.14;
  _alpha3 = 0.69;
  // Mass differences and widths
  _deltaMb = 45.78*MeV;
  _gammaB  = 0.010*MeV; 
  _gammaD0 = 0.060*MeV;
  _gammaDp = 0.088*MeV;
}

void GoityRobertsDecayer::doinitrun() {
  unsigned int ix,iy;
  _current->initrun();
  DecayIntegrator::doinitrun();
  _weights.clear();_wgtloc.clear();_wgtmax.clear();
  for(ix=0;ix<numberModes();++ix) {
    _wgtmax.push_back(mode(ix)->maxWeight());
    _wgtloc.push_back(_weights.size());
    for(iy=0;iy<mode(ix)->numberChannels();++iy) {
      _weights.push_back(mode(ix)->channelWeight(iy));
    }
  }
}

void GoityRobertsDecayer::persistentOutput(PersistentOStream & os) const {
  os << _current << ounit(_GF,1./GeV2) << _includeDstar << _includeDstarstar 
     << ounit(_beta1S,GeV) << ounit(_beta2S,GeV) << ounit(_beta1P,GeV) 
     << ounit(_beta1D,GeV) << ounit(_fpi,GeV) << ounit(_deltaM2S,GeV)
     << ounit(_deltaM1P,GeV) << ounit(_deltaM1D,GeV) << ounit(_gamma2S,GeV)
     << ounit(_gamma1P,GeV) << ounit(_gamma1D,GeV) << ounit(_lambdabar,GeV) 
     << _g << _alpha1 << _alpha2 << _alpha3 << _wgtloc << _wgtmax << _weights 
     << ounit(_deltaMb,GeV) << ounit(_gammaB,GeV) << ounit(_gammaD0,GeV) << ounit(_gammaDp,GeV);
}

void GoityRobertsDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _current >> iunit(_GF,1./GeV2) >> _includeDstar >> _includeDstarstar 
     >> iunit(_beta1S,GeV) >> iunit(_beta2S,GeV) >> iunit(_beta1P,GeV) 
     >> iunit(_beta1D,GeV) >> iunit(_fpi,GeV) >> iunit(_deltaM2S,GeV) 
     >> iunit(_deltaM1P,GeV) >> iunit(_deltaM1D,GeV) >> iunit(_gamma2S,GeV) 
     >> iunit(_gamma1P,GeV) >> iunit(_gamma1D,GeV) >> iunit(_lambdabar,GeV) 
     >> _g >> _alpha1 >> _alpha2 >> _alpha3 >> _wgtloc >> _wgtmax >> _weights 
     >> iunit(_deltaMb,GeV) >> iunit(_gammaB,GeV) >> iunit(_gammaD0,GeV) >> iunit(_gammaDp,GeV);
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
     &GoityRobertsDecayer::_beta1S, GeV, 0.285*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<GoityRobertsDecayer,Energy> interfaceBeta2S
    ("Beta2S",
     "The beta parameter for the 2S states",
     &GoityRobertsDecayer::_beta2S, GeV, 0.285*GeV, 0.0*GeV, 10.0*GeV,
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
     &GoityRobertsDecayer::_lambdabar, GeV, 0.75*GeV, 0.0*GeV, 10.0*GeV,
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

  static Parameter<GoityRobertsDecayer,Energy> interfaceDeltaMB
    ("DeltaMB",
     "The mass difference between the B and B* mesons",
     &GoityRobertsDecayer::_deltaMb, MeV, 45.78*MeV, 0.0*MeV, 1000.0*MeV,
     false, false, Interface::limited);

  static Parameter<GoityRobertsDecayer,Energy> interfaceBstarWidth
    ("BstarWidth",
     "The width of the B*",
     &GoityRobertsDecayer::_gammaB, MeV, 0.01*MeV, 0.0*MeV, 1.0*MeV,
     false, false, Interface::limited);

  static Parameter<GoityRobertsDecayer,Energy> interfaceDstarPlusWidth
    ("DstarPlusWidth",
     "The width of the D*+",
     &GoityRobertsDecayer::_gammaDp, MeV, 0.088*MeV, 0.0*MeV, 1.0*MeV,
     false, false, Interface::limited);

  static Parameter<GoityRobertsDecayer,Energy> interfaceDstar0Width
    ("Dstar0Width",
     "The width of the D*+",
     &GoityRobertsDecayer::_gammaD0, MeV, 0.060*MeV, 0.0*MeV, 1.0*MeV,
     false, false, Interface::limited);
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
  // effective D* resonance parameters
  Energy Mstar = sqrt(extpart[1]->mass()*(2.*D01m->mass()-extpart[1]->mass())
		      +sqr(extpart[2]->mass()));
  Energy Gstar = _gammaD0*extpart[1]->mass()/Mstar;
  for(iy=0;iy<_current->numberOfModes();++iy) {
    channel.clear();extpart.resize(3);
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
    // add the D* diagrams
    channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
    channel.back()->addIntermediate(extpart[0],0,0.0,-1,-2);
    channel.back()->addIntermediate(D01m,0,0.0,1,2);
    // add the D** diagrams
    channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
    channel.back()->addIntermediate(extpart[0],0,0.0,-1,-2);
    channel.back()->addIntermediate(D00p,0,0.0,1,2);
    channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
    channel.back()->addIntermediate(extpart[0],0,0.0,-1,-2);
    channel.back()->addIntermediate(D01p,0,0.0,1,2);
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
    mode->resetIntermediate(D01m,Mstar,Gstar);
    addMode(mode,maxweight,channelwgts);
  }
  // then B- to D0 pi0
  extpart[0]=getParticleData(ParticleID::Bminus);
  extpart[1]=getParticleData(ParticleID::D0);
  extpart[2]=getParticleData(ParticleID::pi0);
  min = extpart[0]->massMax()-extpart[1]->massMin()-extpart[2]->massMin();
  Wcharge = extpart[0]->iCharge()-extpart[1]->iCharge()-extpart[2]->iCharge();
  // effective D* resonance parameters
  Mstar = sqrt(extpart[1]->mass()*(2.*D01m->mass()-extpart[1]->mass())
		      +sqr(extpart[2]->mass()));
  Gstar = _gammaD0*extpart[1]->mass()/Mstar;
  for(iy=0;iy<_current->numberOfModes();++iy) {
    channel.clear();extpart.resize(3);
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
    // add the D* diagrams
    channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
    channel.back()->addIntermediate(extpart[0],0,0.0,-1,-2);
    channel.back()->addIntermediate(D01m,0,0.0,1,2);
    // add the D** diagrams
    channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
    channel.back()->addIntermediate(extpart[0],0,0.0,-1,-2);
    channel.back()->addIntermediate(D00p,0,0.0,1,2);
    channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
    channel.back()->addIntermediate(extpart[0],0,0.0,-1,-2);
    channel.back()->addIntermediate(D01p,0,0.0,1,2);
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
    mode->resetIntermediate(D01m,Mstar,Gstar);
    addMode(mode,maxweight,channelwgts);
  }
  // then B0bar  to D0 pi+
  extpart[0]=getParticleData(ParticleID::Bbar0);
  extpart[1]=getParticleData(ParticleID::D0);
  extpart[2]=getParticleData(ParticleID::piplus);
  min = extpart[0]->massMax()-extpart[1]->massMin()-extpart[2]->massMin();
  Wcharge = extpart[0]->iCharge()-extpart[1]->iCharge()-extpart[2]->iCharge();
  // effective D* resonance parameters
  Mstar =  sqrt(extpart[1]->mass()*(2.*Dp1m->mass()-extpart[1]->mass())
		      +sqr(extpart[2]->mass()));
  Gstar = _gammaDp*extpart[1]->mass()/Mstar;
  for(iy=0;iy<_current->numberOfModes();++iy) {
    channel.clear();extpart.resize(3);
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
    // add the D* diagrams
    channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
    channel.back()->addIntermediate(extpart[0],0,0.0,-1,-2);
    channel.back()->addIntermediate(Dp1m,0,0.0,1,2);
    // add the D** diagrams
    channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
    channel.back()->addIntermediate(extpart[0],0,0.0,-1,-2);
    channel.back()->addIntermediate(Dp0p,0,0.0,1,2);
    channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
    channel.back()->addIntermediate(extpart[0],0,0.0,-1,-2);
    channel.back()->addIntermediate(Dp1p,0,0.0,1,2);
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
    mode->resetIntermediate(Dp1m,Mstar,Gstar);
    addMode(mode,maxweight,channelwgts);
  }
  // first B- to D+ pi0
  extpart[0]=getParticleData(ParticleID::Bbar0);
  extpart[1]=getParticleData(ParticleID::Dplus);
  extpart[2]=getParticleData(ParticleID::pi0);
  min = extpart[0]->massMax()-extpart[1]->massMin()-extpart[2]->massMin();
  Wcharge = extpart[0]->iCharge()-extpart[1]->iCharge()-extpart[2]->iCharge();
  // effective D* resonance parameters
  Mstar =  sqrt(extpart[1]->mass()*(2.*Dp1m->mass()-extpart[1]->mass())
		      +sqr(extpart[2]->mass()));
  Gstar = _gammaDp*extpart[1]->mass()/Mstar;
  for(iy=0;iy<_current->numberOfModes();++iy) {
    channel.clear();extpart.resize(3);
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
    // add the D* diagrams
    channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
    channel.back()->addIntermediate(extpart[0],0,0.0,-1,-2);
    channel.back()->addIntermediate(Dp1m,0,0.0,1,2);
    // add the D** diagrams
    channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
    channel.back()->addIntermediate(extpart[0],0,0.0,-1,-2);
    channel.back()->addIntermediate(Dp0p,0,0.0,1,2);
    channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
    channel.back()->addIntermediate(extpart[0],0,0.0,-1,-2);
    channel.back()->addIntermediate(Dp1p,0,0.0,1,2);
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
    mode->resetIntermediate(Dp1m,Mstar,Gstar);
    addMode(mode,maxweight,channelwgts);
  }
  // B- to D*+ pi-
  extpart[0]=getParticleData(ParticleID::Bminus);
  extpart[1]=getParticleData(ParticleID::Dstarplus);
  extpart[2]=getParticleData(ParticleID::piminus);
  min = extpart[0]->massMax()-extpart[1]->massMin()-extpart[2]->massMin();
  Wcharge = extpart[0]->iCharge()-extpart[1]->iCharge()-extpart[2]->iCharge();
  for(iy=0;iy<_current->numberOfModes();++iy) {
    channel.clear();extpart.resize(3);
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
    // add the D* diagrams
    channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
    channel.back()->addIntermediate(extpart[0],0,0.0,-1,-2);
    channel.back()->addIntermediate(D01m,0,0.0,1,2);
    // add the D** diagrams
    channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
    channel.back()->addIntermediate(extpart[0],0,0.0,-1,-2);
    channel.back()->addIntermediate(D00p,0,0.0,1,2);
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
    channel.clear();extpart.resize(3);
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
    // add the D* diagrams
    channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
    channel.back()->addIntermediate(extpart[0],0,0.0,-1,-2);
    channel.back()->addIntermediate(D01m,0,0.0,1,2);
    // add the D** diagrams
    channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
    channel.back()->addIntermediate(extpart[0],0,0.0,-1,-2);
    channel.back()->addIntermediate(D00p,0,0.0,1,2);
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
    channel.clear();extpart.resize(3);
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
    // add the D* diagrams
    channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
    channel.back()->addIntermediate(extpart[0],0,0.0,-1,-2);
    channel.back()->addIntermediate(Dp1m,0,0.0,1,2);
    // add the D** diagrams
    channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
    channel.back()->addIntermediate(extpart[0],0,0.0,-1,-2);
    channel.back()->addIntermediate(Dp0p,0,0.0,1,2);
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
    channel.clear();extpart.resize(3);
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
    // add the D* diagrams
    channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
    channel.back()->addIntermediate(extpart[0],0,0.0,-1,-2);
    channel.back()->addIntermediate(Dp1m,0,0.0,1,2);
    // add the D** diagrams
    channel.push_back(new_ptr(DecayPhaseSpaceChannel(mode)));
    channel.back()->addIntermediate(extpart[0],0,0.0,-1,-2);
    channel.back()->addIntermediate(Dp0p,0,0.0,1,2);
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
  output << "set " << fullName() << ":Lambdabar " << _lambdabar/GeV << " \n";
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

bool GoityRobertsDecayer::accept(tcPDPtr parent,
				 const PDVector & children) const {
  if(!_current) return false;
  bool allowed(false);
  int id0(parent->id()),idtemp,idD(0),idp(0);
  vector<int> idother;
  // check number of decay products
  if(children.size()!=4) return false;
  PDVector::const_iterator pit(children.begin()),pend(children.end());
  for( ; pit!=pend;++pit) {
    idtemp=(**pit).id();
    if(abs(idtemp)<=16){idother.push_back(idtemp);}
    else if(abs(idtemp/100)==4){idD=idtemp;}
    else {idp=idtemp;}
  }
  if(idother.size()!=2||idD==0||idp==0) return false;
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
  if(!allowed) return false;
  // and the leptons
  return _current->accept(idother);
}

int  GoityRobertsDecayer::modeNumber(bool & cc,tcPDPtr parent,
				     const PDVector & children) const {
  // find the ids of the particles for the decay current
  PDVector::const_iterator pit = children.begin();
  PDVector::const_iterator pend = children.end();
  int id0(parent->id()),idD(0),idtemp;
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
  // gcc-3.3 workaround
  tPPtr temp1 = const_ptr_cast<tPPtr>(&inpart);
  ScalarWaveFunction(temp1,incoming,true,vertex);
  // spin info for the outgoing pion
  tPPtr temp2 = decay[1];
  ScalarWaveFunction(temp2,outgoing,true,vertex);
  // calculate some common variables
  Energy mb(inpart.mass()),md(decay[0]->mass());
  double omega(inpart.momentum()*decay[0]->momentum()/mb/md);
  Complex ii(0.,1.);
  // dot products we will need
  Energy dotv (inpart   .momentum()*decay[1]->momentum()/mb),
         dotvp(decay[0]->momentum()*decay[1]->momentum()/md);
  // delta M parameters
  complex<Energy> dmt1(_deltaM1P-0.5*ii*_gamma1P);
  complex<Energy> dmt2(_deltaM1D-0.5*ii*_gamma1D);
  complex<Energy> dmt3(_deltaM2S-0.5*ii*_gamma2S);
  // calculate the mass splittings for the lightest multiplet
  // calculate the mass differences
  complex<Energy> dmb=_deltaMb-ii*0.5*_gammaB;
  complex<Energy> dmd=-decay[0]->mass();
  double epssign = decay[2]->id()<0 ? 1. : -1.;
  if(decay[1]->id()==ParticleID::pi0) {
    dmd = getParticleData((abs(decay[0]->id())/10)*10+3)->mass()
      -getParticleData((abs(decay[0]->id())/10)*10+1)->mass();
    dmd += ((abs(decay[0]->id())/10)*10+3)==ParticleID::Dstarplus 
      ? -ii*0.5*_gammaDp : -ii*0.5*_gammaD0;
  }
  else {
    int iout   = 400 + ((abs(decay[0]->id())-400)/10)*10;
    int iother = iout==410 ? 420 : 410; 
    if(abs(decay[0]->id())%10==1) {
      dmd = getParticleData(iother+3)->mass()-
	    getParticleData(iout+1  )->mass();
    }
    else {
      dmd = getParticleData(iout+3  )->mass()-
	    getParticleData(iother+1)->mass();
    }
    dmd += iout==410 ? -ii*0.5*_gammaDp : -ii*0.5*_gammaD0;
  }
  // get the IW form factors and calculate some useful prefactors
  double xi,xi1,rho1,rho2;
  calculateFormFactors(omega,xi,xi1,rho1,rho2);
//   InvEnergy xfact  = 0.5*_g/_fpi*xi;
//   InvEnergy a1fact = _alpha1*rho1/2./_fpi;
//   InvEnergy a2fact = _alpha2*rho2/6./_fpi;
//   InvEnergy a3fact = _alpha3*xi1 /2./_fpi;
   InvEnergy xfact  = 0./MeV;
   InvEnergy a1fact = _alpha1*rho1/2./_fpi;
   InvEnergy a2fact = 0./MeV;
   InvEnergy a3fact = 0./MeV;
  complex<Energy> dotvb  = dotv  + dmb;
  complex<Energy> dotvpd = dotvp - dmd;
  // hadronic current
  vector<LorentzVector<complex<InvEnergy> > > hadron;
  // calculate the matrix element
  // for D pi
  if(decay[0]->dataPtr()->iSpin()==PDT::Spin0) {
    // spin info for the outgoing scalar
    // gcc 3.3 workaround
    tPPtr temp = decay[0];
    ScalarWaveFunction(temp,outgoing,true,vertex);
    // non-resonant form factors ( without D*)
    complex<InvEnergy2> hnr  = xfact             /dotvb;
    complex<InvEnergy2> A1nr =-xfact*(1.+omega)  /dotvb;
    complex<InvEnergy > A2nr = xfact*(dotv+dotvp)/dotvb;
    complex<InvEnergy > A3nr = 0./MeV;
    // add D* if needed
    if(_includeDstar) {
      hnr  +=-xfact             /dotvpd;
      A1nr += xfact*(1.+omega)  /dotvpd;
      A3nr +=-xfact*(dotv+dotvp)/dotvpd;
    }
    // resonant pieces
    complex<InvEnergy2> hr  = a2fact*(omega-1.)/(dotv+dmt2)
                             +a3fact           /(dotv+dmt3);
    complex<InvEnergy2> A1r =-a2fact*(sqr(omega)-1.)/(dotv+dmt2)
                             -a3fact*(1.+omega     )/(dotv+dmt3);
    complex<InvEnergy > A2r = 
       a1fact/(dotv+dmt1)*dotv
      +a2fact/(dotv+dmt2)*(omega*dotvp-dotv+2.*(dotvp-omega*dotv))
      +a3fact/(dotv+dmt3)*(dotv+dotvp);
    complex<InvEnergy > A3r =
      -   a1fact/(dotv+dmt1)*dotv
      -3.*a2fact/(dotv+dmt2)*(dotvp-omega*dotv);
    // add D** resonances
    if(_includeDstarstar) {
      hr  +=-a2fact/(dotvp-dmt2)*(omega-1.)
	    -a3fact/(dotvp-dmt3);
      A1r += a2fact/(dotvp-dmt2)*(omega*omega-1.)
	    +a3fact/(dotvp-dmt3)*(1.+omega);
      A2r +=    a1fact/(dotvp-dmt1)*dotvp
	    +3.*a2fact/(dotvp-dmt2)*(dotv-omega*dotvp);
      A3r +=-a1fact/(dotvp-dmt1)*dotvp
	    -a2fact/(dotvp-dmt2)*(omega*dotv-dotvp+2.*(dotv-omega*dotvp))
	    -a3fact/(dotvp-dmt3)*(dotv+dotvp);
    }
    hadron.resize(1);
    // construct the hadron vector
    hadron[0] =
      -ii*epssign*(hnr+hr)/mb/md*Helicity::epsilon(inpart.momentum(),
						   decay[0]->momentum(),
						   decay[1]->momentum())
      +(A1nr+A1r)*decay[1]->momentum()
      +(A2nr+A2r)/mb*inpart.momentum()
      +(A3nr+A3r)/md*decay[0]->momentum();
  }
  // for D* pi
  else {
    // construct the spin info and get the polarization vectors
    vector<LorentzPolarizationVector> eps;
    VectorWaveFunction(eps,decay[0],outgoing,true,false,vertex);
    // non-resonant pieces (D and D* left in)
    complex<InvEnergy > h1nr =-2.*xfact*dotv/dotvb;
    complex<InvEnergy2> h2nr =-2.*xfact     /dotvb;
    complex<InvEnergy2> h3nr =-2.*xfact*(1./dotvb-(1.+omega)/dotvp);
    complex<InvEnergy2> f1nr =-   xfact*(1./dotvb-1./(dotvp+dmd));
    complex<InvEnergy > f5nr =    xfact*(1.+dotv/dotvb);
    complex<InvEnergy2> f6nr =    xfact*(1./dotvb-1./dotvp);
    complex<InvEnergy > knr  =    xfact*(complex<double>((dotvp-omega*dotv)/dotvb)+
					 complex<double>((dotv-omega*dotvp)/dotvp));
    complex<InvEnergy2> g2nr = 0. /MeV/MeV;
    complex<InvEnergy2> g4nr = 2.*xfact/dotvp;
    // resonant part (D** removed)
    complex<InvEnergy > h1r =-2.*a1fact/(dotv+dmt1)*dotv
                             +2.*a2fact/(dotv+dmt2)*(dotv*(1.+2.*omega)-dotvp)
                             -2.*a3fact/(dotv+dmt3)*dotv;
    complex<InvEnergy2> h2r =-2.*a2fact/(dotv+dmt2)*(1.+omega)
                             -2.*a3fact/(dotv+dmt3);
    complex<InvEnergy2> h3r = 2.*a2fact/(dotv+dmt2)*(1.+omega)
                             -2.*a3fact/(dotv+dmt3);
    complex<InvEnergy2> f1r =-   a2fact/(dotv+dmt2)*(omega-1.)
                             -   a3fact/(dotv+dmt3);
    complex<InvEnergy > f5r =    a1fact/(dotv+dmt1)*dotv
                             +   a2fact/(dotv+dmt2)*(3.*dotvp-dotv*(1.+2.*omega))
                             +   a3fact/(dotv+dmt3)*dotv;
    complex<InvEnergy2> f6r =    a2fact/(dotv+dmt2)*(omega-1.)
                             +   a3fact/(dotv+dmt3);
    complex<InvEnergy > kr  =-   a1fact/(dotv+dmt1)*(omega-1.)*dotv              
                             -2.*a2fact/(dotv+dmt2)*(omega-1.)*(dotvp-omega*dotv)
                             +   a3fact/(dotv+dmt3)           *(dotvp-omega*dotv);
    complex<InvEnergy2> g2r = 0. /MeV/MeV;
    complex<InvEnergy2> g4r = 4.*a2fact/(dotv+dmt2);
    // D** part of the resonant part
    if(_includeDstarstar) {
      h1r += -2.*a1fact/(dotvp-dmt1)*dotvp	
	     +3.*a2fact/(dotvp-dmt2)*(omega*dotvp-dotv);
      h3r += -   a2fact/(dotvp-dmt2)*(sqr(omega)-1.)
	     +2.*a3fact/(dotvp-dmt3)*(1.+omega);
      f1r +=     a2fact/(dotvp-dmt2)*(omega-1.)
	     +   a3fact/(dotvp-dmt3);
      f5r +=     a1fact/(dotvp-dmt1)*dotvp
	     +   a2fact/(dotvp-dmt2)*(3.*dotv-dotvp*(1.+2.*omega))
	     +   a3fact/(dotvp-dmt3)*dotvp;
      f6r += -   a2fact/(dotvp-dmt2)*(omega-1.)
	     -   a3fact/(dotvp-dmt3);
      kr  += -   a1fact/(dotvp-dmt1)*(omega-1.)*dotvp             
	     -2.*a2fact/(dotvp-dmt2)*(omega-1.)*(dotv-omega*dotvp)
	     +   a3fact/(dotvp-dmt3)*           (dotv-omega*dotvp);
      g2r += -3.*a2fact/(dotvp-dmt2);
      g4r += -   a2fact/(dotvp-dmt2)*(2.+omega)
	     +2.*a3fact/(dotvp-dmt3);
    }
    complex<Energy> pieps,Beps;
    hadron.resize(3);
    for(unsigned int ix=0;ix<3;++ix) {
      pieps = eps[ix]*decay[1]->momentum();
      Beps  = eps[ix]*inpart.momentum();
      Complex iii=ii*epssign;
      hadron[ix] = 
 	(knr+kr)*eps[ix]+
 	(f1nr+f1r)/mb*pieps*inpart.momentum()+
 	(f1nr+f1r)/md*pieps*decay[0]->momentum()+
 	(f5nr+f5r)/mb/md*Beps *decay[0]->momentum()+
 	(f6nr+f6r)/mb*Beps *decay[1]->momentum()
	  	+iii*0.5*((g2r+g2nr)*Beps)/mb/mb/md*
 	Helicity::epsilon(inpart.momentum(),
 			     decay[0]->momentum(),
 			     decay[1]->momentum())
	+iii*0.5/mb/md*(eps[ix].dot(Helicity::epsilon(inpart.momentum(),
						      decay[0]->momentum(),
						      decay[1]->momentum())))*
	(-(g2r+g2nr)/mb*inpart.momentum()
	 +(g4r+g4nr)/md*decay[0]->momentum())
 	+iii*0.5*
 	((h1r+h1nr)/mb/md*Helicity::epsilon(eps[ix],
 					       inpart.momentum(),
 					       decay[0]->momentum())+
 	 (h2r+h2nr)/mb   *Helicity::epsilon(eps[ix],
 					       inpart.momentum(),
 					       decay[1]->momentum())+
 	 (h3r+h3nr)/md   *Helicity::epsilon(eps[ix],
 					       decay[0]->momentum(),
 					       decay[1]->momentum()));
    }
  }
  // construct the lepton current
  Energy scale;
  ParticleVector leptons;
  leptons.push_back(decay[2]);leptons.push_back(decay[3]);
  int mode=(abs(decay[2]->id())-11)/2;
  vector<LorentzPolarizationVectorE> lepton(_current->current(vertex,mode,
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
	newME(ihel)= lepton[2*ihel[3]+ihel[4]].dot(hadron[ihel[1]]);
      }
    }
  }
  RhoDMatrix temp(PDT::Spin0); temp.average();
  // store the matrix element
  ME(newME);
  double me=0.5*(newME.contract(temp)).real()*md*mb*SM().CKM(1,2)*sqr(inpart.mass()*_GF);
  if(abs(decay[1]->id())==ParticleID::piplus) me*=2.;
  return me;
}
