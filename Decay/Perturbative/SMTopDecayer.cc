// -*- C++ -*-
//
// SMTopDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMTopDecayer class.
//

#include "SMTopDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig/Decay/DecayVertex.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/PDT/ThreeBodyAllOn1IntegralCalculator.h"
#include "Herwig/Shower/RealEmissionProcess.h"
#include "Herwig/Shower/Core/Base/ShowerProgenitor.h"
#include "Herwig/Shower/Core/Base/ShowerParticle.h"
#include "Herwig/Shower/Core/Base/Branching.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

SMTopDecayer::SMTopDecayer() 
  : _wquarkwgt(6,0.),_wleptonwgt(3,0.), _xg_sampling(1.5), 
    _initialenhance(1.),  _finalenhance(2.3), _useMEforT2(true) {
  _wleptonwgt[0] = 0.302583;
  _wleptonwgt[1] = 0.301024;
  _wleptonwgt[2] = 0.299548;
  _wquarkwgt[0]  = 0.851719;
  _wquarkwgt[1]  = 0.0450162;
  _wquarkwgt[2]  = 0.0456962;
  _wquarkwgt[3]  = 0.859839;
  _wquarkwgt[4]  = 3.9704e-06;
  _wquarkwgt[5]  = 0.000489657; 
  generateIntermediates(true);
}
  
bool SMTopDecayer::accept(tcPDPtr parent, const tPDVector & children) const {
  if(abs(parent->id()) != ParticleID::t) return false;
  int id0(0),id1(0),id2(0);
  for(tPDVector::const_iterator it = children.begin();
      it != children.end();++it) {
    int id=(**it).id(),absid(abs(id));
    if(absid==ParticleID::b&&double(id)/double(parent->id())>0) {
      id0=id;
    }
    else {
      switch (absid) {
      case ParticleID::nu_e: 
      case ParticleID::nu_mu:
      case ParticleID::nu_tau:
	id1 = id;
	break;
      case ParticleID::eminus:
      case ParticleID::muminus:
      case ParticleID::tauminus:
	id2 = id;
	break;
      case ParticleID::b:
      case ParticleID::d:
      case ParticleID::s:
	id1 = id;
	break;
      case ParticleID::u:
      case ParticleID::c:
	id2=id;
	break;
      default :
	break;
      }
    }
  }
  if(id0==0||id1==0||id2==0) return false;
  if(double(id1)/double(id2)>0) return false;
  return true;
}
  
ParticleVector SMTopDecayer::decay(const Particle & parent,
				   const tPDVector & children) const {
  int id1(0),id2(0);
  for(tPDVector::const_iterator it = children.begin();
      it != children.end();++it) {
    int id=(**it).id(),absid=abs(id);
    if(absid == ParticleID::b && double(id)/double(parent.id())>0) continue;
    //leptons
    if(absid > 10 && absid%2==0) id1=absid;
    if(absid > 10 && absid%2==1) id2=absid;
    //quarks
    if(absid < 10 && absid%2==0) id2=absid;
    if(absid < 10 && absid%2==1) id1=absid;
  }
  unsigned int imode(0);
  if(id2 >=11 && id2<=16) imode = (id1-12)/2;
  else imode = id1+1+id2/2;
  bool cc = parent.id() == ParticleID::tbar;
  ParticleVector out(generate(true,cc,imode,parent));
  //arrange colour flow
  PPtr pparent=const_ptr_cast<PPtr>(&parent);
  out[1]->incomingColour(pparent,out[1]->id()<0);
  ParticleVector products = out[0]->children();
  if(products[0]->hasColour())
    products[0]->colourNeighbour(products[1],true);
  else if(products[0]->hasAntiColour())
    products[0]->colourNeighbour(products[1],false);
  return out;
}   
 
void SMTopDecayer::persistentOutput(PersistentOStream & os) const {
  os << _wvertex << _wquarkwgt << _wleptonwgt << _wplus << _alpha
     << _initialenhance << _finalenhance << _xg_sampling << _useMEforT2; 
}
  
void SMTopDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _wvertex >> _wquarkwgt >> _wleptonwgt >> _wplus >> _alpha
     >> _initialenhance >> _finalenhance >> _xg_sampling >> _useMEforT2;
}
  
ClassDescription<SMTopDecayer> SMTopDecayer::initSMTopDecayer;
// Definition of the static class description member.
  
void SMTopDecayer::Init() {
    
  static ClassDocumentation<SMTopDecayer> documentation
    ("This is the implementation of the SMTopDecayer which "
     "decays top quarks into bottom quarks and either leptons  "
     "or quark-antiquark pairs including the matrix element for top decay",
     "The matrix element correction for top decay \\cite{Hamilton:2006ms}.",
     "%\\cite{Hamilton:2006ms}\n"
     "\\bibitem{Hamilton:2006ms}\n"
     "  K.~Hamilton and P.~Richardson,\n"
     "  ``A simulation of QCD radiation in top quark decays,''\n"
     "  JHEP {\\bf 0702}, 069 (2007)\n"
     "  [arXiv:hep-ph/0612236].\n"
     "  %%CITATION = JHEPA,0702,069;%%\n");
  
  static ParVector<SMTopDecayer,double> interfaceQuarkWeights
    ("QuarkWeights",
     "Maximum weights for the hadronic decays",
     &SMTopDecayer::_wquarkwgt, 6, 1.0, 0.0, 10.0,
     false, false, Interface::limited);
  
  static ParVector<SMTopDecayer,double> interfaceLeptonWeights
    ("LeptonWeights",
     "Maximum weights for the semi-leptonic decays",
     &SMTopDecayer::_wleptonwgt, 3, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<SMTopDecayer,double> interfaceEnhancementFactor
    ("InitialEnhancementFactor",
     "The enhancement factor for initial-state radiation in the shower to ensure"
     " the weight for the matrix element correction is less than one.",
     &SMTopDecayer::_initialenhance, 1.0, 1.0, 10000.0,
     false, false, Interface::limited);

  static Parameter<SMTopDecayer,double> interfaceFinalEnhancementFactor
    ("FinalEnhancementFactor",
     "The enhancement factor for final-state radiation in the shower to ensure"
     " the weight for the matrix element correction is less than one",
     &SMTopDecayer::_finalenhance, 1.6, 1.0, 1000.0,
     false, false, Interface::limited);

  static Parameter<SMTopDecayer,double> interfaceSamplingTopHardMEC
    ("SamplingTopHardMEC",
     "The importance sampling power for choosing an initial xg, "
     "to sample xg according to xg^-_xg_sampling",
     &SMTopDecayer::_xg_sampling, 1.5, 1.2, 2.0,
     false, false, Interface::limited);

  static Switch<SMTopDecayer,bool> interfaceUseMEForT2
    ("UseMEForT2",
     "Use the matrix element correction, if available to fill the T2"
     " region for the decay shower and don't fill using the shower",
     &SMTopDecayer::_useMEforT2, true, false, false);
  static SwitchOption interfaceUseMEForT2Shower
    (interfaceUseMEForT2,
     "Shower",
     "Use the shower to fill the T2 region",
     false);
  static SwitchOption interfaceUseMEForT2ME
    (interfaceUseMEForT2,
     "ME",
     "Use the Matrix element to fill the T2 region",
     true);

  static Reference<SMTopDecayer,ShowerAlpha> interfaceCoupling
    ("Coupling",
     "Pointer to the object to calculate the coupling for the correction",
     &SMTopDecayer::_alpha, false, false, true, false, false);
}

double SMTopDecayer::me2(const int, const Particle & inpart,
			 const ParticleVector & decay,
			 MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,
					 PDT::Spin1Half,PDT::Spin1Half)));
  // spinors etc for the decaying particle
  if(meopt==Initialize) {
    // spinors and rho
    if(inpart.id()>0)
      SpinorWaveFunction   ::calculateWaveFunctions(_inHalf,_rho,
						    const_ptr_cast<tPPtr>(&inpart),
						    incoming);
    else
      SpinorBarWaveFunction::calculateWaveFunctions(_inHalfBar,_rho,
						    const_ptr_cast<tPPtr>(&inpart),
						    incoming);
  }
  // setup spin info when needed
  if(meopt==Terminate) {
    // for the decaying particle
    if(inpart.id()>0) {
      SpinorWaveFunction::
	constructSpinInfo(_inHalf,const_ptr_cast<tPPtr>(&inpart),incoming,true);
      SpinorBarWaveFunction::constructSpinInfo(_inHalfBar,decay[0],outgoing,true);
      SpinorWaveFunction   ::constructSpinInfo(_outHalf   ,decay[1],outgoing,true);
      SpinorBarWaveFunction::constructSpinInfo(_outHalfBar,decay[2],outgoing,true);
    }
    else {
      SpinorBarWaveFunction::
	constructSpinInfo(_inHalfBar,const_ptr_cast<tPPtr>(&inpart),incoming,true);
      SpinorWaveFunction::constructSpinInfo(_inHalf,decay[0],outgoing,true);
      SpinorBarWaveFunction::constructSpinInfo(_outHalfBar,decay[1],outgoing,true);
      SpinorWaveFunction   ::constructSpinInfo(_outHalf   ,decay[2],outgoing,true);
    }
  }

  if ( ( decay[1]->momentum() + decay[2]->momentum() ).m()
       < decay[1]->data().constituentMass() + decay[2]->data().constituentMass() )
    return 0.0;

  // spinors for the decay product
  if(inpart.id()>0) {
    SpinorBarWaveFunction::calculateWaveFunctions(_inHalfBar ,decay[0],outgoing);
    SpinorWaveFunction   ::calculateWaveFunctions(_outHalf   ,decay[1],outgoing);
    SpinorBarWaveFunction::calculateWaveFunctions(_outHalfBar,decay[2],outgoing);
  }
  else {
    SpinorWaveFunction   ::calculateWaveFunctions(_inHalf    ,decay[0],outgoing);
    SpinorBarWaveFunction::calculateWaveFunctions(_outHalfBar,decay[1],outgoing);
    SpinorWaveFunction   ::calculateWaveFunctions(_outHalf   ,decay[2],outgoing);
  }
  Energy2 scale(sqr(inpart.mass()));
  if(inpart.id() == ParticleID::t) {
    //Define intermediate vector wave-function for Wplus 
    tcPDPtr Wplus(getParticleData(ParticleID::Wplus));
    VectorWaveFunction inter;
    unsigned int thel,bhel,fhel,afhel;
    for(thel = 0;thel<2;++thel){
      for(bhel = 0;bhel<2;++bhel){	  
	inter = _wvertex->evaluate(scale,1,Wplus,_inHalf[thel],
				   _inHalfBar[bhel]);
	for(afhel=0;afhel<2;++afhel){
	  for(fhel=0;fhel<2;++fhel){
	    (*ME())(thel,bhel,afhel,fhel) = 
	      _wvertex->evaluate(scale,_outHalf[afhel],
				 _outHalfBar[fhel],inter);
	  }
	}
      }
    }
  }
  else if(inpart.id() == ParticleID::tbar) {
    VectorWaveFunction inter;
    tcPDPtr Wminus(getParticleData(ParticleID::Wminus));
    unsigned int tbhel,bbhel,afhel,fhel;
    for(tbhel = 0;tbhel<2;++tbhel){
      for(bbhel = 0;bbhel<2;++bbhel){
	inter = _wvertex->
	  evaluate(scale,1,Wminus,_inHalf[bbhel],_inHalfBar[tbhel]);
	for(afhel=0;afhel<2;++afhel){
	  for(fhel=0;fhel<2;++fhel){
	    (*ME())(tbhel,bbhel,fhel,afhel) = 
	      _wvertex->evaluate(scale,_outHalf[afhel],
				 _outHalfBar[fhel],inter);
	  }
	}
      }
    }
  }
  double output = (ME()->contract(_rho)).real();
  if(abs(decay[1]->id())<=6) output *=3.;
  return output;
}

void SMTopDecayer::doinit() {
  DecayIntegrator::doinit();
  //get vertices from SM object
  tcHwSMPtr hwsm = dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  if(!hwsm) throw InitException() << "Must have Herwig::StandardModel in "
				  << "SMTopDecayer::doinit()";
  _wvertex = hwsm->vertexFFW();
  //initialise
  _wvertex->init();
  //set up decay modes
  _wplus = getParticleData(ParticleID::Wplus);
  DecayPhaseSpaceModePtr mode;
  DecayPhaseSpaceChannelPtr Wchannel;
  tPDVector extpart(4);
  vector<double> wgt(1,1.0);
  extpart[0] = getParticleData(ParticleID::t);
  extpart[1] = getParticleData(ParticleID::b);
  //lepton modes
  for(int i=11; i<17;i+=2) {
    extpart[2] = getParticleData(-i);
    extpart[3] = getParticleData(i+1);
    mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
    Wchannel = new_ptr(DecayPhaseSpaceChannel(mode));
    Wchannel->addIntermediate(extpart[0],0,0.0,-1,1);
    Wchannel->addIntermediate(_wplus,0,0.0,2,3);
    Wchannel->init();
    mode->addChannel(Wchannel);
    addMode(mode,_wleptonwgt[(i-11)/2],wgt);
  }
  //quark modes
  unsigned int iz=0;
  for(int ix=1;ix<6;ix+=2) {
    for(int iy=2;iy<6;iy+=2) {
      // check that the combination of particles is allowed
      if(_wvertex->allowed(-ix,iy,ParticleID::Wminus)) {
	extpart[2] = getParticleData(-ix);
	extpart[3] = getParticleData( iy);
	mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
	Wchannel = new_ptr(DecayPhaseSpaceChannel(mode));
	Wchannel->addIntermediate(extpart[0],0,0.0,-1,1);
	Wchannel->addIntermediate(_wplus,0,0.0,2,3);
	Wchannel->init();
	mode->addChannel(Wchannel);
	addMode(mode,_wquarkwgt[iz],wgt);
	++iz;
      }
      else {
	throw InitException() << "SMTopDecayer::doinit() the W vertex" 
			      << "cannot handle all the quark modes" 
			      << Exception::abortnow;
      }
    }
  }
}

void SMTopDecayer::dataBaseOutput(ofstream & os,bool header) const {
  if(header) os << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  for(unsigned int ix=0;ix<_wquarkwgt.size();++ix) {
    os << "newdef " << name() << ":QuarkWeights " << ix << " "
	   << _wquarkwgt[ix] << "\n";
  }
  for(unsigned int ix=0;ix<_wleptonwgt.size();++ix) {
    os << "newdef " << name() << ":LeptonWeights " << ix << " "
	   << _wleptonwgt[ix] << "\n";
  }
  DecayIntegrator::dataBaseOutput(os,false);
  if(header) os << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}

void SMTopDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<numberModes();++ix) {
      if(ix<3) _wleptonwgt[ix  ] = mode(ix)->maxWeight();
      else     _wquarkwgt [ix-3] = mode(ix)->maxWeight();
    }
  }
}

WidthCalculatorBasePtr SMTopDecayer::threeBodyMEIntegrator(const DecayMode & dm) const {
  // identify W decay products
  int sign = dm.parent()->id() > 0 ? 1 : -1;
  int iferm(0),ianti(0);
  for(ParticleMSet::const_iterator pit=dm.products().begin();
      pit!=dm.products().end();++pit) {
    int id = (**pit).id();
    if(id*sign != ParticleID::b) {
      if   (id*sign > 0 ) iferm = id*sign;
      else                ianti = id*sign;
    }
  }
  assert(iferm!=0&&ianti!=0);
  // work out which mode we are doing
  int imode(-1);
  for(unsigned int ix=0;ix<numberModes();++ix) {
    if(mode(ix)->externalParticles(2)->id() == ianti &&
       mode(ix)->externalParticles(3)->id() == iferm ) {
      imode = ix;
      break;
    }
  }
  assert(imode>=0);
  // get the masses we need
  Energy m[3] = {mode(imode)->externalParticles(1)->mass(),
		 mode(imode)->externalParticles(3)->mass(),
		 mode(imode)->externalParticles(2)->mass()};
  return 
    new_ptr(ThreeBodyAllOn1IntegralCalculator<SMTopDecayer>
	    (3,_wplus->mass(),_wplus->width(),0.0,*this,imode,m[0],m[1],m[2]));
}

InvEnergy SMTopDecayer::threeBodydGammads(const int imode, const Energy2 mt2,
					  const Energy2 mffb2, const Energy mb,
					  const Energy mf, const Energy mfb) const {
  Energy mffb(sqrt(mffb2));
  Energy mw(_wplus->mass());
  Energy2 mw2(sqr(mw)),gw2(sqr(_wplus->width()));
  Energy mt(sqrt(mt2));
  Energy Eb  = 0.5*(mt2-mffb2-sqr(mb))/mffb;
  Energy Ef  = 0.5*(mffb2-sqr(mfb)+sqr(mf))/mffb;
  Energy Ebm = sqrt(sqr(Eb)-sqr(mb));
  Energy Efm = sqrt(sqr(Ef)-sqr(mf));
  Energy2 upp = sqr(Eb+Ef)-sqr(Ebm-Efm);
  Energy2 low = sqr(Eb+Ef)-sqr(Ebm+Efm);
  InvEnergy width=(dGammaIntegrand(mffb2,upp,mt,mb,mf,mfb,mw)-
		   dGammaIntegrand(mffb2,low,mt,mb,mf,mfb,mw))
    /32./mt2/mt/8/pow(Constants::pi,3)/(sqr(mffb2-mw2)+mw2*gw2);
  // couplings
  width *= 0.25*sqr(4.*Constants::pi*generator()->standardModel()->alphaEM(mt2)/
		    generator()->standardModel()->sin2ThetaW());
  width *= generator()->standardModel()->CKM(*mode(imode)->externalParticles(0),
					     *mode(imode)->externalParticles(1));
  if(abs(mode(imode)->externalParticles(2)->id())<=6) {
    width *=3.;
    if(abs(mode(imode)->externalParticles(2)->id())%2==0)
      width *=generator()->standardModel()->CKM(*mode(imode)->externalParticles(2),
						*mode(imode)->externalParticles(3));
    else
      width *=generator()->standardModel()->CKM(*mode(imode)->externalParticles(3),
						*mode(imode)->externalParticles(2));
  }
  // final spin average
  assert(!std::isnan(double(width*MeV)));
  return 0.5*width;
}

Energy6 SMTopDecayer::dGammaIntegrand(Energy2 mffb2, Energy2 mbf2, Energy mt,
				      Energy mb, Energy mf, Energy mfb, Energy mw) const {
  Energy2 mt2(sqr(mt)) ,mb2(sqr(mb)) ,mf2(sqr(mf )),mfb2(sqr(mfb )),mw2(sqr(mw ));
  Energy4 mt4(sqr(mt2)),mb4(sqr(mb2)),mf4(sqr(mf2)),mfb4(sqr(mfb2)),mw4(sqr(mw2));
  return -mbf2 * ( + 6 * mb2 * mf2 * mfb2 * mffb2    +   6 * mb2 * mt2 * mfb2 * mffb2 
		   + 6 * mb2 * mt2 * mf2  * mffb2    +  12 * mb2 * mt2 * mf2 * mfb2 
		   - 3  * mb2 * mfb4  * mffb2        +   3 * mb2 * mf2 * mffb2 * mffb2 
		   - 3  * mb2 * mf4   * mffb2        -   6 * mb2 * mt2 * mfb4 
		   - 6  * mb2 * mt2 * mf4            -   3 * mb4 * mfb2 * mffb2 
		   - 3  * mb4 * mf2 * mffb2          -   6 * mb4 * mf2 * mfb2
		   + 3  * mt4 * mf4                  +   3 * mb4 * mfb4 
		   + 3  * mb4 * mf4                  +   3 * mt4 * mfb4
		   + 3  * mb2 * mfb2 * mffb2 * mffb2 +   3 * mt2 * mfb2 * mffb2 * mffb2 
		   - 3  * mt2 * mfb4 * mffb2         +   3 * mt2 * mf2 * mffb2 * mffb2 
		   - 3  * mt2 * mf4 * mffb2          -   3 * mt4 * mfb2 * mffb2 
		   - 3  * mt4 * mf2 * mffb2          -   6 * mt4 * mf2 * mfb2 
		   + 6  * mt2 * mf2 * mfb2 * mffb2   +  12 * mt2 * mf2 * mw4 
		   + 12 * mb2 * mfb2 * mw4           +  12 * mb2 * mt2 * mw4 
		   + 6  * mw2 * mt2 * mfb2 * mbf2    -  12 * mw2 * mt2 * mf2 * mffb2 
		   - 6  * mw2 * mt2 * mf2 * mbf2     -  12 * mw2 * mt2 * mf2 * mfb2 
		   - 12 * mw2 * mb2  * mfb2 * mffb2  -   6 * mw2 * mb2 * mfb2 * mbf2 
		   + 6  * mw2 * mb2  * mf2 * mbf2    -  12 * mw2 * mb2 * mf2 * mfb2 
		   - 12 * mw2 * mb2 * mt2 * mfb2     -  12 * mw2 * mb2 * mt2 * mf2 
		   + 12 * mf2 * mfb2 * mw4           +   4 * mbf2 * mbf2 * mw4 
		   -  6 * mfb2 * mbf2 * mw4          -   6 * mf2 * mbf2 * mw4 
		   -  6 * mt2 * mbf2 * mw4           -   6 * mb2 * mbf2 * mw4 
		   + 12 * mw2 * mt2 * mf4            +  12 * mw2 * mt4 * mf2 
		   + 12 * mw2 * mb2 * mfb4           +  12 * mw2 * mb4 * mfb2) /mw4 / 3.;
}

void SMTopDecayer::initializeMECorrection(RealEmissionProcessPtr born, double & initial,
					  double & final) {
  // check the outgoing particles
  PPtr part[2];
  for(unsigned int ix=0;ix<born->bornOutgoing().size();++ix) {
    part[ix]= born->bornOutgoing()[ix];
  }
  // check the final-state particles and get the masses
  if(abs(part[0]->id())==ParticleID::Wplus&&abs(part[1]->id())==ParticleID::b) {
    _ma=part[0]->mass();
    _mc=part[1]->mass();
  }
  else if(abs(part[1]->id())==ParticleID::Wplus&&abs(part[0]->id())==ParticleID::b) {
    _ma=part[1]->mass();
    _mc=part[0]->mass();
  }
  else {
    return;
  }
  // set the top mass
  _mt=born->bornIncoming()[0]->mass();
  // set the gluon mass
  _mg=getParticleData(ParticleID::g)->constituentMass();
  // set the radiation enhancement factors
  initial = _initialenhance;
  final   = _finalenhance;
  // reduced mass parameters
  _a=sqr(_ma/_mt);
  _g=sqr(_mg/_mt);
  _c=sqr(_mc/_mt);
  double lambda = sqrt(1.+sqr(_a)+sqr(_c)-2.*_a-2.*_c-2.*_a*_c);
  _ktb = 0.5*(3.-_a+_c+lambda);
  _ktc = 0.5*(1.-_a+3.*_c+lambda);
  useMe();
}

RealEmissionProcessPtr SMTopDecayer::applyHardMatrixElementCorrection(RealEmissionProcessPtr born) {
  // Get b and a and put them in particle vector ba in that order...
  ParticleVector ba; 
  for(unsigned int ix=0;ix<born->bornOutgoing().size();++ix)
    ba.push_back(born->bornOutgoing()[ix]);
  if(abs(ba[0]->id())!=5) swap(ba[0],ba[1]);
  assert(born->bornIncoming().size()==1);
  // Now decide if we get an emission into the dead region.
  // If there is an emission newfs stores momenta for a,c,g 
  // according to NLO decay matrix element. 
  vector<Lorentz5Momentum> newfs = applyHard(ba,_ktb,_ktc);
  // If there was no gluon emitted return.
  if(newfs.size()!=3) return RealEmissionProcessPtr();
  // Sanity checks to ensure energy greater than mass etc :)
  bool check = true; 
  tcPDPtr gluondata=getParticleData(ParticleID::g);
  if (newfs[0].e()<ba[0]->data().constituentMass()) check = false;
  if (newfs[1].e()<ba[1]->mass())                   check = false;
  if (newfs[2].e()<gluondata->constituentMass())    check = false;
  // Return if insane:
  if (!check) return RealEmissionProcessPtr();
  // // Set masses in 5-vectors:
  newfs[0].setMass(ba[0]->mass());
  newfs[1].setMass(ba[1]->mass());
  newfs[2].setMass(ZERO);
  // The next part of this routine sets the colour structure.
  // To do this for decays we assume that the gluon comes from c!
  // First create new particle objects for c, a and gluon:
  PPtr newg = gluondata->produceParticle(newfs[2]);
  PPtr newc = ba[0]->data().produceParticle(newfs[0]);
  PPtr newa = ba[1]->data().produceParticle(newfs[1]);
  born->spectator(0);
  born->emitted(3);
  // decaying particle
  born->incoming().push_back(born->bornIncoming()[0]->dataPtr()->
			     produceParticle(born->bornIncoming()[0]->momentum()));
  // colour flow
  newg->incomingColour(born->incoming()[0],ba[0]->id()<0);
  newg->colourConnect(newc                ,ba[0]->id()<0);
  if(born->bornOutgoing()[0]->id()==newc->id()) {
    born->outgoing().push_back(newc);
    born->outgoing().push_back(newa);
    born->emitter(1);
  }
  else {
    born->outgoing().push_back(newa);
    born->outgoing().push_back(newc);
    born->emitter(2);
  }
  born->outgoing().push_back(newg);
  // boost for the W
  LorentzRotation trans(ba[1]->momentum().findBoostToCM());
  trans.boost(newfs[1].boostVector());
  born->transformation(trans);
  if(!inTheDeadRegion(_xg,_xa,_ktb,_ktc)) {
    generator()->log()
      << "SMTopDecayer::applyHardMatrixElementCorrection()\n"
      << "Just found a point that escaped from the dead region!\n"
      << "   _xg: " << _xg << "   _xa: " << _xa 
      << "   newfs.size(): " << newfs.size() << endl;
  }
  born->interaction(ShowerInteraction::QCD);
  return born;
}

vector<Lorentz5Momentum> SMTopDecayer::
applyHard(const ParticleVector &p,double ktb, double ktc) { 
  // ********************************* //
  // First we see if we get a dead     //
  // region event: _xa,_xg             //
  // ********************************* //
  vector<Lorentz5Momentum> fs; 
  // Return if there is no (NLO) gluon emission:
  
  double weight = getHard(ktb,ktc);
  if(weight>1.) {
    generator()->log() << "Weight greater than 1 for hard emission in "
		       << "SMTopDecayer::applyHard xg = " << _xg 
		       << " xa = " << _xa << "\n";
      weight=1.;
  }
  // Accept/Reject
  if (weight<UseRandom::rnd()||p.size()!= 2) return fs; 
   // Drop events if getHard returned a negative weight 
  // as in events that, somehow have escaped from the dead region
  // or, worse, the allowed region.
  if(weight<0.) return fs;
 
  // Calculate xc by momentum conservation:
  _xc = 2.-_xa-_xg;

  // ************************************ //
  // Now we get the boosts & rotations to //
  // go from lab to top rest frame with   //
  // a in the +z direction.               //
  // ************************************ //
  Lorentz5Momentum pa_lab,pb_lab,pc_lab,pg_lab;
  // Calculate momentum of b:
  pb_lab = p[0]->momentum() + p[1]->momentum(); 
   // Define/assign momenta of c,a and the gluon:
  if(abs(p[0]->id())==5) {
    pc_lab = p[0]->momentum(); 
    pa_lab = p[1]->momentum(); 
  } else {
    pc_lab = p[1]->momentum(); 
    pa_lab = p[0]->momentum(); 
  }
  // Calculate the boost to the b rest frame:
  SpinOneLorentzRotation rot0(pb_lab.findBoostToCM());
  // Calculate the rotation matrix to position a along the +z direction
  // in the rest frame of b and does a random rotation about z:
  SpinOneLorentzRotation    rot1 = rotateToZ(rot0*pa_lab);
  // Calculate the boost from the b rest frame back to the lab:
  // and the inverse of the random rotation about the z-axis and the 
  // rotation required to align a with +z:
  SpinOneLorentzRotation invrot = rot0.inverse()*rot1.inverse();

  // ************************************ //
  // Now we construct the momenta in the  //
  // b rest frame using _xa,_xg.          //
  // First we construct b, then c and g,  //
  // finally we generate a by momentum    //
  // conservation.                        //
  // ************************************ //
  Lorentz5Momentum pa_brf, pb_brf(_mt), pc_brf, pg_brf;
  // First we set the top quark to being on-shell and at rest.
  // Second we set the energies of c and g,
  pc_brf.setE(0.5*_mt*(2.-_xa-_xg));
  pg_brf.setE(0.5*_mt*_xg);
  // then their masses,
  pc_brf.setMass(_mc);
  pg_brf.setMass(ZERO);
  // Now set the z-component of c and g. For pg we simply start from
  // _xa and _xg, while for pc we assume it is equal to minus the sum
  // of the z-components of a (assumed to point in the +z direction) and g.
  double root=sqrt(_xa*_xa-4.*_a);
  pg_brf.setZ(_mt*(1.-_xa-_xg+0.5*_xa*_xg-_c+_a)/root);
  pc_brf.setZ(-1.*( pg_brf.z()+_mt*0.5*root));
  // Now set the y-component of c and g's momenta
  pc_brf.setY(ZERO);
  pg_brf.setY(ZERO);
  // Now set the x-component of c and g's momenta
  pg_brf.setX(sqrt(sqr(pg_brf.t())-sqr(pg_brf.z())));
  pc_brf.setX(-pg_brf.x());
  // Momenta b,c,g are now set. Now we obtain a from momentum conservation,
  pa_brf = pb_brf-pc_brf-pg_brf;
  pa_brf.setMass(pa_brf.m());
  pa_brf.rescaleEnergy();
 
  // ************************************ //
  // Now we orient the momenta and boost  //
  // them back to the original lab frame. //
  // ************************************ //
  // As in herwig6507 we assume that, in the rest frame
  // of b, we have aligned the W boson momentum in the 
  // +Z direction by rot1*rot0*pa_lab, therefore
  // we obtain the new pa_lab by applying:
  // invrot*pa_brf.
  pa_lab = invrot*pa_brf;    
  pb_lab = invrot*pb_brf;    
  pc_lab = invrot*pc_brf;    
  pg_lab = invrot*pg_brf;    
  fs.push_back(pc_lab); 
  fs.push_back(pa_lab); 
  fs.push_back(pg_lab); 
  return fs;
}

double SMTopDecayer::getHard(double ktb, double ktc) {
  // zero the variables
  _xg = 0.;    
  _xa = 0.;   
  _xc = 0.;
  // Get a phase space point in the dead region: 
  double volume_factor = deadRegionxgxa(ktb,ktc);
  // if outside region return -1
  if(volume_factor<0) return volume_factor;
  // Compute the weight for this phase space point:
  double weight = volume_factor*me(_xa,_xg)*(1.+_a-_c-_xa); 
  // Alpha_S and colour factors - this hard wired Alpha_S needs removing.
  weight *= (4./3.)/Constants::pi
    *(_alpha->value(_mt*_mt*_xg*(1.-_xa+_a-_c)
		    /(2.-_xg-_xa-_c)));
  return weight; 
}

bool SMTopDecayer::softMatrixElementVeto(ShowerProgenitorPtr initial,
					 ShowerParticlePtr parent,Branching br) {
  // check if we need to apply the full correction
  long id[2]={abs(initial->progenitor()->id()),abs(parent->id())};
  // the initial-state correction
  if(id[0]==ParticleID::t&&id[1]==ParticleID::t) {
    Energy pt=br.kinematics->pT();
    // check if hardest so far
    // if not just need to remove effect of enhancement
    bool veto(false);
    // if not hardest so far
    if(pt<initial->highestpT())
      veto=!UseRandom::rndbool(1./_initialenhance);
    // if hardest so far do calculation
    else {
      // values of kappa and z
      double z(br.kinematics->z()),kappa(sqr(br.kinematics->scale()/_mt));
      // parameters for the translation
      double w(1.-(1.-z)*(kappa-1.)),u(1.+_a-_c-(1.-z)*kappa),v(sqr(u)-4.*_a*w*z);
      // veto if outside phase space
      if(v<0.) 
	veto=true;
      // otherwise calculate the weight
      else {
	v = sqrt(v);
	double xa((0.5*(u+v)/w+0.5*(u-v)/z)),xg((1.-z)*kappa);
	double f(me(xa,xg)),
	  J(0.5*(u+v)/sqr(w)-0.5*(u-v)/sqr(z)+_a*sqr(w-z)/(v*w*z));
	double wgt(f*J*2./kappa/(1.+sqr(z)-2.*z/kappa)/_initialenhance);
	// This next `if' prevents the hardest emission from the 
	// top shower ever entering the so-called T2 region of the
	// phase space if that region is to be populated by the hard MEC.
	if(_useMEforT2&&xg>xgbcut(_ktb)) wgt = 0.;
	if(wgt>1.) {
	  generator()->log() << "Violation of maximum for initial-state "
			     << " soft veto in "
			     << "SMTopDecayer::softMatrixElementVeto"
			     << "xg = " << xg << " xa = " << xa 
			     << "weight =  " << wgt << "\n";
	  wgt=1.;
	}
	// compute veto from weight
	veto = !UseRandom::rndbool(wgt);
      }
      // if not vetoed reset max
      if(!veto) initial->highestpT(pt);
    }
    // if vetoing reset the scale
    if(veto) parent->vetoEmission(br.type,br.kinematics->scale());
    // return the veto
    return veto;
  }
  // final-state correction
  else if(id[0]==ParticleID::b&&id[1]==ParticleID::b) {
    Energy pt=br.kinematics->pT();
    // check if hardest so far
    // if not just need to remove effect of enhancement
    bool veto(false);
    // if not hardest so far
    if(pt<initial->highestpT()) return !UseRandom::rndbool(1./_finalenhance);
    // if hardest so far do calculation
    // values of kappa and z
    double z(br.kinematics->z()),kappa(sqr(br.kinematics->scale()/_mt));
    // momentum fractions
    double xa(1.+_a-_c-z*(1.-z)*kappa),r(0.5*(1.+_c/(1.+_a-xa))),root(sqr(xa)-4.*_a);
    if(root<0.) {
      generator()->log() << "Imaginary root for final-state veto in "
			 << "SMTopDecayer::softMatrixElementVeto"
			 << "\nz =  " << z  << "\nkappa = " << kappa
			 << "\nxa = " << xa 
			 << "\nroot^2= " << root;
      parent->vetoEmission(br.type,br.kinematics->scale());
      return true;
    } 
    root=sqrt(root);
    double xg((2.-xa)*(1.-r)-(z-r)*root);
    // xfact (below) is supposed to equal xg/(1-z). 
    double xfact(z*kappa/2./(z*(1.-z)*kappa+_c)*(2.-xa-root)+root);
    // calculate the full result
    double f(me(xa,xg));
    // jacobian
    double J(z*root);
    double wgt(f*J*2.*kappa/(1.+sqr(z)-2.*_c/kappa/z)/sqr(xfact)/_finalenhance);
    if(wgt>1.) {
      generator()->log() << "Violation of maximum for final-state  soft veto in "
			 << "SMTopDecayer::softMatrixElementVeto"
			 << "xg = " << xg << " xa = " << xa 
			 << "weight =  " << wgt << "\n";
      wgt=1.;
    }
    // compute veto from weight
    veto = !UseRandom::rndbool(wgt);
    // if vetoing reset the scale
    if(veto) parent->vetoEmission(br.type,br.kinematics->scale());
    // return the veto
    return veto;
  }
  // otherwise don't veto
  else return !UseRandom::rndbool(1./_finalenhance);
}

double SMTopDecayer::me(double xw,double xg) {
  double prop(1.+_a-_c-xw),xg2(sqr(xg));
  double lambda=sqrt(1.+_a*_a+_c*_c-2.*_a-2.*_c-2.*_a*_c);
  double denom=(1.-2*_a*_a+_a+_c*_a+_c*_c-2.*_c);
  double wgt=-_c*xg2/prop+(1.-_a+_c)*xg-(prop*(1 - xg)+xg2)
    +(0.5*(1.+2.*_a+_c)*sqr(prop-xg)*xg+2.*_a*prop*xg2)/denom;
  return wgt/(lambda*prop);
}


// This function is auxiliary to the xab function.
double SMTopDecayer::xgbr(int toggle) { 
  return 1.+toggle*sqrt(_a)-_c*(1.-toggle*sqrt(_a))/(1.-_a);
}

// This function is auxiliary to the xab function.
double SMTopDecayer::ktr(double xgb, int toggle) { 
  return 2.*xgb/
    (xgb+toggle*sqrt((1.-1./_a)
		     *(xgb-xgbr( 1))
		     *(xgb-xgbr(-1))));
}

// Function xab determines xa (2*W energy fraction) for a given value
// of xg (2*gluon energy fraction) and kappa tilde (q tilde squared over
// m_top squared). Hence this function allows you to draw 1: the total
// phase space volume in the xa vs xg plane 2: for a given value of 
// kappa tilde (i.e. starting evolution scale) the associated contour 
// in the xa vs xg plane (and hence the regions that either shower can 
// populate). This calculation is done assuming the emission came from
// the top quark i.e. kappa tilde here is the q tilde squared of the TOP
// quark divided by m_top squared. 
double SMTopDecayer::xab(double xgb, double kt, int toggle) { 
  double xab;
  if(toggle==2) {
    // This applies for g==0.&&kt==ktr(a,c,0.,xgb,1).
    xab = -2.*_a*(xgb-2.)/(1.+_a-_c-xgb);
  } else if(toggle==1) {
    // This applies for kt==1&&g==0.
    double lambda = sqrt(sqr(xgb-1.+_a+_c)-4.*_a*_c);
    xab = (0.5/(kt-xgb))*(kt*(1.+_a-_c-xgb)-lambda)
      + (0.5/(kt+xgb*(1.-kt)))*(kt*(1.+_a-_c-xgb)+lambda);
  } else {
    // This is the form of xab FOR _g=0.
    double ktmktrpktmktrm = kt*kt - 4.*_a*(kt-1.)*xgb*xgb
                                  / (sqr(1.-_a-_c-xgb)-4.*_a*_c);
    if(fabs(kt-(2.*xgb-2.*_g)/(xgb-sqrt(xgb*xgb-4.*_g)))/kt>1.e-6) {
      double lambda = sqrt((sqr(1.-_a-_c-xgb)-4.*_a*_c)*ktmktrpktmktrm);
      xab = (0.5/(kt-xgb))*(kt*(1.+_a-_c-xgb)-lambda)
	  + (0.5/(kt+xgb*(1.-kt)))*(kt*(1.+_a-_c-xgb)+lambda);
    }
    else {
      // This is the value of xa as a function of xb when kt->infinity. 
      // Where we take any kt > (2.*xgb-2.*_g)/(xgb-sqrt(xgb*xgb-4.*_g)) 
      // as being effectively infinite. This kt value is actually the 
      // maximum allowed value kt can have if the phase space is calculated  
      // without the approximation of _g=0 (massless gluon). This formula
      // for xab below is then valid for _g=0 AND kt=infinity only.
      xab = ( 2.*_c+_a*(xgb-2.)
	    + 3.*xgb
	    - xgb*(_c+xgb+sqrt(_a*_a-2.*(_c-xgb+1.)*_a+sqr(_c+xgb-1.)))
            - 2.
            )/2./(xgb-1.);
    }
  }
  if(std::isnan(xab)) {
    double ktmktrpktmktrm = ( sqr(xgb*kt-2.*(xgb-_g))
			      -kt*kt*(1.-1./_a)*(xgb-xgbr( 1)-_g/(1.+sqrt(_a)))
			      *(xgb-xgbr(-1)-_g/(1.-sqrt(_a)))
			      )/
      (xgb*xgb-(1.-1./_a)*(xgb-xgbr( 1)-_g/(1.+sqrt(_a)))
       *(xgb-xgbr(-1)-_g/(1.-sqrt(_a)))
       );
    double lambda = sqrt((xgb-1.+sqr(sqrt(_a)+sqrt(_c-_g)))
			 *(xgb-1.+sqr(sqrt(_a)-sqrt(_c-_g)))*
			 ktmktrpktmktrm);
    xab = (0.5/(kt-xgb+_g))*(kt*(1.+_a-_c+_g-xgb)-lambda)
      + (0.5/(kt+xgb*(1.-kt)-_g))*(kt*(1.+_a-_c+_g-xgb)+lambda);
    if(std::isnan(xab)) 
	throw Exception() << "TopMECorrection::xab complex x_a value.\n"
			  << "  xgb    = " << xgb    << "\n"
			  << "  xab    = " << xab    << "\n"
			  << "  toggle = " << toggle << "\n"
			  << "  ktmktrpktmktrm = "   << ktmktrpktmktrm 
			  << Exception::eventerror;
  }
  return xab;
}

// xgbcut is the point along the xg axis where the upper bound on the 
// top quark (i.e. b) emission phase space goes back on itself in the 
// xa vs xg plane i.e. roughly mid-way along the xg axis in
// the xa vs xg Dalitz plot.
double SMTopDecayer::xgbcut(double kt) { 
  double lambda2 = 1.+_a*_a+_c*_c-2.*_a-2.*_c-2.*_a*_c; 
  double num1    = kt*kt*(1.-_a-_c);
  double num2    = 2.*kt*sqrt(_a*(kt*kt*_c+lambda2*(kt-1.)));
  return (num1-num2)/(kt*kt-4.*_a*(kt-1.));
}

double SMTopDecayer::xaccut(double kt) { 
    return 1.+_a-_c-0.25*kt;
}

double SMTopDecayer::z(double xac, double kt, 
		       int toggle1, int toggle2) { 
  double z = -1.0;
  if(toggle2==0) { 
    z = (kt+toggle1*sqrt(kt*(kt-4.*(1.+_a-_c-xac))))/(2.*kt); 
  } else if(toggle2==1) {
    z = ((1.+_a+_c-xac)+toggle1*(1.+_a-_c-xac))
      /(2.*(1.+_a-xac));
  } else if(toggle2==2) {
    z = 0.5;
  } else {
    throw Exception() << "Cannot determine z in SMTopDecayer::z()"
		      << Exception::eventerror;
  }
  return z;
}

double SMTopDecayer::xgc(double xac, double kt, 
			 int toggle1, int toggle2) { 
  double tiny(1.e-6);
  double xaToMinBoundary(xac*xac-4.*_a);
  if(xaToMinBoundary<0) {
    if(fabs(xaToMinBoundary/(1.-_a)/(1.-_a))<tiny)
      xaToMinBoundary *= -1.;
    else
      throw Exception() << "SMTopDecayer::xgc xa not in phase space!"
			<< Exception::eventerror;
  }
  return (2.-xac)*(1.-0.5*(1.+_c/(1.+_a-xac)))
        -(z(xac,kt,toggle1,toggle2)-0.5*(1.+_c/(1.+_a-xac)))
        *sqrt(xaToMinBoundary);
}

double SMTopDecayer::xginvc0(double xg , double kt) { 
  // The function xg(kappa_tilde_c,xa) surely, enough, draws a  
  // line of constant kappa_tilde_c in the xg, xa Dalitz plot. 
  // Such a function can therefore draw the upper and lower 
  // edges of the phase space for emission from c (the b-quark).
  // However, to sample the soft part of the dead zone effectively
  // we want to generate a value of xg first and THEN distribute
  // xa in the associated allowed part of the dead zone. Hence, the 
  // function we want, to define the dead zone in xa for a given
  // xg, is the inverse of xg(kappa_tilde_c,xa). The full expression 
  // for xg(kappa_tilde_c,xa) is complicated and, sure enough, 
  // does not invert. Therefore we try to overestimate the size
  // of the dead zone initially, rejecting events which do not 
  // fall exactly inside it afterwards, with the immediate aim 
  // of getting an approximate version of xg(kappa_tilde_c,xa) 
  // that can  be inverted. We do this by simply setting c=0 i.e.
  // the b-quark mass to zero (and the gluon mass of course), in 
  // the full expression xg(...). The result of inverting this 
  // function is the output of this routine (a value of xa) hence 
  // the name xginvc0. xginvc0 is calculated to be,
  // xginvc0 = (1./3.)*(1.+a+pow((U+sqrt(4.*V*V*V+U*U))/2.,1./3.)
  //                      -V*pow(2./(U+sqrt(4.*V*V*V+U*U)),1./3.)
  //                   )
  // U = 2.*a*a*a - 66.*a*a + 9.*a*kt*xg + 18.*a*kt
  //   - 66.*a + 27.*kt*xg*xg - 45.*kt*xg +18.*kt +2. ;
  // V = -1.-a*a-14.*a-3.kt*xg+3.*kt;
  // This function, as with many functions in this ME correction,
  // is plagued by cuts that have to handled carefully in numerical 
  // implementation. We have analysed the cuts and hence we implement 
  // it in the following way, with a series of 'if' statements. 
  //
  // A useful -definition- to know in deriving the v<0 terms is
  // that tanh^-1(z) = 0.5*(log(1.+z)-log(1.-z)).
  double u,v,output;
  u = 2.*_a*_a*_a-66.*_a*_a
     +9.*xg*kt*_a+18.*kt*_a
     -66.*_a+27.*xg*xg*kt
     -45.*xg*kt+18.*kt+2.;
  v = -_a*_a-14.*_a-3.*xg*kt+3.*kt-1.;
  double u2=u*u,v3=v*v*v;
  if(v<0.) {
    if(u>0.&&(4.*v3+u2)<0.)      output = cos(  atan(sqrt(-4.*v3-u2)/u)/3.);
    else if(u>0.&&(4.*v3+u2)>0.) output = cosh(atanh(sqrt( 4.*v3+u2)/u)/3.);
    else                         output = cos(( atan(sqrt(-4.*v3-u2)/u)
					       +Constants::pi)/3.);
    output *= 2.*sqrt(-v);
  } else {
    output = sinh(log((u+sqrt(4.*v3+u2))/(2.*sqrt(v3)))/3.);
    output *= 2.*sqrt(v);
  }
  if(!isfinite(output)) {
      throw Exception() << "TopMECorrection::xginvc0:\n"
	  << "possible numerical instability detected.\n"
	  << "\n v = " <<  v << "   u = " << u << "\n4.*v3+u2 = " << 4.*v3+u2
	  << "\n_a = " << _a << "  ma = " << sqrt(_a*_mt*_mt/GeV2) 
	  << "\n_c = " << _c << "  mc = " << sqrt(_c*_mt*_mt/GeV2) 
	  << "\n_g = " << _g << "  mg = " << sqrt(_g*_mt*_mt/GeV2) 
	  << Exception::eventerror;
  }
  return ( 1.+_a +output)/3.;
}

double SMTopDecayer::approxDeadMaxxa(double xg,double ktb,double ktc) {
  double maxxa(0.);
  double x  = min(xginvc0(xg,ktc),
		  xab(xg,(2.*xg-2.*_g)/(xg-sqrt(xg*xg-4.*_g)),0));
  double y(-9999999999.);
  if(xg>2.*sqrt(_g)&&xg<=xgbcut(ktb)) {
      y = max(xab(xg,ktb,0),xab(xg,1.,1));
  } else if(xg>=xgbcut(ktb)&&xg<=1.-sqr(sqrt(_a)+sqrt(_c))) {
      y = max(xab(xg,ktr(xg,1),2),xab(xg,1.,1));
  }
  if(xg>2.*sqrt(_g)&&xg<=1.-sqr(sqrt(_a)+sqrt(_c))) {
    if(x>=y) { maxxa =  x       ; }
    else     { maxxa = -9999999.; }
  } else {
    maxxa = -9999999.;
  }
  return maxxa;
}

double SMTopDecayer::approxDeadMinxa(double xg,double ktb,double ktc) {
  double minxa(0.);
  double x  = min(xginvc0(xg,ktc),
		  xab(xg,(2.*xg-2.*_g)/(xg-sqrt(xg*xg-4.*_g)),0));
  double y(-9999999999.);
  if(xg>2.*sqrt(_g)&&xg<=xgbcut(ktb)) {
      y = max(xab(xg,ktb,0),xab(xg,1.,1));
  } else if(xg>=xgbcut(ktb)&&xg<=1.-sqr(sqrt(_a)+sqrt(_c))) {
      if(_useMEforT2) y = xab(xg,1.,1);
      else            y = max(xab(xg,ktr(xg,1),2),xab(xg,1.,1));
  }
  if(xg>2.*sqrt(_g)&&xg<=1.-sqr(sqrt(_a)+sqrt(_c))) {
      if(x>=y) { minxa =  y  ; }
      else     { minxa = 9999999.; }
  } else {
      minxa = 9999999.;
  }
  return minxa;
}

// This function returns true if the phase space point (xg,xa) is in the 
// kinematically allowed phase space.
bool SMTopDecayer::inTheAllowedRegion(double xg , double xa) {
    bool output(true);
    if(xg<2.*sqrt(_g)||xg>1.-sqr(sqrt(_a)+sqrt(_c)))      output = false;
    if(xa<xab(xg,1.,1))                                   output = false;
    if(xa>xab(xg,(2.*xg-2.*_g)/(xg-sqrt(xg*xg-4.*_g)),0)) output = false;
    return output;
}

// This function returns true if the phase space point (xg,xa) is in the 
// approximate (overestimated) dead region.
bool SMTopDecayer::inTheApproxDeadRegion(double xg , double xa,
					 double ktb, double ktc) {
    bool output(true);
    if(!inTheAllowedRegion(xg,xa))       output = false;
    if(xa<approxDeadMinxa(xg,ktb,ktc))   output = false;
    if(xa>approxDeadMaxxa(xg,ktb,ktc))   output = false;
    return output;
}

// This function returns true if the phase space point (xg,xa) is in the 
// dead region.
bool SMTopDecayer::inTheDeadRegion(double xg , double xa,
				   double ktb, double ktc) {
    bool output(true);
    if(!inTheApproxDeadRegion(xg,xa,ktb,ktc)) output = false;
    if(xa>xaccut(ktc)) {
	if(xg<xgc(max(xaccut(ktc),2.*sqrt(_a)),ktc, 1,2)&&
           xg>xgc(xa,ktc, 1,0)) { output = false; } 
	if(xg>xgc(max(xaccut(ktc),2.*sqrt(_a)),ktc,-1,2)&&
           xg<xgc(xa,ktc,-1,0)) { output = false; } 
    } 
    return output;
}

// This function attempts to generate a phase space point in the dead
// region and returns the associated phase space volume factor needed for
// the associated event weight.
double SMTopDecayer::deadRegionxgxa(double ktb,double ktc) {
  _xg=0.;
  _xa=0.;
  // Here we set limits on xg and generate a value inside the bounds.
  double xgmin(2.*sqrt(_g)),xgmax(1.-sqr(sqrt(_a)+sqrt(_c)));
  // Generate _xg.
  if(_xg_sampling==2.) {
      _xg=xgmin*xgmax/(xgmin+UseRandom::rnd()*(xgmax-xgmin));
  } else {
      _xg=xgmin*xgmax/pow((  pow(xgmin,_xg_sampling-1.)
			    + UseRandom::rnd()*(pow(xgmax,_xg_sampling-1.)
					       -pow(xgmin,_xg_sampling-1.))
			   ),1./(_xg_sampling-1.));
  }
  // Here we set the bounds on _xa for given _xg.
  if(_xg<xgmin||xgmin>xgmax) 
      throw Exception() << "TopMECorrection::deadRegionxgxa:\n"
			<< "upper xg bound is less than the lower xg bound.\n"
			<< "\n_xg         = " << _xg 
			<< "\n2.*sqrt(_g) = " << 2.*sqrt(_g) 
			<< "\n_a  = " << _a  << "  ma = " << sqrt(_a*_mt*_mt/GeV2) 
			<< "\n_c  = " << _c  << "  mc = " << sqrt(_c*_mt*_mt/GeV2) 
			<< "\n_g  = " << _g  << "  mg = " << sqrt(_g*_mt*_mt/GeV2) 
			<< Exception::eventerror;
  double xamin(approxDeadMinxa(_xg,ktb,ktc));
  double xamax(approxDeadMaxxa(_xg,ktb,ktc));
  // Are the bounds sensible? If not return.
  if(xamax<=xamin) return -1.;
  _xa=1.+_a-(1.+_a-xamax)*pow((1.+_a-xamin)/(1.+_a-xamax),UseRandom::rnd());
  // If outside the allowed region return -1.
  if(!inTheDeadRegion(_xg,_xa,ktb,ktc)) return -1.;
  // The integration volume for the weight
  double xg_vol,xa_vol; 
  if(_xg_sampling==2.) {
      xg_vol = (xgmax-xgmin)
             / (xgmax*xgmin);
  } else {
      xg_vol = (pow(xgmax,_xg_sampling-1.)-pow(xgmin,_xg_sampling-1.))
	     / ((_xg_sampling-1.)*pow(xgmax*xgmin,_xg_sampling-1.));
  }
  xa_vol = log((1.+_a-xamin)/(1.+_a-xamax));
  // Here we return the integral volume factor multiplied by the part of the 
  // weight left over which is not included in the BRACES function, i.e.
  // the part of _xg^-2 which is not absorbed in the integration measure.
  return xg_vol*xa_vol*pow(_xg,_xg_sampling-2.);
}

LorentzRotation SMTopDecayer::rotateToZ(Lorentz5Momentum v) {
  // compute the rotation matrix
  LorentzRotation trans;
  // rotate so in z-y plane
  trans.rotateZ(-atan2(v.y(),v.x()));
  // rotate so along Z
  trans.rotateY(-acos(v.z()/v.vect().mag()));
  // generate random rotation
  double c,s,cs;
  do
    {
      c = 2.*UseRandom::rnd()-1.;
      s = 2.*UseRandom::rnd()-1.;
      cs = c*c+s*s;
    }
  while(cs>1.||cs==0.);
  double cost=(c*c-s*s)/cs,sint=2.*c*s/cs;
  // apply random azimuthal rotation
  trans.rotateZ(atan2(sint,cost));
  return trans;
}
