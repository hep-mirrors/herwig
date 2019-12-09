// -*- C++ -*-
//
// SMTopDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMTopDecayer class.
//

#include "SMTopDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
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
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

SMTopDecayer::SMTopDecayer() 
  : _wquarkwgt(6,0.),_wleptonwgt(3,0.), _xg_sampling(1.5), 
    _initialenhance(1.),  _finalenhance(2.3) {
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
  os << FFWVertex_ << FFGVertex_ << FFPVertex_ << WWWVertex_
     << _wquarkwgt << _wleptonwgt << _wplus
     << _initialenhance << _finalenhance << _xg_sampling;
}
  
void SMTopDecayer::persistentInput(PersistentIStream & is, int) {
  is >> FFWVertex_ >> FFGVertex_ >> FFPVertex_ >> WWWVertex_
     >> _wquarkwgt >> _wleptonwgt >> _wplus
     >> _initialenhance >> _finalenhance >> _xg_sampling;
}
  
// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SMTopDecayer,PerturbativeDecayer>
describeHerwigSMTopDecayer("Herwig::SMTopDecayer", "HwPerturbativeDecay.so");
  
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
    // fix rho if no correlations
    fixRho(_rho);
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
	inter = FFWVertex_->evaluate(scale,1,Wplus,_inHalf[thel],
				   _inHalfBar[bhel]);
	for(afhel=0;afhel<2;++afhel){
	  for(fhel=0;fhel<2;++fhel){
	    (*ME())(thel,bhel,afhel,fhel) = 
	      FFWVertex_->evaluate(scale,_outHalf[afhel],
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
	inter = FFWVertex_->
	  evaluate(scale,1,Wminus,_inHalf[bbhel],_inHalfBar[tbhel]);
	for(afhel=0;afhel<2;++afhel){
	  for(fhel=0;fhel<2;++fhel){
	    (*ME())(tbhel,bbhel,fhel,afhel) = 
	      FFWVertex_->evaluate(scale,_outHalf[afhel],
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
  PerturbativeDecayer::doinit();
  //get vertices from SM object
  tcHwSMPtr hwsm = dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  if(!hwsm) throw InitException() << "Must have Herwig::StandardModel in "
				  << "SMTopDecayer::doinit()";
  FFWVertex_ = hwsm->vertexFFW();
  FFGVertex_ = hwsm->vertexFFG();
  FFPVertex_ = hwsm->vertexFFP();
  WWWVertex_ = hwsm->vertexWWW();
  //initialise
  FFWVertex_->init();
  FFGVertex_->init();
  FFPVertex_->init();
  WWWVertex_->init();
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
    mode->addChannel(Wchannel);
    addMode(mode,_wleptonwgt[(i-11)/2],wgt);
  }
  //quark modes
  unsigned int iz=0;
  for(int ix=1;ix<6;ix+=2) {
    for(int iy=2;iy<6;iy+=2) {
      // check that the combination of particles is allowed
      if(FFWVertex_->allowed(-ix,iy,ParticleID::Wminus)) {
	extpart[2] = getParticleData(-ix);
	extpart[3] = getParticleData( iy);
	mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
	Wchannel = new_ptr(DecayPhaseSpaceChannel(mode));
	Wchannel->addIntermediate(extpart[0],0,0.0,-1,1);
	Wchannel->addIntermediate(_wplus,0,0.0,2,3);
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
  // parameters for the PerturbativeDecayer base class
  for(unsigned int ix=0;ix<_wquarkwgt.size();++ix) {
    os << "newdef " << name() << ":QuarkWeights " << ix << " "
	   << _wquarkwgt[ix] << "\n";
  }
  for(unsigned int ix=0;ix<_wleptonwgt.size();++ix) {
    os << "newdef " << name() << ":LeptonWeights " << ix << " "
	   << _wleptonwgt[ix] << "\n";
  }
  PerturbativeDecayer::dataBaseOutput(os,false);
  if(header) os << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}

void SMTopDecayer::doinitrun() {
  PerturbativeDecayer::doinitrun();
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
  if(born->bornOutgoing().size()!=2) return;
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


bool SMTopDecayer::softMatrixElementVeto(PPtr parent,
					 PPtr progenitor,
					 const bool & ,
					 const Energy & highestpT,
					 const vector<tcPDPtr> &,
					 const double & z,
					 const Energy & scale,
					 const Energy & pt) {
  // check if we need to apply the full correction
  // the initial-state correction
  if(abs(progenitor->id())==ParticleID::t&&abs(parent->id())==ParticleID::t) {
    // check if hardest so far
    // if not just need to remove effect of enhancement
    bool veto(false);
    // if not hardest so far
    if(pt<highestpT)
      veto=!UseRandom::rndbool(1./_initialenhance);
    // if hardest so far do calculation
    else {
      // values of kappa and z
      double kappa(sqr(scale/_mt));
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
	if(useMEforT2()&&xg>xgbcut(_ktb)) wgt = 0.;
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
    }
    // return the veto
    return veto;
  }
  // final-state correction
  else if(abs(progenitor->id())==ParticleID::b&&abs(parent->id())==ParticleID::b) {
    // check if hardest so far
    // if not just need to remove effect of enhancement
    // if not hardest so far
    if(pt<highestpT) return !UseRandom::rndbool(1./_finalenhance);
    // if hardest so far do calculation
    // values of kappa and z
    double kappa(sqr(scale/_mt));
    // momentum fractions
    double xa(1.+_a-_c-z*(1.-z)*kappa),r(0.5*(1.+_c/(1.+_a-xa))),root(sqr(xa)-4.*_a);
    if(root<0.) {
      generator()->log() << "Imaginary root for final-state veto in "
			 << "SMTopDecayer::softMatrixElementVeto"
			 << "\nz =  " << z  << "\nkappa = " << kappa
			 << "\nxa = " << xa 
			 << "\nroot^2= " << root;
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
    // compute veto from weight and return
    return !UseRandom::rndbool(wgt);
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

double SMTopDecayer::loME(const Particle & inpart, const ParticleVector & decay) {
  // spinors
  vector<SpinorWaveFunction   > swave;
  vector<SpinorBarWaveFunction> awave;
  vector<VectorWaveFunction> vwave;
  tPPtr Wboson = abs(decay[0]->id())==ParticleID::Wplus ? decay[0] : decay[1];
  tPPtr bquark = abs(decay[0]->id())==ParticleID::Wplus ? decay[1] : decay[0];
  // spinors 
  if(inpart.id()>0) {
    SpinorWaveFunction   ::calculateWaveFunctions(swave,const_ptr_cast<tPPtr>(&inpart),
						  incoming);
    SpinorBarWaveFunction::calculateWaveFunctions(awave,bquark,outgoing);
  }
  else {
    SpinorBarWaveFunction::calculateWaveFunctions(awave,const_ptr_cast<tPPtr>(&inpart),
						  incoming);
    SpinorWaveFunction   ::calculateWaveFunctions(swave,bquark,outgoing);
  }
  // polarization vectors
  VectorWaveFunction::calculateWaveFunctions(vwave,Wboson,outgoing,false);
  Energy2 scale(sqr(inpart.mass()));
  double me=0.;
  if(inpart.id() == ParticleID::t) {
    for(unsigned int thel = 0; thel < 2; ++thel) {
      for(unsigned int bhel = 0; bhel < 2; ++bhel) {
	for(unsigned int whel = 0; whel < 3; ++whel) {
	  Complex diag = FFWVertex_->evaluate(scale,swave[thel],awave[bhel],vwave[whel]);
	  me += norm(diag);
	}
      }
    }
  }
  else if(inpart.id() == ParticleID::tbar) {
    for(unsigned int thel = 0; thel < 2; ++thel) {
      for(unsigned int bhel = 0; bhel < 2; ++bhel){ 
	for(unsigned int whel = 0; whel < 3; ++whel) {
	  Complex diag = FFWVertex_->evaluate(scale,swave[bhel],awave[thel],vwave[whel]);
	  me += norm(diag);
 	}
      }
    }
  }
  return me;
}


double SMTopDecayer::realME(const Particle & inpart, const ParticleVector & decay,
			    ShowerInteraction inter) {
  // vertex for emission from fermions
  AbstractFFVVertexPtr vertex = inter==ShowerInteraction::QCD ? FFGVertex_ : FFPVertex_;
  // spinors
  vector<SpinorWaveFunction   > swave;
  vector<SpinorBarWaveFunction> awave;
  vector<VectorWaveFunction> vwave,gwave;
  tPPtr Wboson = abs(decay[0]->id())==ParticleID::Wplus ? decay[0] : decay[1];
  tPPtr bquark = abs(decay[0]->id())==ParticleID::Wplus ? decay[1] : decay[0];
  // spinors 
  if(inpart.id()>0) {
    SpinorWaveFunction   ::calculateWaveFunctions(swave,const_ptr_cast<tPPtr>(&inpart),
						  incoming);
    SpinorBarWaveFunction::calculateWaveFunctions(awave,bquark,outgoing);
  }
  else {
    SpinorBarWaveFunction::calculateWaveFunctions(awave,const_ptr_cast<tPPtr>(&inpart),
						  incoming);
    SpinorWaveFunction   ::calculateWaveFunctions(swave,bquark,outgoing);
  }
  // polarization vectors
  VectorWaveFunction::calculateWaveFunctions(vwave,Wboson,outgoing,false);
  VectorWaveFunction::calculateWaveFunctions(gwave,decay[2],outgoing,true );
  Energy2 scale(sqr(inpart.mass()));
  double me=0.;
  vector<Complex> diag(3,0.);
  if(inpart.id() == ParticleID::t) {
    for(unsigned int thel = 0; thel < 2; ++thel) {
      for(unsigned int bhel = 0; bhel < 2; ++bhel) {
	for(unsigned int whel = 0; whel < 3; ++whel) {
	  for(unsigned int ghel =0; ghel <3; ghel+=2) {
	    // emission from top
	    SpinorWaveFunction interF = vertex->evaluate(scale,3,inpart.dataPtr(),swave[thel],gwave[ghel]);
	    diag[0] = FFWVertex_->evaluate(scale,interF,awave[bhel],vwave[whel]);
	    // emission from bottom
	    SpinorBarWaveFunction  interB = vertex->evaluate(scale,3,bquark->dataPtr()->CC(),awave[bhel],gwave[ghel]);
	    diag[1] = FFWVertex_->evaluate(scale,swave[thel],interB,vwave[whel]);
	    // emission from W
	    if(inter==ShowerInteraction::QED) {
	      VectorWaveFunction interV = WWWVertex_->evaluate(scale,3,Wboson->dataPtr()->CC(),vwave[whel],gwave[ghel]);
	      diag[1] = FFWVertex_->evaluate(scale,swave[thel],awave[bhel],interV);
	    }
	    Complex sum = std::accumulate(diag.begin(),diag.end(),Complex(0.));
	    me += norm(sum);
	  }
	}
      }
    }
  }
  else if(inpart.id() == ParticleID::tbar) {
    for(unsigned int thel = 0; thel < 2; ++thel) {
      for(unsigned int bhel = 0; bhel < 2; ++bhel){ 
	for(unsigned int whel = 0; whel < 3; ++whel) {
	  for(unsigned int ghel =0; ghel <3; ghel+=2) {
	    // emission from top
	    SpinorBarWaveFunction  interB = vertex->evaluate(scale,3,inpart.dataPtr(),awave[thel],gwave[ghel]);
	    diag[1] = FFWVertex_->evaluate(scale,swave[bhel],interB,vwave[whel]);
	    // emission from bottom
	    SpinorWaveFunction interF = vertex->evaluate(scale,3,bquark->dataPtr()->CC(),swave[bhel],gwave[ghel]);
	    diag[0] = FFWVertex_->evaluate(scale,interF,awave[thel],vwave[whel]);
	    // emission from W
	    if(inter==ShowerInteraction::QED) {
	      VectorWaveFunction interV = WWWVertex_->evaluate(scale,3,Wboson->dataPtr()->CC(),vwave[whel],gwave[ghel]);
	      diag[1] = FFWVertex_->evaluate(scale,swave[bhel],awave[thel],interV);
	    }
	    Complex sum = std::accumulate(diag.begin(),diag.end(),Complex(0.));
	    me += norm(sum);
	  }
	}
      }
    }
  }
  // divide out the coupling
  me /= norm(vertex->norm());
  // return the total
  return me;
}

double SMTopDecayer::matrixElementRatio(const Particle & inpart,
					const ParticleVector & decay2,
					const ParticleVector & decay3,
					MEOption ,
					ShowerInteraction inter) {
  double Nc = standardModel()->Nc();
  double Cf = (sqr(Nc) - 1.) / (2.*Nc);  
  // if(inter==ShowerInteraction::QED) return 0.;
  // double f  = (1. + sqr(e2()) - 2.*sqr(s2()) + s2() + s2()*e2() - 2.*e2());
  // 
  // 
  // double B  = f/s2();
  
  // Energy2 PbPg = decay3[0]->momentum()*decay3[2]->momentum();
  // Energy2 PtPg = inpart.momentum()*decay3[2]->momentum();
  // Energy2 PtPb = inpart.momentum()*decay3[0]->momentum();

  // double R = Cf *((-4.*sqr(mb())*f/s2()) * ((sqr(mb())*e2()/sqr(PbPg)) + 
  // 		  (sqr(mb())/sqr(PtPg)) - 2.*(PtPb/(PtPg*PbPg))) +
  // 		  (16. + 8./s2() + 8.*e2()/s2()) * ((PtPg/PbPg) + (PbPg/PtPg)) -
  // 		  (16./s2()) * (1. + e2()));
  // return R/B*Constants::pi;
  double Bnew = loME(inpart,decay2);
  double Rnew = realME(inpart,decay3,inter);
  double output = Rnew/Bnew*4.*Constants::pi*sqr(inpart.mass())*UnitRemoval::InvE2;
  if(inter==ShowerInteraction::QCD) output *= Cf;
  return output;
}
