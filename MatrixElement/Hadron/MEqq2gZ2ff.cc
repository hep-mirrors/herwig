// -*- C++ -*-
//
// MEqq2gZ2ff.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEqq2gZ2ff class.
//

#include "MEqq2gZ2ff.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "Herwig/MatrixElement/HardVertex.h"

using namespace Herwig;

MEqq2gZ2ff::MEqq2gZ2ff() : _minflavour(1), _maxflavour(5), 
			   _gammaZ(0), _process(0), spinCorrelations_(true) {
  massOption(vector<unsigned int>(2,1));
}

void MEqq2gZ2ff::doinit() {
  DrellYanBase::doinit();
  _z0=getParticleData(ThePEG::ParticleID::Z0);
  _gamma=getParticleData(ThePEG::ParticleID::gamma);
  // cast the SM pointer to the Herwig SM pointer
  tcHwSMPtr hwsm=ThePEG::dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  // do the initialisation
  if(hwsm) {
    _theFFZVertex = hwsm->vertexFFZ();
    _theFFPVertex = hwsm->vertexFFP();
  }
  else
    throw InitException() << "Must be the Herwig StandardModel class in "
			  << "MEqq2gZ2ff::doinit" << Exception::abortnow;
}

void MEqq2gZ2ff::getDiagrams() const {
  // which intermediates to include
  bool gamma = _gammaZ==0 || _gammaZ==1;
  bool Z0    = _gammaZ==0 || _gammaZ==2;
  // loop over the processes we need
  for ( int ix=1; ix<17; ++ix ) {
    // increment counter to switch between quarks and leptons
    if(ix==7) ix+=4;
    // is it a valid quark process
    bool quark = ix<=6 && (_process==0 || _process==1 || _process-10==ix);
    // is it a valid lepton process
    bool lepton= ix>=11 && ix<=16 
      && (_process==0
	  || _process==2
	  || (_process==3 && ix%2==1)
	  || (_process==4 && ix%2==0)
	  || (ix%2==0 && (ix-10)/2==_process-7)
	  || (ix%2==1 && (ix-9)/2 ==_process-4)
	  );
    // if not a valid process continue
    if(!(quark||lepton)) continue;
    tcPDPtr lm = getParticleData(ix);
    tcPDPtr lp = lm->CC();
    for(int i = _minflavour; i <= _maxflavour; ++i) {
      tcPDPtr q  = getParticleData(i);
      tcPDPtr qb = q->CC();
      if(Z0)    add(new_ptr((Tree2toNDiagram(2), q, qb, 1, _z0   , 3, lm, 3, lp, -1)));
      if(gamma) add(new_ptr((Tree2toNDiagram(2), q, qb, 1, _gamma, 3, lm, 3, lp, -2)));
    }
  }
}

Energy2 MEqq2gZ2ff::scale() const {
  return sHat();
}

double MEqq2gZ2ff::me2() const {
  vector<SpinorWaveFunction>    fin,aout;
  vector<SpinorBarWaveFunction> ain,fout;
  SpinorWaveFunction       q(meMomenta()[0],mePartonData()[0],incoming);
  SpinorBarWaveFunction qbar(meMomenta()[1],mePartonData()[1],incoming);
  SpinorBarWaveFunction    f(meMomenta()[2],mePartonData()[2],outgoing);
  SpinorWaveFunction    fbar(meMomenta()[3],mePartonData()[3],outgoing);
  for(unsigned int ix=0;ix<2;++ix) {
    q.reset(ix)   ; fin.push_back(q);
    qbar.reset(ix); ain.push_back(qbar);
    f.reset(ix)   ;fout.push_back(f);
    fbar.reset(ix);aout.push_back(fbar);
  }
  return qqbarME(fin,ain,fout,aout,false);
}

unsigned int MEqq2gZ2ff::orderInAlphaS() const {
  return 0;
}

unsigned int MEqq2gZ2ff::orderInAlphaEW() const {
  return 2;
}

Selector<MEBase::DiagramIndex>
MEqq2gZ2ff::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) {
    if ( diags[i]->id() == -1 ) sel.insert(meInfo()[0], i);
    else if ( diags[i]->id() == -2 ) sel.insert(meInfo()[1], i);
  }
  return sel;
}

Selector<const ColourLines *>
MEqq2gZ2ff::colourGeometries(tcDiagPtr) const {
  static const ColourLines c1("1 -2");
  static const ColourLines c2("1 -2,4 -5");
  Selector<const ColourLines *> sel;
  if(abs(mePartonData()[2]->id())<=6) sel.insert(1.0, &c2);
  else                                sel.insert(1.0, &c1);
  return sel;
}

void MEqq2gZ2ff::persistentOutput(PersistentOStream & os) const {
  os << _minflavour << _maxflavour << _gammaZ << _process
     << _theFFZVertex << _theFFPVertex 
     << _gamma << _z0 << spinCorrelations_;
}

void MEqq2gZ2ff::persistentInput(PersistentIStream & is, int) { 
  is >> _minflavour >> _maxflavour >> _gammaZ >> _process
     >> _theFFZVertex >> _theFFPVertex 
     >> _gamma >> _z0 >> spinCorrelations_; 
}

ClassDescription<MEqq2gZ2ff> MEqq2gZ2ff::initMEqq2gZ2ff;
// Definition of the static class description member.

void MEqq2gZ2ff::Init() {

  static ClassDocumentation<MEqq2gZ2ff> documentation
    ("The MEqq2gZ2ff class implements the matrix element for"
     "q qbar to Standard Model fermions via Z and photon exchange using"
     " helicity amplitude techniques");

  static Parameter<MEqq2gZ2ff,int> interfaceMaxFlavour
    ("MaxFlavour",
     "The maximum incoming quark flavour this matrix element is allowed to handle",
     &MEqq2gZ2ff::_maxflavour, 5, 1, 5,
     false, false, Interface::limited);

  static Parameter<MEqq2gZ2ff,int> interfaceMinFlavour
    ("MinFlavour",
     "The minimum incoming quark flavour this matrix element is allowed to handle",
     &MEqq2gZ2ff::_minflavour, 1, 1, 5,
     false, false, Interface::limited);

  static Switch<MEqq2gZ2ff,unsigned int> interfaceGammaZ
    ("GammaZ",
     "Which terms to include",
     &MEqq2gZ2ff::_gammaZ, 0, false, false);
  static SwitchOption interfaceGammaZAll
    (interfaceGammaZ,
     "All",
     "Include both gamma and Z terms",
     0);
  static SwitchOption interfaceGammaZGamma
    (interfaceGammaZ,
     "Gamma",
     "Only include the photon",
     1);
  static SwitchOption interfaceGammaZZ
    (interfaceGammaZ,
     "Z",
     "Only include the Z",
     2);

  static Switch<MEqq2gZ2ff,int> interfaceProcess
    ("Process",
     "Which process to included",
     &MEqq2gZ2ff::_process, 0, false, false);
  static SwitchOption interfaceProcessAll
    (interfaceProcess,
     "All",
     "Include all SM fermions as outgoing particles",
     0);
  static SwitchOption interfaceProcessQuarks
    (interfaceProcess,
     "Quarks",
     "All include the quarks as outgoing particles",
     1);
  static SwitchOption interfaceProcessLeptons
    (interfaceProcess,
     "Leptons",
     "Only include the leptons as outgoing particles",
     2);
  static SwitchOption interfaceProcessChargedLeptons
    (interfaceProcess,
     "ChargedLeptons",
     "Only include the charged leptons as outgoing particles",
     3);
  static SwitchOption interfaceProcessNeutrinos
    (interfaceProcess,
     "Neutrinos",
     "Only include the neutrinos as outgoing particles",
     4);
  static SwitchOption interfaceProcessElectron
    (interfaceProcess,
     "Electron",
     "Only include e+e- as outgoing particles",
     5);
  static SwitchOption interfaceProcessMuon
    (interfaceProcess,
     "Muon",
     "Only include mu+mu- as outgoing particles",
     6);
  static SwitchOption interfaceProcessTau
    (interfaceProcess,
     "Tau",
     "Only include tau+tau- as outgoing particles",
     7);
  static SwitchOption interfaceProcessNu_e
    (interfaceProcess,
     "Nu_e",
     "Only include nu_e ne_ebar as outgoing particles",
     8);
  static SwitchOption interfaceProcessnu_mu
    (interfaceProcess,
     "Nu_mu",
     "Only include nu_mu nu_mubar as outgoing particles",
     9);
  static SwitchOption interfaceProcessnu_tau
    (interfaceProcess,
     "Nu_tau",
     "Only include nu_tau nu_taubar as outgoing particles",
     10);
  static SwitchOption interfaceProcessDown
    (interfaceProcess,
     "Down",
     "Only include d dbar as outgoing particles",
     11);
  static SwitchOption interfaceProcessUp
    (interfaceProcess,
     "Up",
     "Only include u ubar as outgoing particles",
     12);
  static SwitchOption interfaceProcessStrange
    (interfaceProcess,
     "Strange",
     "Only include s sbar as outgoing particles",
     13);
  static SwitchOption interfaceProcessCharm
    (interfaceProcess,
     "Charm",
     "Only include c cbar as outgoing particles",
     14);
  static SwitchOption interfaceProcessBottom
    (interfaceProcess,
     "Bottom",
     "Only include b bbar as outgoing particles",
     15);
  static SwitchOption interfaceProcessTop
    (interfaceProcess,
     "Top",
     "Only include t tbar as outgoing particles",
     16);

  static Switch<MEqq2gZ2ff,bool> interfaceSpinCorrelations
    ("SpinCorrelations",
     "Which on/off spin correlations in the hard process",
     &MEqq2gZ2ff::spinCorrelations_, true, false, false);
  static SwitchOption interfaceSpinCorrelationsYes
    (interfaceSpinCorrelations,
     "Yes",
     "Switch correlations on",
     true);
  static SwitchOption interfaceSpinCorrelationsNo
    (interfaceSpinCorrelations,
     "No",
     "Switch correlations off",
     false);

}

double MEqq2gZ2ff::qqbarME(vector<SpinorWaveFunction>    & fin ,
			   vector<SpinorBarWaveFunction> & ain ,
			   vector<SpinorBarWaveFunction> & fout,
			   vector<SpinorWaveFunction>    & aout,
			   bool calc) const {
  // scale
  Energy2 mb2(scale());
  // matrix element to be stored
  ProductionMatrixElement menew(PDT::Spin1Half,PDT::Spin1Half,
				PDT::Spin1Half,PDT::Spin1Half);
  // which intermediates to include
  bool gamma=_gammaZ==0||_gammaZ==1;
  bool Z0   =_gammaZ==0||_gammaZ==2;
  // declare the variables we need
  unsigned int ihel1,ihel2,ohel1,ohel2;
  VectorWaveFunction inter[2],test;
  double me[3]={0.,0.,0.};
  Complex diag1,diag2;
  // sum over helicities to get the matrix element
  for(ihel1=0;ihel1<2;++ihel1) {
    for(ihel2=0;ihel2<2;++ihel2) {
      // intermediate for Z
      if(Z0)    inter[0]=_theFFZVertex->evaluate(mb2,1,_z0   ,fin[ihel1],ain[ihel2]);
      // intermediate for photon
      if(gamma) inter[1]=_theFFPVertex->evaluate(mb2,1,_gamma,fin[ihel1],ain[ihel2]);
      for(ohel1=0;ohel1<2;++ohel1) {
	for(ohel2=0;ohel2<2;++ohel2) {
	  // first the Z exchange diagram
	  diag1 = Z0 ?
	    _theFFZVertex->evaluate(mb2,aout[ohel2],fout[ohel1],inter[0]) : 0.;
	  // first the photon exchange diagram
	  diag2 = gamma ?
	    _theFFPVertex->evaluate(mb2,aout[ohel2],fout[ohel1],inter[1]) : 0.;
	  // add up squares of individual terms
	  me[1] += real(diag1*conj(diag1));
	  me[2] += real(diag2*conj(diag2));
	  // the full thing including interference
	  diag1 +=diag2;
	  me[0] += real(diag1*conj(diag1));
	  if(calc) menew(ihel1,ihel2,ohel1,ohel2) = diag1;
	}
      }
    }
  }
  // spin and colour factor
  double colspin=1./12.;
  if(abs(fout[0].id())<=6) colspin*=3.;
  // results
  // factor 12 from 4 helicity and 3 colour
  for(int ix=0;ix<3;++ix){me[ix]*=colspin;}
  DVector save;
  save.push_back(me[1]);
  save.push_back(me[2]);
  meInfo(save);
  if(calc) _me.reset(menew);
  return me[0];
}

void MEqq2gZ2ff::constructVertex(tSubProPtr sub) {
  if(!spinCorrelations_) return;
  // extract the particles in the hard process
  ParticleVector hard;
  hard.push_back(sub->incoming().first);hard.push_back(sub->incoming().second);
  hard.push_back(sub->outgoing()[0]);hard.push_back(sub->outgoing()[1]);
  // order of particles
  unsigned int order[4]={0,1,2,3};
  if(hard[0]->id()<0){order[0]=1;order[1]=0;}
  if(hard[2]->id()<0){order[2]=3;order[3]=2;}
  vector<SpinorWaveFunction>    fin,aout;
  vector<SpinorBarWaveFunction> ain,fout;
  SpinorWaveFunction(   fin ,hard[order[0]],incoming,false,true);
  SpinorBarWaveFunction(ain ,hard[order[1]],incoming,false,true);
  SpinorBarWaveFunction(fout,hard[order[2]],outgoing,true ,true);
  SpinorWaveFunction(   aout,hard[order[3]],outgoing,true ,true);
  qqbarME(fin,ain,fout,aout,true);
  // construct the vertex
  HardVertexPtr hardvertex=new_ptr(HardVertex());
  // set the matrix element for the vertex
  hardvertex->ME(_me);
  // set the pointers and to and from the vertex
  for(unsigned int ix=0;ix<4;++ix) 
    hard[order[ix]]->spinInfo()->productionVertex(hardvertex);
}
