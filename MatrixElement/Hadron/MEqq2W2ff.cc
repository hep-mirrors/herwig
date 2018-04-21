// -*- C++ -*-
//
// MEqq2W2ff.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEqq2W2ff class.
//

#include "MEqq2W2ff.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Herwig/MatrixElement/HardVertex.h"
#include <cassert>

using namespace Herwig;

MEqq2W2ff::MEqq2W2ff() : _maxflavour(5), _plusminus(0), _process(0) {
  massOption(vector<unsigned int>(2,1));
}

void MEqq2W2ff::doinit() {
  DrellYanBase::doinit();
  _wp=getParticleData(ThePEG::ParticleID::Wplus);
  _wm=getParticleData(ThePEG::ParticleID::Wminus);
  // cast the SM pointer to the Herwig SM pointer
  tcHwSMPtr hwsm=ThePEG::dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  // do the initialisation
  if(hwsm) 
    _theFFWVertex = hwsm->vertexFFW();
  else
    throw InitException() << "Must be the Herwig StandardModel class in "
			  << "MEqq2W2ff::doinit" << Exception::abortnow;
}

void MEqq2W2ff::getDiagrams() const {
  // which intgermediates to include
  bool wplus =_plusminus==0||_plusminus==1;
  bool wminus=_plusminus==0||_plusminus==2;
  // possible incoming and outgoing particles
  typedef std::vector<pair<long,long> > Pairvector;
  // possible parents
  Pairvector parentpair;
  parentpair.reserve(6);
  // don't even think of putting 'break' in here!
  switch(_maxflavour) {
  case 5:
    parentpair.push_back(make_pair(ParticleID::b, ParticleID::cbar));
    parentpair.push_back(make_pair(ParticleID::b, ParticleID::ubar));
    [[fallthrough]];
  case 4:
    parentpair.push_back(make_pair(ParticleID::s, ParticleID::cbar));
    parentpair.push_back(make_pair(ParticleID::d, ParticleID::cbar));
    [[fallthrough]];
  case 3:
    parentpair.push_back(make_pair(ParticleID::s, ParticleID::ubar));
    [[fallthrough]];
  case 2:
    parentpair.push_back(make_pair(ParticleID::d, ParticleID::ubar));
    [[fallthrough]];
  default:
    ;
  }
  // possible children
  Pairvector childpair;
  childpair.reserve(9);
  childpair.push_back(make_pair(ParticleID::eminus,   ParticleID::nu_ebar));
  childpair.push_back(make_pair(ParticleID::muminus,  ParticleID::nu_mubar));
  childpair.push_back(make_pair(ParticleID::tauminus, ParticleID::nu_taubar));
  childpair.push_back(make_pair(ParticleID::d, ParticleID::ubar));
  childpair.push_back(make_pair(ParticleID::s, ParticleID::ubar));
  childpair.push_back(make_pair(ParticleID::b, ParticleID::ubar));
  childpair.push_back(make_pair(ParticleID::d, ParticleID::cbar));
  childpair.push_back(make_pair(ParticleID::s, ParticleID::cbar));
  childpair.push_back(make_pair(ParticleID::b, ParticleID::cbar));
  // loop over the children
  bool lepton,quark;
  Pairvector::const_iterator child = childpair.begin();
  for (; child != childpair.end(); ++child) {
    assert(child->first > 0 && child->second < 0);
    // allowed leptonic decay
    lepton = child->first > 10 && (_process==0 
				   || _process==2 
				   || (child->first - 5)/2 == int(_process));
    // allowed quark decay
    quark  = child->first < 10 && (_process==0 
				   || _process==1 
				   || (child->second == -2 && (child->first + 11)/2 == int(_process)) 
				   || (child->second == -4 && (child->first + 17)/2 == int(_process))
				   );
    // if decay not allowed skip
    if(!(quark || lepton)) continue;
    // decay products
    tcPDPtr lNeg1 = getParticleData(child->first);
    tcPDPtr lNeg2 = getParticleData(child->second);
    tcPDPtr lPos1 = lNeg2->CC();
    tcPDPtr lPos2 = lNeg1->CC();
    Pairvector::const_iterator parent = parentpair.begin();
    for (; parent != parentpair.end(); ++parent) {
      // parents
      tcPDPtr qNeg1 = getParticleData(parent->first);
      tcPDPtr qNeg2 = getParticleData(parent->second);
      tcPDPtr qPos1 = qNeg2->CC();
      tcPDPtr qPos2 = qNeg1->CC();
      // diagrams
      if(wminus) add(new_ptr((Tree2toNDiagram(2), 
			      qNeg1, qNeg2, 
			      1, _wm, 
			      3, lNeg1, 3, lNeg2, -1)));

      if(wplus)  add(new_ptr((Tree2toNDiagram(2), 
			      qPos1, qPos2, 
			      1, _wp, 
			      3, lPos1, 3, lPos2, -2)));
    }
  }
}

Energy2 MEqq2W2ff::scale() const {
  return sHat();
}

double MEqq2W2ff::me2() const {
  vector<SpinorWaveFunction>    fin,aout;
  vector<SpinorBarWaveFunction> ain,fout;
  SpinorWaveFunction    q   (meMomenta()[0],mePartonData()[0],incoming);
  SpinorBarWaveFunction qbar(meMomenta()[1],mePartonData()[1],incoming);
  SpinorBarWaveFunction f   (meMomenta()[2],mePartonData()[2],outgoing);
  SpinorWaveFunction    fbar(meMomenta()[3],mePartonData()[3],outgoing);
  for(unsigned int ix=0;ix<2;++ix) {
    q.reset(ix)   ; fin.push_back(q);
    qbar.reset(ix); ain.push_back(qbar);
    f.reset(ix)   ;fout.push_back(f);
    fbar.reset(ix);aout.push_back(fbar);
  }
  return qqbarME(fin,ain,fout,aout,false);
}

unsigned int MEqq2W2ff::orderInAlphaS() const {
  return 0;
}

unsigned int MEqq2W2ff::orderInAlphaEW() const {
  return 2;
}

Selector<MEBase::DiagramIndex>
MEqq2W2ff::diagrams(const DiagramVector &) const {
  Selector<DiagramIndex> sel;
  sel.insert(1.0, 0);
  return sel;
}

Selector<const ColourLines *>
MEqq2W2ff::colourGeometries(tcDiagPtr) const {
  static const ColourLines c1("1 -2");
  static const ColourLines c2("1 -2,4 -5");
  Selector<const ColourLines *> sel;
  if(abs(mePartonData()[2]->id())<=6) sel.insert(1.0, &c2);
  else                                sel.insert(1.0, &c1);
  return sel;
}

void MEqq2W2ff::persistentOutput(PersistentOStream & os) const {
  os << _maxflavour << _plusminus << _process << _theFFWVertex << _wp << _wm;
}

void MEqq2W2ff::persistentInput(PersistentIStream & is, int) {
  is >> _maxflavour >> _plusminus >> _process >> _theFFWVertex >> _wp >> _wm;
}

ClassDescription<MEqq2W2ff> MEqq2W2ff::initMEqq2W2ff;
// Definition of the static class description member.

void MEqq2W2ff::Init() {

    static ClassDocumentation<MEqq2W2ff> documentation
    ("The MEqq2W2ff class implements the matrix element for"
     "q qbar to Standard Model fermions via W exchange using helicity amplitude"
     "techniques");

  static Parameter<MEqq2W2ff,unsigned int> interfaceMaxFlavour
    ( "MaxFlavour",
      "The heaviest incoming quark flavour this matrix element is allowed to handle "
      "(if applicable).",
      &MEqq2W2ff::_maxflavour, 5, 0, 5, false, false, true);

  static Switch<MEqq2W2ff,unsigned int> interfacePlusMinus
    ("Wcharge",
     "Which intermediate W bosons to include",
     &MEqq2W2ff::_plusminus, 0, false, false);
  static SwitchOption interfacePlusMinusAll
    (interfacePlusMinus,
     "Both",
     "Include W+ and W-",
     0);
  static SwitchOption interfacePlusMinusPlus
    (interfacePlusMinus,
     "Plus",
     "Only include W+",
     1);
  static SwitchOption interfacePlusMinusMinus
    (interfacePlusMinus,
     "Minus",
     "Only include W-",
     2);

  static Switch<MEqq2W2ff,unsigned int> interfaceProcess
    ("Process",
     "Which processes to include",
     &MEqq2W2ff::_process, 0, false, false);
  static SwitchOption interfaceProcessAll
    (interfaceProcess,
     "All",
     "Include all SM fermions as outgoing particles",
     0);
  static SwitchOption interfaceProcessQuarks
    (interfaceProcess,
     "Quarks",
     "Only include outgoing quarks",
     1);
  static SwitchOption interfaceProcessLeptons
    (interfaceProcess,
     "Leptons",
     "All include outgoing leptons",
     2);
  static SwitchOption interfaceProcessElectron
    (interfaceProcess,
     "Electron",
     "Only include outgoing e nu_e",
     3);
  static SwitchOption interfaceProcessMuon
    (interfaceProcess,
     "Muon",
     "Only include outgoing mu nu_mu",
     4);
  static SwitchOption interfaceProcessTau
    (interfaceProcess,
     "Tau",
     "Only include outgoing tauu nu_tau",
     5);
  static SwitchOption interfaceProcessUpDown
    (interfaceProcess,
     "UpDown",
     "Only include outgoing u dbar/ d ubar",
     6);
  static SwitchOption interfaceProcessUpStrange
    (interfaceProcess,
     "UpStrange",
     "Only include outgoing u sbar/ s ubar",
     7);
  static SwitchOption interfaceProcessUpBottom
    (interfaceProcess,
     "UpBottom",
     "Only include outgoing u bbar/ b ubar",
     8);
  static SwitchOption interfaceProcessCharmDown
    (interfaceProcess,
     "CharmDown",
     "Only include outgoing c dbar/ d cbar",
     9);
  static SwitchOption interfaceProcessCharmStrange
    (interfaceProcess,
     "CharmStrange",
     "Only include outgoing c sbar/ s cbar",
     10);
  static SwitchOption interfaceProcessCharmBottom
    (interfaceProcess,
     "CharmBottom",
     "Only include outgoing c bbar/ b cbar",
     11);

}

double MEqq2W2ff::qqbarME(vector<SpinorWaveFunction>    & fin ,
			  vector<SpinorBarWaveFunction> & ain ,
			  vector<SpinorBarWaveFunction> & fout,
			  vector<SpinorWaveFunction>    & aout,
			  bool calc) const {
  // matrix element to be stored
  ProductionMatrixElement newme(PDT::Spin1Half,PDT::Spin1Half,
				PDT::Spin1Half,PDT::Spin1Half);
  // positive or negative W boson
  bool positive = mePartonData()[0]->iCharge() 
    + mePartonData()[1]->iCharge() > 0;
  unsigned int ihel1,ihel2,ohel1,ohel2;
  // sum over helicities to get the matrix element
  double me = 0.;
  VectorWaveFunction inter;
  for(ihel1=0;ihel1<2;++ihel1) {
    for(ihel2=0;ihel2<2;++ihel2) {
      for(ohel1=0;ohel1<2;++ohel1) {
	for(ohel2=0;ohel2<2;++ohel2) {
	  Complex diag = 0.0;
	  if (positive) {
	    // the Wp exchange
	    inter = _theFFWVertex->evaluate(scale(),1,_wp,fin[ihel1],ain[ihel2]);
	    diag = _theFFWVertex->evaluate(scale(),aout[ohel2],fout[ohel1],inter);
	  } 
	  else {
	    // the Wm exchange
	    inter = _theFFWVertex->evaluate(scale(),1,_wm,fin[ihel1],ain[ihel2]);
	    diag = _theFFWVertex->evaluate(scale(),aout[ohel2],fout[ohel1],inter);
	  }
	  // sum over helicities
	  me += real(diag*conj(diag));
	  if(calc) newme(ihel1,ihel2,ohel1,ohel2) = diag;
	}
      }
    }
  }
  // results
  // spin and colour factor
  double colspin=1./12.;
  if(abs(fout[0].id())<=6) colspin*=3.;
  if(calc) _me.reset(newme);
  return me*colspin;
}

void MEqq2W2ff::constructVertex(tSubProPtr sub) {
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
