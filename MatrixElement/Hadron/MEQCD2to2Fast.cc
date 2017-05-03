// -*- C++ -*-
//
// MEQCD2to2Fast.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEQCD2to2Fast class.
//

#include "MEQCD2to2Fast.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/Cuts/Cuts.h"

using namespace Herwig;

IBPtr MEQCD2to2Fast::clone() const {
  return new_ptr(*this);
}

IBPtr MEQCD2to2Fast::fullclone() const {
  return new_ptr(*this);
}

Energy2 MEQCD2to2Fast::scale() const {
  Energy2 s(sHat()),u(uHat()),t(tHat());
  return 2.*s*t*u/(s*s+t*t+u*u);
}

void MEQCD2to2Fast::persistentOutput(PersistentOStream & os) const {
  os << _maxflavour << _process << _strictFlavourScheme;
}

void MEQCD2to2Fast::persistentInput(PersistentIStream & is, int) {
  is >> _maxflavour >> _process >> _strictFlavourScheme;
}

unsigned int MEQCD2to2Fast::orderInAlphaS() const {
  return 2;
}

unsigned int MEQCD2to2Fast::orderInAlphaEW() const {
  return 0;
}

ClassDescription<MEQCD2to2Fast> MEQCD2to2Fast::initMEQCD2to2Fast;
// Definition of the static class description member.

void MEQCD2to2Fast::Init() {

  static ClassDocumentation<MEQCD2to2Fast> documentation
    ("The MEQCD2to2Fast class implements the QCD 2->2 processes in hadron-hadron"
     " collisions");

  static Parameter<MEQCD2to2Fast,unsigned int> interfaceMaximumFlavour
    ("MaximumFlavour",
     "The maximum flavour of the quarks in the process",
     &MEQCD2to2Fast::_maxflavour, 5, 1, 5,
     false, false, Interface::limited);

  static Switch<MEQCD2to2Fast,unsigned int> interfaceProcess
    ("Process",
     "Which subprocesses to include",
     &MEQCD2to2Fast::_process, 0, false, false);
  static SwitchOption interfaceProcessAll
    (interfaceProcess,
     "All",
     "Include all subprocesses",
     0);
  static SwitchOption interfaceProcess1
    (interfaceProcess,
     "gg2gg",
     "Include only gg->gg subprocesses",
     1);
  static SwitchOption interfaceProcess2
    (interfaceProcess,
     "gg2qqbar",
     "Include only gg -> q qbar processes",
     2);
  static SwitchOption interfaceProcessqqbargg
    (interfaceProcess,
     "qqbar2gg",
     "Include only q qbar -> gg processes",
     3);
  static SwitchOption interfaceProcessqgqg
    (interfaceProcess,
     "qg2qg",
     "Include only q g -> q g processes",
     4);
  static SwitchOption interfaceProcessqbargqbarg
    (interfaceProcess,
     "qbarg2qbarg",
     "Include only qbar g -> qbar g processes",
     5);
  static SwitchOption interfaceProcessqqqq
    (interfaceProcess,
     "qq2qq",
     "Include only q q -> q q processes",
     6);
  static SwitchOption interfaceProcessqbarqbarqbarqbar
    (interfaceProcess,
     "qbarqbar2qbarqbar",
     "Include only qbar qbar -> qbar qbar processes",
     7);
  static SwitchOption interfaceProcessqqbarqqbar
    (interfaceProcess,
     "qqbar2qqbar",
     "Include only q qbar -> q qbar processes",
     8);

  static Switch<MEQCD2to2Fast,bool> interfaceStrictFlavourScheme
    ("StrictFlavourScheme",
     "Switch to not include massive initial state quarks.",
     &MEQCD2to2Fast::_strictFlavourScheme, false, false, false);
  static SwitchOption interfaceStrictFlavourSchemeYes
    (interfaceStrictFlavourScheme,
     "Yes",
     "Do not include massive initial states.",
     true);
  static SwitchOption interfaceStrictFlavourSchemeNo
    (interfaceStrictFlavourScheme,
     "No",
     "Allow massive initial states.",
     false);
}

Selector<MEBase::DiagramIndex>
MEQCD2to2Fast::diagrams(const DiagramVector & diags) const {
  // select the diagram, this is easy for us as we have already done it
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) {
    if(diags[i]->id()==-int(_diagram)) sel.insert(1.0, i);
    else sel.insert(0., i);
  }
  return sel;
}

void MEQCD2to2Fast::getDiagrams() const {
  // get the particle data objects
  PDPtr gluon=getParticleData(ParticleID::g);
  vector<PDPtr> quark,antiquark;
  for(int ix=1;ix<=int(_maxflavour);++ix) {
    quark.push_back(    getParticleData( ix));
    antiquark.push_back(getParticleData(-ix));
  }
  // gg-> gg subprocess
  if(_process==0||_process==1) {
    // s-channel
    add(new_ptr((Tree2toNDiagram(2),gluon,gluon, 1, gluon,
		 3,gluon, 3, gluon, -1)));
    // first  t-channel
    add(new_ptr((Tree2toNDiagram(3),gluon,gluon,gluon,
		 1,gluon, 2,gluon,-2)));
    // second t-channel
    add(new_ptr((Tree2toNDiagram(3),gluon,gluon,gluon,
		 2,gluon, 1,gluon,-3)));
  }
  // processes involving one quark line
  for(unsigned int ix=0;ix<_maxflavour;++ix) {
    // gg -> q qbar subprocesses
    if(_process==0||_process==2) {
      // first t-channel
      add(new_ptr((Tree2toNDiagram(3),gluon,antiquark[ix],gluon,
		   1,quark[ix], 2,antiquark[ix],-4)));
      // interchange
      add(new_ptr((Tree2toNDiagram(3),gluon,antiquark[ix],gluon,
		   2,quark[ix], 1,antiquark[ix],-5)));
      // s-channel
      add(new_ptr((Tree2toNDiagram(2),gluon,gluon, 1, gluon,
		   3,quark[ix], 3, antiquark[ix], -6)));
    }
    // q qbar -> g g subprocesses
    if((_process==0||_process==3) && (quark[ix]->hardProcessMass() == ZERO || !_strictFlavourScheme)) {
      // first t-channel
      add(new_ptr((Tree2toNDiagram(3),quark[ix],antiquark[ix],antiquark[ix],
		   1,gluon, 2,gluon,-7)));
      // second t-channel
      add(new_ptr((Tree2toNDiagram(3),quark[ix],antiquark[ix],antiquark[ix],
		   2,gluon, 1,gluon,-8)));
      // s-channel
      add(new_ptr((Tree2toNDiagram(2),quark[ix], antiquark[ix],
		   1, gluon, 3, gluon, 3, gluon,-9)));
    }
    // q g -> q g subprocesses
    if((_process==0||_process==4) && (quark[ix]->hardProcessMass() == ZERO || !_strictFlavourScheme)) {
      // s-channel
      add(new_ptr((Tree2toNDiagram(2),quark[ix], gluon,
		   1, quark[ix], 3, quark[ix], 3, gluon,-10)));
      // quark t-channel
      add(new_ptr((Tree2toNDiagram(3),quark[ix],quark[ix],gluon,
		   2,quark[ix],1,gluon,-11)));
      // gluon t-channel
      add(new_ptr((Tree2toNDiagram(3),quark[ix],gluon,gluon,
		   1,quark[ix],2,gluon,-12)));
    }
    // qbar g -> qbar g subprocesses
    if((_process==0||_process==5) && (antiquark[ix]->hardProcessMass() == ZERO || !_strictFlavourScheme)) {
      // s-channel
      add(new_ptr((Tree2toNDiagram(2),antiquark[ix], gluon,
		   1, antiquark[ix], 3, antiquark[ix], 3, gluon,-13)));
      // quark t-channel
      add(new_ptr((Tree2toNDiagram(3),antiquark[ix],antiquark[ix],gluon,
		   2,antiquark[ix],1,gluon,-14)));
      // gluon t-channel
      add(new_ptr((Tree2toNDiagram(3),antiquark[ix],gluon,gluon,
		   1,antiquark[ix],2,gluon,-15)));
    }
    // processes involving two quark lines
    for(unsigned int iy=ix;iy<_maxflavour;++iy) {
      // q q -> q q subprocesses
      if((_process==0||_process==6) && 
	 ((quark[ix]->hardProcessMass() == ZERO && quark[iy]->hardProcessMass() == ZERO) || !_strictFlavourScheme)) {
	// gluon t-channel
	add(new_ptr((Tree2toNDiagram(3),quark[ix],gluon,quark[iy],
		     1,quark[ix],2,quark[iy],-16)));
	// exchange for identical quarks
	if(ix==iy)
	  add(new_ptr((Tree2toNDiagram(3),quark[ix],gluon,quark[iy],
		       2,quark[ix],1,quark[iy],-17)));
      }
      // qbar qbar -> qbar qbar subprocesses
      if((_process==0||_process==7) && 
	 ((antiquark[ix]->hardProcessMass() == ZERO && antiquark[iy]->hardProcessMass() == ZERO) || !_strictFlavourScheme)) {
	// gluon t-channel
	add(new_ptr((Tree2toNDiagram(3),antiquark[ix],gluon,antiquark[iy],
		     1,antiquark[ix],2,antiquark[iy],-18)));
	// exchange for identical quarks
	if(ix==iy)
	  add(new_ptr((Tree2toNDiagram(3),antiquark[ix],gluon,antiquark[iy],
		       2,antiquark[ix],1,antiquark[iy],-19)));
      }
    }
    for(unsigned int iy=0;iy<_maxflavour;++iy) {
      // q qbar -> q qbar
      if(_process==0||_process==8) {
	// gluon s-channel
	if ( quark[ix]->hardProcessMass() == ZERO || !_strictFlavourScheme )
	  add(new_ptr((Tree2toNDiagram(2),quark[ix], antiquark[ix],
		       1, gluon, 3, quark[iy], 3, antiquark[iy],-20)));
	// gluon t-channel
	if ( (quark[ix]->hardProcessMass() == ZERO && antiquark[iy]->hardProcessMass() == ZERO) || 
	     !_strictFlavourScheme )
	  add(new_ptr((Tree2toNDiagram(3),quark[ix],gluon,antiquark[iy],
		       1,quark[ix],2,antiquark[iy],-21)));
      }
    }
  }
}

Selector<const ColourLines *>
MEQCD2to2Fast::colourGeometries(tcDiagPtr diag) const {
  // colour lines for gg to gg
  static const ColourLines cgggg[12]={ColourLines("1 -2, -1 -3 -5, 5 -4, 2 3 4"),// A_2 s
				      ColourLines("-1 2, 1 3 5, -5 4, -2 -3 -4"),// A_1 s
				      ColourLines("1 5, -1 -2 3, -3 -4, -5 2 4"),// A_1 u
				      ColourLines("-1 -5, 1 2 -3, 3 4, 5 -2 -4"),// A_2 u
				      ColourLines("1 -2, -1 -3 -4, 4 -5, 2 3 5"),// B_2 s
				      ColourLines("-1 2, 1 3 4, -4 5, -2 -3 -5"),// B_1 s
				      ColourLines("1 4, -1 -2 3, -3 -5, -4 2 5"),// B_1 t
				      ColourLines("-1 -4, 1 2 -3, 3 5, 4 -2 -5"),// B_2 t
				      ColourLines("1 4, -1 -2 -5, 3 5, -3 2 -4"),// C_1 t
				      ColourLines("-1 -4, 1 2 5, -3 -5, 3 -2 4"),// C_2 t
				      ColourLines("1 5, -1 -2 -4, 3 4, -3 2 -5"),// C_1 u
				      ColourLines("-1 -5, 1 2 4, -3 -4, 3 -2 5") // C_2 u
  };
  // colour lines for gg to q qbar
  static const ColourLines cggqq[4]={ColourLines("1  4, -1 -2 3, -3 -5"),
				     ColourLines("3  4, -3 -2 1, -1 -5"),
				     ColourLines("2 -1,  1  3 4, -2 -3 -5"),
				     ColourLines("1 -2, -1 -3 -5, 2 3 4")};
  // colour lines for q qbar to gg
  static const ColourLines cqqgg[4]={ColourLines("1 4, -4 -2 5, -3 -5"),
				     ColourLines("1 5, -3 -4, 4 -2 -5"),
				     ColourLines("1 3 4, -4 5, -2 -3 -5"),
				     ColourLines("1 3 5, -5 4, -2 -3 -4")};
  // colour lines for q g to q g
  static const ColourLines cqgqg[4]={ColourLines("1 -2, 2 3 5, 4 -5"),
				     ColourLines("1 5, 3 4,-3 2 -5 "),
				     ColourLines("1 2 -3, 3 5, -5 -2 4"),
				     ColourLines("1 -2 5,3 2 4,-3 -5")};
  // colour lines for qbar g -> qbar g
  static const ColourLines cqbgqbg[4]={ColourLines("-1 2, -2 -3 -5, -4 5"),
				       ColourLines("-1 -5, -3 -4, 3 -2 5"),
				       ColourLines("-1 -2 3, -3 -5, 5 2 -4"),
				       ColourLines("-1 2 -5,-3 -2 -4, 3 5")};
  // colour lines for q q -> q q 
  static const ColourLines cqqqq[2]={ColourLines("1 2 5,3 -2 4"),
				     ColourLines("1 2 4,3 -2 5")};
  // colour lines for qbar qbar -> qbar qbar
  static const ColourLines cqbqbqbqb[2]={ColourLines("-1 -2 -5,-3 2 -4"),
					 ColourLines("-1 -2 -4,-3 2 -5")};
  // colour lines for q qbar -> q qbar
  static const ColourLines cqqbqqb[2]={ColourLines("1 3 4,-2 -3 -5"),
				       ColourLines("1 2 -3,4 -2 -5")};
  // select the colour flow (as all ready picked just insert answer)
  Selector<const ColourLines *> sel;
  switch(abs(diag->id())) {
    // gg -> gg subprocess
  case 1:
    if(_flow==1) {
      sel.insert(0.5, &cgggg[0]);
      sel.insert(0.5, &cgggg[1]);
    }
    else {
      sel.insert(0.5, &cgggg[4]);
      sel.insert(0.5, &cgggg[5]);
    }
    break;
  case 2: 
    if(_flow==2) {
      sel.insert(0.5, &cgggg[6]);
      sel.insert(0.5, &cgggg[7]);
    }
    else {
      sel.insert(0.5, &cgggg[8]);
      sel.insert(0.5, &cgggg[9]);
    }
    break;
  case 3:
    if(_flow==1) {
      sel.insert(0.5, &cgggg[2]);
      sel.insert(0.5, &cgggg[3]);
    }
    else {
      sel.insert(0.5, &cgggg[10]);
      sel.insert(0.5, &cgggg[11]);
    }
    break;
    // gg -> q qbar subprocess
  case 4: case 5:
    sel.insert(1.0, &cggqq[abs(diag->id())-4]);
    break;
  case 6:
    sel.insert(1.0, &cggqq[1+_flow]);
    break;
    // q qbar -> gg subprocess
  case 7: case 8:
    sel.insert(1.0, &cqqgg[abs(diag->id())-7]);
    break;
  case 9:
    sel.insert(1.0, &cqqgg[1+_flow]);
    break;
    // q g -> q g subprocess
  case 10: case 11:
    sel.insert(1.0, &cqgqg[abs(diag->id())-10]);
    break;
  case 12:
    sel.insert(1.0, &cqgqg[1+_flow]);
    break;
    // q g -> q g subprocess
  case 13: case 14:
    sel.insert(1.0, &cqbgqbg[abs(diag->id())-13]);
    break;
  case 15:
    sel.insert(1.0, &cqbgqbg[1+_flow]);
    break;
    // q q -> q q subprocess
  case 16: case 17:
    sel.insert(1.0, &cqqqq[abs(diag->id())-16]);
    break;
    // qbar qbar -> qbar qbar subprocess
  case 18: case 19:
    sel.insert(1.0, &cqbqbqbqb[abs(diag->id())-18]);
    break;
    // q qbar -> q qbar subprocess
  case 20: case 21:
    sel.insert(1.0, &cqqbqqb[abs(diag->id())-20]);
    break;
  }
  return sel;
}

double MEQCD2to2Fast::me2() const {
  // total matrix element
  double me(0.);
  // gg initiated processes
  if(mePartonData()[0]->id()==ParticleID::g&&mePartonData()[1]->id()==ParticleID::g) {
    // gg -> gg
    if(mePartonData()[2]->id()==ParticleID::g) me = gg2ggME();
    // gg -> q qbar
    else  me=gg2qqbarME();
  }
  // quark first processes
  else if(mePartonData()[0]->id()>0) {
    // q g -> q g
    if(mePartonData()[1]->id()==ParticleID::g) me = qg2qgME();
    else if(mePartonData()[1]->id()<0) {
      // q qbar initiated processes( q qbar -> gg)
      if(mePartonData()[2]->id()==ParticleID::g) me = qqbar2ggME();
      // q qbar to q qbar 
      else me = qqbar2qqbarME();
    }
    // q q -> q q 
    else if(mePartonData()[1]->id()>0) me = qq2qqME();
  }
  // antiquark first processes
  else if(mePartonData()[0]->id()<0) {
    // qbar g -> qbar g
    if(mePartonData()[1]->id()==ParticleID::g) me = qbarg2qbargME();
    // qbar qbar -> qbar qbar
    else if(mePartonData()[1]->id()<0)         me = qbarqbar2qbarqbarME();
  }
  else {
    throw Exception() << "Unknown process in MEQCD2to2Fast::me2()" 
		      << Exception::abortnow;
  }
  // multpliy by alpha_S^2 and return  answer
  double alphas(4.*Constants::pi*SM().alphaS(scale()));
  return me*sqr(alphas);
}

void MEQCD2to2Fast::doinit() {
  if ( _strictFlavourScheme )
    massOption(vector<unsigned int>(2,1));
  HwMEBase::doinit();
}
