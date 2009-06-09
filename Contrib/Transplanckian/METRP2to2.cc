// -*- C++ -*-
//
// METRP2to2.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the METRP2to2 class.
//

#include "METRP2to2.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/Cuts/Cuts.h"
#include <fstream>

using namespace Herwig;


void METRP2to2::persistentOutput(PersistentOStream & os) const {
  os << _maxflavour << _process << _ndim << _planckmass;
}

void METRP2to2::persistentInput(PersistentIStream & is, int) {
  is >> _maxflavour >> _process >> _ndim >> _planckmass;
}


unsigned int METRP2to2::orderInAlphaS() const {
  return 0;
}

unsigned int METRP2to2::orderInAlphaEW() const {
  return 0;
}

Energy2 METRP2to2::scale() const {
  Energy2 s(sHat()),t(tHat());
  double bcsq = 0;
  double qsq = -t;
  bcsq = pow(bccalc(s),2);
  if(qsq > 1/bcsq) { return (1/bcsq); } else { return qsq; }
}


ClassDescription<METRP2to2> METRP2to2::initMETRP2to2;
// Definition of the static class description member.

void METRP2to2::Init() {

  static ClassDocumentation<METRP2to2> documentation
    ("The METRP2to2 class implements the transplanckian 2->2 processes in hadron-hadron"
     " collisions");

  static Parameter<METRP2to2,unsigned int> interfaceMaximumFlavour
    ("MaximumFlavour",
     "The maximum flavour of the quarks in the process",
     &METRP2to2::_maxflavour, 5, 1, 5,
     false, false, Interface::limited);

 static Parameter<METRP2to2, double> interfacePlanckMass
    ("PlanckMass",
     "The Planck Mass",
     &METRP2to2::_planckmass, 2000.0, 200.0,200000.0,
     false, false, Interface::limited);


  
  static Parameter<METRP2to2, unsigned int> interfaceNumberExtraDimensions
    ("NumberExtraDimensions",
     "The number of extra dimensions to consider",
     &METRP2to2::_ndim, 3, 2, 6,
     false, false, Interface::limited);


  static Switch<METRP2to2,unsigned int> interfaceProcess
    ("Process",
     "Which subprocesses to include",
     &METRP2to2::_process, 0, false, false);
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
}

Selector<MEBase::DiagramIndex>
METRP2to2::diagrams(const DiagramVector & diags) const {
  // select the diagram, this is easy for us as we have already done it
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) {
    if(diags[i]->id()==-int(_diagram)) sel.insert(1.0, i);
    else sel.insert(0., i);
  }
  return sel;
}

void METRP2to2::getDiagrams() const {
  // get the particle data objects
  PDPtr gluon=getParticleData(ParticleID::g);
  PDPtr trpon = getParticleData(991);

  vector<PDPtr> quark,antiquark;
  for(int ix=1;ix<=int(_maxflavour);++ix) {
    quark.push_back(    getParticleData( ix));
    antiquark.push_back(getParticleData(-ix));
  }
  // gg-> gg subprocess
  if(_process==0||_process==1) {
   
    add(new_ptr((Tree2toNDiagram(3),gluon,trpon,gluon,
		 1,gluon, 2,gluon,-2)));
   
  }
  // processes involving one quark line
    for(unsigned int ix=0;ix<_maxflavour;++ix) {
    
    // q g -> q g subprocesses
    if(_process==0||_process==4) {
     
     
    
      add(new_ptr((Tree2toNDiagram(3),quark[ix],trpon,gluon,
		   1,quark[ix],2,gluon,-12)));

    }

    // qbar g -> qbar g subprocesses
    if(_process==0||_process==5) {
   
   
      add(new_ptr((Tree2toNDiagram(3),antiquark[ix],trpon,gluon,
		   1,antiquark[ix],2,gluon,-15)));
     
    }
    // processes involving two quark lines
    for(unsigned int iy=0;iy<_maxflavour;++iy) {
      // q q -> q q subprocesses
      if(_process==0||_process==6) {
	// gluon t-channel
	add(new_ptr((Tree2toNDiagram(3),quark[ix],trpon,quark[iy],
		     1,quark[ix],2,quark[iy],-16)));
	// exchange for identical quarks
	if(ix==iy)
	  add(new_ptr((Tree2toNDiagram(3),quark[ix],trpon,quark[iy],
		       2,quark[ix],1,quark[iy],-17)));
      }
      // qbar qbar -> qbar qbar subprocesses
      if(_process==0||_process==7) {
	// gluon t-channel
		add(new_ptr((Tree2toNDiagram(3),antiquark[ix],trpon,antiquark[iy],
		     1,antiquark[ix],2,antiquark[iy],-18)));
	// exchange for identical quarks
	if(ix==iy)
	  add(new_ptr((Tree2toNDiagram(3),antiquark[ix],trpon,antiquark[iy],
		       2,antiquark[ix],1,antiquark[iy],-19)));
      }
      // q qbar -> q qbar
      if(_process==0||_process==8) {
	add(new_ptr((Tree2toNDiagram(3),quark[ix],trpon,antiquark[iy],
		     1,quark[ix],2,antiquark[iy],-21)));
      }
    }
  }
}

Selector<const ColourLines *>
METRP2to2::colourGeometries(tcDiagPtr diag) const {
  // colour lines for gg to gg
  static const ColourLines cgggg[2]={
    ColourLines("1 4, -1 -4, 3 5, -3 -5"),
    ColourLines("1 4, -1 -4, 3 5, -3 -5")};
  // colour lines for q g to q g
  static const ColourLines cqgqg[4]={ColourLines("1 4, 3 5, -3 -5"),ColourLines("1 4, 3 5, -3 -5")};
  // colour lines for qbar g -> qbar g
  static const ColourLines cqbgqbg[4]={ColourLines("-1 -4, -3 -5, 3 5"),
				       ColourLines("-1 -4, -3 -5, 3 5")};

  // colour lines for q q -> q q 
  static const ColourLines cqqqq[2]={ColourLines("1 4,3 5"), ColourLines("1 4,3 5")};

  // colour lines for qbar qbar -> qbar qbar
  static const ColourLines cqbqbqbqb[2]={ColourLines("-1 -4,-3 -5"),
					 ColourLines("-1 -4,-3 -5")};
  // colour lines for q qbar -> q qbar
  static const ColourLines cqqbqqb[2]={ColourLines("1 4,-3 -5"),ColourLines("1 4,-3 -5")};
  // select the colour flow (as all ready picked just insert answer)
  Selector<const ColourLines *> sel;
  switch(abs(diag->id())) {
    //gg -> gg 
  case 2: 
    sel.insert(0.5, &cgggg[0]);
    sel.insert(0.5, &cgggg[0]);
    break;
    // q g -> q g subprocess
  case 12:
    sel.insert(1.0, &cqgqg[0]);
    break;
    // qbar g -> qbar g subprocess
  case 15:
    sel.insert(1.0, &cqbgqbg[0]);
    break;
    // q q -> q q subprocess
  case 16: case 17:
    sel.insert(1.0, &cqqqq[0]);
    break;
    // qbar qbar -> qbar qbar subprocess
  case 18: case 19:
    sel.insert(1.0, &cqbqbqbqb[0]);
    break;
    // q qbar -> q qbar subprocess
  case 21:
    sel.insert(1.0, &cqqbqqb[0]);
    break;
  }
  return sel;
}


double METRP2to2::me2() const {
  double me(0.);
  ofstream outd, outme2;
  // outd.open("qdist.dat", fstream::app);
  //outme2.open("me2.dat", fstream::app);
  me = ME();
  double pifac = 3.1415926536;
  double fac = 16 * sqr(pifac);
  me *= fac;
  // cout << "ecm = " << sqrt(sHat()) << " q = " << sqrt(-tHat()) << " me = " << me << " costh = " << 1 - 2 * fabs(tHat()) / sHat() << endl; 
  //outd << sqrt(-tHat()) << endl;
  // outme2 << me << endl;
  // outme2.close();
  //outd.close();
  return me;
  
}
