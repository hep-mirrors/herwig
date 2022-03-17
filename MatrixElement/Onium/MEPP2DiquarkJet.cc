// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2DiquarkJet class.
//

#include "MEPP2DiquarkJet.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;


void MEPP2DiquarkJet::getDiagrams() const {
  // particle data objects of the quarks in the diquark
  tcPDPtr q1 = getParticleData( state_->id()/1000);
  tcPDPtr q2 = getParticleData((state_->id()%1000)/100);
  tcPDPtr g  = getParticleData(ParticleID::g);
  // first  g q1 -> diquark q2bar
  add(new_ptr((Tree2toNDiagram(2), g, q1, 1, q1 , 3, state_, 3, q2->CC(), -1)));
  add(new_ptr((Tree2toNDiagram(3), g, q1, q1, 2, state_, 1, q2->CC(),     -2)));
  // second g q2 -> diquark q1bar
  if(q1!=q2) {
    add(new_ptr((Tree2toNDiagram(2), g, q2, 1, q2 , 3, state_, 3, q1->CC(), -1)));
    add(new_ptr((Tree2toNDiagram(3), g, q2, q2, 2, state_, 1, q1->CC(),     -2)));
  }
  // q1 q2 -< diquark g (swap automatically included)
  add(new_ptr((Tree2toNDiagram(3), q1, q2, q2, 1, state_, 2, g,     -3)));
  add(new_ptr((Tree2toNDiagram(3), q1, q1, q2, 2, state_, 1, g,     -4)));
  // CC processes
  // first  g q1bar -> antidiquark q2
  add(new_ptr((Tree2toNDiagram(2), g, q1->CC(), 1, q1->CC() , 3, state_->CC(), 3, q2, -5)));
  add(new_ptr((Tree2toNDiagram(3), g, q1->CC(), q1->CC(), 2, state_->CC(), 1, q2,     -6)));
  // second g q2bar -> antidiquark q1
  if(q1!=q2) {
    add(new_ptr((Tree2toNDiagram(2), g, q2->CC(), 1, q2->CC() , 3, state_->CC(), 3, q1, -5)));
    add(new_ptr((Tree2toNDiagram(3), g, q2->CC(), q2->CC(), 2, state_->CC(), 1, q1,     -6)));
  }
  // q1 q2 -< diquark g (swap automatically included)
  add(new_ptr((Tree2toNDiagram(3), q1->CC(), q2->CC(), q2->CC(), 1, state_->CC(), 2, g,     -7)));
  add(new_ptr((Tree2toNDiagram(3), q1->CC(), q2->CC(), q1->CC(), 2, state_->CC(), 1, g,     -8)));
}

Energy2 MEPP2DiquarkJet::scale() const {
  return sHat();
}

Selector<MEBase::DiagramIndex>
MEPP2DiquarkJet::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) 
    if ( abs(diags[i]->id())%2==1 )
      sel.insert(meInfo()[0], i);
    else
      sel.insert(meInfo()[1], i);
  return sel;
}

Selector<const ColourLines *>
MEPP2DiquarkJet::colourGeometries(tcDiagPtr diag) const {
  static ColourLines c[8] = {ColourLines("1 3 -6, 2 -1,-4 -6, -5 -6"),
			     ColourLines("-1 -5, 1 2 -6,3 -6, -4 -6"),
			     ColourLines("1 -6, 3 5, -4 2 -6, -5 -6"),
			     ColourLines("1 5, 3 -6, -4 -6, -5 2 -6"),
			     ColourLines("-1 -3 6, -2 1,4 6, 5 6"),
			     ColourLines("1 5, -1 -2 6,-3 6, 4 6"),
			     ColourLines("-1 6, -3 -5, 4 -2 6, 5 6"),
			     ColourLines("-1 -5, -3 6, 4 6, 5 -2 6")};
  Selector<const ColourLines *> sel;
  sel.insert(1.0, &c[abs(diag->id())-1]);
  return sel;
}

void MEPP2DiquarkJet::persistentOutput(PersistentOStream & os) const {
  os << state_;
}

void MEPP2DiquarkJet::persistentInput(PersistentIStream & is, int) {
  is >> state_;
}

void MEPP2DiquarkJet::doinit() {
  HwMEBase::doinit();
  long pid = id_;
  state_ = getParticleData(pid);
  if(!state_)
    throw Exception() << "No diquark state with pid = " << pid << "in MEPP2DiquarkJet::doinit()" << Exception::runerror;
  setMassGenerator(dynamic_ptr_cast<GenericMassGeneratorPtr>(state_->massGenerator()));
}

// The following static variable is needed for the type
// description system in ThePEG. 
DescribeAbstractClass<MEPP2DiquarkJet,MassiveIncoming>
describeHerwigMEPP2DiquarkJet("Herwig::MEPP2DiquarkJet", "HwOniumParameters.so HwMEHadronOnium.so");

void MEPP2DiquarkJet::Init() {

  static ClassDocumentation<MEPP2DiquarkJet> documentation
    ("Base class for 2->2 diquark processes");

  static Reference<MEPP2DiquarkJet,OniumParameters> interfaceParameters
    ("Parameters",
     "Quarkonium parameters",
     &MEPP2DiquarkJet::params_, false, false, true, false, false);

}

