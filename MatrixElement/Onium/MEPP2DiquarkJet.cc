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
  // Here is an example on how to specify diagrams.

  tcPDPtr g = getParticleData(ParticleID::g);
  for ( int i = 1; i <= 5; ++i ) {
    tcPDPtr q = getParticleData(i);
    tcPDPtr qb = q->CC();

    // For each flavour we add:
    add(new_ptr((Tree2toNDiagram(3), q, qb, qb, 1, g, 2, g, -1)));
    // t-channel q + qbar -> g + g
    add(new_ptr((Tree2toNDiagram(3), q, qb, qb, 2, g, 1, g, -2)));
    // u-channel q + qbar -> g + g
    add(new_ptr((Tree2toNDiagram(2), q, qb, 1, g , 3, g, 3, g, -3)));
    // s-channel q + qbar -> g + g
  }
}

Energy2 MEPP2DiquarkJet::scale() const {
  return sHat();
}

Selector<MEBase::DiagramIndex>
MEPP2DiquarkJet::diagrams(const DiagramVector & diags) const {
  // This example corresponds to the diagrams specified in the example
  // in the getDiagrams() function.

  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) 
    if ( diags[i]->id() == -1 ) sel.insert(1.0, i);
    else if ( diags[i]->id() == -2 )  sel.insert(1.0, i);
    else if ( diags[i]->id() == -3 )  sel.insert(1.0, i);
  // You probably do not want equal weights here...
  return sel;

  // If there is only one possible diagram you can override the
  // MEBase::diagram function instead.

}

Selector<const ColourLines *>
MEPP2DiquarkJet::colourGeometries(tcDiagPtr diag) const {
  // This example corresponds to the diagrams specified in the example
  // in the getDiagrams() function.

  static ColourLines ctST("1 4, -4 -2 5, -5 -3");
  static ColourLines ctSU("1 5, -5 -2 4, -4 -3");

  Selector<const ColourLines *> sel;
  if ( diag->id() == -1 || diag->id() == -3 )
    sel.insert(1.0, &ctST);
  else
    sel.insert(1.0, &ctSU);
  return sel;

  // If there is only one possible colour geometry you can override the
  // MEBase::selectColourGeometry function instead.

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

