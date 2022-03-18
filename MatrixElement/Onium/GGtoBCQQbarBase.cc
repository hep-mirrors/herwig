// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GGtoBCQQbarBase class.
//

#include "GGtoBCQQbarBase.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;


#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"

void GGtoBCQQbarBase::getDiagrams() const {
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

Selector<MEBase::DiagramIndex>
GGtoBCQQbarBase::diagrams(const DiagramVector & diags) const {
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
GGtoBCQQbarBase::colourGeometries(tcDiagPtr diag) const {
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


void GGtoBCQQbarBase::persistentOutput(PersistentOStream & os) const {
  os << n_ << params_ << state_;
}

void GGtoBCQQbarBase::persistentInput(PersistentIStream & is, int) {
  is >> n_ >> params_ >> state_;
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractClass<GGtoBCQQbarBase,TwoToThree>
  describeHerwigGGtoBCQQbarBase("Herwig::GGtoBCQQbarBase", "GGtoBCQQbarBase.so");

void GGtoBCQQbarBase::Init() {

  static ClassDocumentation<GGtoBCQQbarBase> documentation
    ("The GGtoBCQQbarBase class is the base class fpr g g -> B_c b c processes");
  
  static Parameter<GGtoBCQQbarBase,unsigned int> interfacePrincipalQuantumNumber
    ("PrincipalQuantumNumber",
     "The principle quantum number of the states",
     &GGtoBCQQbarBase::n_, 1, 1, 10,
     false, false, Interface::limited);

  static Reference<GGtoBCQQbarBase,OniumParameters> interfaceParameters
    ("Parameters",
     "Quarkonium parameters",
     &GGtoBCQQbarBase::params_, false, false, true, false, false);

}


void GGtoBCQQbarBase::doinit() {
  TwoToThree::doinit();
  long pid = 100000*(n_-1)+id_;
  state_ = getParticleData(pid);
  if(!state_)
    throw Exception() << "No B_c state with pid = " << pid << "in MEGGBCQQbarBase::doinit()" << Exception::runerror;
  setMassGenerator(dynamic_ptr_cast<GenericMassGeneratorPtr>(state_->massGenerator()));
}
