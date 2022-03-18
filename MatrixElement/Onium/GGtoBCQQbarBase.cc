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
  tcPDPtr g = getParticleData(ParticleID::g);
  tcPDPtr b = getParticleData(ParticleID::b);
  tcPDPtr c = getParticleData(ParticleID::c);
  add(new_ptr((Tree2toNDiagram(2), g, g, 1, state_      , 1, c->CC(), 1, b,-1)));
  add(new_ptr((Tree2toNDiagram(2), g, g, 1, state_->CC(), 1, b->CC(), 1, c,-1)));
}

Selector<MEBase::DiagramIndex>
GGtoBCQQbarBase::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) 
    sel.insert(1.0, i);
  return sel;
}

Selector<const ColourLines *>
GGtoBCQQbarBase::colourGeometries(tcDiagPtr) const {
  // {'c2': 'T(aa,ab,if,ie)', 'c3': 'T(ab,aa,if,ie)'}
  static const ColourLines c1[2] = {ColourLines("1 5, -1 2, -2 -4"),
				    ColourLines("1 -2, -1 -4, 2 5")};

  Selector<const ColourLines *> sel;
  for(unsigned int ix=0;ix<2;++ix)
    sel.insert(meInfo()[ix], &c1[ix]);
  return sel;
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
