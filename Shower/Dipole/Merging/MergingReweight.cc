  // -*- C++ -*-
  //
  // MergingReweight.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
  // Copyright (C) 2002-2007 The Herwig Collaboration
  //
  // Herwig is licenced under version 2 of the GPL, see COPYING for details.
  // Please respect the MCnet academic guidelines, see GUIDELINES for details.
  //
#include "MergingReweight.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;


IBPtr MergingReweight::clone() const {
  return new_ptr(*this);
}

IBPtr MergingReweight::fullclone() const {
  return new_ptr(*this);
}

double MergingReweight::weight() const {
  assert(false);
  Energy maxpt = 0.*GeV;
  Energy ht=0*GeV;
  
  for ( int i = 0, N = subProcess()->outgoing().size(); i < N; ++i )
    if ( !onlyColoured || subProcess()->outgoing()[i]->coloured() ){
      maxpt = max(maxpt, subProcess()->outgoing()[i]->momentum().perp());
      ht+=subProcess()->outgoing()[i]->momentum().perp();
    }
  if (maxpt==0*GeV||ht==0*GeV) {
    cout<<"\nMergingReweight: no particles contribute to reweight. Return 1.";
    return 1.;
  }
  return pow(maxpt/scale, MaxPTPower)*pow(ht/scale, HTPower);
}

void MergingReweight::persistentOutput(PersistentOStream & os) const {
  os << HTPower << MaxPTPower << ounit(scale,GeV) << onlyColoured;
}

void MergingReweight::persistentInput(PersistentIStream & is, int) {
  is >> HTPower >> MaxPTPower >> iunit(scale,GeV) >> onlyColoured;
}

// Definition of the static class description member.

  // *** Attention *** The following static variable is needed for the type
  // description system in ThePEG. Please check that the template arguments
  // are correct (the class and its base class), and that the constructor
  // arguments are correct (the class name and the name of the dynamically
  // loadable library where the class implementation can be found).
DescribeClass<MergingReweight, ThePEG::ReweightBase>
describeHerwigMergingReweight("Herwig::MergingReweight", "HwDipoleShower.so");


void MergingReweight::Init() {

  static ClassDocumentation<MergingReweight> documentation
    ("There is no documentation for the ThePEG::MergingReweight class");

  static Parameter<MergingReweight,double> interfacePower
    ("HTPower",
     "Ht power",
     &MergingReweight::HTPower, 4.0, -10.0, 10.0, false, false, true);

  
  static Parameter<MergingReweight,double> interfaceMaxPtPower
  ("MaxPTPower",
   "PT2 power",
   &MergingReweight::MaxPTPower, 4.0, -10.0, 10.0, false, false, true);

  

  
  
  static Parameter<MergingReweight,Energy> interfaceScale
    ("Scale",
     "The reference scale",
     &MergingReweight::scale, GeV, 50.0*GeV, ZERO, ZERO,
     false, false, Interface::lowerlim);


  static Switch<MergingReweight,bool> interfaceOnlyColoured
    ("OnlyColoured",
     "Only consider coloured particles in the SubProcess when finding the minimum transverse momentum.",
     &MergingReweight::onlyColoured, false, true, false);
  static SwitchOption interfaceOnlyColouredTrue
    (interfaceOnlyColoured,
     "True",
     "Use only coloured particles.",
     true);
  static SwitchOption interfaceOnlyColouredFalse
    (interfaceOnlyColoured,
     "False",
     "Use all particles.",
     false);



}

