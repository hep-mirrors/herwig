  // -*- C++ -*-
  //
  // MergingReweight.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
  // Copyright (C) 2002-2019 The Herwig Collaboration
  //
  // Herwig is licenced under version 3 of the GPL, see COPYING for details.
  // Please respect the MCnet academic guidelines, see GUIDELINES for details.
  //
#include "MergingReweight.h"

using namespace Herwig;


IBPtr MergingReweight::clone() const {
  return new_ptr(*this);
}

IBPtr MergingReweight::fullclone() const {
  return new_ptr(*this);
}

double MergingReweight::weight() const {
  Energy maxpt = ZERO;
  Energy ht = ZERO;
  Energy maxmjj = ZERO;
 
  const auto sub=subProcess()->head()?subProcess()->head():subProcess();

  for (auto const & out : sub->outgoing())
    if ( !onlyColoured || out->coloured() ){
      for (auto const & out2 : sub->outgoing())
        if (!onlyColoured || out2->coloured() )
          maxmjj=max(maxmjj,(out->momentum()+out2->momentum()).m());
      maxpt = max(maxpt, out->momentum().perp());
      ht+=out->momentum().perp();
    }
  if (maxpt==ZERO||ht==ZERO) {
    return 1.;
  }
  return pow(maxpt/scale, MaxPTPower)*pow(ht/scale, HTPower)*pow(maxmjj/scale,MaxMjjPower);
}

#include "ThePEG/Persistency/PersistentOStream.h"
void MergingReweight::persistentOutput(PersistentOStream & os) const {
  os << HTPower << MaxPTPower<<MaxMjjPower << ounit(scale,GeV) << onlyColoured;
}

#include "ThePEG/Persistency/PersistentIStream.h"
void MergingReweight::persistentInput(PersistentIStream & is, int) {
  is >> HTPower >> MaxPTPower >>MaxMjjPower>> iunit(scale,GeV) >> onlyColoured;
}

// Definition of the static class description member.

  // *** Attention *** The following static variable is needed for the type
  // description system in ThePEG. Please check that the template arguments
  // are correct (the class and its base class), and that the constructor
  // arguments are correct (the class name and the name of the dynamically
  // loadable library where the class implementation can be found).

#include "ThePEG/Utilities/DescribeClass.h"
DescribeClass<MergingReweight, ThePEG::ReweightBase>
describeHerwigMergingReweight("Herwig::MergingReweight", "HwDipoleShower.so");




#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
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

  
  static Parameter<MergingReweight,double> interfaceMaxMjjPower
  ("MaxMjjPower",
   "Mjj power",
   &MergingReweight::MaxMjjPower, 0.0, 0.0, 10.0, false, false, true);
  
  
  static Parameter<MergingReweight,Energy> interfaceScale
    ("Scale",
     "The reference scale",
     &MergingReweight::scale, GeV, 50.0*GeV, ZERO, ZERO,
     false, false, Interface::lowerlim);


  static Switch<MergingReweight,bool> interfaceOnlyColoured
    ("OnlyColoured",
     "Only consider coloured particles in the SubProcess when finding the minimum transverse momentum.",
     &MergingReweight::onlyColoured, false, true, false);
  static SwitchOption interfaceOnlyColouredYes
    (interfaceOnlyColoured,
     "Yes",
     "Use only coloured particles.",
     true);
  static SwitchOption interfaceOnlyColouredNo
    (interfaceOnlyColoured,
     "No",
     "Use all particles.",
     false);



}

