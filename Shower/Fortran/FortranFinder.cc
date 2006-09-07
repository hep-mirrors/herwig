// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FortranFinder class.
//

#include "FortranFinder.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "Herwig++/Shower/Base/ShowerParticle.h"

using namespace Herwig;

NoPIOClassDescription<FortranFinder> FortranFinder::initFortranFinder;
// Definition of the static class description member.

void FortranFinder::Init() {

  static ClassDocumentation<FortranFinder> documentation
    ("This class sets the initial evolution scales for the FORTRAN HERWIG shower");

}

pair<Energy,Energy> FortranFinder::
calculateInitialFinalScales(const ShowerPPair &ppair, const bool isDecayCase) {
  Energy2 ertxi=ppair.first->momentum()*ppair.second->momentum();
  if(ertxi<0.) ertxi=0.;
  if(ppair.first==ppair.second) ertxi=0.;
  if(isDecayCase) {
    Energy scale=ertxi/ppair.first->mass();
    return make_pair(0.*GeV,scale);
  }
  else {
    Energy scale=sqrt(ertxi);
    if(ertxi==0.) return make_pair(0.*GeV,0.*GeV);
    else          return make_pair(scale,scale);
  }
}

pair<Energy,Energy> FortranFinder::
calculateInitialInitialScales(const ShowerPPair &ppair) {
  Energy2 ertxi=ppair.first->momentum()*ppair.second->momentum();
  if(ertxi<0.) ertxi=0.;
  if(ppair.first==ppair.second) ertxi=0.;
  Energy scale=sqrt(ertxi);
  if(ertxi==0.) return make_pair(0.*GeV,0.*GeV);
  else          return make_pair(scale,scale);
}

pair<Energy,Energy> FortranFinder::
calculateFinalFinalScales(const ShowerPPair & ppair) {
  Energy2 ertxi=ppair.first->momentum()*ppair.second->momentum();
  if(ertxi<0.) ertxi=0.;
  if(ppair.first==ppair.second) ertxi=0.;
  Energy scale=sqrt(ertxi);
  if(ertxi==0.) return make_pair(0.*GeV,0.*GeV);
  else          return make_pair(scale,scale);
}
