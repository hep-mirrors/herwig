// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QtildaShowerKinematics1to2 class.
//

#include "QtildaShowerKinematics1to2.h"
#include "Pythia7/Interface/ClassDocumentation.h"

using namespace Herwig;


QtildaShowerKinematics1to2::~QtildaShowerKinematics1to2() {}


AbstractClassDescription<QtildaShowerKinematics1to2> 
QtildaShowerKinematics1to2::initQtildaShowerKinematics1to2;
// Definition of the static class description member.


void QtildaShowerKinematics1to2::Init() {

  static ClassDocumentation<QtildaShowerKinematics1to2> documentation
    ("This abstract class describes the common features for initial and final ",
     "state radiation kinematics for 1 -> 2 branchings and for the choice of ",
     "Qtilda as evolution variable.");

}


