
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TTbAModel class.
//

#include "TTbAModel.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;
using namespace ThePEG::Helicity;

void TTbAModel::doinit()  {
  addVertex(_theWPTDVertex);

  StandardModel::doinit();
}

TTbAModel::TTbAModel(): _gWPTD_L(1.0), _gWPTD_R(1.0) {}

IBPtr TTbAModel::clone() const {
  return new_ptr(*this);
}
IBPtr TTbAModel::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void TTbAModel::persistentOutput(PersistentOStream & os) const {
  os << _theWPTDVertex
     << _gWPTD_L
     << _gWPTD_R;
}

void TTbAModel::persistentInput(PersistentIStream & is, int) {
  is >> _theWPTDVertex
     >> _gWPTD_L
     >> _gWPTD_R;
    
  
}

ClassDescription<TTbAModel> TTbAModel::initTTbAModel;
// Definition of the static class description member.

void TTbAModel::Init() {
  
  static Reference<TTbAModel,ThePEG::Helicity::AbstractFFVVertex> interfaceVertexWPTD
  ("Vertex/WPTD",
   "Reference to the W prime Top Down vertex",
   &TTbAModel::_theWPTDVertex, false, false, true, false, false);

  static Parameter<TTbAModel, double> interfaceWPTDLCoupling
    ("WPTDLCoupling",
     "The left-handed W prime coupling to top down",
     &TTbAModel::_gWPTD_L, 1.0, 0., 10.0,
     false, false, Interface::limited);

  static Parameter<TTbAModel, double> interfaceWPTDRCoupling
    ("WPTDRCoupling",
     "The right-handed W prime coupling to top down",
     &TTbAModel::_gWPTD_R, 1.0, 0., 10.0,
     false, false, Interface::limited);

 
  static ClassDocumentation<TTbAModel> documentation
    ("There is no documentation for the TTbAModel class");

}

