// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LeptoquarkModel class.
//

#include "LeptoquarkModel.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;
using namespace ThePEG::Helicity;

void LeptoquarkModel::doinit()  {
  addVertex(_theSLQSLQGVertex);
  addVertex(_theSLQSLQGGVertex);
  addVertex(_theSLQFFVertex);
  
  StandardModel::doinit();
}

LeptoquarkModel::LeptoquarkModel() :  _CouplFF(0.312), _leftcoup(1.0), _rightcoup(1.0), _rightcouptilde(1.0), _leftcoup1(1.0)  {}

LeptoquarkModel::~LeptoquarkModel() {}

IBPtr LeptoquarkModel::clone() const {
  return new_ptr(*this);
}
IBPtr LeptoquarkModel::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void LeptoquarkModel::persistentOutput(PersistentOStream & os) const {
  os <<  _theSLQSLQGGVertex
     << _theSLQSLQGVertex
     << _theSLQFFVertex
     << _CouplFF
     << _leftcoup
     << _rightcoup
     << _leftcoup1
     << _rightcouptilde;

    
  
}

void LeptoquarkModel::persistentInput(PersistentIStream & is, int) {
  is >> _theSLQSLQGGVertex
     >> _theSLQSLQGVertex
     >> _theSLQFFVertex
     >> _CouplFF
     >> _leftcoup
     >> _rightcoup
     >> _leftcoup1
     >> _rightcouptilde;
}

ClassDescription<LeptoquarkModel> LeptoquarkModel::initLeptoquarkModel;
// Definition of the static class description member.

void LeptoquarkModel::Init() {
  
  static Reference<LeptoquarkModel,ThePEG::Helicity::AbstractVSSVertex> interfaceVertexSLQSLQG
  ("Vertex/SLQSLQG",
   "Reference to the scalar leptoquark-scalar leptoquark-gluon vertex",
   &LeptoquarkModel::_theSLQSLQGVertex, false, false, true, false, false);

  static Reference<LeptoquarkModel,ThePEG::Helicity::AbstractVVSSVertex> interfaceVertexSLQSLQGG
  ("Vertex/SLQSLQGG",
   "Reference to the scalar leptoquark-scalar leptoquark-gluon-gluon vertex",
   &LeptoquarkModel::_theSLQSLQGGVertex, false, false, true, false, false);

  static Reference<LeptoquarkModel,ThePEG::Helicity::AbstractFFSVertex> interfaceVertexSLQFF
  ("Vertex/SLQFF",
   "Reference to the scalar leptoquark-scalar-quark-lepton",
   &LeptoquarkModel::_theSLQFFVertex, false, false, true, false, false);


  static Parameter<LeptoquarkModel, double> interfaceLQCoupling
    ("LQCoupling",
     "The overall Leptoquark Coupling",
     &LeptoquarkModel::_CouplFF, 0.312, 0., 10.0,
     false, false, Interface::limited);

  static Parameter<LeptoquarkModel, double> interfacegLQ_L
    ("g_S0_L",
     "The leptoquark S0 coupling LQ-lepton_left-quark_right",
     &LeptoquarkModel::_leftcoup, 1.0, 0., 1.0,
     false, false, Interface::limited);

  static Parameter<LeptoquarkModel, double> interfacegLQ_R
    ("g_S0_R",
     "The leptoquark S0 coupling LQ-lepton_right-quark_left",
     &LeptoquarkModel::_rightcoup, 1.0, 0., 1.0,
     false, false, Interface::limited);

  static Parameter<LeptoquarkModel, double> interfacegLQ_Rt
    ("g_S0t_R",
     "The leptoquark ~S0 coupling LQ-lepton_right-quark_left",
     &LeptoquarkModel::_rightcouptilde, 1.0, 0., 1.0,
     false, false, Interface::limited);

  static Parameter<LeptoquarkModel, double> interfacegLQ_L1
    ("g_S1_L",
     "The leptoquark S1 coupling LQ-lepton_left-quark_right",
     &LeptoquarkModel::_leftcoup1, 1.0, 0., 1.0,
     false, false, Interface::limited);
  
  static ClassDocumentation<LeptoquarkModel> documentation
    ("There is no documentation for the LeptoquarkModel class");

}

