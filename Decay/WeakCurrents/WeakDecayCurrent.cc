// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the WeakDecayCurrent class.
//

#include "WeakDecayCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "WeakDecayCurrent.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;

WeakDecayCurrent::~WeakDecayCurrent() {}

void WeakDecayCurrent::persistentOutput(PersistentOStream & os) const {
  os << _quark << _antiquark;
}

void WeakDecayCurrent::persistentInput(PersistentIStream & is, int) {
  is >> _quark >> _antiquark;
}

AbstractClassDescription<WeakDecayCurrent> WeakDecayCurrent::initWeakDecayCurrent;
// Definition of the static class description member.

void WeakDecayCurrent::Init() {

  static ClassDocumentation<WeakDecayCurrent> documentation
    ("The \\classname{WeakDecayCurrent} class is the basse class for the"
     " implementation of hadronic currents in weak decays.");

  static ParVector<WeakDecayCurrent,int> interfaceQuark
    ("Quark",
     "The PDG code for the quark.",
     &WeakDecayCurrent::_quark,
     0, 0, 0, 0, 20, false, false, true);

  static ParVector<WeakDecayCurrent,int> interfaceAntiQuark
    ("AntiQuark",
     "The PDG code for the antiquark.",
     &WeakDecayCurrent::_antiquark,
     0, 0, 0, -20, 0, false, false, true);
}

}
