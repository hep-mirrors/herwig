#include "DiagramContainer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

DiagramContainer::DiagramContainer(const MEPtr me):theME(me){}

DiagramContainer::~DiagramContainer() {}


IBPtr DiagramContainer::clone() const {
  return new_ptr(*this);
}

IBPtr DiagramContainer::fullclone() const {
  return new_ptr(*this);
}

void DiagramContainer::persistentOutput(PersistentOStream &) const {
}

void DiagramContainer::persistentInput(PersistentIStream &, int) {
}

void DiagramContainer::Init() {

  static ClassDocumentation<DiagramContainer> documentation
    ("DiagramContainer");
}
