// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SplitFun class.
//

#include "SplitFun.h"
#include "Pythia7/Interface/ClassDocumentation.h"
#include "Pythia7/Persistency/PersistentOStream.h"
#include "Pythia7/Persistency/PersistentIStream.h"

using namespace Herwig;


SplitFun::~SplitFun() {}


void SplitFun::persistentOutput(PersistentOStream & os) const {
  os << _interaction << _numProducts << _idEmitter << _mEmitter; 
}


void SplitFun::persistentInput(PersistentIStream & is, int) {
  int interactionInt;
  is >> interactionInt >> _numProducts >> _idEmitter >> _mEmitter; 
  _interaction = ShowerIndex::int2Interaction( interactionInt );
}


AbstractClassDescription<SplitFun> SplitFun::initSplitFun;
// Definition of the static class description member.


void SplitFun::Init() {

  static ClassDocumentation<SplitFun> documentation
    ("There is the abstract class from which all other splitting function class inherit");

}

