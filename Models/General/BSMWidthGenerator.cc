// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BSMWidthGenerator class.
//

#include "BSMWidthGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Decay/General/GeneralTwoBodyDecayer.h"

using namespace Herwig;

BSMWidthGenerator::~BSMWidthGenerator() {}

void BSMWidthGenerator::persistentOutput(PersistentOStream & os) const {
  os << theModes;
}

void BSMWidthGenerator::persistentInput(PersistentIStream & is, int) {
  is >> theModes;
}

ClassDescription<BSMWidthGenerator> BSMWidthGenerator::initBSMWidthGenerator;
// Definition of the static class description member.

void BSMWidthGenerator::Init() {

  static ClassDocumentation<BSMWidthGenerator> documentation
    ("A width generator for BSM particles.");

}

void BSMWidthGenerator::setupMode(tcDMPtr mode, tDecayIntegratorPtr decayer, 
				  unsigned int) {
  tcGeneralTwoBodyDecayerPtr dec = 
    dynamic_ptr_cast<tcGeneralTwoBodyDecayerPtr>(decayer);
  theModes.push_back(make_pair(mode, dec));
}

Energy BSMWidthGenerator::partial2BodyWidth(int iloc, Energy m0,
					    Energy m1, Energy m2) const {
  if( m0 < (m1 + m2) ) return Energy();
  //need pointers to particles involved
  tcDMPtr dm = theModes[iloc].first;
  assert(dm->products().size() == 2);
  GeneralTwoBodyDecayer::PMPair p0 = make_pair(dm->parent(), m0);
  GeneralTwoBodyDecayer::PMPair p1 = make_pair(dm->orderedProducts()[0], m1);
  GeneralTwoBodyDecayer::PMPair p2 = make_pair(dm->orderedProducts()[1], m2);
  return theModes[iloc].second->partialWidth(p0, p1, p2);
}
