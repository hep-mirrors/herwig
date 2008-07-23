// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2ZH class.
//

#include "MEPP2ZH.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/PDT/DecayMode.h"

using namespace Herwig;

MEPP2ZH::MEPP2ZH() : _maxflavour(5)
{}

void MEPP2ZH::getDiagrams() const {
  tPDPtr higgs = getParticleData(ParticleID::h0);
  // find possible Z decays
  typedef Selector<tDMPtr> DecaySelector;
  DecaySelector Zdec = Z0()->decaySelector();
  vector<PDPair> Zdecays;
  for(DecaySelector::const_iterator cit=Zdec.begin();cit!=Zdec.end();++cit) {
    if(cit->second->orderedProducts().size()!=2) continue;
    if(cit->second->orderedProducts()[0]->id()>0)
      Zdecays.push_back(make_pair(cit->second->orderedProducts()[0],
				  cit->second->orderedProducts()[1]));
    else
      Zdecays.push_back(make_pair(cit->second->orderedProducts()[1],
				  cit->second->orderedProducts()[0]));
  }
  // create the diagrams
  for(unsigned int ix=1;ix<=_maxflavour;++ix) {
    tcPDPtr q    = getParticleData(ix);
    tcPDPtr qbar = q->CC();
    for(unsigned int iz=0;iz<Zdecays.size();++iz) {
      add(new_ptr((Tree2toNDiagram(2), q, qbar,  
		   1, Z0(), 3, higgs, 3, Z0(), 
		   5, Zdecays[iz].first,5, Zdecays[iz].second,-1)));
    }
  }
}

void MEPP2ZH::persistentOutput(PersistentOStream & os) const {
  os << _maxflavour;
}

void MEPP2ZH::persistentInput(PersistentIStream & is, int) {
  is >> _maxflavour;
}

ClassDescription<MEPP2ZH> MEPP2ZH::initMEPP2ZH;
// Definition of the static class description member.

void MEPP2ZH::Init() {

  static ClassDocumentation<MEPP2ZH> documentation
    ("The MEPP2ZH class implements the matrix element for q qbar -> Z H");

  static Parameter<MEPP2ZH,unsigned int> interfaceMaxFlavour
    ( "MaxFlavour",
      "The heaviest incoming quark flavour this matrix element is allowed to handle "
      "(if applicable).",
      &MEPP2ZH::_maxflavour, 5, 1, 5, false, false, true);
}

