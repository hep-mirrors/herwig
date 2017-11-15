// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEee2ZH class.
//

#include "MEee2ZH.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig/Models/StandardModel/StandardModel.h"

using namespace Herwig;

// The following static variable is needed for the type
// description system in ThePEG.
DescribeNoPIOClass<MEee2ZH,MEfftoVH>
describeHerwigMEee2ZH("Herwig::MEee2ZH", "HwMELepton.so");

void MEee2ZH::Init() {

  static ClassDocumentation<MEee2ZH> documentation
    ("There is no documentation for the MEee2ZH class");

}

void MEee2ZH::getDiagrams() const {
  tcPDPtr eplus  = getParticleData(ParticleID::eplus  );
  tcPDPtr eminus = getParticleData(ParticleID::eminus );
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
  for(unsigned int ix=0;ix<Zdecays.size();++ix) {
    add(new_ptr((Tree2toNDiagram(2), eminus, eplus, 
		 1, Z0(), 3, higgs(), 3, Z0(), 
		 5, Zdecays[ix].first,5, Zdecays[ix].second,-1)));
  }
}

void MEee2ZH::doinit() {
  // get the vedrtex pointers from the SM object
  tcHwSMPtr hwsm= dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  if(!hwsm)
    throw InitException() << "Wrong type of StandardModel object in "
			  << "MEeeto2ZH::doinit() the Herwig"
			  << " version must be used" 
			  << Exception::runerror;
  // set the vertex
  setWWHVertex(hwsm->vertexWWH());
  higgs(getParticleData(ParticleID::h0));
  MEfftoVH::doinit();
}
