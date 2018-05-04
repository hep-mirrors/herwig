// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2WH class.
//

#include "MEPP2WH.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "Herwig/Models/StandardModel/StandardModel.h"

using namespace Herwig;

MEPP2WH::MEPP2WH() : _plusminus(0)
{}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEPP2WH,MEfftoVH>
describeHerwigMEPP2WH("Herwig::MEPP2WH", "HwMEHadron.so");

void MEPP2WH::persistentOutput(PersistentOStream & os) const {
  os << _plusminus;
}

void MEPP2WH::persistentInput(PersistentIStream & is, int) {
  is >> _plusminus;
}

void MEPP2WH::Init() {

  static ClassDocumentation<MEPP2WH> documentation
    ("The MEPP2WH class implements the matrix element for the  Bjorken"
     " process q qbar -> WH");

  static Switch<MEPP2WH,unsigned int> interfacePlusMinus
    ("Wcharge",
     "Which intermediate W bosons to include",
     &MEPP2WH::_plusminus, 0, false, false);
  static SwitchOption interfacePlusMinusAll
    (interfacePlusMinus,
     "Both",
     "Include W+ and W-",
     0);
  static SwitchOption interfacePlusMinusPlus
    (interfacePlusMinus,
     "Plus",
     "Only include W+",
     1);
  static SwitchOption interfacePlusMinusMinus
    (interfacePlusMinus,
     "Minus",
     "Only include W-",
     2);

}

void MEPP2WH::getDiagrams() const {
  // which intgermediates to include
  bool wplus =_plusminus==0||_plusminus==1;
  bool wminus=_plusminus==0||_plusminus==2;
  tPDPtr higgs = getParticleData(ParticleID::h0);
  // possible parents
  vector<PDPair> parentpair;
  parentpair.reserve(6);
  // don't even think of putting 'break' in here!
  switch(maxFlavour()) {
  case 5:
    parentpair.push_back(make_pair(getParticleData(ParticleID::b), 
				   getParticleData(ParticleID::cbar)));
    parentpair.push_back(make_pair(getParticleData(ParticleID::b), 
				   getParticleData(ParticleID::ubar)));
    [[fallthrough]];
  case 4:
    parentpair.push_back(make_pair(getParticleData(ParticleID::s), 
				   getParticleData(ParticleID::cbar)));
    parentpair.push_back(make_pair(getParticleData(ParticleID::d), 
				   getParticleData(ParticleID::cbar)));
    [[fallthrough]];
  case 3:
    parentpair.push_back(make_pair(getParticleData(ParticleID::s), 
				   getParticleData(ParticleID::ubar)));
    [[fallthrough]];
  case 2:
    parentpair.push_back(make_pair(getParticleData(ParticleID::d), 
				   getParticleData(ParticleID::ubar)));
    [[fallthrough]];
  default:
    ;
  }
  // possible children
  typedef Selector<tDMPtr> DecaySelector;
  // for W+
  DecaySelector wpdec = getParticleData(ParticleID::Wplus)->decaySelector();
  vector<PDPair> wpdecays;
  for(DecaySelector::const_iterator cit=wpdec.begin();cit!=wpdec.end();++cit) {
    if(cit->second->orderedProducts().size()!=2) continue;
    if(cit->second->orderedProducts()[0]->id()>0)
      wpdecays.push_back(make_pair(cit->second->orderedProducts()[0],
				   cit->second->orderedProducts()[1]));
    else
      wpdecays.push_back(make_pair(cit->second->orderedProducts()[1],
				   cit->second->orderedProducts()[0]));
  }
  // for W-
  DecaySelector wmdec = getParticleData(ParticleID::Wminus)->decaySelector();
  vector<PDPair> wmdecays;
  for(DecaySelector::const_iterator cit=wmdec.begin();cit!=wmdec.end();++cit) {
    if(cit->second->orderedProducts().size()!=2) continue;
    if(cit->second->orderedProducts()[0]->id()>0)
      wmdecays.push_back(make_pair(cit->second->orderedProducts()[0],
				   cit->second->orderedProducts()[1]));
    else
      wmdecays.push_back(make_pair(cit->second->orderedProducts()[1],
				   cit->second->orderedProducts()[0]));
  }
  vector<PDPair>::const_iterator parent = parentpair.begin();
  for (; parent != parentpair.end(); ++parent) {
    // W- modes
    if(wminus) {
      for(unsigned int ix=0;ix<wmdecays.size();++ix) {
	add(new_ptr((Tree2toNDiagram(2), parent->first, parent->second, 
		     1, WMinus(), 3, higgs, 3, WMinus(), 
		     5, wmdecays[ix].first,5, wmdecays[ix].second,-1)));
      }
    }
    // W+ modes
    if(wplus) {
      for(unsigned int ix=0;ix<wpdecays.size();++ix) {
	add(new_ptr((Tree2toNDiagram(2), parent->second->CC(), parent->first->CC(), 
		     1, WPlus(), 3, higgs, 3, WPlus(), 5, wpdecays[ix].first,
		     5, wpdecays[ix].second, -1)));
      }
    }  
  }
}

void MEPP2WH::doinit() {
  // get the vedrtex pointers from the SM object
  tcHwSMPtr hwsm= dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  if(!hwsm)
    throw InitException() << "Wrong type of StandardModel object in "
			  << "MEPP2WH::doinit() the Herwig"
			  << " version must be used" 
			  << Exception::runerror;
  // set the vertex
  setWWHVertex(hwsm->vertexWWH());
  higgs(getParticleData(ParticleID::h0));
  MEfftoVH::doinit();
}
