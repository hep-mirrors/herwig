// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GeneralfftoVH class.
//

#include "GeneralfftoVH.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include"ThePEG/Utilities/EnumIO.h"
#include "ThePEG/PDT/DecayMode.h"

using namespace Herwig;

GeneralfftoVH::GeneralfftoVH() {}

void GeneralfftoVH::getDiagrams() const {
  if(process_==Lepton) {
    for(long ix=11;ix<=13;ix+=2) {
      tcPDPtr eplus  = getParticleData( ix);
      tcPDPtr eminus = getParticleData(-ix);
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
  }
  else if(process_==HadronZ) {
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
    for ( int ix=1; ix<=int(maxFlavour()); ++ix ) {
      tcPDPtr q    = getParticleData(ix);
      tcPDPtr qbar = q->CC();
      for(unsigned int iz=0;iz<Zdecays.size();++iz) {
	add(new_ptr((Tree2toNDiagram(2), q, qbar,  
		     1, Z0(), 3, higgs(), 3, Z0(), 
		     5, Zdecays[iz].first,5, Zdecays[iz].second,-1)));
      }
    }
  }
  else if(process_==HadronWplus||process_==HadronWminus) {
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
    case 4:
      parentpair.push_back(make_pair(getParticleData(ParticleID::s), 
				     getParticleData(ParticleID::cbar)));
      parentpair.push_back(make_pair(getParticleData(ParticleID::d), 
				     getParticleData(ParticleID::cbar)));
    case 3:
      parentpair.push_back(make_pair(getParticleData(ParticleID::s), 
				     getParticleData(ParticleID::ubar)));
    case 2:
      parentpair.push_back(make_pair(getParticleData(ParticleID::d), 
				     getParticleData(ParticleID::ubar)));
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
      if(process_==HadronWminus) {
	for(unsigned int ix=0;ix<wmdecays.size();++ix) {
	  add(new_ptr((Tree2toNDiagram(2), parent->first, parent->second, 
		       1, WMinus(), 3, higgs(), 3, WMinus(), 
		       5, wmdecays[ix].first,5, wmdecays[ix].second,-1)));
	}
      }
      else {
	// W+ modes
	for(unsigned int ix=0;ix<wpdecays.size();++ix) {
	  add(new_ptr((Tree2toNDiagram(2), parent->second->CC(), parent->first->CC(), 
		       1, WPlus(), 3, higgs(), 3, WPlus(), 5, wpdecays[ix].first,
		       5, wpdecays[ix].second, -1)));
	}
      }
    }
  }
  else assert(false);
}

IBPtr GeneralfftoVH::clone() const {
  return new_ptr(*this);
}

IBPtr GeneralfftoVH::fullclone() const {
  return new_ptr(*this);
}

void GeneralfftoVH::persistentOutput(PersistentOStream & os) const {
  os << oenum(process_);
}

void GeneralfftoVH::persistentInput(PersistentIStream & is, int) {
  is >> ienum(process_);
}

ClassDescription<GeneralfftoVH> GeneralfftoVH::initGeneralfftoVH;
// Definition of the static class description member.

void GeneralfftoVH::Init() {

  static ClassDocumentation<GeneralfftoVH> documentation
    ("The GeneralfftoVH class provides");

}

void GeneralfftoVH::setProcessInfo(Process proc, PDPtr hin,
				   AbstractVVSVertexPtr vertex,
				   unsigned int shapeOpt) {
  higgs(hin);
  process_ = proc;
  setWWHVertex(vertex);
  lineShape(shapeOpt);
}
