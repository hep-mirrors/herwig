// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2HiggsVBF class.
//

#include "MEPP2HiggsVBF.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"

using namespace Herwig;

MEPP2HiggsVBF::MEPP2HiggsVBF() : _maxflavour(5) {}

void MEPP2HiggsVBF::getDiagrams() const {
  // get the quark particle data objects as we'll bew using them
  tcPDPtr q[6],qbar[6];
  for(unsigned int ix=0;ix<5;++ix) {
    q   [ix] = getParticleData(ix+1);
    qbar[ix] = q[ix]->CC();
  }
  // and the higgs
  tcPDPtr higgs(getParticleData(ParticleID::h0));
  // WW processes
  if(process()==0||process()==1) {
    std::vector<pair<tcPDPtr,tcPDPtr> > parentpair;
    parentpair.reserve(6);
    // don't even think of putting 'break' in here!
    switch(_maxflavour) {
    case 5:
      parentpair.push_back(make_pair(getParticleData(ParticleID::b),
				     getParticleData(ParticleID::c)));
      parentpair.push_back(make_pair(getParticleData(ParticleID::b),
				     getParticleData(ParticleID::u)));
    case 4:
      parentpair.push_back(make_pair(getParticleData(ParticleID::s),
				     getParticleData(ParticleID::c)));
      parentpair.push_back(make_pair(getParticleData(ParticleID::d),
				     getParticleData(ParticleID::c)));
    case 3:
      parentpair.push_back(make_pair(getParticleData(ParticleID::s),
				     getParticleData(ParticleID::u)));
    case 2:
      parentpair.push_back(make_pair(getParticleData(ParticleID::d),
				     getParticleData(ParticleID::u)));
    default:
      ;
    }
    for(unsigned int ix=0;ix<parentpair.size();++ix) {
      for(unsigned int iy=0;iy<parentpair.size();++iy) {
 	// q1 q2 -> q1' q2' h
	add(new_ptr((Tree2toNDiagram(4), parentpair[ix].first, WMinus(), WPlus(), 
		     parentpair[iy].second, 1,
		     parentpair[ix].second, 4, parentpair[iy].first, 2, higgs,-1)));
	// q1 qbar2 -> q1' qbar2' h
	add(new_ptr((Tree2toNDiagram(4), parentpair[ix].first, WMinus(), WPlus(), 
		     parentpair[iy].first->CC(), 1,
		     parentpair[ix].second, 4, parentpair[iy].second->CC(),
		     2, higgs,-1)));
	add(new_ptr((Tree2toNDiagram(4), parentpair[ix].second->CC(), WMinus(), WPlus(), 
		     parentpair[iy].second, 1,
		     parentpair[ix].first->CC(), 4, parentpair[iy].first,
		     2, higgs,-1)));
 	// qbar1 qbar2 -> qbar1' qbar2' h
	add(new_ptr((Tree2toNDiagram(4), parentpair[ix].first->CC(), WPlus(), WMinus(), 
		     parentpair[iy].second->CC(), 1,
		     parentpair[ix].second->CC(), 4, parentpair[iy].first->CC(),
		     2, higgs,-1))); 
      }
    }
  }
  // ZZ processes
  if(process()==0||process()==2) {
    for(unsigned int ix=0;ix<_maxflavour;++ix) {
      for(unsigned int iy=ix;iy<_maxflavour;++iy) {
	// q    q    -> q    q    H
	add(new_ptr((Tree2toNDiagram(4), q[ix], Z0(), Z0(), q[iy], 
		     1, q[ix], 4, q[iy], 2, higgs,-2))); 
	// qbar qbar -> qbar qbar H
	add(new_ptr((Tree2toNDiagram(4), qbar[ix], Z0(), Z0(), qbar[iy], 
		     1, qbar[ix], 4, qbar[iy], 2, higgs,-2)));
      }
      // q    qbar -> q    qbar H
      for(unsigned int iy=0;iy<_maxflavour;++iy) {
	add(new_ptr((Tree2toNDiagram(4), q[ix], Z0(), Z0(), qbar[iy], 
		     1, q[ix], 4, qbar[iy], 2, higgs,-2))); 
      }
    }
  }
}

void MEPP2HiggsVBF::persistentOutput(PersistentOStream & os) const {
  os << _maxflavour;
}

void MEPP2HiggsVBF::persistentInput(PersistentIStream & is, int) {
  is >> _maxflavour;
}

ClassDescription<MEPP2HiggsVBF> MEPP2HiggsVBF::initMEPP2HiggsVBF;
// Definition of the static class description member.

void MEPP2HiggsVBF::Init() {

  static ClassDocumentation<MEPP2HiggsVBF> documentation
    ("The MEPP2HiggsVBF class implements Higgs production via vector-boson fusion");

  static Parameter<MEPP2HiggsVBF,unsigned int> interfaceMaxFlavour
    ( "MaxFlavour",
      "The heaviest incoming quark flavour this matrix element is allowed to handle "
      "(if applicable).",
      &MEPP2HiggsVBF::_maxflavour, 5, 0, 5, false, false, true);

}

