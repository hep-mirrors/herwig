// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEee2HiggsVBF class.
//

#include "MEee2HiggsVBF.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"

using namespace Herwig;

void MEee2HiggsVBF::getDiagrams() const {
  // WW processes
  tcPDPtr em(getParticleData(ParticleID::eminus));
  tcPDPtr ep(em->CC());
  tcPDPtr higgs(getParticleData(ParticleID::h0));
  if(process()==0||process()==1) {
    tcPDPtr nue(getParticleData(ParticleID::nu_e));
    tcPDPtr nueb(nue->CC());
    add(new_ptr((Tree2toNDiagram(4), em, WMinus(), WPlus(), ep, 
		 1, nue, 4, nueb, 2, higgs,-1))); 
  }
  // ZZ processes
  if(process()==0||process()==2) {
    add(new_ptr((Tree2toNDiagram(4), em, Z0(), Z0(), ep, 
		 1, em, 4, ep, 2, higgs,-2))); 
  }
}

NoPIOClassDescription<MEee2HiggsVBF> MEee2HiggsVBF::initMEee2HiggsVBF;
// Definition of the static class description member.

void MEee2HiggsVBF::Init() {

  static ClassDocumentation<MEee2HiggsVBF> documentation
    ("The MEee2HiggsVBF class implements VBF type matrix elements for e+e- collisions");

}

