// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEee2HiggsVBF class.
//

#include "MEee2HiggsVBF.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "Herwig/Models/StandardModel/StandardModel.h"

using namespace Herwig;

void MEee2HiggsVBF::getDiagrams() const {
  // WW processes
  tcPDPtr em(getParticleData(ParticleID::eminus));
  tcPDPtr ep(em->CC());
  if(process()==0||process()==1) {
    tcPDPtr nue(getParticleData(ParticleID::nu_e));
    tcPDPtr nueb(nue->CC());
    add(new_ptr((Tree2toNDiagram(4), em, WMinus(), WPlus(), ep, 
		 1, nue, 3, nueb, 2, higgs(),-1))); 
  }
  // ZZ processes
  if(process()==0||process()==2) {
    add(new_ptr((Tree2toNDiagram(4), em, Z0(), Z0(), ep, 
		 1, em, 3, ep, 2, higgs(),-2))); 
  }
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeNoPIOClass<MEee2HiggsVBF,MEfftoffH>
describeHerwigMEee2HiggsVBF("Herwig::MEee2HiggsVBF", "HwMELepton.so");

void MEee2HiggsVBF::Init() {

  static ClassDocumentation<MEee2HiggsVBF> documentation
    ("The MEee2HiggsVBF class implements VBF type matrix elements for e+e- collisions");

}

void MEee2HiggsVBF::doinit() {
  // get the vedrtex pointers from the SM object
  tcHwSMPtr hwsm= dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  if(!hwsm)
    throw InitException() << "Wrong type of StandardModel object in "
			  << "MEPP2HiggsVBF::doinit() the Herwig"
			  << " version must be used" 
			  << Exception::runerror;
  // set the vertex
  setWWHVertex(hwsm->vertexWWH());
  higgs(getParticleData(ParticleID::h0));
  MEfftoffH::doinit();
}

