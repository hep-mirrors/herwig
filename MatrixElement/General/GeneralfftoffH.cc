// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GeneralfftoffH class.
//

#include "GeneralfftoffH.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include"ThePEG/Utilities/EnumIO.h"

using namespace Herwig;

GeneralfftoffH::GeneralfftoffH() {}

IBPtr GeneralfftoffH::clone() const {
  return new_ptr(*this);
}

IBPtr GeneralfftoffH::fullclone() const {
  return new_ptr(*this);
}

void GeneralfftoffH::persistentOutput(PersistentOStream & os) const {
  os << oenum(_proc);
}

void GeneralfftoffH::persistentInput(PersistentIStream & is, int) {
  is >> ienum(_proc);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<GeneralfftoffH,MEfftoffH>
describeHerwigGeneralfftoffH("Herwig::GeneralfftoffH", "Herwig.so");

void GeneralfftoffH::Init() {

  static ClassDocumentation<GeneralfftoffH> documentation
    ("There is no documentation for the GeneralfftoffH class");

}

void GeneralfftoffH::getDiagrams() const {
  if(_proc==Lepton) {
    for(long ix=11;ix<=13;ix+=1) {
      tcPDPtr em(getParticleData(ix));
      tcPDPtr ep(em->CC());
      // WW processes
      if(process()==0||process()==1) {
	tcPDPtr nue(getParticleData(ix+1));
	tcPDPtr nueb(nue->CC());
	add(new_ptr((Tree2toNDiagram(4), em, WMinus(), WPlus(), ep, 
		     1, nue, 4, nueb, 2, higgs(),-1))); 
      }
      // ZZ processes
      if(process()==0||process()==2) {
	add(new_ptr((Tree2toNDiagram(4), em, Z0(), Z0(), ep, 
		     1, em, 4, ep, 2, higgs(),-2))); 
      }
    }
  }
  else {
    // WW processes
    if(process()==0||process()==1) {
      std::vector<pair<tcPDPtr,tcPDPtr> > parentpair;
      parentpair.reserve(6);
      // don't even think of putting 'break' in here!
      switch(maxFlavour()) {
      case 5:
	if (minFlavour()<=4)
	  parentpair.push_back(make_pair(getParticleData(ParticleID::b),
					 getParticleData(ParticleID::c)));
	if (minFlavour()<=2)
	  parentpair.push_back(make_pair(getParticleData(ParticleID::b),
					 getParticleData(ParticleID::u)));
	[[fallthrough]];
      case 4:
	if (minFlavour()<=3)
	  parentpair.push_back(make_pair(getParticleData(ParticleID::s),
					 getParticleData(ParticleID::c)));
	if (minFlavour()<=1)
	  parentpair.push_back(make_pair(getParticleData(ParticleID::d),
					 getParticleData(ParticleID::c)));
	[[fallthrough]];
      case 3:
	if (minFlavour()<=2)
	  parentpair.push_back(make_pair(getParticleData(ParticleID::s),
					 getParticleData(ParticleID::u)));
	[[fallthrough]];
      case 2:
	if (minFlavour()<=1)
	  parentpair.push_back(make_pair(getParticleData(ParticleID::d),
					 getParticleData(ParticleID::u)));
	[[fallthrough]];
      default:
	;
      }
      for(unsigned int ix=0;ix<parentpair.size();++ix) {
	for(unsigned int iy=0;iy<parentpair.size();++iy) {
	  // q1 q2 -> q1' q2' h
	  if(parentpair[ix].first->id()<parentpair[iy].second->id()) {
	    add(new_ptr((Tree2toNDiagram(4), parentpair[ix].first, WMinus(), WPlus(), 
			 parentpair[iy].second, 1, parentpair[ix].second, 4, 
			 parentpair[iy].first, 2, higgs(),-1)));
	  }
	  else {
	    add(new_ptr((Tree2toNDiagram(4), parentpair[iy].second, WPlus(), WMinus(), 
			 parentpair[ix].first, 1, parentpair[iy].first, 4,
			 parentpair[ix].second, 2, higgs(),-1)));
	  }
	  // q1 qbar2 -> q1' qbar2' h
	  add(new_ptr((Tree2toNDiagram(4), parentpair[ix].first, WMinus(), WPlus(), 
		       parentpair[iy].first->CC(), 1,
		       parentpair[ix].second, 4, parentpair[iy].second->CC(),
		       2, higgs(),-1)));
	  add(new_ptr((Tree2toNDiagram(4),parentpair[iy].second, WPlus(), WMinus(),
		       parentpair[ix].second->CC(), 1, parentpair[iy].first,
		       4, parentpair[ix].first->CC(), 
		       2, higgs(),-1)));
	  // qbar1 qbar2 -> qbar1' qbar2' h
	  if(parentpair[ix].first->id()<parentpair[ix].second->id()) {
	    add(new_ptr((Tree2toNDiagram(4), parentpair[ix].first->CC(),
			 WPlus(), WMinus(), 
			 parentpair[iy].second->CC(), 1,
			 parentpair[ix].second->CC(), 4, parentpair[iy].first->CC(),
			 2, higgs(),-1))); 
	  }
	  else {
	    add(new_ptr((Tree2toNDiagram(4), parentpair[iy].second->CC(),
			 WMinus(), WPlus(),
			 parentpair[ix].first->CC(), 1, 
			 parentpair[iy].first->CC(), 4, parentpair[ix].second->CC(),
			 2, higgs(),-1))); 
	  }
	}
      }
    }
    // ZZ processes
    if(process()==0||process()==2) {
      // get the quark particle data objects as we'll be using them
      tcPDPtr q[6],qbar[6];
      for ( int ix=0; ix<5; ++ix ) {
	q   [ix] = getParticleData(ix+1);
	qbar[ix] = q[ix]->CC();
      }
      for(unsigned int ix=minFlavour()-1;ix<maxFlavour();++ix) {
	for(unsigned int iy=ix;iy<maxFlavour();++iy) {
	  // q    q    -> q    q    H
	  add(new_ptr((Tree2toNDiagram(4), q[ix], Z0(), Z0(), q[iy], 
		       1, q[ix], 4, q[iy], 2, higgs(),-2))); 
	  // qbar qbar -> qbar qbar H
	  add(new_ptr((Tree2toNDiagram(4), qbar[ix], Z0(), Z0(), qbar[iy], 
		       1, qbar[ix], 4, qbar[iy], 2, higgs(),-2)));
	}
	// q    qbar -> q    qbar H
	for(unsigned int iy=minFlavour()-1;iy<maxFlavour();++iy) {
	  add(new_ptr((Tree2toNDiagram(4), q[ix], Z0(), Z0(), qbar[iy], 
		       1, q[ix], 4, qbar[iy], 2, higgs(),-2))); 
	}
      }
    }
  }
}

void GeneralfftoffH::setProcessInfo(Process proc, PDPtr hin,
				    AbstractVVSVertexPtr vertex,
				    unsigned int shapeOpt,
				    unsigned int iproc) {
  higgs(hin);
  _proc = proc;
  setWWHVertex(vertex);
  lineShape(shapeOpt);
  process(iproc);
}
