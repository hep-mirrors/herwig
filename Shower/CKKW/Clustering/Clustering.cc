// -*- C++ -*-
//
// Clustering.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Clustering class.
//

#include "Clustering.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef HERWIG_DEBUG_CKKW
#include "ThePEG/Repository/EventGenerator.h"
#endif

#include "Clusterer.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "Clustering.tcc"
#endif


using namespace Herwig;

Clustering::~Clustering() {}

NoPIOClassDescription<Clustering> Clustering::initClustering;
// Definition of the static class description member.

void Clustering::Init() {

  static ClassDocumentation<Clustering> documentation
    ("Clustering is the base class for storing information on"
     " a clustering during reconstruction of a parton shower history.");

}

void Clustering::perform () {

#ifdef HERWIG_DEBUG_CKKW_EXTREME
  generator()->log() << "=== Clustering::perform ; this = " << this << std::endl;
#endif

  _clusterer->doKinematics(this);
  for (vector<tClusteringParticlePtr>::iterator c = _children.begin();
       c != _children.end(); ++c)
    (**c).wasClustered(this);
  unsigned int index = 0;
  for (vector<ClusteringParticlePtr>::iterator p = _parents.begin();
       p != _parents.end(); ++p) {
    (**p).emergedFromClustering(index,this);
    index += 1;
  }

}

vector<ClusteringParticlePtr> Clustering::undo () {

#ifdef HERWIG_DEBUG_CKKW_EXTREME
  generator()->log() << "=== Clustering::undo ; this = " << this << std::endl;
#endif

  for (vector<tClusteringParticlePtr>::iterator c = _children.begin();
       c != _children.end(); ++c)
    (**c).wasUnclustered();  
  return _parents;
}


#ifdef HERWIG_DEBUG_CKKW

void Clustering::debugDump (ostream& os) {

  os << "-- Clustering ------------------------------------------------------------------\n";

  os << this << endl;

  os << "clusterer " << _clusterer << " postclustering " << _postClustering
     << " clustering config " << _clusteringConfiguration << "\n";

  os << "clustering scale / GeV2 = " << _clusteringScale/GeV2
     << " alpha scale / GeV2 = " << _alphaScale/GeV2
     << " z = " << _z << " weight = " << _weight << "\n";

  if (!_parents.empty()) {
    os << "parents : \n";
    for (vector<ClusteringParticlePtr>::iterator p = _parents.begin();
	 p != _parents.end(); ++p)
      os << *p << " ";
  }

  os << "\n";

  if (!_children.empty()) {
    os << "children : \n";
    for (vector<tClusteringParticlePtr>::iterator p = _children.begin();
	 p != _children.end(); ++p)
      os << *p << " ";
  }

  os << endl;

}

#endif
