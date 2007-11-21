// -*- C++ -*-
//
// ClusteringParticle.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ClusteringParticle class.
//

#include "ClusteringParticle.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#include "Clustering.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ClusteringParticle.tcc"
#endif


using namespace Herwig;

ClusteringParticle::~ClusteringParticle() {}

  vector<tClusteringParticlePtr> ClusteringParticle::children () const {
    if (_splitting) return _splitting->children();
    else return vector<tClusteringParticlePtr>();
  }

  vector<ClusteringParticlePtr> ClusteringParticle::parents () const {
    if (_production) return _production->parents();
    else return vector<ClusteringParticlePtr>();
  }

NoPIOClassDescription<ClusteringParticle> ClusteringParticle::initClusteringParticle;
// Definition of the static class description member.

void ClusteringParticle::Init() {

  static ClassDocumentation<ClusteringParticle> documentation
    ("ClusteringParticle stores all information relevant for particles "
     "in a parton shower history.");

}

#ifdef HERWIG_DEBUG_CKKW

void ClusteringParticle::debugDump (ostream& os) {

  os << "-- ClusteringParticle ----------------------------------------------------------\n";
  os << this
     << " PDG " << _data.partonId.PDGId << " state " << _data.partonId.state
     << " c " << _data.colour << " cbar " << _data.antiColour << "\n";

  os << "5-Momentum / GeV = ( "
     << _momentum.t()/GeV << " , "
     << _momentum.x()/GeV << " , "
     << _momentum.y()/GeV << " , "
     << _momentum.z()/GeV << " ; "
     << _momentum.mass2()/GeV2 << " ) "
     << " x = " << _x << "\n";

  os << "production scale / GeV2 = " << _productionScale/GeV2
     << " splitting scale / GeV2 = " << _splittingScale/GeV2
     << " shower scale / GeV2 = " << _showerScale/GeV2 << "\n";

  os << "production clustering " << _production
     << " splitting clustering " << _splitting << "\n";

  if (!_indexStack.empty()) {
    os << "index stack : ";
    for (list<unsigned int>::iterator i = _indexStack.begin(); 
	 i != _indexStack.end(); ++i)
      os << *i << " ";
    os << "\n";
  }

  if (!_momentumBackup.empty()) {
    os << "momentum backup :\n";
    for (list<Lorentz5Momentum>::iterator m = _momentumBackup.begin(); m != _momentumBackup.end(); ++m) {
      os << "5-Momentum / GeV = ( "
	 << (*m).t()/GeV << " , "
	 << (*m).x()/GeV << " , "
	 << (*m).y()/GeV << " , "
	 << (*m).z()/GeV << " ; "
	 << (*m).mass2()/GeV2 << " ) ";
	}
  }

  if (!_xBackup.empty()) {
    os << "x backup : ";
    for (list<double>::iterator i = _xBackup.begin(); 
	 i != _xBackup.end(); ++i)
      os << *i << " ";
    os << "\n";
  }

  os << endl;

}

#endif

