// -*- C++ -*-
//
// ClusteringConfiguration.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ClusteringConfiguration class.
//

#include "ClusteringConfiguration.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ClusteringConfiguration.tcc"
#endif


using namespace Herwig;

ClusteringConfiguration::~ClusteringConfiguration() {}

NoPIOClassDescription<ClusteringConfiguration> ClusteringConfiguration::initClusteringConfiguration;
// Definition of the static class description member.

void ClusteringConfiguration::Init() {

  static ClassDocumentation<ClusteringConfiguration> documentation
    ("ClusterConfiguration is used for internal communication"
     " inside a parton shower history reconstruction.");

}

#ifdef HERWIG_DEBUG_CKKW
#include "Clusterer.h"

void ClusteringConfiguration::debugDump (ostream& os) {
  os << "-- ClusteringConfiguration -----------------------------------------------------\n";
  
  os << this << "\n";

  os << "clustering indices : ";
  for(vector<unsigned int>::iterator i = _clusteringIndices.begin(); i != _clusteringIndices.end(); ++i)
    os << *i << " ";
  os << "\n";

  os << "particles to be clustered : \n";
  for (vector<ClusteringParticleData>::iterator p = _toBeClustered.begin(); p != _toBeClustered.end(); ++p)
     os << " PDG " << (*p).partonId.PDGId << " state " << (*p).partonId.state
	<< " c " << (*p).colour << " cbar " << (*p).antiColour << "\n";

  os << "particles emerging from clustering : \n";
  for (vector<ClusteringParticleData>::iterator p = _emergingFromClustering.begin(); p != _emergingFromClustering.end(); ++p)
     os << " PDG " << (*p).partonId.PDGId << " state " << (*p).partonId.state
	<< " c " << (*p).colour << " cbar " << (*p).antiColour << "\n";

  os << "clusterer : " << _clusterer << " " << _clusterer->name() << " interaction "
     << _interaction << "\n" << endl;

}

#endif
