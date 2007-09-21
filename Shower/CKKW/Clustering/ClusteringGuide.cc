// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ClusteringGuide class.
//

#include "ClusteringGuide.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef HERWIG_DEBUG_CKKW
#include "ThePEG/Repository/EventGenerator.h"
#endif

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ClusteringGuide.tcc"
#endif


using namespace Herwig;

ClusteringGuide::~ClusteringGuide() {}



bool ClusteringGuide::generateGuide () {

#ifdef HERWIG_DEBUG_CKKW_EXTREME
  generator()->log() << "=== ClusteringGuide::generateGuide" << endl;
  generator()->log() << "working on : \n" ;
  for (unsigned int p =0; p< _particles.size(); ++p) {
    generator()->log() << "#" << p
		       << " PDG " << _particles[p].partonId.PDGId << " state " << _particles[p].partonId.state
		       << " c " << _particles[p].colour << " cbar " << _particles[p].antiColour << "\n";
  }
#endif

  string excptId = "CKKW : ClusteringGuide::generateGuide() : ";

  // first do some consistency checking

  if (_clusterers.empty() || _hardProcesses.empty())
    throw Exception () 
      << excptId << "No clusterers or hard processes found."
      << Exception::runerror;

  // stop on hard process

    for (vector<tCKKWHardProcessPtr>::iterator p = _hardProcesses.begin();
	 p != _hardProcesses.end(); ++p) {
      if((**p).reachedHard(_particles)) {
#ifdef HERWIG_DEBUG_CKKW_EXTREME
	generator()->log() << "got hard process " << (**p).name() << endl;
#endif
	hard(*p);
	_hardProcess = _particles;
	return true;
      }
    }

  // get the possible multiplicities of particles to be clustered
  // also check that we have existing clusteres
  list<unsigned int> multiplicities;
  for (vector<tClustererPtr>::iterator c = _clusterers.begin();
       c != _clusterers.end(); ++c) {
    if (*c) {
      if (find(multiplicities.begin(),multiplicities.end(),(**c).toClusterMultiplicity())
	  == multiplicities.end())
	multiplicities.push_back((**c).toClusterMultiplicity());
    }
    else
      throw Exception () 
	<< excptId << "NULL pointer to Clusterer encountered."
	<< Exception::runerror;      
  }
  multiplicities.sort();

  /**
   * For each multiplicity:
   *
   * - First get the possible partitions.
   *
   *   For each partition:
   *
   *   - Check each clusterer to return the configurations.
   *   - Add each configuration together with a new node to
   *     candidates, constructing the node from _clusterers,
   *     _hardProcesses, _partitioner and configuration->emergingFromClustering()
   *     + the complement of the current particle content wrt the partition
   *     considered and setting the _clusteringIndices of the configuration.
   *
   * For each candidate:
   *
   * - If candidate.second.particles is a hard process, add to clusterings.
   * - If not hard, but candidate.second.generateGuide(), add to clusterings.
   *
   */

  map<ClusteringConfigurationPtr,ClusteringGuidePtr> candidateClusterings;
  ClusteringConfigurationPtr candidateConfig;
  ClusteringGuidePtr candidateGuide;

  for (list<unsigned int>::iterator m = multiplicities.begin();
       m != multiplicities.end(); ++m) {

    // get the possible index sets
    vector<vector<unsigned int> > indexSets;
    _partitioner->index_partitions(0,_particles.size(),*m,indexSets);

#ifdef HERWIG_DEBUG_CKKW_EXTREME
    for(vector<vector<unsigned int> >::iterator iset = indexSets.begin();
	iset != indexSets.end(); ++iset) {
      generator()->log() << "{ ";
      for (vector<unsigned int>::iterator index = iset->begin();
	   index != iset->end(); ++index)
	generator()->log() << *index << " ";
      generator()->log() << "} ";
    }
    generator()->log() << endl;
#endif

    // for each partition
    for (vector<vector<unsigned int> >::iterator is = indexSets.begin();
	 is != indexSets.end(); ++is) {

#ifdef HERWIG_DEBUG_CKKW_EXTREME
      generator()->log() << "considering ";
      for (vector<unsigned int>::iterator i = is->begin(); i != is->end(); ++i)
	generator()->log() << *i << " ";
      generator()->log() << endl;
#endif
      
      // get possible configurations
      vector<ClusteringConfigurationPtr> configurations;
      for (vector<tClustererPtr>::iterator c = _clusterers.begin();
	   c != _clusterers.end(); ++c) {
	vector<ClusteringConfigurationPtr> cconfigurations = (**c).configurations(_partitioner->project(_particles,*is));
#ifdef HERWIG_DEBUG_CKKW_EXTREME
	generator()->log() << "from clusterer " << (**c).name() << " found : \n";
	for (vector<ClusteringConfigurationPtr>::iterator dc = cconfigurations.begin();
	     dc != cconfigurations.end(); ++dc)
	  (**dc).debugDump(generator()->log());
#endif
	appendvector(configurations,cconfigurations);
      }

      // add configuration candidates
      for (vector<ClusteringConfigurationPtr>::iterator conf = configurations.begin();
	   conf != configurations.end(); ++conf) {
	// set the index set
	(**conf)._clusteringIndices = *is;
	// construct the nex node
	vector<ClusteringParticleData> newconf = (**conf).emergingFromClustering();
	appendvector(newconf,_partitioner->projectComplement(_particles,*is));
	ClusteringGuidePtr nextnode =
	  new_ptr(ClusteringGuide(_clusterers,_hardProcesses,_partitioner,newconf, this));
	nextnode->eventGenerator(generator());
	candidateClusterings.insert(make_pair(*conf,nextnode));
      }

    } // end for each partition

  } // end for each multiplicity

  // for each candidate
  for (map<ClusteringConfigurationPtr,ClusteringGuidePtr>::iterator cand =
	 candidateClusterings.begin(); cand != candidateClusterings.end(); ++cand) {
    bool useit = false;

    // have we got a hard process ?
    for (vector<tCKKWHardProcessPtr>::iterator p = _hardProcesses.begin();
	 p != _hardProcesses.end(); ++p) {
      if((**p).reachedHard(cand->second->particles())) {
#ifdef HERWIG_DEBUG_CKKW_EXTREME
	generator()->log() << "got hard process (2) " << (**p).name() << endl;
#endif
	cand->second->hard(*p);
	useit = true;
	break;
      }
    }

    // have we reached a hard process after further clustering?
    if (!useit) useit = cand->second->generateGuide();

    if (useit)
      _clusterings.insert(*cand);

  }

  // are there any clusterings to be used at all?
  if (_clusterings.empty()) return false;
  else return true;

}

NoPIOClassDescription<ClusteringGuide> ClusteringGuide::initClusteringGuide;
// Definition of the static class description member.

void ClusteringGuide::Init() {

  static ClassDocumentation<ClusteringGuide> documentation
    ("ClusteringGuide stores information on possible clustering sequences "
     "yielding a valid CKKW subprocess.");

}

#ifdef HERWIG_DEBUG_CKKW

void ClusteringGuide::debugDump (ostream& os, const string& feed) {
  string myfeed = feed;
  if (myfeed == "")
    debugDumpNode(os,myfeed,true);
  else {
    myfeed += "  ";
    debugDumpNode(os,myfeed,false);
  }

  if (!_clusterings.empty()) {
    for (map<ClusteringConfigurationPtr,ClusteringGuidePtr>::iterator c = _clusterings.begin();
	 c != _clusterings.end(); ++c) {
      c->second->debugDump(os,myfeed);
    }
  }

}

void ClusteringGuide::debugDumpNode (ostream& os, const string& feed, bool showc) {

  os << feed << "-- ClusteringGuide -------------------------------------------------------------\n";

  os << feed << this << "\n";

  if (showc) {
    os << feed << "clusterers : \n";
    for (vector<tClustererPtr>::iterator c = _clusterers.begin(); c != _clusterers.end(); ++c) {
      os << feed << *c << " " << (**c).name() << "\n";
    }
    os << feed << "hard processes : \n";
    for (vector<tCKKWHardProcessPtr>::iterator c = _hardProcesses.begin(); c != _hardProcesses.end(); ++c) {
      os << feed << *c << " " << (**c).name() << "\n";
    }
    os << feed << "partitioner : " << _partitioner << "\n";
  }

  os << feed << "Particles at this node : \n";

  for (vector<ClusteringParticleData>::iterator p = _particles.begin(); p != _particles.end(); ++p) {
     os << feed << " PDG " << (*p).partonId.PDGId << " state " << (*p).partonId.state
	<< " c " << (*p).colour << " cbar " << (*p).antiColour << "\n";
  }

  os << feed << "parent clustering guide: " << _parent << "\n";

  if (!_clusterings.empty()) {
    os << feed << "clusterings at this node : \n";
    for (map<ClusteringConfigurationPtr,ClusteringGuidePtr>::iterator c = _clusterings.begin();
	 c != _clusterings.end(); ++c) {
      c->first->debugDump(os);
      os << feed << "clustering guide continuing here: " << c->second << "\n";
    }
  }

  if (_hard)
    os << feed << "hard process : " << _hard << " " << _hard->name() << "\n";

  if(!_hardProcess.empty()) {
    os << feed << "particles at the hard process : \n";
    for (vector<ClusteringParticleData>::iterator p = _hardProcess.begin(); p != _hardProcess.end(); ++p) {
      os << feed << " PDG " << (*p).partonId.PDGId << " state " << (*p).partonId.state
	 << " c " << (*p).colour << " cbar " << (*p).antiColour << "\n";
    }
  }

}

#endif
