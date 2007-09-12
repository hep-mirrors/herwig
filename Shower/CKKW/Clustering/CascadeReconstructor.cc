// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the CascadeReconstructor class.
//

#include "CascadeReconstructor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Utilities/Throw.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "CascadeReconstructor.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

CascadeReconstructor::~CascadeReconstructor() {}

void CascadeReconstructor::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
  os << _clusterers << _hardProcesses << _selector << _resolution << _forceIncreasing << _mayUseUnordered;
}

void CascadeReconstructor::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
  is >> _clusterers >> _hardProcesses >> _selector >> _resolution >> _forceIncreasing >> _mayUseUnordered;
}

ClassDescription<CascadeReconstructor> CascadeReconstructor::initCascadeReconstructor;
// Definition of the static class description member.

void CascadeReconstructor::Init() {

  static ClassDocumentation<CascadeReconstructor> documentation
    ("CascadeReconstructor is the class which performs reconstruction "
     "of parton shower histories for ME/PS merging approaches.");

  static RefVector<CascadeReconstructor,Clusterer> interfaceClusterers
    ("Clusterers",
     "Set the clusterers to be used",
     &CascadeReconstructor::_clusterers, -1, false, false, true, false, false);

  static RefVector<CascadeReconstructor,CKKWHardProcess> interfaceHardProcesses
    ("HardProcesses",
     "Set the hard processes to look at",
     &CascadeReconstructor::_hardProcesses, -1, false, false, true, false, false);


  static Reference<CascadeReconstructor,ClusteringSelector> interfaceSelector
    ("Selector",
     "Set the ClusteringSelector to be used",
     &CascadeReconstructor::_selector, false, false, true, false, false);


  static Switch<CascadeReconstructor,bool> interfaceForceIncreasing
    ("ForceIncreasing",
     "Switch on/off forcing increasing clustering scales.",
     &CascadeReconstructor::_forceIncreasing, false, false, false);
  static SwitchOption interfaceForceIncreasingOn
    (interfaceForceIncreasing,
     "On",
     "Switch on forcing increasing clustering scales.",
     true);
  static SwitchOption interfaceForceIncreasingOff
    (interfaceForceIncreasing,
     "Off",
     "Switch off forcing increasing clustering scales",
     false);


  static Switch<CascadeReconstructor,bool> interfaceMayUseUnordered
    ("MayUseUnordered",
     "Switch on/off use of unordered histories, if ordered one cannot be obtained.",
     &CascadeReconstructor::_mayUseUnordered, false, false, false);
  static SwitchOption interfaceMayUseUnorderedOn
    (interfaceMayUseUnordered,
     "On",
     "Switch on possible use of unordered histories.",
     true);
  static SwitchOption interfaceMayUseUnorderedOff
    (interfaceMayUseUnordered,
     "Off",
     "Switch off possible use of unordered histories.",
     false);

}

void CascadeReconstructor::setup () {

  // check for clusterer, etc. if incomplete, throw init exception

  if (_clusterers.empty())
    Throw<InitException>() << "CKKW : CascadeReconstructor::setup() : At least one clusterer is needed.";

  if (_hardProcesses.empty())
    Throw<InitException>() << "CKKW : CascadeReconstructor::setup() : At least one hard process needs to be specified.";

  if (!_selector)
    Throw<InitException>() << "CKKW : CascadeReconstructor::setup() : No selector given.";

}

bool CascadeReconstructor::reconstruct (PPair in, pair<double,double> x, ParticleVector out) {

#ifdef HERWIG_DEBUG_CKKW
  generator()->log() << "== CascadeReconstructor::reconstruct (ThePEG particles)" << endl;
#endif

  vector<ClusteringParticlePtr> particles = convert (in,x,out);

  for (vector<ClustererPtr>::iterator c = _clusterers.begin(); c != _clusterers.end(); ++c) {
    (**c).reset();
  }

  string excptId = "CKKW : CascadeReconstructor::reconstruct(vector<ClusteringParticlePtr>) : ";

  // clear dynamic stuff
  _clusterings.clear();
  _currentParticles.clear();
  _allClusterings.clear();

  // set _particles
  _particles = particles;

  // look, if we already got a guide for such a configuration,
  // otherwise generate one

#ifdef HERWIG_DEBUG_CKKW_EXTREME
  generator()->log() << "looking up guide ... ";
#endif

  vector<ClusteringParticleData> pdata;
  for (vector<ClusteringParticlePtr>::iterator p = _particles.begin();
       p != _particles.end(); ++p)
    pdata.push_back((**p).pData());
  map<vector<ClusteringParticleData>, ClusteringGuidePtr>::iterator theGuide
    = _guides.find(pdata);

  if (theGuide == _guides.end()) {
#ifdef HERWIG_DEBUG_CKKW_EXTREME
    generator()->log() << "none yet, generating." << endl;
#endif
    vector<tClustererPtr> clusterers; 
    for (vector<ClustererPtr>::iterator it = _clusterers.begin(); it != _clusterers.end(); ++it)
      clusterers.push_back(*it);
    vector<tCKKWHardProcessPtr> hardProcesses; 
    for (vector<CKKWHardProcessPtr>::iterator it = _hardProcesses.begin(); it != _hardProcesses.end(); ++it)
      hardProcesses.push_back(*it);
    ClusteringGuidePtr guide = new_ptr(ClusteringGuide(clusterers,
						       hardProcesses,
						       &_partitioner,
						       pdata,
						       tClusteringGuidePtr()));
    guide->eventGenerator(generator());
    bool ok = guide->generateGuide();

#ifdef HERWIG_DEBUG_CKKW_EXTREME
    generator()->log() << "generated guide ..." << endl;
    guide->debugDump(generator()->log());
#endif

    if (!ok)
      throw Exception() << excptId 
			<< "Unable to generate ClusteringGuide for given configuration."
			<< Exception::eventerror;

    theGuide = _guides.insert(make_pair(pdata,guide)).first;
  }

#ifdef HERWIG_DEBUG_CKKW_EXTREME
  generator()->log() << "found " << theGuide->second << endl;
#endif

  _guide = theGuide->second;

  appendvector(_currentParticles,_particles);

  // if ordered scales and unordered history may
  // be accepted, then first try to get an ordered
  // one. If failed, try to get an unordered one.

  if (forceIncreasing() && mayUseUnordered()) {
#ifdef HERWIG_DEBUG_CKKW_EXTREME
    generator()->log() << "we may use unordered histories. trying first to get an ordered one ..." << endl;
#endif
    bool ordered = reconstruct();
    if (ordered) return true;
    else {
#ifdef HERWIG_DEBUG_CKKW_EXTREME
    generator()->log() << "no ordered one possible, trying without ordering ..." << endl;
#endif
      // now try unordered one
      _forceIncreasing = false;
      // reset to start
      _clusterings.clear();
      _currentParticles.clear();
      _allClusterings.clear();
      _particles = particles;
      _guide = theGuide->second;
      // try again
      bool unordered = reconstruct();
      if (unordered) {
	_forceIncreasing = true;
	return true;
      } else return false;
    }
  }

  return reconstruct();
  
}

CascadeHistory CascadeReconstructor::history () const {
  CascadeHistory result = { _clusterings, _particles, _currentParticles, _hard  };
  return result;
}

bool CascadeReconstructor::reconstruct () {

#ifdef HERWIG_DEBUG_CKKW_EXTREME
  generator()->log() << "== CascadeReconstructor::reconstruct" << endl;

  generator()->log() << "the guide being used:" << endl;
  _guide->debugDumpNode (generator()->log());

  generator()->log() << "the current particles :\n";

  for(vector<tClusteringParticlePtr>::iterator c = _currentParticles.begin();
      c != _currentParticles.end(); ++c)
    (**c).debugDump(generator()->log());

#endif

  // if hard process, we're done (if at least
  // one hard process accepts the configuration)
  if (_guide->reachedHard() && !_guide->hard()->veto(_currentParticles)) {
#ifdef HERWIG_DEBUG_CKKW_EXTREME
    generator()->log() << "got a hard process" << endl;
#endif
    _resolution->hardScales(_currentParticles);
    _hard = _guide->hard();
    return true;
  }

  /**
   * - Create possible clusterings from the current guide,
   *   remove all vetoed ones. If list is empty, return false.
   *
   * - Sort them using _selector and push_front
   *   list to _allClusterings.
   *
   * - while (!_allClusterings.front().empty)
   *   perform ()
   *   while vetoed { perform(), if !vetoed break, undo }
   *   if !reconstruct(), undo
   *   else return true
   *
   */

#ifdef HERWIG_DEBUG_CKKW_EXTREME
    generator()->log() << "starting reconstruction" << endl;
#endif

  /** Create clusterings from current guide */
  map<ClusteringConfigurationPtr,ClusteringGuidePtr> configs =
    _guide->clusterings();

  list<pair<ClusteringPtr, tClusteringGuidePtr> > theClusterings;

  for (map<ClusteringConfigurationPtr,ClusteringGuidePtr>::iterator c
	 = configs.begin(); c != configs.end(); ++c) {
#ifdef HERWIG_DEBUG_CKKW_EXTREME
      generator()->log() << "considering\n";
      c->first->debugDump(generator()->log());
#endif
    // get the particles to be clustered and the emerging ones
    vector<tClusteringParticlePtr> children =
      _partitioner.project(_currentParticles,c->first->clusteringIndices());
    vector<ClusteringParticlePtr> parents;
    vector<ClusteringParticleData> emData = c->first->emergingFromClustering();
    for (vector<ClusteringParticleData>::iterator e = emData.begin();
	 e != emData.end(); ++e) {
      parents.push_back(new_ptr(ClusteringParticle(*e)));
    }
    ClusteringPtr theClustering = c->first->clusterer()->doScale(children,
								 parents,
								 c->first);
    if (theClustering) {

#ifdef HERWIG_DEBUG_CKKW_EXTREME
    generator()->log() << "got a clustering : " << endl;
    theClustering->debugDump(generator()->log());
#endif

      // if increasing scales are forced, veto clusterings
      // being out of order
      if (forceIncreasing()) {
	// only applies, if at least one clustering
	// has been performed
	if (!_clusterings.empty()) {
#ifdef HERWIG_DEBUG_CKKW_EXTREME
	  generator()->log() << "we should use ordered scales ... " << endl;
#endif
	  if (theClustering->scale() < _clusterings.back()->scale()) {
#ifdef HERWIG_DEBUG_CKKW_EXTREME
	    generator()->log() << "odered scales not fullfilled, vetoing the clustering" << endl;
#endif
	    theClustering->veto();
	  }
	}
      }
      if (!theClustering->wasVetoed()) {
#ifdef HERWIG_DEBUG_CKKW_EXTREME
	generator()->log() << "clustering has not been vetoed" << endl;
#endif
	theClusterings.push_back(make_pair(theClustering,c->second));
      } else {
#ifdef HERWIG_DEBUG_CKKW_EXTREME
	generator()->log() << "clustering has been vetoed" << endl;
#endif
      }
    }
  }
  if (theClusterings.empty()) {
#ifdef HERWIG_DEBUG_CKKW_EXTREME
    generator()->log() << "no clusterings" << endl;
#endif
    return false;
  }

  list<pair<ClusteringPtr, tClusteringGuidePtr> > sorted = _selector->sort(theClusterings);

#ifdef HERWIG_DEBUG_CKKW_EXTREME
  generator()->log() << "Clusterings after sorting : clustering -> following guide" << endl;
  for (list<pair<ClusteringPtr, tClusteringGuidePtr> >::iterator s = sorted.begin(); s != sorted.end(); ++s)
    generator()->log() << s->first << " " << s->second << endl;
#endif

  _allClusterings.push_front(sorted);

#ifdef HERWIG_DEBUG_CKKW_EXTREME
  generator()->log() << "starting to continue with next clusterings ..." << endl;
#endif

  while (true) {

    // no more clusterings left -> reconstruction failed
    if (_allClusterings.front().empty()) {
      _allClusterings.pop_front();
#ifdef HERWIG_DEBUG_CKKW_EXTREME
      generator()->log() << "no clusterings left to try" << endl;
#endif
      return false;
    }

    // try the next one
    perform ();

    // if vetoed, undo and try next one
    if (_clusterings.front()->wasVetoed()) {
#ifdef HERWIG_DEBUG_CKKW_EXTREME
      generator()->log() << "clustering has been vetoed" << endl;
#endif
      undo ();
      continue;
    }

    // else, start from this. If reconstruction failed,
    // undo and try next one. If succeeded, return true
    if (!reconstruct()) {
#ifdef HERWIG_DEBUG_CKKW_EXTREME
      generator()->log() << "subsequent reconstruction failed" << endl;
#endif
      undo ();
      continue;
    } else return true;

  }

}

void CascadeReconstructor::perform () {  

#ifdef HERWIG_DEBUG_CKKW_EXTREME
  generator()->log() << "== CascadeReconstructor::perform" << endl;
#endif

  tClusteringPtr clustering = _allClusterings.front().front().first;

#ifdef HERWIG_DEBUG_CKKW_EXTREME
  generator()->log() << "performing " << clustering << endl;
#endif

  _guide = _allClusterings.front().front().second;
  _allClusterings.front().pop_front();
  // push the clustering
  _clusterings.push_back(clustering);
  // perform the clustering
  clustering->perform();
  // append the products to all particles
  appendvector(_particles,clustering->parents());
  // get the ones unaffected by the clustering
  vector<tClusteringParticlePtr> otherOnes =
    _partitioner.projectComplement(_currentParticles,
				   clustering->clusteringConfiguration()->clusteringIndices());
  // update _currentParticles
  _currentParticles.clear();
  // append the emerging particles in order
  appendvector(_currentParticles,clustering->parents());
  // append the other aprticles keeping track of
  // indices in the clustering step to come
  unsigned int otherIndex = clustering->parents().size();
  for (vector<tClusteringParticlePtr>::iterator o = otherOnes.begin();
       o != otherOnes.end(); ++o) {
    (**o).performedClustering(otherIndex);
    _currentParticles.push_back(*o);
    otherIndex += 1;
  }
  // if the clustering has set a postClustering, perform it now
  if (clustering->postClustering()) {
    clustering->postClustering()->initialize(_currentParticles);
    for (vector<tClusteringParticlePtr>::iterator p = _currentParticles.begin();
	 p != _currentParticles.end(); ++p)
      clustering->postClustering()->doTransform(*p);
  }
}

void CascadeReconstructor::undo () {

#ifdef HERWIG_DEBUG_CKKW_EXTREME
  generator()->log() << "== CascadeReconstructor::undo" << endl;
#endif

  // remove the clustering from list
  ClusteringPtr clustering = _clusterings.back();
#ifdef HERWIG_DEBUG_CKKW_EXTREME
  generator()->log() << "undoing " << clustering << endl;
#endif
  _clusterings.pop_back();
  // get the particles emerging from clustering
  vector<ClusteringParticlePtr> parents = clustering->undo();
  // remove parents from all particles
  for (vector<ClusteringParticlePtr>::iterator p = parents.begin();
       p != parents.end(); ++p)
    _particles.erase(find(_particles.begin(),_particles.end(),*p));
  // get the particles which have been clustered
  vector<tClusteringParticlePtr> clustered = clustering->children();
  // get the particles not clustered and call undoneClustering
  vector<tClusteringParticlePtr> notClustered;
  for (vector<tClusteringParticlePtr>::iterator n = _currentParticles.begin();
       n != _currentParticles.end(); ++n)
    if (find(clustered.begin(),clustered.end(),*n) == clustered.end())
      {
	(**n).undoneClustering();
	notClustered.push_back(*n);
      }
  // sort old particles according to indices
  _currentParticles.clear();
  _currentParticles.resize(clustered.size()+notClustered.size());
  for (vector<tClusteringParticlePtr>::iterator c = clustered.begin();
       c != clustered.end(); ++c)
    _currentParticles[(**c).index()]=*c;
  for (vector<tClusteringParticlePtr>::iterator n = notClustered.begin();
       n != notClustered.end(); ++n)
    _currentParticles[(**n).index()]=*n;
  // reset the guide to the previous one
  _guide = _guide->parent();
}


vector<ClusteringParticlePtr> CascadeReconstructor::convert (PPair in, pair<double,double> x, ParticleVector out) {

#ifdef HERWIG_DEBUG_CKKW
  generator()->log() << "== CascadeReconstructor::convert" << endl;

  generator()->log() << "starting with\n";
  generator()->log() << "incoming\n";

  dumpThePEGParticle(generator()->log(),in.first);
  dumpThePEGParticle(generator()->log(),in.second);

  generator()->log() << "incoming x's = " << x.first << " " << x.second << "\n";

  generator()->log() << "outgoing\n";

  for(ParticleVector::iterator p = out.begin(); p != out.end(); ++p)
    dumpThePEGParticle(generator()->log(),*p);

#endif

  _conversionMap.clear();

  vector<ClusteringParticlePtr> temp;

  unsigned int currentColourIndex = 499;

  map<tColinePtr,unsigned int> colourMap;

  map<tColinePtr,unsigned int>::iterator cit;

  ClusteringParticleData data;
  ClusteringParticlePtr cp;

  data.partonId.state = ClusteringParticleState::initial;

  data.partonId.PDGId = in.first->id ();
  data.colour = 0;
  data.antiColour = 0;

  if (in.first->coloured()) {

    if (in.first->hasColour()) {
      cit = colourMap.find(in.first->colourLine());
      if (cit != colourMap.end()) {
	data.colour = cit->second;
      } else {
	currentColourIndex +=1;
	colourMap.insert(make_pair(in.first->colourLine(),currentColourIndex));
	data.colour = currentColourIndex;
      }
    }

    if (in.first->hasAntiColour()) {
      cit = colourMap.find(in.first->antiColourLine());
      if (cit != colourMap.end()) {
	data.antiColour = cit->second;
      } else {
	currentColourIndex +=1;
	colourMap.insert(make_pair(in.first->antiColourLine(),currentColourIndex));
	data.antiColour = currentColourIndex;
      }
    }

  } else {
    data.colour = 0;
    data.antiColour = 0;
  }

  cp = new_ptr(ClusteringParticle(data,in.first->momentum(),x.first));
  cp->splittingScale(_resolution->minResolvableScale(data.partonId.PDGId,true));
  temp.push_back(cp);
  _conversionMap.insert(make_pair(in.first,cp));

  data.partonId.PDGId = in.second->id ();
  data.colour = 0;
  data.antiColour = 0;

  if (in.second->coloured()) {

    if (in.second->hasColour()) {
      cit = colourMap.find(in.second->colourLine());
      if (cit != colourMap.end()) {
	data.colour = cit->second;
      } else {
	currentColourIndex +=1;
	colourMap.insert(make_pair(in.second->colourLine(),currentColourIndex));
	data.colour = currentColourIndex;
      }
    }

    if (in.second->hasAntiColour()) {
      cit = colourMap.find(in.second->antiColourLine());
      if (cit != colourMap.end()) {
	data.antiColour = cit->second;
      } else {
	currentColourIndex +=1;
	colourMap.insert(make_pair(in.second->antiColourLine(),currentColourIndex));
	data.antiColour = currentColourIndex;
      }
    }

  } else {
    data.colour = 0;
    data.antiColour = 0;
  }

  cp = new_ptr(ClusteringParticle(data,in.second->momentum(),x.second));
  cp->splittingScale(_resolution->minResolvableScale(data.partonId.PDGId,true));
  temp.push_back(cp);
  _conversionMap.insert(make_pair(in.second,cp));

  // final state particles

  data.partonId.state = ClusteringParticleState::final;

  for(ParticleVector::iterator p = out.begin(); p != out.end(); ++p) {

    data.partonId.PDGId = (**p).id ();
    data.colour = 0;
    data.antiColour = 0;

    if ((**p).coloured()) {

      if ((**p).hasColour()) {
	cit = colourMap.find((**p).colourLine());
	if (cit != colourMap.end()) {
	  data.colour = cit->second;
	} else {
	  currentColourIndex +=1;
	  colourMap.insert(make_pair((**p).colourLine(),currentColourIndex));
	  data.colour = currentColourIndex;
	}
      }
      
      if ((**p).hasAntiColour()) {
	cit = colourMap.find((**p).antiColourLine());
	if (cit != colourMap.end()) {
	  data.antiColour = cit->second;
	} else {
	  currentColourIndex +=1;
	  colourMap.insert(make_pair((**p).antiColourLine(),currentColourIndex));
	  data.antiColour = currentColourIndex;
	}
      }
      
    } else {
      data.colour = 0;
      data.antiColour = 0;
    }
    
    cp = new_ptr(ClusteringParticle(data,(**p).momentum()));
    cp->splittingScale(_resolution->minResolvableScale(data.partonId.PDGId,true));
    temp.push_back(cp);
    _conversionMap.insert(make_pair(*p,cp));

  }

#ifdef HERWIG_DEBUG_CKKW

  generator()->log() << "-- conversion map --------------------------------------------------------------\n"
		     << "ThePEG::Particle -> ClusteringParticle\n";

  for(map<tPPtr,tClusteringParticlePtr>::iterator con = _conversionMap.begin();
      con != _conversionMap.end(); ++con)
    generator()->log() << con->first << " -> " << con->second << "\n";

  /* // already in debug output of reconstruct
  generator()->log() << "-- clustering particles converted ----------------------------------------------\n";

  for(vector<ClusteringParticlePtr>::iterator c = temp.begin(); c != temp.end(); ++c)
    (**c).debugDump(generator()->log());
  */

#endif

  return temp;

}

#ifdef HERWIG_DEBUG_CKKW

void CascadeReconstructor::dumpThePEGParticle (ostream& os, PPtr particle) {

  os << "-- ThePEG particle -------------------------------------------------------------\n";
  os << "Particle Ptr = " << particle << "\n";
  os << "PDG Id = " << particle->id() << "\n";
  os << "5-Momentum / GeV = ( "
     << particle->momentum().t()/GeV << " , "
     << particle->momentum().x()/GeV << " , "
     << particle->momentum().y()/GeV << " , "
     << particle->momentum().z()/GeV << " ; "
     << particle->momentum().mass2()/GeV2 << " )\n";
  os << "Colour:" << "\n";
  os << "hasColour = " << particle->hasColour() << " hasAntiColour = "
      << particle->hasAntiColour() << "\n";
  os << "Colour line = " << particle->colourLine() << " AntiColourLine = " << particle->antiColourLine() << "\n";
  os << endl;

}

#endif
