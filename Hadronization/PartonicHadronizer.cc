// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PartonicHadronizer class.
//

#include "PartonicHadronizer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include <ThePEG/Interface/Parameter.h> 
#include <ThePEG/Repository/EventGenerator.h>
#include <ThePEG/Interface/Reference.h>
#include <ThePEG/PDT/DecayMode.h>
#include <ThePEG/Interface/Switch.h>
#include "CluHadConfig.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void PartonicHadronizer::persistentOutput(PersistentOStream & os) const {
  os << _partonSplitter << _clusterFinder << _colourReconnector
     << _clusterFissioner << _lightClusterDecayer << _clusterDecayer << _exclusive
     << _partontries << _inter;
}

void PartonicHadronizer::persistentInput(PersistentIStream & is, int) {
  is >> _partonSplitter >> _clusterFinder >> _colourReconnector
     >> _clusterFissioner >> _lightClusterDecayer >> _clusterDecayer >> _exclusive
     >> _partontries >> _inter;
}

ClassDescription<PartonicHadronizer> PartonicHadronizer::initPartonicHadronizer;
// Definition of the static class description member.

void PartonicHadronizer::Init() {

  static ClassDocumentation<PartonicHadronizer> documentation
    ("The PartonicHadronizer is designed to be used by the DecayHandler to"
     " hadronize partonic decays of hadrons.");
  static Reference<PartonicHadronizer,PartonSplitter> 
    interfacePartonSplitter("PartonSplitter", 
		      "A reference to the PartonSplitter object", 
		      &Herwig::PartonicHadronizer::_partonSplitter,
		      false, false, true, false);
  static Reference<PartonicHadronizer,ClusterFinder> 
    interfaceClusterFinder("ClusterFinder", 
		      "A reference to the ClusterFinder object", 
		      &Herwig::PartonicHadronizer::_clusterFinder,
		      false, false, true, false);
  static Reference<PartonicHadronizer,ColourReconnector> 
    interfaceColourReconnector("ColourReconnector", 
		      "A reference to the ColourReconnector object", 
		      &Herwig::PartonicHadronizer::_colourReconnector,
		      false, false, true, false);
  static Reference<PartonicHadronizer,ClusterFissioner> 
    interfaceClusterFissioner("ClusterFissioner", 
		      "A reference to the ClusterFissioner object", 
		      &Herwig::PartonicHadronizer::_clusterFissioner,
		      false, false, true, false);
  static Reference<PartonicHadronizer,LightClusterDecayer> 
    interfaceLightClusterDecayer("LightClusterDecayer", 
		    "A reference to the LightClusterDecayer object", 
		    &Herwig::PartonicHadronizer::_lightClusterDecayer,
		    false, false, true, false);
  static Reference<PartonicHadronizer,ClusterDecayer> 
    interfaceClusterDecayer("ClusterDecayer", 
		       "A reference to the ClusterDecayer object", 
		       &Herwig::PartonicHadronizer::_clusterDecayer,
		       false, false, true, false);

  static Switch<PartonicHadronizer,bool> interface_exclusive
    ("Exclusive",
     "Ensure that the hadrons produced in the partonic decays of bottom"
     " and charm baryons do not duplicated the inclusive modes.",
     &PartonicHadronizer::_exclusive, true, false, false);
  static SwitchOption interface_exclusiveNoDuplication
    (interface_exclusive,
     "NoDuplication",
     "Forbid duplication",
     true);
  static SwitchOption interface_exclusiveDuplication
    (interface_exclusive,
     "Duplication",
     "Duplication allowed",
     false);
  
  static Switch<PartonicHadronizer,bool> interfaceIntermediates
    ("Intermediates",
     "Whether or not to include the intermediate particles produced by the"
     " cluster alogorithm in the event record.",
     &PartonicHadronizer::_inter, false, false, false);
  static SwitchOption interfaceIntermediatesIntermediates
    (interfaceIntermediates,
     "Intermediates",
     "Include the intermediates",
     true);
  static SwitchOption interfaceIntermediatesNoIntermediates
    (interfaceIntermediates,
     "NoIntermediates",
     "Don't include the intermediates.",
     false);

  static Parameter<PartonicHadronizer,unsigned int> interfacePartonic_Tries
    ("Partonic_Tries",
     "Number of attempts to generator the hadronisation if we are vetoing the"
     " reproduction of the inclusive modes in partonic hadron decays.",
     &PartonicHadronizer::_partontries, 10, 1, 1000,
     false, false, Interface::limited);
}

bool PartonicHadronizer::hadronize(tPPtr parent, 
				   const PVector & decaychildren,
				   EventHandler & ch,
				   vector<tPPtr> & outhad) {

  assert(parent->children().size() == decaychildren.size());

  unsigned int ptry(0);
  bool partonicveto(false);

  PVector currentlist = decaychildren;

  do {
    // split the gluons
    _partonSplitter->split(currentlist);

    // form the clusters
    ClusterVector clusters = _clusterFinder->formClusters(currentlist);

    _clusterFinder->reduceToTwoComponents(clusters);

    // perform colour reconnection if needed and then
    // decay the clusters into one hadron
    bool lightOK = false;
    short tried = 0;
    const ClusterVector savedclusters = clusters;
    tPVector finalHadrons;
    while (!lightOK && tried++ < 10) {
      _colourReconnector->rearrange(ch,clusters);

      finalHadrons = _clusterFissioner->fission(clusters,false);

      lightOK = _lightClusterDecayer->decay(clusters,finalHadrons);

      if (!lightOK) {
	clusters = savedclusters;
	for_each(clusters.begin(), 
		 clusters.end(), 
		 mem_fun(&Particle::undecay));
      }
    }

    if (!lightOK) {
      partonicveto=true;
      for_each(currentlist.begin(), 
	       currentlist.end(), 
	       mem_fun(&Particle::undecay));
    }
    else {
      // decay the remaining clusters
      _clusterDecayer->decay(clusters,finalHadrons);

      partonicveto = duplicateMode(parent,finalHadrons);

      if (!partonicveto)
	outhad = finalHadrons;
    }
  }
  while(++ptry < _partontries && partonicveto);
  if (!partonicveto) {
 
    // parent loses quark children
    for (unsigned int ix=0; ix < currentlist.size(); ++ix) {
      if (currentlist[ix]->coloured())
	parent->abandonChild(currentlist[ix]);
    }

    // parent gains hadron children
    // hadrons lose cluster parent
    for (unsigned ix = 0; ix < outhad.size(); ++ix) {
      parent->addChild(outhad[ix]);
      (*outhad[ix]->parents().begin())->abandonChild(outhad[ix]);
    }
  }

  return !partonicveto;
}

bool PartonicHadronizer::duplicateMode(tPPtr parent, const vector<tPPtr> & hadrons) {
  // now find the hadrons and check them
  cParticleMSet hadronsb;
  bool found(false);
  for (unsigned ix = 0; ix < hadrons.size(); ++ix)
    hadronsb.insert(hadrons[ix]->dataPtr());

  // now check particle's decay modes 
  tcPDPtr pdata(parent->dataPtr());
  Selector<tDMPtr> modes = pdata->decaySelector();
  Selector<tDMPtr>::const_iterator modeptr = modes.begin();
  Selector<tDMPtr>::const_iterator end = modes.end();


  // check not a duplicate of a known mode
  for(;modeptr!=end;++modeptr) {
    tcDMPtr mode=(*modeptr).second;
    // check same number of products
    if (mode->products().size() == hadronsb.size()) {
      ParticleMSet::const_iterator dit;
      cParticleMSet::const_iterator pit;
      
      for(dit=mode->products().begin(), pit=hadronsb.begin();
	  dit!=mode->products().end(); ++dit,++pit) {
	if((*dit)!=(*pit)) break;
      }
      if(dit == mode->products().end()) {
	found = true;
	break;
      }
    }
  }
  return found;
}
