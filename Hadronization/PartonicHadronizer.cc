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

bool PartonicHadronizer::hadronize(tPPtr parent,StepPtr pstep,EventHandler & ch,
				   vector<tPPtr> & outhad) {
  unsigned int ptry(0);
  bool partonicveto(false);
  vector<tcPPtr> pclusters;
  tPVector tagged;
  do {
    tagged.clear();
    ClusterVector clusters;
    for(unsigned int ix=0;ix<parent->children().size();++ix) {
      if(parent->children()[ix]->coloured())
	tagged.push_back(parent->children()[ix]);
    }
    // split the gluons
    tPVector partons=_partonSplitter->split(tagged,pstep);
    // form the clusters
    _clusterFinder->formClusters(ch.currentCollision(),pstep,partons,clusters);
    _clusterFinder->reduceToTwoComponents(pstep,clusters);
    // perform colour reconnection if needed and then
    // decay the clusters into one hadron
    bool lightOK = false;
    short tried = 0;
    while (!lightOK && tried++ < 10) {
      _colourReconnector->rearrange(ch,pstep,clusters);
      _clusterFissioner->fission(pstep,false);
      lightOK = _lightClusterDecayer->decay(pstep);
      if (!lightOK) {
	for(unsigned int ix=0;ix<clusters.size();++ix) {
	  for(unsigned int iy=0;iy<clusters[ix]->children().size();++iy) {
	    tPPtr ptemp=clusters[ix]->children()[iy];
	    // remove parents which are free quarks
	    for(unsigned int iz=0;iz<ptemp->parents().size();++iz) {
	      if(ptemp->parents()[iz]->coloured()&&
		 ptemp->parents()[iz]->parents().empty()) {
		pstep->removeParticle(ptemp->parents()[iz]);
	      }
	    }
	    // remove children which aren't clusters
	    if(ptemp->id()!=81) pstep->removeParticle(ptemp);	    
	  }
	}
      }
    }
    if (!lightOK) {
      partonicveto=true;
      for(unsigned int ix=0;ix<clusters.size();++ix) {
	for(unsigned int iy=0;iy<clusters[ix]->children().size();++iy) {
	  tPPtr ptemp=clusters[ix]->children()[iy];
	  // remove parents which are free quarks
	  for(unsigned int iz=0;iz<ptemp->parents().size();++iz) {
	    if(ptemp->parents()[iz]->coloured()&&
	       ptemp->parents()[iz]->parents().empty()) {
	      pstep->removeParticle(ptemp->parents()[iz]);
	    }
	  }
	  if(ptemp->id()!=81) pstep->removeParticle(ptemp);
	}
      }
      // remove the children of the tagged particles
      for(unsigned int ix=0;ix<tagged.size();++ix) {
	for(unsigned int iy=0;iy<tagged[ix]->children().size();++iy)
	  pstep->removeParticle(tagged[ix]->children()[iy]);
      }
    }
    else {
      // decay the remaining clusters
      _clusterDecayer->decay(pstep);
      // find the clusters
      findPartonicClusters(*pstep,parent,pclusters);
      // check if duplicate modes in partonic decays
      partonicveto=duplicateMode(parent,pclusters,outhad);
      // check if vetoing the hadronization
      if(partonicveto) {
	// find the motherless quarks
	ParticleVector quarks;
	for(unsigned int ix=0;ix<outhad.size();++ix)
	  removeQuarks(outhad[ix],quarks);
	// remove the children of the tagged particles
	for(unsigned int ix=0;ix<tagged.size();++ix) {
	  for(int iy=tagged[ix]->children().size()-1;iy>=0;--iy)
	    pstep->removeParticle(tagged[ix]->children()[iy]);
	}
	// remove motherless quarks
	for(unsigned int ix=0;ix<quarks.size();++ix)
	  pstep->removeParticle(quarks[ix]);
      }
    }
    // increment counter
    ++ptry;
  }
  while(ptry<_partontries&&partonicveto);
  // remove the intermediate particles added by the cluster model if needed
  if(!_inter&&!partonicveto) {
      // find the motherless quarks
    ParticleVector quarks;
    for(unsigned int ix=0;ix<outhad.size();++ix)
      removeQuarks(outhad[ix],quarks);
    // remove the tagged particles
    for(unsigned int ix=0;ix<tagged.size();++ix)
      pstep->removeParticle(tagged[ix]);
    // remove motherless quarks
    for(unsigned int ix=0;ix<quarks.size();++ix)
      pstep->removeParticle(quarks[ix]);
    // add the outgoing hadrons as the children of the decaying particle
    for(unsigned int ix=0;ix<outhad.size();++ix)
      pstep->addDecayProduct(parent,outhad[ix]);
  }
  return !partonicveto;
}


void PartonicHadronizer::findPartonicClusters(Step & pstep,
					      tPPtr parent,vector<tcPPtr> & cluout) {
  vector<tPPtr> clusters;
  ParticleSet::iterator pit=pstep.intermediates().begin();
  // first find the clusters which decay to hadrons
  for(;pit!=pstep.intermediates().end();++pit) {
    if((**pit).id()==81&&hadronicCluster(*pit)) clusters.push_back(*pit);
  }
  // now find the decaying hadrons they match with
  for(unsigned int ix=0;ix<clusters.size();++ix) {
    bool found=false;
    tcPPtr mother,part;
    part=clusters[ix];
    do {
      mother=part->parents()[0];
      if(mother) {
	// looking for the quark which produced the cluster
	if(mother->coloured()) {
	  if(part->parents().size()>1) {
	    // make sure quark not from cluster splitting
	    if(part->parents()[1]->id()==81) {
	      mother=part->parents()[1]; 
	    }
	    else if(mother->parents().size()>0) {
	      if(mother->parents()[0]->id()!=ParticleID::g) found=true;
	    }
	    else {
	      found=true;
	    }
	  }
	  else if(mother->parents().size()>0) {
	    if(mother->parents()[0]->id()!=ParticleID::g) found=true;
	  }
	  else {
	    found=true;
	  }
	}
	part=mother;
      }
      else {
	found=true;
      }
    }
    while(!found);
    // add cluster to list if quark has the right parent
    if(part->parents().size()!=0) {
      if(part->parents()[0]==parent) cluout.push_back(clusters[ix]);
    }
  }
}

bool PartonicHadronizer::hadronicCluster(tPPtr part) {
  if(part->children().size()==0&&part->id()==81) return true;
  else if(part->id()==81) {
    unsigned int nclu(0);
    for(unsigned int ix=0;ix<part->children().size();++ix) {
      if(part->children()[ix]->id()==81) ++nclu;
    }
    if(nclu==part->children().size()&&nclu!=1) return false;
    else return true;
  }
  else return false;
}

bool PartonicHadronizer::duplicateMode(tPPtr parent,vector<tcPPtr> & clusters,
				       vector<tPPtr> & hadrons) {
  // now find the hadrons and check them
  cParticleMSet hadronsb;
  bool found(false);
  // find the hadrons produced in the decay
  hadrons.clear();
  for(unsigned int ix=0;ix<clusters.size();++ix) {
    for(unsigned int iy=0;iy<clusters[ix]->children().size();++iy) {
      if(clusters[ix]->children()[iy]->id()!=81) {
	hadronsb.insert(clusters[ix]->children()[iy]->dataPtr());
	hadrons.push_back(clusters[ix]->children()[iy]);
      }
    }
  }
  // now check particle's decay modes 
  tcPDPtr pdata(parent->dataPtr());
  Selector<tDMPtr> modes=pdata->decaySelector();
  Selector<tDMPtr>::const_iterator start=modes.begin();
  Selector<tDMPtr>::const_iterator end=modes.end();
  ParticleMSet::const_iterator dit;
  cParticleMSet::const_iterator pit;
  tcDMPtr mode;
  // check not a duplicate of a known mode
  for(;start!=end;++start) {
    mode=(*start).second;
    // check same number of products
    if(mode->products().size()==hadronsb.size()) {
      for(dit=mode->products().begin(),pit=hadronsb.begin();
	  dit!=mode->products().end();++dit,++pit) {
	if((*dit)!=(*pit)) break;
      }
      if(dit==mode->products().end()) {
	found=true;break;
      }
    }
  }
  return found;
}

void PartonicHadronizer::removeQuarks(tPPtr outhad,ParticleVector & quarks) {
  tPVector parents=outhad->parents();
  while(!parents.empty()) {
    PPtr part=parents.back();
    parents.pop_back();
    if(part->coloured()&&part->parents().empty())
      quarks.push_back(part);
    else if(part->id()==81)                      
      removeQuarks(part,quarks);
  }
}
