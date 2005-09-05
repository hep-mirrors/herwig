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

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "PartonicHadronizer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;

PartonicHadronizer::~PartonicHadronizer() {}

void PartonicHadronizer::persistentOutput(PersistentOStream & os) const {
  os << _globalParameters  << _partonSplitter << _clusterFinder << _colourReconnector
     << _clusterFissioner << _lightClusterDecayer << _clusterDecayer << _exclusive
     << _partontries;
}

void PartonicHadronizer::persistentInput(PersistentIStream & is, int) {
  is >> _globalParameters  >> _partonSplitter >> _clusterFinder >> _colourReconnector
     >> _clusterFissioner >> _lightClusterDecayer >> _clusterDecayer >> _exclusive
     >> _partontries;
}

ClassDescription<PartonicHadronizer> PartonicHadronizer::initPartonicHadronizer;
// Definition of the static class description member.

void PartonicHadronizer::Init() {

  static ClassDocumentation<PartonicHadronizer> documentation
    ("This is the main handler class for the Cluster Hadronization");

  static Reference<PartonicHadronizer,GlobalParameters> 
    interfaceGlobalParameters("GlobalParameters", 
		      "A reference to the GlobalParameters object", 
		      &Herwig::PartonicHadronizer::_globalParameters,
		      false, false, true, false);
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
    ("_exclusive",
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
  
  static Parameter<PartonicHadronizer,unsigned int> interfacePartonic_Tries
    ("Partonic_Tries",
     "Number of attempts to generator the hadronisation if we are vetoing the"
     " reproduction of the inclusive modes in partonic hadron decays.",
     &PartonicHadronizer::_partontries, 10, 1, 1000,
     false, false, Interface::limited);
}

bool PartonicHadronizer::hadronize(tPPtr parent,StepPtr pstep,EventHandler & ch,
				   vector<tPPtr> & outhad)
{
  unsigned int ptry(0);
  bool partonicveto(false);
  do
    {
      ClusterVector clusters;
      tPVector tagged;
      for(unsigned int ix=0;ix<parent->children().size();++ix)
	{if(parent->children()[ix]->coloured())
	    {tagged.push_back(parent->children()[ix]);}}
      // split the gluons
      tPVector partons=_partonSplitter->split(tagged,pstep);
      // form the clusters
      _clusterFinder->formClusters(ch.currentCollision(),pstep,partons,clusters);
      _clusterFinder->reduceToTwoComponents(pstep,clusters);
      // perform colour reconnection if needed and then
      // decay the clusters into one hadron
      bool lightOK = false;
      short tried = 0;
      vector<tcPPtr> pclusters;
      while (!lightOK && tried++ < 10) 
	{
	  _colourReconnector->rearrange(ch,pstep,clusters);
	  _clusterFissioner->fission(pstep);
	  lightOK = _lightClusterDecayer->decay(pstep);
	  if (!lightOK) 
	    {
	      for(unsigned int ix=0;ix<clusters.size();++ix)
		{
		  for(unsigned int iy=0;iy<clusters[ix]->children().size();++iy)
		    {if(clusters[ix]->children()[iy]->id()!=81)
			{pstep->removeParticle(clusters[ix]->children()[iy]);}}
		}
	    }
	}
      if (!lightOK)
	{
	  partonicveto=true;
	  for(unsigned int ix=0;ix<clusters.size();++ix)
	    {
	      for(unsigned int iy=0;iy<clusters[ix]->children().size();++iy)
		{if(clusters[ix]->children()[iy]->id()!=81)
		    {pstep->removeParticle(clusters[ix]->children()[iy]);}}
	    }
	  // remove the children of the tagged particles
	  for(unsigned int ix=0;ix<tagged.size();++ix)
	    {
	      for(unsigned int iy=0;iy<tagged[ix]->children().size();++iy)
		{pstep->removeParticle(tagged[ix]->children()[iy]);}
	    }
	}
      else
	{
	  // decay the remaining clusters
	  _clusterDecayer->decay(pstep);
	  // find the clusters
	  findPartonicClusters(*pstep,parent,pclusters);
	  // check if duplicate modes in partonic decays
	  partonicveto=duplicateMode(parent,pclusters,outhad);
	  // check if vetoing the hadronization
	  if(partonicveto)
	    {
	      // special for the motherless quarks
	      for(unsigned int ix=0;ix<pclusters.size();++ix)
		{
		  for(unsigned int iy=0;iy<pclusters[ix]->parents().size();++iy)
		    {
		      if(pclusters[ix]->parents()[iy]->coloured()&&
			 pclusters[ix]->parents()[iy]->parents().empty())
			{pstep->removeParticle(pclusters[ix]->parents()[iy]);}
		    }
		}
	      // remove the children of the tagged particles
	      for(unsigned int ix=0;ix<tagged.size();++ix)
		{
		  for(unsigned int iy=0;iy<tagged[ix]->children().size();++iy)
		    {pstep->removeParticle(tagged[ix]->children()[iy]);}
		}
	    }
	}
      // increment counter
      ++ptry;
    }
  while(ptry<_partontries&&partonicveto);
  return ptry!=_partontries;
}


void PartonicHadronizer::findPartonicClusters(Step & pstep,
					      tPPtr parent,vector<tcPPtr> & cluout)
{
  vector<tPPtr> clusters;
  ParticleSet::iterator pit=pstep.intermediates().begin();
  // first find the clusters which decay to hadrons
  for(;pit!=pstep.intermediates().end();++pit)
    {if((**pit).id()==81){if(hadronicCluster(*pit)){clusters.push_back(*pit);}}}
  // now find the decaying hadrons they match with
  for(unsigned int ix=0;ix<clusters.size();++ix)
    {
      bool found=false;
      tcPPtr mother,part;
      part=clusters[ix];
      do
	{
	  mother=part->parents()[0];
	  if(mother)
	    {
	      // looking for the quark which produced the cluster
	      if(mother->coloured())
		{
		  if(part->parents().size()>1)
		    {
		      // make sure quark not from cluster splitting
		      if(part->parents()[1]->id()==81)
			{mother=part->parents()[1];}
		      else if(mother->parents().size()>0)
			{if(mother->parents()[0]->id()!=ParticleID::g){found=true;}}
		      else{found=true;}
		    }
		  else if(mother->parents().size()>0)
		    {if(mother->parents()[0]->id()!=ParticleID::g){found=true;}}
		  else{found=true;}
		}
	      part=mother;
	    }
	  else{found=true;}
	}
      while(!found);
      // add cluster to list if quark has the right parent
      if(part->parents().size()!=0)
	{if(part->parents()[0]==parent){cluout.push_back(clusters[ix]);}}
    }
}

bool PartonicHadronizer::hadronicCluster(tPPtr part)
{
  if(part->children().size()==0&&part->id()==81){return true;}
  else if(part->id()==81)
    {
      for(unsigned int ix=0;ix<part->children().size();++ix)
	{
	  if(part->children()[ix]->id()==81){return false;}
	}
      return true;
    }
  else{return false;}
}

bool PartonicHadronizer::duplicateMode(tPPtr parent,vector<tcPPtr> & clusters,
				       vector<tPPtr> & hadrons)
{
  // now find the hadrons and check them
  cParticleMSet hadronsb;
  bool found(false);
  // find the hadrons produced in the decay
  hadrons.resize(0);
  for(unsigned int ix=0;ix<clusters.size();++ix)
    {
      for(unsigned int iy=0;iy<clusters[ix]->children().size();++iy)
	{
	  if(clusters[ix]->children()[iy]->id()!=81)
	    {
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
  for(;start!=end;++start)
    {
      mode=(*start).second;
      // check same number of products
      if(mode->products().size()==hadronsb.size())
	{
	  for(dit=mode->products().begin(),pit=hadronsb.begin();
	      dit!=mode->products().end();++dit,++pit)
	    {if((*dit)!=(*pit)){break;}}
	  if(dit==mode->products().end()){found=true;break;}
	}
    }
  return found;
}

}
