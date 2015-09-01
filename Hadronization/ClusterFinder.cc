// -*- C++ -*-
//
// ClusterFinder.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ClusterFinder class.
//

#include "ClusterFinder.h"
#include <ThePEG/Interface/ClassDocumentation.h>
#include <ThePEG/Interface/Switch.h>
#include <ThePEG/Persistency/PersistentOStream.h>
#include <ThePEG/Persistency/PersistentIStream.h>
#include <ThePEG/PDT/StandardMatchers.h>
#include <ThePEG/PDT/EnumParticles.h>
#include <ThePEG/Repository/EventGenerator.h>
#include <ThePEG/EventRecord/Collision.h>
#include "CheckId.h"
#include "Herwig/Utilities/EnumParticles.h"
#include "Herwig/Utilities/Kinematics.h"
#include "Cluster.h"
#include <ThePEG/Utilities/DescribeClass.h>

using namespace Herwig;

DescribeClass<ClusterFinder,Interfaced>
describeClusterFinder("Herwig::ClusterFinder","");

IBPtr ClusterFinder::clone() const {
  return new_ptr(*this);
}

IBPtr ClusterFinder::fullclone() const {
  return new_ptr(*this);
}
void ClusterFinder::persistentOutput(PersistentOStream & os) const {
  os << heavyDiquarks_ << diQuarkSelection_ << diQuarkOnShell_;
}

void ClusterFinder::persistentInput(PersistentIStream & is, int) {
  is >> heavyDiquarks_ >> diQuarkSelection_ >> diQuarkOnShell_;
}

void ClusterFinder::Init() {

  static ClassDocumentation<ClusterFinder> documentation
    ("This class is responsible of finding clusters.");

  static Switch<ClusterFinder,unsigned int> interfaceHeavyDiquarks
    ("HeavyDiquarks",
     "How to treat heavy quarks in baryon number violating clusters",
     &ClusterFinder::heavyDiquarks_, 2, false, false);
  static SwitchOption interfaceHeavyDiquarksDefault
    (interfaceHeavyDiquarks,
     "Allow",
     "No special treatment, allow both heavy and doubly heavy diquarks",
     0);
  static SwitchOption interfaceHeavyDiquarksNoDoublyHeavy
    (interfaceHeavyDiquarks,
     "NoDoublyHeavy",
     "Avoid having diquarks with twoo heavy quarks",
     1);
  static SwitchOption interfaceHeavyDiquarksNoHeavy
    (interfaceHeavyDiquarks,
     "NoHeavy",
     "Try and avoid diquarks contain c and b altogether",
     2);

  static Switch<ClusterFinder,unsigned int> interfaceDiQuarkSelection
    ("DiQuarkSelection",
     "Option controlling the selection of quarks to merge into a diquark in baryon-number violating clusters",
     &ClusterFinder::diQuarkSelection_, 1, false, false);
  static SwitchOption interfaceDiQuarkSelectionRandom
    (interfaceDiQuarkSelection,
     "Random",
     "Randomly pick a pair to combine",
     0);
  static SwitchOption interfaceDiQuarkSelectionLowestMass
    (interfaceDiQuarkSelection,
     "LowestMass",
     "Combine the lowest mass pair",
     1);

  static Switch<ClusterFinder,bool> interfaceDiQuarkOnShell
    ("DiQuarkOnShell",
     "Force the diquark produced in baryon-number violating clusters to be on-shell",
     &ClusterFinder::diQuarkOnShell_, false, false, false);
  static SwitchOption interfaceDiQuarkOnShellYes
    (interfaceDiQuarkOnShell,
     "Yes",
     "Force to be on-shell",
     true);
  static SwitchOption interfaceDiQuarkOnShellNo
    (interfaceDiQuarkOnShell,
     "No",
     "Leave off-shell",
     false);

}

ClusterVector ClusterFinder::formClusters(const PVector & partons) {

  set<tPPtr> examinedSet;  // colour particles already included in a cluster
  map<tColinePtr, pair<tPPtr,tPPtr> > quarkQuark; // quark quark 
  map<tColinePtr, pair<tPPtr,tPPtr> > aQuarkQuark; // anti quark anti quark
  ParticleSet inputParticles(partons.begin(),partons.end());

  ClusterVector clusters;
  
  // Loop over all current particles.
  for(PVector::const_iterator pit=partons.begin();pit!=partons.end();++pit){
    // Skip to the next particle if it is not coloured or already examined.
    assert(*pit);
    assert((*pit)->dataPtr());
    if(!(**pit).data().coloured() 
       || examinedSet.find(*pit) != examinedSet.end()) {
      continue;  
    }
    // We assume that a cluster is  made of, at most, 3 constituents although
    // in most cases the number will be 2; however, for baryon violating decays
    // (for example in Susy model without R parity conservation) we can have 3
    // constituents. In the latter case, a quark (antiquark) do not have an
    // anticolour (colour) partner as usual, but its colour line either stems 
    // from a colour source, or ends in a colour sink. In the case of double
    // baryon violating decays, but with overall baryon conservation 
    //  ( for instance:
    //       tilde_u_R -> dbar_1 + dbar_2 
    //       tilde_u_R_star -> d1 + d2 
    //    where tilde_u_R and tilde_u_R_star are colour connected )
    // a special treatment is needed, because first we have to process all
    // partons in the current step, and then for each left pair of quarks which
    // stem from a colour source we have to find the corresponding pair of 
    // anti-quarks which ends in a colour sink and is connected with the
    // above colour source. These special pairs are kept into the maps:
    //    spec/CluHadConfig.hialQuarkQuarkMap   and   specialAntiQuarkAntiQuarkMap.

    tParticleVector connected(3);
    int iElement = 0;         
    connected[iElement++] = *pit;
    bool specialCase = false;

    if((*pit)->hasColour()) {
      tPPtr partner = 
	(*pit)->colourLine()->getColouredParticle(partons.begin(),
						  partons.end(),
						  true);

      if(partner) {
	connected[iElement++]= partner;
      }        
      // colour source : baryon-violating process
      else {
	if((*pit)->colourLine()->sourceNeighbours() != tColinePair()) {
	  tColinePair sourcePair = (*pit)->colourLine()->sourceNeighbours();
	  tColinePtr intCL = tColinePtr();
	  for(int i = 0; i < 2; ++i) {
	    tColinePtr pLine = i==0 ? sourcePair.first : sourcePair.second;
            int saveNumElements = iElement;
	    for(tPVector::const_iterator cit = pLine->coloured().begin(); 
		cit != pLine->coloured().end(); ++cit ) {
	      ParticleSet::const_iterator cjt = inputParticles.find(*cit);
	      if(cjt!=inputParticles.end()) connected[iElement++]= (*cit);
	    }      
	    if(iElement == saveNumElements) intCL = pLine;
	  }
	  if(intCL && iElement == 2) {
	    specialCase = true;
	    pair<tPPtr,tPPtr> qp=pair<tPPtr,tPPtr>(connected[0],connected[1]);
	    quarkQuark.insert(pair<tColinePtr,pair<tPPtr,tPPtr> >(intCL,qp)); 
	  }
	  else if(iElement != 3) {
	    throw Exception() << "Colour connections fail in the hadronization for " 
			      << **pit << "in ClusterFinder::formClusters"
			      << " for a coloured particle."
			      << " Failed to find particles from a source"
			      << Exception::runerror;
	  }
	}
	else {
	  throw Exception() << "Colour connections fail in the hadronization for " 
			    << **pit << "in ClusterFinder::formClusters for"
			    << " a coloured particle"
			    << Exception::runerror;
	}
      }
    }

    if((*pit)->hasAntiColour()) {
      tPPtr partner = 
	(*pit)->antiColourLine()->getColouredParticle(partons.begin(),
						      partons.end(),
						      false);
      if(partner) {
	connected[iElement++]=partner; 
      }
      // colour sink : baryon-violating process
      else {
        if((*pit)->antiColourLine()->sinkNeighbours() != tColinePair()) {
	  tColinePair sinkPair = (*pit)->antiColourLine()->sinkNeighbours();
	  tColinePtr intCL = tColinePtr();
	  for(int i = 0; i < 2; ++i) {
	    tColinePtr pLine = i==0 ? sinkPair.first : sinkPair.second;
            int saveNumElements = iElement;
	    for(tPVector::const_iterator cit = pLine->antiColoured().begin(); 
		cit != pLine->antiColoured().end(); ++cit ) {
	      ParticleSet::const_iterator cjt = inputParticles.find(*cit);
	      if(cjt!=inputParticles.end()) connected[iElement++]= (*cit);
	    }
	    if(iElement == saveNumElements) intCL = pLine;
	  }
	  if(intCL && iElement == 2) {
	    specialCase = true;
	    pair<tPPtr,tPPtr> aqp=pair<tPPtr,tPPtr>(connected[0],connected[1]);
	    aQuarkQuark.insert(pair<tColinePtr,pair<tPPtr,tPPtr> >(intCL,aqp));
	  }
	  else if( iElement !=3) {
	    throw Exception() << "Colour connections fail in the hadronization for "
			      << **pit << "in ClusterFinder::formClusters for"
			      << " an anti-coloured particle."
			      << " Failed to find particles from a sink"
			      << Exception::runerror;
	  }
	}
	else {
	  throw Exception() << "Colour connections fail in the hadronization for " 
			    << **pit << "in ClusterFinder::formClusters for"
			    << " an anti-coloured particle"
			    << Exception::runerror;
	}
      }
    }

    if(!specialCase) {
      // Tag the components of the found cluster as already examined.
      for (int i=0; i<iElement; ++i) examinedSet.insert(connected[i]);
      // Create the cluster object with the colour connected particles
      ClusterPtr cluPtr = new_ptr(Cluster(connected[0],connected[1],
					  connected[2]));
      // add to the step
      connected[0]->addChild(cluPtr);
      connected[1]->addChild(cluPtr);
      if(connected[2]) connected[2]->addChild(cluPtr);
      clusters.push_back(cluPtr);
      // Check if any of the components is a beam remnant, and if this
      // is the case then inform the cluster.   
      // this will only work for baryon collisions  
      for (int i=0; i<iElement; ++i) {
	bool fromRemnant = false;
	tPPtr parent=connected[i];
	while(parent) {
	  if(parent->id()==ParticleID::Remnant) {
	    fromRemnant = true;
	    break;
	  }
	  parent = parent->parents().empty() ? tPPtr() : parent->parents()[0];
	}
	if(fromRemnant&&DiquarkMatcher::Check(connected[i]->id()))
	  cluPtr->isBeamCluster(connected[i]);
      }
    }

  }

  // Treat now the special cases, if any. The idea is to find for each pair
  // of quarks coming from a common colour source the corresponding pair of
  // antiquarks coming from a common colour sink, connected to the above
  // colour source via the same colour line. Then, randomly couple one of
  // the two quarks with one of the two antiquarks, and do the same with the
  // quark and antiquark left.
  for(map<tColinePtr, pair<tPPtr,tPPtr> >::const_iterator 
	cit = quarkQuark.begin(); cit != quarkQuark.end(); ++cit ) {
    tColinePtr coline = cit->first;
    pair<tPPtr,tPPtr> quarkPair = cit->second;
    if(aQuarkQuark.find( coline ) != aQuarkQuark.end()) {
      pair<tPPtr,tPPtr> antiQuarkPair = aQuarkQuark.find(coline)->second;
      ClusterPtr cluPtr1, cluPtr2;
      if ( UseRandom::rndbool() ) {      
	cluPtr1 = new_ptr(Cluster(quarkPair.first , antiQuarkPair.first));
	cluPtr2 = new_ptr(Cluster(quarkPair.second , antiQuarkPair.second));
	quarkPair.first->addChild(cluPtr1);
	antiQuarkPair.first->addChild(cluPtr1);
	quarkPair.second->addChild(cluPtr2);
	antiQuarkPair.second->addChild(cluPtr2);
      } else {
	cluPtr1 = new_ptr(Cluster(quarkPair.first , antiQuarkPair.second));
	cluPtr2 = new_ptr(Cluster(quarkPair.second , antiQuarkPair.first));
	quarkPair.second->addChild(cluPtr2);
	antiQuarkPair.first->addChild(cluPtr2);
	quarkPair.first->addChild(cluPtr1);
	antiQuarkPair.second->addChild(cluPtr1);
      }
      clusters.push_back(cluPtr1);
      clusters.push_back(cluPtr2);
    } 
    else {
      throw Exception() << "ClusterFinder::formClusters : " 
			<< "***Skip event: unable to match pairs in "
			<< "Baryon-violating processes***"
			<< Exception::eventerror;
    }
  }
  return clusters;
}

namespace {
  bool PartOrdering(tPPtr p1,tPPtr p2) {
    return abs(p1->id())<abs(p2->id());
  }
}

void ClusterFinder::reduceToTwoComponents(ClusterVector & clusters) {

  // In order to preserve all of the information, we do not modify the 
  // directly the 3-component clusters, but instead we define new clusters,
  // which are related to the original ones by a child-parent relationship,
  // by considering two randomly chosen components as a diquark (or anti-diquark).
  // These new clusters are first added to the vector  vecNewRedefinedCluPtr,
  // and at the end, when all input clusters have been examined, the elements of 
  // this vector will be copied in  collecCluPtr  (the reason is that it is not 
  // allowed to modify a STL container while iterating over it).
  vector<tClusterPtr> redefinedClusters; 
  for(ClusterVector::iterator cluIter = clusters.begin() ; 
      cluIter != clusters.end() ; ++cluIter) {
    tParticleVector vec;

    if ( (*cluIter)->numComponents() != 3 ||
	 ! (*cluIter)->isAvailable() ) continue;
    
    tPPtr other;
    for(int i = 0; i<(*cluIter)->numComponents(); i++) {
      tPPtr part = (*cluIter)->particle(i);
      if(!DiquarkMatcher::Check(*(part->dataPtr())))
	vec.push_back(part);
      else
	other = part;
    }

    if(vec.size()<2) {
      throw Exception() << "Could not make a diquark for a baryonic cluster decay from "
			<< (*cluIter)->particle(0)->PDGName() << " "
			<< (*cluIter)->particle(1)->PDGName() << " "
			<< (*cluIter)->particle(2)->PDGName() << " "
			<< " in ClusterFinder::reduceToTwoComponents()."
			<< Exception::eventerror;
    }

    // order the vector so heaviest at the end
    std::sort(vec.begin(),vec.end(),PartOrdering);

    // Special treatment of heavy quarks
    // avoid doubly heavy diquarks
    if(heavyDiquarks_>=1   && vec.size()>2 &&
       abs(vec[1]->id())>3 && abs(vec[0]->id())<=3) {
      if(UseRandom::rndbool()) swap(vec[1],vec[2]);
      other = vec[2];
      vec.pop_back();
    }
    // avoid singly heavy diquarks
    if(heavyDiquarks_==2   && vec.size()>2 &&
       abs(vec[2]->id())>3 && abs(vec[1]->id())<=3) {
      other = vec[2];
      vec.pop_back();
    }

    // if there's a choice pick the pair to make a diquark from
    if(vec.size()>2) {
      unsigned int ichoice(0);
      // random choice
      if(diQuarkSelection_==0) {
	ichoice = UseRandom::rnd3(1.0, 1.0, 1.0);
      }
      // pick the lightest quark pair
      else if(diQuarkSelection_==1) {
	Energy m12 = (vec[0]->momentum()+vec[1]->momentum()).m();
	Energy m13 = (vec[0]->momentum()+vec[2]->momentum()).m();
	Energy m23 = (vec[1]->momentum()+vec[2]->momentum()).m();
	if     (m13<=m12&&m13<=m23)  ichoice = 2;
	else if(m23<=m12&&m23<=m13)  ichoice = 1;
      }
      else
	assert(false);
      // make the swaps so select pair first
      switch (ichoice) {
      case 0:
	break;
      case 1:
	swap(vec[2],vec[0]);
	break;
      case 2:
	swap(vec[2],vec[1]);
	break;
      }
    }
    // set up
    tcPDPtr temp1  = vec[0]->dataPtr();
    tcPDPtr temp2  = vec[1]->dataPtr();
    if(!other) other = vec[2];

    tcPDPtr dataDiquark  = CheckId::makeDiquark(temp1,temp2);
    
    if(!dataDiquark) 
      throw Exception() << "Could not make a diquark from"
			<< temp1->PDGName() << " and "
			<< temp2->PDGName() 
			<< " in ClusterFinder::reduceToTwoComponents()"
			<< Exception::eventerror;
   

    // Create the new cluster (with two components) and assign to it the same
    // momentum and position of the original (with three components) one.
    // Furthermore, assign to the diquark component a momentum given by the
    // sum of the two original components from which has been formed; for the
    // position, we are assuming, very simply, that the diquark position is
    // the average positions of the two original components.
    // Notice that the mass (5-th component of the 5-momentum) of the diquark
    // is set by hand to the constituent mass of the diquark (which is equal
    // to the sum of the constituent masses of the two quarks which form the
    // diquark) because the sum of 5-component vectors do add only the "normal"
    // 4-components, not the 5-th one. After that, the 5-momentum of the diquark
    // is in an inconsistent state, because the mass (5-th component) is not
    // equal to the invariant mass obtained from the 4-momemtum. This is not
    // unique to this kind of component (all perturbative components are in
    // a similar situation), but it is not harmful.
    
    // construct the diquark
    PPtr diquark = dataDiquark->produceParticle();
    vec[0]->addChild(diquark);
    vec[1]->addChild(diquark);
    diquark->set5Momentum(Lorentz5Momentum(vec[0]->momentum() + vec[1]->momentum(),
					   dataDiquark->constituentMass()));
    // use the same method as for cluster to determine the diquark position
    diquark->setVertex(Cluster::calculateX(vec[0],vec[1]));
    // put on-shell if required
    if(diQuarkOnShell_) {
      Lorentz5Momentum psum = diquark->momentum()+other->momentum();
      psum.rescaleMass();
      Boost boost = psum.boostVector();
      Lorentz5Momentum pother   =   other->momentum();
      Lorentz5Momentum pdiquark = diquark->momentum();
      pother.boost(-boost);
      pdiquark.boost(-boost);
      Energy pcm = Kinematics::pstarTwoBodyDecay(psum.mass(),
						 other->dataPtr()->constituentMass(),
						 diquark->dataPtr()->constituentMass());
      if(pcm>ZERO) {
	double fact = pcm/pother.vect().mag();
	pother   *= fact;
	pdiquark *= fact; 
	pother  .setMass(other->dataPtr()->constituentMass());
	pdiquark.setMass(dataDiquark     ->constituentMass());
	pother  .rescaleEnergy();
	pdiquark.rescaleEnergy();
	pother  .boost(boost);
	pdiquark.boost(boost);
	other->set5Momentum(pother);
	diquark->set5Momentum(pdiquark);
      }
    }
    // make the new cluster
    ClusterPtr nclus = new_ptr(Cluster(other,diquark));
    //vec[0]->addChild(nclus);
    //diquark->addChild(nclus);

    // Set the parent/children relationship between the original cluster 
    // (the one with three components) with the new one (the one with two components)
    // and add the latter to the vector of new redefined clusters.
    (*cluIter)->addChild(nclus);

    redefinedClusters.push_back(nclus);
  }

  // Add to  collecCluPtr  all of the redefined new clusters (indeed the 
  // pointers to them are added) contained in  vecNewRedefinedCluPtr.
  /// \todo why do we keep the original of the redefined clusters?
  for (tClusterVector::const_iterator it = redefinedClusters.begin();
        it != redefinedClusters.end(); ++it) {
    clusters.push_back(*it);
  }

}
