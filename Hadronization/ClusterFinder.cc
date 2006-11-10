// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ClusterFinder class.
//

#include "ClusterFinder.h"
#include <ThePEG/Interface/ClassDocumentation.h>
#include <ThePEG/Persistency/PersistentOStream.h>
#include <ThePEG/Persistency/PersistentIStream.h>
#include <ThePEG/PDT/StandardMatchers.h>
#include <ThePEG/PDT/EnumParticles.h>
#include <ThePEG/Repository/EventGenerator.h>
#include <ThePEG/EventRecord/Collision.h>
#include "CheckId.h"
#include "Herwig++/Utilities/EnumParticles.h"
#include "Cluster.h"

using namespace Herwig;

NoPIOClassDescription<ClusterFinder> ClusterFinder::initClusterFinder;
// Definition of the static class description member.

void ClusterFinder::Init() {

  static ClassDocumentation<ClusterFinder> documentation
    ("This class is responsible of finding clusters.");

}


void ClusterFinder::formClusters(tCollPtr collisionPtr, const StepPtr & pstep,
				 tPVector partons,
				 ClusterVector& clusters) 
  throw(Veto, Stop, Exception) {
  // Get all the beam remnant partons ((anti-)quark or (anti-)diquark) 
  // of the current collision.
  tParticleSet origremnantSet = collisionPtr->getRemnants(),remnantSet;
  for(tParticleSet::const_iterator pit=origremnantSet.begin();
      pit!=origremnantSet.end();++pit) {
    if((**pit).id()==ExtraParticleID::Remnant)
      remnantSet.insert((**pit).children().back());
  }
  set<tPPtr> examinedSet;  // colour particles already included in a cluster
  map<tColinePtr, pair<tPPtr,tPPtr> > quarkQuark; // quark quark 
  map<tColinePtr, pair<tPPtr,tPPtr> > aQuarkQuark; // anti quark anti quark
  ParticleSet inputParticles = pstep->particles();
  // Loop over all current particles.
  for(tPVector::const_iterator pit=partons.begin();pit!=partons.end();++pit){
    // Skip to the next particle if it is not coloured or already examined.
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
    //       tilda_u_R -> dbar_1 + dbar_2 
    //       tilda_u_R_star -> d1 + d2 
    //    where tilda_u_R and tilda_u_R_star are colour connected )
    // a special treatment is needed, because first we have to process all
    // partons in the current step, and then for each left pair of quarks which
    // stem from a colour source we have to find the corresponding pair of 
    // anti-quarks which ends in a colour sink and is connected with the
    // above colour source. These special pairs are kept into the maps:
    //    specialQuarkQuarkMap   and   specialAntiQuarkAntiQuarkMap.

    tParticleVector connected(3);
    int iElement = 0;         
    connected[iElement++] = *pit;
    bool specialCase = false;

    if((*pit)->hasColour()) {
      tPPtr partner=pstep->antiColourNeighbour(*pit);
      if(partner) {
	connected[iElement++]= partner;
      } 
      else {       // colour source : baryon-violating process

	// this is not supported yet exit
	throw Exception() << "Baryon number violating colour flows not yet implemented"
			  << " if you are not trying to simulation a BNV process some"
			  << "thing is probably wrong with the colour connections"
			  << Exception::runerror;


//         if((*pit)->colourLine()->sourceNeighbours() != tColinePair()) {
// 	  tColinePair sourcePair = (*pit)->colourLine()->sourceNeighbours();
// 	  tColinePtr intCL = tColinePtr();
// 	  for(int i = 0; i < 2; ++i) {
// 	    tColinePtr pLine = sourcePair.first; 
// 	    if(i == 1) pLine = sourcePair.second;
//             int saveNumElements = iElement;
// 	    for(tPVector::const_iterator cit = pLine->coloured().begin(); 
// 		cit != pLine->coloured().end(); ++cit )
// 	      for(ParticleSet::const_iterator cjt = inputParticles.begin();
// 		   cjt != inputParticles.end(); ++cjt)
// 		if((*cit) == (*cjt)) connected[iElement++]= (*cit);      
// 	    if(iElement == saveNumElements) intCL = pLine;
// 	  }
// 	  if(intCL && iElement == 2) {
// 	    specialCase = true;
// 	    pair<tPPtr,tPPtr> qp=pair<tPPtr,tPPtr>(connected[0],connected[1]);
// 	    quarkQuark.insert(pair<tColinePtr,pair<tPPtr,tPPtr> >(intCL,qp)); 
// 	  }
//      }	
      }
    }

    if((*pit)->hasAntiColour()) {
      tPPtr partner=pstep->colourNeighbour(*pit);
      if(partner) {
	connected[iElement++]=partner; 
      } 
      else {       // colour sink : baryon-violating process

	// this is not supported yet exit
	throw Exception() << "Baryon number violating colour flows not yet implemented"
			  << " if you are not trying to simulation a BNV process some"
			  << "thing is probably wrong with the colour connections"
			  << Exception::runerror;
//         if((*pit)->antiColourLine()->sinkNeighbours() != tColinePair()) {
// 	  tColinePair sinkPair = (*pit)->antiColourLine()->sinkNeighbours();
// 	  tColinePtr intCL = tColinePtr();
// 	  for(int i = 0; i < 2; ++i) {
// 	    tColinePtr pLine = sinkPair.first; 
// 	    if(i == 1) pLine = sinkPair.second;
//             int saveNumElements = iElement;
// 	    for(tPVector::const_iterator cit = pLine->antiColoured().begin();
// 		cit != pLine->antiColoured().end(); ++cit) {
// 	      for(ParticleSet::const_iterator cjt = inputParticles.begin();
// 		  cjt != inputParticles.end(); ++cjt) {
// 		if((*cit) == (*cjt)) connected[iElement++]= (*cit);      
// 	      }
// 	    }
// 	    if(iElement == saveNumElements) intCL = pLine;
// 	  }
// 	  if(intCL && iElement == 2) {
// 	    specialCase = true;
// 	    pair<tPPtr,tPPtr> aqp=pair<tPPtr,tPPtr>(connected[0],connected[1]);
// 	    aQuarkQuark.insert(pair<tColinePtr,pair<tPPtr,tPPtr> >(intCL,aqp));
// 	  }
// 	} 
      }
    }

    if(!specialCase) {
      // Tag the components of the found cluster as already examined.
      for (int i=0; i<iElement; ++i) {
	examinedSet.insert(connected[i]);        
      }

      // Create the cluster object with the colour connected particles
      ClusterPtr cluPtr = new_ptr(Cluster(connected[0],
					  connected[1],
					  connected[2]));

      clusters.push_back(cluPtr);

      pstep->addDecayProduct(connected[0], cluPtr);
      pstep->addDecayProduct(connected[1], cluPtr);
      if(connected[2]) pstep->addDecayProduct(connected[2], cluPtr);
   

      // Check if any of the components is a beam remnant, and if this
      // is the case then inform the cluster.   
      // this will only work for baryon collisions  
      for (int i=0; i<iElement; ++i) {
	if(connected[i]->parents()[0]->id()==ExtraParticleID::Remnant&&
	   DiquarkMatcher::Check(connected[i]->id()))
	  cluPtr->isBeamCluster(connected[i]);
      }
    }

  } // end main loop over particles in the step.   
  // Sanity checks (normally skipped) to see if the two special maps 
  // are consistent, that is they have the same size, and for each 
  // quark-quark pair there is a corresponding antiQuark-antiQuark pair,
  // and, viceversa, for each antiQuark-antiQuark pair there is a 
  // corresponding quark-quark pair.

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
	pstep->addDecayProduct(quarkPair.first, cluPtr1);
	pstep->addDecayProduct(antiQuarkPair.first, cluPtr1);
	pstep->addDecayProduct(quarkPair.second, cluPtr2);
	pstep->addDecayProduct(antiQuarkPair.second, cluPtr2);
      } else {
	cluPtr1 = new_ptr(Cluster(quarkPair.first , antiQuarkPair.second));
	cluPtr2 = new_ptr(Cluster(quarkPair.second , antiQuarkPair.first));
	pstep->addDecayProduct(quarkPair.second, cluPtr1);
	pstep->addDecayProduct(antiQuarkPair.first, cluPtr1);
	pstep->addDecayProduct(quarkPair.first, cluPtr2);
	pstep->addDecayProduct(antiQuarkPair.second, cluPtr2);
      }
      clusters.push_back(cluPtr1);
      clusters.push_back(cluPtr2);
    } else {
      throw Exception("ClusterFinder::formClusters : " 
      "***Skip event: unable to match pairs in Baryon-violating processes***",
		      Exception::eventerror);
    }
  }
}


void ClusterFinder::reduceToTwoComponents(const StepPtr & pstep,
					  ClusterVector & clusters) 
  throw(Veto, Stop, Exception) {

  // In order to preserve all of the information, we do not modify the 
  // directly the 3-component clusters, but instead we define new clusters,
  // which are related to the original ones by a child-parent relationship,
  // by considering two randomly chosen components as a diquark (or anti-diquark).
  // These new clusters are first added to the vector  vecNewRedefinedCluPtr,
  // and at the end, when all input clusters have been examined, the elements of 
  // this vector will be copied in  collecCluPtr  (the reason is that it is not 
  // allowed to modify a STL container while iterating over it).
  vector<tClusterPtr> redefinedClusters; 
  ParticleVector vec(3);
  for(ClusterVector::iterator cluIter = clusters.begin() ; 
      cluIter != clusters.end() ; ++cluIter) {

    if ( ! (*cluIter)->isAvailable()  ||  (*cluIter)->numComponents() != 3 ) continue;
    
    for(int i = 0; i<(*cluIter)->numComponents(); i++)
	  vec[i] = (*cluIter)->particle(i);

    
    // Randomly selects two components to be considered as a (anti)diquark
    // and place them as the second and third element of  vec.
    double prob0=1.0/3.0, prob1=1.0/3.0, prob2=1.0/3.0;
    int choice = UseRandom::rnd3(prob0, prob1, prob2);
    switch (choice) {
    case 0: 
      break; 
    case 1:
      swap(vec[0],vec[1]);
      break;
    case 2:
      swap(vec[0],vec[2]);
      break;
    }
    
    long idDiquark = CheckId::diquarkId( vec[1]->id(), vec[2]->id() );
   
    
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
    
    PPtr diquark = getParticle(idDiquark);
    pstep->addDecayProduct(vec.begin()+1, vec.end(), diquark);
    ClusterPtr nclus = new_ptr(Cluster(vec[0],diquark));
    tParticleVector pv;
    pv.push_back(vec[0]);
    pv.push_back(diquark);
    pstep->addDecayProduct(pv.begin(), pv.end(), nclus);
    nclus->set5Momentum((*cluIter)->momentum());
    nclus->setVertex((*cluIter)->vertex());
    for(int i = 0; i<nclus->numComponents(); i++) {
      if(nclus->particle(i)->id() == idDiquark) {
	nclus->particle(i)->set5Momentum(Lorentz5Momentum(vec[1]->momentum()
	                     + vec[2]->momentum(), getParticleData(idDiquark)
			     ->constituentMass()));
        nclus->particle(i)->setVertex(0.5*(vec[1]->vertex() 
			     + vec[2]->vertex()));
      }
    }
    // Set the parent/children relationship between the original cluster 
    // (the one with three components) with the new one (the one with two components)
    // and add the latter to the vector of new redefined clusters.
    //(*cluIter)->addChild(nclus);
    redefinedClusters.push_back(nclus);
  }

  // Add to  collecCluPtr  all of the redefined new clusters (indeed the 
  // pointers to them are added) contained in  vecNewRedefinedCluPtr.
  for (tClusterVector::const_iterator it = redefinedClusters.begin();
        it != redefinedClusters.end(); ++it) {
    clusters.push_back(*it);
  }

}
