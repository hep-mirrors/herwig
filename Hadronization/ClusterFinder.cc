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
#include "Herwig++/Utilities/HwDebug.h"
#include "Herwig++/Utilities/CheckId.h"
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
  for(tParticleSet::const_iterator pit=origremnantSet.begin();pit!=origremnantSet.end();
      ++pit)
    {if((**pit).id()==ExtraParticleID::Remnant)
	{remnantSet.insert((**pit).children().back());}}
  //vector<pair<tClusterPtr,tParticleVector> > parent_child;
  set<tPPtr> examinedSet;  // colour particles already included in a cluster
  map<tColinePtr, pair<tPPtr,tPPtr> > quarkQuark; // quark quark 
  map<tColinePtr, pair<tPPtr,tPPtr> > aQuarkQuark; // anti quark anti quark
  ParticleSet inputParticles = pstep->particles();
  // Loop over all current particles.
  //for(ParticleSet::const_iterator pit = inputParticles.begin(); 
  //    pit != inputParticles.end(); ++pit) {
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
      if(pstep->antiColourNeighbour(*pit)) {
	connected[iElement++]= pstep->antiColourNeighbour( *pit );
      } else {       // colour source : baryon-violating process
        if((*pit)->colourLine()->sourceNeighbours() != tColinePair()) {
	  tColinePair sourcePair = (*pit)->colourLine()->sourceNeighbours();
	  tColinePtr intCL = tColinePtr();
	  for(int i = 0; i < 2; ++i) {
	    tColinePtr pLine = sourcePair.first; 
	    if(i == 1) pLine = sourcePair.second;
            int saveNumElements = iElement;
	    for(tPVector::const_iterator cit = pLine->coloured().begin(); 
		cit != pLine->coloured().end(); ++cit )
	      for(ParticleSet::const_iterator cjt = inputParticles.begin();
		   cjt != inputParticles.end(); ++cjt)
		if((*cit) == (*cjt)) connected[iElement++]= (*cit);      
	    if(iElement == saveNumElements) intCL = pLine;
	  }
	  if(intCL && iElement == 2) {
	    specialCase = true;
	    pair<tPPtr,tPPtr> qp=pair<tPPtr,tPPtr>(connected[0],connected[1]);
	    quarkQuark.insert(pair<tColinePtr,pair<tPPtr,tPPtr> >(intCL,qp)); 
	  }
	} 
      }
    }

    if((*pit)->hasAntiColour()) {
      if(pstep->colourNeighbour(*pit)) {
	connected[iElement++]= pstep->colourNeighbour(*pit);
      } else {       // colour sink : baryon-violating process
        if((*pit)->antiColourLine()->sinkNeighbours() != tColinePair()) {
	  tColinePair sinkPair = (*pit)->antiColourLine()->sinkNeighbours();
	  tColinePtr intCL = tColinePtr();
	  for(int i = 0; i < 2; ++i) {
	    tColinePtr pLine = sinkPair.first; 
	    if(i == 1) pLine = sinkPair.second;
            int saveNumElements = iElement;
	    for(tPVector::const_iterator cit = pLine->antiColoured().begin();
		cit != pLine->antiColoured().end(); ++cit) {
	      for(ParticleSet::const_iterator cjt = inputParticles.begin();
		  cjt != inputParticles.end(); ++cjt) {
		if((*cit) == (*cjt)) connected[iElement++]= (*cit);      
	      }
	    }
	    if(iElement == saveNumElements) intCL = pLine;
	  }
	  if(intCL && iElement == 2) {
	    specialCase = true;
	    pair<tPPtr,tPPtr> aqp=pair<tPPtr,tPPtr>(connected[0],connected[1]);
	    aQuarkQuark.insert(pair<tColinePtr,pair<tPPtr,tPPtr> >(intCL,aqp));
	  }
	} 
      }
    }
    if(HERWIG_DEBUG_LEVEL >= HwDebug::extreme) {
       generator()->log() << "Got to here!\n";
    }
    // Sanity checks (normally skipped) to see if the found cluster is 
    // sensible, that is: it has at least 2 components; none of its components
    // is a gluon (because after the gluon splitting step no gluons should be
    // present); all components have consistent colour connections (that is it
    // has colour or anticolour but not both; if it has colour, it must have a
    // colour line but neither have a anticolour line nor a colour neighbour; 
    // viceversa, if it has anticolour, it must have a anticolour line but 
    // neither have a colour line nor a anticolour neighbour); it has not 
    // components already examined (i.e. in examinedSet, which would imply 
    // that such component belong already to another cluster); the cluster must
    // colourless, that is it is a "closed" set under the operation of 
    // "colour/anticolour neighbour"; finally, in the case of sink or source 
    // of colour, but not in the special case (an number of pair of baryon-
    // violating processes which overall baryon number conservation), the sink
    // or source is consistent (that is, if the particle has colour it must 
    // stem from a colour source and it should not end into a colour sink; 
    // viceversa, if the particle has anticolour, it must end into a colour 
    // sink and it should not stem from a colour source). 
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::minimal_Hadronization ) {
      if (iElement < 2) {
	generator()->logWarning(Exception("ClusterFinder::formClusters "
			     "***Cluster with less than 2 components***",
					   Exception::warning) );
	if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization ) {    
	  generator()->log() << "         ===>" << " num elements=" 
			     << iElement << endl << endl;
	}
      }
      for (int i=0; i<iElement; ++i) {
	if ( connected[i]->id() == ThePEG::ParticleID::g ) {
	  generator()->logWarning(Exception("ClusterFinder::formClusters "
				"***Still gluons after gluon splitting***",
					     Exception::warning) );
	} 
      }
#define ci (connected[i])

      for (int i=1; i<iElement; ++i) { 
        if ((!ci->hasColour() &&  !ci->hasAntiColour()) ||
	    (ci->hasColour() && ci->hasAntiColour()) ||
	    (!ci->colourLine() && !ci->antiColourLine()) ||
	    (ci->colourLine() && ci->antiColourLine()) ||
	    (ci->hasColour() && !ci->colourLine()) ||
	    (ci->hasColour() && ci->antiColourLine()) ||
	    (ci->hasAntiColour() && ! ci->antiColourLine()) ||
	    (ci->hasAntiColour() && ci->colourLine()) ||
	    (ci->colourNeighbour() && ci->antiColourNeighbour()) ||
	    (ci->hasColour() && ci->colourNeighbour()) ||
	    (ci->hasAntiColour() && ci->antiColourNeighbour())) {
	  generator()->logWarning(Exception("ClusterFinder::formClusters "
			      "***Inconsistent (anti)Colour Connection***",
					    Exception::warning) );	  
	} else {
	  if((ci->hasColour() && pstep->antiColourNeighbour(ci))
	     || (ci->hasAntiColour() && pstep->colourNeighbour(ci))) {
	    tPPtr c = pstep->colourNeighbour(ci);
	    tPPtr antic = pstep->antiColourNeighbour(ci);
	    if(examinedSet.find(c) != examinedSet.end() ||  
	       examinedSet.find(antic) != examinedSet.end()) {
	      generator()->logWarning(Exception("ClusterFinder::formClusters "
				      "***Found component already examined***",
						Exception::warning) );   	      
	    }
	    bool cFound = false;
	    bool anticFound = false;
	    if ( ! c ) cFound = true;
	    if ( ! antic ) anticFound = true;
	    for (int j=0; j<iElement; ++j) {
	      if (connected[j] == c) cFound = true;
	      if (connected[j] == antic) anticFound = true; 
	    }
	    if ( ! cFound  ||  ! anticFound ) {
	      generator()->logWarning(Exception("ClusterFinder::formClusters "
					     "***Found cluster with colour***",
						Exception::warning) );
	    }
	  } else {  // there is a sink or a source of colour

#define cli (connected[i]->colourLine())
#define acli (connected[i]->antiColourLine())

	    if((cli && cli->sourceNeighbours() == tColinePair()) || 
	       (cli && cli->sinkNeighbours() != tColinePair()) || 
	       (acli && acli->sinkNeighbours() == tColinePair()) || 
	       (acli && acli->sourceNeighbours() != tColinePair())) { 
	      generator()->logWarning(Exception("ClusterFinder::formClusters "
				    "***Inconsistent colour sink or source***",
						Exception::warning) );	    
	    }
	  }
	}
      }
    }  // end sanity checks

    if(HERWIG_DEBUG_LEVEL >= HwDebug::extreme) {
      generator()->log() << "Now past the sanity checks\n";
    }
#undef ci
#undef cli
#undef acli

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
  if(HERWIG_DEBUG_LEVEL >= HwDebug::minimal_Hadronization) {
    if(quarkQuark.size() != aQuarkQuark.size()) {
      generator()->logWarning(Exception("ClusterFinder::formClusters "
		    "***The two special maps have different sizes!***",
					Exception::warning));
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization ) {    
	generator()->log() << "         ===>" << " size QQ = " 
			   << quarkQuark.size() << "  size QbarQbar = " 
			   << aQuarkQuark.size() << endl << endl;
      }
    } else {
      for(map<tColinePtr, pair<tPPtr,tPPtr> >::const_iterator 
	    cit = quarkQuark.begin(); cit != quarkQuark.end(); ++cit) {
	if(aQuarkQuark.find(cit->first) == aQuarkQuark.end()) {
	  generator()->logWarning(Exception("ClusterFinder::formClusters "
		       "***Not found special antiQuark-antiQuark pair ***",
					    Exception::warning) );
	  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization ) {    
	    generator()->log() << "         ===>" << " Q Q " 
			       << cit->second.first->PDGName() 
                               << "  " << cit->second.second->PDGName()
			       << endl << endl;
	  }
	}
      }
      for(map<tColinePtr, pair<tPPtr,tPPtr> >::const_iterator 
	  cit=aQuarkQuark.begin(); cit!=aQuarkQuark.end(); ++cit) {
	if(quarkQuark.find(cit->first) == quarkQuark.end()) {
	  generator()->logWarning( Exception("ClusterFinder::formClusters "
					     "***Not found special quark-quark pair ***",
					     Exception::warning) );
	  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization ) {    
	    generator()->log() << "         ===>" << " Q Q " << cit->second.first->PDGName() 
                               << "  " << cit->second.second->PDGName()
			       << endl << endl;
	  }
	}
      }
    }
  } // end sanity checks

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
  
  // Debugging.
  if ( HERWIG_DEBUG_LEVEL == HwDebug::extreme_Hadronization ) {    
    generator()->log() << "ClusterFinder::formClusters : "
		       << "*** extreme debugging ***" << endl
		       << " \t num clusters=" << clusters.size() << endl
		       << " \t --- Cluster information : begin --- " << endl;
    int i=1;
    for (ClusterVector::const_iterator it = clusters.begin();
	 it != clusters.end(); ++it) {
      generator()->log() << " \t Cluster number " 
			 << (*it)->number() << " ("
			 << i++ << "th)" << endl 
			 << " \t \t numComponents=" << (*it)->numComponents()
			 << " sumConstituentMasses=" 
			 << (*it)->sumConstituentMasses()
			 << " mass=" << (*it)->mass() 
	                 << ((*it)->isBeamCluster() ? "  beamCluster " : " ") 
			 << endl << " \t \t component ids "; 
 
      for(int i = 0; i<(*it)->numComponents(); i++) {
	generator()->log() << (*it)->particle(i)->id()
		     << ((*it)->isBeamRemnant(i) ? " beamRemnant  " : "  ");
      }
      generator()->log() << endl;
    }
    generator()->log() << " \t --- Cluster information : end --- "  << endl 
		       << endl;
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
    
    //    tCompPtr vec[3];  // store here the pointers to the three components.
    //  int i = 0;
    for(int i = 0; i<(*cluIter)->numComponents(); i++)
	  vec[i] = (*cluIter)->particle(i);
    
    // Sanity check (normally skipped) to see if the three components are
    // consistent: i.e. they can form a baryon and all components are pointing 
    // to a corresponding particle.
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::minimal_Hadronization ) {
      if ( ! CheckId::canBeBaryon(vec[0]->id(), vec[1]->id(), vec[2]->id()) ||
	   ! vec[0]  ||  ! vec[1] ||   ! vec[2]) {
	generator()->logWarning( Exception("ClusterFinder::reduceToTwoComponents "
					   "***Cluster with 3 inconsistent components***", 
					   Exception::warning) );
	if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization ) {    
	  generator()->log() << "         ===>" 
			     << " id1=" << vec[0]->id() << " id2=" << vec[1]->id() 
			     << " id3=" << vec[2]->id() << endl << endl;
	}
      }
    }
    
    // Randomly selects two components to be considered as a (anti)diquark
    // and place them as the second and third element of  vec.
    double prob0=1.0/3.0, prob1=1.0/3.0, prob2=1.0/3.0;
    int choice = UseRandom::rnd3(prob0, prob1, prob2);
    switch (choice) {
    case 0: 
      break; 
    case 1: { 
      swap(vec[0],vec[1]);
      break;
    } 
    case 2: {
      swap(vec[0],vec[2]);
      break;
    } 
    }
    
    long idDiquark = CheckId::diquarkId( vec[1]->id(), vec[2]->id() );
    
    // Sanity check (normally skipped) to see if we really got a diquark.
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::minimal_Hadronization ) {
      if ( ! CheckId::isDiquark(idDiquark) ) {
	generator()->logWarning( Exception("ClusterFinder::reduceToTwoComponents "
					   "***Inconsistent Diquark***",
					   Exception::warning) );
	if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization ) {    
	  generator()->log() << "         ===>" << " id=" << idDiquark << endl << endl;
	} 
      }
    }
    
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

    // Debugging.
    if ( HERWIG_DEBUG_LEVEL == HwDebug::extreme_Hadronization ) {    
      generator()->log() << " ClusterFinder::reduceToTwoComponents : *** extreme debugging ***" << endl
			 << " \t 3 component ids: " << vec[0]->id() << " " 
			 << vec[1]->id() << " " << vec[2]->id() << " "
			 << " final 2 component ids: " << vec[0]->id() << " " << idDiquark
			 << endl << endl;
    }
    
  }

  // Add to  collecCluPtr  all of the redefined new clusters (indeed the 
  // pointers to them are added) contained in  vecNewRedefinedCluPtr.
  for (tClusterVector::const_iterator it = redefinedClusters.begin();
        it != redefinedClusters.end(); ++it) {
    clusters.push_back(*it);
  }

}
