// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ClusterFinder class.
//

#include "ClusterFinder.h"
#include "Pythia7/Interface/ClassDocumentation.h"
#include "Pythia7/Persistency/PersistentOStream.h"
#include "Pythia7/Persistency/PersistentIStream.h"
#include "Pythia7/Handlers/PartialCollisionHandler.h"
#include "Pythia7/PDT/EnumParticles.h"
#include "Pythia7/Repository/EventGenerator.h"
#include "Pythia7/EventRecord/Collision.h"
#include "Herwig++/Utilities/HwDebug.h"
#include "Herwig++/Utilities/CheckId.h"
#include "Cluster.h"
#include "Component.h"


using namespace Herwig;
// using namespace Pythia7;


ClusterFinder::~ClusterFinder() {}


void ClusterFinder::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}


void ClusterFinder::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}


ClassDescription<ClusterFinder> ClusterFinder::initClusterFinder;
// Definition of the static class description member.

void ClusterFinder::Init() {

  static ClassDocumentation<ClusterFinder> documentation
    ("This class is responsible of finding clusters.");

}


void ClusterFinder::formClusters(tCollPtr collisionPtr, const StepPtr & pstep, 
				 CollecCluPtr & collecCluPtr) throw(Veto, Stop, Exception) {

  // Get all the beam remnant partons ((anti-)quark or (anti-)diquark) 
  // of the current collision.
  tParticleSet remnantSet = collisionPtr->getRemnants();
 
  set<tPPtr> examinedSet;  // colour particles already included in a cluster
  map<tColinePtr, pair<tPPtr,tPPtr> > specialQuarkQuarkMap;          // see comment below
  map<tColinePtr, pair<tPPtr,tPPtr> > specialAntiQuarkAntiQuarkMap;  //  "     "      "
  
  // Loop over all current particles.
  for( ParticleSet::const_iterator pit = pstep->particles().begin(); 
       pit != pstep->particles().end(); ++pit) {

    // Skip to the next particle if it is not coloured or already examined.
    if ( ! (**pit).data().coloured()  ||
	 examinedSet.find(*pit) != examinedSet.end() ) continue;  
     
    // We assume that a cluster can be made of, at most, 3 constituents although
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

    tPPtr connectedVec[3] = { tPPtr(), tPPtr(), tPPtr() };   
    int iElement = 0;         
    connectedVec[iElement++] = *pit;
    bool specialCase = false;

    if ( (*pit)->hasColour() ) {
      if ( pstep->antiColourNeighbour( *pit ) ) {
	connectedVec[iElement++]= pstep->antiColourNeighbour( *pit );
      } else {       // colour source : baryon-violating process
        if ( (*pit)->colourLine()->sourceNeighbours() != tColinePair() ) {
	  tColinePair sourcePair = (*pit)->colourLine()->sourceNeighbours();
	  tColinePtr internalColine = tColinePtr();
	  for (int i = 0; i < 2; ++i ) {
	    tColinePtr pLine = sourcePair.first; 
	    if ( i == 1 ) pLine = sourcePair.second;
            int saveNumElements = iElement;
	    for ( tPVector::const_iterator cit = pLine->coloured().begin();
		  cit != pLine->coloured().end(); ++cit ) {
	      for( ParticleSet::const_iterator cjt = pstep->particles().begin(); 
		   cjt != pstep->particles().end(); ++cjt) {
		if ( (*cit) == (*cjt) ) connectedVec[iElement++]= (*cit);      
	      }
	    }
	    if ( iElement == saveNumElements ) internalColine = pLine;
	  }
	  if ( internalColine  &&  iElement == 2 ) {
	    specialCase = true;
	    pair<tPPtr,tPPtr> quarkPair = pair<tPPtr,tPPtr>(connectedVec[0],connectedVec[1]);
	    specialQuarkQuarkMap.insert( pair<tColinePtr,pair<tPPtr,tPPtr> >(internalColine,quarkPair) ); 
	  }
	} 
      }
    }

    if ( (*pit)->hasAntiColour() ) {
      if ( pstep->colourNeighbour( *pit ) ) {
	connectedVec[iElement++]= pstep->colourNeighbour( *pit );
      } else {       // colour sink : baryon-violating process
        if ( (*pit)->antiColourLine()->sinkNeighbours() != tColinePair() ) {
	  tColinePair sinkPair = (*pit)->antiColourLine()->sinkNeighbours();
	  tColinePtr internalColine = tColinePtr();
	  for (int i = 0; i < 2; ++i ) {
	    tColinePtr pLine = sinkPair.first; 
	    if ( i == 1 ) pLine = sinkPair.second;
            int saveNumElements = iElement;
	    for ( tPVector::const_iterator cit = pLine->antiColoured().begin();
		  cit != pLine->antiColoured().end(); ++cit ) {
	      for( ParticleSet::const_iterator cjt = pstep->particles().begin(); 
		   cjt != pstep->particles().end(); ++cjt) {
		if ( (*cit) == (*cjt) ) connectedVec[iElement++]= (*cit);      
	      }
	    }
	    if ( iElement == saveNumElements ) internalColine = pLine;
	  }
	  if ( internalColine  &&  iElement == 2 ) {
	    specialCase = true;
	    pair<tPPtr,tPPtr> antiQuarkPair = pair<tPPtr,tPPtr>(connectedVec[0],connectedVec[1]);
	    specialAntiQuarkAntiQuarkMap.insert( pair<tColinePtr,pair<tPPtr,tPPtr> >
						 (internalColine,antiQuarkPair) ); 
	  }
	} 
      }
    }

    // Sanity checks (normally skipped) to see if the found cluster is sensible, 
    // that is: it has at least 2 components; none of its components is a gluon
    // (because after the gluon splitting step no gluons should be present); 
    // all components have consistent colour connections (that is it has colour 
    // or anticolour but not both; if it has colour, it must have a colour line 
    // but neither have a anticolour line nor a colour neighbour; viceversa, 
    // if it has anticolour, it must have a anticolour line but neither have a 
    // colour line nor a anticolour neighbour); it has not components already 
    // examined (i.e. in examinedSet, which would imply that such component belong 
    // already to another cluster); the cluster must colourless, that is it is a 
    // "closed" set under the operation of "colour/anticolour neighbour"; finally, 
    // in the case of sink or source of colour, but not in the special case (an 
    // number of pair of baryon-violating processes which overall baryon number 
    // conservation), the sink or source is consistent (that is, if the particle 
    // has colour it must stem from a colour source and it should not end into a
    // colour sink; viceversa, if the particle has anticolour, it must end into
    // a colour sink and it should not stem from a colour source). 
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::minimal_Hadronization ) {
      if (iElement < 2) {
	generator()->logWarning( Exception("ClusterFinder::formClusters "
					   "***Cluster with less than 2 components***",
					   Exception::warning) );
	if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization ) {    
	  generator()->log() << "         ===>" << " num elements=" << iElement 
			     << endl << endl;
	}
      }
      for (int i=0; i<iElement; ++i) {
	if ( connectedVec[i]->id() == ParticleID::g ) {
	  generator()->logWarning( Exception("ClusterFinder::formClusters "
					     "***Still gluons after gluon splitting***",
					     Exception::warning) );
	} 
      }
      for (int i=1; i<iElement; ++i) { 
        if ( ( ! connectedVec[i]->hasColour()  &&  ! connectedVec[i]->hasAntiColour() ) ||
	     ( connectedVec[i]->hasColour() && connectedVec[i]->hasAntiColour() ) ||
             ( ! connectedVec[i]->colourLine()  &&  ! connectedVec[i]->antiColourLine() ) ||
             ( connectedVec[i]->colourLine() && connectedVec[i]->antiColourLine() ) ||
             ( connectedVec[i]->hasColour()  &&  ! connectedVec[i]->colourLine() ) ||
             ( connectedVec[i]->hasColour() && connectedVec[i]->antiColourLine() ) ||
             ( connectedVec[i]->hasAntiColour()  &&  ! connectedVec[i]->antiColourLine() ) ||
             ( connectedVec[i]->hasAntiColour() && connectedVec[i]->colourLine() ) ||
             ( connectedVec[i]->colourNeighbour() && connectedVec[i]->antiColourNeighbour() ) ||
             ( connectedVec[i]->hasColour() && connectedVec[i]->colourNeighbour() ) ||
             ( connectedVec[i]->hasAntiColour() && connectedVec[i]->antiColourNeighbour() ) ) {
	  generator()->logWarning( Exception("ClusterFinder::formClusters "
					     "***Inconsistent (anti)Colour Connection***",
					     Exception::warning) );	  
	} else {
	  if ( ( connectedVec[i]->hasColour() && pstep->antiColourNeighbour( connectedVec[i] ) ) ||
	       ( connectedVec[i]->hasAntiColour() && pstep->colourNeighbour( connectedVec[i] ) ) ) {
	    tPPtr c = pstep->colourNeighbour( connectedVec[i] );
	    tPPtr antic = pstep->antiColourNeighbour( connectedVec[i] );
	    if ( examinedSet.find(c) != examinedSet.end()  ||  
		 examinedSet.find(antic) != examinedSet.end() ) {
	      generator()->logWarning( Exception("ClusterFinder::formClusters "
						 "***Found component already examined***",
						 Exception::warning) );   	      
	    }
	    bool cFound = false;
	    bool anticFound = false;
	    if ( ! c ) cFound = true;
	    if ( ! antic ) anticFound = true;
	    for (int j=0; j<iElement; ++j) {
	      if (connectedVec[j] == c) cFound = true;
	      if (connectedVec[j] == antic) anticFound = true; 
	    }
	    if ( ! cFound  ||  ! anticFound ) {
	      generator()->logWarning( Exception("ClusterFinder::formClusters "
						 "***Found cluster with colour***",
						 Exception::warning) );
	    }
	  } else {  // there is a sink or a source of colour
	    if ( ( connectedVec[i]->colourLine()  && 
		   connectedVec[i]->colourLine()->sourceNeighbours() == tColinePair() ) || 
		 ( connectedVec[i]->colourLine()  && 
		   connectedVec[i]->colourLine()->sinkNeighbours() != tColinePair() ) || 
		 ( connectedVec[i]->antiColourLine()  && 
		   connectedVec[i]->antiColourLine()->sinkNeighbours() == tColinePair() ) || 
		 ( connectedVec[i]->antiColourLine()  && 
		   connectedVec[i]->antiColourLine()->sourceNeighbours() != tColinePair() ) ) { 
	      generator()->logWarning( Exception("ClusterFinder::formClusters "
						 "***Inconsistent colour sink or source***",
						 Exception::warning) );	    
	    }
	  }
	}
      }
    }  // end sanity checks

    if ( ! specialCase ) {
      // Tag the components of the found cluster as already examined.
      for (int i=0; i<iElement; ++i) {
	examinedSet.insert( connectedVec[i] );        
      }
      // Create the cluster object with the colour connected particles
      CluPtr cluPtr = new_ptr( Cluster(connectedVec[0],connectedVec[1],connectedVec[2]) );
      collecCluPtr.insert( collecCluPtr.end() , cluPtr);
      // Check if any of the components is a beam remnant, and if this
      // is the case then inform the cluster.    
      for (int i=0; i<iElement; ++i) {
	if ( remnantSet.find( connectedVec[i]->original() ) != remnantSet.end() ) {   
	  cluPtr->isBeamCluster( connectedVec[i] );
	}
      }
    }

  } // end main loop over particles in the step.   

  // Sanity checks (normally skipped) to see if the two special maps 
  // are consistent, that is they have the same size, and for each 
  // quark-quark pair there is a corresponding antiQuark-antiQuark pair,
  // and, viceversa, for each antiQuark-antiQuark pair there is a 
  // corresponding quark-quark pair.
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::minimal_Hadronization ) {
    if ( specialQuarkQuarkMap.size() != specialAntiQuarkAntiQuarkMap.size() ) {
      generator()->logWarning( Exception("ClusterFinder::formClusters "
					 "***The two special maps have different sizes!***",
					 Exception::warning) );
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization ) {    
	generator()->log() << "         ===>" << " size QQ = " << specialQuarkQuarkMap.size()
			   << "  size QbarQbar = " << specialAntiQuarkAntiQuarkMap.size()
			   << endl << endl;
      }
    } else {
      for ( map<tColinePtr, pair<tPPtr,tPPtr> >::const_iterator cit = specialQuarkQuarkMap.begin(); 
	    cit != specialQuarkQuarkMap.end(); ++cit ) {
	if ( specialAntiQuarkAntiQuarkMap.find( cit->first ) == specialAntiQuarkAntiQuarkMap.end() ) {
	  generator()->logWarning( Exception("ClusterFinder::formClusters "
					     "***Not found special antiQuark-antiQuark pair ***",
					     Exception::warning) );
	  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization ) {    
	    generator()->log() << "         ===>" << " Q Q " << cit->second.first->PDGName() 
                               << "  " << cit->second.second->PDGName()
			       << endl << endl;
	  }
	}
      }
      for ( map<tColinePtr, pair<tPPtr,tPPtr> >::const_iterator cit = specialAntiQuarkAntiQuarkMap.begin(); 
	    cit != specialAntiQuarkAntiQuarkMap.end(); ++cit ) {
	if ( specialQuarkQuarkMap.find( cit->first ) == specialQuarkQuarkMap.end() ) {
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
  for ( map<tColinePtr, pair<tPPtr,tPPtr> >::const_iterator cit = specialQuarkQuarkMap.begin(); 
	cit != specialQuarkQuarkMap.end(); ++cit ) {
    tColinePtr coline = cit->first;
    pair<tPPtr,tPPtr> quarkPair = cit->second;
    if ( specialAntiQuarkAntiQuarkMap.find( coline ) != specialAntiQuarkAntiQuarkMap.end() ) {
      pair<tPPtr,tPPtr> antiQuarkPair = specialAntiQuarkAntiQuarkMap.find( coline )->second;
      CluPtr cluPtr1, cluPtr2;
      if ( rndbool() ) {      
	cluPtr1 = new_ptr( Cluster( quarkPair.first , antiQuarkPair.first ) );
	cluPtr2 = new_ptr( Cluster( quarkPair.second , antiQuarkPair.second ) );
      } else {
	cluPtr1 = new_ptr( Cluster( quarkPair.first , antiQuarkPair.second ) );
	cluPtr2 = new_ptr( Cluster( quarkPair.second , antiQuarkPair.first ) );
      }
      collecCluPtr.insert( collecCluPtr.end() , cluPtr1);
      collecCluPtr.insert( collecCluPtr.end() , cluPtr2);
    } else {
      throw Exception("ClusterFinder::formClusters : " 
		      "***Skip event: unable to match pairs in Baryon-violating processes***",
		      Exception::eventerror);
    }
  }
  
  // Debugging.
  if ( HERWIG_DEBUG_LEVEL == HwDebug::extreme_Hadronization ) {    
    generator()->log() << "ClusterFinder::formClusters : *** extreme debugging ***" << endl
		       << " \t num clusters=" << collecCluPtr.size() << endl
		       << " \t --- Cluster information : begin --- " << endl;
    int i=1;
    for (CollecCluPtr::const_iterator it = collecCluPtr.begin();
	 it != collecCluPtr.end(); ++it) {
      generator()->log() << " \t Cluster number " << i++ << endl 
			 << " \t \t numComponents=" << (*it)->numComponents()
			 << " sumConstituentMasses=" << (*it)->sumConstituentMasses()
			 << " mass=" << (*it)->mass() 
	                 << ( (*it)->isBeamCluster() ? "  beamCluster " : " " ) << endl
                         << " \t \t component ids "; 
      for (CollecCompPtr::const_iterator jt = ((*it)->components()).begin();
	   jt != ((*it)->components()).end(); ++jt) {
	generator()->log() << (*jt)->id() 
			   << ( (*jt)->isBeamRemnant() ? " beamRemnant  " : "  " );
      }
      generator()->log() << endl;
    }
    generator()->log() << " \t --- Cluster information : end --- " << endl << endl;
  }

}


void ClusterFinder::reduceToTwoComponents(CollecCluPtr & collecCluPtr) 
  throw(Veto, Stop, Exception) {

  // In order to preserve all of the information, we do not modify the 
  // directly the 3-component clusters, but instead we define new clusters,
  // which are related to the original ones by a child-parent relationship,
  // by considering two randomly chosen components as a diquark (or anti-diquark).
  // These new clusters are first added to the vector  vecNewRedefinedCluPtr,
  // and at the end, when all input clusters have been examined, the elements of 
  // this vector will be copied in  collecCluPtr  (the reason is that it is not 
  // allowed to modify a STL container while iterating over it).
  vector<CluPtr> vecNewRedefinedCluPtr; 
  
  for (CollecCluPtr::iterator cluIter = collecCluPtr.begin() ; 
       cluIter != collecCluPtr.end() ; ++cluIter) {

    if ( ! (*cluIter)->isAvailable()  ||  (*cluIter)->numComponents() != 3 ) continue;
    
    tCompPtr vec[3];  // store here the pointers to the three components.
    int i = 0;
    for (CollecCompPtr::const_iterator iter = (*cluIter)->components().begin();
	 iter != (*cluIter)->components().end(); ++iter) {
      vec[i++] = *iter;
    }
    
    // Sanity check (normally skipped) to see if the three components are
    // consistent: i.e. they can form a baryon and all components are pointing 
    // to a corresponding particle.
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::minimal_Hadronization ) {
      if ( ! CheckId::canBeBaryon( vec[0]->id(), vec[1]->id(), vec[2]->id() ) ||  
	   ! vec[0]->pointerParticle()  ||  ! vec[1]->pointerParticle()  ||  
	   ! vec[2]->pointerParticle() ) {
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
    int choice = rnd3(prob0, prob1, prob2);
    switch (choice) {
    case 0: break; 
    case 1: { 
      tCompPtr save = vec[0]; 
      vec[0] = vec[1];
      vec[1] = save;
      break;
    } 
    case 2: {
      tCompPtr save = vec[0]; 
      vec[0] = vec[2];
      vec[2] = save;
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
    
    CluPtr newCluPtr = new_ptr( Cluster( vec[0]->pointerParticle() , idDiquark) );
    newCluPtr->momentum( (*cluIter)->momentum() );
    newCluPtr->position( (*cluIter)->position() );
    for (CollecCompPtr::const_iterator iter = newCluPtr->components().begin();
	 iter != newCluPtr->components().end(); ++iter) {
      if ( (*iter)->id() == idDiquark ) {
	(*iter)->momentum( Lorentz5Momentum( vec[1]->momentum() + vec[2]->momentum() ,
					     getParticleData( idDiquark )->constituentMass() ) );
	(*iter)->position( 0.5 * ( vec[1]->position() + vec[2]->position() ) );
      }
    }
    
    // Set the parent/children relationship between the original cluster 
    // (the one with three components) with the new one (the one with two components)
    // and add the latter to the vector of new redefined clusters.
    (*cluIter)->addChildrenClusters( newCluPtr );
    newCluPtr->parentCluster( *cluIter );
    vecNewRedefinedCluPtr.push_back( newCluPtr );    

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
  for (vector<CluPtr>::const_iterator it = vecNewRedefinedCluPtr.begin();
       it != vecNewRedefinedCluPtr.end(); ++it) {
    collecCluPtr.insert( collecCluPtr.end(), *it );
  }

}




