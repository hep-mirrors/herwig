// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ClusterHadronizationHandler class.
//

#include "ClusterHadronizationHandler.h"
#include "Pythia7/Interface/ClassDocumentation.h"
#include "Pythia7/Persistency/PersistentOStream.h"
#include "Pythia7/Persistency/PersistentIStream.h"
#include "Pythia7/Interface/Parameter.h" 
#include "Pythia7/Interface/Reference.h" 
#include "Pythia7/Handlers/CollisionHandler.h"
#include "Pythia7/Handlers/Hint.h"
#include "Pythia7/PDT/ParticleData.h"
#include "Pythia7/EventRecord/Particle.h"
#include "Pythia7/EventRecord/Step.h"
#include "Pythia7/PDT/PDT.h"   
#include "Pythia7/PDT/EnumParticles.h"
#include "Pythia7/Repository/EventGenerator.h"
#include "Herwig++/Utilities/HwDebug.h"
#include "Herwig++/Utilities/CheckId.h"
#include "Herwig++/Utilities/Smearing.h"
#include "CluHadConfig.h"
#include "Cluster.h" 
#include "Component.h" 


using namespace Herwig;
// using namespace Pythia7;


ClusterHadronizationHandler::~ClusterHadronizationHandler() {}


void ClusterHadronizationHandler::persistentOutput(PersistentOStream & os) const {
  os << _pointerGlobalParameters 
     << _pointerPartonSplitter 
     << _pointerClusterFinder
     << _pointerColourReconnector
     << _pointerClusterFissioner
     << _pointerLightClusterDecayer
     << _pointerClusterDecayer;
}


void ClusterHadronizationHandler::persistentInput(PersistentIStream & is, int) {
  is >> _pointerGlobalParameters 
     >> _pointerPartonSplitter 
     >> _pointerClusterFinder
     >> _pointerColourReconnector
     >> _pointerClusterFissioner
     >> _pointerLightClusterDecayer
     >> _pointerClusterDecayer;
}


ClassDescription<ClusterHadronizationHandler> ClusterHadronizationHandler::initClusterHadronizationHandler;
// Definition of the static class description member.


void ClusterHadronizationHandler::Init() {

  static ClassDocumentation<ClusterHadronizationHandler> documentation
    ("This is the main handler class for the Cluster Hadronization");

  static Reference<ClusterHadronizationHandler,GlobalParameters> 
    interfaceGlobalParameters("GlobalParameters", 
			      "A reference to the GlobalParameters object", 
			      &Herwig::ClusterHadronizationHandler::_pointerGlobalParameters,
			      false, false, true, false);
  static Reference<ClusterHadronizationHandler,PartonSplitter> 
    interfacePartonSplitter("PartonSplitter", 
			    "A reference to the PartonSplitter object", 
			    &Herwig::ClusterHadronizationHandler::_pointerPartonSplitter,
			    false, false, true, false);
  static Reference<ClusterHadronizationHandler,ClusterFinder> 
    interfaceClusterFinder("ClusterFinder", 
			   "A reference to the ClusterFinder object", 
			   &Herwig::ClusterHadronizationHandler::_pointerClusterFinder,
			   false, false, true, false);
  static Reference<ClusterHadronizationHandler,ColourReconnector> 
    interfaceColourReconnector("ColourReconnector", 
			       "A reference to the ColourReconnector object", 
			       &Herwig::ClusterHadronizationHandler::_pointerColourReconnector,
			       false, false, true, false);
  static Reference<ClusterHadronizationHandler,ClusterFissioner> 
    interfaceClusterFissioner("ClusterFissioner", 
			      "A reference to the ClusterFissioner object", 
			      &Herwig::ClusterHadronizationHandler::_pointerClusterFissioner,
			      false, false, true, false);
  static Reference<ClusterHadronizationHandler,LightClusterDecayer> 
    interfaceLightClusterDecayer("LightClusterDecayer", 
				 "A reference to the LightClusterDecayer object", 
				 &Herwig::ClusterHadronizationHandler::_pointerLightClusterDecayer,
				 false, false, true, false);
  static Reference<ClusterHadronizationHandler,ClusterDecayer> 
    interfaceClusterDecayer("ClusterDecayer", 
			    "A reference to the ClusterDecayer object", 
			    &Herwig::ClusterHadronizationHandler::_pointerClusterDecayer,
			    false, false, true, false);
  
}


void ClusterHadronizationHandler::doinitrun() {
  // The run initialization is used here only to allow the non-interfaced and
  // non-persistent classes Component and Cluster to have access to the 
  // GlobalParameters class instance, via a static pointer.
  Component trashComponentObject(1); // any Component instance would be fine.
  trashComponentObject.setPointerGlobalParameters( _pointerGlobalParameters );
  Cluster trashClusterObject(1,-2);  // any Cluster instance would be fine.
  trashClusterObject.setPointerGlobalParameters( _pointerGlobalParameters );
}


void ClusterHadronizationHandler::
handle(PartialCollisionHandler & ch, const tPVector & tagged,
       const Hint & hint) throw(Veto, Stop, Exception) {

  StepPtr pstep = ch.newStep();
 
  if ( HERWIG_DEBUG_LEVEL == HwDebug::extreme_Hadronization ) {   // Occasional debugging
    printStep(pstep,"At the beginning of ClusterHadronizationHandler");
  }
  
  _pointerPartonSplitter->split(tagged,pstep);
  if ( HERWIG_DEBUG_LEVEL == HwDebug::extreme_Hadronization ) {   // Occasional debugging
    printStep(pstep,"After PartonSplitter");
  }
  
  _collecCluPtr.clear();
  
  _pointerClusterFinder->formClusters(ch.currentCollision(),pstep,_collecCluPtr);  
  _pointerClusterFinder->reduceToTwoComponents(_collecCluPtr);  
    if ( HERWIG_DEBUG_LEVEL == HwDebug::extreme_Hadronization ) {   // Occasional debugging
    printStep(pstep,"After ClusterFinder");
  }
  
  _pointerColourReconnector->rearrange(ch,pstep,_collecCluPtr);
  
  _pointerClusterFissioner->fission(_collecCluPtr);
  
  _pointerLightClusterDecayer->decay(_collecCluPtr);
  
  _pointerClusterDecayer->decay(_collecCluPtr);
  
  pstep = ch.newStep();
  
  recordAfterClusterDecays(pstep);    
  if ( HERWIG_DEBUG_LEVEL == HwDebug::extreme_Hadronization ) {   // Occasional debugging
    printStep(pstep,"After recordAfterClusterDecays");
  }
  
  // Debugging
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization ) {    
    debuggingInfo(ch);
  }

  //***LOOKHERE*** In the case of soft underlying event ON, the
  //               beam clusters are as slways contained in the
  //               _collecCluPtr, but have being tagged (by
  //               ClusterFissioner) as not available, and therefore
  //               skipped during the cluster hadronization.
  //               At this point (end of the cluster hadronization)
  //               it should be the responsability of this class
  //               (ClusterHadronizationHandler) to pass these
  //               beam clusters to the class responsible to the
  //               soft underlying event. So, when the latter class
  //               will be implemented, you should write few lines
  //               of code right below here to do it.
  
}


void ClusterHadronizationHandler::printStep(tStepPtr ptrStep, const string & title) {
  generator()->log() << " ^^^ Begin Print Step ^^^ : " << title << endl;
  for ( ParticleSet::const_iterator it = ptrStep->particles().begin();
  	 it != ptrStep->particles().end(); ++it ) {
    generator()->log() << "\t" << (*it)->data().PDGName() 
		       << "\t" << (*it)->number() 
		       << "\t" << ( (*it)->colourNeighbour() ? 
				    (*it)->colourNeighbour()->number() : 0 ) 
		       << "\t" << ( (*it)->antiColourNeighbour() ? 
				    (*it)->antiColourNeighbour()->number() : 0 )
                       << "\t" << (*it)->momentum().m() 
		       << "\t" << (*it)->momentum() 
		       << endl;     
  }
  generator()->log() << " ^^^ End Print Step ^^^ " << endl;
}


void ClusterHadronizationHandler::recordAfterClusterDecays(tStepPtr ptrStep) {

  // Loop over all clusters, but skip the ones that are not final 
  // (the ones, by definition, that do not decay into hadrons). 
  // Then, for each final cluster, it goes up the hierarchy through 
  // the parent cluster pointer, until it is reached the initial
  // cluster, that was formed during the cluster finder phase 
  // (and whose constituents are present in the previous Event Record).
  // Notice that this method is general enough also for the case
  // of cluster decay into three (or even more) hadrons. 
 
  for (CollecCluPtr::const_iterator it = _collecCluPtr.begin();
  	 it != _collecCluPtr.end(); ++it) {
  
    if ( (*it)->isAvailable() && (*it)->isStatusFinal() ) {
  
  	// The parents of the hadrons are all of the original
  	// quarks-antiquarks, present after the non-perturbative splitting,
  	// that are associated, directly or indirectly, with the hadrons.
  	// The algorithm is, de facto, recursive, but for efficiency reason
  	// we avoid the recursion call using a dynamic vector as stack.
  
  	vector<tPPtr> theParents;  // particles present in the last Event Record
  	vector<tCluPtr> theStack;  // to avoid recursion
  	theStack.push_back( *it );
  	do {
  	  tCluPtr ptrClu = theStack.back();
  	  theStack.pop_back();        
  	  if ( ptrClu->isStatusInitial() ) {
  	    for ( CollecCompPtr::const_iterator jt = ptrClu->components().begin();
  		  jt != ptrClu->components().end(); ++jt) {
  	      theParents.push_back( (*jt)->pointerParticle() );
  	    }
  	  } else {
  	    theStack.push_back( ptrClu->parentCluster() );
  	  }
  	} while ( ! theStack.empty() );
  	ptrStep->addDecayProduct( theParents.begin(), theParents.end(),
  				  (*it)->childrenHadrons().begin(), 
				  (*it)->childrenHadrons().end() );  
    }
  }
}


void ClusterHadronizationHandler::debuggingInfo(PartialCollisionHandler & ch) {

  // Define static variables to store statistics information to be 
  // printed out at the end of the last event.
  static double statNumberBeamClusters      = 0;
  static double statNumberInitialClusters   = 0;
  static double statNumberFinalClusters     = 0;
  static double statNumber1HadronClusters   = 0;
  static double statNumberReshuffledClusters= 0;
  static double statNumber3ComponentClusters= 0;
  static double statEventHadronMultiplicity = 0;
  static double statNumberMesons            = 0;  
  static double	statNumberLightHadrons      = 0;  
  static double	statNumberStrangeHadrons    = 0;  
  static double	statNumberCharmHadrons      = 0;  
  static double	statNumberBHadrons          = 0;  
  static double	statNumberPlusPions         = 0;  
  static double	statNumberNeutralPions      = 0;  
  static double	statNumberMinusPions        = 0;  
  static map<long,double> statHadrons;

  // The variables below are related to the positions of components and
  // clusters with respect to the collision vertex or to each other.
  // Notice that the invariant space-time distance is positive/negative for
  // time-like/space-like distance, respectively. Also the pure space
  // distance in the Lab frame is reported (it is never negative).
  // In the case of e+ - e- interactions (i.e. no clusters from the initial
  // state radiation) the invariant distance between a cluster and its
  // children clusters is time-like, because information could flow
  // from the parent to the children. The distance between a final cluster
  // and its hadron children is, in average, of type space-like, because
  // as hadron positions we define the production vertex of the hadron,
  // assuming being gaussian smeared, with width inversely proportional
  // to the parent cluster mass, around the parent cluster position:
  // therefore, being the smearing the same in the four space-time components,
  // in average the displacement is of type space, because there are
  // three space components against the only time one. Other invariant
  // distances, like between a light cluster and its reshuffling partner,
  // could be a priori equally time-like as space-like.
  static Length statInvDistInitialComponent_BeamVtx          = Length();
  static Length statSpaceDistInitialComponent_BeamVtx        = Length();
  static Length statInvDistInitialComponent_Cluster          = Length();
  static Length statSpaceDistInitialComponent_Cluster        = Length();
  static double statNumberInitialComponents                  = 0;  
  static Length statInvDistFinalHadron_BeamVtx               = Length();
  static Length statSpaceDistFinalHadron_BeamVtx             = Length();
  static double statNumberChildHadrons                       = 0;
  static Length statInvDistInitialCluster_BeamVtx            = Length();
  static Length statSpaceDistInitialCluster_BeamVtx          = Length();
  // use statNumberInitialClusters for counting the entries.
  static Length statInvDistFinalCluster_BeamVtx              = Length();
  static Length statSpaceDistFinalCluster_BeamVtx            = Length();
  // use statNumberFinalClusters for counting the entries.
  static Length statInvDistAnyCluster_BeamVtx                = Length();
  static Length statSpaceDistAnyCluster_BeamVtx              = Length();
  static double statNumberClusters                           = 0;
  static Length statInvDistParentCluster_ChildCluster        = Length();
  static Length statSpaceDistParentCluster_ChildCluster      = Length();
  static double statNumberChildClusters                      = 0;
  static Length statInvDistParentCluster_ChildHadron         = Length();
  static Length statSpaceDistParentCluster_ChildHadron       = Length();
  // use count statNumberChildHadrons for counting the entries.
  static Length statInvDistLightCluster_ReshufflingCluster   = Length();
  static Length statSpaceDistLightCluster_ReshufflingCluster = Length();
  static double statNumberReshufflingClusterPartners         = 0;

  // Print some information about masses only once. 
  if ( generator()->currentEventNumber() == 1 ) {
    generator()->log() << "ClusterHadronizationHandler::debuggingInfo "
                       << " %%% BEGIN  PRELIMINARY INFOS %%% " << endl;
    // First print info for g,d,u,s,c,b
    for (int i=0; i<= 5; ++i) {
      int id = i;
      if (i == 0) id = ParticleID::g;
      generator()->log() << "\t id=" << id << "   " 
		         << getParticleData(id)->PDGName() << "  mass: "
                         << " nominal=" << getParticleData(id)->mass()
                         << " constituent=" << getParticleData(id)->constituentMass()   
                         << " iSpin=" << getParticleData(id)->iSpin()
		         << endl;
    }
    // Then print info for diquarks.
    for (int i=1; i<= 5; ++i) {
      for (int j=i; j<= 5; ++j) {
	int id = j*1000 + i*100 + 1;
	if (i == j) id = id + 2;     // if equal quarks the less significant digit is 3 not 1
	generator()->log() << "\t id=" << id << "   " 
			   << getParticleData(id)->PDGName() << "  mass: "
			   << " nominal=" << getParticleData(id)->mass()
			   << " constituent=" << getParticleData(id)->constituentMass()
                           << " iSpin=" << getParticleData(id)->iSpin();
	Energy diff = fabs( getParticleData(id)->constituentMass() - 
			    (  getParticleData(i)->constituentMass() +
			       getParticleData(j)->constituentMass() ) );
        if ( diff > 0.001*GeV ) {
	  generator()->log() << " <--- NOT EQUAL SUM CONSTITUENT QUARK MASSES! "; 
	} 
	generator()->log() << endl; 
      }
    }

    // Print about spin: just four examples, spin 0, 1/2, 1, 3/2
    // The method  ParticleData::spin()  returns the spin in standard units, 
    // whereas  ParticleData::iSpin()  returns 2J+1 in units of hbar/2.
    generator()->log() << " eta:    iSpin=" << getParticleData(ParticleID::eta)->iSpin()
                       << "   spin=" << getParticleData(ParticleID::eta)->spin() << endl
                       << " mu:     iSpin=" << getParticleData(ParticleID::muminus)->iSpin()
                       << "   spin=" << getParticleData(ParticleID::muminus)->spin() << endl
                       << " phi:    iSpin=" << getParticleData(ParticleID::phi)->iSpin()
                       << "   spin=" << getParticleData(ParticleID::phi)->spin() << endl
                       << " Omega-: iSpin=" << getParticleData(ParticleID::Omegaminus)->iSpin()
                       << "   spin=" << getParticleData(ParticleID::Omegaminus)->spin() << endl;
    
    // Print about mixing, but only if the debugging level is 100 (see CheckId)
    CheckId::hasStrangeness( ParticleID::eta, 0.5 );
    CheckId::hasStrangeness( ParticleID::etaprime, 0.5 );
    CheckId::hasStrangeness( ParticleID::phi, 0.5 );
    CheckId::hasStrangeness( ParticleID::omega, 0.5 );
    CheckId::hasStrangeness( ParticleID::hprime_1, 0.5 );
    CheckId::hasStrangeness( ParticleID::h_1, 0.5 );
    CheckId::hasStrangeness( ParticleID::f_0, 0.5 );
    CheckId::hasStrangeness( ParticleID::fprime_1, 0.5 );
    CheckId::hasStrangeness( ParticleID::f_1, 0.5 );
    CheckId::hasStrangeness( ParticleID::fprime_2, 0.5 );
    CheckId::hasStrangeness( ParticleID::f_2, 0.5 );

    generator()->log() << "ClusterHadronizationHandler::debuggingInfo "
                       << " %%% END  PRELIMINARY INFOS %%% " << endl;
  }

  // Loop over all clusters, and print information that allows to
  // check the logical consistency of such clusters. Before doing that,
  // associate to each cluster an integer number which corresponds to 
  // the order in which such cluster appear in the container _collecCluPtr,
  // starting with 1 for the first. This is necessary in order to
  // translate in a readable way the cross references between clusters. 
 
  int count = 0;
  map<tCluPtr,int> orderingMap;
  for (CollecCluPtr::const_iterator it = _collecCluPtr.begin();
  	 it != _collecCluPtr.end(); ++it) {
    orderingMap.insert( orderingMap.end(), pair<tCluPtr,int>(*it,++count) );
  }

  LorentzPoint vertexPosition = ch.currentCollision()->vertex();
  generator()->log() << "ClusterHadronizationHandler::debuggingInfo "
                     << " ===> START DEBUGGING <=== "
		     << "   EventNumber=" << generator()->currentEventNumber() << endl
                     << "  Total number of clusters in _collecCluPtr : " << count << endl
		     << "  Lab vertex of current collision : " << vertexPosition << "  [mm]" 
		     << endl;
                     
  // Now detailed information cluster by cluster. Basically everything
  // that can be checked is print out.

  for ( map<tCluPtr,int>::const_iterator it = orderingMap.begin();
        it != orderingMap.end(); ++it ) {

    tCluPtr ptrClu = it->first;  
    int i = it->second;      
    statInvDistAnyCluster_BeamVtx   += ( ptrClu->position() - vertexPosition ).mag();
    statSpaceDistAnyCluster_BeamVtx += ( ptrClu->position().vect() - vertexPosition.vect() ).mag();
    statNumberClusters++;

    generator()->log() << "  --- Cluster --- " << i << endl;

    // Information about the components, including sum of charges and momenta.
    // Notice that the charge of the cluster is defined (hence no consistency
    // check is possible) as the sum of the charges of its components.
    generator()->log() << "\t numComponents     = " << ptrClu->numComponents() << endl;
    int j=0;
    Lorentz5Momentum sumMomentumComponents = Lorentz5Momentum();
    Charge sumChargeComponents = Charge();        
    for ( CollecCompPtr::const_iterator jt = ptrClu->components().begin();
	  jt != ptrClu->components().end() ; ++jt ) {
      generator()->log() << "\t Component " << ++j << endl 
	                 << "\t \t id = " << (*jt)->id() << "    "
			 << getParticleData( (*jt)->id() )->PDGName() 
			 << "     component mass = " << (*jt)->mass() << endl
	                 << "\t \t momentum        = " << (*jt)->momentum() << endl
	                 << "\t \t position        = " << (*jt)->position() << endl
	                 << "\t \t isPerturbative  = " << (*jt)->isPerturbative() << endl
	                 << "\t \t isBeamRemnant   = " << (*jt)->isBeamRemnant() << endl;
      if ( (*jt)->pointerParticle() ) {
	generator()->log() << "\t \t HAS pointed Particle :   number = " 
			   << (*jt)->pointerParticle()->number() << endl
	                   << "\t \t \t positions: labVertex="
	                   << (*jt)->pointerParticle()->labVertex() 
	                   << "   lifeLength=" 
			   << (*jt)->pointerParticle()->lifeLength() << endl
	                   << "\t \t \t masses:   constituent="  
			   << (*jt)->pointerParticle()->data().constituentMass()
			   << "   current=" << (*jt)->pointerParticle()->mass()
			   << "   invariant=" << (*jt)->pointerParticle()->momentum().m() 
			   << endl;
	if ( ptrClu->isStatusInitial() ) {
	  statInvDistInitialComponent_BeamVtx += ( (*jt)->position() - vertexPosition ).mag();
	  statInvDistInitialComponent_Cluster += ( (*jt)->position() - ptrClu->position() ).mag();
	  statSpaceDistInitialComponent_BeamVtx += 
	    ( (*jt)->position().vect() - vertexPosition.vect() ).mag();
	  statSpaceDistInitialComponent_Cluster += 
	    ( (*jt)->position().vect() - ptrClu->position().vect() ).mag();
	  statNumberInitialComponents++;
	}
      } else {
	generator()->log() << "\t \t NO pointed Particle " << endl;  
      }
      sumMomentumComponents += (*jt)->momentum();
      sumChargeComponents += getParticleData( (*jt)->id() )->charge();      
    }
    generator()->log() << "\t charge (sum components) = " << sumChargeComponents << endl
		       << "\t sumConstituentMasses    = " << ptrClu->sumConstituentMasses() 
		       << endl
		       << "\t mass              = " << ptrClu->mass() << endl
		       << "\t momentum          = " << ptrClu->momentum() << endl
		       << "\t position          = " << ptrClu->position() << endl;    
    Lorentz5Momentum diff = sumMomentumComponents - ptrClu->momentum();
    Energy ediff = fabs( diff.m() );
    if ( ediff > 1e-3*GeV ) {
      generator()->log() << "\t ***ERROR***: MOMENTUM NOT CONSERVED! " << endl;
      generator()->log() << "\t  --> sumMomentumComponents = " << sumMomentumComponents << endl 
			 << "\t      diff = " << diff << endl;
    } 
    generator()->log() << "\t isBeamCluster     = " << ptrClu->isBeamCluster() << endl;
    if ( ptrClu->isBeamCluster() ) statNumberBeamClusters++;

    if ( ! ptrClu->isAvailable() ) {
      generator()->log() << "\t isAvailable     = " << ptrClu->isAvailable() << endl;
    } else {

      generator()->log() << "\t isStatusInitial   = " << ptrClu->isStatusInitial() << endl
			 << "\t isReadyToDecay    = " << ptrClu->isReadyToDecay() << endl
			 << "\t isRedefined       = " << ptrClu->isRedefined() << endl
			 << "\t hasBeenReshuffled = " << ptrClu->hasBeenReshuffled() << endl
			 << "\t isStatusFinal     = " << ptrClu->isStatusFinal() << endl;

      // Statistics information
      if ( ptrClu->isStatusInitial() ) {
	statInvDistInitialCluster_BeamVtx += ( ptrClu->position() - vertexPosition ).mag();	
	statSpaceDistInitialCluster_BeamVtx += 
	  ( ptrClu->position().vect() - vertexPosition.vect() ).mag();	
	statNumberInitialClusters++;
      }
      if ( ptrClu->isStatusFinal() ) {
	statInvDistFinalCluster_BeamVtx += ( ptrClu->position() - vertexPosition ).mag();	
	statSpaceDistFinalCluster_BeamVtx +=
	  ( ptrClu->position().vect() - vertexPosition.vect() ).mag();	
	statNumberFinalClusters++;
      }
      if ( ptrClu->childrenHadrons().size() == 1 ) statNumber1HadronClusters++;
      if ( ptrClu->hasBeenReshuffled() ) {
	statNumberReshuffledClusters++;
	if ( ptrClu->isRedefined() ) {
	  statInvDistLightCluster_ReshufflingCluster += 
	    ( ptrClu->position() - ptrClu->reshufflingPartnerCluster()->position() ).mag();
	  statSpaceDistLightCluster_ReshufflingCluster += 
	    ( ptrClu->position().vect() - ptrClu->reshufflingPartnerCluster()->position().vect() ).mag();
	  statNumberReshufflingClusterPartners++;
	}
      } else if ( ptrClu->isRedefined() ) {
	statNumber3ComponentClusters++;
      }

      // Information about the (eventual) parent cluster. 
      int parent = 0;
      if ( ptrClu->parentCluster() ) {
	if ( orderingMap.find( ptrClu->parentCluster() ) != orderingMap.end() ) {
	  parent = orderingMap.find( ptrClu->parentCluster() )->second;
	} else {
	  parent = -999;  // Error: it shouldn't happen!
	}
      }
      generator()->log() << "\t parentCluster     = " << parent << endl;
      
      // Information about the (eventual) reshuffling partner cluster. 
      int reshuffling = 0;
      if ( ptrClu->reshufflingPartnerCluster() ) {
	if ( orderingMap.find( ptrClu->reshufflingPartnerCluster() ) != orderingMap.end() ) {
	  reshuffling = orderingMap.find( ptrClu->reshufflingPartnerCluster() )->second;
	} else {
	  reshuffling = -999;  // Error: it shouldn't happen!
	}
      }
      generator()->log() << "\t reshufflingPartner= " << reshuffling << endl;
      
      // Information about the (eventual) children clusters, including the check
      // of conservation of the charge and energy-momentum.
      Lorentz5Momentum sumMomentumChildrenClusters = Lorentz5Momentum();
      Charge sumChargeChildrenClusters = Charge();
      if ( ptrClu->childrenClusters() != CollecCluPtr() ) {
	generator()->log() << "\t childrenClusters  = ";
	for ( CollecCluPtr::const_iterator jt = ptrClu->childrenClusters().begin();
	      jt != ptrClu->childrenClusters().end(); ++jt ) {	
	  statInvDistParentCluster_ChildCluster += ( ptrClu->position() - (*jt)->position() ).mag();
	  statSpaceDistParentCluster_ChildCluster += 
	    ( ptrClu->position().vect() - (*jt)->position().vect() ).mag();
	  statNumberChildClusters++;
	  int child = 0;
	  if ( orderingMap.find( *jt ) != orderingMap.end() ) {
	    child = orderingMap.find( *jt )->second;
	  } else {
	    child = -999;  // Error: it shouldn't happen!
	  }
	  generator()->log() << child << "   ";
	  Charge sumChargeChildrenClusterComponents = Charge();
	  // Remind that the charge of a cluster is defined as the sum of the 
	  // charges of its components. 
	  for ( CollecCompPtr::const_iterator kt = (*jt)->components().begin();
		kt != (*jt)->components().end() ; ++kt ) {
	    sumChargeChildrenClusterComponents += getParticleData( (*kt)->id() )->charge();
	  }
	  sumChargeChildrenClusters   += sumChargeChildrenClusterComponents;
	  sumMomentumChildrenClusters += (*jt)->momentum();
	}
	generator()->log() << endl
			   << "\t  --> sumChargeChildrenClusters = " 
			   << sumChargeChildrenClusters << endl
			   << "\t  --> sumMomentumChildrenClusters = " 
			   << sumMomentumChildrenClusters << endl;
      } else {      
	generator()->log() << "\t NO childrenClusters " << endl;
      }
      
      // Information about the (eventual) children hadrons, including the check
      // of conservation of the charge and energy-momentum.
      Lorentz5Momentum sumMomentumChildrenHadrons = Lorentz5Momentum();
      Charge sumChargeChildrenHadrons = Charge();
      if ( ptrClu->childrenHadrons() != PVector() ) {
	generator()->log() << "\t num childrenHadrons  = " 
			   << ptrClu->childrenHadrons().size() << endl;
	int j = 0;
	for ( PVector::const_iterator jt = ptrClu->childrenHadrons().begin();
	      jt != ptrClu->childrenHadrons().end(); ++jt ) {	
	  statInvDistFinalHadron_BeamVtx += ( (*jt)->labVertex() - vertexPosition ).mag();
	  statSpaceDistFinalHadron_BeamVtx += 
	    ( (*jt)->labVertex().vect() - vertexPosition.vect() ).mag();
	  statInvDistParentCluster_ChildHadron += ( ptrClu->position() - (*jt)->labVertex() ).mag();
	  statSpaceDistParentCluster_ChildHadron += 
	    ( ptrClu->position().vect() - (*jt)->labVertex().vect() ).mag();
	  statNumberChildHadrons++;
	  generator()->log() << "\t Hadron " << ++j << endl 
			     << "\t \t id = " << (*jt)->id() << "   " << (*jt)->PDGName() 
			     << "   mass = " << (*jt)->mass() << endl
			     << "\t \t momentum        = " << (*jt)->momentum() << endl
			     << "\t \t labVertex       = " << (*jt)->labVertex() << endl;
	  // Information about the angles, in Cluster CM and in Lab frame, 
	  // between the hadron and the component of the cluster parent from 
	  // which the hadron come from. Although not nice, we use the fact the 
	  // order of cluster's components is the same as cluster's hadrons.
	  int k = 0;
	  Lorentz5Momentum pQLab = Lorentz5Momentum();
	  for ( CollecCompPtr::const_iterator kt = ptrClu->components().begin();
		kt != ptrClu->components().end() ; ++kt ) {
	    if ( ++k == j ) pQLab = (*kt)->momentum();
	  }
	  Lorentz5Momentum pQCm = pQLab;
	  pQCm.boost( - ptrClu->momentum().boostVector() ); 
	  Lorentz5Momentum pHadCm = (*jt)->momentum();
	  pHadCm.boost( - ptrClu->momentum().boostVector() ); 
	  generator()->log() << "\t \t angles:  CM = "
			     << pHadCm.vect().angle( pQCm.vect().unit() ) << "   Lab = " 
			     << (*jt)->momentum().vect().angle( pQLab.vect().unit() ) 
			     << endl;
	  if ( (*jt)->parents() != tParticleVector() ) {
	    generator()->log() << "\t \t num parents = " << (*jt)->parents().size() << endl;
	    int k = 0;
	    for ( tParticleVector::const_iterator kt = (*jt)->parents().begin();
		  kt != (*jt)->parents().end(); ++kt ) {	
	      generator()->log() << "\t \t Parent " << ++k << endl
				 << "\t \t \t id = " << (*kt)->id() << "   " << (*kt)->PDGName() 
				 << endl
				 << "\t \t \t momentum  = " << (*kt)->momentum() << endl
				 << "\t \t \t labVertex = " << (*kt)->labVertex() << endl;
	    }
	  } else {      
	    generator()->log() << "\t \t NO parents for this hadron" << endl;
	  }	
	  sumChargeChildrenHadrons   += (*jt)->data().charge();
	  sumMomentumChildrenHadrons += (*jt)->momentum();
	}
	generator()->log() << endl 
			   << "\t  --> sumChargeChildrenHadrons = " 
			   << sumChargeChildrenHadrons << endl	  
			   << "\t  --> sumMomentumChildrenHadrons = " 
			   << sumMomentumChildrenHadrons << endl;   
      } else {      
	generator()->log() << "\t NO childrenHadrons " << endl;
      }
    
      if ( ptrClu->childrenClusters() != CollecCluPtr()  || 
	   ptrClu->childrenHadrons()  != PVector() ) {
	Charge delta_ch = fabs( getParticleData(ParticleID::d)->charge() ) / 10.0; 
	if ( fabs(sumChargeChildrenClusters + sumChargeChildrenHadrons 
		  - sumChargeComponents) > delta_ch ) {
	  generator()->log() << "\t ***ERROR***: CHARGE NOT CONSERVED! " << endl;	
	} 
	diff = sumMomentumChildrenClusters + sumMomentumChildrenHadrons - ptrClu->momentum();
	ediff = fabs( diff.m() );
	if ( ediff > 1e-3*GeV ) {
	  if ( ptrClu->hasBeenReshuffled() ) {
	    generator()->log() << "\t ***WARNING***: MOMENTUM NOT CONSERVED"
			       << " BUT RESHUFFLING OCCURED!" << endl; 
	  } else {
	    generator()->log() << "\t ***ERROR***: MOMENTUM NOT CONSERVED! " << endl;
	  }
	  generator()->log() << "\t  --> sumMomentumChildren = " 
			     << sumMomentumChildrenClusters + sumMomentumChildrenHadrons << endl 
			     << "\t      diff = " << diff << endl;
	} 
      }

    } // end else part of if ( ! ptrClu->isAvailable() )    
      
  } // end main loop over clusters (indeed, over the elements of the map)

  // We want to focus directly now on the connection between the 
  // original partons (before the cluster hadronization but after the 
  //  nonperturbative gluon splitting) and the final hadrons.
  // To do so, we start as before by looping over all clusters, but
  // then we select only the initial clusters (the ones obtained by
  //  the ClusterFinder from the partons), from which we extract the
  // constituent partons, and then we consider all of the children
  // hadrons belonging to all final clusters (by definition the ones
  //  that decay into hadrons) which derive directly or indirectly 
  // (through a cluster product of a fission of an heavy cluster) 
  // from the original initial cluster. Finally, we check the charge 
  // and momentum conservation between the original parent partons 
  // and the final hadrons. Furthermore, we update all the static
  // variables used for storing statistical information. 
  // Notice that we do *not* use at all the Event Record, therefore
  // these checks could be useful not only to debug the mechanics of
  // Cluster Hadronization algorithm, but also bugs or improper use
  // of the Event Record interface (for example in the use of
  //      ptrStep->addDecayProduct(...)
  // in the method recordAfterClusterDecays).

  generator()->log() << " ---------- Begin analysis of produced Hadrons ---------- " << endl;
  for ( map<tCluPtr,int>::const_iterator it = orderingMap.begin();
        it != orderingMap.end(); ++it ) {
    tCluPtr ptrClu = it->first;  
    if ( ptrClu->isAvailable() && ptrClu->isStatusInitial() ) {
      generator()->log() << "  Cluster (Initial) " << it->second << endl;
      bool reshuffled = false;
      vector<tPPtr> vecPartonParents;
      for ( CollecCompPtr::const_iterator jt = ptrClu->components().begin();
	    jt != ptrClu->components().end(); ++jt) {
	vecPartonParents.push_back( (*jt)->pointerParticle() );
      }
      vector<tPPtr> vecChildrenHadrons;
      vector<tCluPtr> theStack;  // to avoid recursion
      theStack.push_back( ptrClu );
      do {
	tCluPtr ptrClu = theStack.back();
	theStack.pop_back();        
	if ( ptrClu->isStatusFinal() ) {
	  for ( PVector::const_iterator jt = ptrClu->childrenHadrons().begin();
		jt != ptrClu->childrenHadrons().end(); ++jt) {
	    vecChildrenHadrons.push_back( *jt );
	  }
	} 
	for ( CollecCluPtr::const_iterator jt = ptrClu->childrenClusters().begin();
	      jt != ptrClu->childrenClusters().end(); ++jt) {
	  theStack.push_back( *jt );
	}
	if ( ptrClu->hasBeenReshuffled() ) {
	  reshuffled = true;
	}
      } while ( ! theStack.empty() );    
      Lorentz5Momentum sumMomentumPartons = Lorentz5Momentum();
      Charge sumChargePartons = Charge();
      generator()->log() << "\t made of the following  " << vecPartonParents.size()
	                 << "  partons " << endl;
      for ( vector<tPPtr>::const_iterator jt = vecPartonParents.begin();
	    jt != vecPartonParents.end(); ++jt) {
	generator()->log() << "\t \t id = " << (*jt)->id() << "   " << (*jt)->PDGName() 
			   << "   number = " << (*jt)->number() << "   " 
			   << (*jt)->momentum() << endl;
	sumMomentumPartons += (*jt)->momentum();
        sumChargePartons += (*jt)->data().charge();
      }
      Lorentz5Momentum sumMomentumHadrons = Lorentz5Momentum();
      Charge sumChargeHadrons = Charge();
      generator()->log() << "\t has produced the following  " << vecChildrenHadrons.size()
	                 << "  hadrons " << endl;
      for ( vector<tPPtr>::const_iterator jt = vecChildrenHadrons.begin();
	    jt != vecChildrenHadrons.end(); ++jt) {
	generator()->log() << "\t \t id = " << (*jt)->id() << "   " << (*jt)->PDGName() 
			   << "   " << (*jt)->momentum() << endl;
	sumMomentumHadrons += (*jt)->momentum();
        sumChargeHadrons += (*jt)->data().charge();
        // Add here all other statistical information
        statEventHadronMultiplicity++;
        if ( statHadrons.find( (*jt)->id() ) != statHadrons.end() ) {
	  statHadrons.find( (*jt)->id() )->second = statHadrons.find( (*jt)->id() )->second + 1;
	} else {
	  statHadrons.insert( pair<long,double>( (*jt)->id(), 1.0 ) );
	}
        if ( CheckId::isMeson( (*jt)->id() ) ) statNumberMesons++;  
        if ( CheckId::hasBeauty( (*jt)->id() ) ) {
	  statNumberBHadrons++;
	} else if ( CheckId::hasCharm( (*jt)->id() ) ) {
	  statNumberCharmHadrons++;
	} else if ( CheckId::hasStrangeness( (*jt)->id(), rnd() ) ) {
	  statNumberStrangeHadrons++;
	} else {
	  statNumberLightHadrons++;
	}
        if ( (*jt)->id() == ParticleID::piplus )  statNumberPlusPions++;  
        if ( (*jt)->id() == ParticleID::pi0 )     statNumberNeutralPions++;  
        if ( (*jt)->id() == ParticleID::piminus ) statNumberMinusPions++;          
      }
      Charge delta_ch = fabs( getParticleData(ParticleID::d)->charge() ) / 10.0; 
      if ( fabs(sumChargePartons - sumChargeHadrons) > delta_ch ) {
	generator()->log() << "\t ***ERROR***: CHARGE NOT CONSERVED! " << endl;
      } 
      generator()->log() << "\t Total Charge Partons = " << sumChargePartons << endl
			 << "\t Total Charge Hadrons = " << sumChargeHadrons << endl;
      Lorentz5Momentum diff = sumMomentumPartons - sumMomentumHadrons;
      Energy ediff = fabs( diff.m() );
      if ( ediff > 1e-3*GeV ) {
	generator()->log() << "\t ***ERROR***: MOMENTUM NOT CONSERVED!"
			   << "  reshuffled=" << reshuffled << endl;
      }
      generator()->log() << "\t Total Momentum Partons = " << sumMomentumPartons << endl
			 << "\t Total Momentum Hadrons = " << sumMomentumHadrons << endl
			 << "\t   ---> Diff Momentum = " << diff << endl;
    } // end of if isAvailable && isStatusInitial
  } // endl loop over the map
  generator()->log() << " ---------- End analysis of produced Hadrons ---------- " << endl;

  generator()->log() << "ClusterHadronizationHandler::debuggingInfo "
                     << " ===> END DEBUGGING <=== " << endl;

  // At the end of the last event to be generated, print out statistics information.
  if ( generator()->currentEventNumber() == generator()->N() ) {
    generator()->log() << endl << endl << "ClusterHadronizationHandler::debuggingInfo "
		       << " ===> BEGIN ***FINAL STATISTICS*** <=== " << endl << endl;
    double dN = double ( generator()->N() );
    int i=0;
    for ( map<long,double>::const_iterator iter = statHadrons.begin();
	  iter != statHadrons.end(); ++iter ) {
      generator()->log() << "\t" << ++i << ")  id=" << iter->first << "   "
	                 << ( getParticleData( iter->first ) ? getParticleData( iter->first )->PDGName()
			      : "UNKNOWN" ) 
			 << "   " << iter->second / dN << " per event;   fraction="
			 <<  100.0 * iter->second / statEventHadronMultiplicity << " %" << endl;
    }
    generator()->log() << endl
		       << "\t Counts for position statistics " << endl
                       << "\t \t   statNumberInitialComponents =" << statNumberInitialComponents << endl
                       << "\t \t   statNumberInitialClusters   =" << statNumberInitialClusters << endl
                       << "\t \t   statNumberClusters          =" << statNumberClusters << endl
                       << "\t \t   statNumberFinalClusters     =" << statNumberFinalClusters << endl
                       << "\t \t   statNumberChildClusters     =" << statNumberChildClusters << endl
                       << "\t \t   statNumberChildHadrons      =" << statNumberChildHadrons << endl
                       << "\t \t   statNumberReshufflingClusterParters =" 
		       << statNumberReshufflingClusterPartners << endl 
		       << endl
                       << "\t Invariant (space-time) Distance (in millimeters) : " << endl
                       << "\t \t | initial component - beam vertex | = " 
                       << ( statNumberInitialComponents 
			    ? statInvDistInitialComponent_BeamVtx / statNumberInitialComponents  
			    : 0.0 ) << endl
                       << "\t \t | initial component - cluster | = " 
                       << ( statNumberInitialComponents 
			    ? statInvDistInitialComponent_Cluster / statNumberInitialComponents  
			    : 0.0 ) << endl
                       << "\t \t | initial cluster - beam vertex | = " 
                       << ( statNumberInitialClusters 
			    ? statInvDistInitialCluster_BeamVtx / statNumberInitialClusters  
			    : 0.0 ) << endl
                       << "\t \t | any cluster - beam vertex | = " 
                       << ( statNumberClusters 
			    ? statInvDistAnyCluster_BeamVtx / statNumberClusters 
			    : 0.0 ) << endl
                       << "\t \t | final cluster - beam vertex | = " 
                       << ( statNumberFinalClusters 
			    ? statInvDistFinalCluster_BeamVtx / statNumberFinalClusters  
			    : 0.0 ) << endl
                       << "\t \t | final hadron - beam vertex | = " 
                       << ( statNumberChildHadrons 
			    ? statInvDistFinalHadron_BeamVtx / statNumberChildHadrons  
			    : 0.0 ) << endl
                       << "\t \t | parent cluster - child cluster | = " 
                       << ( statNumberChildClusters 
			    ? statInvDistParentCluster_ChildCluster / statNumberChildClusters  
			    : 0.0 ) << endl
                       << "\t \t | parent cluster - child hadron | = " 
                       << ( statNumberChildHadrons 
			    ? statInvDistParentCluster_ChildHadron / statNumberChildHadrons  
			    : 0.0 ) << endl
                       << "\t \t | light cluster - reshuffling partner cluster | = " 
                       << ( statNumberReshufflingClusterPartners
			    ? statInvDistLightCluster_ReshufflingCluster 
			    / statNumberReshufflingClusterPartners : 0.0 ) << endl
		       << endl
                       << "\t Lab space distance (abs) (in millimeters) : " << endl
                       << "\t \t | initial component - beam vertex | = " 
                       << ( statNumberInitialComponents 
			    ? statSpaceDistInitialComponent_BeamVtx / statNumberInitialComponents  
			    : 0.0 ) << endl
                       << "\t \t | initial component - cluster | = " 
                       << ( statNumberInitialComponents 
			    ? statSpaceDistInitialComponent_Cluster / statNumberInitialComponents  
			    : 0.0 ) << endl
                       << "\t \t | initial cluster - beam vertex | = " 
                       << ( statNumberInitialClusters 
			    ? statSpaceDistInitialCluster_BeamVtx / statNumberInitialClusters  
			    : 0.0 ) << endl
                       << "\t \t | any cluster - beam vertex | = " 
                       << ( statNumberClusters 
			    ? statSpaceDistAnyCluster_BeamVtx / statNumberClusters 
			    : 0.0 ) << endl
                       << "\t \t | final cluster - beam vertex | = " 
                       << ( statNumberFinalClusters 
			    ? statSpaceDistFinalCluster_BeamVtx / statNumberFinalClusters  
			    : 0.0 ) << endl
                       << "\t \t | final hadron - beam vertex | = " 
                       << ( statNumberChildHadrons 
			    ? statSpaceDistFinalHadron_BeamVtx / statNumberChildHadrons  
			    : 0.0 ) << endl
                       << "\t \t | parent cluster - child cluster | = " 
                       << ( statNumberChildClusters 
			    ? statSpaceDistParentCluster_ChildCluster / statNumberChildClusters  
			    : 0.0 ) << endl
                       << "\t \t | parent cluster - child hadron | = " 
                       << ( statNumberChildHadrons 
			    ? statSpaceDistParentCluster_ChildHadron / statNumberChildHadrons  
			    : 0.0 ) << endl
                       << "\t \t | light cluster - reshuffling partner cluster | = " 
                       << ( statNumberReshufflingClusterPartners
			    ? statSpaceDistLightCluster_ReshufflingCluster 
			    / statNumberReshufflingClusterPartners : 0.0 ) << endl
                       << endl;
    generator()->log() << endl
                       << "\t Clusters:  BeamCluster  average number  = " 
                       << statNumberBeamClusters / dN << endl
                       << "\t            Initial      average number  = " 
                       << statNumberInitialClusters / dN << endl
                       << "\t            Final        average number  = " 
                       << statNumberFinalClusters / dN << endl
                       << "\t            1-Hadron     average number  = " 
                       << statNumber1HadronClusters / dN << "  fraction="
		       << statNumber1HadronClusters / statNumberFinalClusters << endl
                       << "\t            Reshuffled   average number  = " 
                       << statNumberReshuffledClusters / dN << "  fraction="
		       << statNumberReshuffledClusters / statNumberFinalClusters << endl
                       << "\t            3-Component  average number  = " 
                       << statNumber3ComponentClusters / dN << "  fraction="
		       << statNumber3ComponentClusters / statNumberFinalClusters <<endl<<endl 
                       << "\t Average Hadron Multiplicity = " 
                       << statEventHadronMultiplicity / dN << endl
                       << "\t Mesons:  " << statNumberMesons / dN << " per event;    fraction="
                       << 100.0 * statNumberMesons / statEventHadronMultiplicity << " %" << endl
                       << "\t Baryons: " << (statEventHadronMultiplicity-statNumberMesons) / dN 
		       << " per event;    fraction="
                       << 100.0 * ( 1.0 - statNumberMesons / statEventHadronMultiplicity ) 
                       << " %" << endl
                       << "\t b-Hadrons    = " 
                       << statNumberBHadrons / dN << " per event;    fraction="
                       << 100.0 * statNumberBHadrons / statEventHadronMultiplicity << " %" << endl
                       << "\t c-Hadrons    = " 
                       << statNumberCharmHadrons / dN << " per event;    fraction="
                       << 100.0 * statNumberCharmHadrons / statEventHadronMultiplicity << " %" << endl
                       << "\t s-Hadrons    = " 
                       << statNumberStrangeHadrons / dN << " per event;    fraction="
                       << 100.0 * statNumberStrangeHadrons / statEventHadronMultiplicity << " %" << endl
                       << "\t u,d-Hadrons  = " 
                       << statNumberLightHadrons / dN << " per event;    fraction="
                       << 100.0 * statNumberLightHadrons / statEventHadronMultiplicity << " %" << endl
                       << "\t \t pi+  = " 
                       << statNumberPlusPions / dN << " per event;    fraction="
                       << 100.0 * statNumberPlusPions / statEventHadronMultiplicity << " %" << endl
                       << "\t \t pi0  = " 
                       << statNumberNeutralPions / dN << " per event;    fraction="
                       << 100.0 * statNumberNeutralPions / statEventHadronMultiplicity << " %" << endl
                       << "\t \t pi-  = "
                       << statNumberMinusPions / dN << " per event;    fraction=" 
                       << 100.0 * statNumberMinusPions / statEventHadronMultiplicity << " %" << endl;
    generator()->log() << endl << "ClusterHadronizationHandler::debuggingInfo "
		       << " ===> END ***FINAL STATISTICS*** <=== " << endl << endl;
  }
  
}

