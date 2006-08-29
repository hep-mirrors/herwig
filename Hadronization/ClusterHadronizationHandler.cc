// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ClusterHadronizationHandler class.
//

#include "ClusterHadronizationHandler.h"
#include <ThePEG/Interface/ClassDocumentation.h>
#include <ThePEG/Persistency/PersistentOStream.h>
#include <ThePEG/Persistency/PersistentIStream.h>
#include <ThePEG/Interface/Parameter.h> 
#include <ThePEG/Interface/Reference.h>
#include <ThePEG/Handlers/EventHandler.h>
#include <ThePEG/Handlers/Hint.h>
#include <ThePEG/PDT/ParticleData.h>
#include <ThePEG/EventRecord/Particle.h>
#include <ThePEG/EventRecord/Step.h>
#include <ThePEG/PDT/PDT.h>
#include <ThePEG/PDT/EnumParticles.h>
#include "Remnant.h"
#include "Herwig++/Utilities/EnumParticles.h"
#include <ThePEG/Repository/EventGenerator.h>
#include "Herwig++/Utilities/HwDebug.h"
#include "Herwig++/Utilities/CheckId.h"
#include "Herwig++/Utilities/Smearing.h"
#include "CluHadConfig.h"
#include "Cluster.h"  
#include "Remnant.h"
#include <iostream>

using namespace Herwig;

void ClusterHadronizationHandler::persistentOutput(PersistentOStream & os) 
  const {
  os << _globalParameters 
     << _partonSplitter 
     << _clusterFinder
     << _colourReconnector
     << _clusterFissioner
     << _lightClusterDecayer
     << _clusterDecayer
     << _forcedSplitter;
}


void ClusterHadronizationHandler::persistentInput(PersistentIStream & is, int) {
  is >> _globalParameters 
     >> _partonSplitter 
     >> _clusterFinder
     >> _colourReconnector
     >> _clusterFissioner
     >> _lightClusterDecayer
     >> _clusterDecayer
     >> _forcedSplitter;
}

ClassDescription<ClusterHadronizationHandler> ClusterHadronizationHandler::initClusterHadronizationHandler;
// Definition of the static class description member.


void ClusterHadronizationHandler::Init() {

  static ClassDocumentation<ClusterHadronizationHandler> documentation
    ("This is the main handler class for the Cluster Hadronization");

  static Reference<ClusterHadronizationHandler,GlobalParameters> 
    interfaceGlobalParameters("GlobalParameters", 
		      "A reference to the GlobalParameters object", 
		      &Herwig::ClusterHadronizationHandler::_globalParameters,
		      false, false, true, false);

  static Reference<ClusterHadronizationHandler,PartonSplitter> 
    interfacePartonSplitter("PartonSplitter", 
		      "A reference to the PartonSplitter object", 
		      &Herwig::ClusterHadronizationHandler::_partonSplitter,
		      false, false, true, false);

  static Reference<ClusterHadronizationHandler,ClusterFinder> 
    interfaceClusterFinder("ClusterFinder", 
		      "A reference to the ClusterFinder object", 
		      &Herwig::ClusterHadronizationHandler::_clusterFinder,
		      false, false, true, false);

  static Reference<ClusterHadronizationHandler,ColourReconnector> 
    interfaceColourReconnector("ColourReconnector", 
		      "A reference to the ColourReconnector object", 
		      &Herwig::ClusterHadronizationHandler::_colourReconnector,
		      false, false, true, false);

  static Reference<ClusterHadronizationHandler,ClusterFissioner> 
    interfaceClusterFissioner("ClusterFissioner", 
		      "A reference to the ClusterFissioner object", 
		      &Herwig::ClusterHadronizationHandler::_clusterFissioner,
		      false, false, true, false);

  static Reference<ClusterHadronizationHandler,LightClusterDecayer> 
    interfaceLightClusterDecayer("LightClusterDecayer", 
		    "A reference to the LightClusterDecayer object", 
		    &Herwig::ClusterHadronizationHandler::_lightClusterDecayer,
		    false, false, true, false);

  static Reference<ClusterHadronizationHandler,ClusterDecayer> 
    interfaceClusterDecayer("ClusterDecayer", 
		       "A reference to the ClusterDecayer object", 
		       &Herwig::ClusterHadronizationHandler::_clusterDecayer,
		       false, false, true, false);

  static Reference<ClusterHadronizationHandler,ForcedSplitting> interfaceForcedSplitting
    ("ForcedSplitting",
     "Object responsible for the forced splitting of the Remnant",
     &ClusterHadronizationHandler::_forcedSplitter, false, false, true, false, false);

}


void ClusterHadronizationHandler::doinitrun() {
  // The run initialization is used here to all Cluster to have access to the
  // GlobalParameters class instance, via a static pointer.
  Cluster::setPointerGlobalParameters(_globalParameters);
}


void ClusterHadronizationHandler::
handle(EventHandler & ch, const tPVector & tagged,
       const Hint & hint) throw(Veto, Stop, Exception) {
  ClusterVector clusters;
  StepPtr pstep = ch.newStep();


  if(HERWIG_DEBUG_LEVEL == HwDebug::extreme_Hadronization) { 
    printStep(pstep,"At the beginning of ClusterHadronizationHandler");
  }
  // split the remnants if needed
  tPVector partonsA=_forcedSplitter->split(tagged,pstep);

  // split the gluons
  tPVector partons=_partonSplitter->split(partonsA,pstep);

  if(HERWIG_DEBUG_LEVEL == HwDebug::extreme_Hadronization) { 
    printStep(pstep,"After PartonSplitter");
  }

  // form the clusters
  pstep = ch.newStep();
  _clusterFinder->formClusters(ch.currentCollision(),pstep,partons,clusters); 
  _clusterFinder->reduceToTwoComponents(pstep,clusters); 
  if(HERWIG_DEBUG_LEVEL == HwDebug::extreme_Hadronization) {
    printStep(pstep,"After ClusterFinder");
  }

  // perform colour reconnection if needed and then
  // decay the clusters into one hadron
  bool lightOK = false;
  short tried = 0;
  while (!lightOK && tried++ < 10) {
    pstep = ch.newStep();
    _colourReconnector->rearrange(ch,pstep,clusters);
    _clusterFissioner->fission(pstep);

    if(HERWIG_DEBUG_LEVEL == HwDebug::extreme_Hadronization) {
      printStep(pstep,"After ClusterFissioner");
    }

    pstep = ch.newStep();
    lightOK = _lightClusterDecayer->decay(pstep);
    if (!lightOK) {
      ch.popStep(); ch.popStep(); 
      if(HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization) {    
	generator()->log() << "CluHad::handle(): throw away " 
			   << tried 
			   << ". attempt of LightClusterDecayer!" 
			   << endl;
      }
      generator()->out() << "CluHad::handle(): throw away " 
			 << tried 
			 << ". attempt of LightClusterDecayer!" 
			 << endl 
			 << "  in evt#" << generator()->currentEventNumber() 
			 << endl;	
    }
  } 
  if (!lightOK)
    throw Exception("CluHad::handle(): tried LightClusterDecayer 10 times!", Exception::eventerror);
  
  if (lightOK && tried > 1 && HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization) 
    printStep(pstep,"After LightClusterDecayer");

  // decay the remaining clusters
  _clusterDecayer->decay(pstep);

  // Debugging
  if(HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization) {    
    for (int j=4; j<6; j++) {
      cStepPtr temp = ch.currentCollision()->step(j); 
      for (ParticleSet::iterator it = temp->all().begin();
	   it!= temp->all().end(); it++) { 
	if((*it)->id() == ExtraParticleID::Cluster) 
	  clusters.push_back(dynamic_ptr_cast<ClusterPtr>(*it));
      }
    }
    debuggingInfo(ch, clusters);
  }

  // ***LOOKHERE*** In the case of soft underlying event ON, the
  //                beam clusters are as slways contained in the
  //                _collecCluPtr, but have being tagged (by
  //                ClusterFissioner) as not available, and therefore
  //                skipped during the cluster hadronization.
  //                At this point (end of the cluster hadronization)
  //                it should be the responsability of this class
  //                (ClusterHadronizationHandler) to pass these
  //                beam clusters to the class responsible to the
  //                soft underlying event. So, when the latter class
  //                will be implemented, you should write few lines
  //                of code right below here to do it.
}


/****************************************************************
 * The remaining code is used for debugging purposes. The code  *
 * itself should be of no real interest.                        *
 ***************************************************************/
void ClusterHadronizationHandler::printStep(tStepPtr ptrStep, const string & title) {
  generator()->log() << "############" << title << "##############" << endl;
  for ( ParticleSet::const_iterator it = ptrStep->particles().begin();
  	 it != ptrStep->particles().end(); ++it ) {
    generator()->log() << *(*it);
  }
  generator()->log() << "###########################################" << endl;
}

void ClusterHadronizationHandler::debuggingInfo(EventHandler & ch,
						ClusterVector &clusters) {

  // Define static variables to store statistics information to be 
  // printed out at the end of the last event.
  static double numberBeamClusters      = 0;
  static double numberInitialClusters   = 0;
  static double numberFinalClusters     = 0;
  static double number1HadronClusters   = 0;
  static double numberReshuffledClusters= 0;
  static double number3ComponentClusters= 0;
  static double eventHadronMultiplicity = 0;
  static double numberMesons            = 0;  
  static double	numberLightHadrons      = 0;  
  static double	numberStrangeHadrons    = 0;  
  static double	numberCharmHadrons      = 0;  
  static double	numberBHadrons          = 0;  
  static double	numberPlusPions         = 0;  
  static double	numberNeutralPions      = 0;  
  static double	numberMinusPions        = 0;  
  static map<long,double> hadrons;

  /****************************************************************************
   * The variables below are related to the positions of components and       *
   * clusters with respect to the collision vertex or to each other.          *
   * Notice that the invariant space-time distance is positive/negative for   *
   * time-like/space-like distance, respectively. Also the pure space         *
   * distance in the Lab frame is reported (it is never negative).            *
   * In the case of e+ - e- interactions (i.e. no clusters from the initial   *
   * state radiation) the invariant distance between a cluster and its        *
   * children clusters is time-like, because information could flow           *
   * from the parent to the children. The distance between a final cluster    *
   * and its hadron children is, in average, of type space-like, because      *
   * as hadron positions we define the production vertex of the hadron,       *
   * assuming being gaussian smeared, with width inversely proportional       *
   * to the parent cluster mass, around the parent cluster position:          *
   * therefore, being the smearing the same in the four space-time components,*
   * in average the displacement is of type space, because there are          *
   * three space components against the only time one. Other invariant        *
   * distances, like between a light cluster and its reshuffling partner,     *
   * could be a priori equally time-like as space-like.                       *
   ***************************************************************************/
  static Length invDistInitialComponent_BeamVtx          = Length();
  static Length spaceDistInitialComponent_BeamVtx        = Length();
  static Length invDistInitialComponent_Cluster          = Length();
  static Length spaceDistInitialComponent_Cluster        = Length();
  static int    numberInitialComponents                  = 0;  
  static Length invDistFinalHadron_BeamVtx               = Length();
  static Length spaceDistFinalHadron_BeamVtx             = Length();
  static double numberChildHadrons                       = 0;
  static Length invDistInitialCluster_BeamVtx            = Length();
  static Length spaceDistInitialCluster_BeamVtx          = Length();
  // use numberInitialClusters for counting the entries.
  static Length invDistFinalCluster_BeamVtx              = Length();
  static Length spaceDistFinalCluster_BeamVtx            = Length();
  // use numberFinalClusters for counting the entries.
  static Length invDistAnyCluster_BeamVtx                = Length();
  static Length spaceDistAnyCluster_BeamVtx              = Length();
  static double numberClusters                           = 0;
  static Length invDistParentCluster_ChildCluster        = Length();
  static Length spaceDistParentCluster_ChildCluster      = Length();
  static double numberChildClusters                      = 0;
  static Length invDistParentCluster_ChildHadron         = Length();
  static Length spaceDistParentCluster_ChildHadron       = Length();
  // use count numberChildHadrons for counting the entries.
  static Length invDistLightCluster_ReshufflingCluster   = Length();
  static Length spaceDistLightCluster_ReshufflingCluster = Length();
  static double numberReshufflingClusterPartners         = 0;

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
	generator()->log() 
	  << "\t id=" << id << "   " 
	  << getParticleData(id)->PDGName() << "  mass: "
	  << " nominal=" << getParticleData(id)->mass()
	  << " constituent=" << getParticleData(id)->constituentMass()
	  << " iSpin=" << getParticleData(id)->iSpin();
	Energy diff = fabs( getParticleData(id)->constituentMass() - 
			    (  getParticleData(i)->constituentMass() +
			       getParticleData(j)->constituentMass() ) );
        if ( diff > 0.001*GeV ) {
	  generator()->log() 
	    << " <--- NOT EQUAL SUM CONSTITUENT QUARK MASSES! "; 
	} 
	generator()->log() << endl; 
      }
    }

    // Print about spin: just four examples, spin 0, 1/2, 1, 3/2
    // The method  ParticleData::spin()  returns the spin in standard units, 
    // whereas  ParticleData::iSpin()  returns 2J+1 in units of hbar/2.
    generator()->log() << " eta:    iSpin=" 
		       << getParticleData(ParticleID::eta)->iSpin()
                       << "   spin=" 
		       << getParticleData(ParticleID::eta)->spin() << endl
                       << " mu:     iSpin=" 
		       << getParticleData(ParticleID::muminus)->iSpin()
                       << "   spin=" 
		       << getParticleData(ParticleID::muminus)->spin() << endl
                       << " phi:    iSpin=" 
		       << getParticleData(ParticleID::phi)->iSpin()
                       << "   spin=" 
		       << getParticleData(ParticleID::phi)->spin() << endl
                       << " Omega-: iSpin=" 
		       << getParticleData(ParticleID::Omegaminus)->iSpin()
                       << "   spin=" 
		       << getParticleData(ParticleID::Omegaminus)->spin() 
		       << endl;
    
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

  LorentzPoint vertexPosition = ch.currentCollision()->vertex();
  generator()->log() << "ClusterHadronizationHandler::debuggingInfo "
                     << " ===> START DEBUGGING <=== "
		     << "   EventNumber=" << generator()->currentEventNumber() 
		     << endl
                     << "  Total number of clusters: " << clusters.size()
		     << endl << "  Lab vertex of current collision : "
		     << vertexPosition << "  [mm]" << endl;
   
  // Loop over all clusters, and print information that allows to
  // check the logical consistency of such clusters.
                    
  for(ClusterVector::iterator it=clusters.begin(); it!=clusters.end(); it++) {
    tClusterPtr ptrClu = *it;
    invDistAnyCluster_BeamVtx += (ptrClu->vertex() - vertexPosition).mag();
    spaceDistAnyCluster_BeamVtx += (ptrClu->vertex().vect() - vertexPosition.vect()).mag();
    numberClusters++;

    generator()->log() << "  --- Cluster --- " << ptrClu->number() << endl;

    // Information about the components, including sum of charges and momenta.
    // Notice that the charge of the cluster is defined (hence no consistency
    // check is possible) as the sum of the charges of its components.
    generator()->log() << "\t numComponents     = " << ptrClu->numComponents()
		       << endl;
    //int j=0;
    #define pc(i) (ptrClu->particle(i))
    Lorentz5Momentum sumMomentumComponents = Lorentz5Momentum();
    Charge sumChargeComponents = Charge();        
    for(int i = 0; i<ptrClu->numComponents(); i++) {
      generator()->log() << "\t Component " << i << endl 
	                 << "\t \t id = " << pc(i)->id() << "    "
			 << pc(i)->PDGName() 
			 << "     component mass = " << pc(i)->mass() << endl
	                 << "\t \t momentum        = " << pc(i)->momentum() 
			 << endl
	                 << "\t \t position        = " << pc(i)->vertex() 
			 << endl << "\t \t isPerturbative  = " 
			 << ptrClu->isPerturbative(i) << endl
	                 << "\t \t isBeamRemnant   = " 
			 << ptrClu->isBeamRemnant(i) << endl;
      if(pc(i)) {
	generator()->log() << "\t \t HAS pointed Particle :   number = " 
			   << pc(i)->number() << endl
	                   << "\t \t \t positions: labVertex="
	                   << pc(i)->labVertex() 
	                   << "   lifeLength=" 
			   << pc(i)->lifeLength() << endl
	                   << "\t \t \t masses:   constituent="  
			   << pc(i)->data().constituentMass()
			   << "   current=" << pc(i)->mass()
			   << "   invariant=" << pc(i)->momentum().m() 
			   << endl
			   << *pc(i) << endl;
	if(ptrClu->isStatusInitial()) {
	  invDistInitialComponent_BeamVtx += 
	    (pc(i)->vertex() - vertexPosition).mag();
	  invDistInitialComponent_Cluster += 
	    (pc(i)->vertex() - ptrClu->vertex()).mag();
	  spaceDistInitialComponent_BeamVtx += 
	    (pc(i)->vertex().vect() - vertexPosition.vect()).mag();
	  spaceDistInitialComponent_Cluster += 
	    (pc(i)->vertex().vect() - ptrClu->vertex().vect()).mag();
	  numberInitialComponents++;
	}
      } else {
	generator()->log() << "\t \t NO pointed Particle " << endl;  
      }
      sumMomentumComponents += pc(i)->momentum();
      sumChargeComponents += pc(i)->data().charge();      
    }
    generator()->log() << "\t charge (sum components) = " 
		       << sumChargeComponents << endl
		       << "\t sumConstituentMasses    = " 
		       << ptrClu->sumConstituentMasses() 
		       << endl
		       << "\t mass              = " << ptrClu->mass() << endl
		       << "\t momentum          = " << ptrClu->momentum() 
		       << endl << "\t position          = " 
		       << ptrClu->vertex() << endl;    
    Lorentz5Momentum diff = sumMomentumComponents - ptrClu->momentum();
    Energy ediff = fabs( diff.m() );
    if ( ediff > 1e-3*GeV ) {
      generator()->log() 
	<< "\t ***ERROR***: MOMENTUM NOT CONSERVED! " << endl
	<< "\t  --> sumMomentumComponents = " << sumMomentumComponents << endl 
	<< "\t      diff = " << diff << endl;
    } 
    generator()->log() << "\t isBeamCluster     = " << ptrClu->isBeamCluster()
		       << endl;
    if(ptrClu->isBeamCluster()) numberBeamClusters++;

    if(!ptrClu->isAvailable()) {
      generator()->log() << "\t isAvailable     = " << ptrClu->isAvailable() 
			 << endl;
    } else {
      generator()->log() << "\t isStatusInitial   = " 
			 << ptrClu->isStatusInitial() << endl
			 << "\t isReadyToDecay    = " 
			 << ptrClu->isReadyToDecay() << endl
			 << "\t isRedefined       = " 
			 << ptrClu->isRedefined() << endl
			 << "\t hasBeenReshuffled = " 
			 << ptrClu->hasBeenReshuffled() << endl
			 << "\t isStatusFinal     = " 
			 << ptrClu->isStatusFinal() << endl;

      // Statistics information
      if(ptrClu->isStatusInitial()) {
	invDistInitialCluster_BeamVtx += 
	  (ptrClu->vertex() - vertexPosition).mag();	
	spaceDistInitialCluster_BeamVtx += 
	  (ptrClu->vertex().vect() - vertexPosition.vect()).mag();	
	numberInitialClusters++;
      }
      if(ptrClu->isStatusFinal()) {
	invDistFinalCluster_BeamVtx += (ptrClu->vertex()-vertexPosition).mag();
	spaceDistFinalCluster_BeamVtx +=
	  (ptrClu->vertex().vect() - vertexPosition.vect()).mag();	
	numberFinalClusters++;
      }
      if(ptrClu->children().size() == 1) number1HadronClusters++;
      if(ptrClu->hasBeenReshuffled()) {
	numberReshuffledClusters++;
	if(ptrClu->isRedefined()) {
	  invDistLightCluster_ReshufflingCluster += (ptrClu->vertex() -
	    ptrClu->reshufflingPartnerCluster()->vertex()).mag();
	  spaceDistLightCluster_ReshufflingCluster += (ptrClu->vertex().vect()-
	    ptrClu->reshufflingPartnerCluster()->vertex().vect()).mag();
	  numberReshufflingClusterPartners++;
	}
      } else if(ptrClu->isRedefined()) {
	number3ComponentClusters++;
      }

      // Information about the (eventual) parent cluster. 
      int parent = 0;
      for(unsigned int i = 0; i<ptrClu->parents().size(); i++) {
	if(ptrClu->parents()[i]->PDGName() == "Cluster") {
	  parent = ptrClu->parents()[i]->number();
	  break;
	}
      }
      generator()->log() << "\t parentCluster     = " << parent << endl;
      
      // Information about the (eventual) reshuffling partner cluster. 
      int reshuffling = 0;
      if(ptrClu->reshufflingPartnerCluster()) {
	reshuffling = ptrClu->reshufflingPartnerCluster()->number();
      }
      generator()->log() << "\t reshufflingPartner= " << reshuffling << endl;
      
      // Info about the (eventual) children clusters, including the check
      // of conservation of the charge and energy-momentum.
      Lorentz5Momentum sumMomentumChildrenClusters = Lorentz5Momentum();
      Charge sumChargeChildrenClusters = Charge();
      int numChildClus = 0;
      if(ptrClu->children().size()) {
	generator()->log() << "\t childrenClusters  = ";
	for(ParticleVector::const_iterator jt = ptrClu->children().begin();
	    jt != ptrClu->children().end(); ++jt) {	
	  if((*jt)->PDGName() == "Cluster") {
	    numChildClus++;
	    invDistParentCluster_ChildCluster += 
	      (ptrClu->vertex() - (*jt)->vertex()).mag();
	    spaceDistParentCluster_ChildCluster += 
	      (ptrClu->vertex().vect() - (*jt)->vertex().vect()).mag();
	    numberChildClusters++;
	    generator()->log() << (*jt)->number() << "   ";
	    Charge sumChargeChildClustComp = Charge();
	    // The charge of a cluster is defined as the sum of the 
	    // charges of its components. 
	    ClusterPtr cp = dynamic_ptr_cast<ClusterPtr>(*jt);
	    for(int i = 0; i<cp->numComponents(); i++)
	      sumChargeChildClustComp += cp->particle(i)->data().charge();
	    sumChargeChildrenClusters   += sumChargeChildClustComp;
	    sumMomentumChildrenClusters += (*jt)->momentum();
	  }
	}
      }
      if(numChildClus) {
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
      if(ptrClu->children().size() 
	 != static_cast<unsigned int>(numChildClus)) {
	generator()->log() << "\t num childrenHadrons  = " 
			   << ptrClu->children().size()-numChildClus << endl;
	int j = 0;
	for(PVector::const_iterator jt = ptrClu->children().begin();
	      jt != ptrClu->children().end(); ++jt ) {	
	  if((*jt)->PDGName() != "Cluster") {
	    invDistFinalHadron_BeamVtx += 
	      ((*jt)->labVertex() - vertexPosition).mag();
	    spaceDistFinalHadron_BeamVtx += 
	      ((*jt)->labVertex().vect() - vertexPosition.vect()).mag();
	    invDistParentCluster_ChildHadron += 
	      (ptrClu->labVertex() - (*jt)->labVertex()).mag();
	    spaceDistParentCluster_ChildHadron += 
	      (ptrClu->labVertex().vect() - (*jt)->labVertex().vect()).mag();
	    numberChildHadrons++;
	    generator()->log() << "\t Hadron " << ++j << endl 
			       << "\t \t id = " << (*jt)->id() << "   " 
			       << (*jt)->PDGName() 
			       << "   mass = " << (*jt)->mass() << endl
			       << "\t \t momentum        = " 
			       << (*jt)->momentum() << endl
			       << "\t \t labVertex       = " 
			       << (*jt)->labVertex() << endl;
	  // Information about the angles, in Cluster CM and in Lab frame, 
	  // between the hadron and the component of the cluster from which
	  // the hadron came from. Although not nice, we use the fact the 
	  // order of cluster's components is the same as cluster's hadrons.
      	  int k = 0;
	  Lorentz5Momentum pQLab = Lorentz5Momentum();
	  for(int i = 0; i<ptrClu->numComponents(); i++) {
	    if(++k == j) pQLab = ptrClu->particle(i)->momentum();
	  }
	  Lorentz5Momentum pQCm = pQLab;
	  pQCm.boost( - ptrClu->momentum().boostVector() ); 
	  Lorentz5Momentum pHadCm = (*jt)->momentum();
	  pHadCm.boost( - ptrClu->momentum().boostVector() ); 
	  generator()->log() 
	               << "\t \t angles:  CM = "
		       << pHadCm.vect().angle(pQCm.vect().unit()) 
		       << "   Lab = " 
		       << (*jt)->momentum().vect().angle(pQLab.vect().unit()) 
		       << endl;
	  if((*jt)->parents().size()) {
	    generator()->log() << "\t \t num parents = " 
			       << (*jt)->parents().size() << endl;
	    int k = 0;
	    for(tParticleVector::const_iterator kt = (*jt)->parents().begin();
		  kt != (*jt)->parents().end(); ++kt ) {	
	      generator()->log() << "\t \t Parent " << ++k << endl
				 << "\t \t \t id = " << (*kt)->id() << "   " 
				 << (*kt)->PDGName() 
				 << endl
				 << "\t \t \t momentum  = " 
				 << (*kt)->momentum() << endl
				 << "\t \t \t labVertex = " 
				 << (*kt)->labVertex() << endl;
	    }
	  } else {      
	    generator()->log() << "\t \t NO parents for this hadron" << endl;
	  }	
	  sumChargeChildrenHadrons   += (*jt)->data().charge();
	  sumMomentumChildrenHadrons += (*jt)->momentum();
	  generator()->log() << endl 
			     << "\t  --> sumChargeChildrenHadrons = " 
			     << sumChargeChildrenHadrons << endl	  
			     << "\t  --> sumMomentumChildrenHadrons = " 
			     << sumMomentumChildrenHadrons << endl;   
	  } else {      
	    generator()->log() << "\t NO childrenHadrons " << endl;
	  }
	}
      }

      if(ptrClu->children().size()) {
	Charge delta_ch = fabs(getParticleData(ParticleID::d)->charge())/10.0;
	if(fabs(sumChargeChildrenClusters + sumChargeChildrenHadrons 
		  - sumChargeComponents) > delta_ch ) {
	  generator()->log() <<"\t ***ERROR***: CHARGE NOT CONSERVED! "<< endl;
	} 
	diff = sumMomentumChildrenClusters + sumMomentumChildrenHadrons - 
	  ptrClu->momentum();
	ediff = fabs(diff.m());
	if(ediff > 1e-3*GeV) {
	  if(ptrClu->hasBeenReshuffled()) {
	    generator()->log() << "\t ***WARNING***: MOMENTUM NOT CONSERVED"
			       << " BUT RESHUFFLING OCCURED!" << endl; 
	  } else {
	    generator()->log() << "\t ***ERROR***: MOMENTUM NOT CONSERVED! " 
			       << endl;
	  }
	  generator()->log() 
	    << "\t  --> sumMomentumChildren = " 
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
  
  generator()->log() 
    << " ---------- Begin analysis of produced Hadrons ---------- " << endl;
  /*for ( map<tCluPtr,int>::const_iterator it = orderingMap.begin();
    it != orderingMap.end(); ++it ) {*/
  for(ClusterVector::iterator it=clusters.begin(); it!=clusters.end(); it++) {
    tClusterPtr ptrClu = (*it);  
    if(ptrClu->isAvailable() && ptrClu->isStatusInitial()) {
      generator()->log() << "  Cluster (Initial) " << ptrClu->number() << endl;
      bool reshuffled = false;
      vector<tPPtr> vecPartonParents;
      for(int i = 0; i<ptrClu->numComponents(); i++)
	vecPartonParents.push_back(ptrClu->particle(i));
      tParticleVector vecChildrenHadrons;
      tClusterVector theStack;  // to avoid recursion
      theStack.push_back(ptrClu);
      do {
	tClusterPtr ptrClu = theStack.back();
	theStack.pop_back();        
	for(PVector::const_iterator jt = ptrClu->children().begin();
	    jt != ptrClu->children().end(); ++jt) {
	  if((*jt)->PDGName() == "Cluster") 
	    theStack.push_back(dynamic_ptr_cast<ClusterPtr>(*jt));
	  else vecChildrenHadrons.push_back(*jt);
	}
	if ( ptrClu->hasBeenReshuffled() ) {
	  reshuffled = true;
	}
      } while(!theStack.empty());    
      Lorentz5Momentum sumMomentumPartons = Lorentz5Momentum();
      Charge sumChargePartons = Charge();
      generator()->log() << "\t made of the following  " 
			 << vecPartonParents.size()
	                 << "  partons " << endl;
      for(tParticleVector::const_iterator jt = vecPartonParents.begin();
	  jt != vecPartonParents.end(); ++jt) {
	generator()->log() << "\t \t id = " << (*jt)->id() << "   " 
			   << (*jt)->PDGName() 
			   << "   number = " << (*jt)->number() << "   " 
			   << (*jt)->momentum() << endl;
	sumMomentumPartons += (*jt)->momentum();
        sumChargePartons += (*jt)->data().charge();
      }
      Lorentz5Momentum sumMomentumHadrons = Lorentz5Momentum();
      Charge sumChargeHadrons = Charge();
      generator()->log() << "\t has produced the following  " 
			 << vecChildrenHadrons.size()
	                 << "  hadrons " << endl;
      for(tParticleVector::const_iterator jt = vecChildrenHadrons.begin();
	  jt != vecChildrenHadrons.end(); ++jt) {
	generator()->log() << "\t \t id = " << (*jt)->id() << "   " 
			   << (*jt)->PDGName() 
			   << "   " << (*jt)->momentum() << endl;
	sumMomentumHadrons += (*jt)->momentum();
        sumChargeHadrons += (*jt)->data().charge();

        // Add here all other statistical information
        eventHadronMultiplicity++;
        if(hadrons.find((*jt)->id()) != hadrons.end()) {
	  hadrons.find((*jt)->id())->second = hadrons.find((*jt)->id())->second
	    + 1;
	} else {
	  hadrons.insert(pair<long,double>((*jt)->id(), 1.0));
	}
        if(CheckId::isMeson((*jt)->id())) numberMesons++;  
        if(CheckId::hasBeauty((*jt)->id())) {
	  numberBHadrons++;
	} else if(CheckId::hasCharm((*jt)->id())) {
	  numberCharmHadrons++;
	} else if(CheckId::hasStrangeness((*jt)->id(), rnd())) {
	  numberStrangeHadrons++;
	} else {
	  numberLightHadrons++;
	}
        if((*jt)->id() == ParticleID::piplus )  numberPlusPions++;  
        if((*jt)->id() == ParticleID::pi0 )     numberNeutralPions++;  
        if((*jt)->id() == ParticleID::piminus ) numberMinusPions++; 
      }
      Charge delta_ch = fabs(getParticleData(ParticleID::d)->charge())/10.0; 
      if(fabs(sumChargePartons - sumChargeHadrons) > delta_ch) {
	generator()->log() << "\t ***ERROR***: CHARGE NOT CONSERVED! " << endl;
      } 
      generator()->log() << "\t Total Charge Partons = " << sumChargePartons 
			 << endl
			 << "\t Total Charge Hadrons = " << sumChargeHadrons 
			 << endl;
      Lorentz5Momentum diff = sumMomentumPartons - sumMomentumHadrons;
      Energy ediff = fabs(diff.m());
      if(ediff > 1e-3*GeV) {
	generator()->log() << "\t ***ERROR***: MOMENTUM NOT CONSERVED!"
			   << "  reshuffled=" << reshuffled << endl;
      }
      generator()->log() 
	<< "\t Total Momentum Partons = " << sumMomentumPartons << endl
	<< "\t Total Momentum Hadrons = " << sumMomentumHadrons << endl
	<< "\t   ---> Diff Momentum = " << diff << endl;
    } // end of if isAvailable && isStatusInitial
  } // endl loop over the map
  generator()->log() 
    << " ---------- End analysis of produced Hadrons ---------- " << endl;

  generator()->log() << "ClusterHadronizationHandler::debuggingInfo "
                     << " ===> END DEBUGGING <=== " << endl;

  // At the end of the last event, print out statistics information.
  if(generator()->currentEventNumber() == generator()->N()) {
    generator()->log() << endl << endl 
		       << "ClusterHadronizationHandler::debuggingInfo "
		       << " ===> BEGIN ***FINAL STATISTICS*** <=== " 
		       << endl << endl;
    double dN = double(generator()->N());
    int i=0;
    for(map<long,double>::const_iterator iter = hadrons.begin();
	iter != hadrons.end(); ++iter ) {
      generator()->log() << "\t" << ++i << ")  id=" << iter->first << "   "
	                 << (getParticleData(iter->first) ? 
			     getParticleData(iter->first)->PDGName()
			     : "UNKNOWN") 
			 << "   " << iter->second / dN 
			 << " per event;   fraction="
			 <<  100.0 * iter->second/eventHadronMultiplicity 
			 << " %" << endl;
    }
    generator()->log() 
      << endl
      << "\t Counts for position statistics " << endl
      << "\t \t   numberInitialComponents =" << numberInitialComponents << endl
      << "\t \t   numberInitialClusters   =" << numberInitialClusters << endl
      << "\t \t   numberClusters          =" << numberClusters << endl
      << "\t \t   numberFinalClusters     =" << numberFinalClusters << endl
      << "\t \t   numberChildClusters     =" << numberChildClusters << endl
      << "\t \t   numberChildHadrons      =" << numberChildHadrons << endl
      << "\t \t   numberReshufflingClusterParters =" 
      << numberReshufflingClusterPartners << endl 
      << endl
      << "\t Invariant (space-time) Distance (in millimeters) : " << endl
      << "\t \t | initial component - beam vertex | = " 
      << (numberInitialComponents 
	  ? invDistInitialComponent_BeamVtx / numberInitialComponents
	  : 0.0) << endl
      << "\t \t | initial component - cluster | = " 
      << (numberInitialComponents 
	  ? invDistInitialComponent_Cluster / numberInitialComponents  
	  : 0.0) << endl
      << "\t \t | initial cluster - beam vertex | = " 
      << (numberInitialClusters 
	  ? invDistInitialCluster_BeamVtx / numberInitialClusters  
	  : 0.0) << endl
      << "\t \t | any cluster - beam vertex | = " 
      << (numberClusters 
	  ? invDistAnyCluster_BeamVtx / numberClusters 
	  : 0.0) << endl
      << "\t \t | final cluster - beam vertex | = " 
      << (numberFinalClusters 
	  ? invDistFinalCluster_BeamVtx / numberFinalClusters  
	  : 0.0) << endl
      << "\t \t | final hadron - beam vertex | = " 
      << (numberChildHadrons 
	  ? invDistFinalHadron_BeamVtx / numberChildHadrons  
	  : 0.0) << endl
      << "\t \t | parent cluster - child cluster | = " 
      << (numberChildClusters 
	  ? invDistParentCluster_ChildCluster / numberChildClusters  
	  : 0.0) << endl
      << "\t \t | parent cluster - child hadron | = " 
      << (numberChildHadrons 
	  ? invDistParentCluster_ChildHadron / numberChildHadrons  
	  : 0.0) << endl
      << "\t \t | light cluster - reshuffling partner cluster | = " 
      << (numberReshufflingClusterPartners
	  ? invDistLightCluster_ReshufflingCluster 
	  / numberReshufflingClusterPartners : 0.0) << endl
      << endl
      << "\t Lab space distance (abs) (in millimeters) : " << endl
      << "\t \t | initial component - beam vertex | = " 
      << (numberInitialComponents 
	  ? spaceDistInitialComponent_BeamVtx / numberInitialComponents  
	  : 0.0) << endl
      << "\t \t | initial component - cluster | = " 
      << (numberInitialComponents 
	  ? spaceDistInitialComponent_Cluster / numberInitialComponents  
	  : 0.0) << endl
      << "\t \t | initial cluster - beam vertex | = " 
      << (numberInitialClusters 
	  ? spaceDistInitialCluster_BeamVtx / numberInitialClusters  
	  : 0.0) << endl
      << "\t \t | any cluster - beam vertex | = " 
      << (numberClusters 
	  ? spaceDistAnyCluster_BeamVtx / numberClusters 
	  : 0.0) << endl
      << "\t \t | final cluster - beam vertex | = " 
      << (numberFinalClusters 
	  ? spaceDistFinalCluster_BeamVtx / numberFinalClusters  
	  : 0.0) << endl
      << "\t \t | final hadron - beam vertex | = " 
      << (numberChildHadrons 
	  ? spaceDistFinalHadron_BeamVtx / numberChildHadrons  
	  : 0.0) << endl
      << "\t \t | parent cluster - child cluster | = " 
      << (numberChildClusters 
	  ? spaceDistParentCluster_ChildCluster / numberChildClusters  
	  : 0.0) << endl
      << "\t \t | parent cluster - child hadron | = " 
      << (numberChildHadrons 
	  ? spaceDistParentCluster_ChildHadron / numberChildHadrons  
	  : 0.0) << endl
      << "\t \t | light cluster - reshuffling partner cluster | = " 
      << (numberReshufflingClusterPartners
	  ? spaceDistLightCluster_ReshufflingCluster 
	  / numberReshufflingClusterPartners : 0.0) << endl
      << endl;
    generator()->log() 
      << endl
      << "\t Clusters:  BeamCluster  average number  = " 
      << numberBeamClusters / dN << endl
      << "\t            Initial      average number  = " 
      << numberInitialClusters / dN << endl
      << "\t            Final        average number  = " 
      << numberFinalClusters / dN << endl
      << "\t            1-Hadron     average number  = " 
      << number1HadronClusters / dN << "  fraction="
      << number1HadronClusters / numberFinalClusters << endl
      << "\t            Reshuffled   average number  = " 
      << numberReshuffledClusters / dN << "  fraction="
      << numberReshuffledClusters / numberFinalClusters << endl
      << "\t            3-Component  average number  = " 
      << number3ComponentClusters / dN << "  fraction="
      << number3ComponentClusters / numberFinalClusters <<endl<<endl 
      << "\t Average Hadron Multiplicity = " 
      << eventHadronMultiplicity / dN << endl
      << "\t Mesons:  " << numberMesons / dN << " per event;    fraction="
      << 100.0 * numberMesons / eventHadronMultiplicity << " %" << endl
      << "\t Baryons: " << (eventHadronMultiplicity-numberMesons) / dN 
      << " per event;    fraction="
      << 100.0 * ( 1.0 - numberMesons / eventHadronMultiplicity ) 
      << " %" << endl
      << "\t b-Hadrons    = " 
      << numberBHadrons / dN << " per event;    fraction="
      << 100.0 * numberBHadrons / eventHadronMultiplicity << " %" << endl
      << "\t c-Hadrons    = " 
      << numberCharmHadrons / dN << " per event;    fraction="
      << 100.0 * numberCharmHadrons / eventHadronMultiplicity << " %" << endl
      << "\t s-Hadrons    = " 
      << numberStrangeHadrons / dN << " per event;    fraction="
      << 100.0 * numberStrangeHadrons / eventHadronMultiplicity << " %" << endl
      << "\t u,d-Hadrons  = " 
      << numberLightHadrons / dN << " per event;    fraction="
      << 100.0 * numberLightHadrons / eventHadronMultiplicity << " %" << endl
      << "\t \t pi+  = " 
      << numberPlusPions / dN << " per event;    fraction="
      << 100.0 * numberPlusPions / eventHadronMultiplicity << " %" << endl
      << "\t \t pi0  = " 
      << numberNeutralPions / dN << " per event;    fraction="
      << 100.0 * numberNeutralPions / eventHadronMultiplicity << " %" << endl
      << "\t \t pi-  = "
      << numberMinusPions / dN << " per event;    fraction=" 
      << 100.0 * numberMinusPions / eventHadronMultiplicity << " %" << endl;
    generator()->log() 
      << endl << "ClusterHadronizationHandler::debuggingInfo "
      << " ===> END ***FINAL STATISTICS*** <=== " << endl << endl;
  }
}
