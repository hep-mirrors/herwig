// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ClusterFissioner class.
//

#include "ClusterFissioner.h"
#include "Pythia7/Interface/ClassDocumentation.h"
#include "Pythia7/Interface/Reference.h" 
#include "Pythia7/Interface/Parameter.h"
#include "Pythia7/Persistency/PersistentOStream.h"
#include "Pythia7/Persistency/PersistentIStream.h"
#include "Pythia7/PDT/EnumParticles.h"
#include "Pythia7/Repository/EventGenerator.h"
#include "Pythia7/EventRecord/Collision.h"
#include "Herwig++/Utilities/Kinematics.h"
#include "Herwig++/Utilities/HwDebug.h"
#include "Herwig++/Utilities/CheckId.h"
#include "Cluster.h"


using namespace Herwig;
// using namespace Pythia7;


ClusterFissioner::~ClusterFissioner() {}


void ClusterFissioner::persistentOutput(PersistentOStream & os) const {
  os << _hadronsSelector << _globalParameters
     << _ClMax << _ClPow << _PSplt1 << _PSplt2;
}


void ClusterFissioner::persistentInput(PersistentIStream & is, int) {
  is >> _hadronsSelector >> _globalParameters
     >> _ClMax >> _ClPow >> _PSplt1 >> _PSplt2;
}


ClassDescription<ClusterFissioner> ClusterFissioner::initClusterFissioner;
// Definition of the static class description member.


void ClusterFissioner::Init() {

  static ClassDocumentation<ClusterFissioner> documentation
    ("Class responsibles for chopping up the clusters");

  static Reference<ClusterFissioner,HadronsSelector> 
    interfaceHadronsSelector("HadronsSelector", 
                             "A reference to the HadronsSelector object", 
                             &Herwig::ClusterFissioner::_hadronsSelector,
			     false, false, true, false);

  static Reference<ClusterFissioner,GlobalParameters> 
    interfaceGlobalParameters("GlobalParameters", 
			      "A reference to the GlobalParameters object", 
			      &Herwig::ClusterFissioner::_globalParameters,
			      false, false, true, false);
  
  static Parameter<ClusterFissioner,Energy>
    interfaceClMax ("ClMax","cluster max mass  (unit [GeV])",
                    &ClusterFissioner::_ClMax, GeV, 3.35*GeV, 0.0*GeV, 10.0*GeV);
  static Parameter<ClusterFissioner,double>
    interfaceClPow ("ClPow","cluster mass exponent",
                    &ClusterFissioner::_ClPow, 0, 2.0, 0.0, 10.0);
  static Parameter<ClusterFissioner,double>
    interfacePSplt1 ("PSplt1","cluster mass splitting param for u,d,s,c",
                    &ClusterFissioner::_PSplt1, 0, 1.0, 0.0, 10.0);
  static Parameter<ClusterFissioner,double>
    interfacePSplt2 ("PSplt2","cluster mass splitting param for b",
                    &ClusterFissioner::_PSplt2, 0, 1.0, 0.0, 10.0);

}


void ClusterFissioner::fission(const StepPtr &pstep, ClusterVector &clusters)
{

  // Loop over the (input) collection of cluster pointers, and store in 
  // the vector  vecSplitCluPtr  all the clusters that need to be split
  // (these are beam clusters, if soft underlying event is off, and 
  //  heavy non-beam clusters).
  vector<tClusterPtr> splitClusters; 
  int numBeamClusters = 0;
  for(ClusterVector::iterator it = clusters.begin() ; 
      it != clusters.end() ; ++it) {
 
    // Skip 3-component clusters that have been redefined (as 2-component clusters).
    // or not available clusters. The latter check is indeed redundant now, 
    // but it is used for possible future extensions in which, for some
    // reasons, some of the clusters found by ClusterFinder are tagged
    // straight away as not available.
    if((*it)->isRedefined() || !(*it)->isAvailable()) continue;
    if((*it)->isBeamCluster()) {
      numBeamClusters++;
      // Tagged as not available the beam clusters if soft underlying event if on.
      if ( _globalParameters->isSoftUnderlyingEventON() ) {
	(*it)->isAvailable(false);
      } else {
	splitClusters.push_back(*it);
      }
    } else {      
      // If the cluster is heavy add it to the vector of clusters to be split.
      if(pow((*it)->mass() , _ClPow) > 
	 pow(_ClMax, _ClPow) + pow((*it)->sumConstituentMasses(), _ClPow)) {
	splitClusters.push_back(*it);
      }
    }
  } // end loop over clusters
  
  // Safety check, usually skipped.
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::minimal_Hadronization ) {
    // There should be no more than at most 4 beam-remnant clusters:
    // take the case of gluon-gluon hard scattering, we will have two 
    // beam parton remnants (one (anti-) quark and one (anti-) diquark) 
    // for each of the two beams, which can formed max 4 beam remnant clusters.
    if ( numBeamClusters > 4 ) {
      generator()->logWarning( Exception("ClusterFissioner::fission "
		  "*** numBeamClusters > 4 *** ", Exception::warning));
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization ) {    
	generator()->log() << "         ===>" << " numBeamClusters = " 
			   << numBeamClusters << endl << endl;
      }
    }
    // Now more debugging infos
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization ) {
      generator()->log() << "ClusterFissioner::fission   ===> START DEBUGGING <=== "
                         << "   EventNumber=" << generator()->currentEventNumber() << endl
                         << "   Number Initial Clusters = " << clusters.size() << endl
			 << "   Number of Beam Clusters = " << numBeamClusters << endl;
      //      for(unsigned int i = 0; i<clusters.size(); i++)
      //generator()->log() << *clusters[i] << endl;
      generator()->log() << "ClusterFissioner::fission   ===> END DEBUGGING "
			 << "<=== " << endl; 
    }  
  }
  
  // Loop over the vector of clusters to be split, vecSplitCluPtr,
  // and pass each element to the method  cut. This latter method modifies, 
  // by adding the cluster children of such cluster, the collection of cluster 
  // pointers  collecCluPtr : that's why  we cannot loop directly on the  collecCluPtr  
  // and we need instead the vector  vecHeavyCluPtr. In fact, it is not allowed to 
  // modify a STL container during the iteration over it!  
  for (vector<tClusterPtr>::const_iterator iter = splitClusters.begin() ; 
       iter != splitClusters.end() ; ++iter) {
    //cout << "Calling cut on " << **iter << endl;
    cut(*iter, pstep, clusters);
  }
}

void ClusterFissioner::cut(tClusterPtr cluPtr, const StepPtr &pstep, 
   			   ClusterVector &clusters) {

  // This method does the splitting of the cluster pointed by  cluPtr
  // and "recursively" by all of its cluster children, if heavy. All of these
  // new children clusters are added (indeed the pointers to them) to the
  // collection of cluster pointers  collecCluPtr. The method works as follows.
  // Initially the vector vecCluPtr contains just the input pointer to the
  // cluster to be split. Then it will be filled "recursively" by all
  // of the cluster's children that are heavy enough to require, in their turn,
  // to be split. In each loop, the last element of the vector vecCluPtr is 
  // considered (only once because it is then removed from the vector).
  // This approach is conceptually recursive, but avoid the overhead of
  // a concrete recursive function. Furthermore it requires minimal changes
  // in the case that the fission of an heavy cluster could produce more
  // than two cluster children as assumed now. 
  vector<tClusterPtr> clusterStack; 
  clusterStack.push_back(cluPtr);
  while ( ! clusterStack.empty() ) {

    tClusterPtr iCluPtr = clusterStack.back();  // consider element on the back
    clusterStack.pop_back();                    // delete element on the back 
 
    // We need to require (at least at the moment, maybe in the future we 
    // could change it) that the cluster has exactly two components, 
    // because otherwise we don't know how to deal with the kinematics of 
    // the splitting of a heavy cluster. If this is not the case, then
    // send a warning because it is not suppose to happen, and then
    // do nothing with (ignore) such cluster.
    if ( iCluPtr->numComponents() != 2 ) {
      generator()->logWarning( Exception("ClusterFissioner::cut "
	     "***Still cluster with not exactly 2 components*** ", 
					 Exception::warning) );
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization ) {    
	generator()->log() << "         ===>" << " num components = " 
			   << iCluPtr->numComponents() << endl << endl;
      }
      continue;
    }
    
    // Extract the id and particle pointer of the two components of the cluster.
    //tCompPtr compPtr1 = tCompPtr(), compPtr2 = tCompPtr();
    long idQ1 = 0, idQ2 = 0;
    tPPtr ptrQ1 = tPPtr(), ptrQ2 = tPPtr();
    ptrQ1 = iCluPtr->particle(0);
    if(ptrQ1) idQ1 = ptrQ1->id();
    ptrQ2 = iCluPtr->particle(1);
    if(ptrQ2) idQ2 = ptrQ2->id();

    // Sanity check (normally skipped) to control that the two components of a
    // cluster are consistent, that is they can form a meson or a baryon.
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::minimal_Hadronization ) {
      if (!ptrQ1  ||  !ptrQ2 || !ptrQ1->dataPtr() || !ptrQ2->dataPtr()) {
	generator()->logWarning( Exception("ClusterFissioner::cut "
		       "***Cluster with inconsistent components***", 
					   Exception::warning));
	if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization ) {    
	  generator()->log() << "         ===>" << " idQ1=" << idQ1 
			     << " compPtr1=" << ptrQ1 << " compPtr2=" << ptrQ2 
			     << endl << endl;
	}
      }
      if(!CheckId::canBeMeson(idQ1,idQ2) && !CheckId::canBeBaryon(idQ1,idQ2)) {
	generator()->logWarning( Exception("ClusterFissioner::cut "
	 "***The two components of the cluster are inconsistent***", 
					   Exception::warning) );
	if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization ) {    
	  generator()->log() << "         ===>"  << " idQ1=" << idQ1 
			     << " idQ2=" << idQ2 << " " << CheckId::canBeMeson(idQ1,idQ2) << " " << CheckId::canBeBaryon(idQ1,idQ2) << "iCluPtr->#() "<<iCluPtr->number() << endl << endl;
	}
      }
    }

    // Randomly select a flavor necessary to create two clusters by splitting 
    // the original one. Notice that new id is initially > 0, but then it sign
    // must be consistently defined in order to produce either a meson or baryon.
    long idNew = drawnNewFlavour();      // draw the new flavour (idNew > 0)
    if(!CheckId::canBeMeson(idQ1,-idNew) && !CheckId::canBeBaryon(idQ1,-idNew))
      idNew = -idNew;

    // Determine the masses of the two children clusters.
    // Notice that the exponent for the assumed distribution of cluster masses
    // is different in the case of b (anti-)quark (or (anti-)diquark with b-flavour)
    Energy Mclu = iCluPtr->mass(), Mclu1 = Energy(), Mclu2 = Energy();
    Energy m1 = ptrQ1->data().constituentMass();
    Energy m2 = ptrQ2->data().constituentMass();
    Energy m  = getParticleData(abs(idNew))->constituentMass();

    // Do not split in the case there is no phase space available
    // (it happens sometimes for clusters with b-flavour)
    if(Mclu <  m1+m + m2+m) continue;

    double exponent1=_PSplt1, exponent2=_PSplt1;
    if(CheckId::hasBeauty(idQ1)) exponent1 = _PSplt2;
    if(CheckId::hasBeauty(idQ2)) exponent2 = _PSplt2;

    // Draw the masses: for normal, non-beam clusters a power-like mass distribution
    // is used, whereas for beam clusters a fast-decreasing esponential mass 
    // distribution is used instead (to avoid many iterative splitting which could 
    // produce an unphysical large transverse energy from a supposed soft beam
    // remnant process).
    // More precisely, the choice of which mass distribution should be used for each 
    // of the two cluster children is dictated by the integer  iRemnant, defined
    // as follows:
    // i)   iRemnant == 0  then both Mclu1 and Mclu2 are extracted by the 
    //                     power-like distribution;
    // ii)  iRemnant == 1  then Mclu1 is extracted from the exponential distribution,
    //                     whereas Mclu2 is extracted from the power-like distribution;
    // iii) iRemnant == 2  then Mclu1 is extracted from the power-like distribution,
    //                     whereas Mclu2 is extracted from the exponential distribution;
    // iv)  otherwise      both Mclu1 and Mclu2 are extracted from the exponential
    //                     distribution
    // Case i) is the standard one, always used for non-beam clusters;
    // case ii) and iii) occur for beam clusters with only one beam remnant and
    // when the flag _IOpRem is 1 (which is also its default value) which means 
    // that the cluster child containing the beam remnant should have "soft" mass
    // whereas the other one should have "normal" mass;
    // case iv), finally, occur for beam clusters with two beam remnants 
    // (regardless of the value of _IOpRem) or when _IOpRem is 0 which means 
    // that both cluster children masses should be soft.
    //      NB) Furthermore, although not needed at the moment, in the case iv)
    //          the integer iRemnant carries more information than simply inform 
    //          that the exponential distribution should be used for both children:
    //            iRemnant = 10 or 20  then both components are beam remnants;
    //            iRemnant = 11        component 1 is the beam remnant;
    //            iRemnant = 12        component 2 is the beam remnant.
    // If, during the drawning of candidate masses, too many attempts fail 
    // (because the phase space available is tiny) then give up (the cluster 
    //  is not plit).
    int iRemnant = 0;
    if (cluPtr->isBeamCluster()) {
      if ( cluPtr->isBeamRemnant(0) ) {
        iRemnant = 1;
	if (cluPtr->isBeamRemnant(1)) iRemnant = 10;
      } else if (cluPtr->isBeamRemnant(1)) iRemnant = 2;
      if ( _IOpRem == 0 ) iRemnant += 10;
    }
    if (! drawnChildrenMasses(Mclu,m1,m2,m,Mclu1,Mclu2,
	                      exponent1,exponent2,_BtClM,iRemnant) ) {
      continue;
    } 
 
    // New (not present in Fortran Herwig):
    // check whether the fragment masses  Mclu1  and  Mclu2  are above the 
    // threshold for the production of the lightest pair of hadrons with the 
    // right flavours. If not, then set by hand the mass to the lightest 
    // single hadron with the right flavours, in order to solve correctly
    // the kinematics, and (later in this method) create directly such hadron
    // and add it to the children hadrons of the cluster that undergoes the
    // fission (i.e. the one pointed by iCluPtr). Notice that in this special
    // case, the heavy cluster that undergoes the fission has one single 
    // cluster child and one single hadron child. We prefer this approach,
    // rather than to create a light cluster, with the mass set equal to
    // the lightest hadron, and let then the class LightClusterDecayer to do 
    // the job to decay it to that single hadron, for two reasons: 
    // First, because the sum of the masses of the two constituents can be, 
    // in this case, greater than the mass of that hadron, hence it would
    // be impossible to solve the kinematics for such two components, and
    // therefore we would have a cluster whose components are undefined.
    // Second, the algorithm is faster, because it avoids the reshuffling
    // procedure that would be necessary if we used LightClusterDecayer
    // to decay the light cluster to the lightest hadron.   
    bool decayOneHadronClu1 = false;
    if ( Mclu1 < _hadronsSelector->massLightestHadronsPair(idQ1,-idNew)) { 
      Mclu1 =  _hadronsSelector->massLightestHadron(idQ1,-idNew);          
      decayOneHadronClu1 = true;
    }
    bool decayOneHadronClu2 = false;
    if ( Mclu2 < _hadronsSelector->massLightestHadronsPair(idQ2,idNew)) { 
      Mclu2 =  _hadronsSelector->massLightestHadron(idQ2,idNew);           
      decayOneHadronClu2 = true;
    }
    // Check if the decay kinematics is still possible: if not then 
    // force the one-hadron decay for the other cluster as well.
    if ( Mclu1 + Mclu2  >  Mclu ) {
      if ( ! decayOneHadronClu1 ) {
	Mclu1 =  _hadronsSelector->massLightestHadron(idQ1, -idNew); 
	decayOneHadronClu1 = true;	
      } else if ( ! decayOneHadronClu2 ) {
	Mclu2 =  _hadronsSelector->massLightestHadron(idQ2, idNew);
	decayOneHadronClu2 = true;
      }
      // Sanity check (normally skipped) to see if at this point we still have
      // that the sum of the masses of the two children is above the one of
      // their parent: it should never happen!
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::minimal_Hadronization ) {
	if ( Mclu1 + Mclu2  >  Mclu ) {
	  generator()->logWarning( Exception("ClusterFissioner::cut "
					     "***Impossible too low mass for an heavy cluster***",
					     Exception::warning) );
	  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization ) {    
	    generator()->log() << "         ===>" 
			       << " \t Cluster split : mass=" << Mclu
			       << " components ids= " << idQ1 << " " << idQ2 
			       << endl << " \t Clu1 : mass=" << Mclu1 << " " 
			       << " components ids= " << idQ1 << " " << -idNew 
                               <<"   decayOneHadronClu1=" << decayOneHadronClu1
			       << endl << " \t Clu2 : mass=" << Mclu2 << " " 
			       << " components ids= " << idQ2 << " " << idNew 
                               <<"   decayOneHadronClu2=" << decayOneHadronClu2
			       << endl << endl;
	  }
	}
      }
    }
 
    // Determined the (5-components) momenta (all in the LAB frame)
    Lorentz5Momentum pClu = iCluPtr->momentum(); // known
    Lorentz5Momentum p0Q1 = ptrQ1->momentum();// known (momentum of Q1 before fission)
    Lorentz5Momentum pClu1, pClu2, pQ1, pQone, pQtwo, pQ2; //unknown
    pClu1.setMass(Mclu1);
    pClu2.setMass(Mclu2);
    pQ1.setMass(m1);
    pQ2.setMass(m2);
    pQone.setMass(m); 
    pQtwo.setMass(m);
    
    calculateKinematics(pClu,p0Q1,decayOneHadronClu1,decayOneHadronClu2, // in
			pClu1,pClu2,pQ1,pQone,pQtwo,pQ2);                // out

    // Determine the positions of the two children clusters.
    LorentzPoint positionClu1 = LorentzPoint();
    LorentzPoint positionClu2 = LorentzPoint();
    calculatePositions(pClu, iCluPtr->vertex(), pClu1, pClu2,   // input
		       positionClu1, positionClu2 );            // output
   
    // The previous methods have determined the kinematics and positions
    // of C -> C1 + C2. 
    // In the case that one of the two product is light, that means either
    // decayOneHadronClu1 or decayOneHadronClu2 is true, then the momenta
    // of the components of that light product have not been determined,
    // and a (light) cluster will not be created: the heavy father cluster
    // decays, in this case, into a single (not-light) cluster and a
    // single hadron. In the other, "normal", cases the father cluster
    // decays into two clusters, each of which has well defined components.
    // Notice that, in the case of components which point to particles, the
    // momenta of the components is properly set to the new values, whereas
    // we do not change the momenta of the pointed particles, because we 
    // want to keep all of the information (that is the new momentum of a
    // component after the splitting, which is contained in the _momentum
    // member of the Component class, and the (old) momentum of that component
    // before the splitting, which is contained in the momentum of the
    // pointed particle). Please not make confusion of this only apparent
    // inconsistency!

    if(decayOneHadronClu1) {
      long idhad = _hadronsSelector->lightestHadron(idQ1,-idNew);      
      PPtr ptrhad = getParticle(idhad);  // create the hadron.
      if(!ptrhad) {
	generator()->logWarning(Exception("ClusterFissioner::cut "
		"***Cannot create a particle with specified id***", 
					  Exception::warning));
	if(HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization) {    
	  generator()->log() << "         ===>" 
			     << " idQ1=" << idQ1 << " -idNew=" << -idNew 
			     << " idhad=" << idhad << endl << endl;
	}
      } else {
	ptrhad->set5Momentum(pClu1);
	ptrhad->setLabVertex(positionClu1);
	pstep->addDecayProduct(iCluPtr,ptrhad);
      }

    } else {
      // Create the two new clusters with the proper flavour. Notice the sign 
      // of idNew is necessary to produce consistent clusters (see previous 
      // comment)
      ClusterPtr ptrClu1;
      PPtr ptr2 = getParticle(-idNew);
      pstep->addIntermediate(ptr2);
      if (ptrQ1) {
	ptrClu1 = new_ptr(Cluster(ptrQ1,ptr2));
	ptr2->addChild(ptrClu1);
      } else {
	PPtr ptr1 = getParticle(idQ1);
	ptrClu1 = new_ptr(Cluster(ptr1,ptr2));
	ptr1->addChild(ptrClu1);
	ptr2->addChild(ptrClu1);
	pstep->addIntermediate(ptr1);
      }
      ptrClu1->set5Momentum(pClu1);
      ptrClu1->setVertex(positionClu1);

      for(int i = 0; i<ptrClu1->numComponents(); i++) {
	if(ptrClu1->particle(i)->id() == idQ1) 
	  ptrClu1->particle(i)->set5Momentum(pQ1);
	else if(ptrClu1->particle(i)->id() == -idNew) 
	  ptrClu1->particle(i)->set5Momentum(pQone);
      }

      // Set the parent/children relationship between the clusters, and add
      // the children clusters (indeed the pointers to them) to the collection
      // of cluster pointers  collecCluPtr. Finally, add also in the vector
      // vecCluPtr  those children that are heavy enough to be split.
      //iCluPtr->addChild(ptrClu1);
      pstep->addDecayProduct(iCluPtr, ptrClu1);
      clusters.push_back(ptrClu1);
      if(pow(ptrClu1->mass(), _ClPow) > 
	 pow(_ClMax, _ClPow) + pow(ptrClu1->sumConstituentMasses(), _ClPow)) {
	clusterStack.push_back(ptrClu1);
      } 

    }
  
    if(decayOneHadronClu2) {

      long idhad = _hadronsSelector->lightestHadron(idQ2, idNew);      
      PPtr ptrhad = getParticle( idhad );
      if ( ! ptrhad ) {
	generator()->logWarning( Exception("ClusterFissioner::cut "
					   "***Cannot create a particle with specified id***", 
					   Exception::warning) );
	if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization ) {    
	  generator()->log() << "         ===>" 
			     << " idQ2=" << idQ2 << "  idNew=" <<  idNew 
			     << " idhad=" << idhad << endl << endl;
	}
      } else {
	ptrhad->set5Momentum(pClu2);
	ptrhad->setLabVertex(positionClu2);
	//iCluPtr->addChild(ptrhad);
	pstep->addDecayProduct(iCluPtr,ptrhad);
      }

    } else {
      ClusterPtr ptrClu2 = ClusterPtr();
      PPtr ptr1;
      PPtr ptr2 = getParticle(idNew);
      pstep->addIntermediate(ptr2);
      //ParticleVector dp;
      //dp.push_back(ptr2);
      if(ptrQ2) {
	ptrClu2 = new_ptr(Cluster(ptrQ2,ptr2));
	//dp.push_back(ptrQ2);
	//pstep->addDecayProduct(iCluPtr,ptrQ2);
	ptr2->addChild(ptrClu2);
      } else {
	ptr1 = getParticle(idQ2);
	ptrClu2 = new_ptr(Cluster(ptr1,ptr2));
	///dp.push_back(ptr1);
	pstep->addIntermediate(ptr1);
	ptr1->addChild(ptrClu2);
	ptr2->addChild(ptrClu2);
	pstep->addIntermediate(ptr1);
      }

      ptrClu2->set5Momentum(pClu2);
      ptrClu2->setVertex(positionClu2);
      for(int i = 0; i<ptrClu2->numComponents(); i++) {
	if(ptrClu2->particle(i)->id() == idQ2)
	  ptrClu2->particle(i)->set5Momentum(pQ2);
	else if(ptrClu2->particle(i)->id() == idNew) 
	  ptrClu2->particle(i)->set5Momentum(pQtwo);
      }

      //iCluPtr->addChild(ptrClu2); 
      pstep->addDecayProduct(iCluPtr, ptrClu2);
      clusters.push_back(ptrClu2);
      if(pow(ptrClu2->mass(), _ClPow) > 
	 pow(_ClMax, _ClPow) + pow(ptrClu2->sumConstituentMasses(), _ClPow)) {
	clusterStack.push_back(ptrClu2);
      } 

    }
    //cout << "Sanity check\n";
    // Sanity check (normally skipped) to see if the energy-momentum is conserved.
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::minimal_Hadronization ) {    
      Lorentz5Momentum diff = pClu - ( pClu1 + pClu2 );
      Energy ediff = fabs( diff.m() );
      if ( ediff > 1e-3*GeV ) {
	generator()->logWarning( Exception("ClusterFissioner::cut " 
					   "***Violation of energy-momentum conservation***", 
					   Exception::warning) );    
	if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization ) {    
	  generator()->log() << "         ===> " << endl
			     << " \t Cluster split : mass=" << Mclu
			     << " components ids= " << idQ1 << " " << idQ2 << endl
			     << " \t Clu1 : mass=" << Mclu1 << " " 
			     << " components ids= " << idQ1 << " " << -idNew << endl
			     << " \t Clu2 : mass=" << Mclu2 << " " 
			     << " components ids= " << idQ2 << " " << idNew << endl
			     << "   diff " << diff << endl
			     << "   " << pClu << " ---> " << (pClu1+pClu2) << endl
			     << " = " << pClu1 << " + " << pClu2
			     << endl << endl;
	}      
      }
    }

  }  // end of the main (while) loop over vecCluPtr
  //cout << "done\n";
}


long ClusterFissioner::drawnNewFlavour() const {

  // Flavour is assumed to be only  u, d, s,  with weights
  // (which are not normalized probabilities) given
  // by the same weights as used in HadronsSelector for
  // the decay of clusters into two hadrons. 
  double prob_d = _hadronsSelector->pwtDquark();
  double prob_u = _hadronsSelector->pwtUquark();
  double prob_s = _hadronsSelector->pwtSquark();
  double sum = prob_d + prob_u + prob_s;
  prob_d = prob_d / sum;
  prob_u = prob_u / sum;
  prob_s = prob_s / sum;
  int choice = rnd3(prob_u, prob_d, prob_s);
  long idNew = 0;
  switch (choice) {
  case 0: idNew = Pythia7::ParticleID::u; break;  
  case 1: idNew = Pythia7::ParticleID::d; break;
  case 2: idNew = Pythia7::ParticleID::s; break;
  }
  return idNew;
}


bool ClusterFissioner::drawnChildrenMasses(const Energy Mclu, const Energy m1,    
					   const Energy m2, const Energy m, 
					   Energy & Mclu1, Energy & Mclu2,    // output
					   const double exponent1, const double exponent2, 
					   const Energy average, const int iRemnant) const {

  // This method, given in input the cluster mass Mclu of an heavy cluster C,
  // made of consituents of masses m1 and m2, draws the masses Mclu1 and Mclu2
  // of, respectively, the children cluster C1, made of constituent masses m1 and m, 
  // and cluster C2, of mass Mclu2 and made of constituent masses m2 and m.
  // The mass is extracted from one of the two following mass distributions:
  //   --- power-like ("normal" distribution)
  //                        d(Prob) / d(M^exponent) = const
  //       where the exponent can be different from the two children C1 (exponent1)
  //       and C2 (exponent2).
  //   --- exponential ("soft" distribution) 
  //                        d(Prob) / d(M^2) = exp(-b*M) 
  //       where b = 2.0 / average.
  // Such distributions are limited below by the masses of
  // the constituents quarks, and above from the mass of decaying cluster C. 
  // The choice of which of the two mass distributions to use for each of the
  // two cluster children is dictated by  iRemnant  (see below).
  // If the number of attempts to extract a pair of mass values that are 
  // kinematically acceptable is above some fixed number (max_loop, see below)
  // the method gives up and returns false; otherwise, when it succeeds, it
  // returns true. 
  
  // Initialization for the power-like ("normal") mass distribution
  double Mmin_toExp1 = pow( m1+m , exponent1 );
  double Mmax_toExp1 = pow( Mclu  , exponent1 );
  double Mmin_toExp2 = pow( m2+m , exponent2 );
  double Mmax_toExp2 = pow( Mclu  , exponent2 );

  // Initialization for the exponential ("soft") mass distribution.
  InvEnergy b = 2.0 / average;
  Energy max  = Mclu - m1 - m2 - 2.0*m;
  double rmin = 0.0;
  if ( exp(-b*max) < 50.0 ) { 
    rmin = exp(-b*max);
  } 
  int counter = 0, max_loop = 1000;

  do {    
    // The choice of which of the two mass distributions to use for each of the
    // two cluster children is dictated by  iRemnant, as follows:
    // i)   iRemnant == 0  then both Mclu1 and Mclu2 are extracted by the 
    //                     power-like distribution;
    // ii)  iRemnant == 1  then Mclu1 is extracted from the exponential distribution,
    //                     whereas Mclu2 is extracted from the power-like distribution;
    // iii) iRemnant == 2  then Mclu1 is extracted from the power-like distribution,
    //                     whereas Mclu2 is extracted from the exponential distribution;
    // iv)  otherwise      both Mclu1 and Mclu2 are extracted from the exponential
    //                     distribution,
    if ( iRemnant == 0  ||  iRemnant == 2 ) { 
      Mclu1 = pow( rnd( Mmin_toExp1, Mmax_toExp1) , 1.0/exponent1 );
    } else { 
      double r1 = rnd(rmin, 1.0-rmin) * rnd(rmin, 1.0-rmin);
      if ( r1 > rmin ) {
	Mclu1 = m1 + m - log(r1) / b;
      } else {
	Mclu1 = Energy();
      }
    }
    if ( iRemnant == 0  ||  iRemnant == 1 ) { 
      Mclu2 = pow( rnd( Mmin_toExp2, Mmax_toExp2) , 1.0/exponent2 );
    } else {
      double r2 = rnd(rmin, 1.0-rmin) * rnd(rmin, 1.0-rmin);
      if ( r2 > rmin ) {
	Mclu2 = m + m2 - log(r2) / b;
      } else {
	Mclu2 = Energy();
      }
    }
    counter++;
  } while ( ( Mclu1 < m1+m  ||  Mclu2 < m+m2  ||  Mclu1+Mclu2 > Mclu )  &&  counter < max_loop );
  
  if ( counter < max_loop ) {
    return true;
  } else {
    Mclu1 = Mclu2 = Energy();
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization ) {    
      generator()->logWarning( Exception("ClusterFissioner::drawnChildrenMasses "
					 "***Too many unsuccessful attempts: give up***", 
					 Exception::warning) );
      generator()->log() << "         ===>" 
			 << " Mclu=" << Mclu << " m1=" << m1 << " m2=" << m2 << " m=" << m 
			 << endl << endl;
    }
    return false;
  }
}


void ClusterFissioner::calculateKinematics(const Lorentz5Momentum & pClu, // input
					   const Lorentz5Momentum & p0Q1, // input
					   const bool decayOneHadronClu1, // input
					   const bool decayOneHadronClu2, // input
					   Lorentz5Momentum & pClu1,      // output
					   Lorentz5Momentum & pClu2,      // all 
					   Lorentz5Momentum & pQ1,        // the 
					   Lorentz5Momentum & pQbar,      // rest
					   Lorentz5Momentum & pQ,         
					   Lorentz5Momentum & pQ2bar ) const {

  // This method solves the kinematics of the two body cluster decay:
  //    C (Q1 Q2bar)  --->  C1 (Q1 Qbar)  +  C2 (Q Q2bar)
  // In input we receive the momentum of C, pClu, and the momentum 
  // of the quark Q1 (constituent of C), p0Q1, both in the LAB frame.
  // Furthermore, two boolean variables inform whether the two fission
  // products (C1, C2) decay immediately into a single hadron (in which
  // case the cluster itself is identify with that hadron) and we do 
  // not have to solve the kinematics of the components (Q1,Qbar) for
  // C1 and (Q,Q2bar) for C2.
  // The output is given by the following momenta (all 5-components, 
  // and all in the LAB frame):
  //   pClu1 , pClu2   respectively of   C1 , C2  
  //   pQ1 , pQbar     respectively of   Q1 , Qbar  in  C1
  //   pQ  , pQ2bar    respectively of   Q  , Q2    in  C2
  // The assumption, suggested from the string model, is that, in C frame,
  // C1 and its constituents Q1 and Qbar are collinear, and collinear to
  // the direction of Q1 in C (that is before cluster decay); similarly,
  // (always in the C frame) C2 and its constituents Q and Q2bar are
  // collinear (and therefore anti-collinear with C1,Q1,Qbar).
  // The solution is then obtained by using Lorentz boosts, as follows.
  // The kinematics of C1 and C2 is solved in their parent C frame,
  // and then boosted back in the LAB. The kinematics of Q1 and Qbar
  // is solved in their parent C1 frame and then boosted back in the LAB; 
  // similarly, the kinematics of Q and Q2bar is solved in their parent 
  // C2 frame and then boosted back in the LAB. In each of the three
  // "two-body decay"-like cases, we use the fact that the direction 
  // of the motion of the decay products is known in the rest frame of 
  // their parent. This is obvious for the first case in which the 
  // parent rest frame is C; but it is also true in the other two cases
  // where the rest frames are C1 and C2. This is because C1 and C2 
  // are boosted w.r.t. C in the same direction where their components, 
  // respectively (Q1,Qbar) and (Q,Q2bar) move in C1 and C2 rest frame
  // respectively.      
  // Of course, although the notation used assumed that C = (Q1 Q2bar)
  // where Q1 is a quark and Q2bar an antiquark, indeed everything remain
  // unchanged also in all following cases:
  //  Q1 quark, Q2bar antiquark;           --> Q quark;
  //  Q1 antiquark , Q2bar quark;          --> Q antiquark;  
  //  Q1 quark, Q2bar diquark;             --> Q quark
  //  Q1 antiquark, Q2bar anti-diquark;    --> Q antiquark
  //  Q1 diquark, Q2bar quark              --> Q antiquark
  //  Q1 anti-diquark, Q2bar antiquark;    --> Q quark

  // Calculate the unit three-vector, in the C frame, along which
  // all of the constituents and children clusters move.
  Lorentz5Momentum u(p0Q1);
  u.boost( -pClu.boostVector() );        // boost from LAB to C
  // the unit three-vector is then  u.vect().unit()
						       
  // Calculate the momenta of C1 and C2 in the (parent) C frame first, 
  // where the direction of C1 is u.vect().unit(), and then boost back in the LAB frame.
  Kinematics::twoBodyDecay(pClu, pClu1.mass(), pClu2.mass(), u.vect().unit(), pClu1, pClu2);

  // In the case that cluster1 does not decay immediately into a single hadron,
  // calculate the momenta of Q1 (as constituent of C1) and Qbar in the
  // (parent) C1 frame first, where the direction of Q1 is u.vect().unit(), 
  // and then boost back in the LAB frame. 
  if ( ! decayOneHadronClu1 ) {
    Kinematics::twoBodyDecay(pClu1, pQ1.mass(), pQ.mass(), u.vect().unit(), pQ1, pQbar);
  }

  // In the case that cluster2 does not decay immediately into a single hadron,
  // Calculate the momenta of Q and Q2bar (as constituent of C2) in the
  // (parent) C2 frame first, where the direction of Q is u.vect().unit(), 
  // and then boost back in the LAB frame. 
  if ( ! decayOneHadronClu2 ) {
    Kinematics::twoBodyDecay(pClu2, pQ.mass(), pQ2bar.mass(), u.vect().unit(), pQ, pQ2bar);
  }
}


void ClusterFissioner::calculatePositions( const Lorentz5Momentum & pClu,    // input
					   const LorentzPoint & positionClu, // input
					   const Lorentz5Momentum & pClu1,   // input
					   const Lorentz5Momentum & pClu2,   // input
					   LorentzPoint & positionClu1,      // output
					   LorentzPoint & positionClu2 ) const { // output

  // Determine positions of cluster children.
  // See Marc Smith's thesis, page 127, formulas (4.122) and (4.123).

  Energy Mclu  = pClu.m();
  Energy Mclu1 = pClu1.m();
  Energy Mclu2 = pClu2.m();

  // Calculate the unit three-vector, in the C frame, along which 
  // children clusters move.
  Lorentz5Momentum u(pClu1);
  u.boost( -pClu.boostVector() );        // boost from LAB to C frame
  // the unit three-vector is then  u.vect().unit()

  Energy pstarChild = Kinematics::pstarTwoBodyDecay(Mclu,Mclu1,Mclu2);

  // First, determine the relative positions of the children clusters
  // in the parent cluster reference frame.
  double GeV2mm = _globalParameters->conversionFactorGeVtoMillimeter();
  Length x1 = GeV2mm * (Mclu*0.25 + 0.5
                     *(pstarChild + (sqr(Mclu2) - sqr(Mclu1))/(2.0*Mclu)))/GeV;
  Length t1 = ((Mclu/GeV) * GeV2mm - x1); 
  LorentzDistance distanceClu1( x1 * u.vect().unit(), t1 );
  Length x2 = GeV2mm * (-Mclu*0.25 + 0.5
                     *(-pstarChild + (sqr(Mclu2) - sqr(Mclu))/(2.0*Mclu)))/GeV;
  Length t2 = ((Mclu/GeV) * GeV2mm + x2);
  LorentzDistance distanceClu2(x2 * u.vect().unit(), t2);

  // Debugging
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Hadronization ) {
    generator()->log() << "ClusterFissioner::calculatePositions : *** extreme debugging ***" << endl
                       << "\t distanceClu1 = " << distanceClu1 
                       << "\t invariant length = " << distanceClu1.mag() << "  [mm] " << endl
                       << "\t distanceClu2 = " << distanceClu2 
                       << "\t invariant length = " << distanceClu2.mag() << "  [mm] " << endl;
  }

  // Then, transform such relative positions from the parent cluster
  // reference frame to the Lab frame.
  distanceClu1.boost( pClu.boostVector() );
  distanceClu2.boost( pClu.boostVector() );

  // Finally, determine the absolute positions in the Lab frame.
  positionClu1 = positionClu + distanceClu1;
  positionClu2 = positionClu + distanceClu2;

}


    
