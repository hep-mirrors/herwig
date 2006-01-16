// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LightClusterDecayer class.
//

#include "LightClusterDecayer.h"
#include <ThePEG/Interface/ClassDocumentation.h>
#include <ThePEG/Interface/Parameter.h>
#include <ThePEG/Interface/Reference.h>
#include <ThePEG/Persistency/PersistentOStream.h>
#include <ThePEG/Persistency/PersistentIStream.h>
#include <ThePEG/PDT/EnumParticles.h>
#include <ThePEG/Repository/EventGenerator.h>
#include "Cluster.h"
#include "HadronSelector.h"
#include "Herwig++/Utilities/HwDebug.h"
#include "Herwig++/Utilities/Kinematics.h"
#include "Herwig++/Utilities/CheckId.h"


using namespace Herwig;
// using namespace ThePEG;


LightClusterDecayer::~LightClusterDecayer() {}


void LightClusterDecayer::persistentOutput(PersistentOStream & os) const 
{os << _hadronsSelector << _B1Lim;}

void LightClusterDecayer::persistentInput(PersistentIStream & is, int) 
{is >> _hadronsSelector >> _B1Lim;}


ClassDescription<LightClusterDecayer> LightClusterDecayer::initLightClusterDecayer;
// Definition of the static class description member.


void LightClusterDecayer::Init() {

  static ClassDocumentation<LightClusterDecayer> documentation
    ("There is the class responsible for the one-hadron decay of light clusters");

  static Reference<LightClusterDecayer,HadronSelector> 
    interfaceHadronSelector("HadronSelector", 
			     "A reference to the HadronSelector object", 
			     &Herwig::LightClusterDecayer::_hadronsSelector,
			     false, false, true, false);
  
  static Parameter<LightClusterDecayer,Energy>
    interfaceB1Lim ("B1Lim","one-hadron decay of b-cluster over threshold",
                    &LightClusterDecayer::_B1Lim, 0, 0.0, 0.0, 100.0,false,false,false);

}


bool LightClusterDecayer::decay(const StepPtr &pstep) {
  //  throw (Veto, Stop, Exception) {

  // Loop over all clusters, and for those that were not heavy enough
  // to undergo to fission, check if they are below the threshold
  // for normal two-hadron decays. If this is the case, then the cluster
  // should be decayed into a single hadron: this can happen only if
  // it is possible to reshuffle momenta between the cluster and
  // another one; in the rare occasions in which such exchange of momenta
  // is not possible (because all of the clusters are too light) then 
  // the event is skipped. 
  // Notice that, differently from what happens in Fortran Herwig, 
  // light (that is below the threshold for the production of the lightest
  // pair of hadrons with the proper flavours) fission products, produced 
  // by the fission of heavy clusters in class  ClusterFissioner
  // have been already "decayed" into single hadron (the lightest one
  // with proper flavour) by the same latter class, without requiring
  // any reshuffling. Therefore the light clusters that are treated in
  // this  LightClusterDecayer  class are produced directly
  // (originally) by the class  ClusterFinder. 
    
  // To preserve all of the information, the cluster partner with which 
  // the light cluster (that decays into a single hadron) exchanges
  // momentum in the reshuffling procedure is redefined and inserted 
  // in the vector  vecNewRedefinedCluPtr. Only at the end, when all
  // light clusters have been examined, the elements this vector will be 
  // copied in  collecCluPtr  (the reason is that it is not allowed to 
  // modify a STL container while iterating over it. At the same time,
  // this ensures that a cluster can be redefined only once, which seems
  // sensible although not strictly necessary). 
  // Notice that the cluster reshuffling partner is normally redefined
  // and inserted in the vector  vecNewRedefinedCluPtr, but not always:
  // in the case it is also light, then it is also decayed immediately
  // into a single hadron, without redefining it (the reason being that,
  // otherwise, the would-be redefined cluster could have undefined
  // components). 
  vector<tClusterPtr> redefinedClusters;
  ClusterVector clusters; 
  for (ParticleSet::iterator it = pstep->particles().begin();
       it!= pstep->particles().end(); it++) { 
    if((*it)->id() == ExtraParticleID::Cluster) 
      clusters.push_back(dynamic_ptr_cast<ClusterPtr>(*it));
  }

  for (ClusterVector::iterator it = clusters.begin();
       it != clusters.end(); ++it) {
        
    // Skip the clusters that are not available or that are 
    // heavy, intermediate, clusters that have undergone to fission,
    if ( ! (*it)->isAvailable()  ||  ! (*it)->isReadyToDecay() ) continue; 
                                     
    // We need to require (at least at the moment, maybe in the future we 
    // could change it) that the cluster has exactly two components, 
    // because otherwise we don't know how to deal with the kinematics.
    // If this is not the case, then send a warning because it is not suppose 
    // to happen, and then do nothing with (ignore) such cluster.
    if ( (*it)->numComponents() != 2 ) {
      generator()->logWarning( Exception("LightClusterDecayer::decay "
					 "***Still cluster with not exactly 2 components*** ", 
					 Exception::warning) );
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization ) {    
	generator()->log() << "         ===>" << " num components = " << (*it)->numComponents()
			   << endl << endl;
      }
      continue;
    }
    
    // Extract the id and particle pointer of the two components of the cluster.
    long idQ1 = 0, idQ2 = 0;
    tPPtr ptrQ1 = tPPtr(), ptrQ2 = tPPtr();
    ptrQ1 = (*it)->particle(0);
    idQ1 = ptrQ1->id();
    ptrQ2 = (*it)->particle(1);
    idQ2 = ptrQ2->id();

    // Sanity check (normally skipped) to control that the two components of a
    // cluster are consistent, that is they can form a meson or a baryon.
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::minimal_Hadronization ) {
      if ( ! CheckId::canBeMeson(idQ1,idQ2)  &&  ! CheckId::canBeBaryon(idQ1,idQ2) ) {
	generator()->logWarning( Exception("LightClusterDecayer::decay "
	      "***The two components of the cluster are inconsistent***", 
					   Exception::warning) );
	std::cout << "LightClusterDecayer::decay ****the two components are "
		  << idQ1 << " and " << idQ2 << " ***\n";
	if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization ) {    
	  generator()->log() << "         ===>" 
			     << " idQ1=" << idQ1 << " idQ2=" << idQ2 << endl << endl;
	}
      }
    }
      
    // Determine the sum of the nominal masses of the two lightest hadrons
    // with the right flavour numbers as the cluster under consideration.
    // Notice that we don't need real masses (drawn by a Breit-Wigner 
    // distribution) because the lightest pair of hadrons does not involve
    // any broad resonance.
    Energy threshold = _hadronsSelector->massLightestHadronPair(idQ1,idQ2);
    // Special for b-flavour: it allows one-hadron decays also above threshold.
    if ( CheckId::hasBeauty(idQ1,idQ2) ) {threshold *= (1.0 + rnd()*_B1Lim);} 
    
    //***TRICK***: scale artificially threshold if you want to test
    //             LightClusterDecayer with a huge statistics.
    // // // threshold *= (1.0 + rnd()*2.0);
    
    if ( (*it)->mass() < threshold ) {
      
      long idhad = _hadronsSelector->lightestHadron(idQ1,idQ2);
      
      // We assume that the candidate reshuffling cluster partner, 
      // with whom the light cluster can exchange momenta,
      // is chosen as the closest in space-time between the available
      // clusters. Notice that an alternative, sensible approach
      // could be to consider instead the "closeness" in the colour
      // structure... 
      // Notice that nor a light cluster (which decays into a single hadron) 
      // neither its cluster reshuffling partner (which either has a 
      // redefined cluster or also decays into a single hadron) can be
      // a reshuffling partner of another light cluster.
      // This because we are requiring that the considered candidate cluster
      // reshuffling partner has the status  "isAvailable && isReadyToDecay"  true; 
      // furthermore, the new redefined clusters are not added to the collection 
      // of cluster before the end of the entire reshuffling procedure, avoiding
      // in this way that the redefined cluster of a cluster reshuffling partner 
      // is used again later. Needless to say, this is just an assumption,
      // although reasonable, but nothing more than that!
      
      // Build a multimap of available reshuffling cluster partners,
      // with key given by the module of the invariant space-time distance 
      // w.r.t. the light cluster, so that this new collection is automatically 
      // ordered in increasing distance values. 
      // We use a multimap, rather than a map, just for precaution against not properly 
      // defined cluster positions which could produce all identical (null) distances.
      multimap<Length,tClusterPtr> candidates;
      for ( ClusterVector::iterator jt = clusters.begin();
	    jt != clusters.end(); ++jt ) {
	if ((*jt)->isAvailable() && (*jt)->isReadyToDecay() && jt != it) {
	  Length distance = fabs (((*it)->vertex() - (*jt)->vertex()).mag());
	  candidates.insert(pair<Length,tClusterPtr>(distance,*jt)); 
	}
      }

      // Loop sequentially the multimap.
      multimap<Length,tClusterPtr>::const_iterator mmapIt = candidates.begin();
      bool found = false;
      while (!found && mmapIt  != candidates.end()) {
	found = reshuffling(idhad, *it, (*mmapIt).second, pstep, redefinedClusters);
	if (!found) ++mmapIt;
      }
      
      if (!found) {
	if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization ) {    
	  int numCluAvailable = 0, numCluAlreadyReshuffled = 0;
	  for (ClusterVector::const_iterator iter = clusters.begin();
		iter != clusters.end(); ++iter ) {
	    if ( (*iter)->isAvailable() && (*iter)->isReadyToDecay() ) {
	      numCluAvailable++;
	    } else if((*iter)->isAvailable() && (*iter)->hasBeenReshuffled()) {
	      numCluAlreadyReshuffled++;
	    }
	  }
	  generator()->log() << "LightClusterDecayer::decay "
			     << " (just before throwing the exception) " << endl 
			     << "         ===> num clusters:  all=" << clusters.size()
			     << "  already reshuffled=" << numCluAlreadyReshuffled
			     << "  still available=" << numCluAvailable 
			     << endl << endl;
	}
// 	throw Exception("LightClusterDecayer::decay "
// 			"***Skip event: too light clusters (reshuffling problem)***",
// 			Exception::eventerror);
	// special for partonic b/c decays
	return partonicReshuffle(idhad,*it,pstep);
      } 
      
      // Debugging.
      if ( HERWIG_DEBUG_LEVEL == HwDebug::extreme_Hadronization ) {    
	generator()->log() << "LightClusterDecayer::decay : *** extreme debugging ***" << endl
			   << " \t Cluster : mass=" << (*it)->mass()
			   << " components ids= " << idQ1 << " " << idQ2
			   << " ---> hadron id=" << idhad << endl 
			   << " before " << ( (*it)->momentum() + 
					      ((*mmapIt).second)->momentum() )  
			   << endl 
			   << " = " << (*it)->momentum() 
			   << " + " << ((*mmapIt).second)->momentum() << endl
	                   << "\t multimap of reshuffling candidates : size = " 
                           << candidates.size() << endl;
	for(multimap<Length,tClusterPtr>::const_iterator
	  mmapIt = candidates.begin(); mmapIt != candidates.end(); ++mmapIt ) {
	  generator()->log() << "\t \t distance = " << (*mmapIt).first << "  [mm] " << endl;
	}
      }
	
    } // end if mass < threshold
  } // end loop over collecCluPtr 

  // Add to  collecCluPtr  all of the redefined new clusters (indeed the 
  // pointers to them are added) contained in  vecNewRedefinedCluPtr.
  for (tClusterVector::const_iterator it = redefinedClusters.begin();
       it != redefinedClusters.end(); ++it) {
    clusters.push_back(*it);
  }
  return true;
}


bool LightClusterDecayer::reshuffling(const long idhad1, 
				      tClusterPtr cluPtr1, 
				      tClusterPtr cluPtr2,
				      const StepPtr pstep,
				      tClusterVector & redefinedClusters )
  throw (Veto, Stop, Exception) {
  // don't reshuffle with beam clusters
  if(cluPtr2->isBeamCluster()) return false;

  // This method does the reshuffling of momenta between the cluster "1", 
  // that must decay into a single hadron (with id equal to idhad1), and 
  // the candidate cluster "2". It returns true if the reshuffling succeed,
  // false otherwise.

  PPtr ptrhad1 = getParticle( idhad1 );
  if ( ! ptrhad1 ) {
    generator()->logWarning( Exception("LightClusterDecayer::reshuffling"
				       "***Cannot create a particle with specified id***", 
				       Exception::warning) );
    return false;
  }
  Energy mhad1 = ptrhad1->mass();

  // Let's call "3" and "4" the two constituents of the second cluster
  tPPtr part3 = cluPtr2->particle(0);
  tPPtr part4 = cluPtr2->particle(1);

  // Sanity check (normally skipped) to see if the second cluster
  // is consistently defined.
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::minimal_Hadronization ) {
    if ( cluPtr2->numComponents() != 2  ||  !part3  ||  !part4 ) { 
      generator()->logWarning( Exception("LightClusterDecayer::reshuffling "
					 "***Inconsistent cluster***", 
					 Exception::warning) );
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization ) {    
	generator()->log() << "         ===>" << " n=" << cluPtr2->numComponents()
			   << " compPtr3=" << part3 << " compPtr4=" << part4 
			   << endl << endl;
      }
    }
  }
  
  // Check if the system of the two clusters can kinematically be replaced by
  // an hadron of mass mhad1 (which is the lightest single hadron with the 
  // same flavour numbers as the first cluster) and the second cluster.
  // If not, then try to replace the second cluster with the lightest hadron
  // with the same flavour numbers; if it still fails, then give up!
  Lorentz5Momentum pSystem = cluPtr1->momentum() + cluPtr2->momentum();
  pSystem.rescaleMass();  // set the mass as the invariant of the quadri-vector
  Energy mSystem = pSystem.mass();
  Energy mclu2 = cluPtr2->mass();
  //***TRICK***: uncomment the following line and replace it to the 
  //             if statement below if you want to test  LightClusterDecayer
  //             with a huge statistics.
  // // // if ( mSystem  <  (mhad1 + mclu2)*(1.0 + rnd()*5.0) ) {
  bool singleHadron = false;
  Energy mLHP2 = _hadronsSelector->massLightestHadronPair(part3->id(), part4->id());
  Energy mLH2 = _hadronsSelector->massLightestHadron(part3->id(), part4->id());

  if(mSystem > mhad1 + mclu2 && mclu2 > mLHP2) { singleHadron = false; } 
  else if(mSystem > mhad1 + mLH2) { singleHadron = true; mclu2 = mLH2; }
  else return false;
  // Let's call from now on "Sys" the system of the two clusters, and
  // had1 (of mass mhad1) the lightest hadron in which the first
  // cluster decays, and clu2 (of mass mclu2) either the second
  // cluster or the lightest hadron in which it decays (depending
  // which one is kinematically allowed, see above).
  // The idea behind the reshuffling is to replace the system of the
  // two clusters by the system of the hadron had1 and (cluster or hadron) clu2,
  // but leaving the overall system unchanged. Furthermore, the motion
  // of had1 and clu2 in the Sys frame is assumed to be parallel to, respectively,
  // those of the original cluster1 and cluster2 in the same Sys frame.

  // Calculate the unit three-vector, in the frame "Sys" along which the 
  // two initial clusters move.
  Lorentz5Momentum u( cluPtr1->momentum() );
  u.boost( - pSystem.boostVector() );         // boost from LAB to Sys
							 
  // Calculate the momenta of had1 and clu2 in the Sys frame first, 
  // and then boost back in the LAB frame.
  Lorentz5Momentum phad1, pclu2;
  Kinematics::twoBodyDecay(pSystem, mhad1, mclu2, u.vect().unit(), phad1, pclu2);

  // Sanity check (normally skipped) to see if the energy-momentum is conserved.
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::minimal_Hadronization ) {    
    Lorentz5Momentum diff = pSystem - ( phad1 + pclu2 );
    Energy ediff = fabs( diff.m() );
    if ( ediff > 1e-3*GeV ) {
      generator()->logWarning( Exception("LightClusterDecayer::reshuffling " 
					 "***Energy-momentum NOT conserved: system***", 
					 Exception::warning) );    
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization ) {    
	generator()->log() << "         ===> " << endl
			   << "   diff " << diff << endl
			   << " before " << pSystem << " ---> " << (phad1+pclu2) << endl
			   << "  = " << phad1 << " + " << pclu2 << endl << endl;
      }      
    }
  }

  ptrhad1->set5Momentum( phad1 );               // set momentum of first hadron.
  ptrhad1->setLabVertex(cluPtr1->vertex()); // set hadron vertex position to the
                                                // parent cluster position.
  //cluPtr1->addChild(ptrhad1);
  pstep->addDecayProduct(cluPtr1, ptrhad1);
  cluPtr1->reshufflingPartnerCluster( cluPtr2 );
  cluPtr2->reshufflingPartnerCluster( cluPtr1 );

  if(singleHadron) {  

    // In the case that also the cluster reshuffling partner is light
    // it is decayed into a single hadron, *without* creating the
    // redefined cluster (this choice is justified in order to avoid
    // clusters that could have undefined components).

    long idhad2 = _hadronsSelector->lightestHadron(part3->id(), part4->id());
    PPtr ptrhad2 = getParticle( idhad2 );
    ptrhad2->set5Momentum( pclu2 );            
    ptrhad2->setLabVertex( cluPtr2->vertex() ); // set hadron vertex position to the
                                                  // parent cluster position.
    //    cluPtr2->addChild(ptrhad2);   
    pstep->addDecayProduct(cluPtr2, ptrhad2);
  } else {

    // Create the new cluster which is the redefinitions of the cluster  
    // partner (cluster "2") used in the reshuffling procedure of the
    // light cluster (cluster "1").
    // The rationale of this is to preserve completely all of the information.
    ClusterPtr cluPtr2new = ClusterPtr();
    if(part3 && part4) cluPtr2new = new_ptr(Cluster(part3,part4));

    cluPtr2new->set5Momentum( pclu2 );
    cluPtr2new->setVertex( cluPtr2->vertex() );
    //    cluPtr2->addChild( cluPtr2new );
    pstep->addDecayProduct(cluPtr2, cluPtr2new);
    redefinedClusters.push_back( cluPtr2new );
       
    // Set consistently the momenta of the two components of the second cluster
    // after the reshuffling. To do that we first calculate the momenta of the 
    // constituents in the initial cluster rest frame; then we boost them back 
    // in the lab but using this time the new cluster rest frame. Finally we store 
    // these information in the new cluster. Notice that we do *not* set 
    // consistently also the momenta of the (eventual) particles pointed by the 
    // two components: that's because we do not need to do so, being the momentum
    // an explicit private member of the class Component (which is set equal
    // to the momentum of the eventual particle pointed only in the constructor,
    // but then later should not necessary be the same), and furthermore it allows
    // us not to loose any information, in the sense that we can always, later on,
    // to find the original momenta of the two components before the reshuffling.
    Lorentz5Momentum p3 = part3->momentum(); //p3new->momentum();
    p3.boost( - (cluPtr2->momentum()).boostVector() );  // from LAB to clu2 (old) frame
    p3.boost( pclu2.boostVector() );    // from clu2 (new) to LAB frame
    Lorentz5Momentum p4 = part4->momentum(); //p4new->momentum();
    p4.boost( - (cluPtr2->momentum()).boostVector() );  // from LAB to clu2 (old) frame 
    p4.boost( pclu2.boostVector() );    // from clu2 (new) to LAB frame    
    cluPtr2new->particle(0)->set5Momentum(p3);
    cluPtr2new->particle(1)->set5Momentum(p4);

    // Sanity check (normally skipped) to see if the energy-momentum is conserved.
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::minimal_Hadronization ) {    
      Lorentz5Momentum sumMomentaComponents = Lorentz5Momentum(); 
      for(int i = 0; i<cluPtr2new->numComponents(); i++)
	sumMomentaComponents += cluPtr2new->particle(i)->momentum();

      Lorentz5Momentum diff = sumMomentaComponents - cluPtr2new->momentum(); 
      Energy ediff = fabs( diff.m() );
      if ( ediff > 1e-3*GeV ) {
	generator()->logWarning( Exception("LightClusterDecayer::reshuffling " 
					   "***Energy-momentum NOT conserved: clu2new***", 
					   Exception::warning) );    
	if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Hadronization ) {    
	  generator()->log() << "         ===> " << endl
			     << "   sumMomentaComponents = " << sumMomentaComponents << endl 
			     << "           momentumClu2 = " << cluPtr2new->momentum() << endl
			     << "                   diff = " << diff << endl << endl;
	}      
      }
    }


  } // end of  if (singleHadron) 


  return true;
}

bool LightClusterDecayer::partonicReshuffle(const long idhad,const PPtr cluster,
					    const StepPtr pstep)
{
  tPPtr meson(cluster);
  if(!meson->parents().empty()) meson=meson->parents()[0];
  if(!meson->parents().empty()) meson=meson->parents()[0];
  // check b/c hadron decay
  int ptype(abs(meson->id()/1000));
  if(ptype>5) ptype/=10;
  if(ptype!=4&&ptype!=5) return false;
  // get the leptons
  PPtr leptons[2];
  unsigned int nlep(0);
  Lorentz5Momentum pleptons;
  Energy mLHP2=Energy();
  for(unsigned int ix=0;ix<meson->children().size();++ix)
    {
      if(!(meson->children()[ix]->dataPtr()->coloured()))
	{
	  if(nlep<=2) leptons[nlep]=meson->children()[ix];
	  ++nlep;
	  pleptons+=meson->children()[ix]->momentum();
	  mLHP2+=meson->children()[ix]->mass();
	}
      pleptons.rescaleMass();
    }
  // check leptons
  if(nlep==1)
    {
      mLHP2=Energy();
      for(unsigned int ix=0;ix<leptons[0]->children().size();++ix)
	{mLHP2+=leptons[0]->children()[ix]->mass();}
    }
  else if (nlep>2) return false;
  // check we can do the reshuffling
  PPtr ptrhad = getParticle(idhad);
  Energy mhad1(ptrhad->mass());
  Lorentz5Momentum pSystem = pleptons + cluster->momentum();
  pSystem.rescaleMass();
  Energy mSystem = pSystem.mass();
  Energy mclu2 = pleptons.mass();
  // check if we can reshuffle
  if(!(mSystem > mhad1+mclu2 && mclu2>mLHP2)) return false;
  // Calculate the unit three-vector, in the frame "Sys" along which the 
  // two initial clusters move.
  Lorentz5Momentum u(cluster->momentum());
  u.boost( - pSystem.boostVector() );         // boost from LAB to Sys
  // Calculate the momenta of had1 and clu2 in the Sys frame first, 
  // and then boost back in the LAB frame.
  Lorentz5Momentum phad1, pclu2;
  Kinematics::twoBodyDecay(pSystem, mhad1, mclu2, u.vect().unit(), phad1, pclu2);
  ptrhad->set5Momentum( phad1 );         // set momentum of first hadron.
  ptrhad->setLabVertex(cluster->vertex()); // set hadron vertex position to the
  // parent cluster position.
  pstep->addDecayProduct(cluster, ptrhad);
  // reshuffle the leptons
  // boost the leptons to the rest frame of the system
  Hep3Vector boost1(-pleptons.boostVector());
  Hep3Vector boost2( pclu2.boostVector());
  for(unsigned int ix=0;ix<nlep;++ix)
    {
      leptons[ix]->deepBoost(boost1);
      leptons[ix]->deepBoost(boost2);
    }
  return true;
}
