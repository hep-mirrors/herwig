// -*- C++ -*-
//
// LightClusterDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LightClusterDecayer class.
//

#include "LightClusterDecayer.h"
#include <ThePEG/Interface/ClassDocumentation.h>
#include <ThePEG/Interface/Parameter.h>
#include <ThePEG/Interface/Reference.h>
#include <ThePEG/Interface/Switch.h>
#include <ThePEG/Persistency/PersistentOStream.h>
#include <ThePEG/Persistency/PersistentIStream.h>
#include <ThePEG/PDT/EnumParticles.h>
#include <ThePEG/Repository/EventGenerator.h>
#include "Cluster.h"
#include "CheckId.h"
#include "Herwig/Utilities/Kinematics.h"
#include <ThePEG/Utilities/DescribeClass.h>

using namespace Herwig;

DescribeClass<LightClusterDecayer,Interfaced>
describeLightClusterDecayer("Herwig::LightClusterDecayer","");

IBPtr LightClusterDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr LightClusterDecayer::fullclone() const {
  return new_ptr(*this);
}


void LightClusterDecayer::persistentOutput(PersistentOStream & os) const {
  os << _hadronSelector;
} 

void LightClusterDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _hadronSelector;
}

void LightClusterDecayer::Init() {

  static ClassDocumentation<LightClusterDecayer> documentation
    ("There is the class responsible for the one-hadron decay of light clusters");

  static Reference<LightClusterDecayer,HadronSelector> 
    interfaceHadronSelector("HadronSelector", 
			     "A reference to the HadronSelector object", 
			     &Herwig::LightClusterDecayer::_hadronSelector,
			     false, false, true, false);

}

bool LightClusterDecayer::decay(ClusterVector & clusters, tPVector & finalhadrons) {

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


  for (ClusterVector::const_iterator it = clusters.begin();
       it != clusters.end(); ++it) {
    
    // Skip the clusters that are not available or that are 
    // heavy, intermediate, clusters that have undergone to fission,
    if ( ! (*it)->isAvailable()  ||  ! (*it)->isReadyToDecay() ){
     continue; 
    }                                 
    // We need to require (at least at the moment, maybe in the future we 
    // could change it) that the cluster has exactly two components, 
    // because otherwise we don't know how to deal with the kinematics.
    // If this is not the case, then send a warning because it is not suppose 
    // to happen, and then do nothing with (ignore) such cluster.
    if ( (*it)->numComponents() != 2 ) {
      generator()->logWarning( Exception("LightClusterDecayer::decay "
					 "***Still cluster with not exactly"
					 " 2 components*** ", 
					 Exception::warning) );
      continue;
    }
    
    // select the hadron for single hadron decay
    tcPDPtr hadron = _hadronSelector->chooseSingleHadron((*it)->particle(0)->dataPtr(),
							 (*it)->particle(1)->dataPtr(),
							 (**it).mass());
    // if not single decay continue
    if(!hadron){ 
continue;    
    }
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
	Length distance = abs (((*it)->vertex() - (*jt)->vertex()).m());
	candidates.insert(pair<Length,tClusterPtr>(distance,*jt)); 
      }
    }
    
    // Loop sequentially the multimap.
    multimap<Length,tClusterPtr>::const_iterator mmapIt = candidates.begin();
    bool found = false;
    while (!found && mmapIt  != candidates.end()) {
      found = reshuffling(hadron, *it, (*mmapIt).second, redefinedClusters, finalhadrons);
      if (!found) ++mmapIt;
    }
    
    if (!found) return partonicReshuffle(hadron,*it,finalhadrons);
  } // end loop over collecCluPtr 

  // Add to  collecCluPtr  all of the redefined new clusters (indeed the 
  // pointers to them are added) contained in  vecNewRedefinedCluPtr.
  for (tClusterVector::const_iterator it = redefinedClusters.begin();
       it != redefinedClusters.end(); ++it) {
    clusters.push_back(*it);
  }
  return true;
 }


bool LightClusterDecayer::reshuffling(const tcPDPtr pdata1, 
				      tClusterPtr cluPtr1, 
				      tClusterPtr cluPtr2,
				      tClusterVector & redefinedClusters,
				      tPVector & finalhadrons)
  {
  // don't reshuffle with beam clusters
  if(cluPtr2->isBeamCluster()) return false;

  // This method does the reshuffling of momenta between the cluster "1", 
  // that must decay into a single hadron (with id equal to idhad1), and 
  // the candidate cluster "2". It returns true if the reshuffling succeed,
  // false otherwise.

  PPtr ptrhad1 = pdata1->produceParticle();
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
  
  // Check if the system of the two clusters can kinematically be replaced by
  // an hadron of mass mhad1 (which is the lightest single hadron with the 
  // same flavour numbers as the first cluster) and the second cluster.
  // If not, then try to replace the second cluster with the lightest hadron
  // with the same flavour numbers; if it still fails, then give up!
  Lorentz5Momentum pSystem = cluPtr1->momentum() + cluPtr2->momentum();
  pSystem.rescaleMass();  // set the mass as the invariant of the quadri-vector
  Energy mSystem = pSystem.mass();
  Energy mclu2 = cluPtr2->mass();
  bool singleHadron = false;
  Energy mLHP2 = _hadronSelector->massLightestHadronPair(part3->dataPtr(),part4->dataPtr());
  Energy mLH2 = _hadronSelector->massLightestHadron(part3->dataPtr(),part4->dataPtr());

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

  if (pSystem.m() < mhad1 + mclu2 ) {
    throw Exception() << "Impossible Kinematics in LightClusterDecayer::reshuffling()" 
		      << Exception::eventerror;
  }
  
  Kinematics::twoBodyDecay(pSystem, mhad1, mclu2, u.vect().unit(), phad1, pclu2);

  ptrhad1->set5Momentum( phad1 );        // set momentum of first hadron.
  ptrhad1->setVertex(cluPtr1->vertex()); // set hadron vertex position to the
                                         // parent cluster position.
  cluPtr1->addChild(ptrhad1);
  finalhadrons.push_back(ptrhad1);
  cluPtr1->flagAsReshuffled();
  cluPtr2->flagAsReshuffled();


  if(singleHadron) {  

    // In the case that also the cluster reshuffling partner is light
    // it is decayed into a single hadron, *without* creating the
    // redefined cluster (this choice is justified in order to avoid
    // clusters that could have undefined components).

    PPtr ptrhad2 = _hadronSelector->lightestHadron(part3->dataPtr(),part4->dataPtr())
      ->produceParticle();
    ptrhad2->set5Momentum( pclu2 );            
    ptrhad2->setVertex( cluPtr2->vertex() ); // set hadron vertex position to the
                                             // parent cluster position.
    cluPtr2->addChild(ptrhad2);  
    finalhadrons.push_back(ptrhad2);
  } else {

    // Create the new cluster which is the redefinitions of the cluster  
    // partner (cluster "2") used in the reshuffling procedure of the
    // light cluster (cluster "1").
    // The rationale of this is to preserve completely all of the information.
    ClusterPtr cluPtr2new = ClusterPtr();
    if(part3 && part4) cluPtr2new = new_ptr(Cluster(part3,part4));

    cluPtr2new->set5Momentum( pclu2 );
    cluPtr2new->setVertex( cluPtr2->vertex() );
    cluPtr2->addChild( cluPtr2new );
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
  } // end of  if (singleHadron) 
  return true;
}
  
bool LightClusterDecayer::partonicReshuffle(const tcPDPtr had,
					    const PPtr cluster,
					    tPVector & finalhadrons) {
  tPPtr meson(cluster);
  if(!meson->parents().empty()) meson=meson->parents()[0];
  if(!meson->parents().empty()) meson=meson->parents()[0];
  // check b/c hadron decay
  int ptype(abs(meson->id())%10000);
  bool heavy = (ptype/1000 == 5 || ptype/1000 ==4 );
  heavy |= (ptype/100  == 5 || ptype/100  ==4 );
  heavy |= (ptype/10   == 5 || ptype/10   ==4 );
  if(!heavy) return false;
  // find the leptons
  tPVector leptons;
  for(unsigned int ix=0;ix<meson->children().size();++ix) {
    if(!(meson->children()[ix]->dataPtr()->coloured())) {
      leptons.push_back(meson->children()[ix]);
    }
  }
  if(leptons.size()==1) {
    tPPtr w=leptons[0];
    leptons.pop_back();
    for(unsigned int ix=0;ix<w->children().size();++ix) {
      if(!w->children()[ix]->dataPtr()->coloured()) {
	leptons.push_back(w->children()[ix]);
      }
    }
  }
  if(leptons.size()!=2) return false;
  // get momentum of leptonic system and the its minimum possible mass
  Energy mmin(ZERO);
  Lorentz5Momentum pw;
  for(unsigned int ix=0;ix<leptons.size();++ix) {
    pw+=leptons[ix]->momentum();
    mmin+=leptons[ix]->mass();
  }
  pw.rescaleMass();
  // check we can do the reshuffling
  PPtr ptrhad = had->produceParticle();
  // total momentum fo the system
  Lorentz5Momentum pSystem = pw + cluster->momentum();
  pSystem.rescaleMass();
  // normal case get additional energy by rescaling momentum in rest frame of
  // system
  if(pSystem.mass()>ptrhad->mass()+pw.mass()&&pw.mass()>mmin) {
    // Calculate the unit three-vector, in the frame "Sys" along which the 
    // two initial clusters move.
    Lorentz5Momentum u(cluster->momentum());
    u.boost( - pSystem.boostVector() );
    // Calculate the momenta of had1 and clu2 in the Sys frame first, 
    // and then boost back in the LAB frame.
    Lorentz5Momentum phad1, pclu2;
    Kinematics::twoBodyDecay(pSystem, ptrhad->mass(), pw.mass(), 
			     u.vect().unit(), phad1, pclu2);
    // set momentum of first hadron.
    ptrhad->set5Momentum( phad1 );
    // set hadron vertex position to the parent cluster position.    
    ptrhad->setLabVertex(cluster->vertex());
    // add hadron
    cluster->addChild(ptrhad);
    finalhadrons.push_back(ptrhad);
    // reshuffle the leptons
    // boost the leptons to the rest frame of the system
    Boost boost1(-pw.boostVector());
    Boost boost2( pclu2.boostVector());
    for(unsigned int ix=0;ix<leptons.size();++ix) {
      leptons[ix]->deepBoost(boost1);
      leptons[ix]->deepBoost(boost2);
    }
    return true;
  }
  else {
    return false;
  }
}
