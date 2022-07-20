// -*- C++ -*-
//
// ColourReconnector.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ColourReconnector class.
//

#include "ColourReconnector.h"
#include "Cluster.h"

#include <ThePEG/Utilities/DescribeClass.h>
#include <ThePEG/Repository/UseRandom.h>
#include <ThePEG/PDT/StandardMatchers.h>
#include <ThePEG/Persistency/PersistentOStream.h>
#include <ThePEG/Persistency/PersistentIStream.h>

#include <ThePEG/Interface/Switch.h>
#include <ThePEG/Interface/Parameter.h>

#include "Herwig/Utilities/Maths.h"

using namespace Herwig;

using CluVecIt = ColourReconnector::CluVecIt;
using Constants::pi;
using Constants::twopi;


DescribeClass<ColourReconnector,Interfaced>
describeColourReconnector("Herwig::ColourReconnector","Herwig.so");

IBPtr ColourReconnector::clone() const {
  return new_ptr(*this);
}

IBPtr ColourReconnector::fullclone() const {
  return new_ptr(*this);
}

void ColourReconnector::rearrange(ClusterVector & clusters) {
  if (_clreco == 0) return;

  // need at least two clusters
  if (clusters.size() < 2) return;

  // do the colour reconnection
  switch (_algorithm) {
  case 0:
    _doRecoPlain(clusters);
    break;
  case 1:
    _doRecoStatistical(clusters);
    break;
  case 2:
    _doRecoBaryonic(clusters);
    break;
  case 3:
    _doRecoBaryonicMesonic(clusters);
    break;
  }
}


Energy2 ColourReconnector::_clusterMassSum(const PVector & q, 
                                           const PVector & aq) const {
  const size_t nclusters = q.size();
  assert (aq.size() == nclusters);
  Energy2 sum = ZERO;
  for (size_t i = 0; i < nclusters; i++)
    sum += ( q[i]->momentum() + aq[i]->momentum() ).m2();
  return sum;
}

double ColourReconnector::_displacement(tcPPtr p, tcPPtr q) const {
  double deltaRap = (p->rapidity() - q->rapidity());
  double deltaPhi = (p->momentum().phi() - q->momentum().phi());

  return sqrt(deltaRap * deltaRap + deltaPhi * deltaPhi);
}


double ColourReconnector::_displacementBaryonic(tcPPtr q1, tcPPtr q2, tcPPtr q3) const {
  if (_junctionMBCR) {
    /**
     * Junction-like option i.e. displacement
     * from "junction centre" (mean rapidity/phi)
     */
    double rap1=q1->rapidity();
    double rap2=q2->rapidity();
    double rap3=q3->rapidity();

    double phi1=q1->momentum().phi();
    double phi2=q2->momentum().phi();
    double phi3=q3->momentum().phi();
    double meanRap=(rap1 + rap2 + rap3)/3.0;
    double meanPhi=(phi1 + phi2 + phi3)/3.0;
    double delR;

    delR  = sqrt( (rap1-meanRap)*(rap1-meanRap) + (phi1-meanPhi)*(phi1-meanPhi) );
    delR += sqrt( (rap2-meanRap)*(rap2-meanRap) + (phi2-meanPhi)*(phi2-meanPhi) );
    delR += sqrt( (rap3-meanRap)*(rap3-meanRap) + (phi3-meanPhi)*(phi3-meanPhi) );
    return delR;
  } else {
    /* just summing up all possible 2 quark displacements */
    return _displacement(q1, q2) + _displacement(q1, q3) + _displacement(q2, q3);
  }
}
bool ColourReconnector::_containsColour8(const ClusterVector & cv,
                                         const vector<size_t> & P) const {
  assert (P.size() == cv.size());
  for (size_t i = 0; i < cv.size(); i++) {
    tcPPtr p = cv[i]->colParticle();
    tcPPtr q = cv[P[i]]->antiColParticle();
    if (_isColour8(p, q)) return true;
  }
  return false;
}


void ColourReconnector::_doRecoStatistical(ClusterVector & cv) const {

  const size_t nclusters = cv.size();

  // initially, enumerate (anti)quarks as given in the cluster vector
  ParticleVector q, aq;
  for (size_t i = 0; i < nclusters; i++) {
    q.push_back( cv[i]->colParticle() );
    aq.push_back( cv[i]->antiColParticle() );
  }

  // annealing scheme
  Energy2 t, delta;
  Energy2 lambda = _clusterMassSum(q,aq);
  const unsigned _ntries = _triesPerStepFactor * nclusters;
  
  // find appropriate starting temperature by measuring the largest lambda
  // difference in some dry-run random rearrangements
  {
    vector<Energy2> typical;
    for (int i = 0; i < 10; i++) {
      const pair <int,int> toswap = _shuffle(q,aq,5);
      ParticleVector newaq = aq;
      swap (newaq[toswap.first], newaq[toswap.second]);
      Energy2 newlambda = _clusterMassSum(q,newaq);
      typical.push_back( abs(newlambda - lambda) );
    }
    t = _initTemp * Math::median(typical);
  }

  // anneal in up to _annealingSteps temperature steps
  for (unsigned step = 0; step < _annealingSteps; step++) {

    // For this temperature step, try to reconnect _ntries times. Stop the
    // algorithm if no successful reconnection happens.
    unsigned nSuccess = 0;
    for (unsigned it = 0; it < _ntries; it++) {
      
      // make a random rearrangement
      const unsigned maxtries = 10;
      const pair <int,int> toswap = _shuffle(q,aq,maxtries);
      const int i = toswap.first;
      const int j = toswap.second;

      // stop here if we cannot find any allowed reconfiguration
      if (i == -1) break;

      // create a new antiquark vector with the two partons swapped
      ParticleVector newaq = aq;
      swap (newaq[i], newaq[j]);

      // Check if lambda would decrease. If yes, accept the reconnection. If no,
      // accept it only with a probability given by the current Boltzmann
      // factor. In the latter case we set p = 0 if the temperature is close to
      // 0, to avoid division by 0.
      Energy2 newlambda = _clusterMassSum(q,newaq);
      delta = newlambda - lambda;
      double prob = 1.0;
      if (delta > ZERO) prob = ( abs(t) < 1e-8*MeV2 ) ? 0.0 : exp(-delta/t);
      if (UseRandom::rnd() < prob) {
        lambda = newlambda;
        swap (newaq, aq);
        nSuccess++;
      }
    }
    if (nSuccess == 0) break;

    // reduce temperature
    t *= _annealingFactor;
  }

  // construct the new cluster vector
  ClusterVector newclusters;
  for (size_t i = 0; i < nclusters; i++) {
    ClusterPtr cl = new_ptr( Cluster( q[i], aq[i] ) );
    newclusters.push_back(cl);
  }
  swap(newclusters,cv);
  return;
}


void ColourReconnector::_doRecoPlain(ClusterVector & cv) const {

  ClusterVector newcv = cv;

  // try to avoid systematic errors by randomising the reconnection order
  long (*p_irnd)(long) = UseRandom::irnd;
  random_shuffle( newcv.begin(), newcv.end(), p_irnd );

  // iterate over all clusters
  for (CluVecIt cit = newcv.begin(); cit != newcv.end(); cit++) {
    // find the cluster which, if reconnected with *cit, would result in the
    // smallest sum of cluster masses
    // NB this method returns *cit if no reconnection partner can be found
    CluVecIt candidate = _findRecoPartner(cit, newcv);

    // skip this cluster if no possible reshuffling partner can be found
    if (candidate == cit) continue;

    // accept the reconnection with probability _preco.
    if (UseRandom::rnd() < _preco) {

      pair <ClusterPtr,ClusterPtr> reconnected = _reconnect(*cit, *candidate);

      // Replace the clusters in the ClusterVector. The order of the
      // colour-triplet partons in the cluster vector is retained here.

      // replace *cit by reconnected.first
      *cit = reconnected.first;

      // replace candidate by reconnected.second
      *candidate = reconnected.second;
    }
  }

  swap(cv,newcv);
  return;
}

namespace {
inline bool hasDiquark(CluVecIt cit) {
  for (int i = 0; i<(*cit)->numComponents(); i++) {
    if (DiquarkMatcher::Check(*((*cit)->particle(i)->dataPtr())))
      return true;
  }
  return false;
}
}


// Implementation of the baryonic reconnection algorithm
void ColourReconnector::_doRecoBaryonic(ClusterVector & cv) const {

  ClusterVector newcv = cv;

  ClusterVector deleted; deleted.reserve(cv.size());

  // try to avoid systematic errors by randomising the reconnection order
  long (*p_irnd)(long) = UseRandom::irnd;
  random_shuffle(newcv.begin(), newcv.end(), p_irnd);

  // iterate over all clusters
  for (CluVecIt cit = newcv.begin(); cit != newcv.end(); ++cit) {
    //avoid clusters already containing diuarks
    if (hasDiquark(cit)) continue;

    //skip the cluster to be deleted later 3->2 cluster
    if (find(deleted.begin(), deleted.end(), *cit) != deleted.end())
      continue;

    // Skip all found baryonic clusters, this biases the algorithm but implementing
    // something like re-reconnection is ongoing work
    if ((*cit)->numComponents()==3) continue;

    // Find a candidate suitable for reconnection
    CluVecIt baryonic1, baryonic2;
    bool isBaryonicCandidate = false;
    CluVecIt candidate = _findPartnerBaryonic(cit, newcv,
                         isBaryonicCandidate,
                         deleted,
                         baryonic1, baryonic2);

    // skip this cluster if no possible reconnection partner can be found
    if (!isBaryonicCandidate && candidate==cit)
      continue;

    if (isBaryonicCandidate
        && UseRandom::rnd() < _precoBaryonic) {
      deleted.push_back(*baryonic2);

      // Function that does the reconnection from 3 -> 2 clusters
      ClusterPtr b1, b2;
      _makeBaryonicClusters(*cit, *baryonic1, *baryonic2, b1, b2);

      *cit = b1;
      *baryonic1 = b2;

      // Baryonic2 is easily skipped in the next loop
    }

    // Normal 2->2 Colour reconnection
    if (!isBaryonicCandidate
        && UseRandom::rnd() < _preco) {
      auto reconnected = _reconnect(*cit, *candidate);
      *cit = reconnected.first;
      *candidate = reconnected.second;
    }
  }

  // create a new vector of clusters except for the ones which are "deleted" during
  // baryonic reconnection
  ClusterVector clustervector;
  for (const auto & cluster : newcv)
    if (find(deleted.begin(),
             deleted.end(), cluster) == deleted.end())
      clustervector.push_back(cluster);

  swap(cv, clustervector);



}


bool ColourReconnector::_clustersFarApart( const std::vector<CluVecIt> & clu ) const {
  int Ncl=clu.size();
  assert(Ncl<=3);
  if (Ncl==1) {
    return false;
  } else if (Ncl==2) {
    // TODO: keep turned off all until things are more clear
    // BUG: space-time difference compared to _maxDistance
    // if (((*clu[0])->vertex()-(*clu[1])->vertex()).m() >_maxDistance) return true;
    // Causal selection if desired
    // if (((*clu[0])->vertex()-(*clu[1])->vertex()).m() >ZERO) return true;
  } else if (Ncl==3) {
    // TODO: keep turned off all until things are more clear
    // BUG: space-time difference compared to _maxDistance
    // if (((*clu[0])->vertex()-(*clu[1])->vertex()).m()> _maxDistance) return true;
    // Causal selection if desired
    // if (((*clu[0])->vertex()-(*clu[1])->vertex()).m()> ZERO) return true;
    // if (((*clu[1])->vertex()-(*clu[2])->vertex()).m()> _maxDistance) return true;
    // if (((*clu[1])->vertex()-(*clu[2])->vertex()).m()> ZERO) return true;
    // if (((*clu[0])->vertex()-(*clu[2])->vertex()).m()> _maxDistance) return true;
    // if (((*clu[0])->vertex()-(*clu[2])->vertex()).m()> ZERO) return true;
  }

  return false;
}



void ColourReconnector::_doReco2BeamClusters(ClusterVector & cv) const {
  // try other option
  PPtr p1Di=(cv[0])->colParticle();
  PPtr p2Di=(cv[1])->colParticle();

  PPtr p1Q=(cv[0])->antiColParticle();
  PPtr p2Q=(cv[1])->antiColParticle();

  double min_dist=_displacement(p1Di,p1Q)+_displacement(p2Di,p2Q);

  if ((_displacement(p1Di,p2Q)+_displacement(p1Di,p2Q))<min_dist) {
    _reconnect(cv[0],cv[1]);
  }
  return;
}



void ColourReconnector::_doRecoBaryonicMesonic(ClusterVector & cv) const {
  if (cv.size() < 3) {
    /*
     * if the option _cr2BeamClusters!=0 is chosen then we try to
     * colour reconnect the special case of 2 beam clusters with
     * probability 1.0 if there is a better _displacement
     * */
    if( _cr2BeamClusters && cv.size()==2 ) _doReco2BeamClusters(cv);
    return;
  }

  ClusterVector newcv = cv;
  newcv.reserve(2*cv.size());

  ClusterVector deleted;
  deleted.reserve(cv.size());

  // counters for numbers of mesons and baryons selected
  unsigned num_meson = 0;
  unsigned num_baryon = 0;

  // vector of selected clusters
  std::vector<CluVecIt>  sel;

  unsigned number_of_tries = _stepFactor*cv.size()*cv.size();
  if (number_of_tries<1) number_of_tries=1;

  long (*p_irnd)(long) = UseRandom::irnd;
  for (unsigned reconnections_tries = 0; reconnections_tries < number_of_tries; reconnections_tries++) {
    num_meson = 0;
    num_baryon = 0;

    // flag if we are able to find a suitable combinations of clusters
    bool _found = false;

    // Shuffle list of clusters to avoid systematic bias in cluster selection
    random_shuffle(newcv.begin(), newcv.end(), p_irnd);

    // loop over clustervector to find CR candidates
    for (CluVecIt cit = newcv.begin(); cit != newcv.end(); ++cit) {

      // skip the clusters to be deleted later from 3->2 cluster CR
      if (find(deleted.begin(), deleted.end(), *cit) != deleted.end()) continue;

      // avoid clusters already containing diuarks
      if (hasDiquark(cit)) continue;

      // add to selection
      sel.push_back(cit);

      if (_clustersFarApart(sel)) {
        // reject far appart CR
        // TODO: after discussion maybe to be omitted
        sel.pop_back();
        continue;
      }

      bool isMeson=((*cit)->numComponents() == 2);

      if ( isMeson && (num_meson ==0|| num_meson==1) && num_baryon ==0) {
        num_meson++;
        /**
         * now we habe either 1 or 2 mesonic clusters and have to continue
         */
        continue;
      } else if ( isMeson && (num_baryon == 1 || num_meson ==2)) {
        num_meson++;
        _found = true;
        /**
         * we have either 3 mesonic or 1 mesonic and 1 baryonic cluster
         * and try to colour reconnect
         */
        break;
      } else if (num_baryon ==0 && num_meson==0) {
        num_baryon++;
        /**
         * now we have 1 baryonic cluster and have to continue
         */
        continue;
      } else if (num_meson == 2) {
        /**
         * we already have 2 mesonic clusters and dont want a baryonic one
         * since this is an invalid selection
         */
        // remove previously added cluster
        sel.pop_back();
        continue;
      } else {
        num_baryon++;
        _found = true;
        /**
         * now we have either 2 baryonic clusters or 1 mesonic and 1 baryonic cluster
         * and try to colour reconnect
         */
        break;
      }
    }

    // added for more efficent rejection if some reco probabilities are 0
    if ( _found ) {

      // reject MBtoMB candidates if _precoMB_MB=0
      if ( _precoMB_MB == 0 && (num_baryon == 1 && num_meson == 1) ) {
        _found=false;
      }

      // reject BbarBto3M candidates if _precoBbarB_3M=0
      if ( _precoBbarB_3M== 0 && num_baryon == 2 ) {
        bool isBbarBto3Mcandidate=(
                                    (*sel[0])->particle(0)->hasColour() && (*sel[1])->particle(0)->hasColour(true) )
                                  || ( (*sel[0])->particle(0)->hasColour(true) && (*sel[1])->particle(0)->hasColour() );

        if ( isBbarBto3Mcandidate) _found=false;
      }

      // reject 2Bto2B candidates if _preco2B_2B=0
      if ( _preco2B_2B == 0 && num_baryon == 2 ) {
        bool is2Bto2Bcandidate=(
                                 (*sel[0])->particle(0)->hasColour() && (*sel[1])->particle(0)->hasColour() )
                               || ( (*sel[0])->particle(0)->hasColour(true) && (*sel[1])->particle(0)->hasColour(true) );

        if ( is2Bto2Bcandidate ) _found=false;
      }
    }
    // were we able to find a combination?
    if (_found==false) {
      // clear the selection if we did not find a valid set of  clusters
      sel.erase(sel.begin(), sel.end());
      continue;
    }
    assert(sel.size()<4);
    assert(sel.size()>1);

    string kind_of_reco = "";
    int reco_info[3];

    // find best CR option for the selection
    _findbestreconnectionoption(sel, num_baryon, kind_of_reco, reco_info);

    if (kind_of_reco == "") {
      // no reconnection was found
      sel.erase(sel.begin(), sel.end());
      continue;
    } else if (kind_of_reco == "3Mto3M" && UseRandom::rnd() < _preco3M_3M) {
      // 3Mto3M colour reconnection
      auto reconnected = _reconnect3Mto3M(*sel[0], *sel[1], *sel[2],
                                          reco_info);
      (*sel[0]) = std::get<0>(reconnected);
      (*sel[1]) = std::get<1>(reconnected);
      (*sel[2]) = std::get<2>(reconnected);
    } else if (kind_of_reco=="2Bto3M" && UseRandom::rnd() < _precoBbarB_3M) {
      // antibaryonic and baryonic to 3 mesonic reconnecion
      auto reconnected = _reconnectBbarBto3M(*sel[0], *sel[1],
                                             reco_info[0], reco_info[1], reco_info[2]);
      (*sel[0]) = std::get<0>(reconnected);
      (*sel[1]) = std::get<1>(reconnected);
      newcv.push_back(std::get<2>(reconnected));
    } else if (kind_of_reco=="3Mto2B" && UseRandom::rnd() < _preco3M_BBbar) {
      // 3 mesonic to antibaryonic and baryonic reconnection
      ClusterPtr b1, b2;
      _makeBaryonicClusters(*sel[0], *sel[1], *sel[2], b1, b2);
      (*sel[0]) = b1;
      (*sel[1]) = b2;
      deleted.push_back(*sel[2]);
    } else if (kind_of_reco=="2Bto2B" && UseRandom::rnd() < _preco2B_2B) {
      // 2 (anti)baryonic to 2 (anti)baryonic reconnection
      auto reconnected = _reconnect2Bto2B(*sel[0], *sel[1],
                                          reco_info[0], reco_info[1]);
      (*sel[0]) = reconnected.first;
      (*sel[1]) = reconnected.second;
    } else if (kind_of_reco=="MBtoMB" && UseRandom::rnd() < _precoMB_MB) {
      // (anti)baryonic and mesonic to (anti)baryonic and mesonic reconnection
      auto reconnected = _reconnectMBtoMB(*sel[0], *sel[1],
                                          reco_info[0]);
      (*sel[0]) = reconnected.first;
      (*sel[1]) = reconnected.second;
    }
    // erase the sel-vector
    sel.erase(sel.begin(), sel.end());
  }

  // write to clustervector new CR'd clusters and deleting
  // all deleted clusters
  ClusterVector clustervector;
  for (const auto & cluster : newcv)
    if (find(deleted.begin(), deleted.end(), cluster) == deleted.end())
      clustervector.push_back(cluster);

  swap(cv, clustervector);
}

namespace {

double calculateRapidityRF(const Lorentz5Momentum & q1,
                           const Lorentz5Momentum & p2) {
  //calculate rapidity wrt the direction of q1
  //angle between the particles in the RF of cluster of q1

  // calculate the z component of p2 w.r.t the direction of q1
  if(q1.rho2()==ZERO) return 0.;
  const Energy pz = p2.vect() * q1.vect().unit();
  if ( pz == ZERO ) return 0.;

  // Transverse momentum of p2 w.r.t the direction of q1
  const Energy pt = sqrt(p2.vect().mag2() - sqr(pz));

  // Transverse mass pf p2 w.r.t to the direction of q1
  const Energy mtrans = sqrt(p2.mass()*p2.mass() + (pt*pt));

  // Correct formula
  const double y2 = log((p2.t() + abs(pz))/mtrans);

  return ( pz < ZERO ) ? -y2 : y2;
}

}


void ColourReconnector::_findbestreconnectionoption(std::vector<CluVecIt> & cls, const unsigned & baryonic,
    string & kind_of_reco, int (&reco_info)[3]) const {
  double min_displacement;
  if (baryonic==0) {
    // case with 3 mesonic clusters
    assert(cls.size()==3);

    // calculate the initial displacement sum
    min_displacement  = _mesonToBaryonFactor * _displacement((*cls[0])->particle(0), (*cls[0])->particle(1));
    min_displacement += _mesonToBaryonFactor * _displacement((*cls[1])->particle(0), (*cls[1])->particle(1));
    min_displacement += _mesonToBaryonFactor * _displacement((*cls[2])->particle(0), (*cls[2])->particle(1));

    // find best CR reco_info and kind_of_reco
    _3MtoXreconnectionfinder(cls,
                             reco_info[0], reco_info[1], reco_info[2], min_displacement, kind_of_reco);
    /**
     * kind_of_reco either "3Mto3M" or "3Mto2B" (or "" if no better configuration is found)
     * case 3Mto3M:    the coloured particle of the i-th cluster forms a new cluster with the
     * 				   antiparticle of the reco_info[i]-th cluster
     * case 3MtoBbarB: all 3 (anti)coloured particle form a new (anti)baryonic cluster
     */
  } else if (baryonic == 1) {
    // case 1 baryonic and 1 mesonic cluster
    assert(cls.size()==2);

    // make mesonic cluster always the cls[0]
    if ((*cls[0])->numComponents() == 3) {
      ClusterPtr zw = *cls[0];
      *cls[0] = *cls[1];
      *cls[1] = zw;
    }

    // calculate the initial displacement sum
    min_displacement  = _mesonToBaryonFactor *_displacement((*cls[0])->particle(0), (*cls[0])->particle(1));
    min_displacement += _displacementBaryonic((*cls[1])->particle(0), (*cls[1])->particle(1), (*cls[1])->particle(2));

    // find best CR reco_info and kind_of_reco
    _BMtoBMreconnectionfinder(*cls[0], *cls[1],
                              reco_info[0], min_displacement, kind_of_reco);
    /**
     * reco_info[0] is the index of the (anti)quarks of the baryonic cluster cls[1], which should
     * be swapped with the (anti)quarks of the mesonic cluster cls[0]
     */

  } else {
    assert(baryonic==2);
    assert(cls.size()==2);

    // calculate the initial displacement sum
    min_displacement  = _displacementBaryonic((*cls[0])->particle(0), (*cls[0])->particle(1), (*cls[0])->particle(2));
    min_displacement += _displacementBaryonic((*cls[1])->particle(0), (*cls[1])->particle(1), (*cls[1])->particle(2));

    // case 2 (anti)baryonic clusters to 2 other (anti)baryonic clusters
    if (      ( (*cls[0])->particle(0)->hasColour()     && (*cls[1])->particle(0)->hasColour()     )
              || ( (*cls[0])->particle(0)->hasColour(true) && (*cls[1])->particle(0)->hasColour(true) ) ) {
      // find best CR reco_info and kind_of_reco
      _2Bto2BreconnectionFinder(*cls[0], *cls[1],
                                reco_info[0], reco_info[1], min_displacement, kind_of_reco);
      /**
       * swap the reco_info[0]-th particle of the first cluster in the vector with the
       * reco_info[1]-th particle of the second cluster
       */
    } else {
      // case 1 baryonic and 1 antibaryonic cluster to 3 mesonic clusters

      // find best CR reco_info and kind_of_reco
      _BbarBto3MreconnectionFinder(*cls[0], *cls[1],
                                   reco_info[0], reco_info[1], reco_info[2], min_displacement, kind_of_reco);
      /**
       * the i-th particle of the first cluster form a new mesonic cluster with the
       * reco_info[i]-th particle of the second cluster
       */
    }
  }
  return;
}


CluVecIt ColourReconnector::_findPartnerBaryonic(
  CluVecIt cl, ClusterVector & cv,
  bool & baryonicCand,
  const ClusterVector& deleted,
  CluVecIt &baryonic1,
  CluVecIt &baryonic2 ) const {

  using Constants::pi;
  using Constants::twopi;

  // Returns a candidate for possible reconnection
  CluVecIt candidate = cl;

  bool bcand = false;

  double maxrap = 0.0;
  double minrap = 0.0;
  double maxrapNormal = 0.0;
  double minrapNormal = 0.0;
  double maxsumnormal = 0.0;

  double maxsum = 0.0;
  double secondsum = 0.0;


  // boost into RF of cl
  Lorentz5Momentum cl1 = (*cl)->momentum();
  const Boost boostv(-cl1.boostVector());
  cl1.boost(boostv);
  // boost constituents of cl into RF of cl
  Lorentz5Momentum p1col = (*cl)->colParticle()->momentum();
  Lorentz5Momentum p1anticol = (*cl)->antiColParticle()->momentum();
  p1col.boost(boostv);
  p1anticol.boost(boostv);

  for (CluVecIt cit=cv.begin(); cit != cv.end(); ++cit) {
    //avoid looping over clusters containing diquarks
    if ( hasDiquark(cit) ) continue;
    if ( (*cit)->numComponents()==3 ) continue;
    if ( cit==cl ) continue;

    //skip the cluster to be deleted later 3->2 cluster
    if ( find(deleted.begin(), deleted.end(), *cit) != deleted.end() )
      continue;

    if ( (*cl)->isBeamCluster() && (*cit)->isBeamCluster() )
      continue;

    // stop it putting far apart clusters together
    // BUG: space-time difference compared to _maxDistance
    // if (((**cl).vertex()-(**cit).vertex()).m() >_maxDistance)
    // Causal selection if desired
    // if (((**cl).vertex()-(**cit).vertex()).m() >ZERO) continue;

    const bool Colour8 =
      _isColour8( (*cl)->colParticle(), (*cit)->antiColParticle() )
      ||
      _isColour8( (*cit)->colParticle(), (*cl)->antiColParticle() ) ;
    if ( Colour8 ) continue;


    // boost constituents of cit into RF of cl
    Lorentz5Momentum p2col = (*cit)->colParticle()->momentum();
    Lorentz5Momentum p2anticol = (*cit)->antiColParticle()->momentum();

    p2col.boost(boostv);
    p2anticol.boost(boostv);

    // calculate the rapidity of the other constituents of the clusters
    // w.r.t axis of p1anticol.vect.unit
    const double rapq = calculateRapidityRF(p1anticol,p2col);
    const double rapqbar = calculateRapidityRF(p1anticol,p2anticol);

    // configuration for normal CR
    if ( rapq > 0.0 && rapqbar < 0.0
         && rapq > maxrap
         && rapqbar < minrap ) {
      maxrap = rapq;
      minrap = rapqbar;
      //sum of rapidities of quarks
      const double normalsum = abs(rapq) + abs(rapqbar);
      if ( normalsum > maxsumnormal ) {
        maxsumnormal = normalsum;
        maxrapNormal = rapq;
        minrapNormal = rapqbar;
        bcand = false;
        candidate = cit;
      }
    }

    if ( rapq < 0.0 && rapqbar >0.0
         && rapqbar > maxrapNormal
         && rapq < minrapNormal ) {
      maxrap = rapqbar;
      minrap = rapq;
      const double sumrap = abs(rapqbar) + abs(rapq);
      // first candidate gets here. If second baryonic candidate has higher Ysum than the first
      // one, the second candidate becomes the first one and the first the second.
      if (sumrap > maxsum) {
        if (maxsum != 0) {
          baryonic2 = baryonic1;
          baryonic1 = cit;
          bcand = true;
        } else {
          baryonic1 = cit;
        }
        maxsum = sumrap;
      } else {
        if (sumrap > secondsum && sumrap != maxsum) {
          secondsum = sumrap;
          bcand = true;
          baryonic2 = cit;
        }
      }
    }

  }

  if (bcand == true) {
    baryonicCand = true;
  }

  return candidate;
}


CluVecIt ColourReconnector::_findRecoPartner(CluVecIt cl,
    ClusterVector & cv) const {

  CluVecIt candidate = cl;
  Energy minMass = 1*TeV;
  for (CluVecIt cit=cv.begin(); cit != cv.end(); ++cit) {

    // don't even look at original cluster
    if (cit==cl) continue;

    // don't allow colour octet clusters
    if ( _isColour8( (*cl)->colParticle(),
                     (*cit)->antiColParticle() )  ||
         _isColour8( (*cit)->colParticle(),
                     (*cl)->antiColParticle() ) ) {
      continue;
    }

    // stop it putting beam remnants together
    if ((*cl)->isBeamCluster() && (*cit)->isBeamCluster()) continue;

    // stop it putting far apart clusters together
    // BUG: space-time difference compared to _maxDistance
    // if (((**cl).vertex()-(**cit).vertex()).m()>_maxDistance) continue;
    // Causal selection if desired
    // if (((**cl).vertex()-(**cit).vertex()).m()>ZERO) continue;

    // momenta of the old clusters
    Lorentz5Momentum p1 = (*cl)->colParticle()->momentum() +
                          (*cl)->antiColParticle()->momentum();
    Lorentz5Momentum p2 = (*cit)->colParticle()->momentum() +
                          (*cit)->antiColParticle()->momentum();

    // momenta of the new clusters
    Lorentz5Momentum p3 = (*cl)->colParticle()->momentum() +
                          (*cit)->antiColParticle()->momentum();
    Lorentz5Momentum p4 = (*cit)->colParticle()->momentum() +
                          (*cl)->antiColParticle()->momentum();

    Energy oldMass = abs( p1.m() ) + abs( p2.m() );
    Energy newMass = abs( p3.m() ) + abs( p4.m() );


    if ( newMass < oldMass && newMass < minMass ) {
      minMass = newMass;
      candidate = cit;
    }
  }

  return candidate;
}

// forms two baryonic clusters from three clusters
void ColourReconnector::_makeBaryonicClusters(
  ClusterPtr &c1, ClusterPtr &c2,
  ClusterPtr &c3,
  ClusterPtr &newcluster1,
  ClusterPtr &newcluster2) const {

  //make sure they all have 2 components
  assert(c1->numComponents()==2);
  assert(c2->numComponents()==2);
  assert(c3->numComponents()==2);
  //abandon children
  c1->colParticle()->abandonChild(c1);
  c1->antiColParticle()->abandonChild(c1);
  c2->colParticle()->abandonChild(c2);
  c2->antiColParticle()->abandonChild(c2);
  c3->colParticle()->abandonChild(c3);
  c3->antiColParticle()->abandonChild(c3);

  newcluster1 = new_ptr(Cluster(c1->colParticle(),c2->colParticle(), c3->colParticle()));
  c1->colParticle()->addChild(newcluster1);
  c2->colParticle()->addChild(newcluster1);
  c3->colParticle()->addChild(newcluster1);
  newcluster1->setVertex(LorentzPoint());

  newcluster2 = new_ptr(Cluster(c1->antiColParticle(), c2->antiColParticle(),
                                c3->antiColParticle()));
  c1->antiColParticle()->addChild(newcluster2);
  c2->antiColParticle()->addChild(newcluster2);
  c3->antiColParticle()->addChild(newcluster2);
  newcluster2->setVertex(LorentzPoint());
}

pair <ClusterPtr,ClusterPtr>
ColourReconnector::_reconnect2Bto2B(ClusterPtr &c1, ClusterPtr &c2, const int s1, const int s2) const {

  // form the first new cluster

  // separate the quarks from their original cluster
  c1->particleB((s1+1)%3)->abandonChild(c1);
  c1->particleB((s1+2)%3)->abandonChild(c1);
  c2->particleB(s2)->abandonChild(c2);

  // now the new cluster
  ClusterPtr newCluster1 = new_ptr(Cluster(c1->particleB((s1+1)%3), c1->particleB((s1+2)%3), c2->particleB(s2)));

  c1->particleB((s1+1)%3)->addChild(newCluster1);
  c1->particleB((s1+2)%3)->addChild(newCluster1);
  c2->particleB(s2)->addChild(newCluster1);

  // set new vertex
  newCluster1->setVertex(LorentzPoint());

  // set beam remnants for new cluster
  if (c1->isBeamRemnant((s1+1)%3)) newCluster1->setBeamRemnant(0, true);
  if (c1->isBeamRemnant((s1+2)%3)) newCluster1->setBeamRemnant(1, true);
  if (c2->isBeamRemnant(s2)) newCluster1->setBeamRemnant(2, true);

  // for the second cluster same  procedure
  c2->particleB((s2+1)%3)->abandonChild(c2);
  c2->particleB((s2+2)%3)->abandonChild(c2);
  c1->particleB(s1)->abandonChild(c1);

  ClusterPtr newCluster2 = new_ptr(Cluster(c2->particleB((s2+1)%3), c2->particleB((s2+2)%3), c1->particleB(s1)));

  c2->particleB((s2+1)%3)->addChild(newCluster2);
  c2->particleB((s2+2)%3)->addChild(newCluster2);
  c1->particleB(s1)->addChild(newCluster2);

  newCluster2->setVertex(LorentzPoint());

  if (c2->isBeamRemnant((s2+1)%3)) newCluster2->setBeamRemnant(0, true);
  if (c2->isBeamRemnant((s2+2)%3)) newCluster2->setBeamRemnant(1, true);
  if (c1->isBeamRemnant(s1)) newCluster2->setBeamRemnant(2, true);

  return pair <ClusterPtr, ClusterPtr> (newCluster1, newCluster2);
}


std::tuple  <ClusterPtr, ClusterPtr, ClusterPtr>
ColourReconnector::_reconnectBbarBto3M(ClusterPtr & c1, ClusterPtr & c2, const int s0, const int s1, const int s2) const {
  // make sure they all have 3 components
  assert(c1->numComponents()==3);
  assert(c2->numComponents()==3);

  // first Cluster
  c1->particleB(0)->abandonChild(c1);
  c2->particleB(s0)->abandonChild(c2);

  ClusterPtr newCluster1 = new_ptr(Cluster(c1->particleB(0), c2->particleB(s0)));

  c1->particleB(0)->addChild(newCluster1);
  c2->particleB(s0)->addChild(newCluster1);

  // set new vertex
  newCluster1->setVertex(0.5*(c1->particleB(0)->vertex() + c2->particleB(s0)->vertex()));

  // set beam remnants for new cluster
  if (c1->isBeamRemnant(0)) newCluster1->setBeamRemnant(0, true);
  if (c2->isBeamRemnant(s0)) newCluster1->setBeamRemnant(1, true);

  // same for second cluster
  c1->particleB(1)->abandonChild(c1);
  c2->particleB(s1)->abandonChild(c2);

  ClusterPtr newCluster2 = new_ptr(Cluster(c1->particleB(1), c2->particleB(s1)));

  c1->particleB(1)->addChild(newCluster2);
  c2->particleB(s1)->addChild(newCluster2);

  newCluster2->setVertex(0.5*(c1->particleB(1)->vertex() + c2->particleB(s1)->vertex()));

  if (c1->isBeamRemnant(1)) newCluster2->setBeamRemnant(0, true);
  if (c2->isBeamRemnant(s1)) newCluster2->setBeamRemnant(1, true);

  // same for third cluster
  c1->particleB(2)->abandonChild(c1);
  c2->particleB(s2)->abandonChild(c2);

  ClusterPtr newCluster3 = new_ptr(Cluster(c1->particleB(2), c2->particleB(s2)));

  c1->particleB(2)->addChild(newCluster3);
  c2->particleB(s2)->addChild(newCluster3);

  newCluster3->setVertex(0.5*(c1->particleB(2)->vertex() + c2->particleB(s2)->vertex()));

  if (c1->isBeamRemnant(2)) newCluster3->setBeamRemnant(0, true);
  if (c2->isBeamRemnant(s2)) newCluster3->setBeamRemnant(1, true);

  return std::tuple  <ClusterPtr, ClusterPtr, ClusterPtr> (newCluster1, newCluster2, newCluster3);
}

pair <ClusterPtr,ClusterPtr>
ColourReconnector::_reconnect(ClusterPtr &c1, ClusterPtr &c2) const {

  // choose the other possibility to form two clusters from the given
  // constituents

  assert(c1->numComponents()==2);
  assert(c2->numComponents()==2);
  int c1_col(-1),c1_anti(-1),c2_col(-1),c2_anti(-1);
  for(unsigned int ix=0; ix<2; ++ix) {
    if     (c1->particle(ix)->hasColour(false)) c1_col  = ix;
    else if(c1->particle(ix)->hasColour(true )) c1_anti = ix;
    if     (c2->particle(ix)->hasColour(false)) c2_col  = ix;
    else if(c2->particle(ix)->hasColour(true )) c2_anti = ix;
  }
  assert(c1_col>=0&&c2_col>=0&&c1_anti>=0&&c2_anti>=0);

  c1->colParticle()->abandonChild(c1);
  c2->antiColParticle()->abandonChild(c2);

  ClusterPtr newCluster1
    = new_ptr( Cluster( c1->colParticle(), c2->antiColParticle() ) );

  c1->colParticle()->addChild(newCluster1);
  c2->antiColParticle()->addChild(newCluster1);

  /*
   * TODO: Questionable setting of the vertex
   * */
  newCluster1->setVertex(0.5*(c1->colParticle()->vertex() +
                              c2->antiColParticle()->vertex()));

  if(c1->isBeamRemnant(c1_col )) newCluster1->setBeamRemnant(0,true);
  if(c2->isBeamRemnant(c2_anti)) newCluster1->setBeamRemnant(1,true);

  c1->antiColParticle()->abandonChild(c1);
  c2->colParticle()->abandonChild(c2);

  ClusterPtr newCluster2
    = new_ptr( Cluster( c2->colParticle(), c1->antiColParticle() ) );

  c1->antiColParticle()->addChild(newCluster2);
  c2->colParticle()->addChild(newCluster2);

  /*
   * TODO: Questionable setting of the vertex
   * */
  newCluster2->setVertex(0.5*(c2->colParticle()->vertex() +
                              c1->antiColParticle()->vertex()));

  if(c2->isBeamRemnant(c2_col )) newCluster2->setBeamRemnant(0,true);
  if(c1->isBeamRemnant(c1_anti)) newCluster2->setBeamRemnant(1,true);

  return pair <ClusterPtr,ClusterPtr> (newCluster1, newCluster2);
}

std::tuple  <ClusterPtr, ClusterPtr, ClusterPtr>
ColourReconnector::_reconnect3Mto3M(ClusterPtr & c1, ClusterPtr & c2, ClusterPtr & c3, const int infos [3]) const {
  // check if mesonic clusters
  assert(c1->numComponents()==2);
  assert(c2->numComponents()==2);
  assert(c3->numComponents()==2);

  ClusterVector oldclusters = {c1, c2, c3};
  ClusterVector newclusters;

  for (int i=0; i<3; i++) {

    int c1_col=-1;
    int c2_anticol=-1;

    // get which index is coloured and which anticolour
    for (unsigned int ix=0; ix<2; ++ix) {
      if (oldclusters[i]->particle(ix)->hasColour(false)) c1_col  = ix;
      if (oldclusters[infos[i]]->particle(ix)->hasColour(true)) c2_anticol  = ix;
    }

    assert(c1_col>=0);
    assert(c2_anticol>=0);

    oldclusters[i]->colParticle()->abandonChild(oldclusters[i]);
    oldclusters[infos[i]]->antiColParticle()->abandonChild(oldclusters[infos[i]]);

    // form new cluster
    ClusterPtr newCluster = new_ptr(Cluster(oldclusters[i]->colParticle(), oldclusters[infos[i]]->antiColParticle()));

    oldclusters[i]->colParticle()->addChild(newCluster);
    oldclusters[infos[i]]->antiColParticle()->addChild(newCluster);

    // set new vertex
    newCluster->setVertex(0.5*(oldclusters[i]->colParticle()->vertex() +
                               oldclusters[infos[i]]->antiColParticle()->vertex()));

    // set beam remnants for new cluster
    if (oldclusters[i]->isBeamRemnant(c1_col)) newCluster->setBeamRemnant(0, true);
    if (oldclusters[infos[i]]->isBeamRemnant(c2_anticol)) newCluster->setBeamRemnant(1, true);
    newclusters.push_back(newCluster);
  }
  return std::tuple <ClusterPtr, ClusterPtr, ClusterPtr> (newclusters[0], newclusters[1], newclusters[2]);
}


pair  <ClusterPtr, ClusterPtr>
ColourReconnector::_reconnectMBtoMB(ClusterPtr & c1, ClusterPtr & c2, const int s0) const {
  // make c1 the mesonic cluster
  if (c1->numComponents()==2) {
    assert(c2->numComponents()==3);
  } else {
    return _reconnectMBtoMB(c2,c1,s0);
  }

  int c1_col=-1;
  int c1_anti=-1;
  // get which index is coloured and which anticolour
  for (unsigned int ix=0; ix<2; ++ix) {
    if (c1->particle(ix)->hasColour(false)) c1_col  = ix;
    else if (c1->particle(ix)->hasColour(true)) c1_anti = ix;

  }
  assert(c1_col>=0);
  assert(c1_anti>=0);

  // pointers for the new clusters
  ClusterPtr newCluster1;
  ClusterPtr newCluster2;
  if (c2->particle(0)->hasColour()==true) {
    // first case: we have a baryonic clusters

    // first make the new mesonic cluster
    c1->antiColParticle()->abandonChild(c1);
    c2->particleB(s0)->abandonChild(c2);

    newCluster1 = new_ptr(Cluster(c1->antiColParticle(), c2->particleB(s0)));

    c1->antiColParticle()->addChild(newCluster1);
    c2->particleB(s0)->addChild(newCluster1);

    // set new vertex
    newCluster1->setVertex(0.5*(c1->antiColParticle()->vertex() +
                                c2->particleB(s0)->vertex()));

    // set beam remnants for new cluster
    if (c1->isBeamRemnant(c1_anti)) newCluster1->setBeamRemnant(0, true);
    if (c2->isBeamRemnant(s0)) newCluster1->setBeamRemnant(1, true);

    // then the baryonic one
    c1->colParticle()->abandonChild(c1);
    c2->particleB((s0+1)%3)->abandonChild(c2);
    c2->particleB((s0+2)%3)->abandonChild(c2);

    newCluster2 = new_ptr(Cluster(c1->colParticle(), c2->particleB((s0+1)%3), c2->particleB((s0+2)%3)));

    c1->colParticle()->addChild(newCluster2);
    c2->particleB((s0+1)%3)->addChild(newCluster2);
    c2->particleB((s0+2)%3)->addChild(newCluster2);

    // set new vertex
    newCluster2->setVertex(LorentzPoint());
  } else {
    // second case we have an antibaryonic cluster

    // first make the new mesonic cluster
    c1->colParticle()->abandonChild(c1);
    c2->particleB(s0)->abandonChild(c2);

    newCluster1 = new_ptr(Cluster(c1->colParticle(), c2->particleB(s0)));

    c1->colParticle()->addChild(newCluster1);
    c2->particleB(s0)->addChild(newCluster1);

    // set new vertex
    newCluster1->setVertex(0.5*(c1->colParticle()->vertex() +
                                c2->particleB(s0)->vertex()));

    // set beam remnants for new cluster
    if (c1->isBeamRemnant(c1_col)) newCluster1->setBeamRemnant(0, true);
    if (c2->isBeamRemnant(s0)) newCluster1->setBeamRemnant(1, true);

    // then the baryonic one
    c1->antiColParticle()->abandonChild(c1);
    c2->particleB((s0+1)%3)->abandonChild(c2);
    c2->particleB((s0+2)%3)->abandonChild(c2);

    newCluster2 =  new_ptr(Cluster(c1->antiColParticle(), c2->particleB((s0+1)%3), c2->particleB((s0+2)%3)));

    c1->antiColParticle()->addChild(newCluster2);
    c2->particleB((s0+1)%3)->addChild(newCluster2);
    c2->particleB((s0+2)%3)->addChild(newCluster2);

    // set new vertex
    newCluster2->setVertex(LorentzPoint());
  }
  return pair <ClusterPtr, ClusterPtr> (newCluster1, newCluster2);
}

void ColourReconnector::_2Bto2BreconnectionFinder(ClusterPtr & c1, ClusterPtr & c2,
    int & bswap1, int & bswap2, double min_displ_sum, string & kind_of_reco) const {
  double tmp_delta;
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      // try swapping particle i of c1 with particle j of c2
      tmp_delta  = _displacementBaryonic(c2->particle(j), c1->particle((i+1)%3), c1->particle((i+2)%3));
      tmp_delta += _displacementBaryonic(c1->particle(i), c2->particle((j+1)%3), c2->particle((j+2)%3));

      if (tmp_delta < min_displ_sum) {
        // if minimal displacement select the 2Bto2B CR option
        min_displ_sum = tmp_delta;
        bswap1 = i;
        bswap2 = j;
        kind_of_reco = "2Bto2B";
      }
    }
  }

}

void ColourReconnector::_BbarBto3MreconnectionFinder(ClusterPtr & c1, ClusterPtr & c2, int & mswap0, int & mswap1, int & mswap2,
    double min_displ_sum, string & kind_of_reco) const {
  double pre_tmp_delta;
  double tmp_delta;
  for (int p1=0; p1 <3; p1++) {
    // make sure not to form a mesonic octet
    if (_isColour8(c1->particle(0), c2->particle(p1))) continue;

    pre_tmp_delta = _displacement(c1->particle(0), c2->particle(p1));
    for (int p2=1; p2<3; p2++) {

      // make sure not to form a mesonic octet
      if (_isColour8(c1->particle(1), c2->particle((p1+p2)%3))) continue;
      if (_isColour8(c1->particle(2), c2->particle(3-p1-((p1+p2)%3)))) continue;

      tmp_delta  = pre_tmp_delta + _displacement(c1->particle(1), c2->particle((p1+p2)%3));
      tmp_delta += 				 _displacement(c1->particle(2), c2->particle(3-p1-((p1+p2)%3)));

      // factor _mesonToBaryonFactor to compare Baryonic an mesonic cluster
      tmp_delta *=_mesonToBaryonFactor;

      if (tmp_delta < min_displ_sum) {
        // if minimal displacement select the 2Bto3M CR option
        min_displ_sum = tmp_delta;
        mswap0 = p1;
        mswap1 = (p1+p2)%3;
        mswap2 = 3-p1-((p1+p2)%3);
        kind_of_reco = "2Bto3M";

      }
    }
  }
}

void ColourReconnector::_BMtoBMreconnectionfinder(ClusterPtr & c1, ClusterPtr & c2, int & swap, double min_displ_sum,
    string & kind_of_reco) const {
  assert(c1->numComponents()==2);
  assert(c2->numComponents()==3);
  double tmp_displ = 0;
  for (int i=0; i<3; i++) {
    // Differ if the second cluster is baryonic or antibaryonic
    if (c2->particle(0)->hasColour()) {
      // c2 is baryonic

      // veto mesonic octets
      if (_isColour8(c2->particle(i), c1->antiColParticle())) continue;

      // factor _mesonToBaryonFactor to compare Baryonic an mesonic cluster
      tmp_displ = _mesonToBaryonFactor * _displacement(c2->particle(i), c1->antiColParticle());
      tmp_displ += _displacementBaryonic(c1->colParticle(), c2->particle((i+1)%3), c2->particle((i+2)%3));
    } else {
      // c2 is antibaryonic

      // veto mesonic octets
      if (_isColour8(c2->particle(i), c1->colParticle())) continue;

      // factor _mesonToBaryonFactor to compare Baryonic an mesonic cluster
      tmp_displ = _mesonToBaryonFactor * _displacement(c2->particle(i), c1->colParticle());
      tmp_displ *= _displacementBaryonic(c1->antiColParticle(), c2->particle((i+1)%3), c2->particle((i+2)%3));
    }
    if (tmp_displ < min_displ_sum) {
      // if minimal displacement select the MBtoMB CR option
      min_displ_sum = tmp_displ;
      swap = i;
      kind_of_reco = "MBtoMB";
    }
  }
  return;
}

void ColourReconnector::_3MtoXreconnectionfinder(std::vector<CluVecIt> & cv, int & swap0, int & swap1,
    int & swap2, double min_displ_sum, string & kind_of_reco) const {
  // case of 3M->BbarB CR
  double _tmp_displ;
  _tmp_displ  = _displacementBaryonic((*cv[0])->colParticle(),     (*cv[1])->colParticle(),     (*cv[2])->colParticle());
  _tmp_displ += _displacementBaryonic((*cv[0])->antiColParticle(), (*cv[1])->antiColParticle(), (*cv[2])->antiColParticle());
  if (_tmp_displ < min_displ_sum) {
    // if minimal displacement select the 3Mto2B CR option
    kind_of_reco = "3Mto2B";
    min_displ_sum = _tmp_displ;
  }
  // case for 3M->3M CR
  /**
   * if 3Mto3M reco probability (_preco3M_3M) is 0 we skip this loop
   * since no 3Mto3M CR shall be performed
   */
  int i,j;
  int i1,i2,i3;
  for (i = 0; _preco3M_3M && i<3; i++) {
    // veto mesonic octets
    if (_isColour8((*cv[0])->colParticle(), (*cv[i])->antiColParticle())) continue;

    // factor _mesonToBaryonFactor to compare baryonic an mesonic cluster
    _tmp_displ = _mesonToBaryonFactor * _displacement((*cv[0])->colParticle(), (*cv[i])->antiColParticle());
    for (j=1; j<3; j++) {
      // i1, i2, i3 are pairwise distinct
      i1=i;
      i2=((j+i)%3);
      if (i1==0 && i2==1) continue;
      i3=(3-i-((j+i)%3));

      // veto mesonic octets
      if (_isColour8((*cv[1])->colParticle(), (*cv[i2])->antiColParticle())) continue;
      if (_isColour8((*cv[2])->colParticle(), (*cv[i3])->antiColParticle())) continue;

      _tmp_displ += _mesonToBaryonFactor * _displacement((*cv[1])->colParticle(), (*cv[i2])->antiColParticle());
      _tmp_displ += _mesonToBaryonFactor * _displacement((*cv[2])->colParticle(), (*cv[i3])->antiColParticle());

      if (_tmp_displ < min_displ_sum) {
        // if minimal displacement select the 3Mto3M CR option
        kind_of_reco = "3Mto3M";
        min_displ_sum = _tmp_displ;
        swap0 = i1;
        swap1 = i2;
        swap2 = i3;
      }
    }
  }
}


pair <int,int> ColourReconnector::_shuffle
(const PVector & q, const PVector & aq, unsigned maxtries) const {

  const size_t nclusters = q.size();
  assert (nclusters > 1);
  assert (aq.size() == nclusters);

  int i, j;
  unsigned tries = 0;
  bool octet;

  do {
    // find two different random integers in the range [0, nclusters)
    i = UseRandom::irnd( nclusters );
    do {
      j = UseRandom::irnd( nclusters );
    } while (i == j);

    // check if one of the two potential clusters would be a colour octet state
    octet = _isColour8( q[i], aq[j] ) || _isColour8( q[j], aq[i] ) ;
    tries++;
  } while (octet && tries < maxtries);

  if (octet) i = j = -1;
  return make_pair(i,j);
}



bool ColourReconnector::_isColour8(tcPPtr p, tcPPtr q) const {
  bool octet = false;

  // make sure we have a triplet and an anti-triplet
  if ( ( p->hasColour() && q->hasAntiColour() ) ||
       ( p->hasAntiColour() && q->hasColour() ) ) {

    // true if p and q are originated from a colour octet
    if ( !p->parents().empty() && !q->parents().empty() ) {
      octet = ( p->parents()[0] == q->parents()[0] ) &&
              ( p->parents()[0]->data().iColour() == PDT::Colour8 );
    }

    // (Final) option: check if same colour8 parent
    // or already found an octet.
    if(_octetOption==0||octet) return octet;

    // (All) option handling more octets
    // by browsing particle history/colour lines.
    tColinePtr cline,aline;

    // Get colourlines form final states.
    if(p->hasColour() && q->hasAntiColour()) {
      cline  = p->    colourLine();
      aline  = q->antiColourLine();
    } else {
      cline  = q->    colourLine();
      aline  = p->antiColourLine();
    }

    // Follow the colourline of p.
    if ( !p->parents().empty() ) {
      tPPtr parent = p->parents()[0];
      while (parent) {
        if(parent->data().iColour() == PDT::Colour8) {
          // Coulour8 particles should have a colour
          // and an anticolour line. Currently the
          // remnant has none of those. Since the children
          // of the remnant are not allowed to emit currently,
          // the colour octet remnant is handled by the return
          // statement above. The assert also catches other
          // colour octets without clines. If the children of
          // a remnant should be allowed to emit, the remnant
          // should get appropriate colour lines and
          // colour states.
          //  See Ticket: #407
          //  assert(parent->colourLine()&&parent->antiColourLine());
          octet = (parent->    colourLine()==cline &&
                   parent->antiColourLine()==aline);
        }
        if(octet||parent->parents().empty()) break;
        parent = parent->parents()[0];
      }
    }
  }

  return octet;
}


void ColourReconnector::persistentOutput(PersistentOStream & os) const {
  os
      << _clreco
      << _algorithm
      << _annealingFactor
      << _annealingSteps
      << _triesPerStepFactor
      << _initTemp
      << _preco
      << _precoBaryonic
      << _preco3M_3M
      << _preco3M_BBbar
      << _precoBbarB_3M
      << _preco2B_2B
      << _precoMB_MB
      << _stepFactor
      << _mesonToBaryonFactor
      << ounit(_maxDistance, femtometer)
      << _octetOption
      << _cr2BeamClusters
      << _debug
      << _junctionMBCR
      ;
}

void ColourReconnector::persistentInput(PersistentIStream & is, int) {
  is
      >> _clreco
      >> _algorithm
      >> _annealingFactor
      >> _annealingSteps
      >> _triesPerStepFactor
      >> _initTemp
      >> _preco
      >> _precoBaryonic
      >> _preco3M_3M
      >> _preco3M_BBbar
      >> _precoBbarB_3M
      >> _preco2B_2B
      >> _precoMB_MB
      >> _stepFactor
      >> _mesonToBaryonFactor
      >> iunit(_maxDistance, femtometer)
      >> _octetOption
      >> _cr2BeamClusters
      >> _debug
      >> _junctionMBCR
      ;
}


void ColourReconnector::Init() {
  static ClassDocumentation<ColourReconnector> documentation
  ("This class is responsible of the colour reconnection.");


  static Switch<ColourReconnector,int> interfaceColourReconnection
  ("ColourReconnection",
   "Colour reconnections",
   &ColourReconnector::_clreco, 0, true, false);
  static SwitchOption interfaceColourReconnectionNo
  (interfaceColourReconnection,
   "No",
   "Colour reconnections off",
   0);
  static SwitchOption interfaceColourReconnectionYes
  (interfaceColourReconnection,
   "Yes",
   "Colour reconnections on",
   1);

  // Algorithm interface
  static Switch<ColourReconnector, int> interfaceAlgorithm
  ("Algorithm",
   "Specifies the colour reconnection algorithm",
   &ColourReconnector::_algorithm, 0, true, false);
  static SwitchOption interfaceAlgorithmPlain
  (interfaceAlgorithm,
   "Plain",
   "Plain colour reconnection as in Herwig 2.5.0",
   0);
  static SwitchOption interfaceAlgorithmStatistical
  (interfaceAlgorithm,
   "Statistical",
   "Statistical colour reconnection using simulated annealing",
   1);
  static SwitchOption interfaceAlgorithmBaryonic
  (interfaceAlgorithm,
   "Baryonic",
   "Baryonic cluster reconnection",
   2);
  static SwitchOption interfaceAlgorithmBaryonicMesonic
  (interfaceAlgorithm,
   "BaryonicMesonic",
   "Baryonic cluster reconnection with reconnections to and from Mesonic Clusters",
   3);



  // Statistical CR Parameters:
  static Parameter<ColourReconnector, double> interfaceMtrpAnnealingFactor
  ("AnnealingFactor",
   "The annealing factor is the ratio of the temperatures in two successive "
   "temperature steps.",
   &ColourReconnector::_annealingFactor, 0.9, 0.0, 1.0,
   false, false, Interface::limited);

  static Parameter<ColourReconnector,unsigned> interfaceMtrpAnnealingSteps
  ("AnnealingSteps",
   "Number of temperature steps in the statistical annealing algorithm",
   &ColourReconnector::_annealingSteps, 50, 1, 10000,
   false, false, Interface::limited);

  static Parameter<ColourReconnector,double> interfaceMtrpTriesPerStepFactor
  ("TriesPerStepFactor",
   "The number of reconnection tries per temperature steps is the number of "
   "clusters times this factor.",
   &ColourReconnector::_triesPerStepFactor, 5.0, 0.0, 100.0,
   false, false, Interface::limited);


  static Parameter<ColourReconnector,double> interfaceMtrpInitialTemp
  ("InitialTemperature",
   "Factor used to determine the initial temperature from the median of the "
   "energy change in a few random rearrangements.",
   &ColourReconnector::_initTemp, 0.1, 0.00001, 100.0,
   false, false, Interface::limited);




  // Plain and Baryonic CR Paramters
  static Parameter<ColourReconnector, double> interfaceRecoProb
  ("ReconnectionProbability",
   "Probability that a found two meson to two meson reconnection possibility is actually accepted (used in Plain & Baryonic)",
   &ColourReconnector::_preco, 0.5, 0.0, 1.0,
   false, false, Interface::limited);

  static Parameter<ColourReconnector,double> interfaceRecoProbBaryonic
  ("ReconnectionProbabilityBaryonic",
   "Probability that a found reconnection possibility is actually accepted (used in Baryonic)",
   &ColourReconnector::_precoBaryonic, 0.5, 0.0, 1.0,
   false, false, Interface::limited);


  // BaryonicMesonic CR Paramters
  static Parameter<ColourReconnector, double> interfaceReconnectionProbability3Mto3M
  ("ReconnectionProbability3Mto3M",
   "Probability that a reconnection candidate is accepted for reconnecting 3M -> 3M\'",
   &ColourReconnector::_preco3M_3M, 0.5, 0.0, 1.0,
   false, false, Interface::limited);
  static Parameter<ColourReconnector, double> interfaceReconnectionProbability3MtoBBbar
  ("ReconnectionProbability3MtoBBbar",
   "Probability that a reconnection candidate is accepted for reconnecting 3M -> B,Bbar",
   &ColourReconnector::_preco3M_BBbar, 0.5, 0.0, 1.0,
   false, false, Interface::limited);
  static Parameter<ColourReconnector, double> interfaceReconnectionProbabilityBbarBto3M
  ("ReconnectionProbabilityBbarBto3M",
   "Probability that a reconnection candidate is accepted for reconnecting B,Bbar -> 3M",
   &ColourReconnector::_precoBbarB_3M, 0.5, 0.0, 1.0,
   false, false, Interface::limited);
  static Parameter<ColourReconnector, double> interfaceReconnectionProbability2Bto2B
  ("ReconnectionProbability2Bto2B",
   "Probability that a reconnection candidate is accepted for reconnecting 2B -> 2B\' or 2Bbar -> 2Bbar\'",
   &ColourReconnector::_preco2B_2B, 0.5, 0.0, 1.0,
   false, false, Interface::limited);
  static Parameter<ColourReconnector, double> interfaceReconnectionProbabilityMBtoMB
  ("ReconnectionProbabilityMBtoMB",
   "Probability that a reconnection candidate is accepted for reconnecting M,B -> M\',B\' or M,Bbar -> M\',Bbar\'",
   &ColourReconnector::_precoMB_MB, 0.5, 0.0, 1.0,
   false, false, Interface::limited);

  static Parameter<ColourReconnector, double> interfaceFactorforStep
  ("StepFactor",
   "Factor for how many reconnection-tries are made in the BaryonicMesonic algorithm",
   &ColourReconnector::_stepFactor, 1.0, 0.11111, 10.,
   false, false, Interface::limited);// at least 3 Clusters -> _stepFactorMin=1/9

  static Parameter<ColourReconnector, double> interfaceMesonToBaryonFactor
  ("MesonToBaryonFactor",
   "Factor for comparing mesonic clusters to baryonic clusters in the displacement if BaryonicMesonic CR model is chosen",
   &ColourReconnector::_mesonToBaryonFactor, 2.0, 0.5, 3.0,
   false, false, Interface::limited);



  // General Parameters and switches
  static Parameter<ColourReconnector, Length> interfaceMaxDistance
  ("MaxDistance",
   "Maximum distance between the clusters at which to consider rearrangement"
   " to avoid colour reconneections of displaced vertices (used in all Algorithms). No unit means femtometer",
   &ColourReconnector::_maxDistance, femtometer, 1000.*femtometer, 0.0*femtometer, 1e100*femtometer,
   false, false, Interface::limited);


  static Switch<ColourReconnector, unsigned int> interfaceOctetTreatment
  ("OctetTreatment",
   "Which octets are not allowed to be reconnected (used in all Algorithms)",
   &ColourReconnector::_octetOption, 0, false, false);
  static SwitchOption interfaceOctetTreatmentFinal
  (interfaceOctetTreatment,
   "Final",
   "Only prevent for the final (usuaslly non-perturbative) g -> q qbar splitting",
   0);
  static SwitchOption interfaceOctetTreatmentAll
  (interfaceOctetTreatment,
   "All",
   "Prevent for all octets",
   1);
  static Switch<ColourReconnector, int> interfaceCR2BeamClusters
  ("CR2BeamClusters",
   "Option for colour reconnecting 2 beam remnant clusters if the number of clusters is 2.",
   &ColourReconnector::_cr2BeamClusters, 0, true, false);
  static SwitchOption interfaceCR2BeamClustersYes
  (interfaceCR2BeamClusters,
   "Yes",
   "If possible CR 2 beam clusters",
   1);
  static SwitchOption interfaceCR2BeamClustersNo
  (interfaceCR2BeamClusters,
   "No",
   "If possible do not CR 2 beam clusters",
   0);
  static Switch<ColourReconnector, int> interfaceJunction
  ("Junction",
   "Option for using Junction-like displacement in rapidity-phi plane to compare baryonic cluster "
   "instead of pairwise distance (for BaryonicMesonic model)",
   &ColourReconnector::_junctionMBCR, 1, true, false);
  static SwitchOption interfaceJunctionYes
  (interfaceJunction,
   "Yes",
   "Using junction-like model instead of pairwise distance model",
   1);
  static SwitchOption interfaceJunctionNo
  (interfaceJunction,
   "No",
   "Using pairwise distance model instead of junction-like model",
   0);

  // Debug
  static Switch<ColourReconnector, int> interfaceDebug
  ("Debug",
   "Make a file with some Information of the BaryonicMesonic Algorithm",
   &ColourReconnector::_debug, 0, true, false);
  static SwitchOption interfaceDebugNo
  (interfaceDebug,
   "No",
   "Debug Information for ColourReconnector Off",
   0);
  static SwitchOption interfaceDebugYes
  (interfaceDebug,
   "Yes",
   "Debug Information for ColourReconnector On",
   1);

}

