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
    case 0: _doRecoPlain(clusters); break;
    case 1: _doRecoStatistical(clusters); break;
    case 2: _doRecoBaryonic(clusters); break;
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
    for(int i = 0; i<(*cit)->numComponents(); i++) {
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
    random_shuffle( newcv.begin(), newcv.end(), p_irnd ); 

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
      if ( !isBaryonicCandidate && candidate==cit ) 
        continue;

      if ( isBaryonicCandidate 
           && UseRandom::rnd() < _precoBaryonic ) {
        deleted.push_back(*baryonic2);
     
        // Function that does the reconnection from 3 -> 2 clusters
        ClusterPtr b1, b2;
        _makeBaryonicClusters(*cit,*baryonic1,*baryonic2, b1, b2);

        *cit = b1;
        *baryonic1 = b2;

        // Baryonic2 is easily skipped in the next loop 
      }
  
      // Normal 2->2 Colour reconnection
      if ( !isBaryonicCandidate 
           && UseRandom::rnd() < _preco ) {
        auto reconnected = _reconnectBaryonic(*cit, *candidate);
        *cit = reconnected.first;
        *candidate = reconnected.second;
      }
    }

    // create a new vector of clusters except for the ones which are "deleted" during
    // baryonic reconnection
    ClusterVector clustervector;
    for ( const auto & cluster : newcv )
      if ( find(deleted.begin(),
                deleted.end(), cluster) == deleted.end() )
        clustervector.push_back(cluster);

    swap(cv,clustervector);



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
      if ( ( (**cl).vertex()-(**cit).vertex() ).m() >_maxDistance )
        continue;
  
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
          if(maxsum != 0){
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

    if(bcand == true){
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
    if(cit==cl) continue;

    // don't allow colour octet clusters
    if ( _isColour8( (*cl)->colParticle(),
              (*cit)->antiColParticle() )  ||
         _isColour8( (*cit)->colParticle(),
              (*cl)->antiColParticle() ) ) {
      continue;
    }

    // stop it putting beam remnants together
    if((*cl)->isBeamCluster() && (*cit)->isBeamCluster()) continue;

    // stop it putting far apart clusters together
    if(((**cl).vertex()-(**cit).vertex()).m()>_maxDistance) continue;

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
                ClusterPtr &newcluster2) const{

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
ColourReconnector::_reconnect(ClusterPtr &c1, ClusterPtr &c2) const {

  // choose the other possibility to form two clusters from the given
  // constituents

  assert(c1->numComponents()==2);
  assert(c2->numComponents()==2);
  int c1_col(-1),c1_anti(-1),c2_col(-1),c2_anti(-1);
  for(unsigned int ix=0;ix<2;++ix) {
    if     (c1->particle(ix)->hasColour(false)) c1_col  = ix;
    else if(c1->particle(ix)->hasColour(true )) c1_anti = ix;
    if     (c2->particle(ix)->hasColour(false)) c2_col  = ix;
    else if(c2->particle(ix)->hasColour(true )) c2_anti = ix;
  }
  assert(c1_col>=0&&c2_col>=0&&c1_anti>=0&&c2_anti>=0);

  ClusterPtr newCluster1
    = new_ptr( Cluster( c1->colParticle(), c2->antiColParticle() ) );

  newCluster1->setVertex(0.5*( c1->colParticle()->vertex() + 
                               c2->antiColParticle()->vertex() ));

  if(c1->isBeamRemnant(c1_col )) newCluster1->setBeamRemnant(0,true);
  if(c2->isBeamRemnant(c2_anti)) newCluster1->setBeamRemnant(1,true);

  ClusterPtr newCluster2
    = new_ptr( Cluster( c2->colParticle(), c1->antiColParticle() ) );

  newCluster2->setVertex(0.5*( c2->colParticle()->vertex() + 
                               c1->antiColParticle()->vertex() ));

  if(c2->isBeamRemnant(c2_col )) newCluster2->setBeamRemnant(0,true);
  if(c1->isBeamRemnant(c1_anti)) newCluster2->setBeamRemnant(1,true);

  return pair <ClusterPtr,ClusterPtr> (newCluster1, newCluster2);
}





pair <ClusterPtr,ClusterPtr>
ColourReconnector::_reconnectBaryonic(ClusterPtr &c1, ClusterPtr &c2) const {

  // choose the other possibility to form two clusters from the given
  // constituents

  assert(c1->numComponents()==2);
  assert(c2->numComponents()==2);
  int c1_col(-1),c1_anti(-1),c2_col(-1),c2_anti(-1);
  for(unsigned int ix=0;ix<2;++ix) {
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

  newCluster1->setVertex(0.5*( c1->colParticle()->vertex() + 
             c2->antiColParticle()->vertex() ));

  if(c1->isBeamRemnant(c1_col )) newCluster1->setBeamRemnant(0,true);
  if(c2->isBeamRemnant(c2_anti)) newCluster1->setBeamRemnant(1,true);

  c1->antiColParticle()->abandonChild(c1);
  c2->colParticle()->abandonChild(c2);

  ClusterPtr newCluster2
    = new_ptr( Cluster( c2->colParticle(), c1->antiColParticle() ) );

    c1->antiColParticle()->addChild(newCluster2);
    c2->colParticle()->addChild(newCluster2);

  newCluster2->setVertex(0.5*( c2->colParticle()->vertex() + 
             c1->antiColParticle()->vertex() ));

  if(c2->isBeamRemnant(c2_col )) newCluster2->setBeamRemnant(0,true);
  if(c1->isBeamRemnant(c1_anti)) newCluster2->setBeamRemnant(1,true);

  return pair <ClusterPtr,ClusterPtr> (newCluster1, newCluster2);
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
    do { j = UseRandom::irnd( nclusters ); } while (i == j);

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
    }
    else {
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
  os << _clreco << _preco << _precoBaryonic << _algorithm << _initTemp << _annealingFactor
     << _annealingSteps << _triesPerStepFactor << ounit(_maxDistance,femtometer)
     << _octetOption;
}

void ColourReconnector::persistentInput(PersistentIStream & is, int) {
  is >> _clreco >> _preco >> _precoBaryonic >> _algorithm >> _initTemp >> _annealingFactor
     >> _annealingSteps >> _triesPerStepFactor >> iunit(_maxDistance,femtometer)
     >> _octetOption;
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
  

  static Parameter<ColourReconnector,double> interfaceMtrpAnnealingFactor
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


  static Parameter<ColourReconnector,double> interfaceRecoProb 
    ("ReconnectionProbability",
     "Probability that a found reconnection possibility is actually accepted",
     &ColourReconnector::_preco, 0.5, 0.0, 1.0,
     false, false, Interface::limited);

  static Parameter<ColourReconnector,double> interfaceRecoProbBaryonic
    ("ReconnectionProbabilityBaryonic",
     "Probability that a found reconnection possibility is actually accepted",
     &ColourReconnector::_precoBaryonic, 0.5, 0.0, 1.0,
     false, false, Interface::limited);


  static Switch<ColourReconnector,int> interfaceAlgorithm
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

  static Parameter<ColourReconnector,Length> interfaceMaxDistance
    ("MaxDistance",
     "Maximum distance between the clusters at which to consider rearrangement"
     " to avoid colour reconneections of displaced vertices",
     &ColourReconnector::_maxDistance, femtometer, 1000.*femtometer, 0.0*femtometer, 1e100*femtometer,
     false, false, Interface::limited);


  static Switch<ColourReconnector,unsigned int> interfaceOctetTreatment
    ("OctetTreatment",
     "Which octets are not allowed to be reconnected",
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

}

