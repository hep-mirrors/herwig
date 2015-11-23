// -*- C++ -*-
//
// ColourReconnector.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ColourReconnector class.
//

#include "ColourReconnector.h"
#include "Cluster.h"
#include "Herwig/Utilities/Maths.h"
#include <ThePEG/Interface/Switch.h>
#include "ThePEG/Interface/Parameter.h"
#include <ThePEG/Persistency/PersistentOStream.h>
#include <ThePEG/Persistency/PersistentIStream.h>
#include <ThePEG/Repository/UseRandom.h>
#include <algorithm>
#include <ThePEG/Utilities/DescribeClass.h>
#include <ThePEG/Repository/EventGenerator.h>


using namespace Herwig;

typedef ClusterVector::iterator CluVecIt;

DescribeClass<ColourReconnector,Interfaced>
describeColourReconnector("Herwig::ColourReconnector","");

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
    case 0: _doRecoPlain(clusters);
            break;
    case 1: _doRecoStatistical(clusters);
            break;
  }

  return;
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
    if (isColour8(p, q)) return true;
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


CluVecIt ColourReconnector::_findRecoPartner(CluVecIt cl,
                                             ClusterVector & cv) const {

  CluVecIt candidate = cl;
  Energy minMass = 1*TeV;
  for (CluVecIt cit=cv.begin(); cit != cv.end(); ++cit) {

    // don't even look at original cluster
    if(cit==cl) continue;

    // don't allow colour octet clusters
    if ( isColour8( (*cl)->colParticle(),
	            (*cit)->antiColParticle() )  ||
         isColour8( (*cit)->colParticle(),
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


pair <ClusterPtr,ClusterPtr>
ColourReconnector::_reconnect(ClusterPtr c1, ClusterPtr c2) const {

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
    octet = isColour8( q[i], aq[j] ) || isColour8( q[j], aq[i] ) ;
    tries++;
  } while (octet && tries < maxtries);

  if (octet) i = j = -1;
  return make_pair(i,j);
}


bool ColourReconnector::isColour8(cPPtr p, cPPtr q) {
  bool octet = false;

  // make sure we have a triplet and an anti-triplet
  if ( ( p->hasColour() && q->hasAntiColour() ) ||
       ( p->hasAntiColour() && q->hasColour() ) ) {
    if ( !p->parents().empty() && !q->parents().empty() ) {
      // true if p and q are originated from a colour octet
      octet = ( p->parents()[0] == q->parents()[0] ) &&
              ( p->parents()[0]->data().iColour() == PDT::Colour8 );
    }
  }

  return octet;
}


void ColourReconnector::persistentOutput(PersistentOStream & os) const {
  os << _clreco << _preco << _algorithm << _initTemp << _annealingFactor
     << _annealingSteps << _triesPerStepFactor << ounit(_maxDistance,femtometer);
}

void ColourReconnector::persistentInput(PersistentIStream & is, int) {
  is >> _clreco >> _preco >> _algorithm >> _initTemp >> _annealingFactor
     >> _annealingSteps >> _triesPerStepFactor >> iunit(_maxDistance,femtometer);
}


void ColourReconnector::Init() {

  static ClassDocumentation<ColourReconnector> documentation
    ("This class is responsible of the colour reconnection.");


  static Switch<ColourReconnector,int> interfaceColourReconnection
    ("ColourReconnection",
     "Colour reconnections",
     &ColourReconnector::_clreco, 0, true, false);
  static SwitchOption interfaceColourReconnectionOff
    (interfaceColourReconnection,
     "No",
     "Colour reconnections off",
     0);
  static SwitchOption interfaceColourReconnectionOn
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


  static Parameter<ColourReconnector,Length> interfaceMaxDistance
    ("MaxDistance",
     "Maximum distance between the clusters at which to consider rearrangement"
     " to avoid colour reconneections of displaced vertices",
     &ColourReconnector::_maxDistance, femtometer, 1000.*femtometer, 0.0*femtometer, 1e100*femtometer,
     false, false, Interface::limited);

}
