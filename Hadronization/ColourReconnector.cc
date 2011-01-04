// -*- C++ -*-
//
// ColourReconnector.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ColourReconnector class.
//

#include "ColourReconnector.h"
#include "Cluster.h"

#include <ThePEG/Interface/ClassDocumentation.h>
#include <ThePEG/Interface/Switch.h>
#include "ThePEG/Interface/Parameter.h"
#include <ThePEG/Persistency/PersistentOStream.h>
#include <ThePEG/Persistency/PersistentIStream.h>
#include <ThePEG/Repository/UseRandom.h>

#include <algorithm>

using namespace Herwig;


void ColourReconnector::rearrange(EventHandler & ch, 
				  ClusterVector & clusters) {
  if (_clreco == 0) return;

  ClusterVector newClusters = clusters;

  // try avoiding systematic errors by randomizing the reconnection order
  long (*p_irnd)(long) = UseRandom::irnd;
  random_shuffle( newClusters.begin(), newClusters.end(), p_irnd ); 

  // iterate over all clusters
  for (ClusterVector::iterator currentCl=newClusters.begin();
      currentCl != newClusters.end(); currentCl++) {

    // find cluster, where the new clusters are lighter
    ClusterPtr candidate = _findRecoPartner(*currentCl, newClusters);

    // skip this cluster if no possible reshuffling partner can be found
    if (candidate == *currentCl) {
      continue;
    }

    // accept the reconnection with probability _preco.
    if (UseRandom::rnd() < _preco) {

      pair <ClusterPtr,ClusterPtr> reconnected =
	_reconnect(*currentCl, candidate);

      // reset children-parents relations: 
      _get3(*currentCl)->abandonChild(*currentCl);
      _get3(*currentCl)->addChild(reconnected.first);
      _getAnti3(*currentCl)->abandonChild(*currentCl);
      _getAnti3(*currentCl)->addChild(reconnected.second);
      _get3(candidate)->abandonChild(candidate);
      _get3(candidate)->addChild(reconnected.second);
      _getAnti3(candidate)->abandonChild(candidate);
      _getAnti3(candidate)->addChild(reconnected.first);

      // Replace the clusters in the ClusterVector. The order of the partons in
      // the cluster vector carrying colour (i.e. not anticolour) is preserved.

      // replace *currentCl by reconnected.first
      replace( newClusters.begin(), newClusters.end(), *currentCl,
	  reconnected.first ); 

      // replace candidate by reconnected.second
      replace( newClusters.begin(), newClusters.end(), candidate,
	  reconnected.second ); 
    }
  }

  // override pristine clusters with the new ones
  swap(clusters,newClusters);
}


ClusterPtr ColourReconnector::_findRecoPartner(ClusterPtr cl,
    ClusterVector cv) const {

  ClusterPtr candidate = cl;
  Energy minMass = 1*TeV;
  for (ClusterVector::const_iterator cit=cv.begin(); cit != cv.end(); ++cit) {

    // don't allow colour octet clusters
    if ( _isColour8( _get3(cl),
	             _getAnti3(*cit) )  ||
         _isColour8( _get3(*cit),
	             _getAnti3(cl) ) ) {
      continue;
    }

    // momenta of the old clusters
    Lorentz5Momentum p1 = _get3(cl)->momentum() + 
                          _getAnti3(cl)->momentum();
    Lorentz5Momentum p2 = _get3(*cit)->momentum() + 
                          _getAnti3(*cit)->momentum();

    // momenta of the new clusters
    Lorentz5Momentum p3 = _get3(cl)->momentum() + 
                          _getAnti3(*cit)->momentum();
    Lorentz5Momentum p4 = _get3(*cit)->momentum() + 
                          _getAnti3(cl)->momentum();

    Energy oldMass = abs( p1.m() ) + abs( p2.m() );
    Energy newMass = abs( p3.m() ) + abs( p4.m() );

    if ( newMass < oldMass && newMass < minMass ) {
      minMass = newMass;
      candidate = *cit;
    }
  }
  return candidate;
}


tPPtr ColourReconnector::_getAnti3(ClusterPtr cl) const {
  if ( cl->_original[0]->hasAntiColour() ) {
    return cl->_original[0];
  } else if ( cl->_original[1]->hasAntiColour() ) {
    return cl->_original[1];
  } else {
    return tPPtr();
  }
}


tPPtr ColourReconnector::_get3(ClusterPtr cl) const {
  if ( cl->_original[0]->hasColour() ) {
    return cl->_original[0];
  } else if ( cl->_original[1]->hasColour() ) {
    return cl->_original[1];
  } else {
    return tPPtr();
  }
}


bool ColourReconnector::_isColour8(tPPtr p1, tPPtr p2) const {
  bool octet = false;
  // make sure we have a triplet and an anti-triplet
  if ( ( p1->hasColour() && p2->hasAntiColour() ) ||
       ( p1->hasAntiColour() && p2->hasColour() ) ) {
    if ( p1->parents().size()>0 && p2->parents().size()>0 ) {
      // true if p1 and p2 are originated from the same gluon
      octet = ( p1->parents()[0] == p2->parents()[0] ) &&
              ( p1->parents()[0]->data().iColour() == PDT::Colour8 );
    }
  }
  return octet;
}


pair <ClusterPtr,ClusterPtr> ColourReconnector::_reconnect(ClusterPtr c1,
    ClusterPtr c2) const {
  // choose the other possibility to form two clusters, that have no netto
  // colour
  ClusterPtr newCluster1 = new_ptr( Cluster( _get3(c1),
                                             _getAnti3(c2) ) );
  ClusterPtr newCluster2 = new_ptr( Cluster( _get3(c2),
                                             _getAnti3(c1) ) );
  return pair <ClusterPtr,ClusterPtr> (newCluster1,newCluster2);
}


void ColourReconnector::persistentOutput(PersistentOStream & os) const {
  os << _clreco << _preco;
}

void ColourReconnector::persistentInput(PersistentIStream & is, int) {
  is >> _clreco >> _preco;
}

ClassDescription<ColourReconnector> ColourReconnector::initColourReconnector;
// Definition of the static class description member.


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
  
  static Parameter<ColourReconnector,double> interfaceRecoProb 
    ("ReconnectionProbability",
     "Probability that a found reconnection possibility is actually accepted",
     &ColourReconnector::_preco, 0.5, 0.0, 1.0,
     false, false, Interface::limited);

}
