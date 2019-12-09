// -*- C++ -*-
//
// Cluster.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Cluster class.
//

#include "Cluster.h"
#include <ThePEG/Repository/UseRandom.h>
#include <ThePEG/Repository/CurrentGenerator.h>
#include <ThePEG/PDT/ParticleData.h>
#include "ClusterHadronizationHandler.h"

using namespace Herwig;

PPtr Cluster::clone() const {
  return new_ptr(*this);
}
  
PPtr Cluster::fullclone() const {
  return clone();
}

ClassDescription<Cluster> Cluster::initCluster;

Cluster::Cluster() 
  : Particle(), 
    _isAvailable(true),
    _hasReshuffled(false),
    _component(),
    _original(),
    _isBeamRemnant(),
    _isPerturbative(),
    _numComp(0),
    _id(0) {}

void Cluster::persistentOutput(PersistentOStream & os) const {
  os << _isAvailable << _hasReshuffled << _component << _original
     << _isBeamRemnant << _isPerturbative << _numComp << _id;
}

void Cluster::persistentInput(PersistentIStream & is, int) {
  is >> _isAvailable >> _hasReshuffled >> _component >> _original
     >> _isBeamRemnant >> _isPerturbative >> _numComp >> _id;
}



Cluster::Cluster(tPPtr p1, tPPtr p2, tPPtr p3)
  : Particle(CurrentGenerator::current().
	     getParticleData(long(ParticleID::Cluster))),
    _isAvailable(true), _hasReshuffled(false)
{
  if(!dataPtr()) {
    cerr << "Cluster Particle Data not defined. Cannot complete Hadronization "
	 << "without ParticleData for id " << ParticleID::Cluster << '\n';
  }
  _component.push_back(new_ptr(Particle(*p1))); 
  _component.push_back(new_ptr(Particle(*p2))); 
  if(p3) _component.push_back(new_ptr(Particle(*p3)));
  _original.push_back(p1); 
  _original.push_back(p2); 
  if(p3) _original.push_back(p3);

  _isPerturbative.push_back(initPerturbative(p1));
  _isPerturbative.push_back(initPerturbative(p2));
  if(p3) _isPerturbative.push_back(initPerturbative(p3));
  else _isPerturbative.push_back(false);
  for(int i = 0; i<3; i++) _isBeamRemnant.push_back(false);

  if(p3) {
    _numComp = 3;
    _id = 100*abs(p1->id()) + 10*abs(p2->id()) + abs(p3->id());
  } else {
    _numComp = 2;
    int i1,i2;
    if(p2->id() > 10) {
      i1 = abs(p2->id()/100);
      i2 = abs(p1->id());
    }
    else if(p1->id() > 10) {
      i1 = abs(p1->id()/100);
      i2 = abs(p2->id());
    }
    else {
      i1 = abs(p1->id());
      i2 = abs(p2->id());
    }
    if(i1>i2) swap (i1,i2);
    _id = 10*i1+i2;
  }
  // calculate the momentum
  calculateP();
  // calculate the position
  // Only in the case of two components we have a definition of cluster
  // position in terms of the two components.
  if ( _numComp != 2 ) {
    // take the average
    setVertex(LorentzPoint());
  }
  else {
    setVertex(calculateX(_component[0],_component[1]));
  }
}    

Cluster::Cluster(tcEventPDPtr x) 
  : Particle(x),
    _isAvailable(false),
    _hasReshuffled(false),
    _component(),
    _original(),
    _isBeamRemnant(),
    _isPerturbative(),
    _numComp(0),
    _id(0) {}

Cluster::Cluster(const Particle &x) 
  : Particle(x),
    _isAvailable(false),
    _hasReshuffled(false),
    _component(),
    _original(),
    _isBeamRemnant(),
    _isPerturbative(),
    _numComp(0),
    _id(0) {}

Energy Cluster::sumConstituentMasses() const 
{
  if(_numComp == 3) { 
    return _component[0]->mass() + 
           _component[1]->mass() +
           _component[2]->mass();
  } else if(_numComp == 2) 
    return _component[0]->mass() + _component[1]->mass();
  else return ZERO;
}


void Cluster::calculateP() {
  Lorentz5Momentum m;
  for(int i = 0; i<_numComp; i++) 
    m += _component[i]->momentum();
  m.rescaleMass();
  set5Momentum(m);
}


LorentzPoint Cluster::calculateX(tPPtr q1, tPPtr q2) {
  // Get the needed parameters. 
  Energy2 vmin2 
    = ClusterHadronizationHandler::currentHandler()->minVirtuality2();
  Length dmax 
    = ClusterHadronizationHandler::currentHandler()->maxDisplacement();
  
  // Get the positions and displacements of the two components (Lab frame).
  LorentzPoint pos1 = q1->vertex();
  Lorentz5Momentum p1 = q1->momentum();
  LorentzDistance displace1 = -log( UseRandom::rnd() ) * 
    hbarc * p1 * (1 / sqrt(sqr(p1.m2() - p1.mass2()) + sqr(vmin2)));
  if ( abs( displace1.m() ) > dmax ) {
    displace1 *= dmax / abs( displace1.m() );
  }
  LorentzPoint pos2 = q2->vertex();
  Lorentz5Momentum p2 = q2->momentum();
  LorentzDistance displace2 = -log( UseRandom::rnd() ) * 
    hbarc * p2 * (1 / sqrt(sqr(p2.m2() - p2.mass2()) + sqr(vmin2)));
  if ( abs( displace2.m() ) > dmax ) {
    displace2 *= dmax / abs( displace2.m() );
  }
  double s1 = 0.0, s2 = 0.0;
  Lorentz5Momentum pcl = p1 + p2;
  if ( abs( pcl.vect().dot( displace1.vect() ) ) > 1.0e-20*MeV*mm  &&
       abs( pcl.vect().dot( displace2.vect() ) ) > 1.0e-20*MeV*mm  ) {
    // The displacement with the smallest projection along pcl.vect()
    // is scaled up such that both displacements have equal projections
    // along pcl.vect().
    double ratio = ( abs( pcl.vect().dot( displace1.vect() ) ) / 
		     abs( pcl.vect().dot( displace2.vect() ) ) );
    if ( pcl.vect().dot(displace1.vect()) * 
	 pcl.vect().dot(displace2.vect())  <  0.0*sqr(MeV*mm) ) {
      ratio *= -1;
    }
    if ( abs( ratio ) > 1.0 ) {
      displace2 *= ratio;
    } else {
      displace1 *= ratio;
    }
    // Now determine the s1 and s2 values.
    double s1minusS2 = ( pcl.vect().dot( pos2.vect() - pos1.vect() ) / 
			 pcl.vect().dot( displace1.vect() ) );
    if ( s1minusS2 < 0 ) {
      s1 = 1.0;
      s2 = s1 - s1minusS2;
    } else if ( s1minusS2 > 0 ) {
      s2 = 1;
      s1 = s2 + s1minusS2;
    }
  }
  // Now, finally, determine the cluster position
  LorentzPoint position = 0.5 * (pos1 + pos2 + s1*displace1 + s2*displace2);
  // set the decay vertex of the two particles via the lifeLength
  q1->setLifeLength(position-q1->vertex());
  q2->setLifeLength(position-q2->vertex());
  // return the answer
  return position;
}

bool Cluster::isBeamCluster() const {
  for(int i = 0; i<_numComp; i++)
    if(_isBeamRemnant[i]) return true;
  return false;
}


void Cluster::isBeamCluster(tPPtr part) {
  for(int i = 0; i<_numComp; i++) {
    if(_original[i] == part) {
      _isBeamRemnant[i] = true;
      break;
    }
  }
}

bool Cluster::isStatusFinal() const {
  int s = children().size();
  for(unsigned int i = 0; i<children().size(); i++) 
    if(children()[i]->PDGName() == "Cluster") s--;
  return ( s > 0);
}

tPPtr Cluster::particle(int i) const { 
  return (i < _numComp) ? _component[i] : PPtr(); 
}

tPPtr Cluster::colParticle(bool anti) const {
  if ( _numComp != 2 ) return PPtr();
  if ( _original[0]->hasColour(anti) ) return _original[0];
  else if ( _original[1]->hasColour(anti) ) return _original[1];
  else return PPtr();
}

tPPtr Cluster::antiColParticle() const {
  return colParticle(true);
}

bool Cluster::isPerturbative(int i) const { 
  return _isPerturbative[i]; 
}

bool Cluster::isBeamRemnant(int i) const { 
  return _isBeamRemnant[i]; 
}

void Cluster::setBeamRemnant(int i, bool b) {
  if(i < _numComp)
    _isBeamRemnant[i] = b;
}

bool Cluster::initPerturbative(tPPtr p)
{ 
  Energy mg 
    = CurrentGenerator::current().getParticleData(ParticleID::g)->constituentMass();
  return p->scale() > sqr(mg);
}
