// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Cluster class.
//

#include "Cluster.h"
#include <ThePEG/Repository/UseRandom.h>
#include <ThePEG/Repository/CurrentGenerator.h>
#include <iostream>
#include <iomanip>
#include <cctype>
#include <ThePEG/EventRecord/Step.h>
#include <ThePEG/EventRecord/Event.h>
#include <ThePEG/EventRecord/ColourLine.h>
#include <ThePEG/Utilities/Rebinder.h>
#include <ThePEG/Config/algorithm.h>
#include <ThePEG/EventRecord/ParticleTraits.h>
#include <ThePEG/Persistency/PersistentOStream.h>
#include <ThePEG/Persistency/PersistentIStream.h>
#include <ThePEG/PDT/DecayMode.h>

using namespace Herwig;
// using namespace ThePEG;

GlobParamPtr Cluster::_globalParameters = GlobParamPtr();

ClassDescription<Cluster> Cluster::initCluster;

Cluster::Cluster() 
  : Particle(CurrentGenerator::current().getParticleData(ExtraParticleID::Cluster)), 
    _isAvailable(true),
    _reshufflingPartner(),
    _component(),
    _original(),
    _isBeamRemnant(),
    _isPerturbative(),
    _numComp(0),
    _id(0) {}  

Cluster::Cluster(tPPtr p1, tPPtr p2, tPPtr p3)
  : Particle(CurrentGenerator::current().getParticleData(ExtraParticleID::Cluster)), _isAvailable(true) 
{
  if(!dataPtr()) {
    cerr << "Cluster Particle Data not defined. Cannot complete Hadronization "
	 << "without ParticleData for id " << ExtraParticleID::Cluster << '\n';
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
    if(p2->id() > 10)
      _id = 10*abs(p2->id()/100) + abs(p1->id());
    else if(p1->id() > 10)
      _id = 10*abs(p1->id()/100) + abs(p2->id());
    else
      _id = 10*abs(p1->id()) + abs(p2->id());
  }
  calculateP();
  calculateX();

}    

Cluster::Cluster(tcEventPDPtr x) 
  : Particle(x),
    _isAvailable(false),
    _reshufflingPartner(),
    _component(),
    _original(),
    _isBeamRemnant(),
    _isPerturbative(),
    _numComp(0),
    _id(0) {}
    
Cluster::~Cluster() { /*cout << "Destroying Cluster\n"; _component.clear(); */}    

Cluster::Cluster(const Cluster &x) 
  : Particle(x),
    _isAvailable(x._isAvailable),  
    _reshufflingPartner(x._reshufflingPartner),
    _component(x._component),
    _original(x._original),
    _isBeamRemnant(x._isBeamRemnant),
    _isPerturbative(x._isPerturbative),
    _numComp(x._numComp),
    _id(x._id) {}

Cluster::Cluster(const Particle &x) 
  : Particle(x),
    _isAvailable(false),
    _reshufflingPartner(),
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
  else return 0.0;
}


void Cluster::calculateP() {
  Lorentz5Momentum m;
  for(int i = 0; i<_numComp; i++) 
    m += _component[i]->momentum();
  m.rescaleMass();
  set5Momentum(m);
}


void Cluster::calculateX() {
  if ( _numComp != 2 ) {
    // Only in the case of two components we have a definition of cluster
    // position in terms of the two components.
    setVertex(LorentzPoint());
  } else {
    // Get the needed global parameters. 
    double theConversionFactorGeVtoMillimeter = 0.0;
    Energy2 vmin2 = Energy2();
    Length dmax = Length();
    if ( _globalParameters ) {
      theConversionFactorGeVtoMillimeter = 
    	_globalParameters->conversionFactorGeVtoMillimeter();
      vmin2 = _globalParameters->minVirtuality2();
      dmax = _globalParameters->maxDisplacement();
    }
    // Get the positions and displacements of the two components (Lab frame).
    LorentzPoint pos1 = _component[0]->vertex();
    Lorentz5Momentum p1 = _component[0]->momentum();
    LorentzDistance displace1 = 
      -log( UseRandom::rnd() ) * theConversionFactorGeVtoMillimeter 
      * (p1/GeV) * (1.0/sqrt(sqr(p1.m2()/GeV2 - p1.mass2()/GeV2) + 
			     sqr(vmin2/GeV2)));
    if ( fabs( displace1.mag() ) > dmax ) {
      displace1 *= dmax / fabs( displace1.mag() );
    }
    LorentzPoint pos2 = _component[1]->vertex();
    Lorentz5Momentum p2 = _component[1]->momentum();
    LorentzDistance displace2 = 
      -log( UseRandom::rnd() ) * theConversionFactorGeVtoMillimeter 
      * (p2/GeV) * (1.0/sqrt(sqr(p2.m2()/GeV2 - p2.mass2()/GeV2) + 
			     sqr(vmin2 /GeV2)));
    if ( fabs( displace2.mag() ) > dmax ) {
      displace2 *= dmax / fabs( displace2.mag() );
    }
    double s1 = 0.0, s2 = 0.0;
    Lorentz5Momentum pcl = p1 + p2;
    if ( fabs( pcl.vect().dot( displace1.vect() ) ) > 1.0e-20  &&
	 fabs( pcl.vect().dot( displace2.vect() ) ) > 1.0e-20  ) {
      // The displacement with the smallest projection along pcl.vect()
      // is scaled up such that both displacements have equal projections
      // along pcl.vect().
      double ratio = ( fabs( pcl.vect().dot( displace1.vect() ) ) / 
		       fabs( pcl.vect().dot( displace2.vect() ) ) );
      if ( pcl.vect().dot(displace1.vect()) * 
	   pcl.vect().dot(displace2.vect())  <  0.0 ) {
	ratio *= -1;
      }
      if ( fabs( ratio ) > 1.0 ) {
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
    setVertex(0.5 * (pos1 + pos2 + s1*displace1 + s2*displace2));
  } // end else part of if ( _collecCompPtr.size() != 2 )   
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
  if(i < _numComp)
    return _component[i]; 
  return _component[2];
}

bool Cluster::isPerturbative(int i) const { 
  return _isPerturbative[i]; 
}

void Cluster::setPerturbative(int i, bool b) { 
  if(i < _numComp)
    _isPerturbative[i] = b; 
}

bool Cluster::isBeamRemnant(int i) const { 
  return _isBeamRemnant[i]; 
}

void Cluster::setBeamRemnant(int i, bool b) {
  if(i < _numComp)
    _isBeamRemnant[i] = b;
}

PPtr Cluster::clone() const {
  return dynamic_ptr_cast<PPtr>(ptr_new<ClusterPtr>(*this));
}

PPtr Cluster::fullclone() const {
  return clone();
}

bool Cluster::initPerturbative(tPPtr p) {
  Energy Q0 = CurrentGenerator::current().getParticleData(ParticleID::g)->constituentMass();
  if(p->scale() > Q0*Q0) return true;
  return false;
}
