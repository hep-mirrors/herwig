// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Cluster class.
//

#include "Cluster.h"
#include "Component.h"
#include "Pythia7/Repository/UseRandom.h"


using namespace Herwig;
// using namespace Pythia7;
 

Ptr<GlobalParameters>::pointer Cluster::_pointerGlobalParameters = Ptr<GlobalParameters>::pointer();


Cluster::Cluster(tPPtr part1, tPPtr part2, tPPtr part3) : 
  _isAvailable(true), 
  _collecCompPtr(2) 
{
  _collecCompPtr[0] = new_ptr( Component(part1) );
  _collecCompPtr[1] = new_ptr( Component(part2) );
  if (part3) _collecCompPtr.push_back( new_ptr( Component(part3) ) ); 
  calculateP();
  calculateX();
}    

    
Cluster::Cluster(tPPtr part1, const long id2) :
  _isAvailable(true), _collecCompPtr(2) 
{
  _collecCompPtr[0] = new_ptr( Component(part1) );  
  _collecCompPtr[1] = new_ptr( Component(id2) );  
}
    

Cluster::Cluster(const long id1, const long id2) :
  _isAvailable(true), _collecCompPtr(2) 
{
  _collecCompPtr[0] = new_ptr( Component(id1) );  
  _collecCompPtr[1] = new_ptr( Component(id2) );  
}
    

Energy Cluster::sumConstituentMasses() const {
  Energy sum = 0.0;
  for (CollecCompPtr::const_iterator cit = _collecCompPtr.begin();
       cit != _collecCompPtr.end(); ++cit){
    sum += (*cit)->mass();
  }          
  return sum; 
}


void Cluster::calculateP() {
  _momentum = Lorentz5Momentum();
  for (CollecCompPtr::const_iterator cit = _collecCompPtr.begin();
       cit != _collecCompPtr.end(); ++cit){
    _momentum += (*cit)->momentum();
  }
  _momentum.rescaleMass();
}


void Cluster::calculateX() {
  if ( _collecCompPtr.size() != 2 ) {
    // Only in the case of two components we have a definition of cluster
    // position in terms of the two components.
    _position = LorentzPoint();
  } else {
  
    // Get the needed global parameters. 
    double theConversionFactorGeVtoMillimeter = 0.0;
    Energy2 vmin2 = Energy2();
    Length dmax = Length();
    if ( _pointerGlobalParameters ) {
      theConversionFactorGeVtoMillimeter = 
    	_pointerGlobalParameters->conversionFactorGeVtoMillimeter();
      vmin2 = _pointerGlobalParameters->minVirtuality2();
      dmax = _pointerGlobalParameters->maxDisplacement();
    }
    // Get the positions and displacements of the two components (in Lab frame).
    CollecCompPtr::const_iterator cit = _collecCompPtr.begin();
    LorentzPoint pos1 = (*cit)->position();
    Lorentz5Momentum p1 = (*cit)->momentum();
    LorentzDistance displace1 = - log( UseRandom::rnd() ) * theConversionFactorGeVtoMillimeter 
      * ( p1 / GeV ) * ( 1.0 / sqrt( sqr( (*cit)->momentum().m2() / GeV2 - 
					  (*cit)->momentum().mass2() / GeV2 ) + 
				     sqr( vmin2 / GeV2 ) ) );
    if ( fabs( displace1.mag() ) > dmax ) {
      displace1 *= dmax / fabs( displace1.mag() );
      cout << "Cluster::calculateX ***extreme debugging*** : MAX DISPLACEMENT " << endl;
    }
    ++cit;
    LorentzPoint pos2 = (*cit)->position();
    Lorentz5Momentum p2 = (*cit)->momentum();
    LorentzDistance displace2 = - log( UseRandom::rnd() ) * theConversionFactorGeVtoMillimeter 
      * ( p2 / GeV) * ( 1.0 / sqrt( sqr( (*cit)->momentum().m2() / GeV2 - 
					 (*cit)->momentum().mass2() / GeV2 ) + 
				    sqr( vmin2 / GeV2 ) ) );
    if ( fabs( displace2.mag() ) > dmax ) {
      displace2 *= dmax / fabs( displace2.mag() );
      // cout << "Cluster::calculateX ***extreme debugging*** : MAX DISPLACEMENT " << endl;
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
      if ( pcl.vect().dot( displace1.vect() ) * pcl.vect().dot( displace2.vect() )  <  0.0 ) {
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
      // cout << "Cluster::calculateX ***extreme debugging*** : s1minusS2 = " << s1minusS2 << endl;
      if ( s1minusS2 < 0 ) {
	s1 = 1.0;
	s2 = s1 - s1minusS2;
      } else if ( s1minusS2 > 0 ) {
	s2 = 1;
	s1 = s2 + s1minusS2;
      }
    }
    // Now, finally, determine the cluster position
    _position = 0.5 * ( pos1 + pos2 + s1*displace1 + s2*displace2 );   

    // Debugging
    // cout << "Cluster::calculateX ***extreme debugging***" << endl
    //      << "\t position1 = " << pos1
    //      << "\t invariant length = " << pos1.mag() << "  [mm] " << endl
    //      << "\t displace1 = " << displace1 
    //      << "\t invariant length = " << displace1.mag() << endl 
    //      << "\t position2 = " << pos2
    //      << "\t invariant length = " << pos2.mag() << "  [mm] " << endl
    //      << "\t displace2 = " << displace2 
    //      << "\t invariant length = " << displace2.mag() << endl
    //      << "\t 0.5 * ( pos1 + pos2 ) = " << 0.5 * ( pos1 + pos2 )
    //      << "\t invariant length = " << ( 0.5 * ( pos1 + pos2 ) ).mag() << endl
    //      << "\t s1 = " << s1 << "\t s2 = " << s2 << endl
    //      << "\t _position = " << _position
    //      << "\t invariant length = " << _position.mag() << endl;
    
  } // end else part of if ( _collecCompPtr.size() != 2 )   
}


void Cluster::addChildrenClusters(const tCluPtr child1, const tCluPtr child2) { 
  if (child1) _collecChildCluPtr.push_back(child1);  
  if (child2) _collecChildCluPtr.push_back(child2); 
}


void Cluster::addChildrenHadrons(const tPPtr had1, const tPPtr had2, const tPPtr had3) {
  if (had1) {
    _collecChildHadPtr.push_back(had1);
    if (had2) {
      _collecChildHadPtr.push_back(had2);
      if (had3) _collecChildHadPtr.push_back(had3);
    }
  }
}


bool Cluster::isBeamCluster() const {
  for (CollecCompPtr::const_iterator cit = _collecCompPtr.begin();
       cit != _collecCompPtr.end(); ++cit){
    if ( (*cit)->isBeamRemnant() ) return true;
  }
  return false;
}


void Cluster::isBeamCluster(tPPtr part) {
  for (CollecCompPtr::iterator it = _collecCompPtr.begin();
       it != _collecCompPtr.end(); ++it){
    if ( (*it)->pointerParticle() == part ) {
      (*it)->isBeamRemnant(true);
      break;
    }
  }
}
