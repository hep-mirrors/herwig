// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HardTree class.
//

#include "HardTree.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Repository/CurrentGenerator.h"

using namespace Herwig;

HardTree::HardTree(vector<HardBranchingPtr> branchings,
		   vector<HardBranchingPtr> spacelike) 
  : _branchings(branchings.begin(),branchings.end()),
    _spacelike (spacelike .begin(),spacelike .end())
{}

bool HardTree::findNodes() {
  //clear all containers to be filled
  _theExternals.clear();
  _theNodes.clear();
  _theInternals.clear();
  //call function to recursively fill _theExternals and _theNodes
  for( set< HardBranchingPtr>::const_iterator cit = _branchings.begin();
       cit != _branchings.end(); ++cit ){
    if( (*cit)->incoming() ) continue;
    fillNodes( *cit, HardBranchingPtr() );
  }
  //make _internals from nodes
  //get the intermediates - there is one intermediate for each node
  for( map< HardBranchingPtr, Energy >::const_iterator cit = _theNodes.begin();
       cit != _theNodes.end(); ++cit ) {
    //end scale and intermediate are given by theNodes
    Energy endScale = cit->second;
    Energy startScale;
    if( cit->first->parent() ) {
      HardBranchingPtr intParent = cit->first->parent();
      startScale = intParent->scale() * cit->first->z();
    }
    else startScale = cit->first->branchingParticle()->evolutionScale();
    long intID = cit->first->branchingParticle()->id();
    if ( startScale < endScale ) endScale = startScale;
    _theInternals.insert( make_pair( intID, make_pair(  startScale, endScale  ) ) );
  }
  return true;
}

bool HardTree::fillNodes( HardBranchingPtr branch, HardBranchingPtr parentBranch ){
  if( ! branch->children().empty() ) {
    _theNodes.insert( make_pair( branch, branch->scale() ) );
    //set the parent of branching (used in finding the internal lines)
    branch->parent( parentBranch );
    fillNodes( branch->children()[0], branch );
    fillNodes( branch->children()[1], branch );
  }
  //external branching found
  else  {
    _theExternals.insert( make_pair( branch->branchingParticle(), parentBranch ) );
    if( parentBranch && ( !_lowestPt || branch->pT() < _lowestPt->pT() ) )
      _lowestPt = parentBranch;
  }
  return true;
}

void HardBranching::setMomenta(LorentzRotation R,double aparent,
			       Lorentz5Momentum ptparent,
			       bool setMomentum) {
  if(setMomentum) _original=R*_particle->momentum();
  // compute the shower variables
  Energy2 dot = _n*_p;
  double alpha = (_original*_n)/dot;
  _z=alpha/aparent;
  double beta = ((_original*_p)-alpha*sqr(_p.mass()))/dot;
  _qt = _original - alpha*_p - beta*_n - _z*ptparent;
  _pt=sqrt(max(-_qt*_qt,ZERO));
  // reconstruct children
  for(unsigned int ix=0;ix<_children.size();++ix) {
    _children[ix]->_p=_p;
    _children[ix]->_n=_n;
    _children[ix]->setMomenta( R, alpha, _qt + _z*ptparent, setMomentum);
  }
  // calculate the evolution scale and phi
  if(!_children.empty()) {
    double z = _children[0]->_z;
    Energy pt = _children[0]->_pt;
    IdList ids(3);
    ids[0]=_particle->id();
    ids[1]=_children[0]->_particle->id();
    ids[2]=_children[1]->_particle->id();
    _scale=_sudakov->calculateScale(z,pt,ids,_incoming ? 1 : 0);
    // get the pt vector
    Lorentz5Momentum vect=_children[0]->_qt;
    Boost beta_bb = -(_p+ _n).boostVector();
    Lorentz5Momentum p_bb = _p;
    vect.boost(beta_bb);
    p_bb.boost( beta_bb );
    Axis axis(p_bb.vect().unit());
    LorentzRotation rot;
    if(axis.perp2()>0.) {
      double sinth(sqrt(1.-sqr(axis.z())));
      rot.setRotate(-acos(axis.z()),
		    Axis(-axis.y()/sinth,axis.x()/sinth,0.));
      vect.transform(rot);
    }
    else if(axis.z()<0.) {
      vect.setZ(vect.z());
    }
    _phi= atan2(vect.y(),vect.x());
    if(_phi<0.)                 _phi+=Constants::twopi;
    if(_children[1]->_incoming) _phi+=Constants::pi;
  }
}

bool HardTree::connect(ShowerTreePtr shower) {
  _particles.clear();
  // extract the progenitors from the ShowerTree
  vector<ShowerProgenitorPtr> progenitors = shower->extractProgenitors();
  // connect the trees up
  for(set<HardBranchingPtr>::iterator it=branchings().begin();
      it!=branchings().end();++it) {
    CurrentGenerator::log() << "looknig for match of "
			    << *(**it).branchingParticle() << "\n";
    Energy2 dmin(1e30*GeV2);
    tShowerParticlePtr partner;   
    for(unsigned int ix=0;ix<progenitors.size();++ix) {
      if((**it).branchingParticle()->id()!=progenitors[ix]->progenitor()->id()) continue;
      if((**it).incoming()==progenitors[ix]->progenitor()->isFinalState()) continue;
      Energy2 dtest = 
	sqr(progenitors[ix]->progenitor()->momentum().x()-(**it).showerMomentum().x())+
	sqr(progenitors[ix]->progenitor()->momentum().y()-(**it).showerMomentum().y())+
	sqr(progenitors[ix]->progenitor()->momentum().z()-(**it).showerMomentum().z())+
	sqr(progenitors[ix]->progenitor()->momentum().t()-(**it).showerMomentum().t());
      if(dtest<dmin) {
	partner = progenitors[ix]->progenitor();
	dmin = dtest;
      }
    }
    if(!partner) return false;
    connect(partner,*it);
    if((**it).incoming()) {
      double z((**it).z());
      tHardBranchingPtr parent=(**it).parent();
      while (parent) {
	z *= parent->z();
	parent=parent->parent();
      }
      partner->x(z);
    }
  }
  // return false if not matched
  return particles().size()==progenitors.size();
}

HardBranching::HardBranching(ShowerParticlePtr particle, SudakovPtr sudakov,
			     tHardBranchingPtr parent,bool incoming) 
  : _particle(particle), _incoming(incoming), _parent(parent),
    _sudakov(sudakov)
{}

void HardBranching::fixColours() {
  if(!_sudakov) return;
  if(!_incoming&&_children.empty()) return;
  if(_incoming && !_parent) return;
  if(_incoming)
    _sudakov->splittingFn()->
      colourConnection(_parent->_particle,_particle,
		       _parent->children()[1]->_particle,true);
  else
    _sudakov->splittingFn()->
      colourConnection(_particle,_children[0]->_particle,
		       _children[1]->_particle,false);
}

ostream & Herwig::operator<<(ostream & os, const HardTree & x) {
  os << "Output of HardTree " << &x << "\n";
  for(set<HardBranchingPtr>::const_iterator it=x._branchings.begin();
      it!=x._branchings.end();++it) {
    os << "Hard Particle: " << *(**it).branchingParticle() << " has Sudakov " 
       << (**it).sudakov() << "\n";
    os << "It's colour lines are " << (**it).branchingParticle()->colourLine() << "\t" 
       <<  (**it).branchingParticle()->antiColourLine() << "\n";
    for(unsigned int iy=0;iy<(**it).children().size();++iy) {
      os << "\t Children : " << *(**it).children()[iy]->branchingParticle()
	 << "\n";
      os << "It's colour lines are " 
	 << (**it).children()[iy]->branchingParticle()->colourLine() << "\t" 
	 <<  (**it).children()[iy]->branchingParticle()->antiColourLine() << "\n";
    }
  }
  for(set<HardBranchingPtr>::const_iterator it=x._spacelike.begin();
      it!=x._spacelike.end();++it) {
    os << "SpaceLike: " << *(**it).branchingParticle() << " has Sudakov" 
       << (**it).sudakov() << "\n";
    os << "It's colour lines are " 
       << (**it).branchingParticle()->colourLine() << "\t" 
       << (**it).branchingParticle()->antiColourLine() << "\n";
    for(unsigned int iy=0;iy<(**it).children().size();++iy) {
      os << "\t Children: " << *(**it).children()[iy]->branchingParticle()
	 << "\n";
      os << "It's colour lines are " 
	 << (**it).children()[iy]->branchingParticle()->colourLine() << "\t" 
	 << (**it).children()[iy]->branchingParticle()->antiColourLine() << "\n";
    }
  }
  return os;
}

void HardTree::fillHardScales( HardBranchingPtr branch, vector< pair< Energy, double > > & currentLine ){
  if( branch->children().empty() ) return;
  else{
    //child[0] continues currentline child[1] creates a new line
    //copy contents of old line into newline
    vector< pair< Energy, double > > newHardLine = currentLine;
    currentLine.push_back( make_pair( branch->scale(),
				      branch->children()[0]->z() ) );
    fillHardScales( branch->children()[0], currentLine );
    newHardLine.push_back( make_pair( branch->scale(),
				      branch->children()[1]->z() ) );
    fillHardScales( branch->children()[1], newHardLine );
    _hard_line_scales.push_back( newHardLine );
  }
}

bool HardTree::checkHardOrdering() {
  //this function also caculates sum of pts of all branchings
  _total_pt = 0. * GeV;
  _hard_line_scales.clear();
  for( set<HardBranchingPtr>::const_iterator it = 
	 this->branchings().begin();
       it != this->branchings().end(); ++it)  {
    if( (*it)->incoming() ) continue;
    //if the branching has children then fill qtilde and z of branching and
    //continue recursively on the children 
   
    if( ! (*it)->children().empty() ) {
      vector< pair< Energy, double > > new_hard_line1;
      new_hard_line1.push_back( make_pair( (*it)->scale(), (*it)->children()[0]->z() ) );
      //pts of children are equal so just add once 
      _total_pt += (*it)->children()[0]->pT();

      fillHardScales( (*it)->children()[0], new_hard_line1 );
      vector< pair< Energy, double > > new_hard_line2;
      new_hard_line2.push_back( make_pair( (*it)->scale(), (*it)->children()[1]->z() ) );
      fillHardScales( (*it)->children()[1], new_hard_line2 );
      _hard_line_scales.push_back( new_hard_line1 );
      _hard_line_scales.push_back( new_hard_line2 );
    }    
  }
  //  bool ordered = true;
  for(unsigned int ix = 0; ix < _hard_line_scales.size(); ix++ ){
    for(unsigned int jx = 0; jx < _hard_line_scales[ix].size(); jx++ ){
      if( jx == 0 ) 
	continue;
      //angular ordering condition: z_1*q_1 > q2
      if( _hard_line_scales[ix][jx].first  > _hard_line_scales[ix][jx - 1].first * _hard_line_scales[ix][ jx - 1 ].second )
	return false;
    }
  }
  return true;
}
Energy HardTree::lowestPt( int jetMeasureMode, Energy2 s ){
  //check to see we found the _lowest pt correctly
  if(  _lowestPt->children().size() != 2 ) {
    cerr<< "HardTree: wrong no. children = " << _lowestPt->children().size() << " \n";
    return 0.*GeV;
  }
  if( ! _lowestPt->children()[0] ) {
    cerr<< "HardTree: null child \n";
    return 0.*GeV;
  }

  Energy kt_softest = _lowestPt->children()[0]->pT();
  if( jetMeasureMode != 1 ){
    Energy pt = kt_softest;
    double z = _lowestPt->children()[0]->z();
   
    Energy2 m0 = sqr( _lowestPt->branchingParticle()->nominalMass() );
    Energy2 m1 = sqr( _lowestPt->children()[0]->branchingParticle()->nominalMass() );
    Energy2 m2 = sqr( _lowestPt->children()[1]->branchingParticle()->nominalMass() );

    double lambda = sqrt( 1. - 4.*m0/s );
    double beta1 = 2.*( m1 - sqr(z)*m0 + sqr(pt) )
      / z / lambda / ( lambda + 1. ) / s;
    double beta2 = 2.*( m2 - sqr( 1. - z )*m0 + sqr(pt) )
      / ( 1. - z ) / lambda / ( lambda + 1. ) / s;

    Energy E1 = sqrt(s)/2.*( z + lambda*beta1 );
    Energy E2 = sqrt(s)/2.*( (1.-z) + lambda*beta2 );
    Energy Z1 = sqrt(s)/2.*lambda*( z - beta1 );
    Energy Z2 = sqrt(s)/2.*lambda*( (1.-z) - beta2 );

    double costheta = ( Z1*Z2 - sqr(pt) )
      / sqrt( sqr(Z1)+sqr(pt) ) / sqrt( sqr(Z2)+sqr(pt) );

    Energy2 kt_measure;
    if(  jetMeasureMode == 0 )
      kt_measure = 2.*min( sqr(E1), sqr(E2) )*( 1. - costheta );
    else if( jetMeasureMode == 2 )
      kt_measure = 2.*sqr(E1)*sqr(E2)/sqr(E1+E2)*( 1. - costheta );
    
    kt_softest = sqrt( kt_measure );
  }
  return kt_softest;
}
