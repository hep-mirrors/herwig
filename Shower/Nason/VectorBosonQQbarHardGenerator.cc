// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorBosonQQbarHardGenerator class.
//

#include "VectorBosonQQbarHardGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDF/BeamParticleData.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Shower/Base/ShowerTree.h"
#include "Herwig++/Shower/Base/Evolver.h"
#include "Herwig++/Shower/Base/KinematicsReconstructor.h"
#include "Herwig++/Shower/Base/PartnerFinder.h"
#include "Herwig++/Shower/Base/MECorrectionBase.h"
#include "Herwig++/Utilities/Histogram.h"
#include "ThePEG/Repository/EventGenerator.h"


using namespace Herwig;

void VectorBosonQQbarHardGenerator::persistentOutput(PersistentOStream & os) const {
  os << _alphaS << _prefactor << _power;
}

void VectorBosonQQbarHardGenerator::persistentInput(PersistentIStream & is, int) {
  is >> _alphaS >> _prefactor >> _power;
}

ClassDescription<VectorBosonQQbarHardGenerator> VectorBosonQQbarHardGenerator::initVectorBosonQQbarHardGenerator;
// Definition of the static class description member.

void VectorBosonQQbarHardGenerator::Init() {

  static ClassDocumentation<VectorBosonQQbarHardGenerator> documentation
    ("There is no documentation for the VectorBosonQQbarHardGenerator class");

  static Reference<VectorBosonQQbarHardGenerator,ShowerAlpha> interfaceShowerAlpha
    ("ShowerAlpha",
     "The object calculating the strong coupling constant",
     &VectorBosonQQbarHardGenerator::_alphaS, false, false, true, false, false);
}


NasonTreePtr VectorBosonQQbarHardGenerator::generateHardest(ShowerTreePtr tree) {

   // get the particles to be showered
  vector<tcPDPtr> partons;

  partons.clear();
  // find the incoming particles
  ShowerParticleVector incoming;

  vector<ShowerProgenitorPtr> particlesToShower;

  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cjt;

  for ( cjt=tree->outgoingLines().begin();
	cjt!=tree->outgoingLines().end(); ++cjt ) {
    incoming.push_back( cjt->first->progenitor() );
    partons.push_back( cjt->first->progenitor()->dataPtr() );
    particlesToShower.push_back( cjt->first );
    _quark.push_back( cjt->first->copy()->momentum() );
  }

  PPtr boson = tree->incomingLines().begin()->first->copy();

  // We are assuming quark first, swap order to ensure this
  if( partons[0]->id() < partons[1]->id() ) {
    swap( partons[0] ,partons[1] );
    swap( particlesToShower[0] , particlesToShower[1] );
    swap( _quark[0] , _quark[1] );
  }
  
  // Finds the boost to lab frame that should be applied to particles
  // generated in c.o.m frame by getEvent(). 
  // Final particles are boosted by this in getEvent()
  _eventFrame = getTransf();

  // vector to hold momenta of outgoingparticles
  vector<Lorentz5Momentum> pnew;

  // Generate emission and change _quark[0,1] and _g to momenta
  // of q, qbar and g after the hardest emission.
  getEvent(); 

  pnew.push_back( _quark[0] );
  pnew.push_back( _quark[1] );
  pnew.push_back( _g );
  pnew.push_back( boson->momentum() );

  ShowerParticleVector newparticles;
  tcPDPtr gluon = getParticleData( ParticleID::g );
  
  // make the particles for the NasonTree
  newparticles.push_back( new_ptr( ShowerParticle( partons[0], true ) ) );
  newparticles.push_back( new_ptr( ShowerParticle( partons[1], true ) ) );
  newparticles.push_back( new_ptr( ShowerParticle( gluon, true ) ) );
  newparticles.push_back( new_ptr( ShowerParticle( boson->dataPtr(), false ) ) );

  // set the momenta
  for( unsigned int ix=0; ix < 4; ++ix ) newparticles[ix]->set5Momentum( pnew[ix] );
  
  // create the intermediate off-shell emitting particle
  Lorentz5Momentum poff = pnew[_iemit] + pnew[2];
  poff.rescaleMass();
  newparticles.push_back( new_ptr( ShowerParticle( partons[_iemit], true ) ) );
  newparticles.back()->set5Momentum( poff );

  // find the sudakov for the branching
  BranchingList branchings = 
    evolver()->splittingGenerator()->finalStateBranchings();
  long index = abs( partons[ _iemit]->id() ); 

  IdList br(3);
  // types of particle in the branching
  br[0] = abs( newparticles[_iemit]->id() );
  br[1] = abs( newparticles[4]->id() );
  br[2] = newparticles[2]->id();

  SudakovPtr sudakov;
  // for loop cycles through the Branchinglist - to find the sudakov
  for( BranchingList::const_iterator cit = branchings.lower_bound(index); 
      cit != branchings.upper_bound(index); ++cit ) {
    IdList ids = cit->second.second;
    if( ids[0] == br[0] &&
	ids[1] == br[1] &&
	ids[2] == br[2] ) {
      sudakov = cit->second.first;
      break;
    }
  }
  // check sudakov has been created
  if( ! sudakov ) throw Exception() << "Can't find Sudakov for the hard emission in "
				 << "DrellYanHardGenerator::generateHardest()" 
				 << Exception::runerror;


  // create the vectors of NasonBranchings to create the NasonTree
  vector<NasonBranchingPtr> nasonin, nasonhard;
  // incoming boson
  nasonin.push_back( new_ptr( NasonBranching( newparticles[3], SudakovPtr(),
					    NasonBranchingPtr(), true ) ) );
  // outgoing particles from hard emission
  if( _iemit == 0 ) {
    nasonhard.push_back( new_ptr( NasonBranching( newparticles[4], sudakov,
					        NasonBranchingPtr(), false)));
    nasonhard.push_back( new_ptr( NasonBranching( newparticles[1], SudakovPtr(),
					      NasonBranchingPtr(), false) ) );
  }
  else {
    nasonhard.push_back( new_ptr( NasonBranching( newparticles[0], SudakovPtr(),
					        NasonBranchingPtr(), false ) ) );
    nasonhard.push_back( new_ptr( NasonBranching( newparticles[4], sudakov,
					        NasonBranchingPtr(), false ) ) );
  }
  // add g and q(bar) emitted particles as children of emitting particle
  nasonhard[_iemit]->addChild( new_ptr( NasonBranching( newparticles[_iemit], SudakovPtr(), 
							NasonBranchingPtr(), false) ) );
  nasonhard[_iemit]->addChild( new_ptr ( NasonBranching( newparticles[2], SudakovPtr(),
						     NasonBranchingPtr(), false) ) );
  // incoming boson add to nasonhard
  nasonhard.push_back( nasonin.back() );
  // make the tree
  NasonTreePtr nasontree = new_ptr( NasonTree( nasonhard, nasonin ) );

  // connect the ShowerParticles with the branchings
  // and set the maximum pt for the radiation
  set<NasonBranchingPtr> hard = nasontree->branchings();

  for( unsigned int ix=0 ; ix < particlesToShower.size() ; ++ix ) {
    //set the pt veto on both showers
    particlesToShower[ix]->maximumpT(_pt);

    for( set<NasonBranchingPtr>::const_iterator mit = hard.begin();
	mit != hard.end(); ++mit ) {
      //if the particle in current nasonbranching is to be showered and both ingoing/outgoing
      if( particlesToShower[ix]->progenitor()->id() == ((*mit)->_particle->id()) &&
	 particlesToShower[ix]->progenitor()->isFinalState() !=  ( *mit)->_incoming )
	{
	//connect the particle with that nason branching
	nasontree->connect(particlesToShower[ix]->progenitor(),*mit);

	if((*mit)->_incoming) {
	  (*mit)->_beam = particlesToShower[ix]->original()->parents()[0];
	}

	NasonBranchingPtr parent=(*mit)->_parent;
	while(parent) {
	  parent->_beam = particlesToShower[ix]->original()->parents()[0];
	  parent=parent->_parent;
	};

      }
    }
  }

   //do the colour reconnections
   //create the two colour lines
   ColinePtr blueLine = new_ptr( ColourLine() );
   ColinePtr greenLine = new_ptr( ColourLine() );
   
   //quark emits
   if ( _iemit == 0) {
     blueLine->addColoured( newparticles[0] );
     blueLine->addAntiColoured( newparticles[2] );
     greenLine->addColoured( newparticles[2] );
     greenLine->addAntiColoured( newparticles[1] );
     greenLine->addColoured( newparticles[4] );
   }
   else {
     blueLine->addColoured( newparticles[2] );
     blueLine->addAntiColoured( newparticles[1] );
     greenLine->addAntiColoured( newparticles[2] );
     greenLine->addAntiColoured( newparticles[4] );
     greenLine->addColoured( newparticles[0] );
   }

  // calculate the shower variables
  evolver()->showerModel()->kinematicsReconstructor()->
    reconstructDecayShower(nasontree,evolver());
  return nasontree;

  //calculate partonic event shapes from hard emission event.

  double thrust;
  if(_iemit == 0) thrust = _x1;
  else thrust = _x2;

  Energy mass;

  if( ( _quark[0] )(3) > ( _quark[1] )(3) && ( _quark[0] )(3) > _g(3) )
    mass = ( _quark[1] + _g ).mag();

  else if( ( _quark[1] )(3) > ( _quark[0] )(3) && ( _quark[1] )(3) > _g(3) )
    mass = ((_quark[0])+_g).mag();

  else mass = ( ( _quark[1] ) + ( _quark[0] ) ).mag();

  _ptplot.push_back( _pt / GeV );
  _yplot.push_back(_y);
  _x1plot.push_back(_x1);
  _x2plot.push_back(_x2);
  (*_hthrust) += 1. - thrust;
  (*_hthrustlow) += 1. - thrust;
  (*_hy) += _y;
  (*_hplow) += _pt / GeV;
  (*_hphigh) += _pt / GeV;
}

bool VectorBosonQQbarHardGenerator::canHandle(ShowerTreePtr tree) {
  if( tree->incomingLines().size() !=1 ){
    return false;    
  }
  if( ( tree->incomingLines().begin()->first->id() == 22 )
      && ( tree->incomingLines().begin()->first->progenitor()->id() == 23 ) )
    { 
    return false;
  }

  if( (tree->outgoingLines().size() != 2) )
    return false;

  unsigned int ix(0);
  ShowerParticlePtr part[2];
  ix = 0;
 
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cjt;

  for( cjt=tree->outgoingLines().begin();
       cjt!=tree->outgoingLines().end(); ++cjt ) {
    part[ix] = cjt->first->progenitor();
    ++ix;
  }   

  // outgoing particles check q qbar
  if( !( part[0]->id() > 0 && part[0]->id() < 6 && 
       part[1]->id() < 0 && part[1]->id() > -6 )
     && !( part[1]->id() > 0 && part[1]->id() < 6 && 
	  part[0]->id() < 0 && part[0]->id() > -6 ) ) {
    return false;
  }
  return true;
}

void VectorBosonQQbarHardGenerator::doinit() throw(InitException) {
  HardestEmissionGenerator::doinit();
}

void VectorBosonQQbarHardGenerator::dofinish() {
  HardestEmissionGenerator::dofinish();

  ofstream hist_out("hist3.top");

  _hy->topdrawOutput( hist_out, true, false, false, false,
		      "BLACK",
		      "e+e-  y", 
		      " ", 
		      " ",
		      " ", 
		      "y", 
		      " " );
 _hplow->topdrawOutput( hist_out, true, false, false, false,
		      "BLACK",
		      "e+e- low pT", 
		      " ", 
		      " ",
		      " ", 
		      "pT / GeV", 
		      " " );
 _hphigh->topdrawOutput( hist_out, true, false, false, false,
		      "BLACK",
		      "e+e- high pT", 
		      " ", 
		      " ",
		      " ", 
		      "pT / GeV", 
		      " " );
 _hthrust->topdrawOutput( hist_out, true, false, false, false,
		      "BLACK",
		      "e+e- 1-T", 
		      " ", 
		      " ",
		      " ", 
		      "1-T", 
		      " " );
 _hthrustlow->topdrawOutput( hist_out, true, false, false, false,
		      "BLACK",
		      "e+e- 1-T", 
		      " ", 
		      " ",
		      " ", 
		      "1-T", 
		      " " );

}

void VectorBosonQQbarHardGenerator::doinitrun() {
  _s = sqr( generator()->maximumCMEnergy() );

  _power = 1.0; 
  _prefactor = 1.0;

  _hthrust = new_ptr( Histogram( 0., 0.5, 100) );
  _hthrustlow = new_ptr( Histogram( 0., 0.01, 100) );
  _hy = new_ptr( Histogram( -8., 8., 100 ) );
  _hplow = new_ptr( Histogram( 0., 5., 100 ) );
  _hphigh = new_ptr( Histogram( 0., 100., 100) );
  HardestEmissionGenerator::doinitrun();
}

//private functions-internal workings

Lorentz5Momentum VectorBosonQQbarHardGenerator::getEvent(){
  
  Energy minp = 0.1*GeV;  
  Energy maxp = sqrt(0.5)*generator()->maximumCMEnergy();
  double miny = -8.;
  double maxy = 8.;
  double wgt;
  bool reject;
  Energy last_pt;
  last_pt = maxp;
  
  do {
    // veto algorithm for all C/Pt^pow except pow = 1
    //  _pt = GeV / pow (
    //	       pow ( GeV / last_pt ,_power - 1. ) - log( UseRandom::rnd() )
    //       * ( _power - 1. ) / _prefactor / ( maxy - miny ),
    //       1. / ( _power - 1. ) );

    // veto algorithm for pow = 1
      _pt = last_pt * pow( UseRandom::rnd() ,
    		    1. / ( maxy - miny ) / _prefactor );

     _y = UseRandom::rnd() * ( maxy - miny ) + miny;

     _x1 = 1. - _pt / sqrt( _s ) * exp( -_y );
     _x2 = 1. - _pt / sqrt( _s ) * exp( _y );

     wgt = getResult() / ( _prefactor * pow( GeV / _pt , _power ) );
     reject = UseRandom::rnd() > wgt || ! inRange();
     
     last_pt = _pt;
     if( inRange() ){
       //  last_pt=_pt;
       if ( wgt>1. ){ 
	 cerr << "PROBLEM!!!!"<< endl;
	 cerr<< " res = "<< getResult() << "overfn = " << 
	   _prefactor*pow(GeV/_pt,_power) << endl;
       }
     }
     //no emission event if p goes past p min - basically set to outside
     //of the histogram bounds (hopefully hist object just ignores it)
     if( _pt < minp ){
       _pt = 0. * GeV;//no emission event
       _y = -10;
       reject = false;
     }
  }while ( reject );

  //generate herwig variables (need to choose 1->2 splitting type)
  if( _x2 > _x1 ){
    _iemit = 0;
    _z = ( _x1 + _x2 - 1. ) / _x2;
    _ktild = ( 1. - _x2 ) / _z / ( 1. - _z );
  }
  else{
    _iemit = 1;
    _z = ( _x2 + _x1 - 1. ) / _x1;
    _ktild = ( 1. - _x1 ) / _z / ( 1. - _z ); 
  }

  _k = sqrt( _z * _z * ( 1. - _z ) * ( 1. - _z ) * _ktild );

  //construct vectors in com z frame
  constructVectors();
  azimuthal();
 
  //boost constructed vectors into the event frame
  _quark[0] = _eventFrame * _quark[0];
  _quark[1] = _eventFrame * _quark[1];
  _g = _eventFrame * _g;

  return _g;
  }

 double VectorBosonQQbarHardGenerator::getResult() {
   double res = 4. / 3. / pi * _pt / _s *
     ( sqr ( _x1 ) + sqr( _x2 ) ) / ( 1. - _x1 ) / ( 1. -_x2 ) * GeV;
   res *= _alphaS->value( sqr( _pt ) );
   return res;
 }

double VectorBosonQQbarHardGenerator::getMax(int num){
  
  double res;   
  Energy minp = 0.1 * GeV;  
  Energy maxp = sqrt(0.5) * generator()->maximumCMEnergy();
  double miny = -8.;
  double maxy = 8.;
  _max = 0.;

  for( int i = 0; i < num; i++ ) {
    do{
      _pt = UseRandom::rnd() * ( maxp - minp ) + minp;
      _y = UseRandom::rnd() * ( maxy - miny )+ miny;
      _x1 = 1.-_pt / sqrt(_s) * exp(-_y);
      _x2 = 1.-_pt / sqrt(_s) * exp(_y);
      res = getResult();
    }while( ! inRange() );
    if ( res * pow( ( _pt / GeV ), _power ) > _max ) 
      _max = res * pow( ( _pt / GeV ),_power );
  }

  _prefactor = _max;
  return _max;
}


//momentum construction
LorentzRotation VectorBosonQQbarHardGenerator::getTransf(){
 
  LorentzRotation transf( ( _quark[0] + _quark[1] ).findBoostToCM() );

  Lorentz5Momentum q1 = transf * ( _quark[0] );
  
  if( q1(0) == 0.0 ) transf.rotateZ( -pi / 2. );
  else transf.rotateZ( -atan( q1(1) / q1(0) ) );

  if( q1(2) == 0.0) transf.rotateY( pi / 2. );
  else transf.rotateY( atan ( sqrt( q1(0) * q1(0) + q1(1) * q1(1) ) / q1(2) ) );

  Lorentz5Momentum q2 = transf * ( _quark[0] );

  if( q2(2) < 0. )transf.rotateY( pi );

  transf.invert();

  return transf;
}

void VectorBosonQQbarHardGenerator::azimuthal() {
   if( UseRandom::rnd() < _x1 * _x1 / ( _x1 * _x1 + _x2 * _x2 ) ) {
     _r.setRotate( UseRandom::rnd() * 2. * pi , _quark[0].vect());
     _quark[1] = _r * _quark[1];
   }
   else{
     _r.setRotate( UseRandom::rnd() * 2. * pi, _quark[1].vect());
     _quark[0] = _r * _quark[0];
   }
    _g = _r * _g;
    return;
}

void VectorBosonQQbarHardGenerator::constructVectors(){

  _phi = UseRandom::rnd() * 2.* pi;

  //quark emitted
  if( _iemit == 0 ){

   _quark[0].setT( sqrt(_s) * ( _z + _k * _k / _z ) / 2. );
   _quark[0].setX( sqrt(_s) * _k * cos( _phi ) );
   _quark[0].setY( sqrt(_s) * _k * sin( _phi ) );
   _quark[0].setZ( sqrt(_s) * ( _z - _k * _k / _z ) / 2. );

   _quark[1].setT( sqrt(_s) * ( 1. - _k * _k / _z / ( 1.-_z ) ) / 2.);
   _quark[1].setX(0.);
   _quark[1].setY(0.);
   _quark[1].setZ( sqrt(_s)*( -1. + _k * _k / _z / (1.-_z) ) / 2. );
    
   _g.setT( sqrt(_s) * ( 1. - _z + _k * _k / ( 1.- _z ) ) / 2. );
   _g.setX( sqrt(_s) * -_k * cos ( _phi ) );
   _g.setY( sqrt(_s) * -_k * sin( _phi ) );
   _g.setZ( sqrt(_s) * ( 1. - _z - _k * _k / ( 1. - _z) ) / 2. );

  }
  //antiquark emitted
  else{

   _quark[0].setT( sqrt( _s ) * ( 1. - _k * _k / _z / ( 1. - _z ) ) / 2.);
   _quark[0].setX(0.);
   _quark[0].setY(0.);
   _quark[0].setZ( sqrt(_s) * ( 1. - _k * _k / _z / ( 1. - _z ) ) / 2.);

   _quark[1].setT( sqrt(_s) * ( _z + _k * _k / _z ) / 2. );
   _quark[1].setX( sqrt(_s) * _k * cos( _phi ) );
   _quark[1].setY( sqrt(_s) * _k * sin( _phi ) );
   _quark[1].setZ( sqrt(_s) * ( -_z + _k * _k / _z ) / 2. );
    
   _g.setT( sqrt(_s) * ( ( 1. - _z ) + _k * _k / ( 1. - _z ) ) / 2. );
   _g.setX( sqrt(_s) * -_k * cos( _phi ) );
   _g.setY( sqrt(_s) * -_k * sin( _phi ) );
   _g.setZ( sqrt(_s) * ( -( 1. - _z ) + _k * _k / ( 1. - _z ) ) / 2.);

  }
 
  return;  
}

