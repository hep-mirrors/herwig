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

  // Get the progenitors: Q and Qbar.
  vector<tcPDPtr> partons(2);
  ShowerProgenitorPtr QProgenitor,QbarProgenitor;
  if(tree->outgoingLines().begin()->first->progenitor()->id()>0) {
      QProgenitor=tree->outgoingLines().begin()->first;
      QbarProgenitor=tree->outgoingLines().rbegin()->first;
  } else {
      QProgenitor=tree->outgoingLines().rbegin()->first;
      QbarProgenitor=tree->outgoingLines().begin()->first;
  }
  partons[0]=QProgenitor->progenitor()->dataPtr();
  partons[1]=QbarProgenitor->progenitor()->dataPtr();
  _quark.resize(2);
  _quark[0]=QProgenitor->copy()->momentum();
  _quark[1]=QbarProgenitor->copy()->momentum();

  PPtr boson = tree->incomingLines().begin()->first->copy();

  // Finds the boost to lab frame that should be applied to particles
  // generated in c.o.m frame by getEvent():
  _eventFrame = getTransf();

  // Generate emission and change _quark[0,1] and _g to momenta
  // of q, qbar and g after the hardest emission:
  getEvent(); 

  // Make the particles for the NasonTree:
  tcPDPtr gluon_data = getParticleData(ParticleID::g);
  ShowerParticlePtr emitter(new_ptr(ShowerParticle(partons[_iemitter],true)));
  emitter->set5Momentum(_quark[_iemitter]); 
  ShowerParticlePtr spectator(new_ptr(ShowerParticle(partons[_ispectator],true)));
  spectator->set5Momentum(_quark[_ispectator]);  
  ShowerParticlePtr gluon(new_ptr(ShowerParticle(gluon_data,true)));
  gluon->set5Momentum(_g);  
  ShowerParticlePtr vboson(new_ptr(ShowerParticle(boson->dataPtr(),false)));
  vboson->set5Momentum(boson->momentum());  
  ShowerParticlePtr parent(new_ptr(ShowerParticle(partons[_iemitter],true)));
  Lorentz5Momentum parentMomentum(_quark[_iemitter]+_g);
  parentMomentum.rescaleMass();
  parent->set5Momentum(parentMomentum);

  // Find the sudakov for the branching:
  BranchingList branchings = 
    evolver()->splittingGenerator()->finalStateBranchings();
  long index = abs(emitter->id()); 
		
  SudakovPtr sudakov;
  // For loop cycles through the Branchinglist - to find the sudakov:
  for(BranchingList::const_iterator cit = branchings.lower_bound(index);
      cit != branchings.upper_bound(index); ++cit ) {
    IdList ids = cit->second.second;
    if(ids[0] == abs(emitter->id()) &&
       ids[1] == abs(parent->id())  &&
       ids[2] == gluon->id()) {
	sudakov = cit->second.first;
	break; 	    
    }
  }
  // Check sudakov has been created:
  if(!sudakov ) throw Exception() 
      << "No Sudakov for the hard emission in "
      << "VectorBosonQQbarHardGenerator::generateHardest()" 
      << Exception::runerror;
		
  // Create the vectors of NasonBranchings to create the NasonTree:
  vector<NasonBranchingPtr> oldBranchings, newBranchings;
  // Incoming boson:
  oldBranchings.push_back(new_ptr(NasonBranching(vboson,SudakovPtr(),
						 NasonBranchingPtr(),true)));
  // Outgoing particles from hard emission:
  NasonBranchingPtr spectatorBranch(new_ptr(NasonBranching(spectator,
				    SudakovPtr(),NasonBranchingPtr(),false)));
  NasonBranchingPtr emitterBranch(new_ptr(NasonBranching(parent,
				    sudakov,NasonBranchingPtr(),false)));
  emitterBranch->addChild(new_ptr(NasonBranching(emitter, 
				    SudakovPtr(),NasonBranchingPtr(),false)));
  emitterBranch->addChild(new_ptr(NasonBranching(gluon,
				    SudakovPtr(),NasonBranchingPtr(),false)));
  if(_iemitter== 0) {
    newBranchings.push_back(emitterBranch);
    newBranchings.push_back(spectatorBranch);
  } else {
    newBranchings.push_back(spectatorBranch);
    newBranchings.push_back(emitterBranch);
  }
  // Incoming boson add to newBranchings
  newBranchings.push_back(oldBranchings.back());
  // Make the tree
  NasonTreePtr nasontree = new_ptr(NasonTree(newBranchings,oldBranchings));
	
  // Connect the particles with the branchings
  // and set the maximum pt for the radiation
  set<NasonBranchingPtr> hard = nasontree->branchings();
  set<NasonBranchingPtr>::const_iterator mit;
  QProgenitor->maximumpT(_pt);
  QbarProgenitor->maximumpT(_pt);
  for(mit = hard.begin(); mit != hard.end(); ++mit) {
    // Connect the Q/Qbar progenitors to the Nason branchings.
    if(QProgenitor->progenitor()->id()==(*mit)->_particle->id())
        nasontree->connect(QProgenitor->progenitor(),*mit);
    if(QbarProgenitor->progenitor()->id()==(*mit)->_particle->id())
        nasontree->connect(QbarProgenitor->progenitor(),*mit);
  }

  // Create the two colour lines and connect the particles:
  ColinePtr blueLine  = new_ptr(ColourLine());
  ColinePtr greenLine = new_ptr(ColourLine());
  blueLine->addColoured(emitter,_iemitter);
  blueLine->addColoured(gluon,_ispectator);
  greenLine->addColoured(gluon,_iemitter);
  greenLine->addColoured(parent,_iemitter);
  greenLine->addColoured(spectator,_ispectator);
	
  // Calculate the shower variables:
  evolver()->showerModel()->kinematicsReconstructor()->
      reconstructDecayShower(nasontree,evolver());
  return nasontree;

  // Calculate partonic event shapes from hard emission event?
  double thrust;
  if(_iemitter == 0) thrust=_x1;
  else thrust=_x2;
  _ptplot.push_back(_pt/GeV);
  _yplot.push_back(_y);
  _x1plot.push_back(_x1);
  _x2plot.push_back(_x2);
  (*_hthrust) += 1.-thrust;
  (*_hthrustlow) += 1.-thrust;
  (*_hy) += _y;
  (*_hplow) += _pt/GeV;
  (*_hphigh) += _pt/GeV;
}

bool VectorBosonQQbarHardGenerator::canHandle(ShowerTreePtr tree) {

  if(tree->incomingLines().size()!=1) return false;    
  if((tree->incomingLines().begin()->first->id()==22)&&
  (tree->incomingLines().begin()->first->progenitor()->id()==23)) return false;

  map<ShowerProgenitorPtr,tShowerParticlePtr> outgoing=tree->outgoingLines();
  if(outgoing.size()!=2) return false;
  if(abs(outgoing.begin()->first->progenitor()->id())>6)  return false;
  if(outgoing.begin()->first->progenitor()->id()!=
     -1*outgoing.rbegin()->first->progenitor()->id())     return false;

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
    _iemitter   = 0;
    _ispectator = 1;
    _z = ( _x1 + _x2 - 1. ) / _x2;
    _ktild = ( 1. - _x2 ) / _z / ( 1. - _z );
  }
  else{
    _iemitter   = 1;
    _ispectator = 0;
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
   double res = 4. / 3. / Constants::pi * _pt / _s *
     ( sqr ( _x1 ) + sqr( _x2 ) ) / ( 1. - _x1 ) / ( 1. -_x2 ) * GeV;
   res *= _alphaS->value( sqr( _pt ) );
   return res;
 }

//momentum construction
LorentzRotation VectorBosonQQbarHardGenerator::getTransf(){
 
  LorentzRotation transf( ( _quark[0] + _quark[1] ).findBoostToCM() );

  Lorentz5Momentum q1 = transf * ( _quark[0] );
  using Constants::pi;
  if( q1.x() == 0.0*MeV ) transf.rotateZ( -pi / 2. );
  else transf.rotateZ( -atan( q1.y() / q1.x() ) );

  if( q1.z() == 0.0*MeV ) transf.rotateY( pi / 2. );
  else transf.rotateY( atan ( sqrt( q1.x() * q1.x() + q1.y() * q1.y() ) / q1.z() ) );

  Lorentz5Momentum q2 = transf * ( _quark[0] );

  if( q2.z() < 0.*MeV )transf.rotateY( pi );

  transf.invert();

  return transf;
}

void VectorBosonQQbarHardGenerator::azimuthal() {
  using Constants::pi;
   if( UseRandom::rnd() < _x1 * _x1 / ( _x1 * _x1 + _x2 * _x2 ) ) {
     _r.setRotate( UseRandom::rnd() * Constants::twopi , 
		   _quark[0].vect().unit());
     _quark[1] = _r * _quark[1];
   }
   else{
     _r.setRotate( UseRandom::rnd() * Constants::twopi, 
		   _quark[1].vect().unit());
     _quark[0] = _r * _quark[0];
   }
    _g = _r * _g;
    return;
}

void VectorBosonQQbarHardGenerator::constructVectors(){

  _phi = UseRandom::rnd() * Constants::twopi;
  if( _iemitter == 0 ){
   //quark emitted
   _quark[0].setT( sqrt(_s) * ( _z + _k * _k / _z ) / 2. );
   _quark[0].setX( sqrt(_s) * _k * cos( _phi ) );
   _quark[0].setY( sqrt(_s) * _k * sin( _phi ) );
   _quark[0].setZ( sqrt(_s) * ( _z - _k * _k / _z ) / 2. );

   _quark[1].setT( sqrt(_s) * ( 1. - _k * _k / _z / ( 1.-_z ) ) / 2.);
   _quark[1].setX(0.*MeV);
   _quark[1].setY(0.*MeV);
   _quark[1].setZ( sqrt(_s)*( -1. + _k * _k / _z / (1.-_z) ) / 2. );
  } else{
   //antiquark emitted
   _quark[0].setT( sqrt( _s ) * ( 1. - _k * _k / _z / ( 1. - _z ) ) / 2.);
   _quark[0].setX(0.*MeV);
   _quark[0].setY(0.*MeV);
   _quark[0].setZ( sqrt(_s) * ( 1. - _k * _k / _z / ( 1. - _z ) ) / 2.);

   _quark[1].setT( sqrt(_s) * ( _z + _k * _k / _z ) / 2. );
   _quark[1].setX( sqrt(_s) * _k * cos( _phi ) );
   _quark[1].setY( sqrt(_s) * _k * sin( _phi ) );
   _quark[1].setZ( sqrt(_s) * ( -_z + _k * _k / _z ) / 2. );
  }

  _g=-_quark[0]-_quark[1];
  _g.setT(sqrt(_s)+_g.t());

  return;  
}
