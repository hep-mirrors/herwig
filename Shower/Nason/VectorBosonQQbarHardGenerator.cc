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
  os << _alphaS << _alphaS_max;
}

void VectorBosonQQbarHardGenerator::persistentInput(PersistentIStream & is, int) {
  is >> _alphaS >> _alphaS_max;
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
  ShowerProgenitorPtr 
    QProgenitor   =tree->outgoingLines().begin()->first,
    QbarProgenitor=tree->outgoingLines().rbegin()->first;
  if(QProgenitor->id()<0) swap(QProgenitor   ,QbarProgenitor);
  partons[0]=QProgenitor->progenitor()->dataPtr();
  partons[1]=QbarProgenitor->progenitor()->dataPtr();
  _quark.resize(2);
  _quark[0]=QProgenitor->copy()->momentum();
  _quark[1]=QbarProgenitor->copy()->momentum();
  int q_id(abs(QProgenitor->progenitor()->id()));
  int qbar_id(abs(QbarProgenitor->progenitor()->id()));

  // Get the gauge boson.
  PPtr boson = tree->incomingLines().begin()->first->copy();
  // Get data for the gluon.
  tcPDPtr gluon_data = getParticleData(ParticleID::g);

  // Get the list of possible branchings.
  BranchingList branchings = 
    evolver()->splittingGenerator()->finalStateBranchings();

  // Find the sudakovs for the q/qbar->q/qbarg branchings.
  SudakovPtr q_sudakov,qbar_sudakov;
  long q_index(q_id),qbar_index(qbar_id); 
  for(BranchingList::const_iterator cit = branchings.lower_bound(q_index);
      cit != branchings.upper_bound(q_index); ++cit ) {
    IdList ids = cit->second.second;
    if(ids[0]==q_id&&ids[1]==q_id&&ids[2]==gluon_data->id()) {
      q_sudakov = cit->second.first;
      break; 	    
    }
  }
  for(BranchingList::const_iterator cit = branchings.lower_bound(qbar_index);
      cit != branchings.upper_bound(qbar_index); ++cit ) {
    IdList ids = cit->second.second;
    if(ids[0]==qbar_id&&ids[1]==qbar_id&&ids[2]==gluon_data->id()) {
      qbar_sudakov = cit->second.first;
      break; 	    
    }
  }

  // Check sudakovs got created:
  if(!q_sudakov||!qbar_sudakov) throw Exception() 
    << "No Sudakov for the hard emission in "
    << "VectorBosonQQbarHardGenerator::generateHardest()" 
    << Exception::runerror;
  
  // Get all the gluon mass assuming massless quarks: 
  _Qg = q_sudakov->kinematicCutOff(q_sudakov->kinScale(),0.*GeV);

  // Finds the boost to lab frame that should be applied to particles
  // generated in c.o.m frame by getEvent():
  _eventFrame = getTransf();

  // Generate emission and change _quark[0,1] and _g to momenta
  // of q, qbar and g after the hardest emission:
  getEvent(); 

  // now we need to assign the emitter based on evolution scales
  // rather than for the correlations
  _iemitter   = _quark[0]*_g>_quark[1]*_g ? 1 : 0;
  _ispectator = _iemitter==1              ? 0 : 1; 
  // Set the sudakov for the q/qbar->q/qbarg branching.
  SudakovPtr sudakov = _iemitter==0 ? q_sudakov : qbar_sudakov;

  // Make the particles for the NasonTree:
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

  // Create the vectors of NasonBranchings to create the NasonTree:
  vector<NasonBranchingPtr> spaceBranchings,allBranchings;
  // Incoming boson:
  spaceBranchings.push_back(new_ptr(NasonBranching(vboson,SudakovPtr(),
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
  if(_iemitter==0) {
    allBranchings.push_back(emitterBranch);
    allBranchings.push_back(spectatorBranch);
  } else {
    allBranchings.push_back(spectatorBranch);
    allBranchings.push_back(emitterBranch);
  }
  // Incoming boson add to allBranchings
  allBranchings.push_back(spaceBranchings.back());
  // Make the tree
  NasonTreePtr nasontree = new_ptr(NasonTree(allBranchings,spaceBranchings));
	
  // Connect the particles with the branchings in the NasonTree
  // and set the maximum pt for the radiation
  double xfact2 = _xb>_xc ? sqr(_xb) : sqr(_xc);
  Energy ptveto = _pt *sqrt((_xb+_xc-1.)/xfact2);
  //if no emission event set ptveto to minpt = _Qg
  if( _pt < _Qg ) ptveto = _Qg;
  QProgenitor->maximumpT(ptveto);
  QbarProgenitor->maximumpT(ptveto);
  nasontree->connect(QProgenitor->progenitor(),allBranchings[0]);
  nasontree->connect(QbarProgenitor->progenitor(),allBranchings[1]);

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
  // Calculate partonic event shapes from hard emission event?
  double thrust = _iemitter == 0 ? _xb : _xc;
  _ptplot.push_back(_pt/GeV);
  _yplot.push_back(_y);
  _xbplot.push_back(_xb);
  _xcplot.push_back(_xc);
  (*_hthrust) += 1.-thrust;
  (*_hthrustlow) += 1.-thrust;
  (*_hy) += _y;
  (*_hplow) += _pt/GeV;
  (*_hphigh) += _pt/GeV;
  return nasontree;
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

  ofstream hist_out("hardEmissHists.top");
  using namespace HistogramOptions;

  _hy->topdrawOutput( hist_out, Frame,
		      "BLACK",
		      "e+e-  y", 
		      " ", 
		      " ",
		      " ", 
		      "y", 
		      " " );
 _hplow->topdrawOutput( hist_out, Frame,
		      "BLACK",
		      "e+e- low pT", 
		      " ", 
		      " ",
		      " ", 
		      "pT / GeV", 
		      " " );
 _hphigh->topdrawOutput( hist_out, Frame,
		      "BLACK",
		      "e+e- high pT", 
		      " ", 
		      " ",
		      " ", 
		      "pT / GeV", 
		      " " );
 _hthrust->topdrawOutput( hist_out, Frame,
		      "BLACK",
		      "e+e- 1-T", 
		      " ", 
		      " ",
		      " ", 
		      "1-T", 
		      " " );
 _hthrustlow->topdrawOutput( hist_out, Frame,
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

  _alphaS_max = _alphaS->overestimateValue();

  _hthrust = new_ptr( Histogram( 0., 0.5, 100) );
  _hthrustlow = new_ptr( Histogram( 0., 0.1, 100) );
  _hy = new_ptr( Histogram( -6., 6., 100 ) );
  _hplow = new_ptr( Histogram( 0., 10., 100 ) );
  _hphigh = new_ptr( Histogram( 0., 100., 100) );
  HardestEmissionGenerator::doinitrun();
}

//private functions-internal workings

Lorentz5Momentum VectorBosonQQbarHardGenerator::getEvent(){
  
  Energy pt_min = _Qg;  
  Energy pt_max = 0.5*sqrt(_s);
// Define over valued y_max & y_min according to the associated pt_min cut.
  double y_max  =  acosh(pt_max/pt_min);
  double y_min  = -acosh(pt_max/pt_min);
  double wgt;
  bool reject;
  Energy last_pt(pt_max);
  double prefactor = 2.*_alphaS_max*(4./3.)/Constants::pi;
  
  do {
    reject=true;
    _pt = last_pt*pow(UseRandom::rnd(),1./(prefactor*(y_max-y_min)));
    _y  = UseRandom::rnd()*(y_max-y_min)+y_min;
    
    _xb = 1.-_pt/sqrt(_s)*exp(-_y);
    _xc = 1.-_pt/sqrt(_s)*exp( _y);

    last_pt = _pt;

    wgt    = getResult()/(prefactor*GeV/_pt);
    reject = !inRange() || UseRandom::rnd()>wgt;
    
    if(inRange()&&wgt>1.) { 
      cerr << "VectorBosonQQbarHardGenerator::getEvent() excess weight.\n";
      cerr << "exact = "<< getResult() << "   crude = " << 
	prefactor*GeV/_pt << endl;
    }
    // No emission event if pt goes past pt_min - basically set to outside
    // of the histogram bounds (hopefully hist object just ignores it).
    if(_pt<pt_min) {
      _pt = 0. * GeV; 
      reject = false;
    }
  } while (reject);
  
// _z and _ktild are the z & \tilde{kappa} variables in the new variables
// paper, for the final-final colour connnection with massless partons (b=c=0).
  if(UseRandom::rnd()<sqr(_xb)/(sqr(_xb)+sqr(_xc))) {
    _iemitter   = 1;
    _ispectator = 0;
    _z = (_xc+_xb-1.)/_xb; 
    _ktild = (1.-_xb)/_z/(1.-_z); 
  } else{
    _iemitter   = 0;
    _ispectator = 1;
    _z = (_xb+_xc-1.)/_xc;
    _ktild = (1.-_xc)/_z/(1.-_z);
  }

  _k = _z*(1.-_z)*sqrt(_ktild);
     
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
    ( sqr ( _xb ) + sqr( _xc ) ) / ( 1. - _xb ) / ( 1. -_xc ) * GeV;
  //res *= _alphaS->value( sqr( _pt ) );
  // KH - shouldn't we use the ACTUAL pt inside alphaS instead i.e. 
  double xfact2 = _xb>_xc ? sqr(_xb) : sqr(_xc);
  res *= _alphaS->value( sqr( _pt )*(_xb+_xc-1.)/xfact2 );
  return res;
}

//momentum construction
LorentzRotation VectorBosonQQbarHardGenerator::getTransf() {
 
  LorentzRotation transf((_quark[0]+_quark[1]).findBoostToCM() );
  Lorentz5Momentum quark = transf*_quark[0];
  transf.rotateZ(-quark.phi());
  transf.rotateY(-quark.theta());
  transf.invert();

  return transf;
}

void VectorBosonQQbarHardGenerator::azimuthal() {
  using Constants::pi;
// KH - Check no emission events - they gave nans before when in azimuthal().
  _r.setRotate( UseRandom::rnd() * Constants::twopi, 
	_quark[_ispectator].vect().unit());
  _quark[_iemitter] = _r*_quark[_iemitter];
  _g = _r*_g;
  return;
}

void VectorBosonQQbarHardGenerator::constructVectors(){
  // Construct momenta in boson COM frame with spectator along +/-Z axis: 
  _phi = UseRandom::rnd() * Constants::twopi;

//KH - remind me to check these.
  _quark[_iemitter].setT(sqrt(_s)*(_z+_k*_k/_z)/2.);
  _quark[_iemitter].setX(sqrt(_s)*_k*cos(_phi));
  _quark[_iemitter].setY(sqrt(_s)*_k*sin(_phi));
  _quark[_iemitter].setZ(sqrt(_s)*(_z-_k*_k/_z)/2.);

  _quark[_ispectator].setT(sqrt(_s)*(1.-_k*_k/_z/(1.-_z ))/2.);
  _quark[_ispectator].setX(0.*MeV);
  _quark[_ispectator].setY(0.*MeV);
  _quark[_ispectator].setZ(sqrt(_s)*(-1.+_k*_k/_z/(1.-_z))/2.);
  
  // If Qbar emits then reflect the z components of emitter & spectator.
  if(_iemitter==1) { 
    _quark[_iemitter].setZ(-1.*sqrt(_s)*(_z-_k*_k/_z)/2.);
    _quark[_ispectator].setZ(-1.*sqrt(_s)*(-1.+_k*_k/_z/(1.-_z))/2.);
  }

  _g=-_quark[0]-_quark[1];
  _g.setT(sqrt(_s)+_g.t());

  return;  
}
