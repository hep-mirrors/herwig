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
#include "ThePEG/StandardModel/StandardModelBase.h"


using namespace Herwig;

void VectorBosonQQbarHardGenerator::persistentOutput(PersistentOStream & os) const {
  os << _alphaS << _alphaS_max << ounit( _Qg, GeV );
}

void VectorBosonQQbarHardGenerator::persistentInput(PersistentIStream & is, int) {
  is >> _alphaS >> _alphaS_max >> iunit( _Qg, GeV );
}

ClassDescription<VectorBosonQQbarHardGenerator> VectorBosonQQbarHardGenerator::initVectorBosonQQbarHardGenerator;
// Definition of the static class description member.

void VectorBosonQQbarHardGenerator::Init() {

  static ClassDocumentation<VectorBosonQQbarHardGenerator> documentation
    ("The VectorBosonQQbarHardGenerator class generates the hardest emission for"
     "vector boson decays to fermion-antifermion events in the Nason approach");

  static Reference<VectorBosonQQbarHardGenerator,ShowerAlpha> interfaceShowerAlpha
    ("ShowerAlpha",
     "The object calculating the strong coupling constant",
     &VectorBosonQQbarHardGenerator::_alphaS, false, false, true, false, false);
  
  static Parameter<VectorBosonQQbarHardGenerator, Energy> interfacePtMin
    ("minPt",
     "The pt cut on hardest emision generation",
     &VectorBosonQQbarHardGenerator::_Qg, GeV, 1.*GeV, 0*GeV, 100000.0*GeV,
     false, false, Interface::limited);
}

NasonTreePtr VectorBosonQQbarHardGenerator::generateHardest(ShowerTreePtr tree) {

  // Get the progenitors: Q and Qbar.
  ShowerProgenitorPtr 
    QProgenitor    = tree->outgoingLines().begin()->first,
    QbarProgenitor = tree->outgoingLines().rbegin()->first;
  if(QProgenitor->id()<0) swap(QProgenitor   ,QbarProgenitor);
  _partons.resize(2);
  _partons[0] = QProgenitor->progenitor()->dataPtr();
  _partons[1] = QbarProgenitor->progenitor()->dataPtr();
  // Get data for the gluon.
  tcPDPtr gluon_data = getParticleData(ParticleID::g);

  // momentum of the partons
  _quark.resize(2);
  _quark[0]=QProgenitor->copy()->momentum();
  _quark[1]=QbarProgenitor->copy()->momentum();
  // Set the existing mass entries in partons 5 vectors with the
  // once and for all.
  _quark[0].setMass(_partons[0]->mass());
  _quark[1].setMass(_partons[1]->mass());
  _g.setMass(0.*MeV);

  // PDG codes of the partons
  int q_id   (abs(QProgenitor->progenitor()->id()   ));
  int qbar_id(abs(QbarProgenitor->progenitor()->id()));

  // Get the gauge boson.
  _boson = tree->incomingLines().begin()->first->copy();

  // Get the gauge boson mass.
  _s = (_quark[0]+_quark[1])*(_quark[0]+_quark[1]);

  // Get the list of possible branchings.
  BranchingList branchings = 
    evolver()->splittingGenerator()->finalStateBranchings();

  // Find the sudakovs for the q/qbar->q/qbarg branchings.
  SudakovPtr q_sudakov,qbar_sudakov;
  // quark
  long q_index(q_id),qbar_index(qbar_id); 
  for(BranchingList::const_iterator cit = branchings.lower_bound(q_index);
      cit != branchings.upper_bound(q_index); ++cit ) {
    IdList ids = cit->second.second;
    if(ids[0]==q_id&&ids[1]==q_id&&ids[2]==gluon_data->id()) {
      q_sudakov = cit->second.first;
      break; 	    
    }
  }
  // antiquark
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
  // _Qg = q_sudakov->kinematicCutOff(q_sudakov->kinScale(),0.*GeV);

  // Generate emission and set _quark[0,1] and _g to be the 
  // momenta of q, qbar and g after the hardest emission:
  if(!getEvent()){
    generator()->log() << "generateHardest: NO emission was generated\n";
    QProgenitor->maximumpT(0.*GeV);
    QbarProgenitor->maximumpT(0.*GeV);
    return NasonTreePtr();
  }
  // Ensure the energies are greater than the constituent masses:
  for (int i=0; i<2; i++)
     if (_quark[i].e() < _partons[i]->constituentMass()) return NasonTreePtr();
  if (_g.e() < gluon_data->constituentMass()) return NasonTreePtr();

  // Set masses as done in VectorBosonQQbarMECorrection
  // (note this was already done at the start and the mass
  // component of the five-vector is never changed from that 
  // starting value, but just in case someone rewrites the
  // code above we include the following as in the MECorr code):
  _quark[0].setMass(_partons[0]->mass());
  _quark[1].setMass(_partons[1]->mass());
  _g.setMass(0.*MeV);

  // assign the emitter based on evolution scales
  // rather than for the correlations - we might want to try 
  // making this choice in the same way as VectorBosonQQbarMECorrection
  // (based on relative pT's). 
  _iemitter   = _quark[0]*_g>_quark[1]*_g ? 1 : 0;
  _ispectator = _iemitter==1              ? 0 : 1; 

  // Set the sudakov for the q/qbar->q/qbarg branching.
  SudakovPtr sudakov = _iemitter==0 ? q_sudakov : qbar_sudakov;

  // Make the particles for the NasonTree:
  ShowerParticlePtr emitter(new_ptr(ShowerParticle(_partons[_iemitter],true)));
  ShowerParticlePtr spectator(new_ptr(ShowerParticle(_partons[_ispectator],true)));
  ShowerParticlePtr gluon(new_ptr(ShowerParticle(gluon_data,true)));
  ShowerParticlePtr vboson(new_ptr(ShowerParticle(_boson->dataPtr(),false)));
  ShowerParticlePtr parent(new_ptr(ShowerParticle(_partons[_iemitter],true)));
  emitter->set5Momentum(_quark[_iemitter]); 
  spectator->set5Momentum(_quark[_ispectator]);  
  gluon->set5Momentum(_g);  
  vboson->set5Momentum(_boson->momentum());  
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
  } 
  else {
    allBranchings.push_back(spectatorBranch);
    allBranchings.push_back(emitterBranch);
  }
  // Add incoming boson to allBranchings
  allBranchings.push_back(spaceBranchings.back());

  // Make the NasonTree from the NasonBranching vectors.
  NasonTreePtr nasontree = new_ptr(NasonTree(allBranchings,spaceBranchings));
	
  // Set the maximum pt for all other emissions
  double xfact2 = _xq>_xqb ? sqr(_xq) : sqr(_xqb);
  Energy ptveto = _pt *sqrt((_xq+_xqb-1.)/xfact2);
  QProgenitor   ->maximumpT(ptveto);
  QbarProgenitor->maximumpT(ptveto);

  // Connect the particles with the branchings in the NasonTree
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

  // KMH - why don't we do the next step in reconstructDecayShower? 
  // Reset the momenta to ensure the correct momenta after shower recon
  // if emitter for Kleiss trick and shower are different
  for(map<ShowerParticlePtr,tNasonBranchingPtr>::const_iterator 
  	mit=nasontree->particles().begin();mit!=nasontree->particles().end();++mit)
    mit->first->set5Momentum(mit->second->_shower);

  // Return the NasonTree
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
  _alphaS_max = _alphaS->overestimateValue();

  _hthrust = new_ptr( Histogram( 0., 0.5, 100) );
  _hthrustlow = new_ptr( Histogram( 0., 0.1, 100) );
  _hy = new_ptr( Histogram( -6., 6., 100 ) );
  _hplow = new_ptr( Histogram( 0., 10., 100 ) );
  _hphigh = new_ptr( Histogram( 0., 100., 100) );
  HardestEmissionGenerator::doinitrun();
}

bool VectorBosonQQbarHardGenerator::getEvent(){

  
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
    
    _xq  = 1.-_pt/sqrt(_s)*exp(-_y);
    _xqb = 1.-_pt/sqrt(_s)*exp( _y);
    _xg  = 2.-_xq-_xqb;

    last_pt = _pt;

    wgt    = getResult()/(prefactor*GeV/_pt);
    reject = !inRange() || UseRandom::rnd()>wgt;
    
    if(inRange()&&wgt>1.) { 
      cerr << "VectorBosonQQbarHardGenerator::getEvent() excess weight.\n";
      cerr << "exact = "<< getResult() << "   crude = " << 
	prefactor*GeV/_pt << endl;
    }
    // No emission event if pt goes past pt_min.
    if(_pt<pt_min) {
      _pt = 0. * GeV; 
      return false;
    }
  } while (reject);
  
  // Mike's "Simple Prescription..." paper says that q or qb retains 
  // it's parton model direction with relative prob xq^2 or xqb^2 
  // respectively. _the_spectator_ retains it's parton model direction
  // but the transverse recoil is absorbed by the _the_emitter_. 
  // For xq->1 the gluon and particle qb are collinear, and acollinear 
  // to q, implying qb did the emitting, so for xq>>xqb xq is the 
  // spectator, xqb is the emitter, so select xq as spectator with prob 
  // xq^2/(xq^2+xqb^2)
  UseRandom::rnd()<(sqr(_xq)/(sqr(_xq)+sqr(_xqb))) ? 
    _iemitter = 1: _iemitter = 0 ; 
  _ispectator = !_iemitter;

  //construct vectors in com z frame
  constructVectors();
  return true;
}

double VectorBosonQQbarHardGenerator::getResult() {
  using Constants::pi;
  double CF = 4. / 3.; 
  // factor to get exact pt for argument of alphaS
  double xfact2 = _xq>_xqb ? sqr(_xq) : sqr(_xqb);

  //dimensionless quark mass squared M^2/s
  double rho = sqr( _partons[0]->mass() ) /_s;
  //set to zero for testing
  rho = 0.;

  //Vector and axial couplings- depends on quarks and Vector Boson type
  double sin2ThetaW = generator()->standardModel()->sin2ThetaW();
  
  //the vector and axial coupings
  double Vf, Af;
  double T3f, Qf;
  //down type quark
  if( abs( _partons[0]->id() ) == 1 || abs( _partons[0]->id() ) == 3 ||
      abs( _partons[0]->id() ) == 5 ){
    T3f = - 1. / 2.;
    Qf = - 1. / 3.;
  }
  //up type quark
  else {
    T3f =  1. / 2.;
    Qf =  2. / 3.;
  }

  //sets the couplings according to whether intermediate was a photon or Z
  if( _boson->dataPtr()->id() == 23 ){
    Vf = T3f - 2. * Qf * sin2ThetaW;
    Af = T3f;
  }
  else{
    Vf = 1.;
    Af = 0.;
  }

  //common factors in couplings and flux factors ignored
  //notation from hep-ph/0310083v2
  //born contribution in massive quark case
  double sigB =  sqr( Vf ) * ( 1. + 2. * rho ) + sqr( Af ) * ( 1. - 4. *  rho ) ;

  //Traces for the radiative contribution
  //Vector current trace
  double TrA =  ( sqr( _xq + 2. * rho ) + sqr( _xqb + 2. * rho ) 
		 + 2. * rho * ( sqr( 5. - _xq - _xqb ) - 19. + 4. * rho ) ) 
    / ( 1. - _xq ) / ( 1. - _xqb ) / ( 1. - 4. * rho )
    + 1 / sqr( 1. - _xq ) / sqr( 1. - _xqb ) *
    ( - 2. * rho  * sqr( 1. - _xq ) - 2. * rho  * sqr( 1. - _xqb ) );
  //Vector current trace
		 double TrV =  ( sqr( _xq + 2. * rho ) + sqr( _xqb + 2. * rho ) - 8. * rho * ( 1. + 2. * rho ) ) 
		   / ( 1. - _xq ) / ( 1. - _xqb ) / ( 1. + 2. * rho )
		   + 1 / sqr( 1. - _xq ) / sqr( 1. - _xqb ) *
		   ( - 2. * rho  * sqr( 1. - _xq ) - 2. * rho  * sqr( 1. - _xqb ) );

  //radiative contribution in quark massive case
  double sigR = _alphaS->value( sqr( _pt )*(_xq+_xqb-1.)/xfact2 ) * CF / 2. / pi 
    * ( sqr( Vf ) * TrV + sqr( Af ) * TrA ) 
    * 2. * _pt / _s * GeV;

  double res = sigR / sigB;

  return res;
}

void VectorBosonQQbarHardGenerator::constructVectors(){
  using Constants::twopi;
  using Constants::pi;
  // Find the boost from the lab to the c.o.m with the spectator 
  // along the -z axis, and then invert it.
  LorentzRotation eventFrame((_quark[0]+_quark[1]).findBoostToCM() );
  Lorentz5Momentum spectator = eventFrame*_quark[_ispectator];
  eventFrame.rotateZ(-spectator.phi());
  eventFrame.rotateY(-spectator.theta()-pi);
  eventFrame.invert();

  // Get the COM energy:
  Energy rts(sqrt(_s));

  // Extract the reduced (constituent) masses:
  double mu_e(_quark[_iemitter].mass()/rts);
  double mu_s(_quark[_ispectator].mass()/rts);
  double mu_g(_g.mass()/rts);

  // Get masses to avoid floating point errors later.
  Energy init_g_mass(_g.mass());
  Energy init_e_mass(_quark[_iemitter].mass());
  Energy init_s_mass(_quark[_ispectator].mass());

  // Get the Dalitz variables:
  double xe, xs, xg, b_xe, b_xs, b_xg;
  xe   = _iemitter==0   ? _xq : _xqb;
  xs   = _ispectator==0 ? _xq : _xqb;
  xg   = _xg;
  b_xe = _iemitter==0   ? _b_xq : _b_xqb;
  b_xs = _ispectator==0 ? _b_xq : _b_xqb;
  b_xg = _b_xg;

  // Get the cosines and sines of emitter w.r.t spectator:
  double c_se,s_se;
  c_se  = (xe*xs - 2.*(1.0 - mu_e*mu_e - mu_s*mu_s + mu_g*mu_g - xg))
          /(b_xs*b_xe);
  s_se  = _rt_mlambda/(2.*b_xe*b_xs);

  if(abs(c_se)>1.||abs(s_se)>1.) {
    throw Exception() 
      << "VectorBosonQQBarHardGenerator::constructVectors()"
      << "angle between emitter and spectator not physical." 
      << Exception::abortnow; 
  }

  // Construct momenta in boson COM frame with spectator along +/-Z axis: 
  _phi = UseRandom::rnd() * twopi;  

  // momentum of emitter
  _quark[_iemitter].setT( 0.5*rts * xe);
  _quark[_iemitter].setX( 0.5*rts * b_xe*s_se*cos(_phi));
  _quark[_iemitter].setY( 0.5*rts * b_xe*s_se*sin(_phi));
  _quark[_iemitter].setZ(-0.5*rts * b_xe*c_se);
  _quark[_iemitter].rescaleRho();
  // momentum of spectator
  _quark[_ispectator].setT( 0.5*rts * xs);
  _quark[_ispectator].setX( 0.*MeV);
  _quark[_ispectator].setY( 0.*MeV);
  _quark[_ispectator].setZ(-0.5*rts * b_xs);
  _quark[_ispectator].rescaleRho();
  // momentum of gluon
  _g=-_quark[0]-_quark[1];
  _g.setT(sqrt(_s)+_g.t());
  _g.setMass(init_g_mass);
  _g.rescaleRho();
  // boost constructed vectors into the event frame
  _quark[0] = eventFrame * _quark[0];
  _quark[1] = eventFrame * _quark[1];
  _g        = eventFrame * _g;

  // need to reset masses because for whatever reason the boost  
  // touches the mass component of the five-vector and can make  
  // zero mass objects acquire a floating point negative mass(!).
  _g.setMass(init_g_mass);
  _quark[_iemitter].setMass(init_e_mass);
  _quark[_ispectator].setMass(init_s_mass);
  _g.rescaleRho();
  _quark[_iemitter].rescaleRho();
  _quark[_ispectator].rescaleRho();
}

