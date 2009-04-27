// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PowhegEvolver class.
//

#include "PowhegEvolver.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Shower/Base/KinematicsReconstructor.h"
#include "Herwig++/Shower/Base/MECorrectionBase.h"
#include "Herwig++/Shower/Base/PartnerFinder.h"
#include "Herwig++/Shower/ShowerHandler.h"
#include "Herwig++/Shower/Base/HardTree.h"

using namespace Herwig;

void PowhegEvolver::doinitrun(){
  Evolver::doinitrun();
  //init hists and other variables
  _h_Xdiff = new_ptr( Histogram( -10., 2., 120 ) );
  _h_Ydiff = new_ptr( Histogram( -10., 2., 120 ) );
  _h_Zdiff = new_ptr( Histogram( -10., 2., 120 ) );
  _h_Ediff = new_ptr( Histogram( -10., 2., 120 ) );

  _no_events = 0;
  _mom_fails = 0;
}

void PowhegEvolver::dofinish() {
  Evolver::dofinish();

  ofstream hist_out1("momDiffHists.top");

  using namespace HistogramOptions;

  if( _hardonly ){
    _h_Xdiff->topdrawOutput( hist_out1, Frame| Ylog,
			     "BLACK",
			     "sum of abs(pdiff).x()", 
			     " ", 
			     " ",
			     " ", 
			     "log10( pdiff.x / GeV )", 
			     " " );
    _h_Ydiff->topdrawOutput( hist_out1, Frame| Ylog,
			     "BLACK",
			     "sum of abs(pdiff).y()", 
			     " ", 
			     " ",
			     " ", 
			     "log10( pdiff.y / GeV )", 
			     " " );
    _h_Zdiff->topdrawOutput( hist_out1, Frame| Ylog,
			     "BLACK",
			     "sum of abs(pdiff).z()", 
			     " ", 
			     " ",
			     " ", 
			     "log10( pdiff.z / GeV )", 
			     " " );
    _h_Ediff->topdrawOutput( hist_out1, Frame| Ylog,
			     "BLACK",
			     "sum of abs(pdiff).e()", 
			     " ", 
			     " ",
			     " ", 
			     "log10( pdiff.e / GeV )", 
			     " " );
    cerr<<"\n\n\npercentage of failures = "<< double(_mom_fails)/double(_no_events)*100.<<"\n\n\n";
  }
}

void PowhegEvolver::doinit() {
  Evolver::doinit();
  for(unsigned int ix=0;ix<_hardgenerator.size();++ix)
    _hardgenerator[ix]->setEvolver(this);
}

void PowhegEvolver::persistentOutput(PersistentOStream & os) const {
  os << _hardgenerator << _hardonly << _trunc_Mode;
}

void PowhegEvolver::persistentInput(PersistentIStream & is, int) {
  is >> _hardgenerator >> _hardonly >> _trunc_Mode;
}

ClassDescription<PowhegEvolver> PowhegEvolver::initPowhegEvolver;
// Definition of the static class description member.

void PowhegEvolver::Init() {

  static ClassDocumentation<PowhegEvolver> documentation
    ("The PowhegEvolver implements the POWHEG approach to MC\\@NLO");

  static RefVector<PowhegEvolver,HardestEmissionGenerator> interfaceHardGenerator
    ("HardGenerator",
     "The objects responsible for generating the hardestr emission",
     &PowhegEvolver::_hardgenerator, -1, false, false, true, false, false);

  static Switch<PowhegEvolver,bool> interfaceHardOnly
    ("HardOnly",
     "Only generate the emission supplied by the hardest emission"
     " generator, for testing only.",
     &PowhegEvolver::_hardonly, false, false, false);
  static SwitchOption interfaceHardOnlyNo
    (interfaceHardOnly,
     "No",
     "Generate full shower",
     false);
  static SwitchOption interfaceHardOnlyYes
    (interfaceHardOnly,
     "Yes",
     "Only the hardest emission",
     true);

  static Switch<PowhegEvolver,bool> interfaceTruncMode
    ("TruncatedShower", "Include the truncated shower?", 
     &PowhegEvolver::_trunc_Mode, 1, false, false);
  static SwitchOption interfaceTruncMode0
    (interfaceTruncMode,"No","Truncated Shower is OFF", 0);
  static SwitchOption interfaceTruncMode1
    (interfaceTruncMode,"Yes","Truncated Shower is ON", 1);

}
// here use hardemission stuff rather than me correction as in base class
vector<ShowerProgenitorPtr> PowhegEvolver::setupShower( bool hard ) {
  // generate the hardest emission
  hardestEmission();
  // set the colour partners
  setEvolutionPartners(hard,ShowerInteraction::QCD);
  // get the particles to be showered
  map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator cit;
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cjt;
  vector<ShowerProgenitorPtr> particlesToShower;
  // incoming particles
  for(cit=currentTree()->incomingLines().begin();
      cit!=currentTree()->incomingLines().end();++cit)
    particlesToShower.push_back(((*cit).first));
  assert((particlesToShower.size()==1&&!hard)||(particlesToShower.size()==2&&hard));
  // outgoing particles
  for(cjt=currentTree()->outgoingLines().begin();
      cjt!=currentTree()->outgoingLines().end();++cjt)
    particlesToShower.push_back(((*cjt).first));
  // return the answer
  return particlesToShower;
}

void PowhegEvolver::hardestEmission() {
  // see if there is an appropriate hard emission generator
  HardestEmissionGeneratorPtr currenthard;
  for( unsigned int ix = 0; ix < _hardgenerator.size(); ++ix ) {
    if( !_hardgenerator[ix]->canHandle( currentTree() ) ) continue;
    if( currenthard ) {
      ostringstream output;
      output << "There is more than one possible hard emission generator"
	     << " which could be used for ";
      map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator cit;
      map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cjt;
      for(cit=currentTree()->incomingLines().begin();
	  cit!=currentTree()->incomingLines().end();++cit)
	{output << cit->first->progenitor()->PDGName() << " ";}
      output << " -> ";
      for(cjt=currentTree()->outgoingLines().begin();
	  cjt!=currentTree()->outgoingLines().end();++cjt)
	{output << cjt->first->progenitor()->PDGName() << " ";}
      output << "in PowhegEvolver::hardestEmission()\n";
      throw Exception() << output << Exception::runerror;
    }
    currenthard=_hardgenerator[ix];
  }
  // if no suitable generator return
  _nasontree=HardTreePtr();
  if(!currenthard) return;
  // generate the hardest emission
  _nasontree = currenthard->generateHardest( currentTree() );
}

bool PowhegEvolver::truncatedTimeLikeShower(tShowerParticlePtr particle,
					    HardBranchingPtr branch,
					    ShowerInteraction::Type type) {
  Branching fb;
  unsigned int iout=0;
  tcPDPtr pdata[2];
  while (true) {
    // no truncated shower break
    if(!isTruncatedShowerON()||_hardonly) break;
    // generate emission
    fb=splittingGenerator()->chooseForwardBranching(*particle,1.,type);
    // no emission break
    if(!fb.kinematics) break;
    // check haven't evolved too far
    if(fb.kinematics->scale() < branch->scale()) {
      fb=Branching();
      break;
    }
    // get the particle data objects
    for(unsigned int ix=0;ix<2;++ix) pdata[ix]=getParticleData(fb.ids[ix+1]);
    if(particle->id()!=fb.ids[0]) {
      for(unsigned int ix=0;ix<2;++ix) {
	tPDPtr cc(pdata[ix]->CC());
	if(cc) pdata[ix]=cc;
      }
    }
    // find the truncated line
    iout=0;
    if(pdata[0]->id()!=pdata[1]->id()) {
      if(pdata[0]->id()==particle->id())       iout=1;
      else if (pdata[1]->id()==particle->id()) iout=2;
    }
    else if(pdata[0]->id()==particle->id()) {
      if(fb.kinematics->z()>0.5) iout=1;
      else                       iout=2;
    }
    // apply the vetos for the truncated shower
    // no flavour changing branchings
    if(iout==0) {
      particle->setEvolutionScale(fb.kinematics->scale());
      continue;
    }
    double zsplit = iout==1 ? fb.kinematics->z() : 1-fb.kinematics->z();
    // only if same interaction for forced branching 
    // and evolution
    if(type==branch->sudakov()->interactionType()) {
      if(zsplit < 0.5 || // hardest line veto
	 fb.kinematics->scale()*zsplit < branch->scale() ) { // angular ordering veto
	particle->setEvolutionScale(fb.kinematics->scale());
	continue;
      }
    }
    // pt veto
    if(fb.kinematics->pT() > progenitor()->maximumpT()) {
      particle->setEvolutionScale(fb.kinematics->scale());
      continue;
    }
    // should do base class vetos as well
    if(timeLikeVetoed(fb,particle)) {
      particle->setEvolutionScale(fb.kinematics->scale());
      continue;
    }
    break;
  }
  // if no branching set decay matrix and return
  if(!fb.kinematics) {
    // construct the kinematics for the hard emission
    ShoKinPtr showerKin=
      branch->sudakov()->createFinalStateBranching(branch->scale(),
						   branch->children()[0]->z(),
						   branch->phi(),
						   branch->children()[0]->pT());
    particle->setEvolutionScale(branch->scale() );
    showerKin->initialize( *particle,PPtr() );
    IdList idlist(3);
    idlist[0] = particle->id();
    idlist[1] = branch->children()[0]->branchingParticle()->id();
    idlist[2] = branch->children()[1]->branchingParticle()->id();
    fb = Branching( showerKin, idlist, branch->sudakov() );
    // Assign the shower kinematics to the emitting particle.
    particle->setShowerKinematics( fb.kinematics );
    // Assign the splitting function to the emitting particle. 
    // For the time being we are considering only 1->2 branching
    // Create the ShowerParticle objects for the two children of
    // the emitting particle; set the parent/child relationship
    // if same as definition create particles, otherwise create cc
    ShowerParticleVector theChildren;
    theChildren.push_back(new_ptr(ShowerParticle(branch->children()[0]->
						 branchingParticle()->dataPtr(),true)));
    theChildren.push_back(new_ptr(ShowerParticle(branch->children()[1]->
						 branchingParticle()->dataPtr(),true)));
    particle->showerKinematics()->
      updateChildren(particle, theChildren,type==branch->sudakov()->interactionType());
    // update the history if needed
    if(particle==currentTree()->getFinalStateShowerProduct(progenitor()))
      currentTree()->updateFinalStateShowerProduct(progenitor(),
						   particle,theChildren);
    currentTree()->addFinalStateBranching(particle,theChildren);
    // shower the first  particle
    if( branch->children()[0]->children().empty() ) {
      if( ! _hardonly )
	timeLikeShower(theChildren[0],type);
    }
    else {
      truncatedTimeLikeShower( theChildren[0],branch->children()[0],type);
    } 
    // shower the second particle
    if( branch->children()[1]->children().empty() ) {
      if( ! _hardonly )
	timeLikeShower( theChildren[1] , type);
    }
    else {
      truncatedTimeLikeShower( theChildren[1],branch->children()[1] ,type);
    }
    return true;
  }
  // has emitted
  // Assign the shower kinematics to the emitting particle.
  particle->setShowerKinematics(fb.kinematics);
  // Assign the splitting function to the emitting particle. 
  // For the time being we are considering only 1->2 branching
  // Create the ShowerParticle objects for the two children of
  // the emitting particle; set the parent/child relationship
  // if same as definition create particles, otherwise create cc
  ShowerParticleVector theChildren; 
  theChildren.push_back( new_ptr( ShowerParticle( pdata[0], true ) ) ); 
  theChildren.push_back( new_ptr( ShowerParticle( pdata[1], true ) ) );
  particle->showerKinematics()->
    updateChildren( particle, theChildren , true);
  // update the history if needed
  if( particle == currentTree()->getFinalStateShowerProduct( progenitor() ) )
    currentTree()->updateFinalStateShowerProduct( progenitor(),
						  particle, theChildren );
  currentTree()->addFinalStateBranching( particle, theChildren );
  // shower the first  particle
  if( iout == 1 ) truncatedTimeLikeShower( theChildren[0], branch , type );
  else            timeLikeShower( theChildren[0]  , type);
  // shower the second particle
  if( iout == 2 ) truncatedTimeLikeShower( theChildren[1], branch , type );
  else            timeLikeShower( theChildren[1]  , type);
  // branching has happened
  return true;
}

bool PowhegEvolver::truncatedSpaceLikeShower(tShowerParticlePtr particle, PPtr beam,
					     HardBranchingPtr branch,
					     ShowerInteraction::Type type) {
  Branching bb;
  // generate branching
  tcPDPtr part[2];
  while (true) {
    if( isTruncatedShowerON() || _hardonly ) break;
    bb = splittingGenerator()->chooseBackwardBranching( *particle, 
							beam, 1., beamParticle(), 
							type );
    if( !bb.kinematics || bb.kinematics->scale() < branch->scale() ) {
      bb = Branching();
      break;
    }
    // particles as in Sudakov form factor
    part[0] = getParticleData( bb.ids[0] );
    part[1] = getParticleData( bb.ids[2] );
    
    //is emitter anti-particle
    if( particle->id() != bb.ids[1]) {
      if( part[0]->CC() ) part[0] = part[0]->CC();
      if( part[1]->CC() ) part[1] = part[1]->CC();
    }
    double zsplit = bb.kinematics->z();
    // apply the vetos for the truncated shower
    // if doesn't carry most of momentum
    if(type==branch->sudakov()->interactionType() &&
       zsplit < 0.5) {
      particle->setEvolutionScale(bb.kinematics->scale() );
      continue;
    }
    // others
    if( part[0]->id() != particle->id() || // if particle changes type
	bb.kinematics->pT() > progenitor()->maximumpT() ||   // pt veto
	bb.kinematics->scale() < branch->scale()) { // angular ordering veto
      particle->setEvolutionScale(bb.kinematics->scale() );
      continue;
    }
    // and those from the base class
    if(spaceLikeVetoed(bb,particle)) {
      particle->setEvolutionScale(bb.kinematics->scale() );
      continue;
    }
    break;
  }
  if( !bb.kinematics ) {
    //do the hard emission
    double z(0.);
    HardBranchingPtr timelike;
    for( unsigned int ix = 0; ix < branch->children().size(); ++ix ) {
      if( branch->children()[ix]->status() ==HardBranching::Outgoing) {
	timelike = branch->children()[ix];
      }
      if( branch->children()[ix]->status() ==HardBranching::Incoming )
	z = branch->children()[ix]->z();
    }
    ShoKinPtr kinematics =
      branch->sudakov()->createInitialStateBranching( branch->scale(), z, branch->phi(),
						      branch->children()[0]->pT() );
    kinematics->initialize( *particle, beam );
    // assign the splitting function and shower kinematics
    particle->setShowerKinematics( kinematics );
    // For the time being we are considering only 1->2 branching
    // Now create the actual particles, make the otherChild a final state
    // particle, while the newParent is not
    ShowerParticlePtr newParent = 
      new_ptr( ShowerParticle( branch->branchingParticle()->dataPtr(), false ) );
    ShowerParticlePtr otherChild = 
      new_ptr( ShowerParticle( timelike->branchingParticle()->dataPtr(),
			       true, true ) );
    ShowerParticleVector theChildren;
    theChildren.push_back( particle ); 
    theChildren.push_back( otherChild );
    particle->showerKinematics()->
      updateParent( newParent, theChildren, type==branch->sudakov()->interactionType() );
    // update the history if needed
    currentTree()->updateInitialStateShowerProduct( progenitor(), newParent );
    currentTree()->addInitialStateBranching( particle, newParent, otherChild );
    // for the reconstruction of kinematics, parent/child
    // relationships are according to the branching process:
    // now continue the shower
    bool emitted=false;
    if(!_hardonly) {
      if( branch->parent() ) {
	emitted = truncatedSpaceLikeShower( newParent, beam, branch->parent() , type);
      }
      else {
	emitted = spaceLikeShower( newParent, beam , type);
      }
    }
    if( !emitted ) {
      if( intrinsicpT().find( progenitor() ) == intrinsicpT().end() ) {
	kinematics->updateLast( newParent, ZERO, ZERO );
      }
      else {
	pair<Energy,double> kt = intrinsicpT()[progenitor()];
	kinematics->updateLast( newParent,
				kt.first*cos( kt.second ),
				kt.first*sin( kt.second ) );
      }
    }
    particle->showerKinematics()->
      updateChildren( newParent, theChildren,
		      type==branch->sudakov()->interactionType() );
    if(_hardonly) return true;
    // perform the shower of the final-state particle
    if( timelike->children().empty() ) {
      timeLikeShower( otherChild , type);
    }
    else {
      truncatedTimeLikeShower( otherChild, timelike , type);
    }
    // return the emitted
    return true;
  }
  // assign the splitting function and shower kinematics
  particle->setShowerKinematics( bb.kinematics );
  // For the time being we are considering only 1->2 branching
  // Now create the actual particles, make the otherChild a final state
  // particle, while the newParent is not
  ShowerParticlePtr newParent = new_ptr( ShowerParticle( part[0], false ) );
  ShowerParticlePtr otherChild = new_ptr( ShowerParticle( part[1], true, true ) );
  ShowerParticleVector theChildren; 
  theChildren.push_back( particle ); 
  theChildren.push_back( otherChild );
  particle->showerKinematics()->updateParent( newParent, theChildren , true);
  // update the history if needed
  currentTree()->updateInitialStateShowerProduct( progenitor(), newParent );
  currentTree()->addInitialStateBranching( particle, newParent, otherChild );
  // for the reconstruction of kinematics, parent/child
  // relationships are according to the branching process:
  // now continue the shower
  bool emitted = truncatedSpaceLikeShower( newParent, beam, branch,type);
  // now reconstruct the momentum
  if( !emitted ) {
    if( intrinsicpT().find( progenitor() ) == intrinsicpT().end() ) {
      bb.kinematics->updateLast( newParent, ZERO, ZERO );
    }
    else {
      pair<Energy,double> kt = intrinsicpT()[ progenitor() ];
      bb.kinematics->updateLast( newParent,
				 kt.first*cos( kt.second ),
				 kt.first*sin( kt.second ) );
    }
  }
  particle->showerKinematics()->updateChildren( newParent, theChildren , true);
  // perform the shower of the final-state particle
  timeLikeShower( otherChild , type);
  // return the emitted
  return true;
}

void PowhegEvolver::setEvolutionPartners(bool hard,
					 ShowerInteraction::Type type) {
  // if no hard tree use the methid in the base class
  if(!_nasontree) {
    Evolver::setEvolutionPartners(hard,type);
    return;
  }
  // match the particles in the ShowerTree and NasonTree
  if(!_nasontree->connect(currentTree()))
    throw Exception() << "Can't match trees in "
		      << "PowhegEvolver::setEvolutionPartners()"
		      << Exception::eventerror;
  // sort out the colour partners
  vector<ShowerParticlePtr> particles;
  map<ShowerProgenitorPtr, ShowerParticlePtr>::const_iterator cit;
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cjt;
  for(cit=currentTree()->incomingLines().begin();
      cit!=currentTree()->incomingLines().end();++cit)
    particles.push_back(cit->first->progenitor());
  // outgoing particles
  for(cjt=currentTree()->outgoingLines().begin();
      cjt!=currentTree()->outgoingLines().end();++cjt)
    particles.push_back(cjt->first->progenitor());
  // find the partner
  for(unsigned int ix=0;ix<particles.size();++ix) {
    tHardBranchingPtr partner = 
      _nasontree->particles()[particles[ix]]->colourPartner();
    if(!partner) continue;
    for(map<ShowerParticlePtr,tHardBranchingPtr>::const_iterator
	  it=_nasontree->particles().begin();
	it!=_nasontree->particles().end();++it) {
      if(it->second==partner) particles[ix]->setPartner(it->first);
    }
    if(!particles[ix]->partner()) 
      throw Exception() << "Can't match partners in "
			<< "PowhegEvolver::setEvolutionPartners()"
			<< Exception::eventerror;
  }
  // Set the initial evolution scales
  showerModel()->partnerFinder()->
    setInitialEvolutionScales(particles,!hard,type,false);
}

bool PowhegEvolver::startTimeLikeShower(ShowerInteraction::Type type) {
  if(_nasontree) {
    map<ShowerParticlePtr,tHardBranchingPtr>::const_iterator 
      eit=_nasontree->particles().end(),
      mit = _nasontree->particles().find(progenitor()->progenitor());
    if( mit != eit && !mit->second->children().empty() ) {
      return truncatedTimeLikeShower(progenitor()->progenitor(), mit->second ,type);
    }
  }
  return  _hardonly ? false :
    timeLikeShower(progenitor()->progenitor() ,type) ;
}

bool PowhegEvolver::startSpaceLikeDecayShower(Energy maxscale,Energy minimumMass,
					      ShowerInteraction::Type type) {
  if(_nasontree) {
    map<ShowerParticlePtr,tHardBranchingPtr>::const_iterator 
      eit =_nasontree->particles().end(),
      mit = _nasontree->particles().find(progenitor()->progenitor());
    if( mit != eit && mit->second->parent() ) {
      HardBranchingPtr branch=mit->second;
      while(branch->parent()) branch=branch->parent();
      return truncatedSpaceLikeDecayShower(progenitor()->progenitor(),maxscale,
					   minimumMass, branch ,type);
    }
  }
  return  _hardonly ? false :
    spaceLikeDecayShower(progenitor()->progenitor(),maxscale,minimumMass,type);
}

bool PowhegEvolver::startSpaceLikeShower(PPtr parent,
					 ShowerInteraction::Type type) {
  if(_nasontree) {
    map<ShowerParticlePtr,tHardBranchingPtr>::const_iterator 
      eit =_nasontree->particles().end(),
      mit = _nasontree->particles().find(progenitor()->progenitor());
    if( mit != eit && mit->second->parent() ) {
      return truncatedSpaceLikeShower( progenitor()->progenitor(),
				       parent, mit->second->parent(), type );
    } 
  }
  return  _hardonly ? false :
    spaceLikeShower(progenitor()->progenitor(),parent,type);
}

bool PowhegEvolver::checkShowerMomentum( vector<ShowerProgenitorPtr> particlesToShower ){
  if( _hardonly && _nasontree){
    _no_events ++;
    // extract the particles from end point of the shower
    multimap<long,PPtr> outgoing;
    //loop over all final state particles in particlesToShower
    for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
      if(!particlesToShower[ix]->progenitor()->isFinalState()) continue;
      //insert 3 outgoing particles into outgoing map
      PPtr temp = particlesToShower[ix]->progenitor();
      if(!temp->children().empty()) {
	outgoing.insert(make_pair(temp->children()[0]->id(),temp->children()[0]));
	outgoing.insert(make_pair(temp->children()[1]->id(),temp->children()[1]));
      }
      else {
	outgoing.insert(make_pair(temp->id(),temp));
      }
    }
    // extract the particles from the nason tree
    vector<PPtr> outb;
    for(set<HardBranchingPtr>::const_iterator it = 
	  _nasontree->branchings().begin();
	it != _nasontree->branchings().end(); ++it)  {
      if(!(**it).children().empty()) {
	for(vector<HardBranchingPtr>::const_iterator jt=(**it).children().begin();
	    jt!=(**it).children().end();++jt) {
	  outb.push_back((**jt).branchingParticle());
	}
      }
      else if((**it).branchingParticle()->isFinalState()) {
	outb.push_back((**it).branchingParticle());
      }
    }
    static Energy eps = 0.00001*GeV;

    Lorentz5Momentum diff_tot( 0.*GeV, 0.*GeV,  0.*GeV,  0.*GeV );

    for(unsigned int ix=0;ix<outb.size();++ix) {
      //pointer to the shower particle that matches the powheg particle
      //do by looking at all with matching ids and choosing the one with
      //minimum momentum difference magnitude
      PPtr ShowerMatched;
      Energy magDiff = 10000. * GeV;
      
      //loop over all outgoing particles that match the id of outb
      //choose one that is best match for momentum
      multimap< long, PPtr >::const_iterator cjt;
      for( cjt =  outgoing.lower_bound( outb[ix]->id() );
	   cjt != outgoing.upper_bound( outb[ix]->id() );
	   ++cjt ) {
	Lorentz5Momentum mom_diff = outb[ix]->momentum() 
	  - cjt->second->momentum();
	Energy magComp = abs( mom_diff.mag() );
	if( magComp < magDiff ){
	  magDiff = magComp;
	  ShowerMatched = cjt->second;
	}
      }
      if( !ShowerMatched ) {
	generator()->log()<< "could not find a matching shower particle for: \n"
			  << *outb[ix] << " !!!!\n";
      }

      Lorentz5Momentum pdiff = outb[ix]->momentum();
      pdiff -= ShowerMatched->momentum();

      Lorentz5Momentum abs_pdiff( abs( pdiff.x() ),  abs( pdiff.y() ),
				  abs( pdiff.z() ),  abs( pdiff.t() ) );
      diff_tot += abs_pdiff;
      if( diff_tot.x() > eps || diff_tot.y() > eps || diff_tot.z() > eps || diff_tot.t() > eps ){
	generator()->log() << "show part: " << *ShowerMatched << "\n"
			   << "powh part: " << *outb[ix] << "\n"
			   << "diff    =  " << pdiff / GeV << "\n\n";
      }
      
    }
    if( diff_tot.x() > eps || diff_tot.y() > eps || diff_tot.z() > eps || diff_tot.t() > eps ){
      generator()->log() << "\n total abs momentum difference = "
			 << diff_tot / GeV <<"\n \n \n";
      _mom_fails ++;
    }
    (*_h_Xdiff) += log10( diff_tot.x() / GeV );
    (*_h_Ydiff) += log10( diff_tot.y() / GeV );
    (*_h_Zdiff) += log10( diff_tot.z() / GeV );
    (*_h_Ediff) += log10( diff_tot.t() / GeV );
  }
  return true;
}

bool PowhegEvolver::
truncatedSpaceLikeDecayShower(tShowerParticlePtr particle, Energy maxscale,
			      Energy minmass, HardBranchingPtr branch,
			      ShowerInteraction::Type type) {
  Branching fb;
  unsigned int iout=0;
  tcPDPtr pdata[2];
  while (true) {
    // no truncated shower break
    if(!isTruncatedShowerON()||_hardonly) break;
    fb=splittingGenerator()->chooseDecayBranching(*particle,maxscale,minmass,1.,type);
    // return if no radiation
    if(!fb.kinematics) break;
    // check haven't evolved too far
    if(fb.kinematics->scale() < branch->scale()) {
      fb=Branching();
      break;
    }
    // get the particle data objects
    for(unsigned int ix=0;ix<2;++ix) pdata[ix]=getParticleData(fb.ids[ix+1]);
    if(particle->id()!=fb.ids[0]) {
      for(unsigned int ix=0;ix<2;++ix) {
	tPDPtr cc(pdata[ix]->CC());
	if(cc) pdata[ix]=cc;
      }
    }
    // find the truncated line
    iout=0;
    if(pdata[0]->id()!=pdata[1]->id()) {
      if(pdata[0]->id()==particle->id())       iout=1;
      else if (pdata[1]->id()==particle->id()) iout=2;
    }
    else if(pdata[0]->id()==particle->id()) {
      if(fb.kinematics->z()>0.5) iout=1;
      else                       iout=2;
    }
    // apply the vetos for the truncated shower
    // no flavour changing branchings
    if(iout==0) {
      particle->setEvolutionScale(fb.kinematics->scale());
      continue;
    }
    double zsplit = iout==1 ? fb.kinematics->z() : 1-fb.kinematics->z();
    if(type==branch->sudakov()->interactionType()) {
      if(zsplit < 0.5 || // hardest line veto
	 fb.kinematics->scale()*zsplit < branch->scale() ) { // angular ordering veto
	particle->setEvolutionScale(fb.kinematics->scale());
	continue;
      }
    }
    // pt veto
    if(fb.kinematics->pT() > progenitor()->maximumpT()) {
      particle->setEvolutionScale(fb.kinematics->scale());
      continue;
    }
    // should do base class vetos as well
    // if not vetoed break
    if(!spaceLikeDecayVetoed(fb,particle)) break;
    // otherwise reset scale and continue
    particle->setEvolutionScale(fb.kinematics->scale());
  }
  // if no branching set decay matrix and return
  if(!fb.kinematics) {
    // construct the kinematics for the hard emission
    cerr << "forcing " << branch->children()[0]->z() << "\n";
    ShoKinPtr showerKin=
      branch->sudakov()->createDecayBranching(branch->scale(),
					      branch->children()[0]->z(),
					      branch->phi(),
					      branch->children()[0]->pT());
    particle->setEvolutionScale(branch->scale() );
    showerKin->initialize( *particle,PPtr() );
    IdList idlist(3);
    idlist[0] = particle->id();
    idlist[1] = branch->children()[0]->branchingParticle()->id();
    idlist[2] = branch->children()[1]->branchingParticle()->id();
    fb = Branching( showerKin, idlist, branch->sudakov() );
    // Assign the shower kinematics to the emitting particle.
    particle->setShowerKinematics( fb.kinematics );
    // Assign the splitting function to the emitting particle. 
    // For the time being we are considering only 1->2 branching
    // Create the ShowerParticle objects for the two children of
    // the emitting particle; set the parent/child relationship
    // if same as definition create particles, otherwise create cc
    ShowerParticleVector theChildren;
    theChildren.push_back(new_ptr(ShowerParticle(branch->children()[0]->
						 branchingParticle()->dataPtr(),true)));
    theChildren.push_back(new_ptr(ShowerParticle(branch->children()[1]->
						 branchingParticle()->dataPtr(),true)));
    particle->showerKinematics()->
      updateChildren(particle, theChildren,
		     type==branch->sudakov()->interactionType());
    if(theChildren[0]->id()==particle->id()) {
      // update the history if needed
      currentTree()->updateInitialStateShowerProduct(progenitor(),theChildren[0]);
      currentTree()->addInitialStateBranching(particle,theChildren[0],theChildren[1]);
      // shower the space-like particle
      if( branch->children()[0]->children().empty() ) {
	if( ! _hardonly ) spaceLikeDecayShower(theChildren[0],maxscale,minmass,type);
      }
      else {
	truncatedSpaceLikeDecayShower( theChildren[0],maxscale,minmass,
				       branch->children()[0],type);
      }
      // shower the second particle
      if( branch->children()[1]->children().empty() ) {
	if( ! _hardonly ) timeLikeShower( theChildren[1] , type);
      }
      else {
	truncatedTimeLikeShower( theChildren[1],branch->children()[1] ,type);
      }
    }
    else {
      // update the history if needed
      currentTree()->updateInitialStateShowerProduct(progenitor(),theChildren[1]);
      currentTree()->addInitialStateBranching(particle,theChildren[0],theChildren[1]);
      // shower the space-like particle
      if( branch->children()[1]->children().empty() ) {
	if( ! _hardonly ) spaceLikeDecayShower(theChildren[1],maxscale,minmass,type);
      }
      else {
	truncatedSpaceLikeDecayShower( theChildren[1],maxscale,minmass,
				       branch->children()[1],type);
      }
      // shower the second particle
      if( branch->children()[0]->children().empty() ) {
	if( ! _hardonly ) timeLikeShower( theChildren[0] , type);
      }
      else {
	truncatedTimeLikeShower( theChildren[0],branch->children()[0] ,type);
      }
    }
    return true;
  }
  // has emitted
  // Assign the shower kinematics to the emitting particle.
  particle->setShowerKinematics(fb.kinematics);
  // For the time being we are considering only 1->2 branching
  // Create the ShowerParticle objects for the two children of
  // the emitting particle; set the parent/child relationship
  // if same as definition create particles, otherwise create cc
  ShowerParticleVector theChildren; 
  theChildren.push_back(new_ptr(ShowerParticle(pdata[0],true))); 
  theChildren.push_back(new_ptr(ShowerParticle(pdata[1],true)));
  particle->showerKinematics()->updateChildren(particle, theChildren,true);
  // In the case of splittings which involves coloured particles,
  // set properly the colour flow of the branching.
  // update the history if needed
  currentTree()->updateInitialStateShowerProduct(progenitor(),theChildren[0]);
  currentTree()->addInitialStateBranching(particle,theChildren[0],theChildren[1]);
  // shower the first  particle
  truncatedSpaceLikeDecayShower(theChildren[0],maxscale,minmass,branch,type);
  // shower the second particle
  timeLikeShower(theChildren[1],type);
  // branching has happened
  return true;
}
