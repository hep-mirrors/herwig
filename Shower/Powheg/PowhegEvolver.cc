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
#include "Herwig++/Utilities/Histogram.h"


using namespace Herwig;

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
  setColourPartners(hard);
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

bool PowhegEvolver::truncatedTimeLikeShower( tShowerParticlePtr particle,
					    HardBranchingPtr branch ) {
  Branching fb;
  unsigned int iout=0;
  tcPDPtr pdata[2];
  while (true) {
    // no truncated shower break
    if(!isTruncatedShowerON()||_hardonly) break;
    // generate emission
    fb=splittingGenerator()->chooseForwardBranching(*particle,1.);
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
    if(zsplit < 0.5 || // hardest line veto
       fb.kinematics->scale()*zsplit < branch->scale() || // angular ordering veto
       fb.kinematics->pT() > progenitor()->maximumpT()) { // pt veto
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
    particle->showerKinematics()->updateChildren(particle, theChildren);
    // update the history if needed
    if(particle==currentTree()->getFinalStateShowerProduct(progenitor()))
      currentTree()->updateFinalStateShowerProduct(progenitor(),
						   particle,theChildren);
    currentTree()->addFinalStateBranching(particle,theChildren);

    

    // shower the first  particle
    if( branch->children()[0]->children().empty() ) {
      if( ! _hardonly )
	timeLikeShower(theChildren[0]);
    }
    else {
      truncatedTimeLikeShower( theChildren[0],branch->children()[0] );
    } 
    // shower the second particle
    if( branch->children()[1]->children().empty() ) {
      if( ! _hardonly )
	timeLikeShower( theChildren[1] );
    }
    else {
      truncatedTimeLikeShower( theChildren[1],branch->children()[1] );
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
  particle->showerKinematics()->updateChildren( particle, theChildren );
  // update the history if needed
  if( particle == currentTree()->getFinalStateShowerProduct( progenitor() ) )
    currentTree()->updateFinalStateShowerProduct( progenitor(),
						  particle, theChildren );
  currentTree()->addFinalStateBranching( particle, theChildren );
  // shower the first  particle
  if( iout == 1 ) truncatedTimeLikeShower( theChildren[0], branch );
  else            timeLikeShower( theChildren[0] );
  // shower the second particle
  if( iout == 2 ) truncatedTimeLikeShower( theChildren[1], branch );
  else            timeLikeShower( theChildren[1] );
  // branching has happened
  return true;
}

bool PowhegEvolver::truncatedSpaceLikeShower(tShowerParticlePtr particle, PPtr beam,
					    HardBranchingPtr branch) {
  //  bool vetoed(true);
  Branching bb;
  // generate branching
  tcPDPtr part[2];
  while (true) {
    if( isTruncatedShowerON() || _hardonly ) break;
    bb = splittingGenerator()->chooseBackwardBranching( *particle, 
							beam, 1., beamParticle() );
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
    if( part[0]->id() != particle->id() || // if particle changes type
	zsplit < 0.5 || // if doesn't carry most of momentum
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
      if( !branch->children()[ix]->incoming() ) {
	timelike = branch->children()[ix];
      }
      if( branch->children()[ix]->incoming() ) z = branch->children()[ix]->z();
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
    particle->showerKinematics()->updateParent( newParent, theChildren );
    // update the history if needed
    currentTree()->updateInitialStateShowerProduct( progenitor(), newParent );
    currentTree()->addInitialStateBranching( particle, newParent, otherChild );
    // for the reconstruction of kinematics, parent/child
    // relationships are according to the branching process:
    // now continue the shower
    bool emitted=false;
    if(!_hardonly) {
      if( branch->parent() ) {
	emitted = truncatedSpaceLikeShower( newParent, beam, branch->parent() );
      }
      else {
	emitted = spaceLikeShower( newParent, beam );
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
    particle->showerKinematics()->updateChildren( newParent, theChildren );
    if(_hardonly) return true;
    // perform the shower of the final-state particle
    if( timelike->children().empty() ) {
      timeLikeShower( otherChild );
    }
    else {
      truncatedTimeLikeShower( otherChild, timelike );
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
  particle->showerKinematics()->updateParent( newParent, theChildren );
  // update the history if needed
  currentTree()->updateInitialStateShowerProduct( progenitor(), newParent );
  currentTree()->addInitialStateBranching( particle, newParent, otherChild );
  // for the reconstruction of kinematics, parent/child
  // relationships are according to the branching process:
  // now continue the shower
  bool emitted = truncatedSpaceLikeShower( newParent, beam, branch);
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
  particle->showerKinematics()->updateChildren( newParent, theChildren );
  // perform the shower of the final-state particle
  timeLikeShower( otherChild );
  // return the emitted
  return true;
}

void PowhegEvolver::setColourPartners(bool hard) {
  // if no hard tree use the methid in the base class
  if(!_nasontree) {
    Evolver::setColourPartners(hard);
    return;
  }
  // match the particles in the ShowerTree and NasonTree
  if(!_nasontree->connect(currentTree()))
    throw Exception() << "Can't match trees in "
		      << "PowhegEvolver::setColourPartners()"
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
    if(!particles[ix]->dataPtr()->coloured()) continue;
    tHardBranchingPtr partner = 
      _nasontree->particles()[particles[ix]]->colourPartner();
    for(map<ShowerParticlePtr,tHardBranchingPtr>::const_iterator
	  it=_nasontree->particles().begin();
	it!=_nasontree->particles().end();++it) {
      if(it->second==partner) {
	particles[ix]->setPartner(it->first);
      }
    }
    if(!particles[ix]->partner()) 
      throw Exception() << "Can't match partners in "
			<< "PowhegEvolver::setColourPartners()"
			<< Exception::eventerror;
  }
  // Set the initial evolution scales
  showerModel()->partnerFinder()->
    setInitialEvolutionScales(particles,!hard,false);
}

bool PowhegEvolver::startTimeLikeShower() {
  if(_nasontree) {
    map<ShowerParticlePtr,tHardBranchingPtr>::const_iterator 
      eit=_nasontree->particles().end(),
      mit = _nasontree->particles().find(progenitor()->progenitor());
    if( mit != eit && !mit->second->children().empty() ) {
      return truncatedTimeLikeShower(progenitor()->progenitor(), mit->second );
    }
  }
  return  _hardonly ? false :
    timeLikeShower(progenitor()->progenitor()) ;
}

bool PowhegEvolver::startSpaceLikeShower(PPtr parent) {
  if(_nasontree) {
    map<ShowerParticlePtr,tHardBranchingPtr>::const_iterator 
      eit =_nasontree->particles().end(),
      mit = _nasontree->particles().find(progenitor()->progenitor());
    if( _nasontree && mit != eit && mit->second->parent() ) {
      return truncatedSpaceLikeShower( progenitor()->progenitor(),
				       parent, mit->second->parent() );
    } 
  }
  return  _hardonly ? false :
    spaceLikeShower(progenitor()->progenitor(),parent);
}
