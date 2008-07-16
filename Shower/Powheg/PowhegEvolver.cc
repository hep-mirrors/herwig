// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PowhegEvolver class.
//

#include "PowhegEvolver.h"
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
#include "Herwig++/Shower/Powheg/HardGenerators/DefaultEmissionGenerator.h"
#include "Herwig++/Shower/Default/FS_QtildaShowerKinematics1to2.h"
#include "HardTree.h"
#include "Herwig++/Utilities/Histogram.h"


using namespace Herwig;

void PowhegEvolver::doinit() throw(InitException) {
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
    ("The PowhegEvolver implements the POWHEG approach to MC@NLO");

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


void PowhegEvolver::showerDecay(ShowerTreePtr tree) {
  // set the tree
  currentTree(tree);
  // set up the shower
  vector<ShowerProgenitorPtr> particlesToShower = setupShower(false);
  setupMaximumScales(tree, particlesToShower);
  // main showering loop
  unsigned int ntry(0);
  do {
    // clear results of last attempt
    if(ntry!=0) {
      currentTree()->clear();
      setColourPartners(false);
    }
    // initial-state radiation
//     if(_splittingGenerator->isISRadiationON()) {
//       // compute the minimum mass of the final-state
//       Energy minmass(0.);
//       for(unsigned int ix=0;ix<particlesToShower.size();++ix)
// 	{if(particlesToShower[ix]->progenitor()->isFinalState())
// 	    minmass+=particlesToShower[ix]->progenitor()->mass();}
//       for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
// 	// only consider initial-state particles
// 	if(particlesToShower[ix]->progenitor()->isFinalState()) continue;
// 	// perform shower
// 	_progenitor=particlesToShower[ix];
// 	// set the scales correctly. The current scale is the maximum scale for
// 	// emission not the starting scale
// 	vector<Energy> maxscale=progenitor()->progenitor()->evolutionScales();
// 	Energy startScale=progenitor()->progenitor()->mass();
// 	progenitor()->progenitor()->setEvolutionScale(ShowerIndex::QCD,startScale);
// 	progenitor()->progenitor()->setEvolutionScale(ShowerIndex::QED,startScale);
// 	progenitor()->progenitor()->setEvolutionScale(ShowerIndex::EWK,startScale);
// 	// perform the shower
// 	progenitor()->hasEmitted(spaceLikeDecayShower(progenitor()->progenitor(),
// 						     maxscale,minmass)); 
//       }
//     }
    // final-state radiation
    map<ShowerParticlePtr,tHardBranchingPtr>::const_iterator mit,eit;
    if(_nasontree) eit=_nasontree->particles().end();
    if(isFSRadiationON()) {
      for( unsigned int ix=0; ix < particlesToShower.size(); ++ix ) {
        // only consider final-state particles
	if( !particlesToShower[ix]->progenitor()->isFinalState() ) continue;
	// perform shower
	progenitor( particlesToShower[ix] );
	if( _nasontree ) {
	  mit = _nasontree->particles().find( progenitor()->progenitor() );
	  if( mit != eit && !mit->second->children().empty() ) {
	    progenitor()->hasEmitted( truncatedTimeLikeShower(
			  particlesToShower[ix]->progenitor(), mit->second ) );
	  } 
	  else {
	    progenitor()->hasEmitted( _hardonly ? false :
		timeLikeShower( particlesToShower[ix]->progenitor() ) );
	  }
	}
	else {
	  progenitor()->hasEmitted( _hardonly ? false :
               timeLikeShower( particlesToShower[ix]->progenitor() ) );
	}
      }
    }
  }
  while( !showerModel()->kinematicsReconstructor()->reconstructDecayJets(tree)&&
	maximumTries()>++ntry);
  if(maximumTries()==ntry) 
    throw Exception() << "Failed to generate the shower after "
		      << ntry << " attempts in PowhegEvolver::showerDecay()"
		      << Exception::eventerror;
  currentTree()->hasShowered(true);

  // test that generated pt after IMR and MR  matches original n+1 of hardest
  if(_hardonly && _nasontree) {
    // extract the particles from end point of the shower
    map<int,PPtr> outgoing;
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
    static Energy eps = 1e-5*GeV;
    for(unsigned int ix=0;ix<outb.size();++ix) {
      Lorentz5Momentum pdiff = outb[ix]->momentum();
      pdiff -=outgoing[outb[ix]->id()]->momentum();
     
      if(abs(pdiff.x()) > eps ||
	 abs(pdiff.y()) > eps ||
	 abs(pdiff.z()) > eps ||
	 abs(pdiff.t()) > eps ) {
	generator()->log() << "testing recon problem " << *outgoing[outb[ix]->id()] << "\n"
			   << "                      " << pdiff/GeV << "\n";
	generator()->log() << "testing " << abs(pdiff.x())/eps  << " "
			   << abs(pdiff.y())/eps  << " "
			   << abs(pdiff.z())/eps  << " "
			   << abs(pdiff.t())/eps  << "\n";
      }
    }
  }
}

void PowhegEvolver::showerHardProcess(ShowerTreePtr tree) {
  // set the tree
  currentTree(tree);
  // set up the shower
  vector<ShowerProgenitorPtr> particlesToShower = setupShower(true);
  setupMaximumScales( tree, particlesToShower );
  // generate the intrinsic p_T once and for all
  generateIntrinsicpT(particlesToShower);
  unsigned int ntry(0);
  do {
    // clear results of last attempt and reset colour partners
    if(ntry!=0) {
      currentTree()->clear();
      setColourPartners(true);
    }
    // initial-state radiation
    map<ShowerParticlePtr,tHardBranchingPtr>::const_iterator mit,eit;
    if(_nasontree) eit=_nasontree->particles().end();
    if(isISRadiationON()) {
      for( unsigned int ix=0; ix < particlesToShower.size(); ++ix ) {
	// only consider initial-state particles
	if( particlesToShower[ix]->progenitor()->isFinalState() ) continue;
	// get the PDF
	setBeamParticle( particlesToShower[ix]->beam() );
	if( !beamParticle() ) throw Exception() << "The Beam particle does not have"
						<< " BeamParticleData in Evolver::" 
						<< "showerhardProcess()" 
						<< Exception::runerror;
	// perform the shower
	progenitor( particlesToShower[ix] );
	//vgbbb
	if( _nasontree ) mit = _nasontree->
			  particles().find( progenitor()->progenitor() );
	if( _nasontree && mit != eit && mit->second->parent() ) {
	  progenitor()->
	    hasEmitted( truncatedSpaceLikeShower( particlesToShower[ix]->progenitor(),
						  particlesToShower[ix]->
						  original()->parents()[0],
						  mit->second->parent() ) );
	}
	else {
	  progenitor()->hasEmitted( _hardonly ? false :
				    spaceLikeShower( particlesToShower[ix]->progenitor(),
						     particlesToShower[ix]->
						     original()->parents()[0] ) );
	}
      }
    }
    // final-state radiation
    if( isFSRadiationON() ) {
      for( unsigned int ix = 0; ix < particlesToShower.size(); ++ix ) {
	// only consider final-state particles
	if( !particlesToShower[ix]->progenitor()->isFinalState() ) continue;
	// perform shower
	progenitor( particlesToShower[ix] );
	progenitor()->hasEmitted( _hardonly ? false :
				  timeLikeShower( particlesToShower[ix]->progenitor() ) );
      }
    }
  }
  while(!showerModel()->kinematicsReconstructor()->
	reconstructHardJets( tree, intrinsicpT() ) &&
 	maximumTries() > ++ntry );
  if( maximumTries() == ntry ) 
    throw Exception() << "Failed to generate thPersistency shower after " << ntry 
		      << " attempts in PowhegEvolver::showerHardProcess()"
		      << Exception::eventerror;
  currentTree()->hasShowered( true );
  // test the momenta are the same
  if( _hardonly && _nasontree ) {
    // extract the particles from end point of the shower
    list<ShowerParticlePtr> partons;
    for(unsigned int ix = 0; ix < particlesToShower.size(); ++ix) {
      if(  particlesToShower[ix]->progenitor()->isFinalState() ) continue;
      ShowerParticlePtr incoming=particlesToShower[ix]->progenitor();
      if(!incoming->parents().empty()&&incoming->parents()[0]->dataPtr()->coloured()) {
	ShowerParticlePtr incoming2=dynamic_ptr_cast<ShowerParticlePtr>(incoming->parents()[0]);
	if(incoming2) {
	  partons.push_back(incoming2);
	  for(unsigned int iy=0;iy<incoming2->children().size();++iy) {
	    ShowerParticlePtr child = dynamic_ptr_cast<ShowerParticlePtr>(incoming2->children()[iy]);
	    if(child&&child!=incoming) partons.push_back(child);
	  }
	}
	else {
	  partons.push_back(incoming);
	}
      }
      else {
	partons.push_back(incoming);
      }
    }
    // extract the particles from the nason tree
    list<ShowerParticlePtr> nason;
    for(set<HardBranchingPtr>::const_iterator it = _nasontree->incoming().begin();
	it != _nasontree->incoming().end(); ++it) {
      nason.push_back((*it)->branchingParticle());
      for(unsigned int ix=0;ix<(*it)->children().size();++ix) {
	if((*it)->children()[ix]->incoming()) continue;
	nason.push_back((*it)->children()[ix]->branchingParticle());
      }
    }
    static Energy eps = 1e-4*GeV;
    list<ShowerParticlePtr>::iterator it1,it2,it3;
    for(it1=partons.begin();it1!=partons.end();++it1) {
      Energy2 delta(1e30*GeV2);
      it3=nason.end();
      for(it2=nason.begin();it2!=nason.end();++it2) {
	if((**it1).isFinalState()!=(**it2).isFinalState()) continue;
	if((**it1).id()!=(**it2).id()) continue;
	Lorentz5Momentum pdiff=(**it1).momentum()-(**it2).momentum();
	Energy2 test = sqr(pdiff.x())+sqr(pdiff.y())+sqr(pdiff.z())+sqr(pdiff.t());
	if(test<delta) {
	  delta=test;
	  it3=it2;
	}
      }
      if(it3!=nason.end()) {
	Lorentz5Momentum pdiff=(**it1).momentum()-(**it3).momentum();
	if(abs(pdiff.x()) > eps || abs(pdiff.y()) > eps ||
	   abs(pdiff.z()) > eps || abs(pdiff.t()) > eps ) {
	  generator()->log() << "testing recon problem " << **it1 << "\n"
			     << "                      " << **it3 << "\n"
			     << "                      " << pdiff/GeV << "\n";
	  nason.erase(it3);
	} 
      }
      else {
	generator()->log() << "recon problem - no match for \n"
			   << **it1;
      }
    }
  }
}

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
  // remake the colour partners if needed
  if(currentTree()->hardMatrixElementCorrection()) {
    setColourPartners(hard);
    currentTree()->resetShowerProducts();
  }
  return particlesToShower;
}

void PowhegEvolver::hardestEmission() {
  // see if there is an appropriate hard emission generator
  HardestEmissionGeneratorPtr currenthard = HardestEmissionGeneratorPtr();
  for( unsigned int ix = 0; ix < _hardgenerator.size(); ++ix ) {
    if( !_hardgenerator[ix]->canHandle( currentTree() ) ) continue;
    if( currenthard ) {
      if( dynamic_ptr_cast< DefaultEmissionGeneratorPtr > ( currenthard ) ) {
	currenthard = _hardgenerator[ix];
      }
      else if( !dynamic_ptr_cast< DefaultEmissionGeneratorPtr >( _hardgenerator[ix] ) ) {
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
    }
    else {
      currenthard=_hardgenerator[ix];
    }
  }
  // if no suitable generator return
  _nasontree=HardTreePtr();
  if(!currenthard) return;
  // generate the hardest emission
  _nasontree = currenthard->generateHardest( currentTree() );
}

bool PowhegEvolver::truncatedTimeLikeShower( tShowerParticlePtr particle,
					    HardBranchingPtr branch ) {
  bool vetoed = true;
  Branching fb;
  unsigned int iout;
  tcPDPtr pdata[2];
  while (vetoed) {
    iout=0;
    vetoed = !isTruncatedShowerON();
    fb=splittingGenerator()->chooseForwardBranching(*particle,1.);
    // check haven't evolved too far
    if(!fb.kinematics||fb.kinematics->scale()<branch->scale()) {
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
    if(pdata[0]->id()!=pdata[1]->id()) {
      if(pdata[0]->id()==particle->id())       iout=1;
      else if (pdata[1]->id()==particle->id()) iout=2;
    }
    else if(pdata[0]->id()==particle->id()) {
      if(fb.kinematics->z()>0.5) iout=1;
      else                       iout=2;
    }
    // apply the vetos for the truncated shower
    double zsplit(0.);
    if(iout==0) vetoed=true;
    else {
      zsplit = iout==1 ? fb.kinematics->z() : 1-fb.kinematics->z();
    }
    // hardest line veto
    if(zsplit<0.5) vetoed=true;
    // angular ordering veto
    if(fb.kinematics->scale()<branch->scale()/zsplit) vetoed=true;
    // pt veto
    if(fb.kinematics->pT()>progenitor()->maximumpT()) vetoed = true;
    // no truncated shower if only generating hardest emission
    if(_hardonly) vetoed = true;
    // if vetoed reset scale
    if(vetoed) particle->setEvolutionScale(ShowerIndex::QCD,fb.kinematics->scale());
  }
  // if no branching set decay matrix and return
  if(!fb.kinematics) {
    // construct the kinematics for the hard emission
    ShoKinPtr showerKin=
      branch->sudakov()->createFinalStateBranching(branch->scale(),
						   branch->children()[0]->z(),
						   branch->phi(),
						   branch->children()[0]->pT());
    particle->setEvolutionScale( ShowerIndex::QCD,branch->scale() );
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
    //
    //  need to set rho here
    //
    if( branch->children()[0]->children().empty() ) {
      timeLikeShower(theChildren[0]);
    }
    else {
      truncatedTimeLikeShower( theChildren[0],branch->children()[0] );
    } 
    // shower the second particle
    //
    //   need to set rho here
    //
    if( branch->children()[1]->children().empty() ) {
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
  //
  //  need to set rho here
  //
  if( iout == 1 ) truncatedTimeLikeShower( theChildren[0], branch );
  else        timeLikeShower( theChildren[0] );
  // shower the second particle
  //
  //   need to set rho here
  //
  if( iout == 2 ) truncatedTimeLikeShower( theChildren[1], branch );
  else        timeLikeShower( theChildren[1] );
  // branching has happened
  return true;
}

bool PowhegEvolver::truncatedSpaceLikeShower(tShowerParticlePtr particle, PPtr beam,
					    HardBranchingPtr branch) {
  bool vetoed(true);
  Branching bb;
  // generate branching
  tcPDPtr part[2];
  while (vetoed) {
    if( isTruncatedShowerON() ) vetoed = false;
    bb = splittingGenerator()->chooseBackwardBranching( *particle, beam, 1., beamParticle() );
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
    // if particle changes type
    if( part[0]->id() != particle->id() ) vetoed=true;
    // if doesn't carry most of momentum
    if( zsplit < 0.5 ) vetoed=true;
    // pt veto
    if( bb.kinematics->pT() > progenitor()->maximumpT() ) vetoed = true;
    //hardest only switch
    if( _hardonly ) vetoed = true;
    // angular ordering veto
    if( bb.kinematics->scale() < branch->scale() / zsplit ) vetoed = true;
    //reset scale if vetoed
    if( vetoed ) particle->setEvolutionScale( ShowerIndex::QCD, bb.kinematics->scale() );
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
	kinematics->updateLast( newParent, 0.*MeV, 0.*MeV );
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
      bb.kinematics->updateLast( newParent, 0.*MeV, 0.*MeV );
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
	particles[ix]->setPartner(ShowerIndex::QCD,it->first);
      }
    }
    if(!particles[ix]->partners()[ShowerIndex::QCD]) 
      throw Exception() << "Can't match partners in "
			<< "PowhegEvolver::setColourPartners()"
			<< Exception::eventerror;
  }
  // Set the initial evolution scales
  if(splittingGenerator()->isInteractionON(ShowerIndex::QCD))
    showerModel()->partnerFinder()->
      setQCDInitialEvolutionScales(particles,!hard,false);
  if(splittingGenerator()->isInteractionON(ShowerIndex::QED))
    showerModel()->partnerFinder()->setQEDInitialEvolutionScales(particles,!hard);
  if(splittingGenerator()->isInteractionON(ShowerIndex::EWK))
    showerModel()->partnerFinder()->setEWKInitialEvolutionScales(particles,!hard);
}
