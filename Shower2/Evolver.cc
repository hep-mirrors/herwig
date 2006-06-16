// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Evolver class.
//

#include "Evolver.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/RefVector.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "Evolver.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Shower2/Kinematics/KinematicsReconstructor.h"
#include "Herwig++/Shower2/Kinematics/ShowerParticle.h"
#include "Herwig++/Shower2/Kinematics/QtildaShowerKinematics1to2.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "Herwig++/Shower2/MECorrections/MECorrectionBase.h"
#include "ThePEG/PDF/BeamParticleData.h"
#include "Herwig++/Utilities/EnumParticles.h"
#include "Herwig++/Hadronization/Remnant.h"
#include "ThePEG/Utilities/Timer.h"
#include "ShowerTree.h"
#include "ShowerProgenitor.h"

using namespace Herwig;

void Evolver::persistentOutput(PersistentOStream & os) const {
  os << _splittingGenerator << _partnerFinder << _mecorrections 
     << _showerVariables << _kinematicsrecon;
}

void Evolver::persistentInput(PersistentIStream & is, int) {
  is >> _splittingGenerator >> _partnerFinder >> _mecorrections 
     >> _showerVariables >> _kinematicsrecon;
}

void Evolver::doinit() throw(InitException) {
  Interfaced::doinit();
  _splittingGenerator->setShowerVariables(_showerVariables);
  _partnerFinder->setShowerVariables(_showerVariables);
}

void Evolver::doinitrun() {
  Interfaced::doinitrun();
  for(unsigned int ix=0;ix<_mecorrections.size();++ix)
    {_mecorrections[ix]->showerVariables(_showerVariables);}
  _kinematicsrecon->showerVariables(_showerVariables);
}
ClassDescription<Evolver> Evolver::initEvolver;
// Definition of the static class description member.

void Evolver::Init() {

  static ClassDocumentation<Evolver> documentation
    ("This class is responsible for carrying out the showering,",
     "including the kinematics reconstruction, in a given scale range.");

  static Reference<Evolver,SplittingGenerator> 
    interfaceSplitGen("SplittingGenerator", 
		      "A reference to the SplittingGenerator object", 
		      &Herwig::Evolver::_splittingGenerator,
		      false, false, true, false);

  static Reference<Evolver,PartnerFinder> 
    interfacePartnerFinder("PartnerFinder", 
                           "A reference to the PartnerFinder object", 
                           &Herwig::Evolver::_partnerFinder,
			   false, false, true, false);

  static RefVector<Evolver,MECorrectionBase> interfaceMECorrection
    ("MECorrection",
     "A vector of references to the classes that handle the matrix element"
     " correction for specific processes",
     &Evolver::_mecorrections, -1, false, false, true, false, false);

  static Reference<Evolver,KinematicsReconstructor> interfaceKinematicsReconstructor
    ("KinematicsReconstructor",
     "Reference to the object which performs the kinematic reconstruction",
     &Evolver::_kinematicsrecon, false, false, true, false, false);

  static Reference<Evolver,ShowerVariables> 
    interfaceShowerVariables("ShowerVariables", 
			       "A reference to the ShowerVariables object", 
			       &Herwig::Evolver::_showerVariables,
			       false, false, true, false);

}

void Evolver::showerHardProcess(ShowerTreePtr hard)
{
  Timer<1020> timera("Evolver::showerHardProcess");
  // set the current tree
  _currenttree=hard;
  // generate the hard matrix element correction if needed
  hardMatrixElementCorrection();
  // put all the particles into a data structure which has the particles
  // and the maximum pt for emission from the particle
  vector<ShowerProgenitorPtr> particlesToShower;
  vector<ShowerParticlePtr> particles;
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cit;
  // incoming particles
  for(cit=hard->incomingLines().begin();cit!=hard->incomingLines().end();++cit)
    {
      particlesToShower.push_back(((*cit).first));
      particles.push_back(particlesToShower.back()->progenitor());
    }
  assert(particles.size()==2);
  // outgoing particles
  for(cit=hard->outgoingLines().begin();cit!=hard->outgoingLines().end();++cit)
    {
      particlesToShower.push_back(((*cit).first));
      particles.push_back(particlesToShower.back()->progenitor());
    }
   // Set the initial evolution scales
   if(_splittingGenerator->isInteractionON(ShowerIndex::QCD))
     _partnerFinder->setQCDInitialEvolutionScales(particles);
   if(_splittingGenerator->isInteractionON(ShowerIndex::QED))
     _partnerFinder->setQEDInitialEvolutionScales(particles);
   if(_splittingGenerator->isInteractionON(ShowerIndex::EWK))
     _partnerFinder->setEWKInitialEvolutionScales(particles);
   // initial-state radiation
   if(_splittingGenerator->isISRadiationON())
     {
       for(unsigned int ix=0;ix<particlesToShower.size();++ix)
	 {
	   // only consider initial-state particles
	   if(particlesToShower[ix]->progenitor()->isFinalState()) continue;
  	  // get the PDF
	   Ptr<BeamParticleData>::const_pointer 
	     beam=dynamic_ptr_cast<Ptr<BeamParticleData>::const_pointer>
	     (particlesToShower[ix]->progenitor()->parents()[0]->dataPtr());
	   if(!beam) throw Exception() << "The Beam particle does not have"
				       << " BeamParticleData in Evolver::" 
				       << "showerhardProcess()" 
				       << Exception::runerror;
	   _showerVariables->setBeamParticle(beam);
	   _showerVariables->setCurrentPDF(beam->pdf());
	   // perform the shower
	   _progenitor=particlesToShower[ix];
	   _progenitor->hasEmitted(spaceLikeShower(particlesToShower[ix]->progenitor()));
	 }
     }
   // final-state radiation
   if(_splittingGenerator->isFSRadiationON())
     {
       for(unsigned int ix=0;ix<particlesToShower.size();++ix)
	 {
	   // only consider final-state particles
	   if(!particlesToShower[ix]->progenitor()->isFinalState()) continue;
	   // perform shower
	   _progenitor=particlesToShower[ix];
	   _progenitor->hasEmitted(timeLikeShower(particlesToShower[ix]->progenitor()));
	 }
     }
   // if needed do the kinematic reconstruction
   if(!_kinematicsrecon->reconstructHardJets(hard)) throw Veto();
}

void Evolver::hardMatrixElementCorrection()
{
  // set me correction to null pointer
  _currentme=MECorrectionPtr();
  // set the initial enhancement factors for the soft correction
  _showerVariables->initialStateRadiationEnhancementFactor(1.);
  _showerVariables->finalStateRadiationEnhancementFactor(1.);
  // if hard matrix element switched off return
  if(!_showerVariables->MECOn()) return;
  // see if there is an appropraite matrix element correction
  for(unsigned int ix=0;ix<_mecorrections.size();++ix)
    {
      if(!_mecorrections[ix]->canHandle(_currenttree)) continue;
      if(_currentme)
	{
	  ostringstream output;
	  output << "There is more than one possible matrix"
		 << "element which could be applied for ";
	  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cit;
	  for(cit=_currenttree->incomingLines().begin();
	      cit!=_currenttree->incomingLines().end();++cit)
	    {output << cit->first->progenitor()->PDGName() << " ";}
	  output << " -> ";
	  for(cit=_currenttree->outgoingLines().begin();
	      cit!=_currenttree->outgoingLines().end();++cit)
	    {output << cit->first->progenitor()->PDGName() << " ";}
	  output << "in Evolver::hardMatrixElementCorrection()\n";
	  throw Exception() << output << Exception::runerror;
	}
      else _currentme=_mecorrections[ix];
    }
  // if no suitable me correction
  if(!_currentme) return; 
  // now apply the hard correction
  if(_showerVariables->hardMEC())
    _currentme->applyHardMatrixElementCorrection(_currenttree);
}

bool Evolver::timeLikeShower(tShowerParticlePtr particle)
{
  Timer<1005> timer("Evolver::timeLikeShower");
  bool vetoed = true;
  Branching fb;
  while (vetoed) 
    {
      vetoed = false; 
      fb=_splittingGenerator->chooseForwardBranching(*particle);
      // apply vetos if needed
      if(fb.kinematics && fb.sudakov && _currentme && _showerVariables->softMEC())
	{vetoed=_currentme->softMatrixElementVeto(_progenitor,particle,fb);}
    }
  // if no branching set decay matrix and return
  if(!fb.kinematics||!fb.sudakov)
    {
      //
      // add decay matrix stuff here
      //
      return false;
    }
  // has emitted
  // Assign the shower kinematics to the emitting particle.
  particle->setShowerKinematics(fb.kinematics);
  // Assign the splitting function to the emitting particle.
  tSplittingFnPtr splitF = fb.sudakov->splittingFn();
  assert(splitF);
  particle->setSplittingFn(splitF); 
  // For the time being we are considering only 1->2 branching
  // Create the ShowerParticle objects for the two children of
  // the emitting particle; set the parent/child relationship
  ShowerParticlePtr showerProduct1,showerProduct2;
  // if same as definition create particles, otherwise create cc
  tcPDPtr pdata[2];
  for(unsigned int ix=0;ix<2;++ix) pdata[ix]=getParticleData(fb.ids[ix+1]);
  if(particle->id()!=fb.ids[0])
    {
      for(unsigned int ix=0;ix<2;++ix)
	{
	  tPDPtr cc(pdata[ix]->CC());
	  if(cc) pdata[ix]=cc;
	}
    }
  showerProduct1 = new_ptr(ShowerParticle(pdata[0]));
  showerProduct2 = new_ptr(ShowerParticle(pdata[1]));
  ShowerParticleVector theChildren; 
  theChildren.push_back(showerProduct1); 
  theChildren.push_back(showerProduct2); 
  // some code moved to updateChildren
  particle->showerKinematics()->updateChildren(particle, theChildren);
  // In the case of splittings which involves coloured particles,
  // set properly the colour flow of the branching.
  // Notice that the methods:  ColourLine::addColoured  and
  // ColourLine::addAntiColoured  automatically set also,
  // respectively, the colourLine and antiColourLine of the 
  // Particle  object they received as argument.
  ColinePair parentColinePair = ColinePair(particle->colourLine(), 
					   particle->antiColourLine());
  ColinePair showerProduct1ColinePair = ColinePair();
  ColinePair showerProduct2ColinePair = ColinePair();
  splitF->colourConnection(parentColinePair,
			   showerProduct1ColinePair, 
			   showerProduct2ColinePair);
  if ( showerProduct1ColinePair.first )
    showerProduct1ColinePair.first->addColoured( showerProduct1 );
  if ( showerProduct1ColinePair.second ) 
    showerProduct1ColinePair.second->addAntiColoured( showerProduct1 );
  if ( showerProduct2ColinePair.first ) 
    showerProduct2ColinePair.first->addColoured( showerProduct2 );
  if ( showerProduct2ColinePair.second ) 
    showerProduct2ColinePair.second->addAntiColoured( showerProduct2 );
  particle->addChild(showerProduct1);
  particle->addChild(showerProduct2);
  // update the history if needed
  if(particle==_currenttree->getFinalStateShowerProduct(_progenitor))
    _currenttree->updateFinalStateShowerProduct(_progenitor,
						particle,theChildren);
  _currenttree->addFinalStateBranching(particle,theChildren);
  // shower the first  particle
  //
  //  need to set rho here
  //
  timeLikeShower(showerProduct1);
  // shower the second particle
  //
  //   need to set rho here
  //
  timeLikeShower(showerProduct2);
  // branching has happened
  return true;
}

void Evolver::setBackwardColour(ShowerParticlePtr &newParent,
				ShowerParticlePtr &oldParent,
				ShowerParticlePtr &otherChild) {
  ColinePair parent = ColinePair();
  ColinePair child1 = ColinePair(oldParent->colourLine(),oldParent->antiColourLine());
  ColinePair child2 = ColinePair();
  // We had a colour octet
  if(child1.first && child1.second) 
    {
      if(newParent->id() == ParticleID::g) 
	{
	  if (UseRandom::rndbool()) 
	    {
	      parent.first = child1.first;
	      child2.first = child1.second;
	      parent.second = new_ptr(ColourLine());
	      child2.second = parent.second;
	    } 
	  else 
	    {
	      parent.second = child1.second;
	      child2.second = child1.first;
	      parent.first = new_ptr(ColourLine());
	      child2.first = parent.first;
	    }
	} 
      // a 3 bar state
      else 
	{
	  if(newParent->id() < 0) 
	    {
	      parent.second = child1.second;
	      child2.second = child1.first;
	    } 
	  // colour triplet
	  else 
	    {
	      parent.first = child1.first;
	      child2.first = child1.second;
	    }
	}
    } 
  // The child is a colour triplet
  else if(child1.first) 
    { 
      // colour octet
      if(newParent->hasColour() && newParent->hasAntiColour()) 
	{ 
	  parent.first = child1.first;
	  parent.second = child2.second = new_ptr(ColourLine());
	} 
      // it must be a colour triplet, so child2 is a colour octet
      else 
	{ 
	  child2.second = child1.first;
	  child2.first = parent.first = new_ptr(ColourLine());
	}
    } 
  // The child is a 3 bar state
  else if(child1.second) 
    { 
      // colour octet
      if(newParent->hasColour() && newParent->hasAntiColour()) 
	{
	  parent.second = child1.second;
	  parent.first = child2.first = new_ptr(ColourLine());
	} 
      // it must be a colour triplet, so child2 is a colour octet
      else { 
	child2.first = child1.second;
	child2.second = parent.second = new_ptr(ColourLine());
      }
    } 
  else throw Exception() << "Evolver::setBackwardColour() "
			 << "No colour info for parton!\n"
			 << Exception::eventerror;
  if(parent.first) parent.first->addColoured(newParent);
  if(parent.second) parent.second->addAntiColoured(newParent);
  if(child2.first) child2.first->addColoured(otherChild);
  if(child2.second) child2.second->addAntiColoured(otherChild);
}

bool Evolver::spaceLikeShower(tShowerParticlePtr particle)
{
  Timer<1006> timer("Evolver::spaceLikeShower");
   bool vetoed(true);
   Branching bb;
   // generate branching
   while (vetoed)
     {
       vetoed=false;
       bb=_splittingGenerator->chooseBackwardBranching(*particle);
       // apply the soft correction
       if(bb.kinematics && bb.sudakov && _currentme && _showerVariables->softMEC())
	 {vetoed=_currentme->softMatrixElementVeto(_progenitor,particle,bb);}
     }
   if(!bb.kinematics||!bb.sudakov) return false;
   // assign the splitting function and shower kinematics
   tSplittingFnPtr splitF = bb.sudakov->splittingFn();
   assert(splitF);
   particle->setShowerKinematics(bb.kinematics);
   particle->setSplittingFn(splitF); 
   // For the time being we are considering only 1->2 branching
   // particles as in Sudakov form factor
   tcPDPtr part[2]={getParticleData(bb.ids[0]),
		    getParticleData(bb.ids[2])};
   if(particle->id()!=bb.ids[1])
     {
       if(part[0]->CC()) part[0]=part[0]->CC();
       if(part[1]->CC()) part[1]=part[1]->CC();
     }
   // Now create the actual particles, make the otherChild a final state
   // particle, while the newParent is not
   ShowerParticlePtr newParent=new_ptr(ShowerParticle(part[0]));
   newParent->setFinalState(false);
   ShowerParticlePtr otherChild = new_ptr(ShowerParticle(part[1]));
   otherChild->setFinalState(true);
   otherChild->setInitiatesTLS(true);
   // Set up the colour connections and the parent/child relationships
   createBackwardBranching(particle,newParent,otherChild,
			   particle->showerKinematics()->qtilde(),
			   particle->showerKinematics()->z(),
			   splitF->interactionType());
   // update the history if needed
   _currenttree->updateInitialStateShowerProduct(_progenitor,particle,newParent);
   _currenttree->addInitialStateBranching(particle,newParent,otherChild);
   // now continue the shower
   bool emitted=spaceLikeShower(newParent);
   // now reconstruct the momentum
   if(!emitted)
     bb.kinematics->updateLast(newParent,0);
   // update properties of children needed for branching
   // time-like child
   // the alpha decomposition variable
   double z(bb.kinematics->z());
   double alpha(newParent->sudAlpha());
   otherChild->sudAlpha((1.-z)*alpha);
   // the transverse momentum
   double cphi = cos(bb.kinematics->phi());
   double sphi = sin(bb.kinematics->phi());
   Energy pt = bb.kinematics->pT();
   Energy kx = (1.-z)*newParent->sudPx() - cphi*pt;
   Energy ky = (1.-z)*newParent->sudPy() - sphi*pt; 
   otherChild->sudPx(kx);
   otherChild->sudPy(ky);
   // space-like child
   particle->sudAlpha(newParent->sudAlpha() - otherChild->sudAlpha());
   particle->sudBeta( newParent->sudBeta()  - otherChild->sudBeta() );
   particle->sudPx(   newParent->sudPx()    - otherChild->sudPx()   );
   particle->sudPy(   newParent->sudPy()    - otherChild->sudPy()   );
   // perform the shower of the final-state particle
   timeLikeShower(otherChild);
   // return the emitted
   return true;
}

void Evolver::createBackwardBranching(ShowerParticlePtr part,
				      ShowerParticlePtr newParent,
				      ShowerParticlePtr otherChild,
				      Energy scale, double zz,
				      ShowerIndex::InteractionType inter) 
{
  // no z for angular ordering in backward branchings
  newParent->setEvolutionScale(inter, scale);
  otherChild->setEvolutionScale(inter, (1.-zz)*scale);
  // for the reconstruction of kinematics, parent/child
  // relationships are according to the branching process:
  ShowerParticleVector theChildren; 
  theChildren.push_back(newParent); 
  theChildren.push_back(otherChild); 
  part->showerKinematics()->updateChildren(part, theChildren);
  // *** set proper colour connections
  setBackwardColour(newParent,part,otherChild);
  // *** set proper parent/child relationships
  newParent->addChild(part);
  newParent->addChild(otherChild);
  newParent->x(part->x()/part->showerKinematics()->z()); 
  // Now fix the hadrons connections
  tPPtr hadron;
  if(part->parents().size() == 2) hadron = part->parents()[0];
  else throw Exception() << "Evolver::createBackwardBranching: not one parent!" 
			 << Exception::runerror; 
  hadron->abandonChild(part);
  hadron->addChild(newParent);
}

void Evolver::makeRemnants(ShowerTreePtr hard)
{
  // loop over the jets
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cit;
  for(cit = hard->incomingLines().begin();cit!=hard->incomingLines().end();++cit)
    {
      ShowerParticlePtr particle=(*cit).first->progenitor();
      // if initial-state, from hard process and
      if(!particle->isFinalState()&&particle->perturbative()==1&&
	 isISRadiationON())
	{
	  tShowerParticlePtr current(particle),next;
	  Lorentz5Momentum ptotal(current->momentum());
	  do
	    {
	      // find the parent of the particle
	      if(!current->parents().empty())
		{
		  next=dynamic_ptr_cast<tShowerParticlePtr>(current->parents()[0]);
		  if(next)
		    {
		      tParticleSet siblings=current->siblings();
		      tParticleSet::iterator sib;
		      for(sib=siblings.begin();sib!=siblings.end();++sib)
			{ptotal+=(*sib)->momentum();}
		      current=next;
		    }
		}
	      else
		next=tShowerParticlePtr();
	    }
	  while(next);
	  // now we have the total momentum of the shower and are at the top
	  // of the tree so find the remnant and update it
	  tParticleSet siblings=current->siblings();
	  tParticleSet::iterator sib;
	  for(sib=siblings.begin();sib!=siblings.end();++sib)
	    {
	      // if the remnant perform cast and update
	      if((*sib)->id()==ExtraParticleID::Remnant)
		{
		  tRemnantPtr rem=dynamic_ptr_cast<tRemnantPtr>(*sib);
		  if(rem)
		    {rem->regenerate(current,current->parents()[0]->momentum()-ptotal);}
		}
	    }
	}
    }
  tParticleSet remnantSet
    = generator()->currentEventHandler()->currentCollision()->getRemnants();
  tParticleSet::iterator rem;
}

void Evolver::showerDecay(ShowerTreePtr decay)
{
  Timer<1004> timer("Evolver::showerDecay");
  _currenttree=decay;
  // put all the particles into a data structure which has the particles
   // and the maximum pt for emission from the particle
   vector<ShowerProgenitorPtr> particlesToShower;
   vector<ShowerParticlePtr> particles;
   map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cit;
   // incoming particles
   for(cit=decay->incomingLines().begin();cit!=decay->incomingLines().end();++cit)
     {
       particlesToShower.push_back(((*cit).first));
       particles.push_back(particlesToShower.back()->progenitor());
     };
   assert(particles.size()==1);
   // outgoing particles
   for(cit=decay->outgoingLines().begin();cit!=decay->outgoingLines().end();++cit)
     {
       particlesToShower.push_back(((*cit).first));
       particles.push_back(particlesToShower.back()->progenitor());
     }
   // generate the hard matrix element correction if needed
   hardMatrixElementCorrection();
   // Set the initial evolution scales
   if(_splittingGenerator->isInteractionON(ShowerIndex::QCD))
     _partnerFinder->setQCDInitialEvolutionScales(particles,true);
   if(_splittingGenerator->isInteractionON(ShowerIndex::QED))
     _partnerFinder->setQEDInitialEvolutionScales(particles,true);
   if(_splittingGenerator->isInteractionON(ShowerIndex::EWK))
     _partnerFinder->setEWKInitialEvolutionScales(particles,true);
   // initial-state radiation
   if(_splittingGenerator->isISRadiationON())
     {
       // compute the minimum mass of the final-state
       Energy minmass(0.);
       for(unsigned int ix=0;ix<particlesToShower.size();++ix)
	 {if(particlesToShower[ix]->progenitor()->isFinalState())
	     minmass+=particlesToShower[ix]->progenitor()->mass();}
       for(unsigned int ix=0;ix<particlesToShower.size();++ix)
	 {
	   // only consider initial-state particles
	   if(particlesToShower[ix]->progenitor()->isFinalState()) continue;
	   // perform shower
	   _progenitor=particlesToShower[ix];
	   // set the scales correctly. The current scale is the maximum scale for
	   // emission not the starting scale
	   vector<Energy> maxscale=_progenitor->progenitor()->evolutionScales();
	   Energy startScale=_progenitor->progenitor()->mass();
	   _progenitor->progenitor()->setEvolutionScale(ShowerIndex::QCD,startScale);
	   _progenitor->progenitor()->setEvolutionScale(ShowerIndex::QED,startScale);
	   _progenitor->progenitor()->setEvolutionScale(ShowerIndex::EWK,startScale);
	   // perform the shower
	   _progenitor->hasEmitted(spaceLikeDecayShower(_progenitor->progenitor(),
	   						maxscale,minmass));
	   //_progenitor->hasEmitted(false);
	 }
     }
   // final-state radiation
   if(_splittingGenerator->isFSRadiationON())
     {
       for(unsigned int ix=0;ix<particlesToShower.size();++ix)
	 {
	   // only consider final-state particles
	   if(!particlesToShower[ix]->progenitor()->isFinalState()) continue;
	   // perform shower
	   _progenitor=particlesToShower[ix];
	   _progenitor->hasEmitted(timeLikeShower(particlesToShower[ix]->progenitor()));
	 }
     }
   // if needed do the kinematic reconstruction
   if(!_kinematicsrecon->reconstructDecayJets(decay)) throw Veto();
}

bool Evolver::spaceLikeDecayShower(tShowerParticlePtr particle,vector<Energy> maxscale,
				   Energy minmass)
{
  bool vetoed = true;
  Branching fb;
  while (vetoed) 
    {
      vetoed = false;
      fb=_splittingGenerator->chooseDecayBranching(*particle,maxscale,minmass);
      // apply the soft correction
      if(fb.kinematics && fb.sudakov && _currentme && _showerVariables->softMEC())
	{vetoed=_currentme->softMatrixElementVeto(_progenitor,particle,fb);}
    }
  // if no branching set decay matrix and return
  if(!fb.kinematics||!fb.sudakov)
    {
      //
      // add decay matrix stuff here
      //
      return false;
    }
  // has emitted
  // Assign the shower kinematics to the emitting particle.
  particle->setShowerKinematics(fb.kinematics);
  // Assign the splitting function to the emitting particle.
  tSplittingFnPtr splitF = fb.sudakov->splittingFn();
  assert(splitF);
  particle->setSplittingFn(splitF); 
  // For the time being we are considering only 1->2 branching
  // Create the ShowerParticle objects for the two children of
  // the emitting particle; set the parent/child relationship
  ShowerParticlePtr showerProduct1,showerProduct2;
  // if same as definition create particles, otherwise create cc
  tcPDPtr pdata[2];
  for(unsigned int ix=0;ix<2;++ix) pdata[ix]=getParticleData(fb.ids[ix+1]);
  if(particle->id()!=fb.ids[0])
    {
      for(unsigned int ix=0;ix<2;++ix)
 	{
 	  tPDPtr cc(pdata[ix]->CC());
 	  if(cc) pdata[ix]=cc;
 	}
    }
  showerProduct1 = new_ptr(ShowerParticle(pdata[0]));
  showerProduct2 = new_ptr(ShowerParticle(pdata[1]));
  ShowerParticleVector theChildren; 
  theChildren.push_back(showerProduct1); 
  theChildren.push_back(showerProduct2); 
  // some code moved to updateChildren
  particle->showerKinematics()->updateChildren(particle, theChildren);
  // In the case of splittings which involves coloured particles,
  // set properly the colour flow of the branching.
  // Notice that the methods:  ColourLine::addColoured  and
  // ColourLine::addAntiColoured  automatically set also,
  // respectively, the colourLine and antiColourLine of the 
  // Particle  object they received as argument.
  ColinePair parentColinePair = ColinePair(particle->colourLine(), 
					   particle->antiColourLine());
  ColinePair showerProduct1ColinePair = ColinePair();
  ColinePair showerProduct2ColinePair = ColinePair();
  splitF->colourConnection(parentColinePair,
			   showerProduct1ColinePair, 
			   showerProduct2ColinePair);
  if ( showerProduct1ColinePair.first )
    showerProduct1ColinePair.first->addColoured( showerProduct1 );
  if ( showerProduct1ColinePair.second ) 
    showerProduct1ColinePair.second->addAntiColoured( showerProduct1 );
  if ( showerProduct2ColinePair.first ) 
    showerProduct2ColinePair.first->addColoured( showerProduct2 );
  if ( showerProduct2ColinePair.second ) 
    showerProduct2ColinePair.second->addAntiColoured( showerProduct2 );
  particle->addChild(showerProduct1);
  particle->addChild(showerProduct2);
  // update the history if needed
  _currenttree->updateInitialStateShowerProduct(_progenitor,particle,showerProduct1);
  _currenttree->addInitialStateBranching(particle,showerProduct1,showerProduct2);
  // shower the first  particle
  //
  //  need to set rho here
  //
  spaceLikeDecayShower(showerProduct1,maxscale,minmass);
  // shower the second particle
  //
  //   need to set rho here
  //
  timeLikeShower(showerProduct2);
  // branching has happened
  return true;
}
