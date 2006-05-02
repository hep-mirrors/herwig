// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ForwardEvolver class.
//

#include "ForwardEvolver.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Reference.h" 
#include "Herwig++/Utilities/HwDebug.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ShowerParticle.h"
#include "ShowerKinematics.h"
#include "QtildaShowerKinematics1to2.h"
#include "Herwig++/Shower/SplittingFunctions/SplittingFunction.h"
#include "ThePEG/EventRecord/ColourLine.h"
#include "ShowerConfig.h"
#include "QQbarG.h"

using namespace Herwig;
using namespace ThePEG;

ForwardEvolver::~ForwardEvolver() {}


void ForwardEvolver::persistentOutput(PersistentOStream & os) const {
  os << _splittingGenerator;
}


void ForwardEvolver::persistentInput(PersistentIStream & is, int) {
  is >> _splittingGenerator;
}


ClassDescription<ForwardEvolver> ForwardEvolver::initForwardEvolver;
// Definition of the static class description member.


void ForwardEvolver::Init() {

  static ClassDocumentation<ForwardEvolver> documentation
    ("This class is responsible for the forward showering of time-like particles",
     "It does also the special forward showering for decaying particles",
     "in which the angular ordering is reversed");

  static Reference<ForwardEvolver,SplittingGenerator> 
    interfaceSplitGen("SplittingGenerator", 
		      "A reference to the SplittingGenerator object", 
		      &Herwig::ForwardEvolver::_splittingGenerator,
		      false, false, true, false);

}

//-----------------------------------------------------------------------------

bool ForwardEvolver::timeLikeShower(tEHPtr ch, 
				    const tShowerVarsPtr showerVariables, 
				    tShowerParticlePtr particle, 
				    ShowerParticleVector & collecShoPar,
				    const bool specialDecay )  
  throw (Veto, Stop, Exception) {

  bool hasEmitted = false;
  Energy ptMax=0.0*GeV;
  // maximum pt for emmision
  if(particle->dataPtr()->iColour()==PDT::Colour3) {
    ptMax = showerVariables->largestPtQ();
  } else if(particle->dataPtr()->iColour()==PDT::Colour3bar) {
    ptMax = showerVariables->largestPtQbar();
  }
  // this isn't used at the moment, needed for top etc???
  particle->undecay();
  tShowerParticleVector particlesYetToShower;
  particlesYetToShower.push_back(particle);

  //if(particle->id()==21) return false;
  //return false;
  
  do 
    {
      tShowerParticlePtr part = particlesYetToShower.back();
      particlesYetToShower.pop_back();
      bool vetoed=false;
      Branching fb;
      do 
	{
	  vetoed=false;
	  // select branching
	  fb=_splittingGenerator->chooseForwardBranching(ch, *part, specialDecay);
	  // apply veto to emission for me correction if needed
	  if(showerVariables->softMEC()&&fb.first && fb.second &&
	     fb.second->splittingFn()->interactionType() == ShowerIndex::QCD &&
	     abs(part->id()) < 7&&part->id() == particle->id())
	    {
	      double d_z = fb.first->z();
	      Energy d_qt = fb.first->qtilde();
	      Energy2 d_m2 = part->momentum().m2();
	      Energy pPerp = (1.-d_z)*sqrt( sqr(d_z*d_qt) - d_m2);
	      if(pPerp>ptMax)
		{
		  ptMax = pPerp; 
		  vetoed = MEVeto(particle, d_qt, d_z);
		  if(vetoed)
		    {
		      part->setEvolutionScale(fb.second->splittingFn()->interactionType(),
					      fb.first->qtilde());
		    }
		}
	    }
	}
      while(vetoed);
      // if no branching
      if(fb.first == ShoKinPtr() || fb.second == tSudakovPtr()) continue;

      hasEmitted = true;
      _splittingGenerator->generateBranchingKinematics(ch, *part, fb);
      // Assign the splitting function and the shower kinematics
      // to the emitting particle.
      part->setShowerKinematics(fb.first);
      part->setSplittingFn(fb.second->splittingFn()); 
      // For the time being we are considering only 1->2 branching
      tSplittingFnPtr splitF = fb.second->splittingFn();
      if(!splitF) throw Exception() << "Must have a splitting function !!!!! "
				    << Exception::runerror;
      // Create the ShowerParticle objects for the two children of
      // the emitting particle; set the parent/child relationship;
      // add them to the collecShoPar and particlesYetToShower
      // collections.  Remember that we our splitting functions
      // are associated only with particles (id>0) and not with
      // antiparticles (id<0)
      // (because we are assuming CP-conserving vertices):
      // therefore the signs of the decays products must be set by
      // hand explicitly; to simplify this, we assume that the
      // first product must always have the sign as the parent:
      // q->q+g, qbar->qbar+g, g->q+qbar, etc.
      // Notice that the momenta of the shower products is not
      // set: only at the end of the showering, during the
      // kinematics reconstruction such momenta are calculated and
      // set.
      
      ShowerParticlePtr showerProduct1;
      ShowerParticlePtr showerProduct2;
      // if same as definition create particles, otherwise create cc
      // generalized by PR 
      if(part->id()==fb.third[0]) {
	showerProduct1 = new_ptr(ShowerParticle
				 (getParticleData(fb.third[1])));
	showerProduct2 = new_ptr(ShowerParticle
				 (getParticleData(fb.third[2])));
      } else { 
	tPDPtr cc(getParticleData(fb.third[1])->CC());
	if(cc)showerProduct1 = new_ptr(ShowerParticle(cc));
	else showerProduct1 = new_ptr(ShowerParticle
				      (getParticleData(fb.third[1])));
	cc=getParticleData(fb.third[2])->CC();
	if(cc) showerProduct2 = new_ptr(ShowerParticle(cc));
	else showerProduct2 = new_ptr(ShowerParticle
				      (getParticleData(fb.third[2])));
      }
      
      const ShowerIndex::InteractionType interaction = splitF->interactionType();
      const Energy scale = part->showerKinematics()->qtilde();
      double zz = dynamic_ptr_cast<Ptr<QtildaShowerKinematics1to2>::
	transient_pointer>(part->showerKinematics())->z();
      // note that 1st child gets z, 2nd gets (1-z) by our convention.
      showerProduct1->setEvolutionScale(interaction, zz*scale);
      showerProduct2->setEvolutionScale(interaction, (1.-zz)*scale);
      showerProduct1->setInitiatesTLS(false);
      showerProduct2->setInitiatesTLS(false);
      
      ShowerParticleVector theChildren; 
      theChildren.push_back(showerProduct1); 
      theChildren.push_back(showerProduct2); 
      part->showerKinematics()->updateChildren(part, theChildren); 
      
      // In the case of splittings which involves coloured particles,
      // set properly the colour flow of the branching.
      // Notice that the methods:  ShowerColourLine::addColoured  and
      // ShowerColourLine::addAntiColoured  automatically set also,
      // respectively, the colourLine and antiColourLine of the 
      // ShowerParticle  object they received as argument.
      ShoColinePair parentShoColinePair = 
	ShoColinePair(part->colourLine(), part->antiColourLine());
      ShoColinePair showerProduct1ShoColinePair = ShoColinePair();
      ShoColinePair showerProduct2ShoColinePair = ShoColinePair();
      splitF->colourConnection(parentShoColinePair,
			       showerProduct1ShoColinePair, 
			       showerProduct2ShoColinePair);
      if ( showerProduct1ShoColinePair.first ) {
	showerProduct1ShoColinePair.first->addColoured( showerProduct1 );
      }
      if ( showerProduct1ShoColinePair.second ) {
	showerProduct1ShoColinePair.second->
	  addAntiColoured( showerProduct1 );
      }
      if ( showerProduct2ShoColinePair.first ) {
	showerProduct2ShoColinePair.first->addColoured( showerProduct2 );
      }
      if ( showerProduct2ShoColinePair.second ) {
	showerProduct2ShoColinePair.second->
	  addAntiColoured( showerProduct2 );
      }
      
      part->addChild(showerProduct1);
      part->addChild(showerProduct2);
      collecShoPar.push_back(showerProduct1);
      collecShoPar.push_back(showerProduct2);
      particlesYetToShower.push_back(showerProduct1);
      particlesYetToShower.push_back(showerProduct2);
    } 
  while (!particlesYetToShower.empty());
  return hasEmitted;  
}



bool ForwardEvolver::MEVeto(tcPPtr p, const Energy &q, const double &z) {
  bool veto = true; 
  double weight = 0.;
  QQbarG qqg(2*p->momentum().e(), p->momentum().m());
  if (abs(p->id()) < 7) { 
    if (p->id()>0) weight = qqg.qWeightX(q, z);
    else weight = qqg.qbarWeightX(q, z);
    veto = !rndbool(weight);
  }
  return veto; 
}


