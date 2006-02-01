// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HwDecayHandler class.
//

#include "HwDecayHandler.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "ThePEG/Handlers/Hint.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/PDT/Decayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/Step.h"
#include "ThePEG/EventRecord/Collision.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/Timer.h"
#include "ThePEG/Helicity/SpinInfo.h"
#include "Herwig++/Config/Herwig.h"
#include "Herwig++/Utilities/HwDebug.h"
#include <ThePEG/Repository/EventGenerator.h>
#include "DecayIntegrator.h"
#include "DecayPhaseSpaceMode.h"



namespace Herwig {
using namespace ThePEG;
using ThePEG::Helicity::tcSpinfoPtr;
using ThePEG::Helicity::SpinfoPtr;

HwDecayHandler::~HwDecayHandler() {}

void HwDecayHandler::
handle(EventHandler & ch, const tPVector & tagged,
       const Hint & hint) ThePEG_THROW_SPEC((Veto, Stop, Exception)) {
  // First go through the tagged particles for unstable ones
  Timer<46> timer("HwDecayHandler::handle");
  tPVector parents;
  for(int i = 0, N = tagged.size(); i<N; ++i)
    {
      if(tagged[i])
	{
	  // add to tagged if not stable
	  if(!tagged[i]->data().stable())
	    {parents.push_back(tagged[i]);}
	  // if stable and has spinInfo set the developed flag
	  else if(tagged[i]->spinInfo())
	    {
	      tcSpinfoPtr 
		hwspin=dynamic_ptr_cast<tcSpinfoPtr>(tagged[i]->spinInfo());
	      if(hwspin){hwspin->setDeveloped(true);}
	    }
	}
    }
  if(parents.empty()) return;
  // Create a new step, decay all particles and add their children to the step
  tStepPtr newStep = ch.newStep();
  for(int i = 0, N = parents.size(); i<N; ++i)
    {performDecay(newStep->find(parents[i]), *newStep);}
}


// perform decay method including modifications for spin correlations
// and for the decayer to specify intermediate decay products
void HwDecayHandler::performDecay(tPPtr parent, Step & s) const
  ThePEG_THROW_SPEC((Veto, Exception)) {
  Timer<47> timer("HwDecayHandler::performDecay");
  long ntry = 0;
  tcSpinfoPtr hwspin;
  while ( 1 ) 
    {
      // exit if fails
      if ( ++ntry >= maxLoop() ){throw DecHdlDecayFailed(parent->data(), maxLoop());}
      tDMPtr dm(parent->data().selectMode(*parent));
      if ( !dm ) throw DecHdlNoDecayMode(parent->data());
      if ( !dm->decayer() ) throw DecHdlNoDecayer(parent->data(), *dm);
      try 
	{
	  unsigned int hadronizetries(0);
	  bool hadronized;
	  // start of rejection loop to regenerate same decay mode with
	  // different kinematics this is sometimes needed for the hadronization
	  do
	    {
	      ParticleVector children = dm->decayer()->decay(*dm, *parent);
	      hadronized=true;
	      // decay generates decay products
	      if ( !children.empty() ) 
		{
		  // generate radiation in the decay
		  tDecayIntegratorPtr 
		    hwdec=dynamic_ptr_cast<tDecayIntegratorPtr>(dm->decayer());
		  if(hwdec)
		    {if(hwdec->canGeneratePhotons())
			{children=hwdec->generatePhotons(*parent,children);}}
		  // set up parent
		  parent->decayMode(dm);
		  // add children
		  for ( int i = 0, N = children.size(); i < N; ++i )
		    {
		      children[i]->setLabVertex(parent->labDecayVertex());
		      if ( !s.addDecayProduct(parent, children[i]) )
			throw DecHdlChildFail(parent->data(), children[i]->data());
		    }
		  parent->scale(0.0*GeV2);
		  // loop over the children
		  for ( int i = 0, N = children.size(); i < N; ++i )
		    {
		      // if the child has already been decayed add products to the record
		      if(children[i]->decayed()){addDecayedParticle(children[i],s);}
		      // if not stable decay the child
		      else if (!children[i]->data().stable())
			{performDecay(children[i], s);}
		      // if stable and has spinInfo set up decay matrices etc.
		      else if(children[i]->spinInfo())
			{
			  hwspin=dynamic_ptr_cast<tcSpinfoPtr>(children[i]->spinInfo());
			  if(hwspin){hwspin->setDeveloped(true);}
			}
		    }
		  // sort out the spinInfo for the parent after the decays
		  if(parent->spinInfo())
		    {
		      hwspin=dynamic_ptr_cast<tcSpinfoPtr>(parent->spinInfo());
		      // if the parent has the right kind of spinInfo
		      if(hwspin)
			{
			  // if the parent has been given a decay vertex
			  // calculate the decay matrix for the decay
			  if(hwspin->getDecayVertex()){hwspin->develop();}
			  // if the particle was scalar then it doesn't matter that it
			  // doesn't have a decay vertex as there's no correlations
			  else if(hwspin->iSpin()==PDT::Spin0){hwspin->setDeveloped(true);}
			}
		    }
		  // special for partonic decays
		  // see if colour particles produced
		  bool partonic(false);
		  for(unsigned int ix=0;ix<children.size();++ix)
		    {if(children[ix]->coloured()){partonic=true;break;}}
		  // if coloured particles were produced hadronize them
		  if(partonic&&_partonhad)
		    {
		      vector<tPPtr> hadrons;
		      hadronized=_partonhad->hadronize(parent,StepPtr(&s),
						       *(generator()->currentEventHandler()),
						       hadrons);
		      // if hadronization fails delete decay products and try again
		      if(!hadronized)
			{for(unsigned int ix=0;ix<children.size();++ix)
			    {s.removeParticle(children[ix]);}}
		      // otherwise decay the produced hadrons
		      else
			{
			  for(unsigned int ix=0;ix<hadrons.size();++ix)
			    {if(!hadrons[ix]->data().stable())
				{performDecay(hadrons[ix],s);}}
			}
		    }
		}
	      ++hadronizetries;
	    }
	  while(hadronizetries<_hadtry&&!hadronized);
	  // veto the decay and generate a new one 
	  if(hadronizetries>=_hadtry&&!hadronized)
	    {
	      string mode=dm->parent()->PDGName() + " -> ";
	      for(unsigned int ix=0;ix<dm->orderedProducts().size();++ix)
		{mode+=dm->orderedProducts()[ix]->PDGName() + " ";}
	      generator()->log() << "Failed to hadronize a partonic decay " << mode << "in"
				 << " HwDecayHandler::performDecay()" << endl;
	      throw Veto();
	    }
	  return;
	}
      catch (DecHdlChildFail) 
	{throw;}
      catch (Veto) 
	{}
    }
}

// method to add an intermediate which has already been decayed to the event record
void HwDecayHandler::addDecayedParticle(tPPtr parent, Step & s) const
  ThePEG_THROW_SPEC((Veto, Exception)) 
{
  try {
    for ( int i = 0, N = parent->children().size(); i < N; ++i ) {
      parent->children()[i]->setLabVertex(parent->labDecayVertex());
      s.addDecayProduct(parent->children()[i]);
    }
    parent->scale(0.0*GeV2);
    for ( int i = 0, N = parent->children().size(); i < N; ++i )
      {
	if((parent->children()[i])->decayed())
	  {
	    for(unsigned int ix=0;ix<(parent->children()[i])->children().size();++ix)
	      {addDecayedParticle(parent->children()[i],s);}
	  }
	else if ( !(parent->children()[i])->data().stable() )
	  {performDecay((parent->children()[i]), s);}
	else if(parent->children()[i]->data().stable())
	  {
	    if((parent->children()[i])->spinInfo())
	      {
		tcSpinfoPtr hwspin=
		  ThePEG::dynamic_ptr_cast<tcSpinfoPtr>(parent->children()[i]->spinInfo());
		if(hwspin){hwspin->setDeveloped(true);}
	      }
	  }
      }
    return;
  }
  catch (DecHdlChildFail) {
    throw;
  }
}

void HwDecayHandler::persistentOutput(PersistentOStream & os) const 
{os << _partonhad << _hadtry;}

void HwDecayHandler::persistentInput(PersistentIStream & is, int) 
{is >> _partonhad >> _hadtry;}

ClassDescription<HwDecayHandler> HwDecayHandler::initHwDecayHandler;

void HwDecayHandler::Init() 
{

  static ClassDocumentation<HwDecayHandler> documentation
    ("This is the handler for decays in Herwig++.");

  static Reference<HwDecayHandler,PartonicHadronizer> interfacePartonicHadronizer
    ("PartonicHadronizer",
     "Pointer to the object which hadronizes partonic decays.",
     &HwDecayHandler::_partonhad, false, false, true, true, false);


  static Parameter<HwDecayHandler,unsigned int> interfaceHadronizationTries
    ("HadronizationTries",
     "Number of times to regenerate a partonic decay to try and sucessfully"
     " hadronize it.",
     &HwDecayHandler::_hadtry, 20, 1, 50,
     false, false, Interface::limited);

}

}
