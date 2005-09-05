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



using namespace ThePEG;
using namespace Herwig;
using Helicity::tcSpinfoPtr;
using Helicity::SpinfoPtr;

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
  while ( 1 ) 
    {
      if ( ++ntry >= maxLoop() ){throw DecHdlDecayFailed(parent->data(), maxLoop());}
      tDMPtr dm(parent->data().selectMode(*parent));
      if ( !dm ) throw DecHdlNoDecayMode(parent->data());
      if ( !dm->decayer() ) throw DecHdlNoDecayer(parent->data(), *dm);
      try {
	unsigned int hadronizetries(0);
	bool hadronized;
	do
	  {
	    ParticleVector children = dm->decayer()->decay(*dm, *parent);
	    hadronized=true;
	    if ( !children.empty() ) 
	      {
		parent->decayMode(dm);
		for ( int i = 0, N = children.size(); i < N; ++i )
		  {
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
		    // if stable and has spinInfo
		    else if(children[i]->spinInfo())
		      {
			tcSpinfoPtr 
			  hwspin=dynamic_ptr_cast<tcSpinfoPtr>(children[i]->spinInfo());
			if(hwspin){hwspin->setDeveloped(true);}
		      }
		  }
		// sort out the spinInfo for the parent after the decays
		if(parent->spinInfo())
		  {
		    tcSpinfoPtr hwspin=dynamic_ptr_cast<tcSpinfoPtr>(parent->spinInfo());
		    // if the parent has the right kind of spinInfo
		    if(hwspin)
		      {
			// if the parent has been given a decay vertex
			if(hwspin->getDecayVertex())
			  {
			    // calculate the decay matrix for the decay
			    hwspin->develop();
			    // if debugging check whether all the products were handled
			    if(HERWIG_DEBUG_LEVEL >=HwDebug::full)
			      {
				tcSpinfoPtr hwtemp;
				bool alldeveloped=true; 
				for(unsigned int ix=0,
				      N=hwspin->getDecayVertex()->outgoing().size();
				    ix<N;++ix)
				  {
				    hwtemp = dynamic_ptr_cast<tcSpinfoPtr>
				      (hwspin->getDecayVertex()->outgoing()[ix]);
				    if(hwtemp)
				      {if(!(hwtemp->developed())){alldeveloped=false;}}
				  }
				if(!alldeveloped)
				  {(generator()->log()) << "HwDecayHandler::performDecay " 
							<< "all the decay products of " 
							<< parent->id() 
							<< " where not developed" << endl;}
			      }
			  }
			// if the particle was scalar then it does matter that it does have a 
			// decay vertex as there's no correlations
			else if(hwspin->iSpin()==PDT::Spin0)
			  {hwspin->setDeveloped(true);}
			else if(HERWIG_DEBUG_LEVEL >=HwDebug::full)
			  {generator()->log() << "HwDecayHandler::performDecay " 
					      << parent->id() 
					      << "was not decayed by a decay with " 
					      << "correlations" << endl;}
		      }
		  }
		// special for partonic decays
		bool partonic(false);
		for(unsigned int ix=0;ix<children.size();++ix)
		  {if(children[ix]->coloured()){partonic=true;break;}}
		if(partonic)
		  {
		    vector<tPPtr> hadrons;
		    hadronized=
		      _partonhad->hadronize(parent,StepPtr(&s),
					    *(generator()->currentEventHandler()),
					    hadrons);
		    // if hadronization fails delete decay products and try again
		    if(!hadronized)
		      {
			for(unsigned int ix=0;ix<children.size();++ix)
			  {s.removeParticle(children[ix]);}
		      }
		    // decay the produced hadrons
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
	if(hadronizetries>=_hadtry)
	  {
	    string mode=dm->parent()->PDGName() + " -> ";
	    for(unsigned int ix=0;ix<dm->orderedProducts().size();++ix)
	      {mode+=dm->orderedProducts()[ix]->PDGName() + " ";}
	    throw Exception() << "Failed to hadronize a partonic decay "
			      << mode << "in"
			      << " HwDecayHandler::performDecay()" 
			      << Exception::eventerror;
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
    for ( int i = 0, N = parent->children().size(); i < N; ++i )
	{s.addDecayProduct(parent->children()[i]);}
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
		Helicity::tcSpinfoPtr hwspin=ThePEG::dynamic_ptr_cast<tcSpinfoPtr>(parent->children()[i]->spinInfo());
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
     &HwDecayHandler::_partonhad, false, false, true, false, false);

  static Parameter<HwDecayHandler,unsigned int> interfaceHadronizationTries
    ("HadronizationTries",
     "Number of times to regenerate a partonic decay to try and sucessfully"
     " hadronize it.",
     &HwDecayHandler::_hadtry, 10, 1, 50,
     false, false, Interface::limited);

}

