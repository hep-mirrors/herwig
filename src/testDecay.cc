#include "Pythia7/Config/Pythia7.h"
#include "Pythia7/Repository/Pythia7Initializer.h"
#include "Pythia7/Repository/Repository.h"
#include "Pythia7/Repository/EventGenerator.h"
#include "Pythia7/Persistency/PersistentIStream.h"
#include "Pythia7/PDT/ParticleData.h"
#include "Pythia7/PDT/EnumParticles.h"
#include "Pythia7/EventRecord/Event.h"
#include "Pythia7/Utilities/Debug.h"
#include "Herwig++/Utilities/HwDebug.h"
#include "Herwig++/Decay/ExampleDecayer.h"
#include <fstream>

// using namespace Pythia7;
using namespace Herwig;


void initFirstStep(StepPtr& s, const EGPtr & eg){

  //*** FILL BY HAND THE EVENT RECORD *** 
  // (Below a very simple example)

  PPtr Z0    = eg->getParticle(ParticleID::Z0);
  Momentum3 p_Z0(10.0*GeV, 10.0*GeV, 10.0*GeV);
  Z0->set3Momentum(p_Z0);
  s->addParticle(Z0);
   
}



int main(int argc, char* argv[]){
  
  string run;
  long N = -1;

  for ( int iarg = 1; iarg < argc; ++iarg ) {
    string arg = argv[iarg];
    if ( arg == "-d" ) Debug::level = atoi(argv[++iarg]);
    else if ( arg.substr(0,2) == "-d" && arg.substr(0,4) != "-dHw" )
	Debug::level = atoi(arg.substr(2).c_str());
    else if ( arg == "-dHw" ) Herwig::HwDebug::level = atoi(argv[++iarg]);
    else if ( arg.substr(0,4) == "-dHw" )
	Herwig::HwDebug::level = atoi(arg.substr(4).c_str());
    else if ( arg == "-N" ) N = atoi(argv[++iarg]);
    else if ( arg.substr(0,2) == "-N" ) N = atoi(arg.substr(2).c_str());
    else if ( arg == "-h" ) {
    cerr << "Usage: " << argv[0]
	 << " [-d debug-level] [-dHw herwig-debug-level] run-file"
	 << endl;
      return 3;
    }
    else
      run = arg;
  }

  if ( run.empty() ) {
    cerr << "No run-file specified. Usage: \n " << argv[0]
	 << " [-d debug-level] [-dHw herwig-debug-level] run-file"
	 << endl;
    return 1;
  }

  PersistentIStream is(run);
  EGPtr eg;
  is >> eg;
  eg->initialize();

  //ALB vector<tPPtr> products;  

  breakPythia7();

  long ntry = eg->N();
  if (N > 0) ntry = N;

  for(long itry=0; itry!=ntry; ++itry){

    EventPtr event = new_ptr(Event(PPair()));

    StepPtr firstStep = event->newStep();   // Create an empty step.
     
    //ALB initFirstStep(firstStep, eg);

    PPtr Z0    = eg->getParticle(ParticleID::Z0);
    Momentum3 p_Z0(10.0*GeV, 10.0*GeV, 10.0*GeV);
    Z0->set3Momentum(p_Z0);
    firstStep->addParticle(Z0);

    //ALB if(products.size() ) products.clear();

    event = eg->partialEvent( event );

    tPVector products = event->getFinalState();
    //ALB event->selectFinalState(back_inserter(products) );

    if( itry < 2 ) cout << *event;
    
    int npart = products.size();
    cout << " event = " << itry
         << "\t products.size()= " << npart << endl;

  }
  eg->finish();

  return(0);
}

