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
#include "Herwig++/Shower/ShowerHandler.h"
#include <fstream>

// using namespace Pythia7;
using namespace Herwig;


void initFirstStep(StepPtr& s, const EGPtr & eg){

  //*** FILL BY HAND THE EVENT RECORD *** 
  // (Below a very simple example)

  PPtr u    = eg->getParticle(ParticleID::u);
  PPtr ubar = eg->getParticle(ParticleID::ubar);
  
  Momentum3 p_u(0.0*GeV, 0.0*GeV, 10.0*GeV);
  Momentum3 p_ubar(0.0*GeV, 0.0*GeV, -10.0*GeV);
  
  u->set3Momentum(p_u);
  ubar->set3Momentum(p_ubar);
  
  // Let's suppose that these partons have perturbative origin,
  // for example at scale (100*GeV)^2
  u->scale( sqr(100.0*GeV) );
  ubar->scale( sqr(100.0*GeV) );
 
  u->antiColourNeighbour(ubar);
  ubar->colourNeighbour(u);
 
  s->addParticle(u);
  s->addParticle(ubar);
   
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

  vector<tPPtr> products;  

  breakPythia7();

  long ntry = eg->N();
  if (N > 0) ntry = N;

  for(long itry=0; itry!=ntry; ++itry){

    EventPtr event = new_ptr(Event(PPair()));

    StepPtr firstStep = event->newStep();   // Create an empty step.
     
    initFirstStep(firstStep, eg);

    if(products.size() ) products.clear();

    event = eg->partialEvent( event );

    event->selectFinalState(back_inserter(products) );

    if( itry < 2 ) cout << *event;
    
    int npart = products.size();
    cout << " event = " << itry
         << "\t products.size()= " << npart << endl;

  }
  eg->finish();
  
  cout << "program finished." << endl; 

  return(0);
}

