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
#include "Herwig++/Hadronization/ClusterHadronizationHandler.h"
#include <fstream>

// using namespace Pythia7;
using namespace Herwig;


void initFirstStep(StepPtr& s, const EGPtr & eg) {
  
//   //------------------ TEST GLUON SPLITTING --------------------
//  
//   PPtr u    = eg->getParticle(ParticleID::u);
//   PPtr ubar = eg->getParticle(ParticleID::ubar);
//   PPtr g1   = eg->getParticle(ParticleID::g);
//    
//   Momentum3 p_u(0.0*GeV, 0.0*GeV, 10.0*GeV);
//   Momentum3 p_ubar(0.0*GeV, 0.0*GeV, -10.0*GeV);
//   Momentum3 p_g1(10.0*GeV, 0.0*GeV, 0.0*GeV);
//  
//   u->set3Momentum(p_u);
//   ubar->set3Momentum(p_ubar);
//   g1->set3Momentum(p_g1);
//  
//   // Let's suppose that these partons have perturbative origin,
//   // for example at scale (100*GeV)^2
//   u->scale( sqr(100.0*GeV) );
//   ubar->scale( sqr(100.0*GeV) );
//   g1->scale( sqr(100.0*GeV) );
//   
//   //***LOOKHERE*** To avoid on-shell gluon
//   Lorentz5Momentum gluon5Momentum = Lorentz5Momentum( 0.750*GeV, p_g1);
//   g1->set5Momentum( gluon5Momentum );
// 
//   u->antiColourNeighbour(g1);
//   g1->colourNeighbour(u);
//   g1->antiColourNeighbour(ubar);
//   ubar->colourNeighbour(g1);
//   
//   // Equivalent to what is done above   
//   // tPPtr t_g1 = g1;
//   // tPPtr t_ubar = ubar;
//   // ColinePtr coline1 = ColourLine::create(u,t_g1);
//   // ColinePtr coline2 = ColourLine::create(g1,t_ubar); 
// 
//   s->addParticle(u);
//   s->addParticle(ubar);
//   s->addParticle(g1);
 
//   //------------------ TEST COLOUR CONNECTIONS --------------------
//  
//   PPtr c    = eg->getParticle(ParticleID::c);
//   PPtr cbar = eg->getParticle(ParticleID::cbar);
//   PPtr b    = eg->getParticle(ParticleID::b);
//   PPtr bbar = eg->getParticle(ParticleID::bbar);
//   PPtr mu   = eg->getParticle(ParticleID::muminus);
//   
//   Momentum3 p_c(10.0*GeV, 0.0*GeV, 0.0*GeV);
//   Momentum3 p_cbar(-10.0*GeV, 0.0*GeV, 0.0*GeV);
//   Momentum3 p_b(0.0*GeV, 10.0*GeV, 0.0*GeV);
//   Momentum3 p_bbar(0.0*GeV, -10.0*GeV, 0.0*GeV);
//   Momentum3 p_mu(10.0*GeV, 10.0*GeV, 10.0*GeV);
//   
//   c->set3Momentum(p_c);
//   cbar->set3Momentum(p_cbar);
//   b->set3Momentum(p_b);
//   bbar->set3Momentum(p_bbar);
//   mu->set3Momentum(p_mu);
//   
//   c->scale( sqr(100.0*GeV) );
//   cbar->scale( sqr(100.0*GeV) );
//   b->scale( sqr(100.0*GeV) );
//   bbar->scale( sqr(100.0*GeV) );
//   mu->scale( sqr(100.0*GeV) );
//  
//   c->antiColourNeighbour(bbar);
//   cbar->colourNeighbour(b);
//   b->antiColourNeighbour(cbar);
//   bbar->colourNeighbour(c);
// 
//   // Equivalent to what is done above   
//   // tPPtr t_bbar = bbar;
//   // tPPtr t_cbar = cbar;  
//   // ColinePtr coline1 = ColourLine::create(c,t_bbar);
//   // ColinePtr coline2 = ColourLine::create(b,t_cbar);
// 
//   s->addParticle(c);
//   s->addParticle(cbar);
//   s->addParticle(b);
//   s->addParticle(bbar);
//   s->addParticle(mu);
 
//   //------------------ TEST 3-COMPONENT CLUSTERS AND 
//   //                   (AT THE SAME TIME) HEAVY DIQUARK   --------------------
// 
//   PPtr strange    = eg->getParticle(ParticleID::s);
//   PPtr strangebar = eg->getParticle(ParticleID::sbar);
//   PPtr c          = eg->getParticle(ParticleID::c);
//   PPtr cbar       = eg->getParticle(ParticleID::cbar);
//   PPtr b          = eg->getParticle(ParticleID::b);
//   PPtr bbar       = eg->getParticle(ParticleID::bbar);
//   
//   Momentum3 p_strange(5.0*GeV, 0.0*GeV, 5.0*GeV);
//   Momentum3 p_strangebar(-5.0*GeV, 0.0*GeV, -5.0*GeV);
//   Momentum3 p_c(10.0*GeV, 0.0*GeV, 0.0*GeV);
//   Momentum3 p_cbar(-10.0*GeV, 0.0*GeV, 0.0*GeV);
//   Momentum3 p_b(0.0*GeV, 10.0*GeV, 0.0*GeV);
//   Momentum3 p_bbar(0.0*GeV, -10.0*GeV, 0.0*GeV);
//   
//   strange->set3Momentum(p_strange);
//   strangebar->set3Momentum(p_strangebar);
//   c->set3Momentum(p_c);
//   cbar->set3Momentum(p_cbar);
//   b->set3Momentum(p_b);
//   bbar->set3Momentum(p_bbar);
//   
//   strange->scale( sqr(100.0*GeV) );
//   strangebar->scale( sqr(100.0*GeV) );
//   c->scale( sqr(100.0*GeV) );
//   cbar->scale( sqr(100.0*GeV) );
//   b->scale( sqr(100.0*GeV) );
//   bbar->scale( sqr(100.0*GeV) );
//   
//   ColinePtr l_strange = ColourLine::create(strange);
//   ColinePtr l_c = ColourLine::create(c);
//   ColinePtr l_b = ColourLine::create(b);
//   l_b->setSourceNeighbours(l_strange,l_c);
//   
//   bool anti = true;
//   ColinePtr l_strangebar = ColourLine::create(strangebar,anti);
//   ColinePtr l_cbar = ColourLine::create(cbar,anti);
//   ColinePtr l_bbar = ColourLine::create(bbar,anti);
//   l_bbar->setSinkNeighbours(l_strangebar,l_cbar);
//   
//   s->addParticle(strange);
//   s->addParticle(strangebar);
//   s->addParticle(c);
//   s->addParticle(cbar);
//   s->addParticle(b);
//   s->addParticle(bbar);

  //------------------ TEST SPECIAL 2-COMPONENT CLUSTERS FROM 
  //                   TWO BARYON VIOLATING DECAYS             --------------------

  PPtr tilda_u_r = eg->getParticle(ParticleID::SUSY_u_R);
  PPtr tilda_u_r_star = eg->getParticle(ParticleID::SUSY_u_Rbar);
  tilda_u_r->antiColourNeighbour(tilda_u_r_star);
  tilda_u_r_star->colourNeighbour(tilda_u_r);
  s->addParticle(tilda_u_r);
  s->addParticle(tilda_u_r_star);

  PPtr d    = eg->getParticle(ParticleID::d);
  PPtr dbar = eg->getParticle(ParticleID::dbar);
  PPtr b    = eg->getParticle(ParticleID::b);
  PPtr bbar = eg->getParticle(ParticleID::bbar);  
  Momentum3 p_d(5.0*GeV, 0.0*GeV, 5.0*GeV);
  Momentum3 p_b(-5.0*GeV, 0.0*GeV, -5.0*GeV);
  Momentum3 p_dbar(10.0*GeV, 0.0*GeV, 0.0*GeV);
  Momentum3 p_bbar(-10.0*GeV, 0.0*GeV, 0.0*GeV);  
  d->set3Momentum(p_d);
  b->set3Momentum(p_b);
  dbar->set3Momentum(p_dbar);
  bbar->set3Momentum(p_bbar);  
  d->scale( sqr(100.0*GeV) );
  b->scale( sqr(100.0*GeV) );
  dbar->scale( sqr(100.0*GeV) );
  bbar->scale( sqr(100.0*GeV) );
  
  ColinePtr line_d = ColourLine::create(d);
  ColinePtr line_b = ColourLine::create(b);
  tilda_u_r_star->antiColourLine()->setSourceNeighbours(line_d,line_b);
  
  bool anti = true;
  ColinePtr line_dbar = ColourLine::create(dbar,anti);
  ColinePtr line_bbar = ColourLine::create(bbar,anti);
  tilda_u_r->colourLine()->setSinkNeighbours(line_dbar,line_bbar);
  
  s->addDecayProduct(tilda_u_r,dbar);
  s->addDecayProduct(tilda_u_r,bbar);
  s->addDecayProduct(tilda_u_r_star,d);
  s->addDecayProduct(tilda_u_r_star,b);
  s->fixColourFlow();       	
	
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

  return(0);
}

