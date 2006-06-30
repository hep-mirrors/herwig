// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HeavyDecayer class.
//

#include "HeavyDecayer.h"
#include <ThePEG/PDT/EnumParticles.h>
#include <ThePEG/PDT/DecayMode.h>
#include <ThePEG/Interface/ClassDocumentation.h>
#include <ThePEG/Interface/Parameter.h>
#include "Herwig++/Utilities/Kinematics.h"
#include <ThePEG/Repository/EventGenerator.h>
#include <ThePEG/Repository/Repository.h>
#include <ThePEG/Utilities/Rebinder.h>
#include <ThePEG/Interface/Reference.h>
#include <ThePEG/Interface/Parameter.h>
#include <ThePEG/Persistency/PersistentOStream.h>
#include <ThePEG/Persistency/PersistentIStream.h>
#include <ThePEG/Utilities/Math.h>
#include "Herwig++/Utilities/HwDebug.h"
#include <ThePEG/Handlers/HandlerGroup.h>
#include <ThePEG/Handlers/Hint.h>
#include <ThePEG/EventRecord/Event.h>
#include <ThePEG/Handlers/EventHandler.h>

using namespace Herwig;
using namespace ThePEG;

void HeavyDecayer::Init() {
   static ClassDocumentation<HeavyDecayer> documentation
     ("Class to decay all particles in HERWIG by the algorithms used in HERWIG 6.4");

   static Parameter<HeavyDecayer, int>
      interfaceME("MECode", "The code for the ME type to use in the decay",
                  &Herwig::HeavyDecayer::MECode, 0, 0, 0, true, false, false);
}

ClassDescription<HeavyDecayer> HeavyDecayer::initHeavyDecayer;
long HeavyDecayer::lastAddedNumber = -1;

bool HeavyDecayer::accept(const DecayMode &dm) const { 
  long id = dm.parent()->id();
  int flav1, flav2;
  if((id / 1000)%10) {
    flav1 = (id/1000)%10;
    flav2 = (id/10)%100;
  } else {
    flav1 = id/100;  
    flav2 = (id/10)%10;
  }
  if(!flav1 || !flav2) return false;
  if(dm.products().size() != 4) {
    generator()->log() << "Heavy Decayer can only decay to 4 partons, given\n"
	               << "a mode with too many or too few products." 
		       << dm.tag() << endl;
    return false;
  } else return true;
}

/****** 
 * This function actually decays a particle based on the dm given.
 * The basic idea is that the heavy parton in the heavy meson will decay
 * weakly while the other parton will be more or less unchanged. This produces
 * 4 partons, the spectator plus the result of the weak decay. The W then
 * decays again into two more partons.
 * e.g. a decay of B0 --> d, cbar, dbar, u
 *
 *    bbar  ----> cbar (colour connected to spectator)
 *                W+--> dbar (colour connected to other W+ product)
 *                      u 
 *    d     ----> spectator --> d
 * The resulting partons are then hadronized and decayed again.
 *****/
ParticleVector HeavyDecayer::decay(const DecayMode &dm, const Particle &p) const
{
   if( HERWIG_DEBUG_LEVEL >= HwDebug::full) {
     generator()->log() << "HeavyDecayer::decay called on " << p.PDGName() 
			<< "\n";
   }
   ParticleVector partons = dm.produceProducts();

   if( HERWIG_DEBUG_LEVEL >= HwDebug::full) {
     generator()->log() << "HeavyDecayer::Decaying " << p.PDGName() << " via " 
	                << dm.tag() << "...\n";
   }

   Lorentz5Momentum products[4];
   for(int i = 0; i<4; i++) 
      products[i].setMass(partons[i]->mass());

   // Fraction of momentum to spectator
   double xs = partons[3]->mass()/p.momentum().e();
   // Fraction of momentum to heavy quark (which is weakly decaying)
   double xb = 1.-xs;
   partons[3]->setMomentum(p.momentum()*xs);

   // Get the particle that is decaying
   long idQ, idSpec;
   idSpec = partons[3]->id();
   idQ = (p.id()/1000)%10;
   if(!idQ) idQ = (p.id()/100)%10;

   // Now the odd case of a B_c where the c decays, not the b
   if(idSpec == idQ) idQ = (p.id()/10)%10;

   PPtr inter = getParticleData(idQ)->produceParticle(p.momentum()*xb);

   // Three Body Decay
   // Free Massless (V-A)*(V-A) ME
   if(MECode == 100) {
     double EMwSq, GMwSq, EmLim, EmTest, GamW;
     Energy Mw = getParticleData(ParticleID::Wplus)->mass();
     GamW = getParticleData(ParticleID::Wplus)->width();
     EMwSq = sqr(Mw);
     GMwSq = sqr(Mw*GamW);
     EmLim = GMwSq + sqr(EMwSq - sqr(inter->mass()-p.mass()));
     do {
       Kinematics::threeBodyDecay(inter->momentum(),products[0], products[1], 
				  products[2], &VAWt);
       Lorentz5Momentum pw = products[0]+products[1];
       double pw2 = pw*pw;
       EmTest = sqr(pw2-EMwSq);
     } while((EmTest+GMwSq)*rnd() > EmLim);
   }
   else Kinematics::threeBodyDecay(inter->momentum(),products[0], products[1],
				   products[2]);

   for(int i = 0; i<3; i++) partons[i]->setMomentum(products[i]);

   //cout << "Doing colour connections for " << dm.tag() << endl;
   // Set up colour connections based on the diagram above and input order
   if(partons[0]->coloured()) {
      if(partons[0]->id() > 0)
         partons[0]->antiColourNeighbour(partons[1]);
      else partons[0]->colourNeighbour(partons[1]);
   }

   if(partons[2]->id() > 0) 
      partons[2]->antiColourNeighbour(partons[3]);
   else
      partons[2]->colourNeighbour(partons[3]);

   if( HERWIG_DEBUG_LEVEL >= HwDebug::full) {
     generator()->log() << "HeavyDecayer::Decaying " << "...Done\n";
   }


   // Insure we only add new handler once per event
   /*long evtNum = generator()->currentEventNumber();
   if(evtNum != lastAddedNumber) {
     lastAddedNumber = evtNum; 
     tcCollHdlPtr ch = dynamic_ptr_cast<tcCollHdlPtr>(
		                   generator()->currentEvent()->handler());
     tCollHdlPtr cp = const_ptr_cast<tCollHdlPtr>(ch);
     cp->addStep(ThePEG::Group::main, ThePEG::Group::hadron, 
  	         StepHdlPtr(), Hint::Default());
   }
*/
   //for(ParticleVector::iterator it=partons.begin(); it!=partons.end(); it++) 
   //  cout << *(*it) << endl;

   //cout << "Returning partons" << endl;
   return partons;
}
   
void HeavyDecayer::persistentOutput(PersistentOStream &os) const { os << MECode; }
void HeavyDecayer::persistentInput(PersistentIStream &is, int i) { is >> MECode; }

double HeavyDecayer::VAWt(double *temp) 
{ return (temp[1]-temp[0])*(temp[0]-temp[2])*temp[3]; }
