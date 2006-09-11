// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HeavyDecayer class.
//

#include "QuarkoniumDecayer.h"
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

void QuarkoniumDecayer::Init() {
   static ClassDocumentation<QuarkoniumDecayer> documentation
     ("Class to decay all particles in HERWIG by the algorithms used in HERWIG 6.4");

   static Parameter<QuarkoniumDecayer, int>
      interfaceME("MECode", "The code for the ME type to use in the decay",
                  &Herwig::QuarkoniumDecayer::MECode, 0, 0, 0, true, false, false);
   
}

ClassDescription<QuarkoniumDecayer> QuarkoniumDecayer::initQuarkoniumDecayer;

bool QuarkoniumDecayer::accept(const DecayMode &dm) const { 
  if(dm.products().size() == 3 || dm.products().size() == 2)
    return true;
  return false;
}

/****** 

 *****/
ParticleVector QuarkoniumDecayer::decay(const DecayMode &dm, const Particle &p) const
{
   if( HERWIG_DEBUG_LEVEL >= HwDebug::full) {
     generator()->log() << "Hw64Decay::decay called on " << p.PDGName() 
			<< "\n";
   }

   ParticleVector partons = dm.produceProducts();
   Lorentz5Momentum products[4]; // Shouldn't be any bigger!
   Energy gluMass = getParticleData(ParticleID::g)->constituentMass();
   for(unsigned int i = 0; i<partons.size(); i++) {
     if(partons[i]->id() == ParticleID::g) products[i].setMass(gluMass);
     else                                  products[i].setMass(partons[i]->mass());
   }

   //cout << "generating a quarkonium decay " << dm.tag() << endl;
   if(partons.size() == 3) {
     // 3-gluon or 2-gluon + photon decay
     if(MECode == 130) { // Ore & Powell orthopositronium matrix element
       double x1, x2, x3, test;
       do {
	 Kinematics::threeBodyDecay(p.momentum(), products[0], products[1],
				    products[2]);
	 x1 = 2.*(p.momentum()*products[0])/sqr(p.mass());
	 x2 = 2.*(products[1]*products[2])/sqr(p.mass());
	 x3 = 2. - x1 - x2;
	 test = sqr(x1*(1.-x1)) + sqr(x2*(1.-x2)) + sqr(x3*(1.-x3));
	 test /= sqr(x1*x2*x3);
       } while(test < 2.*UseRandom::rnd());
     } else Kinematics::threeBodyDecay(p.momentum(), products[0], products[1],
				       products[2]);
     
     for(unsigned int i = 0; i<partons.size(); i++)
       partons[i]->set5Momentum(products[i]);

     // Now set colour connections
     if(partons[2]->id() == ParticleID::g) {
       partons[0]->colourNeighbour(partons[1]);
       partons[1]->colourNeighbour(partons[2]);
       partons[2]->colourNeighbour(partons[0]);
       //cout << "Colour connecting : " << partons[0]->PDGName() << "-->"
//	    << partons[1]->PDGName() << "-->" << partons[2]->PDGName() 
//	    << "==>" << partons[0]->PDGName() << endl;
     } else {
       partons[0]->colourNeighbour(partons[1]);
       partons[0]->antiColourNeighbour(partons[1]);
 //      cout << "Colour connecting : "<< partons[0]->PDGName() << "<==>"
//	    << partons[1]->PDGName() << endl;
     }
   } else {
     // 2 gluon or q-qbar decay
     double Theta, Phi;
     Kinematics::generateAngles(Theta, Phi);
     Energy p1, p2;
     p1 = partons[0]->mass();
     p2 = partons[1]->mass();
     if(p1 == 0.0) p1 = gluMass;
     if(p2 == 0.0) p2 = gluMass;
     Kinematics::twoBodyDecay(p.momentum(), p1, p2, Theta, Phi,
			      products[0], products[1]);
     for(unsigned int i = 0; i<partons.size(); i++) {
       partons[i]->set5Momentum(products[i]);
       //cout << *partons[i] << endl;
       //cout << products[i] << endl;
     }

     int first, second; 
     if(partons[0]->id() > 0) { first = 0; second = 1; }
     else { first = 1; second = 0; }
     partons[first]->antiColourNeighbour(partons[second]);
     if(abs(partons[first]->id()) == ParticleID::g) 
       partons[first]->colourNeighbour(partons[second]);
     //cout << "Colour connecting : "<< partons[0]->PDGName() << "<-->"
	  //<< partons[1]->PDGName() << endl;
   }
   if( HERWIG_DEBUG_LEVEL >= HwDebug::full) {
     generator()->log() << "QuarkoniumDecayer::Decaying " << p.PDGName() 
			<< " via " << dm.tag() << "...\n";
   }

   if( HERWIG_DEBUG_LEVEL >= HwDebug::full) {
     generator()->log() << "QuarkoniumDecayer::Decaying " << "...Done\n";
   }

   //cout << "done decaying\n";
   return partons;
}
   
void QuarkoniumDecayer::persistentOutput(PersistentOStream &os) const { 
  os << MECode; 
}

void QuarkoniumDecayer::persistentInput(PersistentIStream &is, int i) { 
  is >> MECode; 
}


