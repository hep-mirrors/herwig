// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Hw64Decayer class.
//

#include "Hw64Decayer.h"
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

using namespace Herwig;
using namespace ThePEG;

Hw64Decayer::~Hw64Decayer() {}

void Hw64Decayer::Init() {
   static ClassDocumentation<Hw64Decayer> documentation
     ("Class to decay all particles in HERWIG by the algorithms used in HERWIG 6.4");

   static Parameter<Hw64Decayer, int>
      interfaceME("MECode", "The code for the ME type to use in the decay",
                  &Herwig::Hw64Decayer::MECode, 0, 0, 0, true, false, false);
}

ClassDescription<Hw64Decayer> Hw64Decayer::initHw64Decayer;

bool Hw64Decayer::accept(const DecayMode &dm) const { return true; }

/****** 
 * This function actually decays a particle based on the dm given.
 *****/
ParticleVector Hw64Decayer::decay(const DecayMode &dm, const Particle &p) const
{
   if( HERWIG_DEBUG_LEVEL >= HwDebug::full) {
     generator()->log() << "Hw64Decay::decay called on " << p.PDGName() << "\n";
   }
   ParticleVector rval;
   ParticleMSet productParticles = dm.products();
   int numProds = productParticles.size();

   // Create a vector of up to 5 momentum
   vector<Lorentz5Momentum> products(5);
   vector<Energy> masses(5);

   if( HERWIG_DEBUG_LEVEL >= HwDebug::full) {
     generator()->log() << "Hw64Decay::Decaying " << p.PDGName() << " via " << dm.tag() << "...\n";
   }

   if(numProds > 5) { 
     generator()->log() << "Number of decay products is greater than 5\n"
                        << "Hw64Decayer cannot handle these decays properly. "
                        << "Will treat as 5 body decay and ignore excess "
                        << "particles.\n";
     numProds = 5;
   }

   int i = 0;
   for(ParticleMSet::iterator it = productParticles.begin(); i<numProds; i++, it++) 
   {
      masses[i] = (*it)->mass();
      products[i].setMass(masses[i]);
   }
   if(p.mass() < 0.000001) {
     generator()->log() << "HwDecayer called on a particle with no mass " 
			<< p.PDGName() << ", " << p.mass() << endl;
     int i = 0;
     for(ParticleMSet::iterator it = dm.products().begin(); it!=dm.products().end(); it++) 
           rval.push_back((*it)->produceParticle(products[i++]));	
     return rval;
   }

   // The K -> KL0 and KS0
   if(numProds == 1) {
      oneBodyDecay(p.momentum(), products[0]);

   // 2 Body Decay
   } else if(numProds == 2) {
      double CosAngle, AzmAngle;
      
      Kinematics::generateAngles(CosAngle, AzmAngle);

      /*******
       * It appears polarized mesons aren't needed currently 
      if(p->id() == ) {
         if() {
            CosAng = Math::absmax<double>(generator()->rnd(), 
		           Math::absmax<double>(generator()->rnd(),  generator()->rnd()))*2.0 - 1.0;
         } else if() {
            CosAng = 2.0*cos((acos(generator()->rnd()*2.0-1.0)+M_PI)/3.0);
         } else {
            CosAng = Math::absmin<double>(generator()->rnd(),
                           Math::absmin<double>(generator()->rnd(),  generator()->rnd()))*2.0 - 1.0;
         }
      }
       *******/
      Kinematics::twoBodyDecay(p.momentum(), masses[0], masses[1],
			       CosAngle, AzmAngle, products[0], products[1]);

      // Now we rotate the result (why?)
      // LorentzRotation R = Kinematics::rotate(p.momentum(), 1.0, 0.0);
      // products[0] *= R; products[1]*=R;
 
   // Three Body Decay
   } else if(numProds == 3) {
      // Free Massless (V-A)*(V-A) ME
      if(MECode == 100) {
         Kinematics::threeBodyDecay(p.momentum(), products[0], products[1], 
				    products[2], &VAWt);
         
      // Bound Massless (V-A)*(V-A) ME
      } else if(MECode == 101) {
         double wtmx, wtmx2, xs, dot1, dot2;
         wtmx = ( (p.mass() - products[2][5])*(p.mass() + products[2][5])
                + (products[1][5] - products[0][5])*
		  (products[1][5] + products[0][5]))/2.0;
         wtmx2 = sqr(wtmx);

         // Find sum of masses of constituent particles 
         int IPDG = abs(p.id());
         double m1, m2, m3;
         if(IPDG >= 1000)
	   m1 = generator()->getParticleData((IPDG/1000)%10)->mass();
	 else
	   m1 = 0.0;
         m2 = generator()->getParticleData((IPDG/100)%10)->mass();
         m3 = generator()->getParticleData((IPDG/10)%10)->mass();
         xs = 1.0 - Math::absmax<double>(m1, Math::absmax<double>(m2, m3))/(m1+m2+m3);

	 // Do decay, repeat until meets condition
         do {
            Kinematics::threeBodyDecay(p.momentum(), products[1], products[2], 
				       products[0], &VAWt);
            dot1 = p.momentum().dot(products[1]);
            dot2 = p.momentum().dot(products[0]);
         } while(dot1*(wtmx-dot1-xs*dot2) < generator()->rnd()*wtmx2);
      } 
 
      // Or use free massless ((V-A)*TB1+(V+A)*CT1)*((V-A)*TB2+(V+A)*CT2)) Matrix Element
      // No particles seem to use this code, we will ignore it currently
      /*else if(MECode == 200) {
         // Sort tan(beta)
         #define IDP (products[0].id())
         #define PID(a) (ParticleID::a)
         if(IDP == PID(u) || IDP == PID(c) || IDP == PID(t) || 
            IDP == PID(ubar) || IDP == PID(cbar) || IDP == PID(tbar) ||
            IDP == PID(nu_e) || IDP == PID(nu_mu) || IDP == PID(nu_tau) ||
            IDP == PID(nu_ebar) || IDP == PID(nu_mubar) || IDP == PID(nu_taubar))
       */

      else Kinematics::threeBodyDecay(p.momentum(), products[0], products[1], 
				      products[2]);

   // Four Body Decay
   } else if(numProds == 4) {
     Kinematics::fourBodyDecay(p.momentum(), products[0], products[1], 
			       products[2], products[3]);

   // Five Body Decay
   } else if(numProds == 5) {
     Kinematics::fiveBodyDecay(p.momentum(), products[0], products[1], 
			       products[2], products[3], products[4]);
   }

   if(products[0] == Lorentz5Momentum()) {
     generator()->warning() << "The Decay mode " << dm.tag() << " cannot "
	                    << "proceed, not enough phase space\n";   
   }
   setParticleMomentum(rval, productParticles, products);

   if( HERWIG_DEBUG_LEVEL >= HwDebug::full) {
     generator()->log() << "Hw64Decay::Decaying " << "...Done\n";
   }

   return rval;
}
   
void Hw64Decayer::persistentOutput(PersistentOStream &os) const { os << MECode; }
void Hw64Decayer::persistentInput(PersistentIStream &is, int i) { is >> MECode; }

/******
 * This function takes the array of momentum generated and sets the momentum
 * to the particles.
 *****/
void Hw64Decayer::setParticleMomentum(ParticleVector &out, ParticleMSet particles, 
                                      vector<Lorentz5Momentum> moms) const {
   PPtr child;
   int i = 0;
   int numProds = particles.size();
   // Can't handle higher decays...warning has already been given
   if(numProds > 5) numProds = 5;

   for(ParticleMSet::iterator it = particles.begin(); i<numProds; i++, it++) {
      child = (*it)->produceParticle(moms[i]);
      out.push_back(child);
   }
}

double Hw64Decayer::VAWt(double *temp) { 
  //double emsq, double a, double b, double c) const {
  // return (a-emsq)*(emsq-b)*c;
  return (temp[1]-temp[0])*(temp[0]-temp[2])*temp[3];
}

//double Hw64Decayer::PhaseSpaceWt() const {
//   return 1.0;
//}

void Hw64Decayer::oneBodyDecay(Lorentz5Momentum p0, Lorentz5Momentum &p1)
{
   p1 = (LorentzVector)p0;
}


