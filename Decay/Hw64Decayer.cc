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
#include "Herwig++/PDT/GenericMassGenerator.h"

using namespace Herwig;
using namespace ThePEG;

Hw64Decayer::~Hw64Decayer() {}

void Hw64Decayer::Init() {
   static ClassDocumentation<Hw64Decayer> documentation
     ("Class to decay all particles in HERWIG by the algorithms used in HERWIG 6.4");

   static Parameter<Hw64Decayer, int>
      interfaceME("MECode", "The code for the ME type to use in the decay",
                  &Herwig::Hw64Decayer::MECode, 0, 0, 0, true, false, false);

  static Parameter<Hw64Decayer,unsigned int> interfaceMassTry
    ("MassTry",
     "The maximum number of attempts to generate the off-shell masses of the"
     " decay products.",
     &Hw64Decayer::_masstry, 50, 1, 1000,
     false, false, Interface::limited);

}

ClassDescription<Hw64Decayer> Hw64Decayer::initHw64Decayer;

bool Hw64Decayer::accept(const DecayMode &dm) const { return true; }

/****** 
 * This function actually decays a particle based on the dm given.
 *****/
ParticleVector Hw64Decayer::decay(const DecayMode &dm, const Particle &p) const
{
  // storage for the decay products and number of decay products
   ParticleVector rval;
   unsigned int numProds(dm.products().size());
   // can't handle more than 5 products throw a veto
   if(numProds > 5) { 
     generator()->log() << "Number of decay products is greater than 5\n"
                        << "Hw64Decayer cannot handle these decays properly. "
                        << "Will veto decay and select a new decay mode\n";
     throw Veto();
   }
   // check that it is possible to kinematically perform the decay
   Energy minmass(0.);
   vector<Energy> minmasses(numProds);
   vector<tcGenericMassGeneratorPtr> massgen(numProds,tcGenericMassGeneratorPtr());
   tcMassGenPtr mtemp;
   for(unsigned int ix=0;ix<numProds;++ix)
     {
       minmasses[ix]=dm.orderedProducts()[ix]->massMin();
       minmass+=minmasses[ix];
       mtemp=dm.orderedProducts()[ix]->massGenerator();
       if(mtemp){massgen[ix]=dynamic_ptr_cast<tcGenericMassGeneratorPtr>(mtemp);}
     }
   // throw a veto if not kinematically possible
   if(minmass>p.mass())
     {
       generator()->log() << "Hw64Decayer cannot perform the decay " << dm.tag() 
			  << " for this particle instance as the minimum mass of "
			  << "the decay products exceeds the mass of the particle"
			  << endl;
       throw Veto();
     }
   // check not decaying a massless particle
   if(p.mass() < 0.000001) {
     generator()->log() << "HwDecayer called on a particle with no mass " 
			<< p.PDGName() << ", " << p.mass() << endl;
     throw Veto();
   }
   // Create a vectors for momenta and masses
   vector<Lorentz5Momentum> products(numProds);
   vector<Energy> masses(numProds);
   // now generate the masses of the particles starting with a random one
   // to avoid bias
   unsigned int ntry(0);
   Energy outmass;
   do
     {
       unsigned int istart=UseRandom::irnd(numProds);
       outmass=0.;
       for(unsigned int ix=istart;ix<numProds;++ix)
	 {
	   if(massgen[ix])
	     {masses[ix]=massgen[ix]->mass(*(dm.orderedProducts()[ix]),minmasses[ix],
					   p.mass()-minmass+minmasses[ix]);}
	   else
	     {masses[ix]=dm.orderedProducts()[ix]->generateMass();}
	   outmass+=masses[ix];
	   if(outmass>p.mass()){break;}
	 }
       for(unsigned int ix=0;ix<istart;++ix)
	 {
	   if(massgen[ix])
	     {masses[ix]=massgen[ix]->mass(*(dm.orderedProducts()[ix]),minmasses[ix],
					   p.mass()-minmass+minmasses[ix]);}
	   else
	     {masses[ix]=dm.orderedProducts()[ix]->generateMass();}
	   outmass+=masses[ix];
	   if(outmass>p.mass()){break;}
	 }
     }
   while(ntry<_masstry&&outmass>p.mass());
   if(outmass>p.mass())
     {
       generator()->log() << "Hw64Decayer failed to generate the masses of the"
			  << " decay products for the decay " << dm.tag() 
			  << " after " << _masstry << " attempts";
       throw Veto();
     }
   for(unsigned int ix=0;ix<numProds;++ix){products[ix].setMass(masses[ix]);}
   // The K -> KL0 and KS0
   if(numProds == 1)
     {oneBodyDecay(p.momentum(), products[0]);}
   // 2 Body Decay
   else if(numProds == 2) 
     {
       double CosAngle, AzmAngle; 
       Kinematics::generateAngles(CosAngle, AzmAngle);
      /*******
       * It appears polarized mesons aren't needed currently 
      if(p->id() == ) {
         if() {
            CosAng = Math::absmax<double>(UseRandom::rnd(), 
		           Math::absmax<double>(UseRandom::rnd(),  UseRandom::rnd()))*2.0 - 1.0;
         } else if() {
            CosAng = 2.0*cos((acos(UseRandom::rnd()*2.0-1.0)+M_PI)/3.0);
         } else {
            CosAng = Math::absmin<double>(UseRandom::rnd(),
                           Math::absmin<double>(UseRandom::rnd(),  UseRandom::rnd()))*2.0 - 1.0;
         }
      }
       *******/
      Kinematics::twoBodyDecay(p.momentum(), masses[0], masses[1],
			       CosAngle, AzmAngle, products[0], products[1]);
      
      // Now we rotate the result (why?)
      // LorentzRotation R = Kinematics::rotate(p.momentum(), 1.0, 0.0);
      // products[0] *= R; products[1]*=R;
     }
   // Three Body Decay
   else if(numProds == 3) 
     {
       // Free Massless (V-A)*(V-A) ME
       if(MECode == 100) 
	 {Kinematics::threeBodyDecay(p.momentum(), products[0], products[1], 
				     products[2], &VAWt);} 
       // Bound Massless (V-A)*(V-A) ME
       else if(MECode == 101) 
	 {
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
	   do 
	     {
	       Kinematics::threeBodyDecay(p.momentum(), products[1], products[2], 
					  products[0], &VAWt);
	       dot1 = p.momentum().dot(products[1]);
	       dot2 = p.momentum().dot(products[0]);
	     } 
	   while(dot1*(wtmx-dot1-xs*dot2) < UseRandom::rnd()*wtmx2);
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
       // Three Body via phase space
       else
	 {Kinematics::threeBodyDecay(p.momentum(), products[0], products[1], 
				     products[2]);}
     }  
   // Four Body Decay
   else if(numProds == 4) 
     {Kinematics::fourBodyDecay(p.momentum(), products[0], products[1], 
				products[2], products[3]);}
   // Five Body Decay
   else if(numProds == 5) 
     {Kinematics::fiveBodyDecay(p.momentum(), products[0], products[1], 
				products[2], products[3], products[4]);}
   if(products[0] == Lorentz5Momentum()) {
     generator()->log() << "The Decay mode " << dm.tag() << " cannot "
	                << "proceed, not enough phase space\n";   
     return ParticleVector();
   }
   cPDVector productParticles(numProds);
   for(unsigned int ix=0;ix<numProds;++ix)
     {productParticles[ix]=dm.orderedProducts()[ix];}
   // set the momenta of the particles and return the answer
   setParticleMomentum(rval, productParticles, products);
   return rval;
}
   
void Hw64Decayer::persistentOutput(PersistentOStream &os) const 
{os << MECode << _masstry;}
void Hw64Decayer::persistentInput(PersistentIStream &is, int i)
{is >> MECode >> _masstry;}

/******
 * This function takes the array of momentum generated and sets the momentum
 * to the particles.
 *****/
void Hw64Decayer::setParticleMomentum(ParticleVector &out, cPDVector particles, 
                                      vector<Lorentz5Momentum> moms) const 
{
   unsigned int numProds = particles.size();
   for(unsigned int ix=0;ix<numProds;++ix)
     {out.push_back(particles[ix]->produceParticle(moms[ix]));}
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


