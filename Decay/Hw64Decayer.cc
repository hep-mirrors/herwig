// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Hw64Decayer class.
//

#include "Hw64Decayer.h"
#include "Pythia7/PDT/EnumParticles.h"
#include "Pythia7/PDT/DecayMode.h"
#include "Pythia7/Interface/ClassDocumentation.h"
#include "Pythia7/Interface/Parameter.h"
#include "Herwig++/Utilities/Kinematics.h"
#include "Pythia7/Repository/EventGenerator.h"
#include "Pythia7/Repository/Repository.h"
#include "Pythia7/Utilities/Rebinder.h"
#include "Pythia7/Interface/Reference.h"
#include "Pythia7/Interface/Parameter.h"
#include "Pythia7/Persistency/PersistentOStream.h"
#include "Pythia7/Persistency/PersistentIStream.h"
#include "Pythia7/Utilities/Math.h"
#include "Herwig++/Utilities/HwDebug.h"
#include "Kinematics.h"

using namespace Herwig;
using namespace Pythia7;

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
     generator()->log() << "Number of particles to decay is greater then 5\n"
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
   // The K -> KL0 and KS0
   if(numProds == 1) {
      oneBodyDecay(p.momentum(), products[0]);

   // 2 Body Decay
   } else if(numProds == 2) {
      double CosAngle, AzmAngle;
      
      MyKinematics::generateAngles(generator(), CosAngle, AzmAngle);

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
      MyKinematics::twoBodyDecay(p.momentum(), masses[0], masses[1],
			         CosAngle, AzmAngle, products[0], products[1]);

      // Now we rotate the result (why?)
      // LorentzRotation R = Kinematics::rotate(p.momentum(), 1.0, 0.0);
      // products[0] *= R; products[1]*=R;
 
   // Three Body Decay
   } else if(numProds == 3) {
      // Free Massless (V-A)*(V-A) ME
      if(MECode == 100) {
         threeBodyDecay(p.momentum(), products[0], products[1], products[2]);
         
      // Bound Massless (V-A)*(V-A) ME
      } else if(MECode == 101) {
         double wtmx, wtmx2, xs, dot1, dot2;
         wtmx = ( (p.mass() - products[2][5])*(p.mass() + products[2][5])
                + (products[1][5] - products[0][5])*(products[1][5] + products[0][5]))/2.0;
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
            threeBodyDecay(p.momentum(), products[1], products[2], products[0]);
            dot1 = p.momentum().dot(products[2]);
            dot2 = p.momentum().dot(products[1]);
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

      else threeBodyDecay(p.momentum(), products[0], products[1], products[2]);

   // Four Body Decay
   } else if(numProds == 4) {
     fourBodyDecay(p.momentum(), products[0], products[1], products[2],
		   products[3]);

   // Five Body Decay
   } else if(numProds == 5) {
     fiveBodyDecay(p.momentum(), products[0], products[1], products[2],
		   products[3], products[4]);
   }

   setParticleMomentum(rval, productParticles, products);
 
   // If there is a gluon in the products, put it first. Do to the color connections inherent
   // in pythia7, this is a simple hack to work most of the time. Needs to be fixed/checked. 
   reorderProducts(rval);

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

void Hw64Decayer::reorderProducts(ParticleVector &out) const {
   ParticleVector in = out;
   out.clear();
   for(unsigned int i = 0; i<in.size(); i++) {
      if(in[i]->id() == ParticleID::g) out.push_back(in[i]);
      else out.insert(out.begin(),in[i]);
   }
}

double Hw64Decayer::EMMasslessWt(double emsq, double a, double b, double c) const {
   return (a-emsq)*(emsq-b)*c;
}

double Hw64Decayer::PhaseSpaceWt() const {
   return 1.0;
}

void Hw64Decayer::oneBodyDecay(Lorentz5Momentum p0, Lorentz5Momentum &p1) const 
{
   p1 = (LorentzVector)p0;
}

/*****
 * This function, as the name implies, performs a three body decay. The decay
 * products are distributed uniformly in all three directions.
 ****/
void Hw64Decayer::threeBodyDecay(Lorentz5Momentum p0, Lorentz5Momentum &p1, 
				 Lorentz5Momentum &p2, Lorentz5Momentum &p3) const {
    // Variables needed in calculation...named same as fortran version
   double a,b,c,d,aa,bb,cc,dd,ee,ff,pp,qq,ww,rr;   

   a = p0.mass() + p1.mass();
   b = p0.mass() - p1.mass();
   c = p2.mass() + p3.mass();

   if(b<c) {
      generator()->log() << "No Phase space available for decay\n";
   }

   d = fabs(p2.mass()-p3.mass());
   aa = sqr(a); bb = sqr(b); cc = sqr(c); dd = sqr(d); ee = (b-c)*(a-d);
   a = 0.5 * (aa+bb);
   b = 0.5 * (cc+dd);
   c = 4./(sqr(a-b));

   // Choose mass of subsystem 23 with prescribed distribution
   do {
      // ff is the mass squared of the 23 subsystem
      ff = generator()->rnd()*(cc-bb)+bb;

      // pp is ((m0+m1)^2 - m23^2)((m0-m1)^2-m23)
      pp = (aa-ff)*(bb-ff);

      // qq is ((m2+m3)^2 - m23^2)(|m2-m3|^2-m23^2)
      qq = (cc-ff)*(dd-ff);

      if(MECode == 100 || MECode == 101) ww = EMMasslessWt(ff,a,b,c);
      else ww =  PhaseSpaceWt();
      ww = sqr(ww);
      rr = ee*ff*generator()->rnd();
   } while(pp*qq*ww < rr*rr);

   // ff is the mass squared of subsystem 23
   // do 2 body decays 0->1+23, 23->2+3
   double CosAngle, AzmAngle;
   Lorentz5Momentum p23;

   p23.setMass(sqrt(ff));

   MyKinematics::generateAngles(generator(),CosAngle,AzmAngle);
   MyKinematics::twoBodyDecay(p0,p1.mass(),p23.mass(),CosAngle,AzmAngle,p1,p23);

   MyKinematics::generateAngles(generator(),CosAngle,AzmAngle);
   MyKinematics::twoBodyDecay(p23,p2.mass(),p3.mass(),CosAngle,AzmAngle,p2,p3);

   // Some code to check distribution of angles...
   // Now lets find dist of angle between 1 and 2
   //double L1 = sqrt((p1.px()*p1.px()) + (p1.py()*p1.py()) + (p1.pz()*p1.pz()));
   //double L2 = sqrt((p2.px()*p2.px()) + (p2.py()*p2.py()) + (p2.pz()*p2.pz()));
   //double L3 = sqrt((p3.px()*p3.px()) + (p3.py()*p3.py()) + (p3.pz()*p3.pz()));
   //double L23 = sqrt((p23.px()*p23.px()) + (p23.py()*p23.py()) + (p23.pz()*p23.pz()));
   //double dot12 = (p1.px()*p2.px()) + (p1.py()*p2.py()) + (p1.pz()*p2.pz());
   //double dot13 = (p1.px()*p3.px()) + (p1.py()*p3.py()) + (p1.pz()*p3.pz());
   //double dot123 = (p1.px()*p23.px()) + (p1.py()*p23.py()) + (p1.pz()*p23.pz());
   //double dot23 = (p2.px()*p3.px()) + (p2.py()*p3.py()) + (p2.pz()*p3.pz());
   //cout << dot23/(L2*L3) << "\t" << "1.0" << endl;
}

/******
 * Again, as the name implies this routine takes a parent particle, p0,
 * and decays it into four particles, given by p1..p4, which are distributed
 * uniformly
 ******/
void Hw64Decayer::fourBodyDecay(Lorentz5Momentum  p0, Lorentz5Momentum &p1,
				Lorentz5Momentum &p2, Lorentz5Momentum &p3,
				Lorentz5Momentum &p4) const {
  // Again, for comparison, we use the same variables as used in fortran code
  double b,c,aa,bb,cc,dd,ee,tt,s1,rs1,ff,s2,pp,qq,rr, temp;

  b = p0.mass() - p1.mass();
  c = p2.mass() + p3.mass() + p4.mass();
  if(b < c) 
    generator()->log() << "Hw64Decayer::fourBodyDecay: no phase space available\n";
  aa = sqr(p0.mass() + p1.mass());
  bb = sqr(b);
  cc = sqr(c);
  dd = sqr(p3.mass()+p4.mass());
  ee = sqr(p3.mass()-p4.mass());
  tt = (b-c)*pow(p0.mass(),7)/16.;
  
  // Select squared masses, S1 and S2 or 234 and 34 subsystems
  do {
    s1 = bb + generator()->rnd()*(cc-bb);
    rs1 = sqrt(s1);
    ff = sqr(rs1-p2.mass());
    s2 = dd + generator()->rnd()*(ff-dd);
    pp = (aa-s1)*(bb-s1);
    qq = (sqr(rs1+p2.mass())-s2)*(ff-s2)/s1;
    rr = (s2 - dd)*(s2-ee)/s2;

    // Since sqr is a macro, need to store rnd() in variable before sqr,
    // otherwise it will use two random numbers
    temp = generator()->rnd();
  } while(pp*qq*rr*sqr(ff-dd) < tt*s1*s2*sqr(temp));

  // Now we have chosen how energy fractions go, do two body decays on subsystems
  // in order to get angles 
  double CosAngle, AzmAngle;
  Lorentz5Momentum p234, p34;

  p234.setMass(rs1);
  p34.setMass(sqrt(s2));

  MyKinematics::generateAngles(generator(),CosAngle,AzmAngle); 
  MyKinematics::twoBodyDecay(p0,p1.mass(),p234.mass(),CosAngle, AzmAngle, p1, p234);

  MyKinematics::generateAngles(generator(),CosAngle,AzmAngle); 
  MyKinematics::twoBodyDecay(p234,p2.mass(),p34.mass(),CosAngle, AzmAngle, p2, p34);

  MyKinematics::generateAngles(generator(),CosAngle,AzmAngle); 
  MyKinematics::twoBodyDecay(p34,p3.mass(),p4.mass(),CosAngle, AzmAngle, p3, p4);
}

/******
 * Lastly, as the name implies this routine takes a parent particle, p0,
 * and decays it into five particles, given by p1..p5, which are distributed
 * uniformly
 ******/
void Hw64Decayer::fiveBodyDecay(Lorentz5Momentum  p0, Lorentz5Momentum &p1,
				Lorentz5Momentum &p2, Lorentz5Momentum &p3,
				Lorentz5Momentum &p4, Lorentz5Momentum &p5) 
const {
  // once more, using the same variable names as fortran code
  double b,c,aa,bb,cc,dd,ee,ff,tt,s1,rs1,gg,s2,rs2,hh,s3,pp,qq,rr,ss,temp;

  b = p0.mass()-p1.mass();
  c = p2.mass()+p3.mass()+p4.mass()+p5.mass();
  if(b<c) 
    generator()->log() << "Hw64Decayer::fiveBodyDecay - No phase space available\n";
  aa = sqr(p0.mass()+p1.mass());
  bb = sqr(b);
  cc = sqr(c);
  dd = sqr(p3.mass()+p4.mass()+p5.mass());
  ee = sqr(p4.mass()+p5.mass());
  ff = sqr(p4.mass()-p5.mass());
  tt = (b-c)*pow(p0.mass(),11)/729.;
  
  // Now lets select the masses of all the subsystems
  do {
    s1 = bb+generator()->rnd()*(cc-bb);
    rs1 = sqrt(s1);
    gg = sqr(rs1-p2.mass());
    s2 = dd + generator()->rnd()*(gg-dd);
    rs2 = sqrt(s2);
    hh = sqr(rs2-p3.mass());
    s3 = ee+generator()->rnd()*(hh-ee);
    pp = (aa-s1)*(bb-s1);
    qq = (sqr(rs1+p2.mass())-s2)*(gg-s2)/s1;
    rr = (sqr(rs2+p3.mass())-s3)*(hh-s3)/s2;
    ss = (s3-ee)*(s3-ff)/s3;
    // Again, since sqr is a macro, the random number must be stored or
    // two different random numbers will be used
    temp = generator()->rnd();
  } while(pp*qq*rr*qq*sqr((gg-dd)*(hh-ee)) < tt*s1*s2*s3*sqr(temp));

  // Now decay the subsystems
  double CosAngle, AzmAngle;
  Lorentz5Momentum p2345, p345, p45;

  p2345.setMass(rs1);
  p345.setMass(rs2);
  p45.setMass(sqrt(s3));

  MyKinematics::generateAngles(generator(),CosAngle,AzmAngle); 
  MyKinematics::twoBodyDecay(p0,p1.mass(),p2345.mass(),CosAngle, AzmAngle, p1, p2345);
  
  MyKinematics::generateAngles(generator(),CosAngle,AzmAngle); 
  MyKinematics::twoBodyDecay(p2345,p2.mass(),p345.mass(),CosAngle, AzmAngle, p2, p345);
  
  MyKinematics::generateAngles(generator(),CosAngle,AzmAngle); 
  MyKinematics::twoBodyDecay(p345,p3.mass(),p45.mass(),CosAngle, AzmAngle, p3, p45);
  
  MyKinematics::generateAngles(generator(),CosAngle,AzmAngle); 
  MyKinematics::twoBodyDecay(p45,p4.mass(),p5.mass(),CosAngle, AzmAngle, p4, p5);
}
