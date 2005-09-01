// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Kinematics class.
//

#include "Kinematics.h"
#include <ThePEG/CLHEPWrap/Lorentz5Vector.h>
#include <ThePEG/CLHEPWrap/LorentzVector.h>
#include <ThePEG/CLHEPWrap/LorentzRotation.h>
#include <ThePEG/Repository/EventGenerator.h>
#include <ThePEG/Repository/CurrentGenerator.h>

using namespace Herwig;
using namespace ThePEG;

Energy Kinematics::CMMomentum(const Energy M, 
			      const Energy m1, 
			      const Energy m2) {
  return ( M <= 0.0  ||  m1 < 0.0  ||  m2 < 0.0  ||  M <= m1+m2  ?  0.0  : 
	   sqrt(( M*M - (m1+m2)*(m1+m2) )*( M*M - (m1-m2)*(m1-m2) ))/(2.0*M)); 
}

void Kinematics::twoBodyDecay(const Lorentz5Momentum & p,
			      const double m1, const double m2,
			      const Vector3 & unitDir1,
			      Lorentz5Momentum & p1, Lorentz5Momentum & p2 ) {
  if ( p.m() >= m1 + m2  &&  m1 >= 0.0  &&  m2 >= 0.0  ) {
    Momentum3 pstarVector = unitDir1;
    pstarVector *= pstarTwoBodyDecay(p.m(),m1,m2);
    p1 = Lorentz5Momentum(m1,pstarVector);
    p2 = Lorentz5Momentum(m2,-pstarVector);
    p1.boost( p.boostVector() );   // boost from CM to LAB
    p2.boost( p.boostVector() );
  }
}

void Kinematics::twoBodyDecay(const Lorentz5Momentum & p,                    
			      const double m1, const double m2,  
			      const double cosThetaStar1, 
			      const double phiStar1, 
			      Lorentz5Momentum & p1, Lorentz5Momentum & p2 ) {
  twoBodyDecay(p,m1,m2,unitDirection(cosThetaStar1,phiStar1),p1,p2); 
}
/* Doesn't work with new LorentzRotation and isn't used anyway PR 1/9/05
     *
     * This returns the Matrix that rotates along the direction p by
     * phi (cp = cos(phi), sp = sin(phi))
     *
    static LorentzRotation rotation(const Lorentz5Momentum p,
				    const double cp, 
				    const double sp);

// Utility to create a rotation matrix to take p to the z axis, rotate it
// by psi about the z axis (cp = cos(psi), sp = sin(psi))
LorentzRotation Kinematics::rotation(const Lorentz5Momentum p,
				     const double cp, 
				     const double sp) 
{
   double pt, pp, ct, st, cf, sf;

#define PTCUT 1E-20
#define sq(a) (a)*(a)

   pt = sq(p.x())+ sq(p.y());
   pp = sq(p.z()) + pt;
   if(pt < pp*PTCUT) {
      ct = p.z() > 0 ? 1.0 : -1.0;
      st = 0.0; cf = 1.0; sf = 0.0;
   } else {
      pp = sqrt(pp); pt = sqrt(pt); 
      ct = p.z()/pp; st = pt/pp; cf = p.x()/pt; sf = p.y()/pt;
   }

#undef PTCUT
#undef sq

   return 
     LorentzRotation(LorentzVector(cp*cf*ct+sp*sf, -cp*sf+sp*cf*ct, cf*st, 0),
		     LorentzVector(cp*sf*ct-sp*cf,  cp*cf+sp*sf*ct, sf*st, 0),
		     LorentzVector(-cp*st,          -sp*st,         ct,    0),
		     LorentzVector(0,               0,              0,     1));
}
*/
void Kinematics::generateAngles(double &ct, double &az) {
  ct = CurrentGenerator::current().rnd()*2.0 - 1.0;  // Flat from -1..1
  az = CurrentGenerator::current().rnd()*2.0*M_PI;   // Flat from 0..2 Pi
}

/*****
 * This function, as the name implies, performs a three body decay. The decay
 * products are distributed uniformly in all three directions.
 ****/
void Kinematics::threeBodyDecay(Lorentz5Momentum p0, Lorentz5Momentum &p1, 
				Lorentz5Momentum &p2, Lorentz5Momentum &p3,
				double (*fcn)(double*)) {
    // Variables needed in calculation...named same as fortran version
   double a,b,c,d,aa,bb,cc,dd,ee,ff,pp,qq,ww,rr;   

   a = p0.mass() + p1.mass();
   b = p0.mass() - p1.mass();
   c = p2.mass() + p3.mass();

   if(b < c) {
     CurrentGenerator::log() << "***Kinematics::threeBodyDecay***"
			     << " not enough phase space!\n";
   }

   d = fabs(p2.mass()-p3.mass());
   aa = sqr(a); bb = sqr(b); cc = sqr(c); dd = sqr(d); ee = (b-c)*(a-d);
   a = 0.5 * (aa+bb);
   b = 0.5 * (cc+dd);
   c = 4./(sqr(a-b));

   // Choose mass of subsystem 23 with prescribed distribution
   do {
      // ff is the mass squared of the 23 subsystem
      ff = CurrentGenerator::current().rnd()*(cc-bb)+bb;

      // pp is ((m0+m1)^2 - m23^2)((m0-m1)^2-m23)
      pp = (aa-ff)*(bb-ff);

      // qq is ((m2+m3)^2 - m23^2)(|m2-m3|^2-m23^2)
      qq = (cc-ff)*(dd-ff);

      if(fcn != NULL) {
	double temp[4] = {ff,a,b,c};
	ww = (*fcn)(temp);
      } else ww = 1.0;
      //if(MECode == 100 || MECode == 101) ww = EMMasslessWt(ff,a,b,c);
      //else ww =  PhaseSpaceWt();
      ww = sqr(ww);
      rr = ee*ff*CurrentGenerator::current().rnd();
   } while(pp*qq*ww < rr*rr);

   // ff is the mass squared of subsystem 23
   // do 2 body decays 0->1+23, 23->2+3
   double CosAngle, AzmAngle;
   Lorentz5Momentum p23;

   p23.setMass(sqrt(ff));

   generateAngles(CosAngle,AzmAngle);
   twoBodyDecay(p0,p1.mass(),p23.mass(),CosAngle,AzmAngle,p1,p23);

   generateAngles(CosAngle,AzmAngle);
   twoBodyDecay(p23,p2.mass(),p3.mass(),CosAngle,AzmAngle,p2,p3);
}


/******
 * Again, as the name implies this routine takes a parent particle, p0,
 * and decays it into four particles, given by p1..p4, which are distributed
 * uniformly
 ******/
void Kinematics::fourBodyDecay(Lorentz5Momentum  p0, Lorentz5Momentum &p1,
			       Lorentz5Momentum &p2, Lorentz5Momentum &p3,
			       Lorentz5Momentum &p4) {
  // Again, for comparison, we use the same variables as used in fortran code
  double b,c,aa,bb,cc,dd,ee,tt,s1,rs1,ff,s2,pp,qq,rr, temp;

  b = p0.mass() - p1.mass();
  c = p2.mass() + p3.mass() + p4.mass();
  if(b < c) 
    CurrentGenerator::log() 
      << "Kinematics::fourBodyDecay: no phase space available\n";
  aa = sqr(p0.mass() + p1.mass());
  bb = sqr(b);
  cc = sqr(c);
  dd = sqr(p3.mass()+p4.mass());
  ee = sqr(p3.mass()-p4.mass());
  tt = (b-c)*pow(p0.mass(),7)/16.;
  
  // Select squared masses, S1 and S2 or 234 and 34 subsystems
  do {
    s1 = bb + CurrentGenerator::current().rnd()*(cc-bb);
    rs1 = sqrt(s1);
    ff = sqr(rs1-p2.mass());
    s2 = dd + CurrentGenerator::current().rnd()*(ff-dd);
    pp = (aa-s1)*(bb-s1);
    qq = (sqr(rs1+p2.mass())-s2)*(ff-s2)/s1;
    rr = (s2 - dd)*(s2-ee)/s2;

    // Since sqr is a macro, need to store rnd() in variable before sqr,
    // otherwise it will use two random numbers
    temp = CurrentGenerator::current().rnd();
  } while(pp*qq*rr*sqr(ff-dd) < tt*s1*s2*sqr(temp));

  // Now we have chosen how energy fractions go, do two body decays on subsystems
  // in order to get angles 
  double CosAngle, AzmAngle;
  Lorentz5Momentum p234, p34;

  p234.setMass(rs1);
  p34.setMass(sqrt(s2));

  generateAngles(CosAngle,AzmAngle); 
  twoBodyDecay(p0,p1.mass(),p234.mass(),CosAngle, AzmAngle, p1, p234);

  generateAngles(CosAngle,AzmAngle); 
  twoBodyDecay(p234,p2.mass(),p34.mass(),CosAngle, AzmAngle, p2, p34);

  generateAngles(CosAngle,AzmAngle); 
  twoBodyDecay(p34,p3.mass(),p4.mass(),CosAngle, AzmAngle, p3, p4);
}

/******
 * Lastly, as the name implies this routine takes a parent particle, p0,
 * and decays it into five particles, given by p1..p5, which are distributed
 * uniformly
 ******/
void Kinematics::fiveBodyDecay(Lorentz5Momentum  p0, Lorentz5Momentum &p1,
			       Lorentz5Momentum &p2, Lorentz5Momentum &p3,
			       Lorentz5Momentum &p4, Lorentz5Momentum &p5) {
  // once more, using the same variable names as fortran code
  double b,c,aa,bb,cc,dd,ee,ff,tt,s1,rs1,gg,s2,rs2,hh,s3,pp,qq,rr,ss,temp;

  b = p0.mass()-p1.mass();
  c = p2.mass()+p3.mass()+p4.mass()+p5.mass();
  if(b<c) 
    CurrentGenerator::log() 
      << "Kinematics::fiveBodyDecay - No phase space available\n";
  aa = sqr(p0.mass()+p1.mass());
  bb = sqr(b);
  cc = sqr(c);
  dd = sqr(p3.mass()+p4.mass()+p5.mass());
  ee = sqr(p4.mass()+p5.mass());
  ff = sqr(p4.mass()-p5.mass());
  tt = (b-c)*pow(p0.mass(),11)/729.;
  
  // Now lets select the masses of all the subsystems
  do {
    s1 = bb+CurrentGenerator::current().rnd()*(cc-bb);
    rs1 = sqrt(s1);
    gg = sqr(rs1-p2.mass());
    s2 = dd + CurrentGenerator::current().rnd()*(gg-dd);
    rs2 = sqrt(s2);
    hh = sqr(rs2-p3.mass());
    s3 = ee+CurrentGenerator::current().rnd()*(hh-ee);
    pp = (aa-s1)*(bb-s1);
    qq = (sqr(rs1+p2.mass())-s2)*(gg-s2)/s1;
    rr = (sqr(rs2+p3.mass())-s3)*(hh-s3)/s2;
    ss = (s3-ee)*(s3-ff)/s3;
    // Again, since sqr is a macro, the random number must be stored or
    // two different random numbers will be used
    temp = CurrentGenerator::current().rnd();
  } while(pp*qq*rr*qq*sqr((gg-dd)*(hh-ee)) < tt*s1*s2*s3*sqr(temp));

  // Now decay the subsystems
  double CosAngle, AzmAngle;
  Lorentz5Momentum p2345, p345, p45;

  p2345.setMass(rs1);
  p345.setMass(rs2);
  p45.setMass(sqrt(s3));

  generateAngles(CosAngle,AzmAngle); 
  twoBodyDecay(p0,p1.mass(),p2345.mass(),CosAngle, AzmAngle, p1, p2345);
  
  generateAngles(CosAngle,AzmAngle); 
  twoBodyDecay(p2345,p2.mass(),p345.mass(),CosAngle, AzmAngle, p2, p345);
  
  generateAngles(CosAngle,AzmAngle); 
  twoBodyDecay(p345,p3.mass(),p45.mass(),CosAngle, AzmAngle, p3, p45);
  
  generateAngles(CosAngle,AzmAngle); 
  twoBodyDecay(p45,p4.mass(),p5.mass(),CosAngle, AzmAngle, p4, p5);
}


