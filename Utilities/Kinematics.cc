// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Kinematics class.
//

#include "Kinematics.h"
#include <ThePEG/Vectors/Lorentz5Vector.h>
#include <ThePEG/Vectors/LorentzVector.h>
#include <ThePEG/Vectors/LorentzRotation.h>
#include <ThePEG/Repository/EventGenerator.h>
#include <ThePEG/Repository/CurrentGenerator.h>
#include <ThePEG/EventRecord/Event.h>

using namespace Herwig;
using namespace ThePEG;

Energy Kinematics::CMMomentum(const Energy M, 
			      const Energy m1, 
			      const Energy m2) {
  return ( M <= Energy()  ||  m1 < Energy()  
	   ||  m2 < Energy()  ||  M <= m1+m2  ?  
	   Energy()  : 
	   Energy(sqrt(( M*M - (m1+m2)*(m1+m2) )*( M*M - (m1-m2)*(m1-m2) ))
	   /(2.0*M))); 
}

bool Kinematics::twoBodyDecay(const Lorentz5Momentum & p,
			      const Energy m1, const Energy m2,
			      const Axis & unitDir1,
			      Lorentz5Momentum & p1, Lorentz5Momentum & p2 ) {
  Energy min=p.m();
  if ( min >= m1 + m2  &&  m1 >= Energy()  &&  m2 >= Energy()  ) {
    Momentum3 pstarVector = unitDir1 * pstarTwoBodyDecay(min,m1,m2);
    p1 = Lorentz5Momentum(m1, pstarVector);
    p2 = Lorentz5Momentum(m2,-pstarVector);
    // boost from CM to LAB
    Boost bv=p.boostVector();
    p1.boost( bv );   
    p2.boost( bv );
    return true;
  } else {
    CurrentGenerator::log() 
      << "Kinematics::twoBodyDecay() phase space problem\n" 
      << "p = " << p / GeV 
      << " p.m() = " << min / GeV
      << " -> " << m1/GeV 
      << ' ' << m2/GeV << '\n';
    return false;
  }
}

bool Kinematics::twoBodyDecay(const Lorentz5Momentum & p,                    
			      const Energy m1, const Energy m2,  
			      const double cosThetaStar1, 
			      const double phiStar1, 
			      Lorentz5Momentum & p1, Lorentz5Momentum & p2 ) {
  return twoBodyDecay(p,m1,m2,unitDirection(cosThetaStar1,phiStar1),p1,p2); 
}

void Kinematics::generateAngles(double &ct, double &az) {
  ct = UseRandom::rnd()*2.0 - 1.0;  // Flat from -1..1
  az = UseRandom::rnd()*2.0*Constants::pi;   // Flat from 0..2 Pi
}

/*****
 * This function, as the name implies, performs a three body decay. The decay
 * products are distributed uniformly in all three directions.
 ****/
bool Kinematics::threeBodyDecay(Lorentz5Momentum p0, Lorentz5Momentum &p1, 
				Lorentz5Momentum &p2, Lorentz5Momentum &p3,
				double (*fcn)(Energy2,Energy2,Energy2,InvEnergy4)) {
  // Variables needed in calculation...named same as fortran version
  Energy a = p0.mass() + p1.mass();
  Energy b = p0.mass() - p1.mass();
  Energy c = p2.mass() + p3.mass();
  
  if(b < c) {
     CurrentGenerator::log() 
       << "Kinematics::threeBodyDecay() phase space problem\n"
       << p0.mass()/GeV << " -> "
       << p1.mass()/GeV << ' '
       << p2.mass()/GeV << ' '
       << p3.mass()/GeV << '\n';
     return false;
  }
  
  Energy d = abs(p2.mass()-p3.mass());
  Energy2 aa = sqr(a); 
  Energy2 bb = sqr(b); 
  Energy2 cc = sqr(c); 
  Energy2 dd = sqr(d); 
  Energy2 ee = (b-c)*(a-d);
  
  Energy2 a1 = 0.5 * (aa+bb);
  Energy2 b1 = 0.5 * (cc+dd);
  InvEnergy4 c1 = 4./(sqr(a1-b1));
  
  Energy2 ff; 
  double ww; 
  Energy4 pp,qq,rr;
  // Choose mass of subsystem 23 with prescribed distribution
  do {
    // ff is the mass squared of the 23 subsystem
    ff = UseRandom::rnd()*(cc-bb)+bb;
    
    // pp is ((m0+m1)^2 - m23^2)((m0-m1)^2-m23)
    pp = (aa-ff)*(bb-ff);
    
    // qq is ((m2+m3)^2 - m23^2)(|m2-m3|^2-m23^2)
    qq = (cc-ff)*(dd-ff);
    
    if(fcn != NULL) {
      ww = (*fcn)(ff,a1,b1,c1);
    } else ww = 1.0;
    //if(MECode == 100 || MECode == 101) ww = EMMasslessWt(ff,a1,b1,c1);
    //else ww =  PhaseSpaceWt();
    ww = sqr(ww);
    rr = ee*ff*UseRandom::rnd();
  } while(pp*qq*ww < rr*rr);
  
  // ff is the mass squared of subsystem 23
  // do 2 body decays 0->1+23, 23->2+3
  double CosAngle, AzmAngle;
  Lorentz5Momentum p23;
  
  p23.setMass(sqrt(ff));
  
  generateAngles(CosAngle,AzmAngle);
  bool status = twoBodyDecay(p0,p1.mass(),p23.mass(),CosAngle,AzmAngle,p1,p23);
  
  generateAngles(CosAngle,AzmAngle);
  status &= twoBodyDecay(p23,p2.mass(),p3.mass(),CosAngle,AzmAngle,p2,p3);
  return status;
}


/******
 * Again, as the name implies this routine takes a parent particle, p0,
 * and decays it into four particles, given by p1..p4, which are distributed
 * uniformly
 ******/
bool Kinematics::fourBodyDecay(Lorentz5Momentum  p0, Lorentz5Momentum &p1,
			       Lorentz5Momentum &p2, Lorentz5Momentum &p3,
			       Lorentz5Momentum &p4) {
  // Again, for comparison, we use the same variables as used in fortran code
  //  double b,c,aa,bb,cc,dd,ee,tt,s1,rs1,ff,s2,pp,qq,rr, temp;

  Energy b = p0.mass() - p1.mass();
  Energy c = p2.mass() + p3.mass() + p4.mass();
  if(b < c) {
    CurrentGenerator::log() 
      << "Kinematics::fourBodyDecay() phase space problem\n";
    return false;
  }
  Energy2 aa = sqr(p0.mass() + p1.mass());
  Energy2 bb = sqr(b);
  Energy2 cc = sqr(c);
  Energy2 dd = sqr(p3.mass()+p4.mass());
  Energy2 ee = sqr(p3.mass()-p4.mass());
  Energy8 tt = (b-c)*pow<7,1>(p0.mass())/16.;
  
  // Select squared masses, S1 and S2 or 234 and 34 subsystems
  Energy2 s1, ff, s2, qq, rr;
  Energy rs1;
  Energy4 pp;
  double temp;
  do {
    s1 = bb + UseRandom::rnd()*(cc-bb);
    rs1 = sqrt(s1);
    ff = sqr(rs1-p2.mass());
    s2 = dd + UseRandom::rnd()*(ff-dd);
    pp = (aa-s1)*(bb-s1);
    qq = (sqr(rs1+p2.mass())-s2)*(ff-s2)/s1;
    rr = (s2 - dd)*(s2-ee)/s2;

    temp = UseRandom::rnd();
  } while(pp*qq*rr*sqr(ff-dd) < tt*s1*s2*sqr(temp));

  // Now we have chosen how energy fractions go, do two body decays on subsystems
  // in order to get angles 
  double CosAngle, AzmAngle;
  Lorentz5Momentum p234, p34;

  p234.setMass(rs1);
  p34.setMass(sqrt(s2));

  generateAngles(CosAngle,AzmAngle); 
  bool status = twoBodyDecay(p0,p1.mass(),p234.mass(),CosAngle, AzmAngle, p1, p234);

  generateAngles(CosAngle,AzmAngle); 
  status &= twoBodyDecay(p234,p2.mass(),p34.mass(),CosAngle, AzmAngle, p2, p34);

  generateAngles(CosAngle,AzmAngle); 
  status &= twoBodyDecay(p34,p3.mass(),p4.mass(),CosAngle, AzmAngle, p3, p4);
  return status;
}

/******
 * Lastly, as the name implies this routine takes a parent particle, p0,
 * and decays it into five particles, given by p1..p5, which are distributed
 * uniformly
 ******/
bool Kinematics::fiveBodyDecay(Lorentz5Momentum  p0, Lorentz5Momentum &p1,
			       Lorentz5Momentum &p2, Lorentz5Momentum &p3,
			       Lorentz5Momentum &p4, Lorentz5Momentum &p5) {
  Energy b(p0.mass()-p1.mass());
  Energy c(p2.mass()+p3.mass()+p4.mass()+p5.mass());
  if(b<c) {
    CurrentGenerator::log() 
      << "Kinematics::fiveBodyDecay() phase space problem\n";
    return false;
  }
  Energy2 aa(sqr(p0.mass()+p1.mass()));
  Energy2 bb(b*b),cc(c*c);
  Energy2 dd(sqr(p3.mass()+p4.mass()+p5.mass()));
  Energy2 ee(sqr(p4.mass()+p5.mass()));
  Energy2 ff(sqr(p4.mass()-p5.mass()));
  Energy12 tt((b-c)*pow<11,1>(p0.mass())/729.);
  // Select squared masses S1, S2 and S3 of 2345, 345 and 45 subsystems
  //double s1,rs1,s2,rs2,gg,hh,s3,pp,qq,rr,ss,temp;
  Energy2 s1, gg, s2, hh, s3, qq, rr, ss;
  Energy rs1, rs2;
  Energy4 pp;
  double temp;
  // protect against infinite loops
  static const unsigned int MAXTRY = 100;
  unsigned int count(0);
  do
    {
      s1=bb+UseRandom::rnd()*(cc-bb);
      rs1=sqrt(s1);
      gg=sqr(rs1-p2.mass());
      s2=dd+UseRandom::rnd()*(gg-dd);
      rs2=sqrt(s2);
      hh=sqr(rs2-p3.mass());
      s3=ee+UseRandom::rnd()*(hh-ee);
      pp=(aa-s1)*(bb-s1);
      qq=(sqr(rs1+p2.mass())-s2)*(gg-s2)/s1;
      rr=(sqr(rs2+p3.mass())-s3)*(hh-s3)/s2;
      ss=(s3-ee)*(s3-ff)/s3;
      temp = UseRandom::rnd();
      ++count;
    }
  while(pp*qq*rr*ss*sqr((gg-dd)*(hh-ee)) < tt*s1*s2*s3*sqr(temp) && count<MAXTRY);
  if(count==MAXTRY)
    {
      CurrentGenerator::log() << "Kinematics::fiveBodyDecay can't generate momenta" 
			      << " after " << MAXTRY << " attempts\n";
      throw Veto();
    }
  // Now decay the subsystems
  double CosAngle, AzmAngle;
  Lorentz5Momentum p2345, p345, p45;

  p2345.setMass(rs1);
  p345.setMass(rs2);
  p45.setMass(sqrt(s3));

  generateAngles(CosAngle,AzmAngle); 
  bool status = twoBodyDecay(p0,p1.mass(),p2345.mass(),CosAngle, AzmAngle, p1, p2345);
  
  generateAngles(CosAngle,AzmAngle); 
  status &= twoBodyDecay(p2345,p2.mass(),p345.mass(),CosAngle, AzmAngle, p2, p345);
  
  generateAngles(CosAngle,AzmAngle); 
  status &= twoBodyDecay(p345,p3.mass(),p45.mass(),CosAngle, AzmAngle, p3, p45);
  
  generateAngles(CosAngle,AzmAngle); 
  status &= twoBodyDecay(p45,p4.mass(),p5.mass(),CosAngle, AzmAngle, p4, p5);
  return status;
}


