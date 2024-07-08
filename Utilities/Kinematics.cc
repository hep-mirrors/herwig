// -*- C++ -*-
//
// Kinematics.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
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

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>


using namespace Herwig;
using namespace ThePEG;

/**
 * Boost consistently the Lorentz5Momenta in momenta (given in the pi and pj COM frame)
 * into the pi pj LAB frame (piLab,pjLab).
 * NOTE: the momenta must be given in the COM frame of piLab+pjLab, where piCOM points
 * 			 in the positive z-Axis and pjCOM in the negative z-Axis
 * */
void Kinematics::BoostIntoTwoParticleFrame(const Energy M, const Lorentz5Momentum & piLab,
		const Lorentz5Momentum & pjLab,
		std::vector<Lorentz5Momentum * > momenta) {
	if (piLab.vect().mag()==ZERO || pjLab.vect().mag()==ZERO) {
		if (piLab.vect().mag()==ZERO && pjLab.vect().mag()==ZERO)
			// even more trivial case where we are already in the correct frame
			return;
		// Trival 1+1 dimensional boost case
		Lorentz5Momentum pClu(M,(piLab+pjLab).vect());
		Boost bv = pClu.boostVector();
		bool PiNonZero = piLab.vect().mag()>ZERO;
		Lorentz5Momentum pConst = PiNonZero ? Lorentz5Momentum(piLab.m(),piLab.vect()):Lorentz5Momentum(pjLab.m(),pjLab.vect());
		Axis u=pConst.boost(-bv).vect().unit();
		for (unsigned int it = 0; it < momenta.size(); it++) {
			momenta[it]->rotateUz(u.unit());
			// mirror momenta if pjLab is considered as positive z-Axis
			if (!PiNonZero) momenta[it]->vect() = - momenta[it]->vect();
			momenta[it]->boost(bv);
		}
		return;
	}
	double cosPhi=piLab.vect().cosTheta(pjLab.vect());
	double Phi=acos(cosPhi);
	double sinPhi=sin(Phi);
	// If Phi==0 use regular 1+1D boost
	double epsilon=std::numeric_limits<double>::epsilon();
	if (fabs(cosPhi-1.0)<=epsilon || fabs(cosPhi+1.0)<=epsilon) {
		Lorentz5Momentum pClu(M,(piLab+pjLab).vect());
		Boost bv = pClu.boostVector();
		Lorentz5Momentum pConst = Lorentz5Momentum(piLab.m(),piLab.vect());
		Axis u=pConst.boost(-bv).vect().unit();
		for (unsigned int it = 0; it < momenta.size(); it++) {
			// rotate positive z-Axis into the original constituent direction of pi
			momenta[it]->rotateUz(u.unit());
			momenta[it]->boost(bv);
		}
		return;
	}
	Energy Ei=piLab.e();
	Energy Ej=pjLab.e();
	if (std::isnan(Phi) || std::isinf(Phi)) throw Exception() << "NAN or INF in Phi in Kinematics::BoostIntoTwoParticleFrame\n"
	<< Exception::runerror;
	Energy mi=piLab.mass();
	Energy mj=pjLab.mass();
	Energy2 mi2=mi*mi;
	Energy2 mj2=mj*mj;
	Energy Pi=piLab.vect().mag();
	Energy Pj=pjLab.vect().mag();
	assert(Pi>ZERO);

	std::vector<boost::numeric::ublas::vector<double>> momentaHat;

	std::vector<Energy> Masses;
	for (unsigned int it = 0; it < momenta.size(); it++) {
		Masses.push_back(momenta[it]->mass());
		momentaHat.push_back(boost::numeric::ublas::vector<double>(3));
		momentaHat[it](0) = momenta[it]->e()/GeV;
		momentaHat[it](1) = momenta[it]->x()/GeV;
		momentaHat[it](2) = momenta[it]->z()/GeV;
	}
	

	// Lorentz Matrix Lambda maps:
	// piHat = (Ecomi,0,0, Pcom) to piLab = (Ei,          0,0,Pi         )
	// pjHat = (Ecomj,0,0,-Pcom) to pjLab = (Ej,Pj*sin(phi),0,Pj*cos(phi))
	// and therefore maps also correctly momentaHat to momentaOut into the Lab frame
	boost::numeric::ublas::matrix<double> Lambda(3,3);
  
	Energy pstar=Kinematics::pstarTwoBodyDecay(M,mi,mj);
	Energy2 A2=pstar*M;
	Energy2 B2=Pi*Pj*sinPhi;
	Energy B2divPi=Pj*sinPhi;
	Energy2 Deltaij=Pi*Ej-Ei*Pj*cosPhi;
	double delta2=mi2*Pj*Pj*sinPhi*sinPhi/(Deltaij*Deltaij);
	double Lambda11=0;

	// better numerics	
	if (delta2<1e-13) Lambda11 = Deltaij>=ZERO ? (1.0-0.5*delta2):-(1.0-0.5*delta2);
	else if (Deltaij!=ZERO) Lambda11= Deltaij>=ZERO ? 1.0/sqrt(1.0+delta2):-1.0/sqrt(1.0+delta2);


	if (std::isnan(A2/GeV2) || std::isinf(A2/GeV2)) throw Exception() << "NAN in A2/GeV2\n"
			  << Exception::runerror;

	Lambda(0,0) = (Ei+Ej)/M;
	Lambda(0,1) = B2/A2;
	Lambda(0,2) = (Ei-Ej)/(2.0*pstar)-((mi2-mj2)*(Ei+Ej))/(2.0*M*A2);

	Lambda(1,0) = B2divPi/M;
	Lambda11 = Pi*Pj*(std::expm1(0.5*std::log1p(mj2/(Pj*Pj)))-std::expm1(0.5*std::log1p(mi2/(Pi*Pi)))+2*pow(sin(Phi/2.0),2)*sqrt(1.0+mi2/(Pi*Pi)))/(A2);
	Lambda(1,1) = Lambda11; // This should be just Deltaij/(M*pstar)
	Lambda(1,2) = -(M*M-(mj2-mi2))*B2divPi/(2.0*M*A2);

	Lambda(2,0) = (Pi+Pj*cosPhi)/M;
	Lambda(2,1) = Ei*B2divPi/A2;
	Lambda(2,2) = (A2*A2-0.5*Ei*(Ej*(M*M-(mj2-mi2))-Ei*(M*M-(mi2-mj2))))/(Pi*M*A2);

	Axis zAxis(0,0,1);
	Axis xAxis(1,0,0);
	Lorentz5Momentum piRes(mi,Pi*zAxis);
	Lorentz5Momentum pjRes(mj,Pj*(xAxis*sinPhi+zAxis*cosPhi));

	std::vector<Lorentz5Momentum> momentaRes;
	Lorentz5Momentum pClu1,pClu2;
	boost::numeric::ublas::vector<double> momentaOut(3);
	unsigned int iter = 0;
	bool isAligned;
	Momentum3 piHat(ZERO, ZERO,  pstar);
	Momentum3 pjHat(ZERO, ZERO, -pstar);
	Momentum3 pAligned(ZERO, ZERO, ZERO);
	// TODO FIX THIS ERROR consistently
	// Horrible fix below but works partially
	// TODO in this case just rescale piLab pjLab correspondingly, but the directions shall not change
	// pClu1 aligned with piLab
	for (auto & pHat : momentaHat)
	{
		isAligned = false;
		if (momenta[iter]->vect().mag()>ZERO) {
			if (
					fabs(1.0 - momenta[iter]->vect().cosTheta(piHat)) < 1.0e-14
				 ) {
				isAligned = true;
				double factor = momenta[iter]->z()/pstar;
				assert(momenta[iter]->z()/pstar > 0);
				double otherFactor = (momenta[iter]->e()-factor*sqrt(mi2+sqr(pstar)))/M;
				// doing the Boost into the Lab frame analytically:
				pAligned = (factor*piRes.vect() + otherFactor*(piRes.vect() + pjRes.vect()));
			} else if (
					fabs(1.0 - momenta[iter]->vect().cosTheta(pjHat)) < 1.0e-14
					) {	
				isAligned = true;
				double factor = fabs(momenta[iter]->z()/pstar);
				double otherFactor = (momenta[iter]->e()-factor*sqrt(mj2+sqr(pstar)))/M;
				// doing the Boost into the Lab frame analytically:
				pAligned = (factor*pjRes.vect() + otherFactor*(piRes.vect() + pjRes.vect()));
			}
		}
		// doing the Boost into the Lab frame:
		if (!isAligned) {
			momentaOut = boost::numeric::ublas::prod(Lambda,pHat);
			momentaRes.push_back(Lorentz5Momentum(Masses[iter],
						GeV*Axis(momentaOut(1), double(momenta[iter]->y()/GeV), momentaOut(2))));
		}
		else {
			momentaRes.push_back(Lorentz5Momentum(Masses[iter], pAligned));
		}
		iter++;
	}
	// Computing the correct rotation, which maps pi/jRes into pi/jLab
	Axis omega1=piRes.vect().unit().cross(piLab.vect().unit());
	double cosAngle1=piRes.vect().unit()*piLab.vect().unit();
	double angle1=acos(cosAngle1);

	if (omega1.mag() > ZERO){
		// Rotate piRes into piLab
		piRes.rotate(angle1, omega1);
		pjRes.rotate(angle1, omega1);
		// Correspondingly do the actual rotation on all momenta
		for(auto & pRes : momentaRes)
			pRes.rotate(angle1, omega1);
	}
	else {
		std::cout << "OMEGA1 == ZERO = " << omega1.mag()<< std::endl;
	}

	Axis omega2=piRes.vect().unit();
	Momentum3 r1dim=(pjLab.vect()-piRes.vect().unit()*(pjLab.vect()*piRes.vect().unit()));
	Momentum3 r2dim=(pjRes.vect()-piRes.vect().unit()*(pjRes.vect()*piRes.vect().unit()));

	if (r1dim.mag()==ZERO || r2dim.mag()==ZERO || fabs(sinPhi)<1e-14) //trivial rotation so we are done
	{
		for (unsigned int i = 0; i < momentaRes.size(); i++) {
			// copy the final momenta
			*(momenta[i]) = momentaRes[i];
		}
		return;
	}
	Axis r1=r1dim.unit();
	Axis r2=r2dim.unit();

	// signs for 2nd rotation
	int signToPi = (piRes.vect()*pjLab.vect())/GeV2 > 0 ? 1:-1;
	int signToR1R2 = signToPi*(r2.cross(r1)*piRes.vect())/GeV> 0 ? 1:-1;
	double angle2=acos(r1*r2);
	if (signToR1R2<0) angle2=-angle2;

	// Rotate pjRes into pjLab
	pjRes.rotate(angle2, signToPi*omega2);
	// Correspondingly do the actual rotation on all momenta
	for (unsigned int i = 0; i < momentaRes.size(); i++) {
		momentaRes[i].rotate(angle2, signToPi*omega2);
		// copy the final momenta
		*(momenta[i]) = momentaRes[i];
	}
}


bool Kinematics::twoBodyDecay(const Lorentz5Momentum & p, 
           const Energy m1, const Energy m2,
           const Axis & unitDir1,
           Lorentz5Momentum & p1, Lorentz5Momentum & p2) {
      Energy min=p.mass();
      if ( min >= m1 + m2  &&  m1 >= ZERO  &&  m2 >= ZERO  ) {
		Momentum3 pstarVector = unitDir1 * Kinematics::pstarTwoBodyDecay(min,m1,m2);
  p1 = Lorentz5Momentum(m1, pstarVector);
  p2 = Lorentz5Momentum(m2,-pstarVector);
  // boost from CM to LAB
  Boost bv = p.boostVector();
  double gammarest = p.e()/p.mass();
  p1.boost( bv, gammarest );
  p2.boost( bv, gammarest );
  return true;
      }
      return false;
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
  const unsigned int MAXTRY = 100;
  unsigned int ntry=0;
  do {
    // ff is the mass squared of the 23 subsystem
    ff = UseRandom::rnd()*(cc-bb)+bb;
    
    // pp is ((m0+m1)^2 - m23^2)((m0-m1)^2-m23)
    pp = (aa-ff)*(bb-ff);
    
    // qq is ((m2+m3)^2 - m23^2)(|m2-m3|^2-m23^2)
    qq = (cc-ff)*(dd-ff);
    
    // weight
    ww = (fcn != NULL) ? (*fcn)(ff,a1,b1,c1) : 1.0;
    ww = sqr(ww);
    rr = ee*ff*UseRandom::rnd();
    ++ntry;
  } 
  while(pp*qq*ww < rr*rr && ntry < MAXTRY );
  if(ntry >= MAXTRY) {
    CurrentGenerator::log() << "Kinematics::threeBodyDecay can't generate momenta" 
			    << " after " << MAXTRY << " attempts\n";
    return false;
  }

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
