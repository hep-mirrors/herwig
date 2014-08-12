// -*- C++ -*-
//
// HalfHalfOneSplitFn.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HalfHalfOneSplitFn class.
//

#include "HalfHalfOneSplitFn.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig++/Shower/Base/ShowerParticle.h"

using namespace Herwig;

DescribeNoPIOClass<HalfHalfOneSplitFn,Herwig::SplittingFunction>
describeHalfHalfOneSplitFn ("Herwig::HalfHalfOneSplitFn","HwShower.so");

void HalfHalfOneSplitFn::Init() {

  static ClassDocumentation<HalfHalfOneSplitFn> documentation
    ("The HalfHalfOneSplitFn class implements the q -> qg splitting function");

}

double HalfHalfOneSplitFn::P(const double z, const Energy2 t,
			     const IdList &ids, const bool mass) const {
  double val = (1. + sqr(z))/(1.-z);
  if(mass) {
    Energy m = getParticleData(ids[0])->mass();  
    val -= 2.*sqr(m)/t;
  }
  return colourFactor(ids)*val;
}

double HalfHalfOneSplitFn::overestimateP(const double z,
					 const IdList & ids) const { 
  return 2.*colourFactor(ids)/(1.-z); 
}

double HalfHalfOneSplitFn::ratioP(const double z, const Energy2 t,
				  const IdList & ids, const bool mass) const {
  double val = 1. + sqr(z);
  if(mass) {
    Energy m = getParticleData(ids[0])->mass();
    val -= 2.*sqr(m)*(1.-z)/t;
  } 
  return 0.5*val;
} 

double HalfHalfOneSplitFn::integOverP(const double z,
				      const IdList & ids,
				      unsigned int PDFfactor) const {
  switch (PDFfactor) {
  case 0:
    return -2.*colourFactor(ids)*Math::log1m(z);
  case 1:
    return  2.*colourFactor(ids)*log(z/(1.-z));
  case 2:
    return  2.*colourFactor(ids)/(1.-z);
  case 3:
  default:
    throw Exception() << "HalfHalfOneSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  } 
}

double HalfHalfOneSplitFn::invIntegOverP(const double r, const IdList & ids,
					 unsigned int PDFfactor) const {
  switch (PDFfactor) {
  case 0:
    return 1. - exp(- 0.5*r/colourFactor(ids)); 
  case 1:
    return 1./(1.-exp(-0.5*r/colourFactor(ids)));
  case 2:
    return 1.-2.*colourFactor(ids)/r;
  case 3:
  default:
    throw Exception() << "HalfHalfOneSplitFn::invIntegOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  } 
}

bool HalfHalfOneSplitFn::accept(const IdList &ids) const {
  // 3 particles and in and out fermion same
  if(ids.size()!=3 || ids[0]!=ids[1]) return false;
  tcPDPtr q=getParticleData(ids[0]);
  tcPDPtr g=getParticleData(ids[2]);
  if(q->iSpin()!=PDT::Spin1Half ||
     g->iSpin()!=PDT::Spin1) return false;
  return checkColours(ids);
}

vector<pair<int, Complex> > 
HalfHalfOneSplitFn::generatePhi(ShowerParticle & ,ShoKinPtr ,
				const double z, const Energy2, const IdList & ,
				const RhoDMatrix &) {
  // no dependence on the spin density matrix, dependence on off-diagonal terms cancels
  // and rest = splitting function for Tr(rho)=1 as required by defn
  return vector<pair<int, Complex> >(1,make_pair(0,1.));
}

DecayMatrixElement HalfHalfOneSplitFn::matrixElement(ShowerParticle & particle,ShoKinPtr showerkin,
						     const double z, const Energy2 t, 
						     const IdList & ids, const double phi) {
  // calculate the kernal
  DecayMatrixElement kernal(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin1);
  Energy m = particle.dataPtr()->mass();
  double mt = m/sqrt(t);
  double root = sqrt(1.-(1.-z)*sqr(m)/z/t);
  double romz = sqrt(1.-z); 
  double rz   = sqrt(z);
  kernal(0,0,0) = -root/romz;
  kernal(0,0,2) =  root/romz*z;
  kernal(0,1,0) =  mt*(1.-z)/rz;
  kernal(0,1,2) =  0.;
  kernal(1,0,0) =  0.;
  kernal(1,0,2) =  mt*(1.-z)/rz;
  kernal(1,1,0) = -root/romz*z;
  kernal(1,1,2) =  root/romz;
  for(unsigned int lamA=0;lamA<2;++lamA) {
    Complex factA = exp(Complex(0.,1.)*double(2*lamA-1)/2.*phi);
    for(unsigned int lamB=0;lamB<2;++lamB) {
      Complex factB = exp(-Complex(0.,1.)*double(2*lamB-1)/2.*phi);
      for(unsigned int lamC=0;lamC<3;lamC+=2) {
	Complex factC = exp(-Complex(0.,1.)*double(lamC-1)*phi);
	kernal(lamA,lamB,lamC/2) *= factA*factB*factC;
      }
    }
  }
  return kernal;
}
