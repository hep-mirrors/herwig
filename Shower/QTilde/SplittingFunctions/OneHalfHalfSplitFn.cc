// -*- C++ -*-
//
// OneHalfHalfSplitFn.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the OneHalfHalfSplitFn class.
//

#include "OneHalfHalfSplitFn.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"
#include "Herwig/Shower/QTilde/QTildeShowerHandler.h"

using namespace Herwig;

DescribeNoPIOClass<OneHalfHalfSplitFn,Herwig::Sudakov1to2FormFactor>
describeOneHalfHalfSplitFn ("Herwig::OneHalfHalfSplitFn","HwShower.so");

void OneHalfHalfSplitFn::Init() {

  static ClassDocumentation<OneHalfHalfSplitFn> documentation
    ("The OneHalfHalfSplitFn class implements the splitting function for g->q qbar");

}

double OneHalfHalfSplitFn::integOverP(const double z, const IdList &,
				      unsigned int PDFfactor) const { 
  switch(PDFfactor) {
  case 0:
    return z; 
  case 1:
    return log(z);
  case 2:
    return -log(1.-z);
  case 3:
    return log(z/(1.-z));
  case 4:
    return 2.*sqrt(z);
  case 5:
    return 2./3.*z*sqrt(z);
  default:
    throw Exception() << "OneHalfHalfSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}

double OneHalfHalfSplitFn::invIntegOverP(const double r,
					 const IdList &,
					 unsigned int PDFfactor) const {
  switch(PDFfactor) {
  case 0:
    return r; 
  case 1:
    return exp(r);
  case 2:
    return 1.-exp(-r);
  case 3:
    return 1./(1.+exp(-r));
  case 4:
    return 0.25*sqr(r);
  case 5:
    return pow(1.5*r,2./3.);
  default:
    throw Exception() << "OneHalfHalfSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}

vector<pair<int, Complex> > 
OneHalfHalfSplitFn::generatePhiForward(const double z, const Energy2 t, const IdList & ids,
				       const RhoDMatrix & rho) { 
  assert(rho.iSpin()==PDT::Spin1);
  double modRho = abs(rho(0,2));
  Energy mq = ids[1]->mass();
  Energy2 mq2 = sqr(mq);
  double fact = z*(1.-z)-mq2/t;
  double max = 1.+2.*fact*(-1.+2.*modRho);
  vector<pair<int, Complex> > output;
  output.reserve(3);
  output.push_back(make_pair( 0,(rho(0,0)+rho(2,2))*(1.-2.*fact)/max));
  output.push_back(make_pair(-2,2.*fact*rho(0,2)/max));
  output.push_back(make_pair( 2,2.*fact*rho(2,0)/max));
  return output;
}

DecayMEPtr OneHalfHalfSplitFn::matrixElement(const double z, const Energy2 t, 
					     const IdList & ids, const double phi,
                                             bool timeLike) {
  static const Complex ii(0.,1.);
  // calculate the kernal
  DecayMEPtr kernal(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1,PDT::Spin1Half,PDT::Spin1Half)));
  double mt = !timeLike ? ZERO : ids[1]->mass()/sqrt(t);
  double root =sqrt(1.-sqr(mt)/z/(1.-z));
  (*kernal)(0,0,0) = mt/sqrt(z*(1.-z));
  (*kernal)(2,1,1) = (*kernal)(0,0,0);
  (*kernal)(0,0,1) = -z*root*exp(-ii*phi);
  (*kernal)(2,1,0) = -conj((*kernal)(0,0,1));
  (*kernal)(0,1,0) = (1.-z)*exp(-ii*phi)*root;
  (*kernal)(2,0,1) = -conj((*kernal)(0,1,0));
  (*kernal)(0,1,1) = 0.;
  (*kernal)(2,0,0) = 0.;
  return kernal;
}

void OneHalfHalfSplitFn::findZero() {
  const double xTest = 0.2;
  vector<Energy2> m2;
  for(int id=4;id<=6;++id) {
    tPDPtr parton = getParticleData(id);
    Energy2 sMin(0.5*GeV2),sMax(4.*sqr(parton->mass()));
    // ensure min and max in correct regions
    while (PDF()->xfx(beam(),parton,sMin,xTest)!=0. && sMin>0.01*GeV2) {
      sMin *=0.5;
    }
    while (PDF()->xfx(beam(),parton,sMax,xTest)==0. && sMax<1e10*GeV2) {
      sMax *=2.;
    }
    if(sMin<0.01*GeV2 || sMax>1e10*GeV2) {
      m2.push_back(ZERO);
      continue;
    }
    // bisect to find zero
    while(sMax-sMin>1e-10*GeV2) {
      Energy2 sMid = 0.5*(sMax+sMin);
      double test = PDF()->xfx(beam(),parton,sMid,xTest);
      if(test!=0.) sMax = sMid;
      else         sMin = sMid;
    }
    m2.push_back(sMin);
  }
  mq_[PDF()] = m2;
}

void OneHalfHalfSplitFn::guesstz(Energy2 t1,unsigned int iopt,
                                 const IdList &ids,
                                 double enhance,bool ident,
                                 double detune,
                                 Energy2 &t_main, double &z_main) {
  mq2_ = ZERO;
  // get the mass to fix
  if(iopt==1 && abs(ids[1]->id())>3) {
    map<cPDFPtr,vector<Energy2>>::const_iterator cit = mq_.find(PDF());
    if(cit== mq_.end()) findZero();
    cit = mq_.find(PDF());
    mq2_ = cit->second[abs(ids[1]->id())-4];
    if(t1-mq2_<1e-10*GeV2) {
      t_main = ZERO;
      return;
    }
  }
  unsigned int pdfopt = iopt!=1 ? 0 : pdfFactor();
  double lower = integOverP(zLimits().first ,ids,pdfopt);
  double upper = integOverP(zLimits().second,ids,pdfopt);
  double c = 1./((upper - lower) * colourFactor()
                 * alpha()->showerOverestimateValue()/Constants::twopi*enhance*detune);
  double r = UseRandom::rnd();
  assert(iopt<=2);
  if(iopt==1) {
    c/=pdfMax();
    //symmetry of FS gluon splitting
    if(ident) c*=0.5;
  }
  else if(iopt==2) c*=-1.;
  // guessing t
  if(mq2_!=ZERO) {
    t_main = (t1-mq2_)*pow(r,c) + mq2_;
  }
  else {
    if(iopt!=2 || c*log(r)<log(Constants::MaxEnergy2/t1)) {
      t_main = t1*pow(r,c);
    }
    else
      t_main = Constants::MaxEnergy2;
  }
  // guessing z
  z_main = invIntegOverP(lower + UseRandom::rnd()
                         *(upper - lower),ids,pdfopt);
}

double OneHalfHalfSplitFn::PDFVetoRatio(const Energy2 t, const double x, const double z,
                                        const tcPDPtr parton0, const tcPDPtr parton1,double factor) const {
  assert(PDF());
  Energy2 theScale = t * sqr(ShowerHandler::currentHandler()->factorizationScaleFactor()*factor);
  if (theScale < sqr(freeze())) theScale = sqr(freeze());
  
  const double newpdf=PDF()->xfx(beam(),parton0,theScale,x/z);
  if(newpdf<=0.) return 0.;

  const double oldpdf=PDF()->xfx(beam(),parton1,theScale,x);
  if(oldpdf<=0.) return 1.;
  
  double ratio = newpdf/oldpdf;
  double maxpdf = pdfMax();
  if(mq2_!=ZERO) {
    ratio *= (t-mq2_)/t;
  }
  switch (pdfFactor()) {
  case 0: break;
  case 1: maxpdf /= z; break;
  case 2: maxpdf /= 1.-z; break;
  case 3: maxpdf /= (z*(1.-z)); break;
  case 4: maxpdf /= sqrt(z); break;
  case 5: maxpdf *= sqrt(z); break;
  default :
    throw Exception() << "OneHalfHalfSplitFn::PDFVetoRatio invalid PDFfactor = "
		      << pdfFactor() << Exception::runerror;
    
  }
  
  if (ratio > maxpdf && (mq2_==ZERO || (mq2_!=ZERO && (t-mq2_)/t >1e-6  ) ) ) 
    generator()->log() << "OneHalfHalfSplitFn::PDFVetoRatio PDFVeto warning: Ratio > " << name()
                       << ":PDFmax (by a factor of " << ratio/maxpdf <<") for "
                       << parton0->PDGName() << " to " << parton1->PDGName() << " " << (t-mq2_)/t << "\n";
  return ratio/maxpdf ;
}
