// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ISGWFormFactor class.
//

#include "ISGWFormFactor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ISGWFormFactor.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;

ISGWFormFactor::~ISGWFormFactor() {}

void ISGWFormFactor::persistentOutput(PersistentOStream & os) const {
  os << _kappa << _mdown << _mup << _mstrange << _mcharm << _mbottom << _mquark
     << _betaSud << _betaSus << _betaSuc << _betaSub << _betaS
     << _betaPud << _betaPus << _betaPuc << _betaP;
}

void ISGWFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> _kappa >> _mdown >> _mup >> _mstrange >> _mcharm >> _mbottom >> _mquark
     >> _betaSud >> _betaSus >> _betaSuc >> _betaSub >> _betaS
     >> _betaPud >> _betaPus >> _betaPuc >> _betaP;
}

ClassDescription<ISGWFormFactor> ISGWFormFactor::initISGWFormFactor;
// Definition of the static class description member.

void ISGWFormFactor::Init() {

  static ClassDocumentation<ISGWFormFactor> documentation
    ("The \\classname{ISGWFormFactor} class implements the ISGW model of"
     "Phys. Rev. D39, 799 (1989) for the scalar meson form-factors.");

  static Parameter<ISGWFormFactor,double> interfaceKappa
    ("Kappa",
     "The relavistic compensation factor of the ISGW model",
     &ISGWFormFactor::_kappa, 0.7, 0.0, 5.0,
     false, false, true);

  static Parameter<ISGWFormFactor,Energy> interfaceDownMass
    ("DownMass",
     "The mass of the down quark in the ISGW model (this is a consituent mass)",
     &ISGWFormFactor::_mdown, GeV, 0.33*GeV, 0.0*GeV, 1.0*GeV,
     false, false, true);

  static Parameter<ISGWFormFactor,Energy> interfaceUpMass
    ("UpMass",
     "The mass of the up quark in the ISGW model (this is a consituent mass)",
     &ISGWFormFactor::_mup, GeV, 0.33*GeV, 0.0*GeV, 1.0*GeV,
     false, false, true);

  static Parameter<ISGWFormFactor,Energy> interfaceStrangeMass
    ("StrangeMass",
     "The mass of the strange quark in the ISGW model (this is a consituent mass)",
     &ISGWFormFactor::_mstrange, GeV, 0.55*GeV, 0.0*GeV, 1.0*GeV,
     false, false, true);

  static Parameter<ISGWFormFactor,Energy> interfaceCharmMass
    ("CharmMass",
     "The mass of the charm quark in the ISGW model (this is a consituent mass)",
     &ISGWFormFactor::_mcharm, GeV, 1.82*GeV, 0.0*GeV, 3.0*GeV,
     false, false, true);

  static Parameter<ISGWFormFactor,Energy> interfaceBottomMass
    ("BottomMass",
     "The mass of the bottom quark in the ISGW model (this is a consituent mass)",
     &ISGWFormFactor::_mbottom, GeV, 5.12*GeV, 3.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ISGWFormFactor,Energy> interfaceBetaSud
    ("BetaSud",
     "The variational parameter for s-wave ud mesons",
     &ISGWFormFactor::_betaSud, GeV, 0.31*GeV, 0.0*GeV, 1.0*GeV,
     false, false, true);

  static Parameter<ISGWFormFactor,Energy> interfaceBetaSus
    ("BetaSus",
     "The variational parameter for s-wave us mesons",
     &ISGWFormFactor::_betaSus, GeV, 0.34*GeV, 0.0*GeV, 1.0*GeV,
     false, false, true);

  static Parameter<ISGWFormFactor,Energy> interfaceBetaSuc
    ("BetaSuc",
     "The variational parameter for s-wave uc mesons",
     &ISGWFormFactor::_betaSuc, GeV, 0.39*GeV, 0.0*GeV, 1.0*GeV,
     false, false, true);

  static Parameter<ISGWFormFactor,Energy> interfaceBetaSub
    ("BetaSub",
     "The variational parameter for s-wave ub mesons",
     &ISGWFormFactor::_betaSub, GeV, 0.41*GeV, 0.0*GeV, 1.0*GeV,
     false, false, true);

  static Parameter<ISGWFormFactor,Energy> interfaceBetaPud
    ("BetaPud",
     "The variational parameter for p-wave ud mesons",
     &ISGWFormFactor::_betaPud, GeV, 0.27*GeV, 0.0*GeV, 1.0*GeV,
     false, false, true);

  static Parameter<ISGWFormFactor,Energy> interfaceBetaPus
    ("BetaPus",
     "The variational parameter for p-wave us mesons",
     &ISGWFormFactor::_betaPus, GeV, 0.30*GeV, 0.0*GeV, 1.0*GeV,
     false, false, true);

  static Parameter<ISGWFormFactor,Energy> interfaceBetaPuc
    ("BetaPuc",
     "The variational parameter for p-wave uc mesons",
     &ISGWFormFactor::_betaPuc, GeV, 0.34*GeV, 0.0*GeV, 1.0*GeV,
     false, false, true);
}


// form-factor for scalar to scalar
void ISGWFormFactor::ScalarScalarFormFactor(Energy2 q2, int iloc,int id0, int id1,
					     Energy mY, Energy mX,
					     Complex & f0,Complex & fp) const
 {
   Complex d1(0.),d2(0.);
   formFactor(q2,iloc,id0,id1,mY,mX,f0,fp,d1,d2);
 }

// form-factor for scalar to vector
void ISGWFormFactor::ScalarVectorFormFactor(Energy2 q2, int iloc, int id0, int id1, 
					     Energy mY, Energy mX,
					     Complex & A0,Complex & A1,
					     Complex & A2,Complex & V) const
 {formFactor(q2,iloc,id0,id1,mY,mX,A0,A1,A2,V);}


// form-factor for scalar to tensor
void ISGWFormFactor::ScalarTensorFormFactor(Energy2 q2, int iloc, int id0, int id1, 
					     Energy mY, Energy mX,
					     Complex & h,Complex & k,
					     Complex & bp,Complex & bm) const
 {formFactor(q2,iloc,id0,id1,mY,mX,h,k,bp,bm);}

// member which does the work
void ISGWFormFactor::formFactor(Energy2 q2, int iloc, int id0, int id1, Energy mY,
				 Energy mX, Complex & f1,Complex & f2,Complex & f3,
				 Complex & f4) const
{
  // work out the flavour of the heavy quarks
  int ifl0=(abs(id0)/100)%10;
  int ifls=(abs(id0)/10)%10;
  int ifl1=(abs(id1)/100)%10;
  int ispin=(abs(id1)/1000);
  int jspin=abs(id1)%10;
  // masses of the quarks
  Energy mQ=_mquark[ifl0-1];
  Energy mq=_mquark[ifl1-1];
  Energy ms=_mquark[ifls-1];
  // of the mesons
  Energy mtildeX=mq+ms;
  Energy mtildeY=mQ+ms;
  // wavefunction parameters for the mesons
  Energy betaX,betaY;
  // spin-0 outgoing mesons
  if(jspin==1)
    {
      // the wavefunction parameter for the incoming meson
      if(ifl0>=2){betaY=_betaS[ifl0-2];}
      else {betaY=_betaS[0];}
      // the wavefunction parameter for the outgoing meson
      if(ispin==0||ispin>50)
	{
	  if(ifl1>=2){betaX=_betaS[ifl1-2];}
	  else {betaX=_betaS[0];}
	}
      else
	{
	  if(ifl1>=2){betaX=_betaP[ifl1-2];}
	  else {betaX=_betaP[0];}
	}
    }
  // spin-1 outgoing mesons
  else if(jspin==3)
    {
      // the wavefunction parameter for the incoming meson
      if(ifl0>=2){betaY=_betaS[ifl0-2];}
      else {betaY=_betaS[0];}
      // the wavefunction parameter for the outgoing meson
      if(ispin==0||ispin>50)
	{
	  if(ifl1>=2){betaX=_betaS[ifl1-2];}
	  else {betaX=_betaS[0];}
	}
      else
	{
	  if(ifl1>=2){betaX=_betaP[ifl1-2];}
	  else {betaX=_betaP[0];}
	}
    }
  // spin-2 outgoing mesons
  else if(jspin==5)
    {
      // the wavefunction parameter for the incoming meson
      if(ifl0>=2){betaY=_betaS[ifl0-2];}
      else {betaY=_betaS[0];}
      // the wavefunction parameter for the outgoing meson
      if(ifl1>=2){betaX=_betaP[ifl1-2];}
      else {betaX=_betaP[0];}
    }
  else
    {throw Exception() << "ISGWFormFactor::formFactor" 
			<< " unknown spin of outgoing meson." << Exception::abortnow;}
  // compute the F_n function we will need
  Energy2 beta2XY = 0.5*(betaX*betaX+betaY*betaY);
  Energy2 tm=(mY-mX)*(mY-mX);
  double betar = betaX*betaY/beta2XY;
  double kappa2=_kappa*_kappa;
  double slope=(tm-q2)/(kappa2*beta2XY);
  Energy mup=mq*mQ/(mQ+mq);
  Energy mum=mq*mQ/(mQ-mq);
  double fn = sqrt(mtildeX/mtildeY)*betar*sqrt(betar)*
    exp(-0.25*ms*ms/(mtildeX*mtildeY)*slope);
  // now we can compute the form-factors
  // for scalars
  if(jspin==1)
    {
      Complex fp,fm;
      // 1 1S0
      if(ispin==0)
	{
	  double yratio=betaY*betaY/beta2XY;
	  fp = fn*(1.+0.5*mQ/mum-0.25*mQ*mq/mup/mum*ms/mtildeX*yratio);
	  fm = fn*(1.-(mtildeX+mtildeY)*(0.5/mq-0.25/mup*ms/mtildeX*yratio));
	}
      // 1 3P0
      else if(ispin<100)
	{
	  // extra power of beta factors
	  fn*=betar;
	  fp = fn*ms*mQ*mq/sqrt(6.)/betaY/mtildeX/mum;
	  fm =-fn*ms/betaY/sqrt(6.)*(mtildeX+mtildeY)/mtildeX;
	}
      // 2 1S0
      else
	{throw Exception() << "ISGWFormFactor::formFactor" 
			   << " 2S not implemented" << Exception::abortnow;}
      // convert to the standard form
      f1 = q2/(mY*mY-mX*mX)*fm+fp;
      f2 = fp;
    }
  // for vectors
  else if(jspin==3)
    {
      Complex f,g,ap,am;
      Energy2 betaX2=betaX*betaX;
      Energy2 betaY2=betaY*betaY;
      //  1 3S1
      if(ispin==0)
	{
	  f  = 2.*mtildeY*fn;
	  g  = 0.5*fn*(1./mq-0.5/mum*ms/mtildeX*betaY2/beta2XY);
	  ap =-0.5*fn/mtildeX*(1.+ms/mQ*(betaY2-betaX2)/(betaX2+betaY2)
			       -0.25*ms*ms/mum/mtildeY*betaX2*betaX2/beta2XY/beta2XY);
	  am = 0.5*fn/mtildeX*(1.+ms/mQ+ms*ms/mq/mQ*betaX2/beta2XY*
			       (1.-0.25*(mtildeX+mtildeY)/mtildeY*betaX2/beta2XY));
	}
      // 1 3P1
      else if(ispin==20)
	{
	  fn*=betar;
	  f  =-fn*mtildeY*betaY*(1./mum+0.5*ms/mtildeY*slope*
				 (1./mq-0.5/mum*ms/mtildeX*betaY2/beta2XY));
	  g  = 0.5*fn*ms/mtildeX/betaY;
	  ap = 0.25*fn*ms*mQ/mtildeY/betaY/mum*(1.-0.5*ms*mq/mtildeX/mum*betaY2/beta2XY);
	  am = -0.25*fn*ms*(mtildeX+mtildeY)/mq/betaY/mtildeY*
	    (1.-0.5*ms*mq/mtildeX/mum*betaY2/beta2XY);
	}
      //  1 1P1
      else if(ispin==10)
	{
	  fn*=betar;
	  double ort=1./sqrt(2.);
	  f  = fn*ort*mtildeY*betaY/mup;
	  g  = 0.25*fn*mtildeY*betaY*ort/mq/mQ/mtildeX;
	  Energy msum=mtildeX+mtildeY;
	  ap = fn*ms*ort/betaY/mtildeY*(1.+0.5*mQ/mum
					-0.25*mq*mQ*ms/mum/mup/mtildeX*betaY2/beta2XY);
	  am = fn*ms*ort/betaY/mtildeY*(1.-0.5/mq*msum
					+0.25*ms*betaY2/mup/beta2XY*msum/mtildeX);
	}
      // 2 1S0
      else
	{throw Exception() << "ISGWScalarVectorDecayer::formFactor" 
			    << " 2S not implemented" << Exception::abortnow;}
      // convert to the standard notation
      Energy msum=mX+mY,mdiff=mY-mX;
      Complex ii(0.,1.);
      f2 = -ii*f/msum;
      f3 = +ii*ap*msum;
      f1 = -ii*0.5/mX*(am*q2-ii*msum*f2+ii*mdiff*f3);
      f4 = -ii*g*msum;
    }
  // for tensors
  else if(jspin==5)
    {
      fn *=fn/sqrt(2.);
      double betaXb2=betaX*betaX/beta2XY;
      // 1 3P2
      if(ispin==0)
	{
	  f1 = 0.5*fn*ms/mtildeY/betaY*(1./mq
				       -0.5*ms/mtildeX/mum*betaY*betaY/beta2XY);
	  f2 = sqrt(2.)*fn*ms/betaY;
	  f3 =-0.5*fn*ms/mtildeX/mQ/betaY*
	    (1.-0.5*ms*mQ/mup/mtildeY*betaXb2
	     +0.25*ms*mQ/mtildeY/mum*betaXb2*(1.-0.5*ms*betaXb2/mtildeY));
	  f4 = 0.5*fn*ms/mtildeX/mQ/betaY*
	    (1.-0.5*ms*mQ/mup/mtildeY*betaXb2+
	     0.25*ms*betaXb2/mq*(mtildeX+mtildeY)/mtildeY*(1.-0.5*ms*betaXb2/mtildeY));
	}
    }
}
}
