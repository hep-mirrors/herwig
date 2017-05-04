// -*- C++ -*-
//
// ISGWFormFactor.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ISGWFormFactor class.
//

#include "ISGWFormFactor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

ISGWFormFactor::ISGWFormFactor() {
  // default values of the parameters
  // fudge factor
  _kappa=0.7;
  // quark masses
  _mdown   =0.33*GeV;
  _mup     =0.33*GeV;
  _mstrange=0.55*GeV;
  _mcharm  =1.82*GeV;
  _mbottom =5.12*GeV;
  // beta values
  _betaSud = 0.31*GeV;
  _betaSus = 0.34*GeV;
  _betaSuc = 0.39*GeV;
  _betaSub = 0.41*GeV;
  _betaPud = 0.27*GeV;
  _betaPus = 0.30*GeV;
  _betaPuc = 0.34*GeV;
  // the mixing for eta/eta'
  _thetaeta=-Constants::pi/18.;
  // B_u decays to d ubar
  addFormFactor(-521,-211  ,0,-2, 5, 1);
  addFormFactor(-521,-213  ,1,-2, 5, 1);
  addFormFactor(-521,-215  ,2,-2, 5, 1);
  addFormFactor(-521,-10213,1,-2, 5, 1);
  addFormFactor(-521,-20213,1,-2, 5, 1);
  addFormFactor(-521,-10211,0,-2, 5, 1);
  // B_u to uu (I=0)
  addFormFactor(-521, 221  ,0,-2, 5, 2);
  addFormFactor(-521, 331  ,0,-2, 5, 2);
  addFormFactor(-521, 223  ,1,-2, 5, 2);
  addFormFactor(-521, 225  ,2,-2, 5, 2);
  addFormFactor(-521, 10223,1,-2, 5, 2);
  addFormFactor(-521, 20223,1,-2, 5, 2);
  addFormFactor(-521, 10221,0,-2, 5, 2);
  // B_u to uu (I=1)
  addFormFactor(-521, 111  ,0,-2, 5, 2);
  addFormFactor(-521, 113  ,1,-2, 5, 2);
  addFormFactor(-521, 115  ,2,-2, 5, 2);
  addFormFactor(-521, 10113,1,-2, 5, 2);
  addFormFactor(-521, 20113,1,-2, 5, 2);
  addFormFactor(-521, 10111,0,-2, 5, 2);
  // B_u decays to d ubar
  addFormFactor(-521,-321  ,0,-2, 5, 3);
  addFormFactor(-521,-323  ,1,-2, 5, 3);
  addFormFactor(-521,-325  ,2,-2, 5, 3);
  addFormFactor(-521,-10323,1,-2, 5, 3);
  addFormFactor(-521,-20323,1,-2, 5, 3);
  addFormFactor(-521,-10321,0,-2, 5, 3);
  // B_u decays to c ubar
  addFormFactor(-521, 421  ,0,-2, 5, 4);
  addFormFactor(-521, 423  ,1,-2, 5, 4);
  addFormFactor(-521, 425  ,2,-2, 5, 4);
  addFormFactor(-521, 10423,1,-2, 5, 4);
  addFormFactor(-521, 20423,1,-2, 5, 4);
  addFormFactor(-521, 10421,0,-2, 5, 4);
  // B_d to d dbar (I=0)
  addFormFactor(-511, 221  ,0, 1,-5,-1);
  addFormFactor(-511, 331  ,0, 1,-5,-1);
  addFormFactor(-511, 223  ,1, 1,-5,-1);
  addFormFactor(-511, 225  ,2, 1,-5,-1);
  addFormFactor(-511, 10223,1, 1,-5,-1);
  addFormFactor(-511, 20223,1, 1,-5,-1);
  addFormFactor(-511, 10221,0, 1,-5,-1);
  // B_d to d dbar (I=1)
  addFormFactor(-511, 111  ,0, 1,-5,-1);
  addFormFactor(-511, 113  ,1, 1,-5,-1);
  addFormFactor(-511, 115  ,2, 1,-5,-1);
  addFormFactor(-511, 10113,1, 1,-5,-1);
  addFormFactor(-511, 20113,1, 1,-5,-1);
  addFormFactor(-511, 10111,0, 1,-5,-1);
  // B_d to u dbar
  addFormFactor(-511, 211  ,0, 1,-5,-2);
  addFormFactor(-511, 213  ,1, 1,-5,-2);
  addFormFactor(-511, 215  ,2, 1,-5,-2);
  addFormFactor(-511, 10213,1, 1,-5,-2);
  addFormFactor(-511, 20213,1, 1,-5,-2);
  addFormFactor(-511, 10211,0, 1,-5,-2);
  // B_d to s dbar 
  addFormFactor(-511, 311  ,0, 1,-5,-3);
  addFormFactor(-511, 313  ,1, 1,-5,-3);
  addFormFactor(-511, 315  ,2, 1,-5,-3);
  addFormFactor(-511, 10313,1, 1,-5,-3);
  addFormFactor(-511, 20313,1, 1,-5,-3);
  addFormFactor(-511, 10311,0, 1,-5,-3);
  // B_d decays to  c dbar
  addFormFactor(-511, 411  ,0, 1,-5,-4);
  addFormFactor(-511, 413  ,1, 1,-5,-4);
  addFormFactor(-511, 415  ,2, 1,-5,-4);
  addFormFactor(-511, 10413,1, 1,-5,-4);
  addFormFactor(-511, 20413,1, 1,-5,-4);
  addFormFactor(-511, 10411,0, 1,-5,-4);
  // D0 to d ubar
  addFormFactor( 421,-211  ,0,-2, 4, 1);
  addFormFactor( 421,-213  ,1,-2, 4, 1);
  addFormFactor( 421,-215  ,2,-2, 4, 1);
  addFormFactor( 421,-10213,1,-2, 4, 1);
  addFormFactor( 421,-20213,1,-2, 4, 1);
  addFormFactor( 421,-10211,0,-2, 4, 1);
  // D0 to d ubar (I=1)
  addFormFactor( 421, 111  ,0,-2, 4, 2);
  addFormFactor( 421, 113  ,1,-2, 4, 2);
  addFormFactor( 421, 115  ,2,-2, 4, 2);
  addFormFactor( 421, 10113,1,-2, 4, 2);
  addFormFactor( 421, 20113,1,-2, 4, 2);
  addFormFactor( 421, 10111,0,-2, 4, 2);
  // D0 to d ubar (I=0)
  addFormFactor( 421, 221  ,0,-2, 4, 2);
  addFormFactor( 421, 331  ,0,-2, 4, 2);
  addFormFactor( 421, 223  ,1,-2, 4, 2);
  addFormFactor( 421, 225  ,2,-2, 4, 2);
  addFormFactor( 421, 10223,1,-2, 4, 2);
  addFormFactor( 421, 20223,1,-2, 4, 2);
  addFormFactor( 421, 10221,0,-2, 4, 2);
  // D0 to s ubar
  addFormFactor( 421,-321  ,0,-2, 4, 3);
  addFormFactor( 421,-323  ,1,-2, 4, 3);
  addFormFactor( 421,-325  ,2,-2, 4, 3);
  addFormFactor( 421,-10323,1,-2, 4, 3);
  addFormFactor( 421,-20323,1,-2, 4, 3);
  addFormFactor( 421,-10321,0,-2, 4, 3);
  // D+ to d dbar I=0
  addFormFactor( 411, 221  ,0,-1, 4, 1);
  addFormFactor( 411, 331  ,0,-1, 4, 1);
  addFormFactor( 411, 223  ,1,-1, 4, 1);
  addFormFactor( 411, 225  ,2,-1, 4, 1); 
  addFormFactor( 411, 10223,1,-1, 4, 1); 
  addFormFactor( 411, 20223,1,-1, 4, 1); 
  addFormFactor( 411, 10221,0,-1, 4, 1);
  // D+ to d dbar I=1
  addFormFactor( 411, 111  ,0,-1, 4, 1);
  addFormFactor( 411, 113  ,1,-1, 4, 1); 
  addFormFactor( 411, 115  ,2,-1, 4, 1); 
  addFormFactor( 411, 10113,1,-1, 4, 1); 
  addFormFactor( 411, 20113,1,-1, 4, 1); 
  addFormFactor( 411, 10111,0,-1, 4, 1);
  // D+ to u dbar
  addFormFactor( 411, 211  ,0,-1, 4, 2);
  addFormFactor( 411, 213  ,1,-1, 4, 2);
  addFormFactor( 411, 215  ,2,-1, 4, 2);
  addFormFactor( 411, 10213,1,-1, 4, 2);
  addFormFactor( 411, 20213,1,-1, 4, 2);
  addFormFactor( 411, 10211,0,-1, 4, 2);
  // D+ to s dbar
  addFormFactor( 411,-311  ,0,-1, 4, 3);
  addFormFactor( 411,-313  ,1,-1, 4, 3);
  addFormFactor( 411,-315  ,2,-1, 4, 3);
  addFormFactor( 411,-10313,1,-1, 4, 3);
  addFormFactor( 411,-20313,1,-1, 4, 3);
  addFormFactor( 411,-10311,0,-1, 4, 3);
  // set the initial number of modes
  initialModes(numberOfFactors());
}

void ISGWFormFactor::doinit() {
  ScalarFormFactor::doinit();
  // set up the quark masses
  _mquark.resize(5);
  _mquark[0]=_mdown;
  _mquark[1]=_mup;
  _mquark[2]=_mstrange;
  _mquark[3]=_mcharm;
  _mquark[4]=_mbottom;
  // and the beta values
  _betaS.resize(5,vector<Energy>(5));
  _betaP.resize(5,vector<Energy>(5));
  _betaS[0][0] = _betaSud;_betaP[0][0] = _betaPud;
  _betaS[1][0] = _betaSud;_betaP[1][0] = _betaPud;
  _betaS[2][0] = _betaSus;_betaP[2][0] = _betaPus;
  _betaS[3][0] = _betaSuc;_betaP[3][0] = _betaPuc;
  _betaS[4][0] = _betaSub;_betaP[4][0] = ZERO  ;
  _betaS[0][1] = _betaSud;_betaP[0][1] = _betaPud;
  _betaS[1][1] = _betaSud;_betaP[1][1] = _betaPud;
  _betaS[2][1] = _betaSus;_betaP[2][1] = _betaPus;
  _betaS[3][1] = _betaSuc;_betaP[3][1] = _betaPuc;
  _betaS[4][1] = _betaSub;_betaP[4][1] = ZERO  ;
  _betaS[0][2] = _betaSus;_betaP[0][2] = _betaPus;
  _betaS[1][2] = _betaSus;_betaP[1][2] = _betaPus;
  _betaS[2][2] = ZERO  ;_betaP[2][2] = ZERO  ;
  _betaS[3][2] = ZERO  ;_betaP[3][2] = ZERO  ;
  _betaS[4][2] = ZERO  ;_betaP[4][2] = ZERO  ;
  _betaS[0][3] = _betaSuc;_betaP[0][3] = _betaPuc;
  _betaS[1][3] = _betaSuc;_betaP[1][3] = _betaPuc;
  _betaS[2][3] = ZERO  ;_betaP[2][3] = ZERO  ;
  _betaS[3][3] = ZERO  ;_betaP[3][3] = ZERO  ;
  _betaS[4][3] = ZERO  ;_betaP[4][3] = ZERO  ;
  _betaS[0][4] = ZERO  ;_betaP[0][4] = ZERO  ;
  _betaS[1][4] = ZERO  ;_betaP[1][4] = ZERO  ;
  _betaS[2][4] = ZERO  ;_betaP[2][4] = ZERO  ;
  _betaS[3][4] = ZERO  ;_betaP[3][4] = ZERO  ;
  _betaS[4][4] = ZERO  ;_betaP[4][4] = ZERO  ;
}

void ISGWFormFactor::persistentOutput(PersistentOStream & os) const {
  os << _kappa << ounit(_mdown,GeV) << ounit(_mup,GeV) << ounit(_mstrange,GeV) 
     << ounit(_mcharm,GeV) << ounit(_mbottom,GeV) << ounit(_betaSud,GeV) 
     << ounit(_betaSus,GeV) << ounit(_betaSuc,GeV) << ounit(_betaSub,GeV) 
     << ounit(_betaPud,GeV) << ounit(_betaPus,GeV) << ounit(_betaPuc,GeV)
     << _thetaeta << ounit(_mquark,GeV) << ounit(_betaS,GeV) << ounit(_betaP,GeV);
}

void ISGWFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> _kappa >> iunit(_mdown,GeV) >> iunit(_mup,GeV) >> iunit(_mstrange,GeV) 
     >> iunit(_mcharm,GeV) >> iunit(_mbottom,GeV) >> iunit(_betaSud,GeV) 
     >> iunit(_betaSus,GeV) >> iunit(_betaSuc,GeV) >> iunit(_betaSub,GeV) 
     >> iunit(_betaPud,GeV) >> iunit(_betaPus,GeV) >> iunit(_betaPuc,GeV)
     >> _thetaeta >> iunit(_mquark,GeV) >> iunit(_betaS,GeV) >> iunit(_betaP,GeV);
}

ClassDescription<ISGWFormFactor> ISGWFormFactor::initISGWFormFactor;
// Definition of the static class description member.

void ISGWFormFactor::Init() {

  static ClassDocumentation<ISGWFormFactor> documentation
    ("The ISGWFormFactor class implements the ISGW model of"
     "Phys. Rev. D39, 799 (1989) for the scalar meson form-factors.",
     "The form factor model of ISGW \\cite{Isgur:1988gb} together with the"
     "form factors for the term which are supressed by the lepton mass from"
     "\\cite{Scora:1989ys,Isgur:1990jf}",
     "\\bibitem{Isgur:1988gb} N.~Isgur, D.~Scora, B.~Grinstein and M.~B.~Wise,\n"
     "Phys.\\ Rev.\\  D {\\bf 39} (1989) 799.\n"
     "%%CITATION = PHRVA,D39,799;%%\n"
     "\\bibitem{Scora:1989ys} D.~Scora and N.~Isgur, \n"
     "Phys.\\ Rev.\\  D {\\bf 40} (1989) 1491.\n"
     "%%CITATION = PHRVA,D40,1491;%%\n"
     "\\bibitem{Isgur:1990jf} N.~Isgur and M.~B.~Wise,\n"
     "Phys.\\ Rev.\\  D {\\bf 43} (1991) 819.\n"
     "%%CITATION = PHRVA,D43,819;%%\n");

  static Parameter<ISGWFormFactor,double> interfaceKappa
    ("Kappa",
     "The relavistic compensation factor of the ISGW model",
     &ISGWFormFactor::_kappa, 0.7, 0.0, 1.0,
     false, false, true);

  static Parameter<ISGWFormFactor,Energy> interfaceDownMass
    ("DownMass",
     "The mass of the down quark in the ISGW model (this is a consituent mass)",
     &ISGWFormFactor::_mdown, GeV, 0.33*GeV, ZERO, 1.0*GeV,
     false, false, true);

  static Parameter<ISGWFormFactor,Energy> interfaceUpMass
    ("UpMass",
     "The mass of the up quark in the ISGW model (this is a consituent mass)",
     &ISGWFormFactor::_mup, GeV, 0.33*GeV, ZERO, 1.0*GeV,
     false, false, true);

  static Parameter<ISGWFormFactor,Energy> interfaceStrangeMass
    ("StrangeMass",
     "The mass of the strange quark in the ISGW model (this is a consituent mass)",
     &ISGWFormFactor::_mstrange, GeV, 0.55*GeV, ZERO, 1.0*GeV,
     false, false, true);

  static Parameter<ISGWFormFactor,Energy> interfaceCharmMass
    ("CharmMass",
     "The mass of the charm quark in the ISGW model (this is a consituent mass)",
     &ISGWFormFactor::_mcharm, GeV, 1.82*GeV, ZERO, 3.0*GeV,
     false, false, true);

  static Parameter<ISGWFormFactor,Energy> interfaceBottomMass
    ("BottomMass",
     "The mass of the bottom quark in the ISGW model (this is a consituent mass)",
     &ISGWFormFactor::_mbottom, GeV, 5.12*GeV, 3.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ISGWFormFactor,Energy> interfaceBetaSud
    ("BetaSud",
     "The variational parameter for s-wave ud mesons",
     &ISGWFormFactor::_betaSud, GeV, 0.31*GeV, ZERO, 1.0*GeV,
     false, false, true);

  static Parameter<ISGWFormFactor,Energy> interfaceBetaSus
    ("BetaSus",
     "The variational parameter for s-wave us mesons",
     &ISGWFormFactor::_betaSus, GeV, 0.34*GeV, ZERO, 1.0*GeV,
     false, false, true);

  static Parameter<ISGWFormFactor,Energy> interfaceBetaSuc
    ("BetaSuc",
     "The variational parameter for s-wave uc mesons",
     &ISGWFormFactor::_betaSuc, GeV, 0.39*GeV, ZERO, 1.0*GeV,
     false, false, true);

  static Parameter<ISGWFormFactor,Energy> interfaceBetaSub
    ("BetaSub",
     "The variational parameter for s-wave ub mesons",
     &ISGWFormFactor::_betaSub, GeV, 0.41*GeV, ZERO, 1.0*GeV,
     false, false, true);

  static Parameter<ISGWFormFactor,Energy> interfaceBetaPud
    ("BetaPud",
     "The variational parameter for p-wave ud mesons",
     &ISGWFormFactor::_betaPud, GeV, 0.27*GeV, ZERO, 1.0*GeV,
     false, false, true);

  static Parameter<ISGWFormFactor,Energy> interfaceBetaPus
    ("BetaPus",
     "The variational parameter for p-wave us mesons",
     &ISGWFormFactor::_betaPus, GeV, 0.30*GeV, ZERO, 1.0*GeV,
     false, false, true);

  static Parameter<ISGWFormFactor,Energy> interfaceBetaPuc
    ("BetaPuc",
     "The variational parameter for p-wave uc mesons",
     &ISGWFormFactor::_betaPuc, GeV, 0.34*GeV, ZERO, 1.0*GeV,
     false, false, true);

  static Parameter<ISGWFormFactor,double> interfaceThetaEtaEtaPrime
    ("ThetaEtaEtaPrime",
     "The eta-eta' mixing angle",
     &ISGWFormFactor::_thetaeta, -Constants::pi/18., -Constants::pi, Constants::pi,
     false, false, true);

}


// form-factor for scalar to scalar
void ISGWFormFactor::ScalarScalarFormFactor(Energy2 q2, unsigned int iloc,int id0,
					    int id1,Energy mY,Energy mX,
					    Complex & f0,Complex & fp) const {
  Complex d1(0.),d2(0.);
  formFactor(q2,iloc,id0,id1,mY,mX,f0,fp,d1,d2);
}

// form-factor for scalar to vector
void ISGWFormFactor::ScalarVectorFormFactor(Energy2 q2, unsigned int iloc, int id0,
					    int id1,Energy mY, Energy mX,
					    Complex & A0,Complex & A1,
					    Complex & A2,Complex & V) const {
  formFactor(q2,iloc,id0,id1,mY,mX,A0,A1,A2,V);
}


// form-factor for scalar to tensor
void ISGWFormFactor::ScalarTensorFormFactor(Energy2 q2, unsigned int iloc, int id0,
					    int id1, Energy mY, Energy mX,
					    complex<InvEnergy2> & h,Complex & k,
					    complex<InvEnergy2> & bp,
					    complex<InvEnergy2> & bm) const {
  Complex f1,f2,f3,f4;
  formFactor(q2,iloc,id0,id1,mY,mX,f1,f2,f3,f4);
  Energy msum(mX+mY);
  h = f1/sqr(msum);
  k = f2;
  bp = f3/sqr(msum);
  bm = f4/sqr(msum);
}

// member which does the work
void ISGWFormFactor::formFactor(Energy2 q2, unsigned int iloc, int, int id1,
				Energy mY,Energy mX, Complex & f1,Complex & f2,
				Complex & f3, Complex & f4) const {
  useMe();
  // work out the flavour of the heavy quarks etc
  int jspin,spect,inquark,outquark;
  formFactorInfo(iloc,jspin,spect,inquark,outquark);
  int ifl0(abs(inquark)),ifl1(abs(outquark)),ifls(abs(spect));
  // determine the multiplet
  int ispin(abs(id1)/1000);
  // masses of the quarks
  Energy mQ(_mquark[ifl0-1]),mq(_mquark[ifl1-1]),ms(_mquark[ifls-1]);
  Energy mtildeX(mq+ms),mtildeY(mQ+ms);
  // wavefunction parameters for the mesons
  Energy betaX(ZERO),betaY(_betaS[ifl0-1][ifls-1]);
  // spin-0 outgoing mesons
  if(ispin==0&&jspin<2) {
    betaX=_betaS[ifl1-1][ifls-1];
  }
  else {
    betaX=_betaP[ifl1-1][ifls-1];
  }
  // compute the F_n function we will need
  Energy2 beta2XY(0.5*(betaX*betaX+betaY*betaY)),tm((mY-mX)*(mY-mX));
  double betar(betaX*betaY/beta2XY),kappa2(_kappa*_kappa),
    slope((tm-q2)/(kappa2*beta2XY));
  Energy mup(mq*mQ/(mQ+mq)),mum(mq*mQ/(mQ-mq));
  double fn(sqrt(mtildeX/mtildeY)*betar*sqrt(betar)*
	    exp(-0.25*ms*ms/(mtildeX*mtildeY)*slope));
  // now we can compute the form-factors
  // for scalars
  if(jspin==0) {
    Complex fp,fm;
    // 1 1S0
    if(ispin==0) {
      double yratio(ms/mtildeX*betaY*betaY/beta2XY);
      fp = fn*(1.+0.5*mQ/mum-0.25*mQ*mq/mup/mum*yratio);
      fm = fn*(1.-(mtildeX+mtildeY)*(0.5/mq-0.25/mup*yratio));
    }
    // 1 3P0
    else if(ispin<100) {
      // extra power of beta factors
      fn*=betar;
      fp = fn*ms*mQ*mq/sqrt(6.)/betaY/mtildeX/mum;
      fm =-fn*ms/betaY/sqrt(6.)*(mtildeX+mtildeY)/mtildeX;
    }
    // 2 1S0
    else {
      throw Exception() << "ISGWFormFactor::formFactor" 
			<< " 2S not implemented" << Exception::abortnow;
    }
    // convert to the standard form
    f1 = Complex(q2/(mY*mY-mX*mX)*fm)+fp;
    f2 = fp;
  }
  // for vectors
  else if(jspin==1) {
    complex<Energy> f;
    complex<InvEnergy> g,ap,am;
    Energy2 betaX2(betaX*betaX),betaY2(betaY*betaY);
    //  1 3S1
    if(ispin==0) {
      f  = 2.*mtildeY*fn;
      g  = 0.5*fn*(1./mq-0.5/mum*ms/mtildeX*betaY2/beta2XY);
      ap =-0.5*fn/mtildeX*(1.+ms/mQ*(betaY2-betaX2)/(betaX2+betaY2)
			   -0.25*ms*ms/mum/mtildeY*betaX2*betaX2/beta2XY/beta2XY);
      am = 0.5*fn/mtildeX*(1.+ms/mQ+ms*ms/mq/mQ*betaX2/beta2XY*
 			   (1.-0.25*(mtildeX+mtildeY)/mtildeY*betaX2/beta2XY));
    }
    // 1 3P1
    else if(ispin==20) {
      fn*=betar;
      f  =-fn*mtildeY*betaY*(1./mum+0.5*ms/mtildeY*beta2XY*slope/betaY2*
			     (1./mq-0.5/mum*ms/mtildeX*betaY2/beta2XY));
      g  = 0.5*fn*ms/mtildeX/betaY;
      ap = 0.25*fn*ms*mQ/mtildeY/betaY/mum*(1.-0.5*ms*mq/mtildeX/mum*betaY2/beta2XY);
      am = -0.25*fn*ms*(mtildeX+mtildeY)/mq/betaY/mtildeY*
 	(1.-0.5*ms*mq/mtildeX/mum*betaY2/beta2XY);
    }
    //  1 1P1
    else if(ispin==10) {
      fn*=betar;
      double ort(1./sqrt(2.));
      f  = fn*ort*mtildeY*betaY/mup;
      g  = 0.25*fn*mtildeY*betaY*ort/mq/mQ/mtildeX;
      Energy msum(mtildeX+mtildeY);
      ap = fn*ms*ort/betaY/mtildeY*(1.+0.5*mQ/mum
				    -0.25*mq*mQ*ms/mum/mup/mtildeX*betaY2/beta2XY);
      am = fn*ms*ort/betaY/mtildeY*(1.-0.5/mq*msum
				    +0.25*ms*betaY2/mup/beta2XY*msum/mtildeX);
    }
    // 2 1S0
    else {
      throw Exception() << "ISGWFormFactor::formFactor" 
			<< " 2S not implemented" << Exception::abortnow;
    }
    // convert to the standard notation
    Energy msum(mX+mY),mdiff(mY-mX);
    Complex ii(0.,1.);
    f2 = -ii*f/msum;
    f3 = +ii*ap*msum;
    f1 = -ii*0.5/mX*(am*q2+ii*msum*f2-ii*mdiff*f3);
    f4 =  ii*g*msum;
  }
  // for tensors
  else if(jspin==2) {
    Energy msum(mX+mY);
    fn *=betar/sqrt(2.);
    double betaXb2(betaX*betaX/beta2XY);
    // 1 3P2
    if(ispin==0) {
      f1 = 0.5*fn*ms/mtildeY/betaY*(1./mq
				    -0.5*ms/mtildeX/mum*betaY*betaY/beta2XY)*sqr(msum);
      f2 = 2.*fn*ms/betaY;
      f3 =-0.5*fn*ms/mtildeX/mQ/betaY*
	(1.-0.5*ms*mQ/mup/mtildeY*betaXb2
	 +0.25*ms*mQ/mtildeY/mum*betaXb2*(1.-0.5*ms*betaXb2/mtildeY))* 
	sqr(msum);
      f4 = 0.5*fn*ms/mtildeX/mQ/betaY*
 	(1.-0.5*ms*mQ/mup/mtildeY*betaXb2+
 	 0.25*ms*betaXb2/mq*(mtildeX+mtildeY)/mtildeY*(1.-0.5*ms*betaXb2/mtildeY))*
	sqr(msum);
    }
  }
  else {
    throw Exception() << "ISGWFormFactor::FormFactor spin = " << jspin 
		      << " but spins higher than 2 not implemented"
		      << Exception::runerror;
  }
  // special for mixing
  double fact(1.);
  if(id1==ParticleID::eta) {
    if(ifl1==3&&ifls==3) fact = -2.*cos(_thetaeta)/sqrt(6.) - sin(_thetaeta)/sqrt(3.);
    else                 fact =     cos(_thetaeta)/sqrt(6.) - sin(_thetaeta)/sqrt(3.);
   
  }
  else if(id1==ParticleID::etaprime) {
    if(ifl1==3&&ifls==3) fact = -2.*sin(_thetaeta)/sqrt(6.) + cos(_thetaeta)/sqrt(3.);
    else                 fact =     sin(_thetaeta)/sqrt(6.) + cos(_thetaeta)/sqrt(3.);
  }
  else if(ifl1==ifls&&ifl1<3) {
    if(abs(ifl1)==1&&int(id1/10)%100==1) fact = -sqrt(0.5);
    else                                 fact =  sqrt(0.5);
  }
  f1*=fact;
  f2*=fact;
  f3*=fact;
  f4*=fact;
}

void ISGWFormFactor::dataBaseOutput(ofstream & output,bool header,bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::ISGWFormFactor " << name() << "\n";
  output << "newdef " << name() << ":Kappa    "    << _kappa        << "\n";
  output << "newdef " << name() << ":DownMass "    << _mdown/GeV    << "\n";
  output << "newdef " << name() << ":UpMass "      << _mup/GeV      << "\n";
  output << "newdef " << name() << ":StrangeMass " << _mstrange/GeV << "\n";
  output << "newdef " << name() << ":CharmMass "   << _mcharm/GeV   << "\n";
  output << "newdef " << name() << ":BottomMass "  << _mbottom/GeV  << "\n";
  output << "newdef " << name() << ":BetaSud "     << _betaSud/GeV  << "\n";
  output << "newdef " << name() << ":BetaSus "     << _betaSus/GeV  << "\n";
  output << "newdef " << name() << ":BetaSuc "     << _betaSuc/GeV  << "\n";
  output << "newdef " << name() << ":BetaSub "     << _betaSub/GeV  << "\n";
  output << "newdef " << name() << ":BetaPud "     << _betaPud/GeV  << "\n";
  output << "newdef " << name() << ":BetaPus "     << _betaPus/GeV  << "\n";
  output << "newdef " << name() << ":BetaPuc "     << _betaPuc/GeV  << "\n";
  output << "newdef " << name() << ":ThetaEtaEtaPrime " << _thetaeta  << "\n";
  ScalarFormFactor::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}
