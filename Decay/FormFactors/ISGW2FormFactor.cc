// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ISGW2FormFactor class.
//

#include "ISGW2FormFactor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ISGW2FormFactor.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;

ISGW2FormFactor::~ISGW2FormFactor() {}

void ISGW2FormFactor::persistentOutput(PersistentOStream & os) const {
  os <<_mdown << _mup << _mstrange << _mcharm << _mbottom << _mquark << _beta1S0ud 
     << _beta1S0us << _beta1S0ss << _beta1S0cu << _beta1S0cs << _beta1S0ub 
     << _beta1S0sb << _beta1S0cc << _beta1S0bc << _beta3S1ud << _beta3S1us 
     << _beta3S1ss << _beta3S1cu << _beta3S1cs << _beta3S1ub << _beta3S1sb 
     << _beta3S1cc << _beta3S1bc << _beta1Pud  << _beta1Pus  << _beta1Pss  
     << _beta1Pcu  << _beta1Pcs  << _beta1Pub  << _beta1Psb  << _beta1Pcc << _beta1Pbc
     << _alphamuQM << _alphaQ;
  for(unsigned int ix=0;ix<5;++ix)
    {
      for(unsigned int iy=0;iy<5;++iy)
	{
	  os << _beta1S0[ix][iy] << _mass1S0[ix][iy] << _beta3S1[ix][iy]
	     << _beta1P[ ix][iy] << _massPoh[ix][iy] << _massPth[ix][iy];
	}
    }
}

void ISGW2FormFactor::persistentInput(PersistentIStream & is, int) {
  is >>_mdown >> _mup >> _mstrange >> _mcharm >> _mbottom >> _mquark >> _beta1S0ud 
     >> _beta1S0us >> _beta1S0ss >> _beta1S0cu >> _beta1S0cs >> _beta1S0ub 
     >> _beta1S0sb >> _beta1S0cc >> _beta1S0bc >> _beta3S1ud >> _beta3S1us 
     >> _beta3S1ss >> _beta3S1cu >> _beta3S1cs >> _beta3S1ub >> _beta3S1sb 
     >> _beta3S1cc >> _beta3S1bc >> _beta1Pud  >> _beta1Pus  >> _beta1Pss  
     >> _beta1Pcu  >> _beta1Pcs  >> _beta1Pub  >> _beta1Psb  >> _beta1Pcc >> _beta1Pbc
     >> _alphamuQM >> _alphaQ;
  for(unsigned int ix=0;ix<5;++ix)
    {
      for(unsigned int iy=0;iy<5;++iy)
	{
	  is >> _beta1S0[ix][iy] >> _mass1S0[ix][iy] >> _beta3S1[ix][iy]
	     >> _beta1P[ ix][iy] >> _massPoh[ix][iy] >> _massPth[ix][iy];
	}
    }
}

ClassDescription<ISGW2FormFactor> ISGW2FormFactor::initISGW2FormFactor;
// Definition of the static class description member.

void ISGW2FormFactor::Init() {

  static Parameter<ISGW2FormFactor,Energy> interfaceDownMass
    ("DownMass",
     "The mass of the down quark in the ISGW model (this is a consituent mass)",
     &ISGW2FormFactor::_mdown, GeV, 0.33*GeV, 0.0*GeV, 1.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceUpMass
    ("UpMass",
     "The mass of the up quark in the ISGW model (this is a consituent mass)",
     &ISGW2FormFactor::_mup, GeV, 0.33*GeV, 0.0*GeV, 1.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceStrangeMass
    ("StrangeMass",
     "The mass of the strange quark in the ISGW model (this is a consituent mass)",
     &ISGW2FormFactor::_mstrange, GeV, 0.55*GeV, 0.0*GeV, 1.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceCharmMass
    ("CharmMass",
     "The mass of the charm quark in the ISGW model (this is a consituent mass)",
     &ISGW2FormFactor::_mcharm, GeV, 1.82*GeV, 0.0*GeV, 3.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBottomMass
    ("BottomMass",
     "The mass of the bottom quark in the ISGW model (this is a consituent mass)",
     &ISGW2FormFactor::_mbottom, GeV, 5.20*GeV, 3.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1S0ud
    ("Beta1S0ud",
     "The beta wavefunction parameter for the ud meson in the 1 1S0 level",
     &ISGW2FormFactor::_beta1S0ud, GeV, 0.41*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1S0us
    ("Beta1S0us",
     "The beta wavefunction parameter for the us meson in the 1 1S0 level",
     &ISGW2FormFactor::_beta1S0us, GeV, 0.44*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1S0ss
    ("Beta1S0ss",
     "The beta wavefunction parameter for the ss meson in the 1 1S0 level",
     &ISGW2FormFactor::_beta1S0ss, GeV, 0.53*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1S0cu
    ("Beta1S0cu",
     "The beta wavefunction parameter for the cu meson in the 1 1S0 level",
     &ISGW2FormFactor::_beta1S0cu, GeV, 0.45*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1S0cs
    ("Beta1S0cs",
     "The beta wavefunction parameter for the cs meson in the 1 1S0 level",
     &ISGW2FormFactor::_beta1S0cs, GeV, 0.56*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1S0ub
    ("Beta1S0ub",
     "The beta wavefunction parameter for the ub meson in the 1 1S0 level",
     &ISGW2FormFactor::_beta1S0ub, GeV, 0.43*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1S0sb
    ("Beta1S0sb",
     "The beta wavefunction parameter for the sb meson in the 1 1S0 level",
     &ISGW2FormFactor::_beta1S0sb, GeV, 0.54*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1S0cc
    ("Beta1S0cc",
     "The beta wavefunction parameter for the cc meson in the 1 1S0 level",
     &ISGW2FormFactor::_beta1S0cc, GeV, 0.88*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1S0bc
    ("Beta1S0bc",
     "The beta wavefunction parameter for the bc meson in the 1 1S0 level",
     &ISGW2FormFactor::_beta1S0bc, GeV, 0.92*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1Pud
    ("Beta1Pud",
     "The beta wavefunction parameter for the ud meson in the 1P level",
     &ISGW2FormFactor::_beta1Pud, GeV, 0.28*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1Pus
    ("Beta1Pus",
     "The beta wavefunction parameter for the us meson in the 1P level",
     &ISGW2FormFactor::_beta1Pus, GeV, 0.30*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1Pss
    ("Beta1Pss",
     "The beta wavefunction parameter for the ss meson in the 1P level",
     &ISGW2FormFactor::_beta1Pss, GeV, 0.33*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1Pcu
    ("Beta1Pcu",
     "The beta wavefunction parameter for the cu meson in the 1P level",
     &ISGW2FormFactor::_beta1Pcu, GeV, 0.33*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1Pcs
    ("Beta1Pcs",
     "The beta wavefunction parameter for the cs meson in the 1P level",
     &ISGW2FormFactor::_beta1Pcs, GeV, 0.38*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1Pub
    ("Beta1Pub",
     "The beta wavefunction parameter for the ub meson in the 1P level",
     &ISGW2FormFactor::_beta1Pub, GeV, 0.35*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1Psb
    ("Beta1Psb",
     "The beta wavefunction parameter for the sb meson in the 1P level",
     &ISGW2FormFactor::_beta1Psb, GeV, 0.41*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1Pcc
    ("Beta1Pcc",
     "The beta wavefunction parameter for the cc meson in the 1P level",
     &ISGW2FormFactor::_beta1Pcc, GeV, 0.52*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1Pbc
    ("Beta1bc",
     "The beta wavefunction parameter for the bc meson in the 1P level",
     &ISGW2FormFactor::_beta1Pbc, GeV, 0.60*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta3S1ud
    ("Beta3S1ud",
     "The beta wavefunction parameter for the ud meson in the 3S1 level",
     &ISGW2FormFactor::_beta3S1ud, GeV, 0.30*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta3S1us
    ("Beta3S1us",
     "The beta wavefunction parameter for the us meson in the 3S1 level",
     &ISGW2FormFactor::_beta3S1us, GeV, 0.33*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta3S1ss
    ("Beta3S1ss",
     "The beta wavefunction parameter for the ss meson in the 3S1 level",
     &ISGW2FormFactor::_beta3S1ss, GeV, 0.37*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta3S1cu
    ("Beta3S1cu",
     "The beta wavefunction parameter for the cu meson in the 3S1 level",
     &ISGW2FormFactor::_beta3S1cu, GeV, 0.38*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta3S1cs
    ("Beta3S1cs",
     "The beta wavefunction parameter for the cs meson in the 3S1 level",
     &ISGW2FormFactor::_beta3S1cs, GeV, 0.44*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta3S1ub
    ("Beta3S1ub",
     "The beta wavefunction parameter for the ub meson in the 3S1 level",
     &ISGW2FormFactor::_beta3S1ub, GeV, 0.40*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta3S1sb
    ("Beta3S1sb",
     "The beta wavefunction parameter for the sb meson in the 3S1 level",
     &ISGW2FormFactor::_beta3S1sb, GeV, 0.49*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta3S1cc
    ("Beta3S1cc",
     "The beta wavefunction parameter for the cc meson in the 3S1 level",
     &ISGW2FormFactor::_beta3S1cc, GeV, 0.62*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta3S1bc
    ("Beta1bc",
     "The beta wavefunction parameter for the bc meson in the 3S1 level",
     &ISGW2FormFactor::_beta3S1bc, GeV, 0.75*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,double> interfaceCfDrho
    ("CfDrho",
     "The relativistic correction factor for D -> rho",
     &ISGW2FormFactor::_CfDrho, 0.889, 0.0, 1.0,
     false, false, true);

  static Parameter<ISGW2FormFactor,double> interfaceCfDKstar
    ("CfDKstar",
     "The relativistic correction factor for D -> Kstar",
     &ISGW2FormFactor::_CfDKstar, 0.928, 0.0, 1.0,
     false, false, true);

  static Parameter<ISGW2FormFactor,double> interfaceCfDsKstar
    ("CfDsKstar",
     "The relativistic correction factor for Ds -> Kstar",
     &ISGW2FormFactor::_CfDsKstar, 0.873, 0.0, 1.0,
     false, false, true);

  static Parameter<ISGW2FormFactor,double> interfaceCfDsphi
    ("CfDsphi",
     "The relativistic correction factor for Ds -> phi",
     &ISGW2FormFactor::_CfDsphi, 0.911, 0.0, 1.0,
     false, false, true);

  static Parameter<ISGW2FormFactor,double> interfaceCfBrho
    ("CfBrho",
     "The relativistic correction factor for B -> rho",
     &ISGW2FormFactor::_CfBrho, 0.905, 0.0, 1.0,
     false, false, true);

  static Parameter<ISGW2FormFactor,double> interfaceCfBDstar
    ("CfBDstar",
     "The relativistic correction factor for B -> Dstar",
     &ISGW2FormFactor::_CfBDstar, 0.989, 0.0, 1.0,
     false, false, true);

  static Parameter<ISGW2FormFactor,double> interfaceCfBsKstar
    ("CfBsKstar",
     "The relativistic correction factor for Bs -> Kstar",
     &ISGW2FormFactor::_CfBsKstar, 0.892, 0.0, 1.0,
     false, false, true);

  static Parameter<ISGW2FormFactor,double> interfaceCfBsDstar
    ("CfBsDstar",
     "The relativistic correction factor for Bs -> Dstar",
     &ISGW2FormFactor::_CfBsDstar, 0.984, 0.0, 1.0,
     false, false, true);

  static Parameter<ISGW2FormFactor,double> interfaceCfBcDstar
    ("CfBcDstar",
     "The relativistic correction factor for Bc -> Dstar",
     &ISGW2FormFactor::_CfBcDstar, 0.868, 0.0, 1.0,
     false, false, true);

  static Parameter<ISGW2FormFactor,double> interfaceCfBcpsi
    ("CfBcpsi",
     "The relativistic correction factor for Bc -> psi",
     &ISGW2FormFactor::_CfBcpsi, 0.967, 0.0, 1.0,
     false, false, true);

  static Parameter<ISGW2FormFactor,double> interfaceCfBcBsstar
    ("CfBcBsstar",
     "The relativistic correction factor for Bc -> Bsstar",
     &ISGW2FormFactor::_CfBcBsstar, 1.000, 0.0, 1.0,
     false, false, true);

  static Parameter<ISGW2FormFactor,double> interfaceCfBcBstar
    ("CfBcBstar",
     "The relativistic correction factor for Bc -> Bstar",
     &ISGW2FormFactor::_CfBcBstar, 1.000, 0.0, 1.0,
     false, false, true);

}

// member which does the work
void ISGW2FormFactor::formFactor(Energy2 q2, int iloc, int id0, int id1, Energy mY,
				 Energy mX, Complex & f1,Complex & f2,Complex & f3,
				 Complex & f4) const
{
  // work out the flavours and spectators
  int  inf1=(abs(id0)/100)%10;
  int  inf2=(abs(id0)/10)%10;
  int outf1=(abs(id1)/100)%10;
  int outf2=(abs(id1)/10)%10;
  int id1f=abs(id1)%1000;
  int ifl0,ifls,ifl1;
  if(id1f==ParticleID::pi0  || id1f==ParticleID::eta   || id1f==ParticleID::etaprime ||
     id1f==ParticleID::rho0 || id1f==ParticleID::omega || id1f==ParticleID::phi)
    {
      ifl0=inf1;
      if((id1f==ParticleID::eta||id1f==ParticleID::etaprime||
	  id1f==ParticleID::phi)&&inf2==3)
	{ifl1=3;}
      else
	{ifl1=inf1-3;}
      ifls=inf2;
    }
  else if(inf1==outf1+1&&inf2==outf2&&outf1>=3)
    {ifl0=inf1;ifl1=outf1;ifls=inf2;}
  else if(inf1==outf1&&inf2==outf2+1&&outf2>=3)
    {ifl0=inf2;ifl1=outf2;ifls=inf1;}
  else if(inf1==outf1&&inf2==outf2+3)
    {ifl0=inf2;ifl1=outf2;ifls=inf1;}
  else if(inf1==outf2+3&&inf2==outf1)
    {ifl0=inf1;ifl1=outf2;ifls=inf2;}
  else if(inf1==outf1+3&&inf2==outf2)
    {ifl0=inf1;ifl1=outf1;ifls=inf2;}
  else
    {throw Exception() << "Unknown flavour combination in "
		       << "ISGW2FormFactor::formFactor"
		       << Exception::abortnow;}
  // determine the multiplet
  int ispin(abs(id1)/1000),jspin(abs(id1)%10);
  // masses of the quarks
  Energy mQ=_mquark[ifl0-1];
  Energy mq=_mquark[ifl1-1];
  Energy ms=_mquark[ifls-1];
  // of the mesons
  Energy mtildeX=mq+ms;
  Energy mtildeY=mQ+ms;
  Energy mup=mq*mQ/(mQ+mq);
  Energy mum=mq*mQ/(mQ-mq);
  // wavefunction parameters for the mesons
  Energy betaX(0.),betaY(0.),mbarX(0.),mbarY(0.);
  betaY = _beta1S0[ifl0-1][ifls-1];
  mbarY = _mass1S0[ifl0-1][ifls-1];
  double Cf(1.);
  // the wavefunction parameter for the outgoing meson
  // 1S0
  if(ispin==0&&jspin==1)
    {
      betaX=_beta1S0[ifl1-1][ifls-1];
      mbarX=_mass1S0[ifl1-1][ifls-1];
    }
  // 1S1
  else if(ispin==0&&jspin==3)
    {
      betaX = _beta3S1[ifl1-1][ifls-1];
      mbarX = _mass1S0[ifl1-1][ifls-1];
      // set the relativistic correction parameter
      // decaying b
      if(ifl0==5)
	{
	  if(ifls<3)
	    {
	      if(ifl1<3){Cf=_CfBrho;}
	      else {Cf=_CfBDstar;}
	    }
	  else if(ifls==3)
	    {
	      if(ifl1==4){Cf=_CfBsDstar;}
	      else{Cf=_CfBsKstar;}
	    }
	  else if(ifls==4)
	    {
	      if(ifl1==4){Cf=_CfBcpsi;}
	      else {Cf=_CfBcDstar;}
	    }
	}
      else if(ifl0==4)
	{
	  if(ifls<3)
	    {
	      if(ifl1<3){Cf=_CfDrho;}
	      else {Cf=_CfDKstar;}
	    }
	  else if(ifls==3)
	    {
	      if(ifl1<3){Cf=_CfDsKstar;}
	      else{Cf=_CfDsphi;}
	    }
	  else if(ifls==5)
	    {
	      if(ifl1<3){Cf=_CfBcBstar;}
	      else{Cf=_CfBcBsstar;}
	    }
	} 
    }
  else if(ispin==10&&jspin==1)
    {
      betaX=_beta1P[ifl1-1][ifls-1];
      mbarX=_massPoh[ifl1-1][ifls-1];
    }
  // 1 3/2 P 1 (1 P1)
  else if(ispin==10&&jspin==3)
    {
      betaX = _beta1P[ifl1-1][ifls-1];
      mbarX=_massPth[ifl1-1][ifls-1];
    }
  // 1 1/2 P1 ( 3 P1) 
  else if(ispin==20&&jspin==3)
    {
      betaX = _beta1P[ifl1-1][ifls-1];
      mbarX=_massPoh[ifl1-1][ifls-1];
    }
  else
    {throw Exception() << "ISGWS2FormFactor::formFactor" 
		       << " unknown multiplet" << Exception::abortnow;
    }
  Energy2 beta2XY = 0.5*(betaX*betaX+betaY*betaY);
  // number of active flavours
  int Nf  = ifl0-1;
  int Nfp = ifl1-1; 
  // first piece of the f_n function
  double betar = betaX*betaY/beta2XY;
  double fn = sqrt(mtildeX/mtildeY)*betar*sqrt(betar);
  // q dependent piece
  Energy2 tm=(mY-mX)*(mY-mX);
  Energy2 tmmt=(tm-q2);
  // radius parameter
  InvEnergy2 r2 = 0.75/mQ/mq+1.5*ms*ms/mbarX/mbarY/beta2XY
    +16./mbarX/mbarY/(33-2.*Nfp)*log(_alphamuQM/_alphaQ[ifl1-1]);
  // the parameters for the form-factor depenedent piece
  double rmbmtY=sqrt(mbarY/mtildeY);
  double rmbmtX=sqrt(mbarX/mtildeX);
  // work out wtilde
  double wt = 1.+0.5*tmmt/mbarX/mbarY;
  // storage of the form factors
  double fpmfm(0.),fppfm(0.),f(0.),g(0.),appam(0.),apmam(0.),
    h(0.),k(0.),bp(0.),bm(0.);
  // scalar and vector from 1S levels
  if(ispin==0)
    {
      // parameters for the beta functions
      double asopi  = min(_alphamuQM,alphaS(mq*mQ))/pi;  
      double w   = 1.+0.5*tmmt/mX/mY;
      double aI  = -6./(33.-2.*Nf);
      double rw  = 1./sqrt(w*w-1)*log(w+sqrt(w*w-1));
      double aLw = 8./(33.-2.*Nfp)*(w*rw-1.); 
      double cji = pow(_alphaQ[ifl0-1]/_alphaQ[ifl1-1],aI)*
	pow(_alphaQ[ifl1-1]/_alphamuQM,aLw);
      double zji   = mq/mQ; 
      double gamji =-2.*zji/(1.-zji)*log(zji)-2.;
      double chiji = -1.-gamji/(1.-zji);
      // scalar
      if(jspin==1)
	{
	  double fact=(1.+1./12.*r2*tmmt);
	  fn/=(fact*fact);
	  fact = (1.-0.5*ms*mq/mup/mtildeX*betaY*betaY/beta2XY);
	  fppfm = fn*rmbmtX/rmbmtY*cji*(1.+asopi*(gamji-2./3.*chiji))*
	    (2.-mtildeX/mq*fact);
	  fpmfm = fn*rmbmtY/rmbmtX*cji*(1.+asopi*(gamji+2./3.*chiji))*mtildeY/mq*fact;
	}
      else if(jspin==3)
	{
	  // factors for the F and R functions
	  double fact=(1.+1./12.*r2*tmmt);
	  fn/=(fact*fact);
	  double betaapmam=1./3.-4./3./(1-zji)-chiji
	    +gamji*(1.-2./3.*(1.+zji)/(1.-zji)/(1.-zji));
	  f     = Cf*fn*rmbmtX*rmbmtY*cji*(1.+asopi*(-2./3.+gamji));
	  g     = fn/rmbmtX/rmbmtY*cji*(1.+asopi*( 2./3.+gamji));
	  appam = fn*rmbmtX/(rmbmtY*rmbmtY*rmbmtY)*cji;
	  apmam = fn/rmbmtX/rmbmtY*cji*(1.+asopi*betaapmam);
	  // rest of the calculation
	  f     *=mtildeY*(1.+wt+0.5*ms*(wt-1.)/mup);
	  g     *=0.5*(1./mq-0.5*ms*betaY*betaY/mum/mtildeX/beta2XY);
	  appam *=(ms*betaX*betaX/(1.+wt)/mq/mQ/beta2XY*
		   (1.-0.5*ms*betaX*betaX/mtildeY/beta2XY)
		   +asopi/mtildeY*(-1.-chiji+4./3./(1.-zji)
				   +2./3.*(1.+zji)/(1.-zji)/(1.-zji)*gamji));
	  apmam *=-1./mtildeX*(mtildeY/mQ
			       -0.5*ms*betaX*betaX/mup/beta2XY
			       +wt*ms*mtildeY*betaX*betaX/(wt+1.)/mq/mQ/beta2XY*
			       (1.-0.5*ms/mtildeY*betaX*betaX/beta2XY)); 
	}
      else if(jspin==5)
	{
	  // factors for the F function
	  double fact=(1.+1./18.*r2*tmmt);
	  fn*=betar/(fact*fact*fact);
	  h = fn/rmbmtX/(rmbmtY*rmbmtY*rmbmtY);
	  k = fn*rmbmtX/rmbmtY;
	  double bppbm = fn*rmbmtX/(rmbmtY*rmbmtY*rmbmtY*rmbmtY*rmbmtY);
	  double bpmbm = fn/rmbmtX/(rmbmtY*rmbmtY*rmbmtY);
	  // functions themselves
	  double or2=sqrt(0.5);
	  h *= 0.5*ms*or2/mtildeY/betaY*(1./mq-0.5*ms*betaY*betaY/mum/mtildeX/beta2XY);
	  k *= or2*ms/betaY*(1.+wt);
	  bppbm *= 0.25*or2*ms*ms/mq/mQ/mtildeY/betaY*betaX*betaX/beta2XY*
	    (1.-0.5*ms/mtildeY*betaX*betaX/beta2XY);
	  bpmbm *= -or2*ms/mQ/mtildeX/betaY*
	    (1.-0.5*ms*mQ/mup/mtildeY*betaX*betaX/beta2XY
	     +0.25*ms/mq*betaX*betaX/beta2XY*(1.-0.5*ms/mtildeY*betaX*betaX/beta2XY));
	  // conversion
	  bp = 0.5*(bppbm+bpmbm);
	  bm = 0.5*(bppbm-bpmbm);
	}
    }
  // 1 3P0
  else if(ispin==10&&jspin==1)
    {
      fn*=betar;
      double fact=(1.+1./18.*r2*tmmt);
      fn/=(fact*fact*fact);
      fact = sqrt(2./3.)*ms/betaY;
      fppfm =-fn*rmbmtX/rmbmtY;
      fpmfm = fn*rmbmtY/rmbmtX*mtildeY/mtildeX;
    }

  // 1 3/2 P1 ( 1 P1) 
  else if(ispin==10&&jspin==3)
    {
      // factors for the F and R functions
      double fact=(1.+1./18.*r2*tmmt);
      fn*=betar/(fact*fact*fact);
      f     = fn*rmbmtX*rmbmtY;
      g     = fn/rmbmtX/rmbmtY;
      appam = fn*rmbmtX/(rmbmtY*rmbmtY*rmbmtY);
      apmam = fn/rmbmtX/rmbmtY;
      // light meson
      if(ifls<3&&ifl1<3)
	{
	  f     *=-mtildeY*betaY*(1./mum
				  +ms*mtildeX*(wt-1.)/betaY/betaY*
				  ((5.+wt)/6./mq-0.5/mum*ms/mtildeX*
				   betaY*betaY/beta2XY)); 
	  g     *=-0.5*ms/mtildeX/betaY*(5.+wt)/6.;
	  appam *=-0.5*ms*mtildeX/mq/mtildeY/betaY*
	    (1.-0.5*ms*mq/mtildeX/mum*betaY*betaY/beta2XY); 
	  apmam *= 0.5*ms/mq/betaY*((wt+2.)/3.
				    -0.5*ms*mq/mtildeX/mum*betaY*betaY/beta2XY);
	}
      // heavy meson
      else
	{
	  double oor3=1./sqrt(3.);
	  f     *=-2.*oor3*mtildeY*betaY*
	    (1.0/mq+0.5*mtildeX*ms*(wt-1.)/betaY/betaY*
	     (0.5*(wt+1.)/mq-0.5*ms*betaY*betaY/mum/mtildeX/beta2XY)); 
	  g     *=-0.5*oor3*(0.5*(1.+wt)+0.5*betaY*betaY*mtildeY/ms/mq/mQ)*ms
	    /betaY/mtildeY;
	  appam *=-0.5/oor3*ms/betaY/mtildeY*(1.-ms/3./mq-ms/3.*betaY*betaY/beta2XY*
					      (0.5/mum-1./mup));
	  apmam *=-0.5*oor3*ms/betaY/mtildeX*((2.-wt)*mtildeX/mq+ms*betaY*betaY/beta2XY*
					      (0.5/mum-1./mup));
	}
    }
  // 1 1/2 P 1 (3 P1)
  else if(ispin==20&&jspin==3)
    {
      // factors for the F and R functions
      double fact=(1.+1./18.*r2*tmmt);
      fn*=betar/(fact*fact*fact);
      f     = fn*rmbmtX*rmbmtY;
      g     = fn/rmbmtX/rmbmtY;
      appam = fn*rmbmtX/(rmbmtY*rmbmtY*rmbmtY);
      apmam = fn/rmbmtX/rmbmtY;
      // light meson
      if(ifls<3&&ifl1<3)
	{
	  double oor2 = 1./sqrt(2.);
	  f     *= oor2*mtildeY*betaY*(1./mup
				      +ms*mtildeX/3./mq/betaY/betaY*(wt-1.)*(wt-1.));
	  g     *= oor2*(0.25*mtildeY*betaY/mQ/mq/mtildeX+(wt-1.)*ms/6./mtildeX/betaY);
	  appam *= oor2*ms/mtildeY/betaY*(1./ms/mq+0.5*ms/mup*betaY*betaY/beta2XY);
	  apmam *= oor2*ms/mq/betaY*((4.-wt)/3.
				     -0.5*ms*mq/mtildeX/mup*betaY*betaY/beta2XY);
	}
      // heavy meson
      else
	{
	  double r2o3=sqrt(2./3.);
	  f     *= r2o3*mtildeY*betaY*(0.5/mq-1.5/mQ+ms*mtildeX*(wt-1.)/betaY/betaY*
				       (1./mq-0.5*ms*betaY*betaY/mum/mtildeX/beta2XY));
	  g     *= 0.5*r2o3*ms/betaY/mtildeX*(1.-0.25*betaY*betaY*mtildeY/ms/mq/mQ);
	  appam *= 0.5*r2o3*ms*ms*betaX*betaX/mtildeY/mq/betaY/beta2XY; 
	  apmam *= -r2o3*ms/mtildeX/betaY*(1.+0.5*ms*betaX*betaX/mq/beta2XY);
	}
    }
  else
    {throw Exception() << "ISGWS2FormFactor::formFactor" 
		       << " unknown multiplet" << Exception::abortnow;}
  // the final manipulations
  if(jspin==1)
    {
      double fp,fm;
      fp = 0.5*(fppfm+fpmfm);
      fm = 0.5*(fppfm-fpmfm);
      // convert to the standard form
      f1 = q2/(mY*mY-mX*mX)*fm+fp;
      f2 = fp;
    }
  else if(jspin==3)
    {
      double ap = 0.5*(appam+apmam);
      double am = 0.5*(appam-apmam);
      // convert to the standard notation
      Energy msum=mX+mY,mdiff=mY-mX;
      Complex ii(0.,1.);
      f2 = -ii*f/msum;
      f3 = +ii*ap*msum;
      f1 = -ii*0.5/mX*(am*q2-ii*msum*f2+ii*mdiff*f3);
      f4 = -ii*g*msum;
    }
  else if(jspin==5)
    {
      f1 = h;
      f2 = k;
      f3 = bp;
      f4 = bm;
    }
}

// form-factor for scalar to scalar
void ISGW2FormFactor::ScalarScalarFormFactor(Energy2 q2, int iloc,int id0, int id1,
					     Energy mY, Energy mX,
					     Complex & f0,Complex & fp) const
 {
   Complex d1(0.),d2(0.);
   formFactor(q2,iloc,id0,id1,mY,mX,f0,fp,d1,d2);
 }

// form-factor for scalar to vector
void ISGW2FormFactor::ScalarVectorFormFactor(Energy2 q2, int iloc, int id0, int id1, 
					     Energy mY, Energy mX,
					     Complex & A0,Complex & A1,
					     Complex & A2,Complex & V) const
 {formFactor(q2,iloc,id0,id1,mY,mX,A0,A1,A2,V);}


// form-factor for scalar to tensor
void ISGW2FormFactor::ScalarTensorFormFactor(Energy2 q2, int iloc, int id0, int id1, 
					     Energy mY, Energy mX,
					     Complex & h,Complex & k,
					     Complex & bp,Complex & bm) const
 {formFactor(q2,iloc,id0,id1,mY,mX,h,k,bp,bm);}

}
