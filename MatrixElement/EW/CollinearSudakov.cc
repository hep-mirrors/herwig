// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the CollinearSudakov class.
//

#include "CollinearSudakov.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "GroupInvariants.h"
#include "ElectroWeakReweighter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

CollinearSudakov::CollinearSudakov() : K_ORDER_(3),
				       B_ORDER_(2)
{}

CollinearSudakov::~CollinearSudakov() {}

IBPtr CollinearSudakov::clone() const {
  return new_ptr(*this);
}

IBPtr CollinearSudakov::fullclone() const {
  return new_ptr(*this);
}

void CollinearSudakov::persistentOutput(PersistentOStream & os) const {
  os << K_ORDER_ << B_ORDER_;
}

void CollinearSudakov::persistentInput(PersistentIStream & is, int) {
  is >> K_ORDER_ >> B_ORDER_;
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<CollinearSudakov,Interfaced>
describeHerwigCollinearSudakov("Herwig::CollinearSudakov", "HwMEEW.so");

void CollinearSudakov::Init() {

  static ClassDocumentation<CollinearSudakov> documentation
    ("The CollinearSudakov class implements the collinear evolution");


}


Complex CollinearSudakov::highScaleIntegral(bool SU3, bool SU2, double Y,
					    Energy2 s, Energy mu_h, Energy mu_l,
					    bool fermion, bool longitudinal,
					    double yukFactor) {
  SU3_          = SU3;
  SU2_          = SU2;
  Y_            = Y;
  s_            = s;
  fermion_      = fermion;
  longitudinal_ = longitudinal;
  yukFactor_    = yukFactor;
  // perform the integral
  Complex result;
  high_ = true;
  // real piece
  real_ = true;
  result.real(integrator_.value(*this,mu_h,mu_l));
  // imaginary piece
  real_ = false;
  result.imag(integrator_.value(*this,mu_h,mu_l));
  // return the answer
  return exp(result);
}

Complex  CollinearSudakov::lowScaleIntegral(bool SU3, double Q, Energy2 s, 
					    Energy mu_h, Energy mu_l, bool fermion, 
					    double boostFactor) {
  SU3_          = SU3;
  Q_            = Q;
  s_            = s;
  fermion_      = fermion;
  boostFactor_  = boostFactor;
  // perform the integral
  Complex result;
  high_ = false;
  // real piece
  real_ = true;
  result.real(integrator_.value(*this,mu_h,mu_l));
  // imaginary piece
  real_ = false;
  result.imag(integrator_.value(*this,mu_h,mu_l));
  // return the answer
  return exp(result);
}

InvEnergy CollinearSudakov::highScaleIntegrand(Energy mu) const {
  using namespace GroupInvariants;
  using Constants::pi;
  Complex gamma = 0.0;
  // Include K-factor Contributions (Cusps):
  GaugeContributions cusp = cuspContributions(mu,K_ORDER_,high_);
  // Include B-factors (B1 for U1, B2 for SU2, B3 for SU3):
  GaugeContributions nonCusp = BContributions(mu,B_ORDER_,fermion_,longitudinal_);
  // common log
  Complex plog = PlusLog(s_/sqr(mu));
  // SU(3)
  if(SU3_) {
    if(fermion_)
      gamma += C_F(3)*cusp.SU3*0.5*plog + 0.5*nonCusp.SU3;
    else
      gamma += (C_A(3))*cusp.SU3*0.5*plog + 0.5*nonCusp.SU3;
  }
  // SU(2)
  if(SU2_) {
    if (fermion_ || longitudinal_ )
      gamma += (C_F(2))*cusp.SU2*0.5*plog + 0.5*nonCusp.SU2;
    else
      gamma += (C_A(2))*cusp.SU2*0.5*plog + 0.5*nonCusp.SU2;
  }
  
  if (fermion_ || longitudinal_ ) {
    gamma += sqr(Y_)*(cusp.U1*0.5*plog + 0.5*nonCusp.U1);
  }
  else {
    // U(1) Gauge boson
    if (!SU3_ && !SU2_ && abs(Y_)<0.001) {
      gamma += 0.5*nonCusp.U1;
    }
  }
  // top Yukawa piece
  double y_t = ElectroWeakReweighter::coupling()->y_t(mu);
  gamma += yukFactor_*sqr(y_t)/(16.0*sqr(pi));
  // return the answer
  return real_ ? gamma.real()/mu : gamma.imag()/mu;
}

InvEnergy CollinearSudakov::lowScaleIntegrand(Energy mu) const {
  using namespace GroupInvariants;
  using Constants::pi;
  Complex gamma = 0.0;
  // Include K-factor Contributions (Cusps):
  GaugeContributions cusp = cuspContributions(mu,K_ORDER_,false);
  // Include B-factors (B1 for U1, B2 for SU2, B3 for SU3):
  GaugeContributions nonCusp = BContributionsLow(mu,B_ORDER_,fermion_,boostFactor_);
  // common log
  Complex plog = PlusLog(s_/sqr(mu));
  Complex blog(0.);
  if(boostFactor_ >0.001) blog = PlusLog(4.0*boostFactor_);
  // SU(3)
  if (SU3_) {
    if (fermion_) {
      // not a bHQET top quark field
      if (abs(boostFactor_)<0.001)
 	gamma += C_F(3)*cusp.SU3*0.5*plog + 0.5*nonCusp.SU3;
      else
  	gamma += C_F(3)*cusp.SU3*0.5*blog + 0.5*nonCusp.SU3;
    }
    else {
      gamma += C_A(3)*cusp.SU3*0.5*plog + 0.5*nonCusp.SU3;
    }
  }
  // fermions
  if (fermion_) {
    // not a bHQET top quark field
    if (abs(boostFactor_)<0.001)
      gamma += sqr(Q_)*(cusp.U1*0.5*plog + 0.5*nonCusp.U1);
    else
      gamma += sqr(Q_)*(cusp.U1*0.5*blog + 0.5*nonCusp.U1);
  }
  else {
    // i.e., not a fermion, not a bHQ, not a gluon => photon
    if (abs(boostFactor_)<0.001 && !SU3_)
      gamma += 0.5*nonCusp.U1;
    // i.e., W treated as a bHQET field
    else if (abs(boostFactor_)>0.001) 
      gamma += sqr(Q_)*(cusp.U1*0.5*blog + 0.5*nonCusp.U1);
  }
  // return the answer
  return real_ ? gamma.real()/mu : gamma.imag()/mu;
}

void CollinearSudakov::evaluateHighScale(Energy highScale, Energy EWScale, Energy2 S) {
  double yCoeffQt(0.), yCoefftR(0.), yCoeffPhi(0.);
  if (K_ORDER_ >= 1) {
    yCoeffQt = 0.5;
    yCoefftR = 1.0;
    yCoeffPhi = 3.0;
  }
  highColW_   = highScaleIntegral(false,true ,0.0   ,S,highScale,EWScale,false,false,0.0);
  highColB_   = highScaleIntegral(false,false,0.0   ,S,highScale,EWScale,false,false,0.0);
  highColG_   = highScaleIntegral(true ,false,0.0   ,S,highScale,EWScale,false,false,0.0);
  highColQ_   = highScaleIntegral(true ,true ,1./6. ,S,highScale,EWScale,true,false,0.0);
  highColQt_  = highScaleIntegral(true ,true ,1./6. ,S,highScale,EWScale,true,false,yCoeffQt);
  highColU_   = highScaleIntegral(true ,false,2./3. ,S,highScale,EWScale,true,false,0.0);
  highColtR_  = highScaleIntegral(true ,false,2./3. ,S,highScale,EWScale,true,false,yCoefftR);
  highColD_   = highScaleIntegral(true ,false,-1./3.,S,highScale,EWScale,true,false,0.0);
  highColL_   = highScaleIntegral(false,true ,-1./2.,S,highScale,EWScale,true,false,0.0);
  highColE_   = highScaleIntegral(false,false,-1.   ,S,highScale,EWScale,true,false,0.0);
  highColPhi_ = highScaleIntegral(false,true ,1./2. ,S,highScale,EWScale,false,true,yCoeffPhi);
}

void CollinearSudakov::evaluateLowScale(Energy EWScale, Energy lowScale, Energy2 S) {
  lowColW_ = lowScaleIntegral(false,1.    ,S,EWScale,lowScale,false,
			      S/sqr(2.*ElectroWeakReweighter::coupling()->mW()));
  lowColA_ = lowScaleIntegral(false,0.    ,S,EWScale,lowScale,false,0.0);
  lowColG_ = lowScaleIntegral(true ,0.    ,S,EWScale,lowScale,false,0.0);
  lowColU_ = lowScaleIntegral(true ,2./3. ,S,EWScale,lowScale,true,0.0);
  lowColt_ = lowScaleIntegral(true ,2./3. ,S,EWScale,lowScale,true,
			      S/sqr(2.*ElectroWeakReweighter::coupling()->mT()));
  lowColD_ = lowScaleIntegral(true ,-1./3.,S,EWScale,lowScale,true,0.0);
  lowColE_ = lowScaleIntegral(false,-1.   ,S,EWScale,lowScale,true,0.0);
}

void CollinearSudakov::evaluateMatching(Energy EWScale,Energy2 s) {
  using Constants::pi;
  using GroupInvariants::PlusLog;
  static const Complex I(0,1.0);
  // wave function corrections  
  WaveFunctionCorrections WFC = waveFunctionCorrections(EWScale);
   
  double aS   = ElectroWeakReweighter::coupling()->a3(EWScale);
  double aEM  = ElectroWeakReweighter::coupling()->aEM(EWScale);
  double aW   = ElectroWeakReweighter::coupling()->aW(EWScale);
  double aZ   = ElectroWeakReweighter::coupling()->aZ(EWScale);
  double cW2  = ElectroWeakReweighter::coupling()->Cos2thW(EWScale);
  double sW2  = ElectroWeakReweighter::coupling()->Sin2thW(EWScale);
  Energy mW   = ElectroWeakReweighter::coupling()->mW();
  Energy mZ   = ElectroWeakReweighter::coupling()->mZ();
  Energy mH   = ElectroWeakReweighter::coupling()->mH();
  Energy mT   = ElectroWeakReweighter::coupling()->mT();
  double gPHI = ElectroWeakReweighter::coupling()->g_phiPlus(EWScale);
  double y_t  = ElectroWeakReweighter::coupling()->y_t(EWScale);
   
  double lz = log(mZ/EWScale);
  double lw = log(mW/EWScale);

   complex<double> F_W = 4.0*lw*0.5*PlusLog(s/EWScale/EWScale)-2.0*lw*lw-2.0*lw-5.0*pi*pi/12.0+1.0;
   complex<double> F_Z = 4.0*lz*0.5*PlusLog(s/EWScale/EWScale)-2.0*lz*lz-2.0*lz-5.0*pi*pi/12.0+1.0;
   
   // Taken from Manohar... along with his formulae for F_tL, F_tR, and F_bL (for 3rd/1st Gen. Wavefunction Differences)
   complex<double> W1 = WFC.fFW0-0.5*WFC.aW0-0.5*WFC.cW0;
   complex<double> W2 = WFC.fF0W-0.5*WFC.a0W;
   complex<double> U1 = ElectroWeakReweighter::coupling()->g_Lu(EWScale)*ElectroWeakReweighter::coupling()->g_Lu(EWScale)*(WFC.fFZZ-0.5*WFC.aZZ) - 
                        0.5*WFC.cZZ*(ElectroWeakReweighter::coupling()->g_Lu(EWScale)*ElectroWeakReweighter::coupling()->g_Lu(EWScale) + 
                                     ElectroWeakReweighter::coupling()->g_Ru(EWScale)*ElectroWeakReweighter::coupling()->g_Ru(EWScale)) + 
                        ElectroWeakReweighter::coupling()->g_Lu(EWScale)*ElectroWeakReweighter::coupling()->g_Ru(EWScale)*WFC.bZZ;
   complex<double> U2 = ElectroWeakReweighter::coupling()->g_Ru(EWScale)*ElectroWeakReweighter::coupling()->g_Ru(EWScale)*(WFC.fFZZ-0.5*WFC.aZZ) -
                        0.5*WFC.cZZ*(ElectroWeakReweighter::coupling()->g_Lu(EWScale)*ElectroWeakReweighter::coupling()->g_Lu(EWScale) + 
                                     ElectroWeakReweighter::coupling()->g_Ru(EWScale)*ElectroWeakReweighter::coupling()->g_Ru(EWScale)) +
                        ElectroWeakReweighter::coupling()->g_Lu(EWScale)*ElectroWeakReweighter::coupling()->g_Ru(EWScale)*WFC.bZZ;
   complex<double> HtL = -0.5*y_t*y_t/(16.0*pi*pi)*(0.25-0.5*log(mH/EWScale)-0.5*lz+
                                                    0.5*WFC.atHH+0.5*WFC.atZZ+WFC.ctHH+WFC.ctZZ+
                                                    WFC.ctW0-WFC.btHH+WFC.btZZ);
   complex<double> HtR = HtL-0.5*y_t*y_t/(16.0*pi*pi)*(0.25-lw+WFC.atW0);
   complex<double> HbL = -0.5*y_t*y_t/(16.0*pi*pi)*(0.25-lw+WFC.at0W);
   
   complex<double> F_tL = (4.0/3.0*aS+4.0/9.0*aEM)/(4.0*pi)*(2.0*log(mT/EWScale)*log(mT/EWScale)-log(mT/EWScale)+
                                                             pi*pi/12.0+2.0) + HtL + 
                           aW/(4.0*pi)*0.5*W1 + aZ/(4.0*pi)*U1;
   
   complex<double> F_tR = (4.0/3.0*aS+4.0/9.0*aEM)/(4.0*pi)*(2.0*log(mT/EWScale)*log(mT/EWScale)-log(mT/EWScale)+
                                                             pi*pi/12.0+2.0) + HtR - 
                           aW/(4.0*pi)*0.25*WFC.cW0 + aZ/(4.0*pi)*U2;
   
   complex<double> F_bL = HbL + aW/(4.0*pi)*0.5*W2;
   
   Complex Dw = CollinearDw(s,EWScale);
   Complex Dz = CollinearDz(s,EWScale);

   complex<double> D_C_UL = ElectroWeakReweighter::coupling()->g_Lu(EWScale)*ElectroWeakReweighter::coupling()->g_Lu(EWScale)*Dz + 0.5*Dw;
   complex<double> D_C_DL = ElectroWeakReweighter::coupling()->g_Ld(EWScale)*ElectroWeakReweighter::coupling()->g_Ld(EWScale)*Dz + 0.5*Dw;
   
   complex<double> D_C_UR = ElectroWeakReweighter::coupling()->g_Ru(EWScale)*ElectroWeakReweighter::coupling()->g_Ru(EWScale)*Dz;
   complex<double> D_C_DR = ElectroWeakReweighter::coupling()->g_Rd(EWScale)*ElectroWeakReweighter::coupling()->g_Rd(EWScale)*Dz;
   
   complex<double> D_C_nuL = ElectroWeakReweighter::coupling()->g_Lnu(EWScale)*ElectroWeakReweighter::coupling()->g_Lnu(EWScale)*Dz + 0.5*Dw;
   complex<double> D_C_EL = ElectroWeakReweighter::coupling()->g_Le(EWScale)*ElectroWeakReweighter::coupling()->g_Le(EWScale)*Dz + 0.5*Dw;
   complex<double> D_C_ER = ElectroWeakReweighter::coupling()->g_Re(EWScale)*ElectroWeakReweighter::coupling()->g_Re(EWScale)*Dz;
   
   complex<double> D_C_WW = aW/(4.0*pi)*cW2*(F_Z+WFC.fsWZWZ) + 
     aW/(4.0*pi)*cW2*(F_W+WFC.fs1ZW) + aW/(4.0*pi)*sW2*(F_W+WFC.fs10) + 
     aW/(4.0*pi)*sW2*(2.0*lw*lw-
		      2.0*lw+pi*pi/12.0+2.0) + 0.5*WFC.RW;
   complex<double> D_C_WZ = aW/(4.0*pi)*2.0*(F_W+WFC.fsZW1) + 0.5*WFC.RZ + sqrt(sW2/cW2)*WFC.RAtoZ;
   complex<double> D_C_WA = aW/(4.0*pi)*2.0*(F_W+WFC.fs01) + 0.5*WFC.RA + sqrt(cW2/sW2)*WFC.RZtoA;
   complex<double> D_C_BZ = 0.5*WFC.RZ - sqrt(cW2/sW2)*WFC.RAtoZ;
   complex<double> D_C_BA = 0.5*WFC.RA - sqrt(sW2/cW2)*WFC.RZtoA;
   
   // The WFC.RW and WFC.RZ are used on purpose (instead of WFC.RPhi and WFC.RPhi3, respectively):
   complex<double> D_C_PhiW = aW/(4.0*pi)*0.25*(F_W+WFC.fs1HW) + 
     aW/(4.0*pi)*0.25*(F_W+WFC.fs1ZW) + aZ/(4.0*pi)*gPHI*gPHI*(F_Z+WFC.fsWZWZ) + 
     aW/(4.0*pi)*sW2*(2.0*lw*lw-2.0*lw+pi*pi/12.0+2.0) + 
     0.5*WFC.RW + WFC.EW;
   complex<double> D_C_PhiZ = aW/(4.0*pi)*0.5*(F_W+WFC.fsZW1) + aZ/(4.0*pi)*0.25*(F_Z+WFC.fs1HZ) + 0.5*WFC.RZ + WFC.EZ;
   complex<double> D_C_PhiH = aW/(4.0*pi)*0.5*(F_W+WFC.fsHW1) + aZ/(4.0*pi)*0.25*(F_Z+WFC.fsHZ1) + 0.5*WFC.RH;
   
   complex<double> D_C_GG = aS/(4.0*pi)*2.0/3.0*log(mT/EWScale);
   
   ULcollinearCorr_     = exp(D_C_UL);
   DLcollinearCorr_     = exp(D_C_DL);
   URcollinearCorr_     = exp(D_C_UR);
   DRcollinearCorr_     = exp(D_C_DR);
   
   tLcollinearCorr_     = exp(D_C_UL+F_tL);
   tRcollinearCorr_     = exp(D_C_UR+F_tR);
   bLcollinearCorr_     = exp(D_C_DL+F_bL);
   
   nuLcollinearCorr_    = exp(D_C_nuL);
   ELcollinearCorr_     = exp(D_C_EL);
   ERcollinearCorr_     = exp(D_C_ER);
   
   WtoWcollinearCorr_   = exp(D_C_WW);
   WtoZcollinearCorr_   = exp(D_C_WZ);
   WtoAcollinearCorr_   = exp(D_C_WA);
   BtoZcollinearCorr_   = exp(D_C_BZ);
   BtoAcollinearCorr_   = exp(D_C_BA);
   
   PhitoWcollinearCorr_ = exp(D_C_PhiW);
   PhitoZcollinearCorr_ = exp(D_C_PhiZ);
   PhitoHcollinearCorr_ = exp(D_C_PhiH);
   
   GcollinearCorr_      = exp(D_C_GG);
}

WaveFunctionCorrections CollinearSudakov::waveFunctionCorrections(Energy EWScale) {
  static const Complex I(0.,1.);
  using Constants::pi;
  double lZ = 2.0*log(ElectroWeakReweighter::coupling()->mZ()/EWScale);
  WaveFunctionCorrections WFC;
  // From Manohar, 2/12/12: (these assume mH=125, mZ=91.1876, mW=80.399)
  WFC.RW    = (0.8410283377963967 - 9.424777961271568*I) + 1.2785863646210789*lZ;
  WFC.RA    = (1.4835982362022198 + 1.855775680704845*pow(10.,-9)*I) - 0.27209907467584415*lZ;
  WFC.RZ    = (1.5114724841549798 - 9.926944419863688*I) + 1.0834802397165764*lZ;
  WFC.RAtoZ = (0.3667485032811715 - 2.2770907130064835*I) - 1.2994544609942593*lZ;
  WFC.RZtoA = -0.2095310079712942 + 0.8320191107808546*lZ;
  WFC.RH    = (12.229832449946716 - 1.7643103462419842*10.0*pow(10.,-12)*I) + 5.309998583664737*lZ;
  WFC.RPhi  = (5.569012418081201 + 1.5439133581417356*0.10*pow(10.,-9)*I) + 5.309998583664737*lZ;
  WFC.RPhi3 = (8.945333042265943 + 5.499309445612249*pow(10.,-12)*I) + 5.309998583664737*lZ;
  WFC.EW    = (3.967645734304811 + 4.712388980384717*I) + 2.238332625165702*lZ;
  WFC.EZ    = (5.916079892937651 + 4.96347220970469*I) + 2.1132591719740788*lZ;
  
  WFC.RW    *= ElectroWeakReweighter::coupling()->a2(EWScale)/(4.0*pi);
  WFC.RA    *= ElectroWeakReweighter::coupling()->a2(EWScale)/(4.0*pi);
  WFC.RZ    *= ElectroWeakReweighter::coupling()->a2(EWScale)/(4.0*pi);
  WFC.RAtoZ *= ElectroWeakReweighter::coupling()->a2(EWScale)/(4.0*pi);
  WFC.RZtoA *= ElectroWeakReweighter::coupling()->a2(EWScale)/(4.0*pi);
  WFC.RPhi  *= ElectroWeakReweighter::coupling()->a2(EWScale)/(4.0*pi);
  WFC.EW    *= ElectroWeakReweighter::coupling()->a2(EWScale)/(4.0*pi);
  WFC.EZ    *= ElectroWeakReweighter::coupling()->a2(EWScale)/(4.0*pi);
  WFC.RPhi3 *= ElectroWeakReweighter::coupling()->a2(EWScale)/(4.0*pi);
  WFC.RH    *= ElectroWeakReweighter::coupling()->a2(EWScale)/(4.0*pi);
  
  
  WFC.RW    = WFC.RW.real();
  WFC.RA    = WFC.RA.real();
  WFC.RZ    = WFC.RZ.real();
  WFC.RAtoZ = WFC.RAtoZ.real();
  WFC.RZtoA = WFC.RZtoA.real();
  WFC.RPhi  = WFC.RPhi.real();
  WFC.RPhi3 = WFC.RPhi3.real();
  WFC.RH    = WFC.RH.real();
  
  // Also from Manohar, 2/10/12
  WFC.fFW0 = -3.7946842553189453 - 4.709019671388589*I;
  WFC.fF0W = 3.8181790485161176;
  WFC.fFZZ = 1.364503250377989 + 0.*I;
  WFC.aHH  = -0.474396977740686 + 0.*I;
  WFC.aZZ  = -0.6806563210877914 + 0.*I;
  WFC.aW0  = 0.49036102506811907 + 1.9323351450971642*I;
  WFC.a0W  = -1.2184776671065072;
  WFC.bHH  = 1.234775745474991 + 0.*I;
  WFC.bZZ  = 1.7303526426747613 + 0.*I;
  WFC.cHH  = 0.33140608274387473 + 0.*I;
  WFC.cZZ  = 0.4961363208382017 + 0.*I;
  WFC.cW0  = -1.005299829777395 + 1.063048500002757*I;
  WFC.atHH = -0.237198488870343 + 0.*I;
  WFC.atZZ = -0.3403281605438957 + 0.*I;
  WFC.atW0 = 0.24518051253405954 + 0.9661675725485821*I;
  WFC.at0W = -0.6092388335532536;
  WFC.ctHH = 0.16570304137193737 + 0.*I;
  WFC.ctZZ = 0.24806816041910085 + 0.*I;
  WFC.ctW0 = -0.5026499148886975 + 0.5315242500013785*I;
  WFC.btHH = -0.30869393636874776 + 0.*I;
  WFC.btZZ = -0.4325881606686903 + 0.*I;
  
  WFC.fs10   = -2.289868133696459;
  WFC.fs1ZW  = 1.8627978596240622;
  WFC.fsWZWZ = 1.1866922529667008;
  WFC.fsZW1  = 1.0840307611156266;
  WFC.fs01   = 2.2898681336964595;
  WFC.fsHW1  = -0.32306745562682404;
  WFC.fsHZ1  = 0.4042992116255279;
  WFC.fs1HW  = 3.330954543719127;
  WFC.fs1HZ  = 2.69768201932412;
  return WFC;
}

Complex CollinearSudakov::CollinearDw(Energy2 s, Energy EWScale) {
  using Constants::pi;
  using GroupInvariants::PlusLog;
  double lw = log(ElectroWeakReweighter::coupling()->mW()/EWScale);
  //s = s/2.; // This should not be here... I think this is a discrepency with Sascha
  return ElectroWeakReweighter::coupling()->aW(EWScale)/(4.0*pi)*
    (4.0*lw*0.5*PlusLog(s/EWScale/EWScale) - 2.0*lw*lw -
     3.0*lw - 5.0/12.0*pi*pi + 9.0/4.0);
}

Complex CollinearSudakov::CollinearDz(Energy2 s, Energy EWScale) {
  using Constants::pi;
  using GroupInvariants::PlusLog;
  double lz = log(ElectroWeakReweighter::coupling()->mZ()/EWScale);
  return ElectroWeakReweighter::coupling()->aZ(EWScale)/(4.0*pi)*
    (4.0*lz*0.5*PlusLog(s/EWScale/EWScale) - 2.0*lz*lz -
     3.0*lz - 5.0/12.0*pi*pi + 9.0/4.0);
}

namespace {

void writeLine(ofstream & file, vector<Energy> x, vector<double> y,
	       string name,string title,
	       string colour, string style) {
  file << "# BEGIN HISTO1D "+name +"\n";
  file << "SmoothLine=1\n";
  file << "ErrorBars=0\n";
  file << "Path="+name+"\n";
  file << "Title="+title+"\n";
  file << "LineColor="+colour+"\n";
  file << "LineStyle="+style +"\n";
  file << "# xlow	 xhigh	 val	 errminus	 errplus\n";
  for(unsigned int ix=0;ix<x.size();++ix)
    file << (x[ix]-MeV)/GeV << "\t" << (x[ix]+MeV)/GeV << "\t"
	 << y[ix] << "\t" << 0. << "\t" << 0. << "\n";
  file << "# END HISTO1D\n";
}

}

void CollinearSudakov::makePlots() {
  vector<Energy> np;
  vector<double> uL,uR,tL,tR,dL,dR,bL,bR,nuL,eL,eR,Wgamma,Bgamma,g,WT,WL,WZT,BZT,ZL,H;
  Energy mZ = ElectroWeakReweighter::coupling()->mZ();
  for(Energy x=200.*GeV;x<5000.*GeV;x*=1.02) {
    Energy2 s(sqr(x));
    np.push_back(x);
    evaluateHighScale(x,mZ,s);
    evaluateMatching(mZ,s);
    uL    .push_back(real(highColQ_   * ULcollinearCorr_    ));
    uR    .push_back(real(highColU_   * URcollinearCorr_    ));
    tL    .push_back(real(highColQt_  * tLcollinearCorr_    ));
    tR    .push_back(real(highColtR_  * tRcollinearCorr_    ));
    dL    .push_back(real(highColQ_   * DLcollinearCorr_    ));
    dR    .push_back(real(highColD_   * DRcollinearCorr_    ));
    bL    .push_back(real(highColQt_  * bLcollinearCorr_    ));
    bR    .push_back(real(highColD_   * DRcollinearCorr_    ));
    nuL   .push_back(real(highColL_   * nuLcollinearCorr_   ));
    eL    .push_back(real(highColL_   * ELcollinearCorr_    ));
    eR    .push_back(real(highColE_   * ERcollinearCorr_    ));
    Wgamma.push_back(real(highColW_   * WtoAcollinearCorr_  ));
    Bgamma.push_back(real(highColB_   * BtoAcollinearCorr_  ));
    g     .push_back(real(highColG_   * GcollinearCorr_     ));
    WT    .push_back(real(highColW_   * WtoWcollinearCorr_  ));
    WL    .push_back(real(highColPhi_ * PhitoWcollinearCorr_));
    WZT   .push_back(real(highColW_   * WtoZcollinearCorr_  ));
    BZT   .push_back(real(highColB_   * BtoZcollinearCorr_  ));
    ZL    .push_back(real(highColPhi_ * PhitoZcollinearCorr_));
    H     .push_back(real(highColPhi_ * PhitoHcollinearCorr_));
  }
  ofstream fig1a("fig1a.dat");
  fig1a << "#BEGIN PLOT\n";
  fig1a << "LegendOnly=/uL /uR /tL /tR\n";
  fig1a << "DrawOnly=/uL /uR /tL /tR\n";
  fig1a << "Title=Figure 1a\n";
  fig1a << "RatioPlot=0\n";
  fig1a << "LogX=1\n";
  fig1a << "XLabel=$\\bar{n}\\cdot p$ [GeV]\n";
  fig1a << "YLabel=u, t\n";
  fig1a << "Legend=1\n";
  fig1a << "XMin=250.\n";
  fig1a << "YMin=0.7\n";
  fig1a << "YMax=1.05\n";
  fig1a << "# END PLOT\n";
  writeLine(fig1a,np,uL,"/uL","$u_L$","green","dotted"   );
  writeLine(fig1a,np,uR,"/uR","$u_R$","cyan" ,"solid"    );
  writeLine(fig1a,np,tL,"/tL","$t_L$","red"  ,"dashed"   );
  writeLine(fig1a,np,tR,"/tR","$t_R$","blue" ,"dotdashed");
  fig1a.close();
  ofstream fig1b("fig1b.dat");
  fig1b << "#BEGIN PLOT\n";
  fig1b << "LegendOnly=/dL /dR /bL /bR\n";
  fig1b << "DrawOnly=/dL /dR /bL /bR\n";
  fig1b << "Title=Figure 1b\n";
  fig1b << "RatioPlot=0\n";
  fig1b << "LogX=1\n";
  fig1b << "XLabel=$\\bar{n}\\cdot p$ [GeV]\n";
  fig1b << "YLabel=d, b\n";
  fig1b << "Legend=1\n";
  fig1b << "XMin=250.\n";
  fig1b << "YMin=0.7\n";
  fig1b << "YMax=1.05\n";
  fig1b << "# END PLOT\n";
  writeLine(fig1b,np,dL,"/dL","$d_L$","green","dotted"   );
  writeLine(fig1b,np,dR,"/dR","$d_R$","cyan" ,"solid"    );
  writeLine(fig1b,np,bL,"/bL","$b_L$","red"  ,"dashed"   );
  writeLine(fig1b,np,bR,"/bR","$b_R$","blue" ,"dotdashed");
  fig1b.close();
  ofstream fig2("fig2.dat");
  fig2 << "#BEGIN PLOT\n";
  fig2 << "LegendOnly=/uL /uR /dL /dR\n";
  fig2 << "DrawOnly=/uL /uR /dL /dR\n";
  fig2 << "Title=Figure 2\n";
  fig2 << "RatioPlot=0\n";
  fig2 << "LogX=1\n";
  fig2 << "XLabel=$\\bar{n}\\cdot p$ [GeV]\n";
  fig2 << "YLabel=u, d\n";
  fig2 << "Legend=1\n";
  fig2 << "XMin=250.\n";
  fig2 << "YMin=0.7\n";
  fig2 << "YMax=1.05\n";
  fig2 << "# END PLOT\n";
  writeLine(fig2,np,uL,"/uL","$u_L$","green","dotted"   );
  writeLine(fig2,np,uR,"/uR","$u_R$","cyan" ,"solid"    );
  writeLine(fig2,np,dL,"/dL","$d_L$","red"  ,"dashed"   );
  writeLine(fig2,np,dR,"/dR","$d_R$","blue" ,"dotdashed");
  fig2.close();
  ofstream fig3("fig3.dat");
  fig3 << "#BEGIN PLOT\n";
  fig3 << "LegendOnly=/nuL /eL /eR\n";
  fig3 << "DrawOnly=/nuL /eL /eR\n";
  fig3 << "Title=Figure 3\n";
  fig3 << "RatioPlot=0\n";
  fig3 << "LogX=1\n";
  fig3 << "XLabel=$\\bar{n}\\cdot p$ [GeV]\n";
  fig3 << "YLabel=$\\nu$, $e$\n";
  fig3 << "Legend=1\n";
  fig3 << "XMin=250.\n";
  fig3 << "YMin=0.9\n";
  fig3 << "YMax=1.05\n";
  fig3 << "# END PLOT\n";
  writeLine(fig3,np,nuL,"/nuL","$\\nu_L$","blue","dashed");
  writeLine(fig3,np, eL,"/eL" ,"$e_L$"   ,"red" ,"dotted");
  writeLine(fig3,np, eR,"/eR" ,"$e_R$"   ,"red" ,"solid" );
  fig3.close();
  ofstream fig5("fig5.dat");
  fig5 << "#BEGIN PLOT\n";
  fig5 << "LegendOnly=/g /Wgamma /Bgamma\n";
  fig5 << "DrawOnly=/g /Wgamma /Bgamma\n";
  fig5 << "Title=Figure 5\n";
  fig5 << "RatioPlot=0\n";
  fig5 << "LogX=1\n";
  fig5 << "XLabel=$\\bar{n}\\cdot p$ [GeV]\n";
  fig5 << "YLabel=$\\gamma$, g\n";
  fig5 << "Legend=1\n";
  fig5 << "XMin=250.\n";
  fig5 << "YMin=0.7\n";
  fig5 << "YMax=1.05\n";
  fig5 << "# END PLOT\n";
  writeLine(fig5,np,g     ,"/g"     ,"$g$"           ,"blue","dashed");
  writeLine(fig5,np,Wgamma,"/Wgamma","$W\\to\\gamma$","red" ,"solid" );
  writeLine(fig5,np,Bgamma,"/Bgamma","$B\\to\\gamma$","red" ,"dotted");
  fig5.close();
  ofstream fig6a("fig6a.dat");
  fig6a << "#BEGIN PLOT\n";
  fig6a << "LegendOnly=/WT /WL\n";
  fig6a << "DrawOnly=/WT /WL\n";
  fig6a << "Title=Figure 6a\n";
  fig6a << "RatioPlot=0\n";
  fig6a << "LogX=1\n";
  fig6a << "XLabel=$\\bar{n}\\cdot p$ [GeV]\n";
  fig6a << "YLabel=$Z_L$, $Z_T$, $H$\n";
  fig6a << "Legend=1\n";
  fig6a << "XMin=250.\n";
  fig6a << "YMin=0.7\n";
  fig6a << "YMax=1.05\n";
  fig6a << "# END PLOT\n";
  writeLine(fig6a,np,WT,"/WT","$W_T$","red"  ,"solid");
  writeLine(fig6a,np,WL,"/WL","$W_L$","blue" ,"dashed" );
  fig6a.close();


  ofstream fig6b("fig6b.dat");
  fig6b << "#BEGIN PLOT\n";
  fig6b << "LegendOnly=/WZT /BZT /ZL /H\n";
  fig6b << "DrawOnly=/WZT /BZT /ZL /H\n";
  fig6b << "Title=Figure 6b\n";
  fig6b << "RatioPlot=0\n";
  fig6b << "LogX=1\n";
  fig6b << "XLabel=$\\bar{n}\\cdot p$ [GeV]\n";
  fig6b << "YLabel=d, b\n";
  fig6b << "Legend=1\n";
  fig6b << "XMin=250.\n";
  fig6b << "YMin=0.7\n";
  fig6b << "YMax=1.05\n";
  fig6b << "# END PLOT\n";
  writeLine(fig6b,np,WZT,"/WZT","$W\\to Z_T$","red"  ,"solid"    );
  writeLine(fig6b,np,BZT,"/BZT","$B\\to Z_T$","red"  ,"dotted"   );
  writeLine(fig6b,np,ZL ,"/ZL ","$Z_L$"      ,"blue" ,"dashed"   );
  writeLine(fig6b,np,H  ,"/H  ","$H$"        ,"green","dotdashed");
  fig6b.close();


  np.clear();
  vector<double> e30,e50,q30,q50,g30,g50;
  for(Energy x=200.*GeV;x<5000.*GeV;x*=1.02) {
    Energy2 s(sqr(x));
    np.push_back(x);
    evaluateLowScale(mZ,30.*GeV,s);
    e30.push_back(real(lowColE_));
    q30.push_back(real(lowColU_));
    g30.push_back(real(lowColG_));
    evaluateLowScale(mZ,50.*GeV,s);
    e50.push_back(real(lowColE_));
    q50.push_back(real(lowColU_));
    g50.push_back(real(lowColG_));
  }
  ofstream fig4a("fig4a.dat");
  fig4a << "#BEGIN PLOT\n";
  fig4a << "LegendOnly=/e30 /e50\n";
  fig4a << "DrawOnly=/e30 /e50\n";
  fig4a << "Title=Figure 4a\n";
  fig4a << "RatioPlot=0\n";
  fig4a << "LogX=1\n";
  fig4a << "XLabel=$\\bar{n}\\cdot p$ [GeV]\n";
  fig4a << "YLabel=e\n";
  fig4a << "Legend=1\n";
  fig4a << "XMin=250.\n";
  fig4a << "YMin=0.7\n";
  fig4a << "YMax=1.05\n";
  fig4a << "# END PLOT\n";
  writeLine(fig4a,np,e30,"/e30","e $(\\mu_f=30\\,\\mathrm{GeV})$","red" ,"solid" );
  writeLine(fig4a,np,e50,"/e50","e $(\\mu_f=50\\,\\mathrm{GeV})$","blue","dashed");
  fig4a.close();
  ofstream fig4b("fig4b.dat");
  fig4b << "#BEGIN PLOT\n";
  fig4b << "LegendOnly=/q30 /q50 /g30 /g50\n";
  fig4b << "DrawOnly=/q30 /q50 /g30 /g50\n";
  fig4b << "Title=Figure 4a\n";
  fig4b << "RatioPlot=0\n";
  fig4b << "LogX=1\n";
  fig4b << "XLabel=$\\bar{n}\\cdot p$ [GeV]\n";
  fig4b << "YLabel=e\n";
  fig4b << "Legend=1\n";
  fig4b << "XMin=250.\n";
  fig4b << "YMin=0.5\n";
  fig4b << "YMax=1.05\n";
  fig4b << "# END PLOT\n";
  writeLine(fig4b,np,q30,"/q30","q $(\\mu_f=30\\,\\mathrm{GeV})$","red" ,"solid" );
  writeLine(fig4b,np,q50,"/q50","q $(\\mu_f=50\\,\\mathrm{GeV})$","blue","dashed");
  writeLine(fig4b,np,g30,"/g30","g $(\\mu_f=30\\,\\mathrm{GeV})$","green" ,"dotted" );
  writeLine(fig4b,np,g50,"/g50","g $(\\mu_f=50\\,\\mathrm{GeV})$","blue","dotdashed");
  fig4b.close();
}
