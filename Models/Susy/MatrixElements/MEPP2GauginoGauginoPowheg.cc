// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2GauginoGauginoPowheg class.
//

#include "MEPP2GauginoGauginoPowheg.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "SusyLoopIntegral.h"

using namespace Herwig;

MEPP2GauginoGauginoPowheg::MEPP2GauginoGauginoPowheg() {}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeAbstractNoPIOClass<MEPP2GauginoGauginoPowheg,NLODrellYanBase>
describeHerwigMEPP2GauginoGauginoPowheg("Herwig::MEPP2GauginoGauginoPowheg", 
					"HwSusy.so HwSusyNLO.so");

void MEPP2GauginoGauginoPowheg::Init() {

  static ClassDocumentation<MEPP2GauginoGauginoPowheg> documentation
    ("There is no documentation for the MEPP2GauginoGauginoPowheg class");

}

double MEPP2GauginoGauginoPowheg::
finiteVirtual(Energy ms, Energy2 mb2,
	      vector<Complex> Cl, vector<Complex> Cr,
	      vector<double > Cs, vector<Complex> Ct) const {
  // cut-off parameter
  Energy2 eps  =  0.1*MeV2;
  Energy  eps2 = 1e-5*MeV;
  // outgoing masses
  Energy m1 = meMomenta()[2].mass();
  Energy m2 = meMomenta()[3].mass();
  Energy2 m12 = sqr(m1), m22 = sqr(m2);
  // average of outgoing masses
  Energy2 mav2 = 0.25*sqr(m1+m2);
  // renormalisation/factorization scale
  Energy2 scale2 = scale();
  // gluino mass
  Energy mg  = getParticleData(ParticleID::SUSY_g)->mass();
  Energy2 mg2(sqr(mg));
  Energy2 ms2(sqr(ms));
  // mandelstam variables
  Energy2 t1 = tHat()-m12, t2 = tHat()-m22;
  Energy2 u1 = uHat()-m12, u2 = uHat()-m22;
  Energy2 ts = tHat()-ms2, us = uHat()-ms2;
  Energy2 tg = tHat()-mg2, ug = uHat()-mg2;
  // mandelstam variables
  Energy2 sz = sHat()-mb2;
  Energy4 kaellen = sqr(sHat()) + pow<4,1>(m1) + pow<4,1>(m2)
    -2.*( sHat()*m12 + sHat()*m22 + sqr(m1*m2 ));
  // vector coupling
  vector<Complex> Cv(4,0.);
  for(unsigned int ix=0;ix<4;++ix) Cv[ix] = conj(Cl[ix])/Cl[ix];
  // now for the LO kinematic pieces
  Complex QBB_ss = 
    16.*Cs[0]*double((t1*t2+u1*u2)/sqr(sHat())) + 
    32.*Cs[1]*double(m1*m2/sHat()) +
    32.*Cs[2]*double((u1*u2- t1*t2)/sHat()/sz);
  Complex QBB_st = 
    16.*double(t1*t2/sHat()/ts) * ( Cl[0]*conj(Cl[1])*Ct[0] + Cr[0]*conj(Cr[1])*Ct[2]) +
    16.*double(m1*m2/ts       ) * ( Cl[0]*conj(Cl[1])*Ct[1] + Cr[0]*conj(Cr[1])*Ct[3]);
  Complex QBB_su =
    -16.*u1*u2/sHat()/us * ( Cl[3]*conj(Cl[2])*Ct[1] + Cr[3]*conj(Cr[2])*Ct[3])
    -16.*m1*m2/us        * ( Cl[3]*conj(Cl[2])*Ct[0] + Cr[3]*conj(Cr[2])*Ct[2]);
  Complex QBB_txs =
    16.*t1*t2/sHat()/ts  * ( Cl[1]*conj(Cl[0]*Ct[0]) + Cr[1]*conj(Cr[0]*Ct[2]));
  Complex QBB_txt =  
    32.*double(t1*t2/sqr(ts))* ( norm(Cl[0])*norm(Cl[1]) + norm(Cr[0])*norm(Cr[1]) );
  Complex QBB_tmu = 
    -32.*m1*m2*sHat()/ts/us  * ( Cl[1]*Cl[3]*conj(Cl[0]*Cl[2]) + 
				 Cr[1]*Cr[3]*conj(Cr[0]*Cr[2]) );
  Complex QBB_tms = 
    16.*m1*m2/ts             * ( Cl[1]*conj(Cl[0]*Ct[1]) + 
				 Cr[1]*conj(Cr[0]*Ct[3]) );
  Complex QBB_uxs =         
    -16.*u1*u2/sHat()/us     * ( Cl[2]*conj(Cl[3]*Ct[1]) + 
				 Cr[2]*conj(Cr[3]*Ct[3]) );
  Complex QBB_uxu =         
    32.*double(u1*u2/sqr(us))* ( norm(Cl[2])*norm(Cl[3]) + 
				 norm(Cr[2])*norm(Cr[3]) );
  Complex QBB_umt =         
    -32.*m1*m2*sHat()/ts/us  * ( Cl[0]*Cl[2]*conj(Cl[1]*Cl[3]) + 
				 Cr[0]*Cr[2]*conj(Cr[1]*Cr[3]) );
  Complex QBB_ums =        
    -16.*m1*m2/us            * ( Cl[2]*conj(Cl[3])*conj(Ct[0]) + 
				 Cr[2]*conj(Cr[3])*conj(Ct[2]) );
  double b_s  = 3.*real( QBB_ss  + QBB_st   + QBB_su );
  double b_tx = 3.*real( QBB_txs + QBB_txt );
  double b_tm = 3.*real( QBB_tmu + QBB_tms );
  double b_ux = 3.*real( QBB_uxs + QBB_uxu );
  double b_um = 3.*real( QBB_umt + QBB_ums );
  double b_t  = b_tx  + b_tm;
  double b_u  = b_ux  + b_um; 


//   cerr << "testing in virtual\n";
//   cerr << "mandelstam" << sHat()/GeV2<< " " <<tHat()/GeV2<< " " <<uHat()/GeV2 << "\n";
//   cerr << "testing chi mass" << m1/GeV<< " " <<m2/GeV<< "\n";
//   cerr << "tesitng mg,ms" << mg/GeV<< " " <<ms/GeV<< "\n";
//   cerr << "testing ts,us" << ts/GeV2<< " " <<us/GeV2<< "\n";
//   cerr << "testing tg,ug" << tg/GeV2<< " " <<ug/GeV2<< "\n";
//   cerr << "testing t1,u1" << t1/GeV2<< " " <<u1/GeV2<< "\n";
//   cerr << "testing t2,u2" << t2/GeV2<< " " <<u2/GeV2<< "\n";

//   cerr << "testing pieces A" << b_s << " " << b_t << " " << b_u<< "\n";
//   cerr << "testing pieces B" << b_tx << " " << b_tm<< "\n";
//   cerr << "testing pieces C" << b_ux << " " << b_um<< "\n";
//   cerr << "testing pieces D" << b_tx << " " << b_tm<< "\n";
//   cerr << "testing pieces E" << b_ux << " " << b_um<< "\n";

 
  // test of the LO
  // //   double test = 0.5*(b_s + b_t + b_u)/16./sqr(Constants::pi);
  //   double test = 0.5*sqr(ee*ee)*(b_s+b_t+b_u)/3./16./sqr(Constants::pi)/12.;
  
  //   cerr << "testing " << test << " " << loWeight() << " "
  //        << test/loWeight() << "\n";
  // now for the finite piece
  Complex values[26];
  // LO like piece
  values[0] =  2.*( b_s + b_t + b_u ) * 
    (-sqr(Constants::pi)/6.- 3.*log(sHat()/mav2)
     +sqr(log(sHat()/mav2))-1. 
     -SusyLoopIntegral::B0P(ZERO,mg,ms,scale2)* ( mg2 - ms2 ));
//   cerr << setprecision(10);
//   cerr << "testing A " << values[0] << "\n";
  // the rest of it
  double B0tg=SusyLoopIntegral::B0 (tHat(),mg,ZERO,scale2);
  values[1] += -4.*B0tg*b_t*(Complex((Cv[1]*m2/t2+Cv[0]*m1/t1)*mg)+double(tg/ts))
    -4.*B0tg*b_um*double(us/sHat())
    +4.*B0tg*b_ux*double((sHat()/u1*(m12/t1+m22/t2)+m22/u1-1.)*us/u2);
  // cerr << "testing B " << values[1] << "\n";
  double B0ug=SusyLoopIntegral::B0 (uHat(),mg,ZERO,scale2);
  values[2] += -4.*B0ug*b_u*(Complex((Cv[2]*m1/u1+Cv[3]*m2/u2)*mg)+double(ug/us))
    -4.*B0ug*b_tm*double(ts/sHat())
    +4.*B0ug*b_tx*double((sHat()/t1*(m12/u1+m22/u2)+m22/t1-1.)*ts/t2);
  // cerr << "testing C " << values[2] << "\n";
  values[3] += -4.*SusyLoopIntegral::B0(ms2,mg,ZERO,scale2)*(b_t/ts+b_u/us)*(mg2-ms2);
  // cerr << "testing D " << values[3] << "\n";
  values[4] +=  4.*SusyLoopIntegral::B0(tHat(),ms,ZERO,scale2)*
    (-b_t*(1.+m12/t1-(tHat()+ms2)/ts+m22/t2)
     -b_tm/sHat()*ts+b_tx*tHat()*ts/t1/t2);
  // cerr << "testing E " << values[4] << "\n";
  values[5] +=  4.*SusyLoopIntegral::B0(uHat(),ms,ZERO,scale2)*
    (-b_u*(1.+m12/u1-(uHat()+ms2)/us+m22/u2)
     -b_um/sHat()*us+b_ux*uHat()*us/u1/u2);
  // cerr << "testing F " << values[5] << "\n";
  values[6] += -8.*ms2*SusyLoopIntegral::B0(ms2,ms,ZERO,scale2)*(b_t/ts + b_u/us);
  // cerr << "testing G " << values[6] << "\n";
  double B012s=SusyLoopIntegral::B0(m12,ms,ZERO,scale2);
  values[7] +=  4.*B012s*m1*mg*(Cv[2]*b_u/u1+Cv[0]*b_t/t1);
  values[7] +=  4.*B012s/kaellen*b_tm*ts/sHat()   *(sHat()*tHat()-m12*u2-m22*t1);
  values[7] +=  4.*B012s/kaellen*b_tx*m12*ts/t1/t2*(sHat()*tHat()-m12*t2-m22*u1);
  values[7] += -4.*B012s/kaellen*b_um*us/sHat()   *(sHat()*tHat()-m12*u2-m22*t1);
  values[7] +=  4.*B012s/kaellen*b_ux*m12*us/u1/u2*(sHat()*uHat()-m12*u2-m22*t1);
  values[7] +=  4.*B012s/kaellen*b_um*us/sHat()   *(sHat()*tHat()-m12*u2-m22*t1);
  values[7] += -4.*B012s/kaellen*b_ux*us*m12/u1/u2*(sHat()*uHat()-m12*u2-m22*t1);
  values[7] += -4.*B012s/kaellen*b_tm/sHat()*ts*   (sHat()*tHat()-m12*u2-m22*t1);
  values[7] += -4.*B012s/kaellen*b_tx*ts*m12/t1/t2*(sHat()*tHat()-m12*t2-m22*u1);
  values[7] +=  4.*B012s*m12*(b_t/t1+b_u/u1);
  values[7] +=  4.*B012s/sHat()*(b_tm*ts+b_um*us);
  values[7] +=  4.*B012s*m12/u1/t1*(b_ux*us+b_tx*ts);
  // cerr << "testing H " << values[7] << "\n";
  double B022s = SusyLoopIntegral::B0(m22,ms,ZERO,scale2);
  values[8] +=  4.*B022s*m2*mg*(Cv[1]*b_t/t2+Cv[3]*b_u/u2);
  values[8] += -4.*B022s/kaellen*b_tm*ts/sHat()*(uHat()*sHat()-m12*u2-m22*t1);
  values[8] +=  4.*B022s/kaellen*b_tx*m22*ts/t1/t2*(sHat()*tHat()-m12*u2-m22*t1);
  values[8] +=  4.*B022s/kaellen*b_um*us/sHat()   *(sHat()*uHat()-m12*u2-m22*t1);
  values[8] +=  4.*B022s/kaellen*b_ux*us*m22/u1/u2*(sHat()*uHat()-m12*t2-m22*u1);
  values[8] += -4.*B022s/kaellen*b_um*us/sHat()*   (uHat()*sHat()-m12*u2-m22*t1);
  values[8] += -4.*B022s/kaellen*b_ux*m22*us/u1/u2*(sHat()*uHat()-m12*t2-m22*u1);
  values[8] += +4.*B022s/kaellen*b_tm*ts/sHat()*   (sHat()*uHat()-m12*u2-m22*t1);
  values[8] += -4.*B022s/kaellen*b_tx*m22*ts/t1/t2*(sHat()*tHat()-m12*u2-m22*t1);
  values[8] +=  4.*B022s*m22*(b_t/t2+b_u/u2);
  values[8] +=  4.*B022s*b_um*us/sHat();
  values[8] +=  4.*B022s*b_ux*m22*us/t2/u2;
  values[8] +=  4.*B022s*b_tm *ts/sHat();
  values[8] +=  4.*B022s*b_tx *m22/t2/u2*ts; 
  // cerr << "testing I " << values[8] << "\n";
  double B0sss = SusyLoopIntegral::B0(sHat(),ms,ms,scale2);
  values[9] +=  4.*B0sss/kaellen*b_um * (uHat()-tHat())*us;
  values[9] +=  4.*B0sss/kaellen*b_ux*us/u1/u2*(2.*m12*m22*sHat()-sqr(m12-m22)*uHat()+(m12+m22)*sHat()*uHat());
  values[9] +=  4.*B0sss/kaellen*b_tm*ts*(tHat()-uHat());
  values[9] +=  4.*B0sss/kaellen*b_tx*ts/t1/t2*(2.*m12*m22*sHat()-sqr(m12-m22)*tHat()+(m12+m22)*sHat()*tHat());
  values[9] +=  2.*B0sss*b_s*(1.+2.*(mg2-ms2)/sHat());
  values[9] +=  4.*B0sss*(b_ux*uHat()*us/u1/u2+b_tx*tHat()*ts/t1/t2);
  // cerr << "testing J " << values[9] << "\n";
  values[10] += -2.*SusyLoopIntegral::B0(ZERO,mg,ms,scale2)*(b_s*(1.+2.*(mg2-ms2)/sHat())+b_t+b_u);
  // cerr << "testing K " << values[10] << "\n";
  double B0s = SusyLoopIntegral::B0(sHat(),ZERO,ZERO,scale2);
  values[11] +=  4.*B0s/kaellen*b_tm*ts*(uHat()-tHat());
  values[11] += -4.*B0s/kaellen*b_tx*ts/t1/t2*(2.*m12*m22*sHat()-sqr(m12-m22)*tHat()+(m12+m22)*sHat()*tHat());
  values[11] +=  4.*B0s/kaellen*b_um*us*(tHat()-uHat());
  values[11] += -4.*B0s/kaellen*b_ux*us/u1/u2*(2.*m12*m22*sHat()-sqr(m12-m22)*uHat()+(m12+m22)*sHat()*uHat());
  values[11] += -6.*B0s*b_s;
  values[11] += -4.*B0s*(b_tx*tHat()*ts/t1/t2+b_ux*uHat()*us/u1/u2);
  // cerr << "testing L " << values[11] << "\n";
  complex<InvEnergy2> test = SusyLoopIntegral::C0(m12,eps,tHat(),eps2,ms  ,mg  );
  double C01tsg = ( sHat()*SusyLoopIntegral::C0(m12,eps,tHat(),eps2,ms  ,mg  )).real();
  values[12] +=  4.*C01tsg*Cv[0]*b_t/sHat()/t1*m1*mg*(mg2-ms2-t1);
  values[12] +=  2.*C01tsg*b_um*us/sqr(sHat())*(2.*(mg2-ms2)-t1);
  values[12] +=  2.*C01tsg*b_ux*us/sHat()/u1/u2/t1*(-(mg2-ms2)*(t1*(tHat()+m12)+2*sHat()*m12)+(ms2-sHat())*t1*t1);
  // cerr << "testing M " << values[12] << "\n";
  double C02tsg = ( sHat()*SusyLoopIntegral::C0(m22,eps,tHat(),eps2,ms  ,mg  ) ).real();
  values[13] +=  4.*C02tsg*Cv[1]*b_t*m2*mg/sHat()/t2*(mg2-ms2-t2);
  values[13] +=  2.*C02tsg*b_um*us/sqr(sHat())*(2.*(mg2-ms2)-t2);
  values[13] +=  2.*C02tsg*b_ux*us/sHat()/u1/u2/t2*(-(mg2-ms2)*(t2*(tHat()+m22)+2*sHat()*m22)+(ms2-sHat())*t2*t2);
  // cerr << "testing N " << values[13] << "\n";
  double C01usg = ( sHat()*SusyLoopIntegral::C0(m12,eps,uHat(),eps2,ms  ,mg  ) ).real();
  values[14] +=  4.*C01usg*Cv[2]*b_u*m1*mg/sHat()/u1*(mg2-ms2-u1);
  values[14] +=  2.*C01usg*b_tm *ts/sqr(sHat())*(2.*(mg2-ms2)-u1);
  values[14] +=  2.*C01usg*b_tx*ts/t1/t2/sHat()/u1*(-(mg2-ms2)*((uHat()+m12)*u1+2*sHat()*m12)+(ms2-sHat())*u1*u1);
  // cerr << "testing O " << values[14] << "\n";
  double C02usg = ( sHat()*SusyLoopIntegral::C0(m22,eps,uHat(),eps2,ms  ,mg  ) ).real();
  values[15] +=  4.*C02usg*Cv[3]*b_u*m2*mg/sHat()/u2*(mg2-ms2-u2);
  values[15] +=  2.*C02usg*b_tm*ts/sqr(sHat())*(2.*(mg2-ms2)-u2);
  values[15] +=  2.*C02usg*b_tx*ts/t1/t2/sHat()/u2*(-(mg2-ms2)*((uHat()+m22)*u2+2*sHat()*m22)+(ms2-sHat())*u2*u2);
  // cerr << "testing P " << values[15] << "\n";
  double C0Dt1 = SusyLoopIntegral::C0div(tHat(),m1,m2,ms,scale2);
  double C0Dt2 = SusyLoopIntegral::C0div(tHat(),m2,m1,ms,scale2);
  values[16] +=  4.*b_t*(C0Dt1/t1*(ms2-m12)+C0Dt2/t2*(ms2-m22));
  values[16] +=  4.*b_tm*ts*(C0Dt1/t1+C0Dt2/t2 );
  values[16] +=  2.*b_tx*ts/t1/t2*(C0Dt1*(tHat()-ms2+m12-m22)+C0Dt2*(tHat()-ms2-m12+m22));
  // cerr << "testing Q " << values[16] << "\n";
  double C0Du1 = SusyLoopIntegral::C0div(uHat(),m1,m2,ms,scale2);
  double C0Du2 = SusyLoopIntegral::C0div(uHat(),m2,m1,ms,scale2);
  values[17] += -4.*b_u*(C0Du1/u1*(m12-ms2)+C0Du2/u2*(m22-ms2));
  values[17] +=  4.*b_um*us*(C0Du1/u1+C0Du2/u2);
  values[17] +=  2.*b_ux*us/u1/u2*(C0Du1*(uHat()-ms2+m12-m22)+C0Du2*(uHat()-ms2-m12+m22));
  // cerr << "testing R " << values[17] << "\n";
  double C000s =( sHat()*SusyLoopIntegral::C0(eps,eps,sHat(),ms  ,mg  ,ms  ) ).real();
  values[18] +=  4.*C000s*b_s * ( mg2*sHat() + sqr(mg2-ms2) )/sqr(sHat());
  values[18] += -2.*C000s*b_um/sHat()*us;
  values[18] +=  2.*C000s*b_ux/u1/u2*us*(2*ms2-mg2-sHat());
  values[18] += -2.*C000s*b_tm/sHat()*ts;
  values[18] +=  2.*C000s*b_tx/t1/t2*ts*(2*ms2-mg2-sHat());
  // cerr << "testing S " << values[18] << "\n";
  double C0div = 0.5*sqr(log(sHat()/scale2))-7.*sqr(Constants::pi)/12.;
  values[19] += -4.*C0div*b_s;
  values[19] +=  2.*C0div*b_tx*ts/t1/t2*(m12+m22-ms2-tHat());
  values[19] +=  2.*C0div*b_ux*us/u1/u2*(m12+m22-ms2-uHat());
  // cerr << "testing T " << values[19] << "\n";
  double C0120s0=( sHat()*SusyLoopIntegral::C0(m12,m22,sHat(),eps2,ms  ,eps2) ).real();
  Energy2 s= sHat(), u=uHat(), t=tHat();


  values[20] += +4.*C0120s0/kaellen*b_tm*ts/sHat()*(2*m12*m22+ms2*(uHat()-tHat())-sqr(uHat())-uHat()*tHat());
  values[20] += +C0120s0/kaellen*b_tx*ts*(  - 4*m12*m22*ms2/t1/t2 - 3*m12*m22*s/t1/t2 - 4*m12*m22/t1 + 4*m12*m22/t2 - 11*m12*m22*m22/t1/t2 - m12*ms2*s/t1/t2 + 2*m12*ms2/t1 + 2*m12*ms2/t2 + m12*s/t1 + 3*m12*s/t2 - 2*m12*u/t1 - 2*m12 + 13*m12*m12*m22/t1/t2 + 2*m12*m12*ms2/t1/t2 + 3*m12*m12*s/t1/t2 + 2*m12*m12/t1 - 5*m12*m12/t2 - 5*m12*m12*m12/t1/t2 - 3*m22*ms2*s/t1/t2 + m22*ms2/t1 + 2*m22*ms2/t2 - 3*m22*s/t1 + 4*m22*s/t2 - m22*u/t1 - 2*m22 + 2*m22*m22*ms2/t1/t2 - 4*m22*m22*s/t1/t2 + 2*m22*m22/t1 - 3*m22*m22/t2 + 3*m22*m22*m22/t1/t2 );
  values[20] += +C0120s0/kaellen*b_tx*ts/t1/t2*(ms2*(t2*m22+sHat()*(m12-m22)-4*sHat()*t2)
					 +3*sHat()*(-u1*m22-sHat()*tHat())-uHat()*m22*t2);
  values[20] += +C0120s0/kaellen*b_um *us* ( 4*m12*ms2/sHat() + 4*m12/sHat()*u + 8*m12 - 4*m12*m12/sHat() + 4*m22*ms2/sHat() + 4*m22/sHat()*u + 8*m22 - 4*m22*m22/sHat() - 8*ms2/sHat()*u - 4*ms2 - 4*s - 4*u );
  values[20] += +C0120s0/kaellen*b_ux *us* (  - 4*m12*m22*ms2/u1/u2 - 9*m12*m22*s/u1/u2 - 3*m12*m22/u1 - 7*m12*m22/u2 + 7*m12*m22*m22/u1/u2 - m12*ms2*s/u1/u2 + 3*m12*ms2/u1 + 2*m12*ms2/u2 - 3*m12*s/u1 + 2*m12*s/u2 + m12*u/u2 - m12 - 5*m12*m12*m22/u1/u2 + 2*m12*m12*ms2/u1/u2 + m12*m12*m12/u1/u2 - 3*m22*ms2*s/u1/u2 + 2*m22*ms2/u1 + m22*ms2/u2 + 3*m22*s/u1 - 3*m22*s/u2 + m22*u/u1 - 2*m22 + 2*m22*m22*ms2/u1/u2 + 5*m22*m22*s/u1/u2 - 2*m22*m22/u1 + m22*m22/u2 - 3*m22*m22*m22/u1/u2 - 3*ms2*s/u1 - ms2*s/u2 - ms2*u/u1 );
  values[20] += +C0120s0/kaellen*b_ux *us* ( ms2*u/u2 + s*u/u1 - s + s*s/u1 - s*s/u2 - u + u*u/u2 );
  values[20] += -4.*C0120s0*b_tm*ts/sHat();
  values[20] += +C0120s0*b_tx*ts*(  - m12/sHat()/t2 + 2*m12/t1/t2 - 4*m22/t1/t2 + 2*ms2/sHat()/t1 + 2*ms2/sHat()/t2 + 2*ms2/t1/t2 - 2/sHat()*u/t1 - 1/s*u/t2 - 3/sHat() - 1./t1 + 3/t2 );
  values[20] += -4.*C0120s0*b_um*us/sHat();
  values[20] +=     C0120s0*b_ux*us/sHat()/u1/u2*(2.*tHat()*(m12+m22)-4.*m12*m22+2.*ms2*(uHat()-tHat())
					   +sHat()*(2.*uHat()-m12+m22));
  // cerr << "testing U " << values[20] << "\n";
  double C012s0s = ( sHat()*SusyLoopIntegral::C0(m12,m22,sHat(),ms  ,eps2,ms  ) ).real();
  values[21] += + C012s0s/kaellen*b_um*us*( 8*m12*m22/s + 4*m12*ms2/s - 4*m12/s*u + 4*m22*ms2/s - 4*m22/s*u - 8*ms2/s*u - 4*ms2 + 4*u );
  values[21] += + C012s0s/kaellen*b_ux*us*(  - 4*m12*m22*ms2/u1/u2 + 6*m12*m22*s/u1/u2 + 6*m12*m22/u1 - 2*m12*m22/u2 + 10*m12*m22*m22/u1/u2 - 2*m12*ms2*s/u1/u2 + 2*m12*ms2/u1 + 2*m12*ms2/u2 - 6*m12*s/u1 - 8*m12*s/u2 - 2*m12 - 8*m12*m12*m22/u1/u2 + 2*m12*m12*ms2/u1/u2 - 8*m12*m12*s/u1/u2 + 2*m12*m12/u2 + 2*m12*m12*m12/u1/u2 - 2*m22*ms2*s/u1/u2 + 2*m22*ms2/u1 + 2*m22*ms2/u2 - 2*m22*s/u1 - 6*m22*s/u2 - 4*m22*u/u1 + 4*m22 + 2*m22*m22*ms2/u1/u2 + 6*m22*m22*s/u1/u2 + 4*m22*m22/u2 - 4*m22*m22*m22/u1/u2 - 2*ms2*s/u1 - 2*ms2*s/u2 + 6*s*u/u1 );
  values[21] += + C012s0s/kaellen*b_ux*us* (  - 6*s + 4*s*s/u1 - 2*u + 2*u*u/u1 );
  values[21] += + C012s0s/kaellen*b_tm*ts* (  - 4*m12*ms2/s + 4*m12/s*u + 8*m12 - 4*m12*m12/s - 4*m22*ms2/s + 4*m22/s*u + 8*m22 - 4*m22*m22/s + 8*ms2/s*u + 4*ms2 - 4*s - 4*u );
  values[21] += + C012s0s/kaellen*b_tx*ts* (  - 4*m12*m22*ms2/t1/t2 + 6*m12*m22*s/t1/t2 - 2*m12*m22/t1 + 6*m12*m22/t2 - 8*m12*m22*m22/t1/t2 - 2*m12*ms2*s/t1/t2 + 2*m12*ms2/t1 + 2*m12*ms2/t2 - 6*m12*s/t1 + 4*m12*s/t2 + 2*m12 + 10*m12*m12*m22/t1/t2 + 2*m12*m12*ms2/t1/t2 + 6*m12*m12*s/t1/t2 + 4*m12*m12/t1 - 2*m12*m12/t2 - 4*m12*m12*m12/t1/t2 - 2*m22*ms2*s/t1/t2 + 2*m22*ms2/t1 + 2*m22*ms2/t2 - 8*m22*s/t1 - 4*m22*s/t2 - 4*m22*u/t2 - 4*m22 + 2*m22*m22*ms2/t1/t2 - 8*m22*m22*s/t1/t2 + 2*m22*m22/t1 + 2*m22*m22/t2 + 2*m22*m22*m22/t1/t2 );
  values[21] += + C012s0s/kaellen*b_tx*ts * (  - 2*ms2*s/t1 - 2*ms2*s/t2 - 2*s*u/t2 - 4*s + 2*u + 2*u*u/t2 );
  values[21] +=  2.*C012s0s*b_um*us/sqr(sHat())*(2.*(ms2-mg2)+tHat()-uHat());
  values[21] += + C012s0s*b_ux*us*( 2*m12*mg2/s/u1/u2 - 2*m12*ms2/s/u1/u2 + 2*m22*mg2/s/u1/u2 - 2*m22*ms2/s/u1/u2 + 4*m22/u1/u2 - 2*mg2/u1/u2 + 2*ms2/s/u1 + 2*ms2/s/u2 + 4*ms2/u1/u2 - 2*s/u1/u2 - 2/u1 - 2/u2 );
  values[21] +=  2.*C012s0s*b_tm*ts/sqr(sHat())*(2.*(ms2-mg2)+uHat()-tHat());
  values[21] +=  2.*C012s0s*b_tx *ts/sHat()/t1/t2*((mg2-ms2)*(uHat()+tHat())+2*m12*sHat()+(sHat()-ms2)*(uHat()-tHat()));
  // cerr << "testing V " << values[21] << "\n";
  double D0ts = SusyLoopIntegral::D_div(tHat(),m1,m2,ms,sHat(),scale2);
  values[22] += -4.*D0ts*b_tm;
  values[22] += + D0ts*b_tx*(  - 2 + 4*m12*ms2/t1/t2 + 2*m12/t1 - 2*m12/t2 - 2*m12*m12/t1/t2 - 2*ms2/t1 + 2*ms2/t2 - 2*ms2*ms2/t1/t2);
  // cerr << "testing W " << values[22] << "\n";
  double D0us = SusyLoopIntegral::D_div(uHat(),m2,m1,ms,sHat(),scale2);
  values[23] += -4.*D0us*b_um;
  values[23] +=     D0us*b_ux * (  - 3 + m12/u2 + 4*m22*ms2/u1/u2 - 2*m22/u1 + 2*m22/u2 - 2*m22*m22/u1/u2 + 2*ms2/u1 - 2*ms2/u2 - 2*ms2*ms2/u1/u2 - s/u2 - t/u2 );
  // cerr << "testing X " << values[23] << "\n";
  double D0st = sHat()*tg*SusyLoopIntegral::D_fin(eps,eps,m12,m22,sHat(),tHat(),ms,mg,ms,eps2);
  values[24] += + D0st*b_um* (  - 2*m12*mg2/s/s/tg*us + 2*m12*ms2/s/s/tg*us - 2*m22*mg2/s/s/tg*us + 2*m22*ms2/s/s/tg*us + 4*mg2*ms2/s/s/tg*us + 4*mg2/s/s*us - 4*ms2/s/s*us - 4*ms2*ms2/s/s/tg*us + 2/s*us );
  values[24] += + D0st*b_ux * (  - 4*m12*mg2*ms2/s/u1/tg/u2*us + 4*m12*mg2/u1/tg/u2*us + 2*m12*mg2*mg2/s/u1/tg/u2*us - 4*m12*ms2/u1/tg/u2*us + 2*m12*ms2*ms2/s/u1/tg/u2*us + 2*m12*s/u1/tg/u2*us - 4*m22*mg2*ms2/s/u1/tg/u2*us + 2*m22*mg2/u1/tg/u2*us + 2*m22*mg2*mg2/s/u1/tg/u2*us + 2*m22*ms2*ms2/s/u1/tg/u2*us + 2*mg2*ms2/s/u1/tg*us + 2*mg2*ms2/s/tg/u2*us + 8*mg2*ms2/u1/tg/u2*us - 4*mg2*s/u1/tg/u2*us - 2*mg2/u1/tg*us - 2*mg2*mg2/u1/tg/u2*us + 8*ms2*s/u1/tg/u2*us + 6*ms2/u1/tg*us );
  values[24] += + D0st*b_ux * ( 2*ms2/tg/u2*us - 2*ms2*ms2/s/u1/tg*us - 2*ms2*ms2/s/tg/u2*us - 8*ms2*ms2/u1/tg/u2*us - 2*s/u1/tg*us - 2*s*s/u1/tg/u2*us );
  // cerr << "testing Y " << values[24] << "\n";
  double D0su = sHat()*ug*SusyLoopIntegral::D_fin(eps,eps,m22,m12,sHat(),uHat(),ms,mg,ms,eps2);
  values[25] += + D0su*b_tm *ts* (  - 2*m12*mg2/s/s/ug + 2*m12*ms2/s/s/ug - 2*m22*mg2/s/s/ug + 2*m22*ms2/s/s/ug + 4*mg2*ms2/s/s/ug + 4*mg2/s/s - 4*ms2/s/s - 4*ms2*ms2/s/s/ug + 2/s );
  values[25] += + D0su*b_tx *ts* (  - 4*m12*mg2*ms2/s/t1/ug/t2 + 2*m12*mg2/t1/ug/t2 + 2*m12*mg2*mg2/s/t1/ug/t2 + 2*m12*ms2*ms2/s/t1/ug/t2 - 4*m22*mg2*ms2/s/t1/ug/t2 + 4*m22*mg2/t1/ug/t2 + 2*m22*mg2*mg2/s/t1/ug/t2 - 4*m22*ms2/t1/ug/t2 + 2*m22*ms2*ms2/s/t1/ug/t2 + 2*m22*s/t1/ug/t2 + 2*mg2*ms2/s/t1/ug + 2*mg2*ms2/s/ug/t2 + 8*mg2*ms2/t1/ug/t2 - 4*mg2*s/t1/ug/t2 - 2*mg2/ug/t2 - 2*mg2*mg2/t1/ug/t2 + 8*ms2*s/t1/ug/t2 + 2*ms2/t1/ug );
  values[25] += + D0su*b_tx *ts/ug/t1/t2/s*( 6*ms2*t1*s- 2*ms2*ms2*(m12+m22+2*(s-u)) + 2*s*s*u2 );
  // cerr << "testing Z " << values[25] << "\n";
  // Drell-Yan only
  // output.finite =-8.+sqr(Constants::pi);
  Complex sum=0;
  for(unsigned int ix=0;ix<26;++ix) sum += values[ix];
  // cerr << "testing total " << QBV << "\n";
  return 0.5*real(sum)/(b_s+b_t+b_u)*loWeight();
}
