// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEGQBC1S0Q class.
//

#include "MEGQBC1S0Q.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"
#include <numeric>

using namespace Herwig;

double MEGQBC1S0Q::me2() const {
  // gluon wavefunction
  VectorWaveFunction      g1w(rescaledMomenta()[0],mePartonData()[0],incoming);
  vector<VectorWaveFunction> g1;
  for(unsigned int ix=0;ix<2;++ix) {
    g1w.reset(2*ix);
    g1.push_back(g1w);
  }
  // invariants
  Energy M  = rescaledMomenta()[2].mass();
  Energy2 M2 = sqr(M);
  double a1 = rescaledMomenta()[1].mass()/M, a2 = rescaledMomenta()[3].mass()/M;
  Energy2 sh((rescaledMomenta()[0]+rescaledMomenta()[1]).m2());
  Energy2 th((rescaledMomenta()[0]-rescaledMomenta()[2]).m2());
  Energy2 uh((rescaledMomenta()[0]-rescaledMomenta()[3]).m2());
  // matrix element
  ProductionMatrixElement me(PDT::Spin1,PDT::Spin1Half,PDT::Spin0,PDT::Spin1Half);
  vector<Complex> amp(5);
  DVector save(2,0.);
  double meSum(0.);
  // g c -> B_c b
  if(mePartonData()[1]->id()>0) {
    SpinorWaveFunction      q2w(rescaledMomenta()[1],mePartonData()[1],incoming);
    SpinorBarWaveFunction   q4w(rescaledMomenta()[3],mePartonData()[3],outgoing);
    vector<SpinorWaveFunction> q2;
    vector<SpinorBarWaveFunction> q4;
    for(unsigned int ix=0;ix<2;++ix) {
      q2w.reset(ix);
      q2.push_back(q2w);
      q4w.reset(ix);
      q4.push_back(q4w);
    }
    for(unsigned int ih2=0;ih2<2;++ih2) {
      for(unsigned int ih4=0;ih4<2;++ih4) {
	complex<Energy> PS = q2[ih2].dimensionedWave().pseudoScalar(q4[ih4].dimensionedWave());
	LorentzPolarizationVectorE fCurrent = q2[ih2].dimensionedWave().generalCurrent(q4[ih4].dimensionedWave(),-1.,1.);
	complex<Energy2> PSp1 = fCurrent*rescaledMomenta()[0];
	auto fTemp = q2[ih2].dimensionedWave().slash(rescaledMomenta()[0]).generalCurrent(q4[ih4].dimensionedWave(),-1.,1.);
	for(unsigned int ih1=0;ih1<2;++ih1) {
	  complex<Energy> PSe1 = fCurrent*g1[ih1].wave();
	  complex<Energy2> PSp1e1 = fTemp*g1[ih1].wave();
	  complex<Energy> de1p2 = rescaledMomenta()[1]*g1[ih1].wave();
	  complex<Energy> de1k2 = rescaledMomenta()[3]*g1[ih1].wave();
	  amp[0] = M*(2*a1*sqr(a1)*M2*PSe1 - 2*a1*(2*de1p2*M*PS + M*PSp1e1 + PSe1*sh))/sqr(sh -sqr(a1)*M2);
	  amp[1] = M*(-(M2*PSe1) - sqr(a1)*M2* PSe1 + 2*de1k2*PSp1 + 2*a1*(de1p2*M*PS + M2*PSe1 + de1p2*PSp1 - de1k2*(M*PS + PSp1)) + M*PSp1e1 + PSe1*uh)/(4.*(sqr(a1)*M2 - sh)*(M2 - th));
	  amp[2] = M*(2*(-1 + a1)*(M* (2*de1k2*PS + sqr(a2)*M* PSe1 + PSp1e1) - PSe1*uh))/sqr(uh-sqr(a2)*M2);
	  amp[3] = M*(-2*(-1 + a1)*de1k2* (M*PS + PSp1) + 2*de1p2*((-1 + a1)*M*PS + a1*PSp1) + M*(sqr(a1)*M*PSe1 + PSp1e1) - PSe1*sh)/(4.*(M2 - th)*(sqr(a2)* M2 - uh));
	  amp[4] = M*(9*(-2*de1p2*M*PS + 2*de1k2*PSp1 + 2*a1*(de1p2*M*PS + M2*PSe1 + de1p2*PSp1 - de1k2*(M*PS + PSp1)) - PSe1*(M2 + sh - uh)))/(4.*(sqr(a1)*M2 - sh)*(sqr(a2)* M2 - uh));
	  // diagram weights
	  save[0]+=norm(amp[3]);
	  save[1]+=norm(amp[1]);
	  Complex aSum = std::accumulate(amp.begin(),amp.end(),Complex(0.));
	  meSum+=norm(aSum);
	  me(2*ih1,ih2,0,ih4)=aSum;
	}
      }
    }
  }
  // g cbar -> B_c- bbar
  else {
    SpinorBarWaveFunction q2w(rescaledMomenta()[1],mePartonData()[1],incoming);
    SpinorWaveFunction    q4w(rescaledMomenta()[3],mePartonData()[3],outgoing);
    vector<SpinorBarWaveFunction> q2;
    vector<SpinorWaveFunction> q4;
    for(unsigned int ix=0;ix<2;++ix) {
      q2w.reset(ix);
      q2.push_back(q2w);
      q4w.reset(ix);
      q4.push_back(q4w);
    }
    for(unsigned int ih2=0;ih2<2;++ih2) {
      for(unsigned int ih4=0;ih4<2;++ih4) {
	complex<Energy> PS = q4[ih4].dimensionedWave().pseudoScalar(q2[ih2].dimensionedWave());
	LorentzPolarizationVectorE fCurrent = q4[ih4].dimensionedWave().generalCurrent(q2[ih2].dimensionedWave(),-1.,1.);
	complex<Energy2> PSp1 = fCurrent*rescaledMomenta()[0];
	auto fTemp = q4[ih4].dimensionedWave().slash(rescaledMomenta()[0]).generalCurrent(q2[ih2].dimensionedWave(),-1.,1.);
	for(unsigned int ih1=0;ih1<2;++ih1) {
	  complex<Energy> PSe1 = fCurrent*g1[ih1].wave();
	  complex<Energy2> PSp1e1 = fTemp*g1[ih1].wave();
	  complex<Energy> de1p2 = rescaledMomenta()[1]*g1[ih1].wave();
	  complex<Energy> de1k2 = rescaledMomenta()[3]*g1[ih1].wave();
	  amp[0] = M*(-2*a1*M*(-2*de1p2*PS + sqr(a1)*M*PSe1 + PSp1e1) + 2*a1*PSe1*sh)/sqr(sh-sqr(a1)*M2);
	  amp[1] = M*(sqr(a1)*M2*PSe1 - 2*de1k2*PSp1 + 2*a1*(-(M*(de1p2*PS + M*PSe1)) - de1p2*PSp1 + de1k2*(M*PS + PSp1)) + M*(M*PSe1 + PSp1e1) - PSe1*uh)/(4.*(sqr(a1)*M2 - sh)*(M2 - th));
	  amp[2] = M*(2*(-1 + a1)*(M*(-2*de1k2*PS -   sqr(a2)*M*PSe1 + PSp1e1) + PSe1*uh))/sqr(uh-sqr(a2)*M2);
	  amp[3] = M*(-(sqr(a1)*M2*PSe1) + 2*(-1 + a1)*de1k2*(M*PS + PSp1) - 2*de1p2*((-1 + a1)*M*PS + a1*PSp1) + M*PSp1e1 + PSe1*sh)/(4.*(M2 - th)*(sqr(a2)*M2 - uh));
	  amp[4] = M*(9*(2*de1p2*M*PS - 2*de1k2*PSp1 + 2*a1*(-(M*(de1p2*PS + M*PSe1)) -   de1p2*PSp1 + de1k2*(M*PS + PSp1)) + PSe1*(M2 + sh - uh)))/(4.*(sqr(a1)*M2 - sh)*(sqr(a2)*M2 - uh));
	  // diagram weights
	  save[0]+=norm(amp[3]);
	  save[1]+=norm(amp[1]);
	  Complex aSum = std::accumulate(amp.begin(),amp.end(),Complex(0.));
	  meSum+=norm(aSum);
	  me(2*ih1,ih2,0,ih4)=aSum;

	}
      }
    }
  }
  // save the diagram weights
  meInfo(save);
  // test vs spin summed version
  // double test =M2*(-128*pow(a1,5)*pow<9,1>(M2)+1025*pow(a1,6)*pow<9,1>(M2)-2560*pow(a1,7)*pow<9,1>(M2)+5779*pow(a1,8)*pow<9,1>(M2)-6544*pow(a1,9)*pow<9,1>(M2)+5411*pow(a1,10)*pow<9,1>(M2)
  // 		   -2208*pow(a1,11)*pow<9,1>(M2)-1391*pow(a1,12)*pow<9,1>(M2)-144*pow(a1,13)*pow<9,1>(M2)-64*sqr(a1)*pow<8,1>(M2)*sh+768*pow(a1,3)*pow<8,1>(M2)*sh-2307*pow(a1,4)*pow<8,1>(M2)*sh
  // 		   +7240*pow(a1,5)*pow<8,1>(M2)*sh-15236*pow(a1,6)*pow<8,1>(M2)*sh+18496*pow(a1,7)*pow<8,1>(M2)*sh-16863*pow(a1,8)*pow<8,1>(M2)*sh+4168*pow(a1,9)*pow<8,1>(M2)*sh
  // 		   +5346*pow(a1,10)*pow<8,1>(M2)*sh+720*pow(a1,11)*pow<8,1>(M2)*sh+64*pow<7,1>(M2)*sqr(sh)-128*a1*pow<7,1>(M2)*sqr(sh)+899*sqr(a1)*pow<7,1>(M2)*sqr(sh)
  // 		   -5136*pow(a1,3)*pow<7,1>(M2)*sqr(sh)+11218*pow(a1,4)*pow<7,1>(M2)*sqr(sh)-17376*pow(a1,5)*pow<7,1>(M2)*sqr(sh)+16142*pow(a1,6)*pow<7,1>(M2)*sqr(sh)+544*pow(a1,7)*pow<7,1>(M2)*sqr(sh)
  // 		   -7473*pow(a1,8)*pow<7,1>(M2)*sqr(sh)-1440*pow(a1,9)*pow<7,1>(M2)*sqr(sh)-129*pow<6,1>(M2)*pow<3,1>(sh)+456*a1*pow<6,1>(M2)*pow<3,1>(sh)-1892*sqr(a1)*pow<6,1>(M2)*pow<3,1>(sh)+5440*pow(a1,3)*pow<6,1>(M2)*pow<3,1>(sh)
  // 		   -3406*pow(a1,4)*pow<6,1>(M2)*pow<3,1>(sh)-4560*pow(a1,5)*pow<6,1>(M2)*pow<3,1>(sh)+4252*pow(a1,6)*pow<6,1>(M2)*pow<3,1>(sh)+1440*pow(a1,7)*pow<6,1>(M2)*pow<3,1>(sh)
  // 		   +131*pow<5,1>(M2)*pow<4,1>(sh)-16*a1*pow<5,1>(M2)*pow<4,1>(sh)-1217*sqr(a1)*pow<5,1>(M2)*pow<4,1>(sh)+1856*pow(a1,3)*pow<5,1>(M2)*pow<4,1>(sh)
  // 		   -513*pow(a1,4)*pow<5,1>(M2)*pow<4,1>(sh)-720*pow(a1,5)*pow<5,1>(M2)*pow<4,1>(sh)-67*pow<4,1>(M2)*pow<5,1>(sh)+200*a1*pow<4,1>(M2)*pow<5,1>(sh)
  // 		   -222*sqr(a1)*pow<4,1>(M2)*pow<5,1>(sh)+144*pow(a1,3)*pow<4,1>(M2)*pow<5,1>(sh)+pow<3,1>(M2)*pow<6,1>(sh)
  // 		   +sqr(M2)*(pow(a1,5)*(768+a1*(-5429+a1*(11168+a1*(-23052+a1*(23224+a1*(-17121+2*a1*(4472+a1*(1751-12*a1*(7+3*a1)))))))))*pow<6,1>(M2)
  // 			     +sqr(a1)*(448+a1*(-4736+a1*(12399+a1*(-32448+a1*(61696+3*a1*(-21408+a1*(18719+4*a1*(-1712+a1*(-1171+4*a1*(16+9*a1))))))))))*pow<5,1>(M2)*sh
  // 			     +(-448+a1*(896+a1*(-4799+2*a1*(11792+a1*(-22380+a1*(29736+a1*(-28925+3*a1*(1392+a1*(3523-20*a1*(11+9*a1))))))))))*pow<4,1>(M2)*sqr(sh)
  // 			     +(901+2*a1*(-1152+a1*(3440+a1*(-9584+a1*(8005+4*a1*(1096+a1*(-1759+60*a1*(2+3*a1))))))))*pow<3,1>(M2)*pow<3,1>(sh)+(-764+a1*(696+a1*(2507-6*a1*(856+a1*(-563+20*a1*(1+9*a1))))))*sqr(M2)*pow<4,1>(sh)
  // 			     +3*(99+4*a1*(-32+a1*(13+4*a1*(-4+9*a1))))*M2*pow<5,1>(sh)-2*(25+36*(-1+a1)*a1)*pow<6,1>(sh))*th
  // 		   +M2*(pow(a1,5)*(-1920+a1*(11834+a1*(-19072+a1*(34834+a1*(-29760+a1*(16457+a1*(-10480+3*a1*(-981+64*a1*(5+a1)))))))))*pow<6,1>(M2)
  // 			+sqr(a1)*(-1344+a1*(12288+a1*(-27518+a1*(57392+a1*(-95944+a1*(78336+a1*(-57949+2*a1*(12396+a1*(6241-576*a1*(4+a1))))))))))*pow<5,1>(M2)*sh
  // 			+(1344+a1*(-2688+a1*(10574+a1*(-43488+a1*(68796+a1*(-69504+a1*(63466+a1*(-12032+45*a1*(-453+64*a1*(3+a1))))))))))*pow<4,1>(M2)*sqr(sh)
  // 			-2*(1285+a1*(-2584+a1*(4788+a1*(-11520+a1*(9781+2*a1*(1964+5*a1*(-779+192*a1*(2+a1))))))))*pow<3,1>(M2)*pow<3,1>(sh)
  // 			+(1890+a1*(-2112+a1*(-1763+5*a1*(1008+a1*(-1021+576*a1*(1+a1))))))*sqr(M2)*pow<4,1>(sh)+(-649+2*a1*(268+129*a1-576*pow(a1,3)))*M2*pow<5,1>(sh)
  // 			+(113+192*(-1+a1)*a1)*pow<6,1>(sh))*sqr(th)+(pow(a1,5)*(2560+a1*(-13610+a1*(15808+a1*(-24108+a1*(16024+a1*(-3387-32*a1*(-97+2*a1*(-9+2*a1*(5+a1)))))))))*pow<6,1>(M2)
  // 								     +sqr(a1)*(2240+a1*(-17280+a1*(32350+a1*(-49888+a1*(70528+a1*(-37216+a1*(15903+32*a1*(-167+2*a1*(-37+12*a1*(4+a1))))))))))*pow<5,1>(M2)*sh
  // 								     -2*(1120+a1*(-2240+a1*(6175+a1*(-20416+a1*(25324+a1*(-14664+a1*(9935+96*a1*(9+a1*(-19+10*a1*(3+a1))))))))))*pow<4,1>(M2)*sqr(sh)
  // 								     +2*(1925+a1*(-3376+a1*(3312+a1*(-5552+a1*(3127+32*a1*(97+2*a1*(-19+20*a1*(2+a1))))))))*pow<3,1>(M2)*pow<3,1>(sh)
  // 								     +(-2396+a1*(2968+a1*(425-32*a1*(43+2*a1*(-7+30*a1*(1+a1))))))*sqr(M2)*pow<4,1>(sh)+3*(225+32*a1*(-9+2*a1+8*pow(a1,3)))*M2*pow<5,1>(sh)
  // 								     -64*(1+2*(-1+a1)*a1)*pow<6,1>(sh))*pow<3,1>(th)+(pow(a1,5)*(-1920+a1*(8805+a1*(-6272+a1*(7155+16*a1*(-184+a1*(-73+8*a1*(5+2*a1)))))))*pow<5,1>(M2)
  // 														      +sqr(a1)*(-2240+a1*(14080+a1*(-21455-4*a1*(-5298+a1*(6085+16*a1*(-60+a1*(-27+48*a1+22*sqr(a1))))))))*pow<4,1>(M2)*sh
  // 														      +(2240+a1*(-4480+a1*(8175+2*a1*(-10184+a1*(9249+16*a1*(-36+a1+152*sqr(a1)+96*pow(a1,3)))))))*pow<3,1>(M2)*sqr(sh)
  // 														      -(3205+4*a1*(-1362+a1*(709+16*a1*(-40+a1*(13+40*a1+52*sqr(a1))))))*sqr(M2)*pow<3,1>(sh)
  // 														      +(1523+16*a1*(-144+a1*(31+8*a1*(-3+14*a1))))*M2*pow<4,1>(sh)-128*(2+a1*(-4+3*a1))*pow<5,1>(sh))*pow<4,1>(th)
  // 		   +(pow(a1,5)*(768-a1*(3185+32*a1*(-29+23*a1+6*pow(a1,3))))*pow<4,1>(M2)+sqr(a1)*(1344+a1*(-6528+a1*(7907+32*a1*(-109+a1*(119+8*a1*(3+4*a1))))))*pow<3,1>(M2)*sh
  // 		     -(1344+a1*(-2688+a1*(3059+96*a1*(-53+a1*(39+4*a1*(2+5*a1))))))*sqr(M2)*sqr(sh)+(1409+32*a1*(-79+3*a1*(11+8*a1*(-1+2*a1))))*M2*pow<3,1>(sh)-64*(6+a1*(-12+7*a1))*pow<4,1>(sh))*pow<5,1>(th)
  // 		   +16*(pow(a1,5)*(-8+39*a1+8*pow(a1,3))*pow<3,1>(M2)-2*sqr(a1)*(14+a1*(-48+47*a1+16*pow(a1,3)))*sqr(M2)*sh+(28+a1*(-56+a1*(39+8*a1*(-4+5*a1))))*M2*sqr(sh)
  // 			-16*sqr(-1+a1)*pow<3,1>(sh))*pow<6,1>(th)-64*(sqr(a1)*M2-sh)*(pow(a1,4)*M2-sqr(-1+a1)*sh)*pow<7,1>(th))/(4.*pow<4,1>(-(sqr(a1)*M2)+sh)*sqr(M2-th)*pow<4,1>(-(sqr(-1+a1)*M2)+uh));
  // cerr << "testing me " << meSum << " "  << test << " " << (meSum-test)/(meSum+test) << " " << meSum/test << "\n";
  // set matrix element
  setME(me);
  // final factors
  return 16.*O1_*pow<3,1>(Constants::pi*standardModel()->alphaS(scale())/M)/(729.*sqr(a1*a2))*meSum;
}

IBPtr MEGQBC1S0Q::clone() const {
  return new_ptr(*this);
}

IBPtr MEGQBC1S0Q::fullclone() const {
  return new_ptr(*this);
}

void MEGQBC1S0Q::doinit() {
  MEGQBCQBase::doinit();
  O1_ = oniumParameters()->singletMEProduction<0>(bcbar,principleQuantumNumber(),0,0);
}

void MEGQBC1S0Q::persistentOutput(PersistentOStream & os) const {
  os << ounit(O1_,GeV*GeV2);
}

void MEGQBC1S0Q::persistentInput(PersistentIStream & is, int) {
  is >> iunit(O1_,GeV*GeV2);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEGQBC1S0Q,MEGQBCQBase>
describeHerwigMEGQBC1S0Q("Herwig::MEGQBC1S0Q", "HwOniumParameters.so HwMEHadronOnium.so");

void MEGQBC1S0Q::Init() {

  static ClassDocumentation<MEGQBC1S0Q> documentation
    ("The MEGQBC1S0Q class implements the matrix element for g c -> B_c(1S0) b and charged conjugate");

}

