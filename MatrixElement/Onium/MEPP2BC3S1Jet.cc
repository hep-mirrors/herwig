// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2BC3S1Jet class.
//

#include "MEPP2BC3S1Jet.h"
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

double MEPP2BC3S1Jet::me2() const {
  // gluon wavefunction
  VectorWaveFunction g1w(rescaledMomenta()[0],mePartonData()[0],incoming);
  vector<VectorWaveFunction> g1;
  for(unsigned int ix=0;ix<2;++ix) {
    g1w.reset(2*ix);
    g1.push_back(g1w);
  }
  // B_C* wavefunction
  VectorWaveFunction Bcw(rescaledMomenta()[2],mePartonData()[2],outgoing);
  vector<VectorWaveFunction> v3;
  for(unsigned int ix=0;ix<3;++ix) {
    Bcw.reset(ix);
    v3.push_back(Bcw);
  }
  // invariants
  Energy M  = rescaledMomenta()[2].mass();
  Energy2 M2 = sqr(M);
  double a1 = rescaledMomenta()[1].mass()/M, a2 = rescaledMomenta()[3].mass()/M;
  Energy2 sh((rescaledMomenta()[0]+rescaledMomenta()[1]).m2());
  Energy2 th((rescaledMomenta()[0]-rescaledMomenta()[2]).m2());
  Energy2 uh((rescaledMomenta()[0]-rescaledMomenta()[3]).m2());
  // matrix element
  ProductionMatrixElement me(PDT::Spin1,PDT::Spin1Half,PDT::Spin1,PDT::Spin1Half);
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
	complex<Energy> SS = q2[ih2].dimensionedWave().scalar(q4[ih4].dimensionedWave());
	LorentzPolarizationVectorE fCurrent = q2[ih2].dimensionedWave().vectorCurrent(q4[ih4].dimensionedWave());
	complex<Energy2> SSp1 = fCurrent*rescaledMomenta()[0];
	auto fTemp = q2[ih2].dimensionedWave().slash(rescaledMomenta()[0]).vectorCurrent(q4[ih4].dimensionedWave());
	for(unsigned int ih3=0;ih3<3;++ih3) {
	  complex<Energy> SSeBc = fCurrent*v3[ih3].wave();
	  complex<Energy2> SSeBcp1 = fTemp*v3[ih3].wave();
	  complex<Energy> deBcp1 = rescaledMomenta()[0]*v3[ih3].wave();
	  complex<Energy> deBcp2 = rescaledMomenta()[1]*v3[ih3].wave();
	  auto fTemp2 = q2[ih2].dimensionedWave().slash(v3[ih3].wave()).vectorCurrent(q4[ih4].dimensionedWave());
	  auto fTemp3 = q2[ih2].dimensionedWave().slash(rescaledMomenta()[0]).slash(v3[ih3].wave()).vectorCurrent(q4[ih4].dimensionedWave());
	  for(unsigned int ih1=0;ih1<2;++ih1) {
	    complex<Energy> SSe1 = fCurrent*g1[ih1].wave();
	    complex<Energy2> SSe1p1 = fTemp*g1[ih1].wave();
	    Complex de1eBc = g1[ih1].wave()*v3[ih3].wave();
	    complex<Energy> de1p2 = rescaledMomenta()[1]*g1[ih1].wave();
	    complex<Energy> de1k2 = rescaledMomenta()[3]*g1[ih1].wave();
	    complex<Energy> SSe1eBc = g1[ih1].wave()*fTemp2;
	    complex<Energy2> SSe1eBcp1 = g1[ih1].wave()*fTemp3;
	    amp[0] = M*(-2.*a1*M*(SSe1eBcp1 + 2.*de1p2*SSeBc - 2.*de1eBc*SSp1))/sqr(sh-sqr(a1)*M2);
	    amp[1] = M*(M*SSe1eBcp1 + 2.*(deBcp1 + 2.*deBcp2)*SSe1p1 + 2.*de1p2*(2.*deBcp1*SS + a1*M*SSeBc - SSeBcp1) - 
			2.*de1k2*(a1*M*SSeBc + SSeBcp1) + SSe1eBc*((-1 + 2.*a1)*M2 - sh + uh) + 
			2.*de1eBc*(sqr(a2)*M2*SS - M*SSp1 - SS*uh))/(4.*(sqr(a1)*M2 - sh)*(M2 - th));
	    amp[2] = M*(2.*(-1 + a1)*M*(-2.*deBcp1*SSe1 + SSe1eBcp1 + 2.*de1k2*SSeBc))/sqr(uh-sqr(a2)*M2);
	    amp[3] = M*(4*de1p2*deBcp1*SS + M*SSe1eBcp1 + 4*deBcp2*SSe1p1 + 2.*deBcp1*(-(M*SSe1) + SSe1p1) - 
			2.*(-1 + a1)*(de1k2 - de1p2)*M*SSeBc - 2.*de1k2*SSeBcp1 - 2.*de1p2*SSeBcp1 +
			2.*de1eBc*SS*(sqr(a2)*M2 - uh) + SSe1eBc*((-1 + 2.*a1)*M2 - sh + uh))/(4.*(M2 - th)*(sqr(a2)*M2 - uh));
	    amp[4] = M*(9*(2.*(2.*de1p2*deBcp1*SS + (deBcp1 + 2.*deBcp2)*SSe1p1 - a1*de1k2*M*SSeBc + (-1 + a1)*de1p2*M*SSeBc - 
			       (de1k2 + de1p2)*SSeBcp1 + de1eBc*SS*(sqr(a2)*M2 - uh)) + 
			   SSe1eBc*((-1 + 2.*a1)*M2 - sh + uh)))/(4.*(sqr(a1)*M2 - sh)*(sqr(a2)*M2 - uh));
	    // diagram weights
	    save[0]+=norm(amp[3]);
	    save[1]+=norm(amp[1]);
	    Complex aSum = std::accumulate(amp.begin(),amp.end(),Complex(0.));
	    meSum+=norm(aSum);
	    me(2*ih1,ih2,ih3,ih4)=aSum;
	  }
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
	complex<Energy> SS = q4[ih4].dimensionedWave().scalar(q2[ih2].dimensionedWave());
	LorentzPolarizationVectorE fCurrent = q4[ih4].dimensionedWave().vectorCurrent(q2[ih2].dimensionedWave());
	complex<Energy2> SSp1 = fCurrent*rescaledMomenta()[0];
	auto fTemp = q4[ih4].dimensionedWave().slash(rescaledMomenta()[0]).vectorCurrent(q2[ih2].dimensionedWave());
	for(unsigned int ih3=0;ih3<3;++ih3) {
	  complex<Energy> SSeBc = fCurrent*v3[ih3].wave();
	  complex<Energy2> SSeBcp1 = fTemp*v3[ih3].wave();
	  complex<Energy> deBcp1 = rescaledMomenta()[0]*v3[ih3].wave();
	  complex<Energy> deBcp2 = rescaledMomenta()[1]*v3[ih3].wave();
	  auto fTemp2 = q4[ih4].dimensionedWave().slash(v3[ih3].wave()).vectorCurrent(q2[ih2].dimensionedWave());
	  auto fTemp3 = q4[ih4].dimensionedWave().slash(rescaledMomenta()[0]).slash(v3[ih3].wave()).vectorCurrent(q2[ih2].dimensionedWave());
	  for(unsigned int ih1=0;ih1<2;++ih1) {
	    complex<Energy> SSe1 = fCurrent*g1[ih1].wave();
	    complex<Energy2> SSe1p1 = fTemp*g1[ih1].wave();
	    Complex de1eBc = g1[ih1].wave()*v3[ih3].wave();
	    complex<Energy> de1p2 = rescaledMomenta()[1]*g1[ih1].wave();
	    complex<Energy> de1k2 = rescaledMomenta()[3]*g1[ih1].wave();
	    complex<Energy> SSe1eBc = g1[ih1].wave()*fTemp2;
	    complex<Energy2> SSe1eBcp1 = g1[ih1].wave()*fTemp3;
	    amp[0] = M*(2.*a1*M*(2.*deBcp1*SSe1 - SSe1eBcp1 + 2.*de1p2*SSeBc))/sqr(sh-sqr(a1)*M2);
	    amp[1] = M*(2.*sqr(a1)*de1eBc*M2*SS - 2.*de1eBc*sh*SS - 2.*deBcp1*M*SSe1 + M2*SSe1eBc + sh*SSe1eBc + 
			M*SSe1eBcp1 - 2.*deBcp1*SSe1p1 - 4.*deBcp2*SSe1p1 - 2.*a1*M*(M*SSe1eBc + de1p2*SSeBc) + 2.*de1p2*SSeBcp1 + 
			2.*de1k2*(-2.*deBcp1*SS + a1*M*SSeBc + SSeBcp1) - SSe1eBc*uh)/(4.*(sqr(a1)*M2 - sh)*(M2 - th));
	    amp[2] = M*(2.*(-1 + a1)*M*(SSe1eBcp1 - 2.*(de1k2*SSeBc + de1eBc*SSp1)))/sqr(uh-sqr(a2)*M2);
	    amp[3] = M*(2.*sqr(a1)*de1eBc*M2*SS - 2.*de1eBc*sh*SS + sh*SSe1eBc - 2.*(deBcp1 + 2.*deBcp2)*SSe1p1 - 
			2.*a1*M*(M*SSe1eBc + de1p2*SSeBc) + M*(M*SSe1eBc + SSe1eBcp1 + 2.*de1p2*SSeBc) + 
			2.*de1p2*SSeBcp1 + 2.*de1k2*(-2.*deBcp1*SS + (-1 + a1)*M*SSeBc + SSeBcp1) - 
			2.*de1eBc*M*SSp1 - SSe1eBc*uh)/(4.*(M2 - th)*(sqr(a2)*M2 - uh));
	    amp[4] = M*(9*(2.*sqr(a1)*de1eBc*M2*SS - 2.*de1eBc*sh*SS + M2*SSe1eBc + sh*SSe1eBc - 
			   2.*(deBcp1 + 2.*deBcp2)*SSe1p1 + 2.*de1p2*M*SSeBc - 
			   2.*a1*M*(M*SSe1eBc + de1p2*SSeBc) + 2.*de1p2*SSeBcp1 + 
			   2.*de1k2*(-2.*deBcp1*SS + a1*M*SSeBc + SSeBcp1) - SSe1eBc*uh))/
	      (4.*(sqr(a1)*M2 - sh)*(sqr(a2)*M2 - uh));
	    // diagram weights
	    save[0]+=norm(amp[3]);
	    save[1]+=norm(amp[1]);
	    Complex aSum = std::accumulate(amp.begin(),amp.end(),Complex(0.));
	    meSum+=norm(aSum);
	    me(2*ih1,ih2,ih3,ih4)=aSum;
	  }
     	}
      }
    }
  }
  // save the diagram weights
  meInfo(save);
  // test vs spin summed version
  // double test =M2*(-384*pow(a1,5)*pow<9,1>(M2)+3011*pow(a1,6)*pow<9,1>(M2)-11812*pow(a1,7)*pow<9,1>(M2)+25697*pow(a1,8)*pow<9,1>(M2)-37020*pow(a1,9)*pow<9,1>(M2)+33957*pow(a1,10)*pow<9,1>(M2)
  // 		   -20364*pow(a1,11)*pow<9,1>(M2)+5431*pow(a1,12)*pow<9,1>(M2)-980*pow(a1,13)*pow<9,1>(M2)+308*pow(a1,14)*pow<9,1>(M2)-64*pow(a1,15)*pow<9,1>(M2)+4*pow(a1,16)*pow<9,1>(M2)
  // 		   -192*sqr(a1)*pow<8,1>(M2)*sh+2048*pow(a1,3)*pow<8,1>(M2)*sh-8585*pow(a1,4)*pow<8,1>(M2)*sh+32804*pow(a1,5)*pow<8,1>(M2)*sh-74568*pow(a1,6)*pow<8,1>(M2)*sh
  // 		   +115184*pow(a1,7)*pow<8,1>(M2)*sh-115261*pow(a1,8)*pow<8,1>(M2)*sh+72340*pow(a1,9)*pow<8,1>(M2)*sh-23814*pow(a1,10)*pow<8,1>(M2)*sh+5384*pow(a1,11)*pow<8,1>(M2)*sh
  // 		   -1896*pow(a1,12)*pow<8,1>(M2)*sh+448*pow(a1,13)*pow<8,1>(M2)*sh-32*pow(a1,14)*pow<8,1>(M2)*sh+192*pow<7,1>(M2)*sqr(sh)-128*a1*pow<7,1>(M2)*sqr(sh)+4169*sqr(a1)*pow<7,1>(M2)*sqr(sh)
  // 		   -23132*pow(a1,3)*pow<7,1>(M2)*sqr(sh)+66234*pow(a1,4)*pow<7,1>(M2)*sqr(sh)-120744*pow(a1,5)*pow<7,1>(M2)*sqr(sh)+136530*pow(a1,6)*pow<7,1>(M2)*sqr(sh)
  // 		   -96216*pow(a1,7)*pow<7,1>(M2)*sqr(sh)+41229*pow(a1,8)*pow<7,1>(M2)*sqr(sh)-12220*pow(a1,9)*pow<7,1>(M2)*sqr(sh)+4908*pow(a1,10)*pow<7,1>(M2)*sqr(sh)
  // 		   -1344*pow(a1,11)*pow<7,1>(M2)*sqr(sh)+112*pow(a1,12)*pow<7,1>(M2)*sqr(sh)-131*pow<6,1>(M2)*pow<3,1>(sh)+2140*a1*pow<6,1>(M2)*pow<3,1>(sh)-17696*sqr(a1)*pow<6,1>(M2)*pow<3,1>(sh)
  // 		   +44016*pow(a1,3)*pow<6,1>(M2)*pow<3,1>(sh)-63738*pow(a1,4)*pow<6,1>(M2)*pow<3,1>(sh)+58248*pow(a1,5)*pow<6,1>(M2)*pow<3,1>(sh)-35396*pow(a1,6)*pow<6,1>(M2)*pow<3,1>(sh)
  // 		   +14640*pow(a1,7)*pow<6,1>(M2)*pow<3,1>(sh)-6880*pow(a1,8)*pow<6,1>(M2)*pow<3,1>(sh)+2240*pow(a1,9)*pow<6,1>(M2)*pow<3,1>(sh)-224*pow(a1,10)*pow<6,1>(M2)*pow<3,1>(sh)
  // 		   +333*pow<5,1>(M2)*pow<4,1>(sh)-1436*a1*pow<5,1>(M2)*pow<4,1>(sh)+9145*sqr(a1)*pow<5,1>(M2)*pow<4,1>(sh)-15388*pow(a1,3)*pow<5,1>(M2)*pow<4,1>(sh)
  // 		   +15489*pow(a1,4)*pow<5,1>(M2)*pow<4,1>(sh)-9740*pow(a1,5)*pow<5,1>(M2)*pow<4,1>(sh)+5580*pow(a1,6)*pow<5,1>(M2)*pow<4,1>(sh)-2240*pow(a1,7)*pow<5,1>(M2)*pow<4,1>(sh)
  // 		   +280*pow(a1,8)*pow<5,1>(M2)*pow<4,1>(sh)-633*pow<4,1>(M2)*pow<5,1>(sh)+1380*a1*pow<4,1>(M2)*pow<5,1>(sh)-3222*sqr(a1)*pow<4,1>(M2)*pow<5,1>(sh)
  // 		   +3400*pow(a1,3)*pow<4,1>(M2)*pow<5,1>(sh)-2568*pow(a1,4)*pow<4,1>(M2)*pow<5,1>(sh)+1344*pow(a1,5)*pow<4,1>(M2)*pow<5,1>(sh)-224*pow(a1,6)*pow<4,1>(M2)*pow<5,1>(sh)
  // 		   +283*pow<3,1>(M2)*pow<6,1>(sh)-484*a1*pow<3,1>(M2)*pow<6,1>(sh)+596*sqr(a1)*pow<3,1>(M2)*pow<6,1>(sh)-448*pow(a1,3)*pow<3,1>(M2)*pow<6,1>(sh)
  // 		   +112*pow(a1,4)*pow<3,1>(M2)*pow<6,1>(sh)-48*sqr(M2)*pow<7,1>(sh)+64*a1*sqr(M2)*pow<7,1>(sh)-32*sqr(a1)*sqr(M2)*pow<7,1>(sh)+4*M2*pow<8,1>(sh)
  // 		   +M2*(pow(a1,5)*(2304+a1*(-15775+a1*(54708+a1*(-107816+a1*(143704+a1*(-124771+2*a1*(37818+a1*(-11893+8*a1*(265+a1*(-49+4*a1))))))))))*pow<7,1>(M2)
  // 			+sqr(a1)*(1344+a1*(-12416+a1*(46925+a1*(-154876+a1*(318876+a1*(-451808+a1*(435239-4*a1*(69277+2*a1*(-12985+2*a1*(1399+a1*(-297+28*a1)))))))))))*pow<6,1>(M2)*sh
  // 			+(-1344+a1*(896+a1*(-23101+2*a1*(54958+a1*(-142498+a1*(240456+a1*(-267095+a1*(191268+a1*(-89121+8*a1*(3020-753*a1+84*sqr(a1)))))))))))*pow<5,1>(M2)*sqr(sh)
  // 			+(1167-2*a1*(4874+a1*(-37674+a1*(90608+a1*(-131711+4*a1*(30165+2*a1*(-9368+5*a1*(678+a1*(-205+28*a1)))))))))*pow<4,1>(M2)*pow<3,1>(sh)
  // 			+(-1412+a1*(8408+a1*(-41375+2*a1*(32642+a1*(-31111+40*a1*(413+a1*(-159+28*a1)))))))*pow<3,1>(M2)*pow<4,1>(sh)+(1675-4*a1*(1257+2*a1*(-1383+2*a1*(635-339*a1+84*sqr(a1)))))*sqr(M2)*pow<5,1>(sh)
  // 			+2*(-291+8*a1*(74+a1*(-67+28*a1)))*M2*pow<6,1>(sh)+16*(3-4*a1)*pow<7,1>(sh))*th
  // 		   +M2*(pow(a1,5)*(-5760+a1*(33774+a1*(-101256+a1*(175542+a1*(-209224+a1*(163379+a1*(-92852+a1*(30687-4708*a1+580*sqr(a1)))))))))*pow<6,1>(M2)
  // 			+sqr(a1)*(-4032+a1*(31488+a1*(-105978+a1*(293576+a1*(-533872+a1*(665376+a1*(-585607+2*a1*(174494+a1*(-66639+4*(3043-435*a1)*a1)))))))))*pow<5,1>(M2)*sh
  // 			+(4032+a1*(-2688+a1*(53034+a1*(-210424+a1*(482428+a1*(-718512+a1*(741390+a1*(-495848+5*a1*(45321+20*a1*(-511+87*a1))))))))))*pow<4,1>(M2)*sqr(sh)
  // 			-2*(1935+a1*(-9052+a1*(63424+a1*(-138896+a1*(190591+2*a1*(-80534+5*a1*(9369+4*a1*(-689+145*a1))))))))*pow<3,1>(M2)*pow<3,1>(sh)
  // 			+(2750+a1*(-15432+a1*(63727+5*a1*(-17684+a1*(14997+4*a1*(-1579+435*a1))))))*sqr(M2)*pow<4,1>(sh)
  // 			+(-1707+2*a1*(2998+a1*(-5991+4*(1091-435*a1)*a1)))*M2*pow<5,1>(sh)+(363-804*a1+580*sqr(a1))*pow<6,1>(sh))*sqr(th)
  // 		   +(pow(a1,5)*(7680+a1*(-37630+a1*(94184+a1*(-137084+a1*(137592+a1*(-88977+4*a1*(10171+a1*(-3211+16*(23-2*a1)*a1))))))))*pow<6,1>(M2)
  // 		     +sqr(a1)*(6720+a1*(-42880+a1*(126490+a1*(-281624+a1*(435064+a1*(-444896+a1*(329101+12*a1*(-13181+a1*(4637+16*a1*(-39+4*a1))))))))))*pow<5,1>(M2)*sh
  // 		     -2*(3360+a1*(-2240+a1*(32285+a1*(-102588+a1*(200272+a1*(-244584+a1*(214893+4*a1*(-29159+15*a1*(785+16*(-8+a1)*a1)))))))))*pow<4,1>(M2)*sqr(sh)
  // 		     +2*(3215+a1*(-8868+a1*(52836+a1*(-97008+a1*(114549+4*a1*(-19703+5*a1*(1931+16*a1*(-25+4*a1))))))))*pow<3,1>(M2)*pow<3,1>(sh)
  // 		     +(-3108+a1*(12152+a1*(-40357-60*a1*(-745+a1*(505+16*a1*(-9+2*a1))))))*sqr(M2)*pow<4,1>(sh)+(921+4*a1*(-715+3*a1*(377+16*a1*(-11+4*a1))))*M2*pow<5,1>(sh)
  // 		     -64*(1+2*(-1+a1)*a1)*pow<6,1>(sh))*pow<3,1>(th)+(pow(a1,5)*(-5760+a1*(22895+a1*(-44916+a1*(50705+4*a1*(-9307+a1*(4279+8*a1*(-97+16*a1)))))))*pow<5,1>(M2)
  // 								      +sqr(a1)*(-6720+a1*(33280+a1*(-84205-4*a1*(-35029+2*a1*(21525+2*a1*(-7811+a1*(4127+8*a1*(-109+19*a1))))))))*pow<4,1>(M2)*sh
  // 								      +(6720+a1*(-4480+a1*(44045+2*a1*(-52582+a1*(81757+4*a1*(-17841+a1*(11189+8*a1*(-371+72*a1))))))))*pow<3,1>(M2)*sqr(sh)
  // 								      +(-5775+4*a1*(2491-4*a1*(2748+a1*(-3715+a1*(3103+8*a1*(-145+34*a1))))))*sqr(M2)*pow<3,1>(sh)
  // 								      +(1949+4*a1*(-1115+a1*(2327+8*a1*(-193+64*a1))))*M2*pow<4,1>(sh)-128*(2+a1*(-4+3*a1))*pow<5,1>(sh))*pow<4,1>(th)
  // 		   +(pow(a1,5)*(2304+a1*(-7251+4*a1*(2409+a1*(-1857+16*(34-11*a1)*a1))))*pow<4,1>(M2)+sqr(a1)*(4032+a1*(-14208+a1*(29865+4*a1*(-8027+a1*(6931+32*a1*(-69+20*a1))))))*pow<3,1>(M2)*sh
  // 		     +(-4032+a1*(2688+a1*(-16089-4*a1*(-6427+a1*(6883+96*a1*(-31+9*a1))))))*sqr(M2)*sqr(sh)+(2691+4*a1*(-809+a1*(1937+32*a1*(-47+16*a1))))*M2*pow<3,1>(sh)
  // 		     -64*(8+a1*(-12+7*a1))*pow<4,1>(sh))*pow<5,1>(th)+16*(pow(a1,5)*(-24+a1*(65-34*a1+24*sqr(a1)))*pow<3,1>(M2)-2*sqr(a1)*(42+a1*(-88+a1*(145-66*a1+32*sqr(a1))))*sqr(M2)*sh
  // 									  +(84+a1*(-56+a1*(161+2*a1*(-65+28*a1))))*M2*sqr(sh)-16*(2+(-2+a1)*a1)*pow<3,1>(sh))*pow<6,1>(th)
  // 		   -64*(sqr(a1)*M2-sh)*(pow(a1,4)*M2-(3+(-2+a1)*a1)*sh)*pow<7,1>(th))/(4.*pow<4,1>(-(sqr(a1)*M2)+sh)*sqr(M2-th)*pow<4,1>(-(sqr(-1+a1)*M2)+uh));
  // cerr << "testing me " << meSum << " "  << test << " " << (meSum-test)/(meSum+test) << " " << meSum/test << "\n";
  // set matrix element
  setME(me);
  // final factors
  return 16.*O1_*pow<3,1>(Constants::pi*standardModel()->alphaS(scale())/M)/(729.*sqr(a1*a2))*meSum;
}

IBPtr MEPP2BC3S1Jet::clone() const {
  return new_ptr(*this);
}

IBPtr MEPP2BC3S1Jet::fullclone() const {
  return new_ptr(*this);
}

void MEPP2BC3S1Jet::doinit() {
  MEPP2BCJetBase::doinit();
  O1_ = oniumParameters()->singletMEProduction<0>(bcbar,principleQuantumNumber(),1,1);
}

void MEPP2BC3S1Jet::persistentOutput(PersistentOStream & os) const {
  os << ounit(O1_,GeV*GeV2);
}

void MEPP2BC3S1Jet::persistentInput(PersistentIStream & is, int) {
  is >> iunit(O1_,GeV*GeV2);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEPP2BC3S1Jet,MEPP2BCJetBase>
describeHerwigMEPP2BC3S1Jet("Herwig::MEPP2BC3S1Jet", "HwOniumParameters.so HwMEHadronOnium.so");

void MEPP2BC3S1Jet::Init() {

  static ClassDocumentation<MEPP2BC3S1Jet> documentation
    ("The MEPP2BC3S1Jet class implements the matrix element for g c -> B_c(3S1) b and charged conjugate");

}

