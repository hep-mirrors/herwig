// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2NeutralinoNeutralinoPowheg class.
//

#include "MEPP2NeutralinoNeutralinoPowheg.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "Herwig++/Models/Susy/SusyBase.h"
#include "Herwig++/MatrixElement/HardVertex.h"
#include "ThePEG/Helicity/Vertex/Scalar/FFSVertex.h"
#include <numeric>

using namespace Herwig;
using ThePEG::Helicity::VectorWaveFunction;
using ThePEG::Helicity::incoming;
using ThePEG::Helicity::outgoing;


MEPP2NeutralinoNeutralinoPowheg::MEPP2NeutralinoNeutralinoPowheg() 
  : process_(0), maxFlavour_(5) {
  vector<unsigned int> mopt(2,1);
  massOption(mopt);
}

void MEPP2NeutralinoNeutralinoPowheg::doinit() {
  NLODrellYanBase::doinit();
  // get the photon and Z ParticleData objects
  Z0_    = getParticleData(ThePEG::ParticleID::Z0);
  // cast the SM pointer to the Herwig SM pointer
  tcSusyBasePtr hwsm=ThePEG::dynamic_ptr_cast<tcSusyBasePtr>(standardModel());
  if(!hwsm)
    throw InitException() << "Must be the Herwig++ SusyBase class in "
			  << "MEPP2NeutralinoNeutralinoPowheg::doinit" 
			  << Exception::abortnow;
  // do the initialisation (see Herwig::SusyBase Class)
  FFZVertex_ = hwsm->vertexFFZ();
  FFGVertex_ = hwsm->vertexFFG();
  NNZVertex_ = hwsm->vertexNNZ();
  NFSVertex_ = hwsm->vertexNFSF();
  GSSVertex_ = hwsm->vertexGSFSF();
}

Selector<const ColourLines *>
MEPP2NeutralinoNeutralinoPowheg::colourGeometries(tcDiagPtr diag) const {
  static const ColourLines c1("1 -2"), c2("1 2 -3");
  Selector<const ColourLines *> sel;
  if(abs(diag->id())==1)
    sel.insert(1.0, &c1);
  else
    sel.insert(1.0, &c2);
  return sel;
}

void MEPP2NeutralinoNeutralinoPowheg::getDiagrams() const {
  // loop over the processes we need
  tcPDPtr chi[4] = {getParticleData(1000022),getParticleData(1000023),
		    getParticleData(1000025),getParticleData(1000035)};
  for(int i = 1; i <= maxFlavour_; ++i) {
    tcPDPtr q  = getParticleData(i);
    tcPDPtr qb = q->CC();
    tcPDPtr qL = getParticleData(1000000+i);
    tcPDPtr qR = getParticleData(2000000+i);
    for(int c1=0;c1<4;++c1) {
      for(int c2=0;c2<=c1;++c2) {
	if(process_==0 || process_ == 4*c2+c1+1) {
	  add(new_ptr((Tree2toNDiagram(2), q, qb, 1, Z0_   ,
		       3, chi[c1], 3, chi[c2], -1)));
	  add(new_ptr((Tree2toNDiagram(3), q, qL, qb,
		       1, chi[c1], 3, chi[c2], -2)));
	  add(new_ptr((Tree2toNDiagram(3), q, qR, qb,
		       1, chi[c1], 3, chi[c2], -3)));
	  add(new_ptr((Tree2toNDiagram(3), q, qL, qb,
		       3, chi[c1], 1, chi[c2], -4)));
	  add(new_ptr((Tree2toNDiagram(3), q, qR, qb,
		       3, chi[c1], 1, chi[c2], -5)));
	}
      }
    }
  }
}

unsigned int MEPP2NeutralinoNeutralinoPowheg::orderInAlphaS() const {
  return 0;
}

unsigned int MEPP2NeutralinoNeutralinoPowheg::orderInAlphaEW() const {
  return 2;
}

IBPtr MEPP2NeutralinoNeutralinoPowheg::clone() const {
  return new_ptr(*this);
}

IBPtr MEPP2NeutralinoNeutralinoPowheg::fullclone() const {
  return new_ptr(*this);
}

void MEPP2NeutralinoNeutralinoPowheg::persistentOutput(PersistentOStream & os) const {
  os << FFZVertex_ << FFGVertex_ << NNZVertex_ << NFSVertex_ << GSSVertex_
     << Z0_ << process_ << maxFlavour_;
}

void MEPP2NeutralinoNeutralinoPowheg::persistentInput(PersistentIStream & is, int) {
  is >> FFZVertex_ >> FFGVertex_ >> NNZVertex_ >> NFSVertex_ >> GSSVertex_
     >> Z0_ >> process_ >> maxFlavour_;
}

ClassDescription<MEPP2NeutralinoNeutralinoPowheg> 
MEPP2NeutralinoNeutralinoPowheg::initMEPP2NeutralinoNeutralinoPowheg;
// Definition of the static class description member.

void MEPP2NeutralinoNeutralinoPowheg::Init() {

  static ClassDocumentation<MEPP2NeutralinoNeutralinoPowheg> documentation
    ("MEPP2NeutralinoNeutralinoPowheg implements the ME calculation"
     " of the fermion-antifermion to chargino-chargino or"
     " neutralino-neutralino hard process.");


  static Switch<MEPP2NeutralinoNeutralinoPowheg,int> interfaceProcess
    ("Process",
     "Which processes to generate",
     &MEPP2NeutralinoNeutralinoPowheg::process_, 0, false, false);
  static SwitchOption interfaceProcessAll
    (interfaceProcess,
     "All",
     "Generate all the processes"
     " (i.e. both chargino and neutralino pair production"
     " of all mass eigenstates.)",
     0);
  static SwitchOption interfaceProcessNeutralino11
    (interfaceProcess,
     "chi01chi01",
     "Only produce chi01, chi01 pairs.",
     1);
  static SwitchOption interfaceProcessNeutralino12
    (interfaceProcess,
     "chi01chi02",
     "Only produce chi01, chi02 pairs.",
     2);
  static SwitchOption interfaceProcessNeutralino13
    (interfaceProcess,
     "chi01chi03",
     "Only produce chi01, chi03 pairs.",
     3);
  static SwitchOption interfaceProcessNeutralino14
    (interfaceProcess,
     "chi01chi04",
     "Only produce chi01, chi04 pairs.",
     4);
  static SwitchOption interfaceProcessNeutralino22
    (interfaceProcess,
     "chi02chi02",
     "Only produce chi02, chi02 pairs.",
     6);
  static SwitchOption interfaceProcessNeutralino23
    (interfaceProcess,
     "chi02chi03",
     "Only produce chi02, chi03 pairs.",
     7);
  static SwitchOption interfaceProcessNeutralino24
    (interfaceProcess,
     "chi02chi04",
     "Only produce chi02, chi04 pairs.",
     8);
  static SwitchOption interfaceProcessNeutralino33
    (interfaceProcess,
     "chi03chi03",
     "Only produce chi03, chi04 pairs.",
     11);
  static SwitchOption interfaceProcessNeutralino34
    (interfaceProcess,
     "chi03chi04",
     "Only produce chi03, chi04 pairs.",
     12);
  static SwitchOption interfaceProcessNeutralino44
    (interfaceProcess,
     "chi04chi04",
     "Only produce chi04, chi04 pairs.",
     16);

  static Parameter<MEPP2NeutralinoNeutralinoPowheg,int> interfaceMaxFlavour
    ("MaxFlavour",
     "The maximum flavour of the incoming quarks",
     &MEPP2NeutralinoNeutralinoPowheg::maxFlavour_, 5, 1, 5,
     false, false, Interface::limited);

}


Selector<MEBase::DiagramIndex>
MEPP2NeutralinoNeutralinoPowheg::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) {
    if ( diags[i]->id() == -1 )      sel.insert(meInfo()[0], i);
    else if ( diags[i]->id() == -2 ) sel.insert(meInfo()[1], i);
    else if ( diags[i]->id() == -3 ) sel.insert(meInfo()[2], i);
    else if ( diags[i]->id() == -4 ) sel.insert(meInfo()[3], i);
    else if ( diags[i]->id() == -5 ) sel.insert(meInfo()[4], i);
  }
  return sel;
}

NLODrellYanBase::Singular MEPP2NeutralinoNeutralinoPowheg::virtualME() const {
  Singular output;
  output.eps2 = -2;
  output.eps1 = -3;
  // cut-off parameter
  Energy2 eps = 0.1*MeV2;
  // outgoing masses
  Energy m1 = meMomenta()[2].mass();
  Energy m2 = meMomenta()[3].mass();
  // average of outgoing masses
  Energy2 mav2 = 0.25*sqr(m1+m2);
  // renormalisation/factorization scale
  Energy2 scale2 = scale();
  // gluino mass
  Energy mg  = getParticleData(ParticleID::SUSY_g)->mass();
  Energy2 mg2(sqr(mg));
  // average of left/right squark masses
  Energy  ms = 0.5*(getParticleData(1000000+mePartonData()[0]->id())->mass()+
		    getParticleData(2000000+mePartonData()[0]->id())->mass());
  Energy2 ms2(sqr(ms));



  Energy mz = Z0_->mass();
  Energy2 t1 = tHat()-sqr(m1), t2 = tHat()-sqr(m2);
  Energy2 u1 = uHat()-sqr(m1), u2 = uHat()-sqr(m2);
  Energy2 sz = sHat()-sqr(mz);
  

//       QBB_ss = 16.*Cs(1)*(t1*t2+u1*u2)/sqr(sHat())+ 32.*Cs(2)*m1*m2/sHat()
//               +32.*Cs(3)*(u1*u2- t1*t2)/sHat()/sz





//       QBB_st =          s**(-1)*ts**(-1) * (
//      +     + 16*Cl(1)*Clc(2)*Ct(1)*t1*t2
//      +     + 16*Cr(1)*Crc(2)*Ct(3)*t1*t2
//      +     )
//       QBB_st = QBB_st + ts**(-1) * (
//      +     + 16*Cl(1)*Clc(2)*Ct(2)*m1*m2
//      +     + 16*Cr(1)*Crc(2)*Ct(4)*m1*m2
//      +     )

//       QBB_su =          s**(-1)*us**(-1) * (
//      +     - 16*Cl(4)*Clc(3)*Ct(2)*u1*u2
//      +     - 16*Cr(4)*Crc(3)*Ct(4)*u1*u2
//      +     )
//       QBB_su = QBB_su + us**(-1) * (
//      +     - 16*Cl(4)*Clc(3)*Ct(1)*m1*m2
//      +     - 16*Cr(4)*Crc(3)*Ct(3)*m1*m2
//      +     )

//       QBB_txs =         s**(-1)*ts**(-1) * (
//      +      + 16*Cl(2)*Clc(1)*Ctc(1)*t1*t2
//      +      + 16*Cr(2)*Crc(1)*Ctc(3)*t1*t2
//      +     )

//       QBB_txt =         ts**(-2) * (
//      +      + 32*Cl(1)*Cl(2)*Clc(1)*Clc(2)*t1*t2
//      +      + 32*Cr(1)*Cr(2)*Crc(1)*Crc(2)*t1*t2
//      +     )

//       QBB_tmu =         s*ts**(-1)*us**(-1) * (
//      +      - 32*Cl(2)*Cl(4)*Clc(1)*Clc(3)*m1*m2
//      +      - 32*Cr(2)*Cr(4)*Crc(1)*Crc(3)*m1*m2
//      +     )

//       QBB_tms =         ts**(-1) * (
//      +      + 16*Cl(2)*Clc(1)*Ctc(2)*m1*m2
//      +      + 16*Cr(2)*Crc(1)*Ctc(4)*m1*m2
//      +     )

//       QBB_uxs =         s**(-1)*us**(-1) * (
//      +      - 16*Cl(3)*Clc(4)*Ctc(2)*u1*u2
//      +      - 16*Cr(3)*Crc(4)*Ctc(4)*u1*u2
//      +     )

//       QBB_uxu =         us**(-2) * (
//      +      + 32*Cl(3)*Cl(4)*Clc(3)*Clc(4)*u1*u2
//      +      + 32*Cr(3)*Cr(4)*Crc(3)*Crc(4)*u1*u2
//      +     )

//       QBB_umt =         s*ts**(-1)*us**(-1) * (
//      +      - 32*Cl(1)*Cl(3)*Clc(2)*Clc(4)*m1*m2
//      +      - 32*Cr(1)*Cr(3)*Crc(2)*Crc(4)*m1*m2
//      +     )

//       QBB_ums =         us**(-1) * (
//      +      - 16*Cl(3)*Clc(4)*Ctc(1)*m1*m2
//      +      - 16*Cr(3)*Crc(4)*Ctc(3)*m1*m2
//      +     )


//          b_s   = Nc*real( QBB_ss  + QBB_st   + QBB_su )
//          b_tx  = Nc*real( QBB_txs + QBB_txt )
//          b_tm  = Nc*real( QBB_tmu + QBB_tms )
//          b_ux  = Nc*real( QBB_uxs + QBB_uxu )
//          b_um  = Nc*real( QBB_umt + QBB_ums )
//          b_tx4 = Nc*real( QBB_txs + QBB_txt )
//          b_tm4 = Nc*real( QBB_tmu + QBB_tms )
//          b_ux2 = Nc*real( QBB_uxs + QBB_uxu )
//          b_um2 = Nc*real( QBB_umt + QBB_ums )
// 





  // bit proportional to LO
  output.finite = -3.*log(sHat()/mav2) + sqr(log(sHat()/mav2)) -1. - Pi**2/6.

    -SusyLoopIntegral::B0P(ZERO  ,mg    ,ms    ,scale2)*(mg2-ms2);
  output.finite *= loWeight();


//       s  = massin(1)
//       t2 = massin(2)
//       u2 = massin(3)
//       t1 = massin(4)
//       u1 = massin(5)

//       Clot = Cv(1) 
//       Cupt = Cv(2) 
//       Cupu = Cv(3) 
//       Clou = Cv(4) 

//       t   = t2 + m2**2 
//       u   = u2 + m2**2 
//       ts  = t2 + m2**2 - ms**2 
//       us  = u2 + m2**2 - ms**2 
//       tg  = t2 + m2**2 - mg**2 
//       ug  = u2 + m2**2 - mg**2 

//       kaellen = s**2 + m1**4 + m2**4
//      &         - 2*( s*m1**2 + s*m2**2 + m1**2*m2**2 )

//       mav2 = (abs(m1)+abs(m2))**2/4.D0

c               finite part from virtual subtraction matrix element 
      finite = - 3.D0/2.D0 * Pi**2/6.D0 
     &         - 3.D0/2.D0 * log(s/mav2) 
     &         + log(s/mav2)**2 /2.D0


  /*

c               the born type structures 
      call BORN_PARTS_NN(iq,iout,massin,Cs,Ct,Cl,Cr,
     &                b_s,b_tx,b_tm,b_ux,b_um,b_tx4,b_tm4,b_ux2,b_um2)
      b_t  = b_tx  + b_tm
      b_u  = b_ux  + b_um

c               numerator: improves stability
c               denominator: only in scalar integrals and as ts,us
c                            numerator should not matter anymore
ctp      if (ms>20.D0*(abs(m1)+abs(m2))/2.D0) ms=1.D0

c               insert the form output
 


      QBV = QBV + SCB(1,1)*Cupt*b_t * (  - 4*m2*mg*t2**(-1) )
      QBV = QBV + SCB(1,1)*Clot*b_t * (  - 4*m1*mg*t1**(-1) )




      QBV = QBV + SCB(1,1)*b_t * (  - 4*m1**2*ts**(-1) - 4*m2**2*
     +    ts**(-1) + 4*mg**2*ts**(-1) + 4*s*ts**(-1) + 4*u*ts**(-1) )
      QBV = QBV + SCB(1,1)*b_um2 * (  - 4*s**(-1)*us )
      QBV = QBV + SCB(1,1)*b_ux2 * ( 4*m1**2*s*t1**(-1)*u1**(-1)*
     +    u2**(-1)*us + 2*m1**2*s*u1**(-1)*t2**(-1)*u2**(-1)*us + 2*
     +    m1**2*u1**(-1)*u2**(-1)*us + 2*m1**2*t2**(-1)*u2**(-1)*us + 2
     +    *m2**2*s*u1**(-1)*t2**(-1)*u2**(-1)*us + 2*m2**2*u1**(-1)*
     +    u2**(-1)*us - 2*m2**2*t2**(-1)*u2**(-1)*us - 4*u2**(-1)*us )
      QBV = QBV + SCB(1,2)*Cupu*b_u * (  - 4*m1*mg*u1**(-1) )
      QBV = QBV + SCB(1,2)*Clou*b_u * (  - 4*m2*mg*u2**(-1) )
      QBV = QBV + SCB(1,2)*b_u * ( 4*mg**2*us**(-1) - 4*u*us**(-1)
     +     )
      QBV = QBV + SCB(1,2)*b_tm4 * (  - 4*s**(-1)*ts )
      QBV = QBV + SCB(1,2)*b_tx4 * ( 2*m1**2*t1**(-1)*t2**(-1)*ts
     +     + 4*m2**2*s*t1**(-1)*u1**(-1)*t2**(-1)*ts + 4*m2**2*s*
     +    t1**(-1)*t2**(-1)*u2**(-1)*ts + 2*m2**2*t1**(-1)*t2**(-1)*ts
     +     + 2*s*t1**(-1)*u1**(-1)*ts - 2*s*t1**(-1)*t2**(-1)*ts - 2*s*
     +    u1**(-1)*t2**(-1)*ts - 2*s**2*t1**(-1)*u1**(-1)*t2**(-1)*ts
     +     - 2*t1**(-1)*ts + 2*u1**(-1)*ts )
      QBV = QBV + SCB(1,3)*b_t * (  - 4*mg**2*ts**(-1) + 4*ms**2*
     +    ts**(-1) )
      QBV = QBV + SCB(1,3)*b_u * (  - 4*mg**2*us**(-1) + 4*ms**2*
     +    us**(-1) )
      QBV = QBV + SCB(3,1)*b_t * (  - 4 - 4*m1**2*t1**(-1) + 4*
     +    m1**2*ts**(-1) - 4*m2**2*t2**(-1) + 4*m2**2*ts**(-1) + 4*
     +    ms**2*ts**(-1) - 4*s*ts**(-1) - 4*u*ts**(-1) )
      QBV = QBV + SCB(3,1)*b_tm * (  - 4*s**(-1)*ts )
      QBV = QBV + SCB(3,1)*b_tx * ( 3*m1**2*t1**(-1)*t2**(-1)*ts
     +     + m2**2*t1**(-1)*t2**(-1)*ts + t1**(-1)*ts + 3*t2**(-1)*ts )
      QBV = QBV + SCB(3,2)*b_u * (  - 4 - 4*m1**2*u1**(-1) - 4*
     +    m2**2*u2**(-1) + 4*ms**2*us**(-1) + 4*u*us**(-1) )
      QBV = QBV + SCB(3,2)*b_um * (  - 4*s**(-1)*us )
      QBV = QBV + SCB(3,2)*b_ux * ( 3*m1**2*u1**(-1)*u2**(-1)*us
     +     + m2**2*u1**(-1)*u2**(-1)*us + u1**(-1)*us + 3*u2**(-1)*us )
      QBV = QBV + SCB(3,3)*b_t * (  - 8*ms**2*ts**(-1) )
      QBV = QBV + SCB(3,3)*b_u * (  - 8*ms**2*us**(-1) )
      QBV = QBV + SCB(3,4)*Cupu*b_u * ( 4*m1*mg*u1**(-1) )
      QBV = QBV + SCB(3,4)*Clot*b_t * ( 4*m1*mg*t1**(-1) )
      QBV = QBV + SCB(3,4)*kaellen**(-1)*b_tm * ( 4*m1**2*m2**2*
     +    s**(-1)*ts - 4*m1**2*s**(-1)*u*ts + 4*m1**2*ts + 4*m2**2*
     +    s**(-1)*u*ts + 8*m2**2*ts - 4*m2**4*s**(-1)*ts - 4*s*ts - 4*u
     +    *ts )
      QBV = QBV + SCB(3,4)*kaellen**(-1)*b_tx * ( 6*m1**2*m2**2*s*
     +    t1**(-1)*t2**(-1)*ts + 4*m1**2*m2**2*t2**(-1)*ts - 4*m1**2*
     +    m2**4*t1**(-1)*t2**(-1)*ts + 2*m1**2*s*t1**(-1)*ts + 2*m1**2*
     +    s*t2**(-1)*ts + 8*m1**4*m2**2*t1**(-1)*t2**(-1)*ts + 2*m1**4*
     +    s*t1**(-1)*t2**(-1)*ts - 4*m1**4*t2**(-1)*ts - 4*m1**6*
     +    t1**(-1)*t2**(-1)*ts )
      QBV = QBV + SCB(3,4)*kaellen**(-1)*b_um * (  - 4*m1**2*m2**2
     +    *s**(-1)*us + 4*m1**2*s**(-1)*u*us - 4*m1**2*us - 4*m2**2*
     +    s**(-1)*u*us - 8*m2**2*us + 4*m2**4*s**(-1)*us + 4*s*us + 4*u
     +    *us )
      QBV = QBV + SCB(3,4)*kaellen**(-1)*b_ux * ( 10*m1**2*m2**2*s
     +    *u1**(-1)*u2**(-1)*us + 2*m1**2*m2**2*u1**(-1)*us + 2*m1**2*
     +    m2**2*u2**(-1)*us - 3*m1**2*m2**4*u1**(-1)*u2**(-1)*us + 4*
     +    m1**2*s*u1**(-1)*us + 3*m1**4*m2**2*u1**(-1)*u2**(-1)*us - 3*
     +    m1**4*u1**(-1)*us - m1**4*u2**(-1)*us - m1**6*u1**(-1)*
     +    u2**(-1)*us - 2*m2**2*s*u1**(-1)*us + 2*m2**2*s*u2**(-1)*us
     +     - 2*m2**4*s*u1**(-1)*u2**(-1)*us + m2**4*u1**(-1)*us - m2**4
     +    *u2**(-1)*us + m2**6*u1**(-1)*u2**(-1)*us )
      QBV = QBV + SCB(3,4)*kaellen**(-1)*b_um2 * ( 4*m1**2*m2**2*
     +    s**(-1)*us - 4*m1**2*s**(-1)*u*us + 4*m1**2*us + 4*m2**2*
     +    s**(-1)*u*us + 8*m2**2*us - 4*m2**4*s**(-1)*us - 4*s*us - 4*u
     +    *us )
      QBV = QBV + SCB(3,4)*kaellen**(-1)*b_ux2 * (  - 4*m1**2*
     +    m2**2*s*u1**(-1)*u2**(-1)*us - 4*m1**2*m2**2*u2**(-1)*us + 6*
     +    m1**2*m2**4*u1**(-1)*u2**(-1)*us + 4*m1**2*s*u1**(-1)*us - 8*
     +    m1**2*s*u2**(-1)*us - 6*m1**4*m2**2*u1**(-1)*u2**(-1)*us - 8*
     +    m1**4*s*u1**(-1)*u2**(-1)*us + 2*m1**4*u1**(-1)*us + 2*m1**4*
     +    u2**(-1)*us + 2*m1**6*u1**(-1)*u2**(-1)*us + 4*m2**2*s*
     +    u1**(-1)*us - 4*m2**2*s*u2**(-1)*us + 4*m2**4*s*u1**(-1)*
     +    u2**(-1)*us - 2*m2**4*u1**(-1)*us + 2*m2**4*u2**(-1)*us - 2*
     +    m2**6*u1**(-1)*u2**(-1)*us )
      QBV = QBV + SCB(3,4)*kaellen**(-1)*b_tm4 * (  - 4*m1**2*
     +    m2**2*s**(-1)*ts + 4*m1**2*s**(-1)*u*ts - 4*m1**2*ts - 4*
     +    m2**2*s**(-1)*u*ts - 8*m2**2*ts + 4*m2**4*s**(-1)*ts + 4*s*ts
     +     + 4*u*ts )
      QBV = QBV + SCB(3,4)*kaellen**(-1)*b_tx4 * (  - 8*m1**2*
     +    m2**2*s*t1**(-1)*t2**(-1)*ts - 4*m1**2*m2**2*t2**(-1)*ts + 4*
     +    m1**2*m2**4*t1**(-1)*t2**(-1)*ts + 2*m1**2*s*t1**(-1)*ts - 6*
     +    m1**2*s*t2**(-1)*ts - 8*m1**4*m2**2*t1**(-1)*t2**(-1)*ts - 6*
     +    m1**4*s*t1**(-1)*t2**(-1)*ts + 4*m1**4*t2**(-1)*ts + 4*m1**6*
     +    t1**(-1)*t2**(-1)*ts + 6*m2**2*s*t1**(-1)*ts - 6*m2**2*s*
     +    t2**(-1)*ts + 6*m2**4*s*t1**(-1)*t2**(-1)*ts )
      QBV = QBV + SCB(3,4)*b_t * ( 4*m1**2*t1**(-1) )
      QBV = QBV + SCB(3,4)*b_u * ( 4*m1**2*u1**(-1) )
      QBV = QBV + SCB(3,4)*b_tm * ( 4*s**(-1)*ts )
      QBV = QBV + SCB(3,4)*b_tx * ( m1**2*t1**(-1)*t2**(-1)*ts - 
     +    m2**2*t1**(-1)*t2**(-1)*ts - t1**(-1)*ts + t2**(-1)*ts )
      QBV = QBV + SCB(3,4)*b_ux * ( m1**2*u1**(-1)*u2**(-1)*us - 
     +    m2**2*u1**(-1)*u2**(-1)*us - u1**(-1)*us + u2**(-1)*us )
      QBV = QBV + SCB(3,4)*b_um2 * ( 4*s**(-1)*us )
      QBV = QBV + SCB(3,4)*b_ux2 * (  - 4*m1**2*s*t1**(-1)*
     +    u1**(-1)*u2**(-1)*us - 4*m1**2*u1**(-1)*u2**(-1)*us )
      QBV = QBV + SCB(3,4)*b_tx4 * ( 2*m1**2*t1**(-1)*u1**(-1)*ts
     +     - 4*m1**2*t1**(-1)*t2**(-1)*ts - 4*m2**2*s*t1**(-1)*u1**(-1)
     +    *t2**(-1)*ts - 2*m2**2*t1**(-1)*u1**(-1)*ts + 2*s*t1**(-1)*
     +    t2**(-1)*ts + 2*s*u1**(-1)*t2**(-1)*ts + 2*s**2*t1**(-1)*
     +    u1**(-1)*t2**(-1)*ts + 2*t1**(-1)*ts - 2*t2**(-1)*ts )
      QBV = QBV + SCB(3,5)*Cupt*b_t * ( 4*m2*mg*t2**(-1) )
      QBV = QBV + SCB(3,5)*Clou*b_u * ( 4*m2*mg*u2**(-1) )
      QBV = QBV + SCB(3,5)*kaellen**(-1)*b_tm * (  - 4*m1**2*m2**2
     +    *s**(-1)*ts + 4*m1**2*s**(-1)*u*ts - 4*m2**2*s**(-1)*u*ts - 4
     +    *m2**2*ts + 4*m2**4*s**(-1)*ts - 4*u*ts )
      QBV = QBV + SCB(3,5)*kaellen**(-1)*b_tx * ( 6*m1**2*m2**2*s*
     +    t1**(-1)*t2**(-1)*ts + 6*m1**2*m2**2*t1**(-1)*ts - 2*m1**2*
     +    m2**2*t2**(-1)*ts + 9*m1**2*m2**4*t1**(-1)*t2**(-1)*ts + 2*
     +    m1**2*s*t1**(-1)*ts - 2*m1**2*s*t2**(-1)*ts - 9*m1**4*m2**2*
     +    t1**(-1)*t2**(-1)*ts - 2*m1**4*s*t1**(-1)*t2**(-1)*ts - 3*
     +    m1**4*t1**(-1)*ts + 3*m1**4*t2**(-1)*ts + 3*m1**6*t1**(-1)*
     +    t2**(-1)*ts + 4*m2**2*s*t1**(-1)*ts + 4*m2**4*s*t1**(-1)*
     +    t2**(-1)*ts - 3*m2**4*t1**(-1)*ts - m2**4*t2**(-1)*ts - 3*
     +    m2**6*t1**(-1)*t2**(-1)*ts )
      QBV = QBV + SCB(3,5)*kaellen**(-1)*b_um * ( 4*m1**2*m2**2*
     +    s**(-1)*us - 4*m1**2*s**(-1)*u*us + 4*m2**2*s**(-1)*u*us + 4*
     +    m2**2*us - 4*m2**4*s**(-1)*us + 4*u*us )
      QBV = QBV + SCB(3,5)*kaellen**(-1)*b_ux * ( 2*m1**2*m2**2*s*
     +    u1**(-1)*u2**(-1)*us + 4*m1**2*m2**2*u1**(-1)*us + 8*m1**2*
     +    m2**4*u1**(-1)*u2**(-1)*us - 4*m1**4*m2**2*u1**(-1)*u2**(-1)*
     +    us + 6*m2**2*s*u1**(-1)*us - 2*m2**2*s*u2**(-1)*us + 6*m2**4*
     +    s*u1**(-1)*u2**(-1)*us - 4*m2**4*u1**(-1)*us - 4*m2**6*
     +    u1**(-1)*u2**(-1)*us )
      QBV = QBV + SCB(3,5)*kaellen**(-1)*b_um2 * (  - 4*m1**2*
     +    m2**2*s**(-1)*us + 4*m1**2*s**(-1)*u*us - 4*m2**2*s**(-1)*u*
     +    us - 4*m2**2*us + 4*m2**4*s**(-1)*us - 4*u*us )
      QBV = QBV + SCB(3,5)*kaellen**(-1)*b_ux2 * (  - 8*m1**2*
     +    m2**2*s*u1**(-1)*u2**(-1)*us - 4*m1**2*m2**2*u1**(-1)*us - 8*
     +    m1**2*m2**4*u1**(-1)*u2**(-1)*us - 6*m1**2*s*u1**(-1)*us + 6*
     +    m1**2*s*u2**(-1)*us + 4*m1**4*m2**2*u1**(-1)*u2**(-1)*us + 6*
     +    m1**4*s*u1**(-1)*u2**(-1)*us - 6*m2**2*s*u1**(-1)*us + 2*
     +    m2**2*s*u2**(-1)*us - 6*m2**4*s*u1**(-1)*u2**(-1)*us + 4*
     +    m2**4*u1**(-1)*us + 4*m2**6*u1**(-1)*u2**(-1)*us )
      QBV = QBV + SCB(3,5)*kaellen**(-1)*b_tm4 * ( 4*m1**2*m2**2*
     +    s**(-1)*ts - 4*m1**2*s**(-1)*u*ts + 4*m2**2*s**(-1)*u*ts + 4*
     +    m2**2*ts - 4*m2**4*s**(-1)*ts + 4*u*ts )
      QBV = QBV + SCB(3,5)*kaellen**(-1)*b_tx4 * (  - 4*m1**2*
     +    m2**2*s*t1**(-1)*t2**(-1)*ts - 4*m1**2*m2**2*t1**(-1)*ts - 6*
     +    m1**2*m2**4*t1**(-1)*t2**(-1)*ts - 4*m1**2*s*t1**(-1)*ts + 4*
     +    m1**2*s*t2**(-1)*ts + 6*m1**4*m2**2*t1**(-1)*t2**(-1)*ts + 4*
     +    m1**4*s*t1**(-1)*t2**(-1)*ts + 2*m1**4*t1**(-1)*ts - 2*m1**4*
     +    t2**(-1)*ts - 2*m1**6*t1**(-1)*t2**(-1)*ts - 8*m2**2*s*
     +    t1**(-1)*ts + 4*m2**2*s*t2**(-1)*ts - 8*m2**4*s*t1**(-1)*
     +    t2**(-1)*ts + 2*m2**4*t1**(-1)*ts + 2*m2**4*t2**(-1)*ts + 2*
     +    m2**6*t1**(-1)*t2**(-1)*ts )
      QBV = QBV + SCB(3,5)*b_t * ( 4*m2**2*t2**(-1) )
      QBV = QBV + SCB(3,5)*b_u * ( 4*m2**2*u2**(-1) )
      QBV = QBV + SCB(3,5)*b_tx * (  - 3*m1**2*t1**(-1)*t2**(-1)*
     +    ts + 3*m2**2*t1**(-1)*t2**(-1)*ts + 3*t1**(-1)*ts - 3*
     +    t2**(-1)*ts )
      QBV = QBV + SCB(3,5)*b_um * ( 4*s**(-1)*us )
      QBV = QBV + SCB(3,5)*b_ux * (  - 3*m1**2*u1**(-1)*u2**(-1)*
     +    us + 3*m2**2*u1**(-1)*u2**(-1)*us + 3*u1**(-1)*us - 3*
     +    u2**(-1)*us )
      QBV = QBV + SCB(3,5)*b_ux2 * (  - 2*m1**2*s*u1**(-1)*
     +    t2**(-1)*u2**(-1)*us - 2*m1**2*t2**(-1)*u2**(-1)*us - 2*m2**2
     +    *s*u1**(-1)*t2**(-1)*u2**(-1)*us - 4*m2**2*u1**(-1)*u2**(-1)*
     +    us + 2*m2**2*t2**(-1)*u2**(-1)*us - 2*u1**(-1)*us + 2*
     +    u2**(-1)*us )
      QBV = QBV + SCB(3,5)*b_tm4 * ( 4*s**(-1)*ts )
      QBV = QBV + SCB(3,5)*b_tx4 * (  - 4*m2**2*s*t1**(-1)*
     +    t2**(-1)*u2**(-1)*ts - 4*m2**2*t1**(-1)*t2**(-1)*ts )
      QBV = QBV + SCB(4,1)*kaellen**(-1)*b_um2 * (  - 4*m1**2*us
     +     - 4*m2**2*us + 4*s*us + 8*u*us )
      QBV = QBV + SCB(4,1)*kaellen**(-1)*b_ux2 * ( 12*m1**2*m2**2*
     +    s*u1**(-1)*u2**(-1)*us + 4*m1**2*m2**2*u1**(-1)*us + 4*m1**2*
     +    m2**2*u2**(-1)*us + 2*m1**2*m2**4*u1**(-1)*u2**(-1)*us + 2*
     +    m1**2*s*u1**(-1)*us + 2*m1**2*s*u2**(-1)*us + 2*m1**4*m2**2*
     +    u1**(-1)*u2**(-1)*us + 2*m1**4*s*u1**(-1)*u2**(-1)*us - 2*
     +    m1**4*u1**(-1)*us - 2*m1**4*u2**(-1)*us - 2*m1**6*u1**(-1)*
     +    u2**(-1)*us + 2*m2**2*s*u1**(-1)*us + 2*m2**2*s*u2**(-1)*us
     +     + 2*m2**4*s*u1**(-1)*u2**(-1)*us - 2*m2**4*u1**(-1)*us - 2*
     +    m2**4*u2**(-1)*us - 2*m2**6*u1**(-1)*u2**(-1)*us )
      QBV = QBV + SCB(4,1)*kaellen**(-1)*b_tm4 * ( 4*m1**2*ts + 4*
     +    m2**2*ts - 4*s*ts - 8*u*ts )
      QBV = QBV + SCB(4,1)*kaellen**(-1)*b_tx4 * ( 12*m1**2*m2**2*
     +    s*t1**(-1)*t2**(-1)*ts + 4*m1**2*m2**2*t1**(-1)*ts + 4*m1**2*
     +    m2**2*t2**(-1)*ts + 2*m1**2*m2**4*t1**(-1)*t2**(-1)*ts + 2*
     +    m1**2*s*t1**(-1)*ts + 2*m1**2*s*t2**(-1)*ts + 2*m1**4*m2**2*
     +    t1**(-1)*t2**(-1)*ts + 2*m1**4*s*t1**(-1)*t2**(-1)*ts - 2*
     +    m1**4*t1**(-1)*ts - 2*m1**4*t2**(-1)*ts - 2*m1**6*t1**(-1)*
     +    t2**(-1)*ts + 2*m2**2*s*t1**(-1)*ts + 2*m2**2*s*t2**(-1)*ts
     +     + 2*m2**4*s*t1**(-1)*t2**(-1)*ts - 2*m2**4*t1**(-1)*ts - 2*
     +    m2**4*t2**(-1)*ts - 2*m2**6*t1**(-1)*t2**(-1)*ts )
      QBV = QBV + SCB(4,1)*b_s * ( 2 + 4*mg**2*s**(-1) - 4*ms**2*
     +    s**(-1) )
      QBV = QBV + SCB(4,1)*b_ux2 * ( 2*m1**2*u1**(-1)*u2**(-1)*us
     +     + 2*m2**2*u1**(-1)*u2**(-1)*us + 2*u1**(-1)*us + 2*u2**(-1)*
     +    us )
      QBV = QBV + SCB(4,1)*b_tx4 * ( 2*m1**2*t1**(-1)*t2**(-1)*ts
     +     + 2*m2**2*t1**(-1)*t2**(-1)*ts + 2*t1**(-1)*ts + 2*t2**(-1)*
     +    ts )
      QBV = QBV + SCB(6,1)*b_s * (  - 2 - 4*(mg**2-ms**2)*s**(-1) )
      QBV = QBV + SCB(6,1)*b_t * (  - 2 )
      QBV = QBV + SCB(6,1)*b_u * (  - 2 )

      QBV = QBV + SCB(6,1)*b_tx4 * (  - 2*m1**2*t1**(-1)*u1**(-1)*
     +    ts + 2*m2**2*t1**(-1)*u1**(-1)*ts - 2*s*t1**(-1)*u1**(-1)*ts
     +     - 2*t1**(-1)*ts - 2*u1**(-1)*ts )
      QBV = QBV + SCB(7,1)*kaellen**(-1)*b_tm * (  - 4*m1**2*ts - 
     +    4*m2**2*ts + 4*s*ts + 8*u*ts )
      QBV = QBV + SCB(7,1)*kaellen**(-1)*b_tx * (  - 12*m1**2*
     +    m2**2*s*t1**(-1)*t2**(-1)*ts - 6*m1**2*m2**2*t1**(-1)*ts - 2*
     +    m1**2*m2**2*t2**(-1)*ts - 5*m1**2*m2**4*t1**(-1)*t2**(-1)*ts
     +     - 4*m1**2*s*t1**(-1)*ts + m1**4*m2**2*t1**(-1)*t2**(-1)*ts
     +     + 3*m1**4*t1**(-1)*ts + m1**4*t2**(-1)*ts + m1**6*t1**(-1)*
     +    t2**(-1)*ts - 4*m2**2*s*t1**(-1)*ts - 4*m2**4*s*t1**(-1)*
     +    t2**(-1)*ts + 3*m2**4*t1**(-1)*ts + m2**4*t2**(-1)*ts + 3*
     +    m2**6*t1**(-1)*t2**(-1)*ts )
      QBV = QBV + SCB(7,1)*kaellen**(-1)*b_um * ( 4*m1**2*us + 4*
     +    m2**2*us - 4*s*us - 8*u*us )
      QBV = QBV + SCB(7,1)*kaellen**(-1)*b_ux * (  - 12*m1**2*
     +    m2**2*s*u1**(-1)*u2**(-1)*us - 6*m1**2*m2**2*u1**(-1)*us - 2*
     +    m1**2*m2**2*u2**(-1)*us - 5*m1**2*m2**4*u1**(-1)*u2**(-1)*us
     +     - 4*m1**2*s*u1**(-1)*us + m1**4*m2**2*u1**(-1)*u2**(-1)*us
     +     + 3*m1**4*u1**(-1)*us + m1**4*u2**(-1)*us + m1**6*u1**(-1)*
     +    u2**(-1)*us - 4*m2**2*s*u1**(-1)*us - 4*m2**4*s*u1**(-1)*
     +    u2**(-1)*us + 3*m2**4*u1**(-1)*us + m2**4*u2**(-1)*us + 3*
     +    m2**6*u1**(-1)*u2**(-1)*us )
      QBV = QBV + SCB(7,1)*b_s * (  - 6 )
      QBV = QBV + SCB(7,1)*b_tx * (  - m1**2*t1**(-1)*t2**(-1)*ts
     +     - 3*m2**2*t1**(-1)*t2**(-1)*ts - 3*t1**(-1)*ts - t2**(-1)*ts
     +     )
      QBV = QBV + SCB(7,1)*b_ux * (  - m1**2*u1**(-1)*u2**(-1)*us
     +     - 3*m2**2*u1**(-1)*u2**(-1)*us - 3*u1**(-1)*us - u2**(-1)*us
     +     )

      QBV = QBV + SCC(1,1)*Clot*b_t * (  - 4*m1*mg*ms**2*s**(-1)*
     +    t1**(-1) - 4*m1*mg*s**(-1) + 4*m1*mg**3*s**(-1)*t1**(-1) )
      QBV = QBV + SCC(1,1)*b_um2 * ( 2*m1**2*s**(-2)*us + 4*mg**2*
     +    s**(-2)*us - 4*ms**2*s**(-2)*us - 2*s**(-2)*t*us )
      QBV = QBV + SCC(1,1)*b_ux2 * (  - 2*m1**2*mg**2*s**(-1)*
     +    u1**(-1)*u2**(-1)*us - 4*m1**2*mg**2*t1**(-1)*u1**(-1)*
     +    u2**(-1)*us + 2*m1**2*ms**2*s**(-1)*u1**(-1)*u2**(-1)*us + 4*
     +    m1**2*ms**2*t1**(-1)*u1**(-1)*u2**(-1)*us - 2*m2**2*mg**2*
     +    s**(-1)*u1**(-1)*u2**(-1)*us + 2*m2**2*ms**2*s**(-1)*u1**(-1)
     +    *u2**(-1)*us + 2*mg**2*s**(-1)*u2**(-1)*us + 2*mg**2*u1**(-1)
     +    *u2**(-1)*us - 2*ms**2*s**(-1)*u1**(-1)*us - 2*ms**2*s**(-1)*
     +    u2**(-1)*us - 4*ms**2*u1**(-1)*u2**(-1)*us + 2*s*u1**(-1)*
     +    u2**(-1)*us + 2*u1**(-1)*us )
      QBV = QBV + SCC(1,2)*Cupt*b_t * (  - 4*m2*mg*ms**2*s**(-1)*
     +    t2**(-1) - 4*m2*mg*s**(-1) + 4*m2*mg**3*s**(-1)*t2**(-1) )
      QBV = QBV + SCC(1,2)*b_um2 * ( 2*m2**2*s**(-2)*us + 4*mg**2*
     +    s**(-2)*us - 4*ms**2*s**(-2)*us - 2*s**(-2)*t*us )
      QBV = QBV + SCC(1,2)*b_ux2 * (  - 2*m1**2*mg**2*s**(-1)*
     +    u1**(-1)*u2**(-1)*us - 2*m1**2*mg**2*s**(-1)*t2**(-1)*
     +    u2**(-1)*us - 2*m1**2*mg**2*u1**(-1)*t2**(-1)*u2**(-1)*us + 2
     +    *m1**2*ms**2*s**(-1)*u1**(-1)*u2**(-1)*us + 2*m1**2*ms**2*
     +    s**(-1)*t2**(-1)*u2**(-1)*us + 2*m1**2*ms**2*u1**(-1)*
     +    t2**(-1)*u2**(-1)*us + 2*m1**2*s**(-1)*u2**(-1)*us - 2*m2**2*
     +    mg**2*s**(-1)*u1**(-1)*u2**(-1)*us + 2*m2**2*mg**2*s**(-1)*
     +    t2**(-1)*u2**(-1)*us - 2*m2**2*mg**2*u1**(-1)*t2**(-1)*
     +    u2**(-1)*us + 2*m2**2*ms**2*s**(-1)*u1**(-1)*u2**(-1)*us - 2*
     +    m2**2*ms**2*s**(-1)*t2**(-1)*u2**(-1)*us + 2*m2**2*ms**2*
     +    u1**(-1)*t2**(-1)*u2**(-1)*us + 2*mg**2*s**(-1)*u2**(-1)*us
     +     + 2*mg**2*u1**(-1)*u2**(-1)*us - 4*ms**2*s**(-1)*u2**(-1)*us
     +     - 4*ms**2*u1**(-1)*u2**(-1)*us - 2*s**(-1)*t*u2**(-1)*us - 2
     +    *s**(-1)*us + 2*s*u1**(-1)*u2**(-1)*us )
      QBV = QBV + SCC(1,3)*Cupu*b_u * (  - 4*m1*mg*ms**2*s**(-1)*
     +    u1**(-1) - 4*m1*mg*s**(-1) + 4*m1*mg**3*s**(-1)*u1**(-1) )
      QBV = QBV + SCC(1,3)*b_tm4 * ( 2*m1**2*s**(-2)*ts + 4*mg**2*
     +    s**(-2)*ts - 4*ms**2*s**(-2)*ts - 2*s**(-2)*u*ts )
      QBV = QBV + SCC(1,3)*b_tx4 * (  - 2*m1**2*m2**2*s**(-1)*
     +    t1**(-1)*u1**(-1)*ts + 2*m1**2*mg**2*s**(-1)*t1**(-1)*
     +    u1**(-1)*ts - 2*m1**2*mg**2*s**(-1)*t1**(-1)*t2**(-1)*ts + 2*
     +    m1**2*ms**2*s**(-1)*t1**(-1)*t2**(-1)*ts + 2*m1**2*s**(-1)*
     +    t1**(-1)*ts + 2*m1**2*s**(-1)*u1**(-1)*ts + 2*m1**2*t1**(-1)*
     +    u1**(-1)*ts + 2*m1**4*s**(-1)*t1**(-1)*u1**(-1)*ts - 2*m2**2*
     +    mg**2*s**(-1)*t1**(-1)*u1**(-1)*ts - 2*m2**2*mg**2*s**(-1)*
     +    t1**(-1)*t2**(-1)*ts - 4*m2**2*mg**2*t1**(-1)*u1**(-1)*
     +    t2**(-1)*ts + 2*m2**2*ms**2*s**(-1)*t1**(-1)*t2**(-1)*ts + 4*
     +    m2**2*ms**2*t1**(-1)*u1**(-1)*t2**(-1)*ts + 2*m2**2*s**(-1)*
     +    t1**(-1)*ts + 2*mg**2*s**(-1)*t1**(-1)*ts + 2*mg**2*s*
     +    t1**(-1)*u1**(-1)*t2**(-1)*ts + 4*mg**2*t1**(-1)*t2**(-1)*ts
     +     + 2*mg**2*u1**(-1)*t2**(-1)*ts - 2*ms**2*s**(-1)*t1**(-1)*ts
     +     + 2*ms**2*s**(-1)*u1**(-1)*ts - 2*ms**2*s*t1**(-1)*u1**(-1)*
     +    t2**(-1)*ts + 2*ms**2*t1**(-1)*u1**(-1)*ts - 6*ms**2*t1**(-1)
     +    *t2**(-1)*ts )
      QBV = QBV + SCC(1,3)*b_tx4 * (  - 2*ms**2*u1**(-1)*t2**(-1)*
     +    ts - 2*s**(-1)*u*t1**(-1)*ts - 2*s**(-1)*ts + 2*s*t1**(-1)*
     +    t2**(-1)*ts )
      QBV = QBV + SCC(1,4)*Clou*b_u * (  - 4*m2*mg*ms**2*s**(-1)*
     +    u2**(-1) - 4*m2*mg*s**(-1) + 4*m2*mg**3*s**(-1)*u2**(-1) )
      QBV = QBV + SCC(1,4)*b_tm4 * ( 2*m2**2*s**(-2)*ts + 4*mg**2*
     +    s**(-2)*ts - 4*ms**2*s**(-2)*ts - 2*s**(-2)*u*ts )
      QBV = QBV + SCC(1,4)*b_tx4 * (  - 2*m1**2*mg**2*s**(-1)*
     +    t1**(-1)*t2**(-1)*ts + 2*m1**2*ms**2*s**(-1)*t1**(-1)*
     +    t2**(-1)*ts - 2*m2**2*mg**2*s**(-1)*t1**(-1)*t2**(-1)*ts - 4*
     +    m2**2*mg**2*t1**(-1)*t2**(-1)*u2**(-1)*ts + 2*m2**2*ms**2*
     +    s**(-1)*t1**(-1)*t2**(-1)*ts + 4*m2**2*ms**2*t1**(-1)*
     +    t2**(-1)*u2**(-1)*ts + 2*mg**2*s**(-1)*t1**(-1)*ts + 2*mg**2*
     +    t1**(-1)*t2**(-1)*ts - 2*ms**2*s**(-1)*t1**(-1)*ts - 2*ms**2*
     +    s**(-1)*t2**(-1)*ts - 4*ms**2*t1**(-1)*t2**(-1)*ts + 2*s*
     +    t1**(-1)*t2**(-1)*ts + 2*t2**(-1)*ts )
      QBV = QBV + SCC(2,1)*b_t * (  - 4*m1**2*t1**(-1) + 4*ms**2*
     +    t1**(-1) )
      QBV = QBV + SCC(2,1)*b_tm * ( 4*t1**(-1)*ts )
      QBV = QBV + SCC(2,1)*b_tx * ( m1**2*t1**(-1)*t2**(-1)*ts + 
     +    m2**2*t1**(-1)*t2**(-1)*ts - 2*ms**2*t1**(-1)*t2**(-1)*ts + 3
     +    *t1**(-1)*ts - t2**(-1)*ts )
      QBV = QBV + SCC(2,2)*b_t * (  - 4*m2**2*t2**(-1) + 4*ms**2*
     +    t2**(-1) )
      QBV = QBV + SCC(2,2)*b_tm * ( 4*t2**(-1)*ts )
      QBV = QBV + SCC(2,2)*b_tx * ( 2*m1**2*t1**(-1)*t2**(-1)*ts
     +     - 2*ms**2*t1**(-1)*t2**(-1)*ts - 2*t1**(-1)*ts + 4*t2**(-1)*
     +    ts )
      QBV = QBV + SCC(2,3)*b_u * (  - 4*m1**2*u1**(-1) + 4*ms**2*
     +    u1**(-1) )
      QBV = QBV + SCC(2,3)*b_um * ( 4*u1**(-1)*us )
      QBV = QBV + SCC(2,3)*b_ux * (  - m1**2*u1**(-1)*u2**(-1)*us
     +     + 3*m2**2*u1**(-1)*u2**(-1)*us - 2*ms**2*u1**(-1)*u2**(-1)*
     +    us + 5*u1**(-1)*us - 3*u2**(-1)*us )
      QBV = QBV + SCC(2,4)*b_u * (  - 4*m2**2*u2**(-1) + 4*ms**2*
     +    u2**(-1) )
      QBV = QBV + SCC(2,4)*b_um * ( 4*u2**(-1)*us )
      QBV = QBV + SCC(2,4)*b_ux * ( 2*m2**2*u1**(-1)*u2**(-1)*us
     +     - 2*ms**2*u1**(-1)*u2**(-1)*us + 2*u2**(-1)*us )
      QBV = QBV + SCC(3,1)*b_s * (  - 8*mg**2*ms**2*s**(-2) + 4*
     +    mg**2*s**(-1) + 4*mg**4*s**(-2) + 4*ms**4*s**(-2) )
      QBV = QBV + SCC(3,1)*b_um2 * (  - 2*s**(-1)*us )
      QBV = QBV + SCC(3,1)*b_ux2 * (  - 2*mg**2*u1**(-1)*u2**(-1)*
     +    us + 4*ms**2*u1**(-1)*u2**(-1)*us - 2*s*u1**(-1)*u2**(-1)*us
     +     )
      QBV = QBV + SCC(3,1)*b_tm4 * (  - 2*s**(-1)*ts )
      QBV = QBV + SCC(3,1)*b_tx4 * (  - 2*mg**2*t1**(-1)*t2**(-1)*
     +    ts + 4*ms**2*t1**(-1)*t2**(-1)*ts - 2*s*t1**(-1)*t2**(-1)*ts
     +     )
      QBV = QBV + SCC(4,1)*b_s * (  - 4 )
      QBV = QBV + SCC(4,1)*b_tx * ( 2*m1**2*t1**(-1)*t2**(-1)*ts
     +     - 2*ms**2*t1**(-1)*t2**(-1)*ts - 2*t1**(-1)*ts )
      QBV = QBV + SCC(4,1)*b_ux * ( 2*m2**2*u1**(-1)*u2**(-1)*us
     +     - 2*ms**2*u1**(-1)*u2**(-1)*us - 2*u2**(-1)*us )
      QBV = QBV + SCC(5,1)*kaellen**(-1)*b_tm * ( 8*m1**2*m2**2*
     +    s**(-1)*ts - 4*m1**2*ms**2*s**(-1)*ts - 4*m1**2*s**(-1)*u*ts
     +     - 4*m2**2*ms**2*s**(-1)*ts - 4*m2**2*s**(-1)*u*ts + 8*ms**2*
     +    s**(-1)*u*ts + 4*ms**2*ts + 4*u*ts )
      QBV = QBV + SCC(5,1)*kaellen**(-1)*b_tx * (  - 4*m1**2*m2**2
     +    *ms**2*t1**(-1)*t2**(-1)*ts - 3*m1**2*m2**2*s*t1**(-1)*
     +    t2**(-1)*ts - 4*m1**2*m2**2*t1**(-1)*ts + 4*m1**2*m2**2*
     +    t2**(-1)*ts - 11*m1**2*m2**4*t1**(-1)*t2**(-1)*ts - m1**2*
     +    ms**2*s*t1**(-1)*t2**(-1)*ts + 2*m1**2*ms**2*t1**(-1)*ts + 2*
     +    m1**2*ms**2*t2**(-1)*ts + m1**2*s*t1**(-1)*ts + 3*m1**2*s*
     +    t2**(-1)*ts - 2*m1**2*u*t1**(-1)*ts - 2*m1**2*ts + 13*m1**4*
     +    m2**2*t1**(-1)*t2**(-1)*ts + 2*m1**4*ms**2*t1**(-1)*t2**(-1)*
     +    ts + 3*m1**4*s*t1**(-1)*t2**(-1)*ts + 2*m1**4*t1**(-1)*ts - 5
     +    *m1**4*t2**(-1)*ts - 5*m1**6*t1**(-1)*t2**(-1)*ts - 3*m2**2*
     +    ms**2*s*t1**(-1)*t2**(-1)*ts + m2**2*ms**2*t1**(-1)*ts + 2*
     +    m2**2*ms**2*t2**(-1)*ts - 3*m2**2*s*t1**(-1)*ts + 4*m2**2*s*
     +    t2**(-1)*ts - m2**2*u*t1**(-1)*ts - 2*m2**2*ts + 2*m2**4*
     +    ms**2*t1**(-1)*t2**(-1)*ts - 4*m2**4*s*t1**(-1)*t2**(-1)*ts
     +     + 2*m2**4*t1**(-1)*ts - 3*m2**4*t2**(-1)*ts + 3*m2**6*
     +    t1**(-1)*t2**(-1)*ts )
      QBV = QBV + SCC(5,1)*kaellen**(-1)*b_tx * (  - 2*ms**2*s*
     +    t1**(-1)*ts - ms**2*s*t2**(-1)*ts + ms**2*u*t1**(-1)*ts + 
     +    ms**2*ts + 2*s*u*t1**(-1)*ts + 3*s*ts - u*ts - u**2*t1**(-1)*
     +    ts )
      QBV = QBV + SCC(5,1)*kaellen**(-1)*b_um * ( 4*m1**2*ms**2*
     +    s**(-1)*us + 4*m1**2*s**(-1)*u*us + 8*m1**2*us - 4*m1**4*
     +    s**(-1)*us + 4*m2**2*ms**2*s**(-1)*us + 4*m2**2*s**(-1)*u*us
     +     + 8*m2**2*us - 4*m2**4*s**(-1)*us - 8*ms**2*s**(-1)*u*us - 4
     +    *ms**2*us - 4*s*us - 4*u*us )
      QBV = QBV + SCC(5,1)*kaellen**(-1)*b_ux * (  - 4*m1**2*m2**2
     +    *ms**2*u1**(-1)*u2**(-1)*us - 9*m1**2*m2**2*s*u1**(-1)*
     +    u2**(-1)*us - 3*m1**2*m2**2*u1**(-1)*us - 7*m1**2*m2**2*
     +    u2**(-1)*us + 7*m1**2*m2**4*u1**(-1)*u2**(-1)*us - m1**2*
     +    ms**2*s*u1**(-1)*u2**(-1)*us + 3*m1**2*ms**2*u1**(-1)*us + 2*
     +    m1**2*ms**2*u2**(-1)*us - 3*m1**2*s*u1**(-1)*us + 2*m1**2*s*
     +    u2**(-1)*us + m1**2*u*u2**(-1)*us - m1**2*us - 5*m1**4*m2**2*
     +    u1**(-1)*u2**(-1)*us + 2*m1**4*ms**2*u1**(-1)*u2**(-1)*us + 
     +    m1**6*u1**(-1)*u2**(-1)*us - 3*m2**2*ms**2*s*u1**(-1)*
     +    u2**(-1)*us + 2*m2**2*ms**2*u1**(-1)*us + m2**2*ms**2*
     +    u2**(-1)*us + 3*m2**2*s*u1**(-1)*us - 3*m2**2*s*u2**(-1)*us
     +     + m2**2*u*u1**(-1)*us - 2*m2**2*us + 2*m2**4*ms**2*u1**(-1)*
     +    u2**(-1)*us + 5*m2**4*s*u1**(-1)*u2**(-1)*us - 2*m2**4*
     +    u1**(-1)*us + m2**4*u2**(-1)*us - 3*m2**6*u1**(-1)*u2**(-1)*
     +    us - 3*ms**2*s*u1**(-1)*us - ms**2*s*u2**(-1)*us - ms**2*u*
     +    u1**(-1)*us )
      QBV = QBV + SCC(5,1)*kaellen**(-1)*b_ux * ( ms**2*u*u2**(-1)
     +    *us + s*u*u1**(-1)*us - s*us + s**2*u1**(-1)*us - s**2*
     +    u2**(-1)*us - u*us + u**2*u2**(-1)*us )
      QBV = QBV + SCC(5,1)*b_tm * (  - 4*s**(-1)*ts )
      QBV = QBV + SCC(5,1)*b_tx * (  - m1**2*s**(-1)*t2**(-1)*ts
     +     + 2*m1**2*t1**(-1)*t2**(-1)*ts - 4*m2**2*t1**(-1)*t2**(-1)*
     +    ts + 2*ms**2*s**(-1)*t1**(-1)*ts + 2*ms**2*s**(-1)*t2**(-1)*
     +    ts + 2*ms**2*t1**(-1)*t2**(-1)*ts - 2*s**(-1)*u*t1**(-1)*ts
     +     - s**(-1)*u*t2**(-1)*ts - 3*s**(-1)*ts - t1**(-1)*ts + 3*
     +    t2**(-1)*ts )
      QBV = QBV + SCC(5,1)*b_um * (  - 4*s**(-1)*us )
      QBV = QBV + SCC(5,1)*b_ux * (  - 2*m1**2*s**(-1)*u2**(-1)*us
     +     - 2*m2**2*s**(-1)*u1**(-1)*us - 3*m2**2*s**(-1)*u2**(-1)*us
     +     - 2*m2**2*u1**(-1)*u2**(-1)*us + 2*ms**2*s**(-1)*u1**(-1)*us
     +     + 2*ms**2*s**(-1)*u2**(-1)*us + 2*ms**2*u1**(-1)*u2**(-1)*us
     +     + 3*s**(-1)*u*u2**(-1)*us - 3*s**(-1)*us - u1**(-1)*us + 3*
     +    u2**(-1)*us )
      QBV = QBV + SCC(6,1)*kaellen**(-1)*b_um2 * ( 8*m1**2*m2**2*
     +    s**(-1)*us + 4*m1**2*ms**2*s**(-1)*us - 4*m1**2*s**(-1)*u*us
     +     + 4*m2**2*ms**2*s**(-1)*us - 4*m2**2*s**(-1)*u*us - 8*ms**2*
     +    s**(-1)*u*us - 4*ms**2*us + 4*u*us )
      QBV = QBV + SCC(6,1)*kaellen**(-1)*b_ux2 * (  - 4*m1**2*
     +    m2**2*ms**2*u1**(-1)*u2**(-1)*us + 6*m1**2*m2**2*s*u1**(-1)*
     +    u2**(-1)*us + 6*m1**2*m2**2*u1**(-1)*us - 2*m1**2*m2**2*
     +    u2**(-1)*us + 10*m1**2*m2**4*u1**(-1)*u2**(-1)*us - 2*m1**2*
     +    ms**2*s*u1**(-1)*u2**(-1)*us + 2*m1**2*ms**2*u1**(-1)*us + 2*
     +    m1**2*ms**2*u2**(-1)*us - 6*m1**2*s*u1**(-1)*us - 8*m1**2*s*
     +    u2**(-1)*us - 2*m1**2*us - 8*m1**4*m2**2*u1**(-1)*u2**(-1)*us
     +     + 2*m1**4*ms**2*u1**(-1)*u2**(-1)*us - 8*m1**4*s*u1**(-1)*
     +    u2**(-1)*us + 2*m1**4*u2**(-1)*us + 2*m1**6*u1**(-1)*u2**(-1)
     +    *us - 2*m2**2*ms**2*s*u1**(-1)*u2**(-1)*us + 2*m2**2*ms**2*
     +    u1**(-1)*us + 2*m2**2*ms**2*u2**(-1)*us - 2*m2**2*s*u1**(-1)*
     +    us - 6*m2**2*s*u2**(-1)*us - 4*m2**2*u*u1**(-1)*us + 4*m2**2*
     +    us + 2*m2**4*ms**2*u1**(-1)*u2**(-1)*us + 6*m2**4*s*u1**(-1)*
     +    u2**(-1)*us + 4*m2**4*u2**(-1)*us - 4*m2**6*u1**(-1)*u2**(-1)
     +    *us - 2*ms**2*s*u1**(-1)*us - 2*ms**2*s*u2**(-1)*us + 6*s*u*
     +    u1**(-1)*us )
      QBV = QBV + SCC(6,1)*kaellen**(-1)*b_ux2 * (  - 6*s*us + 4*
     +    s**2*u1**(-1)*us - 2*u*us + 2*u**2*u1**(-1)*us )
      QBV = QBV + SCC(6,1)*kaellen**(-1)*b_tm4 * (  - 4*m1**2*
     +    ms**2*s**(-1)*ts + 4*m1**2*s**(-1)*u*ts + 8*m1**2*ts - 4*
     +    m1**4*s**(-1)*ts - 4*m2**2*ms**2*s**(-1)*ts + 4*m2**2*s**(-1)
     +    *u*ts + 8*m2**2*ts - 4*m2**4*s**(-1)*ts + 8*ms**2*s**(-1)*u*
     +    ts + 4*ms**2*ts - 4*s*ts - 4*u*ts )
      QBV = QBV + SCC(6,1)*kaellen**(-1)*b_tx4 * (  - 4*m1**2*
     +    m2**2*ms**2*t1**(-1)*t2**(-1)*ts + 6*m1**2*m2**2*s*t1**(-1)*
     +    t2**(-1)*ts - 2*m1**2*m2**2*t1**(-1)*ts + 6*m1**2*m2**2*
     +    t2**(-1)*ts - 8*m1**2*m2**4*t1**(-1)*t2**(-1)*ts - 2*m1**2*
     +    ms**2*s*t1**(-1)*t2**(-1)*ts + 2*m1**2*ms**2*t1**(-1)*ts + 2*
     +    m1**2*ms**2*t2**(-1)*ts - 6*m1**2*s*t1**(-1)*ts + 4*m1**2*s*
     +    t2**(-1)*ts + 2*m1**2*ts + 10*m1**4*m2**2*t1**(-1)*t2**(-1)*
     +    ts + 2*m1**4*ms**2*t1**(-1)*t2**(-1)*ts + 6*m1**4*s*t1**(-1)*
     +    t2**(-1)*ts + 4*m1**4*t1**(-1)*ts - 2*m1**4*t2**(-1)*ts - 4*
     +    m1**6*t1**(-1)*t2**(-1)*ts - 2*m2**2*ms**2*s*t1**(-1)*
     +    t2**(-1)*ts + 2*m2**2*ms**2*t1**(-1)*ts + 2*m2**2*ms**2*
     +    t2**(-1)*ts - 8*m2**2*s*t1**(-1)*ts - 4*m2**2*s*t2**(-1)*ts
     +     - 4*m2**2*u*t2**(-1)*ts - 4*m2**2*ts + 2*m2**4*ms**2*
     +    t1**(-1)*t2**(-1)*ts - 8*m2**4*s*t1**(-1)*t2**(-1)*ts + 2*
     +    m2**4*t1**(-1)*ts + 2*m2**4*t2**(-1)*ts + 2*m2**6*t1**(-1)*
     +    t2**(-1)*ts )
      QBV = QBV + SCC(6,1)*kaellen**(-1)*b_tx4 * (  - 2*ms**2*s*
     +    t1**(-1)*ts - 2*ms**2*s*t2**(-1)*ts - 2*s*u*t2**(-1)*ts - 4*s
     +    *ts + 2*u*ts + 2*u**2*t2**(-1)*ts )
      QBV = QBV + SCC(6,1)*b_um2 * ( 2*m1**2*s**(-2)*us + 2*m2**2*
     +    s**(-2)*us - 4*mg**2*s**(-2)*us + 4*ms**2*s**(-2)*us - 4*
     +    s**(-2)*u*us - 2*s**(-1)*us )
      QBV = QBV + SCC(6,1)*b_ux2 * ( 2*m1**2*mg**2*s**(-1)*
     +    u1**(-1)*u2**(-1)*us - 2*m1**2*ms**2*s**(-1)*u1**(-1)*
     +    u2**(-1)*us + 2*m2**2*mg**2*s**(-1)*u1**(-1)*u2**(-1)*us - 2*
     +    m2**2*ms**2*s**(-1)*u1**(-1)*u2**(-1)*us + 4*m2**2*u1**(-1)*
     +    u2**(-1)*us - 2*mg**2*u1**(-1)*u2**(-1)*us + 2*ms**2*s**(-1)*
     +    u1**(-1)*us + 2*ms**2*s**(-1)*u2**(-1)*us + 4*ms**2*u1**(-1)*
     +    u2**(-1)*us - 2*s*u1**(-1)*u2**(-1)*us - 2*u1**(-1)*us - 2*
     +    u2**(-1)*us )
      QBV = QBV + SCC(6,1)*b_tm4 * (  - 2*m1**2*s**(-2)*ts - 2*
     +    m2**2*s**(-2)*ts - 4*mg**2*s**(-2)*ts + 4*ms**2*s**(-2)*ts + 
     +    4*s**(-2)*u*ts + 2*s**(-1)*ts )
      QBV = QBV + SCC(6,1)*b_tx4 * ( 2*m1**2*mg**2*s**(-1)*
     +    t1**(-1)*t2**(-1)*ts - 2*m1**2*ms**2*s**(-1)*t1**(-1)*
     +    t2**(-1)*ts + 4*m1**2*t1**(-1)*t2**(-1)*ts + 2*m2**2*mg**2*
     +    s**(-1)*t1**(-1)*t2**(-1)*ts - 2*m2**2*ms**2*s**(-1)*t1**(-1)
     +    *t2**(-1)*ts - 2*mg**2*t1**(-1)*t2**(-1)*ts + 2*ms**2*s**(-1)
     +    *t1**(-1)*ts + 2*ms**2*s**(-1)*t2**(-1)*ts + 4*ms**2*t1**(-1)
     +    *t2**(-1)*ts - 2*s*t1**(-1)*t2**(-1)*ts - 2*t1**(-1)*ts - 2*
     +    t2**(-1)*ts )
      QBV = QBV + SCD(1,1)*b_tm * (  - 4 )
      QBV = QBV + SCD(1,1)*b_tx * (  - 2 + 4*m1**2*ms**2*t1**(-1)*
     +    t2**(-1) + 2*m1**2*t1**(-1) - 2*m1**2*t2**(-1) - 2*m1**4*
     +    t1**(-1)*t2**(-1) - 2*ms**2*t1**(-1) + 2*ms**2*t2**(-1) - 2*
     +    ms**4*t1**(-1)*t2**(-1) )
      QBV = QBV + SCD(1,2)*b_um * (  - 4 )
      QBV = QBV + SCD(1,2)*b_ux * (  - 3 + m1**2*u2**(-1) + 4*
     +    m2**2*ms**2*u1**(-1)*u2**(-1) - 2*m2**2*u1**(-1) + 2*m2**2*
     +    u2**(-1) - 2*m2**4*u1**(-1)*u2**(-1) + 2*ms**2*u1**(-1) - 2*
     +    ms**2*u2**(-1) - 2*ms**4*u1**(-1)*u2**(-1) - s*u2**(-1) - t*
     +    u2**(-1) )
      QBV = QBV + SCD(2,1)*b_um2 * (  - 2*m1**2*mg**2*s**(-2)*
     +    tg**(-1)*us + 2*m1**2*ms**2*s**(-2)*tg**(-1)*us - 2*m2**2*
     +    mg**2*s**(-2)*tg**(-1)*us + 2*m2**2*ms**2*s**(-2)*tg**(-1)*us
     +     + 4*mg**2*ms**2*s**(-2)*tg**(-1)*us + 4*mg**2*s**(-2)*us - 4
     +    *ms**2*s**(-2)*us - 4*ms**4*s**(-2)*tg**(-1)*us + 2*s**(-1)*
     +    us )
      QBV = QBV + SCD(2,1)*b_ux2 * (  - 4*m1**2*mg**2*ms**2*
     +    s**(-1)*u1**(-1)*tg**(-1)*u2**(-1)*us + 4*m1**2*mg**2*
     +    u1**(-1)*tg**(-1)*u2**(-1)*us + 2*m1**2*mg**4*s**(-1)*
     +    u1**(-1)*tg**(-1)*u2**(-1)*us - 4*m1**2*ms**2*u1**(-1)*
     +    tg**(-1)*u2**(-1)*us + 2*m1**2*ms**4*s**(-1)*u1**(-1)*
     +    tg**(-1)*u2**(-1)*us + 2*m1**2*s*u1**(-1)*tg**(-1)*u2**(-1)*
     +    us - 4*m2**2*mg**2*ms**2*s**(-1)*u1**(-1)*tg**(-1)*u2**(-1)*
     +    us + 2*m2**2*mg**2*u1**(-1)*tg**(-1)*u2**(-1)*us + 2*m2**2*
     +    mg**4*s**(-1)*u1**(-1)*tg**(-1)*u2**(-1)*us + 2*m2**2*ms**4*
     +    s**(-1)*u1**(-1)*tg**(-1)*u2**(-1)*us + 2*mg**2*ms**2*s**(-1)
     +    *u1**(-1)*tg**(-1)*us + 2*mg**2*ms**2*s**(-1)*tg**(-1)*
     +    u2**(-1)*us + 8*mg**2*ms**2*u1**(-1)*tg**(-1)*u2**(-1)*us - 4
     +    *mg**2*s*u1**(-1)*tg**(-1)*u2**(-1)*us - 2*mg**2*u1**(-1)*
     +    tg**(-1)*us - 2*mg**4*u1**(-1)*tg**(-1)*u2**(-1)*us + 8*ms**2
     +    *s*u1**(-1)*tg**(-1)*u2**(-1)*us + 6*ms**2*u1**(-1)*tg**(-1)*
     +    us )
      QBV = QBV + SCD(2,1)*b_ux2 * ( 2*ms**2*tg**(-1)*u2**(-1)*us
     +     - 2*ms**4*s**(-1)*u1**(-1)*tg**(-1)*us - 2*ms**4*s**(-1)*
     +    tg**(-1)*u2**(-1)*us - 8*ms**4*u1**(-1)*tg**(-1)*u2**(-1)*us
     +     - 2*s*u1**(-1)*tg**(-1)*us - 2*s**2*u1**(-1)*tg**(-1)*
     +    u2**(-1)*us )
      QBV = QBV + SCD(2,2)*b_tm4 * (  - 2*m1**2*mg**2*s**(-2)*
     +    ug**(-1)*ts + 2*m1**2*ms**2*s**(-2)*ug**(-1)*ts - 2*m2**2*
     +    mg**2*s**(-2)*ug**(-1)*ts + 2*m2**2*ms**2*s**(-2)*ug**(-1)*ts
     +     + 4*mg**2*ms**2*s**(-2)*ug**(-1)*ts + 4*mg**2*s**(-2)*ts - 4
     +    *ms**2*s**(-2)*ts - 4*ms**4*s**(-2)*ug**(-1)*ts + 2*s**(-1)*
     +    ts )
      QBV = QBV + SCD(2,2)*b_tx4 * (  - 4*m1**2*mg**2*ms**2*
     +    s**(-1)*t1**(-1)*ug**(-1)*t2**(-1)*ts + 2*m1**2*mg**2*
     +    t1**(-1)*ug**(-1)*t2**(-1)*ts + 2*m1**2*mg**4*s**(-1)*
     +    t1**(-1)*ug**(-1)*t2**(-1)*ts + 2*m1**2*ms**4*s**(-1)*
     +    t1**(-1)*ug**(-1)*t2**(-1)*ts - 4*m2**2*mg**2*ms**2*s**(-1)*
     +    t1**(-1)*ug**(-1)*t2**(-1)*ts + 4*m2**2*mg**2*t1**(-1)*
     +    ug**(-1)*t2**(-1)*ts + 2*m2**2*mg**4*s**(-1)*t1**(-1)*
     +    ug**(-1)*t2**(-1)*ts - 4*m2**2*ms**2*t1**(-1)*ug**(-1)*
     +    t2**(-1)*ts + 2*m2**2*ms**4*s**(-1)*t1**(-1)*ug**(-1)*
     +    t2**(-1)*ts + 2*m2**2*s*t1**(-1)*ug**(-1)*t2**(-1)*ts + 2*
     +    mg**2*ms**2*s**(-1)*t1**(-1)*ug**(-1)*ts + 2*mg**2*ms**2*
     +    s**(-1)*ug**(-1)*t2**(-1)*ts + 8*mg**2*ms**2*t1**(-1)*
     +    ug**(-1)*t2**(-1)*ts - 4*mg**2*s*t1**(-1)*ug**(-1)*t2**(-1)*
     +    ts - 2*mg**2*ug**(-1)*t2**(-1)*ts - 2*mg**4*t1**(-1)*ug**(-1)
     +    *t2**(-1)*ts + 8*ms**2*s*t1**(-1)*ug**(-1)*t2**(-1)*ts + 2*
     +    ms**2*t1**(-1)*ug**(-1)*ts )
      QBV = QBV + SCD(2,2)*b_tx4 * ( 6*ms**2*ug**(-1)*t2**(-1)*ts
     +     - 2*ms**4*s**(-1)*t1**(-1)*ug**(-1)*ts - 2*ms**4*s**(-1)*
     +    ug**(-1)*t2**(-1)*ts - 8*ms**4*t1**(-1)*ug**(-1)*t2**(-1)*ts
     +     - 2*s*ug**(-1)*t2**(-1)*ts - 2*s**2*t1**(-1)*ug**(-1)*
     +    t2**(-1)*ts )

*/










  output.finite =-8.+sqr(Constants::pi);







  
  return output;
}

// ofstream myfile ("Our_MSSM_scale.txt");
// if (myfile.is_open())
//   {
//     myfile << "Our scale " << scale() << endl;
//     myfile.close();
//   }
// else cout << "Unable to open file";

double MEPP2NeutralinoNeutralinoPowheg::
qqbarME(vector<SpinorWaveFunction>    & sp ,
	vector<SpinorBarWaveFunction> & sbar ,
	vector<SpinorWaveFunction>    & spout ,
	vector<SpinorBarWaveFunction> & sbarout ,
	//ScalarWaveFunction & interqL,ScalarWaveFunction & interqR,
	bool first) const {
  // scale for the process
  const Energy2 q2(scale());
  // squarks for the t-channel
  tcPDPtr squark[2]= {getParticleData(1000000+abs(mePartonData()[0]->id())),
		      getParticleData(2000000+abs(mePartonData()[0]->id()))};
  // conjugate spinors for t-channel exchange diagram
  vector<SpinorWaveFunction> sbaroutconj;
  vector<SpinorBarWaveFunction> spoutconj;
  for(unsigned int of1=0;of1<2;++of1) {
    sbaroutconj.push_back(SpinorWaveFunction (-sbarout[of1].momentum(),
					      sbarout[of1].particle(),
					      sbarout[of1].wave().bar().conjugate(),
				 	      sbarout[of1].direction()));
    spoutconj.push_back(SpinorBarWaveFunction(-spout[of1].momentum(),
					      spout[of1].particle(),
					      spout[of1].wave().bar().conjugate(),
					      spout[of1].direction()));
  }
  // storage of the matrix elements for specific diagrams
  vector<double> me(5, 0.);
  double me2(0.);
  // storage of the individual diagrams
  vector<Complex> diag(5, 0.);
  ProductionMatrixElement pme(PDT::Spin1Half, PDT::Spin1Half, 
			      PDT::Spin1Half, PDT::Spin1Half);
  unsigned int propopt;
  if(status()==RealQG || status()==RealQBarG){
    propopt = 7;
  }
  else propopt = 3;
  // loop over the helicities and calculate the matrix elements
  for(unsigned int if1 = 0; if1 < 2; ++if1) {
    for(unsigned int if2 = 0; if2 < 2; ++if2) {
      VectorWaveFunction interV = FFZVertex_->
	evaluate(q2, 1, Z0_, sp[if1],sbar[if2]);
      for(unsigned int of1 = 0; of1 < 2; ++of1) {
	for(unsigned int of2 = 0; of2 < 2; ++of2) {
	  // s-channel Z exchange
	  diag[0] = NNZVertex_->evaluate(q2, spout[of1],  sbarout[of2], interV);
	  // t-channel squark exchanges	  
	  for(unsigned int iq=0;iq<2;++iq) {
	    // 1st t-channel
	    ScalarWaveFunction intersq = NFSVertex_->
		evaluate(q2, propopt, squark[iq], sp[if1], sbarout[of2]);
	    diag[2*iq+1] = 
	      NFSVertex_->evaluate(q2, spout[of1], sbar[if2], intersq);
	    // swapped t-channel
	    intersq = NFSVertex_->
	      evaluate(q2, propopt, squark[iq], sp[if1], spoutconj[of1]);
	    diag[2*iq+2] = 
	      -NFSVertex_->evaluate(q2, sbaroutconj[of2], sbar[if2], intersq);
	  }
	  // individual diagrams
	  for(unsigned int id=0;id<5;++id){
	    me[id] += norm(diag[id]);
	  }
	  // sum up the matrix elements
	  Complex total = std::accumulate(diag.begin(),diag.end(),Complex(0.));
	  me2 += norm(total);
	  pme(if1, if2, of1, of2) = total;
	}
      }
    }
  }
  if(first) {
    DVector save(5);
    for(DVector::size_type ix = 0; ix < 5; ++ix)
      save[ix] = me[ix]/12.;
    meInfo(save);
    me_.reset(pme);
  }
  if(mePartonData()[2]->id() == mePartonData()[3]->id()) me2 *= 0.5;
  return me2/12.;
}

double MEPP2NeutralinoNeutralinoPowheg::loME(const cPDVector & particles,
				    const vector<Lorentz5Momentum> & momenta,
				    bool first) const {
  // wavefunctions for the incoming fermions
  vector<SpinorWaveFunction> sp(2);
  vector<SpinorBarWaveFunction> sbar(2);
  for( unsigned int i = 0; i < 2; ++i ) {
    sp[i]   = SpinorWaveFunction   (momenta[0], particles[0], i,
				    incoming);
    sbar[i] = SpinorBarWaveFunction(momenta[1], particles[1], i,
				    incoming);
  }
  // outgoing neutralino wavefunctions
  vector<SpinorWaveFunction> spout(2);
  vector<SpinorBarWaveFunction> sbarout(2);
  for( unsigned int i = 0; i < 2; ++i ) {
    spout[i] = SpinorWaveFunction(momenta[2], particles[2], i,
				  outgoing);
    sbarout[i] = SpinorBarWaveFunction(momenta[3], particles[3], i,
				       outgoing);
  }
  return qqbarME(sp,sbar,spout,sbarout,first);
}


double MEPP2NeutralinoNeutralinoPowheg::
realME(const cPDVector & particles,
       const vector<Lorentz5Momentum> & momenta) const {
  vector<SpinorWaveFunction> sp(2);
  vector<SpinorBarWaveFunction> sbar(2);
  vector<VectorWaveFunction> gluon(2);
  // squarks for the t-channel
  tcPDPtr squark[2]= {getParticleData(1000000+abs(mePartonData()[0]->id())),
		      getParticleData(2000000+abs(mePartonData()[0]->id()))};
  // wavefunctions for the q qbar -> chi chi g process
  if(particles[0]->id()==-particles[1]->id()) {
    for( unsigned int i = 0; i < 2; ++i ) {
      sp[i]   = SpinorWaveFunction   (momenta[0],particles[0],  i,incoming);
      sbar[i] = SpinorBarWaveFunction(momenta[1],particles[1],  i,incoming);
      gluon[i]= VectorWaveFunction   (momenta[4],particles[4],2*i,outgoing);
    }
  }
  else if(particles[0]->id()==ParticleID::g &&
	  particles[1]->id()<0) {
    for( unsigned int i = 0; i < 2; ++i ) {
      sp[i]   = SpinorWaveFunction   (momenta[4],particles[4],  i,outgoing);
      sbar[i] = SpinorBarWaveFunction(momenta[1],particles[1],  i,incoming);
      gluon[i]= VectorWaveFunction   (momenta[0],particles[0],2*i,incoming);
    }
  }
  else if(particles[0]->id()>0 &&
	  particles[1]->id()==ParticleID::g) {
    for( unsigned int i = 0; i < 2; ++i ) {
      sp[i]   = SpinorWaveFunction   (momenta[0],particles[0],  i,incoming);
      sbar[i] = SpinorBarWaveFunction(momenta[4],particles[4],  i,outgoing);
      gluon[i]= VectorWaveFunction   (momenta[1],particles[1],2*i,incoming);
    }
  }
  else {
    for(unsigned int ix=0;ix<particles.size();++ix) {
      cerr << particles[ix]->PDGName() << " " << momenta[ix]/GeV << "\n";
    }
    assert(false);
  }
  // wavefunctions for the outgoing neutralinos
  vector<SpinorWaveFunction> spout(2),sbaroutconj(2);
  vector<SpinorBarWaveFunction> sbarout(2),spoutconj(2);
  ScalarWaveFunction intersq,intersq2;	  
  for( unsigned int i = 0; i < 2; ++i ) {
    spout[i]       =    SpinorWaveFunction(momenta[2], particles[2], i,
					   outgoing);
    spoutconj[i]   = SpinorBarWaveFunction(momenta[2], particles[2],
					   spout[i].wave().bar().conjugate(),
					   outgoing);
    sbarout[i]     = SpinorBarWaveFunction(momenta[3], particles[3], i,
					   outgoing);
    sbaroutconj[i] = SpinorWaveFunction   (momenta[3], particles[3],
					   sbarout[i].wave().bar().conjugate(),
					   outgoing);
  }
  double output(0.);
  vector<Complex> diag(14,0.);
  Energy2 q2 = scale();
  unsigned int propopt;
  if(status()==RealQG || status()==RealQBarG){
    propopt = 7;
  }
  else propopt = 3; 
  for(unsigned int ihel1=0;ihel1<2;++ihel1) {
    for(unsigned int ihel2=0;ihel2<2;++ihel2) {
      for(unsigned int ohel1=0;ohel1<2;++ohel1) {
	SpinorWaveFunction    inters = FFGVertex_->evaluate(q2, 5, sp[ihel1].particle()->CC(),
							    sp[ihel1], gluon[ohel1]);
	SpinorBarWaveFunction interb = FFGVertex_->evaluate(q2,5,sbar[ihel2].particle()->CC(),
							    sbar[ihel2],gluon[ohel1]);
	VectorWaveFunction interV[2] = 
	  {FFZVertex_->evaluate(q2, 1, Z0_, inters,sbar[ihel2]),
	   FFZVertex_->evaluate(q2, 1, Z0_, sp[ihel1],interb)};
	
	for(unsigned int ohel2=0;ohel2<2;++ohel2) {
	  for(unsigned int ohel3=0;ohel3<2;++ohel3) {
	    // s-channel Z exchange diagrams
	    // first Z diagram (emission from quark)
	    diag[0] = NNZVertex_->evaluate(q2, spout[ohel2],  sbarout[ohel3], interV[0]);
	    // second Z diagram (emission from anti-quark)
	    diag[1] = NNZVertex_->evaluate(q2, spout[ohel2],  sbarout[ohel3], interV[1]);
	    // t-channel squark exchanges
	    for(unsigned int iq=0;iq<2;++iq) {
	      //2.0*ms/(abs(m1)+abs(m2))>20.0
	      if(   (2.*squark[iq]->mass())/(momenta[2].mass()+momenta[3].mass()) > 20.){
		diag[6*iq+2]=0.;
		diag[6*iq+3]=0.;
		diag[6*iq+4]=0.;
		diag[6*iq+5]=0.;
		diag[6*iq+6]=0.;
		diag[6*iq+7]=0.;
	      } else{
		
		// 1st t-channel
		// emission quark
		intersq = NFSVertex_->
		  evaluate(q2, propopt, squark[iq], inters, sbarout[ohel3]);
		diag[6*iq+2] = 
		  NFSVertex_->evaluate(q2, spout[ohel2], sbar[ihel2], intersq);
		// emission antiquark
		intersq = NFSVertex_->
		  evaluate(q2, propopt, squark[iq], sp[ihel1], sbarout[ohel3]);
		diag[6*iq+3] = 
		  NFSVertex_->evaluate(q2, spout[ohel2], interb, intersq);
		// emission from intermediate
		intersq2 = GSSVertex_->evaluate(q2,propopt,squark[iq],gluon[ohel1],intersq);
		diag[6*iq+4] = 
		  NFSVertex_->evaluate(q2, spout[ohel2], sbar[ihel2], intersq2);
		// swapped t-channel
		// emission quark
		intersq = NFSVertex_->
		  evaluate(q2, propopt, squark[iq], inters, spoutconj[ohel2]);
		diag[6*iq+5] = 
		  -NFSVertex_->evaluate(q2, sbaroutconj[ohel3], sbar[ihel2], intersq);
		// emission antiquark
		intersq = NFSVertex_->
		  evaluate(q2, propopt, squark[iq], sp[ihel1], spoutconj[ohel2]);
		diag[6*iq+6] = 
		  -NFSVertex_->evaluate(q2, sbaroutconj[ohel3], interb, intersq);
		// emission from intermediate
		intersq2 = GSSVertex_->evaluate(q2,propopt,squark[iq],gluon[ohel1],intersq);
		diag[6*iq+7] = 
		  -NFSVertex_->evaluate(q2, sbaroutconj[ohel3], sbar[ihel2], intersq2);}
	    }
	    // add them up
	    Complex total = std::accumulate(diag.begin(),diag.end(),Complex(0.));
	    output += norm(total);
	  }
	}
      }
    }
  }
  // strong coupling
  double gs2 = norm(FFGVertex_->norm());
  
  
  double totcount = 0.;
  //   bool first = 0.;
  //   bool second = 0.;
  Energy sqmom2,sqmass;


  // subtract the on-shell squark decay if neccessary
  if(abs(particles[4]->id())<=6) {
    Energy roots = (momenta[0]+momenta[1]).m();
    // off-shell masses of the squarks
    Energy2 msq1 = (momenta[2]+momenta[4]).m2();
    Energy2 msq2 = (momenta[3]+momenta[4]).m2();
    Lorentz5Momentum pgluon = 
      particles[0]->id() == ParticleID::g ? momenta[0] : momenta[1];
    tcPDPtr quark = particles[0]->id() == ParticleID::g ? particles[1] : particles[0];
    FFSVertexPtr vertex = dynamic_ptr_cast<FFSVertexPtr>(NFSVertex_);
    double Cf = 4./3.;
    double Nc = 3.;
    Energy2 sh = (momenta[0]+momenta[1]).m2();
    for(unsigned int ix=0;ix<2;++ix) {
      Energy2 sqwidth2 = sqr(squark[ix]->width());
      // is on-shell mass allowed for first neutralino and squark
      //      if( roots >= squark[ix]->mass() + momenta[2].mass() ) {
      if(  (roots >= squark[ix]->mass() + momenta[2].mass() )
	   && (squark[ix]->mass() > momenta[3].mass())){
	// Create a counterterm for the squark decay pole.
	Energy2 mneut2 = sqr(momenta[2].mass()); 
	Energy2 t3 = (pgluon-momenta[3]-momenta[4]).m2()-msq2;
	Energy2 u4 = (pgluon-momenta[2]).m2()-mneut2;
	vertex->setCoupling(q2,quark,particles[2],squark[ix]);
 	double a2 = norm(vertex->norm())*
 	  (norm(vertex->left())+norm(vertex->right()));
	double sqprod2 = (1./8.) * gs2 * Cf * Nc * a2 *
	  (-u4/sh - (2*(msq2-mneut2)*u4)/sh/t3 * (1+mneut2/u4+msq2/t3));
	vertex->setCoupling(q2,quark,particles[3],squark[ix]);
 	a2 = norm(vertex->norm())*
 	  (norm(vertex->left())+norm(vertex->right()));
	Energy2 sqdecay2 = 4. * a2 * (msq2-sqr(particles[3]->mass()));
	Energy4 denom = sqr(msq2-sqr(squark[ix]->mass())) + 
	  sqr(squark[ix]->mass())*sqwidth2;
	double sqcounter = sqprod2 * sqdecay2 * UnitRemoval::E2 / denom;

	if( ((2.*squark[ix]->mass())/(momenta[2].mass()+momenta[3].mass()) > 100.)
	    && (squark[ix]->mass()>1.e4*GeV))
	  sqcounter = 0;

	totcount += sqcounter;
      }

      // is on-shell mass allowed for second neutralino and squark
      //else if( roots >= squark[ix]->mass() + momenta[3].mass() ) {
      if((roots >= squark[ix]->mass() + momenta[3].mass() )
	 && (squark[ix]->mass() > momenta[2].mass())){
	Energy2 mneut2 = sqr(momenta[3].mass());
	Energy2 t3 = (pgluon-momenta[2]-momenta[4]).m2()-msq1;
	Energy2 u4 = (pgluon-momenta[3]).m2()-mneut2;
	vertex->setCoupling(q2,quark,particles[3],squark[ix]);
 	double a2 = norm(vertex->norm())*
 	  (norm(vertex->left())+norm(vertex->right()));
	double sqprod2 = (1./8.) * gs2 * Cf * Nc * a2 *
	  (-u4/sh - (2*(msq1-mneut2)*u4)/sh/t3 * (1+mneut2/u4+msq1/t3));
	vertex->setCoupling(q2,quark,particles[2],squark[ix]);
	a2 = norm(vertex->norm())*
	  (norm(vertex->left())+norm(vertex->right()));
	Energy2 sqdecay2 = 4. * a2 * (msq1-sqr(particles[2]->mass()));
	Energy4 denom = sqr(msq1-sqr(squark[ix]->mass())) + 
	  sqr(squark[ix]->mass())*sqwidth2;
	double sqcounter = sqprod2 * sqdecay2 * UnitRemoval::E2 / denom;

	if( ((2.*squark[ix]->mass())/(momenta[2].mass()+momenta[3].mass()) > 100.)
	    && (squark[ix]->mass()>1.e4*GeV))
	  sqcounter = 0;

	totcount += sqcounter;
      }
    }
  }


//   if( first == 1. || second == 1.){
//     cout << "\n \n" << particles[0]->id() << "\t" << particles[1]->id() << endl;
//     cout << particles[2]->id() << "\t" << particles[3]->id() << "\t" << particles[4]->id() << endl;
//     cout << sqmom2/GeV << "\t" << sqmass/GeV << endl;
//     cout << "ratios: " << output << "\t" << totcount << "\t" << output/totcount << "\t" << (output-totcount)/totcount << "\n \n" << endl;
//   }


  output -= totcount;



  // Colour and spin factors
  if(particles[0]->id()==-particles[1]->id()) {
    output *= 1./9.;
  }
  else  {
    output *= 1./24.;
  }
  if(mePartonData()[2]->id() == mePartonData()[3]->id()) output *= 0.5;
  // divided by 2 g_S^2
  return 0.5*output/gs2;
}

void MEPP2NeutralinoNeutralinoPowheg::constructVertex(tSubProPtr sub) {
//   //get particles
//   ParticleVector ext(4);
//   ext[0] = sub->incoming().first;
//   ext[1] = sub->incoming().second;
//   ext[2] = sub->outgoing()[0];
//   ext[3] = sub->outgoing()[1];
//   if( ext[0]->id() != mePartonData()[0]->id() ) swap(ext[0], ext[1]);
//   if( ext[2]->id() != mePartonData()[2]->id() ) swap(ext[2], ext[3]);

//   //First calculate wave functions with off-shell momenta
//   //to calculate correct spin information
//   vector<SpinorWaveFunction> sp;
//   SpinorWaveFunction(sp, ext[0], incoming, false);
//   vector<SpinorBarWaveFunction> sbar;
//   SpinorBarWaveFunction(sbar, ext[1], incoming, false);
//   vector<SpinorWaveFunction> spout;
//   SpinorWaveFunction(spout, ext[2], outgoing, false);
//   vector<SpinorBarWaveFunction> sbarout;
//   SpinorBarWaveFunction(sbarout, ext[3], outgoing, false);

//   //Need to use rescale momenta to calculate matrix element
//   cPDVector data(4);
//   vector<Lorentz5Momentum> momenta(4);
//   for( size_t i = 0; i < 4; ++i ) {
//     data[i] = ext[i]->dataPtr();
//     momenta[i] = ext[i]->momentum();
//   }
//   rescaleMomenta(momenta, data);
//   SpinorWaveFunction spr(rescaledMomenta()[0], data[0], incoming);
//   SpinorBarWaveFunction sbr(rescaledMomenta()[1], data[1],incoming);
//   SpinorWaveFunction spout1(rescaledMomenta()[2], data[2], outgoing);
//   SpinorBarWaveFunction sbarout2(rescaledMomenta()[3], data[3], outgoing);
//   for( unsigned int ihel = 0; ihel < 2; ++ihel ) {  
//     spr.reset(ihel);
//     sp[ihel] = spr;
//     sbr.reset(ihel);
//     sbar[ihel] = sbr;
//     spout1.reset(ihel);
//     spout[ihel] = spout1;
//     sbarout2.reset(ihel);
//     sbarout[ihel] = sbarout2;
//   }


//   qqbarME(sp, sbar, spout, sbarout, true);
//   HardVertexPtr hv = new_ptr(HardVertex());
//   hv->ME(me_);
//   for(unsigned int i = 0; i < 4; ++i )
//     dynamic_ptr_cast<SpinfoPtr>(ext[i]->spinInfo())->setProductionVertex(hv);  
}





