// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2HiggsVBFPowheg class.
//

#include "MEPP2HiggsVBFPowheg.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "Herwig/PDT/StandardMatchers.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.h"
//ACDC test
//#include "ThePEG/Repository/UseRandom.h"

using namespace Herwig;

MEPP2HiggsVBFPowheg::MEPP2HiggsVBFPowheg() 
  : scaleOpt_(1),  muF_(100.*GeV), scaleFact_(1.), contrib_(1), power_(0.1)
{}

int MEPP2HiggsVBFPowheg::nDim() const {
  //ACDC test 
  // return 0;
    return MEPP2HiggsVBF::nDim()+3;
}
//ACDC test (delete r)
bool MEPP2HiggsVBFPowheg::generateKinematics(const double * r) {
  //bool MEPP2HiggsVBFPowheg::generateKinematics(const double * ) {
  int ndim = nDim();
  //ACDC test
   double r3 = r[ndim-3];
   //double r3 = UseRandom::rnd();

  // Born kinematics
  //ACDC test (comment next line out)
   if(!MEPP2HiggsVBF::generateKinematics(r)) return false;

  // hadron and momentum fraction
  // set Q2 process momenta

        if(r3 > 0.5) {
        r3 = 2.*(r3-0.5);
    if(lastPartons().first ->dataPtr()==mePartonData()[0]&&
       lastPartons().second->dataPtr()==mePartonData()[1]) {
      _hadron = dynamic_ptr_cast<tcBeamPtr>(lastParticles().first->dataPtr());
      _xB = lastX1();
          }
    else {
      _hadron = dynamic_ptr_cast<tcBeamPtr>(lastParticles().second->dataPtr());
      _xB = lastX2();
      }
    _partons[0] = mePartonData()[0]; 
    _partons[1] = mePartonData()[1]; 
    _partons[4] = mePartonData()[4];
    if(!swapOrder()) {
      _pb = meMomenta()[0];
      _pc = meMomenta()[2];
      _pbother = meMomenta()[1];
      _pcother = meMomenta()[3];
      _partons[2] = mePartonData()[2];
      _partons[3] = mePartonData()[3];
    }
    else {
      _pb = meMomenta()[0];
      _pc = meMomenta()[3];
      _pbother = meMomenta()[1];
      _pcother = meMomenta()[2];
      _partons[2] = mePartonData()[3];
      _partons[3] = mePartonData()[2];
    }
      }
else {
      r3 = 2.*r3;
    if(lastPartons().first ->dataPtr()==mePartonData()[0]&&
       lastPartons().second->dataPtr()==mePartonData()[1]) {
      _hadron = dynamic_ptr_cast<tcBeamPtr>(lastParticles().second->dataPtr());
      _xB = lastX2();
    }
    else {
      _hadron = dynamic_ptr_cast<tcBeamPtr>(lastParticles().first->dataPtr());
      _xB = lastX1();
    }
    _partons[0] = mePartonData()[1];
    _partons[1] = mePartonData()[0];
    _partons[4] = mePartonData()[4];
    if(!swapOrder()) {
      _pb = meMomenta()[1];
      _pc = meMomenta()[3];
      _pbother = meMomenta()[0];
      _pcother = meMomenta()[2];
      _partons[2] = mePartonData()[3];
      _partons[3] = mePartonData()[2];
    }
    else {
      _pb = meMomenta()[1];
      _pc = meMomenta()[2];
      _pbother = meMomenta()[0];
      _pcother = meMomenta()[3];
      _partons[2] = mePartonData()[2];
      _partons[3] = mePartonData()[3];
    }
    }
  // LO Momenta assignment
  _loMomenta[0] = _pb;
  _loMomenta[1] = _pbother;
  _loMomenta[2] = _pc;
  _loMomenta[3] = _pcother;
  _pa =  _pc-_pb;
  // xp
  double rhomin = pow(1.-_xB,1.-power_);

  //ACDC test 
  double rho = r[ndim-1]*rhomin;
  //double rho = UseRandom::rnd()*rhomin;

  _xp = 1.-pow(rho,1./(1.-power_));
  // zp
  //ACDC test  
    _zp = r[ndim-2];
    //_zp = UseRandom::rnd();

  // phi
  _phi = r3*Constants::twopi;
  jac_  = rhomin/(1.-power_)*pow(1.-_xp,power_);
  return true;
}

IBPtr MEPP2HiggsVBFPowheg::clone() const {
  return new_ptr(*this);
}

IBPtr MEPP2HiggsVBFPowheg::fullclone() const {
  return new_ptr(*this);
}


void MEPP2HiggsVBFPowheg::persistentOutput(PersistentOStream & os) const {
  os << ounit(muF_,GeV) << scaleFact_ << scaleOpt_ << contrib_
     << ounit(_mz2,GeV2) << ounit(_mw2,GeV2) << power_;
}

void MEPP2HiggsVBFPowheg::persistentInput(PersistentIStream & is, int) {
  is >> iunit(muF_,GeV) >> scaleFact_ >> scaleOpt_ >> contrib_
     >> iunit(_mz2,GeV2) >> iunit(_mw2,GeV2) >> power_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEPP2HiggsVBFPowheg,MEPP2HiggsVBF>
describeHerwigMEPP2HiggsVBFPowheg("Herwig::MEPP2HiggsVBFPowheg", "HwMEHadron.so HwPowhegMEHadron.so");

void MEPP2HiggsVBFPowheg::Init() {

  static ClassDocumentation<MEPP2HiggsVBFPowheg> documentation
    ("The MENeutralCurrentDISPowheg class implements the NLO matrix element"
     " for neutral current DIS in the Powheg scheme.");

  static Switch<MEPP2HiggsVBFPowheg,unsigned int> interfaceContribution
    ("Contribution",
     "Which contributions to the cross section to include",
     &MEPP2HiggsVBFPowheg::contrib_, 1, false, false);
  static SwitchOption interfaceContributionLeadingOrder
    (interfaceContribution,
     "LeadingOrder",
     "Just generate the leading order cross section",
     0);
  static SwitchOption interfaceContributionPositiveNLO
    (interfaceContribution,
     "PositiveNLO",
     "Generate the positive contribution to the full NLO cross section",
     1);
  static SwitchOption interfaceContributionNegativeNLO
    (interfaceContribution,
     "NegativeNLO",
     "Generate the negative contribution to the full NLO cross section",
     2);

  static Switch<MEPP2HiggsVBFPowheg,unsigned int> interfaceScaleOption
    ("ScaleOption",
     "Option for the choice of factorization (and renormalization) scale",
     &MEPP2HiggsVBFPowheg::scaleOpt_, 1, false, false);
  static SwitchOption interfaceDynamic
    (interfaceScaleOption,
     "Dynamic",
     "Dynamic factorization scale equal to the current sqrt(sHat())",
     1);
  static SwitchOption interfaceFixed
    (interfaceScaleOption,
     "Fixed",
     "Use a fixed factorization scale set with FactorizationScaleValue",
     2);

  static Parameter<MEPP2HiggsVBFPowheg,Energy> interfaceFactorizationScale
    ("FactorizationScale",
     "Value to use in the event of a fixed factorization scale",
     &MEPP2HiggsVBFPowheg::muF_, GeV, 100.0*GeV, 1.0*GeV, 500.0*GeV,
     true, false, Interface::limited);

  static Parameter<MEPP2HiggsVBFPowheg,double> interfaceScaleFactor
    ("ScaleFactor",
     "The factor used before Q2 if using a running scale",
     &MEPP2HiggsVBFPowheg::scaleFact_, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<MEPP2HiggsVBFPowheg,double> interfaceSamplingPower
    ("SamplingPower",
     "Power for the sampling of xp",
     &MEPP2HiggsVBFPowheg::power_, 0.6, 0.0, 1.,
     false, false, Interface::limited);

}

Energy2 MEPP2HiggsVBFPowheg::scale() const {
  return scaleOpt_ == 1 ? 
    sqr(scaleFact_)*MEPP2HiggsVBF::scale() : sqr(scaleFact_*muF_);
}

CrossSection MEPP2HiggsVBFPowheg::dSigHatDR() const {
  return NLOWeight()*MEPP2HiggsVBF::dSigHatDR();
}

double MEPP2HiggsVBFPowheg::NLOWeight() const {
  // If only leading order is required return 1:
  if(contrib_==0) return 1.;
  // Boost
  Axis axis(_pa.vect().unit());
  LorentzRotation rot;
  double sinth(sqrt(sqr(axis.x())+sqr(axis.y())));
  rot = LorentzRotation();
  if(axis.perp2()>1e-20) {
    rot.setRotate(-acos(axis.z()),Axis(-axis.y()/sinth,axis.x()/sinth,0.));
    rot.rotateX(Constants::pi);
  }
  if(abs(1.-_pa.e()/_pa.vect().mag())>1e-6) 
    rot.boostZ(_pa.e()/_pa.vect().mag());
  _pb *= rot;
  if(_pb.perp2()/GeV2>1e-20) {
    Boost trans = -1./_pb.e()*_pb.vect();
    trans.setZ(0.);
    rot.boost(trans);
  }
  // momenta of particles
  Lorentz5Momentum p1, p2,p1other,p2other;
  _pa *= rot;
  _q2  = -_pa.m2();
  p1 = rot*_loMomenta[0];
  p2 = rot*_loMomenta[2];
  p1other = rot*_loMomenta[1];
  p2other = rot*_loMomenta[3];

  // scale and prefactors
  Energy2 mu2 = scale();
  Energy Q = sqrt(_q2);
  //  double aS = 2.*SM().alphaS(mu2);
  double aS = 2*0.113291076960184;
  double CFfact = 4./3.*aS/Constants::twopi;
  double TRfact = 1./2.*aS/Constants::twopi; 
  double wgt = 0.;
  // Breit frame variables
  double x1 = -1./_xp;
  double x2 = 1.-(1.-_zp)/_xp;
  double x3 = 2.+x1-x2;
  double xT = 2*sqrt((1.-_xp)*(1.-_zp)*_zp/_xp);
  
  vector<Lorentz5Momentum>  nloMomenta;
  nloMomenta.resize(3);

  nloMomenta[0] = Lorentz5Momentum(ZERO,ZERO,-0.5*Q*x1,-0.5*Q*x1);
  nloMomenta[1] = Lorentz5Momentum( 0.5*Q*xT*cos(_phi),  0.5*Q*xT*sin(_phi),
				   -0.5*Q*x2, 0.5*Q*sqrt(sqr(xT)+sqr(x2)));
  nloMomenta[2] = Lorentz5Momentum(-0.5*Q*xT*cos(_phi), -0.5*Q*xT*sin(_phi),
				   -0.5*Q*x3, 0.5*Q*sqrt(sqr(xT)+sqr(x3)));
  Lorentz5Momentum qnlo = nloMomenta[2]+nloMomenta[1]-nloMomenta[0];
  Lorentz5Momentum r1 = -nloMomenta[0]/x1;
  Lorentz5Momentum r2 =  nloMomenta[1]/x2;
  Lorentz5Momentum r3 = -nloMomenta[2]/x3;

  // LO + dipole subtracted virtual + collinear quark bit with LO pdf
  double virt = 1.+CFfact*(-4.5-1./3.*sqr(Constants::pi)+
			   1.5*log(_q2/mu2/(1.-_xB))+
			   2.*log(1.-_xB)*log(_q2/mu2)+
			   sqr(log(1.-_xB)));
  // PDF from leading-order
  double loPDF = 
    _hadron->pdf()->xfx(_hadron,_partons[0],mu2,_xB)/_xB;
  // NLO gluon PDF

  tcPDPtr gluon = getParticleData(ParticleID::g);
  double gPDF   = 
    _hadron->pdf()->xfx(_hadron,gluon,mu2,_xB/_xp)*_xp/_xB;
  // NLO quark PDF
  double qPDF   = 
    _hadron->pdf()->xfx(_hadron,_partons[0],mu2,_xB/_xp)*_xp/_xB;
  // collinear counterterms
  // gluon
  double collg = 
    TRfact/_xp*gPDF*(2.*_xp*(1.-_xp)+(sqr(_xp)+sqr(1.-_xp))*
    log((1.-_xp)*_q2/_xp/mu2));
  // quark
  double collq = 
    CFfact/_xp*qPDF*(1-_xp-2./(1.-_xp)*log(_xp)-(1.+_xp)*
    log((1.-_xp)/_xp*_q2/mu2))+
    CFfact/_xp*(qPDF-_xp*loPDF)*(2./(1.-_xp)*
    log(_q2*(1.-_xp)/mu2)-1.5/(1.-_xp));
  // Electroweak coefficients
  double c0L,c1L,c0R,c1R;
  //Energy2 mb2;
  // W
  if(_partons[0]->id()!=_partons[2]->id()) {
    //mb2 = _mw2;
    c0L = sqrt(0.5);
    c0R = 0;
    c1L = sqrt(0.5);
    c1R = 0;
  }
  // Z
  else {
    //mb2 = _mz2;
    if(abs(_partons[0]->id())%2==0) {
      c0L = 
	generator()->standardModel()->vu()+
	generator()->standardModel()->au();
      c0R =
	generator()->standardModel()->vu()-
	generator()->standardModel()->au();
    }
    else {
      c0L = 
	generator()->standardModel()->vd()+
	generator()->standardModel()->ad();
      c0R =
	generator()->standardModel()->vd()-
	generator()->standardModel()->ad();
    }
    if(abs(_partons[1]->id())%2==0) {
      c1L = 
	generator()->standardModel()->vu()+
	generator()->standardModel()->au();
      c1R =
	generator()->standardModel()->vu()-
	generator()->standardModel()->au();
    }
    else {
      c1L = 
	generator()->standardModel()->vd()+
	generator()->standardModel()->ad();
      c1R =
	generator()->standardModel()->vd()-
	generator()->standardModel()->ad();
    }
    c0L *= 0.25;
    c0R *= 0.25;
    c1L *= 0.25;
    c1R *= 0.25;
  }
  // Matrix element variables
  double G1 = sqr(c0L*c1L)+sqr(c0R*c1R);
  double G2 = sqr(c0L*c1R)+sqr(c0R*c1L);
  Energy4 term1,term2,term3,loME;
  if(_partons[0]->id()>0) {
    if(_partons[1]->id()>0) {
      term1 = loMatrixElement(r1     ,p1other,qnlo+r1,p2other,G1,G2);
      term2 = loMatrixElement(r2-qnlo,p1other,r2     ,p2other,G1,G2);
      term3 = loMatrixElement(r3     ,p1other,qnlo+r3,p2other,G1,G2);
      loME  = loMatrixElement(p1     ,p1other,p2     ,p2other,G1,G2);
    }
    else {
      term1 = loMatrixElement(r1     ,p2other,qnlo+r1,p1other,G1,G2);
      term2 = loMatrixElement(r2-qnlo,p2other,r2     ,p1other,G1,G2);
      term3 = loMatrixElement(r3     ,p2other,qnlo+r3,p1other,G1,G2);
      loME  = loMatrixElement(p1     ,p2other,p2     ,p1other,G1,G2);
    }
  }
  else {
    if(_partons[1]->id()>0) {
      term1 = loMatrixElement(qnlo+r1,p1other,r1     ,p2other,G1,G2);
      term2 = loMatrixElement(r2     ,p1other,r2-qnlo,p2other,G1,G2);
      term3 = loMatrixElement(qnlo+r3,p1other,r3     ,p2other,G1,G2);
      loME  = loMatrixElement(p2     ,p1other,p1     ,p2other,G1,G2);
    }
    else {
      term1 = loMatrixElement(qnlo+r1,p2other,r1     ,p1other,G1,G2);
      term2 = loMatrixElement(r2     ,p2other,r2-qnlo,p1other,G1,G2);
      term3 = loMatrixElement(qnlo+r3,p2other,r3     ,p1other,G1,G2);
      loME  = loMatrixElement(p2     ,p2other,p1     ,p1other,G1,G2);
    }
  }
  if(1-_xp > 1e-10 && 1.-_zp > 1e-10){
  // q -> qg term
  double real1   = (term1+sqr(_xp)*sqr(x2)*term2)/loME;
  double dipole1 = (sqr(_xp)+sqr(_zp));
  double realq   = 
    CFfact*qPDF/loPDF/_xp/((1.-_xp)*(1.-_zp))*(real1-dipole1);
  
  // g -> q qbar term
  double real2 = sqr(_xp)/loME*(sqr(x2)*term2+sqr(x3)*term3);
  double dipole2 = sqr(_xp)+sqr(1.-_xp);
  double realg = TRfact*gPDF/loPDF/_xp/(1.-_zp)*(real2-dipole2);
     
  // Return The Full Result
  wgt = virt+jac_*((collq+collg)/loPDF+realq+realg);
  }


  return contrib_ == 1 ? max(0.,wgt) : max(0.,-wgt);
}

void MEPP2HiggsVBFPowheg::doinit() {
  MEPP2HiggsVBF::doinit();
  // electroweak parameters
  _mz2 = sqr(getParticleData(ParticleID::Z0)->mass());
  _mw2 = sqr(getParticleData(ParticleID::Wplus)->mass());
  tcPDPtr gluon = getParticleData(ParticleID::g);
}

Energy4 MEPP2HiggsVBFPowheg::
loMatrixElement(const Lorentz5Momentum &p1,
		const Lorentz5Momentum &p2,
		const Lorentz5Momentum &q1,
		const Lorentz5Momentum &q2,
		double G1, double G2) const {
  return G1*(p1*p2)*(q1*q2) + G2*(p1*q2)*(q1*p2);
}



























//   double g;
//     g = 1./sqrt(SM().sin2ThetaW());
//     g = 1./sqrt((1.-SM().sin2ThetaW())*SM().sin2ThetaW());


//   vector<SpinorWaveFunction> f1,f2;
//   vector<SpinorBarWaveFunction> a1,a2;
//   Lorentz5Momentum phiggs = rot*meMomenta()[4];
//   ScalarWaveFunction higgs(phiggs,mePartonData()[4],1.,outgoing);
//   SpinorWaveFunction fin1,fin2;
//   SpinorBarWaveFunction ain1,ain2;
//   if(_partons[0]->id()>0) {
//     fin1 =    SpinorWaveFunction(nloMomenta[0],_partons[0],incoming);
//     ain1 = SpinorBarWaveFunction(nloMomenta[1],_partons[2],outgoing);
//   }
//   else {
//     fin1 =     SpinorWaveFunction(nloMomenta[1],_partons[2],outgoing);
//     ain1 =  SpinorBarWaveFunction(nloMomenta[0],_partons[0],incoming);
//   }
//   if(_partons[1]->id()>0) {
//     fin2 =    SpinorWaveFunction(p1other,_partons[1],incoming);
//     ain2 = SpinorBarWaveFunction(p2other,_partons[3],outgoing);
//   }
//   else {
//     fin2 =     SpinorWaveFunction(p2other,_partons[3],outgoing);
//     ain2 =  SpinorBarWaveFunction(p1other,_partons[1],incoming);
//   }
//   VectorWaveFunction gwave(nloMomenta[2],gluon,outgoing);
//   vector<VectorWaveFunction> g1;
//   for(unsigned int ix=0;ix<2;++ix) {
//     fin1.reset(ix); f1.push_back(fin1);
//     fin2.reset(ix); f2.push_back(fin2);
//     ain1.reset(ix); a1.push_back(ain1);
//     ain2.reset(ix); a2.push_back(ain2);
//     gwave.reset(2*ix); g1.push_back(gwave);
//   }
//   AbstractFFVVertexPtr vertex[2];
//   tcPDPtr vec[2];
//   for(unsigned int ix=0;ix<2;++ix) {
//     int icharge;
//     icharge = _partons[ix]->iCharge()-_partons[ix+2]->iCharge();
//     if(icharge==0)     vec[ix] = _z0;
//     else if(icharge>0) vec[ix] = _wplus;
//     else               vec[ix] = _wminus;
//     vertex[ix] = vec[ix]==_z0 ? _vertexFFZ : _vertexFFW;
//   }
//   tcHwSMPtr hwsm= dynamic_ptr_cast<tcHwSMPtr>(standardModel());
//   AbstractFFVVertexPtr gvertex = hwsm->vertexFFG();
//   VectorWaveFunction inter[2];
//   Complex diag[2];
//   double me(0.);
//   nloMomenta[0].setMass(ZERO);
//   nloMomenta[1].setMass(ZERO);
//   nloMomenta[2].setMass(ZERO);
//   for(unsigned int i1=0;i1<2;++i1) {
//     unsigned int i2=i1;
//     // wavefunction for the 1st intermediate vector
//     inter[0] = vertex[1]->evaluate(scale(),1,vec[1],f2[i1],a2[i2]);
//     for(unsigned int i3=0;i3<2;++i3) {
//       unsigned int i4=i3;
//       for(unsigned int ig=0;ig<2;++ig) {
// 	SpinorWaveFunction soff = 
// 	  gvertex->evaluate(scale(),5,f1[i3].getParticle(),
// 			    f1[i3],g1[ig]);
// 	SpinorBarWaveFunction aoff = 
// 	  gvertex->evaluate(scale(),5,a1[i4].getParticle(),
// 			    a1[i4],g1[ig]);
// 	inter[1] = vertex[0]->evaluate(scale(),1,vec[0],soff,a1[i4]);
// 	diag[0] = _vertexWWH->evaluate(scale(),inter[0],inter[1],higgs);
// 	inter[1] = vertex[0]->evaluate(scale(),1,vec[0],f1[i3],aoff);
// 	diag[1] = _vertexWWH->evaluate(scale(),inter[0],inter[1],higgs);
// // 	cerr << "testing helicity "
// // 	     << i1 << " " << i2 << " " << i3 << " " << i4 << " " 
// // 	     << ig << " " << (diag[0]+diag[1])/(abs(diag[0])+abs(diag[1]))
// // 	     << "\n";
	
// 	me += norm(diag[0]+diag[1]);
//       }
//     }
//   }
//   // spin factor
//   me *=0.25;






//   cerr << "testing the quark emission "
//        << fact*8.*Constants::pi*SM().alphaS(scale())/(-1.-x1)/(1.-x2)/_q2
//     *(sqr(x1)*term1+sqr(x2)*term2)*sqr(MeV2)/me << "\n";










//   vector<SpinorWaveFunction> f1,f2;
//   vector<SpinorBarWaveFunction> a1,a2;
//   Lorentz5Momentum phiggs = rot*meMomenta()[4];
//   ScalarWaveFunction higgs(phiggs,mePartonData()[4],1.,outgoing);
//   SpinorWaveFunction fin1,fin2;
//   SpinorBarWaveFunction ain1,ain2;
//   if(_partons[0]->id()>0) {
//     fin1 =    SpinorWaveFunction(nloMomenta[2],_partons[0]->CC(),outgoing);
//     ain1 = SpinorBarWaveFunction(nloMomenta[1],_partons[2],outgoing);
//   }
//   else {
//     fin1 =     SpinorWaveFunction(nloMomenta[1],_partons[2],outgoing);
//     ain1 =  SpinorBarWaveFunction(nloMomenta[2],_partons[0]->CC(),outgoing);
//   }
//   if(_partons[1]->id()>0) {
//     fin2 =    SpinorWaveFunction(p1other,_partons[1],incoming);
//     ain2 = SpinorBarWaveFunction(p2other,_partons[3],outgoing);
//   }
//   else {
//     fin2 =     SpinorWaveFunction(p2other,_partons[3],outgoing);
//     ain2 =  SpinorBarWaveFunction(p1other,_partons[1],incoming);
//   }
//   VectorWaveFunction gwave(nloMomenta[0],gluon,incoming);
//   vector<VectorWaveFunction> g1;
//   for(unsigned int ix=0;ix<2;++ix) {
//     fin1.reset(ix); f1.push_back(fin1);
//     fin2.reset(ix); f2.push_back(fin2);
//     ain1.reset(ix); a1.push_back(ain1);
//     ain2.reset(ix); a2.push_back(ain2);
//     gwave.reset(2*ix);
//     g1.push_back(gwave);
//   }
//   AbstractFFVVertexPtr vertex[2];
//   tcPDPtr vec[2];
//   for(unsigned int ix=0;ix<2;++ix) {
//     int icharge;
//     icharge = _partons[ix]->iCharge()-_partons[ix+2]->iCharge();
//     if(icharge==0)     vec[ix] = _z0;
//     else if(icharge>0) vec[ix] = _wplus;
//     else               vec[ix] = _wminus;
//     vertex[ix] = vec[ix]==_z0 ? _vertexFFZ : _vertexFFW;
//   }
//   tcHwSMPtr hwsm= dynamic_ptr_cast<tcHwSMPtr>(standardModel());
//   AbstractFFVVertexPtr gvertex = hwsm->vertexFFG();
//   VectorWaveFunction inter[2];
//   Complex diag[2];
//   double me(0.);
//   nloMomenta[0].setMass(ZERO);
//   nloMomenta[1].setMass(ZERO);
//   nloMomenta[2].setMass(ZERO);
//   unsigned int o[2]={1,0};
//   for(unsigned int i1=0;i1<2;++i1) {
//     unsigned int i2=i1;
//     // wavefunction for the 1st intermediate vector
//     inter[0] = vertex[1]->evaluate(scale(),1,vec[1],f2[i1],a2[i2]);
//     for(unsigned int i3=0;i3<2;++i3) {
//       unsigned int i4=o[i3];
//       for(unsigned int ig=0;ig<2;++ig) {
// 	SpinorWaveFunction soff = 
// 	  gvertex->evaluate(scale(),5,f1[i3].getParticle(),
// 			    f1[i3],g1[ig]);
// 	SpinorBarWaveFunction aoff = 
// 	  gvertex->evaluate(scale(),5,a1[i4].getParticle(),
// 			    a1[i4],g1[ig]);
// 	inter[1] = vertex[0]->evaluate(scale(),1,vec[0],soff,a1[i4]);
// 	diag[0] = _vertexWWH->evaluate(scale(),inter[0],inter[1],higgs);
// 	inter[1] = vertex[0]->evaluate(scale(),1,vec[0],f1[i3],aoff);
// 	diag[1] = _vertexWWH->evaluate(scale(),inter[0],inter[1],higgs);
// // 	cerr << "testing helicity "
// // 	     << i1 << " " << i2 << " " << i3 << " " << i4 << " " 
// // 	     << ig << " " << (diag[0]+diag[1])/(abs(diag[0])+abs(diag[1]))
// // 	     << "\n";
	
// 	me += norm(diag[0]+diag[1]);
//       }
//     }
//   }
//   // spin factor
//   me *=0.25;
//   cerr << "testing the gluon emission "
//        << fact*8.*Constants::pi*SM().alphaS(scale())/(1.-x2)/(1.-x3)/_q2
//     *(sqr(x3)*term3+sqr(x2)*term2)*sqr(MeV2)/me << "\n";













//   Energy2 D1 = -_q2-mb2;
//   Energy2 D2 = (p1other-p2other).m2()-mb2;
//   double e = sqrt(4.*Constants::pi*SM().alphaEM(scale()));

//   InvEnergy6 fact = 4.*pow(e*g,6)*mb2/sqr(D1)/sqr(D2);

//   cerr << "testing LO ME in NLO code "
//        << fact*loME*MeV2/_mestore << "\n";



//   vector<SpinorWaveFunction> f1,f2;
//   vector<SpinorBarWaveFunction> a1,a2;
//   Lorentz5Momentum phiggs = rot*meMomenta()[4];
//   ScalarWaveFunction higgs(phiggs,mePartonData()[4],1.,outgoing);
//   SpinorWaveFunction fin1,fin2;
//   SpinorBarWaveFunction ain1,ain2;
//   if(_partons[0]->id()>0) {
//     fin1 =    SpinorWaveFunction(p1,_partons[0],incoming);
//     ain1 = SpinorBarWaveFunction(p2,_partons[2],outgoing);
//   }
//   else {
//     fin1 =     SpinorWaveFunction(p2,_partons[2],outgoing);
//     ain1 =  SpinorBarWaveFunction(p1,_partons[0],incoming);
//   }
//   if(_partons[1]->id()>0) {
//     fin2 =    SpinorWaveFunction(p1other,_partons[1],incoming);
//     ain2 = SpinorBarWaveFunction(p2other,_partons[3],outgoing);
//   }
//   else {
//     fin2 =     SpinorWaveFunction(p2other,_partons[3],outgoing);
//     ain2 =  SpinorBarWaveFunction(p1other,_partons[1],incoming);
//   }
//   for(unsigned int ix=0;ix<2;++ix) {
//     fin1.reset(ix); f1.push_back(fin1);
//     fin2.reset(ix); f2.push_back(fin2);
//     ain1.reset(ix); a1.push_back(ain1);
//     ain2.reset(ix); a2.push_back(ain2);
//   }
//   AbstractFFVVertexPtr vertex[2];
//   tcPDPtr vec[2];
//   for(unsigned int ix=0;ix<2;++ix) {
//     int icharge;
//     icharge = _partons[ix]->iCharge()-_partons[ix+2]->iCharge();
//     if(icharge==0)     vec[ix] = _z0;
//     else if(icharge>0) vec[ix] = _wplus;
//     else               vec[ix] = _wminus;
//     vertex[ix] = vec[ix]==_z0 ? _vertexFFZ : _vertexFFW;
//   }
//   VectorWaveFunction inter[2];
//   Complex diag;
//   double me(0.);
//   for(unsigned int i1=0;i1<2;++i1) {
//     for(unsigned int i2=0;i2<2;++i2) {
//       // wavefunction for the 1st intermediate vector
//       inter[0] = vertex[0]->evaluate(scale(),1,vec[1],f2[i1],a2[i2]);
//       for(unsigned int i3=0;i3<2;++i3) {
// 	for(unsigned int i4=0;i4<2;++i4) {
// 	  // wavefunction for the 2nd intermediate vector
// 	  inter[1] = vertex[1]->evaluate(scale(),1,vec[0],f1[i3],a1[i4]);
// 	  // matrix element
// 	  diag = _vertexWWH->evaluate(scale(),inter[0],inter[1],higgs);
// 	  me += norm(diag);
// 	}
//       }
//     }
//   }
//   // spin factor
//   me *=0.25;

//   cerr << "testing helicity computation " << me/_mestore << "\n";
