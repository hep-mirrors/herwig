// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2HiggsVBFPowheg class.
//

#include "MEPP2HiggsVBFPowheg.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "Herwig++/PDT/StandardMatchers.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.h"

using namespace Herwig;

MEPP2HiggsVBFPowheg::MEPP2HiggsVBFPowheg() 
  : scaleOpt_(1),  muF_(100.*GeV), scaleFact_(1.), contrib_(1), power_(0.1)
{}

int MEPP2HiggsVBFPowheg::nDim() const {
  return MEPP2HiggsVBF::nDim()+3;
}

bool MEPP2HiggsVBFPowheg::generateKinematics(const double * r) {
  // Born kinematics
  if(!MEPP2HiggsVBF::generateKinematics(r)) return false;
  // hadron and momentum fraction
  // set Q2 process momenta
  if(UseRandom::rnd()> 0.5) {
    _hadron = dynamic_ptr_cast<tcBeamPtr>(lastParticles().first->dataPtr());
    _xB = lastX1();
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
    _hadron = dynamic_ptr_cast<tcBeamPtr>(lastParticles().second->dataPtr());
    _xB = lastX2();
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
  int ndim=nDim();
  double rhomin = pow(1.-_xB,1.-power_); 
  double rho = r[ndim-1]*rhomin;

  _xp = 1.-pow(rho,1./(1.-power_));

  // zp 
  //double lambda = r[ndim-2]*rhomin;
  //_zp = 1.-pow(lambda,1./(1.-power_));
    _zp = r[ndim-2];
  // phi
  _phi = r[ndim-3]*Constants::twopi;
  jac_  = rhomin/(1.-power_)*pow(1.-_xp,power_);
  //jac_ *= pow(1.-_zp,power_)/(1.-power_);
  // NO as generating phi between 0 and 2pi and have a factor 1/2pi in the phase space
  // these cancel
  //jac_ *= 1./Constants::twopi;
  jacobian(jacobian()*jac_);
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

ClassDescription<MEPP2HiggsVBFPowheg>
 MEPP2HiggsVBFPowheg::initMEPP2HiggsVBFPowheg;
// Definition of the static class description member.

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
  double sinth(sqr(axis.x())+sqr(axis.y()));
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
  /*
  cerr << "Momenta prima del boost  _pb: "
                   << "(" << _pb.x()/GeV
                   << "," << _pb.y()/GeV
		   << "," << _pb.z()/GeV
		   << "," << _pb.t()/GeV
	           << ")" << 
          " _pc: " << "(" << _pc.x()/GeV
                   << "," << _pc.y()/GeV
		   << "," << _pc.z()/GeV
		   << "," << _pc.t()/GeV
	           << ")" <<
     " _pbother: " << "(" << _pbother.x()/GeV
                   << "," << _pbother.y()/GeV
		   << "," << _pbother.z()/GeV
		   << "," << _pbother.t()/GeV
	           << ")" <<
     " _pcother: " << "(" << _pcother.x()/GeV
                   << "," << _pcother.y()/GeV
		   << "," << _pcother.z()/GeV
                   << "," << _pcother.t()/GeV <<")" << endl; 
*/
  _pa *= rot;
  _q2  = -_pa.m2();
  
//   for(int i = 0; i < 4; i++){
//     if (_partons[i]->id() < 0){ 
//       _loMomenta[i].x() = - _loMomenta[i].x();
//       _loMomenta[i].y() = - _loMomenta[i].y();
//       _loMomenta[i].z() = - _loMomenta[i].z();
//       _loMomenta[i].e() = - _loMomenta[i].e(); 
//     }
//   }
  p1 = rot*_loMomenta[0];
  p2 = rot*_loMomenta[2];
  p1other = rot*_loMomenta[1];
  p2other = rot*_loMomenta[3];
  /*
  cerr << "Momenta dopo il boost - p1: "
                   << "(" << p1.x()/GeV
                   << "," << p1.y()/GeV
		   << "," << p1.z()/GeV
		   << "," << p1.t()/GeV
	           << ")" << 
           " p2: " << "(" << p2.x()/GeV
                   << "," << p2.y()/GeV
		   << "," << p2.z()/GeV
		   << "," << p2.t()/GeV
	           << ")" <<
      " p1other: " << "(" << p1other.x()/GeV
                   << "," << p1other.y()/GeV
		   << "," << p1other.z()/GeV
		   << "," << p1other.t()/GeV
	           << ")" <<
      " p2other: " << "(" << p2other.x()/GeV
                   << "," << p2other.y()/GeV
		   << "," << p2other.z()/GeV
                   << "," << p2other.t()/GeV << ")" <<  endl;
*/



  // scale and prefactors
  Energy2 mu2 = _q2;
  Energy Q = sqrt(_q2);
  double aS = 2.*SM().alphaS(mu2);
  double CFfact = 4./3.*aS/Constants::twopi;
  double TRfact = 1./2.*aS/Constants::twopi; 

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
  
//     if (_partons[0]->id() < 0)
//       nloMomenta[0] = - nloMomenta[0];
//     if (_partons[2]->id() < 0) 
//       nloMomenta[1] = - nloMomenta[1];
//     if (_partons[0]->id() > 0)
//       nloMomenta[2] = - nloMomenta[2];


  Lorentz5Momentum qnlo = nloMomenta[2]+nloMomenta[1]-nloMomenta[0];
  Lorentz5Momentum r1 = -nloMomenta[0]/x1;
  Lorentz5Momentum r2 =  nloMomenta[1]/x2;
  Lorentz5Momentum r3 = -nloMomenta[2]/x3;

  
  /* cerr << "Momenta dopo il boost - p1: "
                   << "(" << p1.x()/GeV
                   << "," << p1.y()/GeV
		   << "," << p1.z()/GeV
		   << "," << p1.t()/GeV
	           << ")" <<
" nloMomenta[0]: " << "(" << nloMomenta[0].x()/GeV
                   << "," << nloMomenta[0].y()/GeV
		   << "," << nloMomenta[0].z()/GeV
		   << "," << nloMomenta[0].t()/GeV
                   << ")" << 
           " r1: " << "(" << r1.x()/GeV
                   << "," << r1.y()/GeV
		   << "," << r1.z()/GeV
		   << "," << r1.t()/GeV
	           << ")" <<
      " p1other: " << "(" << p1other.x()/GeV
                   << "," << p1other.y()/GeV
		   << "," << p1other.z()/GeV
		   << "," << p1other.t()/GeV
	           << ")" << 
           " x1: " << x1 << endl;*/

  // LO + dipole subtracted virtual + collinear quark bit with LO pdf
  double virt = 1.+CFfact*(-4.5-1./3.*sqr(Constants::pi)+
			   1.5*log(_q2/mu2/(1.-_xB))+
			   2.*log(1.-_xB)*log(_q2/mu2)+
			   sqr(log(1.-_xB)));
  virt /= jac_;
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
  Energy2 mb2;
  // W
  if(_partons[0]->id()!=_partons[2]->id()) {
    mb2 = _mw2;
    c0L = 1;
    c0R = 0;
    c1L = 1;
    c1R = 0;
  }
  // Z
  else {
    mb2 = _mz2;
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
  }

  // Matrix element variables
  double G1 = sqr(c0L*c1L)+sqr(c0R*c1R);
  double G2 = sqr(c0L*c1R)+sqr(c0R*c1L);

  Energy4 term1,term2,term3,loME;
  if(_partons[0]->id()>0) {
    if(_partons[1]->id()>0) {
      term1 = loME(r1     ,p1other,qnlo+r1,p2other,G1,G2);
      term2 = loME(r2-qnlo,p1other,r2     ,p2other,G1,G2);
      term3 = loME(r3     ,p1other,qnlo+r3,p2other,G1,G2);
      loME  = loME(p1     ,p1other,p2     ,p2other,G1,G2);
    }
    else {
      term1 = loME(r1     ,p2other,qnlo+r1,p1other,G1,G2);
      term2 = loME(r2-qnlo,p2other,r2     ,p1other,G1,G2);
      term3 = loME(r3     ,p2other,qnlo+r3,p1other,G1,G2);
      loME  = loME(p1     ,p2other,p2     ,p1other,G1,G2);
    }
  }
  else {
    if(_partons[1]->id()>0) {
      term1 = loME(qnlo+r1,p1other,r1     ,p2other,G1,G2);
      term2 = loME(r2     ,p1other,r2-qnlo,p2other,G1,G2);
      term3 = loME(qnlo+r3,p1other,r3     ,p2other,G1,G2);
      loME  = loME(p2     ,p1other,p1     ,p2other,G1,G2);
    }
    else {
      term1 = loME(qnlo+r1,p2other,r1     ,p1other,G1,G2);
      term2 = loME(r2     ,p2other,r2-qnlo,p1other,G1,G2);
      term3 = loME(qnlo+r3,p2other,r3     ,p1other,G1,G2);
      loME  = loME(p2     ,p2other,p1     ,p1other,G1,G2);
    }
  }






  /*
  cerr << "Testing term2: " << (r2-qnlo)*p1other/GeV2 << " " << r2*p2other/GeV2 << " " << (r2-qnlo)*p2other/GeV2 << " " << r2*p1other/GeV2 << endl; */



  /*  cerr << "p1tilde: " << "(" << p1tilde.x()/GeV
                   << "," << p1tilde.y()/GeV
		   << "," << p1tilde.z()/GeV
		   << "," << p1tilde.t()/GeV
	           << ")" << 
           " p2tilde: " << "(" << p2tilde.x()/GeV
                   << "," << p2tilde.y()/GeV
		   << "," << p2tilde.z()/GeV
		   << "," << p2tilde.t()/GeV
	           << ")" <<
           " p1other: " << "(" << p1other.x()/GeV
                   << "," << p1other.y()/GeV
		   << "," << p1other.z()/GeV
		   << "," << p1other.t()/GeV
	           << ")" <<
           " p2other: " << "(" << p2other.x()/GeV
                   << "," << p2other.y()/GeV
		   << "," << p2other.z()/GeV
		   << "," << p2other.t()/GeV
                   << ")" << 
         " qnlo: " << "(" << qnlo.x()/GeV
                   << "," << qnlo.y()/GeV
		   << "," << qnlo.z()/GeV
		   << "," << qnlo.t()/GeV
                   << ")" << endl;*/

  /*
  cerr << "loME: " << loME/sqr(GeV2) << " G1: " << G1 << " G2: " << G2 << " p1*p1other: " << p1*p1other/GeV2 << " p2*p2other: " << p2*p2other/GeV2 << " p1*p2other: " << p1*p2other/GeV2 << " p2*p1other: " << p2*p1other/GeV2 << endl;
  */
  /* cerr << "loMetilde: " <<  loMEtilde/sqr(GeV2)  << " G1: " << G1 << " G2: " << G2 << " p1tilde*p1other: " << p1tilde*p1other/GeV2 << " p2tilde*p2other: " << p2tilde*p2other/GeV2 << " p1tilde*p2other: " << p1tilde*p2other/GeV2 << " p2tilde*p1other: " << p2tilde*p1other/GeV2 << endl;*/

  // dipoles kinematc variables
  /*  double x =0.;
  double z=1.;*/
    
  double x = 1.-nloMomenta[2]*nloMomenta[1]/
            ((nloMomenta[2]+nloMomenta[1])*
            nloMomenta[0]);
  double z = nloMomenta[1]*nloMomenta[0]/
            ((nloMomenta[2]+nloMomenta[1])*
            nloMomenta[0]);

  // q -> qg term
  double R1 = term1/loME;
  double R2 = sqr(x2)/(sqr(x2)+sqr(xT))*(term2/loME);
  double real1   = CFfact*qPDF/loPDF/_xp*1./((1.-_xp)*(1.-_zp))*
                    (R1+sqr(_xp)*(sqr(x2)+sqr(xT))*R2);

  double dipole1 = CFfact*
                   qPDF/loPDF/_xp*(sqr(x)+sqr(z))/
                   ((1.-x)*(1.-z)); 


  double realq   = (real1-dipole1);
   
   if(realq >100000000){
   cerr << "term1: " << term1/sqr(GeV2) << " term2: " << term2/sqr(GeV2) << " xp: " << _xp << " zp: " << _zp << " commfact: " << loME/sqr(GeV2) << " x: " << x << " z: " << z <<  " real1: " << real1 << " dipole1: " << dipole1 << " realq: "<< realq << endl;

    /*cerr << "term2: " <<term2/sqr(GeV2) << " (r2-qnlo)*p1other: " << (r2-qnlo)*p1other/GeV2 << " (r2*p2other): " << (r2*p2other)/GeV2 << " (r2-qnlo)*p2other: " << (r2-qnlo)*p2other/GeV2 << " r2*p1other/GeV2: " << r2*p1other/GeV2 << " loME: " << loME/sqr(GeV2) << " (p1*p1other): " << (p1*p1other)/GeV2 << " (p2*p2other): " << (p2*p2other)/GeV2 << "(p1*p2other): " << (p1*p2other)/GeV2 << " (p2*p1other):  " << (p2*p1other)/GeV2 << endl;
   
    cerr << " r2: " <<  r2/GeV
    <<" p2: " <<  p2/GeV
	           <<  
" nloMomenta[1]: " <<  nloMomenta[1].x()/GeV
                   <<  nloMomenta[1].y()/GeV
		   <<  nloMomenta[1].z()/GeV
		   <<  nloMomenta[1].t()/GeV
	           << 
" _pc: " << "(" << _pc.x()/GeV
                   << "," << _pc.y()/GeV
		   << "," << _pc.z()/GeV
		   << "," << _pc.t()/GeV
	           << ")" << 
           " zp: " << _zp <<             
           " xp: " << _xp << endl;*/


   }


   /*   
    cerr << "realq: " << realq << " x: " << x << " z: " << z << " term1: " << term1/loME << " term2: " << term2/loME  << " TERM1: " << ((r1*p1other)*((qnlo+r1)*p2other))/((p1*p1other)*(p2*p2other))<< " TEMR2: " <<(((r2-qnlo)*p1other)*(r2*p2other))/((p1*p1other)*(p2*p2other))<< "Estremo1: " << (r1*p1other)/(p1*p1other) << " " << ((qnlo+r1)*p2other)/(p2*p2other) << '\n';
*/   

    
   /* cerr << "qnlo: " << "(" << qnlo.x()/GeV
                   << "," << qnlo.y()/GeV
		   << "," << qnlo.z()/GeV
		   << "," << qnlo.t()/GeV
	           << ")" << 
           " r1: " << r1/GeV
	           << ")" <<
           " p1: " << "(" << p1.x()/GeV
                   << "," << p1.y()/GeV
		   << "," << p1.z()/GeV
		   << "," << p1.t()/GeV
	           << ")" <<
           " p2: " << "(" << p2.x()/GeV
                   << "," << p2.y()/GeV
		   << "," << p2.z()/GeV
		   << "," << p2.t()/GeV
	           << ")" <<
      " p1other: " << "(" << p1other.x()/GeV
                   << "," << p1other.y()/GeV
		   << "," << p1other.z()/GeV
		   << "," << p1other.t()/GeV
	           << ")" <<
      " p2other: " << "(" << p2other.x()/GeV
                   << "," << p2other.y()/GeV
		   << "," << p2other.z()/GeV
		   << "," << p2other.t()/GeV
                   << ")" << 

    " " << ((qnlo+r1)*p2other)/(p2*p2other) <<   endl;*/ 

  // g -> q qbar term
   
  double R3 = sqr(x3)/(sqr(x3)+sqr(xT))*(term3/loME);
         R2 = sqr(x2)/(sqr(x2)+sqr(xT))*(term2/loME);
  double real2 = TRfact*gPDF/loPDF/_xp*(1./(1.-_zp)+1./_zp)*
                 (sqr(_xp)*(sqr(x3)+sqr(xT))*R3+
                  sqr(_xp)*(sqr(x2)+sqr(xT))*R2);
   double dipole2 = TRfact*gPDF/loPDF/_xp*
                    ((sqr(x)+sqr(1.-x))/(1.-z)
		    +(sqr(x)+sqr(1.-x))/z);
   
   double realg = real2-dipole2;
      if(realg > 100000000){    
	cerr << "term3: " << term3/sqr(GeV2) << " term2: " << term2/sqr(GeV2) << " xp: " << _xp << " zp: " << _zp << " x: " << x << " z: " << z <<  " real2: " << real2 << " dipole2: " << dipole2 << " realg: "<< realg << endl;

	/* cerr << "term2: " <<term2/sqr(GeV2) << " (r2-qnlo)*p1other: " << (r2-qnlo)*p1other/GeV2 << " (r2*p2other): " << (r2*p2other)/GeV2 << " (r2-qnlo)*p2other: " << (r2-qnlo)*p2other/GeV2 << " r2*p1other/GeV2: " << r2*p1other/GeV2 << " loMEtilde: " << loMEtilde/sqr(GeV2) << " (p1tilde*p1other): " << (p1tilde*p1other)/GeV2 << " (p2tilde*p2other): " << (p2tilde*p2other)/GeV2 << "(p1tilde*p2other): " << (p1tilde*p2other)/GeV2 << " (p2tilde*p1other):  " << (p2tilde*p1other)/GeV2 << endl;

    cerr << 
           " r2: " << "(" << r2.x()/GeV
                   << "," << r2.y()/GeV
		   << "," << r2.z()/GeV
		   << "," << r2.t()/GeV
	           << ")" <<
      " p2tilde: " << "(" << p2tilde.x()/GeV
                   << "," << p2tilde.y()/GeV
		   << "," << p2tilde.z()/GeV
		   << "," << p2tilde.t()/GeV
	           << ")" << 
" nloMomenta[1]: " << "(" << nloMomenta[1].x()/GeV
                   << "," << nloMomenta[1].y()/GeV
		   << "," << nloMomenta[1].z()/GeV
		   << "," << nloMomenta[1].t()/GeV
 	           << ")" <<
          " _pc: " << "(" << _pc.x()/GeV
                   << "," << _pc.y()/GeV
		   << "," << _pc.z()/GeV
		   << "," << _pc.t()/GeV
	           << ")" <<  
           " zp: " << _zp <<  
           " xp: " << _xp << endl;*/

   }
   //cerr << "Testing cancellation realg: " << realg << " x: " << x << " z: " << z << '\n';

  // return the full result
  double wgt = virt+((collq+collg)/loPDF+realq+realg);
  return contrib_ == 1 ? max(0.,wgt) : max(0.,-wgt);
}

void MEPP2HiggsVBFPowheg::doinit() {
  MEPP2HiggsVBF::doinit();
  // electroweak parameters
  _mz2 = sqr(getParticleData(ParticleID::Z0)->mass());
  _mw2 = sqr(getParticleData(ParticleID::Wplus)->mass());
  tcPDPtr gluon = getParticleData(ParticleID::g);
}

Energy4 MEPP2HiggsVBFPowheg::loME(const Lorentz5Momentum &p1,
				  const Lorentz5Momentum &p2,
				  const Lorentz5Momentum &q1,
				  const Lorentz5Momentum &q2,
				  double G1, double G2) const {
  return G1*(p1*p2)*(q1*q2) + G2*(p1*q2)*(q1*p2);
}

