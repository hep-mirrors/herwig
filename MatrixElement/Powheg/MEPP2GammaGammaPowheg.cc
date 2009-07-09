// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2GammaGammaPowheg class.
//

#include "MEPP2GammaGammaPowheg.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDF/BeamParticleData.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include <numeric>

using namespace Herwig;

/**
 *  Typedef for BeamParticleData
 */
typedef Ptr<BeamParticleData>::transient_const_pointer tcBeamPtr;

MEPP2GammaGammaPowheg::MEPP2GammaGammaPowheg()  
  : contrib_(1), scaleopt_(0),
    fixedScale_(100.*GeV), scaleFact_(1.), power_(0.1)
{}

IBPtr MEPP2GammaGammaPowheg::clone() const {
  return new_ptr(*this);
}

IBPtr MEPP2GammaGammaPowheg::fullclone() const {
  return new_ptr(*this);
}

void MEPP2GammaGammaPowheg::persistentOutput(PersistentOStream & os) const {
  os << contrib_ << scaleopt_ << ounit(fixedScale_,GeV) << scaleFact_
     << gluon_ << QEDVertex_ << QCDVertex_ << power_;
}

void MEPP2GammaGammaPowheg::persistentInput(PersistentIStream & is, int) {
  is >> contrib_ >> scaleopt_ >> iunit(fixedScale_,GeV) >> scaleFact_
     >> gluon_ >> QEDVertex_ >> QCDVertex_ >> power_;
}

ClassDescription<MEPP2GammaGammaPowheg> MEPP2GammaGammaPowheg::initMEPP2GammaGammaPowheg;
// Definition of the static class description member.

void MEPP2GammaGammaPowheg::Init() {

  static ClassDocumentation<MEPP2GammaGammaPowheg> documentation
    ("The MEPP2GammaGammaPowheg class implements the NLO matrix elements"
     "for photon pair production");

   static Switch<MEPP2GammaGammaPowheg,unsigned int> interfaceContribution
    ("Contribution",
     "Which contributions to the cross section to include",
     &MEPP2GammaGammaPowheg::contrib_, 1, false, false);
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

  static Switch<MEPP2GammaGammaPowheg,unsigned int> interfaceFactorizationScaleOption
    ("FactorizationScaleOption",
     "Option for the scale to be used",
     &MEPP2GammaGammaPowheg::scaleopt_, 0, false, false);
  static SwitchOption interfaceScaleOptionFixed
    (interfaceFactorizationScaleOption,
     "Fixed",
     "Use a fixed scale",
     0);
  static SwitchOption interfaceScaleOptionsHat
    (interfaceFactorizationScaleOption,
     "Dynamic",
     "Used sHat as the scale",
     1);

  static Parameter<MEPP2GammaGammaPowheg,Energy> interfaceFactorizationScaleValue
    ("FactorizationScaleValue",
     "The fixed scale to use if required",
     &MEPP2GammaGammaPowheg::fixedScale_, GeV, 100.0*GeV, 10.0*GeV, 1000.0*GeV,
     false, false, Interface::limited);

  static Parameter<MEPP2GammaGammaPowheg,double> interfaceScaleFactor
    ("ScaleFactor",
     "The factor used before sHat if using a running scale",
     &MEPP2GammaGammaPowheg::scaleFact_, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<MEPP2GammaGammaPowheg,double> interfaceSamplingPower
    ("SamplingPower",
     "Power for the sampling of xp",
     &MEPP2GammaGammaPowheg::power_, 0.6, 0.0, 1.,
     false, false, Interface::limited);

}

void MEPP2GammaGammaPowheg::getDiagrams() const {
  // call base class
  MEPP2GammaGamma::getDiagrams();
}

Energy2 MEPP2GammaGammaPowheg::scale() const {
  return scaleopt_ == 0 ? sqr(fixedScale_) : scaleFact_*sHat();
}

int MEPP2GammaGammaPowheg::nDim() const {
  return MEPP2GammaGamma::nDim()+3;
}

bool MEPP2GammaGammaPowheg::generateKinematics(const double * r) {
  // generate the ztilde variable
  int ndim=nDim();
//   ztilde_ = 1. - pow(r[ndim-1],1./(1.-power_));
//   jac_.first  = pow(1.-ztilde_,power_)/(1.-power_);
  ztilde_ = r[ndim-1];
  // vtilde
  vtilde_ = r[ndim-2];
//   jac_.second = 1.;
  // generate the azimuth
  phi_ = r[ndim-3]*Constants::twopi;
  // call base class and return
  return MEPP2GammaGamma::generateKinematics(r);
}

double MEPP2GammaGammaPowheg::me2() const {
  return MEPP2GammaGamma::me2();
}

CrossSection MEPP2GammaGammaPowheg::dSigHatDR() const {
  return MEPP2GammaGamma::dSigHatDR()*NLOweight();
}

Selector<MEBase::DiagramIndex>
MEPP2GammaGammaPowheg::diagrams(const DiagramVector & diags) const {
  return MEPP2GammaGamma::diagrams(diags);
}

Selector<const ColourLines *>
MEPP2GammaGammaPowheg::colourGeometries(tcDiagPtr diag) const {
  return MEPP2GammaGamma::colourGeometries(diag);
}

double MEPP2GammaGammaPowheg::NLOweight() const {
  using Constants::pi;
  // If only leading order is required return 1:
  if(contrib_==0) return 1.;
  // if g g -> gamma gamma no correction
  if(mePartonData()[0]->id()==ParticleID::g&&
     mePartonData()[1]->id()==ParticleID::g) return 1.;
  // scales and prefactors
  Energy2 mu2 = scale();
  alphaS_ = SM().alphaS(mu2);
  double CFfact = 4./3.*alphaS_/Constants::twopi;
  double TRfact = 1./2.*alphaS_/Constants::twopi; 
  // first the LO + virtual after subtraction
  double v=1.+tHat()/sHat();
  double lokin = (1.-v)/v+v/(1.-v);
  double virt = 1.
    +CFfact*(3.+sqr(log(v))+sqr(log(1.-v))+3.*log(1.-v)
	     +1./lokin*(2.*log(v)+2.*log(1.-v)
			+3.*(1.-v)/v*log(v/(1.-v))
			+(2.+v/(1.-v))*sqr(log(v))
			+(2.+(1.-v)/v)*sqr(log(1.-v))));
  // extract the partons and stuff for the real emission and
  // collinear counter terms
  pair<tcBeamPtr,tcBeamPtr> hadrons= 
    make_pair(dynamic_ptr_cast<tcBeamPtr>(lastParticles().first->dataPtr()),
	      dynamic_ptr_cast<tcBeamPtr>(lastParticles().second->dataPtr()));
  pair<double,double> x=make_pair(lastX1(),lastX2());
  pair<tcPDPtr,tcPDPtr> partons = make_pair(mePartonData()[0],mePartonData()[1]);
  // If necessary swap the particle data objectss so that 
  // first beam gives the incoming quark
  if(lastPartons().first ->dataPtr()!=partons.first) {
    swap(x.first,x.second);
    swap(hadrons.first,hadrons.second);
  }
  // z for different terms
  pair<double,double> z;
  pair<double,double> zJac;
//   z.first  = x.first  + (1.-x.first )*ztilde_;
//   zJac.first  = (1.-x.first );
//   z.second = x.second + (1.-x.second)*ztilde_; 
//   zJac.second = (1.-x.second);
  double rhomax(pow(1.-x.first,1.-power_));
  double rho = ztilde_*rhomax;
  z.first = 1.-pow(rho,1./(1.-power_));
  zJac.first = rhomax*pow(1.-z.first,power_)/(1.-power_);
  rhomax = pow(1.-x.second,1.-power_);
  rho = ztilde_*rhomax; 
  z.second = 1.-pow(rho,1./(1.-power_));
  zJac.second = rhomax*pow(1.-z.second,power_)/(1.-power_);
  // PDFS
  // LO for q/qbar
  pair<double,double> oldqPDF = 
    make_pair(hadrons.first ->pdf()->xfx(hadrons.first ,partons.first ,scale(),
					 x.first )/x.first ,
	      hadrons.second->pdf()->xfx(hadrons.second,partons.second,scale(),
					 x.second)/x.second);
  // real/coll q/qbar
  pair<double,double> newqPDF = 
    make_pair(hadrons.first ->pdf()->xfx(hadrons.first ,partons.first ,scale(),
					 x.first /z.first )*z.first /x.first ,
	      hadrons.second->pdf()->xfx(hadrons.second,partons.second,scale(),
					 x.second/z.second)*z.second/x.second);
  // real/coll gluon
  pair<double,double> newgPDF =  
    make_pair(hadrons.first ->pdf()->xfx(hadrons.first ,gluon_,scale(),
					 x.first /z.first )*z.first /x.first ,
	      hadrons.second->pdf()->xfx(hadrons.second,gluon_,scale(),
					 x.second/z.second)*z.second/x.second);
  // coll terms
  // g -> q
  double collGQ    = zJac.first /z.first *newgPDF.first /oldqPDF.first *TRfact*
    ((sqr(z.first )+sqr(1.-z.first ))*log(sqr(1.-z.first )*sHat()/z.first /mu2)
     +2.*z.first *(1.-z.first ));
  // g -> qbar
  double collGQbar = zJac.second/z.second*newgPDF.second/oldqPDF.second*TRfact*
    ((sqr(z.second)+sqr(1.-z.second))*log(sqr(1.-z.second)*sHat()/z.second/mu2)
     +2.*z.second*(1.-z.second));
  // q -> q
  double collQQ =
    CFfact*(
	    // this bit is multiplied by LO PDF
	    sqr(pi)/3.-5.+2.*sqr(log(1.-x.first ))
	    +(1.5+2.*log(1.-x.first ))*log(sHat()/mu2)
	    // NLO PDF bit
	    +zJac.first /z.first * newqPDF.first /oldqPDF.first *
	    (1.-z.first -(1.+z.first )*log(sqr(1.-z.first )/z.first )
	     -(1.+z.first )*log(sHat()/mu2)-2.*log(z.first )/(1.-z.first ))
	    // + function bit
	    +zJac.first /z.first *(newqPDF.first /oldqPDF.first -z.first )*
	    2./(1.-z.first )*(log(sHat()*sqr(1.-z.first )/mu2)));
  // qbar -> qbar
  double collQbarQbar =
    CFfact*(
	    // this bit is multiplied by LO PDF
	    sqr(pi)/3.-5.+2.*sqr(log(1.-x.second))
	    +(1.5+2.*log(1.-x.second))*log(sHat()/mu2)
	    // NLO PDF bit
	    +zJac.second/z.second*newqPDF.second/oldqPDF.second*
	    (1.-z.second-(1.+z.second)*log(sqr(1.-z.second)/z.second)
	     -(1.+z.second)*log(sHat()/mu2)-2.*log(z.second)/(1.-z.second))
	    // + function bit
	    +zJac.second/z.second*(newqPDF.second/oldqPDF.second-z.second)*
	    2./(1.-z.second)*(log(sHat()*sqr(1.-z.second)/mu2)));
  // sum
//   double coll = collQQ+collQbarQbar+collGQ+collGQbar;
  double coll = collQQ+collQbarQbar;
  // real emission stuff
  // the vtildes
  pair<double,double> vt = make_pair(vtilde_*(1.-z.first),
				     vtilde_*(1.-z.second));
  pair<double,double> vJac = make_pair(1.-z.first ,
				       1.-z.second);
  // pTs
  pair<Energy,Energy> pT = 
    make_pair(sqrt(sHat()*vt.first *(1.-vt.first -z.first )/z.first ),
	      sqrt(sHat()*vt.second*(1.-vt.second-z.second)/z.second));
  // rapidities
  pair<double,double> rapidity = 
    make_pair(-log(x.second*sqrt(lastS())/pT.first *vt.first ),
	       log(x.first *sqrt(lastS())/pT.second*vt.second));
  // CMS system
  Energy rs=sqrt(lastS());
  Lorentz5Momentum pcmf = Lorentz5Momentum(ZERO,ZERO,0.5*rs*(x.first-x.second),
					   0.5*rs*(x.first+x.second));
  pcmf.rescaleMass();
  Boost blab(pcmf.boostVector());
  // emission from the quark radiation
  vector<Lorentz5Momentum> pfirst(5);
  pfirst [0] = Lorentz5Momentum(ZERO,ZERO,0.5*rs*x.first/z.first,
				0.5*rs*x.first/z.first,ZERO);
  pfirst [1] = Lorentz5Momentum(ZERO,ZERO,-0.5*rs*x.second,
				0.5*rs*x.second,ZERO) ;
  pfirst [2] = meMomenta()[2];
  pfirst [3] = meMomenta()[3];
  pfirst [4] = Lorentz5Momentum(pT.first*cos(phi_),pT.first*sin(phi_),
				pT.first*sinh(rapidity.first),
				pT.first*cosh(rapidity.first), ZERO);
  Lorentz5Momentum K  = pfirst [0]+pfirst [1]-pfirst [4];
  Lorentz5Momentum Kt = pcmf;
  Lorentz5Momentum Ksum = K+Kt;
  Energy2 K2 = K.m2();
  Energy2 Ksum2 = Ksum.m2();
  for(unsigned int ix=2;ix<4;++ix) {
    pfirst [ix].boost(blab);
    pfirst [ix] = pfirst [ix] - 2.*Ksum*(Ksum*pfirst [ix])/Ksum2
      +2*K*(Kt*pfirst [ix])/K2;
  }
  // emission from the antiquark
  vector<Lorentz5Momentum> psecond(5);
  psecond[0] = Lorentz5Momentum(ZERO,ZERO,0.5*rs*x.first,
				0.5*rs*x.first,ZERO);
  psecond[1] = Lorentz5Momentum(ZERO,ZERO,-0.5*rs*x.second/z.second,
				0.5*rs*x.second/z.second,ZERO) ;
  psecond[2] = meMomenta()[2];
  psecond[3] = meMomenta()[3];
  psecond[4] = Lorentz5Momentum(pT.second*cos(phi_),pT.second*sin(phi_),
				pT.second*sinh(rapidity.second),
				pT.second*cosh(rapidity.second), ZERO);
  K  = psecond[0]+psecond[1]-psecond[4];
  Ksum = K+Kt;
  K2 = K.m2();
  Ksum2 = Ksum.m2();
  for(unsigned int ix=2;ix<4;++ix) {
    psecond[ix].boost(blab);
    psecond[ix] = psecond[ix] - 2.*Ksum*(Ksum*psecond[ix])/Ksum2
      +2*K*(Kt*psecond[ix])/K2;
  }
  // phase-space prefactors
  pair<double,double> phase = 
    make_pair(zJac.first *vJac.first /sqr(z.first ),
 	      zJac.second*vJac.second/sqr(z.second));
  // real emission q qbar -> gamma gamma g
  // emission from q
  double realQQ       = CFfact*phase.first *MEqqbarg(pfirst ,true )*
    newqPDF.first /oldqPDF.first ;
  // emission from qbar
  double realQbarQbar = CFfact*phase.second*MEqqbarg(psecond,false)*
    newqPDF.second/oldqPDF.second;
  // real emission q g -> gamma gamma q
  double realGQ       = TRfact*phase.second*MEqgq(psecond)*
    newgPDF.second/oldqPDF.second;
  // real emission g qbar -> gamma gamma qbar
  double realGQbar    = TRfact*phase.first *MEqbargqbar(pfirst)*
    newgPDF.first /oldqPDF.first ;
  // return the answer
  double real = realQQ+realQbarQbar+realGQ+realGQbar;
  // return the answer
  double wgt = virt+coll+real;
//   if(wgt<0) {
//   cerr << "testing wgt   " << wgt << "\n"; 
//     cerr << "testing x1,x2 " << x.first << " " << x.second << "\n";
//     cerr << "testing z1,z2 " << z.first << " " << z.second << "\n";
//     cerr << "testing pieces " << virt << " " << coll << " " << real << "\n";
//     cerr << "testing coll A " << collQQ << " " << collQbarQbar << "\n";
//     cerr << "testing coll B " 
// 	 << CFfact*(sqr(pi)/3.-5.+2.*sqr(log(1.-x.first ))
// 		    +(1.5+2.*log(1.-x.first ))*log(sHat()/mu2)) << " "
// 	 << CFfact*(+zJac.first /z.first * newqPDF.first /oldqPDF.first *
// 		    (1.-z.first -(1.+z.first )*log(sqr(1.-z.first )/z.first )
// 		     -(1.+z.first )*log(sHat()/mu2)-2.*log(z.first )/(1.-z.first ))) << " "
// 	 << CFfact*(+zJac.first /z.first *(newqPDF.first /oldqPDF.first -z.first )*
// 		    2./(1.-z.first )*(log(sHat()*sqr(1.-z.first )/mu2)))
// 	 << "\n";
//     cerr << "testing coll C " 
// 	 << CFfact*(sqr(pi)/3.-5.+2.*sqr(log(1.-x.second))
// 		    +(1.5+2.*log(1.-x.second))*log(sHat()/mu2)) << " "
// 	 << CFfact*(+zJac.second/z.second*newqPDF.second/oldqPDF.second*
// 		    (1.-z.second-(1.+z.second)*log(sqr(1.-z.second)/z.second)
// 		     -(1.+z.second)*log(sHat()/mu2)-2.*log(z.second)/(1.-z.second))) 
// 	 << " "
// 	 << CFfact*(+zJac.second/z.second*(newqPDF.second/oldqPDF.second-z.second)*
// 		    2./(1.-z.second)*(log(sHat()*sqr(1.-z.second)/mu2))) << "\n";
//   }
  return contrib_==1 ? max(0.,wgt) : max(0.,-wgt);
}

void MEPP2GammaGammaPowheg::doinit() {
  MEPP2GammaGamma::doinit();
  gluon_ = getParticleData(ParticleID::g);
  // get the vedrtex pointers from the SM object
  tcHwSMPtr hwsm= dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  // do the initialisation
  if(!hwsm) throw InitException() << "Wrong type of StandardModel object in "
				  << "MEPP2GammaGammaPowheg::doinit() the Herwig++"
				  << " version must be used" 
				  << Exception::runerror;
  QEDVertex_ = hwsm->vertexFFP();
  QCDVertex_ = hwsm->vertexFFG();
}

double MEPP2GammaGammaPowheg::MEqqbarg(const vector<Lorentz5Momentum> & p,
				       bool order) const {
  using namespace ThePEG::Helicity;
  double sum(0.);
  vector<SpinorWaveFunction> qin;
  vector<SpinorBarWaveFunction> qbarin;
  vector<VectorWaveFunction> pout1,pout2,gout;
  SpinorWaveFunction    q_in   (p[0],mePartonData()[0],incoming);
  SpinorBarWaveFunction qbar_in(p[1],mePartonData()[1],incoming);
  VectorWaveFunction     p_out1(p[2],mePartonData()[2],outgoing);
  VectorWaveFunction     p_out2(p[3],mePartonData()[3],outgoing);
  VectorWaveFunction     g_out (p[4],gluon_           ,outgoing);
  for(unsigned int ix=0;ix<2;++ix) {
    q_in.reset(ix);               qin.push_back(q_in   );
    qbar_in.reset(ix);         qbarin.push_back(qbar_in);
    g_out.reset(2*ix);           gout.push_back(g_out  );
    p_out1.reset(2*ix);         pout1.push_back(p_out1 );
    p_out2.reset(2*ix);         pout2.push_back(p_out2 );
  }
  vector<Complex> diag(6);
  for(unsigned int ihel1=0;ihel1<2;++ihel1) {
    for(unsigned int ihel2=0;ihel2<2;++ihel2) {
      for(unsigned int phel1=0;phel1<2;++phel1) {
	for(unsigned int phel2=0;phel2<2;++phel2) {
	  for(unsigned int ghel=0;ghel<2;++ghel) {
 	    // first diagram
	    SpinorWaveFunction inters1 = 
	      QEDVertex_->evaluate(ZERO,5,mePartonData()[0],
				   qin[ihel1],pout1[phel1]);
	    SpinorBarWaveFunction inters2 = 
	      QEDVertex_->evaluate(ZERO,5,mePartonData()[1],
				   qbarin[ihel2],pout2[phel2]);
	    diag[0] = QCDVertex_->evaluate(scale(),inters1,inters2,gout[ghel]);
	    // second diagram
	    SpinorWaveFunction inters3 = 
	      QCDVertex_->evaluate(scale(),5,mePartonData()[0],qin[ihel1],gout[ghel]);
	    SpinorBarWaveFunction inters4 = 
	      QEDVertex_->evaluate(ZERO,5,mePartonData()[1],qbarin[ihel2],pout1[phel1]);
	    diag[1] = QEDVertex_->evaluate(ZERO,inters3,inters4,pout2[phel2]);
	    // fourth diagram
	    diag[2] = QEDVertex_->evaluate(ZERO,inters3,inters2,pout1[phel1]);
	    // fifth diagram
	    SpinorBarWaveFunction inters5 = 
	      QCDVertex_->evaluate(scale(),5,mePartonData()[1],qbarin[ihel2],gout[ghel]);
	    diag[3] = 
	      QEDVertex_->evaluate(ZERO,inters1,inters5,pout2[phel2]);
	    // sixth diagram
	    SpinorWaveFunction inters6 = 
	      QEDVertex_->evaluate(ZERO,5,mePartonData()[0],qin[ihel1],pout2[phel2]);
	    diag[4] = QCDVertex_->evaluate(scale(),inters6,inters4,gout[ghel]);
	    // eighth diagram
	    diag[5] = QEDVertex_->evaluate(ZERO,inters6,inters5,pout1[phel1]);
	    // sum
	    Complex dsum = std::accumulate(diag.begin(),diag.end(),Complex(0.));
	    sum += norm(dsum);
	  }
	}
      }
    }
  }
  // remove some coupling factors
  sum /= (8.*Constants::pi*alphaS_);
  // compute the two dipole terms
  double x = (p[0]*p[1]-p[4]*p[1]-p[4]*p[0])/(p[0]*p[1]);
  Lorentz5Momentum Kt = p[0]+p[1]-p[4];
  Lorentz5Momentum pa[4],pb[4];
  // momenta for q -> q g emission
  pa[0] = x*p[0];
  pa[1] = p[1];
  Lorentz5Momentum K = pa[0]+pa[1];
  Lorentz5Momentum Ksum = K+Kt;
  Energy2 K2 = K.m2();
  Energy2 Ksum2 = Ksum.m2();
  for(unsigned int ix=2;ix<4;++ix)
    pa[ix] = p[ix]-2.*Ksum*(Ksum*p[ix])/Ksum2+2*K*(Kt*p[ix])/K2;
  // momenta
  pb[0] = p[0];
  pb[1] = x*p[1];
  K = pb[0]+pb[1];
  Ksum = K+Kt;
  K2 = K.m2();
  Ksum2 = Ksum.m2();
  for(unsigned int ix=2;ix<4;++ix)
    pb[ix] = p[ix]-2.*Ksum*(Ksum*p[ix])/Ksum2+2*K*(Kt*p[ix])/K2;
  // first LO matrix element
  Energy2 s,t,u;
  s = (pa[0]+pa[1]).m2();
  t = (pa[0]-pa[2]).m2();
  u = (pa[0]-pa[3]).m2();
  double lo1 = 8.*sqr(4.*Constants::pi*SM().alphaEM(ZERO))*(t/u+u/t)*
    pow(double(mePartonData()[0]->iCharge())/3.,4);
  // second LO matrix element
  s = (pb[0]+pb[1]).m2();
  t = (pb[0]-pb[2]).m2();
  u = (pb[0]-pb[3]).m2();
  double lo2 = 8.*sqr(4.*Constants::pi*SM().alphaEM(ZERO))*(t/u+u/t)*
    pow(double(mePartonData()[0]->iCharge())/3.,4);
  // first dipole
  InvEnergy2 D1 = 0.5/(p[0]*p[4])/x*(2./(1.-x)-(1.+x));
  // second dipole
  InvEnergy2 D2 = 0.5/(p[1]*p[4])/x*(2./(1.-x)-(1.+x));
  // result
  double me;
  if(order) {
    me = sHat()*(abs(D1)/(abs(D1)*lo1+abs(D2)*lo2)*UnitRemoval::InvE2*sum-D1);
  }
  else {
    me = sHat()*(abs(D2)/(abs(D1)*lo1+abs(D2)*lo2)*UnitRemoval::InvE2*sum-D2);
  }
  return me;
}

double MEPP2GammaGammaPowheg::MEqgq(const vector<Lorentz5Momentum> & p) const {
  using namespace ThePEG::Helicity;
  double sum(0.);
  vector<SpinorWaveFunction> qin;
  vector<SpinorBarWaveFunction> qout;
  vector<VectorWaveFunction> pout1,pout2,gin;
  SpinorWaveFunction q_in    (p[0],mePartonData()[0]      ,incoming);
  VectorWaveFunction g_in    (p[1],gluon_                 ,incoming);
  VectorWaveFunction p_out1  (p[2],mePartonData()[2]      ,outgoing);
  VectorWaveFunction p_out2  (p[3],mePartonData()[3]      ,outgoing);
  SpinorBarWaveFunction q_out(p[4],mePartonData()[1]->CC(),outgoing);
  for(unsigned int ix=0;ix<2;++ix) {
    q_in .reset(ix);         qin  .push_back(q_in  );
    q_out.reset(ix);         qout .push_back(q_out );
    g_in .reset(2*ix);       gin  .push_back(g_in  );
    p_out1.reset(2*ix);      pout1.push_back(p_out1);
    p_out2.reset(2*ix);      pout2.push_back(p_out2);
  }
  vector<Complex> diag(6);
  for(unsigned int ihel1=0;ihel1<2;++ihel1) {
    for(unsigned int ihel2=0;ihel2<2;++ihel2) {
      for(unsigned int phel1=0;phel1<2;++phel1) {
	for(unsigned int phel2=0;phel2<2;++phel2) {
	  for(unsigned int ohel=0;ohel<2;++ohel) {
 	    // first diagram
	    SpinorWaveFunction inters1 = 
	      QEDVertex_->evaluate(ZERO,5,qin[ihel1].getParticle()->CC(),
				   qin[ihel1],pout1[phel1]);
	    SpinorBarWaveFunction inters2 = 
	      QEDVertex_->evaluate(ZERO,5,qout[ohel].getParticle(),
				   qout[ohel],pout2[phel2]);
	    diag[0] = QCDVertex_->evaluate(scale(),inters1,inters2,gin[ihel2]);
	    // second diagram
	    SpinorWaveFunction inters3 = 
	      QCDVertex_->evaluate(scale(),5,qin[ihel1].getParticle()->CC(),
				   qin[ihel1],gin[ihel2]);
	    SpinorBarWaveFunction inters4 = 
	      QEDVertex_->evaluate(ZERO,5,qout[ohel].getParticle(),
				   qout[ohel],pout1[phel1]);
	    diag[1] = QEDVertex_->evaluate(ZERO,inters3,inters4,pout2[phel2]);
	    // fourth diagram
	    diag[2] = QEDVertex_->evaluate(ZERO,inters3,inters2,pout1[phel1]);
	    // fifth diagram
	    SpinorBarWaveFunction inters5 = 
	      QCDVertex_->evaluate(scale(),5,qout[ohel].getParticle(),
				   qout[ohel],gin[ihel2]);
	    diag[3] = 
	      QEDVertex_->evaluate(ZERO,inters1,inters5,pout2[phel2]);
	    // sixth diagram
	    SpinorWaveFunction inters6 = 
	      QEDVertex_->evaluate(ZERO,5,qin[ihel1].getParticle()->CC(),
				   qin[ihel1],pout2[phel2]);
	    diag[4] = QCDVertex_->evaluate(scale(),inters6,inters4,gin[ihel2]);
	    // eighth diagram
	    diag[5] = QEDVertex_->evaluate(ZERO,inters6,inters5,pout1[phel1]);
	    // sum
	    Complex dsum = std::accumulate(diag.begin(),diag.end(),Complex(0.));
	    sum += norm(dsum);
	  }
	}
      }
    }
  }
  // remove some coupling factors
  sum /= (8.*Constants::pi*alphaS_);
  // compute the two dipole terms
  double x = 1.-(p[4]*p[1]+p[4]*p[0])/(p[0]*p[1]);
  Lorentz5Momentum Kt = p[0]+p[1]-p[4];
  Lorentz5Momentum pa[4],pb[4],pc[4];
  // momenta for ISR
  pa[0] = p[0];
  pa[1] = x*p[1];
  Lorentz5Momentum K = pb[0]+pb[1];
  Lorentz5Momentum Ksum = K+Kt;
  Energy2 K2 = K.m2();
  Energy2 Ksum2 = Ksum.m2();
  for(unsigned int ix=2;ix<4;++ix)
    pa[ix] = p[ix]-2.*Ksum*(Ksum*p[ix])/Ksum2+2*K*(Kt*p[ix])/K2;
  // first LO matrix element
  Energy2 s,t,u;
  s = (pa[0]+pa[1]).m2();
  t = (pa[0]-pa[2]).m2();
  u = (pa[0]-pa[3]).m2();
  double coupling = 
    sqr(4.*Constants::pi*generator()->standardModel()->alphaEM(ZERO))*
    pow(double(mePartonData()[0]->iCharge())/3.,4);
  double lo1 = 8.*coupling*(t/u+u/t);
  // momenta for first FS emission
  double xFS[2],zFS[2];
  for(unsigned int ix=0;ix<2;++ix) {
    xFS[ix] = 1.-(p[ix+2]*p[4])/((p[ix+2]+p[4])*p[0]);
    zFS[ix] = (p[0]*p[4])/((p[ix+2]+p[4])*p[0]);
  }
  // first set of momenta
  pb[0] = xFS[0]*p[0];
  pb[1] = p[1];
  pb[2] = p[3];
  pb[3] = p[2]+p[4]-(1.-xFS[0])*p[0];
  // second set of momenta
  pc[0] = xFS[1]*p[0];
  pc[1] = p[1];
  pc[2] = p[2];
  pc[3] = p[3]+p[4]-(1.-xFS[1])*p[0];
  s = (pb[0]+pb[1]).m2();
  t = (pb[0]-pb[2]).m2();
  u = (pb[0]-pb[3]).m2();
  double lo2 = -8./s/t*(s*s+t*t+2.*u*(s+t+u))*coupling;
  // third  LO matrix element
  s = (pc[0]+pc[1]).m2();
  t = (pc[0]-pc[2]).m2();
  u = (pc[0]-pc[3]).m2();
  double lo3 = -8./s/t*(s*s+t*t+2.*u*(s+t+u))*coupling;
  // first dipole
  InvEnergy2 D1 =  0.5/(p[1]*p[4])/x*(1.-2.*x*(1.-x));
  // second dipole
  InvEnergy2 D2 =  0.5/(p[4]*p[2])*(2./(2.-zFS[0]-xFS[0])-1.-zFS[0]);
  // third  dipole
  InvEnergy2 D3 =  0.5/(p[4]*p[3])*(2./(2.-zFS[1]-xFS[1])-1.-zFS[1]);
  return sHat()*(abs(D1)/(abs(D1)*lo1+abs(D2)*lo2+abs(D3)*lo3)*UnitRemoval::InvE2*sum
		 -D1);
}

double MEPP2GammaGammaPowheg::MEqbargqbar(const vector<Lorentz5Momentum> & p) const {
  return 0.;
}
