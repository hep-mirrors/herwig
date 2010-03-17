// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NLODrellYanBase class.
//

#include "NLODrellYanBase.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/PDF/BeamParticleData.h"
#include <numeric>

using namespace Herwig;

/**
 *  Typedef for BeamParticleData
 */
typedef Ptr<BeamParticleData>::transient_const_pointer tcBeamPtr;

NLODrellYanBase::NLODrellYanBase() : contrib_(1), power_(0.1),
				     alphaS_(0.), fixedAlphaS_(false),
				     supressionFunction_(0), 
				     lambda2_(10000.*GeV2)
{}

Energy2 NLODrellYanBase::scale() const {
  //return sHat();
  return sqr(0.5*(mePartonData()[2]->mass()+mePartonData()[3]->mass()));
}

int NLODrellYanBase::nDim() const {
  return HwMEBase::nDim() + ( contrib_>=1 && contrib_<=3 ? 3 : 0 );
}

bool NLODrellYanBase::generateKinematics(const double * r) {
  if(contrib_>=1&&contrib_<=3) {
    zTilde_ = r[nDim()-1];
    vTilde_ = r[nDim()-2];
    phi_    = Constants::twopi*r[nDim()-3];
  }
  jacobian(1.0);
  return HwMEBase::generateKinematics(r);
}

CrossSection NLODrellYanBase::dSigHatDR() const {
  // old technique
  CrossSection preFactor = 
    jacobian()/(16.0*sqr(Constants::pi)*sHat())*sqr(hbarc);
  loME_ = me2();
  if(contrib_<4) return NLOWeight()*preFactor;
  // folding technique to ensure positive
  double wgt(0.);
  unsigned int ntry(0);
  do {
    // radiative variables
    zTilde_ = UseRandom::rnd();
    vTilde_ = UseRandom::rnd();
    phi_    = Constants::twopi*UseRandom::rnd();
    wgt += NLOWeight();
    ++ntry;
  }
  while (wgt<0.&&ntry<100);
  if(wgt<0.) return ZERO;
  return wgt*preFactor/double(ntry);
}

double NLODrellYanBase::me2() const {
  return loME(mePartonData(),meMomenta(),true);
}

void NLODrellYanBase::persistentOutput(PersistentOStream & os) const {
  os << contrib_ << power_ << gluon_ << fixedAlphaS_ << alphaS_
     << supressionFunction_ << ounit(lambda2_,GeV2);
}

void NLODrellYanBase::persistentInput(PersistentIStream & is, int) {
  is >> contrib_ >> power_ >> gluon_ >> fixedAlphaS_ >> alphaS_
     >> supressionFunction_ >> iunit(lambda2_,GeV2);
}

void NLODrellYanBase::doinit() {
  HwMEBase::doinit();
  gluon_ = getParticleData(ParticleID::g);
}

AbstractClassDescription<NLODrellYanBase> NLODrellYanBase::initNLODrellYanBase;
// Definition of the static class description member.

void NLODrellYanBase::Init() {

  static ClassDocumentation<NLODrellYanBase> documentation
    ("The NLODrellYanBase class provides a base class for the"
     " implementation of the next-to-leading order corrections"
     " to Drell Yan type processes in the POWHEG approach");

   static Switch<NLODrellYanBase,unsigned int> interfaceContribution
    ("Contribution",
     "Which contributions to the cross section to include",
     &NLODrellYanBase::contrib_, 1, false, false);
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
  static SwitchOption interfaceContributionReal
    (interfaceContribution,
     "Real",
     "Generate the high pT real emission only",
     3);
  static SwitchOption interfaceContributionTotalNLO
    (interfaceContribution,
     "TotalNLO",
     "Generate the full NLO cross section using folding",
     4);

  static Parameter<NLODrellYanBase,double> interfaceSamplingPower
    ("SamplingPower",
     "Power for the sampling of xp",
     &NLODrellYanBase::power_, 0.6, 0.0, 1.,
     false, false, Interface::limited);

  static Switch<NLODrellYanBase,bool> interfaceFixedAlphaS
    ("FixedAlphaS",
     "Use a fixed value of alpha_S",
     &NLODrellYanBase::fixedAlphaS_, false, false, false);
  static SwitchOption interfaceFixedAlphaSYes
    (interfaceFixedAlphaS,
     "Yes",
     "Use fixed alpha_S",
     true);
  static SwitchOption interfaceFixedAlphaSNo
    (interfaceFixedAlphaS,
     "No",
     "Use running alpha_S",
     false);

  static Parameter<NLODrellYanBase,double> interfaceAlphaS
    ("AlphaS",
     "The fixed value of alpha_S to use",
     &NLODrellYanBase::alphaS_, 0., 0., 1.,
     false, false, Interface::limited);


  static Switch<NLODrellYanBase,unsigned int> interfaceSupressionFunction
    ("SupressionFunction",
     "Choice of the supression function",
     &NLODrellYanBase::supressionFunction_, 0, false, false);
  static SwitchOption interfaceSupressionFunctionNone
    (interfaceSupressionFunction,
     "None",
     "Default POWHEG approach",
     0);
  static SwitchOption interfaceSupressionFunctionThetaFunction
    (interfaceSupressionFunction,
     "ThetaFunction",
     "Use theta functions at scale Lambda",
     1);
  static SwitchOption interfaceSupressionFunctionSmooth
    (interfaceSupressionFunction,
     "Smooth",
     "Supress high pT by pt^2/(pt^2+lambda^2)",
     2);

  static Parameter<NLODrellYanBase,Energy2> interfaceSupressionScale
    ("SupressionScale",
     "The square of the scale for the supression function",
     &NLODrellYanBase::lambda2_, GeV2, 10000.0*GeV2, 0.0*GeV2, 0*GeV2,
     false, false, Interface::lowerlim);

}

double NLODrellYanBase::subtractedVirtual() const {
  Singular virt = virtualME();
  Singular Ieps;
  Ieps.eps2 = 2;
  Ieps.eps1 = 3;
  Ieps.finite = 10.-sqr(Constants::pi);
  // check the singular pieces cancel
  assert(Ieps.eps2==-virt.eps2 && Ieps.eps1 == -virt.eps1 );
  return virt.finite+Ieps.finite;
}

double NLODrellYanBase::NLOWeight() const {
  // if leading-order return 1
  if(contrib_==0) return loME_;
  // strong coupling
  Energy2 mu2(scale());
  if(!fixedAlphaS_) alphaS_ = SM().alphaS(mu2);
  // prefactors
  CFfact_ = 4./3.*alphaS_/Constants::twopi;
  TRfact_ = 1./2.*alphaS_/Constants::twopi;
  // virtual pieces
  double virt = CFfact_*subtractedVirtual();
  // extract the partons and stuff for the real emission
  //  and collinear counter terms
  // hadrons
  pair<tcBeamPtr,tcBeamPtr> hadrons= 
    make_pair(dynamic_ptr_cast<tcBeamPtr>(lastParticles().first->dataPtr() ),
	      dynamic_ptr_cast<tcBeamPtr>(lastParticles().second->dataPtr()));
  // momentum fractions
  pair<double,double> x = make_pair(lastX1(),lastX2());
  // partons
  pair<tcPDPtr,tcPDPtr> partons = make_pair(mePartonData()[0],mePartonData()[1]);
  // If necessary swap the particle data objects so that 
  // first beam gives the incoming quark
  if(lastPartons().first ->dataPtr()!=partons.first) {
    swap(x.first,x.second);
    swap(hadrons.first,hadrons.second);
  }
  // convert the values of z tilde to z
  pair<double,double> z;
  pair<double,double> zJac;
  double rhomax(pow(1.-x.first,1.-power_));
  double rho = zTilde_*rhomax;
  z.first = 1.-pow(rho,1./(1.-power_));
  zJac.first = rhomax*pow(1.-z.first,power_)/(1.-power_);
  rhomax = pow(1.-x.second,1.-power_);
  rho = zTilde_*rhomax; 
  z.second = 1.-pow(rho,1./(1.-power_));
  zJac.second = rhomax*pow(1.-z.second,power_)/(1.-power_);
  // calculate the PDFs
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
  double collGQ    = collinearGluon(mu2,zJac.first,z.first,
				    oldqPDF.first,newgPDF.first);
  // g -> qbar
  double collGQbar = collinearGluon(mu2,zJac.second,z.second,
				    oldqPDF.second,newgPDF.second);
  // q -> q
  double collQQ       = collinearQuark(x.first ,mu2,zJac.first ,z.first ,
				       oldqPDF.first ,newqPDF.first );
  // qbar -> qbar
  double collQbarQbar = collinearQuark(x.second,mu2,zJac.second,z.second,
				       oldqPDF.second,newqPDF.second);
  // collinear remnants 
  double coll = collQQ+collQbarQbar+collGQ+collGQbar;
  vector<double> real1 = 
    subtractedReal(x,z. first,zJac. first,
		   oldqPDF. first,newqPDF. first,newgPDF. first, true);
  vector<double> real2 = 
    subtractedReal(x,z.second,zJac.second, 
		   oldqPDF.second,newqPDF.second,newgPDF.second,false);
  // add up all the terms and return the answer
  double wgt = loME_*( 1. + virt + coll ) +real1[0] + real1[2] +real2[0] +real2[2];
  if(isnan(wgt)||isinf(wgt)) {
    generator()->log() << "testing bad weight "
	 << collQQ << " " << collQbarQbar << " "
	 << collGQ << " " << collGQbar << " "
	 << virt << " " << coll << " "
	 << real1[0] << " " << real1[2] << " "
	 << real2[0] << " " << real2[2] << "\n";
    generator()->log() << "testing z " << z.first << " " << z.second << "\n";
    generator()->log() << "testing z " << 1.-z.first << " " << 1.-z.second << "\n";
    assert(false);
  }
  weights_.resize(5,0.);
  weights_[1] = real1[1];
  weights_[2] = real1[3];
  weights_[3] = real2[1];
  weights_[4] = real2[3];
  if( contrib_ < 3 ) {
    weights_[0] = wgt;
    wgt += real1[1]+real1[3]+real2[1]+real2[3];
    return contrib_ == 1 ? max(0.,wgt) : max(0.,-wgt);
  }
  else if(contrib_==3) {
    weights_[0] = 0.;
    return real1[1]+real1[3]+real2[1]+real2[3];
  }
  else 
    return wgt;
}

double NLODrellYanBase::collinearQuark(double x, Energy2 mu2, 
				       double jac, double z,
				       double oldPDF, double newPDF) const {
  if(1.-z < 1.e-8) return 0.;
  return CFfact_*(
		  // this bit is multiplied by LO PDF
		  sqr(Constants::pi)/3.-5.+2.*sqr(log(1.-x ))
		  +(1.5+2.*log(1.-x ))*log(sHat()/mu2)
		  // NLO PDF bit
		  +jac /z * newPDF /oldPDF *
		  (1.-z -(1.+z )*log(sqr(1.-z )/z )
		   -(1.+z )*log(sHat()/mu2)-2.*log(z )/(1.-z ))
		  // + function bit
		  +jac /z *(newPDF /oldPDF -z )*
		  2./(1.-z )*log(sHat()*sqr(1.-z )/mu2));
}

double NLODrellYanBase::collinearGluon(Energy2 mu2,
				       double jac, double z,
				       double oldPDF, double newPDF) const {
  if(1.-z < 1.e-8) return 0.;
  return TRfact_*jac/z*newPDF/oldPDF*
    ((sqr(z)+sqr(1.-z))*log(sqr(1.-z)*sHat()/z/mu2)
     +2.*z*(1.-z));
}

vector<double> NLODrellYanBase::subtractedReal(pair<double,double> x, double z,
					       double zJac,
					       double oldqPDF,double newqPDF,
					       double newgPDF, bool order) const {
  double vt   = vTilde_*(1.-z);
  double vJac = 1.-z;
  Energy pT   = sqrt(sHat()*vt*(1.-vt-z)/z);
  // rapidities
  double rapidity;
  if(order) {
    rapidity = -log(x.second*sqrt(lastS())/pT*vt);
  }
  else {
    rapidity =  log(x.first *sqrt(lastS())/pT*vt);
  }
  // CMS system
  Energy rs=sqrt(lastS());
  Lorentz5Momentum pcmf = Lorentz5Momentum(ZERO,ZERO,0.5*rs*(x.first-x.second),
					   0.5*rs*(x.first+x.second));
  pcmf.rescaleMass();
  Boost blab(pcmf.boostVector());
  // emission from the quark radiation
  vector<Lorentz5Momentum> pnew(5);
  if(order) {
    pnew [0] = Lorentz5Momentum(ZERO,ZERO,0.5*rs*x.first/z,
				0.5*rs*x.first/z,ZERO);
    pnew [1] = Lorentz5Momentum(ZERO,ZERO,-0.5*rs*x.second,
				0.5*rs*x.second,ZERO) ;
  }
  else {
    pnew[0] = Lorentz5Momentum(ZERO,ZERO,0.5*rs*x.first,
			       0.5*rs*x.first,ZERO);
    pnew[1] = Lorentz5Momentum(ZERO,ZERO,-0.5*rs*x.second/z,
			       0.5*rs*x.second/z,ZERO) ;
  }
  pnew [2] = meMomenta()[2];
  pnew [3] = meMomenta()[3];
  pnew [4] = Lorentz5Momentum(pT*cos(phi_),pT*sin(phi_),
			      pT*sinh(rapidity),
			      pT*cosh(rapidity), ZERO);
  Lorentz5Momentum K  = pnew [0]+pnew [1]-pnew [4];
  Lorentz5Momentum Kt = pcmf;
  Lorentz5Momentum Ksum = K+Kt;
  Energy2 K2 = K.m2();
  Energy2 Ksum2 = Ksum.m2();
  for(unsigned int ix=2;ix<4;++ix) {
    pnew [ix].boost(blab);
    pnew [ix] = pnew [ix] - 2.*Ksum*(Ksum*pnew [ix])/Ksum2
      +2*K*(Kt*pnew [ix])/K2;
  }
  // phase-space prefactors
  // double phase = zJac*vJac/sqr(z);
  double phase = zJac*vJac/z;
  // real emission q qbar
  vector<double> output(4,0.);
  if(order) {
    realEmissionGluon1_ = pnew;
    realEmissionQuark1_ = pnew;
  }
  else {
    realEmissionGluon2_ = pnew;
    realEmissionQuark2_ = pnew;
  }
  if(!(zTilde_<1e-7 || vt<1e-7 || 1.-z-vt < 1e-7 )) {
    pair<double,double> realQQ = subtractedMEqqbar(pnew,order);
    double fact1 = CFfact_*phase*newqPDF/oldqPDF;
    pair<double,double> realGQ = subtractedMEgqbar(pnew,order);
    double fact2 = TRfact_*phase*newgPDF/oldqPDF;
    output[0] = realQQ.first *fact1;
    output[1] = realQQ.second*fact1;
    output[2] = realGQ.first *fact2;
    output[3] = realGQ.second*fact2;
  }
  // return the answer
  return output;
}

pair<double,double> 
NLODrellYanBase::subtractedMEqqbar(const vector<Lorentz5Momentum> & p,
				   bool order) const {
  // use the inheriting class to calculate the matrix element
  cPDVector particles(mePartonData());
  particles.push_back(gluon_);
  double me = 0.75*realME(particles,p);
  // compute the two dipole terms
  double x = (p[0]*p[1]-p[4]*p[1]-p[4]*p[0])/(p[0]*p[1]);
  Lorentz5Momentum Kt = p[0]+p[1]-p[4];
  vector<Lorentz5Momentum> pa(4),pb(4);
  // momenta for q -> q g emission
  pa[0] = x*p[0];
  pa[1] =   p[1];
  Lorentz5Momentum K = pa[0]+pa[1];
  Lorentz5Momentum Ksum = K+Kt;
  Energy2 K2 = K.m2();
  Energy2 Ksum2 = Ksum.m2();
  for(unsigned int ix=2;ix<4;++ix)
    pa[ix] = p[ix]-2.*Ksum*(Ksum*p[ix])/Ksum2+2*K*(Kt*p[ix])/K2;
  // first LO matrix element
  double lo1 = loME(mePartonData(),pa,false);
  // momenta for qbar -> qbar g emission
  pb[0] =   p[0];
  pb[1] = x*p[1];
  K = pb[0]+pb[1];
  Ksum = K+Kt;
  K2 = K.m2();
  Ksum2 = Ksum.m2();
  for(unsigned int ix=2;ix<4;++ix)
    pb[ix] = p[ix]-2.*Ksum*(Ksum*p[ix])/Ksum2+2*K*(Kt*p[ix])/K2;
  // second LO matrix element
  double lo2 = loME(mePartonData(),pb,false);
  // first dipole
  InvEnergy2 D1 = 0.5/(p[0]*p[4])/x*(2./(1.-x)-(1.+x));
  // second dipole
  InvEnergy2 D2 = 0.5/(p[1]*p[4])/x*(2./(1.-x)-(1.+x));
  // results
  pair<double,double> supressionFactor = supressionFunction(sqr(p[4].x())+sqr(p[4].y()));
  pair<double,double> output = make_pair(0.,0.);
  if(lo1>0.&&lo2>0.) {
    if(order) {
      me *= abs(D1)*lo1/(abs(D1)*lo1+abs(D2)*lo2);
      output.first  = sHat()*(UnitRemoval::InvE2*me*supressionFactor.first -D1*lo1);
      output.second = sHat()*(UnitRemoval::InvE2*me*supressionFactor.second);
    }
    else {
      me *= abs(D2)*lo2/(abs(D1)*lo1+abs(D2)*lo2);
      output.first  = sHat()*(UnitRemoval::InvE2*me*supressionFactor.first -D2*lo2);
      output.second = sHat()*(UnitRemoval::InvE2*me*supressionFactor.second);
    }
  }
  return output;
}

pair<double,double>
NLODrellYanBase::subtractedMEgqbar(const vector<Lorentz5Momentum> & p,
				   bool order) const {
  // use the inheriting class to calculate the matrix element
  cPDVector particles(mePartonData());
  if(order) {
    particles.push_back(particles[0]->CC());
    particles[0] = gluon_;
  }
  else {
    particles.push_back(particles[1]->CC());
    particles[1] = gluon_;
  }
  double me = 2.*realME(particles,p);
  // compute the two dipole terms
  double x = 1.-(p[4]*p[1]+p[4]*p[0])/(p[0]*p[1]);
  Lorentz5Momentum Kt = p[0]+p[1]-p[4];
  vector<Lorentz5Momentum> pa(4);
  // momenta for ISR
  if(order) {
    pa[0] = x*p[0];
    pa[1] =   p[1];
  }
  else {
    pa[0] =   p[0];
    pa[1] = x*p[1];
  }
  Lorentz5Momentum K = pa[0]+pa[1];
  Lorentz5Momentum Ksum = K+Kt;
  Energy2 K2 = K.m2();
  Energy2 Ksum2 = Ksum.m2();
  for(unsigned int ix=2;ix<4;++ix)
    pa[ix] = p[ix]-2.*Ksum*(Ksum*p[ix])/Ksum2+2*K*(Kt*p[ix])/K2; 
  // first LO matrix element 
  double lo1 = loME(mePartonData(),pa,false); 
  // dipole
  InvEnergy2 D1;
  if(order) {
    D1 =  0.5/(p[0]*p[4])/x*(1.-2.*x*(1.-x));
  }
  else {
    D1 =  0.5/(p[1]*p[4])/x*(1.-2.*x*(1.-x));
  }
  pair<double,double> supressionFactor = 
    supressionFunction(sqr(p[4].x())+sqr(p[4].y()));
  return make_pair(sHat()*(UnitRemoval::InvE2*me*supressionFactor.first -D1*lo1),
		   sHat()*(UnitRemoval::InvE2*me*supressionFactor.second));
}

HardTreePtr NLODrellYanBase::generateHardest(ShowerTreePtr tree) {
  // get the particles to be showered
  _beams.clear();
  _partons.clear();
  // find the incoming particles
  ShowerParticleVector incoming;
  map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator cit;
  _quarkplus = true;
  vector<ShowerProgenitorPtr> particlesToShower;
  //progenitor particles are produced in z direction.
  for( cit = tree->incomingLines().begin(); cit != tree->incomingLines().end(); ++cit ) {
    incoming.push_back( cit->first->progenitor() );
    _beams.push_back( cit->first->beam() );
    _partons.push_back( cit->first->progenitor()->dataPtr() );
    // check that quark is along +ve z direction
    if(cit->first->progenitor()->id() > 0 &&
       cit->first->progenitor()->momentum().z() < ZERO ) 
      _quarkplus = false;
    particlesToShower.push_back( cit->first );
  } 
  // we are assuming quark first, swap order to ensure this
  // if antiquark first
  if(_partons[0]->id()<_partons[1]->id()) {
    swap(_partons[0],_partons[1]);
    swap(_beams[0],_beams[1]);
  }
  double wgtb = std::accumulate(++weights_.begin(),weights_.end(),0.);
  int emission_type(-1);
  if(wgtb>UseRandom::rnd()*(weights_[0]+wgtb)) {
    wgtb *= UseRandom::rnd();
    unsigned int itype=1;
    for(;itype<weights_.size();++itype) {
      if(weights_[itype]>=wgtb) break;
      wgtb-=weights_[itype];
    }
    emission_type = itype;
  }
  else {
//   vector<Lorentz5Momentum> pnew;
//   
//   // generate the hard emission and return if no emission
//   if(!getEvent(pnew,emission_type)) {
//     for(unsigned int ix=0;ix<particlesToShower.size();++ix)
//       particlesToShower[ix]->maximumpT(_min_pt);
//     return HardTreePtr();
//   }
    exit(0);
  }



  // construct the HardTree object needed to perform the showers
  ShowerParticleVector newparticles;
  // make the particles for the HardTree
  tcPDPtr gluon=getParticleData(ParticleID::g);
  // create the partons
  // q qbar -> g X
  vector<Lorentz5Momentum> pnew;
  if(emission_type==1||emission_type==3) {
    newparticles.push_back(new_ptr(ShowerParticle(_partons[0]      ,false)));
    newparticles.push_back(new_ptr(ShowerParticle(_partons[1]      ,false)));
    newparticles.push_back(new_ptr(ShowerParticle(gluon            , true)));
    if(emission_type ==1) pnew = realEmissionGluon1_;
    else                  pnew = realEmissionGluon2_;
  }
  // q g    -> q X
  else if(emission_type==4) {
    newparticles.push_back(new_ptr(ShowerParticle(_partons[0]      ,false)));
    newparticles.push_back(new_ptr(ShowerParticle(gluon            ,false)));
    newparticles.push_back(new_ptr(ShowerParticle(_partons[1]->CC(), true)));
    pnew = realEmissionQuark2_;
  }
  // g qbar -> qbar X
  else if(emission_type==2) {
    newparticles.push_back(new_ptr(ShowerParticle(gluon            ,false)));
    newparticles.push_back(new_ptr(ShowerParticle(_partons[1]      ,false)));
    newparticles.push_back(new_ptr(ShowerParticle(_partons[0]->CC(), true)));
    pnew = realEmissionQuark1_;
  }
  else
    assert(false);
  // transform to get directions right
  LorentzRotation R;
  if(!_quarkplus) {
    PPair partons = make_pair(particlesToShower[0]->progenitor(),
			      particlesToShower[1]->progenitor());
    if(partons.first->id()<0) swap(partons.first,partons.second);
    Boost bv = (partons.first->momentum()+
		partons.second->momentum()).boostVector();
    R = LorentzRotation(-bv);
    R.rotateY(-partons.first->momentum().theta());
    R.boost(bv);
    for(unsigned int ix=0;ix<pnew.size();++ix)
      pnew[ix].transform(R);
  }
  // set the momenta
  for(unsigned int ix=0;ix<2;++ix) newparticles[ix]->set5Momentum(pnew[ix]);
  newparticles[2]->set5Momentum(pnew[4]);
  // create the off-shell particle
  Lorentz5Momentum poff=pnew[emission_type>2]-pnew[4];
  poff.rescaleMass();
  newparticles.push_back(new_ptr(ShowerParticle(_partons[emission_type>2],false)));
  newparticles.back()->set5Momentum(poff);
//   // compute the boost 
//   LorentzRotation boost(pboson.findBoostToCM());
//   boost.boost(pnew[3].boostVector());
  for( map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
	 cjt= tree->outgoingLines().begin();
       cjt != tree->outgoingLines().end();++cjt ) {
    newparticles.push_back(new_ptr(ShowerParticle(cjt->first->original()->dataPtr(),
						  true)));
  }
  if(newparticles[4]->id()==mePartonData()[2]->id()) {
    newparticles[4]->set5Momentum(pnew[2]);
    newparticles[5]->set5Momentum(pnew[3]);
  }
  else {
    newparticles[4]->set5Momentum(pnew[3]);
    newparticles[5]->set5Momentum(pnew[2]);
  }
  vector<HardBranchingPtr> nasonin,nasonhard;
  // create the branchings for the incoming particles
  nasonin.push_back(new_ptr(HardBranching(newparticles[0],SudakovPtr(),
					  HardBranchingPtr(),HardBranching::Incoming)));
  nasonin.push_back(new_ptr(HardBranching(newparticles[1],SudakovPtr(),
					  HardBranchingPtr(),HardBranching::Incoming)));
  // intermediate IS particle
  nasonhard.push_back(new_ptr(HardBranching(newparticles[3],SudakovPtr(),
					    nasonin[emission_type>2],HardBranching::Incoming)));
  nasonin[emission_type>2]->addChild(nasonhard.back());
  // create the branching for the emitted jet
  nasonin[emission_type>2]->addChild(new_ptr(HardBranching(newparticles[2],SudakovPtr(),
							   nasonin[emission_type>2],
							   HardBranching::Outgoing)));
  // set the colour partners
  nasonhard.back()->colourPartner(nasonin[emission_type<=2]);
  nasonin[emission_type<=2]->colourPartner(nasonhard.back());
  // add other particle
  nasonhard.push_back(nasonin[emission_type<=2]);
  // outgoing particles
  for(unsigned int ix=4;ix<newparticles.size();++ix) {
    nasonhard.push_back(new_ptr(HardBranching(newparticles[ix],SudakovPtr(),
					      HardBranchingPtr(),HardBranching::Outgoing)));
  }
  // make the tree
  HardTreePtr nasontree=new_ptr(HardTree(nasonhard,nasonin,ShowerInteraction::QCD));
  // connect the ShowerParticles with the branchings
  // and set the maximum pt for the radiation
  Energy pt = pnew[5].perp();
  set<HardBranchingPtr> hard=nasontree->branchings();
  for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
    particlesToShower[ix]->maximumpT(pt);
    for(set<HardBranchingPtr>::const_iterator mit=hard.begin();
	mit!=hard.end();++mit) {
      if(particlesToShower[ix]->progenitor()->id()==(*mit)->branchingParticle()->id()&&
	 (( particlesToShower[ix]->progenitor()->isFinalState()&&
	    (**mit).status()==HardBranching::Outgoing)||
	  (!particlesToShower[ix]->progenitor()->isFinalState()&&
	   (**mit).status()==HardBranching::Incoming))) {
	nasontree->connect(particlesToShower[ix]->progenitor(),*mit);
	if((**mit).status()==HardBranching::Incoming) {
	  (*mit)->beam(particlesToShower[ix]->original()->parents()[0]);
	}
	HardBranchingPtr parent=(*mit)->parent();
	while(parent) {
	  parent->beam(particlesToShower[ix]->original()->parents()[0]);
	  parent=parent->parent();
	};
      }
    }
  }
  ColinePtr newline=new_ptr(ColourLine());
  for(set<HardBranchingPtr>::const_iterator cit=nasontree->branchings().begin();
      cit!=nasontree->branchings().end();++cit) {
    if((**cit).branchingParticle()->dataPtr()->iColour()==PDT::Colour3)
      newline->addColoured((**cit).branchingParticle());
    else if((**cit).branchingParticle()->dataPtr()->iColour()==PDT::Colour3bar)
      newline->addAntiColoured((**cit).branchingParticle());
  }
  // return the tree
  return nasontree;
}
