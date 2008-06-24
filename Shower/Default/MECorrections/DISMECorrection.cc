// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DISMECorrection class.
//

#include "DISMECorrection.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include <numeric>
#include "Herwig++/Utilities/Maths.h"

using namespace Herwig;

void DISMECorrection::persistentOutput(PersistentOStream & os) const {
  os << _meopt << _comptonint << _bgfint << _procprob << _sinW << _cosW 
     << ounit(_mz2,GeV2) << _initial << _final;
}

void DISMECorrection::persistentInput(PersistentIStream & is, int) {
  is >> _meopt >> _comptonint >> _bgfint >> _procprob >> _sinW >> _cosW 
     >> iunit(_mz2,GeV2) >> _initial >> _final;
}

ClassDescription<DISMECorrection> DISMECorrection::initDISMECorrection;
// Definition of the static class description member.

void DISMECorrection::Init() {

  static ClassDocumentation<DISMECorrection> documentation
    ("The DISMECorrection class implements the matrix element correction for DIS");

  static Parameter<DISMECorrection,double> interfaceProcessProbability
    ("ProcessProbability",
     "The probabilty of the QCD compton process for the process selection",
     &DISMECorrection::_procprob, 0.3, 0.0, 1.,
     false, false, Interface::limited);

  static Switch<DISMECorrection,bool> interfaceMEOption
    ("MEOption",
     "Option for the treatment of the matrix element",
     &DISMECorrection::_meopt, false, false);
  static SwitchOption interfaceMEOptionFull
    (interfaceMEOption,
     "Full",
     "Use the full matrix element",
     true);
  static SwitchOption interfaceMEOptionProcessIndependent
    (interfaceMEOption,
     "ProcessIndependent",
     "Only use the process independent part, as in FORTRAN HERWIG",
     false);

}

void DISMECorrection::doinit() throw(InitException) {
  QTildeMECorrection::doinit();
  // electroweak parameters
  _sinW = generator()->standardModel()->sin2ThetaW();
  _cosW = sqrt(1.-_sinW);
  _sinW = sqrt(_sinW);
  _mz2 = sqr(getParticleData(ParticleID::Z0)->mass());
  // integrals of me over phase space
  double r5=sqrt(5.),darg((r5-1.)/(r5+1.)),ath(0.5*log((1.+1./r5)/(1.-1./r5)));
  _comptonint = 2.*(-21./20.-6.*r5/25.*ath+sqr(Constants::pi)/3.
		    -2.*Math::ReLi2(1.-darg)-2.*Math::ReLi2(1.-1./darg));
  _bgfint = 121./9.-56./5.*r5*ath;
}

void DISMECorrection::dofinish() {
  MECorrectionBase::dofinish();
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream outfile(fname.c_str());
  outfile << "SET FONT DUPLEX\n";
  outfile << "SET ORDER X Y\n";
  outfile << "SET WINDOW X 2 9 Y 2 9\n";
  outfile << "TITLE BOTTOM \"x0p1\"\n";
  outfile << "CASE         \" X X\"\n";
  outfile << "TITLE LEFT \"z0p1\"\n";
  outfile << "CASE       \" X X\"\n";
  outfile << "SET LIMITS X 0 1 Y 0 1\n";
  int nfortran(0);
  for(unsigned int ix=0;ix<_compton.size();++ix) {
    double xp = _compton[ix].first, zp = _compton[ix].second;
    double xpmax = (1.+4.*zp*(1.-zp))/(1.+6.*zp*(1.-zp));
    Complex a = 10.-45.*xp+18.*sqr(xp)+3.*Complex(0.,1.)*
      sqrt(3.*(9.+66.*xp-93.*sqr(xp)+12.*pow(xp,3)-8.*pow(xp,4)
	       +24.*pow(xp,5)-8.*pow(xp,6)));
    Complex b = pow(a,1./3.)*(0.5+sqrt(3.)/2.*Complex(0.,1.));
    double zpmax = 1.-2./3.*xp*(1.+b.real());
    outfile << xp << "\t" << zp << "\n";
    if(xp<xpmax&&zp<zpmax) ++nfortran;
    if(ix%50000==0&&ix>0) outfile << "PLOT\n";
  }
  if(!_compton.empty()) outfile << "PLOT\n";
  cerr << nfortran << " would have been accepted by FORTRAN HERWIG for the "
       << "compton process of " << _compton.size() << " accepted by the C++\n";
  for(double xp=0;xp<=1.001;xp+=0.01) {
    outfile << xp << "\t" << 1./(1.+xp-sqr(xp)) << "\n";
  }
  outfile << "JOIN BLUE\n";
  for(double z=0;z<=1.001;z+=0.01) {
    outfile << 1./(1+z*(1.-z)) << "\t" << z << "\n";
  }
  outfile << "JOIN BLUE\n";
  for(double zp=0.;zp<=1.001;zp+=0.01) {
    outfile << (1.+4.*zp*(1.-zp))/(1.+6.*zp*(1.-zp)) << "\t" << zp << "\n";
  }
  outfile << "JOIN GREEN\n";
  for(double xp=0.;xp<=1.001;xp+=0.01) {
    Complex a = 10.-45.*xp+18.*sqr(xp)+3.*Complex(0.,1.)*
      sqrt(3.*(9.+66.*xp-93.*sqr(xp)+12.*pow(xp,3)-8.*pow(xp,4)
	       +24.*pow(xp,5)-8.*pow(xp,6)));
    Complex b = pow(a,1./3.)*(0.5+sqrt(3.)/2.*Complex(0.,1.));
    outfile << xp << "\t" << 1.-2./3.*xp*(1.+b.real()) << "\n";
  }
  outfile << "JOIN GREEN\n";
  outfile << "NEW FRAME\n";
  outfile << "SET FONT DUPLEX\n";
  outfile << "SET ORDER X Y\n";
  outfile << "SET WINDOW X 2 9 Y 2 9\n";
  outfile << "TITLE BOTTOM \"x0p1\"\n";
  outfile << "CASE         \" X X\"\n";
  outfile << "TITLE LEFT \"z0p1\"\n";
  outfile << "CASE       \" X X\"\n";
  outfile << "SET LIMITS X 0 1 Y 0 1\n";
  for(unsigned int ix=0;ix<_comptonover.size();++ix) {
    outfile << _comptonover[ix].first << "\t" << _comptonover[ix].second << "\n";
    if(ix%50000==0&&ix>0) outfile << "PLOT RED\n";
  }
  if(!_comptonover.empty()) outfile << "PLOT RED\n";
  for(double xp=0;xp<=1.001;xp+=0.01) {
    outfile << xp << "\t" << 1./(1.+xp-sqr(xp)) << "\n";
  }
  outfile << "JOIN BLUE\n";
  for(double z=0;z<=1.001;z+=0.01) {
    outfile << 1./(1+z*(1.-z)) << "\t" << z << "\n";
  }
  outfile << "JOIN BLUE\n";
  for(double zp=0.;zp<=1.001;zp+=0.01) {
    outfile << (1.+4.*zp*(1.-zp))/(1.+6.*zp*(1.-zp)) << "\t" << zp << "\n";
  }
  outfile << "JOIN GREEN\n";
  for(double xp=0.;xp<=1.001;xp+=0.01) {
    Complex a = 10.-45.*xp+18.*sqr(xp)+3.*Complex(0.,1.)*
      sqrt(3.*(9.+66.*xp-93.*sqr(xp)+12.*pow(xp,3)-8.*pow(xp,4)
	       +24.*pow(xp,5)-8.*pow(xp,6)));
    Complex b = pow(a,1./3.)*(0.5+sqrt(3.)/2.*Complex(0.,1.));
    outfile << xp << "\t" << 1.-2./3.*xp*(1.+b.real()) << "\n";
  }
  outfile << "JOIN GREEN\n";
  outfile << "NEW FRAME\n";
  outfile << "SET FONT DUPLEX\n";
  outfile << "SET ORDER X Y\n";
  outfile << "SET WINDOW X 2 9 Y 2 9\n";
  outfile << "TITLE BOTTOM \"x0p1\"\n";
  outfile << "CASE         \" X X\"\n";
  outfile << "TITLE LEFT \"z0p1\"\n";
  outfile << "CASE       \" X X\"\n";
  outfile << "SET LIMITS X 0 1 Y 0 1\n";
  nfortran=0; 
  for(unsigned int ix=0;ix<_bgf.size();++ix) {
    double xp = _bgf[ix].first, zp = _bgf[ix].second;
    Complex a = 10.-45.*xp+18.*sqr(xp)+3.*Complex(0.,1.)*
      sqrt(3.*(9.+66.*xp-93.*sqr(xp)+12.*pow(xp,3)-8.*pow(xp,4)
	       +24.*pow(xp,5)-8.*pow(xp,6)));
    Complex b = pow(a,1./3.)*(0.5+sqrt(3.)/2.*Complex(0.,1.));
    double zpmax = 1.-2./3.*xp*(1.+b.real()),zpmin=1.-zpmax;
    if(zp>zpmin&&zp<zpmax) ++nfortran;
    outfile << _bgf[ix].first << "\t" << _bgf[ix].second << "\n";
    if(ix%50000==0&&ix>0) outfile << "PLOT\n";
  }
  cerr << nfortran << " would have been accepted by FORTRAN HERWIG for the "
       << "BGF process of " << _bgf.size() << " accepted by the C++\n";
  if(!_bgf.empty()) outfile << "PLOT\n";
  for(double z=0;z<=1.001;z+=0.01) {
    double xp = 2.*z/(1.+(1.-z)+sqrt(sqr(1.+(1.-z))-4.*z*(1.-z)));
    double zp = 0.5* (1.-(1.-z)+sqrt(sqr(1.+(1.-z))-4.*z*(1.-z)));
    outfile << xp << "\t" << zp << "\n";
  }
  outfile << "JOIN BLUE\n";
  for(double z=0;z<=1.001;z+=0.01) {
    double xp = 2.*z/(1.+(1.-z)+sqrt(sqr(1.+(1.-z))-4.*z*(1.-z)));
    double zp = 0.5* (1.-(1.-z)+sqrt(sqr(1.+(1.-z))-4.*z*(1.-z)));
    outfile << xp << "\t" << 1.-zp << "\n";
  }
  outfile << "JOIN BLUE\n";
  for(double xp=0.;xp<=1.001;xp+=0.01) {
    Complex a = 10.-45.*xp+18.*sqr(xp)+3.*Complex(0.,1.)*
      sqrt(3.*(9.+66.*xp-93.*sqr(xp)+12.*pow(xp,3)-8.*pow(xp,4)
	       +24.*pow(xp,5)-8.*pow(xp,6)));
    Complex b = pow(a,1./3.)*(0.5+sqrt(3.)/2.*Complex(0.,1.));
    outfile << xp << "\t" << 1.-2./3.*xp*(1.+b.real()) << "\n";
  }
  outfile << "JOIN GREEN\n";
  for(double xp=0.;xp<=1.001;xp+=0.01) {
    Complex a = 10.-45.*xp+18.*sqr(xp)+3.*Complex(0.,1.)*
      sqrt(3.*(9.+66.*xp-93.*sqr(xp)+12.*pow(xp,3)-8.*pow(xp,4)
	       +24.*pow(xp,5)-8.*pow(xp,6)));
    Complex b = pow(a,1./3.)*(0.5+sqrt(3.)/2.*Complex(0.,1.));
    outfile << xp << "\t" << +2./3.*xp*(1.+b.real()) << "\n";
  }
  outfile << "JOIN GREEN\n";
  outfile << "NEW FRAME\n";
  outfile << "SET FONT DUPLEX\n";
  outfile << "SET ORDER X Y\n";
  outfile << "SET WINDOW X 2 9 Y 2 9\n";
  outfile << "TITLE BOTTOM \"x0p1\"\n";
  outfile << "CASE         \" X X\"\n";
  outfile << "TITLE LEFT \"z0p1\"\n";
  outfile << "CASE       \" X X\"\n";
  outfile << "SET LIMITS X 0 1 Y 0 1\n";
  for(unsigned int ix=0;ix<_bgfover.size();++ix) {
    outfile << _bgfover[ix].first << "\t" << _bgfover[ix].second << "\n";
    if(ix%50000==0&&ix>0) outfile << "PLOT RED\n";
  }
  if(!_bgfover.empty()) outfile << "PLOT RED\n";
  for(double z=0;z<=1.001;z+=0.01) {
    double xp = 2.*z/(1.+(1.-z)+sqrt(sqr(1.+(1.-z))-4.*z*(1.-z)));
    double zp = 0.5* (1.-(1.-z)+sqrt(sqr(1.+(1.-z))-4.*z*(1.-z)));
    outfile << xp << "\t" << zp << "\n";
  }
  outfile << "JOIN BLUE\n";
  for(double z=0;z<=1.001;z+=0.01) {
    double xp = 2.*z/(1.+(1.-z)+sqrt(sqr(1.+(1.-z))-4.*z*(1.-z)));
    double zp = 0.5* (1.-(1.-z)+sqrt(sqr(1.+(1.-z))-4.*z*(1.-z)));
    outfile << xp << "\t" << 1.-zp << "\n";
  }
  outfile << "JOIN BLUE\n";
  for(double xp=0.;xp<=1.001;xp+=0.01) {
    Complex a = 10.-45.*xp+18.*sqr(xp)+3.*Complex(0.,1.)*
      sqrt(3.*(9.+66.*xp-93.*sqr(xp)+12.*pow(xp,3)-8.*pow(xp,4)
	       +24.*pow(xp,5)-8.*pow(xp,6)));
    Complex b = pow(a,1./3.)*(0.5+sqrt(3.)/2.*Complex(0.,1.));
    outfile << xp << "\t" << 1.-2./3.*xp*(1.+b.real()) << "\n";
  }
  outfile << "JOIN GREEN\n";
  for(double xp=0.;xp<=1.001;xp+=0.01) {
    Complex a = 10.-45.*xp+18.*sqr(xp)+3.*Complex(0.,1.)*
      sqrt(3.*(9.+66.*xp-93.*sqr(xp)+12.*pow(xp,3)-8.*pow(xp,4)
	       +24.*pow(xp,5)-8.*pow(xp,6)));
    Complex b = pow(a,1./3.)*(0.5+sqrt(3.)/2.*Complex(0.,1.));
    outfile << xp << "\t" << +2./3.*xp*(1.+b.real()) << "\n";
  }
  outfile << "JOIN GREEN\n";
  outfile << "NEW FRAME\n";
  double xB=-0.005;
  for(unsigned int ix=0;ix<101;++ix) {
    xB+=0.01;
    if(_comptonxb[ix].second!=0.) 
      outfile << xB << "\t" << _comptonxb[ix].first/_comptonxb[ix].second << "\n";
    else
      outfile << xB << "\t" << 0. << "\n";
  }
  outfile << "HIST RED\n";
  xB=-0.005;
  for(unsigned int ix=0;ix<101;++ix) {
    xB+=0.01;
    if(_bgfxb[ix].second!=0.) 
      outfile << xB << "\t" << _bgfxb[ix].first/_bgfxb[ix].second << "\n";
    else
      outfile << xB << "\t" << 0. << "\n";
  }
  outfile << "HIST BLUE\n";


  outfile.close();
  if(_ntry==0) return;
  generator()->log() << "DISMECorrection when applying the hard correction "
		     << "generated " << _ntry << " trial emissions of which "
		     << _ngen << " were accepted\n";
  if(_nover==0) return;
  generator()->log() << "DISMECorrection when applying the hard correction " 
		     << _nover << " weights larger than one were generated of which"
		     << " the largest was " << _maxwgt.first << " for the QCD compton"
		     << " processes and " << _maxwgt.second << " for the BGF process\n";
}

bool DISMECorrection::canHandle(ShowerTreePtr tree,double & initial,
				double & final,EvolverPtr) {
  // two incoming particles
  if(tree->incomingLines().size()!=2) return false;
  // two outgoing particles
  if(tree->outgoingLines().size()!=2) return false;
  // check incoming quark and lepton
  bool quark(false),lepton(false);
  for(map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator 
	cit=tree->incomingLines().begin();cit!=tree->incomingLines().end();++cit) {
    quark  |=  QuarkMatcher::Check(cit->first->progenitor()->data());
    lepton |= LeptonMatcher::Check(cit->first->progenitor()->data());
  }
  if(!quark||!lepton) return false;
  // check outgoing quark and lepton
  quark = false;
  lepton = false;
  for(map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
	cjt=tree->outgoingLines().begin();cjt!=tree->outgoingLines().end();++cjt) {
    quark  |=  QuarkMatcher::Check(cjt->first->progenitor()->data());
    lepton |= LeptonMatcher::Check(cjt->first->progenitor()->data());
  }
  if(!quark||!lepton) return false;
  // can handle it
  initial = _initial;
  final   = _final;
  return true;
}

void DISMECorrection::applyHardMatrixElementCorrection(ShowerTreePtr tree) {
  ++_ntry;
  // find the incoming and outgoing quarks and leptons
  ShowerParticlePtr quark[2],lepton[2];
  PPtr hadron;
  tcBeamPtr beam;
  // incoming particles
  for(map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator 
	cit=tree->incomingLines().begin();cit!=tree->incomingLines().end();++cit) {
    if(QuarkMatcher::Check(cit->first->progenitor()->data())) {
      hadron = cit->first->original()->parents()[0];
      quark [0] = cit->first->progenitor();
      beam = cit->first->beam();
    }
    else if(LeptonMatcher::Check(cit->first->progenitor()->data())) {
      lepton[0] = cit->first->progenitor();
    }
  }
  tcPDFPtr pdf=beam->pdf();
  assert(beam&&pdf&&quark[0]&&lepton[0]);
  // outgoing particles
  for(map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
	cit=tree->outgoingLines().begin();cit!=tree->outgoingLines().end();++cit) {
    if(QuarkMatcher::Check(cit->first->progenitor()->data()))
      quark [1] = cit->first->progenitor();
    else if(LeptonMatcher::Check(cit->first->progenitor()->data()))
      lepton[1] = cit->first->progenitor();
  }
  assert(quark[1]&&lepton[1]);
  // extract the born variables
  Lorentz5Momentum q =lepton[0]->momentum()-lepton[1]->momentum();
  _q2 = -q.m2();
  double  xB = quark[0]->x();
  double  yB = 
    (                    q*quark[0]->momentum())/
    (lepton[0]->momentum()*quark[0]->momentum()); 
  _l = 2./yB-1.;
  // calculate the A coefficient for the correlations
  _acoeff = A(lepton[0]->dataPtr(),lepton[1]->dataPtr(),
	      quark [0]->dataPtr(),quark [1]->dataPtr());
  vector<double> azicoeff;
  double xp,zp,wgt,x1,x2,x3,xperp;
  tcPDPtr gluon = getParticleData(ParticleID::g);
  if(abs(quark[0]->id())>2) {
    _procprob = 0.34;
  }
  else if(abs(quark[0]->id())==1) {
    _procprob = 0.36;
  }
  else if(abs(quark[0]->id())==2) {
    if(quark[0]->id()*hadron->id()>0) {
      _procprob = 0.40;
    }
    else {
      _procprob = 0.35;
    }
  }
  // select the type of process
  bool BGF = UseRandom::rnd()>_procprob;
  //bool BGF=true;
  // generate a QCD compton process
  if(!BGF) {
    wgt = generateComptonPoint(xp,zp);
    if(xp<1e-6) return;
    // common pieces
    Energy2 scale = _q2*((1.-xp)*(1-zp)*zp/xp+1);
    wgt *= 2./3./Constants::pi*coupling()->value(scale)/_procprob;
    // PDF piece
    wgt *= pdf->xfx(beam,quark[0]->dataPtr(),scale,xB/xp)/
           pdf->xfx(beam,quark[0]->dataPtr(),_q2  ,xB);
    // other bits
    xperp = sqrt(4.*(1.-xp)*(1.-zp)*zp/xp);
    x1 = -1./xp;
    x2 = 1.-(1.-zp)/xp;
    x3 = 2.+x1-x2;
    // matrix element pieces
    azicoeff = ComptonME(xp,x2,xperp,_acoeff,_l,true);
  }
  // generate a BGF process
  else {
    wgt = generateBGFPoint(xp,zp);
    if(xp<1e-6) return;
    // common pieces 
    Energy2 scale = _q2*((1.-xp)*(1-zp)*zp/xp+1);
    wgt *= 0.25/Constants::pi*coupling()->value(scale)/(1.-_procprob);
    // PDF piece
    wgt *= pdf->xfx(beam,gluon              ,scale,xB/xp)/
           pdf->xfx(beam,quark[0]->dataPtr(),_q2  ,xB);
    // other bits
    xperp = sqrt(4.*(1.-xp)*(1.-zp)*zp/xp);
    x1 = -1./xp;
    x2 = 1.-(1.-zp)/xp;
    x3 = 2.+x1-x2;
    // matrix element pieces
    azicoeff = BGFME(xp,x2,x3,xperp,_acoeff,_l,true);
  }
  // compute the azimuthal average of the weight
  wgt *= (azicoeff[0]+0.5*azicoeff[2]);
  // decide whether or not to accept the weight
  if(UseRandom::rnd()>wgt) return;
  // if generate generate phi
  unsigned int itry(0);
  double phimax = std::accumulate(azicoeff.begin(),azicoeff.end(),0.);
  double phiwgt,phi;
  do {
    phi = UseRandom::rnd()*Constants::twopi;
    double cphi(cos(phi));
    phiwgt = azicoeff[0]+azicoeff[1]*cphi+azicoeff[2]*sqr(cphi);
    ++itry;
  }
  while (phimax*UseRandom::rnd() > phiwgt && itry<200);
  if(itry==200) throw Exception() << "Too many tries in DISMECorrection"
				  << "::applyHardMatrixElementCorrection() to"
				  << " generate phi" << Exception::eventerror;
  // compute the new incoming and outgoing momenta
  Energy Q(sqrt(_q2));
  Lorentz5Momentum p1 = Lorentz5Momentum( 0.5*Q*xperp*cos(phi), 0.5*Q*xperp*sin(phi),
					  -0.5*Q*x2,0.*GeV,0.*GeV);
  p1.rescaleEnergy();
  Lorentz5Momentum p2 = Lorentz5Momentum(-0.5*Q*xperp*cos(phi),-0.5*Q*xperp*sin(phi),
					 -0.5*Q*x3,0.*GeV,0.*GeV);
  p2.rescaleEnergy();
  Lorentz5Momentum pin(0.*GeV,0.*GeV,-0.5*x1*Q,-0.5*x1*Q,0.*GeV);
  // construct lorentz transform from lab to breit frame
  Lorentz5Momentum phadron =  hadron->momentum();
  phadron.setMass(0.*GeV);
  phadron.rescaleEnergy();
  Lorentz5Momentum pcmf = phadron+0.5/xB*q;
  pcmf.rescaleMass();
  LorentzRotation rot(-pcmf.boostVector());
  Lorentz5Momentum pbeam = rot*phadron;
  Axis axis(pbeam.vect().unit());
  double sinth(sqrt(1.-sqr(axis.z())));
  rot.rotate(-acos(axis.z()),Axis(-axis.y()/sinth,axis.x()/sinth,0.));
  Lorentz5Momentum pl    = rot*lepton[0]->momentum();
  rot.rotateZ(-atan2(pl.y(),pl.x()));
  // we need the Lorentz transform back to the lab
  rot.invert();
  // transform the momenta to lab frame
  pin *= rot;
  p1  *= rot;
  p2  *= rot;
  // test to ensure outgoing particles can be put on-shell
  if(!BGF) {
    if(p1.e()<quark[1]->dataPtr()->constituentMass()) return;
    if(p2.e()<gluon              ->constituentMass()) return;
  }
  else {
    if(p1.e()<quark[1]->dataPtr()      ->constituentMass()) return;
    if(p2.e()<quark[0]->dataPtr()->CC()->constituentMass()) return;
  }
  // stats for weights > 1
  if(wgt>1.) {
    ++_nover;
    if(!BGF) {
      _maxwgt.first  = max(_maxwgt.first ,wgt);
      _comptonover.push_back(make_pair(xp,zp));
    }
    else {
      _maxwgt.second = max(_maxwgt.second,wgt);
      _bgfover.push_back(make_pair(xp,zp));
    }
  }
  if(!BGF) {
    if(_comptonxb.empty()) _comptonxb.resize(101,make_pair(0.,0.));
    int iloc = int(xB/0.01);
    _comptonxb[iloc].first  += wgt;
    _comptonxb[iloc].second += 1. ;
  }
  else {
    if(_bgfxb.empty()) _bgfxb    .resize(101,make_pair(0.,0.));
    int iloc = int(xB/0.01);
    _bgfxb[iloc].first  += wgt;
    _bgfxb[iloc].second += 1. ;
  }
  // points into hist
  ++_ngen;
  if(!BGF) _compton.push_back(make_pair(xp,zp));
  else     _bgf    .push_back(make_pair(xp,zp));
  // create the new particles and add to ShowerTree
  bool isquark = quark[0]->colourLine();
  if(!BGF) {
    PPtr newin  = new_ptr(Particle(*quark[0]));
    newin->set5Momentum(pin);
    PPtr newg   = gluon              ->produceParticle(p2 );
    PPtr newout = quark[1]->dataPtr()->produceParticle(p1 ); 
    ColinePtr col=isquark ? 
      quark[0]->colourLine() : quark[0]->antiColourLine();
    ColinePtr newline=new_ptr(ColourLine());
    // final-state emission
    if(xp>zp) {
      col->removeColoured(newout,!isquark);
      col->addColoured(newin,!isquark);
      col->addColoured(newg,!isquark);
      newline->addColoured(newg,isquark);
      newline->addColoured(newout,!isquark);
    }
    // initial-state emission
    else {
      col->removeColoured(newin ,!isquark);
      col->addColoured(newout,!isquark);
      col->addColoured(newg,isquark);
      newline->addColoured(newg,!isquark);
      newline->addColoured(newin,!isquark);
    }
    PPtr orig;
    for(map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator 
	  cit=tree->incomingLines().begin();cit!=tree->incomingLines().end();++cit) {
      if(cit->first->progenitor()!=quark[0]) continue;
      // remove old particles from colour line
      col->removeColoured(cit->first->copy(),!isquark);
      col->removeColoured(cit->first->progenitor(),!isquark);
      // insert new particles
      cit->first->copy(newin);
      ShowerParticlePtr sp(new_ptr(ShowerParticle(*newin,1,false)));
      cit->first->progenitor(sp);
      tree->incomingLines()[cit->first]=sp;
      sp->x(xB/xp);
      cit->first->perturbative(xp>zp);
      if(xp<=zp) orig=cit->first->original();
    }
    for(map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
	  cit=tree->outgoingLines().begin();cit!=tree->outgoingLines().end();++cit) {
      if(cit->first->progenitor()!=quark[1]) continue;
      // remove old particles from colour line
      col->removeColoured(cit->first->copy(),!isquark);
      col->removeColoured(cit->first->progenitor(),!isquark);
      // insert new particles
      cit->first->copy(newout);
      ShowerParticlePtr sp(new_ptr(ShowerParticle(*newout,1,true)));
      cit->first->progenitor(sp);
      tree->outgoingLines()[cit->first]=sp;
      cit->first->perturbative(xp<=zp);
      if(xp>zp) orig=cit->first->original();
    }
    assert(orig);
    // add the gluon
    ShowerParticlePtr sg=new_ptr(ShowerParticle(*newg,1,true));
    ShowerProgenitorPtr gluon=new_ptr(ShowerProgenitor(orig,newg,sg));
    gluon->perturbative(false);
    tree->outgoingLines().insert(make_pair(gluon,sg));
    tree->hardMatrixElementCorrection(true);
  }
  else {
    PPtr newin   = gluon                    ->produceParticle(pin);
    PPtr newqbar = quark[0]->dataPtr()->CC()->produceParticle(p2 );
    PPtr newout  = quark[1]->dataPtr()      ->produceParticle(p1 );
    ColinePtr col=isquark ? quark[0]->colourLine() : quark[0]->antiColourLine();
    ColinePtr newline=new_ptr(ColourLine()); 
    col    ->addColoured(newin  ,!isquark);
    newline->addColoured(newin  , isquark);
    col    ->addColoured(newout ,!isquark);
    newline->addColoured(newqbar, isquark);
    PPtr orig;
    for(map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator 
	  cit=tree->incomingLines().begin();cit!=tree->incomingLines().end();++cit) {
      if(cit->first->progenitor()!=quark[0]) continue;
      // remove old particles from colour line
      col->removeColoured(cit->first->copy(),!isquark);
      col->removeColoured(cit->first->progenitor(),!isquark);
      // insert new particles
      cit->first->copy(newin);
      ShowerParticlePtr sp(new_ptr(ShowerParticle(*newin,1,false)));
      cit->first->progenitor(sp);
      tree->incomingLines()[cit->first]=sp;
      sp->x(xB/xp);
      cit->first->perturbative(false);
      orig=cit->first->original();
    }
    for(map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
	  cit=tree->outgoingLines().begin();cit!=tree->outgoingLines().end();++cit) {
      if(cit->first->progenitor()!=quark[1]) continue;
      // remove old particles from colour line
      col->removeColoured(cit->first->copy(),!isquark);
      col->removeColoured(cit->first->progenitor(),!isquark);
      // insert new particles
      cit->first->copy(newout);
      ShowerParticlePtr sp(new_ptr(ShowerParticle(*newout,1,true)));
      cit->first->progenitor(sp);
      tree->outgoingLines()[cit->first]=sp;
      cit->first->perturbative(true);
    }
    assert(orig);
    // add the (anti)quark
    ShowerParticlePtr sqbar=new_ptr(ShowerParticle(*newqbar,1,true));
    ShowerProgenitorPtr qbar=new_ptr(ShowerProgenitor(orig,newqbar,sqbar));
    qbar->perturbative(false);
    tree->outgoingLines().insert(make_pair(qbar,sqbar));
    tree->hardMatrixElementCorrection(true);
  }
}

bool DISMECorrection::softMatrixElementVeto(ShowerProgenitorPtr initial,
					    ShowerParticlePtr parent,Branching br) {
  bool veto = !UseRandom::rndbool(parent->isFinalState() ? 1./_final : 1./_initial);
  // check if me correction should be applied
  long id[2]={initial->id(),parent->id()};
  if(id[0]!=id[1]||id[1]==ParticleID::g) return veto;
  // get the pT
  Energy pT=br.kinematics->pT();
  // check if hardest so far
  if(pT<initial->highestpT()) return veto;
  double kappa(sqr(br.kinematics->scale())/_q2),z(br.kinematics->z());
  double zk((1.-z)*kappa);
  // final-state
  double wgt(0.);
  if(parent->isFinalState()) {
    double zp=z,xp=1./(1.+z*zk);
    double xperp = sqrt(4.*(1.-xp)*(1.-zp)*zp/xp);
    double x2 = 1.-(1.-zp)/xp;
    vector<double> azicoeff = ComptonME(xp,x2,xperp,_acoeff,_l,false);
    wgt = (azicoeff[0]+0.5*azicoeff[2])*xp/(1.+sqr(z))/_final;
    if(wgt<.0||wgt>1.) 
      generator()->log() << "Soft ME correction weight too large or "
			 << "negative for FSR in DISMECorrection::"
			 << "softMatrixElementVeto() soft weight " 
			 << " xp = " << xp
			 << " zp = " << zp
			 << " weight = " << wgt << "\n";
  }
  else {
    double xp = 2.*z/(1.+zk+sqrt(sqr(1.+zk)-4.*z*zk));
    double zp = 0.5* (1.-zk+sqrt(sqr(1.+zk)-4.*z*zk));
    // compton
    if(br.ids[0]!=ParticleID::g) {
      double xperp = sqrt(4.*(1.-xp)*(1.-zp)*zp/xp);
      double x2 = 1.-(1.-zp)/xp;
      vector<double> azicoeff = ComptonME(xp,x2,xperp,_acoeff,_l,false);
      wgt = (azicoeff[0]+0.5*azicoeff[2])*xp*(1.-z)/(1.-xp)/(1.+sqr(z))/
	(1.-zp+xp-2.*xp*(1.-zp));
    }
    // BGF
    else {
      double xperp = sqrt(4.*(1.-xp)*(1.-zp)*zp/xp);
      double x1 = -1./xp, x2 = 1.-(1.-zp)/xp, x3 = 2.+x1-x2;    
      vector<double> azicoeff = BGFME(xp,x2,x3,xperp,_acoeff,_l,true);
      wgt = (azicoeff[0]+0.5*azicoeff[2])*xp/(1.-zp+xp-2.*xp*(1.-zp))/(sqr(z)+sqr(1.-z));
    }
    wgt /=_initial;
    if(wgt<.0||wgt>1.) 
      generator()->log() << "Soft ME correction weight too large or "
			 << "negative for ISR in DISMECorrection::"
			 << "softMatrixElementVeto() soft weight " 
			 << " xp = " << xp
			 << " zp = " << zp
			 << " weight = " << wgt << "\n";
  }
  // if not vetoed
  if(UseRandom::rndbool(wgt)) {
    initial->highestpT(pT);
    return false;
  }
  // otherwise
  parent->setEvolutionScale(ShowerIndex::QCD,br.kinematics->scale());
  return true;
}
