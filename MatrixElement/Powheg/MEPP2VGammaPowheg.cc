// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2VGammaPowheg class.
//

#include "MEPP2VGammaPowheg.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "Herwig++/Utilities/Maths.h"
#include <numeric>

using namespace Herwig;
using Herwig::Math::ReLi2;

MEPP2VGammaPowheg::MEPP2VGammaPowheg() 
  :  _contrib(1), _scaleopt(0),
     _fixedScale(100.*GeV), _scaleFact(1.)
{}

void MEPP2VGammaPowheg::doinit() {
  // colour factors
  CF_ = 4./3.; 
  TR_ = 0.5;
  static const tcHwSMPtr hwsm
    = dynamic_ptr_cast<tcHwSMPtr>(generator()->standardModel());
  if (!hwsm) throw InitException() << "hwsm pointer is null in"
				   << " MEPP2VGamma::doinit()"
				   << Exception::abortnow;
  // get pointers to all required Vertex objects
  FFZvertex_ = hwsm->vertexFFZ();
  FFPvertex_ = hwsm->vertexFFP();
  WWWvertex_ = hwsm->vertexWWW();
  FFWvertex_ = hwsm->vertexFFW();
  FFGvertex_ = hwsm->vertexFFG();

  MEPP2VGamma::doinit();
}


IBPtr MEPP2VGammaPowheg::clone() const {
  return new_ptr(*this);
}

IBPtr MEPP2VGammaPowheg::fullclone() const {
  return new_ptr(*this);
}

void MEPP2VGammaPowheg::persistentOutput(PersistentOStream & os) const {
  os << _contrib << _scaleopt << ounit(_fixedScale,GeV) << _scaleFact     
     << FFPvertex_ << FFWvertex_ << FFZvertex_ << WWWvertex_ << FFGvertex_
     << TR_ << CF_;
}

void MEPP2VGammaPowheg::persistentInput(PersistentIStream & is, int) {
  is >> _contrib >> _scaleopt >> iunit(_fixedScale,GeV) >> _scaleFact     
     >> FFPvertex_ >> FFWvertex_ >> FFZvertex_ >> WWWvertex_ >> FFGvertex_
     >> TR_ >> CF_;
}

ClassDescription<MEPP2VGammaPowheg> MEPP2VGammaPowheg::initMEPP2VGammaPowheg;
// Definition of the static class description member.

void MEPP2VGammaPowheg::Init() {

  static ClassDocumentation<MEPP2VGammaPowheg> documentation
    ("The MEPP2VGammaPowheg class implements the NLO matrix"
     " elements for q qbar -> W/Z gamma in the POWHEG scheme.");

   static Switch<MEPP2VGammaPowheg,unsigned int> interfaceContribution
    ("Contribution",
     "Which contributions to the cross section to include",
     &MEPP2VGammaPowheg::_contrib, 1, false, false);
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

  static Switch<MEPP2VGammaPowheg,unsigned int> interfaceFactorizationScaleOption
    ("FactorizationScaleOption",
     "Option for the scale to be used",
     &MEPP2VGammaPowheg::_scaleopt, 0, false, false);
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

  static Parameter<MEPP2VGammaPowheg,Energy> interfaceFactorizationScaleValue
    ("FactorizationScaleValue",
     "The fixed scale to use if required",
     &MEPP2VGammaPowheg::_fixedScale, GeV, 100.0*GeV, 10.0*GeV, 1000.0*GeV,
     false, false, Interface::limited);

  static Parameter<MEPP2VGammaPowheg,double> interfaceScaleFactor
    ("ScaleFactor",
     "The factor used before sHat if using a running scale",
     &MEPP2VGammaPowheg::_scaleFact, 1.0, 0.0, 10.0,
     false, false, Interface::limited);
}

Energy2 MEPP2VGammaPowheg::scale() const {
  return _scaleopt == 0 ? sqr(_fixedScale) : _scaleFact*sHat();
}

int MEPP2VGammaPowheg::nDim() const {
  return MEPP2VGamma::nDim()+3;
}

bool MEPP2VGammaPowheg::generateKinematics(const double * r) {
  _xtil = *(r +nDim() -1);
  _y1 = *(r +nDim() -3);
  _y2 = *(r +nDim() -2);
  _y3 = *(r +nDim() -1);
  return MEPP2VGamma::generateKinematics(r);
}

CrossSection MEPP2VGammaPowheg::dSigHatDR() const {
  // momentum fractions of the incoming partons
  _xa=  lastX1();
  _xb=  lastX2();
  //incoming parton ParticleData object and momenta
  _partona= mePartonData()[0];
  _partonb= mePartonData()[1];
  _p_partona= meMomenta()[0];
  _p_partonb= meMomenta()[1];
  // outgoing ParticleData object and momenta
  _p_photon= meMomenta()[3];
  _p_boson=  meMomenta()[2];
  _photon = mePartonData()[3];
  _boson= mePartonData()[2];
  // gluon ParticleData object
  _gluon = getParticleData(ParticleID::g);
  //_alphas= SM().alphaS(scale());
  _alphas= SM().alphaS(_p_boson.m2());
  sin2w = SM().sin2ThetaW();
  //alphae = SM().alphaEM(_p_boson.m2());
  //alphae = SM().alphaEM(sHat());
  alphae = SM().alphaEM(scale());
  // BeamParticleData objects for PDF's
  _hadron_A= dynamic_ptr_cast<Ptr<BeamParticleData>::transient_const_pointer>
    (lastParticles().first->dataPtr());  
  //assert(_hadron_A);
  _hadron_B= dynamic_ptr_cast<Ptr<BeamParticleData>::transient_const_pointer>
    (lastParticles().second->dataPtr()); 
  //assert(_hadron_B);

  
  if (_partona->id()<0) {
    swap(_partona,_partonb);
    swap(_p_partona,_p_partonb);
    swap(_hadron_A,_hadron_B);}

  // select quark/anti-quark for quark radiation:
  _iflagcq = 0;
  if(UseRandom::rnd()>0.5) _iflagcq = 1;
  if (_iflagcq==0) {
  _quark = _partonb;
  _quarkpi = _partona->CC();}
  else {
  _quark = _partona;
  _quarkpi = _partonb->CC();}
  

  // calculate CKM matrix square
  int iu,id;
  iu = abs(_partona->id());
  id = abs(_partonb->id());
  if(iu%2!=0) swap(iu,id);
  iu = (iu-2)/2;
  id = (id-1)/2;
  ckm_ = SM().CKM(iu,id);

  _muF2= scale();
  CrossSection lo = MEPP2VGamma::dSigHatDR();
  double weight = NLOweight();
  return lo*weight;
}

// Transformation from Born+rad phase space configuration to m+1 phase space of final states
Lorentz5Momentum MEPP2VGammaPowheg::InvLortr(Lorentz5Momentum kbar_a, Lorentz5Momentum kbar_b, double xi,
					     Lorentz5Momentum k_i, Lorentz5Momentum kbar_j, int fi) const{
  Lorentz5Momentum KK, KKbar;
  Energy2 K2, Kp2, Kdotkj, Kpdotkj;
  KKbar = kbar_a + kbar_b;
  if (fi == 0){
    KK = kbar_a/xi + kbar_b -k_i;}
  else{
    KK = kbar_a + kbar_b/xi -k_i;}
  K2 = KK.m2();
  Kp2 = (KK + KKbar).m2();
  Kdotkj = kbar_j.dot(KKbar);
  Kpdotkj = kbar_j.dot(KK + KKbar);
  return kbar_j -2.0*(Kpdotkj/Kp2)*(KK + KKbar) +2.0*(Kdotkj/K2)*KK;
} 

//Transformation from m+1 phase space to Born phase space configuration of final states
Lorentz5Momentum MEPP2VGammaPowheg::Lortr(Lorentz5Momentum k_a, Lorentz5Momentum k_b, double xi,
					     Lorentz5Momentum k_i, Lorentz5Momentum k_j, int fi) const{
  Lorentz5Momentum KK, KKbar;
  Energy2 K2, Kp2, Kdotkj, Kpdotkj;
  KK = k_a + k_b -k_i;
  if (fi == 0){
    KKbar = xi*k_a + k_b;}
  else{
    KKbar = k_a + xi*k_b;}
  K2 = KK.m2();
  Kp2 = (KK + KKbar).m2();
  Kdotkj = k_j.dot(KK);
  Kpdotkj = k_j.dot(KK + KKbar);
  return k_j -2.0*(Kpdotkj/Kp2)*(KK + KKbar) +2.0*(Kdotkj/K2)*KKbar;
}  

// momentum of radiated gluon transformed from Born phase space
Lorentz5Momentum MEPP2VGammaPowheg::radk(Lorentz5Momentum kbar_a, Lorentz5Momentum kbar_b, double xi,
					 double vi, double phi, int fi) const{
  Energy Epi; 
  double costhei, sinthei, ga;  
  Lorentz5Momentum CMSk, kpi, ki;  
  Epi = sqrt(kbar_a.dot(kbar_b)/(2.0*xi))*(1.0-xi); 
  if (fi == 0){
    costhei = (1.0 -2.0*vi -xi)/(1.0-xi);
    CMSk = kbar_a/xi + kbar_b;}
  else{
    costhei = -(1.0 -2.0*vi -xi)/(1.0-xi);
    CMSk = kbar_a + kbar_b/xi;}
  sinthei = sqrt(1.0 - sqr(costhei));
  kpi = Lorentz5Momentum(Epi*sinthei, ZERO, Epi*costhei, Epi, ZERO);
  ga = CMSk.e()/CMSk.m();
  ki = kpi.boost(0.0, 0.0, tanh(CMSk.rapidity()), ga);
  return ki.rotateZ(phi);
} 

// momentum of radiated quark transformed from gq->qV Born phase space
Lorentz5Momentum MEPP2VGammaPowheg::radkqr(Lorentz5Momentum kbar_ga, Lorentz5Momentum kbar_V, Energy2 MV2, 
					   double zi, double ui, double phi) const{
  Energy Epi; 
  Energy2 kgadotV;
  double costhei, sinthei, zu, phiga, thega;  
  Lorentz5Momentum CMSk, kpi, ki;  
  zu = (1.0 - zi)*(1.0 - ui);
  kgadotV = kbar_ga.dot(kbar_V);
  Epi = (zu + ui)*kgadotV/sqrt(MV2 + kgadotV); 
  costhei = (zu - ui)/(zu + ui)*(MV2 + kgadotV)/(2.0*kgadotV) - MV2/(2.0*kgadotV);
  sinthei = sqrt(1.0 - sqr(costhei));
  kpi = Lorentz5Momentum(Epi*sinthei, ZERO, Epi*costhei, Epi, ZERO);
  kpi = kpi.rotateZ(phi);
  kpi = Lorentz5Momentum(ZERO, ZERO, kbar_ga.e(), kbar_ga.e(), ZERO);
  phiga = kbar_ga.phi();
  thega = kbar_ga.theta();
  ki = kpi.rotateY(thega);
  ki = ki.rotateZ(phiga);
  return ki;
} 

// Born matrix element
double MEPP2VGammaPowheg::MatrBorn(Lorentz5Momentum k_a ,  Lorentz5Momentum k_b,
				   Lorentz5Momentum k_ga,  Lorentz5Momentum k_V,
				   double char0, double char1, Energy2 MV2) const {
  using Constants::pi;
  // kinematic invariants
  Energy2 bs = (k_a + k_b ).m2();
  Energy2 bt = (k_a - k_V ).m2();
  Energy2 bu = (k_a - k_ga).m2();
  double coeff  = sqr(4.0*pi*alphae)* ckm_/sin2w;
  double scfact = 1.0/12.0;
  double matr   = 4.0*sqr(char0*bt + char1*bu)/(bt*bu*sqr(bt+bu))*(bs*MV2 -bt*bu +sqr(bt+bu)/2.0);
  return matr*coeff*scfact;
}


// Tree level matrix element for gq->qV
double MEPP2VGammaPowheg::MatrQV(Lorentz5Momentum p_a, Lorentz5Momentum p_b, Lorentz5Momentum k_q,
			         Lorentz5Momentum k_V, Energy2 MV2) const {
  using Constants::pi;
  Energy2 bs, bt, bu;
  double coeff, scfact, matr;
  Lorentz5Momentum k_a, k_b;
  if (_iflagcq==0) {
    k_a = p_a;
    k_b = p_b;}
  else {
    k_a = p_b;
    k_b = p_a;}

  bs = (k_a + k_b).m2();
  bt = (k_a - k_q).m2();
  bu = (k_a - k_V).m2();

  coeff = sqr(4.0*pi)*alphae*_alphas* ckm_/sin2w;
  // spin = 1/4 colour = 4/(3*8)
  scfact = 1.0/24.0;
  matr = 4.0*(bu*MV2 +sqr(bt) +sqr(bs))/(bt*bs);
  return matr*coeff*scfact;
}



double MEPP2VGammaPowheg::test2to3(Lorentz5Momentum k_a, Lorentz5Momentum k_b, Lorentz5Momentum k_ga,
				    Lorentz5Momentum k_1, Lorentz5Momentum k_2) const{
  using namespace ThePEG::Helicity;
  using namespace ThePEG;
  double sum(0.), scfact;
  Complex diag;
  vector<SpinorWaveFunction> qin,qbarout;
  vector<SpinorBarWaveFunction> qbarin,qout;
  vector<VectorWaveFunction> pout;
  SpinorWaveFunction q_in(k_a, _partona, incoming);
  SpinorBarWaveFunction qbar_in(k_b, _partonb, incoming);
  SpinorBarWaveFunction q_out(k_1, _partona, outgoing);
  SpinorWaveFunction qbar_out(k_2, _partonb, outgoing);
  VectorWaveFunction p_out(k_ga, _photon, outgoing);
  for(unsigned int ix=0;ix<2;++ix) {
      q_in.reset(ix);
      qin.push_back(q_in);
      qbar_in.reset(ix);
      qbarin.push_back(qbar_in);
      q_out.reset(ix);
      qout.push_back(q_out);
      qbar_out.reset(ix);
      qbarout.push_back(qbar_out);
      p_out.reset(2*ix);
      pout.push_back(p_out);
  }
  Energy2 scaleEM;
  scaleEM = scale();
  for(unsigned int ihel1=0;ihel1<2;++ihel1) {
    for(unsigned int ihel2=0;ihel2<2;++ihel2) {
      for(unsigned int ihela=0;ihela<2;++ihela) {
	for(unsigned int phel=0;phel<2;++phel) {
	  for(unsigned int ihelb=0;ihelb<2;++ihelb) {
             SpinorWaveFunction inters1 =
	       FFPvertex_->evaluate(scaleEM,5,_partona,qin[ihela],pout[phel]);
             VectorWaveFunction intersw =
	       FFWvertex_->evaluate(scaleEM,3,_boson->CC(),qbarout[ihel2],qout[ihel1]);
             diag = FFWvertex_->evaluate(scaleEM,inters1,qbarin[ihelb],intersw);
             sum += norm(diag);
	  }
	}
      }
    }
  }
  // final spin and colour factors
  // spin = 1/4 colour = 3/9
  scfact = 1.0/12.0;  
  return sum*scfact;
}
double MEPP2VGammaPowheg::test2to3am(Lorentz5Momentum k_a, Lorentz5Momentum k_b, Lorentz5Momentum k_ga,
				     Lorentz5Momentum k_1, Lorentz5Momentum k_2, Energy2 MV2) const{
  using Constants::pi;
  double scfact, coupfact, part;
  coupfact = sqr(4.0*pi*alphae)*(4.0*pi*alphae)* sqr(ckm_/sin2w)/4.0;
  part = (k_ga.dot(k_1)*k_b.dot(k_2) +k_ga.dot(k_2)*k_b.dot(k_1))/sqr((k_1+k_2).m2()-MV2);
   // final spin and colour factors
  // spin = 1/4 colour = 3/9
  scfact = 1.0/12.0;  
  return part/(2.0*k_a.dot(k_ga))*coupfact*scfact*sqr(GeV);
}

// Gluon Real radiation matrix element
InvEnergy2 MEPP2VGammaPowheg::MatrRealGluon(Lorentz5Momentum k_a, Lorentz5Momentum k_b,
					    Lorentz5Momentum k_ga,
					    Lorentz5Momentum k_V, Lorentz5Momentum k_glu) const{
  using namespace ThePEG::Helicity;
  using namespace ThePEG;
  double sum(0.), scfact;
  Energy2 bosonMass2_;
  vector<SpinorWaveFunction> qin;
  vector<SpinorBarWaveFunction> qbarin;
  vector<VectorWaveFunction> wout,pout,gout;
  SpinorWaveFunction    q_in   (k_a  , _partona, incoming);
  SpinorBarWaveFunction qbar_in(k_b  , _partonb, incoming);
  VectorWaveFunction    v_out  (k_V  , _boson  , outgoing);
  VectorWaveFunction    p_out  (k_ga , _photon , outgoing);
  VectorWaveFunction    g_out  (k_glu, _gluon  , outgoing);
  bosonMass2_ = k_V.m2();
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix<2) {
      q_in.reset(ix);
      qin.push_back(q_in);
      qbar_in.reset(ix);
      qbarin.push_back(qbar_in);
      g_out.reset(2*ix);
      //g_out.reset(10);
      gout.push_back(g_out);
      p_out.reset(2*ix);
      //p_out.reset(10);
      pout.push_back(p_out);
    }
    v_out.reset(ix);
    wout.push_back(v_out);
  }
  Energy2 scaleEM, scaleS;
  scaleS  = bosonMass2_;
  scaleEM = scale();
  vector<Complex> diag(_boson->id()==ParticleID::Z0 ? 6 : 8, 0.);
  AbstractFFVVertexPtr vertex = 
    _boson->id()==ParticleID::Z0 ? FFZvertex_ : FFWvertex_;
  for(unsigned int ihel1=0;ihel1<2;++ihel1) {
    for(unsigned int ihel2=0;ihel2<2;++ihel2) {
      for(unsigned int whel=0;whel<3;++whel) {
	for(unsigned int phel=0;phel<2;++phel) {
	  for(unsigned int ghel=0;ghel<2;++ghel) {
	    //
	    //  Diagrams which are the same for W/Z gamma
	    //
	    // first diagram
	    SpinorWaveFunction inters1 = 
	      FFPvertex_->evaluate(scaleEM,5,_partona,qin[ihel1],pout[phel]);
	    SpinorBarWaveFunction inters2 = 
	      vertex->evaluate(scaleEM,5,_partona->CC(),qbarin[ihel2],wout[whel]);
	    diag[0] = FFGvertex_->evaluate(scaleS,inters1,inters2,gout[ghel]);
	    // second diagram
	    SpinorWaveFunction inters3 = 
	      FFGvertex_->evaluate(scaleS,5,_partona,qin[ihel1],gout[ghel]);
	    SpinorBarWaveFunction inters4 = 
	      FFPvertex_->evaluate(scaleEM,5,_partonb,qbarin[ihel2],pout[phel]);
	    diag[1] = vertex->evaluate(scaleEM,inters3,inters4,wout[whel]);
	    // fourth diagram
	    diag[2] = FFPvertex_->evaluate(scaleEM,inters3,inters2,pout[phel]);
	    // fifth diagram
	    SpinorBarWaveFunction inters5 = 
	      FFGvertex_->evaluate(scaleS,5,_partonb,qbarin[ihel2],gout[ghel]);
	    diag[3] = 
	      vertex->evaluate(scaleEM,inters1,inters5,wout[whel]);
	    // sixth diagram
	    SpinorWaveFunction inters6 = 
	      vertex->evaluate(scaleEM,5,_partonb->CC(),qin[ihel1],wout[whel]);
	    diag[4] = FFGvertex_->evaluate(scaleS,inters6,inters4,gout[ghel]);
	    // eighth diagram
	    diag[5] = FFPvertex_->evaluate(scaleEM,inters6,inters5,pout[phel]);
	    //
	    //  Diagrams only for W gamma
	    //
	    if(_boson->id()!=ParticleID::Z0) {
	      // third diagram
	      VectorWaveFunction interv = 
		WWWvertex_->evaluate(scaleEM,3,_boson->CC(),pout[phel],wout[whel]);

	      diag[6] = vertex->evaluate(scaleEM,inters3,qbarin[ihel2],interv);
	      // seventh diagram
	      diag[7] = vertex->evaluate(scaleEM,qin[ihel1],inters5,interv);
	    }
	    // sum
	    Complex dsum = std::accumulate(diag.begin(),diag.end(),Complex(0.));
// 	    cerr << "testing" 
// 		 << ihel1 << " " << ihel2 << " " 
// 		 << ghel  << " " << phel  << " " << whel << " " 
// 		 << dsum/sqrt(norm(diag[0])+norm(diag[1])+norm(diag[2])+norm(diag[3])+
// 			      norm(diag[4])+norm(diag[5])+norm(diag[6])+norm(diag[7]))<< "\n";
	    sum += norm(dsum);
	  }
	}
      }
    }
  }
  // final spin and colour factors
  // spin = 1/4 colour = 4/9
  scfact = 1./9.;  
  return sum*scfact*UnitRemoval::InvE2;
}


// Quark Real radiation matrix element
double MEPP2VGammaPowheg::MatrRealQuark(Lorentz5Momentum p_a, Lorentz5Momentum p_b, Lorentz5Momentum k_ga,
					Lorentz5Momentum k_V, Lorentz5Momentum k_q) const{
  using namespace ThePEG::Helicity;
  using namespace ThePEG;
  double sum(0.), scfact;
  Energy bosonMass_;
  vector<SpinorWaveFunction> qin;
  vector<SpinorBarWaveFunction> qbarin;
  vector<SpinorBarWaveFunction> qout;
  vector<SpinorWaveFunction> qbarout;
  vector<VectorWaveFunction> vout,pout,gin;
  Lorentz5Momentum k_a, k_b;
  if (_iflagcq==0) {
    k_a = p_a;
    k_b = p_b;}
  else {
    k_a = p_b;
    k_b = p_a;}
  SpinorWaveFunction q_in(k_b, _quark, incoming);
  SpinorBarWaveFunction qbar_in(k_b, _quark, incoming);
  SpinorBarWaveFunction q_out(k_q, _quarkpi, outgoing); //modify
  SpinorWaveFunction qbar_out(k_q, _quarkpi, outgoing); //modify
  VectorWaveFunction v_out(k_V  ,_boson, outgoing);
  VectorWaveFunction p_out(k_ga, _photon, outgoing);
  VectorWaveFunction g_in(k_a, _gluon, incoming);

  bosonMass_ = k_V.mass();
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix<2) {
      q_in.reset(ix);
      qin.push_back(q_in);
      qbar_in.reset(ix);
      qbarin.push_back(qbar_in);
      q_out.reset(ix);
      qout.push_back(q_out);
      qbar_out.reset(ix);
      qbarout.push_back(qbar_out);
      //g_in.reset(10);
      g_in.reset(2*ix);
      gin.push_back(g_in);
      //p_out.reset(10);
      p_out.reset(2*ix);
      pout.push_back(p_out);
    }
    v_out.reset(ix);
    vout.push_back(v_out);
  }
  Energy2 scaleEM, scaleS;
  scaleS = sqr(bosonMass_);
  scaleEM = scale();

  vector<Complex> diag(_boson->id()==ParticleID::Z0 ? 6 : 8, 0.);
  AbstractFFVVertexPtr vertex = 
    _boson->id()==ParticleID::Z0 ? FFZvertex_ : FFWvertex_;
  for(unsigned int ihel1=0;ihel1<2;++ihel1) { //incoming quark
    for(unsigned int ihel2=0;ihel2<2;++ihel2) { //outgoing quark
      for(unsigned int vhel=0;vhel<3;++vhel) {
	for(unsigned int phel=0;phel<2;++phel) {
	  for(unsigned int ghel=0;ghel<2;++ghel) {
	    //
         // suppose parton_b is quark:
         if (_quark->id()>0){
	    //  Diagrams which are the same for W/Z gamma
	    //
	    // first diagram
	    SpinorWaveFunction inters1 = 
	      FFPvertex_->evaluate(scaleEM,5,_quark,qin[ihel1],pout[phel]);
	    SpinorBarWaveFunction inters2 = 
	      FFGvertex_->evaluate(scaleS,5,_quarkpi->CC(),qout[ihel2],gin[ghel]);
	    diag[0] = vertex->evaluate(scaleEM,inters1,inters2,vout[vhel]);
	    // second diagram
	    SpinorWaveFunction inters3 = 
	      vertex->evaluate(scaleEM,5,_quarkpi,qin[ihel1],vout[vhel]);
	    SpinorBarWaveFunction inters4 = 
	      FFPvertex_->evaluate(scaleEM,5,_quarkpi->CC(),qout[ihel2],pout[phel]);
	    diag[1] = FFGvertex_->evaluate(scaleS,inters3,inters4,gin[ghel]);
	    // third diagram
	    diag[2] = FFPvertex_->evaluate(scaleEM,inters3,inters2,pout[phel]);
	    // fourth diagram
	    SpinorWaveFunction inters5 = 
	      FFGvertex_->evaluate(scaleS,5,_quark,qin[ihel1],gin[ghel]);
	    diag[3] = 
	      vertex->evaluate(scaleEM,inters5,inters4,vout[vhel]);
	    // fifth diagram
	    SpinorBarWaveFunction inters6 = 
	      vertex->evaluate(scaleEM,5,_quark->CC(),qout[ihel2],vout[vhel]);
	    diag[4] = FFGvertex_->evaluate(scaleS,inters1,inters6,gin[ghel]);
	    // sixth diagram
	    diag[5] = FFPvertex_->evaluate(scaleEM,inters5,inters6,pout[phel]);
	    //
	    //  Diagrams only for W gamma
	    //
	    if(_boson->id()!=ParticleID::Z0) {
	      // seventh diagram
	      VectorWaveFunction interv = 
		WWWvertex_->evaluate(scaleEM,3,_boson->CC(),pout[phel],vout[vhel]);
	      diag[6] = vertex->evaluate(scaleEM,inters5,qout[ihel2],interv);
	      // eighth diagram
	      diag[7] = vertex->evaluate(scaleEM,qin[ihel1],inters2,interv);
	    }
	  }
	 //
      //suppose parton_b is anti-quark:
	  else {
	    //  Diagrams which are the same for W/Z gamma
	    //
	    // first diagram
	    SpinorBarWaveFunction inters1 = 
	      FFPvertex_->evaluate(scaleEM,5,_quark,qbarin[ihel1],pout[phel]);
	    SpinorWaveFunction inters2 = 
	      FFGvertex_->evaluate(scaleS,5,_quarkpi->CC(),qbarout[ihel2],gin[ghel]);
	    diag[0] = vertex->evaluate(scaleEM,inters2,inters1,vout[vhel]);
	    // second diagram
	    SpinorBarWaveFunction inters3 = 
	      vertex->evaluate(scaleEM,5,_quarkpi,qbarin[ihel1],vout[vhel]);
	    SpinorWaveFunction inters4 = 
	      FFPvertex_->evaluate(scaleEM,5,_quarkpi->CC(),qbarout[ihel2],pout[phel]);
	    diag[1] = FFGvertex_->evaluate(scaleS,inters4,inters3,gin[ghel]);
	    // third diagram
	    diag[2] = FFPvertex_->evaluate(scaleEM,inters2,inters3,pout[phel]);
	    // fourth diagram
	    SpinorBarWaveFunction inters5 = 
	      FFGvertex_->evaluate(scaleS,5,_quark,qbarin[ihel1],gin[ghel]);
	    diag[3] = 
	      vertex->evaluate(scaleEM,inters4,inters5,vout[vhel]);
	    // fifth diagram
	    SpinorWaveFunction inters6 = 
	      vertex->evaluate(scaleEM,5,_quark->CC(),qbarout[ihel2],vout[vhel]);
	    diag[4] = FFGvertex_->evaluate(scaleS,inters6,inters1,gin[ghel]);
	    // sixth diagram
	    diag[5] = FFPvertex_->evaluate(scaleEM,inters6,inters5,pout[phel]);
	    //
	    //  Diagrams only for W gamma
	    //
	    if(_boson->id()!=ParticleID::Z0) {
	      // seventh diagram
	      VectorWaveFunction interv = 
		WWWvertex_->evaluate(scaleEM,3,_boson->CC(),pout[phel],vout[vhel]);
	      diag[6] = vertex->evaluate(scaleEM,qbarout[ihel2],inters5,interv);
	      // eighth diagram
	      diag[7] = vertex->evaluate(scaleEM,inters2,qbarin[ihel1],interv);
	    }
	  }

	    // sum
	    Complex dsum = std::accumulate(diag.begin(),diag.end(),Complex(0.));
	    /*
	    cerr << "testing" 
		 << ihel1 << " " << ihel2 << " " 
		 << ghel  << " " << phel  << " " << whel << " " 
		 << dsum/sqrt(norm(diag[0])+norm(diag[1])+norm(diag[2])+norm(diag[3])+
			      norm(diag[4])+norm(diag[5])+norm(diag[6])+norm(diag[7]))<< "\n";
	    */
	    sum += norm(dsum);
	  }
	}
      }
    }
  }
  // final spin and colour factors
  // spin = 1/4 colour = 4/(3*8)
  scfact = 1.0/24.0;  
  return sum*scfact;
}

// dipole of gluon radiation
InvEnergy2 MEPP2VGammaPowheg::Dipole(Lorentz5Momentum k_a, Lorentz5Momentum k_b,
				     Lorentz5Momentum k_ga,
				     Lorentz5Momentum k_V, double char0,double char1,
				     Energy2 MV2,Energy2 kig,double xi) const {
  using Constants::pi;
  double coeff = 8.*pi*_alphas*CF_;
  InvEnergy2 fact = 1./(2.0*xi *kig) * ( 1. + sqr(xi) )/(1. - xi);
  return MatrBorn(k_a,k_b,k_ga,k_V,char0,char1,MV2) *coeff*fact;
}

// parton distribution function of beam partID
double MEPP2VGammaPowheg::PDFab(double z, int partID)  const {
  if(partID==0){
    //assert(_hadron_A);
    return _hadron_A->pdf()->xfx(_hadron_A,_partona,_muF2,z);}  
  else{
    //assert(_hadron_B);
    return _hadron_B->pdf()->xfx(_hadron_B,_partonb,_muF2,z);}
}

// parton distribution function ratio
double MEPP2VGammaPowheg::PDFratio(double z, double x, int partID)  const {
  double pdfz, pdfx;
  if (partID <3){
  pdfz = PDFab(z, partID);
  pdfx = PDFab(x, partID);}

  if (partID ==3){
    if (_iflagcq ==0){
      pdfx = PDFab(z, 0);
      pdfz = _hadron_A->pdf()->xfx(_hadron_A,_gluon,_muF2,z);}
    else {
      pdfx = PDFab(x, 1);
      pdfz = _hadron_A->pdf()->xfx(_hadron_B,_gluon,_muF2,x);}}

  if (pdfx == 0.0) {
    if (pdfz == 0.0){
      return 1.0;}
    else{
      return 0.0;}}
  else{
    return pdfz/pdfx;}

}


// functionals of PDFs and energy fractions in the collinear remnants:

double MEPP2VGammaPowheg::KP1(double zz, Energy2 shat, Energy2 muF2, int partID) const{
  double x, hpart;
  x = 1.0 -_xtil +zz*_xtil;
  //pdfzx= PDFab(zz/x, partID);
  //pdfz=  PDFab(zz, partID);
  hpart= log(sqr(1.0-x)*shat/(x*muF2))*
    (2.0/(1.0-x) -(1.0+sqr(x))/(1.0-x)*PDFratio(zz/x,zz,partID)/x) +2.0/(1.0-x)*log(x);
  return (1.0-zz)*hpart;
}

double MEPP2VGammaPowheg::KP2(double zz, Energy2 shat, Energy2 muF2) const{
  double x, hpart;
  x = zz*_xtil;
  hpart= log(sqr(1.0-x)*shat/muF2)*(2.0/(1.0-x));
  return zz*hpart;
}

double MEPP2VGammaPowheg::KP3(double zz,  int partID) const{
  double x, hpart;
  x = 1.0 -_xtil +zz*_xtil;
  //pdfzx= PDFab(zz/x, partID);
  //pdfz=  PDFab(zz, partID);
  hpart= (1.0-x)*PDFratio(zz/x,zz,partID)/x;
  return (1.0-zz)*hpart;
}

// definitions of H, F^V:
double MEPP2VGammaPowheg::Hfunc(Energy2 t, Energy2 s, Energy2 m2) const{  
  using Constants::pi;
  return sqr(pi)-sqr(log(s/m2))+sqr(log(-t/s))-sqr(log(-t/m2))
    -2.0*ReLi2(1.0-s/m2)-2.0*ReLi2(1.0-t/m2);
} 


double MEPP2VGammaPowheg::FWfunc(Energy2 t, Energy2 u, Energy2 s, Energy2 m2) const {
  using Constants::pi;
    double y1,y2,y3,y4,y5;
    y1= 4.0*(2.0*s*s/(t*u) +2.0*s/u +t/u)* Hfunc(u,s,m2);
    y2= -8.0/3.0*sqr(pi)*s*(2.0*s/t +t/s -u/s)/(t+u);
    y3= 4.0*(6.0-10.0*u/(t+u)-10.0*s*s/(u*(t+u))-11.0*s/u-5.0*t/u+2.0*s/(t+u)+s/(s+t));
    y4= -4.0*log(s/m2)*(3.0*t/u+2.0*s/u+4.0*s*(t+s)/(u*(t+u))+2.0*t/u*sqr(s/(t+u)));
    y5= 4.0*log(-u/m2)*((4.0*s+u)/(s+t)+s*u/sqr(s+t));
    return y1+y2+y3+y4+y5;
}

double MEPP2VGammaPowheg::FZfunc(Energy2 t, Energy2 u, Energy2 s, Energy2 m2) const{
  using Constants::pi;
    double y1,y2,y3,y4,y5;
    y1= 4.0*(2.0*s*s/(t*u) +2.0*s/u +t/u)* Hfunc(u,s,m2);
    y2= -8.0/3.0*sqr(pi)*s*s/(t*u);
    y3= 4.0*(1.0-5.0*s*s/(t*u)-11.0*s/u-5.0*t/u+2.0*s/(t+u)+s/(s+t));
    y4= 4.0*log(s/m2)*(4.0*s/(t+u)+2.0*sqr(s/(t+u))-3.0*sqr(s+t)/(t*u));
    y5= 4.0*log(-u/m2)*((4.0*s+u)/(s+t)+s*u/sqr(s+t));
    return y1+y2+y3+y4+y5;
}

double MEPP2VGammaPowheg::FVfunc(Energy2 t, Energy2 u, Energy2 s, Energy2 m2) const {
  if(mePartonData()[2]->id()==ParticleID::Z0) {
    return FZfunc(t,u,s,m2);
  }
  else {
    return FWfunc(t,u,s,m2);
  }
}

// ratio of NLO/LO 
double MEPP2VGammaPowheg::NLOweight() const {
  using Constants::pi;
  //  double FV, FW, FZ;
  // If only leading order is required return 1:
  if(_contrib==0) return 1.;
  // kinematic invariants
  _ss= sHat();
  _tt= ( _p_partona - _p_boson ).m2();
  _uu= ( _p_partona - _p_photon).m2();
  Energy2 MV2= _p_boson.m2();
  // strong coupling piece
  double alfsfact = _alphas*CF_/(2.0*pi);
  // ids of the incoming partons
  int ida = _partona->id();
  int idb = _partonb->id();
  // charges of the quarks for the incoming partons
  double charge0 = _partona->iCharge()/3.*ida/abs(ida);
  double charge1 = _partonb->iCharge()/3.*idb/abs(idb);
  // prefactor for the born matrix element
  InvEnergy4 smbfact= 4.0 * sqr( charge0*_tt + charge1*_uu )/(_tt*_uu*sqr(_tt+_uu));
  // eps^0 term in the expansion of the born ME
  Energy4 smborn0 =  (_ss*MV2 -_tt*_uu + 0.5*sqr(_tt+_uu));
  // eps^1 term in the expansion of the born ME
  Energy4 smborn1 = -(_ss*MV2 -_tt*_uu +sqr(_tt+_uu));
  // eps^2 term in the expansion of the born ME
  Energy4 smborn2 = 0.5*sqr(_tt+_uu);
  // divergent Virtual - I(eps) piece
  double partep= 2.0*(5.0-sqr(pi)/3.0) +3.0*smborn1/smborn0 +2.0*smborn2/smborn0;
  // finite virtual piece
  double smloop = 0.5*(charge0*_tt+ charge1*_uu)* 
    ( charge0 * FVfunc(_tt,_uu,_ss,MV2) + 
      charge1 * FVfunc(_uu,_tt,_ss,MV2) )/(_tt+_uu);
  double partloop = smloop/(smborn0*smbfact);
  // PR checked above here 26/8/09
  // collinear remnant (needs checking)
  double partc = -(KP1(_xa,_ss,_muF2,0)+KP2(_xa,_ss,_muF2)-KP3(_xa,0))
    -(KP1(_xb,_ss,_muF2,1)+KP2(_xb,_ss,_muF2)-KP3(_xb,1))-2.0*(5.0-sqr(pi)/3.0);
  // real radiation
  // azimuthal angle
  double phi = Constants::twopi*_y3;
  //  for singularities with parton_a:
  double xiab = 1.0 - _y1*(1.-_xa);
  double vi   = (1.0 - xiab )*_y2;
  Lorentz5Momentum k_a   = _p_partona/xiab;
  Lorentz5Momentum k_b   = _p_partonb;
  Lorentz5Momentum k_glu =  radk(_p_partona, _p_partonb, xiab, vi, phi, 0); 
  Lorentz5Momentum k_ga  = InvLortr(_p_partona, _p_partonb, xiab, k_glu, _p_photon, 0);
  Lorentz5Momentum k_V   = InvLortr(_p_partona, _p_partonb, xiab, k_glu, _p_boson , 0);
  // first dipole term for emission from emitter
  Energy2 kadotg = k_a.dot(k_glu);
  InvEnergy2 dipole1 = Dipole(_p_partona, _p_partonb, _p_photon, _p_boson, 
			      charge0, charge1, MV2, kadotg, xiab);
  // second dipole term for emission from spectator
  Energy2 kbdotg = k_b.dot(k_glu);
  Lorentz5Momentum kbarp_ga = Lortr(k_a, k_b, xiab, k_glu, k_ga, 1);
  Lorentz5Momentum kbarp_V  = Lortr(k_a, k_b, xiab, k_glu, k_V , 1);
  InvEnergy2 dipole2 = Dipole(k_a, xiab*k_b, kbarp_ga, kbarp_V,
			      charge0, charge1, MV2, kbdotg, xiab);
  // dipole subtracted sum
  InvEnergy2 sumdipole = abs(dipole1) + abs(dipole2);
  double partreal = MatrRealGluon(k_a, k_b, k_ga, k_V, k_glu)/sumdipole 
    - dipole1/abs(dipole1);
  Energy2 jacobfact = 2.0*k_a.dot(k_b)*(1.0 - xiab)*(1.0 -_xa)/(16.0*sqr(pi)*xiab);
  double  pdffact = PDFratio(_xa/xiab, _xa, 0);
  double fluxfact = _ss/(k_a+k_b).m2();
  double smreala = partreal * abs(dipole1) * jacobfact * pdffact * fluxfact; 
  // for singularities with parton_b:
  xiab = 1.0 -_y1*(1.-_xb);
  vi = (1.0 -xiab)*_y2;
  k_a = _p_partona;
  k_b = _p_partonb/xiab;
  k_glu =  radk(_p_partona, _p_partonb, xiab, vi, phi, 1); 
  k_ga = InvLortr(_p_partona, _p_partonb, xiab, k_glu, _p_photon, 1);
  k_V = InvLortr(_p_partona, _p_partonb, xiab, k_glu, _p_boson, 1);
  kadotg = k_a.dot(k_glu);
  kbdotg = k_b.dot(k_glu);
  // first dipole term for emission from emitter
  dipole1 = Dipole(_p_partona, _p_partonb, _p_photon, _p_boson,
		   charge0, charge1, MV2, kbdotg, xiab);
  // second dipole term for emission from the spectator
  kbarp_ga = Lortr(k_a, k_b, xiab, k_glu, k_ga, 0);
  kbarp_V  = Lortr(k_a, k_b, xiab, k_glu, k_V, 0);
  dipole2 = Dipole(xiab*k_a, k_b, kbarp_ga, kbarp_V,
		   charge0, charge1, MV2, kadotg, xiab);
  // dipole subtracted sum
  sumdipole = abs(dipole1) + abs(dipole2);
  partreal = MatrRealGluon(k_a, k_b, k_ga, k_V, k_glu)/sumdipole 
    - dipole1/abs(dipole1);
  jacobfact = 2.0*k_a.dot(k_b)*(1.0 -xiab)*(1.0 -_xb)/(16.0*sqr(pi)*xiab);
  pdffact = PDFratio(_xb/xiab, _xb, 1);
  fluxfact = _ss/(k_a+k_b).m2();
  double smrealb = partreal * abs(dipole1) * jacobfact * pdffact * fluxfact;
  // sum over two piece of real parts
  double smreal = (smreala + smrealb)/MatrBorn(_p_partona, _p_partonb, _p_photon,
					       _p_boson, charge0, charge1, MV2);



  /*
  // d\sigma^A+d\sigma^C quark radiation:
  zz = _y1;
  pn = _p_partona+ _p_partonb -zz*_p_photon;
  smbornqr = MatrQV(_p_partona, _p_partonb, zz*_p_photon, _p_boson);
  pz = (1.0+sqr(1.0-zz))/zz;
  hfs = ;
  chargeqr= ;
  lnz = log(sqr(2.0*(1.0-zz)*zz*_p_photon.dot(pn))/(_muF2*pn.m2()));
  partqr =(sqr(chargeqr)*(pz*lnz +zz) -hfs)/sqr(zz);
  pdffactqr = PDFratio(_xa, _xb, 3);
  factqr = alphae/(2.0*pi);
  smqr = smbornqr*partqr*factqr*pdffactqr/MatrBorn(_p_partona, _p_partonb, _p_photon, _p_boson, charge0, charge1, MV2);
  */

//   Lorentz5Momentum k_1, k_2, kquark;
//   Energy E2, EV, PV;
//   double testmn, testma;
//   EV = _p_boson.e();
//   PV = sqrt(sqr(EV)-MV2);
//   E2 = (EV + PV);
//   k_1 = Lorentz5Momentum((E2/PV)*_p_boson.vect(), E2);
//   k_2 = _p_boson - k_1;
//   testmn = test2to3(_p_partona, _p_partonb, k_1, _p_photon, k_2);
//   testma = test2to3am(_p_partona, _p_partonb, k_1, _p_photon, k_2, MV2);
//   kquark = radkqr(_p_photon, _p_boson, MV2, 0.5,0.5,1.);

  //double borm1, borm2;
  //borm2= MEPP2VGamma::me2p();
  //borm1= MEPP2VGamma::me2();
  //Lorentz5Momentum loparta, lopartb, lophoton, loboson;
  //loparta= MEPP2VGamma::mompartona();
  //lopartb= MEPP2VGamma::mompartonb();
  //lophoton= MEPP2VGamma::momphoton();
  //loboson =  MEPP2VGamma::momboson();
  //return 1.+ alfsfact*(partep+partloop+partc)+ smreal;
  //return 1.+ alfsfact*(partep+partloop+partc);
  return 1.0+ smreal;
}
