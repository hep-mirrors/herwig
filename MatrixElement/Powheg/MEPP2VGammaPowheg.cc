// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2VGammaPowheg class.
//
 
#include "MEPP2VGammaPowheg.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h" // from VGammaHardGenerator
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Shower/Base/Evolver.h" // Should have it? from VGammaHardGenerator
#include "Herwig++/Shower/Base/KinematicsReconstructor.h" //?
#include "Herwig++/Shower/Base/PartnerFinder.h" //?
#include "ThePEG/PDT/StandardMatchers.h" //?
#include "ThePEG/Repository/EventGenerator.h" //?
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
//#include "ThePEG/Utilities/SimplePhaseSpace.h"  // Should we include it like in MEPP2VGammaNewPowheg?
#include "ThePEG/Cuts/Cuts.h"
#include "Herwig++/Utilities/GSLIntegrator.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "Herwig++/Utilities/Maths.h"
#include <numeric>

using namespace Herwig;
using Herwig::Math::ReLi2;

  //merge
MEPP2VGammaPowheg::MEPP2VGammaPowheg() 
  :  _contrib(1), _scaleopt(0),
     _fixedScale(100.*GeV), _scaleFact(1.),
     //pTmin_(2.*GeV), qqgFactor_(1.), qgFactor_g(2.),
     //qgFactor_p(0.1), gqbarFactor_g(2.), gqbarFactor_p(0.1),  power_(2.), power_photon(2.5)
     pTmin_(2.*GeV), qqgFactor_(800.), qgFactor_g(800.),
     qgFactor_p(2.3), gqbarFactor_g(800.), gqbarFactor_p(2.3),  power_(1.6), power_photon(2.4)
{} 

  //merge
void MEPP2VGammaPowheg::doinit() {
  // gluon ParticleData object
  _gluon = getParticleData(ParticleID::g); // for HardestEmission ?
  // colour factors
  //  _CF = 4./3.; 
  // _TR = 0.5;
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

//merge
void MEPP2VGammaPowheg::persistentOutput(PersistentOStream & os) const {
  os << _contrib << _scaleopt << ounit(_fixedScale,GeV) << _scaleFact  // ounit exists in both classes!
     << alphaS_ << _gluon <<alphaQCD_<< ounit(pTmin_,GeV)<< fraccut
     << FFPvertex_ << FFWvertex_ << FFZvertex_ << WWWvertex_ << FFGvertex_;
}

//merge
void MEPP2VGammaPowheg::persistentInput(PersistentIStream & is, int) {
  is >> _contrib >> _scaleopt >> iunit(_fixedScale,GeV) >> _scaleFact  // iunit exists in both classes!
     >> alphaS_ >> _gluon >> alphaQCD_>> iunit(pTmin_,GeV) >>fraccut
     >> FFPvertex_ >> FFWvertex_ >> FFZvertex_ >> WWWvertex_ >> FFGvertex_;
}

ClassDescription<MEPP2VGammaPowheg> MEPP2VGammaPowheg::initMEPP2VGammaPowheg;
// Definition of the static class description member.

  //merge
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

  static Parameter<MEPP2VGammaPowheg,double> interfaceAlphaQCD
    ("AlphaQCD",
     "The value of alphaS to use if using a fixed alphaS",
     &MEPP2VGammaPowheg::alphaQCD_, 0.118, 0.0, 0.2,
     false, false, Interface::limited);

  static Reference<MEPP2VGammaPowheg,ShowerAlpha> interfaceShowerAlphaS
    ("ShowerAlphaS",
     "Reference to the object calculating the QCD coupling for the shower",
     &MEPP2VGammaPowheg::alphaS_, false, false, true, false, false);

  static Parameter<MEPP2VGammaPowheg,double> interfaceFraccut
    ("Fraccut",
     "The photon energy fraction cut in the photon cone",
     &MEPP2VGammaPowheg::fraccut, 0.4, 0.0, 1.0,
     false, false, Interface::limited);
}

Energy2 MEPP2VGammaPowheg::scale() const {
  //return _scaleopt == 0 ? sqr(_fixedScale) : _scaleFact*sHat();
  return sHat();
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
  
   Pi= Constants::pi;
  _CF= 4.0/3.0;
  _TR= 1.0/2.0;
  _xa=  lastX1();
  _xb=  lastX2();
  gamnum = 0.5772;

  //incoming parton ParticleData object and momenta
  _partona= mePartonData()[0];
  _partonb= mePartonData()[1];
  _p_partona= meMomenta()[0];
  _p_partonb= meMomenta()[1];

  // gluon ParticleData object
  _gluon = getParticleData(ParticleID::g);

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
    swap(_hadron_A,_hadron_B);
    swap(_xa, _xb);}

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
  ckm = SM().CKM(iu,id);
  // the third componet of isospin of the quark when vector boson is Z0
  I3quark = -double(abs(_partona->id())%2) + 0.5;
//  ParticleData object and momenta of photon and vector boson
  if(mePartonData()[3]->id()==ParticleID::gamma){
      _idboson= 2;
      _p_photon= meMomenta()[3];
      _p_boson=  meMomenta()[2];
      _photon = mePartonData()[3];
      _boson= mePartonData()[2];}
  else{
      _idboson= 3;
      _p_photon= meMomenta()[2];
      _p_boson=  meMomenta()[3];
      _photon = mePartonData()[2];
      _boson= mePartonData()[3];}

  _alphas= SM().alphaS(_p_boson.m2());
  sin2w = SM().sin2ThetaW();
  alphae = SM().alphaEM(scale());
  _muF2= scale();
  double weight = NLOweight();
  return MEPP2VGamma::dSigHatDR()*weight;
}

unsigned int MEPP2VGammaPowheg::orderInAlphaS() const {
  return 0;
}

unsigned int MEPP2VGammaPowheg::orderInAlphaEW() const {
  return 2;
}


// Transformation from Born + CS phase space configuration to m+1 phase space of final states
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
 
// the kinematics of quark radiation from initial gluon is just like the gluon radiation,
// so we can use the functions radk, InvLortr and Lortr.

// momentum of radiated quark transformed from gq->qV Born phase space
Lorentz5Momentum MEPP2VGammaPowheg::radkqr(Lorentz5Momentum kbar_ga, Lorentz5Momentum kbar_V, Energy2 MassV2, 
					   double zi, double ui, double phi, double zcut) const{
  Energy Epi; 
  Energy2 kgadotV;
  double costhei, sinthei, zu, phiga, thega, beta;  
  Lorentz5Momentum CMSk, kpi, ki, kga;  
  zu = (1.0 - zi)*(1.0 - ui);
  if (zi >= zcut){
  kga = kbar_ga*zi;
  //  kga = kbar_ga;
  kgadotV = kga.dot(kbar_V);
  //  Epi = (zu + ui)*kgadotV/sqrt(MV2 + 2.0*kgadotV); 
  Epi = (zu + ui*zi)*kgadotV/sqrt(MassV2 + 2.0*kgadotV)/zi;
  //  costhei = (zu - ui)/(zu + ui)*(MV2 + 2.0*kgadotV)/(2.0*kgadotV) - MV2/(2.0*kgadotV);
  //  costhei = ((zu - ui)*kgadotV -ui*MV2)/((zu + ui)*kgadotV);
  costhei = ((zu - ui*zi)*kgadotV -ui*zi*MassV2)/((zu + ui*zi)*kgadotV); }
  else {
    kga = kbar_ga;
    kgadotV = kga.dot(kbar_V);
    Epi = (zu + ui)*kgadotV/sqrt(MassV2 + 2.0*kgadotV);
    costhei = ((zu - ui)*kgadotV -ui*MassV2)/((zu + ui)*kgadotV);}
  sinthei = sqrt(1.0 - sqr(costhei));
  kpi = Lorentz5Momentum(Epi*sinthei, ZERO, Epi*costhei, Epi, ZERO);
  kpi = kpi.rotateZ(phi);
  beta = -(1.0-zi)*kgadotV/((1.0+zi)*kgadotV + zi*MassV2);
  phiga = kbar_ga.phi();
  thega = kbar_ga.theta();
  if (zi >= zcut)  kpi.boost(0., 0., beta);

  ki = kpi.rotateY(thega);
  ki = ki.rotateZ(phiga);
  //  ki.setMass(ZERO);
  //  ki.rescaleEnergy();
  ki.rescaleMass();
  return ki;
} 

// Born matrix element
double MEPP2VGammaPowheg::MatrBorn(Lorentz5Momentum k_a, Lorentz5Momentum k_b, Lorentz5Momentum k_ga,
				   Lorentz5Momentum k_V, double char0, double char1, Energy2 MassV2) const{
  Energy2 bs, bt, bu;
  double coeff, scfact, matr, sinw, cosw;

  bs = (k_a + k_b).m2();
  bt = (k_a - k_V).m2();
  bu = (k_a - k_ga).m2();

  if(_boson->id()!=ParticleID::Z0) {
    coeff = sqr(4.0*Pi*alphae)* ckm/sin2w;}
  else {
    sinw = sqrt(sin2w);
    cosw = sqrt(1.0-sin2w);
    coeff = sqr(4.0*Pi*alphae)* 2.0*(sqr(I3quark/(sinw*cosw) -char0*sinw/cosw) +sqr(char0*sinw/cosw));}
  scfact = 1.0/12.0;
  matr = 4.0*sqr(char0*bt + char1*bu)/(bt*bu*sqr(bt+bu))*(bs*MassV2 -bt*bu +sqr(bt+bu)/2.0);
  return matr*coeff*scfact;
}

 
// Tree level matrix element for gq->qV
  // merge from VGammaHardGenerator
double MEPP2VGammaPowheg::MatrQV(Lorentz5Momentum p_a, Lorentz5Momentum p_b, Lorentz5Momentum k_q,
				   Lorentz5Momentum k_V, double charqr, Energy2 MassV2, double alphas, int fi) const{
  Energy2 bs, bt, bu;
  double coeff, scfact, matr, sinw, cosw;
  Lorentz5Momentum k_a, k_b;
  if (fi==0) {
    k_a = p_a;
    k_b = p_b;}
  else {
    k_a = p_b;
    k_b = p_a;}

  bs = (k_a + k_b).m2();
  bt = (k_a - k_q).m2();
  bu = (k_a - k_V).m2();

  if(_boson->id()!=ParticleID::Z0) {
    coeff = sqr(4.0*Pi)*alphae*alphas* ckm/(2.0*sin2w);}
  else {
    sinw = sqrt(sin2w);
    cosw = sqrt(1.0-sin2w);
    coeff = sqr(4.0*Pi)*alphae*alphas* (sqr(I3quark/(sinw*cosw) -charqr*sinw/cosw) +sqr(charqr*sinw/cosw));}
  // spin = 1/4 colour = 4/(3*8)
  scfact = 1.0/24.0;
  matr = -4.0*(2.0*bu*MassV2 +sqr(bt) +sqr(bs))/(bt*bs);
  return matr*coeff*scfact;
}


double MEPP2VGammaPowheg::qgME(Lorentz5Momentum p_a, Lorentz5Momentum p_b, Lorentz5Momentum k_q,
			       Lorentz5Momentum k_V, Energy2 MassV2, bool calc)  const {
  using namespace ThePEG::Helicity;
  using namespace ThePEG;
  vector<SpinorWaveFunction> qin;
  vector<SpinorBarWaveFunction> qbarin;
  vector<SpinorBarWaveFunction> qout;
  vector<SpinorWaveFunction> qbarout;
  vector<VectorWaveFunction> vout,gin;
  Lorentz5Momentum k_a, k_b;
  if (_iflagcq==0) {
    k_a = p_a;
    k_b = p_b;}
  else {
    k_a = p_b;
    k_b = p_a;}                                                                                                    
  SpinorWaveFunction q_in(k_b, _quark, incoming);
  SpinorBarWaveFunction qbar_in(k_b, _quark, incoming);
  SpinorBarWaveFunction q_out(k_q, _quarkpi, outgoing);
  SpinorWaveFunction qbar_out(k_q, _quarkpi, outgoing);                                                       
  VectorWaveFunction v_out(k_V, _boson, outgoing);                                                
  VectorWaveFunction g_in(k_a, _gluon, incoming);
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
      g_in.reset(2*ix);
      gin.push_back(g_in);                                                      
    }

    v_out.reset(ix);                                                     
    vout.push_back(v_out);                                                     
  }
  Energy2 scaleS, scaleEM;
  scaleS = MassV2;
  scaleEM = scale();

  vector<Complex> diag(2,0.0);
  double output(0.), colspin;
  AbstractFFVVertexPtr vertex =
    _boson->id()==ParticleID::Z0 ? FFZvertex_ : FFWvertex_;

  for(unsigned int ihel1=0;ihel1<2;++ihel1) {
    for(unsigned int ihel2=0;ihel2<2;++ihel2) {
      for(unsigned int ohel1=0;ohel1<2;++ohel1) {
        for(unsigned int ohel2=0;ohel2<3;++ohel2) {
	  // intermediates for the diagrams                                                                                                      
	  // suppose parton_b is quark:                                                                                                        
	  if (_quark->id()>0){
                                                                                                  
	  SpinorBarWaveFunction interb=FFGvertex_->evaluate(scaleS,5,_quarkpi->CC(),
				      qout[ohel1],gin[ihel2]);
	  SpinorWaveFunction inters=FFGvertex_->evaluate(scaleS,5,_quark,
				      qin[ihel1],gin[ihel2]);
	  diag[0]=vertex->evaluate(scaleEM,qin[ihel1],interb,
				       vout[ohel2]  );
	  diag[1]=vertex->evaluate(scaleEM,inters,qout[ohel1],
				       vout[ohel2]  );
	  }
	  else{
	  SpinorWaveFunction interb=FFGvertex_->evaluate(scaleS,5,_quarkpi->CC(),
					qbarout[ohel1],gin[ihel2]);
	  SpinorBarWaveFunction inters=FFGvertex_->evaluate(scaleS,5,_quark,
					qbarin[ihel1],gin[ihel2]);
	  diag[0]=vertex->evaluate(scaleEM,interb,qbarin[ihel1],
					 vout[ohel2]  );
	  diag[1]=vertex->evaluate(scaleEM,qbarout[ohel1],inters,
					 vout[ohel2]  );
          }
	  // individual diagrams                                                                                                             
	  //for (size_t ii=0; ii<2; ++ii) me[ii] += std::norm(diag[ii]);
	  // full matrix element                                                                         
	  diag[0] += diag[1];
	  output += std::norm(diag[0]);
      // storage of the matrix element for spin correlations                                                     
	  //  if(calc) me_(ihel1,ihel2,ohel1,ohel2) = diag[0];
	}
      }
    }
  }
  //DVector save(3);                                        
  colspin=1./24./4.;                                                                                                                        
  colspin *=4.;
  //for(size_t ix=0;ix<3;++ix) {
  //  save[ix] = colspin * me[ix];
  //}
  //meInfo(save);
  output *= colspin;
  return output;
}


// Gluon Real radiation matrix element
double MEPP2VGammaPowheg::MatrRealGluon(Lorentz5Momentum k_a, Lorentz5Momentum k_b, Lorentz5Momentum k_ga,
					Lorentz5Momentum k_V, Lorentz5Momentum k_glu, Energy2 scaleS) const{
  using namespace ThePEG::Helicity;
  using namespace ThePEG;
  double sum(0.), scfact;
  //Energy2 bosonMass2_;
  vector<SpinorWaveFunction> qin;
  vector<SpinorBarWaveFunction> qbarin;
  vector<VectorWaveFunction> wout,pout,gout;
  SpinorWaveFunction q_in(k_a, _partona, incoming);
  SpinorBarWaveFunction qbar_in(k_b, _partonb, incoming);
  VectorWaveFunction v_out(k_V  ,_boson, outgoing);
  VectorWaveFunction p_out(k_ga, _photon, outgoing);
  VectorWaveFunction g_out(k_glu, _gluon, outgoing);

  //bosonMass2_ = k_V.m2();
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix<2) {
      q_in.reset(ix);
      qin.push_back(q_in);
      qbar_in.reset(ix);
      qbarin.push_back(qbar_in);
      g_out.reset(2*ix);
      gout.push_back(g_out);
      p_out.reset(2*ix);
      pout.push_back(p_out);
    }
    v_out.reset(ix);
    wout.push_back(v_out);
  }
  Energy2 scaleEM;
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

	    sum += norm(dsum);
	  }
	}
      }
    }
  }
  // final spin and colour factors
  // spin = 1/4 colour = 4/9
  scfact = 1.0/9.0;  
  return sum*scfact*1000000.0;
}


// Quark Real radiation matrix element
double MEPP2VGammaPowheg::MatrRealQuark(Lorentz5Momentum p_a, Lorentz5Momentum p_b, Lorentz5Momentum k_ga,
					Lorentz5Momentum k_V, Lorentz5Momentum k_q, Energy2 scaleS, int flavorflag) const{
  using namespace ThePEG::Helicity;
  using namespace ThePEG;
  double sum(0.), scfact;
  //Energy bosonMass_;
  vector<SpinorWaveFunction> qin;
  vector<SpinorBarWaveFunction> qbarin;
  vector<SpinorBarWaveFunction> qout;
  vector<SpinorWaveFunction> qbarout;
  vector<VectorWaveFunction> vout,pout,gin;
  Lorentz5Momentum k_a, k_b;
  if (flavorflag%2==0) {
    k_a = p_a;
    k_b = p_b;}
  else {
    k_a = p_b;
    k_b = p_a;}
  SpinorWaveFunction q_in(k_b, _quark, incoming);
  SpinorBarWaveFunction qbar_in(k_b, _quark, incoming);
  SpinorBarWaveFunction q_out(k_q, _quarkpi, outgoing); 
  SpinorWaveFunction qbar_out(k_q, _quarkpi, outgoing); 
  VectorWaveFunction v_out(k_V, _boson, outgoing);
  VectorWaveFunction p_out(k_ga, _photon, outgoing);
  VectorWaveFunction g_in(k_a, _gluon, incoming);

  //bosonMass_ = k_V.mass();
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
      g_in.reset(2*ix);
      gin.push_back(g_in);
      p_out.reset(2*ix);
      pout.push_back(p_out);
    }

    v_out.reset(ix);
    vout.push_back(v_out);
  }
 
  Energy2 scaleEM;
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
         //if (_quark->id()>0){
         if (flavorflag%2==1) {
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
	    sum += norm(dsum);
	  }
	}
      }
    }
  }
  // final spin and colour factors
  // spin = 1/4 colour = 4/(3*8)
  scfact = 1.0/24.0;  
  return sum*scfact*1000000.;
}

// dipole of gluon radiation
  // merge from VGammaHardGenerator
double MEPP2VGammaPowheg::Dipole(Lorentz5Momentum k_a, Lorentz5Momentum k_b, Lorentz5Momentum k_ga,
				 Lorentz5Momentum k_V, double char0, double char1, Energy2 MassV2, Energy2 kig, double alphas, double xi) const{
  double coeff, fact;
  coeff = 8.0*Pi*alphas*_CF;
  fact = (1.0 + sqr(xi))/(1.0 - xi)/(2.0*xi *kig)*sqr(GeV);
  return MatrBorn(k_a,k_b,k_ga,k_V,char0,char1,MassV2) *coeff*fact;
}


// dipole of quark radiation in the singular region collinear with final state photon
  // merge from VGammaHardGenerator
double MEPP2VGammaPowheg::Dipolepqr(Lorentz5Momentum k_a, Lorentz5Momentum k_b, Lorentz5Momentum k_ga,
				 Lorentz5Momentum k_V, double chargr, Energy2 MassV2, Energy2 kig, double zi, double alphas, int fi) const{
  double coeff, fact;
  coeff = 8.0*Pi*alphae*sqr(chargr);
  fact = (1.0 + sqr(1.0-zi))/zi/(2.0 *kig)*sqr(GeV);
  //cerr<<"kig="<<kig/sqr(GeV)<<" coeff="<<coeff<<" fact="<<fact<<"\n";
  return MatrQV(k_a,k_b,k_ga,k_V,chargr,MassV2,alphas,fi) *coeff*fact;
}


// dipole of quark radiation in the singular region collinear with initial state gluon
  // merge from VGammaHardGenerator
double MEPP2VGammaPowheg::Dipolegluqr(Lorentz5Momentum k_a, Lorentz5Momentum k_b, Lorentz5Momentum k_ga,
				 Lorentz5Momentum k_V, double char0, double char1, Energy2 MassV2, Energy2 kig, double alphas, double xi) const{
  double coeff, fact;
  coeff = 8.0*Pi*alphas*_TR;
  fact = (1.0 - 2.0*xi*(1.0 - xi))/(2.0*xi *kig)*sqr(GeV);
  return  MatrBorn(k_a,k_b,k_ga,k_V,char0,char1,MassV2) *coeff*fact;
}


// parton distribution function of beam partID
double MEPP2VGammaPowheg::PDFab(double z, int partID, Energy2 scale, tcBeamPtr hadron_A, tcBeamPtr hadron_B)  const {
  if(partID==0){
    //assert(_hadron_A);
    return hadron_A->pdf()->xfx(_hadron_A,_partona,scale,z);}  
  else{
    //assert(_hadron_B);
    return hadron_B->pdf()->xfx(_hadron_B,_partonb,scale,z);}
}


// parton distribution function ratio
double MEPP2VGammaPowheg::PDFratio(double z, double x, Energy2 scale, int partID, int fi, tcBeamPtr hadron_A, tcBeamPtr hadron_B)  const {
  double pdfz, pdfx;
  //quark PDF ratio: z is the new fraction, and x the old fraction
  if (partID <3){
  pdfz = PDFab(z, partID, scale, hadron_A, hadron_B)/z;
  pdfx = PDFab(x, partID, scale, hadron_A, hadron_B)/x;}
  //gluon vs. quark PDF ratio with the same fraction: z and x is the fractions of the two initial states
  if (partID ==3){
    if (fi ==0){
      pdfx = PDFab(z, 0, scale, hadron_A, hadron_B);
      pdfz = hadron_A->pdf()->xfx(_hadron_A,_gluon,scale,z);}
    else {
      pdfx = PDFab(x, 1, scale, hadron_A, hadron_B);
      pdfz = hadron_B->pdf()->xfx(_hadron_B,_gluon,scale,x);}}
  //gluon vs. quark PDF ratio with the different fraction: z is the new fraction of gluon and x is the old one of quark
  if (partID ==4){
    if (fi ==0){
      pdfx = PDFab(x, 0, scale, hadron_A, hadron_B)/x;
      pdfz = hadron_A->pdf()->xfx(_hadron_A,_gluon,scale,z)/z;}
    else {
      pdfx = PDFab(x, 1, scale, hadron_A, hadron_B)/x;
      pdfz = hadron_B->pdf()->xfx(_hadron_B,_gluon,scale,z)/z;}}

  if (pdfx == 0.0) {
    if (pdfz == 0.0){
      return 1.0;}
    else{
      return 0.0;}}
  else{
    return pdfz/pdfx;}

}



// functionals of PDFs and energy fractions in the collinear remnants for gluon radiation:

double MEPP2VGammaPowheg::KP1(double zz, Energy2 shat, Energy2 muF2, int partID) const{
  double x, hpart;
     x = 1.0 -_xtil +zz*_xtil;
     hpart= log(sqr(1.0-x)*shat/(x*muF2))*
       (2.0/(1.0-x) -(1.0+sqr(x))/(1.0-x)*PDFratio(zz/x,zz,muF2,partID,_iflagcq,_hadron_A,_hadron_B)/x) +2.0/(1.0-x)*log(x);
     hpart = log(x*muF2/(sqr(1.0-x)*shat))*(PDFratio(zz/x,zz,muF2,partID,_iflagcq,_hadron_A,_hadron_B)/x*(1.0+x*x)/(1.0-x) - 2.0/(1.0-x))
       -2.0*log(sqr(1.0-x));
     return (1.0-zz)*hpart;
}

double MEPP2VGammaPowheg::KP2(double zz, Energy2 shat, Energy2 muF2) const{
  double x, hpart;
  /*
  x = zz*_xtil;
     hpart= log(sqr(1.0-x)*shat/muF2)*(2.0/(1.0-x));
     return zz*hpart;
  */
  hpart = (3.0/2.0 - (zz + zz*zz/2.0))*log(muF2/shat) +4.0*(1.0-zz)*(log(1.0-zz) - 1.0);
  return hpart;
}

double MEPP2VGammaPowheg::KP3(double zz,  int partID) const{
  double x, hpart;
     x = 1.0 -_xtil +zz*_xtil;
     hpart= (1.0-x)*PDFratio(zz/x,zz,_muF2,partID,_iflagcq,_hadron_A,_hadron_B)/x;
     return (1.0-zz)*hpart;
}

// the collinear remnants for quark radiation:
double MEPP2VGammaPowheg::KPpr(double xx, Energy2 shat, Energy2 muF2) const{
  double x1, x2, hpart, px1, px2, lnx1, lnx2, part1, part2, part3, pdffact;
     x1 = 1.0 -_xtil +xx*_xtil;
     //x2 = xx*_xtil;
     pdffact= PDFratio(xx/x1,xx,muF2,4,_iflagcq,_hadron_A,_hadron_B)/x1;
     px1= 1.0 - 2.0*x1*(1.0-x1);
     /*
     px2= 1.0 - 2.0*x2*(1.0-x2);
     lnx1= log(shat*sqr(1.0-x1)/muF2);
     lnx2= log(shat*sqr(1.0-x2)/muF2);
     part1= (_TR*(1.0-px1)- gamnum*px1)*pdffact;
     part2= px1*lnx1*(pdffact -1.0);
     part3= px2*lnx2;
     hpart= (1.0-xx)*(part1+part2) - xx*part3;
     */
     lnx1 = log(shat*sqr(1.0-x1)/x1/muF2);
     hpart= (1.0-xx)*_TR*((px1*lnx1 + (1.0-px1))*pdffact - 2.0*log(1.0-x1)*PDFratio(xx,xx,muF2,4,_iflagcq,_hadron_A,_hadron_B)) + _TR*2.0*(1.0-xx)*(log(1.0-xx)-1.0)*PDFratio(xx,xx,muF2,4,_iflagcq,_hadron_A,_hadron_B);
     return hpart;
}



// definitions of H, F^V:
double MEPP2VGammaPowheg::Hfunc(Energy2 t, Energy2 s, Energy2 m2) const{  
  return sqr(Pi)-sqr(log(s/m2))+sqr(log(-t/s))-sqr(log(-t/m2))
    -2.0*ReLi2(1.0-s/m2)-2.0*ReLi2(1.0-t/m2);
} 


double MEPP2VGammaPowheg::FWfunc(Energy2 t, Energy2 u, Energy2 s, Energy2 m2) const{
    double y1,y2,y3,y4,y5;
    y1= 4.0*(2.0*s*s/(t*u) +2.0*s/u +t/u)* Hfunc(u,s,m2);
    y2= -8.0/3.0*sqr(Pi)*s*(2.0*s/t +t/s -u/s)/(t+u);
    y3= 4.0*(6.0-10.0*u/(t+u)-10.0*s*s/(u*(t+u))-11.0*s/u-5.0*t/u+2.0*s/(t+u)+s/(s+t));
    y4= -4.0*log(s/m2)*(3.0*t/u+2.0*s/u+4.0*s*(t+s)/(u*(t+u))+2.0*t/u*sqr(s/(t+u)));
    y5= 4.0*log(-u/m2)*((4.0*s+u)/(s+t)+s*u/sqr(s+t));
    return y1+y2+y3+y4+y5;
}

double MEPP2VGammaPowheg::FZfunc(Energy2 t, Energy2 u, Energy2 s, Energy2 m2) const{
    double y1,y2,y3,y4,y5;
    y1= 4.0*(2.0*s*s/(t*u) +2.0*s/u +t/u)* Hfunc(u,s,m2);
    y2= -8.0/3.0*sqr(Pi)*s*s/(t*u);
    y3= 4.0*(1.0-5.0*s*s/(t*u)-11.0*s/u-5.0*t/u+2.0*s/(t+u)+s/(s+t));
    y4= 4.0*log(s/m2)*(4.0*s/(t+u)+2.0*sqr(s/(t+u))-3.0*sqr(s+t)/(t*u));
    y5= 4.0*log(-u/m2)*((4.0*s+u)/(s+t)+s*u/sqr(s+t));
    return y1+y2+y3+y4+y5;
}

double MEPP2VGammaPowheg::FVfunc(Energy2 t, Energy2 u, Energy2 s, Energy2 m2) const {
  if(mePartonData()[_idboson]->id()==ParticleID::Z0) {
    return FZfunc(t,u,s,m2);
  }
  else {
    return FWfunc(t,u,s,m2);
  }
}

// ratio of NLO/LO 
double MEPP2VGammaPowheg::NLOweight() const {
  //  double Pi(3.1415926);
  double partep, partloop, partc, alfsfact;
  double smborn0, smborn1, smborn2, smbfact, smloop;
  Energy2 kadotg, kbdotg, scale_S; //MV2,
  double sh,shr,ratea,rateb,va,vb,txa,zx,tva,tvb; //charge0,charge1,
  Lorentz5Momentum k_a, k_b, k_ga, k_V, k_glu, kbarp_ga, kbarp_V;
  double xiab, vi, phi;
  double sumdipole, partreal, partdipole, smreala, smrealb, smreal, jacobfact, pdffact, fluxfact;
  int fi,ida, idb;
  //  double FV, FW, FZ;
  // If only leading order is required return 1:
   if(_contrib==0) return 1.;
   
  _ss= sHat();
  _tt= ( _p_partona - _p_boson).m2();
  _uu= ( _p_partona - _p_photon).m2();
  MV2= _p_boson.m2();
 
  ida = _partona->id();
  idb = _partonb->id();
  charge0= _partona->iCharge()/3.*ida/abs(ida);
  charge1= _partonb->iCharge()/3.*idb/abs(idb);

  // finite contribution from gluon virtual loop:
  smbfact= 4.0*sqr(charge0*_tt + charge1*_uu)/(_tt*_uu*sqr(_tt+_uu))*GeV*GeV*GeV*GeV;
  smborn0= (_ss*MV2 -_tt*_uu +sqr(_tt+_uu)/2.0)/(GeV*GeV*GeV*GeV);
  smborn1= -(_ss*MV2 -_tt*_uu +sqr(_tt+_uu))/(GeV*GeV*GeV*GeV);
  smborn2= (sqr(_tt+_uu)/2.0)/(GeV*GeV*GeV*GeV);
  partep= 2.0*(5.0-sqr(Pi)/3.0)*1.0 +3.0*smborn1/smborn0 +2.0*smborn2/smborn0;

  smloop= (charge0*_tt+ charge1*_uu)*
       (charge0*FVfunc(_tt,_uu,_ss,MV2) +charge1*FVfunc(_uu,_tt,_ss,MV2))/(_tt+_uu)/2.0;
  partloop= smloop/(smborn0*smbfact);  

  // collinear remnants for gluon radiation:
  partc= -(KP1(_xa,_ss,_muF2,0)+KP2(_xa,_ss,_muF2)-KP3(_xa,0))
    -(KP1(_xb,_ss,_muF2,1)+KP2(_xb,_ss,_muF2)-KP3(_xb,1))-2.0*(5.0-2.0*sqr(Pi)/3.0);
  // coupling factor
  alfsfact= _alphas*_CF/(2.0*Pi);
  scale_S = _p_boson.m2();

  //contributions fron real gluon radiation:
  phi = _y3*(2.0*Pi);
  //  in the singular for gluon soft or collinear with parton_a:
  fi = 0;
  if (_xa<=(1.0 - 0.0000001)) {
  xiab = (1.0 - 0.0000001)*(1.0-_y1) + _xa*_y1;
  vi = (1.0 -xiab)*_y2 + 0.0000001*(1.0 - _y2);
  // map to real phase space:
  k_a = _p_partona/xiab;
  k_b = _p_partonb;
  k_glu =  radk(_p_partona, _p_partonb, xiab, vi, phi, fi); 
  k_ga = InvLortr(_p_partona, _p_partonb, xiab, k_glu, _p_photon, fi);
  k_V = InvLortr(_p_partona, _p_partonb, xiab, k_glu, _p_boson, fi);

  kadotg = k_a.dot(k_glu);
  kbdotg = k_b.dot(k_glu);
 
  // map to Born phase space for the other (not singular region) dipole in the denominator dipole sum
  kbarp_ga = Lortr(k_a, k_b, xiab, k_glu, k_ga, 1);
  kbarp_V  = Lortr(k_a, k_b, xiab, k_glu, k_V, 1);

  partdipole = Dipole(_p_partona, _p_partonb, _p_photon, _p_boson, charge0, charge1, MV2, kadotg, _alphas, xiab);
  sumdipole = partdipole +Dipole(k_a, xiab*k_b, kbarp_ga, kbarp_V, charge0, charge1, MV2, kbdotg, _alphas, xiab);

  partreal = MatrRealGluon(k_a, k_b, k_ga, k_V, k_glu, scale_S)/sumdipole - 1.0;
  jacobfact = 2.0*k_a.dot(k_b)*(1.0 -xiab)*(1.0 -_xa)/(16.0*sqr(Pi))/sqr(GeV);
  pdffact = PDFratio(_xa/xiab, _xa, _muF2,0,_iflagcq,_hadron_A,_hadron_B);
  fluxfact = _ss/(k_a+k_b).m2();
  smreala = partreal *partdipole*jacobfact*pdffact*fluxfact;
    }
   else smreala = 0.0; 

 //  in the singular for gluon soft or collinear with parton_b:
  fi = 1;
    if (_xb<=(1.0 - 0.0000001)) {
  xiab = (1.0 - 0.0000001)*(1.0-_y1) + _xb*_y1;
  vi = (1.0 -xiab)*_y2 + 0.0000001*(1.0 - _y2);
  // map to real phase space:  
  k_a = _p_partona;
  k_b = _p_partonb/xiab;
  k_glu =  radk(_p_partona, _p_partonb, xiab, vi, phi, fi); 
  k_ga = InvLortr(_p_partona, _p_partonb, xiab, k_glu, _p_photon, fi);
  k_V = InvLortr(_p_partona, _p_partonb, xiab, k_glu, _p_boson, fi);

  kadotg = k_a.dot(k_glu);
  kbdotg = k_b.dot(k_glu);

  // map to Born phase space for the other (not singular region) dipole in the denominator dipole sum
  kbarp_ga = Lortr(k_a, k_b, xiab, k_glu, k_ga, 0);
  kbarp_V  = Lortr(k_a, k_b, xiab, k_glu, k_V, 0);

  partdipole = Dipole(_p_partona, _p_partonb, _p_photon, _p_boson, charge0, charge1, MV2, kbdotg, _alphas, xiab);
  sumdipole = partdipole +Dipole(xiab*k_a, k_b, kbarp_ga, kbarp_V, charge0, charge1, MV2, kadotg, _alphas, xiab);

  partreal = MatrRealGluon(k_a, k_b, k_ga, k_V, k_glu, scale_S)/sumdipole -1.0;
  jacobfact = 2.0*k_a.dot(k_b)*(1.0 -xiab)*(1.0 -_xb)/(16.0*sqr(Pi))/sqr(GeV);
  pdffact = PDFratio(_xb/xiab, _xb, _muF2, 1, _iflagcq, _hadron_A, _hadron_B);
  fluxfact = _ss/(k_a+k_b).m2();
  smrealb = partreal *partdipole*jacobfact*pdffact*fluxfact;
    }
   else smrealb = 0.0;
  // the total subtracted gluon radiation piece
  smreal = (smreala + smrealb)/MatrBorn(_p_partona, _p_partonb, _p_photon, _p_boson, charge0, charge1, MV2);

  
  // collinear remnants for quark radiation: d\sigma^A+d\sigma^C 
  double zi, zq, uq, ulim, phiq, zlim, shqr, partinitqr, xip, xiabqr, viqr, sumdipoleqr, smrealqrp, smrealqrg;  //chargeqr, 
  double smbornqr, hfs, lnz, pz, lnmf, partqr, smrealqrz_0, pdffactqr, factqr, smcrqr, partrealqr, partdipoleqr, jacobfactqr, smrealqr, conu;
  double zcut, RDz_0, Rz_0, Dz_0, CAz_0, charge_i, charge_V, charge_q;
  Lorentz5Momentum plp, plv, plq, kquark, kphoton, kboson, kgluon, kb_photon, kb_boson, kq_z0, kV_z0;
  Energy QCDcut;
  Energy2 kpdotq, kgdotq;
  int idq;
  double thetacut, Rcut;
 
  // photon fragmentation function collinear remnant
  //fraccut = 0.4;
  //fraccut = 0.3;
  //thetacut = 1.0 - 0.6;
  thetacut = 1.0 - 0.4;
  Rcut = 0.01;

  QCDcut = 0.14*GeV;  
  //QCDcut = 0.2*GeV;
  //lnmf = log(_muF2/sqr(QCDcut));
  //zlim = abs(25.0*GeV/_p_photon.perp());
  //zlim = (_ss - MV2)/(_ss/(_xa*_xb) - MV2);
  // boost from Born CM frame to lab frame to calculate the photon fraction constraint
  tva = (_xa-_xb)/(_xa+_xb);  
  //plp = _p_photon;
  //plp.boost(0.,0.,tva);
  //plv = _p_boson;
  //plv.boost(0.,0.,tva);

  //*******************************************************************************
  //**********  Here is the photon fraction constraint ****************************
  //**********  I impose in fragmantation and real mapping ************************
  //**********  where quark radiate from photon  **********************************
  //**********  The first line of zlim is from photon energy constraint ***********
  //**********  the second line is from experimental Kt cut ***********************
  //*******************************************************************************  
  //  zlim = min(max(max(plp.perp()*max(exp(plp.rapidity())+exp(plv.rapidity()),exp(-plp.rapidity())+exp(-plv.rapidity()))/(_ss/(_xa *_xb) - MV2)*
  //	     sqrt(_ss/(_xa *_xb)), lastCuts().minKT(_photon)*cosh(plp.rapidity())/plp.e()), 0.5), 1.0-0.0000001);
  //zlim = 0.1*fraccut;
  zlim = 0.1;
  //zlim = min(QCDcut/_p_photon.e(),1.0);
  //zlim = max(max(plp.perp()*max(exp(plp.rapidity())+exp(plv.rapidity()),exp(-plp.rapidity())+exp(-plv.rapidity()))/(_ss/(_xa *_xb) - MV2)*
  //		 sqrt(_ss/(_xa *_xb)), lastCuts().minKT(_photon)*cosh(plp.rapidity())/plp.e()), zcut);
  //  zlim = abs(zlim);
  //zlim = max(plp.perp()*max(exp(plp.rapidity())+exp(plv.rapidity()),exp(-plp.rapidity())+exp(-plv.rapidity()))/(_ss/(_xa *_xb) - MV2)*sqrt(_ss/(_xa *_xb)), lastCuts().minKT(_photon)/plp.e() );
  //zlim = min(max(0.3, lastCuts().minKT(_photon)*cosh(plp.rapidity())/plp.e()), 1.0-0.0000001);
  //zlim = 0.1;
  //zlim =  (_ss - MV2)/_ss;

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //+++++++++++++  Here we use FKS subtraction to single out the soft photon ++++++++++++
  //+++++++++++++  divergence to make calculaton safe from singularities     ++++++++++++
  //+++++++++++++  We have aux photon fracton cut z_cut, but the total cross ++++++++++++
  //+++++++++++++  section is independent of z_cut. RDz_0 and CAz_0 is defined+++++++++++
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //  zlim = 0.0;
  zcut = 1.0;
  //zlim = 0.5;
  //************ photon fraction in fragmantation function where zlim is applied ***********************
  zi = 1.0 - _y1 + zlim*_y1;
  //zi = 1.0 - _y1; 
  idq = _quarkpi->id();
  chargeqr = _quarkpi->iCharge()/3.*idq/abs(idq);
  charge_q = _quarkpi->iCharge()/3.;
  charge_i = _quark->iCharge()/3.;
  charge_V = charge_i - charge_q;
  smbornqr = MatrQV(_p_partona, _p_partonb, _p_photon, _p_boson, chargeqr, MV2, _alphas, _iflagcq);
  // the fraction after exchanging the integration
  //zx = zlim/zi;
  zx = zi;
  // another fragmantation function in Baur's paper
  //hfs = (sqr(chargeqr)*(2.21-1.28*zx+1.29*sqr(zx))/(1.0-1.63*log(1.0-zx))*pow(zx,0.049)+0.002*sqr(1.0-zx)*pow(zx,-1.54))*log(_muF2/sqr(QCDcut));
  // A-P splitting function
  pz = (1.0+sqr(1.0-zx))/zx;
  // the aux function I impose to have right mapping in soft photon region
  //lnmf= log(1.0+zx)/log(2.0);
  lnmf = 1.0;
  // u variable upper limit
  ulim = 2.0*(1.0-zx)*_p_photon.dot(_p_boson)*lnmf/(MV2+2.0*_p_photon.dot(_p_boson) -2.0*zx*lnmf*_p_photon.dot(_p_boson));
  // our fragmantation function
  //lnz = log(2.0*_p_photon.dot(_p_boson)*zx*zx*(1.0-zx)*ulim/_muF2);
  lnz = log(2.0*_p_photon.dot(_p_boson)*(1.0-zx)*ulim/_muF2);
  //hfs = pz*log(_muF2/sqr(QCDcut*(1.0-zx))) -13.26;
  if (zi>0.7) hfs = pz*log(_muF2/sqr(QCDcut*(1.0-zx))) -13.26;
  //else hfs = pz*log(1.0/sqr(1.0-zx)) -13.26 + pz*log(_muF2/sqr(QCDcut))*(3.0*sqr(zi/0.7) - 2.0*pow(zi/0.7,3)); 
  else hfs = pz*log(_muF2/sqr(1.0-zx)/(2.0*_p_photon.dot(_p_boson))) -13.26 + pz*log((2.0*_p_photon.dot(_p_boson))/sqr(QCDcut))*(3.0*sqr(zi/0.7) - 2.0*pow(zi/0.7,3));
  //hfs = pz*log(_muF2/sqr(QCDcut*(1.0-zx)));

  /*
  // CA(z=0) for cancelling the soft photon singularity
  CAz_0 = 2.0* log(2.0*_p_photon.dot(_p_boson)*lnmf*2.0*_p_photon.dot(_p_boson)/(sqr(QCDcut)*(MV2+2.0*_p_photon.dot(_p_boson))));
  // photon fragmentation function collinear remnant piece
  //partqr =sqr(chargeqr)*(pz*lnz +zx +hfs);
  if (zi<zcut) {
    partqr =sqr(chargeqr)*((pz*lnz +zx +hfs)*zi -CAz_0);}
  else {
    partqr =sqr(chargeqr)*((pz*lnz +zx +hfs)*zi);}
  */

  // CA(z=0) for cancelling the soft photon singularity
  CAz_0 = 2.0*log(2.0*_p_photon.dot(_p_boson)*lnmf/(MV2+2.0*_p_photon.dot(_p_boson)));
   if (zi >= zlim) 
     partqr = sqr(chargeqr)*((pz*lnz +zx +hfs)*zi - CAz_0);
     //partqr = sqr(chargeqr)*((pz*log((1.0-zx)*ulim) +zx -13.26)*zi - CAz_0);
   //partqr = sqr(chargeqr)*(-13.26)*zi;
    //    cout <<" zi > zlim: "<<zi<<" partqr= "<<partqr << "\n";}
  // else if (zi < zlim && zi >=zcut) 
  //  partqr = sqr(chargeqr)*((pz*log((1.0-zx)*ulim) +zx)*zi);
    //    cout <<" zi <= zlim && zi >zcut: "<<zi<<" partqr= "<<partqr << "\n";}
   else   
    partqr = sqr(chargeqr)*((pz*log((1.0-zx)*ulim) +zx)*zi - CAz_0);
    //    cout <<" zi <zcut: "<<zi<<" partqr= "<<partqr << "\n";}

  //separatation cut:
    if (zi < fraccut) {
      //cout<<" separatation cut the photon FF remnant: z_frac="<<zi<<" partqr="<<partqr<<"\n";
      partqr = 0.0;
    }

  //partqr =sqr(chargeqr)*((pz*lnz +zx +hfs)*zi);
  pdffactqr = PDFratio(_xa, _xb, _muF2, 3, _iflagcq, _hadron_A, _hadron_B);
  factqr = alphae/(2.0*Pi)*(1.0 - zlim)/zi;
  //factqr = alphae/(2.0*Pi)/zi;

  //  PDF collinear remnants for quark radiation from initial quark:
  if (_iflagcq==0){
    partinitqr = KPpr(_xa,_ss,_muF2)*_alphas/(2.0*Pi);}
  else {
    partinitqr = KPpr(_xb,_ss,_muF2)*_alphas/(2.0*Pi);}
  // sum of collinear remnants for quark radiation
  //smcrqr = smbornqr*(partqr + sqr(chargeqr)*log(zcut)*CAz_0*zi)*factqr*
  //          pdffactqr/MatrBorn(_p_partona, _p_partonb, _p_photon, _p_boson, charge0, charge1, MV2)+ partinitqr;
  smcrqr = smbornqr*partqr*factqr*pdffactqr/MatrBorn(_p_partona, _p_partonb, _p_photon, _p_boson, charge0, charge1, MV2)+ partinitqr;
 
  //smcrqr = partinitqr;  
  
  // contribution from real quark radiation:
  phiq = 2.0*Pi*_y3;


  // in the singular region that quark radiation collinear with final state photon (Dpq):

  
  //************ photon fraction in quark-photon collinear real piece where zlim is applied ***********************
  //zq = 1.0 - _y1 + zlim*_y1;
  //zlim = 0.1;
  zq = (1.0 - MV2/2.0/_p_photon.dot(_p_boson)*1.0e-7/(1.0-1.0e-7))*(1.0 - _y1) + zlim*_y1;
  //zq = (1.0 - MV2/2.0/_p_photon.dot(_p_boson)*1.0e-7/(1.0-1.0e-7))*(1.0 - _y1);
  // the aux function I impose to have right mapping in soft photon region
  //  conu = log(1.0+zq)/log(2.0);
  conu = 1.0;
  // u variable upper limit
  ulim = 2.0*(1.0-zq)*_p_photon.dot(_p_boson)*conu/(MV2+2.0*_p_photon.dot(_p_boson) -2.0*zq*conu*_p_photon.dot(_p_boson));
  //  uq = ulim*_y2;
  //  if (ulim > 0.0000001)
    uq = ulim - (ulim-1.0e-7)*_y2;
  //  else uq = 0.0000001;
  //  cout << " zlim="<<zlim<< " zq="<<zq<<" uq="<<uq<<"\n";
  // map to real phase space: 
  kquark = radkqr(_p_photon, _p_boson, MV2, zq, uq, phiq, zlim);
  kphoton = _p_photon*zq;
  //  kphoton.setMass(ZERO);
  //  kphoton.rescaleEnergy();
  kboson = _p_boson -kquark +(1.0-zq)/zq*kphoton;
  kboson.rescaleMass();
  //  kboson.setMass(sqrt(MV2));
  //  kboson.rescaleEnergy();
  kpdotq = kphoton.dot(kquark);
  //  plp = kphoton + kquark;
  //  kpdotq = plp.m2()/2.0;

  //re-map the Born phase space so that no numerical errers sensitivity:
  zq = kphoton.dot(kquark+kboson)/(kphoton.dot(kquark+kboson)+kquark.dot(kboson));
  uq = kpdotq/kphoton.dot(kquark+kboson);
  kb_photon = kphoton/zq;
  kb_boson = kquark + kboson - (1.0 - zq)/zq*kphoton;
  kb_boson.rescaleMass();
  //  kb_boson.setMass(sqrt(MV2));
  //  kb_boson.rescaleEnergy();
  kq_z0 = radkqr(_p_photon, _p_boson, MV2, 0.0, uq, phiq, zlim);
  kV_z0 = _p_boson -kq_z0 + _p_photon;
  kV_z0.rescaleMass();
  //if (uq<0.0000001)
  //cout <<" zq="<<zq<<" uq="<<uq<<" kqdot="<<kq_z0.dot(_p_photon)/sqr(GeV)<<" kVdot="<<(kV_z0.dot(_p_boson)-MV2)/sqr(GeV)<<" conserve?="<<
  //  _p_partona.isNear(kV_z0+kq_z0-_p_partonb, 0.00000001)<<"\n";
   //  cout <<" ratio of k_q.k_p="<<kpdotq/(uq*zq*_p_boson.dot(_p_photon)) <<"\n";
   //  if (uq!=kpdotq/kphoton.dot(kquark+kboson)) {cout<<" !!!!!!!!!!!!!!!!!!!!!!!!!uq is wrong!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<"  uq="<<uq<<"  ?="<<kpdotq/kphoton.dot(kquark+kboson)<<"\n";}
   //  if (zq!=kphoton.dot(kquark+kboson)/(kphoton.dot(kquark+kboson)+kquark.dot(kboson))) 
   //       {cout<<" !!!!!!!!!!!!!!!!!!!!!!!!!zq is wrong!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<"  zq="<<zq<<"  ?="
   //     <<kphoton.dot(kquark+kboson)/(kphoton.dot(kquark+kboson)+kquark.dot(kboson))<<"\n";}

   //Born phase space for Dgq in the denominator dipole sum:
   xip = 1.0- kquark.dot(_p_partona+ _p_partonb)/_p_partona.dot(_p_partonb);
   if (_iflagcq==0){
     fi = 0;
     kgdotq = _p_partona.dot(kquark); 
     }
   else {
     fi = 1;
     kgdotq = _p_partonb.dot(kquark);
     }
   kbarp_ga = Lortr(_p_partona, _p_partonb, xip, kquark, kphoton, fi);
   kbarp_V  = Lortr(_p_partona, _p_partonb, xip, kquark, kboson, fi); 

   smbornqr = MatrQV(_p_partona, _p_partonb, kq_z0, kV_z0, chargeqr, MV2, _alphas, _iflagcq);
   partdipoleqr = Dipolepqr(_p_partona, _p_partonb, kb_photon, kb_boson, chargeqr, MV2, kpdotq, zq, _alphas, _iflagcq);
   
   if (_iflagcq==0){
     sumdipoleqr = partdipoleqr + Dipolegluqr(xip*_p_partona,_p_partonb,kbarp_ga,kbarp_V,charge0,charge1,MV2,kgdotq,_alphas,xip);
     RDz_0 = 8.0*Pi*alphae* (((charge_i*charge_q)*_p_partonb.dot(kq_z0)/_p_partonb.dot(kb_photon)/kq_z0.dot(kb_photon) 
			     + (charge_V*charge_i)*_p_partonb.dot(kV_z0)/_p_partonb.dot(kb_photon)/kV_z0.dot(kb_photon)
                             - (charge_V*charge_q)*kq_z0.dot(kV_z0)/kV_z0.dot(kb_photon)/kq_z0.dot(kb_photon))*smbornqr 
                             - sqr(charge_q)/kq_z0.dot(kb_photon)*MatrQV(_p_partona, _p_partonb, kb_photon, kb_boson, chargeqr, MV2, _alphas, _iflagcq) )
                               *sqr(GeV);
     Dz_0 = 8.0*Pi*alphae*sqr(charge_q)/kq_z0.dot(kb_photon)*MatrQV(_p_partona, _p_partonb, kb_photon, kb_boson, chargeqr, MV2, _alphas, _iflagcq)
             *sqr(GeV);
     Rz_0 = RDz_0 + Dz_0;
   }
   else {
     sumdipoleqr = partdipoleqr + Dipolegluqr(_p_partona,xip*_p_partonb,kbarp_ga,kbarp_V,charge0,charge1,MV2,kgdotq,_alphas,xip);
     RDz_0 = 8.0*Pi*alphae* (((charge_i*charge_q)*_p_partona.dot(kq_z0)/_p_partona.dot(kb_photon)/kq_z0.dot(kb_photon)
			      + (charge_V*charge_i)*_p_partona.dot(kV_z0)/_p_partona.dot(kb_photon)/kV_z0.dot(kb_photon)
			      - (charge_V*charge_q)*kq_z0.dot(kV_z0)/kV_z0.dot(kb_photon)/kq_z0.dot(kb_photon))*smbornqr
                             - sqr(charge_q)/kq_z0.dot(kb_photon)*MatrQV(_p_partona, _p_partonb, kb_photon, kb_boson, chargeqr, MV2, _alphas, _iflagcq) )
                                *sqr(GeV);
     Dz_0 = 8.0*Pi*alphae*sqr(charge_q)/kq_z0.dot(kb_photon)*MatrQV(_p_partona, _p_partonb, kb_photon, kb_boson, chargeqr, MV2, _alphas, _iflagcq)
             *sqr(GeV);
     Rz_0 = RDz_0 + Dz_0;
   }
   
   jacobfactqr = 2.0*kphoton.dot(_p_boson)*ulim*(1.0-zlim)/(16.0*sqr(Pi))/sqr(GeV);
   //jacobfactqr = 2.0*kb_photon.dot(kb_boson)*ulim/(16.0*sqr(Pi))/sqr(GeV);
   partrealqr = MatrRealQuark(_p_partona, _p_partonb, kphoton, kboson, kquark, scale_S, _iflagcq)/sumdipoleqr - 1.0;
   //txa = MatrRealQuark(_p_partona, _p_partonb, kphoton, kboson, kquark);
   pdffactqr = PDFratio(_xa, _xb, _muF2, 3, _iflagcq, _hadron_A, _hadron_B);
   //zcut = 1.0;
   // real piece for quark radiation collinear with final state photon
   if (zq<zcut) {
      smrealqrp = pdffactqr*jacobfactqr*(partrealqr*partdipoleqr*zq*zq -RDz_0)/zq;}
    else {
      smrealqrp = pdffactqr*jacobfactqr*partrealqr*partdipoleqr*zq;}
     //cout <<" zq="<<zq<<" uq="<<uq<<" zlim="<<zlim<<" kq_z0="<<kq_z0.e()/GeV<<" smrealqrp="<<smrealqrp<<"\n";

   //separatation cut:
   plp = kphoton;
   plv = kboson;
   plq = kquark;
   plp.boost(0.,0.,tva);
   plq.boost(0.,0.,tva);
   if (sqrt(sqr(plp.rapidity()-plq.rapidity()) + min(min(sqr(plp.phi()-plq.phi()), sqr(plp.phi()-plq.phi()+2.0*Pi)), sqr(plp.phi()-plq.phi()-2.0*Pi)))<Rcut) {
     if (plp.e()/(plp.e()+plq.e()) <fraccut) {

       //if ((kphoton.dot(kquark)/(kphoton.e()*kquark.e()))<thetacut) {
       //if (kphoton.e()/(kphoton.e()+kquark.e())<fraccut) {
       if (zq<zcut) {
         smrealqrp = smrealqrp - pdffactqr*jacobfactqr*((partrealqr+1.0)*partdipoleqr*zq - Rz_0/zq);
         //cout<<" separatation cut the real piece: E_frac="<<kphoton.e()/(kphoton.e()+kquark.e())<<" smrealqrp="<<smrealqrp<<"\n";
       }
       else {
	 smrealqrp = smrealqrp - pdffactqr*jacobfactqr*(partrealqr+1.0)*partdipoleqr*zq;
       }
     }
   }
   if (zi < fraccut) {
     if (zq<zcut) {
       smrealqrp = smrealqrp + pdffactqr*jacobfactqr*(partdipoleqr*zq - Dz_0/zq);
       //cout<<" separatation cut the subtracted term: z_frac="<<zq<<" smrealqrp="<<smrealqrp<<"\n";
     }
     else {
       smrealqrp = smrealqrp + pdffactqr*jacobfactqr*partdipoleqr*zq;
     }
   }

   // the zcut remnant term in realqr piece after extracting the soft photon singularity 
   ulim = 2.0*_p_photon.dot(_p_boson)*conu/(MV2+2.0*_p_photon.dot(_p_boson));
   uq = ulim - (ulim-1.0e-7)*_y2;
   kq_z0 = radkqr(_p_photon, _p_boson, MV2, 0.0, uq, phiq, zlim);
   kV_z0 = _p_boson -kq_z0 + _p_photon;
   kV_z0.rescaleMass();
   smbornqr = MatrQV(_p_partona, _p_partonb, kq_z0, kV_z0, chargeqr, MV2, _alphas, _iflagcq);
   if (_iflagcq==0){
     RDz_0 = 8.0*Pi*alphae* (((charge_i*charge_q)*_p_partonb.dot(kq_z0)/_p_partonb.dot(kb_photon)/kq_z0.dot(kb_photon)
			      + (charge_V*charge_i)*_p_partonb.dot(kV_z0)/_p_partonb.dot(kb_photon)/kV_z0.dot(kb_photon)
			      - (charge_V*charge_q)*kq_z0.dot(kV_z0)/kV_z0.dot(kb_photon)/kq_z0.dot(kb_photon))*smbornqr
                             - sqr(charge_q)/kq_z0.dot(kb_photon)*MatrQV(_p_partona, _p_partonb, kb_photon, kb_boson, chargeqr, MV2, _alphas, _iflagcq) )
                                *sqr(GeV);
   }
   else {
     RDz_0 = 8.0*Pi*alphae* (((charge_i*charge_q)*_p_partona.dot(kq_z0)/_p_partona.dot(kb_photon)/kq_z0.dot(kb_photon)
                              + (charge_V*charge_i)*_p_partona.dot(kV_z0)/_p_partona.dot(kb_photon)/kV_z0.dot(kb_photon)
                              - (charge_V*charge_q)*kq_z0.dot(kV_z0)/kV_z0.dot(kb_photon)/kq_z0.dot(kb_photon))*smbornqr
                             - sqr(charge_q)/kq_z0.dot(kb_photon)*MatrQV(_p_partona, _p_partonb, kb_photon, kb_boson, chargeqr, MV2, _alphas, _iflagcq) )
                                *sqr(GeV);
   }
   jacobfactqr = 2.0*kb_photon.dot(kb_boson)*ulim/(16.0*sqr(Pi))/sqr(GeV);
   smrealqrz_0 = pdffactqr*jacobfactqr*log(zcut)*RDz_0;   
   //smrealqrz_0 = 0.0;
   //separatation cut:
   //   if ((kphoton.dot(kquark)/(kphoton.e()*kquark.e()))<thetacut) {
   //     if (kphoton.e()/(kphoton.e()+kquark.e())<fraccut) {
   //      smrealqrz_0 = 0.0;
   //  }
   //}
     
   // in the singular region that quark radiation collinear with initial state gluon (Dgq):
   if (_iflagcq==0){
     fi = 0;
     xiabqr = (1.0 - 0.0000001)*(1.0 - _y1) +_xa*_y1;
     k_a = _p_partona/xiabqr;
     k_b = _p_partonb;
     kgluon = k_a;}
   else {
     fi = 1;
     xiabqr = (1.0 - 0.0000001)*(1.0 - _y1) +_xb*_y1;
     k_a = _p_partona;
     k_b = _p_partonb/xiabqr;
     kgluon = k_b;}
   viqr = (1.0-xiabqr)*_y2 + 0.0000001*(1.0 - _y2);
   // map to real phase space:  
   kquark =  radk(_p_partona, _p_partonb, xiabqr, viqr, phiq, fi); 
   kphoton = InvLortr(_p_partona, _p_partonb, xiabqr, kquark, _p_photon, fi);
   kboson = InvLortr(_p_partona, _p_partonb, xiabqr, kquark, _p_boson, fi);

   kgdotq = kgluon.dot(kquark);
   kpdotq = kphoton.dot(kquark);
   //Born phase space for Dpq in the denominator dipole sum:
   zi = 1.0/(1.0 + kquark.dot(kboson)/kphoton.dot(kquark+kboson));
   kbarp_ga = kphoton/zi;
   kbarp_V  = kboson + kquark -(1.0-zi)*kbarp_ga; 
   kbarp_V.rescaleMass();
   ulim = 2.0*(1.0-zi)*kbarp_ga.dot(kbarp_V)/(MV2+ 2.0*(1.0-zi)*kbarp_ga.dot(kbarp_V));

   partdipoleqr = Dipolegluqr(_p_partona, _p_partonb, _p_photon, _p_boson,charge0,charge1,MV2,kgdotq,_alphas,xiabqr);
   if (_iflagcq==0){
     jacobfactqr = 2.0*k_a.dot(k_b)*(1.0-xiabqr)*(1.0-_xa)/(16.0*sqr(Pi))/sqr(GeV);
     pdffactqr = PDFratio(_xa/xiabqr, _xa, _muF2, 4, _iflagcq, _hadron_A, _hadron_B);
     }
   else {
     jacobfactqr = 2.0*k_a.dot(k_b)*(1.0-xiabqr)*(1.0-_xb)/(16.0*sqr(Pi))/sqr(GeV);
     pdffactqr = PDFratio(_xb/xiabqr, _xb, _muF2, 4, _iflagcq, _hadron_A, _hadron_B);
      }
   sumdipoleqr = partdipoleqr + Dipolepqr(k_a, k_b, kbarp_ga, kbarp_V, chargeqr, MV2, kpdotq, zi, _alphas, _iflagcq);
   fluxfact = _ss/(k_a+k_b).m2();
   partrealqr = MatrRealQuark(k_a, k_b, kphoton, kboson, kquark, scale_S, _iflagcq)/sumdipoleqr - 1.0;
   //  real piece for quark radiation collinear with initial state quark
   smrealqrg = pdffactqr*jacobfactqr*partrealqr*partdipoleqr*fluxfact;
   /*  
   plp = kphoton;
   plp.boost(0.,0.,tva);
   if ((_iflagcq==0 && _xa>(1.0 - 1.0e-7)) || (_iflagcq==1 && _xb>(1.0 - 1.0e-7)) || plp.e()/(lastCuts().minKT(_photon)*cosh(plp.rapidity()))<1.0)
     smrealqrg = 0.0;
   */
   //   if ((_iflagcq==0 && _xa>(1.0 - 1.0e-7)) || (_iflagcq==1 && _xb>(1.0 - 1.0e-7)) || kphoton.e()/(lastCuts().minKT(_photon)*cosh(kphoton.rapidity()))<1.0)
   //if ((_iflagcq==0 && _xa>(1.0 - 1.0e-7)) || (_iflagcq==1 && _xb>(1.0 - 1.0e-7)) || kphoton.e()/GeV<20.0)
     //if ((_iflagcq==0 && _xa>(1.0 - 1.0e-7)) || (_iflagcq==1 && _xb>(1.0 - 1.0e-7)))
     //  { cout<< " xa="<<_xa<<" xb="<<_xb<<" xiab="<<xiab<<" vi="<<vi<<" xiabqr="<<xiabqr<<" viqr="<<viqr<<" smreal="<<smreal<<" smrealqr="
     //   <<(smrealqrg)/MatrBorn(_p_partona, _p_partonb, _p_photon, _p_boson, charge0, charge1, MV2)<<" cut="<<lastCuts().minKT(_photon)/GeV<<"\n";

   //separatation cut:
   //if ((kphoton.dot(kquark)/(kphoton.e()*kquark.e()))<thetacut) {
   //  if (kphoton.e()/(kphoton.e()+kquark.e())<fraccut) {
   plp = kphoton;
   plp.boost(0.,0.,tva);
   plq = kquark;
   plq.boost(0.,0.,tva);
   if (sqrt(sqr(plp.rapidity()-plq.rapidity())+min(sqr(plp.phi()-plq.phi()),sqr(2.0*Pi-abs(plp.phi()-plq.phi()))))<Rcut) {
     if (plp.e()/(plp.e()+plq.e())<fraccut) {
       //cout<<" cut in Dgq"<<"\n";
       smrealqrg = 0.0;
     }
   }

   //	 smrealqrg = 0.0;
   // the total subtracted quark radiation piece
   //smrealqr = (smrealqrp + smrealqrg + smrealqrz_0)/MatrBorn(_p_partona, _p_partonb, _p_photon, _p_boson, charge0, charge1, MV2);
   smrealqr = (smrealqrp + smrealqrg)/MatrBorn(_p_partona, _p_partonb, _p_photon, _p_boson, charge0, charge1, MV2); 
   
   //cout <<" zq="<<zq<<" uq="<<uq<<" zlim="<<zlim<<" kq_z0="<<kq_z0.e()/GeV<<" smrealqrp="<<smrealqrp<<" smrealqr="<<smrealqr<<" smcrqr="<<smcrqr<<" RDz_0="<<RDz_0<<"\n";
   //smrealqr = (smrealqrg)/MatrBorn(_p_partona, _p_partonb, _p_photon, _p_boson, charge0, charge1, MV2);
   //cout<< " smreal="<<smreal<< " smreala="<<smreala<< " smrealb="<<smrealb<< " smrealqr="<<smrealqr<<" smcrqr="<<smcrqr<<" xiab="<<xiab<<" vi="<<vi<<" smrealqrg="<<smrealqrg<< " smrealqrp="<<smrealqrp<<" whole NLO="<<alfsfact*(partep+partloop+partc)+ smreal+ smcrqr+ smrealqr<<"\n";
   return  1.0 +alfsfact*(partep+partloop+partc)+ smreal+ smcrqr+ smrealqr;
   //  if (abs(kpdotq/(uq*zq*_p_boson.dot(_p_photon)) - 1.0)>=0.000000001) {
   //   if (zq>= zcut ) {
   //  if (uq<=0.000000001 && abs((partdipoleqr - txa)/txa)>=0.01 )
   //    cout <<" zq="<<zq<< "  uq = "<< uq <<" dipole="<<partdipoleqr<<" sumdipole="<<sumdipoleqr<<" realME="<<txa<<" partrealqr="<<partrealqr
   //	 <<" ratio of kq.kp="<<plp.m2()/(2.0*uq*zq*kb_boson.dot(kb_photon)) <<" D/R="<<partdipoleqr/txa<<"\n";
  //  return  1.0;}
  //  else {
  //return  1.0 + alfsfact*(partep+partloop+partc) + smreal + smrealqr;
  // return  1.0 +alfsfact*(partep+partloop);
}

  //merge from VGammaHardGenerator
HardTreePtr MEPP2VGammaPowheg::
generateHardest(ShowerTreePtr tree, vector<ShowerInteraction::Type>) {
  Energy2 s_;
  Energy rs_;
  // get a pointer to the standard model object in the run
  static const tcHwSMPtr hwsm = dynamic_ptr_cast<tcHwSMPtr>(generator()->standardModel());
  // get the particles to be showered
  beams_.clear();
  partons_.clear();
  // find the incoming particles
  ShowerParticleVector incoming;
  quarkplus_ = true;
  vector<ShowerProgenitorPtr> particlesToShower;
  //generator()->log() << " Here is generateHardest \n";
  //progenitor particles are produced in z direction.
  for( map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator 
	 cit = tree->incomingLines().begin(); 
       cit != tree->incomingLines().end(); ++cit ) {
    incoming.push_back( cit->first->progenitor() );
    beams_.push_back( cit->first->beam() );
    partons_.push_back( cit->first->progenitor()->dataPtr() );
    // check that quark is along +ve z direction
    //if(cit->first->progenitor()->id() > 0 && 
    //   cit->first->progenitor()->momentum().z() < ZERO ) 
    //  quarkplus_ = false;
    particlesToShower.push_back( cit->first );
  }
  // find the outgoing particles
  tShowerParticlePtr boson,photon;
  for(map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
 	cjt=tree->outgoingLines().begin();cjt!=tree->outgoingLines().end();++cjt) {
    if(cjt->first->progenitor()->id()==ParticleID::gamma)
      photon = cjt->first->progenitor();
    else if(abs(cjt->first->progenitor()->id())==ParticleID::Wplus||
	    cjt->first->progenitor()->id()==ParticleID::Z0)
      boson  =  cjt->first->progenitor();
    particlesToShower.push_back(cjt->first);
  }
  assert(boson&&photon);
  // we are assuming quark first, swap order to ensure this
  // if antiquark first
  if(partons_[0]->id()<partons_[1]->id()) {
    swap(partons_[0],partons_[1]);
    swap(beams_[0],beams_[1]);
    swap(particlesToShower[0],particlesToShower[1]);
  }
      quarkplus_ = particlesToShower[0]->progenitor()->momentum().z() > ZERO;
  // Born variables
  // Rapidity of the photon
  _p_photon = photon->momentum();
  _p_boson = boson->momentum();
  photonRapidity_ = photon->momentum().rapidity();
  // Rapidity of the gauge boson
  bosonRapidity_  =  boson->momentum().rapidity();
  // pT of the photon
  photonpT_       = photon->momentum().perp();
  // Azimuth of the photon
  photonAzimuth_  = photon->momentum().phi();
  // gauge boson mass
  //bosonMass_      = boson->mass();
  MV2 = sqr(boson->mass()); 
  // mass of the boson/photon system
  systemMass_     = (_p_photon+_p_boson).m();
  // direction flip if needed
  /*
  if(!quarkplus_) {
    photonRapidity_ *= -1.;
    bosonRapidity_  *= -1.;
    photonAzimuth_  *= -1.;
    _p_photon.setZ((-1.)*_p_photon.z()); 
    _p_boson.setZ((-1.)*_p_boson.z());
  }
  */
  // ParticleData objects
  _photon = photon->dataPtr();
  _boson  = boson->dataPtr();
  // CMS energy of the hadron collision
  s_  = generator()->currentEvent()->primaryCollision()->m2();
  rs_ = sqrt(s_);
  // Borm momentum fractions
  double systemRapidity = (_p_photon + _p_boson).rapidity();
  //if(!quarkplus_) systemRapidity *= -1.;
  x_[0] = systemMass_*exp( systemRapidity)/rs_;          
  x_[1] = systemMass_*exp(-systemRapidity)/rs_;
  //generator()->log()<<"systemMass_="<<systemMass_/GeV<<" rs_*sqrt(x1x2)="<<rs_*sqrt(x_[0]*x_[1])/GeV;
  //incoming parton ParticleData object and momenta
  _partona= partons_[0];
  _partonb= partons_[1];
  _p_partona= Lorentz5Momentum(ZERO, ZERO, x_[0]*rs_/2.0, x_[0]*rs_/2.0, ZERO);
  _p_partonb= Lorentz5Momentum(ZERO, ZERO, -x_[1]*rs_/2.0, x_[1]*rs_/2.0, ZERO);
  if(!quarkplus_) {
    swap(_p_partona,_p_partonb);
    swap(x_[0],x_[1]);}

  //calculate the CKM matrix square
  int iu, id;
  iu= abs(_partona->id());
  id= abs(_partonb->id());
  if(iu%2!=0) swap(iu,id);
  iu = (iu-2)/2;
  id = (id-1)/2;

  ckm = hwsm->CKM(iu,id);
  // the third componet of isospin of the quark when vector boson is Z0
  I3quark = -double(abs(_partona->id())%2) + 0.5;
  // The charges of the quark and anti-quark
  charge0 = _partona->iCharge()/3.*_partona->id()/abs(_partona->id());
  charge1 = _partonb->iCharge()/3.*_partonb->id()/abs(_partonb->id());
  // constants
  Pi= Constants::pi;
  _CF= 4.0/3.0;
  _TR= 1.0/2.0;
  sin2w = hwsm->sin2ThetaW();
  alphae = hwsm->alphaEM(sqr(systemMass_));
  _muF2 = MV2;
  // generate the new configuration
  pTqqbar_ = -GeV;
  pTqg_    = -GeV;
  pTgqbar_ = -GeV;
  //generator()->log()   
  //throw InitException()<< " conserve?="<<(sqr((_p_partona+_p_partonb-_p_photon-_p_boson).e())+sqr((_p_partona+_p_partonb-_p_photon-_p_boson).z()))
  //   /sqr(GeV<<" (_p_photon+_p_boson).z=" <<(_p_photon+_p_boson).z()/GeV<<" (_p_photon+_p_boson).e=" <<(_p_photon+_p_boson).e()/GeV     
  //      << " x_[0]="<<x_[0]<< " x_[1]="<<x_[1]  <<"\n";     
  //    <<<<<<<<<<<<<<<<<<< Here we begin to generate the additional partonic radiation  >>>>>>>>>>>>>>
  // generate qq_bar to gluon Vgamma events
  
  //generator()->log() << " generate qq_bar to gluon Vgamma events" << "\n";

  //throw InitException() << " generate qq_bar to gluon Vgamma events";
  //			<< Exception::abortnow;
  generateQQbarG();
  // generate qg to q Vgamma events (1) and gq_bar to q_bar Vgamma events (2)

  //generator()->log() << " generate qg to q Vgamma events" << "\n";
  generateQGQ(1);

  //generator()->log() << " generate gq_bar to q_bar Vgamma events" << "\n";
  generateQGQ(2);

  if(pTqqbar_<ZERO && pTqg_<ZERO && pTgqbar_<ZERO) {
    for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
      particlesToShower[ix]->maximumpT(pTmin_,ShowerInteraction::QED);
      particlesToShower[ix]->maximumpT(pTmin_,ShowerInteraction::QCD);
    }
    return HardTreePtr();
  }
  // select the type of emission
  int emissionType=0;
  if      (pTqqbar_>pTqg_    && pTqqbar_>pTgqbar_ ) emissionType = 1;
  else if (pTqg_>pTqqbar_    && pTqg_>pTgqbar_    ) emissionType = 2;
  else if (pTgqbar_>pTqqbar_ && pTgqbar_>pTqg_    ) emissionType = 3;
  int iemit;
  ShowerParticleVector newparticles;
  Lorentz5Momentum pboson,pphoton;
  Energy pTEmit;
  // q qbar -> V gamma g
  if(emissionType==1) {
    newparticles.push_back(new_ptr(ShowerParticle(partons_[0]      ,false)));
    newparticles.push_back(new_ptr(ShowerParticle(partons_[1]      ,false)));
    newparticles.push_back(new_ptr(ShowerParticle(_gluon            , true)));
    newparticles[0]->set5Momentum(pQqqbar_);
    newparticles[1]->set5Momentum(pQbarqqbar_);
    newparticles[2]->set5Momentum(pGqqbar_);
    pboson  = pVqqbar_;
    pphoton = pGammaqqbar_; 
    iemit = pQqqbar_.z()/pGqqbar_.rapidity()>ZERO ? 0 : 1;
    pTEmit = pTqqbar_;
  }
  // q g    -> q V gamma
  else if(emissionType==2) {
    iemit=1;
    newparticles.push_back(new_ptr(ShowerParticle(partons_[0]      ,false)));
    newparticles.push_back(new_ptr(ShowerParticle(_gluon           ,false)));
    newparticles.push_back(new_ptr(ShowerParticle(partons_[1]->CC(), true)));
    newparticles[0]->set5Momentum(pQinqg_);
    newparticles[1]->set5Momentum(pGqg_);
    newparticles[2]->set5Momentum(pQoutqg_);
    pboson  = pVqg_;
    pphoton = pGammaqg_;
    pTEmit = pTqg_;
    QGQISR = QGQISR_qg;
  }
  // g qbar -> qbar V gamma
  else {
    iemit=0;
    newparticles.push_back(new_ptr(ShowerParticle(_gluon           ,false)));
    newparticles.push_back(new_ptr(ShowerParticle(partons_[1]      ,false)));
    newparticles.push_back(new_ptr(ShowerParticle(partons_[0]->CC(), true)));
    newparticles[0]->set5Momentum(pGgqbar_);
    newparticles[1]->set5Momentum(pQingqbar_);
    newparticles[2]->set5Momentum(pQoutgqbar_);
    pboson  = pVgqbar_;
    pphoton = pGammagqbar_;
    pTEmit  = pTgqbar_;
    QGQISR = QGQISR_gqbar;
  }

  double labbe;
  labbe = (x_[0]-x_[1])/(x_[0]+x_[1]);
  if(!quarkplus_) labbe *= -1.;
  Lorentz5Momentum ka, kb, kph, kbo, kemit, bpa, bpb, bpV, bpp;
  //look into the events cut out by the photon pT and rapidity cuts:
  ka = newparticles[0]->momentum();
  kb = newparticles[1]->momentum();
  kemit = newparticles[2]->momentum();
  kph = pphoton;
  kbo = pboson;
  ka.boost(0.,0.,-labbe);
  kemit.boost(0.,0.,-labbe);
  kb.boost(0.,0.,-labbe);
  kph.boost(0.,0.,-labbe);
  kbo.boost(0.,0.,-labbe);
  bpa = meMomenta()[0];
  bpb = meMomenta()[1];
  bpV = meMomenta()[2];
  bpp = meMomenta()[3];
  bpa.boost(0.,0.,labbe);
  bpb.boost(0.,0.,labbe);
  bpV.boost(0.,0.,labbe);
  bpp.boost(0.,0.,labbe);
  /*
  if (pphoton.perp()/(20.*GeV)<1. || abs(pphoton.rapidity())>2.7)
    generator()->log()<<" Investigate events by photon pT cut: "<<"     emission="<<emissionType<<"\n  pa = "<<newparticles[0]->momentum()/GeV<<"\n"
                      <<"  pb = "<<newparticles[1]->momentum()/GeV<<"\n"
                      <<" pemit="<<newparticles[2]->momentum()/GeV<<"\n"<<"  pV = "<<pboson/GeV<<"\n"<<" pph = "<<pphoton/GeV<<"\n"
                      <<" x0="<<x_[0]<<" x1="<<x_[1]<<" xa="<<lastX1()<<" quark="<<quarkplus_<<"   photon pT="<<pphoton.perp()/GeV
                      <<"  in lab="<<kph.perp()/GeV<<" y="<<pphoton.rapidity()<<"  Born pT="<<bpp.perp()/GeV<<" Born y="<<bpp.rapidity()<<"\n"
                      <<"   BornPa="<<meMomenta()[0]/GeV<<"       in lab"<<bpa/GeV<<"\n"<<"   BornPb="<<meMomenta()[1]/GeV<<"       in lab"<<bpb/GeV<<"\n"
                      <<"   BornPV="<<meMomenta()[2]/GeV<<"       in lab"<<bpV/GeV<<"\n"<<"  BornPph="<<meMomenta()[3]/GeV<<"       in lab"<<bpp/GeV<<"\n"
                      <<"  ka = "<<ka/GeV<<"\n"<<"  kb = "<<kb/GeV<<"\n"<<" kemit="<<kemit/GeV<<"\n"<<"  kV = "<<kbo/GeV<<"\n"<<" kph = "<<kph/GeV<<"\n";
  */
  /*
  if (abs(pphoton.rapidity())>2.7)
    generator()->log()<<" Investigate events by photon rapidity cut:\n  pa = "<<newparticles[0]->momentum()/GeV<<"\n"
                      <<"  pb = "<<newparticles[1]->momentum()/GeV<<"\n"
                      <<" pemit="<<newparticles[2]->momentum()/GeV<<"\n"<<"  pV = "<<pboson/GeV<<"\n"<<" pph = "<<pphoton/GeV<<"\n"
                      <<" x0="<<x_[0]<<" x1="<<x_[1]<<" xa="<<lastX1()<<" quark="<<quarkplus_<<"   photon pT="<<pphoton.perp()/GeV
                      <<"  in lab="<<kph.perp()/GeV<<" y="<<pphoton.rapidity()<<"  Born pT="<<bpp.perp()/GeV<<" Born y="<<bpp.rapidity()<<"\n"
                      <<"   BornPa="<<meMomenta()[0]/GeV<<"       in lab"<<bpa/GeV<<"\n"<<"   BornPb="<<meMomenta()[1]/GeV<<"       in lab"<<bpb/GeV<<"\n"
                      <<"   BornPV="<<meMomenta()[2]/GeV<<"       in lab"<<bpV/GeV<<"\n"<<"  BornPph="<<meMomenta()[3]/GeV<<"       in lab"<<bpp/GeV<<"\n"
                      <<"  ka = "<<ka/GeV<<"\n"<<"  kb = "<<kb/GeV<<"\n"<<" kemit="<<kemit/GeV<<"\n"<<"  kV = "<<kbo/GeV<<"\n"<<" kph = "<<kph/GeV<<"\n";
  */

  //  QGQISR =1;
  generator()->log()<<" check the emission type: emissionType="<<emissionType<<" QGQISR="<<QGQISR<<"\n";
  // create the photon
  newparticles.push_back(new_ptr(ShowerParticle(photon->dataPtr(),true)));
  newparticles.back()->set5Momentum(pphoton);
  // create the boson
  newparticles.push_back(new_ptr(ShowerParticle(boson ->dataPtr(),true)));
  newparticles.back()->set5Momentum(pboson);
  // newparticle[5]
  Lorentz5Momentum poff;
  // ISR:
  if (emissionType==1 || QGQISR==1) {
    poff = newparticles[iemit]->momentum()-newparticles[2]->momentum();
    poff.rescaleMass();
    newparticles.push_back(new_ptr(ShowerParticle(partons_[iemit],false)));
    newparticles.back()->set5Momentum(poff);}
  //FSR
  else {
    poff = newparticles[2]->momentum() + newparticles[3]->momentum();
    poff.rescaleMass();
    newparticles.push_back(new_ptr(ShowerParticle(partons_[iemit]->CC(),false)));
    newparticles.back()->set5Momentum(poff);}


  
  // find the sudakov for the branching
  //not apply when we merge the Powheg shower class into the matrix element class in the new version of Herwig++ ??
  /* 
  SudakovPtr sudakov;
  tEvolverPtr Evolver_;
  Evolver_ = Evolver();
  // **********should treat it for different types of emission***************
  // try to consider ISR in quark/antiquark radiation first:
  // ISR:
  if (emissionType==1 || QGQISR==1) {
  BranchingList branchings=Evolver_->splittingGenerator()->initialStateBranchings();
  long index = abs(partons_[iemit]->id());
  IdList br(3);
  // types of particle in the branching  
  // *********are these types of branching br[] assigned correctly?  Only apply for the ISR****************
  // gqbar to qar V gamma
  br[0]=newparticles[iemit]->id();
  br[1]=newparticles[  5  ]->id();
  br[2]=newparticles[  2  ]->id();
  if(emissionType==1) { //q qbar to g v gamma
    br[0]=abs(br[0]);
    br[1]=abs(br[1]);
  }
  else if(emissionType==2) { //q g to q V gamma
    br[1]=-br[1];
    br[2]=-br[2];
  }
  // **********should treat it for different types of emission?***************
  for(BranchingList::const_iterator cit = branchings.lower_bound(index); 
      cit != branchings.upper_bound(index); ++cit ) {
    IdList ids = cit->second.second;
    if(ids[0]==br[0]&&ids[1]==br[1]&&ids[2]==br[2]) { 
      sudakov=cit->second.first;
      break;
    }
  }
  }
  // FSR:
  else {
  BranchingList branchings=Evolver_->splittingGenerator()->finalStateBranchings();
  long index = abs(partons_[iemit]->CC()->id());
  IdList brp(3);
  brp[0]=index;
  brp[1]=index;
  brp[2]=newparticles[3]->id();
   for(BranchingList::const_iterator cit = branchings.lower_bound(index); 
      cit != branchings.upper_bound(index); ++cit ) {
    IdList ids = cit->second.second;
    if(ids[0]==brp[0]&&ids[1]==brp[1]&&ids[2]==brp[2]) { 
      sudakov=cit->second.first;
      break;
    }
  }
  }

  if(!sudakov) throw Exception() << "Can't find Sudakov for the hard emission in "
				 << "VGammaHardGenerator::generateHardest()" 
				 << Exception::runerror;
 
  */ 

  vector<HardBranchingPtr> nasonin,nasonhard;
  // ISR:
  if (emissionType==1 || QGQISR==1){  
  // create the branchings for the incoming particles
  nasonin.push_back(new_ptr(HardBranching(newparticles[0], SudakovPtr(),
					  //iemit==0 ? sudakov : SudakovPtr(),
					  HardBranchingPtr(),HardBranching::Incoming)));
  nasonin.push_back(new_ptr(HardBranching(newparticles[1], SudakovPtr(),
					  //iemit==1 ? sudakov : SudakovPtr(),
					  HardBranchingPtr(),HardBranching::Incoming)));
  // create the branching for the emitted jet
  nasonin[iemit]->addChild(new_ptr(HardBranching(newparticles[2],SudakovPtr(),
						 nasonin[iemit],
						 HardBranching::Outgoing)));
  // intermediate IS particle
  nasonhard.push_back(new_ptr(HardBranching(newparticles[5],SudakovPtr(),
					    nasonin[iemit],HardBranching::Incoming))); //   Only apply for the ISR****************
  nasonin[iemit]->addChild(nasonhard.back());
  // set the colour partners
  nasonhard.back()->colourPartner(nasonin[iemit==0 ? 1 : 0]);
  nasonin[iemit==0 ? 1 : 0]->colourPartner(nasonhard.back());
  // add other particle
  nasonhard.push_back(nasonin[iemit==0 ? 1 : 0]);
  // outgoing photon
  nasonhard.push_back(new_ptr(HardBranching(newparticles[3],SudakovPtr(),
					    HardBranchingPtr(),HardBranching::Outgoing)));
  //outgoing boson
  nasonhard.push_back(new_ptr(HardBranching(newparticles[4],SudakovPtr(),
					    HardBranchingPtr(),HardBranching::Outgoing)));
  }
  // FSR:
  else {
  // create the branchings for the incoming particles
  nasonin.push_back(new_ptr(HardBranching(newparticles[0], SudakovPtr(),
					  HardBranchingPtr(),HardBranching::Incoming)));
  nasonin.push_back(new_ptr(HardBranching(newparticles[1], SudakovPtr(),
					  HardBranchingPtr(),HardBranching::Incoming)));

  // try to add the initial state hard tree to nansonhard first:
  nasonhard.push_back(nasonin[0]); // should we add the initial state partons in nasonhard?
  nasonhard.push_back(nasonin[1]); // should we add the initial state partons in nasonhard?
  //nasonhard.push_back(new_ptr(HardBranching(newparticles[4],SudakovPtr(),
  //                                          HardBranchingPtr(),HardBranching::Outgoing)));

  // intermediate FS particle
  nasonhard.push_back(new_ptr(HardBranching(newparticles[5], SudakovPtr(),  //sudakov,
					    HardBranchingPtr(),HardBranching::Outgoing))); 
  nasonhard.back()->addChild(new_ptr(HardBranching(newparticles[2],SudakovPtr(),
						   nasonhard[2],HardBranching::Outgoing)));
  /*
  nasonhard.push_back(new_ptr(HardBranching(newparticles[3],SudakovPtr(),
					    nasonhard[0],HardBranching::Outgoing)));
  nasonhard[0]->addChild(nasonhard.back());
  */
  nasonhard.back()->addChild(new_ptr(HardBranching(newparticles[3],SudakovPtr(),
						   nasonhard[2],HardBranching::Outgoing)));
  //try to set the colour partners, suspect it relates to "Failed to make colour connections in PartnerFinder::setQCDInitialEvolutionScales":                                                
  nasonhard.back()->colourPartner(nasonhard[iemit==0 ? 1 : 0]);
  //nasonhard[iemit==0 ? 1 : 0]->colourPartner(nasonhard[iemit]);
  nasonhard[iemit==0 ? 1 : 0]->colourPartner(nasonhard.back());
  //nasonhard.back()->colourPartner(nasonhard[iemit]);
  //nasonhard[iemit]->colourPartner(nasonhard.back());
  nasonhard[iemit==0 ? 1 : 0]->colourPartner(nasonhard[iemit]);
  nasonhard[iemit]->colourPartner(nasonhard[iemit==0 ? 1 : 0]);

  
  //nasonhard.push_back(nasonin[0]); // should we add the initial state partons in nasonhard? 
  //nasonhard.push_back(nasonin[1]); // should we add the initial state partons in nasonhard? 
  nasonhard.push_back(new_ptr(HardBranching(newparticles[4],SudakovPtr(),
  				    HardBranchingPtr(),HardBranching::Outgoing))); 
  
  // we set the colour lines later: 
  /* 
  // set the colour line: two color lines:   initial q/qbar-gluon and gluon-q/qbar radiation
    ColinePtr cline1(new_ptr(ColourLine())), cline2(new_ptr(ColourLine()));
    bool antiq=newparticles[iemit==0 ? 1 : 0]->dataPtr()->iColour()!=PDT::Colour3;
    cline1->addColoured(newparticles[iemit==0 ? 1 : 0], antiq);
    cline1->addColoured(newparticles[iemit], !antiq);
    cline2->addColoured(newparticles[iemit], antiq); 
    cline2->addColoured(newparticles[2], antiq);
  */
  }

  // make the tree
  // try to separate the ISR (QCD) emission and FSR (QED) cases:
  HardTreePtr hardTree;
  if (emissionType==1 || QGQISR==1)
    hardTree=new_ptr(HardTree(nasonhard,nasonin,ShowerInteraction::QCD));
  else 
    hardTree=new_ptr(HardTree(nasonhard,nasonin,ShowerInteraction::QED));
  // connect the ShowerParticles with the branchings
  // and set the maximum pt for the radiation
  set<HardBranchingPtr> hard=hardTree->branchings();
  for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
    if( pTEmit < pTmin_ ) {
      particlesToShower[ix]->maximumpT(pTmin_,ShowerInteraction::QED);
      particlesToShower[ix]->maximumpT(pTmin_,ShowerInteraction::QCD);
    }
    else {
      particlesToShower[ix]->maximumpT(pTEmit,ShowerInteraction::QED);
      particlesToShower[ix]->maximumpT(pTEmit,ShowerInteraction::QCD);
    }
    for(set<HardBranchingPtr>::const_iterator mit=hard.begin();
	mit!=hard.end();++mit) {
      if(particlesToShower[ix]->progenitor()->id()==(*mit)->branchingParticle()->id()&&
	 (( (**mit).status()==HardBranching::Incoming &&
	    !particlesToShower[ix]->progenitor()->isFinalState())||
	  ( (**mit).status()==HardBranching::Outgoing&&
	    particlesToShower[ix]->progenitor()->isFinalState()))) {
	hardTree->connect(particlesToShower[ix]->progenitor(),*mit);
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

  //generator()->log()<<" ****after connecting the ShowerParticles with the branchings\n";
  // set the colour line:color lines:   initial q/qbar-gluon newline  and gluon-q/qbar radiation cline1 & cline2                                                                     
  ColinePtr cline1(new_ptr(ColourLine())), cline2(new_ptr(ColourLine()));
  ColinePtr newline=new_ptr(ColourLine());
  if (emissionType==1 || QGQISR==1){
    for(set<HardBranchingPtr>::const_iterator cit=hardTree->branchings().begin();
      cit!=hardTree->branchings().end();++cit) {
    if((**cit).branchingParticle()->dataPtr()->iColour()==PDT::Colour3)
      newline->addColoured((**cit).branchingParticle());
    else if((**cit).branchingParticle()->dataPtr()->iColour()==PDT::Colour3bar)
      newline->addAntiColoured((**cit).branchingParticle());
    }
  }
  else {
    for(set<HardBranchingPtr>::const_iterator cit=hardTree->branchings().begin();
	cit!=hardTree->branchings().end();++cit) {
      if((**cit).branchingParticle()->dataPtr()->iColour()==PDT::Colour8) {
	if((**cit).status()==HardBranching::Incoming) {
	  cline1->addColoured    ((**cit).branchingParticle());
	  cline2->addAntiColoured((**cit).branchingParticle());
	}
       
      }
      else if((**cit).branchingParticle()->dataPtr()->iColour()==PDT::Colour3) {
	if((**cit).status()==HardBranching::Incoming) 
	  cline2->addColoured((**cit).branchingParticle());
	else 
	  cline1->addColoured((**cit).branchingParticle());
      }
      else if((**cit).branchingParticle()->dataPtr()->iColour()==PDT::Colour3bar) {
	if((**cit).status()==HardBranching::Incoming) 
	  cline1->addAntiColoured((**cit).branchingParticle());
	else 
	  cline2->addAntiColoured((**cit).branchingParticle());
      }
    }
  }

  ShowerParticleVector particles;
  for(set<HardBranchingPtr>::iterator cit=hardTree->branchings().begin();
      cit!=hardTree->branchings().end();++cit) {
    particles.push_back((*cit)->branchingParticle());
    if (particles.back()->partner())
      generator()->log() << " ***hardTree particle: p="<<particles.back()->data().id()<<"  partner="<<particles.back()->partner()->data().id()<<" cline="<<particles.back()->colourLine()<<" InitScale="<<particles.back()->evolutionScale()/GeV<<"\n";
    else 
      generator()->log() << " ***hardTree particle: p="<<particles.back()->data().id()<<" no partner"<<" cline="<<particles.back()->colourLine()<<" InitScale="<<particles.back()->evolutionScale()/GeV<<"\n";
  }

  generator()->log() <<"  check the nasonin scales: "<<nasonin[0]->scale()/GeV<<" "<< nasonin[1]->scale()/GeV<<" nasonhard: "<<nasonhard[0]->scale()/GeV<<
    " "<<nasonhard[1]->scale()/GeV<<" "<<nasonhard[2]->scale()/GeV<<" "<<nasonhard[3]->scale()/GeV<<"\n";
  // here we try to set the initial evolution scale of the initial partons for the FSR by hand:
  pair<Energy,Energy> initscalepair;
  initscalepair = pair<Energy,Energy>(ZERO,ZERO);
  ShowerParticleVector initpart;
  if  (!(emissionType==1 || QGQISR==1)){
  for(set<HardBranchingPtr>::iterator cit=hardTree->branchings().begin();
      cit!=hardTree->branchings().end();++cit) {
    if ((**cit).branchingParticle()->dataPtr()->iColour()==PDT::Colour8 && (**cit).status()==HardBranching::Incoming) {
	initpart.push_back((*cit)->branchingParticle());
      }
    else if(((**cit).branchingParticle()->dataPtr()->iColour()==PDT::Colour3 || (**cit).branchingParticle()->dataPtr()->iColour()==PDT::Colour3bar) && (**cit).status()==HardBranching::Incoming) {
      initpart.push_back((*cit)->branchingParticle());
    }
  }
    Lorentz5Momentum psum;
    psum = initpart[0]->momentum() + initpart[1]->momentum();
    psum.boost(psum.findBoostToCM());
    int InitConditions (3);
    Energy Qinit = sqrt(psum.m2());
    if(InitConditions==1) {
      initscalepair = pair<Energy,Energy>(sqrt(2.0)*Qinit,sqrt(0.5)*Qinit);
    } else if(InitConditions==2) {
      initscalepair = pair<Energy,Energy>(sqrt(0.5)*Qinit,sqrt(2.0)*Qinit);
    } else {
      initscalepair = pair<Energy,Energy>(Qinit,Qinit);
    }
    initpart[0]->setEvolutionScale(initscalepair.first);
    initpart[1]->setEvolutionScale(initscalepair.second);
  }
  
  generator()->log() << " ***hardTree particle: p0="<<particles[0]->data().id()<<" scale="<<particles[0]->evolutionScale()/GeV<<//"\n";<<particles[0]->colourLine()<<
    ";  p1="<<particles[1]->data().id()<<" scale="<<particles[1]->evolutionScale()/GeV<<" initscale="<<initscalepair.second/GeV<<"\n";
  //";  p2="<<particles[2]->data().id()<<" "<<particles[2]->partner()->data().id()<<particles[2]->colourLine()<<
  //";  p3="<<particles[3]->data().id()<<" "<<particles[3]->partner()->data().id()<<particles[3]->colourLine()<<
    //";  p4="<<particles[4]->data().id()<<" "<<particles[4]->partner()->data().id()<<particles[4]->colourLine()<<
    //";  p5="<<particles[5]->data().id()<<" "<<particles[5]->partner()->data().id()<<//particles[5]->colourLine()<<
    //"\n";
  
  //not apply when we merge the Powheg shower class into the matrix element class in the new version of Herwig++ ??
  /*
  Herwig::Evolver::Evolver().showerModel()->partnerFinder()->
     setInitialEvolutionScales(particles,true,ShowerInteraction::QCD,true);
  
  // calculate the shower variables
  Evolver_->showerModel()->kinematicsReconstructor()->
    deconstructHardJets(hardTree,Evolver_,ShowerInteraction::QCD);
  for(set<HardBranchingPtr>::const_iterator cit=hardTree->branchings().begin();
      cit!=hardTree->branchings().end();++cit) {
    for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
      if((((**cit).status()==HardBranching::Incoming && 
	   !particlesToShower[ix]->progenitor()->isFinalState())||
	  ((**cit).status()==HardBranching::Outgoing &&
	   particlesToShower[ix]->progenitor()->isFinalState()))&&
	 particlesToShower[ix]->progenitor()->id()==(**cit).branchingParticle()->id()) {
	particlesToShower[ix]->progenitor()->set5Momentum((**cit).showerMomentum());
      }
    }
  }
  */

  if(hardTree->interaction()==ShowerInteraction::QCD)
    generator()->log() << "INTERACTION = QCD\n";
  else
    generator()->log() << "INTERACTION = QED\n";
  for(set<HardBranchingPtr>::const_iterator cit=hardTree->branchings().begin();
      cit!=hardTree->branchings().end();++cit) {
    generator()->log() << "testing partners"
		       << *cit << " " << (**cit).colourPartner() << "\n";
  }
  generator()->log() << *hardTree << "\n";




  return hardTree;
}

// generate qq_bar to gluon Vgamma events
  // merge from VGammaHardGenerator
void MEPP2VGammaPowheg::generateQQbarG() {
  Energy2 sHat, scale, kadotg, kbdotg;
  Energy rsHat, rsub;
  // storage of pt and rapidity of the gluon
  Energy pTa, pTb, pTV, mTV;
  double yj,x1,x2,wgta, wgtb, estimatefact, pdffact, jacobfact, partdipole, sumdipole;
  // radiation phase space variables:
  double phi, xiab, vi, rho, rhomx;
  //  double charge0,charge1,sh,shr,ratea,rateb,va,vb,txa,txb,tva,tvb;
  Lorentz5Momentum k_a_a, k_b_a, k_ga_a, k_V_a, k_glu_a, k_a_b, k_b_b, k_ga_b, k_V_b, k_glu_b, kbarp_ga, kbarp_V;
  // limits on the rapidity of the jet
  double minyj = -8.0,maxyj = 8.0, myjpr, exymin, exymax;
  //cout<<" start generateQQbarG \n";
  // PDFs for the Born cross section
  double pdf[4];
  int emitflag= 0;
  pdf[0]=beams_[0]->pdf()->xfx(beams_[0],partons_[0],_muF2,x_[0]);
  pdf[1]=beams_[1]->pdf()->xfx(beams_[1],partons_[1],_muF2,x_[1]);
  // prefactor for veto algorithm
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //double c = alphaS_->overestimateValue()/Constants::twopi*
  //    qqgFactor_*(maxyj-minyj)/(power_-1.);
  double c = alphaS_->overestimateValue()/Constants::twopi*qqgFactor_*(maxyj-minyj);  

  //   for the initial state singularity with parton_a:
  rhomx = (1.0 -x_[0])/(2.0*sqrt(x_[0]));
  rho = rhomx;
  pTa = rho*systemMass_;
  rsub = systemMass_;
  //throw InitException() <<" for the initial state singularity with parton_a \n";
  do {
    // generate pT
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //pTa = rsub*pow(pow(double(pTa/rsub),1.-power_)-log(UseRandom::rnd())/c,1./(1.-power_));
    pTa  *= pow(UseRandom::rnd(),1./c);
    // generate rapidity of the jet
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    minyj = -8.0;
    maxyj = 8.0;
    rho = pTa/systemMass_;
    /*
    exymin= (rhomx/rho - sqrt(sqr(rhomx/rho)-1.0))/sqrt(x_[0]);
    //exymin= (rhomx/rho - sqrt(sqr(rhomx/rho)-1.0))*sqrt(x_[0]);
    if (exymin<rho) exymin = rho;
    myjpr = log(exymin*sqrt(x_[0]/x_[1]));
    //myjpr = log(exymin/sqrt(x_[0]/x_[1]));
    if (myjpr>minyj) minyj = myjpr;
    exymax= (rhomx/rho + sqrt(sqr(rhomx/rho)-1.0))/sqrt(x_[0]);
    //exymax= (rhomx/rho + sqrt(sqr(rhomx/rho)-1.0))*sqrt(x_[0]);
    myjpr = log(exymax*sqrt(x_[0]/x_[1]));
    //myjpr = log(exymax/sqrt(x_[0]/x_[1]));
    if (myjpr<maxyj) maxyj = myjpr; 
    */
    yj = UseRandom::rnd()*(maxyj-minyj)+ minyj;
    // generate phi
    phi = 2.0*Pi*UseRandom::rnd();
    // transform to CS variables:
    vi = rho *exp(-yj) *sqrt(x_[0]/x_[1]);
    xiab = vi*(1.0-vi)/(sqr(rho)+vi);

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    if (xiab>1.|| xiab<x_[0]) continue;
    if (vi<0. || vi>(1.-xiab)) continue;

    //generator()->log() <<" the rapidity of gluon_a: y="<<yj<<" exymin="<<exymin<<" rho="<<rho<<" exymax="<<exymax<<" myjpr="<<myjpr<<" pTa="<<pTa/GeV<<"\n";    
    /*
    // pT of the W/Z
    pTV = sqrt(sqr(photonpT_*sin(photonAzimuth_)+pT*sin(phi))+
	       sqr(photonpT_*cos(photonAzimuth_)+pT*cos(phi)));
    mTV = sqrt(sqr(pTV)+sqr(bosonMass_));
    // azimuth of W/Z
    phiV = Constants::pi+atan2(photonpT_*sin(photonAzimuth_)+pT*sin(phi),
			       photonpT_*cos(photonAzimuth_)+pT*cos(phi));
    // calculate x_1 and x_2
    x1 = 0.5/rs_*(photonpT_*exp( photonRapidity_)+mTV*exp( bosonRapidity_)+pT*exp( yj));   //why there is a 0.5 factor here?
    x2 = 0.5/rs_*(photonpT_*exp(-photonRapidity_)+mTV*exp(-bosonRapidity_)+pT*exp(-yj));
    */
    // calculate x_1 and x_2
    x1 = x_[0]/xiab;
    x2 = x_[1];
    //    throw InitException() << " x1="<< x1<< " xiab="<<xiab<<" vi="<<vi<<Exception::abortnow;
    // sHat
    sHat  = sqr(systemMass_)/xiab; //sHat  = x1*x2*s_;
    rsHat = sqrt(sHat);

    if(!quarkplus_) yj *= -1.;
    // real radiation phase space:
    k_glu_a=     Lorentz5Momentum (pTa *cos(phi ),pTa *sin(phi ),
				   pTa *sinh(       yj     ),
				   pTa *cosh(      yj      ),ZERO);
    //throw InitException() << "x_[0]="<<x_[0]<<" x1="<< x1<<" vi="<<vi<<" xiab="<<xiab<<" yj="<<yj<<" rho="<<rho<<" rho_0="<<(1.0 -x_[0])/(2.0*sqrt(x_[0]))
    //  <<" k_glu_a.e="<<k_glu_a.e()/GeV<<" k_glu_a.z="<<k_glu_a.z()/GeV<<"\n"
    //                      << Exception::abortnow;
    // Suppose we have construct the pT and C-S variables, and then we can generate the n+1 events according to the R/B ratios 
    // of the two singular regions using the highest-kT bid technique: 
    k_a_a = _p_partona/xiab;
    k_b_a = _p_partonb;
    //throw InitException() << " x1="<< x1<<"x_[0]="<<x_[0]<<" vi="<<vi<<" xiab="<<xiab<<" k_glu_a.e="<<k_glu_a.e()/GeV<<" k_a_a.e="<<k_a_a.e()/GeV<<
    //  " k_b_a.e="<<k_b_a.e()/GeV<< Exception::abortnow;
    k_ga_a = InvLortr(_p_partona, _p_partonb, xiab, k_glu_a, _p_photon, 0);
    k_V_a =  InvLortr(_p_partona, _p_partonb, xiab, k_glu_a, _p_boson, 0);

    if (abs(xiab -(1.-k_glu_a.dot(k_a_a+k_b_a)/k_a_a.dot(k_b_a)))>0.0001 || abs(vi -(k_glu_a.dot(k_a_a)/k_a_a.dot(k_b_a)))> 0.0001) {
      generator()->log() <<" At gluon_a singular region CS-mapping is wrong:  "<<" xiab="<<1.-k_glu_a.dot(k_a_a+k_b_a)/k_a_a.dot(k_b_a)<<" ?= "<<xiab<<" vi="<<k_glu_a.dot(k_a_a)/k_a_a.dot(k_b_a)<<" ?= "<<vi<<" xlim="<<x_[0]
			 <<" vlim="<<k_glu_a.dot(k_a_a+k_b_a)/k_a_a.dot(k_b_a)<<" quarkplus="<<quarkplus_<<"\n";
      continue;}

    //    throw InitException() << " x1="<< x1<<" k_ga_a.e="<<k_ga_a.e()/GeV<<" k_V_a.e="<<k_V_a.e()/GeV
    //              <<" conserve?="<<k_glu_a.isNear(k_a_a+k_b_a-k_ga_a-k_V_a, 0.000000000001)<<"\n"
    //                     << Exception::abortnow;    
    kadotg = k_a_a.dot(k_glu_a);
    //  kadotg = vi*k_a_a.dot(k_b_a);
    kbdotg = k_b_a.dot(k_glu_a);

    kbarp_ga = Lortr(k_a_a, k_b_a, xiab, k_glu_a, k_ga_a, 1);
    kbarp_V  = Lortr(k_a_a, k_b_a, xiab, k_glu_a, k_V_a, 1);
    //throw InitException() << " x1="<< x1<<" k_ga_a.e="<<k_ga_a.e()/GeV<<" k_V_a.e="<<k_V_a.e()/GeV<<" kbarp_ga.e="<<kbarp_ga.e()/GeV<<
    //  " kbarp_V.e="<<kbarp_V.e()/GeV<<" conserve?="<<k_glu_a.isNear(k_a_a+k_b_a-k_ga_a-k_V_a, 0.000000000001)<<"\n"
    //    << Exception::abortnow;
    scale = sqr(systemMass_)+ sqr(pTa);
    partdipole = Dipole(_p_partona, _p_partonb, _p_photon, _p_boson, charge0, charge1, MV2, kadotg, alphaS_->value(scale), xiab);
    sumdipole = partdipole +Dipole(k_a_a, xiab*k_b_a, kbarp_ga, kbarp_V, charge0, charge1, MV2, kbdotg, alphaS_->value(scale), xiab);

    //jacobfact = sHat/(16.0*sqr(Pi)*xiab)/(1.0-vi)*sqr(xiab)/sqr(systemMass_);
    //jacobfact = sHat/(16.0*sqr(Pi))/(1.0-vi)*sqr(xiab)/systemMass_;
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    jacobfact = xiab/(16.0*sqr(Pi))/(1.0-vi);
    //estimatefact = alphaS_->overestimateValue()/(4.0*Pi)*
    //     qqgFactor_/sqr(pTa/GeV)*pow(rsub/pTa,(power_-1.0))*(2.*8./(maxyj-minyj));
    estimatefact = _alphas/alphaS_->ratio(scale)/(4.0*Pi)*qqgFactor_/sqr(pTa/GeV);
    pdffact = PDFratio(x1, x_[0], scale, 0, 0, beams_[0], beams_[1]);
    // now for the weight
    // first kinematical factors
    //wgt = 0.5*pow<4,1>(systemMass_)/sqr(sHat);
    // pdf bit
    //Energy2 scale = sqr(systemMass_)+sqr(pT);
    pdf[2]=beams_[0]->pdf()->xfx(beams_[0],partons_[0],scale,x1);
    pdf[3]=beams_[1]->pdf()->xfx(beams_[1],partons_[1],scale,x2);
    if(pdf[0]<=0.||pdf[1]<=0.||pdf[2]<=0.||pdf[3]<=0.) {
      wgta=0.;
      continue;
    }
    //throw InitException() <<" about to QQbarGratio() at gluon_a singular region\n"
    //		  << Exception::abortnow;
    // final bit
    wgta = QQbarGratio(k_a_a, k_b_a, k_ga_a, k_V_a, k_glu_a, scale)*partdipole/sumdipole*pdffact/estimatefact*jacobfact;
    // correct the y's Jacobian:
    //wgta *= 2.*8.0/(maxyj-minyj)*1000000.;
    //generator()->log() <<" QQbarG weight at gluon_a singular region:  "<< wgta<<"\n";
    if(wgta>1.) generator()->log() << "Weight greater than one for gluon_a emission, Weight = " << wgta << "\n";
    if (UseRandom::rnd()<wgta) break;
  }
  //  while(UseRandom::rnd()>wgta && pTa>pTmin_);
  while(pTa>pTmin_);
  
  if(pTa<pTmin_) pTa=-GeV;

  //   for the initial state singularity with parton_b:
  rhomx = (1.0 -x_[1])/(2.0*sqrt(x_[1]));
  rho = rhomx;
  pTb = rho*systemMass_;
  rsub = systemMass_;
  do {
    // generate pT
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //pTb = rsub*pow(pow(double(pTb/rsub),1.-power_)-log(UseRandom::rnd())/c,1./(1.-power_));
    pTb  *= pow(UseRandom::rnd(),1./c);
    // generate rapidity of the jet
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    minyj = -8.0;
    maxyj = 8.0;
    rho = pTb/systemMass_;
    /*
    exymin= (rhomx/rho - sqrt(sqr(rhomx/rho)-1.0))*sqrt(x_[1]);
    myjpr = log(exymin*sqrt(x_[0]/x_[1]));
    if (myjpr>minyj) minyj = myjpr;
    exymax= (rhomx/rho + sqrt(sqr(rhomx/rho)-1.0))*sqrt(x_[1]);
    if (exymax>(1.0/rho)) exymax = 1.0/rho;
    myjpr = log(exymax*sqrt(x_[0]/x_[1]));
    if (myjpr<maxyj) maxyj = myjpr;
    */   
    yj = UseRandom::rnd()*(maxyj-minyj)+ minyj;
    // generate phi
    phi = 2.0*Pi*UseRandom::rnd();
    // transform to CS variables:
    //rho = pTb/systemMass_;
    vi = rho *exp(yj) *sqrt(x_[1]/x_[0]);
    xiab = vi*(1.0-vi)/(sqr(rho)+vi);

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&                                                                                                           
    if (xiab>1.|| xiab<x_[0]) continue;
    if (vi<0. || vi>(1.-xiab)) continue;

    //generator()->log() <<" the rapidity of gluon_b: y="<<yj<<" exymin="<<exymin<<" 1/rho="<<1.0/rho<<" exymax="<<exymax<<" myjpr="<<myjpr<<" pTb="<<pTb/GeV<<"\n";    
    /*
    // pT of the W/Z
    pTV = sqrt(sqr(photonpT_*sin(photonAzimuth_)+pT*sin(phi))+
	       sqr(photonpT_*cos(photonAzimuth_)+pT*cos(phi)));
    mTV = sqrt(sqr(pTV)+sqr(bosonMass_));
    // azimuth of W/Z
    phiV = Constants::pi+atan2(photonpT_*sin(photonAzimuth_)+pT*sin(phi),
			       photonpT_*cos(photonAzimuth_)+pT*cos(phi));
    // calculate x_1 and x_2
    x1 = 0.5/rs_*(photonpT_*exp( photonRapidity_)+mTV*exp( bosonRapidity_)+pT*exp( yj));   //why there is a 0.5 factor here?
    x2 = 0.5/rs_*(photonpT_*exp(-photonRapidity_)+mTV*exp(-bosonRapidity_)+pT*exp(-yj));
    */
    // calculate x_1 and x_2
    x1 = x_[0];
    x2 = x_[1]/xiab;
    // sHat
    sHat  = sqr(systemMass_)/xiab; //sHat  = x1*x2*s_;
    rsHat = sqrt(sHat);

    if(!quarkplus_) yj *= -1.;
    // real radiation phase space:
    k_glu_b=     Lorentz5Momentum (pTb *cos(phi ),pTb *sin(phi ),
				   pTb *sinh(       yj     ),
				   pTb *cosh(      yj      ),ZERO);
    // Suppose we have construct the pT and C-S variables, and then we can generate the n+1 events according to the R/B ratios 
    // of the two singular regions using the highest-kT bid technique: 
    k_a_b = _p_partona;  
    k_b_b = _p_partonb/xiab;
    k_ga_b = InvLortr(_p_partona, _p_partonb, xiab, k_glu_b, _p_photon, 1);
    k_V_b =  InvLortr(_p_partona, _p_partonb, xiab, k_glu_b, _p_boson, 1);

    if(abs(xiab-(1.-k_glu_b.dot(k_a_b+k_b_b)/k_a_b.dot(k_b_b)))>0.0001 || abs(vi-(k_glu_b.dot(k_b_b)/k_a_b.dot(k_b_b)))>0.0001) {
      generator()->log() <<" At gluon_b singular region CS-mapping is wrong:  "<<" xiab="<<1.-k_glu_b.dot(k_a_b+k_b_b)/k_a_b.dot(k_b_b)<<" ?= "<<xiab<<" vi="<<k_glu_b.dot(k_b_b)/k_a_b.dot(k_b_b)<<" ?= "<<vi
			 <<" xlim="<<x_[1]<<" vlim="<<1.-xiab<<" quarkplus="<<quarkplus_<<"\n";
      continue;}

    kadotg = k_a_b.dot(k_glu_b);
    //  kadotg = vi*k_a_a.dot(k_b_a);
    kbdotg = k_b_b.dot(k_glu_b);

    kbarp_ga = Lortr(k_a_b, k_b_b, xiab, k_glu_b, k_ga_b, 0);
    kbarp_V  = Lortr(k_a_b, k_b_b, xiab, k_glu_b, k_V_b, 0);

    scale = sqr(systemMass_)+ sqr(pTb);
    partdipole = Dipole(_p_partona, _p_partonb, _p_photon, _p_boson, charge0, charge1, MV2, kbdotg, alphaS_->value(scale), xiab);
    sumdipole = partdipole +Dipole(xiab*k_a_b, k_b_b, kbarp_ga, kbarp_V, charge0, charge1, MV2, kadotg, alphaS_->value(scale), xiab);

    //jacobfact = sHat/(16.0*sqr(Pi)*xiab)/(1.0-vi)*sqr(xiab)/sqr(systemMass_);
    //jacobfact = sHat/(16.0*sqr(Pi))/(1.0-vi)*sqr(xiab)/systemMass_;
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    jacobfact = xiab/(16.0*sqr(Pi))/(1.0-vi);
    //estimatefact = alphaS_->overestimateValue()/(4.0*Pi)*
    //   qqgFactor_/sqr(pTb/GeV)*pow(rsub/pTb,(power_-1.0))*(2.*8./(maxyj-minyj));
    estimatefact = _alphas/alphaS_->ratio(scale)/(4.0*Pi)*qqgFactor_/sqr(pTb/GeV);
    pdffact = PDFratio(x2, x_[1], scale, 1, 0, beams_[0], beams_[1]);
    // now for the weight
    // first kinematical factors
    //wgt = 0.5*pow<4,1>(systemMass_)/sqr(sHat);
    // pdf bit
    //Energy2 scale = sqr(systemMass_)+sqr(pT);
    pdf[2]=beams_[0]->pdf()->xfx(beams_[0],partons_[0],scale,x1);
    pdf[3]=beams_[1]->pdf()->xfx(beams_[1],partons_[1],scale,x2);
    if(pdf[0]<=0.||pdf[1]<=0.||pdf[2]<=0.||pdf[3]<=0.) {
      wgtb=0.;
      continue;
    }
    //throw InitException() <<" about to QQbarGratio() at gluon_b singular region\n"
    //                     << Exception::abortnow; 
    //generator()->log() <<" about to QQbarGratio() at gluon_b singular region";
    // final bit
    wgtb = QQbarGratio(k_a_b, k_b_b, k_ga_b, k_V_b, k_glu_b, scale)*partdipole/sumdipole*pdffact/estimatefact*jacobfact;
    // correct the y's Jacobian:
    //wgtb *= 2.*8.0/(maxyj-minyj)*1000000.;
    //generator()->log() <<" QQbarG weight at gluon_b singular region:  "<< wgtb<<"\n";
    if(wgtb>1.) generator()->log() << "Weight greater than one for gluon_b emission, Weight = " << wgtb << "\n";

    if (UseRandom::rnd()<wgtb) break;
  }
  //  while(UseRandom::rnd()>wgtb && pTb>pTmin_);
  while(pTb>pTmin_);
  
  if(pTb<pTmin_) pTb=-GeV;
  //generator()->log() << " gluon radiation with Highest-pT-bid method: wgta="<<wgta<<" pTa="<<pTa/GeV<<" wgtb="<<wgtb<<" pTb="<<pTb/GeV<<
  //" sbar="<<systemMass_/GeV<<"\n";

  // select the real radiation event with Highest-pT-bid method
  if(pTa > pTb){
    pTqqbar_ = pTa;
    pGammaqqbar_ = k_ga_a;
    pVqqbar_     = k_V_a;
    pGqqbar_     = k_glu_a;
    pQqqbar_     = k_a_a;
    pQbarqqbar_  = k_b_a;
    emitflag = 0;
    generator()->log() <<" QQbarG weight at gluon_a singular region:  "<<wgta<<"\n";
    //generator()->log() <<" QQbarG weight at gluon_a singular region:  "<<wgta<<" xiab="<<1.-k_glu_a.dot(k_a_a+k_b_a)/k_a_a.dot(k_b_a)<<" vi="<<k_glu_a.dot(k_a_a)/k_a_a.dot(k_b_a)<<" xlim="<<x_[0]
    //		       <<" vlim="<<k_glu_a.dot(k_a_a+k_b_a)/k_a_a.dot(k_b_a)<<" quarkplus="<<quarkplus_<<"\n";
  }
  else{
    pTqqbar_ = pTb;
    pGammaqqbar_ = k_ga_b;
    pVqqbar_     = k_V_b;
    pGqqbar_     = k_glu_b;
    pQqqbar_     = k_a_b;
    pQbarqqbar_  = k_b_b;
    emitflag = 1;
    generator()->log() <<" QQbarG weight at gluon_b singular region:  "<<wgtb<<"\n";
    //generator()->log() <<" QQbarG weight at gluon_b singular region:  "<<wgtb<<" xiab="<<1.-k_glu_b.dot(k_a_b+k_b_b)/k_a_b.dot(k_b_b)<<" ?= "<<xiab<<" vi="<<k_glu_b.dot(k_b_b)/k_a_b.dot(k_b_b)<<" ?= "<<vi
    //		       <<" xlim="<<x_[1]<<" vlim="<<1.-xiab<<" quarkplus="<<quarkplus_<<"\n";
  }

  /*
  if(!quarkplus_) {
    LorentzRotation trans;
    trans.rotateX(Constants::pi);
    pGammaqqbar_.transform(trans);
    pVqqbar_    .transform(trans);
    pGqqbar_    .transform(trans);
    pQqqbar_    .transform(trans);
    pQbarqqbar_ .transform(trans);
  }
  */
  
  generator()->log() << " q qbar -> g V Gamma: testing new \n"
		     << "Photon    = " << pGammaqqbar_/GeV << "\n"
		     << "Boson     = " << pVqqbar_/GeV << "\n"
		     << "Gluon     = " << pGqqbar_/GeV << "\n"
		     << "Quark     = " << pQqqbar_/GeV << "\n"
		     << "Antiquark = " << pQbarqqbar_/GeV << "\n"
		     << "sum " 
		     << (pGammaqqbar_+pVqqbar_+pGqqbar_-pQqqbar_-pQbarqqbar_)/GeV << "\n";
  
  //ssHat_ = sHat;
  //maTV_ = mTV;
  pTparton_ = pTqqbar_;
  //phiparton_ = phi;
  //yparton_ = yj;
  //generator()->log() <<"                                emitflag"<<emitflag<<"\n";
}

  // merge from VGammaHardGenerator
double MEPP2VGammaPowheg::QQbarGratio(Lorentz5Momentum k_a, Lorentz5Momentum  k_b, Lorentz5Momentum k_ga, Lorentz5Momentum k_V, 
                                        Lorentz5Momentum k_glu, Energy2 scale) {
  double mreal, mborn, ratio;
  //Gluon radiation Real Matrix Element:
  mreal = MatrRealGluon(k_a, k_b, k_ga, k_V, k_glu, scale);
  //Born Matrix Element:
  mborn = MatrBorn(_p_partona, _p_partonb, _p_photon, _p_boson, charge0, charge1, MV2);

  ratio= mreal/mborn;
 
  return ratio;
}






  // generate qg to q Vgamma events (qgflag=1) and gq_bar to q_bar Vgamma events (qgflag=2)
  // merge from VGammaHardGenerator
void MEPP2VGammaPowheg::generateQGQ(int qgflag) {
  Energy2 sHat, scale, kgdotq, kpdotq;
  Energy rsHat, rsubg;
  // storage of pt and rapidity of the q/qbar radiation
  Energy pTg, pTp, pTV, mTV;
  double yj,x1,x2,wgtg, wgtp, rho, rhomx, estimatefact, pdffact, jacobfact, partdipole, sumdipole;
  // radiation phase space variables:
  double phi, xiabqr, viqr, zq, uq, zi, xi, Nfactorg, Nfactorp, cg, cp, labbeta;
  Lorentz5Momentum k_a_g, k_b_g, k_ga_g, k_V_g, k_quk_g, k_a_p, k_b_p, k_ga_p, k_V_p, k_quk_p, kbarp_ga, kbarp_V, k_glu, kpgt, testk, kvgt;
  // the charge of radiated q/qbar
  int idq, fi; 
  if (qgflag==1) {
    fi = 1;
    _quark = partons_[0];
    _quarkpi = partons_[1]->CC();
    Nfactorg = qgFactor_g;
    Nfactorp = qgFactor_p;}
  else {
    fi = 0;
    _quark = partons_[1];
    _quarkpi = partons_[0]->CC();
    Nfactorg = gqbarFactor_g;
    Nfactorp = gqbarFactor_p;}
  idq = _quarkpi->id();
  chargeqr= _quarkpi->iCharge()/3.*idq/abs(idq);
  labbeta = (x_[0]-x_[1])/(x_[0]+x_[1]);
  if(!quarkplus_) labbeta *= -1.;
  // limits on the rapidity of the jet
  double minyj = -8.0,maxyj = 8.0, myjpr, exymin, exymax;
  // PDFs for the Born cross section
  double pdf[4];
  pdf[0]=beams_[0]->pdf()->xfx(beams_[0],partons_[0],_muF2,x_[0]);
  pdf[1]=beams_[1]->pdf()->xfx(beams_[1],partons_[1],_muF2,x_[1]);
  
  // in the singular region that quark radiation collinear with initial state gluon (Dgq):
  rhomx = (1.0 -x_[fi])/(2.0*sqrt(x_[fi]));
  rho = rhomx;
  pTg = rho*systemMass_;
  rsubg = systemMass_;
  // prefactor for veto algorithm
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //cg = alphaS_->overestimateValue()/Constants::twopi*
  //Nfactorg*(maxyj-minyj)/(power_-1.);
  cg = alphaS_->overestimateValue()/Constants::twopi*Nfactorg*(maxyj-minyj);
  do {
    UseRandom::rnd();
    // generate pT
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //pTg = rsubg*pow(pow(double(pTg/rsubg),1.-power_)-log(UseRandom::rnd())/cg,1./(1.-power_));
    pTg *= pow(UseRandom::rnd(),1./cg);
    // generate rapidity of the jet
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    minyj = -8.0;
    maxyj = 8.0;
    rho = pTg/systemMass_;
    /*
    exymin= (rhomx/rho - sqrt(sqr(rhomx/rho)-1.0));
    exymax= (rhomx/rho + sqrt(sqr(rhomx/rho)-1.0));
    if (fi==0 && exymin<rho*sqrt(x_[fi])) exymin = rho*sqrt(x_[fi]);
    if (fi==1 && exymax>1.0/(rho*sqrt(x_[fi]))) exymax = 1.0/(rho*sqrt(x_[fi]));
    myjpr = log(exymin) + (double(fi)-0.5)*log(x_[1-fi]);
    if (myjpr>minyj) minyj = myjpr;
    myjpr = log(exymax) + (double(fi)-0.5)*log(x_[1-fi]);
    if (myjpr<maxyj) maxyj = myjpr;
    */     
    yj = UseRandom::rnd()*(maxyj-minyj)+ minyj;
    // generate phi
    phi = UseRandom::rnd()*2.0*Pi;
    // transform to CS variables:
    rho = pTg/systemMass_;
  if (qgflag==2) {
    viqr = rho *exp(-yj) *sqrt(x_[0]/x_[1]);}
  else {
    viqr = rho *exp( yj) *sqrt(x_[1]/x_[0]);}
    xiabqr = viqr*(1.0-viqr)/(sqr(rho)+viqr);

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&                                                                                                                                                           
    if (xiabqr>1.|| xiabqr<x_[fi]) continue;
    if (viqr<0. || viqr>(1.-xiabqr)) continue;

  if (qgflag==2) {
    x1 = x_[0]/xiabqr;
    x2 = x_[1];}
  else {
    x1 = x_[0];
    x2 = x_[1]/xiabqr;
  }
  // sHat
  sHat  = sqr(systemMass_)/xiabqr; //sHat  = x1*x2*s_;
  rsHat = sqrt(sHat);

    if(!quarkplus_) yj *= -1.;
    // real radiation phase space:
    k_quk_g=     Lorentz5Momentum (pTg *cos(phi ),pTg *sin(phi ),
				   pTg *sinh(       yj     ),
				   pTg *cosh(      yj      ),ZERO);
    // Suppose we have construct the pT and C-S variables, and then we can generate the n+1 events according to the R/B ratios 
    // of the two singular regions using the highest-kT bid technique: 
    scale = sqr(systemMass_)+ sqr(pTg);
    if (qgflag==2){
      k_a_g = _p_partona/xiabqr;
      k_b_g = _p_partonb;
      k_glu = k_a_g;
      pdffact = PDFratio(x1, x_[0], scale, 4, fi, beams_[0], beams_[1]);}
    else {
      k_a_g = _p_partona;
      k_b_g = _p_partonb/xiabqr;
      k_glu = k_b_g;
      pdffact = PDFratio(x2, x_[1], scale, 4, fi, beams_[0], beams_[1]);}  
  
    k_ga_g = InvLortr(_p_partona, _p_partonb, xiabqr, k_quk_g, _p_photon, fi);
    k_V_g =  InvLortr(_p_partona, _p_partonb, xiabqr, k_quk_g, _p_boson, fi);
    zi = 1.0/(1.0 + k_quk_g.dot(k_V_g)/k_ga_g.dot(k_quk_g+k_V_g));
    kgdotq = k_glu.dot(k_quk_g);
    kpdotq = k_ga_g.dot(k_quk_g);

    //Born phase space for Dpq:
    kbarp_ga = k_ga_g/zi;
    kbarp_V  = k_V_g + k_quk_g -(1.0-zi)*kbarp_ga;
    //kbarp_V.rescaleMass();
    kbarp_V.setMass(sqrt(MV2));
    kbarp_V.rescaleEnergy();
    //generator()->log() <<"fi="<<fi<<" x2="<<x2<< " x1="<< x1<<" viqr="<<viqr<<" xiabqr="<<xiabqr<<" pTg="<<pTg/GeV<<" k_quk_g.e="<<k_quk_g.e()/GeV
    //		       <<" conserve?="<<k_quk_g.isNear(k_a_g+k_b_g-k_ga_g-k_V_g, 0.000001)<<"\n"; 
    partdipole = Dipolegluqr(_p_partona, _p_partonb, _p_photon, _p_boson, charge0, charge1, MV2, kgdotq, alphaS_->value(scale), xiabqr);
    sumdipole = partdipole +Dipolepqr(k_a_g, k_b_g, kbarp_ga, kbarp_V, chargeqr, MV2, kpdotq, zi, alphaS_->value(scale), fi);

    //    generator()->log() <<" xiabqr="<<xiabqr<<" correct?="<<xiabqr-(1.0-k_quk_g.dot(k_a_g+k_b_g)/k_a_g.dot(k_b_g))<<"  viqr="<<viqr<<" correct?="
    //                   <<viqr-kgdotq/k_a_g.dot(k_b_g)<<" zi="<<zi<<" kgdotq="<<kgdotq/sqr(GeV)<<" kpdotq= "<<kpdotq/sqr(GeV)<<
    //                   " underlyingBorn="<<MatrQV(k_a_g, k_b_g, kbarp_ga, kbarp_V, chargeqr, MV2, alphaS_->value(scale), fi)<<"\n";


    //jacobfact = sHat/(16.0*sqr(Pi)*xiabqr)/(1.0-viqr)*sqr(xiabqr)/sqr(systemMass_);
    //jacobfact = sHat/(16.0*sqr(Pi))/(1.0-viqr)*sqr(xiabqr)/systemMass_;
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    jacobfact = xiabqr/(16.0*sqr(Pi))/(1.0-viqr);
    //estimatefact = alphaS_->overestimateValue()/(4.0*Pi)*
    //  Nfactorg/sqr(pTg/GeV)*pow(rsubg/pTg,(power_-1.0))*(2.*8./(maxyj-minyj));
    estimatefact = _alphas/alphaS_->ratio(scale)/(4.0*Pi)*Nfactorg/sqr(pTg/GeV);

    // now for the weight
    // pdf bit
    if(pdffact<=0.) {
      wgtg=0.;
      continue;
    }

    // final bit
    wgtg = QGQratio(k_a_g, k_b_g, k_ga_g, k_V_g, k_quk_g, scale, qgflag)*partdipole/sumdipole*pdffact/estimatefact*jacobfact;
    //generator()->log() <<"pTg ="<<pTg/GeV<<" partdipole ="<<partdipole<<" sumdipole ="<<sumdipole<<" estimatefact ="<<estimatefact<<" jacobfact ="<<jacobfact<<" pdffact ="<<pdffact<<"  Weight= "<<wgtg<<"\n";
    //wgtg = QGQratio(k_a_g, k_b_g, k_ga_g, k_V_g, k_quk_g, scale, qgflag)*pdffact/estimatefact*jacobfact;
    if(wgtg>1.) generator()->log() << "Weight greater than one for quark emission from "<<qgflag <<"gluon, Weight = " << wgtg << "\n";
    //wgtg = 1.1;
    if (UseRandom::rnd()<wgtg) break;
  }
  //  while(UseRandom::rnd()>wgtg && pTg>pTmin_);
  while(pTg>pTmin_);

  if(pTg<pTmin_) pTg=-GeV;


  //pTg=-GeV;     


  //   in the singular region that quark radiation collinear with final state photon (Dpq):
  Energy2 smm2, mm1, kt2, s2, pkypzu, mm3, mm4, mm5, mm6, tepT;
  Energy4 testi;
  s2 =  sqr(systemMass_);
  smm2 = s2 - MV2;
  Energy smm, pT1, ephoton, pTmax, pTte;
  smm = sqrt(smm2);
  double chy, shy, mm2, emyjpq, ymaxpq, uqlim, tech, tesh, mm2te;
  //generator()->log()<<"Here begin to generate quark radiation collinear with final state photon (Dpq):\n";
  int noy; 
  double Iuz;
  pTmax = (systemMass_-sqrt(MV2))/2.0;
  pT1 = smm2/(2.0*systemMass_);
  pTp = pTmax;
  // prefactor for veto algorithm
  minyj = -8.0;
  maxyj = 8.0;
  cp = alphaS_->overestimateValue()/(2.*Pi)*Nfactorp*(maxyj-minyj)/(power_photon-1.);
  do {
    //  do {
    // generate pT
    pTp = smm*pow(pow(double(pTp/smm),1.-power_photon)-log(UseRandom::rnd())/cp,1./(1.-power_photon));
    pTte = pTp;
    // generate rapidity of the jet
    noy = 0;
    Iuz = 1.0;
    //minyj = -8.0;
    //maxyj = 8.0;
    emyjpq = pT1/pTp - sqrt(MV2)/systemMass_;
    //if (emyjpq>1.0) {
      ymaxpq = log(sqrt(sqr(emyjpq)-1.0)+emyjpq )- 1.0e-7;
      //  if (ymaxpq<maxyj){
        maxyj = ymaxpq;
        minyj = -ymaxpq;
    yj = UseRandom::rnd()*(maxyj-minyj)+ minyj;
    // generate phi
    phi = UseRandom::rnd()*2.*Pi;
    // transform to CS variables:
    kt2 = sqr(pTp);
    chy = (exp(yj)+exp(-yj))/2.;
    shy = (exp(yj)-exp(-yj))/2.;
    mm1 = sqrt(sqr(smm2 - 2.0*pTp*systemMass_*chy)-4.0*kt2*MV2);
    zq = systemMass_* ((smm2 -2.0*pTp*systemMass_*chy +2.0*kt2)*(systemMass_ -pTp*chy)
		       - Iuz*pTp*shy*mm1)/(smm2*(s2+kt2-2.0*pTp*systemMass_*chy));
    uq = 1.0 -(1.0 -2.0*pTp*systemMass_*chy/smm2)/zq;
    testi = (smm2*(1.0+zq*uq) -2.0*uq*s2)*(smm2*(1.0-uq)*(1.0-zq*uq) -2.0*uq*MV2);
    if (testi/sqr(sqr(GeV))<0.0) {
      Iuz = -1.0;
      yj = -yj;
      shy = -shy;
      zq = systemMass_* ((smm2 -2.0*pTp*systemMass_*chy +2.0*kt2)*(systemMass_ -pTp*chy)
    			 - Iuz*pTp*shy*mm1)/(smm2*(s2+kt2-2.0*pTp*systemMass_*chy));
      uq = 1.0 -(1.0 -2.0*pTp*systemMass_*chy/smm2)/zq;
    }
    uqlim = (1.0-zq)*smm2/((1.0-zq)*smm2 + MV2);
    //}
    //else noy = 1;
    //}
    // if (noy==1 || zq<=0.0 || zq>=1.0 || uq>=(1.0-zq)*smm2/((1.0-zq)*smm2+MV2) || uq<=0.0)
    //   generator()->log()<<"!!!!!!!!z and u are outside the phyiscal range!!!!!!!!!  \n";    

    // calculate x_1 and x_2
    x1 = x_[0];
    x2 = x_[1];
   
    tepT = smm2*uq*zq*zq*(smm2*(1.0-uq)*(1.0-zq) -uq*MV2)/(smm2*sqr(1.0-zq*uq) -4.0*zq*uq*MV2);
    tech = (1.0-zq*(1.0-uq))*smm2/(2.0*sqrt(tepT*s2));    
    tesh = (2.0*tepT*s2-zq*((s2+MV2)*uq-smm2*(1.0-uq)*(1.0-zq))*smm2/2.0)/(sqrt(tepT*s2)*sqrt(zq*zq*sqr(smm2)-4.0*s2*tepT));
    //if (abs((tepT-kt2)/kt2)>0.0001)
    //generator()->log()<<"!!!!!!!!test the C-S variables:   kt2="<<kt2/sqr(GeV)<<" tepT2="<<tepT/sqr(GeV)<<" chy="<<chy<<" techy="<<tech<<" shy="<<shy<<" teshy="<<tesh<<" shy not identical?="
    //           <<(abs(shy-tesh)<0.0000001)<<"  Iuz= "<<Iuz<<" testi= "<<testi/sqr(sqr(GeV))<<"\n";
    //generator()->log()<<"@@@@@@check the frame: Born conserve?=  "<<_p_photon.isNear(_p_partona+_p_partonb-_p_boson, 0.000000000001)<<" (_p_partona+_p_partonb).z= "<<(_p_partona.z()+_p_partonb.z())/GeV<<"\n";

    if(pTp/(pT1*zq)>(1.0-1.0e-6) || uq/uqlim>(1.0-1.0e-6)) {
      wgtp=0.;
      continue;
    }

    // real radiation phase space:
    k_quk_p=     Lorentz5Momentum (pTp *cos(phi ),pTp *sin(phi ),
				   pTp *sinh(       yj     ),
				   pTp *cosh(      yj      ),ZERO);
    ephoton = smm2/(2.0*systemMass_);
    kpgt =       Lorentz5Momentum (-pTp*cos(phi)/zq, -pTp*sin(phi)/zq,
				   sqrt(sqr(ephoton)-kt2/(zq*zq)), ephoton, ZERO);
    k_quk_p.rotateZ(-kpgt.phi());
    k_quk_p.rotateY(-kpgt.theta());
    k_quk_p.rotateZ(kpgt.phi());
    kpgt = _p_photon;    
    kpgt.boost(0.,0.,-labbeta);
    testk = Lorentz5Momentum (ZERO, ZERO, ephoton, ephoton, ZERO);
    k_quk_p.rotateY(kpgt.theta());
    k_quk_p.rotateZ(kpgt.phi());
    k_quk_p.boost(0.,0.,labbeta);
    testk.rotateY(kpgt.theta());
    testk.rotateZ(kpgt.phi());
    kvgt = _p_boson;
    kvgt.boost(0.,0.,-labbeta);
    //generator()->log()<<"~~~~~~ test  rotation in CM frame: photon.z= "<<kpgt.z()/GeV<<" boson.z="<<kvgt.z()/GeV<<" photon.e="<<kpgt.e()/GeV<<" smm2/(2.0*systemMass_="
    //                  <<smm2/(2.0*systemMass_)/GeV<<" px?  "<<(kpgt.x()+kvgt.x())/GeV
    //	      <<" near?  "<<testk.isNear(kpgt, 0.0001)<<" boson.e= "<<(s2+MV2)/(2.0*systemMass_)/GeV<<" or? = "<<kvgt.e()/GeV<<"\n";
    testk.boost(0.,0.,labbeta);

    // Suppose we have construct the pT and C-S variables, and then we can generate the n+1 events according to the R/B ratios 
    // of the two singular regions using the highest-kT bid technique: 
    k_a_p = _p_partona;  
    k_b_p = _p_partonb;
    k_ga_p = _p_photon*zq;
    k_V_p =  _p_boson -k_quk_p +(1.0-zq)*_p_photon;
    //k_V_p.rescaleMass();
    k_V_p.setMass(sqrt(MV2));
    k_V_p.rescaleEnergy();
    kpdotq = k_ga_p.dot(k_quk_p);
    //generator()->log()<<"##########check the rotation: uq=k_ga.k_quk/(k_quk+k_V).k_ga? "<<uq-kpdotq/k_ga_p.dot(k_quk_p+k_V_p)<<" zq=z_qVga? "<<zq-1.0/(1.0+k_quk_p.dot(k_V_p)/k_ga_p.dot(k_quk_p+k_V_p))<<
    //  "         the momentum of photon correct? ="<<_p_photon.isNear(testk, 0.00000001)<<"\n";
    //generator()->log() <<" x2="<<x2<< " x1="<< x1<<" uq="<<uq<<" zq="<<zq<<" pTp="<<pTp/GeV<<" k_quk_p.e="<<k_quk_p.e()/GeV
    //                 <<" Is photon soft? photon.e= "<<k_ga_p.e()/GeV<<" conserve?="<<k_quk_p.isNear(k_a_p+k_b_p-k_ga_p-k_V_p, 0.000000001)<<"\n";

    xi = 1.0- k_quk_p.dot(_p_partona+ _p_partonb)/_p_partona.dot(_p_partonb);
    if (qgflag==2){      
      kgdotq = k_a_g.dot(k_quk_p);}
    else {      
      kgdotq = k_b_g.dot(k_quk_p);}

    kbarp_ga = Lortr(k_a_p, k_b_p, xi, k_quk_p, k_ga_p, fi);
    kbarp_V  = Lortr(k_a_p, k_b_p, xi, k_quk_p, k_V_p, fi);

    scale = sqr(systemMass_)+ sqr(pTp);
    partdipole = Dipolepqr(_p_partona, _p_partonb, _p_photon, _p_boson, chargeqr, MV2, kpdotq, zq, alphaS_->value(scale), fi);
    if (qgflag==2){
      sumdipole = partdipole +Dipolegluqr(xi*k_a_p, k_b_p, kbarp_ga, kbarp_V, charge0, charge1, MV2, kgdotq, alphaS_->value(scale), xi);}
    else {
      sumdipole = partdipole +Dipolegluqr(k_a_p, xi*k_b_p, kbarp_ga, kbarp_V, charge0, charge1, MV2, kgdotq, alphaS_->value(scale), xi);}

    //mm2 = 4.0*s2*uq*(smm2*(1.0-uq)*(1.0-zq) -uq*MV2)/(smm2*(smm2*sqr(1.0-zq*uq) -4.0*uq*zq*MV2));
    // modify:
    //mm2 = 4.0*s2*uq*(smm2*(1.0-uq)*(1.0-zq) -uq*MV2)/(smm2*(s2*sqr(1.0-zq*uq) - MV2*sqr(1.0+zq*uq)));
    //if ((1.0-mm2)<0.0 || abs(4.0*s2*kt2/(sqr(zq)*sqr(smm2)) - mm2)>0.01) { 
    //  generator()->log()<<" Jacobian is nan: zq="<<zq<<" uq="<<uq<<" ulim="<<(1.0-zq)*smm2/((1.0-zq)*smm2 + MV2)<<" Iuz="<<Iuz<<" z^2...="<<1.0-4.0*s2*kt2/(sqr(zq)*sqr(smm2))
    //			<<" kT/kT1="<<2.0*pTp*systemMass_/smm2/zq<<" pT="<<pTp/GeV<<" yj="<<yj<<" chy="<<chy<<" shy="<<shy<<" 1.0-mm2="<<1.0-mm2<<"\n";
      //mm2 = 0.0;
    //}
    //if (abs(4.0*s2*kt2/(sqr(zq)*sqr(smm2)) - mm2)>0.01)
    mm2 = 4.0*s2*kt2/(sqr(zq)*sqr(smm2));
    //" chy="<<chy<<" techy="<<tech<<" shy="<<shy<<" teshy="<<tesh
    mm4 = smm2*(1.0+zq*uq)-2.0*uq*s2;
    mm5 = smm2*sqr(1.0+zq*uq)-4.0*uq*zq*s2;
    mm6 = smm2*(1.0-uq)*(1.0-zq*uq)-2.0*uq*MV2;
    mm3 = mm4*mm5/mm6;
    jacobfact = 2.0*_p_photon.dot(_p_boson)*zq/(16.0*sqr(Pi)) *abs(1.0/(sqr(smm2)*zq*zq*sqrt(1.0 -mm2))*mm3);
    mm2te = 4.0*s2*uq*(smm2*(1.0-uq)*(1.0-zq) -uq*MV2)/(smm2*(smm2*sqr(1.0-zq*uq) -4.0*uq*zq*MV2));

    /*
    if (qgflag==2)
      generator()->log() <<" uq="<<uq<<" zq="<<zq<<" xi="<<xi<<" kpdotq="<<kpdotq/sqr(GeV)<<" kgdotq= "<<kgdotq/sqr(GeV)<<" underlyingBorn="<<
	                 MatrBorn(xi*k_a_p, k_b_p, kbarp_ga, kbarp_V, charge0, charge1,MV2)<<" conserve?="
                         <<kbarp_ga.isNear(xi*k_a_p+k_b_p-kbarp_V, 0.0000000001)<<" mm2="<<mm2<<" mm3="<<mm3/sqr(GeV)<<" Iuz="<<Iuz<<"\n";
    else
      generator()->log() <<" uq="<<uq<<" zq="<<zq<<" xi="<<xi<<" kpdotq="<<kpdotq/sqr(GeV)<<" kgdotq= "<<kgdotq/sqr(GeV)<<" underlyingBorn="<<
	                 MatrBorn(k_a_p, xi*k_b_p, kbarp_ga, kbarp_V, charge0, charge1,MV2)<<" conserve?="
                         <<kbarp_ga.isNear(k_a_p+xi*k_b_p-kbarp_V, 0.0000000001)<<" mm2="<<mm2<<" mm3="<<mm3/sqr(GeV)<<" Iuz="<<Iuz<<"\n";
    */

    estimatefact = alphaS_->overestimateValue()/(4.0*Pi)*
      Nfactorp/sqr(pTp/GeV)*pow(smm/pTp,(power_photon-1.0))*(2.*8./(maxyj-minyj));
    pdffact = PDFratio(x1, x2, scale, 3, fi, beams_[0], beams_[1]);
    
    if(pdffact<=0.) {
      wgtp=0.;
      continue;
    }

    // final bit
    pTp =abs(k_quk_p.perp((k_ga_p+k_quk_p).vect())/GeV)*GeV;
    //if (k_quk_p.isNear(k_a_p+k_b_p-k_ga_p-k_V_p, 0.000000000001)==1) {
    wgtp = QGQratio(k_a_p, k_b_p, k_ga_p, k_V_p, k_quk_p, scale, qgflag)*partdipole/sumdipole*pdffact/estimatefact*jacobfact;

    if(wgtp>1.) generator()->log() << "Weight greater than one for quark emission from photon, Weight = " << wgtp << "\n";
    //generator()->log() <<"pTp ="<<pTp/GeV<<" or?="<<abs(k_quk_p.perp((k_ga_p+k_quk_p).vect())/GeV)<<" partdipole ="<<partdipole<<" sumdipole ="<<sumdipole<<" estimatefact ="<<estimatefact<<" jacobfact ="<<jacobfact<<" pdffact ="<<pdffact<<"  Weight= "<<wgtp<<"\n";
    //}
    //pTp =abs(k_quk_p.perp((k_ga_p+k_quk_p).vect())/GeV)*GeV;
    //wgtp = wgtp*10000.0;
    //if(wgtp>1.) generator()->log() << "Weight greater than one for quark emission from photon, Weight = " << wgtp << "\n";
    //wgtp = 1.1;
    if (UseRandom::rnd()<wgtp) break;
  }
  //while(UseRandom::rnd()>wgtp && pTp>pTmin_);
  while(pTp>pTmin_);
  

  if(pTp<pTmin_) pTp=-GeV;

  //generator()->log() << " quark radiation with Highest-pT-bid method: wgtg="<<wgtg<<" pTg="<<pTg/GeV<<" wgtp="<<wgtp<<" pTp="<<pTp/GeV<<
  //  " sbar="<<systemMass_/GeV<<"\n";  
  //pTp=-GeV;
  
  // select the real radiation event with Highest-pT-bid method
  if(pTg > pTp){
    if (qgflag==1) {
      QGQISR_qg = 1;
      pTqg_ = pTg;
      pGammaqg_ = k_ga_g;
      pVqg_     = k_V_g;
      pQoutqg_  = k_quk_g;
      pQinqg_   = k_a_g;
      pGqg_     = k_b_g;
      if (std::isnan(wgtg)) {
	generator()->log() <<"NaNWeight wgtg_a : pV="<<pVqg_/GeV<<" pT="<<pTqg_/GeV<<" real="<<QGQratio(k_a_g,k_b_g,k_ga_g,k_V_g,k_quk_g,sqr(systemMass_)+sqr(pTg),qgflag)
			   <<" jacobfact="<<sHat/(16.0*sqr(Pi)*xiabqr)/(1.0-viqr)*sqr(xiabqr)/sqr(systemMass_)<<" estimatefact="<<alphaS_->overestimateValue()/(4.0*Pi)*
	  Nfactorg/sqr(pTg/GeV)*pow(rsubg/pTg,(power_-1.0))<<"\n";
      }
      generator()->log() <<" QGQ weight at gluon collinear singular region:  "<<wgtg<<"\n";
    }
    else if (qgflag==2) {
      QGQISR_gqbar = 1;
      pTgqbar_     = pTg;
      pGammagqbar_ = k_ga_g;
      pVgqbar_     = k_V_g;
      pQoutgqbar_  = k_quk_g;
      pGgqbar_     = k_a_g;
      pQingqbar_   = k_b_g;
      if (std::isnan(wgtg)) {
	generator()->log() <<"NaNWeight wgtg_b : pV="<<pVgqbar_/GeV<<" pT="<<pTgqbar_/GeV<<" real="<<QGQratio(k_a_g,k_b_g,k_ga_g,k_V_g,k_quk_g,sqr(systemMass_)
	       +sqr(pTg),qgflag)<<" jacobfact="<<sHat/(16.0*sqr(Pi)*xiabqr)/(1.0-viqr)*sqr(xiabqr)/sqr(systemMass_)<<" estimatefact="<<alphaS_->overestimateValue()/(4.0*Pi)*
	  Nfactorg/sqr(pTg/GeV)*pow(rsubg/pTg,(power_-1.0))<<"\n";
      }
      generator()->log() <<" GQarQar weight at gluon collinear singular region:  "<<wgtg<<"\n";
    }
  }
  else{
    if (qgflag==1) {
      QGQISR_qg = 0;
      pTqg_ = pTp;
      pGammaqg_ = k_ga_p;
      pVqg_     = k_V_p;
      pQoutqg_  = k_quk_p;
      pQinqg_   = k_a_p;
      pGqg_     = k_b_p;
      if (std::isnan(wgtp)) {
	generator()->log() <<"NaNWeight wgtp_a : pV="<<k_V_p/GeV<<" pT="<<pTp/GeV<<" pTte="<<pTte/GeV<<" real="<<QGQratio(k_a_p,k_b_p,k_ga_p,k_V_p,k_quk_p,sqr(systemMass_)+sqr(pTp),qgflag)
			   <<" jacobfact="<<2.0*_p_photon.dot(_p_boson)*zq/(16.0*sqr(Pi)) *abs(1.0/(sqr(smm2)*zq*zq*sqrt(1.0 -mm2))*mm3)<<" 1-mm2="<<1.0-mm2<<" 1-mm2te="<<1.0-mm2te<<" estimatefact="
			   <<alphaS_->overestimateValue()/(4.0*Pi)*Nfactorp/sqr(pTp/GeV)*pow(smm/pTp,(power_photon-1.0))<<" smm2="<<smm2/sqr(GeV)<<"\n"<<" zq="<<zq<<" uq="<<uq<<" ulim="<<(1.0-zq)*smm2/((1.0-zq)*smm2 + MV2)<<" Iuz="<<Iuz<<" z^2...="<<1.0-4.0*s2*kt2/(sqr(zq)*sqr(smm2))<<" kT/kT1="
                           <<2.0*pTte*systemMass_/smm2/zq<<" yj="<<yj<<" chy="<<chy<<" shy="<<shy<<" beta="<<labbeta<<" pV="<<kvgt.perp()/GeV<<"\n";
      }
      generator()->log() <<" QGQ weight at photon collinear singular region:  "<<wgtp<<"\n";
    }
    else if (qgflag==2) {
      QGQISR_gqbar = 0;
      pTgqbar_     = pTp;
      pGammagqbar_ = k_ga_p;
      pVgqbar_     = k_V_p;
      pQoutgqbar_  = k_quk_p;
      pGgqbar_     = k_a_p;
      pQingqbar_   = k_b_p;
      if (std::isnan(wgtp)) {
	generator()->log() <<"NaNWeight wgtp_b : pV="<<k_V_p/GeV<<" pT="<<pTp/GeV<<" pTte="<<pTte/GeV<<" real="<<QGQratio(k_a_p,k_b_p,k_ga_p,k_V_p,k_quk_p,sqr(systemMass_)+sqr(pTp),qgflag)
			   <<" jacobfact="<<2.0*_p_photon.dot(_p_boson)*zq/(16.0*sqr(Pi)) *abs(1.0/(sqr(smm2)*zq*zq*sqrt(1.0 -mm2))*mm3)<<" 1-mm2="<<1.0-mm2<<" 1-mm2te="<<1.0-mm2te<<" estimatefact="
			   <<alphaS_->overestimateValue()/(4.0*Pi)*Nfactorp/sqr(pTp/GeV)*pow(smm/pTp,(power_photon-1.0))<<" smm2="<<smm2/sqr(GeV)<<"\n"<<" zq="<<zq<<" uq="<<uq<<" ulim="<<(1.0-zq)*smm2/((1.0-zq)*smm2 + MV2)<<" Iuz="<<Iuz<<" z^2...="<<1.0-4.0*s2*kt2/(sqr(zq)*sqr(smm2))<<" kT/kT1="
                           <<2.0*pTte*systemMass_/smm2/zq<<" yj="<<yj<<" chy="<<chy<<" shy="<<shy<<" beta="<<labbeta<<" pV="<<kvgt.perp()/GeV<<"\n";
      }
      generator()->log() <<" GQarQar weight at photon collinear singular region:  "<<wgtp<<"\n";
    }
  }

 
  /*
  if(!quarkplus_) {
    LorentzRotation trans;
    trans.rotateX(Pi);
    if (qgflag==1) {
    pGammaqg_.transform(trans);
    pVqg_    .transform(trans);
    pQoutqg_ .transform(trans);
    pQinqg_  .transform(trans);
    pGqg_    .transform(trans);
    }
    else if (qgflag==2) {
    pGammagqbar_.transform(trans);
    pVgqbar_    .transform(trans);
    pQoutgqbar_ .transform(trans);
    pQingqbar_  .transform(trans);
    pGgqbar_    .transform(trans);
    }
  }
  */
    if (qgflag==1) {
  generator()->log() << " q g -> q V Gamma: testing new \n"
		     << "Photon    = " << pGammaqg_/GeV << "\n"
		     << "Boson     = " << pVqg_/GeV << "\n"
		     << "Gluon     = " << pGqg_/GeV << "\n"
		     << "IncomgQuark/Antiquark     = " << pQinqg_/GeV << "\n"
		     << "OutgoingQuark/Antiquark = " << pQoutqg_/GeV << "\n"
		     << "sum " 
		     << (pGammaqg_+pVqg_+pQoutqg_-pGqg_-pQinqg_)/GeV << "\n";
    pTparton_ = pTqg_;
    }
    else if (qgflag==2) {
  generator()->log() << " g qbar -> qbar V Gamma: testing new \n"
		     << "Photon    = " << pGammagqbar_/GeV << "\n"
		     << "Boson     = " << pVgqbar_/GeV << "\n"
		     << "Gluon     = " << pGgqbar_/GeV << "\n"
		     << "IncomgQuark/Antiquark     = " << pQingqbar_/GeV << "\n"
		     << "OutgoingQuark/Antiquark = " << pQoutgqbar_/GeV << "\n"
		     << "sum " 
		     << (pGammagqbar_+pVgqbar_+pQoutgqbar_-pGgqbar_-pQingqbar_)/GeV << "\n";
    pTparton_ = pTgqbar_;
    }
    //ssHat_ = s2;
    //maTV_ = mTV;

    //phiparton_ = phi;
    //yparton_ = yj;
}

  // merge from VGammaHardGenerator
double MEPP2VGammaPowheg::QGQratio(Lorentz5Momentum k_a, Lorentz5Momentum  k_b, Lorentz5Momentum k_ga, Lorentz5Momentum k_V, 
                                        Lorentz5Momentum k_q, Energy2 scale, int flavorflag) {
  double mreal,mborn,ratio;
  //Quark radiation Real Matrix Element:
  mreal = MatrRealQuark(k_a, k_b, k_ga, k_V, k_q, scale, flavorflag);

  //Born Matrix Element:
  mborn = MatrBorn(_p_partona, _p_partonb, _p_photon, _p_boson, charge0, charge1, MV2);

  ratio= mreal/mborn;

  return ratio;
}
