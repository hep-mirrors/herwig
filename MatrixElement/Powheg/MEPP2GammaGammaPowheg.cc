// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2GammaGammaPowheg class.
//

#include "MEPP2GammaGammaPowheg.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "Herwig/Utilities/Maths.h"
#include "Herwig/Shower/RealEmissionProcess.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Herwig;

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEPP2GammaGammaPowheg,Herwig::HwMEBase>
describeMEPP2GammaGammaPowheg("Herwig::MEPP2GammaGammaPowheg",
			      "HwMEHadron.so HwPowhegMEHadron.so");

unsigned int MEPP2GammaGammaPowheg::orderInAlphaS() const {
  return 0;
}

unsigned int MEPP2GammaGammaPowheg::orderInAlphaEW() const {
  return 2;
}

IBPtr MEPP2GammaGammaPowheg::clone() const {
  return new_ptr(*this);
}

IBPtr MEPP2GammaGammaPowheg::fullclone() const {
  return new_ptr(*this);
}

MEPP2GammaGammaPowheg::MEPP2GammaGammaPowheg() 
  : contrib_(1), power_(0.1), process_(0), threeBodyProcess_(0),
    maxflavour_(5), alphaS_(0.), fixedAlphaS_(false),
    supressionFunction_(0), supressionScale_(0), lambda_(20.*GeV),
    preQCDqqbarq_(5.), preQCDqqbarqbar_(0.5), preQCDqg_(50.), preQCDgqbar_(50.),
    preQEDqqbarq_(40.), preQEDqqbarqbar_(0.5), preQEDqgq_(1.), preQEDgqbarqbar_(1.),
    minpT_(2.*GeV), scaleChoice_(0), scalePreFactor_(1.)
{}

void MEPP2GammaGammaPowheg::getDiagrams() const {
  tcPDPtr gamma  = getParticleData(ParticleID::gamma);
  tcPDPtr g      = getParticleData(ParticleID::g);
  for(int ix=1;ix<=maxflavour_;++ix) {
    tcPDPtr qk = getParticleData(ix);
    tcPDPtr qb = qk->CC();
    // gamma gamma
    if(process_==0 || process_ == 1) {
      add(new_ptr((Tree2toNDiagram(3), qk, qk, qb, 1, gamma, 2, gamma, -1)));
      add(new_ptr((Tree2toNDiagram(3), qk, qk, qb, 2, gamma, 1, gamma, -2)));
    }
    // gamma +jet
    if(process_==0 || process_ == 2) {
      add(new_ptr((Tree2toNDiagram(3), qk, qb, qb, 1,    gamma,
		   2, g, -4)));
      add(new_ptr((Tree2toNDiagram(3), qk, qk, qb, 2,    gamma,
		   1, g, -5)));
      add(new_ptr((Tree2toNDiagram(3), qk, qk, g,    1,    gamma,
		   2, qk, -6)));
      add(new_ptr((Tree2toNDiagram(2), qk, g, 1, qk, 3,    gamma,
		   3, qk, -7)));
      add(new_ptr((Tree2toNDiagram(3), g, qb, qb,     2,    gamma,
		   1, qb, -8)));
      add(new_ptr((Tree2toNDiagram(2), g, qb,  1, qb, 3,    gamma,
		   3, qb, -9)));
    }
    // gamma + jet + gamma
    if((process_==0 && contrib_==1) || process_ == 3) {
      // gamma + g + gamma
      if(threeBodyProcess_==0 || threeBodyProcess_==1) {
	add(new_ptr((Tree2toNDiagram(4), qk, qk, qk, qb, 1,    gamma,
		     2, gamma, 3, g, -10)));
	add(new_ptr((Tree2toNDiagram(4), qk, qk, qk, qb, 3,    gamma,
		     2, gamma, 1, g, -12)));
      }
      // Z + q + gamma
      if(threeBodyProcess_==0 || threeBodyProcess_==2) {
	add(new_ptr((Tree2toNDiagram(4),qk,qk,qk,g,1,gamma,2,gamma,3,qk, -20)));
	add(new_ptr((Tree2toNDiagram(4),qk,qk,qk,g,2,gamma,1,gamma,3,qk, -21)));
	add(new_ptr((Tree2toNDiagram(3),qk,qk,g,1,gamma,2,qk,5,gamma,5,qk,-22)));
      }
      // Z + qbar + gamma
      if(threeBodyProcess_==0 || threeBodyProcess_==3) {
	add(new_ptr((Tree2toNDiagram(4),g,qb,qb,qb,3,gamma,2,gamma,1,qb     ,-30)));
	add(new_ptr((Tree2toNDiagram(4),g,qb,qb,qb,2,gamma,3,gamma,1,qb     ,-31)));
	add(new_ptr((Tree2toNDiagram(3),g,qb,qb   ,2,gamma,1,qb,5,gamma,5,qb,-32)));
      }
    }
  }
}

Energy2 MEPP2GammaGammaPowheg::scale() const {
  Energy2 scale;
  if(scaleChoice_==0) {
    Energy pt;
    if(meMomenta()[2].perp(meMomenta()[0].vect())>=
       meMomenta()[3].perp(meMomenta()[0].vect())){
      pt = meMomenta()[2].perp(meMomenta()[0].vect());
    } else {
      pt = meMomenta()[3].perp(meMomenta()[0].vect());
    }
    scale = sqr(pt);
  }
  else if(scaleChoice_==1) {
    scale = sHat();
  }
  return scalePreFactor_*scale;
}

int MEPP2GammaGammaPowheg::nDim() const {
  return HwMEBase::nDim() + ( contrib_>=1 ? 3 : 0 );
}

bool MEPP2GammaGammaPowheg::generateKinematics(const double * r) {
  // radiative variables
  if(contrib_>=1) {
    zTilde_ = r[nDim()-1];
    vTilde_ = r[nDim()-2];
    phi_    = Constants::twopi*r[nDim()-3];
  }
  // set the jacobian
  jacobian(1.0);
  // set up the momenta
  for ( int i = 2, N = meMomenta().size(); i < N; ++i )
    meMomenta()[i] = Lorentz5Momentum(ZERO);
  // generate sHat
  Energy2 shat(sHat());
  if(mePartonData().size()==5) {
    double eps = sqr(meMomenta()[2].mass())/shat;
    jacobian(jacobian()*(1.-eps));
    shat *= eps+zTilde_*(1.-eps);
  }
  // momenta of the core process
  double ctmin = -1.0, ctmax = 1.0;
  Energy q = ZERO;
  try {
    q = SimplePhaseSpace::
      getMagnitude(shat, meMomenta()[2].mass(), ZERO);
  } 
  catch ( ImpossibleKinematics & e ) {
    return false;
  }
  Energy e = 0.5*sqrt(shat);
  Energy2 m22 = meMomenta()[2].mass2();
  Energy2 e0e2 = 2.0*e*sqrt(sqr(q) + m22);
  Energy2 e1e2 = 2.0*e*sqrt(sqr(q) + m22);
  Energy2 e0e3 = 2.0*e*sqrt(sqr(q));
  Energy2 e1e3 = 2.0*e*sqrt(sqr(q));
  Energy2 pq = 2.0*e*q;
  if(mePartonData().size()==4) {
    Energy2 thmin = lastCuts().minTij(mePartonData()[0], mePartonData()[2]);
    if ( thmin > ZERO ) ctmax = min(ctmax, (e0e2 - m22 - thmin)/pq);
    thmin = lastCuts().minTij(mePartonData()[1], mePartonData()[2]);
    if ( thmin > ZERO ) ctmin = max(ctmin, (thmin + m22 - e1e2)/pq);
    thmin = lastCuts().minTij(mePartonData()[1], mePartonData()[3]);
    if ( thmin > ZERO ) ctmax = min(ctmax, (e1e3 - thmin)/pq);
    thmin = lastCuts().minTij(mePartonData()[0], mePartonData()[3]);
    if ( thmin > ZERO ) ctmin = max(ctmin, (thmin - e0e3)/pq);
    Energy ptmin = max(lastCuts().minKT(mePartonData()[2]),
		       lastCuts().minKT(mePartonData()[3]));
    if ( ptmin > ZERO ) {
      double ctm = 1.0 - sqr(ptmin/q);
      if ( ctm <= 0.0 ) return false;
      ctmin = max(ctmin, -sqrt(ctm));
      ctmax = min(ctmax, sqrt(ctm));
    }
    double ymin2 = lastCuts().minYStar(mePartonData()[2]);
    double ymax2 = lastCuts().maxYStar(mePartonData()[2]);
    double ymin3 = lastCuts().minYStar(mePartonData()[3]);
    double ymax3 = lastCuts().maxYStar(mePartonData()[3]);
    double ytot = lastCuts().Y() + lastCuts().currentYHat();
    if ( ymin2 + ytot > -0.9*Constants::MaxRapidity )
      ctmin = max(ctmin, sqrt(sqr(q) +  m22)*tanh(ymin2)/q);
    if ( ymax2 + ytot < 0.9*Constants::MaxRapidity )
      ctmax = min(ctmax, sqrt(sqr(q) +  m22)*tanh(ymax2)/q);
    if ( ymin3 + ytot > -0.9*Constants::MaxRapidity )
      ctmax = min(ctmax,                     tanh(-ymin3));
    if ( ymax3 + ytot < 0.9*Constants::MaxRapidity )
      ctmin = max(ctmin,                     tanh(-ymax3));
    if ( ctmin >= ctmax ) return false;
  }
  double cth = getCosTheta(ctmin, ctmax, r[0]);
  Energy pt = q*sqrt(1.0-sqr(cth));
  phi(rnd(2.0*Constants::pi));
  meMomenta()[2].setVect(Momentum3( pt*sin(phi()),  pt*cos(phi()),  q*cth));
  meMomenta()[3].setVect(Momentum3(-pt*sin(phi()), -pt*cos(phi()), -q*cth));
  meMomenta()[2].rescaleEnergy();
  meMomenta()[3].rescaleEnergy();
  // jacobian
  tHat(pq*cth + m22 - e0e2);
  uHat(m22 - shat - tHat());
  jacobian(pq/shat*Constants::pi*jacobian());
  // end for 2->2 processes
  if(mePartonData().size()==4) {
    vector<LorentzMomentum> out(2);
    out[0] = meMomenta()[2];
    out[1] = meMomenta()[3];
    tcPDVector tout(2);
    tout[0] = mePartonData()[2];
    tout[1] = mePartonData()[3];
    if ( !lastCuts().passCuts(tout, out, mePartonData()[0], mePartonData()[1]) )
      return false;
    return true;
  }
  // special for 2-3 processes
  pair<double,double> x = make_pair(lastX1(),lastX2());
  // partons
  pair<tcPDPtr,tcPDPtr> partons = make_pair(mePartonData()[0],mePartonData()[1]);
  // If necessary swap the particle data objects so that 
  // first beam gives the incoming quark
  if(lastPartons().first ->dataPtr()!=partons.first) {
    swap(x.first,x.second);
  }
  // use vTilde to select the dipole for emission
  // gamma gamma g processes
  if(mePartonData()[4]->id()==ParticleID::g) {
    if(vTilde_<=0.5) {
      dipole_ = IIQCD1;
      vTilde_ = 4.*vTilde_;
    }
    else {
      dipole_ = IIQCD2;
      vTilde_ = 4.*(vTilde_-0.25);
    }
    jacobian(2.*jacobian());
  }
  // gamma gamma q processes
  else if(mePartonData()[4]->id()>0&&mePartonData()[4]->id()<6) {
    if(vTilde_<=1./3.) {
      dipole_ = IIQCD2;
      vTilde_ = 3.*vTilde_;
    }
    else if(vTilde_<=2./3.) {
      dipole_ = IFQED1;
      vTilde_ = 3.*vTilde_-1.;
    }
    else {
      dipole_ = FIQED1;
      vTilde_ = 3.*vTilde_-2.;
    }
    jacobian(3.*jacobian());
  }
  // gamma gamma qbar processes
  else if(mePartonData()[4]->id()<0&&mePartonData()[4]->id()>-6) {
    if(vTilde_<=1./3.) {
      dipole_ = IIQCD1;
      vTilde_ = 3.*vTilde_;
    }
    else if(vTilde_<=2./3.) {
      dipole_ = IFQED2;
      vTilde_ = 3.*vTilde_-1.;
    }
    else {
      dipole_ = FIQED2;
      vTilde_ = 3.*vTilde_-2.;
    }
    jacobian(3.*jacobian());
  }
  else {
    assert(false);
  }
  // initial-initial dipoles
  if(dipole_<=4) {
    double z    = shat/sHat();
    double vt   = vTilde_*(1.-z);
    double vJac = 1.-z;
    Energy pT   = sqrt(shat*vt*(1.-vt-z)/z);
    if(pT<MeV) return false;
    double rapidity;
    Energy rs=sqrt(lastS());
    Lorentz5Momentum pcmf;
    // emission from first beam
    if(dipole_<=2) {
      rapidity = -log(x.second*sqrt(lastS())/pT*vt);
      pcmf = Lorentz5Momentum(ZERO,ZERO,
			      0.5*rs*(x.first*z-x.second),
			      0.5*rs*(x.first*z+x.second));
    }
    // emission from second beam
    else {
      rapidity =  log(x.first *sqrt(lastS())/pT*vt);
      pcmf = Lorentz5Momentum(ZERO,ZERO,
			      0.5*rs*(x.first-x.second*z),
			      0.5*rs*(x.first+x.second*z));
    }
    pcmf.rescaleMass();
    Boost blab(pcmf.boostVector());
    // emission from the quark radiation
    vector<Lorentz5Momentum> pnew(5);
    pnew [0] = Lorentz5Momentum(ZERO,ZERO,0.5*rs*x.first,
				0.5*rs*x.first,ZERO);
    pnew [1] = Lorentz5Momentum(ZERO,ZERO,-0.5*rs*x.second,
				0.5*rs*x.second,ZERO) ;
    pnew [2] = meMomenta()[2];
    pnew [3] = meMomenta()[3];
    pnew [4] = Lorentz5Momentum(pT*cos(phi_),pT*sin(phi_),
				pT*sinh(rapidity),
				pT*cosh(rapidity), ZERO);
    pnew[4].rescaleEnergy();
    Lorentz5Momentum K  = pnew [0]+pnew [1]-pnew [4]; 
    Lorentz5Momentum Kt = pcmf;
    Lorentz5Momentum Ksum = K+Kt;
    Energy2 K2 = K.m2();
    Energy2 Ksum2 = Ksum.m2();
    for(unsigned int ix=2;ix<4;++ix) {
      pnew [ix].boost(blab);
      pnew [ix] = pnew [ix] - 2.*Ksum*(Ksum*pnew [ix])/Ksum2
	+2*K*(Kt*pnew [ix])/K2;
      pnew[ix].rescaleEnergy();
    }
    pcmf = Lorentz5Momentum(ZERO,ZERO,
			    0.5*rs*(x.first-x.second),
			    0.5*rs*(x.first+x.second));
    pcmf.rescaleMass();
    blab = pcmf.boostVector();
    for(unsigned int ix=0;ix<pnew.size();++ix)
      pnew[ix].boost(-blab);
    // phase-space prefactors
    jacobian(jacobian()*vJac);
    if(dipole_%2!=0) swap(pnew[3],pnew[4]);
    for(unsigned int ix=2;ix<meMomenta().size();++ix)
      meMomenta()[ix] = pnew[ix];
  }
  else if(dipole_<=8) {
    double x  = shat/sHat();
    double z = vTilde_;
    double x1 = -1./x;
    double x3 = 1.-z/x;
    double x2 = 2.+x1-x3;
    double xT = sqrt(4.*(1-x)*(1-z)*z/x);
    // rotate the momenta into the Breit-frame
    Lorentz5Momentum pin,pcmf;
    if(dipole_<=6) {
      pin = x*meMomenta()[0];
      pcmf = pin+meMomenta()[1];
    }
    else {
      pin = x*meMomenta()[1];
      pcmf = pin+meMomenta()[0];
    }
    Boost bv = pcmf.boostVector();
    meMomenta()[2].boost(bv);
    meMomenta()[3].boost(bv);
    Lorentz5Momentum q = meMomenta()[3]-pin;
    Axis axis(q.vect().unit());
    LorentzRotation rot;
    double sinth(sqrt(sqr(axis.x())+sqr(axis.y())));
    rot = LorentzRotation();
    if(axis.perp2()>1e-20) {
      rot.setRotate(-acos(axis.z()),Axis(-axis.y()/sinth,axis.x()/sinth,0.));
      rot.rotateX(Constants::pi);
    }
    if(abs(1.-q.e()/q.vect().mag())>1e-6) 
      rot.boostZ(q.e()/q.vect().mag());
    pin *= rot;
    if(pin.perp2()/GeV2>1e-20) {
      Boost trans = -1./pin.e()*pin.vect();
      trans.setZ(0.);
      rot.boost(trans);
    }
    rot.invert();
    Energy Q = sqrt(-q.m2());
    meMomenta()[4] = rot*Lorentz5Momentum( 0.5*Q*xT*cos(phi_), 0.5*Q*xT*sin(phi_),
					  -0.5*Q*x2,0.5*Q*sqrt(sqr(x2)+sqr(xT)));
    meMomenta()[3] = rot*Lorentz5Momentum(-0.5*Q*xT*cos(phi_),-0.5*Q*xT*sin(phi_),
					  -0.5*Q*x3,0.5*Q*sqrt(sqr(x3)+sqr(xT)));

    double ratio;
    if(dipole_<=6) {
      ratio =  2.*((meMomenta()[3]+meMomenta()[4])*meMomenta()[0])/sHat();
    }
    else {
      ratio =  2.*((meMomenta()[3]+meMomenta()[4])*meMomenta()[1])/sHat();
    }
    jacobian(jacobian()*ratio);
  }
  else {
    assert(false);
  }
  vector<LorentzMomentum> out(3);
  tcPDVector tout(3);
  for(unsigned int ix=0;ix<3;++ix) {
    out[ix]  = meMomenta()   [2+ix];
    tout[ix] = mePartonData()[2+ix];
  }
  return lastCuts().passCuts(tout, out, mePartonData()[0], mePartonData()[1]);
}

double MEPP2GammaGammaPowheg::me2() const {
  // Born configurations
  if(mePartonData().size()==4) {
    // gamma gamma core process
    if(mePartonData()[3]->id()==ParticleID::gamma) {
      return 2.*Constants::twopi*alphaEM_*
	loGammaGammaME(mePartonData(),meMomenta(),true);
    }
    // V jet core process
    else if(mePartonData()[3]->id()==ParticleID::g) {
      return 2.*Constants::twopi*alphaS_*
	loGammagME(mePartonData(),meMomenta(),true);
    }
    else if(mePartonData()[3]->id()>0) {
      return 2.*Constants::twopi*alphaS_*
	loGammaqME(mePartonData(),meMomenta(),true);
    }
    else if(mePartonData()[3]->id()<0) {
      return 2.*Constants::twopi*alphaS_*
	loGammaqbarME(mePartonData(),meMomenta(),true);
    }
    else {
      assert(false);
      return 0.;
    }
  }
  // hard emission configurations
  else {
    if(mePartonData()[4]->id()==ParticleID::g)
      return sHat()*realGammaGammagME   (mePartonData(),meMomenta(),dipole_,Hard,true);
    else if(mePartonData()[4]->id()>0&&mePartonData()[4]->id()<6)
      return sHat()*realGammaGammaqME   (mePartonData(),meMomenta(),dipole_,Hard,true);
    else if(mePartonData()[4]->id()<0&&mePartonData()[4]->id()>-6)
      return sHat()*realGammaGammaqbarME(mePartonData(),meMomenta(),dipole_,Hard,true);
    else {
      assert(false);
      return 0.;
    }
  }
}

CrossSection MEPP2GammaGammaPowheg::dSigHatDR() const {
  // couplings
  if(!fixedAlphaS_) alphaS_ = SM().alphaS(scale());
  alphaEM_ = SM().alphaEM();
    // cross section
  CrossSection preFactor = 
    jacobian()/(16.0*sqr(Constants::pi)*sHat())*sqr(hbarc);
  loME_ = me2();
 
  if( contrib_== 0 || mePartonData().size()==5 || 
      (mePartonData().size()==4&& mePartonData()[3]->coloured())) 
    return loME_*preFactor;
  else
    return NLOWeight()*preFactor;
}


Selector<MEBase::DiagramIndex>
MEPP2GammaGammaPowheg::diagrams(const DiagramVector & diags) const {
  if(mePartonData().size()==4) {
    if(mePartonData()[3]->id()==ParticleID::gamma) {
 
      Selector<DiagramIndex> sel;
      for ( DiagramIndex i = 0; i < diags.size(); ++i ){
	sel.insert(meInfo()[abs(diags[i]->id())], i);
      }
      return sel;
    }
    else {
      Selector<DiagramIndex> sel;
      for ( DiagramIndex i = 0; i < diags.size(); ++i ){ 
	sel.insert(meInfo()[abs(diags[i]->id())%2], i);
      }
      return sel;
    }
  }
  else {
    Selector<DiagramIndex> sel;
    for ( DiagramIndex i = 0; i < diags.size(); ++i ) { 
      if(abs(diags[i]->id()) == 10 && dipole_ == IIQCD2 ) 
	sel.insert(1., i);
      else if(abs(diags[i]->id()) == 12 && dipole_ == IIQCD1 ) 
	sel.insert(1., i);
      else if(abs(diags[i]->id()) == 20 && dipole_ == IIQCD2 ) 
      sel.insert(1., i);
      else if(abs(diags[i]->id()) == 21 && dipole_ == IFQED1 ) 
	sel.insert(1., i);
      else if(abs(diags[i]->id()) == 22 && dipole_ == FIQED1 ) 
	sel.insert(1., i);
      else 
	sel.insert(0., i);
    }
    return sel;
  }
}

Selector<const ColourLines *>
MEPP2GammaGammaPowheg::colourGeometries(tcDiagPtr diag) const {
  // colour lines for V gamma
  static ColourLines cs("1 -2");
  static ColourLines ct("1 2 -3");
  // colour lines for q qbar -> V g
  static const ColourLines cqqbar[2]={ColourLines("1 -2 5,-3 -5"),
				      ColourLines("1 5,-5 2 -3")};
  // colour lines for q g -> V q
  static const ColourLines cqg   [2]={ColourLines("1 2 -3,3 5"),
				      ColourLines("1 -2,2 3 5")};
  // colour lines for g qbar -> V qbar
  static const ColourLines cgqbar[2]={ColourLines("-3 -2 1,-1 -5"),
				      ColourLines("-2 1,-1 -3 -5")};
  // colour lines for q qbar -> V gamma g
  static const ColourLines cqqbarg[4]={ColourLines("1 2 3 7,-4 -7"),
				       ColourLines("1 2 7,-4 3 -7"),
				       ColourLines("1 7,-4 3 2 -7"),
				       ColourLines("1 2 7,-4 3 -7")};
  // colour lines for q g -> V gamma q
  static const ColourLines cqgq  [3]={ColourLines("1 2 3 -4,4 7"),
				      ColourLines("1 2 3 -4,4 7"),
				      ColourLines("1 2 -3,3 5 7")};
  // colour lines for gbar -> V gamma qbar
  static const ColourLines cqbargqbar[3]={ColourLines("1 -2 -3 -4,-1 -7"),
					  ColourLines("1 -2 -3 -4,-1 -7"),
					  ColourLines("1 -2 -3,-1 -5 -7")};
  Selector<const ColourLines *> sel;
  switch(abs(diag->id())) {
  case 1 :case 2 :
    sel.insert(1.0, &ct); 
    break;
  case 3 : 
    sel.insert(1.0, &cs);
    break;
  case 4 : 
    sel.insert(1.0, &cqqbar[0]);
    break;
  case 5:
    sel.insert(1.0, &cqqbar[1]);
    break;
  case 6: 
    sel.insert(1.0, &cqg[0]);
    break;
  case 7:
    sel.insert(1.0, &cqg[1]);
    break;
  case 8: 
    sel.insert(1.0, &cgqbar[0]);
    break;
  case 9:
    sel.insert(1.0, &cgqbar[1]);
    break;
  case 10: case 11: case 12: case 13:
    sel.insert(1.0, &cqqbarg[abs(diag->id())-10]);
    break;
  case 20: case 21: case 22:
    sel.insert(1.0, &cqgq[abs(diag->id())-20]);
    break;
  case 30: case 31: case 32:
    sel.insert(1.0, &cqbargqbar[abs(diag->id())-30]);
    break;
  default:
    assert(false);
  }
  return sel;
}

void MEPP2GammaGammaPowheg::persistentOutput(PersistentOStream & os) const {
  os << FFPvertex_ << FFGvertex_ 
     << contrib_ << power_ << gluon_ << prefactor_
     << process_ << threeBodyProcess_<< maxflavour_ 
     << alphaS_ << fixedAlphaS_
     << supressionFunction_ << supressionScale_ << ounit(lambda_,GeV)
     << alphaQCD_ << alphaQED_ << ounit(minpT_,GeV)
     << preQCDqqbarq_ << preQCDqqbarqbar_ << preQCDqg_ << preQCDgqbar_
     << preQEDqqbarq_ << preQEDqqbarqbar_ << preQEDqgq_ << preQEDgqbarqbar_
     << scaleChoice_ << scalePreFactor_;
}

void MEPP2GammaGammaPowheg::persistentInput(PersistentIStream & is, int) {
  is >> FFPvertex_ >> FFGvertex_ 
     >> contrib_ >> power_ >> gluon_ >> prefactor_
     >> process_ >> threeBodyProcess_ >> maxflavour_ 
     >> alphaS_ >> fixedAlphaS_
     >> supressionFunction_ >> supressionScale_ >> iunit(lambda_,GeV)
     >> alphaQCD_ >> alphaQED_ >> iunit(minpT_,GeV)
     >> preQCDqqbarq_ >> preQCDqqbarqbar_ >> preQCDqg_ >> preQCDgqbar_
     >> preQEDqqbarq_ >> preQEDqqbarqbar_ >> preQEDqgq_ >> preQEDgqbarqbar_
     >> scaleChoice_ >> scalePreFactor_;
}

void MEPP2GammaGammaPowheg::Init() {

  static ClassDocumentation<MEPP2GammaGammaPowheg> documentation
    ("TheMEPP2GammaGammaPowheg class implements gamma gamma production at NLO");

  static Switch<MEPP2GammaGammaPowheg,unsigned int> interfaceProcess
    ("Process",
     "Which processes to include",
     &MEPP2GammaGammaPowheg::process_, 0, false, false);
  static SwitchOption interfaceProcessAll
    (interfaceProcess,
     "All",
     "Include all the processes",
     0);
  static SwitchOption interfaceProcessGammaGamma
    (interfaceProcess,
     "GammaGamma",
     "Only include gamma gamma",
     1);
  static SwitchOption interfaceProcessVJet
    (interfaceProcess,
     "VJet",
     "Only include gamma + jet",
     2);
  static SwitchOption interfaceProcessHard
    (interfaceProcess,
     "Hard",
     "Only include hard radiation contributions",
     3);

  static Switch<MEPP2GammaGammaPowheg,unsigned int> interfaceThreeBodyProcess
    ("ThreeBodyProcess",
     "The possible three body processes to include",
     &MEPP2GammaGammaPowheg::threeBodyProcess_, 0, false, false);
  static SwitchOption interfaceThreeBodyProcessAll
    (interfaceThreeBodyProcess,
     "All",
     "Include all processes",
     0);
  static SwitchOption interfaceThreeBodyProcessqqbar
    (interfaceThreeBodyProcess,
     "qqbar",
     "Only include q qbar -> gamma gamma g processes",
     1);
  static SwitchOption interfaceThreeBodyProcessqg
    (interfaceThreeBodyProcess,
     "qg",
     "Only include q g -> gamma gamma q processes",
     2);
  static SwitchOption interfaceThreeBodyProcessgqbar
    (interfaceThreeBodyProcess,
     "gqbar",
     "Only include g qbar -> gamma gamma qbar processes",
     3);

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

  static Parameter<MEPP2GammaGammaPowheg,int> interfaceMaximumFlavour
    ("MaximumFlavour",
     "The maximum flavour allowed for the incoming quarks",
     &MEPP2GammaGammaPowheg::maxflavour_, 5, 1, 5,
     false, false, Interface::limited);

  static Parameter<MEPP2GammaGammaPowheg,double> interfaceAlphaS
    ("AlphaS",
     "The value of alphaS to use if using a fixed alphaS",
     &MEPP2GammaGammaPowheg::alphaS_, 0.118, 0.0, 0.2,
     false, false, Interface::limited);

  static Switch<MEPP2GammaGammaPowheg,bool> interfaceFixedAlphaS
    ("FixedAlphaS",
     "Use a fixed value of alphaS",
     &MEPP2GammaGammaPowheg::fixedAlphaS_, false, false, false);
  static SwitchOption interfaceFixedAlphaSYes
    (interfaceFixedAlphaS,
     "Yes",
     "Use a fixed alphaS",
     true);
  static SwitchOption interfaceFixedAlphaSNo
    (interfaceFixedAlphaS,
     "No",
     "Use a running alphaS",
     false);

  static Switch<MEPP2GammaGammaPowheg,unsigned int> interfaceSupressionFunction
    ("SupressionFunction",
     "Choice of the supression function",
     &MEPP2GammaGammaPowheg::supressionFunction_, 0, false, false);
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

  static Parameter<MEPP2GammaGammaPowheg,Energy> interfaceSupressionScale
    ("SupressionScale",
     "The square of the scale for the supression function",
     &MEPP2GammaGammaPowheg::lambda_, GeV, 20.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

  static Switch<MEPP2GammaGammaPowheg,unsigned int> interfaceSupressionScaleChoice
    ("SupressionScaleChoice",
     "Choice of the supression scale",
     &MEPP2GammaGammaPowheg::supressionScale_, 0, false, false);
  static SwitchOption interfaceSupressionScaleChoiceFixed
    (interfaceSupressionScaleChoice,
     "Fixed",
     "Use a fixed scale",
     0);
  static SwitchOption interfaceSupressionScaleChoiceVariable
    (interfaceSupressionScaleChoice,
     "Variable",
     "Use the pT of the hard process as the scale",
     1);

  static Reference<MEPP2GammaGammaPowheg,ShowerAlpha> interfaceShowerAlphaQCD
    ("ShowerAlphaQCD",
     "Reference to the object calculating the QCD coupling for the shower",
     &MEPP2GammaGammaPowheg::alphaQCD_, false, false, true, false, false);

  static Reference<MEPP2GammaGammaPowheg,ShowerAlpha> interfaceShowerAlphaQED
    ("ShowerAlphaQED",
     "Reference to the object calculating the QED coupling for the shower",
     &MEPP2GammaGammaPowheg::alphaQED_, false, false, true, false, false);

  static Parameter<MEPP2GammaGammaPowheg,double> interfacepreQCDqqbarq
    ("preQCDqqbarq",
     "The constant for the Sudakov overestimate for the "
     "q qbar -> V Gamma +g with emission from the q",
     &MEPP2GammaGammaPowheg::preQCDqqbarq_, 23.0, 0.0, 1000.0,
     false, false, Interface::limited);

  static Parameter<MEPP2GammaGammaPowheg,double> interfacepreQCDqqbarqbar
    ("preQCDqqbarqbar",
     "The constant for the Sudakov overestimate for the "
     "q qbar -> V Gamma +g with emission from the qbar",
     &MEPP2GammaGammaPowheg::preQCDqqbarqbar_, 23.0, 0.0, 1000.0,
     false, false, Interface::limited);

  static Switch<MEPP2GammaGammaPowheg,unsigned int> interfaceScaleChoice
    ("ScaleChoice",
     "The scale choice to use",
     &MEPP2GammaGammaPowheg::scaleChoice_, 0, false, false);
  static SwitchOption interfaceScaleChoicepT
    (interfaceScaleChoice,
     "pT",
     "Use the pT of the photons",
     0);
  static SwitchOption interfaceScaleChoiceMGammaGamma
    (interfaceScaleChoice,
     "MGammaGamma",
     "Use the mass of the photon pair",
     1);

  static Parameter<MEPP2GammaGammaPowheg,double> interfaceScalePreFactor
    ("ScalePreFactor",
     "Prefactor to change factorization/renormalisation scale",
     &MEPP2GammaGammaPowheg::scalePreFactor_, 1.0, 0.1, 10.0,
     false, false, Interface::limited);


//   prefactor_.push_back(preQCDqg_);
//   prefactor_.push_back(preQCDgqbar_);
//   prefactor_.push_back(preQEDqqbarq_);
//   prefactor_.push_back(preQEDqqbarqbar_);
//   prefactor_.push_back(preQEDqgq_);
//   prefactor_.push_back(preQEDgqbarqbar_);
}

double MEPP2GammaGammaPowheg::NLOWeight() const {
  // if leading-order return 
  if(contrib_==0) return loME_;
  // prefactors
  CFfact_ = 4./3.*alphaS_/Constants::twopi;
  TRfact_ = 1./2.*alphaS_/Constants::twopi;
  // scale
  Energy2 mu2 = scale();
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
  // real emission contribution
  double real1 = subtractedReal(x,z. first,zJac. first,oldqPDF. first,
				newqPDF. first,newgPDF. first, true);
  double real2 = subtractedReal(x,z.second,zJac.second,oldqPDF.second,
				newqPDF.second,newgPDF.second,false);
  // the total weight
  double wgt = loME_ + loME_*virt + loME_*coll + real1 + real2;
  return contrib_ == 1 ? max(0.,wgt) : max(0.,-wgt);
}

double MEPP2GammaGammaPowheg::loGammaGammaME(const cPDVector & particles,
					     const vector<Lorentz5Momentum> & momenta,
					     bool first) const {
  double output(0.);
  // analytic formula for speed
  if(!first) {
    Energy2 th = (momenta[0]-momenta[2]).m2();
    Energy2 uh = (momenta[0]-momenta[3]).m2();
    output = 4./3.*Constants::pi*SM().alphaEM(ZERO)*(th/uh+uh/th)*
      pow(double(particles[0]->iCharge())/3.,4);
  }
  // HE code result
  else {
    // wavefunctions for the incoming fermions
    SpinorWaveFunction    em_in( momenta[0],particles[0],incoming);
    SpinorBarWaveFunction ep_in( momenta[1],particles[1],incoming);
    // wavefunctions for the outgoing bosons
    VectorWaveFunction v1_out(momenta[2],particles[2],outgoing);
    VectorWaveFunction v2_out(momenta[3],particles[3],outgoing);
    vector<SpinorWaveFunction> f1;
    vector<SpinorBarWaveFunction> a1;
    vector<VectorWaveFunction> v1,v2;
    // calculate the wavefunctions
    for(unsigned int ix=0;ix<2;++ix) {
      em_in.reset(ix);
      f1.push_back(em_in);
      ep_in.reset(ix);
      a1.push_back(ep_in);
      v1_out.reset(2*ix);
      v1.push_back(v1_out);
      v2_out.reset(2*ix);
      v2.push_back(v2_out);
    }
    vector<double> me(4,0.0);
    me_.reset(ProductionMatrixElement(PDT::Spin1Half,PDT::Spin1Half,
				      PDT::Spin1,PDT::Spin1));
    vector<Complex> diag(2,0.0);
    SpinorWaveFunction inter;
    for(unsigned int ihel1=0;ihel1<2;++ihel1) {
      for(unsigned int ihel2=0;ihel2<2;++ihel2) {
	for(unsigned int ohel1=0;ohel1<2;++ohel1) {
	  for(unsigned int ohel2=0;ohel2<2;++ohel2) {
	    inter   = FFPvertex_->evaluate(ZERO,5,f1[ihel1].particle()->CC(),
					   f1[ihel1],v1[ohel1]);
	    diag[0] = FFPvertex_->evaluate(ZERO,inter,a1[ihel2],v2[ohel2]);
	    inter   = FFPvertex_->evaluate(ZERO,5,f1[ihel1].particle()->CC(),
					   f1[ihel1] ,v2[ohel2]);
	    diag[1] = FFPvertex_->evaluate(ZERO,inter,a1[ihel2],v1[ohel1]);
	    // individual diagrams
	    for (size_t ii=0; ii<2; ++ii) me[ii] += std::norm(diag[ii]);
	    // full matrix element
	    diag[0] += diag[1];
	    output += std::norm(diag[0]);
	    // storage of the matrix element for spin correlations
	    me_(ihel1,ihel2,2*ohel1,2*ohel2) = diag[0];
	  }
	}
      }
    }
    // store diagram info, etc.
    DVector save(3);
    for (size_t i = 0; i < 3; ++i) save[i] = 0.25 * me[i];
    meInfo(save);
    // spin and colour factors
    output *= 0.125/3./norm(FFPvertex_->norm());
  }
  return output;
}

double MEPP2GammaGammaPowheg::loGammaqME(const cPDVector & particles,
					 const vector<Lorentz5Momentum> & momenta,
					 bool first) const {
  double output(0.);
  // analytic formula for speed
  if(!first) {
    Energy2 sh = (momenta[0]+momenta[1]).m2();
    Energy2 th = (momenta[0]-momenta[2]).m2();
    Energy2 uh = (momenta[0]-momenta[3]).m2();
    output = -1./3./sh/th*(sh*sh+th*th+2.*uh*(sh+th+uh))*
      4.*Constants::pi*SM().alphaEM(ZERO)*
      sqr(particles[0]->iCharge()/3.);
  }
  // HE result
  else {
    vector<SpinorWaveFunction> fin;
    vector<VectorWaveFunction> gin;
    vector<SpinorBarWaveFunction> fout;
    vector<VectorWaveFunction> vout;
    SpinorWaveFunction    qin (momenta[0],particles[0],incoming);
    VectorWaveFunction    glin(momenta[1],particles[1],incoming);
    VectorWaveFunction    wout(momenta[2],particles[2],outgoing);
    SpinorBarWaveFunction qout(momenta[3],particles[3],outgoing);
    // polarization states for the particles
    for(unsigned int ix=0;ix<2;++ix) {
      qin.reset(ix) ;
      fin.push_back(qin);
      qout.reset(ix);
      fout.push_back(qout);
      glin.reset(2*ix);
      gin.push_back(glin);
      wout.reset(2*ix);
      vout.push_back(wout);
    }
    me_.reset(ProductionMatrixElement(PDT::Spin1Half,PDT::Spin1,
				      PDT::Spin1,PDT::Spin1Half));
    // compute the matrix elements
    double me[3]={0.,0.,0.};
    Complex diag[2];
    SpinorWaveFunction inters;
    SpinorBarWaveFunction interb;
    for(unsigned int ihel1=0;ihel1<2;++ihel1) {
      for(unsigned int ihel2=0;ihel2<2;++ihel2) {
	for(unsigned int ohel1=0;ohel1<2;++ohel1) {
	  // intermediates for the diagrams
	  interb= FFGvertex_->evaluate(scale(),5,particles[3]->CC(),
				       fout[ohel1],gin[ihel2]);
	  inters= FFGvertex_->evaluate(scale(),5,particles[0],
				       fin[ihel1],gin[ihel2]);
	  for(unsigned int vhel=0;vhel<2;++vhel) {
	    diag[0] = FFPvertex_->evaluate(ZERO,fin[ihel1],interb,vout[vhel]);
	    diag[1] = FFPvertex_->evaluate(ZERO,inters,fout[ohel1],vout[vhel]);
	    // diagram contributions
	    me[1] += norm(diag[0]);
	    me[2] += norm(diag[1]);
	    // total
	    diag[0] += diag[1];
	    me[0]   += norm(diag[0]);
	    me_(ihel1,2*ihel2,2*vhel,ohel1) = diag[0];
	  }
	}
      }
    }
    // results
    // initial state spin and colour average
    double colspin = 1./24./4.;
    // and C_F N_c from matrix element
    colspin *= 4.;
    DVector save;
    for(unsigned int ix=0;ix<3;++ix) {
      me[ix] *= colspin;
      if(ix>0) save.push_back(me[ix]);
    }
    meInfo(save);
    output = me[0]/norm(FFGvertex_->norm());
  }
  return output;
}

double MEPP2GammaGammaPowheg::loGammaqbarME(const cPDVector & particles,
					    const vector<Lorentz5Momentum> & momenta,
					    bool first) const {
  double output(0.);
  // analytic formula for speed
  if(!first) {
    Energy2 sh = (momenta[0]+momenta[1]).m2();
    Energy2 uh = (momenta[0]-momenta[2]).m2();
    Energy2 th = (momenta[0]-momenta[3]).m2();
    output = -1./3./sh/th*(sh*sh+th*th+2.*uh*(sh+th+uh))*
      4.*Constants::pi*SM().alphaEM()*
      sqr(particles[1]->iCharge()/3.);
  }
  // HE result
  else {
    vector<SpinorBarWaveFunction>  ain;
    vector<VectorWaveFunction> gin;
    vector<SpinorWaveFunction> aout;
    vector<VectorWaveFunction> vout;
    VectorWaveFunction    glin (momenta[0],particles[0],incoming);
    SpinorBarWaveFunction qbin (momenta[1],particles[1],incoming);
    VectorWaveFunction    wout (momenta[2],particles[2],outgoing);
    SpinorWaveFunction    qbout(momenta[3],particles[3],outgoing);
    // polarization states for the particles
    for(unsigned int ix=0;ix<2;++ix) {
      qbin .reset(ix  );
      ain .push_back(qbin );
      qbout.reset(ix  );
      aout.push_back(qbout);
      glin.reset(2*ix);
      gin.push_back(glin);
      wout.reset(2*ix);
      vout.push_back(wout);
    }
    // if calculation spin corrections construct the me
    me_.reset(ProductionMatrixElement(PDT::Spin1,PDT::Spin1Half,
				      PDT::Spin1,PDT::Spin1Half));
    // compute the matrix elements
    double me[3]={0.,0.,0.};
    Complex diag[2];
    SpinorWaveFunction inters;
    SpinorBarWaveFunction interb;
    for(unsigned int ihel1=0;ihel1<2;++ihel1) {
      for(unsigned int ihel2=0;ihel2<2;++ihel2) {
	for(unsigned int ohel1=0;ohel1<2;++ohel1) {
	  // intermediates for the diagrams
	  inters= FFGvertex_->evaluate(scale(),5,particles[3]->CC(),
				       aout[ohel1],gin[ihel1]);
	  interb= FFGvertex_->evaluate(scale(),5,particles[1],
				       ain[ihel2],gin[ihel1]);
	  for(unsigned int vhel=0;vhel<2;++vhel) {
	    diag[0]= FFPvertex_->evaluate(ZERO,inters,ain[ihel2],vout[vhel]);
	    diag[1]= FFPvertex_->evaluate(ZERO,aout[ohel1],interb,vout[vhel]);
	    // diagram contributions
	    me[1] += norm(diag[0]);
	    me[2] += norm(diag[1]);
	    // total
	    diag[0] += diag[1];
	    me[0]   += norm(diag[0]);
	    me_(2*ihel1,ihel2,2*vhel,ohel1) = diag[0];
	  }
	}
      }
    }
    // results
    // initial state spin and colour average
    double colspin = 1./24./4.;
    // and C_F N_c from matrix element
    colspin *= 4.;
    DVector save;
    for(unsigned int ix=0;ix<3;++ix) {
      me[ix] *= colspin;
      if(ix>0) save.push_back(me[ix]);
    }
    meInfo(save);
    output = me[0]/norm(FFGvertex_->norm());
  }
  return output;
}  

double MEPP2GammaGammaPowheg::loGammagME(const cPDVector & particles,
					 const vector<Lorentz5Momentum> & momenta,
					 bool first) const {
  double output(0.);
  // analytic formula for speed
  if(!first) {
    Energy2 uh = (momenta[0]-momenta[2]).m2();
    Energy2 th = (momenta[0]-momenta[3]).m2();
   output = 8./9.*double((th*th+uh*uh)/uh/th)*
     4.*Constants::pi*SM().alphaEM(ZERO)*
     sqr(particles[0]->iCharge()/3.);
  }
  else {
    vector<SpinorWaveFunction>     fin;
    vector<SpinorBarWaveFunction>  ain;
    vector<VectorWaveFunction> gout;
    vector<VectorWaveFunction> vout;
    SpinorWaveFunction    qin (momenta[0],particles[0],incoming);
    SpinorBarWaveFunction qbin(momenta[1],particles[1],incoming);
    VectorWaveFunction    wout(momenta[2],particles[2],outgoing);
    VectorWaveFunction   glout(momenta[3],particles[3],outgoing);
    // polarization states for the particles
    for(unsigned int ix=0;ix<2;++ix) {
      qin.reset(ix)  ;
      fin.push_back(qin);
      qbin.reset(ix) ;
      ain.push_back(qbin);
      glout.reset(2*ix);
      gout.push_back(glout);
      wout.reset(2*ix);
      vout.push_back(wout);
    }
    // if calculation spin corrections construct the me
    if(first) me_.reset(ProductionMatrixElement(PDT::Spin1Half,PDT::Spin1Half,
						PDT::Spin1,PDT::Spin1));
    // compute the matrix elements
    double me[3]={0.,0.,0.};
    Complex diag[2];
    SpinorWaveFunction inters;
    SpinorBarWaveFunction interb;
    for(unsigned int ihel1=0;ihel1<2;++ihel1) {
      for(unsigned int ihel2=0;ihel2<2;++ihel2) {
	for(unsigned int ohel1=0;ohel1<2;++ohel1) {
	  // intermediates for the diagrams
	  inters= FFGvertex_->evaluate(scale(),5,particles[0],
				       fin[ihel1],gout[ohel1]);
	  interb= FFGvertex_->evaluate(scale(),5,particles[1],
				       ain[ihel2],gout[ohel1]);
	  for(unsigned int vhel=0;vhel<2;++vhel) {
	    diag[0]= FFPvertex_->evaluate(ZERO,fin[ihel1],interb,vout[vhel]);
	    diag[1]= FFPvertex_->evaluate(ZERO,inters,ain[ihel2],vout[vhel]);
	    // diagram contributions
	    me[1] += norm(diag[0]);
	    me[2] += norm(diag[1]);
	    // total
	    diag[0] += diag[1];
	    me[0]   += norm(diag[0]);
	    if(first) me_(ihel1,ihel2,vhel,2*ohel1) = diag[0];
	  }
	}
      }
    }
    // results
    // initial state spin and colour average
    double colspin = 1./9./4.;
    // and C_F N_c from matrix element
    colspin *= 4.;
    DVector save;
    for(unsigned int ix=0;ix<3;++ix) {
      me[ix] *= colspin;
      if(ix>0) save.push_back(me[ix]);
    }
    meInfo(save);
    output = me[0]/norm(FFGvertex_->norm());
  }
  return output;
}

InvEnergy2 MEPP2GammaGammaPowheg::
realGammaGammagME(const cPDVector & particles,
		  const vector<Lorentz5Momentum> & momenta,
		  DipoleType dipole, RadiationType rad,
		  bool ) const {
  // matrix element
  double sum = realME(particles,momenta);
  // loop over the QCD and QCD dipoles
  InvEnergy2 dipoles[2];
  pair<double,double> supress[2];
  // compute the two dipole terms
  unsigned int iemit = 4, ihard = 3;
  double x = (momenta[0]*momenta[1]-momenta[iemit]*momenta[1]-
	      momenta[iemit]*momenta[0])/(momenta[0]*momenta[1]);
  Lorentz5Momentum Kt = momenta[0]+momenta[1]-momenta[iemit];
  vector<Lorentz5Momentum> pa(4),pb(4);
  // momenta for q -> q g/gamma emission
  pa[0] = x*momenta[0];
  pa[1] =   momenta[1];
  Lorentz5Momentum K = pa[0]+pa[1];
  Lorentz5Momentum Ksum = K+Kt;
  Energy2 K2 = K.m2();
  Energy2 Ksum2 = Ksum.m2();
  pa[2] = momenta[2]-2.*Ksum*(Ksum*momenta[2])/Ksum2+2*K*(Kt*momenta[2])/K2;
  pa[2].setMass(momenta[2].mass());
  pa[3] = momenta[ihard]
    -2.*Ksum*(Ksum*momenta[ihard])/Ksum2+2*K*(Kt*momenta[ihard])/K2;
  pa[3].setMass(ZERO);
  cPDVector part(particles.begin(),--particles.end());
  part[3] = particles[ihard];
  // first leading-order matrix element
  double lo1 = loGammaGammaME(part,pa);
  // first dipole
  dipoles[0] = 1./(momenta[0]*momenta[iemit])/x*(2./(1.-x)-(1.+x))*lo1;
  supress[0] = supressionFunction(momenta[iemit].perp(),pa[3].perp());
  // momenta for qbar -> qbar g/gamma emission
  pb[0] =   momenta[0];
  pb[1] = x*momenta[1];
  K = pb[0]+pb[1];
  Ksum = K+Kt;
  K2 = K.m2();
  Ksum2 = Ksum.m2();
  pb[2] = momenta[2]-2.*Ksum*(Ksum*momenta[2])/Ksum2+2*K*(Kt*momenta[2])/K2;
  pb[2].setMass(momenta[2].mass());
  pb[3] = momenta[ihard]
    -2.*Ksum*(Ksum*momenta[ihard])/Ksum2+2*K*(Kt*momenta[ihard])/K2;
  pb[3].setMass(ZERO);
  // second LO matrix element
  double lo2 = loGammaGammaME(part,pb);
  // second dipole
  dipoles[1] = 1./(momenta[1]*momenta[iemit])/x*(2./(1.-x)-(1.+x))*lo2;
  supress[1] = supressionFunction(momenta[iemit].perp(),pb[3].perp());
  for(unsigned int ix=0;ix<2;++ix) dipoles[ix] *= 4./3.;
  // denominator for the matrix element
  InvEnergy2 denom = abs(dipoles[0]) + abs(dipoles[1]);
  // contribution
  if( denom==ZERO || dipoles[(dipole-1)/2]==ZERO ) return ZERO;
  sum *= abs(dipoles[(dipole-1)/2])/denom;
  // final coupling factors
  InvEnergy2 output;
  if(rad==Subtraction) {
    output = alphaS_*alphaEM_*
      (sum*UnitRemoval::InvE2*supress[(dipole-1)/2].first 
       - dipoles[(dipole-1)/2]);
  }
  else {
    output = alphaS_*alphaEM_*sum*UnitRemoval::InvE2;
    if(rad==Hard)        output *=supress[(dipole-1)/2].second;
    else if(rad==Shower) output *=supress[(dipole-1)/2].first ;
  }
  return output;
}

InvEnergy2 MEPP2GammaGammaPowheg::realGammaGammaqME(const cPDVector & particles,
						const vector<Lorentz5Momentum> & momenta,
						DipoleType dipole, RadiationType rad,
						bool ) const {
  double sum = realME(particles,momenta);
  // initial-state QCD dipole
  double x = (momenta[0]*momenta[1]-momenta[4]*momenta[1]-
	      momenta[4]*momenta[0])/(momenta[0]*momenta[1]);
  Lorentz5Momentum Kt = momenta[0]+momenta[1]-momenta[4];
  vector<Lorentz5Momentum> pa(4);
  pa[0] =   momenta[0];
  pa[1] = x*momenta[1];
  Lorentz5Momentum K = pa[0]+pa[1];
  Lorentz5Momentum Ksum = K+Kt;
  Energy2 K2 = K.m2();
  Energy2 Ksum2 = Ksum.m2();
  pa[2] = momenta[2]-2.*Ksum*(Ksum*momenta[2])/Ksum2+2*K*(Kt*momenta[2])/K2;
  pa[2].setMass(momenta[2].mass());
  pa[3] = momenta[3]
    -2.*Ksum*(Ksum*momenta[3])/Ksum2+2*K*(Kt*momenta[3])/K2;
  pa[3].setMass(ZERO);
  cPDVector part(particles.begin(),--particles.end());
  part[1] = particles[4]->CC();
  double lo1 = loGammaGammaME(part,pa);
  InvEnergy2 D1 = 0.5/(momenta[1]*momenta[4])/x*(1.-2.*x*(1.-x))*lo1;
  // initial-final QED dipole
  vector<Lorentz5Momentum> pb(4);
  x = 1.-(momenta[3]*momenta[4])/(momenta[4]*momenta[0]+momenta[0]*momenta[3]);
  pb[3] = momenta[4]+momenta[3]-(1.-x)*momenta[0];
  pb[0] = x*momenta[0];
  pb[1] = momenta[1];
  pb[2] = momenta[2];
  double z = momenta[0]*momenta[3]/(momenta[0]*momenta[3]+momenta[0]*momenta[4]);
  part[1] = particles[1];
  part[3] = particles[4];
  double lo2 = loGammaqME(part,pb);
  Energy pT = sqrt(-(pb[0]-pb[3]).m2()*(1.-x)*(1.-z)*z/x);
  InvEnergy2 DF = 1./(momenta[4]*momenta[3])/x*(1./(1.-x+z)-2.+z)*lo2;
  InvEnergy2 DI = 1./(momenta[0]*momenta[3])/x*(1./(1.-x+z)-1.-x)*lo2;
  DI *= sqr(double(particles[0]->iCharge())/3.);  
  DF *= sqr(double(particles[0]->iCharge())/3.);
  InvEnergy2 denom = abs(D1)+abs(DI)+abs(DF);
  pair<double,double> supress;
  InvEnergy2 term;
  if     ( dipole == IFQED1 ) {
    term  = DI;
    supress = supressionFunction(pT,pb[3].perp());
  }
  else if( dipole == FIQED1 ) {
    term  = DF;
    supress = supressionFunction(pT,pb[3].perp());
  }
  else                        {
    term  = D1;
    supress = supressionFunction(momenta[4].perp(),pa[3].perp()); 
  }
  if( denom==ZERO || term == ZERO ) return ZERO;
  sum *= abs(term)/denom;
  // final coupling factors
  InvEnergy2 output;
  if(rad==Subtraction) {
    output = alphaS_*alphaEM_*
      (sum*UnitRemoval::InvE2*supress.first - term);
  }
  else {
    output = alphaS_*alphaEM_*sum*UnitRemoval::InvE2;
    if(rad==Hard)        output *= supress.second;
    else if(rad==Shower) output *= supress.first ;
  }
  // final coupling factors
  return output;
}

InvEnergy2 MEPP2GammaGammaPowheg::
realGammaGammaqbarME(const cPDVector & particles,
		     const vector<Lorentz5Momentum> & momenta,
		     DipoleType dipole, RadiationType rad,
		     bool) const {
  double sum = realME(particles,momenta);
  // initial-state QCD dipole
  double x = (momenta[0]*momenta[1]-momenta[4]*momenta[1]-momenta[4]*momenta[0])/
    (momenta[0]*momenta[1]);
  Lorentz5Momentum Kt = momenta[0]+momenta[1]-momenta[4];
  vector<Lorentz5Momentum> pa(4);
  pa[0] = x*momenta[0];
  pa[1] =   momenta[1];
  Lorentz5Momentum K = pa[0]+pa[1];
  Lorentz5Momentum Ksum = K+Kt;
  Energy2 K2 = K.m2();
  Energy2 Ksum2 = Ksum.m2();
  pa[2] = momenta[2]-2.*Ksum*(Ksum*momenta[2])/Ksum2+2*K*(Kt*momenta[2])/K2;
  pa[2].setMass(momenta[2].mass());
  pa[3] = momenta[3]
    -2.*Ksum*(Ksum*momenta[3])/Ksum2+2*K*(Kt*momenta[3])/K2;
  pa[3].setMass(ZERO);
  cPDVector part(particles.begin(),--particles.end());
  part[0] = particles[4]->CC();
  double lo1 = loGammaGammaME(part,pa);
  InvEnergy2 D1 = 0.5/(momenta[0]*momenta[4])/x*(1.-2.*x*(1.-x))*lo1;
  // initial-final QED dipole
  vector<Lorentz5Momentum> pb(4);
  x = 1.-(momenta[3]*momenta[4])/(momenta[4]*momenta[1]+momenta[1]*momenta[3]);
  pb[3] = momenta[4]+momenta[3]-(1.-x)*momenta[1];
  pb[0] = momenta[0];
  pb[1] = x*momenta[1];
  pb[2] = momenta[2];
  double z = momenta[1]*momenta[3]/(momenta[1]*momenta[3]+momenta[1]*momenta[4]);
  part[0] = particles[0];
  part[3] = particles[4];
  double lo2 = loGammaqbarME(part,pb);
  Energy pT = sqrt(-(pb[1]-pb[3]).m2()*(1.-x)*(1.-z)*z/x);
  InvEnergy2 DF = 1./(momenta[4]*momenta[3])/x*(2./(1.-x+z)-2.+z)*lo2;
  InvEnergy2 DI = 1./(momenta[0]*momenta[3])/x*(2./(1.-x+z)-1.-x)*lo2;
  InvEnergy2 term;
  DI *= sqr(double(particles[1]->iCharge())/3.);
  DF *= sqr(double(particles[1]->iCharge())/3.);
  InvEnergy2 denom = abs(D1)+abs(DI)+abs(DF);
  pair<double,double> supress;
  if     ( dipole == IFQED2 ) {
    term  = DI;
    supress = supressionFunction(pT,pb[3].perp());
  }
  else if( dipole == FIQED2 ) {
    term  = DF;
    supress = supressionFunction(pT,pb[3].perp());
  }
  else                        {
    term  = D1;
    supress = supressionFunction(momenta[4].perp(),pa[3].perp()); 
  }
  if( denom==ZERO || dipole==ZERO ) return ZERO;
  sum *= abs(term)/denom;
  // final coupling factors
  InvEnergy2 output;
  if(rad==Subtraction) {
    output = alphaS_*alphaEM_*
      (sum*UnitRemoval::InvE2*supress.first - term);
  }
  else {
    output = alphaS_*alphaEM_*sum*UnitRemoval::InvE2;
    if(rad==Hard)        output *= supress.second;
    else if(rad==Shower) output *= supress.first ;
  }
  // final coupling factors
  return output;
}

double MEPP2GammaGammaPowheg::
realME(const cPDVector & particles,
       const vector<Lorentz5Momentum> & momenta) const {
  vector<SpinorWaveFunction> qin;
  vector<SpinorBarWaveFunction> qbarin;
  vector<VectorWaveFunction> wout,pout,gout;
  SpinorWaveFunction    q_in;  
  SpinorBarWaveFunction qbar_in;
  VectorWaveFunction    g_out;  
  VectorWaveFunction    v_out  (momenta[2],particles[2],outgoing);
  VectorWaveFunction    p_out  (momenta[3],particles[3],outgoing);
  // q qbar -> gamma gamma g
  if(particles[4]->id()==ParticleID::g) {
    q_in    = SpinorWaveFunction    (momenta[0],particles[0],incoming);
    qbar_in = SpinorBarWaveFunction (momenta[1],particles[1],incoming);
    g_out   = VectorWaveFunction    (momenta[4],particles[4],outgoing);
  }
  // q g -> gamma gamma q
  else if(particles[4]->id()>0) {
    q_in    = SpinorWaveFunction    (momenta[0],particles[0],incoming);
    qbar_in = SpinorBarWaveFunction (momenta[4],particles[4],outgoing);
    g_out   = VectorWaveFunction    (momenta[1],particles[1],incoming);
  }
  else if(particles[4]->id()<0) {
    q_in    = SpinorWaveFunction    (momenta[4],particles[4],outgoing);
    qbar_in = SpinorBarWaveFunction (momenta[1],particles[1],incoming);
    g_out   = VectorWaveFunction    (momenta[0],particles[0],incoming);
  }
  else assert(false);
  for(unsigned int ix=0;ix<2;++ix) {
    q_in.reset(ix);
    qin.push_back(q_in);
    qbar_in.reset(ix);
    qbarin.push_back(qbar_in);
    g_out.reset(2*ix);
    gout.push_back(g_out);
    p_out.reset(2*ix);
    pout.push_back(p_out);
    v_out.reset(2*ix);
    wout.push_back(v_out);
  }
  vector<Complex> diag(6 , 0.);
  Energy2 mu2 = scale();
  double sum(0.);
  for(unsigned int ihel1=0;ihel1<2;++ihel1) {
    for(unsigned int ihel2=0;ihel2<2;++ihel2) {
      for(unsigned int whel=0;whel<2;++whel) {
	for(unsigned int phel=0;phel<2;++phel) {
	  for(unsigned int ghel=0;ghel<2;++ghel) {
	    // first diagram
	    SpinorWaveFunction inters1 = 
	      FFPvertex_->evaluate(ZERO,5,qin[ihel1].particle()->CC(),qin[ihel1],pout[phel]);
	    SpinorBarWaveFunction inters2 = 
	      FFPvertex_->evaluate(ZERO,5,qin[ihel1].particle(),
			       qbarin[ihel2],wout[whel]);
	    diag[0] = FFGvertex_->evaluate(mu2,inters1,inters2,gout[ghel]);
	    // second diagram
	    SpinorWaveFunction inters3 = 
	      FFGvertex_->evaluate(mu2,5,qin[ihel1].particle()->CC(),qin[ihel1],gout[ghel]);
	    SpinorBarWaveFunction inters4 = 
	      FFPvertex_->evaluate(ZERO,5,qbarin[ihel2].particle()->CC(),
				   qbarin[ihel2],pout[phel]);
	    diag[1] = FFPvertex_->evaluate(ZERO,inters3,inters4,wout[whel]);
	    // fourth diagram
	    diag[2] = FFPvertex_->evaluate(ZERO,inters3,inters2,pout[phel]);
	    // fifth diagram
	    SpinorBarWaveFunction inters5 = 
	      FFGvertex_->evaluate(mu2,5,qbarin[ihel2].particle()->CC(),
				   qbarin[ihel2],gout[ghel]);
	    diag[3] = 
	      FFPvertex_->evaluate(ZERO,inters1,inters5,wout[whel]);
	    // sixth diagram
	    SpinorWaveFunction inters6 = 
	      FFPvertex_->evaluate(ZERO,5,qbarin[ihel2].particle(),
				   qin[ihel1],wout[whel]);
	    diag[4] = FFGvertex_->evaluate(mu2,inters6,inters4,gout[ghel]);
	    // eighth diagram
	    diag[5] = FFPvertex_->evaluate(ZERO,inters6,inters5,pout[phel]);
	    // sum
	    Complex dsum = std::accumulate(diag.begin(),diag.end(),Complex(0.));
	    sum += norm(dsum);
	  }
	}
      }
    }
  }
  // divide out the em and strong couplings
  sum /= norm(FFGvertex_->norm()*FFPvertex_->norm());
  // final spin and colour factors spin = 1/4 colour = 4/9
  if(particles[4]->id()==ParticleID::g) sum /= 9.;
  // final spin and colour factors spin = 1/4 colour = 4/(3*8)
  else                                  sum /= 24.;
  // finally identical particle factor
  return 0.5*sum;
}

double MEPP2GammaGammaPowheg::subtractedVirtual() const {
  double v = 1+tHat()/sHat();
  double born = (1-v)/v+v/(1-v);
  double finite_term = born*
    (2./3.*sqr(Constants::pi)-3.+sqr(log(v))+sqr(log(1-v))+3.*log(1-v))+
    2.+2.*log(v)+2.*log(1-v)+3.*(1-v)/v*(log(v)-log(1-v))+
    (2.+v/(1-v))*sqr(log(v))+(2.+(1-v)/v)*sqr(log(1-v));

  double virt = ((6.-(2./3.)*sqr(Constants::pi))*
		 born-2.+finite_term);

    return virt/born;
}

double MEPP2GammaGammaPowheg::subtractedReal(pair<double,double> x, double z,
					     double zJac, double oldqPDF, double newqPDF,
					     double newgPDF,bool order) const {
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
  double phase = zJac*vJac/z;
  // real emission q qbar
  vector<double> output(4,0.);
  double realQQ(0.),realGQ(0.);
  if(!(zTilde_<1e-7 || vt<1e-7 || 1.-z-vt < 1e-7 )) {
    cPDVector particles(mePartonData());
    particles.push_back(gluon_);
    // calculate the full 2->3 matrix element
    realQQ = sHat()*phase*newqPDF/oldqPDF*
      realGammaGammagME(particles,pnew,order ? IIQCD1 : IIQCD2,Subtraction,false);
    if(order) {
      particles[0] = gluon_;
      particles[4] = mePartonData()[0]->CC();
      realGQ = sHat()*phase*newgPDF/oldqPDF*
	realGammaGammaqbarME(particles,pnew,IIQCD2,Subtraction,false);
    }
    else {
      particles[1] = gluon_;
      particles[4] = mePartonData()[1]->CC();
      realGQ = sHat()*phase*newgPDF/oldqPDF*
	realGammaGammaqME   (particles,pnew,IIQCD1,Subtraction,false);
    }
  }
  // return the answer
  return realQQ+realGQ;
}

double MEPP2GammaGammaPowheg::collinearQuark(double x, Energy2 mu2, double jac, double z,
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

double MEPP2GammaGammaPowheg::collinearGluon(Energy2 mu2, double jac, double z,
					     double oldPDF, double newPDF) const {
  if(1.-z < 1.e-8) return 0.;
  return TRfact_*jac/z*newPDF/oldPDF*
    ((sqr(z)+sqr(1.-z))*log(sqr(1.-z)*sHat()/z/mu2)
     +2.*z*(1.-z));
}

void MEPP2GammaGammaPowheg::doinit() {
  HwMEBase::doinit();
  vector<unsigned int> mopt(2,1);
  massOption(mopt);
  // get the vertices we need
  // get a pointer to the standard model object in the run
  static const tcHwSMPtr hwsm
    = dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  if (!hwsm) throw InitException() << "hwsm pointer is null in"
				   << " MEPP2GammaGamma::doinit()"
				   << Exception::abortnow;
  // get pointers to all required Vertex objects
  FFPvertex_ = hwsm->vertexFFP();
  FFGvertex_ = hwsm->vertexFFG();
  gluon_ = getParticleData(ParticleID::g);
  // sampling factors
  prefactor_.push_back(preQCDqqbarq_);
  prefactor_.push_back(preQCDqqbarqbar_);
  prefactor_.push_back(preQCDqg_);
  prefactor_.push_back(preQCDgqbar_);
  prefactor_.push_back(preQEDqqbarq_);
  prefactor_.push_back(preQEDqqbarqbar_);
  prefactor_.push_back(preQEDqgq_);
  prefactor_.push_back(preQEDgqbarqbar_);
}

RealEmissionProcessPtr MEPP2GammaGammaPowheg::
generateHardest(RealEmissionProcessPtr born,
		ShowerInteraction inter) {
  beams_.clear();
  partons_.clear();
  bool QCDAllowed = inter !=ShowerInteraction::QED;
  bool QEDAllowed = inter !=ShowerInteraction::QCD;
  // find the incoming particles
  // and get the particles to be showered
  ParticleVector incoming,particlesToShower;
  pair<double,double> x;
  for(unsigned int ix=0;ix<born->bornIncoming().size();++ix) {
    incoming.push_back( born->bornIncoming()[ix] );
    beams_.push_back(dynamic_ptr_cast<tcBeamPtr>(born->hadrons()[ix]->dataPtr()));
    partons_.push_back( born->bornIncoming()[ix]->dataPtr() );
    particlesToShower.push_back( born->bornIncoming()[ix] );
    if(ix==0) x.first  = incoming.back()->momentum().rho()/born->hadrons()[ix]->momentum().rho();
    else      x.second = incoming.back()->momentum().rho()/born->hadrons()[ix]->momentum().rho();
  }
  // find the parton which should be first
  if( ( particlesToShower[1]->id() > 0 && particlesToShower[0]->id() < 0 ) || 
      ( particlesToShower[0]->id() == ParticleID::g &&
  	particlesToShower[1]->id() < 6 && particlesToShower[1]->id() > 0 ) ) {
    swap(particlesToShower[0],particlesToShower[1]);
    swap(partons_[0],partons_[1]);
    swap(beams_  [0],beams_  [1]);
    swap(x.first    ,x.second   );
  }
  // check that quark is along +ve z direction
  quarkplus_ = particlesToShower[0]->momentum().z() > ZERO;
  // outgoing partons
  for(unsigned int ix=0;ix<born->bornOutgoing().size();++ix) {
    particlesToShower.push_back( born->bornOutgoing()[ix] );
  }
  if(particlesToShower.size()!=4) return RealEmissionProcessPtr();
  if(particlesToShower[2]->id()!=ParticleID::gamma)
    swap(particlesToShower[2],particlesToShower[3]);
  if(particlesToShower[3]->id()==ParticleID::gamma) {
    if(QCDAllowed) return hardQCDEmission(born,particlesToShower,x);
  }
  else {
    if(QEDAllowed) return hardQEDEmission(born,particlesToShower,x);
  }
  return born;
}

RealEmissionProcessPtr MEPP2GammaGammaPowheg::
hardQCDEmission(RealEmissionProcessPtr born,
		ParticleVector particlesToShower,
		pair<double,double> x) {
  Energy rootS = sqrt(lastS());
  // limits on the rapidity of the jet
  double minyj = -8.0,maxyj = 8.0;
  // generate the hard emission
  vector<Energy> pT;
  Energy pTmax(-GeV);
  cPDVector selectedParticles;
  vector<Lorentz5Momentum> selectedMomenta;
  int iemit(-1);
  for(unsigned int ix=0;ix<4;++ix) {
    pT.push_back(0.5*generator()->maximumCMEnergy());
    double a = alphaQCD_->overestimateValue()/Constants::twopi*
      prefactor_[ix]*(maxyj-minyj);
    cPDVector particles;
    for(unsigned int iy=0;iy<particlesToShower.size();++iy) 
      particles.push_back(particlesToShower[iy]->dataPtr());
    if(ix<2) particles.push_back(gluon_);
    else if(ix==2) {
      particles.push_back(particles[0]->CC());
      particles[0] = gluon_;
    }
    else {
      particles.push_back(particles[1]->CC());
      particles[1] = gluon_;
    }
    vector<Lorentz5Momentum> momenta(5);
    do {
      pT[ix] *= pow(UseRandom::rnd(),1./a);
      double y = UseRandom::rnd()*(maxyj-minyj)+ minyj;
      double vt,z;
      if(ix%2==0) {
  	vt = pT[ix]*exp(-y)/rootS/x.second;
  	z  = (1.-pT[ix]*exp(-y)/rootS/x.second)/(1.+pT[ix]*exp( y)/rootS/x.first );
  	if(z>1.||z<x.first) continue;
      }
      else {
  	vt = pT[ix]*exp( y)/rootS/x.first ;
  	z  = (1.-pT[ix]*exp( y)/rootS/x.first )/(1.+pT[ix]*exp(-y)/rootS/x.second );
  	if(z>1.||z<x.second) continue;
      }
      if(vt>1.-z || vt<0.) continue;
      if(ix%2==0) {
  	momenta[0] = particlesToShower[0]->momentum()/z;
  	momenta[1] = particlesToShower[1]->momentum();
      }
      else {
  	momenta[0] = particlesToShower[0]->momentum();
  	momenta[1] = particlesToShower[1]->momentum()/z;
      }
      double phi = Constants::twopi*UseRandom::rnd();
      momenta[2] = particlesToShower[2]->momentum();
      momenta[3] = particlesToShower[3]->momentum();
      if(!quarkplus_) y *= -1.;
      momenta[4] = Lorentz5Momentum(pT[ix]*cos(phi),pT[ix]*sin(phi),
  				    pT[ix]*sinh(y),pT[ix]*cosh(y), ZERO);
      Lorentz5Momentum K = momenta[0] + momenta[1] - momenta[4]; 
      Lorentz5Momentum Kt = momenta[2]+momenta[3];
      Lorentz5Momentum Ksum = K+Kt;
      Energy2 K2 = K.m2(), Ksum2 = Ksum.m2();
      for(unsigned int iy=2;iy<4;++iy) {
  	momenta [iy] = momenta [iy] - 2.*Ksum*(Ksum*momenta [iy])/Ksum2
  	  +2*K*(Kt*momenta [iy])/K2;
      }
      // matrix element piece
      double wgt = alphaQCD_->ratio(sqr(pT[ix]))*z/(1.-vt)/prefactor_[ix]/loME_;
      if(ix==0)
  	wgt *= sqr(pT[ix])*realGammaGammagME(particles,momenta,IIQCD1,Shower,false);
      else if(ix==1)
  	wgt *= sqr(pT[ix])*realGammaGammagME(particles,momenta,IIQCD2,Shower,false);
      else if(ix==2)
  	wgt *= sqr(pT[ix])*realGammaGammaqbarME(particles,momenta,IIQCD1,Shower,false);
      else if(ix==3)
  	wgt *= sqr(pT[ix])*realGammaGammaqME(particles,momenta,IIQCD2,Shower,false);
      wgt *= 4.*Constants::pi/alphaS_;
      // pdf piece
      double pdf[2];
      if(ix%2==0) {
  	pdf[0] = beams_[0]->pdf()->xfx(beams_[0],partons_ [0],
  				       scale(),            x.first   )  /x.first;
  	pdf[1] = beams_[0]->pdf()->xfx(beams_[0],particles[0],
  				       scale()+sqr(pT[ix]),x.first /z)*z/x.first;
      }
      else {
  	pdf[0] = beams_[1]->pdf()->xfx(beams_[1],partons_ [1],
  				       scale()            ,x.second  )  /x.second;
  	pdf[1] = beams_[1]->pdf()->xfx(beams_[1],particles[1],
  				       scale()+sqr(pT[ix]),x.second/z)*z/x.second;
      }
      if(pdf[0]<=0.||pdf[1]<=0.) continue;
      wgt *= pdf[1]/pdf[0];
      if(wgt>1.) generator()->log() << "Weight greater than one in "
  				    << "MEPP2GammaGammaPowheg::hardQCDEmission() "
  				    << "for channel " << ix 
  				    << " Weight = " << wgt << "\n";
      if(UseRandom::rnd()<wgt) break;
    }
    while(pT[ix]>minpT_);
    if(pT[ix]>minpT_ && pT[ix]>pTmax) {
      pTmax = pT[ix];
      selectedParticles = particles;
      selectedMomenta = momenta;
      iemit=ix;
    }
  }
  // if no emission
  if(pTmax<ZERO) {
    born->pT()[ShowerInteraction::QCD] = minpT_;
    return born;
  }
  // construct the HardTree object needed to perform the showers
  // create the partons
  ParticleVector newparticles;
  newparticles.push_back(selectedParticles[0]->produceParticle(selectedMomenta[0]));
  newparticles.push_back(selectedParticles[1]->produceParticle(selectedMomenta[1]));
  for(unsigned int ix=2;ix<particlesToShower.size();++ix) {
    newparticles.push_back(particlesToShower[ix]->dataPtr()->
			   produceParticle(selectedMomenta[ix]));
  }
  newparticles.push_back(selectedParticles[4]->produceParticle(selectedMomenta[4]));
  // identify the type of process
  // gluon emission
  if(newparticles.back()->id()==ParticleID::g) {
    newparticles[4]->incomingColour(newparticles[0]);
    newparticles[4]->incomingColour(newparticles[1],true);
  }
  // quark
  else if(newparticles.back()->id()>0) {
    iemit=1;
    newparticles[4]->incomingColour(newparticles[1]);
    newparticles[1]-> colourConnect(newparticles[0]);
  }
  // antiquark
  else  {
    iemit=0;
    newparticles[4]->incomingColour(newparticles[0],true);
    newparticles[1]-> colourConnect(newparticles[0]);
  }
  // add incoming
  int ispect = iemit==0 ? 1 : 0;
  if(particlesToShower[0]==born->bornIncoming()[0]) {
    born->incoming().push_back(newparticles[0]);
    born->incoming().push_back(newparticles[1]);
  }
  else {
    born->incoming().push_back(newparticles[1]);
    born->incoming().push_back(newparticles[0]);
    swap(iemit,ispect);
  }
  // add the outgoing
  for(unsigned int ix=2;ix<newparticles.size();++ix)
    born->outgoing().push_back(newparticles[ix]);
  // emitter spectator etc
  born->emitter  (iemit );
  born->spectator(ispect);
  born->emitted(4);
  // x values
  pair<double,double> xnew;
  for(unsigned int ix=0;ix<2;++ix) {
    double x = born->incoming()[ix]->momentum().rho()/born->hadrons()[ix]->momentum().rho();
    if(ix==0) xnew.first  = x;
    else      xnew.second = x;
  }
  born->x(xnew);
  // max pT
  born->pT()[ShowerInteraction::QCD] = pTmax;
  born->interaction(ShowerInteraction::QCD);
  // return the process
  return born;
}

RealEmissionProcessPtr MEPP2GammaGammaPowheg::
hardQEDEmission(RealEmissionProcessPtr born,
		ParticleVector particlesToShower,
		pair<double,double> x) {
  // return if not emission from quark
  if(particlesToShower[0]->id()!=ParticleID::g &&
     particlesToShower[1]->id()!=ParticleID::g )
    return RealEmissionProcessPtr();
  // generate the hard emission
  vector<Energy> pT;
  Energy pTmax(-GeV);
  cPDVector selectedParticles;
  vector<Lorentz5Momentum> selectedMomenta;
  int iemit(-1);
  pair<double,double> mewgt(make_pair(0.,0.));
  for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
    selectedParticles.push_back(particlesToShower[ix]->dataPtr());
    selectedMomenta.push_back(particlesToShower[ix]->momentum());
  }
  selectedParticles.push_back(getParticleData(ParticleID::gamma));
  swap(selectedParticles[3],selectedParticles[4]);
  selectedMomenta.push_back(Lorentz5Momentum());
  swap(selectedMomenta[3],selectedMomenta[4]);
  Lorentz5Momentum pin,pout;
  double xB;
  unsigned int iloc;
  if(particlesToShower[0]->dataPtr()->charged()) {
    pin = particlesToShower[0]->momentum();
    xB = x.first;
    iloc = 6;
  }
  else {
    pin = particlesToShower[1]->momentum();
    xB = x.second;
    iloc = 7;
  }
  pout = particlesToShower[3]->momentum();
  Lorentz5Momentum q = pout-pin;
  Axis axis(q.vect().unit());
  LorentzRotation rot;
  double sinth(sqrt(sqr(axis.x())+sqr(axis.y())));
  rot = LorentzRotation();
  if(axis.perp2()>1e-20) {
    rot.setRotate(-acos(axis.z()),Axis(-axis.y()/sinth,axis.x()/sinth,0.));
    rot.rotateX(Constants::pi);
  }
  if(abs(1.-q.e()/q.vect().mag())>1e-6) 
    rot.boostZ(q.e()/q.vect().mag());
  Lorentz5Momentum ptemp = rot*pin;
  if(ptemp.perp2()/GeV2>1e-20) {
    Boost trans = -1./ptemp.e()*ptemp.vect();
    trans.setZ(0.);
    rot.boost(trans);
  }
  rot.invert();
  Energy Q = sqrt(-q.m2());
  double xT = sqrt((1.-xB)/xB);
  double xTMin = 2.*minpT_/Q;
  double wgt(0.);
  double a = alphaQED_->overestimateValue()*prefactor_[iloc]/Constants::twopi;
  Lorentz5Momentum p1,p2,p3;
  do {
    wgt = 0.;
    // intergration variables dxT/xT^3
    xT *= 1./sqrt(1.-2.*log(UseRandom::rnd())/a*sqr(xT));
    // dz
    double zp = UseRandom::rnd();
    double xp = 1./(1.+0.25*sqr(xT)/zp/(1.-zp));
    if(xT<xTMin) break;
    // check allowed
    if(xp<xB||xp>1.) continue;
    // phase-space piece of the weight
    wgt = 4.*sqr(1.-xp)*(1.-zp)*zp/prefactor_[iloc]/loME_;
    // coupling
    Energy2 pT2 = 0.25*sqr(Q*xT);
    wgt *= alphaQED_->ratio(pT2);
    // matrix element
    wgt *= 4.*Constants::pi/alphaEM_;
    // PDF
    double pdf[2];
    if(iloc==6) {
      pdf[0] = beams_[0]->pdf()->
  	xfx(beams_[0],partons_[0],scale()    ,x.first    );
      pdf[1] = beams_[0]->pdf()->
  	xfx(beams_[0],partons_[0],scale()+pT2,x.first /xp);
    }
    else {
      pdf[0] = beams_[1]->pdf()->
  	xfx(beams_[1],partons_[1],scale()    ,x.second   );
      pdf[1] = beams_[1]->pdf()->
  	xfx(beams_[1],partons_[1],scale()+pT2,x.second/xp);
    }
    if(pdf[0]<=0.||pdf[1]<=0.) {
      wgt = 0.;
      continue;
    }
    wgt *= pdf[1]/pdf[0];
    // matrix element piece
    double phi = Constants::twopi*UseRandom::rnd();
    double x2 = 1.-(1.-zp)/xp;
    double x3 = 2.-1./xp-x2;
    p1=Lorentz5Momentum(ZERO,ZERO,0.5*Q/xp,0.5*Q/xp,ZERO);
    p2=Lorentz5Momentum( 0.5*Q*xT*cos(phi), 0.5*Q*xT*sin(phi),
  			 -0.5*Q*x2,0.5*Q*sqrt(sqr(xT)+sqr(x2)));
    p3=Lorentz5Momentum(-0.5*Q*xT*cos(phi),-0.5*Q*xT*sin(phi),
  			-0.5*Q*x3,0.5*Q*sqrt(sqr(xT)+sqr(x3)));
    selectedMomenta[iloc-6] = rot*p1;
    selectedMomenta[3] = rot*p3;
    selectedMomenta[4] = rot*p2;
    if(iloc==6) {
      mewgt.first  = 
  	sqr(Q)*realGammaGammaqME(selectedParticles,selectedMomenta,IFQED1,Shower,false);
      mewgt.second = 
  	sqr(Q)*realGammaGammaqME(selectedParticles,selectedMomenta,FIQED1,Shower,false);
      wgt *= mewgt.first+mewgt.second;
    }
    else {
      mewgt.first  = 
  	sqr(Q)*realGammaGammaqbarME(selectedParticles,selectedMomenta,IFQED2,Shower,false);
      mewgt.second = 
  	sqr(Q)*realGammaGammaqbarME(selectedParticles,selectedMomenta,FIQED2,Shower,false); 
      wgt *= mewgt.first+mewgt.second;
    }
    if(wgt>1.) generator()->log() << "Weight greater than one in "
  				  << "MEPP2GammaGammaPowheg::hardQEDEmission() "
  				  << "for IF channel "
  				  << " Weight = " << wgt << "\n";
  }
  while(xT>xTMin&&UseRandom::rnd()>wgt);
  // if no emission
  if(xT<xTMin) {
    born->pT()[ShowerInteraction::QED] = minpT_;
    return born;
  }
  pTmax = 0.5*xT*Q;
  iemit = mewgt.first>mewgt.second ? 2 : 3;
  // construct the object needed to perform the showers
  // create the partons
  ParticleVector newparticles;
  newparticles.push_back(selectedParticles[0]->produceParticle(selectedMomenta[0]));
  newparticles.push_back(selectedParticles[1]->produceParticle(selectedMomenta[1]));
  for(unsigned int ix=2;ix<particlesToShower.size();++ix) {
    newparticles.push_back(particlesToShower[ix]->
			   dataPtr()->produceParticle(selectedMomenta[ix==2 ? 2 : 4 ]));
  }
  newparticles.push_back(selectedParticles[3]->produceParticle(selectedMomenta[3]));
  // make the colour connections
  bool col = newparticles[3]->id()<0;
  if(particlesToShower[0]->dataPtr()->charged()) {
    newparticles[3]->incomingColour(newparticles[1],col);
    newparticles[1]->colourConnect (newparticles[0],col);
  }
  else {
    newparticles[3]->incomingColour(newparticles[0],col);
    newparticles[0]->colourConnect (newparticles[1],col);
  }
  bool FSR = iemit==3;
  // add incoming particles
  if(particlesToShower[0]==born->bornIncoming()[0]) {
    born->incoming().push_back(newparticles[0]);
    born->incoming().push_back(newparticles[1]);
  }
  else {
    born->incoming().push_back(newparticles[1]);
    born->incoming().push_back(newparticles[0]);
  }
  // IS radiatng particle
  unsigned int iemitter = born->incoming()[0]->dataPtr()->charged() ? 0 : 1;
  // add outgoing particles
  if(particlesToShower[2]==born->bornOutgoing()[0]) {
    born->outgoing().push_back(newparticles[2]);
    born->outgoing().push_back(newparticles[3]);
  }
  else {
    born->outgoing().push_back(newparticles[3]);
    born->outgoing().push_back(newparticles[2]);
  }
  born->outgoing().push_back(newparticles[4]);
  // outgoing radiating particle
  unsigned int ispectator = born->outgoing()[0]->dataPtr()->charged() ? 2 : 3;
  // get emitter and spectator right
  if(FSR) swap(iemitter,ispectator);
  born->emitter  (iemitter  );
  born->spectator(ispectator);
  born->emitted(4);
  // x values
  pair<double,double> xnew;
  for(unsigned int ix=0;ix<2;++ix) {
    double x = born->incoming()[ix]->momentum().rho()/born->hadrons()[ix]->momentum().rho();
    if(ix==0) xnew.first  = x;
    else      xnew.second = x;
  }
  born->x(xnew);
  // max pT
  born->pT()[ShowerInteraction::QED] = pTmax;
  // return the process
  born->interaction(ShowerInteraction::QED);
  return born;
}
