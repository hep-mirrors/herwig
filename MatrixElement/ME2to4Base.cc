// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ME2to4Base class.
//

#include "ME2to4Base.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "Herwig++/Utilities/Kinematics.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"

using namespace Herwig;

ME2to4Base::ME2to4Base() : _prob(0.5), _onShell(false), _weightOpt(1) {}

int ME2to4Base::nDim() const {
  return 7;
}

bool ME2to4Base::generateKinematics(const double * r) {
  jacobian(1.);
  // identify the bosons
  pair<int,int> iCharge 
    = make_pair(mePartonData()[2]->iCharge()+mePartonData()[3]->iCharge(),
		mePartonData()[4]->iCharge()+mePartonData()[5]->iCharge());
  pair<tcPDPtr,tcPDPtr> bosons=
    make_pair(iCharge.first ==0 ? _z0 : _wPlus,
	      iCharge.second==0 ? _z0 : _wPlus);
  if(iCharge.first <0) bosons.first  = bosons.first ->CC();
  if(iCharge.second<0) bosons.second = bosons.second->CC();
  // generate the boson masses
  pair<Energy,Energy> mb;
  Energy ecm = sqrt(sHat());
  pair<Energy,Energy> mmin = make_pair(bosons.first ->massMin(),
				       bosons.second->massMin());
  // not kinematically possible return
  if(ecm<mmin.first+mmin.second) return false;
  // maximum masses of the outgoing particles, including kinematic limit
  pair<Energy,Energy> mmax = 
    make_pair(min(bosons.first ->massMax(),ecm-mmin.second),
	      min(bosons.second->massMax(),ecm-mmin.first ));
  // generate the mass of the first particle
  if(!bosonMass(bosons.first ,mb.first,r[0],mmin.first,mmax.first)) return false;
  // generate the mass of the second particle
  mmax.second = min(mmax.second,ecm-mb.first);
  if(!bosonMass(bosons.second,mb.second,r[1],mmin.second,mmax.second)) return false;
  // generate the momenta of the bosons
  pair<Lorentz5Momentum,Lorentz5Momentum> pBoson = 
    make_pair(Lorentz5Momentum(mb.first),Lorentz5Momentum(mb.second));
  double ctmin = -1.0, ctmax = 1.0;
  Energy q = ZERO;
  try {
    q = SimplePhaseSpace::
      getMagnitude(sHat(), pBoson.first.mass(), pBoson.second.mass());
  } catch ( ImpossibleKinematics ) {
    return false;
  }
  Energy e = sqrt(sHat())/2.0;
  Energy2 m22 = pBoson.first .mass2();
  Energy2 m32 = pBoson.second.mass2();
  Energy2 e0e2 = 2.0*e*sqrt(sqr(q) + m22);
  Energy2 e1e2 = 2.0*e*sqrt(sqr(q) + m22);
  Energy2 e0e3 = 2.0*e*sqrt(sqr(q) + m32);
  Energy2 e1e3 = 2.0*e*sqrt(sqr(q) + m32);
  Energy2 pq = 2.0*e*q;

  Energy2 thmin = lastCuts().minTij(mePartonData()[0], bosons.first);
  if ( thmin > ZERO ) ctmax = min(ctmax, (e0e2 - m22 - thmin)/pq);

  thmin = lastCuts().minTij(mePartonData()[1], bosons.first);
  if ( thmin > ZERO ) ctmin = max(ctmin, (thmin + m22 - e1e2)/pq);

  thmin = lastCuts().minTij(mePartonData()[1], bosons.second);
  if ( thmin > ZERO ) ctmax = min(ctmax, (e1e3 - m32 - thmin)/pq);

  thmin = lastCuts().minTij(mePartonData()[0], bosons.second);
  if ( thmin > ZERO ) ctmin = max(ctmin, (thmin + m32 - e0e3)/pq);

  Energy ptmin = max(lastCuts().minKT(bosons.first),
   		     lastCuts().minKT(bosons.second));
  if ( ptmin > ZERO ) {
    double ctm = 1.0 - sqr(ptmin/q);
    if ( ctm <= 0.0 ) return false;
    ctmin = max(ctmin, -sqrt(ctm));
    ctmax = min(ctmax, sqrt(ctm));
  }

  double ymin2 = lastCuts().minYStar(bosons.first);
  double ymax2 = lastCuts().maxYStar(bosons.first);
  double ymin3 = lastCuts().minYStar(bosons.second);
  double ymax3 = lastCuts().maxYStar(bosons.second);
  double ytot = lastCuts().Y() + lastCuts().currentYHat();
  if ( ymin2 + ytot > -0.9*Constants::MaxRapidity )
    ctmin = max(ctmin, sqrt(sqr(q) +  m22)*tanh(ymin2)/q);
  if ( ymax2 + ytot < 0.9*Constants::MaxRapidity )
    ctmax = min(ctmax, sqrt(sqr(q) +  m22)*tanh(ymax2)/q);
  if ( ymin3 + ytot > -0.9*Constants::MaxRapidity )
    ctmax = min(ctmax, sqrt(sqr(q) +  m32)*tanh(-ymin3)/q);
  if ( ymax3 + ytot < 0.9*Constants::MaxRapidity )
    ctmin = max(ctmin, sqrt(sqr(q) +  m32)*tanh(-ymax3)/q);
  
  if ( ctmin >= ctmax ) return false;
    
  double cth = getCosTheta(sHat(),sqr(mb.first),sqr(mb.second),
			   ctmin, ctmax, r[2]);
  double phi = Constants::twopi*UseRandom::rnd();
  Energy pt = q*sqrt(1.0-sqr(cth));
  pBoson.first .setVect(Momentum3( pt*sin(phi),  pt*cos(phi),  q*cth));
  pBoson.second.setVect(Momentum3(-pt*sin(phi), -pt*cos(phi), -q*cth));
  pBoson.first .rescaleEnergy();
  pBoson.second.rescaleEnergy();
  jacobian((pq/sHat())*Constants::pi*jacobian());
  // decay of the 1st vector boson
  bool test=Kinematics::twoBodyDecay(pBoson.first,meMomenta()[2].mass(),
				     meMomenta()[3].mass(),
				     -1.+2*r[3],r[4]*Constants::twopi,
				     meMomenta()[2],meMomenta()[3]);
  if(!test) return false;
  jacobian(0.5*jacobian()/sqr(Constants::twopi)*
	   Kinematics::pstarTwoBodyDecay(mb.first,meMomenta()[2].mass(),
					 meMomenta()[3].mass())/mb.first);
  // decay of the 2nd vector boson
  test=Kinematics::twoBodyDecay(pBoson.second,meMomenta()[4].mass(),
				meMomenta()[5].mass(),
				-1.+2*r[5],r[6]*Constants::twopi,
				meMomenta()[4],meMomenta()[5]);
  if(!test) return false;
  jacobian(0.5*jacobian()/sqr(Constants::twopi)*
	   Kinematics::pstarTwoBodyDecay(mb.second,meMomenta()[4].mass(),
					 meMomenta()[5].mass())/mb.second);
  // check cuts
  vector<LorentzMomentum> out;
  tcPDVector tout;
  for(unsigned int ix=2;ix<6;++ix) {
    out .push_back(meMomenta()[ix]);
    tout.push_back(mePartonData()[ix]);
  }
  return lastCuts().passCuts(tout, out, mePartonData()[0], mePartonData()[1]);
}

CrossSection ME2to4Base::dSigHatDR() const {
  double me=me2();
  CrossSection output = me*jacobian()/(16.0*sqr(Constants::pi)*sHat())*sqr(hbarc);
  return output;
}

void ME2to4Base::doinit() {
  MEBase::doinit();
  // get the particle data objects we need
  _wPlus  = getParticleData(ParticleID::Wplus );
  _wMinus = getParticleData(ParticleID::Wminus); 
  _z0     = getParticleData(ParticleID::Z0    );
  _gamma  = getParticleData(ParticleID::gamma);
  // get the vertices we need
  // get a pointer to the standard model object in the run
  static const tcHwSMPtr hwsm
    = dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  if (!hwsm) throw InitException() << "hwsm pointer is null in"
				   << " ME2to4Base::doinit()"
				   << Exception::abortnow;
  // get pointers to all required Vertex objects
  _vertexFFZ = hwsm->vertexFFZ();
  _vertexFFP = hwsm->vertexFFP();
  _vertexWWW = hwsm->vertexWWW();
  _vertexFFW = hwsm->vertexFFW();
}

void ME2to4Base::persistentOutput(PersistentOStream & os) const {
  os << _onShell << _prob << _weightOpt
     << _wPlus << _wMinus << _z0 << _gamma
     << _vertexFFP << _vertexFFW << _vertexFFZ << _vertexWWW;
}

void ME2to4Base::persistentInput(PersistentIStream & is, int) {
  is >> _onShell >> _prob >> _weightOpt
     >> _wPlus >> _wMinus >> _z0 >> _gamma
     >> _vertexFFP >> _vertexFFW >> _vertexFFZ >> _vertexWWW;
}

AbstractClassDescription<ME2to4Base> ME2to4Base::initME2to4Base;
// Definition of the static class description member.

void ME2to4Base::Init() {

  static ClassDocumentation<ME2to4Base> documentation
    ("There is no documentation for the ME2to4Base class");


  static Switch<ME2to4Base,bool> interfaceOnShell
    ("OnShell",
     "Whether or not the intermediate gauge bosons are on-shell",
     &ME2to4Base::_onShell, false, false, false);
  static SwitchOption interfaceOnShellYes
    (interfaceOnShell,
     "Yes",
     "Bosons on-shell",
     true);
  static SwitchOption interfaceOnShellNo
    (interfaceOnShell,
     "No",
     "Bosons can be off-shell",
     false);

}

bool ME2to4Base::bosonMass(tcPDPtr boson,Energy & mb, double r,
			   Energy minMass, Energy maxMass) {
  if(maxMass<minMass) return false;
  // boson mass2 and mass*width
  Energy2 M2 = sqr(boson->mass()), GM = boson->mass()*boson->width();
  // integration limits
  pair<Energy2,Energy2> limitsA = make_pair(sqr(minMass),sqr(maxMass));
  pair<double,double>   limitsB = make_pair(atan((limitsA.first -M2)/GM),
					    atan((limitsA.second-M2)/GM));
  // generate mass
  // on-shell case
  Energy2 mb2;
  if(_onShell) {
    mb2 = M2;
    jacobian(GM*Constants::pi/sHat()*jacobian());
  }
  // off-shell case
  else {
    // according to BW
    if(boson->id()!=ParticleID::Z0||r>_prob) {
      if(boson->id()==ParticleID::Z0) r = (r-_prob)/(1.-_prob);
      double rho = limitsB.second+r*(limitsB.first-limitsB.second);
      mb2 = GM*tan(rho)+M2;
    }
    // according to 1/m^2
    else {
      r   = r/_prob;
      mb2 = limitsA.first*limitsA.second/
	(limitsA.first+r*(limitsA.second-limitsA.first));
    }
    // calculate the jacobian
    InvEnergy2 den = GM/(limitsB.second-limitsB.first)
      /(sqr(mb2-M2)+sqr(GM));
    if(boson->id()==ParticleID::Z0) {
      den *=(1.-_prob);
      den += _prob*limitsA.first*limitsB.second/(limitsA.second-limitsA.first)/mb2;
    }
    jacobian(jacobian()/den/sHat());
  }
  // calculate the mass
  mb = sqrt(mb2);
  // return true if in allowed phase-space
  return mb2>limitsA.first&&mb2<limitsA.second;
}

double ME2to4Base::getCosTheta(Energy2 sHat, Energy2 m12, Energy2 m22,
			       double ctmin, double ctmax, double r) {
  Energy2 D1 = sHat-m12-m22;
  Energy4 lambda = sqr(D1) - 4*m12*m22;
  if (lambda < ZERO)  throw ImpossibleKinematics();
  double D =  D1 / sqrt(lambda);
  if(_weightOpt==1) {
    double fraction = (D-ctmax)/(D-ctmin);
    double costh = D - (D - ctmin) * pow(fraction, r);
    jacobian(jacobian() * (costh - D) * log(fraction));
    return costh;
  }
  else if(_weightOpt==2) {
    double prob = 0.5;
    double costh;
    double fraction1 = (D-ctmax)/(D-ctmin);
    double fraction2 = (D+ctmin)/(D+ctmax);
    if(r<=prob) {
      r /=prob;
      costh = D - (D - ctmin) * pow(fraction1, r);
    }
    else {
      r = (r-prob)/(1.-prob);
      costh =-D + (D + ctmax) * pow(fraction2, r);
    }
    jacobian(jacobian()/(prob     /((costh - D) * log(fraction1))-
			 (1.-prob)/((costh + D) * log(fraction2))));
    return costh;
  }
  else {
    jacobian(jacobian()*(ctmax-ctmin));
    return ctmin+r*(ctmax-ctmin);
  }
}
