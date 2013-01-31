// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMZDecayerPOWHEG class.
//

#include "SMZDecayerPOWHEG.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Utilities/Maths.h"
#include "Herwig++/Shower/Base/ShowerTree.h"
#include "Herwig++/Shower/Base/ShowerProgenitor.h"
#include "Herwig++/Shower/Base/ShowerParticle.h"
#include "Herwig++/Shower/Base/Branching.h"
#include "Herwig++/Shower/Base/HardTree.h"

using namespace Herwig;

IBPtr SMZDecayerPOWHEG::clone() const {
  return new_ptr(*this);
}

IBPtr SMZDecayerPOWHEG::fullclone() const {
  return new_ptr(*this);
}

void SMZDecayerPOWHEG::persistentOutput(PersistentOStream & os) const {
  os << contrib_ << yPow_ << zPow_;
}

void SMZDecayerPOWHEG::persistentInput(PersistentIStream & is, int) {
  is >> contrib_ >> yPow_ >> zPow_;
}

ClassDescription<SMZDecayerPOWHEG> SMZDecayerPOWHEG::initSMZDecayerPOWHEG;
// Definition of the static class description member.

void SMZDecayerPOWHEG::Init() {

  static ClassDocumentation<SMZDecayerPOWHEG> documentation
    ("The SMZDecayerPOWHEG class implements the decay of the Z boson"
     " including the NLO QCD cross sections");

  static Switch<SMZDecayerPOWHEG,unsigned int> interfaceContribution
    ("Contribution",
     "Which contributions to the cross section to include",
     &SMZDecayerPOWHEG::contrib_, 1, false, false);
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

  static Parameter<SMZDecayerPOWHEG,double> interfacezPower
    ("zPower",
     "The sampling power for z",
     &SMZDecayerPOWHEG::zPow_, 0.5, 0.0, 1.0,
     false, false, Interface::limited);

  static Parameter<SMZDecayerPOWHEG,double> interfaceyPower
    ("yPower",
     "The sampling power for y",
     &SMZDecayerPOWHEG::yPow_, 0.9, 0.0, 1.0,
     false, false, Interface::limited);

}

// return the matrix element squared
double SMZDecayerPOWHEG::me2(const int, const Particle & inpart,
			     const ParticleVector & decay,
			     MEOption meopt) const {
  if(abs(decay[0]->id())!=5&&
     abs(decay[0]->id())!=4) return 1.;
  // cast the vertices
  tcFFVVertexPtr Zvertex = dynamic_ptr_cast<tcFFVVertexPtr>(FFZVertex());
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(vectors(),rho(),
					       const_ptr_cast<tPPtr>(&inpart),
					       incoming,false);
    ME(DecayMatrixElement(PDT::Spin1,PDT::Spin1Half,PDT::Spin1Half));
  }
  if(meopt==Terminate) {
    VectorWaveFunction::constructSpinInfo(vectors(),const_ptr_cast<tPPtr>(&inpart),
					  incoming,true,false);
    SpinorBarWaveFunction::
      constructSpinInfo(wavebar(),decay[0],outgoing,true);
    SpinorWaveFunction::
      constructSpinInfo(wave()   ,decay[1],outgoing,true);
    return 0.;
  }
  SpinorBarWaveFunction::
    calculateWaveFunctions(wavebar(),decay[0],outgoing);
  SpinorWaveFunction::
    calculateWaveFunctions(wave()   ,decay[1],outgoing);
  // store the NLO matrix element
  DecayMatrixElement    nloME(PDT::Spin1,PDT::Spin1Half,PDT::Spin1Half);
  LorentzPolarizationVector momDiff = 
    (decay[0]->momentum()-decay[1]->momentum())/2./
    (decay[0]->mass()+decay[1]->mass());
  // compute the matrix element
  Energy2 scale(sqr(inpart.mass()));
  for(unsigned int ifm=0;ifm<2;++ifm) {
    for(unsigned int ia=0;ia<2;++ia) {
      // extra stuff for NLO
      LorentzPolarizationVector left  = 
	wave()[ia].wave(). leftCurrent(wavebar()[ifm].wave());
      LorentzPolarizationVector right = 
	wave()[ia].wave().rightCurrent(wavebar()[ifm].wave());
      Complex scalar = 
	wave()[ia].wave().scalar(wavebar()[ifm].wave());
      for(unsigned int vhel=0;vhel<3;++vhel) {
	ME()(vhel,ifm,ia)=
	  FFZVertex()->evaluate(scale,wave()[ia],wavebar()[ifm],vectors()[vhel]);
	Complex scalarV = vectors()[vhel].wave().dot(momDiff);
	nloME(vhel,ifm,ia)= Complex(0.,1.)*Zvertex->norm()*
	  (Zvertex->right()*( left.dot(vectors()[vhel].wave())) +
	   Zvertex-> left()*(right.dot(vectors()[vhel].wave())) -
	   ( Zvertex-> left()+Zvertex->right())*scalarV*scalar);
      }
    }
  }
  // leading order piece
  double output=(ME().contract(rho())).real()*UnitRemoval::E2/scale;
  // colour connections
  if(decay[0]->hasColour())      decay[0]->antiColourNeighbour(decay[1]);
  else if(decay[1]->hasColour()) decay[1]->antiColourNeighbour(decay[0]);
  if(!decay[0]->dataPtr()->coloured()) return output;
  // NLO piece with different structure
  double nlo   =real(ME().contract(nloME ,rho()) +
		     nloME .contract(ME(),rho()))*UnitRemoval::E2/scale;
  // colour factor
  if(abs(decay[0]->id())<=6) {
    output *= 3.;
    nlo    *= 3.;
  }
  // now for the NLO bit
  double mu2 = 0.25*sqr(decay[0]->mass()+decay[1]->mass())/scale;
  double mu = sqrt(mu2);
  double mu4 = sqr(mu2);
  double lmu = log(mu);
  double v = sqrt(1.-4.*mu2),v2(sqr(v));
  double omv = 4.*mu2/(1.+v);
  double f1,f2,fNS,VNS;
  double CF=4./3.;
  double aS = SM().alphaS(scale);
  double r = omv/(1.+v),lr(log(r));
  // normal form
  if(mu>1e-4) {
    f1 = CF*aS/Constants::pi*
      ( +1. + 3.*log(0.5*(1.+v)) - 1.5*log(0.5*(1.+v2)) + sqr(Constants::pi)/6.
	- 0.5*sqr(lr) - (1.+v2)/v*(lr*log(1.+v2) + sqr(Constants::pi)/12. 
				       -0.5*log(4.*mu2)*lr + 0.25*sqr(lr)));
    fNS =  -0.5*(1.+2.*v2)*lr/v + 1.5*lr - 2./3.*sqr(Constants::pi) + 0.5*sqr(lr)
      + (1.+v2)/v*(Herwig::Math::ReLi2(r) + sqr(Constants::pi)/3. - 0.25*sqr(lr) + lr*log((2.*v/ (1.+v))));
    VNS = 1.5*log(0.5*(1.+v2)) 
      + 0.5*(1.+v2)/v*( 2.*lr*log(2.*(1.+v2)/sqr(1.+v))  + 2.*Herwig::Math::ReLi2(sqr(r)) 
			    - 2.*Herwig::Math::ReLi2(2.*v/(1.+v)) - sqr(Constants::pi)/6.)
      + log(1.-mu) - 2.*log(1.-2.*mu) - 4.*mu2/(1.+v2)*log(mu/(1.-mu)) - mu/(1.-mu)
      + 4.*(2.*mu2-mu)/(1.+v2) + 0.5*sqr(Constants::pi); 
    f2 = CF*aS/Constants::pi*mu2*lr/v;
  }
  // small mass limit
  else {
    f1 = -CF*aS/Constants::pi/6.*
      ( - 6. - 24.*lmu*mu2 - 15.*mu4 - 12.*mu4*lmu - 24.*mu4*sqr(lmu) 
	+ 2.*mu4*sqr(Constants::pi) - 12.*mu2*mu4 - 96.*mu2*mu4*sqr(lmu) 
	+ 8.*mu2*mu4*sqr(Constants::pi) - 80.*mu2*mu4*lmu);
    fNS = - mu2/18.*( + 36.*lmu - 36. - 45.*mu2 + 216.*lmu*mu2 - 24.*mu2*sqr(Constants::pi) 
		      + 72.*mu2*sqr(lmu) - 22.*mu4 + 1032.*mu4 * lmu
		      - 96.*mu4*sqr(Constants::pi) + 288.*mu4*sqr(lmu));
    VNS = - mu2/1260.*(-6930. + 7560.*lmu + 2520.*mu - 16695.*mu2 + 1260.*mu2*sqr(Constants::pi) 
		       + 12600.*lmu*mu2 + 1344.*mu*mu2 - 52780.*mu4 + 36960.*mu4*lmu 
		       + 5040.*mu4*sqr(Constants::pi) - 12216.*mu*mu4);
    f2 = CF*aS*mu2/Constants::pi*( 2.*lmu + 4.*mu2*lmu + 2.*mu2 + 12.*mu4*lmu + 7.*mu4);
  }
  // add up bits for f1
  f1 += CF*aS/Constants::pi*(fNS+VNS);
  // now for the real correction
  // generate y
  double jac = 1.;
  double yminus = 0.; 
  double yplus  = 1.-2.*mu*(1.-mu)/(1.-2*mu2);
  double rhoymax = pow(yplus-yminus,1.-yPow_);
  double rhoy = UseRandom::rnd()*rhoymax;
  double y = yminus+pow(rhoy,1./(1.-yPow_));
  jac *= pow(y-yminus,yPow_)*rhoymax/(1.-yPow_);
  // generate z 
  double vt = sqrt(sqr(2.*mu2+(1.-2.*mu2)*(1.-y))-4.*mu2)/(1.-2.*mu2)/(1.-y);
  double zplus  = (1.+vt)*(1.-2.*mu2)*y/2./(mu2 +(1.-2.*mu2)*y);
  double zminus = (1.-vt)*(1.-2.*mu2)*y/2./(mu2 +(1.-2.*mu2)*y);
  double rhozmax = pow(zplus-zminus,1.-zPow_);
  double rhoz = UseRandom::rnd()*rhozmax;
  double z = zminus+pow(rhoz,1./(1.-zPow_));
  jac *= pow(z-zminus,zPow_)*rhozmax/(1.-zPow_);
  // calculate x1,x2,x3 and xT 
  double x2 = 1. - y*(1.-2.*mu2);
  double x1 = 1. - z*(x2-2.*mu2);
  double x3 = 2.-x1-x2;
  double xT = sqrt(sqr(x3) -0.25*sqr(sqr(x2)+sqr(x3)-sqr(x1))/(sqr(x2)-4.*mu2));
  // calculate the momenta
  Energy M = inpart.mass();
  Lorentz5Momentum pspect(ZERO,ZERO,-0.5*M*sqrt(max(sqr(x2)-4.*mu2,0.)),0.5*M*x2,M*mu); 
  double phi = UseRandom::rnd()*Constants::twopi;
  Lorentz5Momentum pemit (-0.5*M*xT*cos(phi),-0.5*M*xT*sin(phi),
			  0.5*M*sqrt(max(sqr(x1)-sqr(xT)-4.*mu2,0.)),0.5*M*x1,M*mu);
  Lorentz5Momentum pgluon( 0.5*M*xT*cos(phi), 0.5*M*xT*sin(phi),
			   0.5*M*sqrt(max(sqr(x3)-sqr(xT),0.)),0.5*M*x3,ZERO);
  if(abs(pspect.z()+pemit.z()-pgluon.z())/M<1e-6) 
    pgluon.setZ(-pgluon.z());
  else if(abs(pspect.z()-pemit.z()+pgluon.z())/M<1e-6) 
    pemit .setZ(- pemit.z());
  // loop over the possible emitting partons
  vector<cPDPtr> partons;
  partons.push_back(inpart.dataPtr());
  partons.push_back(decay[0]->dataPtr());
  partons.push_back(decay[1]->dataPtr());
  partons.push_back(gluon());
  double realwgt(0.);
  for(unsigned int iemit=0;iemit<2;++iemit) {
    // boost and rotate momenta
    LorentzRotation eventFrame( inpart.momentum().findBoostToCM() );
    Lorentz5Momentum spectator = eventFrame*decay[iemit]->momentum();
    eventFrame.rotateZ( -spectator.phi() );
    eventFrame.rotateY( -spectator.theta()  );
    eventFrame.invert();
    vector<Lorentz5Momentum> momenta;
    momenta.push_back(inpart.momentum());
    momenta.push_back(decay[0]->momentum());
    momenta.push_back(decay[1]->momentum());
    if(iemit==0) {
      momenta[2] = eventFrame*pspect;
      momenta[1] = eventFrame*pemit ;
    }
    else {
      momenta[1] = eventFrame*pspect;
      momenta[2] = eventFrame*pemit ;
    }
    momenta.push_back(eventFrame*pgluon);
    // calculate the weight
    if(1.-x1>1e-5 && 1.-x2>1e-5) 
      realwgt += meRatio(partons,momenta,iemit,true);
  }
  // total real emission contribution
  double realFact = 0.25*(1.-y)*jac*sqr(1.-2.*mu2)/
    sqrt(1.-4.*mu2)/Constants::twopi*2.*CF*aS*realwgt;
  // the born + virtual + real
  output = output*(1. + f1 + realFact) + f2*nlo;
  if(contrib_==2) output *=-1.;
  return max(output,0.);
}
