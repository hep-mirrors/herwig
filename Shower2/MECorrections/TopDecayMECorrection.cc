// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TopDecayMECorrection class.
//

#include "TopDecayMECorrection.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "TopDecayMECorrection.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

TopDecayMECorrection::~TopDecayMECorrection() {}

void TopDecayMECorrection::persistentOutput(PersistentOStream & os) const {
  os << _initialenhance << _finalenhance;
}

void TopDecayMECorrection::persistentInput(PersistentIStream & is, int) {
  is >> _initialenhance >> _finalenhance;
}

ClassDescription<TopDecayMECorrection> TopDecayMECorrection::initTopDecayMECorrection;
// Definition of the static class description member.

void TopDecayMECorrection::Init() {

  static ClassDocumentation<TopDecayMECorrection> documentation
    ("The TopDecayMECorrection class implements the matrix element correction for"
     " top decay");

  static Parameter<TopDecayMECorrection,double> interfaceEnhancementFactor
    ("InitialEnhancementFactor",
     "The enhancement factor for initial-state radiation in the shower to ensure"
     " the weight for the matrix element correction is less than one.",
     &TopDecayMECorrection::_initialenhance, 700.0, 1.0, 10000.0,
     false, false, Interface::limited);

  static Parameter<TopDecayMECorrection,double> interfaceFinalEnhancementFactor
    ("FinalEnhancementFactor",
     "The enhancement factor for final-state radiation in the shower to ensure"
     " the weight for the matrix element correction is less than one",
     &TopDecayMECorrection::_finalenhance, 1.2, 1.0, 1000.0,
     false, false, Interface::limited);

}

bool TopDecayMECorrection::canHandle(ShowerTreePtr tree)
{
  // one incoming particle
  if(tree->incomingLines().size()!=1) return false;
  // it should be top
  if(abs(tree->incomingLines().begin()->first->progenitor()->id())!=ParticleID::t)
    return false;
  // two outgoing particles
  if(tree->outgoingLines().size()!=2) return false;
  // check the outgoing particles
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cit;
  ShowerParticlePtr part[2];
  unsigned int ix(0);
  for(cit=tree->outgoingLines().begin();cit!=tree->outgoingLines().end();++cit)
    {part[ix]=cit->first->progenitor();++ix;}
  // check the final-state particles and get the masses
  if(abs(part[0]->id())==ParticleID::Wplus&&abs(part[1]->id())==ParticleID::b)
    {
      _mw=part[0]->mass();
      _mb=part[1]->mass();
    }
  else if(abs(part[1]->id())==ParticleID::Wplus&&abs(part[0]->id())==ParticleID::b)
    {
      _mw=part[1]->mass();
      _mb=part[0]->mass();
    }
  else 
    return false;
  // set the top mass
  _mt=tree->incomingLines().begin()->first->progenitor()->mass();
  // set the radiation enhancement factors
  showerVariables()->initialStateRadiationEnhancementFactor(_initialenhance);
  showerVariables()->finalStateRadiationEnhancementFactor(_finalenhance);
  // parameters
  _a=sqr(_mw/_mt);
  _c=sqr(_mb/_mt);
  //_c=0.;
  return true;
}

void TopDecayMECorrection::applyHardMatrixElementCorrection(ShowerTreePtr)
{
}

bool TopDecayMECorrection::softMatrixElementVeto(ShowerProgenitorPtr initial,
						 ShowerParticlePtr parent,Branching br)
{
  //cerr << "testing in top decay soft correction " 
  //     << *initial->progenitor() << "\n" << *parent << "\n";
  // check if we need to apply the full correction
  unsigned int id[2]={abs(initial->progenitor()->id()),abs(parent->id())};
  // the initial-state correction
  if(id[0]==ParticleID::t&&id[1]==ParticleID::t)
    {
      Energy pt=br.kinematics->pT();
      // check if hardest so far
      // if not just need to remove effect of enhancement
      bool veto(false);
      // if not hardest so far
      if(pt<initial->pT())
	veto=!UseRandom::rndbool(1./_initialenhance);
      // if hardest so far do calculation
      else
	{
	  // values of kappa and z
	  double z(br.kinematics->z()),kappa(sqr(br.kinematics->qtilde()/_mt));
	  // parameters for the translation
	  double w(1.-(1.-z)*(kappa-1.)),u(1.+_a-_c-(1.-z)*kappa),v(sqr(u)-4.*_a*w*z);
	  // veto if outside phase space
	  if(v<0.) 
	    veto=true;
	  // otherwise calculate the weight
	  else
	    {
	      v = sqrt(v);
	      double xa((0.5*(u+v)/w+0.5*(u-v)/z)),xg((1.-z)*kappa);
	      double f(me(xa,xg)),
		J(0.5*(u+v)/sqr(w)-0.5*(u-v)/sqr(z)+_a*sqr(w-z)/(v*w*z));
	      double wgt(f*J*2./kappa/(1.+sqr(z)-2.*z/kappa)/_initialenhance);
	      if(wgt>1.)
		{
		  cerr << "testing violates max initial " << xg 
		       << " " << xa << " " << wgt << " " << _initialenhance << endl;
		  wgt=1.;
		}
	      // compute veto from weight
	      veto = !UseRandom::rndbool(wgt);
	    }
	  // if not vetoed reset max
	  if(!veto) initial->pT(pt);
	}
      // if vetoing reset the scale
      if(veto) parent->setEvolutionScale(ShowerIndex::QCD,br.kinematics->qtilde());
      // return the veto
      return veto;
    }
  // final-state correction
  else if(id[0]==ParticleID::b&&id[1]==ParticleID::b)
    {
      Energy pt=br.kinematics->pT();
      // check if hardest so far
      // if not just need to remove effect of enhancement
      bool veto(false);
      // if not hardest so far
      if(pt<initial->pT()) return !UseRandom::rndbool(1./_finalenhance);
      // if hardest so far do calculation
      // values of kappa and z
      double z(br.kinematics->z()),kappa(sqr(br.kinematics->qtilde()/_mt));
      // momentum fractions
      double xa(1.+_a-_c-z*(1.-z)*kappa),r(0.5*(1.+_c/(1.+_a-xa))),
	xg((2.-xa)*(1.-r)-(z-r)*(sqr(xa)-4.*_a));
      double root(sqr(xa)-4.*_a),xfact(z*kappa/2./(z*(1.-z)*kappa+_c)*(2.-xa-root)+root);
      // calculate the full result
      double f(me(xa,xg));
      // jacobian
      double J(z*root);
      double wgt(f*J*2.*kappa/(1.+sqr(z)-2.*_c/kappa/z)/sqr(xfact)/_finalenhance);
      if(wgt>1.)
	{
	  cerr << "testing violates max final " << xg 
	       << " " << xa << " " << wgt << endl;
	  wgt=1.;
	}
      // compute veto from weight
      veto = !UseRandom::rndbool(wgt);
      // if not vetoed reset max
      if(!veto) initial->pT(pt);
      // if vetoing reset the scale
      if(veto) parent->setEvolutionScale(ShowerIndex::QCD,br.kinematics->qtilde());
      // return the veto
      return veto;
    }
  // otherwise don't veto
  else return !UseRandom::rndbool(1./_finalenhance);
}

double TopDecayMECorrection::me(double xw,double xg)
{
  double prop(1.+_a-_c-xw),xg2(sqr(xg));
  double lambda=sqrt(1.+_a*_a+_c*_c-2.*_a-2.*_c-2.*_a*_c);
  double denom=(1.-2*_a*_a+_a+_c*_a+_c*_c-2.*_c);
  double wgt=-_c*xg2/prop+(1.-_a+_c)*xg-(prop*(1 - xg)+xg2)
    +(0.5*(1.+2.*_a+_c)*sqr(prop-xg)*xg+2.*_a*prop*xg2)/denom;
  return wgt/(lambda*prop);
}
