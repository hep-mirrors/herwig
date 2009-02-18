// -*- C++ -*-
//
// TopDecayMECorrection.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TopDecayMECorrection class.
//

#include "TopDecayMECorrection.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Shower/Base/Evolver.h"
#include "Herwig++/Shower/Base/KinematicsReconstructor.h"
#include "Herwig++/Shower/Base/PartnerFinder.h"
#include "ThePEG/Repository/EventGenerator.h"

extern "C" int isnan(double) throw();

using namespace Herwig;

TopDecayMECorrection::TopDecayMECorrection() 
  : _xg_sampling(1.5), 
    _initialenhance(1.),  _finalenhance(2.3), _useMEforT2(true) {}

IBPtr TopDecayMECorrection::clone() const {
  return new_ptr(*this);
}

IBPtr TopDecayMECorrection::fullclone() const {
  return new_ptr(*this);
}

void TopDecayMECorrection::persistentOutput(PersistentOStream & os) const {
    os << _initialenhance << _finalenhance << _xg_sampling << _useMEforT2;
}

void TopDecayMECorrection::persistentInput(PersistentIStream & is, int) {
    is >> _initialenhance >> _finalenhance >> _xg_sampling >> _useMEforT2;
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
     &TopDecayMECorrection::_initialenhance, 1.0, 1.0, 10000.0,
     false, false, Interface::limited);

  static Parameter<TopDecayMECorrection,double> interfaceFinalEnhancementFactor
    ("FinalEnhancementFactor",
     "The enhancement factor for final-state radiation in the shower to ensure"
     " the weight for the matrix element correction is less than one",
     &TopDecayMECorrection::_finalenhance, 1.6, 1.0, 1000.0,
     false, false, Interface::limited);

  static Parameter<TopDecayMECorrection,double> interfaceSamplingTopHardMEC
    ("SamplingTopHardMEC",
     "The importance sampling power for choosing an initial xg, "
     "to sample xg according to xg^-_xg_sampling",
     &TopDecayMECorrection::_xg_sampling, 1.5, 1.2, 2.0,
     false, false, Interface::limited);

  static Switch<TopDecayMECorrection,bool> interfaceUseMEForT2
    ("UseMEForT2",
     "Use the matrix element correction, if available to fill the T2"
     " region for the decay shower and don't fill using the shower",
     &TopDecayMECorrection::_useMEforT2, true, false, false);
  static SwitchOption interfaceUseMEForT2Shower
    (interfaceUseMEForT2,
     "Shower",
     "Use the shower to fill the T2 region",
     false);
  static SwitchOption interfaceUseMEForT2ME
    (interfaceUseMEForT2,
     "ME",
     "Use the Matrix element to fill the T2 region",
     true);

}

bool TopDecayMECorrection::canHandle(ShowerTreePtr tree, double & initial, 
				     double & final,EvolverPtr evolver) {
  // check radiation on
  if(!evolver->isFSRadiationON()) return false;
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
  if(abs(part[0]->id())==ParticleID::Wplus&&abs(part[1]->id())==ParticleID::b) {
    _ma=part[0]->mass();
    _mc=part[1]->mass();
  }
  else if(abs(part[1]->id())==ParticleID::Wplus&&abs(part[0]->id())==ParticleID::b) {
    _ma=part[1]->mass();
    _mc=part[0]->mass();
  }
  else 
    return false;
  // set the top mass
  _mt=tree->incomingLines().begin()->first->progenitor()->mass();
  // set the gluon mass
  _mg=getParticleData(ParticleID::g)->constituentMass();
  // set the radiation enhancement factors
  initial = _initialenhance;
  final   = _finalenhance;
  // reduced mass parameters
  _a=sqr(_ma/_mt);
  _g=sqr(_mg/_mt);
  _c=sqr(_mc/_mt);
  return true;
}

void TopDecayMECorrection::applyHardMatrixElementCorrection(ShowerTreePtr tree) {
  // Get b and a and put them in particle vector ba in that order...
  ParticleVector ba; 
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cit;
  for(cit=tree->outgoingLines().begin();cit!=tree->outgoingLines().end();++cit)
    ba.push_back(cit->first->copy());
  PPtr temp;
  if(abs(ba[0]->id())!=5) swap(ba[0],ba[1]);
  // Get the starting scales for the showers $\tilde{\kappa}_{b}$
  // $\tilde{\kappa}_{c}$. These are needed in order to know the boundary
  // of the dead zone.
  double ktb(0.),ktc(0.);
  map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator cjt;
  for(cjt = tree->incomingLines().begin();
      cjt!= tree->incomingLines().end();++cjt) {
    if(abs(cjt->first->progenitor()->id())==6)
      ktb=sqr(cjt->first->progenitor()->evolutionScale()/_mt); 
  }
  for(cit = tree->outgoingLines().begin();
      cit!= tree->outgoingLines().end();++cit) {
    if(abs(cit->first->progenitor()->id())==5)
      ktc=sqr(cit->first->progenitor()->evolutionScale()/_mt); 
  }
  if (ktb<=0.||ktc<=0.) {
    throw Exception() 
      << "TopDecayMECorrection::applyHardMatrixElementCorrection()"
      << " did not set ktb,ktc" 
      << Exception::abortnow; 
  }
  _ktb = ktb;
  _ktc = ktc;
  // Now decide if we get an emission into the dead region.
  // If there is an emission newfs stores momenta for a,c,g 
  // according to NLO decay matrix element. 
  vector<Lorentz5Momentum> newfs = applyHard(ba,ktb,ktc);
  // If there was no gluon emitted return.
  if(newfs.size()!=3) return;
  // Sanity checks to ensure energy greater than mass etc :)
  bool check = true; 
  tcPDPtr gluondata=getParticleData(ParticleID::g);
  if (newfs[0].e()<ba[0]->data().constituentMass()) check = false;
  if (newfs[1].e()<ba[1]->mass())                   check = false;
  if (newfs[2].e()<gluondata->constituentMass())    check = false;
  // Return if insane:
  if (!check) return;
  // Set masses in 5-vectors:
  newfs[0].setMass(ba[0]->mass());
  newfs[1].setMass(ba[1]->mass());
  newfs[2].setMass(ZERO);
  // The next part of this routine sets the colour structure.
  // To do this for decays we assume that the gluon comes from c!
  // First create new particle objects for c, a and gluon:
  PPtr newg = gluondata->produceParticle(newfs[2]);
  PPtr newc,newa;
  newc = ba[0]->data().produceParticle(newfs[0]);
  newa = ba[1];
  newa->set5Momentum(newfs[1]);
  // set the colour lines
  ColinePtr col;
  if(ba[0]->id()>0) {
    col=ba[0]->colourLine();
    col->addColoured(newg);
    newg->colourNeighbour(newc);
  }
  else {     
    col=ba[0]->antiColourLine();
    col->addAntiColoured(newg);
    newg->antiColourNeighbour(newc);
  }
  // change the existing quark and antiquark
  PPtr orig;
  for(cit=tree->outgoingLines().begin();cit!=tree->outgoingLines().end();++cit) {
    if(cit->first->progenitor()->id()==newc->id()) {
      // remove old particles from colour line
      if(newc->id()>0) {
	col->removeColoured(cit->first->copy());
	col->removeColoured(cit->first->progenitor());
      }
      else {
	col->removeAntiColoured(cit->first->copy());
	col->removeAntiColoured(cit->first->progenitor());
      }
      // insert new particles
      cit->first->copy(newc);
      ShowerParticlePtr sp(new_ptr(ShowerParticle(*newc,2,true)));
      cit->first->progenitor(sp);
      tree->outgoingLines()[cit->first]=sp;
      cit->first->perturbative(false);
      orig=cit->first->original();
    }
    else {
      cit->first->copy(newa);
      ShowerParticlePtr sp(new_ptr(ShowerParticle(*newa,2,true)));
      cit->first->progenitor(sp);
      tree->outgoingLines()[cit->first]=sp;
      cit->first->perturbative(true);
    }
  }
  // Add the gluon to the shower:
  ShowerParticlePtr   sg   =new_ptr(ShowerParticle(*newg,2,true));
  ShowerProgenitorPtr gluon=new_ptr(ShowerProgenitor(orig,newg,sg));
  gluon->perturbative(false);
  tree->outgoingLines().insert(make_pair(gluon,sg));
  if(!inTheDeadRegion(_xg,_xa,_ktb,_ktc)) {
    generator()->log()
      << "TopDecayMECorrection::applyHardMatrixElementCorrection()\n"
      << "Just found a point that escaped from the dead region!\n"
      << "   _xg: " << _xg << "   _xa: " << _xa 
      << "   newfs.size(): " << newfs.size() << endl;
  }
  tree->hardMatrixElementCorrection(true);
}

vector<Lorentz5Momentum> TopDecayMECorrection::
applyHard(const ParticleVector &p,double ktb, double ktc) 
{ 
  // ********************************* //
  // First we see if we get a dead     //
  // region event: _xa,_xg             //
  // ********************************* //
  vector<Lorentz5Momentum> fs; 
  // Return if there is no (NLO) gluon emission:
  
  double weight = getHard(ktb,ktc);
  if(weight>1.) {
    generator()->log() << "Weight greater than 1 for hard emission in "
		       << "TopDecayMECorrection::applyHard xg = " << _xg 
		       << " xa = " << _xa << "\n";
      weight=1.;
  }
  // Accept/Reject
  if (weight<UseRandom::rnd()||p.size()!= 2) return fs; 
   // Drop events if getHard returned a negative weight 
  // as in events that, somehow have escaped from the dead region
  // or, worse, the allowed region.
  if(weight<0.) return fs;
 
  // Calculate xc by momentum conservation:
  _xc = 2.-_xa-_xg;

  // ************************************ //
  // Now we get the boosts & rotations to //
  // go from lab to top rest frame with   //
  // a in the +z direction.               //
  // ************************************ //
  Lorentz5Momentum pa_lab,pb_lab,pc_lab,pg_lab;
  // Calculate momentum of b:
  pb_lab = p[0]->momentum() + p[1]->momentum(); 
   // Define/assign momenta of c,a and the gluon:
  if(abs(p[0]->id())==5) {
    pc_lab = p[0]->momentum(); 
    pa_lab = p[1]->momentum(); 
  } else {
    pc_lab = p[1]->momentum(); 
    pa_lab = p[0]->momentum(); 
  }
  // Calculate the boost to the b rest frame:
  SpinOneLorentzRotation rot0(pb_lab.findBoostToCM());
  // Calculate the rotation matrix to position a along the +z direction
  // in the rest frame of b and does a random rotation about z:
  SpinOneLorentzRotation    rot1 = rotateToZ(rot0*pa_lab);
  // Calculate the boost from the b rest frame back to the lab:
  // and the inverse of the random rotation about the z-axis and the 
  // rotation required to align a with +z:
  SpinOneLorentzRotation invrot = rot0.inverse()*rot1.inverse();

  // ************************************ //
  // Now we construct the momenta in the  //
  // b rest frame using _xa,_xg.          //
  // First we construct b, then c and g,  //
  // finally we generate a by momentum    //
  // conservation.                        //
  // ************************************ //
  Lorentz5Momentum pa_brf, pb_brf(_mt), pc_brf, pg_brf;
  // First we set the top quark to being on-shell and at rest.
  // Second we set the energies of c and g,
  pc_brf.setE(0.5*_mt*(2.-_xa-_xg));
  pg_brf.setE(0.5*_mt*_xg);
  // then their masses,
  pc_brf.setMass(_mc);
  pg_brf.setMass(ZERO);
  // Now set the z-component of c and g. For pg we simply start from
  // _xa and _xg, while for pc we assume it is equal to minus the sum
  // of the z-components of a (assumed to point in the +z direction) and g.
  double root=sqrt(_xa*_xa-4.*_a);
  pg_brf.setZ(_mt*(1.-_xa-_xg+0.5*_xa*_xg-_c+_a)/root);
  pc_brf.setZ(-1.*( pg_brf.z()+_mt*0.5*root));
  // Now set the y-component of c and g's momenta
  pc_brf.setY(ZERO);
  pg_brf.setY(ZERO);
  // Now set the x-component of c and g's momenta
  pg_brf.setX(sqrt(sqr(pg_brf.t())-sqr(pg_brf.z())));
  pc_brf.setX(-pg_brf.x());
  // Momenta b,c,g are now set. Now we obtain a from momentum conservation,
  pa_brf = pb_brf-pc_brf-pg_brf;
  pa_brf.setMass(pa_brf.m());
  pa_brf.rescaleEnergy();
 
  // ************************************ //
  // Now we orient the momenta and boost  //
  // them back to the original lab frame. //
  // ************************************ //
  // As in herwig6507 we assume that, in the rest frame
  // of b, we have aligned the W boson momentum in the 
  // +Z direction by rot1*rot0*pa_lab, therefore
  // we obtain the new pa_lab by applying:
  // invrot*pa_brf.
  pa_lab = invrot*pa_brf;    
  pb_lab = invrot*pb_brf;    
  pc_lab = invrot*pc_brf;    
  pg_lab = invrot*pg_brf;    
  fs.push_back(pc_lab); 
  fs.push_back(pa_lab); 
  fs.push_back(pg_lab); 
  return fs;
}

double TopDecayMECorrection::getHard(double ktb, double ktc) 
{
  // zero the variables
  _xg = 0.;    
  _xa = 0.;   
  _xc = 0.;
  // Get a phase space point in the dead region: 
  double volume_factor = deadRegionxgxa(ktb,ktc);
  // if outside region return -1
  if(volume_factor<0) return volume_factor;
  // Compute the weight for this phase space point:
  double weight = volume_factor*me(_xa,_xg)*(1.+_a-_c-_xa); 
  // Alpha_S and colour factors - this hard wired Alpha_S needs removing.
  weight *= (4./3.)/Constants::pi
    *(coupling()->value(_mt*_mt*_xg*(1.-_xa+_a-_c)
			/(2.-_xg-_xa-_c)));
  return weight; 
}

bool TopDecayMECorrection::softMatrixElementVeto(ShowerProgenitorPtr initial,
						 ShowerParticlePtr parent,Branching br)
{
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
      if(pt<initial->highestpT())
	veto=!UseRandom::rndbool(1./_initialenhance);
      // if hardest so far do calculation
      else
	{
	  // values of kappa and z
	  double z(br.kinematics->z()),kappa(sqr(br.kinematics->scale()/_mt));
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
              // This next `if' prevents the hardest emission from the 
              // top shower ever entering the so-called T2 region of the
              // phase space if that region is to be populated by the hard MEC.
              if(_useMEforT2&&xg>xgbcut(_ktb)) wgt = 0.;
	      if(wgt>1.) {
		generator()->log() << "Violation of maximum for initial-state "
				   << " soft veto in "
				   << "TopDecayMECorrection::softMatrixElementVeto"
				   << "xg = " << xg << " xa = " << xa 
				   << "weight =  " << wgt << "\n";
		wgt=1.;
	      }
	      // compute veto from weight
	      veto = !UseRandom::rndbool(wgt);
	    }
	  // if not vetoed reset max
	  if(!veto) initial->highestpT(pt);
	}
      // if vetoing reset the scale
      if(veto) parent->setEvolutionScale(br.kinematics->scale());
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
      if(pt<initial->highestpT()) return !UseRandom::rndbool(1./_finalenhance);
      // if hardest so far do calculation
      // values of kappa and z
      double z(br.kinematics->z()),kappa(sqr(br.kinematics->scale()/_mt));
      // momentum fractions
      double xa(1.+_a-_c-z*(1.-z)*kappa),r(0.5*(1.+_c/(1.+_a-xa))),root(sqr(xa)-4.*_a);
      if(root<0.) {
	  generator()->log() << "Imaginary root for final-state veto in "
			     << "TopDecayMECorrection::softMatrixElementVeto"
			     << "\nz =  " << z  << "\nkappa = " << kappa
			     << "\nxa = " << xa 
			     << "\nroot^2= " << root;
	  parent->setEvolutionScale(br.kinematics->scale());
	  return true;
      } 
      root=sqrt(root);
      double xg((2.-xa)*(1.-r)-(z-r)*root);
      // xfact (below) is supposed to equal xg/(1-z). 
      double xfact(z*kappa/2./(z*(1.-z)*kappa+_c)*(2.-xa-root)+root);
      // calculate the full result
      double f(me(xa,xg));
      // jacobian
      double J(z*root);
      double wgt(f*J*2.*kappa/(1.+sqr(z)-2.*_c/kappa/z)/sqr(xfact)/_finalenhance);
      if(wgt>1.) {
	generator()->log() << "Violation of maximum for final-state  soft veto in "
			   << "TopDecayMECorrection::softMatrixElementVeto"
			   << "xg = " << xg << " xa = " << xa 
			   << "weight =  " << wgt << "\n";
	wgt=1.;
      }
      // compute veto from weight
      veto = !UseRandom::rndbool(wgt);
      // if not vetoed reset max
      if(!veto) initial->highestpT(pt);
      // if vetoing reset the scale
      if(veto) parent->setEvolutionScale(br.kinematics->scale());
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


// This function is auxiliary to the xab function.
double TopDecayMECorrection::xgbr(int toggle) 
{ 
  return 1.+toggle*sqrt(_a)-_c*(1.-toggle*sqrt(_a))/(1.-_a);
}

// This function is auxiliary to the xab function.
double TopDecayMECorrection::ktr(double xgb, int toggle) 
{ 
  return 2.*xgb/
    (xgb+toggle*sqrt((1.-1./_a)
		     *(xgb-xgbr( 1))
		     *(xgb-xgbr(-1))));
}

// Function xab determines xa (2*W energy fraction) for a given value
// of xg (2*gluon energy fraction) and kappa tilde (q tilde squared over
// m_top squared). Hence this function allows you to draw 1: the total
// phase space volume in the xa vs xg plane 2: for a given value of 
// kappa tilde (i.e. starting evolution scale) the associated contour 
// in the xa vs xg plane (and hence the regions that either shower can 
// populate). This calculation is done assuming the emission came from
// the top quark i.e. kappa tilde here is the q tilde squared of the TOP
// quark divided by m_top squared. 
double TopDecayMECorrection::xab(double xgb,
                                        double kt, int toggle) 
{ 
  double xab;
  if(toggle==2) {
    // This applies for g==0.&&kt==ktr(a,c,0.,xgb,1).
    xab = -2.*_a*(xgb-2.)/(1.+_a-_c-xgb);
  } else if(toggle==1) {
    // This applies for kt==1&&g==0.
    double lambda = sqrt(sqr(xgb-1.+_a+_c)-4.*_a*_c);
    xab = (0.5/(kt-xgb))*(kt*(1.+_a-_c-xgb)-lambda)
      + (0.5/(kt+xgb*(1.-kt)))*(kt*(1.+_a-_c-xgb)+lambda);
  } else {
    // This is the form of xab FOR _g=0.
    double ktmktrpktmktrm = ( sqr(xgb*kt-2.*xgb)
			      -kt*kt*(1.-1./_a)*(xgb-xgbr( 1))
			      *(xgb-xgbr(-1))
			      )/
      (xgb*xgb-(1.-1./_a)*(xgb-xgbr( 1))
       *(xgb-xgbr(-1))
       );
    double lambda = sqrt((sqr(1.-_a-_c-xgb)-4.*_a*_c)*
			 ktmktrpktmktrm);
    xab = (0.5/(kt-xgb))*(kt*(1.+_a-_c-xgb)-lambda)
      + (0.5/(kt+xgb*(1.-kt)))*(kt*(1.+_a-_c-xgb)+lambda);
  }
  if(isnan(xab)) {
    double ktmktrpktmktrm = ( sqr(xgb*kt-2.*(xgb-_g))
			      -kt*kt*(1.-1./_a)*(xgb-xgbr( 1)-_g/(1.+sqrt(_a)))
			      *(xgb-xgbr(-1)-_g/(1.-sqrt(_a)))
			      )/
      (xgb*xgb-(1.-1./_a)*(xgb-xgbr( 1)-_g/(1.+sqrt(_a)))
       *(xgb-xgbr(-1)-_g/(1.-sqrt(_a)))
       );
    double lambda = sqrt((xgb-1.+sqr(sqrt(_a)+sqrt(_c-_g)))
			 *(xgb-1.+sqr(sqrt(_a)-sqrt(_c-_g)))*
			 ktmktrpktmktrm);
    xab = (0.5/(kt-xgb+_g))*(kt*(1.+_a-_c+_g-xgb)-lambda)
      + (0.5/(kt+xgb*(1.-kt)-_g))*(kt*(1.+_a-_c+_g-xgb)+lambda);
    if(isnan(xab)) 
	throw Exception() << "TopMECorrection::xab complex x_a value.\n"
			  << "  xgb    = " << xgb    << "\n"
			  << "  xab    = " << xab    << "\n"
			  << "  toggle = " << toggle << "\n"
			  << "  ktmktrpktmktrm = "   << ktmktrpktmktrm 
			  << Exception::eventerror;
  }
  return xab;
}

// xgbcut is the point along the xg axis where the upper bound on the 
// top quark (i.e. b) emission phase space goes back on itself in the 
// xa vs xg plane i.e. roughly mid-way along the xg axis in
// the xa vs xg Dalitz plot.
double TopDecayMECorrection::xgbcut(double kt) 
{ 
  double lambda2 = 1.+_a*_a+_c*_c-2.*_a-2.*_c-2.*_a*_c; 
  double num1    = kt*kt*(1.-_a-_c);
  double num2    = 2.*kt*sqrt(_a*(kt*kt*_c+lambda2*(kt-1.)));
  return (num1-num2)/(kt*kt-4.*_a*(kt-1.));
}

double TopDecayMECorrection::xaccut(double kt) 
{ 
    return 1.+_a-_c-0.25*kt;
}

double TopDecayMECorrection::z(double xac, double kt, 
                                      int toggle1, int toggle2) 
{ 
  double z = -1.0;
  if(toggle2==0) { 
    z = (kt+toggle1*sqrt(kt*(kt-4.*(1.+_a-_c-xac))))/(2.*kt); 
  } else if(toggle2==1) {
    z = ((1.+_a+_c-xac)+toggle1*(1.+_a-_c-xac))
      /(2.*(1.+_a-xac));
  } else if(toggle2==2) {
    z = 0.5;
  } else {
    throw Exception() << "Cannot determine z in TopDecayMECorrection::z()"
		      << Exception::eventerror;
  }
  return z;
}

double TopDecayMECorrection::xgc(double xac, double kt, 
                                        int toggle1, int toggle2) 
{ 
    return (2.-xac)*(1.-0.5*(1.+_c/(1.+_a-xac)))
          -(z(xac,kt,toggle1,toggle2)-0.5*(1.+_c/(1.+_a-xac)))
          *sqrt(xac*xac-4.*_a);
}

double TopDecayMECorrection::xginvc0(double xg , double kt) 
{ 
  // The function xg(kappa_tilde_c,xa) surely, enough, draws a  
  // line of constant kappa_tilde_c in the xg, xa Dalitz plot. 
  // Such a function can therefore draw the upper and lower 
  // edges of the phase space for emission from c (the b-quark).
  // However, to sample the soft part of the dead zone effectively
  // we want to generate a value of xg first and THEN distribute
  // xa in the associated allowed part of the dead zone. Hence, the 
  // function we want, to define the dead zone in xa for a given
  // xg, is the inverse of xg(kappa_tilde_c,xa). The full expression 
  // for xg(kappa_tilde_c,xa) is complicated and, sure enough, 
  // does not invert. Therefore we try to overestimate the size
  // of the dead zone initially, rejecting events which do not 
  // fall exactly inside it afterwards, with the immediate aim 
  // of getting an approximate version of xg(kappa_tilde_c,xa) 
  // that can  be inverted. We do this by simply setting c=0 i.e.
  // the b-quark mass to zero (and the gluon mass of course), in 
  // the full expression xg(...). The result of inverting this 
  // function is the output of this routine (a value of xa) hence 
  // the name xginvc0. xginvc0 is calculated to be,
  // xginvc0 = (1./3.)*(1.+a+pow((U+sqrt(4.*V*V*V+U*U))/2.,1./3.)
  //                      -V*pow(2./(U+sqrt(4.*V*V*V+U*U)),1./3.)
  //                   )
  // U = 2.*a*a*a - 66.*a*a + 9.*a*kt*xg + 18.*a*kt
  //   - 66.*a + 27.*kt*xg*xg - 45.*kt*xg +18.*kt +2. ;
  // V = -1.-a*a-14.*a-3.kt*xg+3.*kt;
  // This function, as with many functions in this ME correction,
  // is plagued by cuts that have to handled carefully in numerical 
  // implementation. We have analysed the cuts and hence we implement 
  // it in the following way, with a series of 'if' statements. 
  //
  // A useful -definition- to know in deriving the v<0 terms is
  // that tanh^-1(z) = 0.5*(log(1.+z)-log(1.-z)).
  double u,v,output;
  u = 2.*_a*_a*_a-66.*_a*_a
     +9.*xg*kt*_a+18.*kt*_a
     -66.*_a+27.*xg*xg*kt
     -45.*xg*kt+18.*kt+2.;
  v = -_a*_a-14.*_a-3.*xg*kt+3.*kt-1.;
  double u2=u*u,v3=v*v*v;
  if(v<0.) {
    if(u>0.&&(4.*v3+u2)<0.)      output = cos(  atan(sqrt(-4.*v3-u2)/u)/3.);
    else if(u>0.&&(4.*v3+u2)>0.) output = cosh(atanh(sqrt( 4.*v3+u2)/u)/3.);
    else                         output = cos(( atan(sqrt(-4.*v3-u2)/u)
					       +Constants::pi)/3.);
    output *= 2.*sqrt(-v);
  } else {
    output = sinh(log((u+sqrt(4.*v3+u2))/(2.*sqrt(v3)))/3.);
    output *= 2.*sqrt(v);
  }
  if(isnan(output)||isinf(output)) {
      throw Exception() << "TopMECorrection::xginvc0:\n"
	  << "possible numerical instability detected.\n"
	  << "\n v = " <<  v << "   u = " << u << "\n4.*v3+u2 = " << 4.*v3+u2
	  << "\n_a = " << _a << "  ma = " << sqrt(_a*_mt*_mt/GeV2) 
	  << "\n_c = " << _c << "  mc = " << sqrt(_c*_mt*_mt/GeV2) 
	  << "\n_g = " << _g << "  mg = " << sqrt(_g*_mt*_mt/GeV2) 
	  << Exception::eventerror;
  }
  return ( 1.+_a +output)/3.;
}

double TopDecayMECorrection::approxDeadMaxxa(double xg,double ktb,double ktc) 
{
  double maxxa(0.);
  double x  = min(xginvc0(xg,ktc),
		  xab(xg,(2.*xg-2.*_g)/(xg-sqrt(xg*xg-4.*_g)),0));
  double y(-9999999999.);
  if(xg>2.*sqrt(_g)&&xg<=xgbcut(ktb)) {
      y = max(xab(xg,ktb,0),xab(xg,1.,1));
  } else if(xg>=xgbcut(ktb)&&xg<=1.-sqr(sqrt(_a)+sqrt(_c))) {
      y = max(xab(xg,ktr(xg,1),2),xab(xg,1.,1));
  }
  if(xg>2.*sqrt(_g)&&xg<=1.-sqr(sqrt(_a)+sqrt(_c))) {
    if(x>=y) { maxxa =  x       ; }
    else     { maxxa = -9999999.; }
  } else {
    maxxa = -9999999.;
  }
  return maxxa;
}

double TopDecayMECorrection::approxDeadMinxa(double xg,double ktb,double ktc) 
{
  double minxa(0.);
  double x  = min(xginvc0(xg,ktc),
		  xab(xg,(2.*xg-2.*_g)/(xg-sqrt(xg*xg-4.*_g)),0));
  double y(-9999999999.);
  if(xg>2.*sqrt(_g)&&xg<=xgbcut(ktb)) {
      y = max(xab(xg,ktb,0),xab(xg,1.,1));
  } else if(xg>=xgbcut(ktb)&&xg<=1.-sqr(sqrt(_a)+sqrt(_c))) {
      if(_useMEforT2) y = xab(xg,1.,1);
      else            y = max(xab(xg,ktr(xg,1),2),xab(xg,1.,1));
  }
  if(xg>2.*sqrt(_g)&&xg<=1.-sqr(sqrt(_a)+sqrt(_c))) {
      if(x>=y) { minxa =  y  ; }
      else     { minxa = 9999999.; }
  } else {
      minxa = 9999999.;
  }
  return minxa;
}

// This function returns true if the phase space point (xg,xa) is in the 
// kinematically allowed phase space.
bool TopDecayMECorrection::inTheAllowedRegion(double xg , double xa)
{
    bool output(true);
    if(xg<2.*sqrt(_g)||xg>1.-sqr(sqrt(_a)+sqrt(_c)))      output = false;
    if(xa<xab(xg,1.,1))                                   output = false;
    if(xa>xab(xg,(2.*xg-2.*_g)/(xg-sqrt(xg*xg-4.*_g)),0)) output = false;
    return output;
}

// This function returns true if the phase space point (xg,xa) is in the 
// approximate (overestimated) dead region.
bool TopDecayMECorrection::inTheApproxDeadRegion(double xg , double xa,
                                                        double ktb, double ktc)
{
    bool output(true);
    if(!inTheAllowedRegion(xg,xa))       output = false;
    if(xa<approxDeadMinxa(xg,ktb,ktc))   output = false;
    if(xa>approxDeadMaxxa(xg,ktb,ktc))   output = false;
    return output;
}

// This function returns true if the phase space point (xg,xa) is in the 
// dead region.
bool TopDecayMECorrection::inTheDeadRegion(double xg , double xa,
                                                  double ktb, double ktc)
{
    bool output(true);
    if(!inTheApproxDeadRegion(xg,xa,ktb,ktc)) output = false;
    if(xa>xaccut(ktc)) {
	if(xg<xgc(max(xaccut(ktc),2.*sqrt(_a)),ktc, 1,2)&&
           xg>xgc(xa,ktc, 1,0)) { output = false; } 
	if(xg>xgc(max(xaccut(ktc),2.*sqrt(_a)),ktc,-1,2)&&
           xg<xgc(xa,ktc,-1,0)) { output = false; } 
    } 
    return output;
}

// This function attempts to generate a phase space point in the dead
// region and returns the associated phase space volume factor needed for
// the associated event weight.
double TopDecayMECorrection::deadRegionxgxa(double ktb,double ktc)
{
  _xg=0.;
  _xa=0.;
  // Here we set limits on xg and generate a value inside the bounds.
  double xgmin(2.*sqrt(_g)),xgmax(1.-sqr(sqrt(_a)+sqrt(_c)));
  // Generate _xg.
  if(_xg_sampling==2.) {
      _xg=xgmin*xgmax/(xgmin+UseRandom::rnd()*(xgmax-xgmin));
  } else {
      _xg=xgmin*xgmax/pow((  pow(xgmin,_xg_sampling-1.)
			    + UseRandom::rnd()*(pow(xgmax,_xg_sampling-1.)
					       -pow(xgmin,_xg_sampling-1.))
			   ),1./(_xg_sampling-1.));
  }
  // Here we set the bounds on _xa for given _xg.
  if(_xg<xgmin||xgmin>xgmax) 
      throw Exception() << "TopMECorrection::deadRegionxgxa:\n"
			<< "upper xg bound is less than the lower xg bound.\n"
			<< "\n_xg         = " << _xg 
			<< "\n2.*sqrt(_g) = " << 2.*sqrt(_g) 
			<< "\n_a  = " << _a  << "  ma = " << sqrt(_a*_mt*_mt/GeV2) 
			<< "\n_c  = " << _c  << "  mc = " << sqrt(_c*_mt*_mt/GeV2) 
			<< "\n_g  = " << _g  << "  mg = " << sqrt(_g*_mt*_mt/GeV2) 
			<< Exception::eventerror;
  double xamin(approxDeadMinxa(_xg,ktb,ktc));
  double xamax(approxDeadMaxxa(_xg,ktb,ktc));
  // Are the bounds sensible? If not return.
  if(xamax<=xamin) return -1.;
  _xa=1.+_a-(1.+_a-xamax)*pow((1.+_a-xamin)/(1.+_a-xamax),UseRandom::rnd());
  // If outside the allowed region return -1.
  if(!inTheDeadRegion(_xg,_xa,ktb,ktc)) return -1.;
  // The integration volume for the weight
  double xg_vol,xa_vol; 
  if(_xg_sampling==2.) {
      xg_vol = (xgmax-xgmin)
             / (xgmax*xgmin);
  } else {
      xg_vol = (pow(xgmax,_xg_sampling-1.)-pow(xgmin,_xg_sampling-1.))
	     / ((_xg_sampling-1.)*pow(xgmax*xgmin,_xg_sampling-1.));
  }
  xa_vol = log((1.+_a-xamin)/(1.+_a-xamax));
  // Here we return the integral volume factor multiplied by the part of the 
  // weight left over which is not included in the BRACES function, i.e.
  // the part of _xg^-2 which is not absorbed in the integration measure.
  return xg_vol*xa_vol*pow(_xg,_xg_sampling-2.);
}

LorentzRotation TopDecayMECorrection::rotateToZ(Lorentz5Momentum v) {
  // compute the rotation matrix
  LorentzRotation trans;
  // rotate so in z-y plane
  trans.rotateZ(-atan2(v.y(),v.x()));
  // rotate so along Z
  trans.rotateY(-acos(v.z()/v.vect().mag()));
  // generate random rotation
  double c,s,cs;
  do
    {
      c = 2.*UseRandom::rnd()-1.;
      s = 2.*UseRandom::rnd()-1.;
      cs = c*c+s*s;
    }
  while(cs>1.||cs==0.);
  double cost=(c*c-s*s)/cs,sint=2.*c*s/cs;
  // apply random azimuthal rotation
  trans.rotateZ(atan2(sint,cost));
  return trans;
}
