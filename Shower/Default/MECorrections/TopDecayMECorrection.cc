// -*- C++ -*-
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

using namespace Herwig;

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
      ktb=sqr(cjt->first->progenitor()->evolutionScales()[0]/_mt); 
  }
  for(cit = tree->outgoingLines().begin();
      cit!= tree->outgoingLines().end();++cit) {
    if(abs(cit->first->progenitor()->id())==5)
      ktc=sqr(cit->first->progenitor()->evolutionScales()[0]/_mt); 
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
  newfs[2].setMass(0.*MeV);
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
  pg_brf.setMass(0.*MeV);
  // Now set the z-component of c and g. For pg we simply start from
  // _xa and _xg, while for pc we assume it is equal to minus the sum
  // of the z-components of a (assumed to point in the +z direction) and g.
  double root=sqrt(_xa*_xa-4.*_a);
  pg_brf.setZ(_mt*(1.-_xa-_xg+0.5*_xa*_xg-_c+_a)/root);
  pc_brf.setZ(-1.*( pg_brf.z()+_mt*0.5*root));
  // Now set the y-component of c and g's momenta
  pc_brf.setY(0.*MeV);
  pg_brf.setY(0.*MeV);
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
      if(veto) parent->setEvolutionScale(ShowerIndex::QCD,br.kinematics->scale());
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
	  parent->setEvolutionScale(ShowerIndex::QCD,br.kinematics->scale());
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
      if(veto) parent->setEvolutionScale(ShowerIndex::QCD,br.kinematics->scale());
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
