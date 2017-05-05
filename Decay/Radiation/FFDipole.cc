// -*- C++ -*-
//
// FFDipole.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFDipole class.
//

#include "FFDipole.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Utilities/Debug.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "YFSFormFactors.h"
#include "Herwig/Decay/DecayPhaseSpaceMode.h"
#include "Herwig/Decay/DecayIntegrator.h"

using namespace Herwig;

void FFDipole::persistentOutput(PersistentOStream & os) const {
  os << ounit(_emin,GeV) << ounit(_eminrest,GeV) << ounit(_eminlab,GeV) 
     << _maxwgt << _weightOutput
     << _mode << _maxtry << _energyopt << _betaopt << _dipoleopt;
}

void FFDipole::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_emin,GeV) >> iunit(_eminrest,GeV) >> iunit(_eminlab,GeV) 
     >> _maxwgt >> _weightOutput
     >> _mode >> _maxtry >> _energyopt >> _betaopt >> _dipoleopt;
}

FFDipole::~FFDipole() {}

ClassDescription<FFDipole> FFDipole::initFFDipole;
// Definition of the static class description member.

void FFDipole::Init() {

  static ClassDocumentation<FFDipole> documentation
    ("The FFDipole class implements the final-final dipole for the SOPTHY algorithm");

  static Switch<FFDipole,unsigned int> interfaceUnWeight
    ("UnWeight",
     "Control the type of unweighting to perform, only one should be used the"
     " other options are for debugging purposes.",
     &FFDipole::_mode, 1, false, false);
  static SwitchOption interfaceUnWeightNoUnweighting
    (interfaceUnWeight,
     "NoUnweighting",
     "Perform no unweighting",
     0);
  static SwitchOption interfaceUnWeightAllWeights
    (interfaceUnWeight,
     "AllWeights",
     "Include all the weights",
     1);
  static SwitchOption interfaceUnWeightNoJacobian
    (interfaceUnWeight,
     "NoJacobian",
     "Only include the dipole and YFS weights",
     2);
  static SwitchOption interfaceUnWeightDipole
    (interfaceUnWeight,
     "Dipole",
     "Only include the dipole weight",
     3);
  static SwitchOption interfaceUnWeightYFS
    (interfaceUnWeight,
     "YFS",
     "Only include the YFS weight",
     4);
  static SwitchOption interfaceUnWeightNLO
    (interfaceUnWeight,
     "NLO",
     "Weight to get the stict NLO rate",
     5);

  static Parameter<FFDipole,unsigned int> interfaceMaximumTries
    ("MaximumTries",
     "Maximum number of attempts to unweight",
     &FFDipole::_maxtry, 500, 10, 100000,
     false, false, Interface::limited);

  static Parameter<FFDipole,Energy> interfaceMinimumEnergyBoosted
    ("MinimumEnergyBoosted",
     "The minimum energy of the photons in the boosted frame in which"
     " they are generated.",
     &FFDipole::_emin, MeV, 1.e-6*MeV, ZERO, 100.0*MeV,
     false, false, Interface::limited);

  static Parameter<FFDipole,Energy> interfaceMinimumEnergyRest
    ("MinimumEnergyRest",
     "The minimum energy of the photons in the rest frame of the decaying particle",
     &FFDipole::_eminrest, MeV, 100.0*MeV, 1.0*MeV, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<FFDipole,Energy> interfaceMinimumEnergyLab
    ("MinimumEnergyLab",
     "The minimum energy of the photons in the lab frame",
     &FFDipole::_eminlab, MeV, 100.0*MeV, 1.0*MeV, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<FFDipole,double> interfaceMaximumWeight
    ("MaximumWeight",
     "The maximum weight for unweighting",
     &FFDipole::_maxwgt, 7.0, 0.0, 100.0,
     false, false, Interface::limited);

  static Switch<FFDipole,unsigned int> interfaceEnergyCutOff
    ("EnergyCutOff",
     "The type of cut-off on the photon energy to apply",
     &FFDipole::_energyopt, 1, false, false);
  static SwitchOption interfaceEnergyCutOffBoostedFrame
    (interfaceEnergyCutOff,
     "BoostedFrame",
     "Only apply cut-off in boosted frame",
     0);
  static SwitchOption interfaceEnergyCutOffRestFrame
    (interfaceEnergyCutOff,
     "RestFrame",
     "Apply cut-off in rest frame",
     1);
  static SwitchOption interfaceEnergyCutOff2
    (interfaceEnergyCutOff,
     "LabFrame",
     "Apply cut-off in lab frame",
     2);

  static Switch<FFDipole,unsigned int> interfaceBetaOption
    ("BetaOption",
     "Option for the inclusive of the higher beta coefficients",
     &FFDipole::_betaopt, 4, false, false);
  static SwitchOption interfaceBetaOptionNone
    (interfaceBetaOption,
     "None",
     "No higher betas included",
     0);
  static SwitchOption interfaceBetaOptionCollinear
    (interfaceBetaOption,
     "Collinear",
     "Include the collinear approx",
     1);
  static SwitchOption interfaceBetaOptionCollinearVirtA
    (interfaceBetaOption,
     "CollinearVirtualA",
     "Include the collinear approx with virtual corrections",
     2);
  static SwitchOption interfaceBetaOptionCollinearVirtB
    (interfaceBetaOption,
     "CollinearVirtualB",
     "Include the collinear approx with virtual corrections",
     3);
  static SwitchOption interfaceBetaOptionExact
    (interfaceBetaOption,
     "Exact",
     "Include the exact higher order terms if available",
     4);

  static Switch<FFDipole,unsigned int> interfaceDipoleOption
    ("DipoleOption",
     "Option for generating the primary dipole distribution",
     &FFDipole::_dipoleopt, 0, false, false);
  static SwitchOption interfaceDipoleOptionNoMass
    (interfaceDipoleOption,
     "NoMass",
     "Don't include the mass terms in the primary distribution",
     0);
  static SwitchOption interfaceDipoleOptionMass
    (interfaceDipoleOption,
     "Mass",
     "Include the mass terms in the primary distribution",
     1);

  static Switch<FFDipole,bool> interfaceWeightOutput
    ("WeightOutput",
     "Whether or not to output the average weight for testing",
     &FFDipole::_weightOutput, false, false, false);
  static SwitchOption interfaceWeightOutputNo
    (interfaceWeightOutput,
     "No",
     "Don't output the average",
     false);
  static SwitchOption interfaceWeightOutputYes
    (interfaceWeightOutput,
     "Yes",
     "Output the average",
     true);

}

void FFDipole::printDebugInfo(const Particle & p,
			      const ParticleVector & children,
			      double wgt) const {
  generator()->log() << "Input masses " 
		     << p.mass()/GeV << " -> " 
		     << children[0]->mass()/GeV << " " 
		     << children[1]->mass()/GeV << '\n'; 
  generator()->log() << "Momenta\n";
  generator()->log() << "parent " << p.momentum()/GeV << '\n';
  for(unsigned int ix=0;ix<2;++ix)
    generator()->log() << "charged " << ix << " " 
		       << _qnewlab[ix]/GeV << " " 
		       << children[ix]->momentum()/GeV << '\n';
  for(unsigned int ix=0;ix<_multiplicity;++ix) {
    generator()->log() << "photons " << ix << " "
		       << "phocut " << _photcut[ix] << ' '
		       << _llab[ix]/GeV << '\n';
  }
  generator()->log() << "wgt         : " << wgt          << '\n';
  generator()->log() << "_mewgt      : " << _mewgt       << '\n';
  generator()->log() << "_jacobianwgt: " << _jacobianwgt << '\n';
  generator()->log() << "_yfswgt     : " << _yfswgt      << '\n';
  generator()->log() << "_dipolewgt  : " << _dipolewgt   << '\n';
  generator()->log() << "dipoleopt   : " << _dipoleopt   << '\n';
}

ParticleVector FFDipole::generatePhotons(const Particle & p,
					 ParticleVector children,
					 tDecayIntegratorPtr decayer) {
  _parent = const_ptr_cast<tPPtr>(&p);
  // set the decayer
  _decayer=decayer;
  // set parameters which won't change in the event loop
  // masses of the particles
  _m[0] = p.mass(); 
  _m[1] = children[0]->mass();
  _m[2] = children[1]->mass();
  // set the maximum photon energy (exact - no approximations here).
  _emax=(0.5*(_m[0]-sqr(_m[1]+_m[2])/_m[0]))*_m[0]/(_m[1]+_m[2]);
  // check masses non-zero
  for(unsigned int ix=0;ix<2;++ix) {
    if(children[ix]->mass()<1e-4*GeV) { 
      ostringstream message;
      message << "FFDipole::generatePhotons() trying to generate QED radiation from "
	      << children[ix]->dataPtr()->PDGName() << "\n with mass " << children[ix]->mass()/GeV
	      << "which is much smaller than the mass of the electron.\n"
	      << "This is probably due to reading events from a LHEF,\nskipping radiation in this case.\n";
      generator()->logWarning( Exception(message.str(), Exception::warning));
      return children;
    }
  }
  // momenta before radiation in lab
  for(unsigned int ix=0;ix<2;++ix)
    _qlab[ix]=children[ix]->momentum();
  // get the charges of the particles in units of the positron charge
  _charge=children[0]->dataPtr()->iCharge()*children[1]->dataPtr()->iCharge()/9.;
  // boost the momenta to the rest frame
  Boost boostv(-p.momentum().boostVector());
  // boost the particles to the parent rest frame
  // and set the initial momenta of the charged particles 
  // in the dipole rest frame: currently this is the same 
  // as the boson rest frame...
  for(unsigned int ix=0;ix<2;++ix) {
    children[ix]->deepBoost(boostv);
    _qdrf[ix]=children[ix]->momentum();
    _qprf[ix]=children[ix]->momentum();
  }
  _parent->boost(boostv);
  // perform the unweighting
  double wgt;
  unsigned int ntry(0);
  do {
    ++ntry;
    wgt = makePhotons(-boostv,children);

    // Error checks
    if ( std::isnan(wgt) ) {
      generator()->log() << "Infinite weight for decay " 
			 << p.PDGName() << " " 
			 << children[0]->PDGName() 
			 << " " << children[1]->PDGName()
			 << '\n';
      wgt = 0.0;
    }
    else if ( wgt < 0.0 && _mode != 5 ) {
      generator()->log() << "Negative weight for decay " 
			 << p.PDGName() << " " 
			 << children[0]->PDGName() 
			 << " " << children[1]->PDGName()
			 << "in FFDipole: Weight = " << wgt << '\n';
      if ( Debug::level ) 
	printDebugInfo(p,children,wgt);
    }
    else if ( wgt > _maxwgt ) {
      generator()->log() << "Weight "<< wgt<<" exceeds maximum for decay " 
			 << p.PDGName() << ' '
			 << children[0]->PDGName() 
			 << " " << children[1]->PDGName()
			 << " in FFDipole:\nresetting maximum weight.\n"
			 << "Old Maximum = " << _maxwgt;

      _maxwgt = min(1.1 * wgt, 10.0);
      generator()->log() << " New Maximum = " << wgt << '\n';
      if ( Debug::level && _mode!=5 ) 
	printDebugInfo(p,children,wgt);
    }
    // End of error checks

    _wgtsum += wgt;
    _wgtsq  += sqr(wgt);
    ++_nweight;

  }
  while ( wgt<(_maxwgt*UseRandom::rnd()) && ntry<_maxtry );
  if(ntry>=_maxtry) {
    generator()->log() << "FFDipole failed to generate QED radiation for the decay " 
		       << p.PDGName() << " -> " 
		       << children[0]->PDGName() << " "
		       << children[1]->PDGName() << '\n';
    _parent->boost(-boostv);
    for(unsigned int ix=0;ix<2;++ix)
      children[ix]->deepBoost(-boostv);
    return children;
  }
  // produce products after radiation if needed
  if(_multiplicity>0) {
    // change the momenta of the children, they are currently
    // in original rest frame 
    for(unsigned int ix=0;ix<2;++ix) {
      // unit vector along direction
      Boost br = children[ix]->momentum().vect().unit();
      // calculate the boost vector using expression accurate for beta->1
      double beta(sqrt((_qdrf[ix].e()+_m[ix+1])*(_qdrf[ix].e()-_m[ix+1]))/
		  _qdrf[ix].e());
      double ombeta(sqr(_m[ix+1]/_qdrf[ix].e())/(1.+beta));
      double betap(sqrt((_qnewdrf[ix].e()+_m[ix+1])*(_qnewdrf[ix].e()-_m[ix+1]))
		   /_qnewdrf[ix].e());
      double ombetap(sqr(_m[ix+1]/_qnewdrf[ix].e())/(1.+betap));
      // boost to get correct momentum in dipole rest frame
      double bv = -(ombetap-ombeta)/(beta*ombetap + ombeta);
      br *= bv;
      children[ix]->deepBoost(br);
      // boost to the parent rest frame
      Lorentz5Momentum pnew(_bigLdrf);
      pnew.setMass(_m[0]);
      pnew.rescaleEnergy();
      br = pnew.findBoostToCM();
      children[ix]->deepBoost(br);
      // boost back to the lab
      children[ix]->deepBoost(-boostv);
    }
    // add the photons to the event record
    tcPDPtr photon=getParticleData(ParticleID::gamma);
    for(unsigned int ix=0;ix<_multiplicity;++ix) {
      // add if not removed because energy too low
      if(!_photcut[ix]) {
	PPtr newphoton=new_ptr(Particle(photon));
	newphoton->set5Momentum(_llab[ix]);
	children.push_back(newphoton);
      }
    }
    _parent->boost(-boostv);

    //printDebugInfo(p, children, wgt);

    return children;
  }
  // otherwise just return the original particles
  else {
    for(unsigned int ix=0;ix<2;++ix)
      children[ix]->deepBoost(-boostv);
    _parent->boost(-boostv);
    return children;
  }
}

// member which generates the photons
double FFDipole::makePhotons(const Boost & boostv,
			     const ParticleVector & children) {
  // set the initial parameters
  // number of photons (zero)
  _multiplicity=0;
  // zero size of photon vectors
  _ldrf.clear();
  _lprf.clear();
  _llab.clear();
  // zero size of angle storage
  _sinphot.clear();
  _cosphot.clear();
  _photcut.clear();
  _photonwgt.clear();
  // zero total momenta of the photons
  _bigLdrf=Lorentz5Momentum();
  _bigLprf=Lorentz5Momentum();
  // set the initial values of the reweighting factors to one
  _dipolewgt   = 1.0;
  _yfswgt      = 1.0;
  _jacobianwgt = 1.0;
  _mewgt       = 1.0;
  // calculate the velocities of the charged particles (crude/overvalued)
  double beta1(sqrt((_qdrf[0].e()+_m[1])*(_qdrf[0].e()-_m[1]))/_qdrf[0].e());
  double beta2(sqrt((_qdrf[1].e()+_m[2])*(_qdrf[1].e()-_m[2]))/_qdrf[1].e());
  // calculate 1-beta to avoid numerical problems
  double ombeta1(sqr(_m[1]/_qdrf[0].e())/(1.+beta1));
  double ombeta2(sqr(_m[2]/_qdrf[1].e())/(1.+beta2));
  // calculate the average photon multiplicity
  double aver(YFSFormFactors::nbarFF(beta1,ombeta1,beta2,ombeta2,_charge,
				     _emax,_emin,_dipoleopt==1));
  // calculate the number of photons using the poisson
  _multiplicity = _mode !=5 ? UseRandom::rndPoisson(aver) : 1;
  // calculate the first part of the YFS factor
  // (N.B. crude form factor is just exp(-aver) to get a poisson)
  _yfswgt *= exp(aver);
  // if photons produced
  if(_multiplicity>0) {
    _photonwgt.resize(_multiplicity);
    // generate the photon momenta with respect to q1 
    // keeping track of the weight
    for(unsigned int ix=0;ix<_multiplicity;++ix)
      _photonwgt[ix] = photon(beta1,ombeta1,beta2,ombeta2);
    // rotate the photons so in dipole rest frame rather 
    // than angle measured w.r.t q1 first work out the rotation
    SpinOneLorentzRotation rotation;
    rotation.setRotateZ(-_qdrf[0].phi());
    rotation.rotateY(_qdrf[0].theta());
    rotation.rotateZ(_qdrf[0].phi());
    // rotate the total
    _bigLdrf *= rotation;
    // rotate the photons
    for(unsigned int ix=0;ix<_multiplicity;++ix)
      _ldrf[ix]*=rotation;
    // boost the momenta without any removal of low energy photons
    // resize arrays
    _photcut.resize(_multiplicity,false);
    _lprf.resize(_multiplicity);
    _llab.resize(_multiplicity);
    // perform the boost
    if(!boostMomenta(boostv)){return 0.;}
    // apply the cut on the photon energy if needed
    unsigned int nremoved(removePhotons());
    // redo the boost if we have removed photons
    if(nremoved!=0){if(!boostMomenta(boostv)){return 0.;}}
    // form factor part of the removal term to remove existing cut
    if(_energyopt!=0) _dipolewgt *= 
      YFSFormFactors::exponentialYFSFF(beta1,ombeta1,beta2,ombeta2,
				       _qdrf[0].e(),_qdrf[1].e(),
				       _m[1],_m[2],_m[0]*_m[0],
				       _charge,_emin);
    // calculate the new dipole weight
    // calculate velocities and 1-velocites
    beta1=sqrt((_qnewdrf[0].e()+_m[1])*(_qnewdrf[0].e()-_m[1]))/_qnewdrf[0].e();
    beta2=sqrt((_qnewdrf[1].e()+_m[2])*(_qnewdrf[1].e()-_m[2]))/_qnewdrf[1].e();
    ombeta1=sqr(_m[1]/_qnewdrf[0].e())/(1.+beta1);
    ombeta2=sqr(_m[2]/_qnewdrf[1].e())/(1.+beta2);
    for(unsigned int ix=0;ix<_multiplicity;++ix) {
      if(!_photcut[ix])
	_dipolewgt *= exactDipoleWeight(beta1,ombeta1,beta2,ombeta2,ix)/
	  _photonwgt[ix];
    }
    // calculate the weight for the photon removal
    Energy2 s((_qnewdrf[0]+_qnewdrf[1]).m2());
    // calculate the second part of the yfs form factor
    // this is different for the different photon removal options
    // option with no removal
    if(_energyopt==0) {
      _yfswgt *= 
	YFSFormFactors::exponentialYFSFF(beta1,ombeta1,beta2,ombeta2,
					 _qnewdrf[0].e(),_qnewdrf[1].e(),
					 _m[1],_m[2],s,_charge,_emin);
    }
    // weight for option with cut in the rest frame
    else if(_energyopt==1) {
      // yfs piece
      double nbeta1(sqrt( (_qnewprf[0].e()+_m[1])*(_qnewprf[0].e()-_m[1]))
		    /_qnewprf[0].e());
      double nbeta2(sqrt( (_qnewprf[1].e()+_m[2])*(_qnewprf[1].e()-_m[2]))
		    /_qnewprf[1].e());
      double nomb1 (sqr(_m[1]/_qnewprf[0].e())/(1.+nbeta1));
      double nomb2 (sqr(_m[2]/_qnewprf[1].e())/(1.+nbeta2));
      _yfswgt *= 
	YFSFormFactors::exponentialYFSFF(nbeta1,nomb1,nbeta2,nomb2,
					 _qnewprf[0].e(),_qnewprf[1].e(),
					 _m[1],_m[2],s,_charge,_eminrest);
      // dipole piece
      // Find the momenta of the particles of original particles in new rest frame 
      Lorentz5Momentum pnew(_bigLdrf.x(),_bigLdrf.y(),
			    _bigLdrf.z(),_bigLdrf.e(),_m[0]);
      pnew.rescaleEnergy();
      SpinOneLorentzRotation boost(pnew.findBoostToCM());
      Lorentz5Momentum q1=boost*_qdrf[0];
      Lorentz5Momentum q2=boost*_qdrf[1];
      // use this to calculate the form factor
      nbeta1=sqrt( (q1.e()+_m[1])*(q1.e()-_m[1]))/q1.e();
      nbeta2=sqrt( (q2.e()+_m[2])*(q2.e()-_m[2]))/q2.e();
      nomb1 =sqr(_m[1]/q1.e())/(1.+nbeta1);
      nomb2 =sqr(_m[2]/q2.e())/(1.+nbeta2);
      _dipolewgt /=YFSFormFactors::exponentialYFSFF(nbeta1,nomb1,nbeta2,nomb2,
						    q1.e(),q2.e(),
						    _m[1],_m[2],_m[0]*_m[0],
						    _charge,_eminrest);
    }
    // weight for option with cut in the rest frame
    else if(_energyopt==2) {
      // yfs piece
      double nbeta1(sqrt( (_qnewlab[0].e()+_m[1])*(_qnewlab[0].e()-_m[1]))
		    /_qnewlab[0].e());
      double nbeta2(sqrt( (_qnewlab[1].e()+_m[2])*(_qnewlab[1].e()-_m[2]))
		    /_qnewlab[1].e());
      double nomb1 (sqr(_m[1]/_qnewlab[0].e())/(1.+nbeta1));
      double nomb2 (sqr(_m[2]/_qnewlab[1].e())/(1.+nbeta2));
      _yfswgt *= 
	YFSFormFactors::exponentialYFSFF(nbeta1,nomb1,nbeta2,nomb2,
					 _qnewlab[0].e(),_qnewlab[1].e(),
					 _m[1],_m[2],s,_charge,_eminlab);
      // dipole piece
      // Find the momenta of the particles of original particles in new rest frame 
      Lorentz5Momentum pnew(_bigLdrf.x(),_bigLdrf.y(),
			    _bigLdrf.z(),_bigLdrf.e(),_m[0]);
      pnew.rescaleEnergy();
      SpinOneLorentzRotation boost(pnew.findBoostToCM());
      Lorentz5Momentum q1=boost*_qdrf[0];
      Lorentz5Momentum q2=boost*_qdrf[1];
      // then boost to the lab
      boost.setBoost(boostv);
      q1 *=boost;
      q2 *=boost;
      // use this to calculate the form factor
      nbeta1=sqrt( (q1.e()+_m[1])*(q1.e()-_m[1]))
	/q1.e();
      nbeta2=sqrt( (q2.e()+_m[2])*(q2.e()-_m[2]))
	/q2.e();
      nomb1 =sqr(_m[1]/q1.e())/(1.+nbeta1);
      nomb2 =sqr(_m[2]/q2.e())/(1.+nbeta2);
      _dipolewgt /=YFSFormFactors::exponentialYFSFF(nbeta1,nomb1,nbeta2,nomb2,
						    q1.e(),q2.e(),_m[1],_m[2],
						    _m[0]*_m[0],_charge,_eminlab);
    }
    // Calculating jacobian weight
    _jacobianwgt = jacobianWeight();
    // Calculate the weight for the corrections
    _mewgt = meWeight(children);
  }
  // otherwise copy momenta
  else {
    for(unsigned int ix=0;ix<2;++ix) {
      _qnewdrf[ix]=_qdrf[ix];
      _qnewprf[ix]=_qprf[ix]; 
      _qnewlab[ix]=_qlab[ix]; 
    }
    _jacobianwgt = 1.0;
    _yfswgt*=YFSFormFactors::exponentialYFSFF(beta1,ombeta1,beta2,ombeta2,
					      _qdrf[0].e(),_qdrf[1].e(),
					      _m[1],_m[2],_m[0]*_m[0],
					      _charge,_emin);
    _dipolewgt   = 1.0;
  }
  double wgt;
  if(_mode!=5) {
    // virtual corrections
    _mewgt += virtualWeight(children);
    // calculate the weight depending on the option
    if(_mode==0)      wgt = _maxwgt;
    else if(_mode==1) wgt = _mewgt*_yfswgt*_jacobianwgt*_dipolewgt;
    else if(_mode==2) wgt = _jacobianwgt*_yfswgt*_dipolewgt;
    else if(_mode==3) wgt = _yfswgt*_dipolewgt;
    else              wgt = _yfswgt;
  }
  // special to test NLO results
  else {
    double beta1   = sqrt((_qdrf[0].e()+_m[1])*(_qdrf[0].e()-_m[1]))/_qdrf[0].e();
    double beta2   = sqrt((_qdrf[1].e()+_m[2])*(_qdrf[1].e()-_m[2]))/_qdrf[1].e();
    double ombeta1 = sqr(_m[1]/_qdrf[0].e())/(1.+beta1);
    double ombeta2 = sqr(_m[2]/_qdrf[1].e())/(1.+beta2);
    double yfs = YFSFormFactors::YFSFF(beta1,ombeta1,beta2,ombeta2,
				       _qdrf[0].e(),_qdrf[1].e(),
				       _m[1],_m[2],_m[0]*_m[0],
				       _charge,_emin);
    double nbar = YFSFormFactors::nbarFF(beta1,ombeta1,beta2,ombeta2,_charge,
					 _emax,_emin,_dipoleopt==1);
    wgt = 1.+virtualWeight(children)+yfs+nbar*_dipolewgt*_mewgt*_jacobianwgt;
  }
  return wgt;
}

double FFDipole::photon(double beta1,double ombeta1,
                        double beta2,double ombeta2) {
  // generate the polar angle
  double r1,r2,costh,sinth,opbc,ombc;
  // relative weights for the two terms
  double Pp(log((1+beta2)/ombeta2));
  double Pm(log((1+beta1)/ombeta1));
  Pp/=(Pp+Pm);
  // generate the angle
  double wgt=1.;
  do {
    r1=UseRandom::rnd();
    r2=UseRandom::rnd();
    // 1/(1+bc) branch
    if(r1<=Pp) {
      opbc  = pow(1.+beta2,r2)*pow(ombeta2,1.-r2);
      costh = -1./beta2*(1.-opbc);
      ombc  = 1.-beta1*costh;
      sinth = sqrt(opbc*(2.-opbc)-(1.+beta2)*ombeta2*sqr(costh));
    }
    // 1/(1-bc) branch
    else {
      ombc  = pow(1.+beta1,1.-r2)*pow(ombeta1,r2);
      costh = 1./beta1*(1.-ombc);
      opbc  = 1.+beta2*costh;
      sinth = sqrt(ombc*(2.-ombc)-(1.+beta1)*ombeta1*sqr(costh));
    }
    // wgt for rejection
    if(_dipoleopt==1)
      wgt = 1.-0.5/(1.+beta1*beta2)*(ombeta1*(1.+beta1)*opbc/ombc+
				     ombeta2*(1.+beta2)*ombc/opbc);
  }
  while(UseRandom::rnd()>wgt);
  // generate the polar angle randomly in -pi->+pi
  double phi(-pi+UseRandom::rnd()*2.*pi);
  // generate the ln(energy) uniformly in ln(_emin)->ln(_emax)
  Energy en(pow(_emax/_emin,UseRandom::rnd())*_emin);
  // calculate the weight (omit the pre and energy factors
  //                       which would cancel later anyway)
  if(_dipoleopt==0)
    wgt = 0.5*(1.+beta1*beta2)/opbc/ombc;
  else
    wgt = 0.25*(2.*(1.+beta1*beta2)/opbc/ombc
		-ombeta1*(1.+beta1)/sqr(ombc)
		-ombeta2*(1.+beta2)/sqr(opbc));
  // store the angles
  _cosphot.push_back(costh);
  _sinphot.push_back(sinth);
  // store the four vector for the photon
  _ldrf.push_back(Lorentz5Momentum(en*sinth*cos(phi),
				   en*sinth*sin(phi),
				   en*costh,en,
				   ZERO));
  // add the photon momentum to the total
  _bigLdrf+=_ldrf.back();
  // return the weight
  return wgt;
}

double FFDipole::meWeight(const ParticleVector & children) {
  if(_multiplicity==0) return 1.;
  // option which does nothing
  if(_betaopt==0) {
    return 1.;
  }
  // collinear approx
  else if(_betaopt <= 3) {
    return collinearWeight(children);
  }
  else if (_betaopt == 4 ) {
    if(_decayer&&_decayer->hasRealEmissionME()) {
      double outwgt=1.;
      // values of beta etc to evaluate the dipole
      double beta1(sqrt( (_qnewdrf[0].e()+_m[1])*(_qnewdrf[0].e()-_m[1]))/
 		   _qnewdrf[0].e());
      double beta2(sqrt( (_qnewdrf[1].e()+_m[2])*(_qnewdrf[1].e()-_m[2]))/
 		   _qnewdrf[1].e());
      double ombeta1(sqr(_m[1]/_qnewdrf[0].e())/(1.+beta1));
      double ombeta2(sqr(_m[2]/_qnewdrf[1].e())/(1.+beta2));
      // storage of the weights
      ParticleVector ptemp;
      for(unsigned int ix=0;ix<children.size();++ix)
	ptemp.push_back(new_ptr(Particle(children[ix]->dataPtr())));
      ptemp.push_back(new_ptr(Particle(getParticleData(ParticleID::gamma))));
      for(unsigned int i=0;i<_multiplicity;++i) {
	PPtr new_parent = new_ptr(Particle(*_parent));
	if(_photcut[i]) continue;
	// compute the angle terms
	// if cos is greater than zero use result accurate as cos->1
	double opbc,ombc;
	if(_cosphot[i]>0) {
	  opbc=1.+beta2*_cosphot[i];
	  ombc=ombeta1+beta1*sqr(_sinphot[i])/(1.+_cosphot[i]);
	}
	// if cos is less    than zero use result accurate as cos->-1
	else {
	  opbc=ombeta2+beta2*sqr(_sinphot[i])/(1.-_cosphot[i]);
	  ombc=1.-beta1*_cosphot[i];
	}
	// dipole factor for denominator
	double dipole = 2./opbc/ombc*(1.+beta1*beta2
				      -0.5*ombeta1*(1.+beta1)*opbc/ombc		 
				      -0.5*ombeta2*(1.+beta2)*ombc/opbc); 
	// energy and momentum of the photon
	Energy L0(_ldrf[i].e()),modL(_ldrf[i].rho());
	// 3-momenta of charged particles
	Energy modq(_qdrf[0].rho());
	// calculate the energy of the fermion pair
	Energy newE12(-L0+sqrt(sqr(_m[0])+sqr(modL)));
	// 3-momentum rescaling factor (NOT energy rescaling).
	double kappa(Kinematics::pstarTwoBodyDecay(newE12,_m[1],_m[2])/modq);
	// calculate the rescaled momenta
	Lorentz5Momentum porig[3];
	for(unsigned int ix=0;ix<2;++ix) {
	  porig[ix] = kappa*_qdrf[ix];
	  porig[ix].setMass(_m[ix+1]);
	  porig[ix].rescaleEnergy();
	}
	porig[2] = _ldrf[i];
	// calculate the momentum of the decaying particle in dipole rest frame
	Lorentz5Momentum pnew(_ldrf[i].x(),_ldrf[i].y(),
			      _ldrf[i].z(),_ldrf[i].e(),_m[0]);
	pnew.rescaleEnergy();
	// Find the momenta of the particles in the rest frame of the parent...
	// First get the boost from the parent particle
	Boost boost = pnew.findBoostToCM();
	LorentzRotation rot1(-boost, pnew.e()/pnew.mass());
	// check the photon energy
	Lorentz5Momentum ptest = _ldrf[i];
	ptest.boost(boost);
	if(_energyopt==1&&ptest.e()<_eminrest) continue;
	new_parent->transform(rot1);
	// rotation to put the emitter along the z axis
	// first particle emits
	unsigned int iemit = _cosphot[i]>0. ? 0 : 1;
	LorentzRotation rot2;
	rot2.setRotateZ(-porig[iemit].phi());
	rot2.rotateY(porig[iemit].theta());
	rot2.rotateZ(porig[iemit].phi());
	rot2.invert();
	// Boost the momenta of the charged particles
	for(unsigned int ix=0;ix<3;++ix) {
 	  porig[ix].transform(rot2);
	  ptemp[ix]->set5Momentum(porig[ix]);
	}
 	new_parent->transform(rot2);
	if(_cosphot[i]>0.) {
	  outwgt -= _decayer->
	    realEmissionME(_decayer->imode(),*new_parent,ptemp,
			   0,_cosphot[i],_sinphot[i],rot1,rot2)/
	    (_charge/sqr(_ldrf[i].e())*dipole);
	}
	else {
	  outwgt -= _decayer->
	    realEmissionME(_decayer->imode(),*new_parent,ptemp,
			   1,-_cosphot[i],_sinphot[i],rot1,rot2)/
	    (_charge/sqr(_ldrf[i].e())*dipole);
	}
	rot1.invert();
 	rot2.invert();
 	new_parent->transform(rot2);
 	new_parent->transform(rot1);
      }
      return outwgt;
    }
    else
      return collinearWeight(children);
  }
  return 1.;
}

double FFDipole::collinearWeight(const ParticleVector & children) {
  double outwgt=1.;
  // spins of the decay products
  PDT::Spin spin1(children[0]->dataPtr()->iSpin());
  PDT::Spin spin2(children[1]->dataPtr()->iSpin());
  // values of beta etc to evaluate the dipole
  double beta1(sqrt( (_qnewdrf[0].e()+_m[1])*(_qnewdrf[0].e()-_m[1]))/
	       _qnewdrf[0].e());
  double beta2(sqrt( (_qnewdrf[1].e()+_m[2])*(_qnewdrf[1].e()-_m[2]))/
	       _qnewdrf[1].e());
  double ombeta1(sqr(_m[1]/_qnewdrf[0].e())/(1.+beta1));
  double ombeta2(sqr(_m[2]/_qnewdrf[1].e())/(1.+beta2));
  // storage of the weights
  double twgt,dipole;
  double opbc,ombc;
  // compute the collinear approx
  for(unsigned int i=0;i<_multiplicity;++i) {
    if(_photcut[i]) continue;
    // compute the angle terms
    // if cos is greater than zero use result accurate as cos->1
    if(_cosphot[i]>0) {
      opbc=1.+beta2*_cosphot[i];
      ombc=ombeta1+beta1*sqr(_sinphot[i])/(1.+_cosphot[i]);
    }
    // if cos is less    than zero use result accurate as cos->-1
    else {
      opbc=ombeta2+beta2*sqr(_sinphot[i])/(1.-_cosphot[i]);
      ombc=1.-beta1*_cosphot[i];
    }
    // dipole factor for denominator
    dipole = 2.*(1.+beta1*beta2
		 -0.5*ombeta1*(1.+beta1)*opbc/ombc		 
		 -0.5*ombeta2*(1.+beta2)*ombc/opbc); 
    twgt=0.;
    // correction for the first particle
    double ratio(_ldrf[i].e()/_qnewdrf[0].e());
    if(spin1==PDT::Spin0)      twgt += 0.;
    else if(spin1==PDT::Spin1Half) 
      twgt += opbc*ratio/(1.+(1.+beta1*beta2)/ratio/opbc);
    else              
      twgt += 2.*sqr(opbc*ratio) *
	(+1./(1+beta1*beta2+_ldrf[i].e()/_qnewdrf[1].e()*ombc)
	 +(1.+beta1*beta2)/sqr(1.+beta1*beta2
			       +_ldrf[i].e()/_qnewdrf[0].e()*opbc));
    // correction for the second particle
    ratio =_ldrf[i].e()/_qnewdrf[1].e();
    if(spin2==PDT::Spin0)      twgt += 0.;
    else if(spin2==PDT::Spin1Half) 
      twgt += ombc*ratio/(1.+(1.+beta1*beta2)/ratio/ombc);
    else       
      twgt += 2.*sqr(ombc*ratio) *
	(1./(1. + beta1*beta2 + _ldrf[i].e()/_qnewdrf[0].e()*opbc)
	 + (1.+beta1*beta2) / sqr(1. + beta1*beta2
				  + _ldrf[i].e()/_qnewdrf[1].e()*ombc));
    twgt/=dipole;
    outwgt+=twgt;
  }
  return outwgt;
}

bool FFDipole::boostMomenta(const Boost & boostv) {
  // total energy  and momentum of photons
  Energy L0(_bigLdrf.e()),modL(_bigLdrf.rho());
  // 3-momenta of charged particles
  Energy modq(_qdrf[0].rho());
  // calculate the energy of the fermion pair
  Energy newE12(-L0+sqrt(_m[0]*_m[0]+modL*modL));
  // check this is allowed
  if(newE12<_m[1]+_m[2]){return false;}
  // 3-momentum rescaling factor (NOT energy rescaling).
  double kappa(Kinematics::pstarTwoBodyDecay(newE12,_m[1],_m[2])/modq);
  // calculate the rescaled momenta
  for(unsigned int ix=0;ix<2;++ix) {
    _qnewdrf[ix] = kappa*_qdrf[ix];
    _qnewdrf[ix].setMass(_m[ix+1]);
    _qnewdrf[ix].rescaleEnergy();
  }
  // calculate the momentum of the decaying particle in dipole rest frame
  Lorentz5Momentum pnew(_bigLdrf.x(),_bigLdrf.y(),
			_bigLdrf.z(),_bigLdrf.e(),_m[0]);
  pnew.rescaleEnergy();
  // Find the momenta of the particles in the rest frame 
  // of the parent...
  // First get the boost from the parent particle
  SpinOneLorentzRotation boost(pnew.findBoostToCM());
  // Boost the momenta of the charged particles
  for(unsigned int ix=0;ix<2;++ix) _qnewprf[ix]=boost*_qnewdrf[ix];
  // Boost the total photon momentum
  _bigLprf=boost*_bigLdrf;
  // Boost the individual photon momenta
  for(unsigned int ix=0;ix<_multiplicity;++ix){_lprf[ix]=boost*_ldrf[ix];}
  // Now boost from the parent rest frame to the lab frame
  boost.setBoost(boostv);
  // Boosting charged particles
  for(unsigned int ix=0;ix<2;++ix){_qnewlab[ix]=boost*_qnewprf[ix];}
  // Boosting total photon momentum
  _bigLlab=boost*_bigLprf;
  // Boosting individual photon momenta
  for(unsigned int ix=0;ix<_multiplicity;++ix){_llab[ix]=boost*_lprf[ix];}
  return true;
}

unsigned int FFDipole::removePhotons() {
  unsigned int nremoved(0);
  // apply the cut in the rest frame
  if(_energyopt==1) {
    for(unsigned int ix=0;ix<_multiplicity;++ix) {
      if(_lprf[ix].e()<_eminrest) {
	++nremoved;
	_photcut[ix]=true;
	_bigLdrf-=_ldrf[ix];
	_ldrf[ix]=Lorentz5Momentum();
      }
    }
  }
  // apply the cut in the lab frame
  else if(_energyopt==2) {
    for(unsigned int ix=0;ix<_multiplicity;++ix) {
      if(_llab[ix].e()<_eminlab) {
	++nremoved;
	_photcut[ix]=true;
	_bigLdrf-=_ldrf[ix];
	_ldrf[ix]=Lorentz5Momentum();
      }
    }
  }
  // correction factor for dipoles if needed
  if(_dipoleopt==0&&nremoved!=0) {
    // calculate the velocities of the charged particles (crude/overvalued)
    double beta1(sqrt((_qdrf[0].e()+_m[1])*(_qdrf[0].e()-_m[1]))/_qdrf[0].e());
    double beta2(sqrt((_qdrf[1].e()+_m[2])*(_qdrf[1].e()-_m[2]))/_qdrf[1].e());
    // calculate 1-beta to avoid numerical problems
    double ombeta1(sqr(_m[1]/_qdrf[0].e())/(1.+beta1));
    double ombeta2(sqr(_m[2]/_qdrf[1].e())/(1.+beta2));
    // calculate the weights
    for(unsigned int ix=0;ix<_multiplicity;++ix) {
      if(_photcut[ix]) _dipolewgt *= 
	exactDipoleWeight(beta1,ombeta1,beta2,ombeta2,ix)/_photonwgt[ix];
    }
  }
  // return number of remove photons
  return nremoved;
}

double FFDipole::virtualWeight(const ParticleVector & children) {
  double output = 0.;
  // Virtual corrections for beta_0:
  // These should be zero for the scalar case as there is no
  // collinear singularity going by the dipoles above...
  // Use mass of decaying particle...
  if(_betaopt==2) { 
    if((children[0]->dataPtr()->iSpin())==2&&
       (children[1]->dataPtr()->iSpin())==2
       ) {
      output += (1.0*YFSFormFactors::_alpha/pi)
	* log(sqr(_m[0]/_m[1]));
    }
  }
  // OR Use invariant mass of final state children...
  else if(_betaopt==3) { 
    if((children[0]->dataPtr()->iSpin())==2&&
       (children[1]->dataPtr()->iSpin())==2
       ) {
      output += (1.0*YFSFormFactors::_alpha/pi)
	* log((_qnewprf[0]+_qnewprf[1]).m2()/sqr(_m[1]));
    }
  }
  else if (_betaopt==4) {
    if(_decayer&&_decayer->hasOneLoopME()) {
      output += 
	_decayer->oneLoopVirtualME(_decayer->imode(),*_parent,
				   children);
    }
    else {
      output += (1.0*YFSFormFactors::_alpha/pi)
	* log(sqr(_m[0]/_m[1]));
    }
  }
  return output;
}

void FFDipole::dofinish() {
  Interfaced::dofinish();
  if(_weightOutput) {
    _wgtsum /= double(_nweight);
    _wgtsq  /= double(_nweight);
    _wgtsq = max(_wgtsq - sqr(_wgtsum),0.);
    _wgtsq /= double(_nweight);
    _wgtsq = sqrt(_wgtsq);
    generator()->log() << "The average weight for QED Radiation in " << fullName() 
		       << " was " << _wgtsum << " +/- " << _wgtsq << '\n';
  }
}
