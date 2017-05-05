// -*- C++ -*-
//
// IFDipole.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IFDipole class.
//

#include "IFDipole.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"

using namespace ThePEG::Helicity;

using namespace Herwig;

void IFDipole::persistentOutput(PersistentOStream & os) const {
  os << _alpha << ounit(_emin,GeV) << _maxwgt
     << _mode  << _maxtry << _energyopt  << _betaopt;
}

void IFDipole::persistentInput(PersistentIStream & is, int) {
  is >> _alpha >> iunit(_emin,GeV) >> _maxwgt
     >> _mode  >> _maxtry >> _energyopt  >> _betaopt;
}

ClassDescription<IFDipole> IFDipole::initIFDipole;
// Definition of the static class description member.

void IFDipole::Init() {
  static ClassDocumentation<IFDipole> documentation
    ("The IFDipole class implements the initial-final dipole for the SOPTHY algorithm");

  static Switch<IFDipole,unsigned int> interfaceUnWeight
    ("UnWeight",
     "Control the type of unweighting to perform, only one should be used the"
     " other options are for debugging purposes.",
     &IFDipole::_mode, 1, false, false);
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

  static Parameter<IFDipole,unsigned int> interfaceMaximumTries
    ("MaximumTries",
     "Maximum number of attempts to unweight",
     &IFDipole::_maxtry, 500, 10, 100000,
     false, false, Interface::limited);

  static Parameter<IFDipole,Energy> interfaceMinimumEnergyRest
    ("MinimumEnergyRest",
     "The minimum energy of the photons in the rest frame of the decaying particle",
     &IFDipole::_emin, MeV, 1.*MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<IFDipole,double> interfaceMaximumWeight
    ("MaximumWeight",
     "The maximum weight for unweighting",
     &IFDipole::_maxwgt, 2.0, 0.0, 100.0,
     false, false, Interface::limited);

  static Switch<IFDipole,unsigned int> interfaceEnergyCutOff
    ("EnergyCutOff",
     "The type of cut-off on the photon energy to apply",
     &IFDipole::_energyopt, 1, false, false);
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

  static Switch<IFDipole,unsigned int> interfaceBetaOption
    ("BetaOption",
     "Option for the inclusive of the higher beta coefficients",
     &IFDipole::_betaopt, 1, false, false);
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

}

ParticleVector IFDipole::generatePhotons(const Particle & p,ParticleVector children) {
  // set parameters which won't change in the event loop
  // masses of the particles
  _m[0] = p.mass(); 
  _m[1] = children[0]->mass();
  _m[2] = children[1]->mass();
  // momenta before radiation in lab
  for(unsigned int ix=0;ix<2;++ix){_qlab[ix]=children[ix]->momentum();}
  // get the charges of the particles in units of the positron charge
  // chrg1 is the charge of the parent and chrg2 is the charge of the
  // charged child. Also we create a map between the arguments of 
  // _q???[X] _m[X] etc so that 
  // _q???[_map[0]] and _m[_map[0]] are the momenta and masses of 
  // the charged child while
  // _q???[_map[1]] and _m[_map[1]] are the momenta and masses of 
  // the neutral child.
  _chrg1 = p.dataPtr()->iCharge()/3.0;
  if(children[1]->dataPtr()->iCharge()/3.0==0.0) {
    _chrg2  = children[0]->dataPtr()->iCharge()/3.0; 
    _map[0] = 0; _map[1] = 1;
  }
  else if(children[0]->dataPtr()->iCharge()/3.0==0.0) {
    _chrg2  = children[1]->dataPtr()->iCharge()/3.0; 
    _map[0] = 1; _map[1] = 0; 
  }
  // check the radiating particle is not massless
  // if(children[1]->mass()<
  if(children[_map[0]]->mass()<1e-4*GeV) { 
    ostringstream message;
    message << "IFDipole::generatePhotons() trying to generate QED radiation from "
	    << children[_map[0]]->dataPtr()->PDGName() << "\n with mass " << children[_map[0]]->mass()/GeV
	    << "which is much smaller than the mass of the electron.\n"
	    << "This is probably due to reading events from a LHEF,\nskipping radiation in this case.\n";
    generator()->logWarning( Exception(message.str(), Exception::warning));
    return children;
  } 
  // boost the momenta to the rest frame
  Boost boostv(p.momentum().boostVector());
  // boost the particles to the parent rest frame
  // and set the initial momenta of the charged particles 
  // in the dipole rest frame: currently this is the same 
  // as the boson rest frame...
  for(unsigned int ix=0;ix<2;++ix) {
    // KMH - 08/11/05 - This used to be boostv instead of -boostv
    // -boostv is the boost from the lab to the parent rest frame
    // whereas boostv goes the other way!!!
    children[ix]->deepBoost(-boostv);
    _qprf[ix]=children[ix]->momentum();
  }
  // perform the unweighting
  double wgt;
  unsigned int ntry(0);
  do {
    wgt =makePhotons(boostv,children);
    ++ntry;
    // Record warnings about large and weird weights in the .log file.
    if(wgt>_maxwgt||wgt<0.0||std::isnan(wgt)) {
      generator()->log() << "IFDipole.cc:\n";
      if(wgt>_maxwgt) {
	generator()->log() << "Weight exceeds maximum for decay!\n"; 
      } 
      if(wgt<0.0) {
	generator()->log() << "Weight is negative! \n"; 
      }
      if(std::isnan(wgt)) {
	generator()->log() << "Weight is NAN! \n";
	wgt = 0.;
      }
      generator()->log() << p.PDGName() << " " 
			 << children[0]->PDGName() << " " 
			 << children[1]->PDGName()
			 << endl 
			 << " Current Maximum = " << _maxwgt 
			 << endl
			 << " Current Weight  = " << wgt  
			 << endl;
      generator()->log() << "Photon Multiplicity      : " 
			 << _multiplicity                          << endl
			 << "Original Parent rest frame momenta: " << endl
			 << "charged child: " << ounit(_qprf[_map[0]],GeV) << endl
			 << "neutral child: " << ounit(_qprf[_map[1]],GeV) << endl
			 << "Parent rest frame momenta: "          << endl
			 << "charged child: " << ounit(_qnewprf[_map[0]],GeV)<< endl
			 << "neutral child: " << ounit(_qnewprf[_map[1]],GeV)<< endl
			 << "photons      : " << ounit(_bigLprf,GeV) << endl
			 << "Weights      : "                      << endl
			 << "_dipolewgt   : " << _dipolewgt        << endl
			 << "_yfswgt      : " << _yfswgt           << endl
			 << "_jacobianwgt : " << _jacobianwgt      << endl
			 << "_mewgt       : " << _mewgt            << endl;
      for(unsigned int ct=0;ct<_multiplicity;ct++) {
	generator()->log() << "_cosphot[" << ct << "]: " << _cosphot[ct] << endl;
	generator()->log() << "_sinphot[" << ct << "]: " << _sinphot[ct] << endl;
      }
      if(wgt>_maxwgt) {
	if(wgt<15.0) { 
	  generator()->log() << "Resetting maximum weight" 
			     << endl << " New Maximum = " << wgt  << endl;
	  _maxwgt=wgt; 
	} else {
	  generator()->log() << "Maximum weight set to limit (15)" << endl;
	  _maxwgt=15.0; 
	}
      }
    }
  } while (wgt<(_maxwgt*UseRandom::rnd()) && ntry<_maxtry);
  if(ntry>=_maxtry) {
    generator()->log() << "IFDipole Failed to generate QED radiation for the decay " 
		       << p.PDGName() << " -> " 
		       << children[0]->PDGName() << " "
		       << children[1]->PDGName() << endl;
    return children;
  }
  // produce products after radiation if needed
  if(_multiplicity>0) {
    // change the momenta of the children, they are currently
    // in parent rest frame
    for(unsigned int ix=0;ix<2;++ix) {
      LorentzRotation boost(solveBoost(_qnewprf[ix],children[ix]->momentum()));
      children[ix]->deepTransform(boost);
      // boost back to the lab
      // KMH - 08/11/05 - This used to be -boostv instead of boostv
      // -boostv is the boost from the lab to the parent rest frame
      // whereas boostv goes the other way!!!
      children[ix]->deepBoost(boostv);
    }
    // add the photons to the event record
    tcPDPtr photon=getParticleData(ParticleID::gamma);
    for(unsigned int ix=0;ix<_multiplicity;++ix) {
      PPtr newphoton=new_ptr(Particle(photon));
      newphoton->set5Momentum(_llab[ix]);
      children.push_back(newphoton);
    }
    return children;
  }
  // otherwise just return the orginial particles
  // boosted back to lab
  else {
    for(unsigned int ix=0;ix<children.size();++ix)
      children[ix]->deepBoost(boostv);
    return children;
  }
}

// member which generates the photons
double IFDipole::makePhotons(Boost boostv,ParticleVector children) {
  // set the initial parameters
  // number of photons (zero)
  _multiplicity=0;
  // zero size of photon vectors
  _lprf.clear();
  _llab.clear();
  // zero size of angle storage
  _sinphot.clear();
  _cosphot.clear();
  // zero total momenta of the photons
  _bigLprf=Lorentz5Momentum();
  // set the initial values of the reweighting factors to one
  _dipolewgt   = 1.0;
  _yfswgt      = 1.0;
  _jacobianwgt = 1.0;
  _mewgt       = 1.0;
  // set the maximum photon energy (exact - no approximations here).
  double boost_factor = 1.0;
  _emax=(0.5*(_m[0]-sqr(_m[1]+_m[2])/_m[0]))*boost_factor;
  // calculate the velocities of the children (crude/overvalued)
  double beta1(sqrt( (_qprf[_map[0]].e()+_m[_map[0]+1])
                    *(_qprf[_map[0]].e()-_m[_map[0]+1])
                   )
                   /_qprf[_map[0]].e());
  double beta2(sqrt( (_qprf[_map[1]].e()+_m[_map[1]+1]) 
                    *(_qprf[_map[1]].e()-_m[_map[1]+1])
                   )
                   /_qprf[_map[1]].e());
  // calculate 1-beta to avoid numerical problems
  double ombeta1(sqr(_m[_map[0]+1]/_qprf[_map[0]].e())/(1.+beta1));
  double ombeta2(sqr(_m[_map[1]+1]/_qprf[_map[1]].e())/(1.+beta2));
  // calculate the average photon multiplicity
  double aver(nbar(beta1,ombeta1));
  // calculate the number of photons using the poisson
  _multiplicity = UseRandom::rndPoisson(aver);
  // calculate the first part of the YFS factor
  _yfswgt/=crudeYFSFormFactor(beta1,ombeta1); 
  // generate the photon momenta with respect to q1 
  // keeping track of the weight
  double dipoles(1.);
  for(unsigned int ix=0;ix<_multiplicity;++ix)
  { dipoles *= photon(beta1,ombeta1); }
  // calculate contributions to the dipole weights so far
  _dipolewgt /=dipoles;

  // now do the momentum reshuffling
  Lorentz5Momentum pmom(ZERO,ZERO,ZERO,_m[0],_m[0]);
  if(_multiplicity>0) {
      // total energy  and momentum of photons
      Energy L0(_bigLprf.e()),modL(_bigLprf.rho());
      // squared invariant mass of final state fermions...
      Energy2 m122 = sqr(_m[0]-L0)-sqr(modL);
      if(m122<sqr(_m[1]+_m[2])) return 0.;
      // 3-momenta of charged particles
      Energy modq(_qprf[_map[0]].rho());
      // total photon momentum perpendicular to charged child...
      Energy LT(_bigLprf.perp());
      // kallen function...
      Energy4 kallen = ( m122 - sqr(_m[1]+_m[2]) )
	             * ( m122 - sqr(_m[1]-_m[2]) );
      // discriminant of rho...
      Energy4 droot = kallen-4.*sqr(_m[_map[0]+1]*LT);
      if(droot<ZERO) return 0.;
      double disc = (_m[0]-L0) *  sqrt(droot) / (2.*modq*(m122+LT*LT));
      // calculate the energy rescaling factor
      double rho  = disc-_bigLprf.z()
	          * (m122+sqr(_m[_map[0]+1])-sqr(_m[_map[1]+1]))
  	          / (2.*modq*(m122+LT*LT));
      // calculate the rescaled charged child momentum
      _qnewprf[_map[0]]=rho*_qprf[_map[0]];
      _qnewprf[_map[0]].setMass(_m[_map[0]+1]);
      _qnewprf[_map[0]].rescaleEnergy();
      // rotate the photons so in parent rest frame rather 
      // than angle measured w.r.t q1 first work out the rotation
      SpinOneLorentzRotation rotation;
      rotation.setRotateZ(-_qprf[_map[0]].phi());
      rotation.rotateY(_qprf[_map[0]].theta());
      rotation.rotateZ(_qprf[_map[0]].phi());
      // rotate the total
      _bigLprf*=rotation;
      // rotate the photons
      for(unsigned int ix=0;ix<_multiplicity;++ix){_lprf[ix]*=rotation;}
      // calculate the rescaled neutral child momentum
      _qnewprf[_map[1]]=pmom-_qnewprf[_map[0]]-_bigLprf;
      _qnewprf[_map[1]].setMass(_m[_map[1]+1]);
      _qnewprf[_map[1]].rescaleEnergy();
      // calculate the new dipole weight
      // Note this (weight) is Lorentz invariant
      // calculate velocities and 1-velocites
      beta1=sqrt( (_qnewprf[_map[0]].e()+_m[_map[0]+1])
                 *(_qnewprf[_map[0]].e()-_m[_map[0]+1]))
                /_qnewprf[_map[0]].e();
      beta2=sqrt( (_qnewprf[_map[1]].e()+_m[_map[1]+1])
                 *(_qnewprf[_map[1]].e()-_m[_map[1]+1]))
                /_qnewprf[_map[1]].e();
      ombeta1=sqr(_m[_map[0]+1]/_qnewprf[_map[0]].e())/(1.+beta1);
      ombeta2=sqr(_m[_map[1]+1]/_qnewprf[_map[1]].e())/(1.+beta2);
      for(unsigned int ix=0;ix<_multiplicity;++ix)
	{_dipolewgt*=exactDipoleWeight(beta1,ombeta1,ix);}
      // calculate the second part of the yfs form factor
      _yfswgt*=exactYFSFormFactor(beta1,ombeta1,beta2,ombeta2);
      // Now boost from the parent rest frame to the lab frame
      SpinOneLorentzRotation boost(boostv);
      // Boosting charged particles
      for(unsigned int ix=0;ix<2;++ix){_qnewlab[ix]=boost*_qnewprf[ix];}
      // Boosting total photon momentum
      _bigLlab=boost*_bigLprf;
      // Boosting individual photon momenta
      for(unsigned int ix=0;ix<_multiplicity;++ix)
        {_llab.push_back(boost*_lprf[ix]);}
      // Calculating jacobian weight
      _jacobianwgt = jacobianWeight();
      // Calculating beta^1  weight
      _mewgt = meWeight(children);
      // Apply phase space vetos...
      if(kallen<(4.*sqr(_m[_map[0]+1]*LT))||m122<sqr(_m[1]+_m[2])||rho<0.0) {  
//           generator()->log() << "Outside Phase Space" << endl;
//           generator()->log() << "Photon Multiplicity: " 
//                              << _multiplicity                          << endl
//                              << "Original Parent rest frame momenta: " << endl
//                              << "charged child: " << _qprf[_map[0]]    << endl
//                              << "neutral child: " << _qprf[_map[1]]    << endl
//                              << "rescaling    : " << rho               << endl
//                              << "Parent rest frame momenta: "          << endl
//                              << "charged child: " << _qnewprf[_map[0]] << endl
//                              << "neutral child: " << _qnewprf[_map[1]] << endl
//                              << "photons      : " << _bigLprf          << endl
//                              << endl;
	_dipolewgt   = 0.0 ;
	_yfswgt      = 0.0 ;
	_jacobianwgt = 0.0 ;
	_mewgt       = 0.0 ;
      }
      _qprf[_map[0]].rescaleEnergy();
      _qprf[_map[1]].rescaleEnergy();
      _qnewprf[_map[0]].rescaleEnergy();
      _qnewprf[_map[1]].rescaleEnergy();
      if( ((abs(_m[0]-_bigLprf.e()-_qnewprf[0].e()-_qnewprf[1].e())>0.00001*MeV)||
           (abs(      _bigLprf.x()+_qnewprf[0].x()+_qnewprf[1].x())>0.00001*MeV)||
           (abs(      _bigLprf.y()+_qnewprf[0].y()+_qnewprf[1].y())>0.00001*MeV)||
           (abs(      _bigLprf.z()+_qnewprf[0].z()+_qnewprf[1].z())>0.00001*MeV))
         &&(_dipolewgt*_jacobianwgt*_yfswgt*_mewgt>0.0)) {
	Lorentz5Momentum ptotal = _bigLprf+_qnewprf[0]+_qnewprf[1];
	ptotal.setE(ptotal.e()-_m[0]);
	generator()->log() 
	  <<   "Warning! Energy Not Conserved! tol = 0.00001 MeV"
	  << "\nwgt               = " << _dipolewgt*_yfswgt*_jacobianwgt*_mewgt
	  << "\nrho               = " << rho
	  << "\nmultiplicity      = " << _multiplicity
	  << "\n_qprf[_map[0]]    = " << _qprf[_map[0]]/GeV
	  << "\n_qprf[_map[1]]    = " << _qprf[_map[1]]/GeV
	  << "\n_qnewprf[_map[0]] = " << _qnewprf[_map[0]]/GeV << " " 
	  << _qnewprf[_map[0]].m()/GeV << " " << _m[_map[0]+1]/GeV
	  << "\n_qnewprf[_map[1]] = " << _qnewprf[_map[1]]/GeV << " " 
	  << _qnewprf[_map[1]].m()/GeV << " " << _m[_map[1]+1]/GeV
	  << "\n_bigLprf          = " << _bigLprf/GeV
	  << "\n_bigLprf.m2()     = " << _bigLprf.m2()/GeV2
	  << "\n_total out -in    = " << ptotal/GeV
	  << "\nRejecting Event.    " << "\n";
	_dipolewgt   = 0.0 ;
	_yfswgt      = 0.0 ;
	_jacobianwgt = 0.0 ;
	_mewgt       = 0.0 ;
      }
    }
  // otherwise copy momenta
  else
    { for(unsigned int ix=0;ix<2;++ix) {
        _qnewprf[ix]=_qprf[ix];
        _qnewlab[ix]=_qlab[ix]; 
      }
      _jacobianwgt = 1.0;
      // calculate the second part of the yfs form factor
      _yfswgt*=exactYFSFormFactor(beta1,ombeta1,beta2,ombeta2);
      _dipolewgt   = 1.0;
    }
// Virtual corrections for beta_0:
// These should be zero for the scalar case as there is no
// collinear singularity going by the dipoles above...
// Use mass of decaying particle...
  if(_betaopt==2) { 
    if((children[_map[0]]->dataPtr()->iSpin())==2) {
      _mewgt += (0.5*_alpha/pi) 
              * log(sqr(_m[0]
                       /_m[_map[0]+1])
                   );
    }
  }
// OR Use invariant mass of final state children...
  if(_betaopt==3) { 
    if((children[_map[0]]->dataPtr()->iSpin())==2) {
      _mewgt += (0.5*_alpha/pi)
              * log((_qnewprf[0]+_qnewprf[1]).m2()
                   /sqr(_m[_map[0]+1])
                   );
    }
  }
  // calculate the weight depending on the option
  double wgt;
  if(_mode==0){wgt=_maxwgt;}
  else if(_mode==1){wgt=_mewgt*_jacobianwgt*_yfswgt*_dipolewgt;}
  else if(_mode==2){wgt=_jacobianwgt*_yfswgt*_dipolewgt;}
  else if(_mode==3){wgt=_yfswgt*_dipolewgt;}
  else             {wgt=_yfswgt;}
  return wgt;
}

double IFDipole::photon(double beta1,double ombeta1)
{
  // generate the azimuthal angle randomly in -pi->+pi
  double phi(-pi+UseRandom::rnd()*2.*pi);
  // generate the polar angle
  double r(UseRandom::rnd());
  double costh,sinth,ombc;
  ombc  = pow(1.+beta1,1.-r)*pow(ombeta1,r);
  costh = 1./beta1*(1.-ombc);
  sinth = sqrt(ombc*(2.-ombc)-(1.+beta1)*ombeta1*sqr(costh));
  // generate the ln(energy) uniformly in ln(_emin)->ln(_emax)
  Energy energy   = pow(_emax/_emin,UseRandom::rnd())*_emin;
  // calculate the weight (omit the pre and energy factors
  //                       which would cancel later anyway)
  double wgt = 2./ombc;
  // store the angles
  _cosphot.push_back(costh);
  _sinphot.push_back(sinth);
  // store the four vector for the photon
  _lprf.push_back(Lorentz5Momentum(energy*sinth*cos(phi),energy*sinth*sin(phi),
				   energy*costh,energy,ZERO));
  // add the photon momentum to the total
  _bigLprf+=_lprf.back();
  // return the weight
  return wgt;
}

double IFDipole::meWeight(ParticleVector children)
{
  unsigned int spin = children[_map[0]]->dataPtr()->iSpin();
  double mewgt = 1.0;
  double beta1=sqrt( (_qnewprf[_map[0]].e()+_m[_map[0]+1])
                    *(_qnewprf[_map[0]].e()-_m[_map[0]+1]))
                   /_qnewprf[_map[0]].e();
  double ombeta1=sqr(_m[_map[0]+1]/_qnewprf[_map[0]].e())/(1.+beta1);
  // option which does nothing
  if(_betaopt==0){mewgt=1.;}
  // collinear approx
  else if(_betaopt==1||_betaopt==2||_betaopt==3)
    {
      double ombc;
      InvEnergy2 dipole;
      for(unsigned int i=0;i<_multiplicity;++i) {
	double opbc;
        if(_cosphot[i]<0.0)
          { opbc=ombeta1+beta1*sqr(_sinphot[i])/(1.-_cosphot[i]); }
        // if cos is greater than zero use result accurate as cos->-1
        else
          { opbc=1.+beta1*_cosphot[i]; }
        // if cos is greater than zero use result accurate as cos->1
        if(_cosphot[i]>0.0)
          { ombc=ombeta1+beta1*sqr(_sinphot[i])/(1.+_cosphot[i]); }
        // if cos is less    than zero use result accurate as cos->-1
        else
          { ombc=1.-beta1*_cosphot[i]; }
        if(((_qnewprf[_map[0]].z()>ZERO)&&(_qprf[_map[0]].z()<ZERO))||
           ((_qnewprf[_map[0]].z()<ZERO)&&(_qprf[_map[0]].z()>ZERO))) {
          dipole = sqr(beta1*_sinphot[i]/(opbc*_lprf[i].e()));
        } else {
          dipole = sqr(beta1*_sinphot[i]/(ombc*_lprf[i].e()));
	}
	// here "dipole" is the exact dipole function divided by alpha/4pi^2.

        if(spin==2) {
	  Energy magpi= sqrt( sqr(_qnewprf[_map[0]].x())
			      + sqr(_qnewprf[_map[0]].y())
			      + sqr(_qnewprf[_map[0]].z())
			      );

	  mewgt += sqr(_lprf[i].e())*_qnewprf[_map[0]].e()*ombc
	         / (sqr(magpi*_sinphot[i])*(_qnewprf[_map[0]].e()+_lprf[i].e()));
        }
        else if(spin==3) {
	  Energy2 pik  = _qnewprf[_map[0]].e()*_lprf[i].e()
	               - _qnewprf[_map[0]].x()*_lprf[i].x()
                       - _qnewprf[_map[0]].y()*_lprf[i].y()
	               - _qnewprf[_map[0]].z()*_lprf[i].z();

	  Energy2 pjk = _m[0]*_lprf[i].e();

	  Energy2 pipj = _m[0]*_qnewprf[_map[0]].e();

          mewgt += (2.*pjk*pipj/(pik*sqr(pipj+pjk))
		   +2.*pjk/(pik*(pipj+pik))
	           )/dipole;
        }
        else {
          mewgt = 1.0;
        }
      }
    }
  return mewgt;
}

double IFDipole::exactYFSFormFactor(double beta1,double ombeta1,
					   double beta2,double ombeta2) {
  double Y    = 0.0    ;
  double b    = beta1  ;
  double omb  = ombeta1;
  double c    = beta2  ;
  double omc  = ombeta2;
  double arg1 = -omc/(2.*c);
  double arg2 = -omb*omc/(2.*(b+c));
  double arg3 = 2.*b/(1.+b);
  if(_m[_map[1]+1]!=ZERO) {
    Y = _chrg1*_chrg2*(_alpha/(2.*pi))*(
         log(_m[0]*_m[_map[1]+1]/sqr(2.*_emin))
        +log(_m[_map[0]+1]*_m[_map[1]+1]/sqr(2.*_emin))
        -(1./b )*log((1.+b)/omb)*log(sqr(_m[_map[1]+1]/(2.*_emin)))
        -(1./b )*log(omb/(1.+b))
        -(0.5/b )*sqr(log(omb/(1.+b)))
        +((b+c  )/(b*omc))*log((b+c  )/(b*omc))
        -((c+b*c)/(b*omc))*log((c+b*c)/(b*omc))
        +((b+c  )/(b+b*c))*log((b+c  )/(b+b*c))
        -((c*omb)/(b+b*c))*log((c*omb)/(b+b*c))
        +(0.5/b)*(   sqr(log(  (b+c)/(b*omc)))-sqr(log((c+b*c)/(b*omc)))
       	           + sqr(log((c*omb)/(b+b*c)))-sqr(log((b+  c)/(b+b*c)))
	         )  
        +(2./b )*(   real(Math::Li2(arg1))
                   - real(Math::Li2(arg2))
                   - real(Math::Li2(arg3))
  	         )
        +(1./b )*log((b+c)/(b+b*c))*log((1.+c)/(2.*c))
        -(1./b )*log((c*omb)/(b*(1.+c)))*log((1.+b)*(1.+c)/(2.*(b+c)))
        -(1./b )*log((2.*c/b)*((b+c)/(omc*(1.+c))))*log((b+c)/(c*omb))
        );
  }
  else if(_m[_map[1]+1]==ZERO) {
    Y = _chrg1*_chrg2*(_alpha/(2.*pi))*(
         log(sqr(_m[0]/(2.*_emin)))
        +log(sqr(_m[_map[0]+1]/(2.*_emin)))
        -(1./b )*log((1.+b)/omb)
                *log((sqr(_m[0])-sqr(_m[_map[0]+1]))/sqr(2.*_emin))
        -0.5*log(omb*(1.+b)/sqr(2.*b))
        +((1.+b)/(2.*b))*log((1.+b)/(2.*b))
        -(   omb/(2.*b))*log(   omb/(2.*b))
        -(1./b )*log((1.-b)/(1.+b))
        +1.
        +(0.5/b)*sqr(log(   omb/(2.*b)))
        -(0.5/b)*sqr(log((1.+b)/(2.*b)))  
        -(0.5/b)*sqr(log((1.-b)/(1.+b)))
        -(2. /b)*real(Math::Li2(arg3))
        );
  }
  return exp(Y);
}

double IFDipole::jacobianWeight() {
  // calculate the velocities of the children (crude/overvalued)
  Energy mag1old  = sqrt( (_qprf[_map[0]].e()   +_m[_map[0]+1])
                         *(_qprf[_map[0]].e()   -_m[_map[0]+1])
                        );
  Energy mag1new  = sqrt( (_qnewprf[_map[0]].e()+_m[_map[0]+1])
                         *(_qnewprf[_map[0]].e()-_m[_map[0]+1])
			);
  Energy magL     = sqrt( sqr(_bigLprf.x())
			  + sqr(_bigLprf.y())
			  + sqr(_bigLprf.z())
			  );

// 14/12/05 - KMH - This was another mistake. This is supposed to be 
// the angel between _qnewprf[_map[0]] and _bigLprf instead of
// between _qnewprf[0] and _bigLprf. Stupid. Hopefully this weight
// is correct now.
//  double cos1L    = (_qnewprf[0].x()*_bigLprf.x()
//                    +_qnewprf[0].y()*_bigLprf.y()
//                    +_qnewprf[0].z()*_bigLprf.z()
//                    )
//                    /(mag1new*magL);
  double cos1L    = (_qnewprf[_map[0]].x()*_bigLprf.x()
                    +_qnewprf[_map[0]].y()*_bigLprf.y()
                    +_qnewprf[_map[0]].z()*_bigLprf.z()
                    )
                    /(mag1new*magL);

  return abs(  (_m[0]*sqr(mag1new)/mag1old)
             / (  mag1new*(_m[0]-_bigLprf.e())
                +_qnewprf[_map[0]].e()*magL*cos1L
	       )
            );
}

LorentzRotation IFDipole::solveBoost(const Lorentz5Momentum & q, 
				     const Lorentz5Momentum & p ) const {
  Energy modp = p.vect().mag();
  Energy modq = q.vect().mag();
  double betam = (p.e()*modp-q.e()*modq)/(sqr(modq)+sqr(modp)+p.mass2());
  Boost beta = -betam*q.vect().unit();
  ThreeVector<Energy2> ax = p.vect().cross( q.vect() ); 
  double delta = p.vect().angle( q.vect() );
  LorentzRotation R;
  using Constants::pi;
  if ( ax.mag2()/GeV2/MeV2 > 1e-16 ) {
    R.rotate( delta, unitVector(ax) ).boost( beta );
  } 
  else {
    if(p.mass()>ZERO) {
      R.boost(p.findBoostToCM(),p.e()/p.mass());
      R.boost(q.boostVector(),q.e()/q.mass());
    }
    else {
      if(modp>modq) beta = -betam*p.vect().unit();
      R.boost( beta );
    }
  } 
  return R;
}

void IFDipole::doinit() {
  Interfaced::doinit();
  // get the value fo alpha from the Standard Model object
  _alpha=generator()->standardModel()->alphaEM();
}
