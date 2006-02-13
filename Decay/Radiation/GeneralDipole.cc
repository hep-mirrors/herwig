// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GeneralDipole class.
//

#include "GeneralDipole.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "GeneralDipole.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

GeneralDipole::~GeneralDipole() {}

void GeneralDipole::persistentOutput(PersistentOStream & os) const {
  os << _emin << _eminrest << _eminlab << _maxwgt << _maxtry << _nphotonmax 
     << _mode << _energyopt << _betaopt;
}

void GeneralDipole::persistentInput(PersistentIStream & is, int) {
  is >> _emin >> _eminrest >> _eminlab >> _maxwgt >> _maxtry >> _nphotonmax 
     >> _mode >> _energyopt >> _betaopt;
}

ClassDescription<GeneralDipole> GeneralDipole::initGeneralDipole;
// Definition of the static class description member.

void GeneralDipole::Init() {

  static ClassDocumentation<GeneralDipole> documentation
    ("There is no documentation for the GeneralDipole class");


  static Parameter<GeneralDipole,Energy> interfaceMinimumEnergyBoosted
    ("MinimumEnergyBoosted",
     "The minimum energy of the photons in the boosted frame in which"
     " they are generated.",
     &GeneralDipole::_emin, MeV, 1.e-6*MeV, 0.0*MeV, 100.0*MeV,
     false, false, Interface::limited);

  static Parameter<GeneralDipole,Energy> interfaceMinimumEnergyRest
    ("MinimumEnergyRest",
     "The minimum energy of the photons in the rest frame of the decaying particle",
     &GeneralDipole::_eminrest, MeV, 100.0*MeV, 1.0*MeV, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<GeneralDipole,Energy> interfaceMinimumEnergyLab
    ("MinimumEnergyLab",
     "The minimum energy of the photons in the lab frame",
     &GeneralDipole::_eminlab, MeV, 100.0*MeV, 1.0*MeV, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<GeneralDipole,double> interfaceMaximumWeight
    ("MaximumWeight",
     "The maximum weight for unweighting",
     &GeneralDipole::_maxwgt, 2.0, 0.0, 100.0,
     false, false, Interface::limited);

  static Parameter<GeneralDipole,unsigned int> interfaceMaximumTries
    ("MaximumTries",
     "Maximum number of attempts to unweight",
     &GeneralDipole::_maxtry, 500, 10, 100000,
     false, false, Interface::limited);

  static Parameter<GeneralDipole,unsigned int> interfaceMaximumNumberOfPhotons
    ("MaximumNumberOfPhotons",
     "The maximum number of photons to produce",
     &GeneralDipole::_nphotonmax, 20, 1, 1000,
     false, false, Interface::limited);

  static Switch<GeneralDipole,unsigned int> interfaceUnWeight
    ("UnWeight",
     "Control the type of unweighting to perform, only one should be used the"
     " other options are for debugging purposes.",
     &GeneralDipole::_mode, 1, false, false);
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

  static Switch<GeneralDipole,unsigned int> interfaceEnergyCutOff
    ("EnergyCutOff",
     "The type of cut-off on the photon energy to apply",
     &GeneralDipole::_energyopt, 1, false, false);
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

  static Switch<GeneralDipole,unsigned int> interfaceBetaOption
    ("BetaOption",
     "Option for the inclusive of the higher beta coefficients",
     &GeneralDipole::_betaopt, 1, false, false);
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
  static SwitchOption interfaceBetaOptionExact
    (interfaceBetaOption,
     "Exact",
     "Include the exact higher order terms if available",
     2);
}

ParticleVector GeneralDipole::generatePhotons(const Particle & p,ParticleVector children)
{
  // obtain a number of physical properties of the particles
  // boost from the rest frame to the lab
  _boosttolab=p.momentum().boostVector();
  // number of decay products
  _nprod=children.size();
  // resize the vectors containing the masses of the particles
  _m.resize(_nprod+1,0.);_m2.resize(_nprod+1,0.);
  // vectors containing the velocities
  _beta.resize(   _nprod+1,0.);_ombeta.resize(   _nprod+1,0.);
  _betanew.resize(_nprod+1,0.);_ombetanew.resize(_nprod+1,0.);
  // resize the vectors containing the momenta
  _qlab.resize(_nprod+1);_qnewlab.resize(_nprod+1);
  _qprf.resize(_nprod+1);_qnewprf.resize(_nprod+1);
  _qdrf.resize(_nprod+1);_qnewdrf.resize(_nprod+1);
  // properties of the decaying particle
  // mass
  _m[0]=p.mass();_m2[0]=_m[0]*_m[0];
  // orginial momentum
  _qlab[0]=p.momentum();
  _qprf[0]=p.momentum();_qprf[0].boost(-_boosttolab);
  _qdrf[0]=_qprf[0];
  // properties of the decay products
  _msum=0.;
  for(unsigned int ix=0;ix<_nprod;++ix)
    {
      // masses of the decay products
      _m[ix+1]=children[ix]->mass();
      _m2[ix+1]=_m[ix+1]*_m[ix+1];
      // sum of masses
      _msum+=_m[ix+1];
      // momenta before radiation in the lab
      _qlab[ix+1]=children[ix]->momentum();
      // boost to the rest frame
      children[ix]->deepBoost(-_boosttolab);
      _qdrf[ix+1]=children[ix]->momentum();
      _qprf[ix+1]=children[ix]->momentum();
      // velocities and 1-velocity
      _beta[ix+1]   = sqrt((_qdrf[ix+1].e()+_m[ix+1])*(_qdrf[ix+1].e()-_m[ix+1]))/
	_qdrf[ix+1].e();
      _ombeta[ix+1] = sqr(_m[ix+1]/_qdrf[ix+1].e())/(1.+_beta[ix+1]);
    }
  // maximum photon energy
  _emax = 0.5*_m[0]*(_m[0]/_msum-_msum/_m[0]);
  // properties of the dipoles
  // resize the vector containing the charges
  _zij.resize(_nprod+1,vector<int>(_nprod+1,0));
  // resize the vector containing the average multiplicities
  _nbar.resize(_nprod+1,vector<double>(_nprod+1,0.));
  // dipole opening angles
  _cosij.resize(_nprod+1,vector<double>(_nprod+1,0.));
  _sinij.resize(_nprod+1,vector<double>(_nprod+1,0.));
  // dipole masses
  _mdipole.resize(_nprod+1,vector<double>(_nprod+1,0.));
  _mdipolenew.resize(_nprod+1,vector<double>(_nprod+1,0.));
  // charge of the decaying particle
  int Zin(p.dataPtr()->iCharge()),Ztemp;
  for(unsigned int ix=1;ix<=_nprod;++ix)
    {
      // charges of initial-final dipoles
      Ztemp = children[ix-1]->dataPtr()->iCharge();
      _zij[ 0][ix]=-Zin*Ztemp;
      _zij[ix][ 0]=_zij[0][ix];
      // angle of the inital-final dipoles
      _cosij[ 0][ix] = _qdrf[ix].cosTheta();
      _cosij[ix][ 0] = _cosij[ 0][ix];
      _sinij[ 0][ix] = sqrt(1.-sqr(_cosij[ 0][ix]));
      _sinij[ix][ 0] = _sinij[ 0][ix];
      // mass of the initial-final dipole
      _mdipole[0 ][ix] = (_qdrf[0]-_qdrf[ix]).m2();
      _mdipole[ix][ 0] = _mdipole[0 ][ix];
      // crude photon multiplcity for initial-final dipoles
      _nbar[ 0][ix] = nbar(0,ix);
      _nbar[ix][0 ] = _nbar[0][ix];
      // properties of the final-final dipoles
      for(unsigned int iy=ix+1;iy<=_nprod;++iy)
	{
	  // charges of the final-final dipoles
	  _zij[ix][iy] = Ztemp*children[iy-1]->dataPtr()->iCharge();
	  _zij[iy][ix] = _zij[ix][iy];
	  // opening angles of the final-final dipoles
	  _cosij[ix][iy] = (_qdrf[ix].vect()).cosTheta(_qdrf[iy].vect());
	  _cosij[iy][ix] = _cosij[ix][iy];
	  _sinij[ix][iy] = sqrt(1.-sqr(_cosij[ix][iy]));
	  _sinij[iy][ix] = _sinij[ix][iy];
	  // mass of the final-final dipole
	  _mdipole[ix][iy] = (_qdrf[ix]+_qdrf[iy]).m2();
	  _mdipole[iy][ix] = _mdipole[ix][iy];
	  // crude photon multiplicity for final-final dipoles
	  _nbar[ix][iy] = nbar(ix,iy);
	  _nbar[iy][ix] = _nbar[ix][iy];
	}
    }
  //  cout << "testing the decay " << p.PDGName() << " -> ";
  //for(unsigned int ix=0;ix<children.size();++ix)
  //  {cout << children[ix]->PDGName() << " ";}
  //cout << endl;
  //cout << "testing the boost " << _boosttolab << endl;
  //cout << "testing number of decay products " << _nprod << endl;
  //cout << "testing the masses and velocities " << endl;
  //for(unsigned int ix=0;ix<=_nprod;++ix)
  //  {
  //    cout << _m[ix]      << " " << _m2[ix] << " " 
  //	   << _beta[ix]   << " " << _qdrf[ix].rho()/_qdrf[ix].e() << " "
  //	   << _ombeta[ix] << " " << 1.-_beta[ix] << endl; 
  // }
  //cout << "testing sum of masses " << _msum << endl;
  //cout << "testing maximum photon energy " << _emax << endl;
  Lorentz5Momentum ptemp;
  //cout << "testing the lab momenta " << endl;
  //for(unsigned int ix=0;ix<=_nprod;++ix)
  //  {
  //   cout << _qlab[ix] << endl;
  //    if(ix==0){ptemp=_qlab[ix];}
  //    else{ptemp-=_qlab[ix];}
  //  }
  //cout << "testing sum " << ptemp << endl;
  //cout << "testing the prf momenta " << endl;
  //for(unsigned int ix=0;ix<=_nprod;++ix)
  //  {
  //    cout << _qprf[ix] << endl;
  //    if(ix==0){ptemp=_qprf[ix];}
  //    else{ptemp-=_qprf[ix];}
  //  }
  //cout << "testing sum " << ptemp << endl;
  //cout << "testing the drf momenta " << endl;
  //for(unsigned int ix=0;ix<=_nprod;++ix)
  //  {
  //    cout << _qdrf[ix] << endl;
  //    if(ix==0){ptemp=_qdrf[ix];}
  //    else{ptemp-=_qdrf[ix];}
  //  }
  //cout << "testing sum " << ptemp << endl;
  //cout << "testing the dipole charges " << endl;
  //for(unsigned int ix=0;ix<=_nprod;++ix)
  //  {
  //    for(unsigned int iy=0;iy<=_nprod;++iy)
  //      {
  //        cout << _zij[ix][iy] << " ";
  //      }
  //    cout << endl;
  //  }
  //cout << "testing the dipole masses " << endl;
  //for(unsigned int ix=0;ix<=_nprod;++ix)
  //  {
  //    for(unsigned int iy=0;iy<=_nprod;++iy)
//	{
//	  cout << _mdipole[ix][iy] << " ";
//	}
//      cout << endl;
//    }
//  cout << "testing the dipole cos " << endl;
//  for(unsigned int ix=0;ix<=_nprod;++ix)
//    {
//      for(unsigned int iy=0;iy<=_nprod;++iy)
//	{
//	  cout << _cosij[ix][iy] << " ";
//	}
//      cout << endl;
//    }
//  cout << "testing the dipole sin " << endl;
//  for(unsigned int ix=0;ix<=_nprod;++ix)
//    {
//      for(unsigned int iy=0;iy<=_nprod;++iy)
//	{
//	  cout << _sinij[ix][iy] << " ";
//	}
//      cout << endl;
//    }
//  cout << "testing the dipole multiplicity " << endl;
//  for(unsigned int ix=0;ix<=_nprod;++ix)
//    {
//      for(unsigned int iy=0;iy<=_nprod;++iy)
//	{
//	  cout << _nbar[ix][iy] << " ";
//	}
//      cout << endl;
//    }
  // perform the unweighting
  double wgt;
  unsigned int ntry(0);
  do
    {
      //     cout << "testing called make photons " << endl;
      wgt = makePhotons();
      //cout << "testing got back ? " << wgt << endl;
      ++ntry;
      if(wgt>_maxwgt)
	{
	  generator()->log() << "Weight exceeds maximum for decay " 
			     << p.PDGName() << " ";
	  for(unsigned int ix=0;ix<_nprod;++ix)
	    {generator()->log() << children[ix]->PDGName() << " ";}
	  generator()->log() << "in GeneralDipole: resetting maximum weight." << endl
			     << " Old Maximum = " << _maxwgt 
			     << " New Maximum = " << wgt << endl;
	  _maxwgt=wgt;
	}
    }
  while(wgt<(_maxwgt*UseRandom::rnd())&&ntry<_maxtry);
  //  cout << "testing got to the end " << ntry << endl;
  if(ntry>=_maxtry)
    {
      generator()->log() << "GeneralDipole failed to generate QED " 
			 << "radiation for the decay " 
			 << p.PDGName() << " -> ";
      for(unsigned int ix=0;ix<_nprod;++ix)
	{generator()->log() << children[ix]->PDGName() << " ";}
      generator()->log() << endl;
      return children;
    }
  // produce products after radiation if needed
  if(_multiplicity>0)
    {
      // change the momenta of the children, they are currently
      // in original rest frame 
      for(unsigned int ix=1;ix<=_nprod;++ix)
	{
	  // unit vector along direction
	  Hep3Vector br(children[ix-1]->momentum().vect().unit());
	  // calculate the boost vector using expression accurate for beta->1
	  double bv(-(_ombetanew[ix]-_ombeta[ix])
		    /((1.+_beta[ix])*(1.-_betanew[ix])+_ombeta[ix]-_ombetanew[ix]));
	  br *=bv;
	  children[ix-1]->deepBoost(br);
	  br = _qnewdrf[0].findBoostToCM();
	  children[ix-1]->deepBoost(br);
	  // boost back to the lab
	  children[ix-1]->deepBoost(_boosttolab);
	}
      // add the photons to the event record
      tcPDPtr photon=getParticleData(ParticleID::gamma);
      for(unsigned int ix=0;ix<_multiplicity;++ix)
	{
	  // add if not removed because energy too low
	  if(!_photcut[ix])
	    {
	      PPtr newphoton=new_ptr(Particle(photon));
	      newphoton->set5Momentum(_llab[ix]);
	      children.push_back(newphoton);
	    }
	}
      /*
      ptemp=p.momentum();
      cout << "testing the parent " << ptemp << endl;
      for(unsigned int ix=0;ix<children.size();++ix)
	{
	  cout << "testing children " << ix << " " 
	       << children[ix]->PDGName() << " " << children[ix]->momentum() << endl;
	  ptemp-=children[ix]->momentum();
	}
      if(children.back()->id()==ParticleID::gamma){cout << "made a photon!!!" << endl;}
      cout << "testing total " << ptemp << endl;
      cout << "testing got to end ??? " << endl;
      */
      return children;
    }
  // otherwise just return the orginial particles
  else
    {
      for(unsigned int ix=0;ix<_nprod;++ix){children[ix]->deepBoost(_boosttolab);}
      ptemp=p.momentum();
      //cout << "testing the parent " << ptemp << endl;
      //for(unsigned int ix=0;ix<children.size();++ix)
      //	{
      //	  cout << "testing children " << ix << " " 
      //	       << children[ix]->PDGName() << " " << children[ix]->momentum() << endl;
      //	  ptemp-=children[ix]->momentum();
      //	}
      //cout << "testing total " << ptemp << endl;
      //cout << "testing got to end ??? " << endl;
      return children;
    }
}



// member which generates the photons
double GeneralDipole::makePhotons()
{
  // set the initial parameters
  // number of photons (zero)
  _multiplicity=0;
  // zero size of photon vectors
  _ldrf.resize(0);
  _lprf.resize(0);
  _llab.resize(0);
  // zero size of angle storage
  _photonwgt.resize(0);
  _photonemit.resize(0);
  _photonspect.resize(0);
  _cosphot.resize(0);
  _sinphot.resize(0);
  _photcut.resize(0);
  // zero total momenta of the photons
  _bigLdrf=Lorentz5Momentum();
  _bigLprf=Lorentz5Momentum();
  // set the initial values of the reweighting factors to one
  _dipolewgt   = 1.0;
  _yfswgt      = 1.0;
  _jacobianwgt = 1.0;
  _mewgt       = 1.0;
  // loop over the dipoles and generate the photons
  unsigned int nemit;
  for(unsigned int ix=0;ix<=_nprod;++ix)
    {
      for(unsigned int iy=ix+1;iy<=_nprod;++iy)
	{
	  // check if the dipole exists and has right charge
	  if(_zij[ix][iy]<0)
	    {
	      nemit=poisson(_nbar[ix][iy]);
	      // generate multiplicity using poisson
	      _yfswgt*=exp(_nbar[ix][iy]);
	      // if any photons produce them
	      if(nemit)
		{
		  // increment the multiplicity
		  _multiplicity+=nemit;
		  // generate the photons
		  photon(ix,iy,nemit);
		}
	    }
	}
    }
  // if photons produced
  //cout << "testing the multiplicity" << endl;
  if(_multiplicity>0)
    {
      // boost the momenta without any removal of low energy photons
      // resize arrays
      _photcut.resize(_multiplicity,false);
      _lprf.resize(_multiplicity);
      _llab.resize(_multiplicity);
      // perform the boost
      //cout << "testing first call " << endl;
      if(!boostMomenta(true)){return 0.;}
      // apply the cut on the photon energy if needed
      unsigned int nremoved(removePhotons());
      // redo the boost if we have removed photons
      //cout << "testing second call " << endl;
      if(nremoved!=0){if(!boostMomenta(nremoved!=_multiplicity)){return 0.;}}
      // reweight
      // masses of dipoles after radiation
      for(unsigned int ix=1;ix<=_nprod;++ix)
	{
	  _mdipolenew[0][ix]=(_qnewdrf[0]-_qnewdrf[ix]).m2();
	  _mdipolenew[ix][0]=_mdipolenew[0][ix];
	  for(unsigned int iy=ix+1;iy<=_nprod;++iy)
	    {
	      _mdipolenew[ix][iy]=(_qnewdrf[ix]+_qnewdrf[iy]).m2();
	      _mdipolenew[iy][ix]=_mdipolenew[ix][iy];
	    }
	}
      // calculate the velocities and 1-velocities for new momenta
      for(unsigned int ix=0;ix<=_nprod;++ix)
	{
	  if(ix==0&&_qnewdrf[ix].e()<_m[ix]){_betanew[ix]=0.;_ombetanew[ix]=1.;}
	  else
	    {
	      _betanew[ix]   = sqrt((_qnewdrf[ix].e()+_m[ix])*
				    (_qnewdrf[ix].e()-_m[ix]))/_qnewdrf[ix].e();
	      _ombetanew[ix] = sqr(_m[ix]/_qnewdrf[ix].e())/(1.+_betanew[ix]);
	    }
	}
      // angles between the decaying particle and the decay products have changed
      for(unsigned int ix=1;ix<=_nprod;++ix)
	{
	  _cosij[0][ix] = (_qdrf[0].vect()).cosTheta(_qdrf[ix].vect());
	  _cosij[ix][0] = _cosij[0][ix];
	}
      // calculate the new dipole weight
      reweightDipole();
      // remove the effect of the old cut
      _dipolewgt*=YFSFormFactor(_qdrf,_mdipole,_emin,false);
      //cout << "testing dipole weight A1 " 
      //	   << YFSFormFactor(_qdrf,_mdipole,_emin,false) << endl;
      //cout << "testing dipole weight A2 " << _dipolewgt << endl;
      // calculate the rest of the weight for photon removal
      if(_energyopt==1)
	{
	  //cout << "testing form factors B " << endl;
	  // yfs wgt
	  _yfswgt*=YFSFormFactor(_qnewprf,_mdipolenew,_eminrest,true);
	  //cout << "testing form factors C " << endl;
	  // find the momenta of the original particles in the new rest frame
	  const Lorentz5Momentum tmp = _qnewdrf[0];
	  HepLorentzRotation boost(HepBoost(tmp.findBoostToCM()));
	  vector<Lorentz5Momentum> ptemp;
	  for(unsigned int ix=0;ix<=_nprod;++ix){ptemp.push_back(boost*_qdrf[ix]);}
	  //cout << "testing form factors D " << endl;
	  // final part of the dipole wgt
	  _dipolewgt /=YFSFormFactor(ptemp,_mdipole,_eminrest,false);
	  //cout << "testing dipole wgt A3 " 
	  //    << YFSFormFactor(ptemp,_mdipole,_eminrest,false)
	  //     << endl;
	  //cout << "testing form factors E " << endl;
	}
      else 
	{
	  // yfs wgt
	  _yfswgt*=YFSFormFactor(_qnewlab,_mdipolenew,_eminlab,true);
	  // find the momenta of the original particles in the lab
	  const Lorentz5Momentum tmp = _qnewdrf[0];
	  HepLorentzRotation boost(HepBoost(tmp.findBoostToCM()));
	  boost = HepBoost(_boosttolab)*boost;
	  vector<Lorentz5Momentum> ptemp;
	  for(unsigned int ix=0;ix<=_nprod;++ix){ptemp.push_back(boost*_qdrf[ix]);}
	  // final part of the dipole wgt
	  _dipolewgt /=YFSFormFactor(ptemp,_mdipole,_eminlab,false);
	}
      //cout << "testing jacobian" << endl;
      // Calculating jacobian weight
      _jacobianwgt = jacobianWeight();
      //cout << "testing full " << endl;
      // different configuration part of the dipole weight
      _dipolewgt *=fullDipoleWeight();
      //cout << "testing got to end ??" << endl;
    }
  // otherwise copy momenta
  else
    {
      for(unsigned int ix=0;ix<=_nprod;++ix) 
	{
	  _qnewdrf[ix]=_qdrf[ix];
	  _qnewprf[ix]=_qprf[ix]; 
	  _qnewlab[ix]=_qlab[ix];
	}
      _jacobianwgt = 1.0;
      _dipolewgt   = 1.0;
      _yfswgt*=YFSFormFactor(_qdrf,_mdipole,_emin,true);
    }
  //cout << "testing end of make photons " << endl;
  // calculate the weight depending on the option
  double wgt;
  int nreal(0);
  for(unsigned int ix=0;ix<_multiplicity;++ix){if(!_photcut[ix]){++nreal;}}
  //cout << "testing the weights A " << nreal << " " << _multiplicity << " "  
  //     << _mewgt << " " << _yfswgt << " " << _jacobianwgt << " " 
  //     << _dipolewgt << endl;
  if(_mode==0){wgt=_maxwgt;}
  else if(_mode==1){wgt=_mewgt*_yfswgt*_jacobianwgt*_dipolewgt;}
  else if(_mode==2){wgt=_jacobianwgt*_yfswgt*_dipolewgt;}
  else if(_mode==3){wgt=_yfswgt*_dipolewgt;}
  else             {wgt=_yfswgt;}
  return wgt;
}


void GeneralDipole::photon(unsigned int i,unsigned int j,unsigned int nphoton)
{
  if(i==0)
    {
      // variables
      double phi,r,denom,costh,wgt,sinth;
      Energy en;
      // rotation matrix we'll need
      HepLorentzRotation rotation(HepRotationZ(-_qdrf[j].phi()));
      rotation.rotateY(-_qdrf[j].theta());
      rotation.invert();

      // generate the momenta
      for(unsigned int iphot=0;iphot<nphoton;++iphot)
	{
	  // generate the azimuthal angle
	  phi=2.*pi*UseRandom::rnd();
	  // generate the polar angle
	  r=UseRandom::rnd();
	  denom=pow(1.+_beta[j],1.-r)*pow(_ombeta[j],r);
	  costh=1./_beta[j]*(1.-denom);
	  sinth=sqrt(1.-sqr(costh));
	  // generate the ln(energy) uniformly in ln(_emin)->ln(_emax)
	  en=pow(_emax/_emin,UseRandom::rnd())*_emin;
	  // generate the momenta (this is wrt to the emitting particle)
	  _ldrf.push_back(Lorentz5Momentum(en*sinth*cos(phi),en*sinth*sin(phi),
					   en*costh,en,0.));
	  // perform the rotation
	  _ldrf.back() *=rotation;
	  // add the momentum to the total
	  _bigLdrf +=_ldrf.back();
	  // sort out the weights
	  // weight for the generation
	  wgt = 2./denom;
	  _dipolewgt /=wgt;
	  // including mass terms
	  wgt -= 1.+_ombeta[j]*(1.+_beta[j])/sqr(denom);
	  _photonwgt.push_back(wgt);
	  // emitter and spectator
	  _photonemit.push_back(j);
	  _photonspect.push_back(i);
	  // angles with respect to the emitter
	  _cosphot.push_back(costh);
	  _sinphot.push_back(sinth);
	  //cout << "testing the angles initial" << endl;
	  //cout << costh << " " << sinth << endl;
	  //cout << (_qdrf[j].vect()).cosTheta(_ldrf.back().vect()) << " " 
	  //     << sqrt(1.-sqr( (_qdrf[j].vect()).cosTheta(_ldrf.back().vect()))) << endl;
	  //cout << "testing store " << _cosphot.back() << " " << _sinphot.back() << endl;
	}

    }
  else if(j==0)
    {throw Exception() << "The lowest numbered particle in the dipole must be first"
		       << " in GeneralDipole::photon()" << Exception::abortnow;}
  else
    {
      // relative weights for the two terms
      double Pp(log((1+_beta[j])/_ombeta[j]));
      double Pm(log((1+_beta[i])/_ombeta[i]));
      Pm/=(Pp+Pm);
      // some variables we will need a lot
      double r1,r2,costh[2],sinth[2],denom1,denom2,phi,cphi,sphi,wgt;
      unsigned int ntry,ibranch,mtry(1000);
      Energy en;
      // storage for the rotation matrices
      HepLorentzRotation rotation[2]={HepLorentzRotation(),HepLorentzRotation()};
      // maximum weight for the unweighting
      double maxwgt(1./(_beta[i]+_beta[j]
			-_beta[i]*_beta[j]*sqrt(2.*(1.+_cosij[i][j]))));
      for(unsigned int iphot=0;iphot<nphoton;++iphot)
	{
	  ntry=0;
	  // generate the polar angle and azimuthal angles
	  do
	    {
	      ++ntry;
	      r1=UseRandom::rnd();
	      r2=UseRandom::rnd();
	      phi=2.*pi*UseRandom::rnd();
	      cphi = cos(phi);
	      sphi = sin(phi);
	      // 1/(1-beta_1c_1) branch
	      if(r1<=Pm)
		{
		  denom1   = pow(1.+_beta[i],1.-r2)*pow(_ombeta[i],r2);
		  costh[0] = 1./_beta[i]*(1.-denom1);
		  sinth[0] = sqrt(denom1*(2.-denom1)-(1.+_beta[i])*_ombeta[i]*sqr(costh[0]));
		  costh[1] = cphi*sinth[0]*_sinij[i][j]+costh[0]*_cosij[i][j];
		  sinth[1] = sqrt(1.-sqr(costh[1]));
		  denom2   = 1.-_beta[j]*costh[1];
		  ibranch=0;
		}
	      // 1/(1-beta_2c_2) branch 
	      else
		{
		  denom2   = pow(1.+_beta[j],1.-r2)*pow(_ombeta[j],r2);
		  costh[1] = 1./_beta[j]*(1.-denom2);
		  sinth[1] = sqrt(denom2*(2.-denom2)-(1.+_beta[j])*_ombeta[j]*sqr(costh[1]));
		  costh[0] = cphi*sinth[1]*_sinij[i][j]+costh[1]*_cosij[i][j];
		  sinth[0] = sqrt(1.-sqr(costh[0]));
		  denom1   = 1.-_beta[i]*costh[0];
		  ibranch=1;
		}
	      wgt = 1./(_beta[i]+_beta[j]-_beta[i]*_beta[j]*(costh[0]+costh[1]));
	      if(wgt>maxwgt){cerr << "testing violates maxweight" << endl;}
	    }
	  while(wgt<maxwgt*UseRandom::rnd()&&ntry<mtry);
	  if(ntry>=mtry)
	    {cerr << "testing too many tries " << endl;}
	  // generate the ln(energy) uniformly in ln(_emin)->ln(_emax)
	  en=pow(_emax/_emin,UseRandom::rnd())*_emin;
	  // generate the momenta (this is wrt to the emitting particle)
	  _ldrf.push_back(Lorentz5Momentum(en*sinth[ibranch]*cphi,en*sinth[ibranch]*sphi,
					   en*costh[ibranch],en,0.));
	  // identify emitter and spectator
	  unsigned int iemit(i),ispect(j);
	  if(ibranch==1){iemit=j;ispect=i;}
	  // calculate the rotation matrix if needed
	  if(rotation[ibranch].isIdentity())
	    {
	      rotation[ibranch]=HepRotationZ(-_qdrf[iemit].phi());
	      rotation[ibranch].rotateY(-_qdrf[iemit].theta());
	      double phiy((rotation[ibranch]*_qdrf[ispect]).phi());
	      rotation[ibranch].rotateZ(-phiy);
	      // the matrix is the inverse of this
	      rotation[ibranch].invert();
	    }
	  // perform the rotation
	  _ldrf.back() *=rotation[ibranch];
	  // add the momentum to the total
	  _bigLdrf +=_ldrf.back();
	  // sort out the weights
	  // weight for the generation
	  wgt = 2.*(1.-_beta[i]*_beta[j]*_cosij[i][j])/denom1/denom2;
	  _dipolewgt /=wgt;
	  // including mass terms
	  wgt -= _ombeta[i]*(1.+_beta[i])/sqr(denom1)+_ombeta[j]*(1.+_beta[j])/sqr(denom2);
	  _photonwgt.push_back(wgt);
	  // emitter and spectator
	  _photonemit.push_back(iemit);
	  _photonspect.push_back(ispect);
	  // angles with respect to the emitter
	  _cosphot.push_back(costh[ibranch]);
	  _sinphot.push_back(sinth[ibranch]);
	  //cout << "testing the angles final " << ibranch << endl;
	  //cout << costh[0] << " " << sinth[0] << endl;
	  //cout << costh[1] << " " << sinth[1] << endl;
	  //cout << (_qdrf[i].vect()).cosTheta(_ldrf.back().vect()) << " " 
	  //     << sqrt(1.-sqr( (_qdrf[i].vect()).cosTheta(_ldrf.back().vect()))) << endl;
	  //cout << (_qdrf[j].vect()).cosTheta(_ldrf.back().vect()) << " " 
	  //     << sqrt(1.-sqr( (_qdrf[j].vect()).cosTheta(_ldrf.back().vect()))) << endl;
	  //cout << "testing store " << _cosphot.back() << " " << _sinphot.back() << endl;
	}
    }
}

bool GeneralDipole::boostMomenta(bool photons)
{
  // if photons need to do all the boosts
  if(photons)
    {
      // total energy  and momentum of photons
      Energy L0(_bigLdrf.e()),modL(_bigLdrf.rho());
      // invariant mass after radiation
      _roots=-L0+sqrt(_m[0]*_m[0]+modL*modL);
      //      cout << "testing the masses " << _roots << " " << _m[0] << endl;
      // check this is allowed
      if(_roots<_msum){return false;}
      // compute all the magnitudes as a vector
      vector<long double> absq2(_nprod+1);
      for(unsigned int ix=1;ix<=_nprod;++ix)
	{absq2[ix]=_qdrf[ix].rho();absq2[ix]*=absq2[ix];}
      // find the rescaling factor using Newton-Raphson
      long double uold(1.),u(0.),f,fprime,err;
      vector<long double> en(_nprod+1);
      const long double tol(1e-12);
      do
	{
	  f = -_roots;
	  fprime=0.0;
	  for(unsigned int ix=1;ix<=_nprod;++ix)
	    {
	      en[ix]=sqrt(uold*uold*absq2[ix]+_m2[ix]);
	      f      +=en[ix];
	      fprime +=absq2[ix]/en[ix];
	    }
	  u = uold- f/(uold*fprime);
	  err = abs(u-uold);
	  uold=u;
	}
      while(err>tol);
      _rescalingfactor=u;
      //cout << "testing rescaling factor " << _rescalingfactor << endl;
      // calculate the rescaled momenta
      for(unsigned int ix=1;ix<=_nprod;++ix)
	{
	  _qnewdrf[ix]=_rescalingfactor*_qdrf[ix];
	  _qnewdrf[ix].setMass(_m[ix]);
	  _qnewdrf[ix].rescaleEnergy();
	}
      _qnewdrf[0]=Lorentz5Momentum(_bigLdrf.x(),_bigLdrf.y(),
				   _bigLdrf.z(),_bigLdrf.e(),_m[0]);
      _qnewdrf[0].rescaleEnergy();
      // Find the momenta of the particles in the rest frame of the parent.
      const Lorentz5Momentum tmp = _qnewdrf[0];
      HepLorentzRotation boost(HepBoost(tmp.findBoostToCM()));
      // Boost the momenta of the charged particles
      for(unsigned int ix=0;ix<=_nprod;++ix){_qnewprf[ix]=boost*_qnewdrf[ix];}
      // Boost the total photon momentum
      _bigLprf=boost*_bigLdrf;
      // Boost the individual photon momenta
      for(unsigned int ix=0;ix<_multiplicity;++ix){_lprf[ix]=boost*_ldrf[ix];}
      // Now boost from the parent rest frame to the lab frame
      boost = HepBoost(_boosttolab);
      // Boosting charged particles
      for(unsigned int ix=0;ix<=_nprod;++ix){_qnewlab[ix]=boost*_qnewprf[ix];}
      // Boosting total photon momentum
      _bigLlab=boost*_bigLprf;
      // Boosting individual photon momenta
      for(unsigned int ix=0;ix<_multiplicity;++ix){_llab[ix]=boost*_lprf[ix];}
      // check everything is O.K.
      //Lorentz5Momentum ptemp;
      //cout << "testing the new momenta in dipole frame" << endl;
      //cout << "particles" << endl;
      //for(unsigned int ix=0;ix<=_nprod;++ix)
      //	{
      //	  cout << _qnewdrf[ix] << endl;
      //	  if(ix>0){ptemp+=_qnewdrf[ix];}
      //	}
      //cout << "testing A " << ptemp.m2() << endl;
      //cout << "photons" << endl;
      //for(unsigned int ix=0;ix<_multiplicity;++ix)
      //	{cout << _ldrf[ix] << endl;ptemp+=_ldrf[ix];}
      //cout << "testing B " << ptemp.m2() << endl;
      //cout << "testing total " << ptemp << endl;;
      //cout << "testing the new momenta in parent frame" << endl;
      //cout << _qprf[0] << endl;
      //ptemp=Lorentz5Momentum();
      //cout << "particles" << endl;
      //for(unsigned int ix=0;ix<=_nprod;++ix)
      //	{
      //	  cout << _qnewprf[ix] << endl;
      //	  if(ix>0){ptemp+=_qnewprf[ix];}
      //	}
      //cout << "testing A " << ptemp.m2() << endl;
      //cout << "photons" << endl;
      //for(unsigned int ix=0;ix<_multiplicity;++ix)
      //	{cout << _lprf[ix] << endl;ptemp+=_lprf[ix];}
      //cout << "testing B " << ptemp.m2() << endl;
      //cout << "testing total " << ptemp << endl;
      //cout << "testing the new momenta in lab frame" << endl;
      //ptemp=Lorentz5Momentum();
      //cout << _qlab[0] << endl;
      //cout << "particles" << endl;
      //for(unsigned int ix=0;ix<=_nprod;++ix)
      //	{
      //	  cout << _qnewlab[ix] << endl;
      //	  if(ix>0){ptemp+=_qnewlab[ix];}
      //	}
      //cout << "photons" << endl;
      //for(unsigned int ix=0;ix<_multiplicity;++ix)
      //	{cout << _llab[ix] << endl;ptemp+=_llab[ix];}
      //cout << "testing total " << ptemp << endl;
    }
  else
    {
      // copy momenta
      for(unsigned int ix=0;ix<=_nprod;++ix)
	{
	  _qnewdrf[ix]=_qdrf[ix];
	  _qnewprf[ix]=_qprf[ix];
	  _qnewlab[ix]=_qlab[ix];
	}
      // relevant variable
      _roots=_m[0];
      _rescalingfactor=1.;
      // total photon energy
      _bigLprf=Lorentz5Momentum();
      _bigLdrf=Lorentz5Momentum();
      _bigLlab=Lorentz5Momentum();
    }
  return true;
}

unsigned int GeneralDipole::removePhotons()
{
  unsigned int nremoved(0);
  // apply the cut in the rest frame
  if(_energyopt==1)
    {
      for(unsigned int ix=0;ix<_multiplicity;++ix)
	{
	  if(_lprf[ix].e()<_eminrest)
	    {
	      ++nremoved;
	      _photcut[ix]=true;
	      _bigLdrf-=_ldrf[ix];
	      _ldrf[ix]=Lorentz5Momentum();
	    }
	}
    }
  // apply the cut in the lab frame
  else if(_energyopt==2)
    {
      for(unsigned int ix=0;ix<_multiplicity;++ix)
	{
	  if(_llab[ix].e()<_eminlab)
	    {
	      ++nremoved;
	      _photcut[ix]=true;
	      _bigLdrf-=_ldrf[ix];
	      _ldrf[ix]=Lorentz5Momentum();
	    }
	}
    }
  // correction factor for dipoles
  if(nremoved!=0)
    {for(unsigned int ix=0;ix<_multiplicity;++ix)
	{if(_photcut[ix]){_dipolewgt*=_photonwgt[ix];}}}
  // return number of removed photons
  return nremoved;
}

double GeneralDipole::YFSFormFactor(const vector<Lorentz5Momentum> & mom,
				    const vector<vector<Energy2> > mdipole,
				    const Energy & ecut,
				    bool poscharge) const
{
  vector<double> btemp(_nprod+1),ombtemp(_nprod+1);
  double output(1.);
  // calculate the velocities
  for(unsigned int ix=0;ix<=_nprod;++ix)
    {
      if(ix==0&&mom[ix].e()<_m[ix]){btemp[ix]=0.;ombtemp[ix]=1.;}
      else
	{
	  btemp[ix]   = sqrt((mom[ix].e()+_m[ix])*(mom[ix].e()-_m[ix]))/mom[ix].e();
	  ombtemp[ix] = sqr(_m[ix]/mom[ix].e())/(1.+btemp[ix]);
	}
    }
  for(unsigned int ix=0;ix<=_nprod;++ix)
    {
      for(unsigned int iy=ix+1;iy<=_nprod;++iy)
	{
	  if((_zij[ix][iy]&&poscharge)||(_zij[ix][iy]<0&&!poscharge))
	    {
	      if(ix==0)
		{output*=YFSFormFactors::exponentialYFSIF(btemp[ix],ombtemp[ix],
							  btemp[iy],ombtemp[iy],
							  mom[ix].e(),mom[iy].e(),
							  _m[ix],_m[iy],mdipole[ix][iy],
							  _zij[ix][iy]/9.,ecut);}
	      else
		{output*=YFSFormFactors::exponentialYFSFF(btemp[ix],ombtemp[ix],
							  btemp[iy],ombtemp[iy],
							  mom[ix].e(),mom[iy].e(),
							  _m[ix],_m[iy],mdipole[ix][iy],
							  _zij[ix][iy]/9.,ecut);}
	    }
	}
    }
  return output;
}

double GeneralDipole::fullDipoleWeight()
{
  // compute the weights
  vector<vector<double> > dipole;
  double wgt,wtemp,wtotal(1.);
  // calculate the dipole weight for the different dipoles
  for(unsigned int iphot=0;iphot<_multiplicity;++iphot)
    {
      if(!_photcut[iphot])
	{
	  wgt=0.;
	  dipole.push_back(vector<double>(0));
	  for(unsigned int ix=0;ix<=_nprod;++ix)
	    {
	      for(unsigned int iy=ix+1;iy<=_nprod;++iy)
		{
		  // calculate the weight
		  if(_zij[ix][iy])
		    {
		      wtemp=dipoleWeight(ix,iy,iphot);
		      // store the answer if needed
		      if(_zij[ix][iy]<0){dipole.back().push_back(wtemp);}
		      // add to the total
		      wgt+=wtemp;
		    }
		}
	    }
	  // multiply the total by the weight
	  wtotal*=wgt;
	}
    }
  if(dipole.size()==0){return 1.;}
  // now we have the weights for all the dipoles
  vector<unsigned int> iphot(dipole.size(),0),nmult(dipole[0].size());
  unsigned int ndipole(dipole[0].size()),ix,iy;
  int ipart;
  // the factorial we will need
  double nfact(1.);
  for(ix=1;ix<=dipole.size();++ix){nfact*=ix;}
  // calculate the denominator
  double denom(0.),pre;
  bool contin;
  do
    {
      // compute the weight for this piece
      // zero multiplicities
      for(ix=0;ix<nmult.size();++ix){nmult[ix]=0;}
      // dipole part of the weight
      wgt=1.;
      for(ix=0;ix<dipole.size();++ix)
	{
	  wgt*=dipole[ix][iphot[ix]];
	  ++nmult[iphot[ix]];
	}
      // combinatoric prefactor
      pre=1.;
      for(iy=0;iy<nmult.size();++iy){for(ix=1;ix<=nmult[iy];++ix){pre*=ix;}}
      pre /=nfact;
      // total weight for configuration
      wgt *=pre;
      // add to the total
      denom+=wgt;
      // work out which one to do next
      ipart=iphot.size()-1;
      do
	{
	  // got to end
	  if(iphot[ipart]==ndipole-1)
	    {
	      iphot[ipart]=0;
	      --ipart;
	      contin=true;
	    }
	  else
	    {
	      ++iphot[ipart];
	      contin=false;
	    }
	}
      while(contin&&ipart>=0);
    }
  while(ipart>=0);
  // return the answer
  return wtotal/denom;
}
