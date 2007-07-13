// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FortranSudakov class.
//

#include "FortranSudakov.h"
#include "FortranShowerKinematics.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "Herwig++/Shower/Couplings/FortranAlphaQCD.h"
#include "ThePEG/Handlers/LuminosityFunction.h"
#include "Herwig++/Utilities/GaussianIntegrator.h"

using namespace Herwig;

IBPtr FortranSudakov::clone() const {
  return new_ptr(*this);
}

IBPtr FortranSudakov::fullclone() const {
  return new_ptr(*this);
}

FortranSudakov::FortranSudakov() : _sudord(1), _inter(3), 
				   _vqcut(0.48*GeV),_vgcut(0.10*GeV),
				   _vpcut(0.40*GeV),_nqev(1024)
{}

void FortranSudakov::persistentOutput(PersistentOStream & os) const {
  os << _sudord << _inter << ounit(_vqcut,GeV) << ounit(_vgcut,GeV) << ounit(_vpcut,GeV) << _nqev 
     << _sudakovQ << _sudakovP << _linearQ << _linearP << _fortranQCD;
}

void FortranSudakov::persistentInput(PersistentIStream & is, int) {
  is >> _sudord >> _inter >> iunit(_vqcut,GeV) >> iunit(_vgcut,GeV) >> iunit(_vpcut,GeV) >> _nqev 
     >> _sudakovQ >> _sudakovP >> _linearQ >> _linearP >> _fortranQCD;
}

ClassDescription<FortranSudakov> FortranSudakov::initFortranSudakov;
// Definition of the static class description member.

void FortranSudakov::Init() {

  static ClassDocumentation<FortranSudakov> documentation
    ("The FortranSudakov class is designed to implement the same Sudakov"
     "form factors as in FORTRAN HERWIG in exactly the same way");

  static Switch<FortranSudakov,unsigned int> interfaceSudakovOrder
    ("SudakovOrder",
     "The order of alpha_S to use for the Sudakov tables. In both cases 2nd"
     " order alpha_S is used but the switch effects how this is achieved.",
     &FortranSudakov::_sudord, 1, false, false);
  static SwitchOption interfaceSudakovOrder1stOrder
    (interfaceSudakovOrder,
     "1stOrder",
     "Use first order in the tables and obtain second order using the veto algorithm",
     1);
  static SwitchOption interfaceSudakovOrder2ndOrder
    (interfaceSudakovOrder,
     "2ndOrder",
     "Use second order in the tables",
     2);

  static Parameter<FortranSudakov,unsigned int> interfaceInterpolationOrder
    ("InterpolationOrder",
     "The order of interpolation to use for the Sudakov tables.",
     &FortranSudakov::_inter, 3, 1, 10,
     false, false, Interface::limited);

  static Parameter<FortranSudakov,Energy> interfaceVQCut
    ("VQCut",
     "Cut-off on quark virtual mass",
     &FortranSudakov::_vqcut, GeV, 0.48*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<FortranSudakov,Energy> interfaceVGCut
    ("VFCut",
     "Cut-off on gluon virtual mass",
     &FortranSudakov::_vgcut, GeV, 0.10*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<FortranSudakov,Energy> interfaceVPCut
    ("VPCut",
     "Cut-off on photon virtual mass",
     &FortranSudakov::_vpcut, GeV, 0.40*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<FortranSudakov,unsigned int> interfaceNumberEntries
    ("NumberEntries",
     "The nmber of entries in the lookup tables",
     &FortranSudakov::_nqev, 1024, 10, 10000,
     false, false, Interface::limited);

}

void FortranSudakov::doinit() throw(InitException) {
  SudakovFormFactor::doinit();
  // get the maximum energy from the event handler
  Energy qlim=generator()->eventHandler()->lumiFnPtr()->maximumCMEnergy();
  // power for interpolation
  double power=1./(_nqev-1);
  // get the constituent masses of the quarks
  vector<Energy> qmass;
  for(unsigned int ix=1;ix<=6;++ix) {
    qmass.push_back(getParticleData(ix)->constituentMass());
  }
  _fortranQCD=dynamic_ptr_cast<FortranAlphaQCDPtr>(alpha());
  if(!_fortranQCD) throw InitException() << "Must be using FortranAlphaQCD in "
					 << "FortranSudakov::doinit()";
  // ensure strong coupling is initialised
  _fortranQCD->init();
  // get lambda
  Energy lam3=_fortranQCD->lambdaQCDThree();
  // integrator
  GaussianIntegrator integrator;
  // compute the sudakovs and construct the interpolation tables
  for(unsigned int ix=0;ix<particles().size();++ix) {
    vector<double> sud,qev;
    Energy q1=cutOff(particles()[ix][1]);
    Energy q2=cutOff(particles()[ix][2]);
    Energy qmin=q1+q2;
    double qfac=pow(1.1*qlim/qmin,power);
    sud.push_back(1.);
    qev.push_back(qmin/GeV);
    Energy qnow;
    double qlam,qrat,zmin,zmax;
    // create the integrands for the two regions
    FortranSudakovIntegrand ilower(lam3,splittingFn(),particles()[ix],_sudord,true,
				   qmass,_fortranQCD);
    FortranSudakovIntegrand iupper(lam3,splittingFn(),particles()[ix],_sudord,false,
				   qmass,_fortranQCD);
    // compute the table
    do {
      qnow = qfac*qev.back()*GeV;
      qlam=qnow/lam3;
      // upper part of the integral
      zmin = q2/qnow;
      qrat = 1./zmin;
      zmax = q2/qmin;
      iupper.setScales(qrat,qlam);
      double g1=0.;
      for(unsigned int ix=3;ix<7;++ix) {
	double zlo = zmin;
	double zhi = zmax;
	if(ix!=6) zlo = max(zlo,q2/qmass[ix]);
	if(ix!=3) zhi = min(zhi,q2/qmass[ix-1]);
	if(zhi>zlo) g1+=integrator.value(iupper,log(zlo),log(zhi));
      }
      // lower part of the the integral      
      zmin = q1/qnow;
      qrat = 1./zmin;
      zmax = q1/qmin;
      ilower.setScales(qrat,qlam);
      double g2=0.;
      for(unsigned int ix=3;ix<7;++ix) {
	double zlo = zmin;
	double zhi = zmax;
	if(ix!=6) zlo = max(zlo,q1/qmass[ix]);
	if(ix!=3) zhi = min(zhi,q1/qmass[ix-1]);
	if(zhi>zlo) g2+=integrator.value(ilower,log(zlo),log(zhi));
      }
      sud.push_back(exp(_fortranQCD->scaleFactor()*(g1+g2)));
      qev.push_back(qnow/GeV);
      //cerr << qev.back() << " " << sud.back() << "\n";
    }
    while (qnow<1.1*qlim);
    //cerr << "join red\n";
    // now construct the interpolators
    Interpolator<double,double>::Ptr newint = make_InterpolatorPtr(sud,qev,_inter);
    _sudakovQ.push_back(newint);
    newint = make_InterpolatorPtr(qev,sud,_inter);
    _sudakovP.push_back(newint);
    newint = make_InterpolatorPtr(sud,qev,1);
    _linearQ.push_back(newint);
    newint = make_InterpolatorPtr(qev,sud,1);
    _linearP.push_back(newint);
  }
}

ShoKinPtr FortranSudakov::generateNextTimeBranching(const Energy startingScale,
						    const IdList &ids,const bool cc,
						    double enhance) {
  // can only enhance by integer factor here
  unsigned int nrej=int(enhance);
  if(enhance-double(nrej)>0.) ++nrej;
  // get the minimum scale
  Energy qmin=cutOff(ids[1])+cutOff(ids[2]);
  cerr << "testing ids " << ids[1] << " " << ids[2] << "\n";
  cerr << "testing ids " << cutOff(ids[1])/MeV << " " << cutOff(ids[2])/MeV << "\n";
  cerr << "testing scales " << startingScale/MeV << " " << qmin/MeV << "\n";
  // return if no allowed branching
  if(qmin>startingScale) return ShoKinPtr();
  // find which sudakov table to use
  unsigned int iloc=findSudakov(ids,cc);
  cerr << "testing sudakov" << iloc << "\n";
  cerr << "testing starting scale " << startingScale/MeV << " " << qmin/MeV << "\n"; 
  Energy qnow=-1.*GeV;
  for(unsigned int irej=0;irej<nrej;++irej) {
    double rn=UseRandom::rnd();
    double slst=(*_sudakovQ[iloc])(startingScale/GeV);
    cerr << "testing slst " << slst <<"\n";
    double snow = rn==0. ? 2. : slst/rn;
    cerr << "testing snow " << snow << " " << rn << "\n";
    Energy qsud=-1.*GeV;
    if(snow<1.) {
      qsud=(*_sudakovP[iloc])(snow)*GeV;
      // if normal interpolation fail use linear
      if(qsud>startingScale) {
	snow=(*_linearQ[iloc])(startingScale/GeV)/rn;
	qsud=(*_linearP[iloc])(snow)*GeV;
	if(qsud>startingScale) {
	  generator()->log() << "Failed to find scale for " 
			     << ids[0] << " " << ids[1] << " " << ids[2] 
			     << " in FortranSudakov::generateNextTimeBranching\n";
	  qsud=-1.*GeV;
	}
      }
    }
    if(qsud>qmin&&qsud>qnow) qnow=qsud;
  }
  cerr << "testing selected" << qnow/GeV << "\n";
  // if no branching selected return
  if(qnow<0.*MeV) return ShoKinPtr();
  // limits on z
  double zmin=cutOff(ids[1])/qnow;
  double zmax=1.-cutOff(ids[2])/qnow;
  cerr << "testing limits " << zmin << " " << zmax << "\n";
  double fmax(splittingFn()->integOverP(zmax)),fmin(splittingFn()->integOverP(zmin));
  double z;
  do {
    z=splittingFn()->invIntegOverP(fmin+UseRandom::rnd()*(fmax-fmin));
  }
  while(splittingFn()->ratioP(z,sqr(qnow),ids,false)<=UseRandom::rnd());
  cerr << "testing selected " << z << endl;
  // create the ShowerKinematics object
  FortranShowerKinematicsPtr showerKin = new_ptr(FortranShowerKinematics());
  showerKin->scale(qnow);
  showerKin->z(z);
  showerKin->phi(2.*pi*UseRandom::rnd());
  //showerKin->pT(pT());
  showerKin->splittingFn(splittingFn());
  // return it
  return showerKin;
}

ShoKinPtr FortranSudakov::generateNextDecayBranching(const Energy startingScale,
						     const Energy stoppingScale,
						     const Energy minmass,
						     const IdList &ids,
						     const bool cc,
						     double enhance) {
  throw Exception() << "FortranSudakov::generateNextDecayBranching() not yet "
		    << "implemented" << Exception::runerror;
  return ShoKinPtr();
}

ShoKinPtr FortranSudakov::
generateNextSpaceBranching(const Energy startingScale,
			   const IdList &ids,double x,
			   const bool cc,
			   double enhance,
			   Ptr<BeamParticleData>::transient_const_pointer beam) {
  throw Exception() << "FortranSudakov::generateNextSpaceBranching() not yet "
		    << "implemented" << Exception::runerror;
  return ShoKinPtr();
}

FortranSudakovIntegrand::FortranSudakovIntegrand(Energy qcdlam, 
						 SplittingFnPtr split, 
						 IdList ids,
						 unsigned int sudord, bool zord,
						 vector<Energy> qmass,
						 FortranAlphaQCDPtr alpha)
  : _qcdlam(qcdlam), _split(split), _ids(ids), _sudord(sudord), _zord(zord), 
    _alpha(alpha)
{
  // compute the beta function coefficients
  for(unsigned int nf=3;nf<=6;++nf) {
    _bet.push_back((11.*3.-2.*nf)/(12.*pi));
    _bep.push_back((17*9.-(5*3.+4.)*nf)/(24.*sqr(pi))/_bet.back());
    if(nf==3) {
      _mumi.push_back(0.*MeV);
      _almi.push_back(1e30);
    }
    else {
      _mumi.push_back(qmass[nf-1]);
      _almi.push_back(_alpha->value(sqr(_mumi.back())));
    }
    if(nf==6) {
      _muma.push_back(1.e30*GeV);
      _alma.push_back(0.);
    }
    else {
      _muma.push_back(qmass[nf  ]);
      _alma.push_back(_alpha->value(sqr(_muma.back())));
    }
    if(nf!=3&&nf!=6) {
      _fint.push_back(alphaIntegral(_almi[nf-3],_alma[nf-3],nf));
    }
    else             _fint.push_back(0.);
  }
}

double FortranSudakovIntegrand::operator() (double zlog) const {
  // convert log(z) to z
  double z = exp(zlog);
  double u = 1.-z;
  // compute the alpha_S piece
  double output=0.;
  // leading order
  if(_sudord==1) {
    double al=log(_qrat*z);
    double bl=log(_qlam*u*z);
    output=log(1.-al/bl)/_bet[2];
  }
  // 2 loop
  else if(_sudord==2) {
    Energy qnow=_qlam*_qcdlam;
    Energy qmin=qnow/_qrat;
    Energy mumin=u*qmin;
    Energy mumax=z*u*qnow;
    if(mumax<=mumin) return 0.;
    double almin=_alpha->value(sqr(mumin));
    double almax=_alpha->value(sqr(mumax));
    unsigned int nf=3;
    while(mumin>_muma[nf-3]) ++nf;
    if(mumax<_muma[nf-3]) output=alphaIntegral(almin,almax,nf);
    else {
      output=alphaIntegral(almin,_alma[nf-3],nf);
      ++nf;
      while (mumax>_muma[nf-3]) {
	output+=_fint[nf-3];
	++nf;
      }
      output+=alphaIntegral(_almi[nf-3],almax,nf);
    }
  }
  // compute the splitting function piece
  if(_zord) output*=z*_split->P(   z,0.*MeV2,_ids,false);
  else      output*=z*_split->P(1.-z,0.*MeV2,_ids,false);
  return output/2./pi;
}

Energy FortranSudakov::calculateScale(double z, Energy pt, IdList ids,
					 unsigned int iopt) {
  throw Exception() << "Base class udakovFormFactor::calculateScale() called "
		    << "this should be overidden in the inheriting class"
		    << Exception::runerror;
}

ShoKinPtr FortranSudakov::createFinalStateBranching(Energy scale,double z,
				    double phi, Energy pt) {
  throw Exception() << " FortranSudakov::createFinalStateBranching() not implemented"
		    << Exception::runerror;
}

ShoKinPtr FortranSudakov::createInitialStateBranching(Energy scale,double z,
						 double phi, Energy pt) {
  throw Exception() << " FortranSudakov::createInitialStateBranching() not implemented"
		    << Exception::runerror;
}
