// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EtaPiPiGammaDecayer class.
//

#include "EtaPiPiGammaDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "EtaPiPiGammaDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/ScalarSpinInfo.h"
#include "ThePEG/Helicity/VectorSpinInfo.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig++/Helicity/EpsFunction.h"
#include "Herwig++/PDT/ThreeBodyAllOnCalculator.h"

namespace Herwig{
using namespace ThePEG;
using namespace ThePEG::Helicity;
using ThePEG::Helicity::RhoDMatrix;
using ThePEG::Helicity::LorentzPolarizationVector;
using ThePEG::Helicity::tcScalarSpinPtr;
using ThePEG::Helicity::VectorSpinPtr;
using ThePEG::Helicity::ScalarSpinInfo;
using ThePEG::Helicity::VectorSpinInfo;
using Helicity::Direction;
using Helicity::incoming;
using Helicity::outgoing;
using Helicity::VectorWaveFunction;
using Helicity::EpsFunction;

  EtaPiPiGammaDecayer::~EtaPiPiGammaDecayer() 
  {
    /*
    if(_Oreal){delete _Oreal;}
    if(_Oimag){delete _Oimag;}
    */
  }

bool EtaPiPiGammaDecayer::accept(const DecayMode & dm) const {
  // check number of external particles
  if(dm.products().size()!=3){return false;}
  // check the outgoing particles
  unsigned int npip=0,npim=0,ngamma=0;
  ParticleMSet::const_iterator pit = dm.products().begin();
  int id;
  for(;pit!=dm.products().end();++pit)
    {
      id=(**pit).id();
      if(id==ParticleID::piplus){++npip;}
      else if(id==ParticleID::piminus){++npim;}
      else if(id==ParticleID::gamma){++ngamma;}
    }
  if(!(npip==1&&npim==1&&ngamma==1)){return false;}
  // and the incoming particle
  id=dm.parent()->id();
  for(unsigned int ix=0;ix<_incoming.size();++ix)
    {if(id==_incoming[ix]){return true;}}
  return false;
}

ParticleVector EtaPiPiGammaDecayer::decay(const DecayMode & dm,
				  const Particle & parent) const {
  int imode=-1; unsigned int ix=0;
  int id=parent.id();
  do
    {
      if(id==_incoming[ix]){imode=ix;}
      ++ix;
    }
  while(imode<0&&ix<_incoming.size());
  bool cc=false;
  return generate(true,cc,imode,parent);
}


void EtaPiPiGammaDecayer::persistentOutput(PersistentOStream & os) const {
  os << _fpi << _incoming << _coupling << _maxweight << _option << _aconst 
     << _cconst <<_mrho << _rhowidth << _rhoconst << _mpi << _localparameters
     << _energy << _Omnesenergy 
     << _phase << _Omnesfunctionreal << _Omnesfunctionimag << _initialize
     << _npoints;
}

void EtaPiPiGammaDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _fpi >> _incoming >> _coupling >> _maxweight >> _option >> _aconst 
     >> _cconst >>_mrho >> _rhowidth >> _rhoconst >> _mpi >> _localparameters
     >> _energy >> _Omnesenergy 
     >> _phase >> _Omnesfunctionreal >> _Omnesfunctionimag >> _initialize
     >> _npoints;
}

ClassDescription<EtaPiPiGammaDecayer> EtaPiPiGammaDecayer::initEtaPiPiGammaDecayer;
// Definition of the static class description member.

void EtaPiPiGammaDecayer::Init() {

  static ClassDocumentation<EtaPiPiGammaDecayer> documentation
    ("The \\classname{EtaPiPiGammaDecayer} class is design for the decay of"
     " the eta and eta prime to pi+pi-gamma");


  static Parameter<EtaPiPiGammaDecayer,Energy> interfacefpi
    ("fpi",
     "The pion decay constant",
     &EtaPiPiGammaDecayer::_fpi, MeV, 130.7*MeV, -1.0e12*MeV, 1.0e12*MeV,
     false, false, false); 

  static ParVector<EtaPiPiGammaDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &EtaPiPiGammaDecayer::_incoming,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<EtaPiPiGammaDecayer,double> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &EtaPiPiGammaDecayer::_coupling,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<EtaPiPiGammaDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &EtaPiPiGammaDecayer::_maxweight,
     0, 0, 0, -10000, 10000, false, false, true);

  static Parameter<EtaPiPiGammaDecayer,Energy> interfaceRhoMass
    ("RhoMass",
     "The mass of the rho",
     &EtaPiPiGammaDecayer::_mrho, MeV, 771.1*MeV, -1.0e12*MeV, 1.0e12*MeV,
     false, false, false);

  static Parameter<EtaPiPiGammaDecayer,Energy> interfaceRhoWidth
    ("RhoWidth",
     "The width of the rho",
     &EtaPiPiGammaDecayer::_rhowidth, MeV, 149.2*MeV, -1.0e12*MeV, 1.0e12*MeV,
     false, false, false);

  static Switch<EtaPiPiGammaDecayer,bool> interfaceLocalParameters
    ("LocalParameters",
     "Use local values of the rho mass and width",
     &EtaPiPiGammaDecayer::_localparameters, true, false, false);
  static SwitchOption interfaceLocalParametersLocal
    (interfaceLocalParameters,
     "Local",
     "Use local parameters",
     true);
  static SwitchOption interfaceLocalParametersParticleData
    (interfaceLocalParameters,
     "ParticleData",
     "Use values from the particle data objects",
     false);

  static Parameter<EtaPiPiGammaDecayer,double> interfaceOmnesC
    ("OmnesC",
     "The constant c for the Omnes form of the prefactor",
     &EtaPiPiGammaDecayer::_cconst, 1.0, -1.0e12, 1.0e12,
     false, false, false);

  static Parameter<EtaPiPiGammaDecayer,InvEnergy2> interfaceOmnesA
    ("OmnesA",
     "The constant a for the Omnes form of the prefactor",
     &EtaPiPiGammaDecayer::_aconst, 1./GeV2, 0.8409082/GeV2, -1.0e12*1./GeV2,
     1.0e12*1./GeV2,
     false, false, false);

  static ParVector<EtaPiPiGammaDecayer,int> interfaceOption
    ("Option",
     "The form of the prefactor 0 is a VMD model using M Gamma for the width term,"
     "1 is a VMD molde using q Gamma for the width term,"
     "2. analytic form of the Omnes function,"
     "3. experimental form of the Omnes function.",
     &EtaPiPiGammaDecayer::_option,
     0, 0, 0, 0, 4, false, false, true);

}

double EtaPiPiGammaDecayer::me2(bool vertex,const int,const Particle & inpart,
				 const ParticleVector & decay) const
{
  // check if the decay particle has spin info 
  tcScalarSpinPtr inspin;
  if(inpart.spinInfo())
    {inspin = dynamic_ptr_cast<tcScalarSpinPtr>(inpart.spinInfo());}
  // if the spin info object exists use it
  if(inspin)
    {inspin->decayed(true);}
  else if(inpart.spinInfo())
    {throw DecayIntegratorError() << "Wrong type of spin info for th incoming particle"
				  << " in EtaPiPiGammaDecayer::me2()" 
				  << Exception::abortnow;}
  else
    {
      SpinPtr newspin=new_ptr(ScalarSpinInfo(inpart.momentum(),true));
      inspin = dynamic_ptr_cast<tcScalarSpinPtr>(newspin);
      inspin->decayed(true);
      const_ptr_cast<tPPtr>(&inpart)->spinInfo(newspin);
    }
  // create the spin info pointers of the outgoing particles
  VectorSpinPtr vspin;
  if(vertex)
    {
      decay[0]->spinInfo(new_ptr(ScalarSpinInfo(decay[0]->momentum(),true)));
      decay[1]->spinInfo(new_ptr(ScalarSpinInfo(decay[1]->momentum(),true)));
      vspin=new_ptr(VectorSpinInfo(decay[2]->momentum(),true));
      decay[2]->spinInfo(vspin);
    }
  // prefactor for the matrix element
  Complex pre=_coupling[imode()]*2.*sqrt(2.)/(_fpi*_fpi*_fpi);
  Lorentz5Momentum ppipi=decay[0]->momentum()+decay[1]->momentum();ppipi.rescaleMass();
  Energy q=ppipi.mass();
  Energy2 q2=q*q;
  Complex ii(0.,1.);
  // first VMD option
  if(_option[imode()]==0)
    {
      Energy pcm = Kinematics::pstarTwoBodyDecay(q,_mpi,_mpi);
      Complex resfact = q2/(_mrho*_mrho-q2-ii*_mrho*pcm*pcm*pcm*_rhoconst/q2);
      pre*=(1.+1.5*resfact);
    }
  // second VMD option
  else if(_option[imode()]==1)
    {
      Energy pcm = Kinematics::pstarTwoBodyDecay(q,_mpi,_mpi);
      Complex resfact = q2/(_mrho*_mrho-q2-ii*pcm*pcm*pcm*_rhoconst/q);
      pre*=(1.+1.5*resfact);
    }
  // analytic omne function
  else if(_option[imode()]==2)
    {pre*=(1.-_cconst+_cconst*(1.+_aconst*q2)/analyticOmnes(q2));}
  // experimental omnes function
  else if(_option[imode()]==3)
    {pre*=(1.-_cconst+_cconst*(1.+_aconst*q2)/experimentalOmnes(q2));}
  // wave function for the vector
  VectorWaveFunction vwave=VectorWaveFunction(decay[2]->momentum(),
					      decay[2]->dataPtr(),outgoing);
  LorentzPolarizationVector ptemp;
  //pre=_coupling[imode()]*2.*sqrt(2.)/(_fpi*_fpi*_fpi);
  LorentzPolarizationVector epstemp=pre*EpsFunction::product(decay[0]->momentum(),
							     decay[1]->momentum(),
							     decay[2]->momentum());
  // compute the matrix element
  DecayMatrixElement newME(1,1,1,3);
  vector<int> ispin(4,0);
  for(int ix=-1;ix<2;++ix)
    {
      // compute the polarization vector for this helicity
      if(ix!=0){vwave.reset(ix);ptemp=vwave.Wave();}
      else{ptemp=LorentzPolarizationVector();}
      if(vertex){vspin->setBasisState(ix,ptemp);}
      // and the matrix element
      ispin[3]=ix;
      newME(ispin)=epstemp*ptemp;
    }
  // contract the whole thing
  ME(newME);
  RhoDMatrix rhoin=RhoDMatrix(1);rhoin.average();
  double me=newME.contract(rhoin).real();
  return me;
}
 
double EtaPiPiGammaDecayer::threeBodyMatrixElement(int imodeb,Energy2 q2, Energy2 s3,
						   Energy2 s2,Energy2 s1,
						   Energy m1,Energy m2,Energy m3)
{
  Complex pre=_coupling[imodeb]*2.*sqrt(2.)/(_fpi*_fpi*_fpi);
  Energy q=sqrt(s3);
  Complex ii(0.,1.);
  // first VMD option
  if(_option[imodeb]==0)
    {
      Energy pcm = Kinematics::pstarTwoBodyDecay(q,_mpi,_mpi);
      Complex resfact = s3/(_mrho*_mrho-s3-ii*_mrho*pcm*pcm*pcm*_rhoconst/s3);
      pre*=(1.+1.5*resfact);
    }
  // second VMD option
  else if(_option[imodeb]==1)
    {
      Energy pcm = Kinematics::pstarTwoBodyDecay(q,_mpi,_mpi);
      Complex resfact = s3/(_mrho*_mrho-s3-ii*pcm*pcm*pcm*_rhoconst/q);
      pre*=(1.+1.5*resfact);
    }
  // analytic omne function
  else if(_option[imodeb]==2)
    {pre*=(1.-_cconst+_cconst*(1.+_aconst*s3)/analyticOmnes(s3));}
  // experimental omnes function
  else if(_option[imodeb]==3)
    {pre*=(1.-_cconst+_cconst*(1.+_aconst*s3)/experimentalOmnes(s3));}
  double factor=(pre*conj(pre)).real();
  Energy mpi2=_mpi*_mpi;
  return factor*((-mpi2*(-2*mpi2+s1+s2)*(-2*mpi2+s1+s2)+(mpi2-s1)*(mpi2-s2)*s3)/4.);
}

WidthCalculatorBasePtr 
EtaPiPiGammaDecayer::threeBodyMEIntegrator(const DecayMode & dm) const
{
  // workout which mode we are doing
  int id=dm.parent()->id(),imode=1;
  if(id==ParticleID::eta){imode=0;}
  // construct the integrator
  vector<double> inweights(1,1.);
  vector<double> inmass(1,getParticleData(ParticleID::rho0)->mass());
  vector<double> inwidth(1,getParticleData(ParticleID::rho0)->width());
  vector<int> intype(1,1);
  tcDecayIntegratorPtr decayer=this;
  WidthCalculatorBasePtr 
    output(new_ptr(ThreeBodyAllOnCalculator(inweights,intype,inmass,inwidth,
					    const_ptr_cast<tDecayIntegratorPtr>(decayer),
					    imode,_mpi,_mpi,0.)));
  return output;
}
}


namespace Herwig{
using namespace Genfun;

FUNCTION_OBJECT_IMP(OmnesFunction)
  
OmnesFunction::OmnesFunction(const OmnesFunction & right) 
{  }

OmnesFunction::OmnesFunction(Interpolator * in,Energy2 eps)
 {
   _interpolator=in;
   _precision=eps;
 }

// destructor
OmnesFunction::~OmnesFunction() {}

  void OmnesFunction::setScale(Energy2 in){_s=in;}
   
double OmnesFunction::operator ()(double xpoint) const
{
  double output;
  Energy q=sqrt(xpoint);
  if(abs(xpoint-_s)>_precision)
    {output= (*_interpolator)(q)/xpoint/(xpoint-_s);}
  else
    {output=0.;}
  return output;
}

}
