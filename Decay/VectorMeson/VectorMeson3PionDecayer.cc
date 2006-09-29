// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorMeson3PionDecayer class.
//

#include "VectorMeson3PionDecayer.h"
#include "Herwig++/Utilities/Kinematics.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "VectorMeson3PionDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig++/Helicity/EpsFunction.h"
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig++/Decay/ThreeBodyIntegrator.h"
#include "Herwig++/PDT/ThreeBodyAllOnCalculator.h"

namespace Herwig {
using namespace ThePEG;
using Herwig::Helicity::outgoing;
using Herwig::Helicity::incoming;
using Herwig::Helicity::EpsFunction;
using Herwig::Helicity::ScalarWaveFunction;
using Herwig::Helicity::VectorWaveFunction;
using ThePEG::Helicity::RhoDMatrix;
using ThePEG::Helicity::LorentzPolarizationVector;

VectorMeson3PionDecayer::VectorMeson3PionDecayer() {
  // generation of intermediates
  generateIntermediates(true);
  // omega decay
  _incoming.push_back(223);
  _coupling.push_back(178.71/GeV);
  _directcoupling.push_back(0.);_directphase.push_back(0.);
  _rho2coupling.push_back(0.);_rho2phase.push_back(0.);
  _rho3coupling.push_back(0.);_rho3phase.push_back(0.);
  _maxwgt.push_back(4.52975);
  _rho1wgt.push_back( 1.0);
  _rho2wgt.push_back(-1.0); 
  _rho3wgt.push_back(-1.0);
  _rho1mass.push_back(0.7758*GeV); 
  _rho2mass.push_back(1.4650*GeV);
  _rho3mass.push_back(1.7000*GeV); 
  _rho1width.push_back(0.1503*GeV); 
  _rho2width.push_back(0.3100*GeV);
  _rho3width.push_back(0.2400*GeV);
  _defaultmass.push_back(true);
  // phi decay
  _incoming.push_back(333);
  _coupling.push_back(9.424029/GeV);
  _directcoupling.push_back(0.78);_directphase.push_back(2.47);
  _rho2coupling.push_back(0.);_rho2phase.push_back(0.);
  _rho3coupling.push_back(0.);_rho3phase.push_back(0.);
  _maxwgt.push_back(5.62103);
  _rho1wgt.push_back( 1.0);
  _rho2wgt.push_back(-1.0); 
  _rho3wgt.push_back(-1.0);
  _rho1mass.push_back(0.7758*GeV); 
  _rho2mass.push_back(1.4500*GeV);
  _rho3mass.push_back(1.7000*GeV); 
  _rho1width.push_back(0.1439*GeV); 
  _rho2width.push_back(0.3100*GeV);
  _rho3width.push_back(0.2400*GeV);
  _defaultmass.push_back(false);
  _mpi0=0.;_mpic=0.;
  // initial size of the arrays
  _initsize=_coupling.size();
}

void VectorMeson3PionDecayer::doinit() throw(InitException) {
  DecayIntegrator::doinit();
  // check the consistence of the decay modes
  unsigned int isize=_incoming.size();
  if(isize!=_coupling.size()       || 
     isize!=_directcoupling.size() || isize!=_directphase.size()  ||
     isize!=_rho2coupling.size()   || isize!=_rho2phase.size()    ||
     isize!=_rho3coupling.size()   || isize!=_rho3phase.size()    ||
     isize!=_maxwgt.size()         || isize!=_rho1wgt.size()      ||
     isize!=_rho2wgt.size()        || isize!=_rho3wgt.size()      ||
     isize!=_rho1mass.size()       || isize!=_rho2mass.size()     ||
     isize!=_rho3mass.size()       || isize!=_rho1width.size()    ||
     isize!=_rho2width.size()      || isize!=_rho3width.size())
    {throw InitException() << "Inconsistent parameters in " 
			   << "VectorMeson3PionDecayer::doinit()" 
			   << Exception::abortnow;}
  // calculate the parameters 
  // set the external particles
  PDVector extpart(4);
  extpart[1]=getParticleData(ParticleID::pi0);
  extpart[2]=getParticleData(ParticleID::piplus);
  extpart[3]=getParticleData(ParticleID::piminus);
  // pointer to the different rho resonances
  // the rho0 resonances
  tPDPtr rho0[3]={getParticleData(113),getParticleData(100113),getParticleData(30113)};
  // the charged rho resonance
  tPDPtr rhom[3]={getParticleData(-213),getParticleData(-100213),
		  getParticleData(-30213)};
  tPDPtr rhop[3]={getParticleData(213),getParticleData(100213),getParticleData(30213)};
  // create the integration channels for the decay
  DecayPhaseSpaceModePtr mode;
  DecayPhaseSpaceChannelPtr newchannel;
  unsigned int iy,iz;
  for(unsigned int ix=0;ix<_incoming.size();++ix)
    {
      extpart[0]=getParticleData(int(_incoming[ix]));
      mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
      // decide which rho resonances to add
      double  temp[3] = {_rho1wgt[ix]  ,_rho2wgt[ix]  ,_rho3wgt[ix]  };
      Energy  mass[3] = {_rho1mass[ix] ,_rho2mass[ix] ,_rho3mass[ix] };
      Energy width[3] = {_rho1width[ix],_rho2width[ix],_rho3width[ix]};
      vector<double> wgt;
      // set the mass parameters to the default if needed
      if(_defaultmass[ix])
	{
	  _rho1mass[ix] = rhom[0]->mass(); _rho1width[ix] = rhom[0]->width();
	  _rho2mass[ix] = rhom[1]->mass(); _rho2width[ix] = rhom[1]->width();
	  _rho3mass[ix] = rhom[2]->mass(); _rho3width[ix] = rhom[2]->width();
	}
      double sumwgt(0);
      for(iy=0;iy<3;++iy){if(temp[iy]>0){sumwgt+=temp[iy];}}
      for(iy=0;iy<3;++iy)
	{
	  if(temp[iy]>0)
	    {
	      // set the weights for the channels
	      for(iz=0;iz<3;++iz){wgt.push_back(temp[iy]/3./sumwgt);}
	      // rho0 channel
	      newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
	      newchannel->addIntermediate(extpart[0],0,0.0,-1,1);
	      newchannel->addIntermediate(rho0[iy]  ,0,0.0, 2,3);
	      mode->addChannel(newchannel);
	      if(!_defaultmass[ix])
		{resetIntermediate(rho0[iy],mass[iy],width[iy]);}
	      // rho+ channel
	      newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
	      newchannel->addIntermediate(extpart[0],0,0.0,-1,3);
	      newchannel->addIntermediate(rhop[iy]  ,0,0.0, 1,2);
	      mode->addChannel(newchannel);
	      if(!_defaultmass[ix])
		{resetIntermediate(rhop[iy],mass[iy],width[iy]);}
	      // rho- channel
	      newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
	      newchannel->addIntermediate(extpart[0],0,0.0,-1,2);
	      newchannel->addIntermediate(rhom[iy]  ,0,0.0, 1,3);
	      mode->addChannel(newchannel);
	      if(!_defaultmass[ix])
		{mode->resetIntermediate(rhom[iy],mass[iy],width[iy]);}
	      addMode(mode,_maxwgt[ix],wgt);
	    }
	}
    }
  // work out the masses and constants for the running widths
  Energy pcm;
  _mpi0=getParticleData(ParticleID::pi0)->mass();
  _mpic=getParticleData(ParticleID::piplus)->mass();
  Complex ii(0.,1.);
  _rhomass.resize(_incoming.size()); _rhomass2.resize(_incoming.size());
  _rho0const.resize(_incoming.size()); _rhocconst.resize( _incoming.size());
  _ccoupling.resize(_incoming.size());
  for(unsigned int ix=0;ix<_incoming.size();++ix)
    {
      // set the masses of the rho resonances
      _rhomass[ix].push_back(_rho1mass[ix]);
      _rhomass[ix].push_back(_rho2mass[ix]);
      _rhomass[ix].push_back(_rho3mass[ix]);
      _rhomass2[ix].push_back(_rho1mass[ix]*_rho1mass[ix]);
      _rhomass2[ix].push_back(_rho2mass[ix]*_rho2mass[ix]);
      _rhomass2[ix].push_back(_rho3mass[ix]*_rho3mass[ix]);
      // set up the constants for the running width
      pcm=Kinematics::pstarTwoBodyDecay(_rho1mass[ix],_mpic,_mpic);
      _rho0const[ix].push_back(_rho1mass[ix]*_rho1mass[ix]*_rho1width[ix]/(pcm*pcm*pcm));
      pcm=Kinematics::pstarTwoBodyDecay(_rho2mass[ix],_mpic,_mpic);
      _rho0const[ix].push_back(_rho2mass[ix]*_rho2mass[ix]*_rho2width[ix]/(pcm*pcm*pcm));
      pcm=Kinematics::pstarTwoBodyDecay(_rho3mass[ix],_mpic,_mpic);
      _rho0const[ix].push_back(_rho3mass[ix]*_rho3mass[ix]*_rho3width[ix]/(pcm*pcm*pcm));
      pcm=Kinematics::pstarTwoBodyDecay(_rho1mass[ix],_mpi0,_mpic);
      _rhocconst[ix].push_back(_rho1mass[ix]*_rho1mass[ix]*_rho1width[ix]/(pcm*pcm*pcm));
      pcm=Kinematics::pstarTwoBodyDecay(_rho2mass[ix],_mpi0,_mpic);
      _rhocconst[ix].push_back(_rho2mass[ix]*_rho2mass[ix]*_rho2width[ix]/(pcm*pcm*pcm));
      pcm=Kinematics::pstarTwoBodyDecay(_rho3mass[ix],_mpi0,_mpic);
      _rhocconst[ix].push_back(_rho3mass[ix]*_rho3mass[ix]*_rho3width[ix]/(pcm*pcm*pcm));
      // set the complex coupling constants
      _ccoupling[ix].push_back(1./_rhomass2[ix][0]);
      _ccoupling[ix].push_back( _rho2coupling[ix]/_rhomass2[ix][0]*
				(cos(_rho2phase[ix])+ii*sin(_rho2phase[ix])));
      _ccoupling[ix].push_back( _rho3coupling[ix]/_rhomass2[ix][0]*
				(cos(_rho3phase[ix])+ii*sin(_rho3phase[ix])));
      _ccoupling[ix].push_back(_directcoupling[ix]/_rhomass2[ix][0]*
			       (cos(_directphase[ix])+ii*sin(_directphase[ix])));
    }
}

VectorMeson3PionDecayer::~VectorMeson3PionDecayer() {}

int VectorMeson3PionDecayer::modeNumber(bool & cc,const DecayMode & dm) const
{
  cc=false;
  int imode(-1),id;
  // must be three outgoing particles
  if(dm.products().size()!=3){return imode;}
  // check the id's of the outgoing particles
  unsigned int npi0(0),npip(0),npim(0);
  ParticleMSet::const_iterator pit(dm.products().begin());
  for(;pit!=dm.products().end();++pit)
    {
      id = (*pit)->id();
      if(id==ParticleID::pi0){++npi0;}
      else if(id==ParticleID::piplus){++npip;}
      else if(id==ParticleID::piminus){++npim;}
    }
  if(!(npi0==1&&npip==1&&npim==1)){return imode;}
  unsigned int ix(0);
  id=dm.parent()->id();
  do{if(_incoming[ix]==id){imode=ix;}++ix;}
  while(imode<0&&ix<_incoming.size());
  return imode;
}

void VectorMeson3PionDecayer::persistentOutput(PersistentOStream & os) const {
  os << _incoming << _coupling << _directcoupling << _rho2coupling << _rho3coupling 
     << _directphase << _rho2phase << _rho3phase 
     << _maxwgt << _rho1wgt << _rho2wgt << _rho3wgt << _rho1mass << _rho2mass
     << _rho3mass << _rho1width << _rho2width << _rho3width << _defaultmass 
     << _rho0const << _rhocconst << _rhomass << _rhomass2 << _ccoupling
     << _mpi0 << _mpic;
}

void VectorMeson3PionDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _incoming >> _coupling >> _directcoupling >> _rho2coupling >> _rho3coupling 
     >> _directphase >> _rho2phase >> _rho3phase 
     >> _maxwgt >> _rho1wgt >> _rho2wgt >> _rho3wgt >> _rho1mass >> _rho2mass
     >> _rho3mass >> _rho1width >> _rho2width >> _rho3width >> _defaultmass
     >> _rho0const >> _rhocconst >> _rhomass >> _rhomass2 >> _ccoupling
     >> _mpi0 >> _mpic;
}

ClassDescription<VectorMeson3PionDecayer> VectorMeson3PionDecayer::initVectorMeson3PionDecayer;
// Definition of the static class description member.

void VectorMeson3PionDecayer::Init() {

  static ClassDocumentation<VectorMeson3PionDecayer> documentation
    ("The VectorMeson3PionDecayer class is designed for the decay "
     "of I=0 vector mesons to three pions via a current taking into account the "
     "rho and a possible direct term");
  
  static ParVector<VectorMeson3PionDecayer,double> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &VectorMeson3PionDecayer::_incoming,
     0, 0, 0, 0, 1000000, false, false, true);
  
  static ParVector<VectorMeson3PionDecayer,InvEnergy> interfaceCoupling
    ("Coupling",
     "The overall coupling for the decay, this is the coupling of the decaying "
     "particle to the lowest lying rho multiplet.",
     &VectorMeson3PionDecayer::_coupling,
     1./GeV, -1, 1./GeV,-1000./GeV, 1000./GeV, false, false, true);
  
  static ParVector<VectorMeson3PionDecayer,double> interfaceDirectCoupling
    ("DirectCoupling",
     "The magnitude of the coupling of the direct term with respect to the "
     "coupling of the lowest lying rho multiplet",
     &VectorMeson3PionDecayer::_directcoupling,
     0, 0, 0, -10., 10., false, false, true);
  
  static ParVector<VectorMeson3PionDecayer,double> interfaceRho2Coupling
    ("Rho2Coupling",
     "The magnitude of the coupling of the second rho multiplet with respect to the "
     "lowest lying  multiplet",
     &VectorMeson3PionDecayer::_rho2coupling,
     0, 0, 0, -10., 10., false, false, true);
  
  static ParVector<VectorMeson3PionDecayer,double> interfaceRho3Coupling
    ("Rho3Coupling",
     "The magntiude of the coupling of the third rho multiplet with respect to the "
     "lowest lying multiplet",
     &VectorMeson3PionDecayer::_rho3coupling,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<VectorMeson3PionDecayer,double> interfaceDirectPhase
    ("DirectPhase",
     "The phase of the coupling of the direct term with respect to the "
     "coupling of the lowest lying rho multiplet",
     &VectorMeson3PionDecayer::_directphase,
     0, 0, 0, -2.*pi, 2.*pi, false, false, true);
  
  static ParVector<VectorMeson3PionDecayer,double> interfaceRho2Phase
    ("Rho2Phase",
     "The phasee of the coupling of the second rho multiplet with respect to the "
     "lowest lying  multiplet",
     &VectorMeson3PionDecayer::_rho2phase,
     0, 0, 0, -2.*pi, 2.*pi, false, false, true);
  
  static ParVector<VectorMeson3PionDecayer,double> interfaceRho3Phase
    ("Rho3Phase",
     "The phase of the coupling of the third rho multiplet with respect to the "
     "lowest lying multiplet",
     &VectorMeson3PionDecayer::_rho3phase,
     0, 0, 0, -2.*pi, 2.*pi, false, false, true);
  
  static ParVector<VectorMeson3PionDecayer,double> interfaceMaxWgt
    ("MaxWeight",
     "The maximum weight for the integration of the channel",
     &VectorMeson3PionDecayer::_maxwgt,
     0, 0, 0, 0., 100., false, false, true);
  
  static ParVector<VectorMeson3PionDecayer,double> interfaceRho1Wgt
    ("Rho1Weight",
     "The weight for the lowest lying rho multiplet's in the integration",
     &VectorMeson3PionDecayer::_rho1wgt,
     0, 0, 0, 0, 1., false, false, true);

  static ParVector<VectorMeson3PionDecayer,double> interfaceRho2Wgt
    ("Rho2Weight",
     "The weight for the second rho multiplet's in the integration",
     &VectorMeson3PionDecayer::_rho2wgt,
     0, 0, 0, -2., 1., false, false, true);

  static ParVector<VectorMeson3PionDecayer,double> interfaceRho3Wgt
    ("Rho3Weight",
     "The weight for the third rho multiplet's in the integration",
     &VectorMeson3PionDecayer::_rho3wgt,
     0, 0, 0, -2., 1., false, false, true);

  static ParVector<VectorMeson3PionDecayer,double> interfaceRho1Mass
    ("Rho1Mass",
     "The mass of the lowest lying rho multiplet",
     &VectorMeson3PionDecayer::_rho1mass,
     GeV, -1, 0.77*GeV, 0.*GeV, 5.*GeV, false, false, true);

  static ParVector<VectorMeson3PionDecayer,double> interfaceRho2Mass
    ("Rho2Mass",
     "The mass of the second rho multiplet",
     &VectorMeson3PionDecayer::_rho2mass,
     GeV, -1, 0.77*GeV, 0.*GeV, 5.*GeV, false, false, true);

  static ParVector<VectorMeson3PionDecayer,double> interfaceRho3Mass
    ("Rho3Mass",
     "The mass of the third rho multiplet",
     &VectorMeson3PionDecayer::_rho3mass,
     GeV, -1, 0.77*GeV, 0.*GeV, 5.*GeV, false, false, true);

  static ParVector<VectorMeson3PionDecayer,double> interfaceRho1Width
    ("Rho1Width",
     "The width of the lowest lying rho multiplet",
     &VectorMeson3PionDecayer::_rho1width,
     GeV, -1, 0.15*GeV, 0.*GeV, 1.*GeV, false, false, true);

  static ParVector<VectorMeson3PionDecayer,double> interfaceRho2Width
    ("Rho2Width",
     "The width of the second rho multiplet",
     &VectorMeson3PionDecayer::_rho2width,
     GeV, -1, 0.15*GeV, 0.*GeV, 1.*GeV, false, false, true);

  static ParVector<VectorMeson3PionDecayer,double> interfaceRho3Width
    ("Rho3Width",
     "The width of the third rho multiplet",
     &VectorMeson3PionDecayer::_rho3width,
     GeV, -1, 0.15*GeV, 0.*GeV, 1.*GeV, false, false, true);

  static ParVector<VectorMeson3PionDecayer,bool> interfaceDefaultParam
    ("DefaultParameters",
     "If true the default rho masses and widths are used, otherwise "
     "the values specified in the rhomass and width arrays are used,",
     &VectorMeson3PionDecayer::_defaultmass,
     0, 0, 1, 0, 1, false, false, true);

}
double VectorMeson3PionDecayer::me2(bool vertex, const int ichan,
				    const Particle & inpart,
				    const ParticleVector & decay) const
{
  // wavefunctions for the decaying particle
  RhoDMatrix rhoin(PDT::Spin1);rhoin.average();
  vector<LorentzPolarizationVector> invec;
  VectorWaveFunction(invec,rhoin,const_ptr_cast<tPPtr>(&inpart),
		     incoming,true,false,vertex);
  // create the spin information for the decay products if needed
  unsigned int ix;
  if(vertex)
    {for(ix=0;ix<decay.size();++ix)
	// workaround for gcc 3.2.3 bug
	//ALB {ScalarWaveFunction(decay[ix],outgoing,true,vertex);}}
	{PPtr mytemp = decay[ix] ; ScalarWaveFunction(mytemp,outgoing,true,vertex);}}
  // compute the matrix element
  DecayMatrixElement newME(PDT::Spin1,PDT::Spin0,PDT::Spin0,PDT::Spin0);
  // work out the prefactor
  Complex pre(0.),resfact,ii(0.,1.);
  if(ichan<0){pre=_ccoupling[imode()][3];}
  Energy pcm;
  // work out the direct invariant masses needed
  Lorentz5Momentum temp(decay[1]->momentum()+decay[2]->momentum());temp.rescaleMass();
  Energy mrho0(temp.mass());
  temp = decay[1]->momentum()+decay[0]->momentum();temp.rescaleMass();
  Energy mrhop(temp.mass());
  temp = decay[2]->momentum()+decay[0]->momentum();temp.rescaleMass();
  Energy mrhom(temp.mass());
  // contribution of the resonances
  int ichannow(-3);
  for(ix=0;ix<3;++ix)
    {
      ichannow+=3;
      if((ix==0 && _rho1wgt[imode()]>0.) || 
	 (ix==1 && _rho2wgt[imode()]>0.) ||
	 (ix==2 && _rho3wgt[imode()]>0.))
	{
	  if(ichan<0)
	    {
	      // rho0 contribution
	      pcm = Kinematics::pstarTwoBodyDecay(mrho0,_mpic,_mpic);
	      resfact = _rhomass2[imode()][ix]/
		(mrho0*mrho0-_rhomass2[imode()][ix]
		 +ii*pcm*pcm*pcm*_rho0const[imode()][ix]/mrho0);
	      // rho+ contribution
	      pcm = Kinematics::pstarTwoBodyDecay(mrhop,_mpic,_mpi0);
	      resfact+= _rhomass2[imode()][ix]/
		(mrhop*mrhop-_rhomass2[imode()][ix]
		 +ii*pcm*pcm*pcm*_rhocconst[imode()][ix]/mrhop);
	      // rho- contribution
	      pcm = Kinematics::pstarTwoBodyDecay(mrhom,_mpic,_mpi0);
	      resfact+= _rhomass2[imode()][ix]/
		(mrhom*mrhom-_rhomass2[imode()][ix]
		 +ii*pcm*pcm*pcm*_rhocconst[imode()][ix]/mrhom);
	      resfact*=_ccoupling[imode()][ix];
	      // add the contribution
	    }
	  else if(ichan==ichannow)
	    {
	      pcm = Kinematics::pstarTwoBodyDecay(mrho0,_mpic,_mpic);
	      resfact = _rhomass2[imode()][ix]/
		(mrho0*mrho0-_rhomass2[imode()][ix]
		 +ii*pcm*pcm*pcm*_rho0const[imode()][ix]/mrho0);
	      resfact*=_ccoupling[imode()][ix];
	    }
	  else if(ichan==ichannow+1)
	    {
	      pcm = Kinematics::pstarTwoBodyDecay(mrhop,_mpic,_mpi0);
	      resfact+= _rhomass2[imode()][ix]/
		(mrhop*mrhop-_rhomass2[imode()][ix]
		 +ii*pcm*pcm*pcm*_rhocconst[imode()][ix]/mrhop);
	      resfact*=_ccoupling[imode()][ix];
	    }
	  else if(ichan==ichannow+2)
	    {
	      pcm = Kinematics::pstarTwoBodyDecay(mrhom,_mpic,_mpi0);
	      resfact+= _rhomass2[imode()][ix]/
		(mrhom*mrhom-_rhomass2[imode()][ix]
		 +ii*pcm*pcm*pcm*_rhocconst[imode()][ix]/mrhom);
	      resfact*=_ccoupling[imode()][ix];
	    }
	  pre+=resfact;
	  ichannow+=3;
	}
    }
  // overall coupling
  pre *=_coupling[imode()];
  // polarization vector piece
  LorentzPolarizationVector 
    scalar=pre*EpsFunction::product(decay[0]->momentum(),decay[1]->momentum(),
				    decay[2]->momentum());
  // compute the matrix element
  for(ix=0;ix<3;++ix)
    {newME(ix,0,0,0)=scalar*invec[ix];}
  ME(newME);
  // return the answer
  return newME.contract(rhoin).real();
}

double VectorMeson3PionDecayer::threeBodyMatrixElement(int imode,Energy2 q2, Energy2 s3,
						       Energy2 s2,Energy2 s1,
						       Energy,Energy,Energy)
{
  Lorentz5Momentum p1,p2,p3; Energy2 ee1,ee2,ee3;Energy pp1,pp2,pp3;
  Energy q(sqrt(q2));
  Energy2 mpi2c(_mpic*_mpic),mpi20(_mpi0*_mpi0);
  p1.setE(0.5*(q2+mpi20-s1)/q); ee1=p1.e()*p1.e(); pp1=sqrt(ee1-mpi20);
  p2.setE(0.5*(q2+mpi2c-s2)/q); ee2=p2.e()*p2.e(); pp2=sqrt(ee2-mpi2c);
  p3.setE(0.5*(q2+mpi2c-s3)/q); ee3=p3.e()*p3.e(); pp3=sqrt(ee3-mpi2c);
  // take momentum of 1 parallel to z axis
  p1.setPx(0.);p1.setPy(0.);p1.setPz(pp1);
  // construct 2 
  double cos2(0.5*(ee1+ee2-ee3-mpi20)/pp1/pp2);
  p2.setPx(pp2*sqrt(1.-cos2*cos2)); p2.setPy(0.); p2.setPz(-pp2*cos2);
  // construct 3
  double cos3(0.5*(ee1-ee2+ee3-mpi20)/pp1/pp3);
  p3.setPx(-pp3*sqrt(1.-cos3*cos3)); p3.setPy(0.); p3.setPz(-pp3*cos3); 
  // compute the prefactor
  Complex pre(_ccoupling[imode][3]),resfact,ii(0.,1.);
  // rho0 contribution
  Energy pcm,mrho1(sqrt(s1)),mrho2(sqrt(s2)),mrho3(sqrt(s3));
  for(unsigned int ix=0;ix<3;++ix)
    {
      // rho0 contribution
      pcm = Kinematics::pstarTwoBodyDecay(mrho1,_mpic,_mpic);
      resfact = _rhomass2[imode][ix]/(mrho1*mrho1-_rhomass2[imode][ix]
				      +ii*pcm*pcm*pcm*_rho0const[imode][ix]/mrho1);
      // rho+ contribution
      pcm = Kinematics::pstarTwoBodyDecay(mrho2,_mpic,_mpi0);
      resfact+= _rhomass2[imode][ix]/(mrho2*mrho3-_rhomass2[imode][ix]
				      +ii*pcm*pcm*pcm*_rhocconst[imode][ix]/mrho2);
      // rho- contribution
      pcm = Kinematics::pstarTwoBodyDecay(mrho3,_mpic,_mpi0);
      resfact+= _rhomass2[imode][ix]/(mrho3*mrho3-_rhomass2[imode][ix]
				      +ii*pcm*pcm*pcm*_rhocconst[imode][ix]/mrho3);
      resfact*=_ccoupling[imode][ix];
      // add the contribution
      pre+=resfact;
    }
  pre*=_coupling[imode];
  LorentzPolarizationVector current(pre*EpsFunction::product(p1,p2,p3));
  Complex temp(current*(current.conjugate()));
  return -temp.real()/3.;
} 

WidthCalculatorBasePtr 
VectorMeson3PionDecayer::threeBodyMEIntegrator(const DecayMode & dm) const
{
  // workout which mode we are doing
  int imode(-1),id(dm.parent()->id());
  unsigned int ix=0;
  do{if(_incoming[ix]==id){imode=ix;}++ix;}
  while(imode<0&&ix<_incoming.size());
  // construct the integrator
  vector<double> inweights(3,1./3.);
  vector<int> intype;intype.push_back(1);intype.push_back(2);intype.push_back(3);
  Energy mrho(getParticleData(ParticleID::rhoplus)->mass());
  Energy wrho(getParticleData(ParticleID::rhoplus)->width());
  vector<double> inmass(3,mrho);
  vector<double> inwidth(3,wrho);
  //tcDecayIntegratorPtr decayer(this);
  WidthCalculatorBasePtr output(
    new_ptr(ThreeBodyAllOnCalculator(inweights,intype,inmass,inwidth,
				     const_ptr_cast<tDecayIntegratorPtr>(this),
				     imode,_mpi0,_mpic,_mpic)));
  return output;
}

void VectorMeson3PionDecayer::dataBaseOutput(ofstream & output,
					     bool header) const
{
  if(header){output << "update decayers set parameters=\"";}
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  for(unsigned int ix=0;ix<_incoming.size();++ix)
    {
      if(ix<_initsize)
	{
	  output << "set " << fullName() << ":Incoming " 
		 << ix << " " << _incoming[ix] << endl;
	  output << "set " << fullName() << ":Coupling " 
		 << ix << " " << _coupling[ix]*GeV << endl;
	  output << "set " << fullName() << ":DirectCoupling " 
		 << ix << " " << _directcoupling[ix] << endl;
	  output << "set " << fullName() << ":Rho2Coupling " 
		 << ix << " " << _rho2coupling[ix] << endl;
	  output << "set " << fullName() << ":Rho3Coupling " 
		 << ix << " " << _rho3coupling[ix] << endl;
	  output << "set " << fullName() << ":DirectPhase " 
		 << ix << " " << _directphase[ix] << endl;
	  output << "set " << fullName() << ":Rho2Phase " 
		 << ix << " " << _rho2phase[ix] << endl;
	  output << "set " << fullName() << ":Rho3Phase " 
		 << ix << " " << _rho3phase[ix] << endl;
	  output << "set " << fullName() << ":MaxWeight " 
		 << ix << " " << _maxwgt[ix] << endl;
	  output << "set " << fullName() << ":Rho1Weight " 
		 << ix << " " << _rho1wgt[ix] << endl;
	  output << "set " << fullName() << ":Rho2Weight " 
		 << ix << " " << _rho2wgt[ix] << endl;
	  output << "set " << fullName() << ":Rho3Weight " 
		 << ix << " " << _rho3wgt[ix] << endl;
	  output << "set " << fullName() << ":Rho1Mass " 
		 << ix << " " << _rho1mass[ix]/GeV << endl;
	  output << "set " << fullName() << ":Rho2Mass " 
		 << ix << " " << _rho2mass[ix]/GeV<< endl;
	  output << "set " << fullName() << ":Rho3Mass " 
		 << ix << " " << _rho3mass[ix]/GeV<< endl;
	  output << "set " << fullName() << ":Rho1Width " 
		 << ix << " " << _rho1width[ix]/GeV << endl;
	  output << "set " << fullName() << ":Rho2Width " 
		 << ix << " " << _rho2width[ix]/GeV << endl;
	  output << "set " << fullName() << ":Rho3Width " 
		 << ix << " " << _rho3width[ix]/GeV << endl;
	  output << "set " << fullName() << ":DefaultParameters " 
		 << ix << " " << _defaultmass[ix] << endl;
	}
      else
	{
	  output << "insert " << fullName() << ":Incoming " 
		 << ix << " " << _incoming[ix] << endl;
	  output << "insert " << fullName() << ":Coupling " 
		 << ix << " " << _coupling[ix]*GeV << endl;
	  output << "insert " << fullName() << ":DirectCoupling " 
		 << ix << " " << _directcoupling[ix] << endl;
	  output << "insert " << fullName() << ":Rho2Coupling " 
		 << ix << " " << _rho2coupling[ix] << endl;
	  output << "insert " << fullName() << ":Rho3Coupling " 
		 << ix << " " << _rho3coupling[ix] << endl;
	  output << "insert " << fullName() << ":DirectPhase " 
		 << ix << " " << _directphase[ix] << endl;
	  output << "insert " << fullName() << ":Rho2Phase " 
		 << ix << " " << _rho2phase[ix] << endl;
	  output << "insert " << fullName() << ":Rho3Phase " 
		 << ix << " " << _rho3phase[ix] << endl;
	  output << "insert " << fullName() << ":MaxWeight " 
		 << ix << " " << _maxwgt[ix] << endl;
	  output << "insert " << fullName() << ":Rho1Weight " 
		 << ix << " " << _rho1wgt[ix] << endl;
	  output << "insert " << fullName() << ":Rho2Weight " 
		 << ix << " " << _rho2wgt[ix] << endl;
	  output << "insert " << fullName() << ":Rho3Weight " 
		 << ix << " " << _rho3wgt[ix] << endl;
	  output << "insert " << fullName() << ":Rho1Mass " 
		 << ix << " " << _rho1mass[ix]/GeV << endl;
	  output << "insert " << fullName() << ":Rho2Mass " 
		 << ix << " " << _rho2mass[ix]/GeV<< endl;
	  output << "insert " << fullName() << ":Rho3Mass " 
		 << ix << " " << _rho3mass[ix]/GeV<< endl;
	  output << "insert " << fullName() << ":Rho1Width " 
		 << ix << " " << _rho1width[ix]/GeV << endl;
	  output << "insert " << fullName() << ":Rho2Width " 
		 << ix << " " << _rho2width[ix]/GeV << endl;
	  output << "insert " << fullName() << ":Rho3Width " 
		 << ix << " " << _rho3width[ix]/GeV << endl;
	  output << "insert " << fullName() << ":DefaultParameters " 
		 << ix << " " << _defaultmass[ix] << endl;
	}
    }
  if(header){output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;}
}
}

