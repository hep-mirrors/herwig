// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorMeson3PionDecayer class.
//

#include "VectorMeson3PionDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "VectorMeson3PionDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig++/Helicity/EpsFunction.h"
#include "ThePEG/Helicity/ScalarSpinInfo.h"
#include "Herwig++/Decay/ThreeBodyIntegrator.h"
#include "Herwig++/PDT/ThreeBodyAllOnCalculator.h"

namespace Herwig {
using namespace ThePEG;
using ThePEG::Helicity::ScalarSpinInfo;

VectorMeson3PionDecayer::~VectorMeson3PionDecayer() {}

bool VectorMeson3PionDecayer::accept(const DecayMode & dm) const {
  // is this mode allowed
  bool allowed=false;
  // must be two outgoing particles
  if(dm.products().size()!=3){return allowed;}
  // check the ids of the inital particles
  int id0=dm.parent()->id();
  for(unsigned int ix=0;ix<_incoming.size();++ix)
    {if(int(_incoming[ix])==id0){allowed=true;}}
  if(!allowed){return allowed;}
  // check the id's of the outgoing particles
  int npi0=0,npip=0,npim=0,id;
  ParticleMSet::const_iterator pit = dm.products().begin();
  for(;pit!=dm.products().end();++pit)
    {
      id = (*pit)->id();
      if(id==ParticleID::pi0){++npi0;}
      else if(id==ParticleID::piplus){++npip;}
      else if(id==ParticleID::piminus){++npim;}
    }
  if(npi0==1&&npip==1&&npim==1){allowed=true;}
  return allowed;
}

ParticleVector VectorMeson3PionDecayer::decay(const DecayMode & dm,
				  const Particle & parent) const {
  // workout which mode we are doing
  int imode=-1;
  int id=parent.id();
  unsigned int ix=0;
  do
    {
      if(_incoming[ix]==id){imode=ix;}
      ++ix;
    }
  while(imode<0&&ix<_incoming.size());
  // perform the decay
  bool cc=false;
  return generate(true,cc,imode,parent);
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
    ("The \\classname{VectorMeson3PionDecayer} class is designed for the decay "
     "of I=0 vector mesons to three pions via a current taking into account the "
     "rho and a possible direct term");
  
  static ParVector<VectorMeson3PionDecayer,double> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &VectorMeson3PionDecayer::_incoming,
     0, 0, 0, -10000, 10000, false, false, true);
  
  static ParVector<VectorMeson3PionDecayer,double> interfaceCoupling
    ("Coupling",
     "The overall coupling for the decay, this is the coupling of the decaying "
     "particle to the lowest lying rho multiplet.",
     &VectorMeson3PionDecayer::_coupling,
     0, 0, 0, -10000, 10000, false, false, true);
  
  static ParVector<VectorMeson3PionDecayer,double> interfaceDirectCoupling
    ("DirectCoupling",
     "The magnitude of the coupling of the direct term with respect to the "
     "coupling of the lowest lying rho multiplet",
     &VectorMeson3PionDecayer::_directcoupling,
     0, 0, 0, -10000, 10000, false, false, true);
  
  static ParVector<VectorMeson3PionDecayer,double> interfaceRho2Coupling
    ("Rho2Coupling",
     "The magnitude of the coupling of the second rho multiplet with respect to the "
     "lowest lying  multiplet",
     &VectorMeson3PionDecayer::_rho2coupling,
     0, 0, 0, -10000, 10000, false, false, true);
  
  static ParVector<VectorMeson3PionDecayer,double> interfaceRho3Coupling
    ("Rho3Coupling",
     "The magntiude of the coupling of the third rho multiplet with respect to the "
     "lowest lying multiplet",
     &VectorMeson3PionDecayer::_rho3coupling,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<VectorMeson3PionDecayer,double> interfaceDirectPhase
    ("DirectPhase",
     "The phase of the coupling of the direct term with respect to the "
     "coupling of the lowest lying rho multiplet",
     &VectorMeson3PionDecayer::_directcoupling,
     0, 0, 0, -10000, 10000, false, false, true);
  
  static ParVector<VectorMeson3PionDecayer,double> interfaceRho2Phase
    ("Rho2Phase",
     "The phasee of the coupling of the second rho multiplet with respect to the "
     "lowest lying  multiplet",
     &VectorMeson3PionDecayer::_rho2coupling,
     0, 0, 0, -10000, 10000, false, false, true);
  
  static ParVector<VectorMeson3PionDecayer,double> interfaceRho3Phase
    ("Rho3Phase",
     "The phase of the coupling of the third rho multiplet with respect to the "
     "lowest lying multiplet",
     &VectorMeson3PionDecayer::_rho3coupling,
     0, 0, 0, -10000, 10000, false, false, true);
  
  static ParVector<VectorMeson3PionDecayer,double> interfaceMaxWgt
    ("MaxWeight",
     "The maximum weight for the integration of the channel",
     &VectorMeson3PionDecayer::_maxwgt,
     0, 0, 0, -10000, 10000, false, false, true);
  
  static ParVector<VectorMeson3PionDecayer,double> interfaceRho1Wgt
    ("Rho1Weight",
     "The weight for the lowest lying rho multiplet's in the integration",
     &VectorMeson3PionDecayer::_rho1wgt,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<VectorMeson3PionDecayer,double> interfaceRho2Wgt
    ("Rho2Weight",
     "The weight for the second rho multiplet's in the integration",
     &VectorMeson3PionDecayer::_rho2wgt,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<VectorMeson3PionDecayer,double> interfaceRho3Wgt
    ("Rho3Weight",
     "The weight for the third rho multiplet's in the integration",
     &VectorMeson3PionDecayer::_rho3wgt,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<VectorMeson3PionDecayer,double> interfaceRho1Mass
    ("Rho1Mass",
     "The mass of the lowest lying rho multiplet",
     &VectorMeson3PionDecayer::_rho1mass,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<VectorMeson3PionDecayer,double> interfaceRho2Mass
    ("Rho2Mass",
     "The mass of the second rho multiplet",
     &VectorMeson3PionDecayer::_rho2mass,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<VectorMeson3PionDecayer,double> interfaceRho3Mass
    ("Rho3Mass",
     "The mass of the third rho multiplet",
     &VectorMeson3PionDecayer::_rho3mass,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<VectorMeson3PionDecayer,double> interfaceRho1Width
    ("Rho1Width",
     "The width of the lowest lying rho multiplet",
     &VectorMeson3PionDecayer::_rho1width,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<VectorMeson3PionDecayer,double> interfaceRho2Width
    ("Rho2Width",
     "The width of the second rho multiplet",
     &VectorMeson3PionDecayer::_rho2width,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<VectorMeson3PionDecayer,double> interfaceRho3Width
    ("Rho3Width",
     "The width of the third rho multiplet",
     &VectorMeson3PionDecayer::_rho3width,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<VectorMeson3PionDecayer,bool> interfaceDefaultParam
    ("DefaultParameters",
     "If less than zero the default rho masses and widths are used, otherwise "
     "the values specified in the rhomass and width arrays are used,",
     &VectorMeson3PionDecayer::_defaultmass,
     0, 0, 0, -10000, 10000, false, false, true);

}

// the hadronic currents   
vector<LorentzPolarizationVector> 
VectorMeson3PionDecayer::decayCurrent(const bool vertex,const int ichan, 
				      const Particle & inpart,
				      const ParticleVector & decay) const
{
  // create the spin information if needed
  if(vertex)
    {
      for(unsigned int ix=0;ix<decay.size();++ix)
	{
	  SpinPtr temp=new_ptr(ScalarSpinInfo(decay[ix]->momentum(),true));
	  decay[ix]->spinInfo(temp);
	}
    }
  // identify the three pions
  int ipi0=-1,ipim=-1,ipip=-1,id;
  for(unsigned int ix=0;ix<decay.size();++ix)
    {
      id=decay[ix]->id();
      if(id==ParticleID::pi0){ipi0=ix;}
      else if(id==ParticleID::piplus){ipip=ix;}
      else if(id==ParticleID::piminus){ipim=ix;}
    }
  // work out the prefactor
  Complex pre=0,resfact,ii(0.,1.);
  if(ichan<0){pre=_ccoupling[imode()][3];}
  Energy pcm;
  // work out the direct invariant masses needed
  Lorentz5Momentum temp=decay[ipip]->momentum()+decay[ipim]->momentum();
  temp.rescaleMass();
  Energy mrho0 = temp.mass();
  temp = decay[ipip]->momentum()+decay[ipi0]->momentum();
  temp.rescaleMass();
  Energy mrhop=temp.mass();
  temp = decay[ipim]->momentum()+decay[ipi0]->momentum();
  temp.rescaleMass();
  Energy mrhom=temp.mass();
  // contribution of the resonances
  int ichannow=-3;
  for(unsigned int ix=0;ix<3;++ix)
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
	      resfact = 1./(mrho0*mrho0-_rhomass2[imode()][ix]
			    +ii*pcm*pcm*pcm*_rho0const[imode()][ix]/mrho0);
	      // rho+ contribution
	      pcm = Kinematics::pstarTwoBodyDecay(mrhop,_mpic,_mpi0);
	      resfact+= 1./(mrhop*mrhop-_rhomass2[imode()][ix]
			    +ii*pcm*pcm*pcm*_rhocconst[imode()][ix]/mrhop);
	      // rho- contribution
	      pcm = Kinematics::pstarTwoBodyDecay(mrhom,_mpic,_mpi0);
	      resfact+= 1./(mrhom*mrhom-_rhomass2[imode()][ix]
			    +ii*pcm*pcm*pcm*_rhocconst[imode()][ix]/mrhom);
	      resfact*=_ccoupling[imode()][ix];
	      // add the contribution
	    }
	  else if(ichan==ichannow)
	    {
	      pcm = Kinematics::pstarTwoBodyDecay(mrho0,_mpic,_mpic);
	      resfact = 1./(mrho0*mrho0-_rhomass2[imode()][ix]
			    +ii*pcm*pcm*pcm*_rho0const[imode()][ix]/mrho0);
	      resfact*=_ccoupling[imode()][ix];
	    }
	  else if(ichan==ichannow+1)
	    {
	      pcm = Kinematics::pstarTwoBodyDecay(mrhop,_mpic,_mpi0);
	      resfact+= 1./(mrhop*mrhop-_rhomass2[imode()][ix]
			    +ii*pcm*pcm*pcm*_rhocconst[imode()][ix]/mrhop);
	      resfact*=_ccoupling[imode()][ix];
	    }
	  else if(ichan==ichannow+2)
	    {
	      pcm = Kinematics::pstarTwoBodyDecay(mrhom,_mpic,_mpi0);
	      resfact+= 1./(mrhom*mrhom-_rhomass2[imode()][ix]
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
  LorentzPolarizationVector pol=Helicity::EpsFunction::product(decay[ipi0]->momentum(),
							       decay[ipip]->momentum(),
							       decay[ipim]->momentum());
  pol *=pre;
  return vector<LorentzPolarizationVector>(1,pol);
}

double VectorMeson3PionDecayer::threeBodyMatrixElement(int imode,Energy2 q2, Energy2 s3,
						       Energy2 s2,Energy2 s1,
						       Energy m1,Energy m2,Energy m3)
{
  Lorentz5Momentum p1,p2,p3; Energy2 ee1,ee2,ee3;Energy pp1,pp2,pp3;
  Energy q=sqrt(q2);
  Energy2 mpi2c=_mpic*_mpic;
  Energy2 mpi20=_mpi0*_mpi0;
  p1.setE(0.5*(q2+mpi20-s1)/q); ee1=p1.e()*p1.e(); pp1=sqrt(ee1-mpi20);
  p2.setE(0.5*(q2+mpi2c-s2)/q); ee2=p2.e()*p2.e(); pp2=sqrt(ee2-mpi2c);
  p3.setE(0.5*(q2+mpi2c-s3)/q); ee3=p3.e()*p3.e(); pp3=sqrt(ee3-mpi2c);
  // take momentum of 1 parallel to z axis
  p1.setPx(0.);p1.setPy(0.);p1.setPz(pp1);
  // construct 2 
  double cos2 = 0.5*(ee1+ee2-ee3-mpi20)/pp1/pp2;
  p2.setPx(pp2*sqrt(1.-cos2*cos2)); p2.setPy(0.); p2.setPz(-pp2*cos2);
  // construct 3
  double cos3 = 0.5*(ee1-ee2+ee3-mpi20)/pp1/pp3;
  p3.setPx(-pp3*sqrt(1.-cos3*cos3)); p3.setPy(0.); p3.setPz(-pp3*cos3); 
  // compute the prefactor
  Complex pre=_ccoupling[imode][3],resfact,ii(0.,1.);
  // rho0 contribution
  Energy pcm,mrho1=sqrt(s1),mrho2=sqrt(s2),mrho3=sqrt(s3);
  for(unsigned int ix=0;ix<3;++ix)
    {
      // rho0 contribution
      pcm = Kinematics::pstarTwoBodyDecay(mrho1,_mpic,_mpic);
      resfact = 1./(mrho1*mrho1-_rhomass2[imode][ix]
		    +ii*pcm*pcm*pcm*_rho0const[imode][ix]/mrho1);
      // rho+ contribution
      pcm = Kinematics::pstarTwoBodyDecay(mrho2,_mpic,_mpi0);
      resfact+= 1./(mrho2*mrho3-_rhomass2[imode][ix]
		    +ii*pcm*pcm*pcm*_rhocconst[imode][ix]/mrho2);
      // rho- contribution
      pcm = Kinematics::pstarTwoBodyDecay(mrho3,_mpic,_mpi0);
      resfact+= 1./(mrho3*mrho3-_rhomass2[imode][ix]
		    +ii*pcm*pcm*pcm*_rhocconst[imode][ix]/mrho3);
      resfact*=_ccoupling[imode][ix];
      // add the contribution
      pre+=resfact;
    }
  pre*=_coupling[imode];
  ThePEG::Helicity::LorentzPolarizationVector current=pre*
    Helicity::EpsFunction::product(p1,p2,p3);
  Complex temp=current*(current.conjugate());
  return -temp.real()/3.;
} 

WidthCalculatorBasePtr 
VectorMeson3PionDecayer::threeBodyMEIntegrator(const DecayMode & dm) const
{
  // workout which mode we are doing
  int imode=-1;
  int id=dm.parent()->id();
  unsigned int ix=0;
  do
    {
      if(_incoming[ix]==id){imode=ix;}
      ++ix;
    }
  while(imode<0&&ix<_incoming.size());
  // construct the integrator
  vector<double> inweights;inweights.push_back(1./3.);inweights.push_back(1./3.);
  inweights.push_back(1./3.);
  vector<int> intype;intype.push_back(1);intype.push_back(2);intype.push_back(3);
  Energy mrho=getParticleData(ParticleID::rhoplus)->mass();
  Energy wrho=getParticleData(ParticleID::rhoplus)->width();
  vector<double> inmass;inmass.push_back(mrho);
  inmass.push_back(mrho);inmass.push_back(mrho);
  vector<double> inwidth;inwidth.push_back(wrho);
  inwidth.push_back(wrho);inwidth.push_back(wrho);
  tcDecayIntegratorPtr decayer=this;
  WidthCalculatorBasePtr output(
    new_ptr(ThreeBodyAllOnCalculator(inweights,intype,inmass,inwidth,
				     const_ptr_cast<tDecayIntegratorPtr>(decayer),
				     imode,_mpi0,_mpic,_mpic)));
  return output;
}
}

