// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TauDecayerBase class.
//
//  Author: Peter Richardson
//

#include "TauDecayerBase.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig++/Helicity/Correlations/DecayVertex.h"
#include "ThePEG/Utilities/Timer.h"

namespace Herwig {

using namespace ThePEG;
using Helicity::SpinorWaveFunction;
using Helicity::SpinorBarWaveFunction;
using ThePEG::Helicity::LorentzPolarizationVector;
using ThePEG::Helicity::LorentzSpinor;
using ThePEG::Helicity::LorentzSpinorBar;
using ThePEG::Helicity::FermionSpinInfo;
using ThePEG::Helicity::tcFermionSpinPtr;
using ThePEG::Helicity::FermionSpinPtr;
using ThePEG::Helicity::SpinInfo;
using ThePEG::Helicity::tcSpinInfoPtr;
using ThePEG::Helicity::RhoDMatrix;
using ThePEG::Helicity::VertexPtr;
using ThePEG::Helicity::DiracRep;
using ThePEG::Helicity::HaberDRep;
using ThePEG::Helicity::HELASDRep;
using Helicity::DVertexPtr;
using Helicity::DecayVertex;
using Helicity::Direction;
using Helicity::incoming;
using Helicity::outgoing;

TauDecayerBase::~TauDecayerBase() {}
  
bool TauDecayerBase::accept(const DecayMode & dm) const {
  return false;
}
  
ParticleVector TauDecayerBase::decay(const DecayMode & dm,
				     const Particle & parent) const {
  ParticleVector children = dm.produceProducts();
  return children;
}
  
  
void TauDecayerBase::persistentOutput(PersistentOStream & os) const 
{os << _GF;}    
  
void TauDecayerBase::persistentInput(PersistentIStream & is, int) 
{is >> _GF;}
  
AbstractClassDescription<TauDecayerBase> TauDecayerBase::initTauDecayerBase;
// Definition of the static class description member.
  
void TauDecayerBase::Init() {
  
  static ClassDocumentation<TauDecayerBase> documentation
    ("The \\classname{TauDecayerBase} class is the base class for the"
     " implementation of tau decays in Herwig++. All actual decays should "
     "inherit from it and implement the calculation of the hadronic current.");
  
  static Parameter<TauDecayerBase,InvEnergy2> interfaceGFermi
    ("GFermi",
     "The Fermi coupling constant",
     &TauDecayerBase::_GF, 1./GeV2, 1.16639E-5/GeV2, -1.0e12*1./GeV2, 1.0e12*1./GeV2,
     false, false, false);
  
}

// the hadronic currents (this is a dummy)
vector<LorentzPolarizationVector> 
TauDecayerBase::hadronCurrent(bool vertex, const int imode, const int ichan,
			      const Particle & inpart,
			      const ParticleVector & decay) const
{vector<LorentzPolarizationVector> temp; return temp;}

// combine the currents to give the matrix element
double TauDecayerBase::me2(bool vertex, 
			   const int imode, const int ichan,
			   const Particle & inpart,
			   const ParticleVector & decay) const
{
  // calculate the lepton and hadron currents
  vector<LorentzPolarizationVector> hadron=hadronCurrent(vertex,imode,ichan,
							 inpart,decay);
  vector<LorentzPolarizationVector> lepton=leptonCurrent(vertex,inpart,decay);
  // get the spin density matrix for the decaying particle
  RhoDMatrix temp(2); temp.average();
  if(inpart.spinInfo())
    {
      tcSpinInfoPtr tauspin=dynamic_ptr_cast<tcSpinInfoPtr>(inpart.spinInfo());
      if(tauspin)
	{
	  tauspin->decay();
	  temp=tauspin->rhoMatrix();
	}
    }
  // work out the mapping for the hadron vector
  vector<int> constants(decay.size()+1), ispin(decay.size()),ihel(decay.size()+1);
  int itemp=1; unsigned int inu=0;
  for(int ix=int(decay.size()-1);ix>=0;--ix)
    {
      ispin[ix]=decay[ix]->data().iSpin();
      if(abs(decay[ix]->id())!=16){itemp*=ispin[ix];constants[ix]=itemp;}
      else{inu=ix;}
    }
  constants[decay.size()]=1;
  constants[inu]=constants[inu+1];
  // compute the matrix element
  DecayMatrixElement newME(2,ispin);
  for(unsigned int hhel=0;hhel<hadron.size();++hhel)
    {
      // map the index for the hadrons to a helicity state
      for(unsigned int ix=decay.size();ix>0;--ix)
	{
	  if(ix-1!=inu)
	    {
	      ihel[ix]=(hhel%constants[ix-1])/constants[ix]-int(ispin[ix-1]/2);
	      if(ispin[ix-1]%2==0&&ihel[ix]>-0&&ispin[ix-1]!=0){++ihel[ix];}
	    }
	}
      // loop over the helicities of the tau and neutrino and set up the matrix 
      // element
      unsigned int ix=0;
      for(int nhel=-1;nhel<2;nhel+=2)
	{
	  ihel[inu+1]=nhel;
	  for(int thel=-1;thel<2;thel+=2)
	    {
	      ihel[0]=thel;
	      ix = (thel+1)+(nhel+1)/2;
	      newME(ihel)= lepton[ix]*hadron[hhel];
	    }
	}
    }
  // store the matrix element
  ME(newME);
  // return the answer
  double me= 2.*(newME.contract(temp)).real()*_GF*_GF;
  return me;  
}

// the lepton currents for the different lepton helicities
vector<LorentzPolarizationVector> 
TauDecayerBase::leptonCurrent(bool vertex, const Particle & inpart,
			      const ParticleVector & decay) const
{
  // storage of the currents (first two components -1 for tau, second two +1 for tau)
  vector<LorentzPolarizationVector> temp;
  // locate the tau neutrino
  unsigned int inu=0;
  for(unsigned int ix=0;ix<decay.size();++ix){if(abs(decay[ix]->id())==16){inu=ix;}}
  // spininfo object for the outgoing neutrino
  SpinPtr nuspin,newtau;tcFermionSpinPtr hwtau,hwnewtau;
  FermionSpinPtr hwnuspin;
  if(vertex)
    {
      nuspin= new_ptr(FermionSpinInfo(decay[inu]->momentum(),true));
      hwnuspin=dynamic_ptr_cast<FermionSpinPtr>(nuspin);
    }
  if(inpart.spinInfo())
    {hwtau= dynamic_ptr_cast<tcFermionSpinPtr>(inpart.spinInfo());}
  // storage for the wavefunctions of the tau and neutrino
  vector<LorentzSpinor> wave;
  vector<LorentzSpinorBar> wavebar;
  // create a spinInfo object for the taudecay if neeed
  if(!hwtau&&vertex)
    {
      newtau   = new_ptr(FermionSpinInfo(inpart.momentum(),true));
      hwnewtau = dynamic_ptr_cast<tcFermionSpinPtr>(newtau);
      hwnewtau->decayed(true);
    }
  // incoming tau-
  if(inpart.id()==ParticleID::tauminus)
    {
      // extract or calculate the spinors for the tau
      if(hwtau)
	{
	  wave.push_back(hwtau->getDecayBasisState(-1));
	  wave.push_back(hwtau->getDecayBasisState( 1));
	}
      else
	{
	  SpinorWaveFunction temp=SpinorWaveFunction(inpart.momentum(),
						     inpart.dataPtr(),-1,incoming);
	  wave.push_back(temp.Wave());
	  temp.reset(1);wave.push_back(temp.Wave());
	  if(vertex)
	    {for(int ix=-1;ix<2;ix+=2){hwnewtau->setDecayState(ix,wave[(ix+1)/2]);}}
	}
      // calculate the wavefunctions of the neutrino
      DiracRep dirac=wave[0].Rep();
      SpinorBarWaveFunction temp=SpinorBarWaveFunction(decay[inu]->momentum(),
						       decay[inu]->dataPtr(),
						       -1,outgoing,dirac);
      wavebar.push_back(temp.Wave());
      temp.reset(1);wavebar.push_back(temp.Wave());
      if(vertex)
	{for(int ix=-1;ix<2;ix+=2){hwnuspin->setBasisState(ix,wavebar[(ix+1)/2].bar());}}
    }
  // incoming tau+
  else
    {
      // extract or calculate the spinors for the tau
      if(hwtau)
	{
	  wavebar.push_back(hwtau->getDecayBasisState(-1).bar());
	  wavebar.push_back(hwtau->getDecayBasisState( 1).bar());
	}
      else
	{
	  SpinorBarWaveFunction temp=SpinorBarWaveFunction(inpart.momentum(),
							   inpart.dataPtr(),
							   -1,incoming);
	  wavebar.push_back(temp.Wave());
	  temp.reset(1);wavebar.push_back(temp.Wave());
	  if(vertex)
	    {for(int ix=-1;ix<2;ix+=2)
	      {hwnewtau->setDecayState(ix,wavebar[(ix+1)/2].bar());}}
	}
      // calculate the wavefunctions of the neutrino
      DiracRep dirac=wavebar[0].Rep();
      SpinorWaveFunction temp = SpinorWaveFunction(decay[inu]->momentum(),
						   decay[inu]->dataPtr(),-1,
						   outgoing,dirac);
      wave.push_back(temp.Wave());
      temp.reset(1);wave.push_back(temp.Wave());
      if(vertex)
	{for(int ix=-1;ix<2;ix+=2){hwnuspin->setBasisState(ix,wave[(ix+1)/2]);}}
    }
  // set the spinInfo's if needed
  if(vertex)
    {
      decay[inu]->spinInfo(nuspin);
      if(!hwtau){const_ptr_cast<tPPtr>(&inpart)->spinInfo(newtau);}
    }
  // now compute the currents
  Complex vec[4], ii(0.,1.);
  temp.resize(4);
  for(unsigned int ix=0;ix<2;++ix)
    {
      for(unsigned int iy=0;iy<2;++iy)
	{
	  // calculate the current
	  if(wave[ix].Rep()==HaberDRep&&wavebar[iy].Rep()==HaberDRep)
	    {
	      Complex s2m4=wave[ix].s2()-wave[ix].s4();
	      Complex s1m3=wave[ix].s1()-wave[ix].s3();
	      vec[0] =   0.5*(-wavebar[iy].s1()*s2m4-wavebar[iy].s2()*s1m3
			      -wavebar[iy].s3()*s2m4-wavebar[iy].s4()*s1m3);
	      vec[1] =ii*0.5*(+wavebar[iy].s1()*s2m4-wavebar[iy].s2()*s1m3
			      +wavebar[iy].s3()*s2m4-wavebar[iy].s4()*s1m3);
	      vec[2] =   0.5*(-wavebar[iy].s1()*s1m3+wavebar[iy].s2()*s2m4
			      -wavebar[iy].s3()*s1m3+wavebar[iy].s4()*s2m4);
	      vec[3] =   0.5*(+wavebar[iy].s1()*s1m3+wavebar[iy].s2()*s2m4
			      +wavebar[iy].s3()*s1m3+wavebar[iy].s4()*s2m4);
	    }
	  else if(wave[ix].Rep()==HELASDRep&&wavebar[iy].Rep()==HELASDRep)
	    {
	      Complex s3s1=wavebar[iy].s3()*wave[ix].s1();
	      Complex s3s2=wavebar[iy].s3()*wave[ix].s2();
	      Complex s4s1=wavebar[iy].s4()*wave[ix].s1();
	      Complex s4s2=wavebar[iy].s4()*wave[ix].s2();
	      vec[0] =   -(s3s2+s4s1);
	      vec[1] = ii*(s3s2-s4s1);
	      vec[2] =   -(s3s1-s4s2);
	      vec[3] =    (s3s1+s4s2);
	    }
	  // add it to the vector
	  if(inpart.id()==15)
	    {temp[2*ix+iy]=LorentzPolarizationVector(vec[0],vec[1],vec[2],vec[3]);}
	  else
	    {temp[2*iy+ix]=LorentzPolarizationVector(vec[0],vec[1],vec[2],vec[3]);}
	}
    }
  // return the answer
  return temp;
}

}
