// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TauDecayer class.
//
//  Author: Peter Richardson
//

#include "TauDecayer.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Reference.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "TauDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Parameter.h"
#include "Herwig++/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig++/Helicity/Correlations/DecayVertex.h"
#include "ThePEG/Helicity/FermionSpinInfo.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

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
using ThePEG::Helicity::tcSpinfoPtr;
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

TauDecayer::~TauDecayer() {}

bool TauDecayer::accept(const DecayMode & dm) const {
  bool allowed=false;
  // find the neutrino 
  int idnu(0),idtemp,idin(dm.parent()->id());
  vector<int> idother;
  ParticleMSet::const_iterator pit  = dm.products().begin();
  ParticleMSet::const_iterator pend = dm.products().end();
  for( ; pit!=pend;++pit)
    {
      idtemp=(**pit).id();
      if(abs(idtemp)==16)
	{idnu=idtemp;}
      else
	{idother.push_back(idtemp);}
    }
  if((idnu==ParticleID::nu_tau    && idin==ParticleID::tauminus)||
     (idnu==ParticleID::nu_taubar && idin==ParticleID::tauplus ))
    {allowed=_current->accept(idother);}
  return allowed;
}

ParticleVector TauDecayer::decay(const DecayMode & dm,
				  const Particle & parent) const {
  // find the ids of the particles for the decay current
  ParticleMSet::const_iterator pit = dm.products().begin();
  ParticleMSet::const_iterator pend = dm.products().end();
  int idtemp;vector<int> idother;
  for( ; pit!=pend;++pit)
    {idtemp=(**pit).id();if(abs(idtemp)!=16){idother.push_back(idtemp);}}
  unsigned int itemp=_current->decayMode(idother),imode(0);
  for(unsigned int ix=0;ix<_modemap.size();++ix)
    {if(_modemap[ix]==itemp){imode=ix;}}
  // perform the decay
  bool cc=parent.id()==ParticleID::tauplus; 
  return generate(true,cc,imode,parent);
}


void TauDecayer::persistentOutput(PersistentOStream & os) const {
  os << _GF << _modemap << _current << _wgtloc << _wgtmax << _weights;
}

void TauDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _GF >> _modemap >> _current >> _wgtloc >> _wgtmax >> _weights;
}

ClassDescription<TauDecayer> TauDecayer::initTauDecayer;
// Definition of the static class description member.

void TauDecayer::Init() {

  static ClassDocumentation<TauDecayer> documentation
    ("The \\classname{TauDecayer} class is designed to use a weak current"
     " to perform the decay of the tau.");

  static Parameter<TauDecayer,InvEnergy2> interfaceGFermi
    ("GFermi",
     "The Fermi coupling constant",
     &TauDecayer::_GF, 1./GeV2, 1.16639E-5/GeV2, -1.0e12*1./GeV2, 1.0e12*1./GeV2,
     false, false, false);

  static Reference<TauDecayer,WeakDecayCurrent> interfaceWeakCurrent
    ("WeakCurrent",
     "The reference for the decay current to be used.",
     &TauDecayer::_current, false, false, true, false, false);

  static ParVector<TauDecayer,int> interfaceWeightLocation
    ("WeightLocation",
     "The locations of the weights for a given channel in the vector",
     &TauDecayer::_wgtloc,
     0, 0, 0, 0, 10000, false, false, true);

  static ParVector<TauDecayer,double> interfaceWeightMax
    ("MaximumWeight",
     "The maximum weight for a given channel.",
     &TauDecayer::_wgtmax,
     0, 0, 0, 0., 10000., false, false, true);

  static ParVector<TauDecayer,double> interfaceWeights
    ("Weights",
     "The weights for the integration.",
     &TauDecayer::_weights,
     0, 0, 0, 0., 10000., false, false, true);
}


// the lepton currents for the different lepton helicities
vector<LorentzPolarizationVector> 
TauDecayer::leptonCurrent(bool vertex, const Particle & inpart,
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
  LorentzPolarizationVector vec;
  Complex ii(0.,1.);
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
	      vec[0] =   (-wavebar[iy].s1()*s2m4-wavebar[iy].s2()*s1m3
			  -wavebar[iy].s3()*s2m4-wavebar[iy].s4()*s1m3);
	      vec[1] =ii*(+wavebar[iy].s1()*s2m4-wavebar[iy].s2()*s1m3
			      +wavebar[iy].s3()*s2m4-wavebar[iy].s4()*s1m3);
	      vec[2] =   (-wavebar[iy].s1()*s1m3+wavebar[iy].s2()*s2m4
			  -wavebar[iy].s3()*s1m3+wavebar[iy].s4()*s2m4);
	      vec[3] =   (+wavebar[iy].s1()*s1m3+wavebar[iy].s2()*s2m4
			  +wavebar[iy].s3()*s1m3+wavebar[iy].s4()*s2m4);
	    }
	  else if(wave[ix].Rep()==HELASDRep&&wavebar[iy].Rep()==HELASDRep)
	    {
	      Complex s3s1=wavebar[iy].s3()*wave[ix].s1();
	      Complex s3s2=wavebar[iy].s3()*wave[ix].s2();
	      Complex s4s1=wavebar[iy].s4()*wave[ix].s1();
	      Complex s4s2=wavebar[iy].s4()*wave[ix].s2();
	      vec[0] =   -2.*(s3s2+s4s1);
	      vec[1] = ii*2.*(s3s2-s4s1);
	      vec[2] =   -2.*(s3s1-s4s2);
	      vec[3] =    2.*(s3s1+s4s2);
	    }
	  // add it to the vector
	  if(inpart.id()==15)
	    {temp[2*ix+iy]=vec;}
	  else
	    {temp[2*iy+ix]=vec;}
	}
    }
  // return the answer
  return temp;
}

// combine the currents to give the matrix element
double TauDecayer::me2(bool vertex, const int ichan,
			   const Particle & inpart,
			   const ParticleVector & decay) const
{
  // map the mode to those in the current
  int mode=_modemap[imode()];
  // calculate the lepton and hadron currents
  vector<LorentzPolarizationVector> hadron(_current->current(vertex,mode,ichan,
							     inpart,decay));
  vector<LorentzPolarizationVector> lepton(leptonCurrent(vertex,inpart,decay));
  // get the spin density matrix for the decaying particle
  RhoDMatrix temp(2); temp.average();
  if(inpart.spinInfo())
    {
      tcSpinfoPtr tauspin=dynamic_ptr_cast<tcSpinfoPtr>(inpart.spinInfo());
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
	      if(ispin[ix-1]%2==0&&ihel[ix]>=0&&ispin[ix-1]!=0){++ihel[ix];}
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
  // multiply by the CKM element
  int iq,ia;
  _current->decayModeInfo(mode,iq,ia);
  double ckm(1.);
  if(iq<=6)
    {
      if(iq%2==0){ckm = SM().CKM(iq/2-1,(abs(ia)-1)/2);}
      else{ckm = SM().CKM(abs(ia)/2-1,(iq-1)/2);}
    }
  // return the answer
  double me= 0.5*ckm*(newME.contract(temp)).real()*_GF*_GF;
  return me;  
}

// output the setup information for the particle database
void TauDecayer::dataBaseOutput(ofstream & output)
{
  output << "update decayers set parameters=\"";
  output << "set " << fullName() << ":Iteration " << _niter << "\n";
  output << "set " << fullName() << ":Ntry " << _ntry << "\n";
  output << "set " << fullName() << ":GFermi "   << _GF*GeV2 << "\n";
  for(unsigned int ix=0;ix<_wgtloc.size();++ix)
    {output << "insert " << fullName() << ":WeightLocation " << ix << " " 
	    << _wgtloc[ix] << "\n";}
  for(unsigned int ix=0;ix<_wgtmax.size();++ix)
    {output << "insert " << fullName() << ":MaximumWeight " << ix << " " 
	    << _wgtmax[ix] << "\n";}
  for(unsigned int ix=0;ix<_weights.size();++ix)
    {output << "insert " << fullName() << ":Weights " << ix << " " 
	    << _weights[ix] << "\n";}
  _current->dataBaseOutput(output);
  output << "set " << fullName() << ":WeakCurrent " << _current->fullName() << " \n";
  output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
} 
}
