// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Tau3MesonDecayerBase class.
//
// Author: Peter Richardson

#include "Tau3MesonDecayerBase.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/DecayMode.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "Tau3MesonDecayerBase.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Helicity/EpsFunction.h"
#include "ThePEG/Helicity/ScalarSpinInfo.h"

namespace Herwig {

using namespace ThePEG;
using ThePEG::Helicity::ScalarSpinInfo;

Tau3MesonDecayerBase::~Tau3MesonDecayerBase() {}

bool Tau3MesonDecayerBase::accept(const DecayMode & dm) const {
  bool allowed=false;
  int idtau=dm.parent()->id();
  if(abs(idtau)==15 && dm.products().size() == 4)
    {
      // identify the decay products
      ParticleMSet::const_iterator pit  = dm.products().begin();
      ParticleMSet::const_iterator pend = dm.products().end();
      int idtemp;
      int npiplus(0),npiminus(0),nkplus(0),nkminus(0),npi0(0),nk0(0),nk0bar(0),neta(0);
      for( ; pit!=pend;++pit)
	{
	  idtemp=(**pit).id();
	  if(idtemp==ParticleID::piplus){++npiplus;}
	  else if(idtemp==ParticleID::piminus){++npiminus;}
	  else if(idtemp==ParticleID::Kplus){++nkplus;}
	  else if(idtemp==ParticleID::Kminus){++nkminus;}
	  else if(idtemp==ParticleID::pi0){++npi0;}
	  else if(idtemp==ParticleID::K0){++nk0;}
	  else if(idtemp==ParticleID::Kbar0){++nk0bar;}
	  else if(idtemp==ParticleID::eta){++neta;}
	}
      // work out which decay mode
      int imode=0;
      if(idtau==ParticleID::tauplus)
	{
	  if(npiplus==2&&npiminus==1){imode=1;}
	  else if(npiplus==1&&npi0==2){imode=2;}
	  else if(nkplus==1&&nkminus==1&&npiplus==1){imode=3;}
	  else if(nk0==1&&nk0bar==1&&npiplus==1){imode=4;}
	  else if(nkplus==1&&nk0bar==1&&npi0==1){imode=5;}
	  else if(nkplus==1&&npi0==2){imode=6;}
	  else if(npiplus==1&&npiminus==1&nkplus==1){imode=7;}
	  else if(nk0==1&&npiplus==1&&npi0==1){imode=8;}
	  else if(npiplus==1&&npi0==1&neta==1){imode=9;}
	}
      // decaying tau-
      else
	{
	  if(npiminus==2&&npiplus==1){imode=1;}
	  else if(npiminus==1&npi0==2){imode=2;}
	  else if(nkplus==1&nkminus==1&npiminus==1){imode=3;}
	  else if(nk0==1&&nk0bar==1&npiminus==1){imode=4;}
	  else if(nkminus==1&npi0==1&&nk0==1){imode=5;}
	  else if(npi0==2&nkminus==1){imode=6;}
	  else if(nkminus==1&npiminus==1&npiplus==1){imode=7;}
	  else if (npiminus==1&nk0bar==1&npi0==1){imode=8;}
	  else if (npiminus==1&npi0==1&neta==1){imode=9;}
	}
      if(imode==0){allowed=false;}
      else{allowed = acceptMode(imode);}
    }
  return allowed;
}
  
ParticleVector Tau3MesonDecayerBase::decay(const DecayMode & dm,
					   const Particle & parent) const {
  // produce the particles
  ParticleVector children = dm.produceProducts();
  // work out the order of the particles for the hadronic current
  int npiplus(0),npiminus(0),nkplus(0),nkminus(0),npi0(0),nk0(0),nk0bar(0),neta(0);
  for(unsigned int ix=0;ix<children.size();++ix)
    {
      if(children[ix]->id()==ParticleID::piplus){++npiplus;}
      else if(children[ix]->id()==ParticleID::piminus){++npiminus;}
      else if(children[ix]->id()==ParticleID::Kplus){++nkplus;}
      else if(children[ix]->id()==ParticleID::Kminus){++nkminus;}
      else if(children[ix]->id()==ParticleID::pi0){++npi0;}
      else if(children[ix]->id()==ParticleID::K0){++nk0;}
      else if(children[ix]->id()==ParticleID::Kbar0){++nk0bar;}
      else if(children[ix]->id()==ParticleID::eta){++neta;}
    }
  // decaying tau+
  int imodeb =0; 
  if(parent.id()==ParticleID::tauplus)
    {
      if(npiplus==2&&npiminus==1){imodeb=1;}
      else if(npiplus==1&&npi0==2){imodeb=2;}
      else if(nkplus==1&&nkminus==1&&npiplus==1){imodeb=3;}
      else if(nk0==1&&nk0bar==1&&npiplus==1){imodeb=4;}
      else if(nkplus==1&&nk0bar==1&&npi0==1){imodeb=5;}
      else if(nkplus==1&&npi0==2){imodeb=6;}
      else if(npiplus==1&&npiminus==1&nkplus==1){imodeb=7;}
      else if(nk0==1&&npiplus==1&&npi0==1){imodeb=8;}
      else if(npiplus==1&&npi0==1&neta==1){imodeb=9;}
    }
  // decaying tau-
  else
    {
      if(npiminus==2&&npiplus==1){imodeb=1;}
      else if(npiminus==1&npi0==2){imodeb=2;}
      else if(nkplus==1&nkminus==1&npiminus==1){imodeb=3;}
      else if(nk0==1&&nk0bar==1&npiminus==1){imodeb=4;}
      else if(nkminus==1&npi0==1&&nk0==1){imodeb=5;}
      else if(npi0==2&nkminus==1){imodeb=6;}
      else if(nkminus==1&npiminus==1&npiplus==1){imodeb=7;}
      else if (npiminus==1&nk0bar==1&npi0==1){imodeb=8;}
      else if (npiminus==1&npi0==1&neta==1){imodeb=9;}
    }
  // map this to the phase space integration channel
  int imode = phaseSpaceMode(imodeb);
  // generate the decay
  generate(true,imode,parent,children);
  return children;
}
  
void Tau3MesonDecayerBase::persistentOutput(PersistentOStream & os) const {
}
  
void Tau3MesonDecayerBase::persistentInput(PersistentIStream & is, int) {
}
  
ClassDescription<Tau3MesonDecayerBase> Tau3MesonDecayerBase::initTau3MesonDecayerBase;
// Definition of the static class description member.
  
void Tau3MesonDecayerBase::Init() {
    
  static ClassDocumentation<Tau3MesonDecayerBase> documentation
    ("The \\classname{Tau3MesonDecayerBase} class is designed to be the "
     "base class for "
     "the three meson decays of the tau, ie pi- pi- pi+, pi0 pi0 pi-, " 
     "K- pi- K+, K0 pi- Kbar0, K- pi0 K0,pi0 pi0 K-, K- pi- pi+, "
     "pi- Kbar0 pi0, pi- pi0 eta.");
  
}
  
// the hadronic current for this decay mode
vector<LorentzPolarizationVector> 
Tau3MesonDecayerBase::hadronCurrent(bool vertex,const int imode, const int ichan,
				    const Particle & inpart,
				    const ParticleVector & decay) const
{    
  // storage for the currents
  vector<LorentzPolarizationVector> temp;
  // work out the order of the particles for the hadronic current 
  int ipiplus[5],ipiminus[5],ikplus[5],ikminus[5],ipi0[5],ik0[5],ik0bar[5],ieta[5];
  int npiplus(0),npiminus(0),nkplus(0),nkminus(0),npi0(0),nk0(0),nk0bar(0),neta(0);
  unsigned int inu=0; int id;
  for(unsigned int ix=0;ix<decay.size();++ix)
    {
      id=decay[ix]->id();
      if(abs(id)==ParticleID::nu_tau){inu=ix;}
      else if(id==ParticleID::piplus){++npiplus;ipiplus[npiplus]=ix;}
      else if(id==ParticleID::piminus){++npiminus;ipiminus[npiminus]=ix;}
      else if(id==ParticleID::Kplus){++nkplus;ikplus[nkplus]=ix;}
      else if(id==ParticleID::Kminus){++nkminus;ikminus[nkminus]=ix;}
      else if(id==ParticleID::pi0){++npi0;ipi0[npi0]=ix;}
      else if(id==ParticleID::K0){++nk0;ik0[nk0]=ix;}
      else if(id==ParticleID::Kbar0){++nk0bar;ik0bar[nk0bar]=ix;}
      else if(id==ParticleID::eta){++neta;ieta[neta]=ix;}
    }
  // spininfo for the particles
  if(vertex)
    {
      for(unsigned int ix=0;ix<decay.size();++ix)
	{if(ix!=inu)
	    {decay[ix]->spinInfo(new_ptr(ScalarSpinInfo(decay[ix]->momentum(),true)));}}
    }
  // decaying tau+
  int imodeb =0; int iord[3]={0,0,0};
  if(inpart.id()==ParticleID::tauplus)
    {
      if(npiplus==2&&npiminus==1)
	{imodeb=1;iord[0]=ipiplus[1];iord[1]=ipiplus[2];iord[2]=ipiminus[1];}
      else if(npiplus==1&&npi0==2)
	{imodeb=2;iord[0]=ipi0[1];iord[1]=ipi0[2];iord[2]=ipiplus[1];}
      else if(nkplus==1&&nkminus==1&&npiplus==1)
	{imodeb=3;;iord[0]=ikplus[1];iord[1]=ipiplus[1];iord[2]=ikminus[1];}
      else if(nk0==1&&nk0bar==1&&npiplus==1)
	{imodeb=4;iord[0]=ik0bar[1];iord[1]=ipiplus[1];iord[2]=ik0[1];}
      else if(nkplus==1&&nk0bar==1&&npi0==1)
	{imodeb=5;iord[0]=ikplus[1];iord[1]=ipi0[1];iord[2]=ik0bar[1];}
      else if(nkplus==1&&npi0==2)
	  {imodeb=6;iord[0]=ipi0[1];iord[1]=ipi0[2];iord[2]=ikplus[1];}
      else if(npiplus==1&&npiminus==1&nkplus==1)
	{imodeb=7;iord[0]=ikplus[1];iord[1]=ipiplus[1];iord[2]=ipiminus[1];}
      else if(nk0==1&&npiplus==1&&npi0==1)
	{imodeb=8;iord[0]=ipiplus[1];iord[1]=ik0[1];iord[2]=ipi0[1];}
      else if(npiplus==1&&npi0==1&neta==1)
	  {imodeb=9;iord[0]=ipiplus[1];iord[1]=ipi0[1];iord[2]=ieta[1];}
    }
  // decaying tau-
  else
    {
      if(npiminus==2&&npiplus==1)
	{imodeb=1;iord[0]=ipiminus[1];iord[1]=ipiminus[2];iord[2]=ipiplus[1];}
      else if(npiminus==1&npi0==2)
	{imodeb=2;iord[0]=ipi0[1];iord[1]=ipi0[2];iord[2]=ipiminus[1];}
      else if(nkplus==1&nkminus==1&npiminus==1)
	{imodeb=3;iord[0]=ikminus[1];iord[1]=ipiminus[1];iord[2]=ikplus[1];}
      else if(nk0==1&&nk0bar==1&npiminus==1)
	{imodeb=4;iord[0]=ik0[1];iord[1]=ipiminus[1];iord[2]=ik0bar[1];}
      else if(nkminus==1&npi0==1&&nk0==1)
	{imodeb=5;iord[0]=ikminus[1];iord[1]=ipi0[1];iord[2]=ik0[1];}
      else if(npi0==2&nkminus==1)
	{imodeb=6;iord[0]=ipi0[1];iord[1]=ipi0[2];iord[2]=ikminus[1];}
      else if(nkminus==1&npiminus==1&npiplus==1)
	{imodeb=7;iord[0]=ikminus[1];iord[1]=ipiminus[1];iord[2]=ipiplus[1];}
      else if (npiminus==1&nk0bar==1&npi0==1)
	{imodeb=8;iord[0]=ipiminus[1];iord[1]=ik0bar[1];iord[2]=ipi0[1];}
      else if (npiminus==1&npi0==1&neta==1)
	{imodeb=9;iord[0]=ipiminus[1];iord[1]=ipi0[1];iord[2]=ieta[1];}
    }
  // calculate q2,s1,s2,s3
  Lorentz5Momentum q=0;
  for(unsigned int ix=0;ix<3;++ix){q+=decay[iord[ix]]->momentum();}
  Energy2 q2=q.m2();
  Energy2 s1 = decay[iord[1]]->momentum().m2()+decay[iord[2]]->momentum().m2()+
    2.*decay[iord[1]]->momentum().dot(decay[iord[2]]->momentum());
  Energy2 s2 = decay[iord[0]]->momentum().m2()+decay[iord[2]]->momentum().m2()+
    2.*decay[iord[0]]->momentum().dot(decay[iord[2]]->momentum());
  Energy2 s3 = decay[iord[0]]->momentum().m2()+decay[iord[1]]->momentum().m2()+
    2.*decay[iord[0]]->momentum().dot(decay[iord[1]]->momentum());
  complex<double> F1,F2,F3,F4,F5;
  calculateFormFactors(imodeb,ichan,q2,s1,s2,s3,F1,F2,F3,F4,F5);
  if(inpart.id()==ParticleID::tauplus){F5=conj(F5);}
  // the first three form-factors
  complex<double> vect[4];
  for(unsigned int ix=0;ix<4;++ix)
    {vect[ix]=
       +F1*(decay[iord[1]]->momentum()[ix]-decay[iord[2]]->momentum()[ix])
       +F2*(decay[iord[2]]->momentum()[ix]-decay[iord[0]]->momentum()[ix])
       +F3*(decay[iord[0]]->momentum()[ix]-decay[iord[1]]->momentum()[ix]);}
  // multiply by the transverse projection operator
  complex<double> dot=vect[3]*q[3]-vect[0]*q[0]-vect[1]*q[1]-vect[2]*q[2];
  dot=dot/q2;
  for(unsigned int ix=0;ix<4;++ix){vect[ix]-=dot*q[ix];}
  // the scalar form-factor
  for(unsigned int ix=0;ix<4;++ix){vect[ix]+=q[ix]*F4;}
  // the parity violating form-factor
  LorentzPolarizationVector 
    v5=Helicity::EpsFunction::product(decay[iord[0]]->momentum(),
				      decay[iord[1]]->momentum(),
				      decay[iord[2]]->momentum());
  // add this term to the result
  for(unsigned int ix=0;ix<4;++ix){vect[ix]+=v5[ix]*F5;}
  // return the answer
  temp.push_back(LorentzPolarizationVector(vect[0],vect[1],vect[2],vect[3]));
  return temp;
}
  
// modes handled by this class
bool Tau3MesonDecayerBase::acceptMode(int imode) const{return false;}

// calculate the form-factors
void Tau3MesonDecayerBase::calculateFormFactors(const int imode,const int ichan,
						Energy2 q2, Energy2 s1,
						Energy2 s2, Energy2 s3,
						complex<double> & F1,
						complex<double> & F2,
						complex<double> & F3,
						complex<double> & F4,
						complex<double> & F5) const
{F1 = 0.;F2 = 0.;F3 = 0.;F4 = 0.;F5 = 0.;}

// mapping of mode to integration channel
int Tau3MesonDecayerBase::phaseSpaceMode(int in) const {return 0;}

}
