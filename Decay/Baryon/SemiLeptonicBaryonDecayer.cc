// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SemiLeptonicBaryonDecayer class.
//

#include "SemiLeptonicBaryonDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/PDT/DecayMode.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "SemiLeptonicBaryonDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/RSSpinorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/RSSpinorBarWaveFunction.h"
#include "ThePEG/Helicity/LorentzPolarizationVector.h"
#include "ThePEG/Helicity/FermionSpinInfo.h"
#include "ThePEG/Helicity/RSFermionSpinInfo.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;
using ThePEG::Helicity::RhoDMatrix;
using ThePEG::Helicity::LorentzPolarizationVector;
using ThePEG::Helicity::tcFermionSpinPtr;
using ThePEG::Helicity::FermionSpinInfo;
using ThePEG::Helicity::tcRSFermionSpinPtr;
using ThePEG::Helicity::RSFermionSpinInfo;
using Helicity::SpinorWaveFunction;
using Helicity::SpinorBarWaveFunction;
using Helicity::RSSpinorWaveFunction;
using Helicity::RSSpinorBarWaveFunction;
using Helicity::Direction;
using Helicity::incoming;
using Helicity::outgoing;

SemiLeptonicBaryonDecayer::~SemiLeptonicBaryonDecayer() {}

void SemiLeptonicBaryonDecayer::doinit() throw(InitException) {
  DecayIntegrator::doinit();
  // make sure the current got initialised
  _current->init();
  // and the form factors
  _form->init();
  // the channels
  PDVector extpart,ptemp;
  _modemap.resize(0);
  double maxweight;
  vector<double> channelwgts(1,1.);
  int id0(0),id1(0),Wcharge(0),inspin,spect1,spect2,inquark,outquark,outspin;
  Energy min;
  DecayPhaseSpaceModePtr mode;
  DecayPhaseSpaceChannelPtr channel;
  int iq,ia;
  for(unsigned int ix=0;ix<_form->numberOfFactors();++ix)
    {
      // get the external particles for this mode
      extpart.resize(2);
      _form->particleID(ix,id0,id1);
      _form->formFactorInfo(ix,inspin,outspin,spect1,spect2,inquark,outquark);
      extpart[0]=getParticleData(id0);
      extpart[1]=getParticleData(id1);
      Wcharge =(extpart[0]->iCharge()-extpart[1]->iCharge());
      min = extpart[0]->mass()+extpart[0]->widthUpCut()
	-extpart[1]->mass()+extpart[1]->widthLoCut();
      _modemap.push_back(numberModes());
      for(unsigned int iy=0;iy<_current->numberOfModes();++iy)
	{
	  extpart.resize(2);
	  ptemp=_current->particles(Wcharge,iy,iq,ia);
	  for(unsigned int iz=0;iz<ptemp.size();++iz){extpart.push_back(ptemp[iz]);}
	  // create the mode
	  mode=new_ptr(DecayPhaseSpaceMode(extpart,this));
	  // create the first piece of the channel
	  channel = new_ptr(DecayPhaseSpaceChannel(mode));
	  channel->addIntermediate(extpart[0],0,0.0,-1,1);
	  bool done=_current->createMode(Wcharge,iy,mode,2,1,channel,min);
	  if(done&&abs(Wcharge)==3&&inspin==2&&(outspin==2||outspin==4))
	    {
	      // the maximum weight
	      if(_maxwgt.size()>numberModes()){maxweight=_maxwgt[numberModes()];}
	      else{maxweight=2.;}
	      channelwgts.resize(mode->numberChannels(),1./(mode->numberChannels()));
	      addMode(mode,maxweight,channelwgts);
	    }
	  
	}
    }
}

bool SemiLeptonicBaryonDecayer::accept(const DecayMode & dm) const {
  // find the non-lepton
  int ibar(0),idtemp,idin(dm.parent()->id());
  vector<int> idother; bool dummy;
  ParticleMSet::const_iterator pit  = dm.products().begin();
  ParticleMSet::const_iterator pend = dm.products().end();
  for( ; pit!=pend;++pit)
    {
      idtemp=(**pit).id();
      if(abs(idtemp)>16){ibar=idtemp;}
      else{idother.push_back(idtemp);}
    }
  // check that the form factor exists
  if(_form->formFactorNumber(idin,ibar,dummy)<0){return false;}
  // and the current
  return _current->accept(idother);
}

ParticleVector SemiLeptonicBaryonDecayer::decay(const DecayMode & dm,
				  const Particle & parent) const {
  // find the ids of the particles for the decay current
  ParticleMSet::const_iterator pit = dm.products().begin();
  ParticleMSet::const_iterator pend = dm.products().end();
  int idtemp,ibar(0),idin(dm.parent()->id());
  vector<int> idother;
  bool cc(false);
  for( ; pit!=pend;++pit)
    {
      idtemp=(**pit).id();
      if(abs(idtemp)>16){ibar=idtemp;}
      else{idother.push_back(idtemp);}
    }
  int imode = _modemap[_form->formFactorNumber(idin,ibar,cc)]
    +_current->decayMode(idother);
  // perform the decay
  return generate(true,cc,imode,parent);
}


void SemiLeptonicBaryonDecayer::persistentOutput(PersistentOStream & os) const {
  os << _current << _form << _maxwgt << _modemap << _GF;}

void SemiLeptonicBaryonDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _current >> _form >> _maxwgt >> _modemap >> _GF;}

ClassDescription<SemiLeptonicBaryonDecayer> SemiLeptonicBaryonDecayer::initSemiLeptonicBaryonDecayer;
// Definition of the static class description member.

void SemiLeptonicBaryonDecayer::Init() {

  static ClassDocumentation<SemiLeptonicBaryonDecayer> documentation
    ("The \\classname{SemiLeptonicBaryonDecayer} class is designed for"
     " the semi-leptonic decay of the baryons.");

  static Parameter<SemiLeptonicBaryonDecayer,InvEnergy2> interfaceGFermi
    ("GFermi",
     "The Fermi coupling constant",
     &SemiLeptonicBaryonDecayer::_GF, 1./GeV2, 1.16639E-5/GeV2,
     -1.0e12*1./GeV2, 1.0e12*1./GeV2,
     false, false, false);

  static Reference<SemiLeptonicBaryonDecayer,LeptonNeutrinoCurrent> interfaceCurrent
    ("Current",
     "The current for the leptons produced in the decay.",
     &SemiLeptonicBaryonDecayer::_current, true, true, true, false, false);


  static Reference<SemiLeptonicBaryonDecayer,BaryonFormFactor> interfaceFormFactor
    ("FormFactor",
     "The form factor",
     &SemiLeptonicBaryonDecayer::_form, true, true, true, false, false);

  static ParVector<SemiLeptonicBaryonDecayer,double> interfaceMaximumWeight
    ("MaximumWeight",
     "The maximum weights for the decays",
     &SemiLeptonicBaryonDecayer::_maxwgt,
     0, 0, 0, 0, 10000, false, false, true);

}

// combine the currents and form-factors to give the matrix element
double SemiLeptonicBaryonDecayer::me2(bool vertex, const int ichan,
				      const Particle & inpart,
				      const ParticleVector & decay) const
{
  double me(0.);
  if(inpart.dataPtr()->iSpin()==2)
    {
      if(decay[0]->dataPtr()->iSpin()==2)
	{me=halfHalf(vertex,ichan,inpart,decay);}
      else if(decay[0]->dataPtr()->iSpin()==4)
	{me=halfThreeHalf(vertex,ichan,inpart,decay);}
      else
	{throw DecayIntegratorError() << "Invalid spins "
				  << " in SemiLeptonicBaryonDecayer::me2()" 
				  << Exception::abortnow;}
    }
  else
    {throw DecayIntegratorError() << "Invalid spins "
				  << " in SemiLeptonicBaryonDecayer::me2()" 
				  << Exception::abortnow;}
  return me;
}

// matrix element for a 1/2 -> 1/2 semi-leptonic decay
double SemiLeptonicBaryonDecayer::halfHalf(bool vertex, const int ichan,
					   const Particle & inpart,
					   const ParticleVector & decay) const
{
  RhoDMatrix temp(inpart.dataPtr()->iSpin()); temp.average();
  // check if the decay particle has spin info 
   tcFermionSpinPtr inspin;
  if(inpart.spinInfo())
    {inspin = dynamic_ptr_cast<tcFermionSpinPtr>(inpart.spinInfo());}
  // if the spin info object exists use it
  vector<LorentzSpinor> sp;
  vector<LorentzSpinorBar> sbar;
  if(inspin)
    {
      inspin->decay();
      temp=inspin->rhoMatrix();
      if(inpart.id()>0)
	{
	  sp.push_back(inspin->getDecayBasisState(-1));
	  sp.push_back(inspin->getDecayBasisState( 1));
	}
      else
	{
	  sbar.push_back(inspin->getDecayBasisState(-1).bar());
	  sbar.push_back(inspin->getDecayBasisState( 1).bar());
	}
    }
  else if(inpart.spinInfo())
    {throw DecayIntegratorError() << "Wrong type of spin info for the incoming particle"
				  << " in SemiLeptonicBaryonDecayer::me2()" 
				  << Exception::abortnow;}
  // spininfo does not exist create it
  else
    {
      SpinPtr newspin=new_ptr(FermionSpinInfo(inpart.momentum(),true));
      inspin = dynamic_ptr_cast<tcFermionSpinPtr>(newspin);
      inspin->decayed(true);
      const_ptr_cast<tPPtr>(&inpart)->spinInfo(newspin);
      if(inpart.id()>0)
	{
	  SpinorWaveFunction stemp(inpart.momentum(),inpart.dataPtr(),-1,incoming);
	  sp.push_back(stemp.Wave());
	  stemp.reset(1);sp.push_back(stemp.Wave());
	  if(vertex)
	    {for(int ix=-1;ix<2;ix+=2)
		{inspin->setDecayState(ix,sp[(ix+1)/2]);}}
	  
	}
      else
	{
	  SpinorBarWaveFunction stemp(inpart.momentum(),inpart.dataPtr(),-1,incoming);
	  sbar.push_back(stemp.Wave());
	  stemp.reset(1);sbar.push_back(stemp.Wave());
	  if(vertex)
	    {for(int ix=-1;ix<2;ix+=2)
	      {inspin->setDecayState(ix,sbar[(ix+1)/2].bar());}}
	}
    }
  // get the information on the form-factor
  int spinin(0),spinout(0),spect1,spect2,inquark,outquark;
  int id0=inpart.id(),id1=decay[0]->id();
  bool cc; int iloc(_form->formFactorNumber(id0,id1,cc));
  _form->formFactorInfo(iloc,spinin,spinout,spect1,spect2,inquark,outquark);
  tcFermionSpinPtr outspin;
  // construct the spin info for the outgoing baryon
  if(vertex)
    {
      SpinPtr stemp=new_ptr(FermionSpinInfo(decay[0]->momentum(),true));
      outspin=dynamic_ptr_cast<tcFermionSpinPtr>(stemp);
      decay[0]->spinInfo(stemp);
    }
  // work out the value of q and calculate the form factors
  Lorentz5Momentum q=inpart.momentum()-decay[0]->momentum();q.rescaleMass();
  Energy m0(inpart.mass()),m1(decay[0]->mass());
  Energy2 q2=q.mass2();
  Lorentz5Momentum sum=inpart.momentum()+decay[0]->momentum();
  // calculate the hadronic current for the decay
  vector<LorentzPolarizationVector> hadron;
  // calculate the wavefunction for the outgoing particles
  if(decay[0]->id()>0)
    {
      DiracRep dirac=sp[0].Rep();
      SpinorBarWaveFunction stemp=SpinorBarWaveFunction(decay[0]->momentum(),
							decay[0]->dataPtr(),
							-1,outgoing,dirac);
      sbar.push_back(stemp.Wave());
      stemp.reset(1);sbar.push_back(stemp.Wave());
      if(vertex)
	{for(int ix=-1;ix<2;ix+=2)
	    {outspin->setBasisState(ix,sbar[(ix+1)/2].bar());}}
    }
  else
    {
      DiracRep dirac=sbar[0].Rep();
      SpinorWaveFunction stemp = SpinorWaveFunction(decay[0]->momentum(),
						    decay[0]->dataPtr(),-1,
						    outgoing,dirac);
      sp.push_back(stemp.Wave());
      stemp.reset(1);sp.push_back(stemp.Wave());
      if(vertex)
	{for(int ix=-1;ix<2;ix+=2)
	    {outspin->setBasisState(ix,sp[(ix+1)/2]);}}
    }
  // calculate the form factors
  Complex f1v,f2v,f3v,f1a,f2a,f3a;
  _form->SpinHalfSpinHalfFormFactor(q2,iloc,id0,id1,m0,m1,
				    f1v,f2v,f3v,f1a,f2a,f3a);
  // now we need to construct the current
  LorentzPolarizationVector vtemp;
  hadron.resize(4);
  Complex ii(0.,1.),left,right;
  left =  f1v-f1a-f2v-(m0-m1)/(m0+m1)*f2a;
  right = f1v+f1a-f2v+(m0-m1)/(m0+m1)*f2a;
  for(unsigned int ix=0;ix<2;++ix)
    {
      for(unsigned int iy=0;iy<2;++iy)
	{
	  Complex vspin,aspin;
	  if(sp[ix].Rep()==HaberDRep)
	    {
	      // the vector and axial vector pieces
	      Complex s2m4(sp[ix].s2()-sp[ix].s4()),s1m3(sp[ix].s1()-sp[ix].s3());
	      Complex s1p3(sp[ix].s1()+sp[ix].s3()),s2p4(sp[ix].s2()+sp[ix].s4());
	      // the left hand like piece
	      vtemp[0] =   0.5*left*(-sbar[iy].s1()*s2m4-sbar[iy].s2()*s1m3
				     -sbar[iy].s3()*s2m4-sbar[iy].s4()*s1m3);
	      vtemp[1] =ii*0.5*left*(+sbar[iy].s1()*s2m4-sbar[iy].s2()*s1m3
				     +sbar[iy].s3()*s2m4-sbar[iy].s4()*s1m3);
	      vtemp[2] =   0.5*left*(-sbar[iy].s1()*s1m3+sbar[iy].s2()*s2m4
				     -sbar[iy].s3()*s1m3+sbar[iy].s4()*s2m4);
	      vtemp[3] =   0.5*left*(+sbar[iy].s1()*s1m3+sbar[iy].s2()*s2m4
				     +sbar[iy].s3()*s1m3+sbar[iy].s4()*s2m4);
	      // the right hand like piece
	      vtemp[0] +=    +0.5*right*(+sbar[iy].s1()*s2p4+sbar[iy].s2()*s1p3
					 -sbar[iy].s3()*s2p4-sbar[iy].s4()*s1p3);
	      vtemp[1] += +ii*0.5*right*(-sbar[iy].s1()*s2p4+sbar[iy].s2()*s1p3
					 +sbar[iy].s3()*s2p4-sbar[iy].s4()*s1p3);
	      vtemp[2] +=    +0.5*right*(+sbar[iy].s1()*s1p3-sbar[iy].s2()*s2p4
					 -sbar[iy].s3()*s1p3+sbar[iy].s4()*s2p4);
	      vtemp[3] +=    +0.5*right*(+sbar[iy].s1()*s1p3+sbar[iy].s2()*s2p4
					 -sbar[iy].s3()*s1p3-sbar[iy].s4()*s2p4);
	      // product for the momentum terms
	      vspin = 
		sbar[iy].s1()*sp[ix].s1()+sbar[iy].s2()*sp[ix].s2()
		+sbar[iy].s3()*sp[ix].s3()+sbar[iy].s4()*sp[ix].s4();
	      aspin = 
		sbar[iy].s1()*sp[ix].s3()+sbar[iy].s2()*sp[ix].s4()
		+sbar[iy].s3()*sp[ix].s1()+sbar[iy].s4()*sp[ix].s2();
	    }
	  else
	    {
	      Complex 
		s3s1(sbar[iy].s3()*sp[ix].s1()),s3s2(sbar[iy].s3()*sp[ix].s2()),
		s4s1(sbar[iy].s4()*sp[ix].s1()),s4s2(sbar[iy].s4()*sp[ix].s2()),
		s1s4(sbar[iy].s1()*sp[ix].s4()),s1s3(sbar[iy].s1()*sp[ix].s3()),
		s2s3(sbar[iy].s2()*sp[ix].s3()),s2s4(sbar[iy].s2()*sp[ix].s4());
	      vtemp[0] =    -left*(s3s2+s4s1)+right*(s1s4+s2s3);
	      vtemp[1] = ii*(left*(s3s2-s4s1)-right*(s1s4-s2s3));
	      vtemp[2] =    -left*(s3s1-s4s2)+right*(s1s3-s2s4);
	      vtemp[3] =     left*(s3s1+s4s2)+right*(s1s3+s2s4);
	      // product for the momentum terms
	      vspin = 
		sbar[iy].s1()*sp[ix].s1()+sbar[iy].s2()*sp[ix].s2()
		+sbar[iy].s3()*sp[ix].s3()+sbar[iy].s4()*sp[ix].s4();
	      aspin = 
		-sbar[iy].s1()*sp[ix].s1()-sbar[iy].s2()*sp[ix].s2()
		+sbar[iy].s3()*sp[ix].s3()+sbar[iy].s4()*sp[ix].s4();
	    }
	  // the momentum like pieces
	  if(inpart.id()>0)
	    {
	      vtemp+= (f2v*vspin+f2a*aspin)/(m0+m1)*sum;
	      vtemp+= (f3v*vspin+f3a*aspin)/(m0+m1)*q;
	    }
	  else
	    {
	      vtemp+= (f2v*vspin-f2a*aspin)/(m0+m1)*sum;
	      vtemp+= (f3v*vspin-f3a*aspin)/(m0+m1)*q;
	    }
	  if(inpart.id()>0){hadron[2*ix+iy]=vtemp;}
	  else{hadron[2*iy+ix]=vtemp;}
	}
    }
  // construct the lepton current
  int mode=(abs(decay[1]->id())-11)/12;
  Energy scale;
  ParticleVector leptons;
  leptons.push_back(decay[decay.size()-2]);
  leptons.push_back(decay[decay.size()-1]);
  vector<LorentzPolarizationVector> lepton(_current->current(vertex,mode,ichan,
							     scale,leptons));
  // work out the mapping for the lepton vector
  vector<int> constants(decay.size()+1), ispin(decay.size()),ihel(decay.size()+1);
  int itemp=1; unsigned int ibar=0;
  for(int ix=int(decay.size()-1);ix>=0;--ix)
    {
      ispin[ix]=decay[ix]->data().iSpin();
      if(abs(decay[ix]->id())<=16){itemp*=ispin[ix];constants[ix]=itemp;}
      else{ibar=ix;}
    }
  constants[decay.size()]=1;
  constants[ibar]=constants[ibar+1];
  DecayMatrixElement newME(spinin,ispin);
  for(unsigned int mhel=0;mhel<hadron.size();++mhel)
    {
      ihel[0     ]=2*(mhel/spinout)-spinin/2;
      ihel[ibar+1]=2*(mhel%spinout)-spinout/2;
      for(unsigned int lhel=0;lhel<lepton.size();++lhel)
	{
	  // map the index for the leptons to a helicity state
	  for(unsigned int ix=decay.size();ix>0;--ix)
	    {
	      if(ix-1!=ibar)
		{
		  ihel[ix]=(lhel%constants[ix-1])/constants[ix]-int(ispin[ix-1]/2);
		  if(ispin[ix-1]%2==0&&ihel[ix]>=0&&ispin[ix-1]!=0){++ihel[ix];}
		}
	    }
	  newME(ihel)= lepton[lhel]*hadron[mhel];
	}
    }  
  // store the matrix element
  ME(newME);
  // ckm factor
  double ckm(1.);
  if(inquark<=6)
    {
      if(inquark%2==0){ckm = SM().CKM(inquark/2-1,(abs(outquark)-1)/2);}
      else{ckm = SM().CKM(abs(outquark)/2-1,(inquark-1)/2);}
    }
  // return the answer
  return 0.5*(newME.contract(temp)).real()*_GF*_GF*ckm;
}

// matrix element for a 1/2 -> 3/2 semi-leptonic decay
double SemiLeptonicBaryonDecayer::halfThreeHalf(bool vertex, const int ichan,
						const Particle & inpart,
						const ParticleVector & decay) const
{
  RhoDMatrix temp(inpart.dataPtr()->iSpin()); temp.average();
  // check if the decay particle has spin info
  tcFermionSpinPtr inspin;
  if(inpart.spinInfo())
    {inspin = dynamic_ptr_cast<tcFermionSpinPtr>(inpart.spinInfo());}
  // if the spin info object exists use it
  vector<LorentzSpinor> sp;
  vector<LorentzSpinorBar> sbar;
  vector<LorentzRSSpinor> RSsp;
  vector<LorentzRSSpinorBar> RSsbar;
  if(inspin)
    {
      inspin->decay();
      temp=inspin->rhoMatrix();
      if(inpart.id()>0)
	{
	  sp.push_back(inspin->getDecayBasisState(-1));
	  sp.push_back(inspin->getDecayBasisState( 1));
	}
      else
	{
	  sbar.push_back(inspin->getDecayBasisState(-1).bar());
	  sbar.push_back(inspin->getDecayBasisState( 1).bar());
	}
    }
  else if(inpart.spinInfo())
    {throw DecayIntegratorError() << "Wrong type of spin info for the incoming particle"
				  << " in SemiLeptonicBaryonDecayer::halfThreeHalf()" 
				  << Exception::abortnow;}
  // spininfo does not exist create it
  else
    {
      SpinPtr newspin=new_ptr(FermionSpinInfo(inpart.momentum(),true));
      inspin = dynamic_ptr_cast<tcFermionSpinPtr>(newspin);
      inspin->decayed(true);
      const_ptr_cast<tPPtr>(&inpart)->spinInfo(newspin);
      if(inpart.id()>0)
	{
	  SpinorWaveFunction stemp(inpart.momentum(),inpart.dataPtr(),-1,incoming);
	  sp.push_back(stemp.Wave());
	  stemp.reset(1);sp.push_back(stemp.Wave());
	  if(vertex)
	    {for(int ix=-1;ix<2;ix+=2)
		{inspin->setDecayState(ix,sp[(ix+1)/2]);}}
	  
	}
      else
	{
	  SpinorBarWaveFunction stemp(inpart.momentum(),inpart.dataPtr(),-1,incoming);
	  sbar.push_back(stemp.Wave());
	  stemp.reset(1);sbar.push_back(stemp.Wave());
	  if(vertex)
	    {for(int ix=-1;ix<2;ix+=2)
	      {inspin->setDecayState(ix,sbar[(ix+1)/2].bar());}}
	}
    }
  // get the information on the form-factor
  int spinin(0),spinout(0),inquark,outquark,spect1,spect2;
  int id0=inpart.id(),id1=decay[0]->id();
  bool cc; int iloc(_form->formFactorNumber(id0,id1,cc));
  _form->formFactorInfo(iloc,spinin,spinout,spect1,spect2,inquark,outquark);
  tcRSFermionSpinPtr outspin;
  // construct the spin info for the outgoing baryon
  if(vertex)
    {
      SpinPtr stemp=new_ptr(RSFermionSpinInfo(decay[0]->momentum(),true));
      outspin=dynamic_ptr_cast<tcRSFermionSpinPtr>(stemp);
      decay[0]->spinInfo(stemp);
    }
  // work out the value of q and calculate the form factors
  Lorentz5Momentum q=inpart.momentum()-decay[0]->momentum();q.rescaleMass();
  Energy m0(inpart.mass()),m1(decay[0]->mass());
  Energy2 q2=q.mass2();
  Lorentz5Momentum sum=inpart.momentum()+decay[0]->momentum();
  // calculate the hadronic current for the decay
  LorentzPolarizationVector hadron[4][2];
  // calculate the wavefunction for the outgoing particles
  if(decay[0]->id()>0)
    {
      DiracRep dirac=sp[0].Rep();
      RSSpinorBarWaveFunction temp(decay[0]->momentum(),decay[0]->dataPtr(),-2,
				   outgoing,dirac);
      for(int ix=-2;ix<3;++ix)
	{
	  if(ix!=0)
	    {
	      temp.reset(ix);
	      sbar.push_back(temp.Wave().dot(inpart.momentum()));
	      RSsbar.push_back(temp.Wave());
	      if(vertex){outspin->setBasisState(ix,temp.Wave().bar());}
	    }
	}
    }
  else
    {
      DiracRep dirac=sbar[0].Rep();
      RSSpinorWaveFunction temp(decay[0]->momentum(),decay[0]->dataPtr(),-2,
				outgoing,dirac);
      for(int ix=-2;ix<3;++ix)
	{
	  if(ix!=0)
	    {
	      temp.reset(ix);
	      sp.push_back(temp.Wave().dot(inpart.momentum()));
	      RSsp.push_back(temp.Wave());
	      if(vertex){outspin->setBasisState(ix,temp.Wave());}
	    }
	}
    }
  // calculate the form factors
  Complex f1v,f2v,f3v,f4v,f1a,f2a,f3a,f4a;
  _form->SpinHalfSpinThreeHalfFormFactor(q2,iloc,id0,id1,m0,m1,
					 f1v,f2v,f3v,f4v,f1a,f2a,f3a,f4a);
  // now we need to construct the current
  LorentzPolarizationVector vtemp;
  Complex ii(0.,1.),lS1,lS2,rS1,rS2,left,right,lV,rV;
  InvEnergy ms(1./(m0+m1)),ms2(ms*ms);
  if(inpart.id()>0)
    {
      left  = f1a-f1v;
      right = f1a+f1v; 
      lS1   = ms2*(f3a-f4a-f3v+f4v);
      rS1   = ms2*(f3a-f4a+f3v-f4v);
      lS2   = ms2*(f4a-f4v);
      rS2   = ms2*(f4a+f4v);
      lV    = ms*(f2a-f2v);
      rV    = ms*(f2a+f2v);
    }
  else
    {
      left  = conj(f1a+f1v);
      right = conj(f1a-f1v); 
      lS1   = ms2*conj(f3a-f4a+f3v-f4v);
      rS1   = ms2*conj(f3a-f4a-f3v+f4v);
      lS2   = ms2*conj(f4a-f4v);
      rS2   = ms2*conj(f4a+f4v);
      lV    = ms *conj(f2a-f2v);
      rV    = ms *conj(f2a+f2v);
    }
  unsigned int ix,iy,ixa,iya,iz;
  // construct the vectors for the decay
  Complex scalar1,scalar2,lfact,rfact,s2m4,s1m3,s1p3,s2p4,s3s2,s1s4,s3s1,s1s3,
    s4s1,s2s3,s4s2,s2s4;
  LorentzPolarizationVector svec,tvec;
  for(iya=0;iya<4;++iya)
    {
      for(ixa=0;ixa<2;++ixa)
	{
	  if(decay[0]->id()>0){ix=iya;iy=ixa;}
	  else{ix=ixa;iy=iya;}
	  // low energy
	  if(sp[iy].Rep()==HaberDRep)
	    {
	      // scalar like terms
	      // left and right pieces
	      lfact = (sbar[ix].s1()-sbar[ix].s3())*(sp[iy].s1()-sp[iy].s3())
		+(sbar[ix].s2()-sbar[ix].s4())*(sp[iy].s2()-sp[iy].s4());
	      rfact = (sbar[ix].s1()+sbar[ix].s3())*(sp[iy].s1()+sp[iy].s3())
		+(sbar[ix].s2()+sbar[ix].s4())*(sp[iy].s2()+sp[iy].s4());
	      scalar1 = 0.5*(lS1*lfact+rS1*rfact);
	      scalar2 = 0.5*(lS2*lfact+rS2*rfact);
	      // vector like term
	      s2m4=sp[iy].s2()-sp[iy].s4();
	      s1m3=sp[iy].s1()-sp[iy].s3();
	      s1p3=sp[iy].s1()+sp[iy].s3();
	      s2p4=sp[iy].s2()+sp[iy].s4();
	      // vector like term
	      svec[0] =   0.5*lV*(-sbar[ix].s1()*s2m4-sbar[ix].s2()*s1m3
				  -sbar[ix].s3()*s2m4-sbar[ix].s4()*s1m3)
		+0.5*rV*(+sbar[ix].s1()*s2p4+sbar[ix].s2()*s1p3
			 -sbar[ix].s3()*s2p4-sbar[ix].s4()*s1p3);
	      svec[1] =ii*0.5*lV*(+sbar[ix].s1()*s2m4-sbar[ix].s2()*s1m3
				  +sbar[ix].s3()*s2m4-sbar[ix].s4()*s1m3)
		+ii*0.5*rV*(-sbar[ix].s1()*s2p4+sbar[ix].s2()*s1p3
			    +sbar[ix].s3()*s2p4-sbar[ix].s4()*s1p3);
	      svec[2] =   0.5*lV*(-sbar[ix].s1()*s1m3+sbar[ix].s2()*s2m4
				  -sbar[ix].s3()*s1m3+sbar[ix].s4()*s2m4)
		+0.5*rV*(+sbar[ix].s1()*s1p3-sbar[ix].s2()*s2p4
			 -sbar[ix].s3()*s1p3+sbar[ix].s4()*s2p4);
	      svec[3] =   0.5*lV*(+sbar[ix].s1()*s1m3+sbar[ix].s2()*s2m4
				  +sbar[ix].s3()*s1m3+sbar[ix].s4()*s2m4)
		+0.5*rV*(+sbar[ix].s1()*s1p3+sbar[ix].s2()*s2p4
			 -sbar[ix].s3()*s1p3-sbar[ix].s4()*s2p4);
	      if(inpart.id()>0)
		{
		  for(iz=0;iz<4;++iz)
		    {
		      tvec[iz] = 
			0.5*left*((sbar[ix].s1()-sbar[ix].s3())*
			      (RSsp[iy](iz,0)-RSsp[iy](iz,2))
			      +(sbar[ix].s2()-sbar[ix].s4())*
			      (RSsp[iy](iz,1)-RSsp[iy](iz,3)))
			+0.5*right*( (sbar[ix].s1()+sbar[ix].s3())*
				 (RSsp[iy](iz,0)+RSsp[iy](iz,2))
				 +(sbar[ix].s2()+sbar[ix].s4())*
				 (RSsp[iy](iz,1)+RSsp[iy](iz,3)));
		    }
		}
	      else
		{
		  for(iz=0;iz<4;++iz)
		    {
		      tvec[iz] = 
			0.5*left*((RSsbar[ix](iz,0)-RSsbar[ix](iz,2))*
			      (sp[iy].s1()-sp[iy].s3())
			      +(RSsbar[ix](iz,1)-RSsbar[ix](iz,3))*
			      (sp[iy].s2()-sp[iy].s4()))
			+0.5*right*( (RSsbar[ix](iz,0)+RSsbar[ix](iz,2))*
				 (sp[iy].s1()+sp[iy].s3())
				 +(RSsbar[ix](iz,1)+RSsbar[ix](iz,3))*
				 (sp[iy].s2()+sp[iy].s4()));
		    }
		}
	    }
	  // high energy conventions
	  else if(sp[iy].Rep()==HELASDRep&&sbar[ix].Rep()==HELASDRep)
	    { 
	      // scalar like term
	      lfact = sbar[ix].s1()*sp[iy].s1()+sbar[ix].s2()*sp[iy].s2();
	      rfact = sbar[ix].s3()*sp[iy].s3()+sbar[ix].s4()*sp[iy].s4();
	      scalar1 = lS1*lfact+rS1*rfact;
	      scalar2 = lS2*lfact+rS2*rfact;
	      // vector like term
	      s3s2=sbar[ix].s3()*sp[iy].s2();s4s1=sbar[ix].s4()*sp[iy].s1();
	      s1s4=sbar[ix].s1()*sp[iy].s4();s2s3=sbar[ix].s2()*sp[iy].s3();
	      s3s1=sbar[ix].s3()*sp[iy].s1();s4s2=sbar[ix].s4()*sp[iy].s2();
	      s1s3=sbar[ix].s1()*sp[iy].s3();s2s4=sbar[ix].s2()*sp[iy].s4();
	      svec[0] =    -lV*(s3s2+s4s1)+rV*(s1s4+s2s3);
	      svec[1] = ii*(lV*(s3s2-s4s1)-rV*(s1s4-s2s3));
	      svec[2] =    -lV*(s3s1-s4s2)+rV*(s1s3-s2s4);
	      svec[3] =     lV*(s3s1+s4s2)+rV*(s1s3+s2s4);
	      if(inpart.id()>0)
		{
		  for(iz=0;iz<4;++iz)
		    {
		      tvec[iz]=  left*(RSsbar[ix](iz,0)*sp[iy].s1()
				       +RSsbar[ix](iz,1)*sp[iy].s2())
			+right*(RSsbar[ix](iz,2)*sp[iy].s3()
				+RSsbar[ix](iz,3)*sp[iy].s4());
		    }
		}
	      else
		{
		  for(iz=0;iz<4;++iz)
		    {
		      tvec[iz]=  
			left*(sbar[ix].s1()*RSsp[iy](iz,0)+sbar[ix].s2()*RSsp[iy](iz,1))
			+right*(sbar[ix].s3()*RSsp[iy](iz,2)+sbar[ix].s4()*RSsp[iy](iz,3));
		    }
		}
	    }
	  // mixed conventions
	  else
	    {
	      // scalar like term
	      lfact = sbar[ix].s1()*sp[iy].s1()+sbar[ix].s2()*sp[iy].s2();
	      rfact = sbar[ix].s3()*sp[iy].s3()+sbar[ix].s4()*sp[iy].s4();
	      scalar1 = lS1*lfact+rS1*rfact;
	      scalar2 = lS2*lfact+rS2*rfact;
	      // vector like term
	      s3s2=sbar[ix].s3()*sp[iy].s2();s4s1=sbar[ix].s4()*sp[iy].s1();
	      s1s4=sbar[ix].s1()*sp[iy].s4();s2s3=sbar[ix].s2()*sp[iy].s3();
	      s3s1=sbar[ix].s3()*sp[iy].s1();s4s2=sbar[ix].s4()*sp[iy].s2();
	      s1s3=sbar[ix].s1()*sp[iy].s3();s2s4=sbar[ix].s2()*sp[iy].s4();
	      svec[0] =    -lV*(s3s2+s4s1)+rV*(s1s4+s2s3);
	      svec[1] = ii*(lV*(s3s2-s4s1)-rV*(s1s4-s2s3));
	      svec[2] =    -lV*(s3s1-s4s2)+rV*(s1s3-s2s4);
	      svec[3] =     lV*(s3s1+s4s2)+rV*(s1s3+s2s4);
	      if(inpart.id()>0)
		{
		  for(iz=0;iz<4;++iz)
		    {
		      tvec[iz] = 
			0.5*left*((RSsbar[ix](iz,0)-RSsbar[ix](iz,2))*
			      (sp[iy].s1()-sp[iy].s3())
			      +(RSsbar[ix](iz,1)-RSsbar[ix](iz,3))*
			      (sp[iy].s2()-sp[iy].s4()))
			+0.5*right*( (RSsbar[ix](iz,0)+RSsbar[ix](iz,2))*
				 (sp[iy].s1()+sp[iy].s3())
				 +(RSsbar[ix](iz,1)+RSsbar[ix](iz,3))*
				 (sp[iy].s2()+sp[iy].s4()));
		    }
		}
	      else
		{
		  for(iz=0;iz<4;++iz)
		    {
		      tvec[iz] = 
			0.5*left*((sbar[ix].s1()-sbar[ix].s3())*
			      (RSsp[iy](iz,0)-RSsp[iy](iz,2))
			      +(sbar[ix].s2()-sbar[ix].s4())*
			      (RSsp[iy](iz,1)-RSsp[iy](iz,3)))
			+0.5*right*( (sbar[ix].s1()+sbar[ix].s3())*
				 (RSsp[iy](iz,0)+RSsp[iy](iz,2))
				 +(sbar[ix].s2()+sbar[ix].s4())*
				 (RSsp[iy](iz,1)+RSsp[iy](iz,3)));
		    }
		}
	    }
	  hadron[iya][ixa] = tvec+svec+scalar1*decay[0]->momentum()
	    +scalar2*inpart.momentum();
	}
    }
  // construct the lepton current
  int mode=(abs(decay[1]->id())-11)/12;
  Energy scale;
  ParticleVector leptons;
  leptons.push_back(decay[decay.size()-2]);
  leptons.push_back(decay[decay.size()-1]);
  vector<LorentzPolarizationVector> lepton(_current->current(vertex,mode,ichan,
							     scale,leptons));
  // work out the mapping for the lepton vector
  vector<int> ispin(decay.size()),ihel(decay.size()+1);
  for(int ix=int(decay.size()-1);ix>=0;--ix){ispin[ix]=decay[ix]->data().iSpin();}
  DecayMatrixElement newME(spinin,ispin);
  for(unsigned int iya=0;iya<4;++iya)
    {
      ihel[1]=iya-2;if(ihel[1]>=0){++ihel[1];}
      for(unsigned int ixa=0;ixa<2;++ixa)
	{
	  ihel[0]=2*ixa-1;
	  for(unsigned int lhel=0;lhel<lepton.size();++lhel)
	    {
	      ihel[2] = 2*int(lhel/2)-1;
	      ihel[3] = 2*int(lhel%2)-1;
	      newME(ihel) = lepton[lhel]*hadron[iya][ixa];
	    }
	}  
    }
  // store the matrix element
  ME(newME);
  // ckm factor
  double ckm(1.);
  if(inquark<=6)
    {
      if(inquark%2==0){ckm = SM().CKM(inquark/2-1,(abs(outquark)-1)/2);}
      else{ckm = SM().CKM(abs(outquark)/2-1,(inquark-1)/2);}
    }
  // return the answer
  return 0.5*(newME.contract(temp)).real()*_GF*_GF*ckm;
}

// output the setup information for the particle database
void SemiLeptonicBaryonDecayer::dataBaseOutput(ofstream & output)
{
  output << "update decayers set parameters=\"";
  output << "set " << fullName() << ":Iteration " << _niter << " \n";
  output << "set " << fullName() << ":Ntry " << _ntry << "\n";
  output << "set " << fullName() << ":GFermi "   << _GF*GeV2 << " \n";
  for(unsigned int ix=0;ix<_maxwgt.size();++ix)
    {output << "insert " << fullName() << ":MaximumWeight " << ix << " " 
	    << _maxwgt[ix] << " \n";}
  _current->dataBaseOutput(output);
  output << "set " << fullName() << ":Current " << _current->fullName() << " \n";
  _form->dataBaseOutput(output);
  output << "set " << fullName() << ":FormFactor " << _form->fullName() << " \n";
  output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}

}
