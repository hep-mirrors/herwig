// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BaryonFactorizedDecayer class.
//

#include "BaryonFactorizedDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/PDT/DecayMode.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "BaryonFactorizedDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Parameter.h"
#include "Herwig++/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/RSSpinorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/RSSpinorBarWaveFunction.h"
#include "Herwig++/Helicity/Correlations/DecayVertex.h"
#include "ThePEG/Helicity/FermionSpinInfo.h"
#include "ThePEG/Helicity/RSFermionSpinInfo.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;
using namespace Herwig::Helicity;

BaryonFactorizedDecayer::~BaryonFactorizedDecayer() {}

void BaryonFactorizedDecayer::doinit() throw(InitException) {
  DecayIntegrator::doinit();
  // make sure the current and form factor got initialised
  _current->init();
  _form->init();
  // set up the phase-space channels
  DecayPhaseSpaceModePtr mode,modeb;
  DecayPhaseSpaceChannelPtr channel,channelb;
  PDVector extpart(2),extpartb(2),ptemp;
  double maxweight;
  vector<double> channelwgts;
  _formmap.resize(0);
  _currentmap.resize(0);
  int iq,ia,spect1,spect2,inquark,outquark,id0,id1,ispin,ospin,icharge;
  unsigned int iy,iform,icurrent;
  Energy mmax;
  bool done,doneb;
  vector<double>::iterator start,end;
  double ckm;
  for(iform=0;iform<_form->numberOfFactors();++iform)
    {
      // get the particles nad information for the form factor
      _form->particleID (iform,id0,id1);
      extpart[0]=getParticleData(id0);extpartb[0]=extpart[0];
      extpart[1]=getParticleData(id1);extpartb[1]=extpart[1];
      _form->formFactorInfo(iform,ispin,ospin,spect1,spect2,inquark,outquark);
      // maximum mass for the particles in the current
      mmax = extpart[0]->mass()+extpart[0]->widthUpCut()
	-(extpart[1]->mass()-extpart[1]->widthUpCut());
      // the charge of the decay products
      icharge = extpart[0]->iCharge()-extpart[1]->iCharge();
      // loop over the modes in the current
      for(icurrent=0;icurrent<_current->numberOfModes();++icurrent)
	{
	  done=false;doneb=false;
	  _current->decayModeInfo(icurrent,iq,ia);
	  extpart.resize(2);
	  extpartb.resize(2);
	  // get the external particles for this mode
	  ptemp=_current->particles(icharge,icurrent,iq,ia);
	  for(iy=0;iy<ptemp.size();++iy){extpart.push_back(ptemp[iy]);}
	  // create the mode
	  mode=new_ptr(DecayPhaseSpaceMode(extpart,this));
	  // create the first piece of the channel
	  channel = new_ptr(DecayPhaseSpaceChannel(mode));
	  channel->addIntermediate(extpart[0],0,0.0,-1,1);
	  done=_current->createMode(icharge,icurrent,mode,2,1,channel,mmax);
	  if(icharge==0&&iq!=-ia)
	    {
	      // create the CC mode
	      // get the external particles for this mode
	      ptemp=_current->particles(icharge,icurrent,-ia,-iq);
	      for(iy=0;iy<ptemp.size();++iy){extpartb.push_back(ptemp[iy]);}
	      // create the mode
	      modeb=new_ptr(DecayPhaseSpaceMode(extpartb,this));
	      // create the first piece of the channel
	      channelb = new_ptr(DecayPhaseSpaceChannel(modeb));
	      channelb->addIntermediate(extpartb[0],0,0.0,-1,1);
	      doneb=_current->createMode(icharge,icurrent,modeb,2,1,channelb,mmax);
	    }
	  if(done)
	    {
	      // the maximum weight and the channel weights
	      // the maximum
	      if(_wgtmax.size()>numberModes()){maxweight=_wgtmax[numberModes()];}
	      else{maxweight=0.;}
	      // the weights for the channel
	      if(_wgtloc.size()>numberModes()&&
		 _wgtloc[numberModes()]+mode->numberChannels()<=_weights.size())
		{
		  start=_weights.begin()+_wgtloc[numberModes()];
		  end  =start+mode->numberChannels();
		  channelwgts=vector<double>(start,end);
		}
	      else
		{channelwgts.resize(mode->numberChannels(),1./(mode->numberChannels()));}
	      _formmap.push_back(iform);
	      _currentmap.push_back(icurrent);
	      // special for the two body modes
	      if(extpart.size()==3)
		{
		  channelwgts.resize(0);
		  mode=new_ptr(DecayPhaseSpaceMode(extpart,this));
		}
	      addMode(mode,maxweight,channelwgts);
	      if(abs(icharge)==3){_modetype.push_back(0);}
	      else{_modetype.push_back(1);}
	      // the ckm factor
	      ckm=1.;
	      if(abs(icharge)==3)
		{
		  // CKM factor for the baryon transition
		  if(abs(inquark)%2==0){ckm*=SM().CKM(abs(inquark )/2-1,
						      (abs(outquark)-1)/2);}
		  else                 {ckm*=SM().CKM(abs(outquark)/2-1,
						      (abs(inquark )-1)/2);}
		  // CKM factor for the current transition
		  if(iq%2==0){ckm*= SM().CKM(abs(iq)/2-1,(abs(ia)-1)/2);}
		  else       {ckm*= SM().CKM(abs(ia)/2-1,(abs(iq)-1)/2);}
		}
	      else
		{
		  if(abs(inquark)%2==0)
		    {
		      ckm*=SM().CKM(abs(inquark )/2-1,(abs(iq)-1)/2);
		      ckm*=SM().CKM(abs(outquark)/2-1,(abs(ia)-1)/2);
		    }
		  else
		    {
		      ckm*=SM().CKM(abs(iq)/2-1,(abs(inquark )-1)/2);
		      ckm*=SM().CKM(abs(ia)/2-1,(abs(outquark)-1)/2);
		    }
		}
	      _CKMfact.push_back(ckm);
	    }
	  if(doneb)
	    {
	      // the maximum weight and the channel weights
	      // the maximum
	      if(_wgtmax.size()>numberModes()){maxweight=_wgtmax[numberModes()];}
	      else{maxweight=0.;}
	      // the weights for the channel
	      if(_wgtloc.size()>numberModes()&&
		 _wgtloc[numberModes()]+modeb->numberChannels()<=_weights.size())
		{
		  start=_weights.begin()+_wgtloc[numberModes()];
		  end  =start+modeb->numberChannels();
		  channelwgts=vector<double>(start,end);
		}
	      else
		{channelwgts.resize(modeb->numberChannels(),1./(modeb->numberChannels()));}
	      _formmap.push_back(iform);
	      _currentmap.push_back(icurrent);
	      // special for the two body modes
	      if(extpart.size()==3)
		{
		  channelwgts.resize(0);
		  modeb=new_ptr(DecayPhaseSpaceMode(extpartb,this));
		}
	      addMode(modeb,maxweight,channelwgts);
	      if(abs(icharge)==3){_modetype.push_back(0);}
	      else{_modetype.push_back(1);}
	      if(abs(inquark)%2==0)
		{
		  ckm*=SM().CKM(abs(inquark )/2-1,(abs(ia)-1)/2);
		  ckm*=SM().CKM(abs(outquark)/2-1,(abs(iq)-1)/2);
		}
	      else
		{
		  ckm*=SM().CKM(abs(ia)/2-1,(abs(inquark )-1)/2);
		  ckm*=SM().CKM(abs(iq)/2-1,(abs(outquark)-1)/2);
		}
	      _CKMfact.push_back(ckm);
	    }
	}
    }
}

bool BaryonFactorizedDecayer::accept(const DecayMode & dm) const {
  bool allowed=false;
  unsigned int iform(0),ix;
  int idin(dm.parent()->id()),ibaryon,foundb,id0,id1;
  vector<int> idall,idother;
  ParticleMSet::const_iterator pit  = dm.products().begin();
  ParticleMSet::const_iterator pend = dm.products().end();
  for( ; pit!=pend;++pit){idall.push_back((**pit).id());}
  // loop over the particles in the form factor
  do
    {
      _form->particleID(iform,id0,id1);
      ibaryon=0;
      if(id0==idin){ibaryon=id1;}
      else if(id0==-idin){ibaryon=-id1;}
      if(ibaryon!=0)
	{
	  foundb=false;
	  idother.resize(0);
	  for(ix=0;ix<idall.size();++ix)
	    {
	      if(idall[ix]==ibaryon){foundb=true;}
	      else{idother.push_back(idall[ix]);}
	    }
	  if(foundb){allowed=_current->accept(idother);}
	}
      ++iform;
    }
  while(!allowed&&iform<_form->numberOfFactors());
  return allowed;
}

ParticleVector BaryonFactorizedDecayer::decay(const DecayMode & dm,
					      const Particle & parent) const {
  unsigned int ix;
  int idin(parent.id()),ibaryon,foundb,id0,id1,icurrent(-1),iform(0);
  vector<int> idall,idother;
  ParticleMSet::const_iterator pit  = dm.products().begin();
  ParticleMSet::const_iterator pend = dm.products().end();
  for( ; pit!=pend;++pit){idall.push_back((**pit).id());}
  // loop over the particles in the form factor
  do
    {
      _form->particleID(iform,id0,id1);
      ibaryon=0;
      if(id0==idin){ibaryon=id1;}
      else if(id0==-idin){ibaryon=-id1;}
      ++iform;
      foundb=false;
      idother.resize(0);
      for(ix=0;ix<idall.size();++ix)
	{
	  if(idall[ix]==ibaryon){foundb=true;}
	  else{idother.push_back(idall[ix]);}
	}
      if(foundb){icurrent=_current->decayMode(idother);}
    }
  while(icurrent<0&&iform<int(_form->numberOfFactors()));
  // now find the mode
  int imode=-1;;
  do
    {
      if(_currentmap[ix]==icurrent&&_formmap[ix]==iform){imode=ix;}
      ++ix;
    }
  while(imode<0&&ix<numberModes());
  // generate the mode
  bool cc(id0!=idin);
  return generate(true,cc,imode,parent);
}


void BaryonFactorizedDecayer::persistentOutput(PersistentOStream & os) const {
  os << _current << _form << _GF << _a1b << _a2b <<_a1c <<_a2c << _currentmap 
     << _formmap << _modetype << _CKMfact << _wgtloc << _wgtmax << _weights;
}

void BaryonFactorizedDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _current >> _form >> _GF >> _a1b >> _a2b >>_a1c >>_a2c >> _currentmap 
     >> _formmap >> _modetype >> _CKMfact >> _wgtloc >> _wgtmax >> _weights;
}

ClassDescription<BaryonFactorizedDecayer> BaryonFactorizedDecayer::initBaryonFactorizedDecayer;
// Definition of the static class description member.

void BaryonFactorizedDecayer::Init() {

  static ClassDocumentation<BaryonFactorizedDecayer> documentation
    ("The BaryonFactorizedDecayer class combines the baryon form factor and a"
     " weak current to perform a decay in the naive factorization approximation.");

  static Parameter<BaryonFactorizedDecayer,InvEnergy2> interfaceGFermi
    ("GFermi",
     "The Fermi coupling constant",
     &BaryonFactorizedDecayer::_GF, 1./GeV2, 1.16639E-5/GeV2, -1.0e12*1./GeV2, 1.0e12*1./GeV2,
     false, false, false);

  static Reference<BaryonFactorizedDecayer,WeakDecayCurrent> interfaceWeakCurrent
    ("Current",
     "The reference for the decay current to be used.",
     &BaryonFactorizedDecayer::_current, false, false, true, false, false);

  static ParVector<BaryonFactorizedDecayer,int> interfaceWeightLocation
    ("WeightLocation",
     "The locations of the weights for a given channel in the vector",
     &BaryonFactorizedDecayer::_wgtloc,
     0, 0, 0, 0, 10000, false, false, true);

  static ParVector<BaryonFactorizedDecayer,double> interfaceWeightMax
    ("MaximumWeight",
     "The maximum weight for a given channel.",
     &BaryonFactorizedDecayer::_wgtmax,
     0, 0, 0, 0., 10000., false, false, true);

  static ParVector<BaryonFactorizedDecayer,double> interfaceWeights
    ("Weights",
     "The weights for the integration.",
     &BaryonFactorizedDecayer::_weights,
     0, 0, 0, 0., 10000., false, false, true);

  static Reference<BaryonFactorizedDecayer,BaryonFormFactor> interfaceFormFactor
    ("FormFactor",
     "The form-factor",
     &BaryonFactorizedDecayer::_form, true, true, true, false, false);
  static Parameter<BaryonFactorizedDecayer,double> interfacea1Bottom
    ("a1Bottom",
     "The factorization paramter a_1 for decays of bottom baryons",
     &BaryonFactorizedDecayer::_a1b, 1.23, -10.0, 10.0,
     false, false, true);

  static Parameter<BaryonFactorizedDecayer,double> interfacea2Bottom
    ("a2Bottom",
     "The factorization paramter a_2 for decays of bottom baryons",
     &BaryonFactorizedDecayer::_a2b, 0.33, -10.0, 10.0,
     false, false, true);

  static Parameter<BaryonFactorizedDecayer,double> interfacea1Charm
    ("a1Charm",
     "The factorization paramter a_1 for decays of charm baryons",
     &BaryonFactorizedDecayer::_a1c, 1.1, -10.0, 10.0,
     false, false, true);

  static Parameter<BaryonFactorizedDecayer,double> interfacea2Charm
    ("a2Charm",
     "The factorization paramter a_2 for decays of charm baryons",
     &BaryonFactorizedDecayer::_a2c, -0.5, -10.0, 10.0,
     false, false, true);
}

double BaryonFactorizedDecayer::me2(bool vertex, const int ichan,
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
				  << " in BaryonFactorizedDecayer::me2()" 
				  << Exception::abortnow;}
    }
  else
    {throw DecayIntegratorError() << "Invalid spins "
				  << " in BaryonFactorizedDecayer::me2()" 
				  << Exception::abortnow;}
  return me;
}

// matrix element for a 1/2 -> 1/2 decay
double BaryonFactorizedDecayer::halfHalf(bool vertex, const int ichan,
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
				  << " in BaryonFactorizedDecayer::halfHalf()" 
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
  _form->formFactorInfo(_formmap[imode()],spinin,spinout,spect1,spect2,inquark,outquark);
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
  // calculate the baryon part of the current for the decay
  vector<LorentzPolarizationVector> baryon;
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
  _form->SpinHalfSpinHalfFormFactor(q2,_formmap[imode()],id0,id1,m0,m1,
				    f1v,f2v,f3v,f1a,f2a,f3a);
  // now we need to construct the current
  LorentzPolarizationVector vtemp;
  baryon.resize(4);
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
	  if(inpart.id()>0){baryon[2*ix+iy]=vtemp;}
	  else{baryon[2*iy+ix]=vtemp;}
	}
    }
  // construct the weak current
  Energy scale;
  ParticleVector::const_iterator start=decay.begin()+1;
  ParticleVector::const_iterator end  =decay.end();
  ParticleVector hadpart(start,end);
  vector<LorentzPolarizationVector> hadron(_current->current(vertex,_currentmap[imode()],
							     ichan,scale,hadpart));
  // prefactor
  double pre = pow(inpart.mass()/scale,int(hadpart.size()-2));pre*=pre;
  // work out the mapping for the hadron vector
  vector<int> constants(decay.size()+1), ispin(decay.size()),ihel(decay.size()+1);
  int itemp=1; unsigned int ibar=0;
  for(int ix=int(decay.size()-1);ix>=0;--ix)
    {
      ispin[ix]=decay[ix]->data().iSpin();
      if(abs(decay[ix]->id())!=id1){itemp*=ispin[ix];constants[ix]=itemp;}
      else{ibar=ix;}
    }
  constants[decay.size()]=1;
  constants[ibar]=constants[ibar+1];
  DecayMatrixElement newME(spinin,ispin);
  for(unsigned int mhel=0;mhel<baryon.size();++mhel)
    {
      ihel[0     ]=2*(mhel/spinout)-spinin/2;
      ihel[ibar+1]=2*(mhel%spinout)-spinout/2;
      for(unsigned int lhel=0;lhel<hadron.size();++lhel)
	{
	  // map the index for the hadrons to a helicity state
	  for(unsigned int ix=decay.size();ix>0;--ix)
	    {
	      if(ix-1!=ibar)
		{
		  ihel[ix]=(lhel%constants[ix-1])/constants[ix]-int(ispin[ix-1]/2);
		  if(ispin[ix-1]%2==0&&ihel[ix]>=0&&ispin[ix-1]!=0){++ihel[ix];}
		}
	    }
	  newME(ihel)= hadron[lhel]*baryon[mhel];
	}
    }  
  // store the matrix element
  ME(newME);
  // return the answer
  return 0.5*pre*(newME.contract(temp)).real()*_GF*_GF*_CKMfact[imode()];
}

// matrix element for a 1/2 -> 3/2 decay
double BaryonFactorizedDecayer::halfThreeHalf(bool vertex, const int ichan,
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
  _form->formFactorInfo(_formmap[imode()],spinin,spinout,spect1,spect2,inquark,outquark);
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
  // calculate the baryonic current for the decay
  LorentzPolarizationVector baryon[4][2];
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
  _form->SpinHalfSpinThreeHalfFormFactor(q2,_formmap[imode()],id0,id1,m0,m1,
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
	  baryon[iya][ixa] = tvec+svec+scalar1*decay[0]->momentum()
	    +scalar2*inpart.momentum();
	}
    }
  // construct the hadron current
  Energy scale;
  ParticleVector::const_iterator start=decay.begin()+1;
  ParticleVector::const_iterator end  =decay.end();
  ParticleVector hadpart(start,end);
  vector<LorentzPolarizationVector> hadron(_current->current(vertex,_currentmap[imode()],
							     ichan,scale,hadpart));
  // prefactor
  double pre = pow(inpart.mass()/scale,int(hadpart.size()-2));pre*=pre;
  // work out the mapping for the hadron vector
  vector<int> constants(decay.size()+1), ispin(decay.size()),ihel(decay.size()+1);
  int itemp=1; unsigned int ibar=0;
  for(int ix=int(decay.size()-1);ix>=0;--ix)
    {
      ispin[ix]=decay[ix]->data().iSpin();
      if(abs(decay[ix]->id())!=id1){itemp*=ispin[ix];constants[ix]=itemp;}
      else{ibar=ix;}
    }
  constants[decay.size()]=1;
  constants[ibar]=constants[ibar+1];
  DecayMatrixElement newME(spinin,ispin);
  for(unsigned int iya=0;iya<4;++iya)
    {
      ihel[1]=iya-2;if(ihel[1]>=0){++ihel[1];}
      for(unsigned int ixa=0;ixa<2;++ixa)
	{
	  ihel[0]=2*ixa-1;
	  for(unsigned int lhel=0;lhel<hadron.size();++lhel)
	    {
	      // map the index for the hadrons to a helicity state
	      for(unsigned int ix=decay.size();ix>0;--ix)
		{
		  if(ix-1!=ibar)
		    {
		      ihel[ix]=(lhel%constants[ix-1])/constants[ix]-int(ispin[ix-1]/2);
		      if(ispin[ix-1]%2==0&&ihel[ix]>=0&&ispin[ix-1]!=0){++ihel[ix];}
		    }
		}
	      newME(ihel) = hadron[lhel]*baryon[iya][ixa];
	    }
	}
    }
  // store the matrix element
  ME(newME);
  // return the answer
  double me=0.5*pre*(newME.contract(temp)).real()*_GF*_GF*_CKMfact[imode()];
  // testing code
  /*
  Energy m2(decay[1]->mass());
  Energy ef=0.5/m0*(m0*m0-m2*m2+m1*m1);
  complex<Energy> fP= inpart.mass()*hadron[0][0]/decay[1]->momentum()[0];
  Energy msum(m0+m1);
  //    cout << "testing A " 
  //     << f1v << " " << f2v << " " << f3v << " " << f4v 
  //     << f1a << " " << f2a << " " << f3a << " " << f4a << endl; 
  Complex C =-fP*(f1a+(m0-m1)*f2a/msum+(m0*ef-m1*m1)*f3a/msum/msum);
  Complex D = fP*(f1v-(m0+m1)*f2v/msum+(m0*ef-m1*m1)*f3v/msum/msum);
  // calculate the matrix element using the KK results
  Energy2 qplus  = (m0+m1)*(m0+m1)-m2*m2;
  Energy2 qminus = (m0-m1)*(m0-m1)-m2*m2;
  Energy rqplus=sqrt(qplus),rqminus=sqrt(qminus); 
  Complex h1,h2;
  Energy pcm = sqrt((m0*m0-(m1+m2)*(m1+m2))*(m0*m0-(m1-m2)*(m1-m2)))/2./m0;
  // the amplitudes
  h1 =-2.*sqrt(2./3.)*pcm*m0/m1*rqminus*D;
  h2 = 2.*sqrt(2./3.)*pcm*m0/m1*rqplus *C;
  double kappa = sqrt((ef-m1)/(ef+m1));
  cout << "testing " << C << " " << D << endl;
  cout << "testing alpha " << 2.*kappa*(conj(C)*D).real()/((kappa*kappa*C*conj(C)+D*conj(D)).real()) << endl;
  cout << "testing C " 
       << 0.5*_GF*_GF*_CKMfact[imode()]*0.25*(h1*conj(h1)+h2*conj(h2))/(me*m0*m0) 
       << endl;
  */
  return me;
}


// output the setup information for the particle database
void BaryonFactorizedDecayer::dataBaseOutput(ofstream & output)
{
  /*
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
  */
}

}
