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
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/RSSpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/RSSpinorBarWaveFunction.h"
#include "ThePEG/Helicity/LorentzPolarizationVector.h"
#include "ThePEG/Helicity/FermionSpinInfo.h"
#include "ThePEG/Helicity/RSFermionSpinInfo.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;
using ThePEG::Helicity::RhoDMatrix;
using ThePEG::Helicity::LorentzPolarizationVector;
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
  tPDVector extpart,ptemp;
  _modemap.clear();
  double maxweight;
  vector<double> channelwgts(1,1.);
  int id0(0),id1(0),Wcharge(0),inspin,spect1,spect2,inquark,outquark,outspin;
  Energy min;
  DecayPhaseSpaceModePtr mode;
  DecayPhaseSpaceChannelPtr channel;
  int iq(0),ia(0);
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

bool SemiLeptonicBaryonDecayer::accept(tcPDPtr parent, const tPDVector & children) const {
  // find the non-lepton
  int ibar(0),idtemp,idin(parent->id());
  vector<int> idother; bool dummy;
  tPDVector::const_iterator pit  = children.begin();
  tPDVector::const_iterator pend = children.end();
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

int SemiLeptonicBaryonDecayer::modeNumber(bool & cc,tcPDPtr parent,
					  const tPDVector & children) const {
  // find the ids of the particles for the decay current
  tPDVector::const_iterator pit  = children.begin();
  tPDVector::const_iterator pend = children.end();
  int idtemp,ibar(0),idin(parent->id());
  vector<int> idother;
  cc=false;
  for( ; pit!=pend;++pit)
    {
      idtemp=(**pit).id();
      if(abs(idtemp)>16){ibar=idtemp;}
      else{idother.push_back(idtemp);}
    }
  return _modemap[_form->formFactorNumber(idin,ibar,cc)]
    +_current->decayMode(idother);
}


void SemiLeptonicBaryonDecayer::persistentOutput(PersistentOStream & os) const {
  os << _current << _form << _maxwgt << _modemap << ounit(_GF,1./GeV2);
}

void SemiLeptonicBaryonDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _current >> _form >> _maxwgt >> _modemap >> iunit(_GF,1./GeV2);
}

ClassDescription<SemiLeptonicBaryonDecayer> SemiLeptonicBaryonDecayer::initSemiLeptonicBaryonDecayer;
// Definition of the static class description member.

void SemiLeptonicBaryonDecayer::Init() {

  static ClassDocumentation<SemiLeptonicBaryonDecayer> documentation
    ("The SemiLeptonicBaryonDecayer class is designed for"
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
  RhoDMatrix temp;
  vector<LorentzSpinor<SqrtEnergy> > sp;
  vector<LorentzSpinorBar<SqrtEnergy> > sbar;
  if(inpart.id()>0)
    {
      SpinorWaveFunction(sp,temp,const_ptr_cast<tPPtr>(&inpart),incoming,true,vertex);
      SpinorBarWaveFunction(sbar,decay[0],outgoing,true,vertex);
    }
  else
    {
      SpinorBarWaveFunction(sbar,temp,const_ptr_cast<tPPtr>(&inpart),incoming,true,
			    vertex);
      SpinorWaveFunction(sp,decay[0],outgoing,true,vertex);
    }
  // get the information on the form-factor
  int spinin(0),spinout(0),spect1,spect2,inquark,outquark;
  int id0(inpart.id()),id1(decay[0]->id());
  bool cc; int iloc(_form->formFactorNumber(id0,id1,cc));
  _form->formFactorInfo(iloc,spinin,spinout,spect1,spect2,inquark,outquark);
  // work out the value of q and calculate the form factors
  Lorentz5Momentum q(inpart.momentum()-decay[0]->momentum());q.rescaleMass();
  Energy m0(inpart.mass()),m1(decay[0]->mass());
  Energy2 q2(q.mass2());
  Lorentz5Momentum sum(inpart.momentum()+decay[0]->momentum());
  // calculate the hadronic current for the decay
  vector<LorentzPolarizationVectorE> hadron;
  // calculate the form factors
  Complex f1v,f2v,f3v,f1a,f2a,f3a;
  _form->SpinHalfSpinHalfFormFactor(q2,iloc,id0,id1,m0,m1,
				    f1v,f2v,f3v,f1a,f2a,f3a);
  // now we need to construct the current
  LorentzPolarizationVectorE vtemp;
  hadron.resize(4);
  Complex left  =f1v-f1a-f2v-double((m0-m1)/(m0+m1))*f2a;
  Complex right =f1v+f1a-f2v+double((m0-m1)/(m0+m1))*f2a;
  for(unsigned int ix=0;ix<2;++ix)
    {
      for(unsigned int iy=0;iy<2;++iy)
	{
	  vtemp = sp[ix].generalCurrent(sbar[iy],left,right);
	  complex<Energy> vspin = sp[ix].scalar(sbar[iy]);
	  complex<Energy> aspin = sp[ix].pseudoScalar(sbar[iy]);
	  // the momentum like pieces
	  if(inpart.id()>0)
	    {
	      vtemp+= (f2v*vspin+f2a*aspin)/(m0+m1)*sum;
	      vtemp+= (f3v*vspin+f3a*aspin)/(m0+m1)*q;
	    }
	  else
	    {
	      vtemp-= (f2v*vspin-f2a*aspin)/(m0+m1)*sum;
	      vtemp+= (f3v*vspin-f3a*aspin)/(m0+m1)*q;
	    }
	  if(inpart.id()>0){hadron[2*ix+iy]=vtemp;}
	  else{hadron[2*iy+ix]=vtemp;}
	}
    }
  // construct the lepton current
  int mode((abs(decay[1]->id())-11)/12);
  Energy scale;
  ParticleVector leptons;
  leptons.push_back(decay[decay.size()-2]);
  leptons.push_back(decay[decay.size()-1]);
  vector<LorentzPolarizationVectorE> lepton(_current->current(vertex,mode,ichan,
							      scale,leptons));
  // work out the mapping for the lepton vector
  vector<unsigned int> constants(decay.size()+1), ihel(decay.size()+1);
  vector<PDT::Spin> ispin(decay.size());
  int itemp(1); unsigned int ibar(0);
  for(int ix=int(decay.size()-1);ix>=0;--ix)
    {
      ispin[ix]=decay[ix]->data().iSpin();
      if(abs(decay[ix]->id())<=16){itemp*=ispin[ix];constants[ix]=itemp;}
      else{ibar=ix;}
    }
  constants[decay.size()]=1;
  constants[ibar]=constants[ibar+1];
  DecayMatrixElement newME(PDT::Spin1Half,ispin);
  unsigned int mhel,ix,lhel;
  for(mhel=0;mhel<hadron.size();++mhel)
    {
      ihel[0     ]=mhel/spinout;
      ihel[ibar+1]=mhel%spinout;
      for(lhel=0;lhel<lepton.size();++lhel)
	{
	  // map the index for the leptons to a helicity state
	  for(ix=decay.size();ix>0;--ix)
	    {if(ix-1!=ibar){ihel[ix]=(lhel%constants[ix-1])/constants[ix];}}
	  newME(ihel)= lepton[lhel].dot(hadron[mhel])*_GF;
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
  return 0.5*(newME.contract(temp)).real()*ckm;
}

// matrix element for a 1/2 -> 3/2 semi-leptonic decay
double SemiLeptonicBaryonDecayer::halfThreeHalf(bool vertex, const int ichan,
						const Particle & inpart,
						const ParticleVector & decay) const
{
  // set up the spins and calculate the spinors
  RhoDMatrix temp;
  vector<LorentzSpinor<SqrtEnergy> > sp;
  vector<LorentzSpinorBar<SqrtEnergy> > sbar;
  vector<LorentzRSSpinor<SqrtEnergy> > RSsp;
  vector<LorentzRSSpinorBar<SqrtEnergy> > RSsbar;
  vector<LorentzPolarizationVector> eps;
  if(inpart.id()>0)
    {
      SpinorWaveFunction(sp,temp,const_ptr_cast<tPPtr>(&inpart),incoming,true,vertex);
      RSSpinorBarWaveFunction(RSsbar,decay[0],outgoing,true,vertex);
      sbar.resize(RSsbar.size());
      for(unsigned int ix=0;ix<RSsbar.size();++ix)
	{sbar[ix]=RSsbar[ix].dot(inpart.momentum());}
    }
  else
    {
      SpinorBarWaveFunction(sbar,temp,const_ptr_cast<tPPtr>(&inpart),incoming,true,
			    vertex);
      RSSpinorWaveFunction(RSsp,decay[0],outgoing,true,vertex);
      sp.resize(RSsp.size());
      for(unsigned int ix=0;ix<RSsp.size();++ix)
	{sp[ix]=RSsp[ix].dot(inpart.momentum());}
    }
  // get the information on the form-factor
  int spinin(0),spinout(0),inquark,outquark,spect1,spect2;
  int id0(inpart.id()),id1(decay[0]->id());
  bool cc; int iloc(_form->formFactorNumber(id0,id1,cc));
  _form->formFactorInfo(iloc,spinin,spinout,spect1,spect2,inquark,outquark);
  // work out the value of q and calculate the form factors
  Lorentz5Momentum q(inpart.momentum()-decay[0]->momentum());q.rescaleMass();
  Energy m0(inpart.mass()),m1(decay[0]->mass());
  Energy2 q2(q.mass2());
  Lorentz5Momentum sum(inpart.momentum()+decay[0]->momentum());
  // calculate the hadronic current for the decay
  LorentzPolarizationVectorE hadron[4][2];
  // calculate the form factors
  Complex f1v,f2v,f3v,f4v,f1a,f2a,f3a,f4a;
  _form->SpinHalfSpinThreeHalfFormFactor(q2,iloc,id0,id1,m0,m1,
					 f1v,f2v,f3v,f4v,f1a,f2a,f3a,f4a);
  // now we need to construct the current
  LorentzPolarizationVector vtemp;
  complex<InvEnergy2> lS1,lS2,rS1,rS2;
  complex<InvEnergy> lV,rV;
  Complex left,right;
  InvEnergy ms(1./(m0+m1));
  InvEnergy2 ms2(ms*ms);
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
  unsigned int ix,iy;
  // construct the vectors for the decay
  Complex scalar1,scalar2;
  complex<Energy> lfact,rfact;
  LorentzPolarizationVectorE tvec;
  LorentzPolarizationVector svec;
  for(unsigned int iya=0;iya<4;++iya)
    {
      for(unsigned int ixa=0;ixa<2;++ixa)
	{
	  if(decay[0]->id()>0){ix=iya;iy=ixa;}
	  else{ix=ixa;iy=iya;}
	  // scalar like terms
	  lfact = sp[iy].leftScalar(sbar[ix]);
	  rfact = sp[iy].rightScalar(sbar[ix]);
	  scalar1 = (lS1*lfact+rS1*rfact)*UnitRemoval::E;
	  scalar2 = (lS2*lfact+rS2*rfact)*UnitRemoval::E;
	  svec = sp[iy].generalCurrent(sbar[ix],lV/ms,rV/ms)*ms;
	  if(inpart.id()>0) tvec=RSsbar[ix].generalCurrent(sp[iy],left,right);
	  else              tvec=RSsp[iy].generalCurrent(sbar[ix],left,right);
 	  hadron[iya][ixa] = tvec+svec*UnitRemoval::E+scalar1*decay[0]->momentum()
	    +scalar2*inpart.momentum();
	}
    }
  // construct the lepton current
  int mode=(abs(decay[1]->id())-11)/12;
  Energy scale;
  ParticleVector leptons;
  leptons.push_back(decay[decay.size()-2]);
  leptons.push_back(decay[decay.size()-1]);
  vector<LorentzPolarizationVectorE> lepton(_current->current(vertex,mode,ichan,
							     scale,leptons));
  // work out the mapping for the lepton vector
  vector<unsigned int> ihel(decay.size()+1);
  vector<PDT::Spin> ispin(decay.size());
  for(int ii=int(decay.size()-1);ii>=0;--ii){ispin[ii]=decay[ii]->data().iSpin();}
  DecayMatrixElement newME(PDT::Spin1Half,ispin);
  for(unsigned int iya=0;iya<4;++iya)
    {
      ihel[1]=iya;
      for(unsigned int ixa=0;ixa<2;++ixa)
	{
	  ihel[0]=ixa;
	  for(unsigned int lhel=0;lhel<lepton.size();++lhel)
	    {
	      ihel[2] = lhel/2;
	      ihel[3] = lhel%2;
	      newME(ihel) = lepton[lhel].dot(hadron[iya][ixa])*_GF;
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
  return 0.5*(newME.contract(temp)).real()*ckm;
}

// output the setup information for the particle database
void SemiLeptonicBaryonDecayer::dataBaseOutput(ofstream & output,bool header) const
{
  if(header){output << "update decayers set parameters=\"";}
  DecayIntegrator::dataBaseOutput(output,false);
  output << "set " << fullName() << ":GFermi "   << _GF*GeV2 << " \n";
  for(unsigned int ix=0;ix<_maxwgt.size();++ix)
    {output << "insert " << fullName() << ":MaximumWeight " << ix << " " 
	    << _maxwgt[ix] << " \n";}
  _current->dataBaseOutput(output,false,true);
  output << "set " << fullName() << ":Current " << _current->fullName() << " \n";
  _form->dataBaseOutput(output,false,true);
  output << "set " << fullName() << ":FormFactor " << _form->fullName() << " \n";
  if(header){output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;}
}

}
