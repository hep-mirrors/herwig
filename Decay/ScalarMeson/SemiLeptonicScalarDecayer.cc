// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SemiLeptonicScalarDecayer class.
//

#include "SemiLeptonicScalarDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "SemiLeptonicScalarDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/LorentzPolarizationVector.h"
#include "Herwig++/Helicity/EpsFunction.h"
#include "ThePEG/Helicity/LorentzTensor.h"
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/TensorWaveFunction.h"

namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;
using ThePEG::Helicity::RhoDMatrix;
using ThePEG::Helicity::LorentzPolarizationVector;
using Helicity::ScalarWaveFunction;
using Helicity::VectorWaveFunction;
using Helicity::TensorWaveFunction;
using Helicity::EpsFunction;
using Helicity::Direction;
using Helicity::incoming;
using Helicity::outgoing;

void SemiLeptonicScalarDecayer::doinit() throw(InitException) {
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
  int id0(0),id1(0),Wcharge(0);
  Energy min;
  DecayPhaseSpaceModePtr mode;
  DecayPhaseSpaceChannelPtr channel;
  int iq,ia; unsigned int ix,iy,iz;
  bool done;
  for(ix=0;ix<_form->numberOfFactors();++ix)
    {
      // get the external particles for this mode
      extpart.resize(2);
      _form->particleID(ix,id0,id1);
      extpart[0]=getParticleData(id0);
      extpart[1]=getParticleData(id1);
      Wcharge =(extpart[0]->iCharge()-extpart[1]->iCharge());
      min = extpart[0]->mass()+extpart[0]->widthUpCut()
	-extpart[1]->mass()+extpart[1]->widthLoCut();
      _modemap.push_back(numberModes());
      for(iy=0;iy<_current->numberOfModes();++iy)
	{
	  extpart.resize(2); 	
	  _current->decayModeInfo(iy,iq,ia);
	  ptemp=_current->particles(Wcharge,iy,iq,ia);
	  for(iz=0;iz<ptemp.size();++iz){extpart.push_back(ptemp[iz]);}
	  // create the mode
	  mode=new_ptr(DecayPhaseSpaceMode(extpart,this));
	  // create the first piece of the channel
	  channel = new_ptr(DecayPhaseSpaceChannel(mode));
	  channel->addIntermediate(extpart[0],0,0.0,1,-1);
	  done=_current->createMode(Wcharge,iy,mode,2,1,channel,min);
	  if(done)
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
  
SemiLeptonicScalarDecayer::~SemiLeptonicScalarDecayer() {}

bool SemiLeptonicScalarDecayer::accept(const DecayMode & dm) const {
  // find the non-lepton
  int imes(0),idtemp,idin(dm.parent()->id());
  vector<int> idother; bool dummy;
  ParticleMSet::const_iterator pit  = dm.products().begin();
  ParticleMSet::const_iterator pend = dm.products().end();
  for( ; pit!=pend;++pit)
    {
      idtemp=(**pit).id();
      if(abs(idtemp)>16){imes=idtemp;}
      else{idother.push_back(idtemp);}
    }
  // check that the form factor exists
  if(_form->formFactorNumber(idin,imes,dummy)<0){return false;}
  // and the current
  return _current->accept(idother);
}

ParticleVector SemiLeptonicScalarDecayer::decay(const DecayMode & dm,
				  const Particle & parent) const {
  // find the ids of the particles for the decay current
  ParticleMSet::const_iterator pit = dm.products().begin();
  ParticleMSet::const_iterator pend = dm.products().end();
  int idtemp,imes(0),idin(dm.parent()->id());
  vector<int> idother;
  bool cc(false);
  for( ; pit!=pend;++pit)
    {
      idtemp=(**pit).id();
      if(abs(idtemp)>16){imes=idtemp;}
      else{idother.push_back(idtemp);}
    }
  int imode = _modemap[_form->formFactorNumber(idin,imes,cc)]
    +_current->decayMode(idother);
  // perform the decay
  return generate(true,cc,imode,parent);
}


void SemiLeptonicScalarDecayer::persistentOutput(PersistentOStream & os) const {
  os << _current << _form << _maxwgt << _modemap << _GF;
}

void SemiLeptonicScalarDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _current >> _form >> _maxwgt >> _modemap >> _GF;
}

ClassDescription<SemiLeptonicScalarDecayer> SemiLeptonicScalarDecayer::initSemiLeptonicScalarDecayer;
// Definition of the static class description member.

void SemiLeptonicScalarDecayer::Init() {

  static ClassDocumentation<SemiLeptonicScalarDecayer> documentation
    ("The SemiLeptonicScalarDecayer class is designed for the"
    "semi-leptonic decay of a (pseudo)-scalar meson.");

  static Parameter<SemiLeptonicScalarDecayer,InvEnergy2> interfaceGFermi
    ("GFermi",
     "The Fermi coupling constant",
     &SemiLeptonicScalarDecayer::_GF, 1./GeV2, 1.16639E-5/GeV2,
     0./GeV2, 1.e-4/GeV2,
     false, false, false);

  static Reference<SemiLeptonicScalarDecayer,LeptonNeutrinoCurrent> interfaceCurrent
    ("Current",
     "The current for the leptons produced in the decay.",
     &SemiLeptonicScalarDecayer::_current, true, true, true, false, false);

  static Reference<SemiLeptonicScalarDecayer,ScalarFormFactor> interfaceFormFactor
    ("FormFactor",
     "The form factor",
     &SemiLeptonicScalarDecayer::_form, true, true, true, false, false);

  static ParVector<SemiLeptonicScalarDecayer,double> interfaceMaximumWeight
    ("MaximumWeight",
     "The maximum weights for the decays",
     &SemiLeptonicScalarDecayer::_maxwgt,
     0, 0, 0, 0, 100., false, false, true);

}

// combine the currents and form-factors to give the matrix element
double SemiLeptonicScalarDecayer::me2(bool vertex, const int ichan,
				      const Particle & inpart,
				      const ParticleVector & decay) const
{
  // workaround for gcc 3.2.3 bug
  // spin info for the decaying particle
  //ALB ScalarWaveFunction(const_ptr_cast<tPPtr>(&inpart),incoming,true,vertex);
  tPPtr mytempInpart = const_ptr_cast<tPPtr>(&inpart);
  ScalarWaveFunction(mytempInpart,incoming,true,vertex);

  // get the information on the form-factor
  int jspin(0),id0(inpart.id()),id1(decay[0]->id());
  bool cc;
  unsigned int iloc(_form->formFactorNumber(id0,id1,cc));
  int spect,iq,ia;
  _form->formFactorInfo(iloc,jspin,spect,iq,ia);
  // work out the value of q and calculate the form factors
  Lorentz5Momentum q(inpart.momentum()-decay[0]->momentum());q.rescaleMass();
  Energy2 q2(q.mass2());
  Lorentz5Momentum sum(inpart.momentum()+decay[0]->momentum());
  // calculate the hadronic current for the decay
  Complex ii(0.,1.);
  vector<LorentzPolarizationVector> hadron;
  if(jspin==0)
    {
      // workaround for gcc 3.2.3 bug
      //ALB ScalarWaveFunction(decay[0],outgoing,true,vertex);
      PPtr mytemp = decay[0];
      ScalarWaveFunction(mytemp,outgoing,true,vertex);

      Complex fp,f0;
      _form->ScalarScalarFormFactor(q2,iloc,id0,id1,inpart.mass(),decay[0]->mass(),
				    f0,fp);
      double pre((inpart.mass()*inpart.mass()-decay[0]->mass()*decay[0]->mass())/q2);
      hadron.push_back(fp*sum+pre*(f0-fp)*q);
    }
  else if(jspin==1)
    {
      vector<LorentzPolarizationVector> vwave;
      VectorWaveFunction(vwave,decay[0],outgoing,true,false,vertex);
      Complex A0,A1,A2,A3,V,dot;
      Energy MP(inpart.mass()),MV(decay[0]->mass()),msum(MP+MV),mdiff(MP-MV);
      _form->ScalarVectorFormFactor(q2,iloc,id0,id1,MP,MV,A0,A1,A2,V);
      A3 = 0.5/MV*(msum*A1-mdiff*A2);
      if(cc){V=-V;}
      // compute the hadron currents
      for(unsigned int ix=0;ix<3;++ix)
	{
	  // dot product
	  dot = vwave[ix]*inpart.momentum();
	  // current
	  hadron.push_back(-ii*msum*A1*vwave[ix]
			   +ii*A2/msum*dot*sum
			   +2.*ii*MV/q2*(A3-A0)*dot*q
			   +2.*V/msum*EpsFunction::product(vwave[ix],inpart.momentum(),
							   decay[0]->momentum()));
	}
    }
  else if(jspin==2)
    {
      vector<LorentzTensor> twave;
      TensorWaveFunction(twave,decay[0],outgoing,true,false,vertex);
      Complex h,k,bp,bm,dot;
      _form->ScalarTensorFormFactor(q2,iloc,id0,id1,inpart.mass(),decay[0]->mass(),
				    h,k,bp,bm);
      if(cc){h=-h;}
      LorentzPolarizationVector dotv;
      // compute the hadron currents
      for(unsigned int ix=0;ix<5;++ix)
	{
	  dotv = twave[ix]*inpart.momentum();
	  dot = dotv*inpart.momentum();
	  hadron.push_back(ii*h*EpsFunction::product(dotv,sum,q)
			   -k*dotv-bp*dot*sum-bm*dot*q);
	}
    }
  int mode=(abs(decay[1]->id())-11)/2;
  // construct the lepton current
  Energy scale;
  ParticleVector leptons;
  leptons.push_back(decay[decay.size()-2]);
  leptons.push_back(decay[decay.size()-1]);
  vector<LorentzPolarizationVector> lepton(_current->current(vertex,mode,
							     ichan,scale,leptons));
  // work out the mapping for the lepton vector
  vector<unsigned int> constants(decay.size()+1),ihel(decay.size()+1);
  vector<PDT::Spin> ispin(decay.size());
  unsigned int itemp(1),imes(0);
  for(int ix=int(decay.size()-1);ix>=0;--ix)
    {
      ispin[ix]=decay[ix]->data().iSpin();
      if(abs(decay[ix]->id())<=16){itemp*=ispin[ix];constants[ix]=itemp;}
      else{imes=ix;}
    }
  constants[decay.size()]=1;
  constants[imes]=constants[imes+1];
  DecayMatrixElement newME(PDT::Spin0,ispin);
  for(unsigned int mhel=0;mhel<hadron.size();++mhel)
    {
      for(unsigned int lhel=0;lhel<lepton.size();++lhel)
	{
	  // map the index for the leptons to a helicity state
	  for(unsigned int ix=decay.size();ix>0;--ix)
	    {if(ix-1!=imes)
		{ihel[ix]=(lhel%constants[ix-1])/constants[ix];}}
	  // helicities of mesons
	  ihel[0]=0;
	  ihel[imes+1]=mhel;
	  newME(ihel)= lepton[lhel]*hadron[mhel];
	}
    }
  RhoDMatrix temp(PDT::Spin0); temp.average();
  // store the matrix element
  ME(newME);
  double ckm(1.);
  if(iq<=6)
    {
      if(iq%2==0){ckm = SM().CKM(abs(iq)/2-1,(abs(ia)-1)/2);}
      else{ckm = SM().CKM(abs(ia)/2-1,(abs(iq)-1)/2);}
    }
  // return the answer
  return 0.5*(newME.contract(temp)).real()*_GF*_GF*ckm; 
}
 
// output the setup information for the particle database
  void SemiLeptonicScalarDecayer::dataBaseOutput(ofstream & output,
						 bool header) const
{
  if(header){output << "update decayers set parameters=\"";}
  DecayIntegrator::dataBaseOutput(output,false);
  output << "set " << fullName() << ":GFermi "   << _GF*GeV2 << "\n";
  for(unsigned int ix=0;ix<_maxwgt.size();++ix)
    {output << "insert " << fullName() << ":MaximumWeight " << ix << " " 
	    << _maxwgt[ix] << "\n";}
  _current->dataBaseOutput(output,false,true);
  output << "set " << fullName() << ":Current " << _current->fullName() << " \n";
  _form->dataBaseOutput(output,false,true);
  output << "set " << fullName() << ":FormFactor " << _form->fullName() << " \n";
  if(header){output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;}
}

}
