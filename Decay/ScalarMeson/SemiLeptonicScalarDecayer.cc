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
#include "ThePEG/Helicity/ScalarSpinInfo.h"
#include "ThePEG/Helicity/LorentzPolarizationVector.h"
#include "ThePEG/Helicity/VectorSpinInfo.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig++/Helicity/EpsFunction.h"
#include "ThePEG/Helicity/TensorSpinInfo.h"
#include "ThePEG/Helicity/LorentzTensor.h"
#include "Herwig++/Helicity/WaveFunction/TensorWaveFunction.h"

namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;
using ThePEG::Helicity::RhoDMatrix;
using ThePEG::Helicity::LorentzPolarizationVector;
using ThePEG::Helicity::tcScalarSpinPtr;
using ThePEG::Helicity::ScalarSpinInfo;
using ThePEG::Helicity::tcVectorSpinPtr;
using ThePEG::Helicity::VectorSpinInfo;
using Helicity::VectorWaveFunction;
using ThePEG::Helicity::tcTensorSpinPtr;
using ThePEG::Helicity::TensorSpinInfo;
using Helicity::TensorWaveFunction;
using Helicity::EpsFunction;
using Helicity::Direction;
using Helicity::incoming;
using Helicity::outgoing;


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
    ("The \\classname{SemiLeptonicScalarDecayer} class is designed for the"
    "semi-leptonic decay of a (pseudo)-scalar meson.");

  static Parameter<SemiLeptonicScalarDecayer,InvEnergy2> interfaceGFermi
    ("GFermi",
     "The Fermi coupling constant",
     &SemiLeptonicScalarDecayer::_GF, 1./GeV2, 1.16639E-5/GeV2,
     -1.0e12*1./GeV2, 1.0e12*1./GeV2,
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
     0, 0, 0, 0, 10000, false, false, true);

}

// combine the currents and form-factors to give the matrix element
double SemiLeptonicScalarDecayer::me2(bool vertex, const int ichan,
				      const Particle & inpart,
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
				  << " in SemiLeptonicScalarDecayer::me2()" 
				  << Exception::abortnow;}
  else
    {
      SpinPtr newspin=new_ptr(ScalarSpinInfo(inpart.momentum(),true));
      inspin = dynamic_ptr_cast<tcScalarSpinPtr>(newspin);
      inspin->decayed(true);
      const_ptr_cast<tPPtr>(&inpart)->spinInfo(newspin);
    }
  // get the information on the form-factor
  int jspin(0);
  int id0=inpart.id(),id1=decay[0]->id();
  bool dummy;
  int iloc(_form->formFactorNumber(id0,id1,dummy));
  int spect,iq,ia;
  _form->formFactorInfo(iloc,jspin,spect,iq,ia);
  SpinPtr smeson; VectorSpinPtr vspin; TensorSpinPtr tspin;
  // construct the spin info for the outgoing meson
  if(vertex)
    {
      if(jspin==0)
	{
	  smeson = new_ptr(ScalarSpinInfo(decay[0]->momentum(),true));
	  decay[0]->spinInfo(smeson);
	}
      else if(jspin==1)
	{
	  smeson=new_ptr(VectorSpinInfo(decay[0]->momentum(),true));
	  decay[0]->spinInfo(smeson);
	  vspin=dynamic_ptr_cast<VectorSpinPtr>(smeson);
	}
      else if(jspin==2)
	{
	  smeson=new_ptr(TensorSpinInfo(decay[0]->momentum(),true));
	  decay[0]->spinInfo(smeson);
	  tspin=dynamic_ptr_cast<TensorSpinPtr>(smeson);
	}
    }
  // work out the value of q and calculate the form factors
  Lorentz5Momentum q=inpart.momentum()-decay[0]->momentum();q.rescaleMass();
  Energy2 q2=q.mass2();
  Lorentz5Momentum sum=inpart.momentum()+decay[0]->momentum();
  // calculate the hadronic current for the decay
  vector<LorentzPolarizationVector> hadron;
  if(jspin==0)
    {
      Complex fp,f0;
      _form->ScalarScalarFormFactor(q2,iloc,id0,id1,inpart.mass(),decay[0]->mass(),
				    f0,fp);
      double pre((inpart.mass()*inpart.mass()-decay[0]->mass()*decay[0]->mass())/q2);
      hadron.push_back(fp*sum+pre*(f0-fp)*q);
    }
  else if(jspin==1)
    {
      Complex A0,A1,A2,A3,V;
      Energy MP=inpart.mass();
      Energy MV=decay[0]->mass();
      Energy msum  = MP+MV;
      Energy mdiff = MP-MV;
      _form->ScalarVectorFormFactor(q2,iloc,id0,id1,MP,MV,A0,A1,A2,V);
      A3 = 0.5/MV*(msum*A1-mdiff*A2);
      // wave function for the vector
      VectorWaveFunction vwave=VectorWaveFunction(decay[0]->momentum(),
						  decay[0]->dataPtr(),outgoing);
      LorentzPolarizationVector ptemp;
      // compute the hadron currents
      Complex dot,ii(0.,1.);
      for(int ix=-1;ix<2;++ix)
	{
	  // compute the polarization vector for this helicity
	  vwave.reset(ix);ptemp=vwave.Wave();
	  if(vertex){vspin->setBasisState(ix,ptemp);}
	      // dot product
	  dot = ptemp*inpart.momentum();
	  hadron.push_back(-ii*msum*A1*ptemp
			   +ii*A2/msum*dot*sum
			   +2.*ii*MV/q2*(A3-ii*A0)*q*dot
			   +2.*V/msum*EpsFunction::product(ptemp,inpart.momentum(),
							   decay[0]->momentum()));
	}
    }
  else if(jspin==2)
    {
      Complex h,k,bp,bm;
      _form->ScalarTensorFormFactor(q2,iloc,id0,id1,inpart.mass(),decay[0]->mass(),
				    h,k,bp,bm);
      // wave function for the tensor
      TensorWaveFunction twave=TensorWaveFunction(decay[0]->momentum(),
						  decay[0]->dataPtr(),outgoing);
      LorentzPolarizationVector dotv;
      LorentzTensor ptemp;
      // compute the hadron currents
      Complex dot,ii(0.,1.);
      for(int ix=-2;ix<3;++ix)
	{
	  // compute the tensor for this helicity
	  twave.reset(ix);ptemp=twave.Wave();
	  if(vertex){tspin->setBasisState(ix,ptemp);}
	  dotv = ptemp*inpart.momentum();
	  dot = dotv*inpart.momentum();
	  hadron.push_back(ii*h*EpsFunction::product(dotv,sum,q)
			   -k*dotv-bp*dot*sum-bm*dot*q);
	}
    }
  int mode=(abs(decay[1]->id())-11)/12;
  // construct the lepton current
  vector<LorentzPolarizationVector> lepton(_current->current(vertex,mode,ichan,
							      inpart,decay));
  // work out the mapping for the lepton vector
  vector<int> constants(decay.size()+1), ispin(decay.size()),ihel(decay.size()+1);
  int itemp=1; unsigned int imes=0;
  for(int ix=int(decay.size()-1);ix>=0;--ix)
    {
      ispin[ix]=decay[ix]->data().iSpin();
      if(abs(decay[ix]->id())<=16){itemp*=ispin[ix];constants[ix]=itemp;}
      else{imes=ix;}
    }
  constants[decay.size()]=1;
  constants[imes]=constants[imes+1];
  DecayMatrixElement newME(1,ispin);
  for(unsigned int mhel=0;mhel<hadron.size();++mhel)
    {
      for(unsigned int lhel=0;lhel<lepton.size();++lhel)
	{
	  // map the index for the leptons to a helicity state
	  for(unsigned int ix=decay.size();ix>0;--ix)
	    {
	      if(ix-1!=imes)
		{
		  ihel[ix]=(lhel%constants[ix-1])/constants[ix]-int(ispin[ix-1]/2);
		  if(ispin[ix-1]%2==0&&ihel[ix]>=0&&ispin[ix-1]!=0){++ihel[ix];}
		}
	    }
	  // helicities of mesons
	  ihel[0]=0;
	  ihel[imes+1]=-jspin+mhel;
	  newME(ihel)= lepton[lhel]*hadron[mhel];
	}
    }
  RhoDMatrix temp(1); temp.average();
  // store the matrix element
  ME(newME);
  double ckm(1.);
  if(iq<=6)
    {
      if(iq%2==0){ckm = SM().CKM(iq/2-1,(abs(ia)-1)/2);}
      else{ckm = SM().CKM(abs(ia)/2-1,(iq-1)/2);}
      //      cout << "testing the CKM factor " << iq << " " << ia << " " << ckm << endl;
    }
  // return the answer
  double me= 0.5*(newME.contract(temp)).real()*_GF*_GF;
  return me;  
}
 
// output the setup information for the particle database
void SemiLeptonicScalarDecayer::dataBaseOutput(ofstream & output)
{
  output << "update decayers set parameters=\"";
  output << "set " << fullName() << ":Iteration " << _niter << "\n";
  output << "set " << fullName() << ":Ntry " << _ntry << "\n";
  output << "set " << fullName() << ":GFermi "   << _GF*GeV2 << "\n";
  for(unsigned int ix=0;ix<_maxwgt.size();++ix)
    {output << "insert " << fullName() << ":MaximumWeight " << ix << " " 
	    << _maxwgt[ix] << "\n";}
  output << "set " << fullName() << ":Current " << _current->fullName() << " \n";
  output << "set " << fullName() << ":FormFactor " << _form->fullName() << " \n";
  output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}

}
