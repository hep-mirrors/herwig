// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LeptonNeutrinoCurrent class.
//

#include "LeptonNeutrinoCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "Herwig++/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/FermionSpinInfo.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "LeptonNeutrinoCurrent.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;
using Helicity::SpinorWaveFunction;
using Helicity::SpinorBarWaveFunction;
using ThePEG::Helicity::LorentzPolarizationVector;
using ThePEG::Helicity::FermionSpinPtr;
using ThePEG::Helicity::FermionSpinInfo;
using ThePEG::Helicity::DiracRep;
using ThePEG::Helicity::HaberDRep;
using ThePEG::Helicity::HELASDRep;
using Helicity::Direction;
using Helicity::incoming;
using Helicity::outgoing;

LeptonNeutrinoCurrent::~LeptonNeutrinoCurrent() {}

void LeptonNeutrinoCurrent::persistentOutput(PersistentOStream & os) const {
  // no data
}

void LeptonNeutrinoCurrent::persistentInput(PersistentIStream & is, int) {
  // no data
}

ClassDescription<LeptonNeutrinoCurrent> LeptonNeutrinoCurrent::initLeptonNeutrinoCurrent;
// Definition of the static class description member.

void LeptonNeutrinoCurrent::Init() {

  static ClassDocumentation<LeptonNeutrinoCurrent> documentation
    ("The \\classname{LeptonNeutrinoCurrent} class is designed to handle the "
     "leptonic decay of the weak current.");
}


// complete the construction of the decay mode for integration
bool LeptonNeutrinoCurrent::createMode(int icharge, unsigned int imode,
				       DecayPhaseSpaceModePtr mode,
				       unsigned int iloc,unsigned int ires,
				       DecayPhaseSpaceChannelPtr phase,Energy upp)
{
  // make sure the the decays are kinematically allowed
  bool kineallowed=true;
  Energy min = getParticleData(11+2*imode)->mass()+getParticleData(12+2*imode)->mass();
  if(min>=upp){kineallowed=false;}
  if(kineallowed==false){return kineallowed;}
  DecayPhaseSpaceChannelPtr newchannel;
  // set the resonances
  tPDPtr res;
  if(icharge==3)
    {res=getParticleData(ParticleID::Wplus);}
  else if(icharge==-3)
    {res=getParticleData(ParticleID::Wminus);}
  else
    {return false;}
  // create the channel
  newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
  newchannel->addIntermediate(res,0,0.0,iloc,iloc+1);
  newchannel->init();
  mode->addChannel(newchannel);
  // return if successful
  return kineallowed;
}

// the particles produced by the current
PDVector LeptonNeutrinoCurrent::particles(int icharge, unsigned int imode,
					  int iq,int ia)
{
  PDVector output(2);
  if(icharge==3)
    {
      output[0]=getParticleData(-11-2*imode);
      output[1]=getParticleData(12+2*imode);
    }
  else if(icharge==-3)
    {
      output[0]=getParticleData(11+2*imode);
      output[1]=getParticleData(-12-2*imode);
    }
  return output;
}

// hadronic current   
vector<LorentzPolarizationVector> 
LeptonNeutrinoCurrent::current(bool vertex, const int imode, const int ichan,
			       Energy & scale,const ParticleVector & outpart) const
{
  Lorentz5Momentum q=outpart[0]->momentum()+outpart[1]->momentum();q.rescaleMass();
  scale=q.mass();
  // storage for the currents
  vector<LorentzPolarizationVector> temp;
  // construct the spin information objects for the  decay products
  FermionSpinPtr lepspin,nuspin;
  if(vertex)
    {
      SpinPtr slep=new_ptr(FermionSpinInfo(outpart[0]->momentum(),true));
      outpart[0]->spinInfo(slep);
      lepspin=dynamic_ptr_cast<FermionSpinPtr>(slep);
      SpinPtr snu =new_ptr(FermionSpinInfo(outpart[1 ]->momentum(),true));
      outpart[1 ]->spinInfo(snu );
      nuspin =dynamic_ptr_cast<FermionSpinPtr>(snu);
    }
  // lepton wavefunctions for the different helicities
  vector<LorentzSpinor> wave;
  vector<LorentzSpinorBar> wavebar;
  if(outpart[0]->id()>0)
    {
      SpinorWaveFunction nu =SpinorWaveFunction(outpart[1]->momentum(),
						outpart[1]->dataPtr(),outgoing);
      SpinorBarWaveFunction lep=SpinorBarWaveFunction(outpart[0]->momentum(),
						      outpart[0]->dataPtr(),outgoing);
      for(int ix=-1;ix<2;ix+=2)
	{
	  nu.reset(ix);wave.push_back(nu.Wave());
	  if(vertex){nuspin->setBasisState(ix,wave[(ix+1)/2]);}
	  lep.reset(ix);wavebar.push_back(lep.Wave());
	  if(vertex){lepspin->setBasisState(ix,wavebar[(ix+1)/2].bar());}
	}
    }
  else
    {
      SpinorWaveFunction lep=SpinorWaveFunction(outpart[0]->momentum(),
						outpart[0]->dataPtr(),outgoing);
      SpinorBarWaveFunction nu=SpinorBarWaveFunction(outpart[1]->momentum(),
						     outpart[1]->dataPtr(),outgoing);
      for(int ix=-1;ix<2;ix+=2)
	{
	  lep.reset(ix);wave.push_back(lep.Wave());
	  if(vertex){lepspin->setBasisState(ix,wave[(ix+1)/2]);}
	  nu.reset(ix);wavebar.push_back(nu.Wave());
	  if(vertex){nuspin->setBasisState(ix,wavebar[(ix+1)/2].bar());}
	}
    }
  // now compute the currents
  LorentzPolarizationVector vec;
  Complex ii(0.,1.);
  temp.resize(4); int iloc=0;
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
	  // location in the vector
	  if(outpart[0]->id()>0){iloc=2*iy+ix;}
	  else{iloc=2*ix+iy;}
	  // add it to the vector
	  temp[iloc]=vec;
	}
    }
  // return the answer
  return temp;
}

bool LeptonNeutrinoCurrent::accept(vector<int> id)
{
  bool allowed(false);
  if(id.size()!=2){return false;}
  if(abs(id[0])%2==0)
    {if((id[0]> 10&&id[0]< 18&&id[1]==-id[0]+1)||
	(id[0]<-10&&id[0]>-18&&id[1]==-id[0]-1)){allowed=true;}}
  else
    {if((id[1]> 10&&id[1]< 18&&id[0]==-id[1]+1)||
	(id[1]<-10&&id[1]>-18&&id[0]==-id[1]-1)){allowed=true;}}
  return allowed;
}

// the decay mode
unsigned int LeptonNeutrinoCurrent::decayMode(vector<int> idout)
{
  unsigned int imode=((abs(idout[0])+abs(idout[0])%2)-12)/2;
  return imode;
}

// output the information for the database
void LeptonNeutrinoCurrent::dataBaseOutput(ofstream & output)
{output << "create /Herwig++/LeptonNeutrinoCurrent " << fullName() << " \n";}

}



