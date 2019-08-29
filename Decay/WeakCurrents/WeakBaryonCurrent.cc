// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the WeakBaryonCurrent class.
//

#include "WeakBaryonCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/HelicityFunctions.h"

using namespace Herwig;

WeakBaryonCurrent::WeakBaryonCurrent() {}

IBPtr WeakBaryonCurrent::clone() const {
  return new_ptr(*this);
}

IBPtr WeakBaryonCurrent::fullclone() const {
  return new_ptr(*this);
}

void WeakBaryonCurrent::persistentOutput(PersistentOStream & os) const {
  os << formFactor_;
}

void WeakBaryonCurrent::persistentInput(PersistentIStream & is, int) {
  is >> formFactor_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<WeakBaryonCurrent,WeakCurrent>
describeHerwigWeakBaryonCurrent("Herwig::WeakBaryonCurrent", "Herwig.so");

void WeakBaryonCurrent::Init() {

  static ClassDocumentation<WeakBaryonCurrent> documentation
    ("The WeakBaryonCurrent class is a wrapper for the BaryonFormFactor"
     " so it can be used as a WeakCurrent");

  static Reference<WeakBaryonCurrent,BaryonFormFactor> interfaceFormFactor
    ("FormFactor",
     "The baryon form factor",
     &WeakBaryonCurrent::formFactor_, false, false, true, false, false);

}

void WeakBaryonCurrent::doinit() {
  // initialize the form factor
  formFactor_->init();
  // set the modes
  for(unsigned int iloc=0;iloc<formFactor_->numberOfFactors();++iloc) {
    int ispin(0), ospin(0), spect1(0), spect2(0), inquark(0), outquark(0);
    formFactor_->formFactorInfo(iloc,ispin,ospin,spect1,spect2,inquark,outquark);
    addDecayMode(outquark,-inquark);
  }
  setInitialModes(formFactor_->numberOfFactors());
  WeakCurrent::doinit();
}

// complete the construction of the decay mode for integration
bool WeakBaryonCurrent::createMode(int icharge, tcPDPtr ,
				   IsoSpin::IsoSpin Itotal, IsoSpin::I3 i3, Strangeness::Strange S,
				   unsigned int imode,PhaseSpaceModePtr mode,
				   unsigned int iloc,int ires,
				   PhaseSpaceChannel phase, Energy upp ) {
  // todo isospin in the form factors
  // no isospin here
  if(Itotal!=IsoSpin::IUnknown || i3 !=IsoSpin::I3Unknown) return false;
  unsigned int iq(0),ia(0);
  tPDVector out = particles(icharge,imode,iq,ia);
  // make sure the the decays are kinematically allowed
  Energy min =out[0]->massMin()+out[1]->massMin();
  if(min>=upp) return false;
  // set the resonances and check charge
  tPDPtr res;
  if(icharge==3)       res=getParticleData(ParticleID::Wplus );
  else if(icharge==-3) res=getParticleData(ParticleID::Wminus);
  else                 res=getParticleData(ParticleID::gamma );
  // create the channel
  mode->addChannel((PhaseSpaceChannel(phase),ires,res,ires+1,iloc+1,ires+1,iloc+2));
  // return if successful
  return true;
}

// the particles produced by the current
tPDVector WeakBaryonCurrent::particles(int icharge, unsigned int imode, int , int ) {
  tPDVector extpart(2);
  int id0(0),id1(0);
  formFactor_->particleID(imode,id0,id1);
  extpart[0] = getParticleData(id0);
  if(extpart[0]->CC()) extpart[0]=extpart[0]->CC();
  extpart[1] = getParticleData(id1);
  int charge = extpart[0]->iCharge()+extpart[1]->iCharge();
  if(charge==icharge)
    return extpart;
  else if(charge==-icharge) {
    for(unsigned int ix=0;ix<2;++ix)
      if(extpart[ix]->CC()) extpart[ix]=extpart[ix]->CC();
    return extpart;
  }
  else
    return tPDVector();
}

void WeakBaryonCurrent::constructSpinInfo(ParticleVector decay) const {
  if(decay[0]->id()>0) {
    SpinorWaveFunction   ::constructSpinInfo(wave_   ,decay[1],outgoing,true);
    SpinorBarWaveFunction::constructSpinInfo(wavebar_,decay[0],outgoing,true);
  }
  else {
    SpinorWaveFunction   ::constructSpinInfo(   wave_,decay[0],outgoing,true);
    SpinorBarWaveFunction::constructSpinInfo(wavebar_,decay[1],outgoing,true);
  }
}

// hadronic current   
vector<LorentzPolarizationVectorE> 
WeakBaryonCurrent::current(tcPDPtr ,
			       IsoSpin::IsoSpin Itotal, IsoSpin::I3 i3, Strangeness::Strange S,
			       const int, const int, Energy & scale, 
			       const tPDVector & outgoing,
			       const vector<Lorentz5Momentum> & momenta,
			       DecayIntegrator::MEOption) const {
  // no isospin here
  if(Itotal!=IsoSpin::IUnknown || i3 !=IsoSpin::I3Unknown) return vector<LorentzPolarizationVectorE>();
  useMe();
  Lorentz5Momentum q = momenta[0]+momenta[1];
  q.rescaleMass();
  scale=q.mass();
  int in  = abs(outgoing[0]->id());
  int out = abs(outgoing[1]->id());
  Energy m1 = outgoing[0]->mass();
  Energy m2 = outgoing[1]->mass();
  bool cc = false;
  unsigned int imode = formFactor_->formFactorNumber(in,out,cc);
  // todo generalize to spin != 1/2
  assert(outgoing[0]->iSpin()==PDT::Spin1Half &&
	 outgoing[1]->iSpin()==PDT::Spin1Half );
  wave_.resize(2);
  wavebar_.resize(2);
  for(unsigned int ix=0;ix<2;++ix) {
    wavebar_[ix] = HelicityFunctions::dimensionedSpinorBar(-momenta[0],ix,Helicity::outgoing);
    wave_[ix] = HelicityFunctions::dimensionedSpinor   (-momenta[1],ix,Helicity::outgoing);
  }
  // get the form factors
  Complex f1v(0.),f2v(0.),f3v(0.),f1a(0.),f2a(0.),f3a(0.);
  formFactor_->SpinHalfSpinHalfFormFactor(sqr(scale),imode,in,out,m1,m2,
					  f1v,f2v,f3v,f1a,f2a,f3a,
					  BaryonFormFactor::TimeLike);
  Complex left  = f1v - f1a + f2v -double((m1-m2)/(m1+m2))*f2a;
  Complex right = f1v + f1a + f2v +double((m1-m2)/(m1+m2))*f2a;
  vector<LorentzPolarizationVectorE> baryon;
  Lorentz5Momentum diff = momenta[0]-momenta[1];
  for(unsigned int ohel1=0;ohel1<2;++ohel1) {
    for(unsigned int ohel2=0;ohel2<2;++ohel2) {
      LorentzPolarizationVectorE 
	vtemp = wave_[ohel2].generalCurrent(wavebar_[ohel1],left,right);
      complex<Energy> vspin=wave_[ohel2].scalar      (wavebar_[ohel1]);      
      complex<Energy> aspin=wave_[ohel2].pseudoScalar(wavebar_[ohel1]);
      vtemp-= (f2v*vspin+f2a*aspin)/(m1+m2)*diff;
      vtemp+= (f3v*vspin+f3a*aspin)/(m1+m2)*q;
      baryon.push_back(vtemp);
    }
  }
  // return the answer
  return baryon;
}

bool WeakBaryonCurrent::accept(vector<int> id) {
  assert(id.size()==2);
  int itemp[2] = {id[0],id[1]};
  for(unsigned int ix=0;ix<2;++ix)
    if(itemp[ix]<0) itemp[ix]=-itemp[ix];
  bool cc = false;
  return formFactor_->formFactorNumber(itemp[0],itemp[1],cc)>=0;
}

// the decay mode
unsigned int WeakBaryonCurrent::decayMode(vector<int> idout) {
  assert(idout.size()==2);
  int itemp[2] = {idout[0],idout[1]};
  for(unsigned int ix=0;ix<2;++ix)
    if(itemp[ix]<0) itemp[ix]=-itemp[ix];
  bool cc = false;
  return formFactor_->formFactorNumber(itemp[0],itemp[1],cc);
}

// output the information for the database
void WeakBaryonCurrent::dataBaseOutput(ofstream & output,bool header,
					   bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::WeakBaryonCurrent " << name() << "\n";
  output << "newdef    " << name() << ":FormFactor " << formFactor_  << "\n";
  WeakCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}
