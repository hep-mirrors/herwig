// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFSVertex class.
//

#include "FFSVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
namespace Helicity {

using namespace ThePEG;
    
FFSVertex::~FFSVertex() {}
    
void FFSVertex::persistentOutput(PersistentOStream & os) const {}
    
void FFSVertex::persistentInput(PersistentIStream & is, int) {}
  
AbstractClassDescription<FFSVertex> FFSVertex::initFFSVertex;
// Definition of the static class description member.
    
void FFSVertex::Init() {

  static ClassDocumentation<FFSVertex> documentation
    ("The FFSVertex class is the implementation of the FFS"
     "vertex. All such vertices shoud inherit from it");
  
}
  
// set coupling member
void FFSVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b,tcPDPtr c){;}

// evaluate the full vertex
Complex FFSVertex::evaluate(Energy2 q2, const SpinorWaveFunction & sp,
				    const SpinorBarWaveFunction & sbar,
				    const ScalarWaveFunction & sca)
{
  // extract the pointers to the particle data objects
  tcPDPtr Psp=sp.getParticle();
  tcPDPtr Psca=sca.getParticle();
  tcPDPtr Psbar=sbar.getParticle();
  // calculate the couplings
  setCoupling(q2,Psp,Psca,Psbar);
  Complex norm=getNorm();
  Complex ii(0.,1.);
  Complex vertex(0.);
  // low energy conventions
  if(sp.Wave().Rep()==HaberDRep&&sbar.Wave().Rep()==HaberDRep)
    {
      vertex = _left*( (sbar.s1()-sbar.s3())*(sp.s1()-sp.s3())
		       +(sbar.s2()-sbar.s4())*(sp.s2()-sp.s4()))
	+_right*( (sbar.s1()+sbar.s3())*(sp.s1()+sp.s3())
		  +(sbar.s2()+sbar.s4())*(sp.s2()+sp.s4()));
      vertex=0.5*vertex;
    }
  // high energy conventions
  else if(sp.Wave().Rep()==HELASDRep&&sbar.Wave().Rep()==HELASDRep)
    {
      vertex=  _left*(sbar.s1()*sp.s1()+sbar.s2()*sp.s2())
	+_right*(sbar.s3()*sp.s3()+sbar.s4()*sp.s4());
    }
  // mixing conventions
  else
    {
      sp.Wave().changeRep(HELASDRep);
      sbar.Wave().changeRep(HELASDRep);
      vertex=  _left*(sbar.s1()*sp.s1()+sbar.s2()*sp.s2())
	+_right*(sbar.s3()*sp.s3()+sbar.s4()*sp.s4());
    }
  // final factors
  vertex = ii*norm*vertex*sca.Wave();
  return vertex;
}

// off-shell scalar
ScalarWaveFunction FFSVertex::evaluate(Energy2 q2,int iopt, tcPDPtr out, 
				       const SpinorWaveFunction & sp,
				       const SpinorBarWaveFunction & sbar)
{
  // work out the momentum of the off-shell particle
  Lorentz5Momentum pout = Lorentz5Momentum(sbar.px()+sp.px(),sbar.py()+sp.py(),
					   sbar.pz()+sp.pz(),sbar.e() +sp.e());
  // extract the pointers to the particle data objects
  tcPDPtr Psbar=sbar.getParticle();
  tcPDPtr Psp=sp.getParticle();
  // first calculate the couplings
  setCoupling(q2,Psp,out,Psbar);
  Energy2 p2 = pout.m2();
  Complex fact=getNorm()*propagator(iopt,p2,out);
  Complex ii(0.,1.); 
  Complex output;
  // low energy conventions
  if(sp.Wave().Rep()==HaberDRep&&sbar.Wave().Rep()==HaberDRep)
    {
      output = _left*( (sbar.s1()-sbar.s3())*(sp.s1()-sp.s3())
		       +(sbar.s2()-sbar.s4())*(sp.s2()-sp.s4()))
	+_right*( (sbar.s1()+sbar.s3())*(sp.s1()+sp.s3())
		  +(sbar.s2()+sbar.s4())*(sp.s2()+sp.s4()));
      output=0.5*output;
    }
  // high energy conventions
  else if(sp.Wave().Rep()==HELASDRep&&sbar.Wave().Rep()==HELASDRep)
    {
      output =  _left*(sbar.s1()*sp.s1()+sbar.s2()*sp.s2())
	+_right*(sbar.s3()*sp.s3()+sbar.s4()*sp.s4());
    }
  else
    {
      sp.Wave().changeRep(HELASDRep);
      sbar.Wave().changeRep(HELASDRep);
      output =  _left*(sbar.s1()*sp.s1()+sbar.s2()*sp.s2())
	+_right*(sbar.s3()*sp.s3()+sbar.s4()*sp.s4());
    }
  // final factors and output
  output=output*fact;
  return ScalarWaveFunction(pout,out,output);
}
    
// off-shell spinor
SpinorWaveFunction FFSVertex::evaluate(Energy2 q2, int iopt,tcPDPtr out,
				       const SpinorWaveFunction & sp,
				       const ScalarWaveFunction & sca, DiracRep dirac)
{
  // work out the momentum of the off-shell particle
  Lorentz5Momentum pout = Lorentz5Momentum(sp.px()+sca.px(),sp.py()+sca.py(),
					   sp.pz()+sca.pz(),sp.e() +sca.e() );
  // extract the pointers to the particle data objects
  tcPDPtr Psca=sca.getParticle();
  tcPDPtr Psp=sp.getParticle();
  // first calculate the couplings
  setCoupling(q2,Psp,Psca,out);
  double p2 = pout.m2();
  Complex fact=-getNorm()*sca.Wave()*propagator(iopt,p2,out);
  Complex ii(0.,1.);
  // useful combinations of the momenta
  Energy  mass  = out->mass();
  Complex p1p2=pout.px()+ii*pout.py();
  Complex p1m2=pout.px()-ii*pout.py();
  Complex s1(0.),s2(0.),s3(0.),s4(0.);
  // ensure the spinor is in the correct dirac representation
  sp.Wave().changeRep(dirac);
  // low energy convention
  if(dirac==HaberDRep)
    {
      fact = 0.5*fact;
      Complex p0pm=pout.e()+mass;
      Complex p0mm=pout.e()-mass;
      Complex pz=pout.pz();
      Complex lpr=_left+_right;
      Complex lmr=_left-_right;
      s1 = fact*( lpr*(sp.s1()*p0pm-  pz*sp.s3()-p1m2*sp.s4())
		  +lmr*(  pz*sp.s1()+p1m2*sp.s2()-p0pm*sp.s3()));
      s2 = fact*( lpr*(sp.s2()*p0pm-p1p2*sp.s3()+  pz*sp.s4())
		  +lmr*(p1p2*sp.s1()-  pz*sp.s2()-p0pm*sp.s4()));
      s3 =-fact*( lpr*(sp.s3()*p0mm-  pz*sp.s1()-p1m2*sp.s2())
		  +lmr*(  pz*sp.s3()+p1m2*sp.s4()-p0mm*sp.s1()));
      s4 =-fact*( lpr*(sp.s4()*p0mm-p1p2*sp.s1()+  pz*sp.s2())
		  +lmr*(p1p2*sp.s3()-  pz*sp.s4()-p0mm*sp.s2()));
    }
  // high energy convention
  else if(dirac==HELASDRep)
    {
      Complex p0p3=pout.e()+pout.pz();
      Complex p0m3=pout.e()-pout.pz();
      s1 = fact*( _left*mass*sp.s1()+_right*(p0m3*sp.s3()-p1m2*sp.s4()));
      s2 = fact*( _left*mass*sp.s2()+_right*(p0p3*sp.s4()-p1p2*sp.s3()));
      s3 = fact*(_right*mass*sp.s3()+ _left*(p0p3*sp.s1()+p1m2*sp.s2()));
      s4 = fact*(_right*mass*sp.s4()+ _left*(p0m3*sp.s2()+p1p2*sp.s1()));
    }
  return SpinorWaveFunction(pout,out,s1,s2,s3,s4,dirac);
}

// off-shell SpinorBar
SpinorBarWaveFunction FFSVertex::evaluate(Energy q2,int iopt,tcPDPtr out,
					  const SpinorBarWaveFunction & sbar,
					  const ScalarWaveFunction & sca,
					  DiracRep dirac)
{
  // work out the momentum of the off-shell particle
  Lorentz5Momentum pout = Lorentz5Momentum(sbar.px()+sca.px(),sbar.py()+sca.py(),
					   sbar.pz()+sca.pz(),sbar.e() +sca.e());
  // extract the pointers to the particle data objects
  tcPDPtr  Psbar=sbar.getParticle();
  tcPDPtr  Psca =sca.getParticle();
  // first calculate the couplings
  setCoupling(q2,out,Psca,Psbar);
  Energy2 p2=pout.m2();
  Complex fact=-getNorm()*sca.Wave()*propagator(iopt,p2,out);
  Complex ii(0.,1.);
  // momentum components
  Energy mass  = out->mass();
  Complex p1p2=pout.px()+ii*pout.py();
  Complex p1m2=pout.px()-ii*pout.py();
  // complex numbers for the spinor
  Complex s1(0.),s2(0.),s3(0.),s4(0.);
  // ensure the spinorbar is in the correct dirac representation
  sbar.Wave().changeRep(dirac);
  // low energy convention
  if(dirac==HaberDRep)
    {
      fact = 0.5*fact;
      Complex p0pm=pout.e()+mass;
      Complex p0mm=pout.e()-mass;
      Complex pz=pout.pz();
      Complex lpr=_left+_right;
      Complex lmr=_left-_right;
      s1 = fact*(-lpr*(p0mm*sbar.s1()+  pz*sbar.s3()+p1p2*sbar.s4())
		 +lmr*(  pz*sbar.s1()+p1p2*sbar.s2()+p0mm*sbar.s3()));
      s2 = fact*(-lpr*(p0mm*sbar.s2()+p1m2*sbar.s3()-  pz*sbar.s4())
		 +lmr*(p1m2*sbar.s1()-  pz*sbar.s2()+p0mm*sbar.s4()));
      s3 = fact*(+lpr*(p0pm*sbar.s3()+  pz*sbar.s1()+p1p2*sbar.s2())
		 -lmr*(  pz*sbar.s3()+p1p2*sbar.s4()-p0pm*sbar.s1()));
      s4 = fact*(+lpr*(p0pm*sbar.s4()+p1m2*sbar.s1()-  pz*sbar.s2())
		 -lmr*(p1m2*sbar.s3()-  pz*sbar.s4()-p0pm*sbar.s2()));
    }
  // high energy convention
  else if(dirac==HELASDRep)
    {
      Complex p0p3=pout.e() +   pout.pz();
      Complex p0m3=pout.e() -   pout.pz();
      s1 = fact*( mass*_left*sbar.s1()-_right*(p0p3*sbar.s3()+p1p2*sbar.s4()));
      s2 = fact*( mass*_left*sbar.s2()-_right*(p1m2*sbar.s3()+p0m3*sbar.s4()));
      s3 = fact*(mass*_right*sbar.s3()- _left*(p0m3*sbar.s1()-p1p2*sbar.s2()));
      s4 = fact*(mass*_right*sbar.s4()+ _left*(p1m2*sbar.s1()-p0p3*sbar.s2()));
    }
  return SpinorBarWaveFunction(pout,out,s1,s2,s3,s4,dirac);
}    

}
}
