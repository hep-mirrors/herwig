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
  
AbstractNoPIOClassDescription<FFSVertex> FFSVertex::initFFSVertex;
// Definition of the static class description member.
    
void FFSVertex::Init() {

  static ClassDocumentation<FFSVertex> documentation
    ("The FFSVertex class is the implementation of the FFS"
     "vertex. All such vertices shoud inherit from it");
  
}
  
// evaluate the full vertex
Complex FFSVertex::evaluate(Energy2 q2, const SpinorWaveFunction & sp,
			    const SpinorBarWaveFunction & sbar,
			    const ScalarWaveFunction & sca) {
  // extract the pointers to the particle data objects
  tcPDPtr Psp=sp.getParticle();
  tcPDPtr Psca=sca.getParticle();
  tcPDPtr Psbar=sbar.getParticle();
  // calculate the couplings
  //integer to determine incoming
  int iint(0);
  if(sp.direction() == Helicity::incoming || 
     sp.direction() == Helicity::intermediate) iint = 1;
  else if(sca.direction() == Helicity::incoming ||
	  sca.direction() == Helicity::intermediate) iint = 2;
  else if(sbar.direction() == Helicity::incoming ||
	  sbar.direction() == Helicity::intermediate) iint = 3;
  else 
    throw HelicityLogicalError() << "There is no incoming particle in "
				 << fullName() << Exception::runerror;
  setCoupling(q2,Psp,Psca,Psbar,iint);
  Complex norm=getNorm();
  Complex ii(0.,1.);
  Complex vertex(0.);
  // low energy conventions
  if(sp.wave().Rep()==HaberDRep&&sbar.wave().Rep()==HaberDRep) {
    vertex = _left*( (sbar.s1()-sbar.s3())*(sp.s1()-sp.s3())
		    +(sbar.s2()-sbar.s4())*(sp.s2()-sp.s4()))
           +_right*( (sbar.s1()+sbar.s3())*(sp.s1()+sp.s3())
		    +(sbar.s2()+sbar.s4())*(sp.s2()+sp.s4()));
    vertex=0.5*vertex;
  }
  // high energy conventions
  else if(sp.wave().Rep()==HELASDRep&&sbar.wave().Rep()==HELASDRep) {
    vertex=  _left*(sbar.s1()*sp.s1()+sbar.s2()*sp.s2())
      +_right*(sbar.s3()*sp.s3()+sbar.s4()*sp.s4());
  }
  // mixing conventions
  else {
    LorentzSpinorBar<double> sbart=sbar.wave().transformRep(HELASDRep);
    LorentzSpinor<double>    spt  =sp  .wave().transformRep(HELASDRep);
    vertex=  _left*(sbart.s1()*spt.s1()+sbart.s2()*spt.s2())
           +_right*(sbart.s3()*spt.s3()+sbart.s4()*spt.s4());
  }
  // final factors
  vertex *= ii*norm*sca.wave();
  return vertex;
}

// off-shell scalar
ScalarWaveFunction FFSVertex::evaluate(Energy2 q2,int iopt, tcPDPtr out, 
				       const SpinorWaveFunction & sp,
				       const SpinorBarWaveFunction & sbar) {
  // work out the momentum of the off-shell particle
  Lorentz5Momentum pout = Lorentz5Momentum(sbar.px()+sp.px(),sbar.py()+sp.py(),
					   sbar.pz()+sp.pz(),sbar.e() +sp.e());
  // extract the pointers to the particle data objects
  tcPDPtr Psbar=sbar.getParticle();
  tcPDPtr Psp=sp.getParticle();
  int iint(0);
  if(sp.direction() == Helicity::incoming) iint = 1;
  else if(sbar.direction() == Helicity::incoming) iint = 3;
  else iint = 2;
  // first calculate the couplings
  setCoupling(q2,Psp,out,Psbar,iint);
  Energy2 p2 = pout.m2();
  Complex fact=getNorm()*propagator(iopt,p2,out);
  Complex output;
  // low energy conventions
  if(sp.wave().Rep()==HaberDRep&&sbar.wave().Rep()==HaberDRep) {
    output = _left*( (sbar.s1()-sbar.s3())*(sp.s1()-sp.s3())
          	    +(sbar.s2()-sbar.s4())*(sp.s2()-sp.s4()))
           +_right*( (sbar.s1()+sbar.s3())*(sp.s1()+sp.s3())
		    +(sbar.s2()+sbar.s4())*(sp.s2()+sp.s4()));
    output*=0.5;
  }
  // high energy conventions
  else if(sp.wave().Rep()==HELASDRep&&sbar.wave().Rep()==HELASDRep) {
    output =  _left*(sbar.s1()*sp.s1()+sbar.s2()*sp.s2())
	    +_right*(sbar.s3()*sp.s3()+sbar.s4()*sp.s4());
  }
  else {
    LorentzSpinor<double>    spt  =sp  .wave().transformRep(HELASDRep);
    LorentzSpinorBar<double> sbart=sbar.wave().transformRep(HELASDRep);
    output =  _left*(sbart.s1()*spt.s1()+sbart.s2()*spt.s2())
            +_right*(sbart.s3()*spt.s3()+sbart.s4()*spt.s4());
  }
  // final factors and output
  output*=fact;
  return ScalarWaveFunction(pout,out,output);
}
    
// off-shell spinor
SpinorWaveFunction FFSVertex::evaluate(Energy2 q2, int iopt,tcPDPtr out,
				       const SpinorWaveFunction & sp,
				       const ScalarWaveFunction & sca,
				       DiracRep dirac) {
  // work out the momentum of the off-shell particle
  Lorentz5Momentum pout = Lorentz5Momentum(sp.px()+sca.px(),sp.py()+sca.py(),
					   sp.pz()+sca.pz(),sp.e() +sca.e() );
  // extract the pointers to the particle data objects
  tcPDPtr Psca=sca.getParticle();
  tcPDPtr Psp=sp.getParticle();
  int iint(0);
  if(sp.direction() == Helicity::incoming) iint = 1;
  else if(sca.direction() == Helicity::incoming) iint = 2;
  else iint = 3;
  // first calculate the couplings
  setCoupling(q2,Psp,Psca,out, iint);
  Energy2 p2 = pout.m2();
  Complex fact=-getNorm()*sca.wave()*propagator(iopt,p2,out);
  Complex ii(0.,1.);
  // useful combinations of the momenta
  Energy  mass  = out->mass();
  complex<Energy> p1p2=pout.x()+ii*pout.y();
  complex<Energy> p1m2=pout.x()-ii*pout.y();
  Complex s1(0.),s2(0.),s3(0.),s4(0.);
  // ensure the spinor is in the correct dirac representation
  LorentzSpinor<double>    spt  =sp.wave().transformRep(dirac);
  // low energy convention
  if(dirac==HaberDRep)
    {
      fact = 0.5*fact;
      complex<Energy> p0pm=pout.e()+mass;
      complex<Energy> p0mm=pout.e()-mass;
      complex<Energy> pz=pout.z();
      Complex lpr=_left+_right;
      Complex lmr=_left-_right;
      s1 = UnitRemoval::InvE * 
	fact*( lpr*(spt.s1()*p0pm-  pz*spt.s3()-p1m2*spt.s4())
	       +lmr*(  pz*spt.s1()+p1m2*spt.s2()-p0pm*spt.s3()));
      s2 = UnitRemoval::InvE * 
	fact*( lpr*(spt.s2()*p0pm-p1p2*spt.s3()+  pz*spt.s4())
	       +lmr*(p1p2*spt.s1()-  pz*spt.s2()-p0pm*spt.s4()));
      s3 = UnitRemoval::InvE * 
	-fact*( lpr*(spt.s3()*p0mm-  pz*spt.s1()-p1m2*spt.s2())
		+lmr*(  pz*spt.s3()+p1m2*spt.s4()-p0mm*spt.s1()));
      s4 = UnitRemoval::InvE * 
	-fact*( lpr*(spt.s4()*p0mm-p1p2*spt.s1()+  pz*spt.s2())
		+lmr*(p1p2*spt.s3()-  pz*spt.s4()-p0mm*spt.s2()));
    }
  // high energy convention
  else if(dirac==HELASDRep)
    {
      complex<Energy> p0p3=pout.e()+pout.z();
      complex<Energy> p0m3=pout.e()-pout.z();
      s1 = UnitRemoval::InvE * 
	fact*( _left*mass*spt.s1()+_right*(p0m3*spt.s3()-p1m2*spt.s4()));
      s2 = UnitRemoval::InvE * 
	fact*( _left*mass*spt.s2()+_right*(p0p3*spt.s4()-p1p2*spt.s3()));
      s3 = UnitRemoval::InvE * 
	fact*(_right*mass*spt.s3()+ _left*(p0p3*spt.s1()+p1m2*spt.s2()));
      s4 = UnitRemoval::InvE * 
	fact*(_right*mass*spt.s4()+ _left*(p0m3*spt.s2()+p1p2*spt.s1()));
    }
  return SpinorWaveFunction(pout,out,s1,s2,s3,s4,dirac);
}

// off-shell SpinorBar
SpinorBarWaveFunction FFSVertex::evaluate(Energy2 q2,int iopt,tcPDPtr out,
					  const SpinorBarWaveFunction & sbar,
					  const ScalarWaveFunction & sca,
					  DiracRep dirac) {
  // work out the momentum of the off-shell particle
  Lorentz5Momentum pout = Lorentz5Momentum(sbar.px()+sca.px(),sbar.py()+sca.py(),
					   sbar.pz()+sca.pz(),sbar.e() +sca.e());
  // extract the pointers to the particle data objects
  tcPDPtr  Psbar=sbar.getParticle();
  tcPDPtr  Psca =sca.getParticle();
  int iint(0);
  if(sca.direction() == Helicity::incoming) iint = 2;
  else if(sbar.direction() == Helicity::incoming) iint = 3;
  else iint = 1;
  // first calculate the couplings
  setCoupling(q2,out,Psca,Psbar, iint);
  Energy2 p2=pout.m2();
  Complex fact=-getNorm()*sca.wave()*propagator(iopt,p2,out);
  Complex ii(0.,1.);
  // momentum components
  Energy mass  = out->mass();
  complex<Energy> p1p2=pout.x()+ii*pout.y();
  complex<Energy> p1m2=pout.x()-ii*pout.y();
  // complex numbers for the spinor
  Complex s1(0.),s2(0.),s3(0.),s4(0.);
  // ensure the spinorbar is in the correct dirac representation
  LorentzSpinorBar<double> sbart=sbar.wave().transformRep(dirac);
  // low energy convention
  if(dirac==HaberDRep)
    {
      fact = 0.5*fact;
      complex<Energy> p0pm=pout.e()+mass;
      complex<Energy> p0mm=pout.e()-mass;
      complex<Energy> pz=pout.z();
      Complex lpr=_left+_right;
      Complex lmr=_left-_right;
      s1 = UnitRemoval::InvE * 
	fact*(-lpr*(p0mm*sbart.s1()+  pz*sbart.s3()+p1p2*sbart.s4())
	      +lmr*(  pz*sbart.s1()+p1p2*sbart.s2()+p0mm*sbart.s3()));
      s2 = UnitRemoval::InvE * 
	fact*(-lpr*(p0mm*sbart.s2()+p1m2*sbart.s3()-  pz*sbart.s4())
	      +lmr*(p1m2*sbart.s1()-  pz*sbart.s2()+p0mm*sbart.s4()));
      s3 = UnitRemoval::InvE * 
	fact*(+lpr*(p0pm*sbart.s3()+  pz*sbart.s1()+p1p2*sbart.s2())
	      -lmr*(  pz*sbart.s3()+p1p2*sbart.s4()-p0pm*sbart.s1()));
      s4 = UnitRemoval::InvE * 
	fact*(+lpr*(p0pm*sbart.s4()+p1m2*sbart.s1()-  pz*sbart.s2())
	      -lmr*(p1m2*sbart.s3()-  pz*sbart.s4()-p0pm*sbart.s2()));
    }
  // high energy convention
  else if(dirac==HELASDRep)
    {
      complex<Energy> p0p3=pout.e() +   pout.z();
      complex<Energy> p0m3=pout.e() -   pout.z();
      s1 = UnitRemoval::InvE * 
	fact*( mass*_left*sbart.s1()-_right*(p0p3*sbart.s3()+p1p2*sbart.s4()));
      s2 = UnitRemoval::InvE * 
	fact*( mass*_left*sbart.s2()-_right*(p1m2*sbart.s3()+p0m3*sbart.s4()));
      s3 = UnitRemoval::InvE * 
	fact*(mass*_right*sbart.s3()- _left*(p0m3*sbart.s1()-p1p2*sbart.s2()));
      s4 = UnitRemoval::InvE * 
	fact*(mass*_right*sbart.s4()+ _left*(p1m2*sbart.s1()-p0p3*sbart.s2()));
    }
  return SpinorBarWaveFunction(pout,out,s1,s2,s3,s4,dirac);
}    

}
}
