// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFTVertex class.
//

#include "FFTVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;
    
FFTVertex::~FFTVertex() {}
    
void FFTVertex::persistentOutput(PersistentOStream & os) const { }
    
void FFTVertex::persistentInput(PersistentIStream & is, int) { }
    
// Definition of the static class description member
ClassDescription<FFTVertex> FFTVertex::initFFTVertex;
    
void FFTVertex::Init() {
  
  static ClassDocumentation<FFTVertex> documentation
    ("The \\classname{FFTVertex} class is the implementation of"
     "the fermion-antifermion tensor vertices for helicity "
     "amplitude calculations. All such vertices should inherit"
     "from it.");
  
}

// coupling
void FFTVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b,tcPDPtr c){;}

// function to evaluate the vertex
Complex FFTVertex::evaluate(Energy2 q2,const SpinorWaveFunction & sp,
				    const SpinorBarWaveFunction & sbar,
				    const TensorWaveFunction & ten)
{
  // pointers to the particles
  tcPDPtr Psp = sp.getParticle();
  tcPDPtr Psbar = sbar.getParticle();
  tcPDPtr Pten = ten.getParticle();
  // set the couplings
  setCoupling(q2,Psp,Psbar,Pten);
  Complex norm=getNorm();
  // first calculate the spinor vector 
  // low energy convention
  Complex aspin[4],ii(0.,1.);
  if(sp.Wave().Rep()==HaberDRep&&sbar.Wave().Rep()==HaberDRep)
    {
      aspin[3] = sbar.s1()*sp.s1()+sbar.s2()*sp.s2()
	-sbar.s3()*sp.s3()-sbar.s4()*sp.s4();
    }
  // high energy convention
  else if(sp.Wave().Rep()==HELASDRep&&sbar.Wave().Rep()==HELASDRep)
    {
      aspin[3] = sbar.s1()*sp.s3()+sbar.s2()*sp.s4()
	+sbar.s3()*sp.s1()+sbar.s4()*sp.s2();
    }
  else
    {
      sp.Wave().changeRep(HELASDRep);
      sbar.Wave().changeRep(HELASDRep);
      aspin[3] = sbar.s1()*sp.s3()+sbar.s2()*sp.s4()
	+sbar.s3()*sp.s1()+sbar.s4()*sp.s2();
    }
  // spatial components are the same in both conventions
  aspin[0] =     +sbar.s1()*sp.s4()+sbar.s2()*sp.s3()
    -sbar.s3()*sp.s2()-sbar.s4()*sp.s1();
  aspin[1] = ii*(-sbar.s1()*sp.s4()+sbar.s2()*sp.s3()
		 +sbar.s3()*sp.s2()-sbar.s4()*sp.s1());
  aspin[2] =     +sbar.s1()*sp.s3()-sbar.s2()*sp.s4()
    -sbar.s3()*sp.s1()+sbar.s4()*sp.s2();
  // difference of spinor momenta
  Energy diff[4]={sp.px()-sbar.px(),sp.py()-sbar.py(),
		  sp.pz()-sbar.pz(),sp.e() -sbar.e()};
  // trace of polarization tensor
  Complex trace = ten.tt()-ten.xx()-ten.yy()-ten.zz();
  // dot products with polarization tensor
  Complex dot[4];
  for(int ix=0;ix<4;++ix)
    {
      dot[ix] = 
	+ten(ix,3)*diff[3]-ten(ix,0)*diff[0]
	-ten(ix,1)*diff[1]-ten(ix,2)*diff[2]
	+ten(3,ix)*diff[3]-ten(0,ix)*diff[0]
	-ten(1,ix)*diff[1]-ten(2,ix)*diff[2]
	-2.*trace*diff[ix];
    }
  // product of spinors
  Complex ffbar= 
    sbar.s1()*sp.s1()+sbar.s2()*sp.s2()
    +sbar.s3()*sp.s3()+sbar.s4()*sp.s4();
  // put everything together
  Complex vertex=-0.125*ii*norm*( aspin[3]*dot[3]-aspin[0]*dot[0]
					  -aspin[1]*dot[1]-aspin[2]*dot[2]
					  +4.0*(Psp->mass())*trace*ffbar);
  return vertex;
}

// member function to evaluate an off-shell spinor
SpinorWaveFunction FFTVertex::evaluate(Energy2 q2, int iopt, tcPDPtr out,
				       const SpinorWaveFunction & sp,
				       const TensorWaveFunction & ten,
				       DiracRep dirac)
{
  // pointers to the particle data objects
  tcPDPtr Psp=sp.getParticle();
  tcPDPtr Pten = ten.getParticle();
  // momentum of the outgoing fermion
  Lorentz5Momentum pout = Lorentz5Momentum(ten.px()+sp.px(),ten.py()+sp.py(),
					   ten.pz()+sp.pz(),ten.e() +sp.e());   
  // set the couplings
  setCoupling(q2,Psp,out,Pten);
  Complex ii(0.,1.);
  // trace of the tensor
  Complex trace = ten.tt()-ten.xx()-ten.yy()-ten.zz();
  // mass of the fermion
  Energy mass = out->mass();
  // overall factor
  Energy2 p2 = pout.m2();
  Complex fact = 0.125*getNorm()*propagator(iopt,p2,out);
  // compute the vector we need
  Complex vec[4],dot;
  for(int ix=0;ix<4;++ix)
    {
      // evaluate the products we need
      dot = (ten(ix,3)+ten(3,ix))*(pout.e()+sp.e());
      for(int iy=0;iy<3;++iy){dot=dot-(ten(ix,iy)+ten(iy,ix))*(pout[iy]+sp[iy]);}
      vec[ix] = dot-2.*trace*(pout[ix]+sp[ix]);
    }
  // combinations of the vector
  Complex a1p2=vec[0]+ii*vec[1];
  Complex a1m2=vec[0]-ii*vec[1];
  // ensure the correct Dirac representation for the spinor
  sp.Wave().changeRep(dirac);
  // now compute the first stage of the spinor wavefunction
  // low energy
  if(dirac==HaberDRep)
    {
      Complex a0=vec[3];
      Complex a3=vec[2];
      vec[0] = a0*sp.s1()-a3*sp.s3()-a1m2*sp.s4(); 
      vec[1] = a0*sp.s2()+a3*sp.s4()-a1p2*sp.s3();
      vec[2] =-a0*sp.s3()+a3*sp.s1()+a1m2*sp.s2();
      vec[3] =-a0*sp.s4()-a3*sp.s2()+a1p2*sp.s1();
    }
  else if(dirac==HELASDRep)
    // high energy
    {
      Complex a0p3=vec[3]+vec[2];
      Complex a0m3=vec[3]-vec[2];
      vec[0] = a0m3*sp.s3()-a1m2*sp.s4(); 
      vec[1] = a0p3*sp.s4()-a1p2*sp.s3();
      vec[2] = a0p3*sp.s1()+a1m2*sp.s2();
      vec[3] = a0m3*sp.s2()+a1p2*sp.s1();
    }
  if(mass!=0.)
    {
      dot = 4.*mass*trace;
      vec[0] = vec[0] + dot*sp.s1(); 
      vec[1] = vec[1] + dot*sp.s2();
      vec[2] = vec[2] + dot*sp.s3();
      vec[3] = vec[3] + dot*sp.s4();
    }
  // combinations of the momentum
  Complex p1p2=pout.px()+ii*pout.py();
  Complex p1m2=pout.px()-ii*pout.py();
  // finally put everything together as the spinor
  Complex ferm[4];
  // low energy
  if(dirac==HaberDRep)
    {
      double p0 = pout.e();
      double p3 = pout.pz();
      ferm[0] = fact*(+p0*vec[0]-  p3*vec[2]-p1m2*vec[3]);
      ferm[1] = fact*(+p0*vec[1]-p1p2*vec[2]+  p3*vec[3]);
      ferm[2] = fact*(-p0*vec[2]+  p3*vec[0]+p1m2*vec[1]);
      ferm[3] = fact*(-p0*vec[3]+p1p2*vec[0]-  p3*vec[1]);
    }
  // high energy
  else if(dirac==HELASDRep)
    {
      Complex p0p3=pout.e() +   pout.pz();
      Complex p0m3=pout.e() -   pout.pz();
      ferm[0] = fact*(+p0m3*vec[2]-p1m2*vec[3]);
      ferm[1] = fact*(-p1p2*vec[2]+p0p3*vec[3]);
      ferm[2] = fact*(+p0p3*vec[0]+p1m2*vec[1]);
      ferm[3] = fact*(+p1p2*vec[0]+p0m3*vec[1]);
    }
  if(mass!=0.)
    {
      ferm[0] = ferm[0] + fact*(mass*vec[0]);
      ferm[1] = ferm[1] + fact*(mass*vec[1]);
      ferm[2] = ferm[2] + fact*(mass*vec[2]);
      ferm[3] = ferm[3] + fact*(mass*vec[3]);
    }
  // return the wavefunction
  return SpinorWaveFunction(pout,out,ferm[0],ferm[1],ferm[2],ferm[3],dirac);
}


// member function to evaluate an off-shell spinor bar
SpinorBarWaveFunction FFTVertex::evaluate(Energy2 q2, int iopt, tcPDPtr out,
					  const SpinorBarWaveFunction & sbar,
					  const TensorWaveFunction & ten,
					  DiracRep dirac)
{
  // pointers to the particle data objects
  tcPDPtr Pten = ten.getParticle();
  tcPDPtr Psbar= sbar.getParticle();
  // momentum of the outgoing fermion
  Lorentz5Momentum pout = Lorentz5Momentum(ten.px()+sbar.px(),ten.py()+sbar.py(),
					   ten.pz()+sbar.pz(),ten.e() +sbar.e());   
  // set the couplings
  setCoupling(q2,out,Psbar,Pten);
  Complex ii(0.,1.);
  // trace of the tensor
  Complex trace = ten.tt()-ten.xx()-ten.yy()-ten.zz();
  // mass of the fermion
  Energy mass = out->mass();
      // overall factor
  Energy2 p2 = pout.m2();
  Complex fact=0.125*getNorm()*propagator(iopt,p2,out);
  // vector
  Complex vec[4],dot;
  for(int ix=0;ix<4;++ix)
    {
      // evaluate the products we need
      dot = -(ten(ix,3)+ten(3,ix))*(pout.e()+sbar.e());
      for(int iy=0;iy<3;++iy){dot=dot+(ten(ix,iy)+ten(iy,ix))*(pout[iy]+sbar[iy]);}
      vec[ix] = dot+2.*trace*(pout[ix]+sbar[ix]);
    }
  // combinations of the vector
  Complex a1p2=vec[0]+ii*vec[1];
  Complex a1m2=vec[0]-ii*vec[1];
  // ensure the correct Dirac representation for the spinor
  sbar.Wave().changeRep(dirac);
  // now compute the first stage of the spinorbar wavefunction
  // low energy
  if(dirac==HaberDRep)
    {
      Complex a0=vec[3];
      Complex a3 = vec[2];
      vec[0] = +a0*sbar.s1()+a3*sbar.s3()+a1p2*sbar.s4(); 
      vec[1] = +a0*sbar.s2()-a3*sbar.s4()+a1m2*sbar.s3();
      vec[2] = -a0*sbar.s3()-a3*sbar.s1()-a1p2*sbar.s2();
      vec[3] = -a0*sbar.s4()+a3*sbar.s2()-a1m2*sbar.s1();
    }
  // high energy
  else if(dirac==HELASDRep)
    {
      Complex a0p3=vec[3]+vec[2];
      Complex a0m3=vec[3]-vec[2];
      vec[0] = a0p3*sbar.s3()+a1p2*sbar.s4(); 
      vec[1] = a0m3*sbar.s4()+a1m2*sbar.s3();
      vec[2] = a0m3*sbar.s1()-a1p2*sbar.s2();
      vec[3] = a0p3*sbar.s2()-a1m2*sbar.s1();
    }
  if(mass!=0.)
    {
      dot =+4.*mass*trace;
      vec[0] = vec[0] + dot*sbar.s1(); 
      vec[1] = vec[1] + dot*sbar.s2();
      vec[2] = vec[2] + dot*sbar.s3();
      vec[3] = vec[3] + dot*sbar.s4();
    }
  // combinations of the momentum
  Complex p1p2=pout.px()+ii*pout.py();
  Complex p1m2=pout.px()-ii*pout.py();
  // finally put everything together as the spinor
  Complex ferm[4];
  // low energy
  if(dirac==HaberDRep)
    {
      double p3=pout.pz();
      double p0=pout.e();
      ferm[0] = fact*(-p0*vec[0]-  p3*vec[2]-p1p2*vec[3]);
      ferm[1] = fact*(-p0*vec[1]-p1m2*vec[2]+  p3*vec[3]);
      ferm[2] = fact*(+p0*vec[2]+  p3*vec[0]+p1p2*vec[1]);
      ferm[3] = fact*(+p0*vec[3]+p1m2*vec[0]-  p3*vec[1]);
    }
  // high energy
  else if(dirac==HELASDRep)
    {
      Complex p0p3=pout.e() +   pout.pz();
      Complex p0m3=pout.e() -   pout.pz();
      ferm[0] = fact*(-p0p3*vec[2]-p1p2*vec[3]);
      ferm[1] = fact*(-p1m2*vec[2]-p0m3*vec[3]);
      ferm[2] = fact*(-p0m3*vec[0]+p1p2*vec[1]);
      ferm[3] = fact*(+p1m2*vec[0]-p0p3*vec[1]);
    }
  if(mass!=0.)
    {
      ferm[0] = ferm[0] + fact*mass*vec[0];
      ferm[1] = ferm[1] + fact*mass*vec[1];
      ferm[2] = ferm[2] + fact*mass*vec[2];
      ferm[3] = ferm[3] + fact*mass*vec[3];
    }
  // return the wavefunction
  return SpinorBarWaveFunction(pout,out,ferm[0],ferm[1],ferm[2],ferm[3],dirac);
}

// member function to evaluate an off-shell tensor
TensorWaveFunction FFTVertex::evaluate(Energy2 q2, int iopt, tcPDPtr out,
				       const SpinorWaveFunction & sp,
				       const SpinorBarWaveFunction & sbar)
{
  // pointers to the particle data objects
  tcPDPtr Psp = sp.getParticle();
  tcPDPtr Psbar = sbar.getParticle();
  // calculating the couplings
  setCoupling(q2,Psp,Psbar,out);
  Complex norm=getNorm();
  Complex ii(0.,1.);
  // momentum of the outgoing tensor
  Lorentz5Momentum pout = Lorentz5Momentum(sp.px()+sbar.px(),sp.py()+sbar.py(),
					   sp.pz()+sbar.pz(),sp.e() +sbar.e());   
  // calculate the prefactor
  Energy2 p2=pout.m2();
  Complex fact=0.125*norm*propagator(iopt,p2,out);
  Energy mass = out->mass();
  Energy2 mass2=mass*mass;
  // spinor vector
  Complex aspin[4];
  if(sp.Wave().Rep()==HaberDRep&&sbar.Wave().Rep()==HaberDRep)
    {
      aspin[3] = sbar.s1()*sp.s1()+sbar.s2()*sp.s2()
	-sbar.s3()*sp.s3()-sbar.s4()*sp.s4();
    }
  // high energy convention
  else if(sp.Wave().Rep()==HELASDRep&&sbar.Wave().Rep()==HELASDRep)
    {
      aspin[3] = sbar.s1()*sp.s3()+sbar.s2()*sp.s4()
	+sbar.s3()*sp.s1()+sbar.s4()*sp.s2();
    }
  else
    {
      sp.Wave().changeRep(HELASDRep);
      sbar.Wave().changeRep(HELASDRep);
      aspin[3] = sbar.s1()*sp.s3()+sbar.s2()*sp.s4()
	+sbar.s3()*sp.s1()+sbar.s4()*sp.s2();
    }
  // spatial components are the same in both conventions
  aspin[0] =     +sbar.s1()*sp.s4()+sbar.s2()*sp.s3()
    -sbar.s3()*sp.s2()-sbar.s4()*sp.s1();
  aspin[1] = ii*(-sbar.s1()*sp.s4()+sbar.s2()*sp.s3()
		 +sbar.s3()*sp.s2()-sbar.s4()*sp.s1());
  aspin[2] =     +sbar.s1()*sp.s3()-sbar.s2()*sp.s4()
    -sbar.s3()*sp.s1()+sbar.s4()*sp.s2();
  // mass dependent term
  Complex ffbar;
  if(Psp->mass()!=0.)
    {
      ffbar = (Psp->mass())*
	(sp.s1()*sbar.s1()+sp.s2()*sbar.s2()+sp.s3()*sbar.s3()+sp.s4()*sbar.s4());
    }
  else
    {
      ffbar = 0.;
    }
  // dot products for the calculation
  Complex dotka = 
    +aspin[3]*pout.e() -aspin[0]*pout.px()
    -aspin[1]*pout.py()-aspin[2]*pout.pz();
  Complex dot12a = 
    +aspin[3]*(sp.e() -sbar.e() )-aspin[0]*(sp.px()-sbar.px())
    -aspin[1]*(sp.py()-sbar.py())-aspin[2]*(sp.pz()-sbar.pz());
  Complex diff=sp.m2()-sbar.m2();
  Complex dotkam=dotka/mass2;
  Complex diffm =diff/mass2;
  Complex p2m = p2/mass2;
  // construct the vectors for the first two terms
  Complex veca[4],vecb[4];
  for(int ix=0;ix<4;++ix)
    {
      veca[ix] = aspin[ix]-dotka*pout[ix];
      vecb[ix] = sp[ix]-sbar[ix]-diffm*pout[ix];
    }
  // coefficients fr hte second two terms
  Complex coeff1 = -4./3.*(2.*ffbar*(1.-p2m)+p2m*dot12a-dotkam*diff);
  Complex coeff2 = -4./3./mass2*( 4.*ffbar*(1.-p2m)
					  -3.*dot12a+2.*p2m*dot12a+diffm*dotka);
  // construct the tensor
  Complex ten[4][4];
  for(int ix=0;ix<4;++ix)
    {
      for(int iy=0;iy<4;++iy)
	{
	  ten[ix][iy] = 2.*(veca[ix]*vecb[iy]+veca[iy]*vecb[ix])
	    +coeff2*pout[ix]*pout[iy];
	}
    }
  ten[0][0]=ten[0][0]-coeff1;
  ten[1][1]=ten[1][1]-coeff1;
  ten[2][2]=ten[2][2]-coeff1;
  ten[3][3]=ten[3][3]+coeff1;
  // multiply by final prefactor
  for(int ix=0;ix<4;++ix)
    {for(int iy=0;iy<4;++iy){ten[ix][iy] = fact*ten[ix][iy];}}
  // return the wavefunction
  return TensorWaveFunction(pout,out,
			    ten[0][0],ten[0][1],ten[0][2],ten[0][3],
			    ten[1][0],ten[1][1],ten[1][2],ten[1][3],
			    ten[2][0],ten[2][1],ten[2][2],ten[2][3],
			    ten[3][0],ten[3][1],ten[3][2],ten[3][3]);
}
    
}
}

