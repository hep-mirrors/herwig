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
    
// Definition of the static class description member
AbstractNoPIOClassDescription<FFTVertex> FFTVertex::initFFTVertex;
    
void FFTVertex::Init() {
  
  static ClassDocumentation<FFTVertex> documentation
    ("The FFTVertex class is the implementation of"
     "the fermion-antifermion tensor vertices for helicity "
     "amplitude calculations. All such vertices should inherit"
     "from it.");
  
}

// function to evaluate the vertex
Complex FFTVertex::evaluate(Energy2 q2,const SpinorWaveFunction & sp,
				    const SpinorBarWaveFunction & sbar,
				    const TensorWaveFunction & ten) {
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
  LorentzSpinorBar<double> sbart=sbar.wave();
  LorentzSpinor<double>    spt  =sp.wave();
  if(sp.wave().Rep()==HaberDRep&&sbar.wave().Rep()==HaberDRep) {
    aspin[3] = sbart.s1()*spt.s1()+sbart.s2()*spt.s2()
              -sbart.s3()*spt.s3()-sbart.s4()*spt.s4();
  }
  // high energy convention
  else if(sp.wave().Rep()==HELASDRep&&sbar.wave().Rep()==HELASDRep) {
    aspin[3] = sbart.s1()*spt.s3()+sbart.s2()*spt.s4()
              +sbart.s3()*spt.s1()+sbart.s4()*spt.s2();
  }
  else {
    sbart.changeRep(HELASDRep);
    spt.changeRep(HELASDRep);
    aspin[3] = sbart.s1()*spt.s3()+sbart.s2()*spt.s4()
              +sbart.s3()*spt.s1()+sbart.s4()*spt.s2();
  }
  // spatial components are the same in both conventions
  aspin[0] =     +sbart.s1()*spt.s4()+sbart.s2()*spt.s3()
                 -sbart.s3()*spt.s2()-sbart.s4()*spt.s1();
  aspin[1] = ii*(-sbart.s1()*spt.s4()+sbart.s2()*spt.s3()
		 +sbart.s3()*spt.s2()-sbart.s4()*spt.s1());
  aspin[2] =     +sbart.s1()*spt.s3()-sbart.s2()*spt.s4()
                 -sbart.s3()*spt.s1()+sbart.s4()*spt.s2();
  // difference of spinor momenta
  Energy diff[4]={sp.px()-sbar.px(),sp.py()-sbar.py(),
		  sp.pz()-sbar.pz(),sp.e() -sbar.e()};
  // trace of polarization tensor
  Complex trace = ten.tt()-ten.xx()-ten.yy()-ten.zz();
  // dot products with polarization tensor
  complex<Energy> dot[4];
  for(int ix=0;ix<4;++ix) {
    dot[ix] = 
      +ten(ix,3)*diff[3]-ten(ix,0)*diff[0]
      -ten(ix,1)*diff[1]-ten(ix,2)*diff[2]
      +ten(3,ix)*diff[3]-ten(0,ix)*diff[0]
      -ten(1,ix)*diff[1]-ten(2,ix)*diff[2]
      -2.*trace*diff[ix];
  }
  // product of spinors
  Complex ffbar=  sbar.s1()*sp.s1()+sbar.s2()*sp.s2()
                 +sbar.s3()*sp.s3()+sbar.s4()*sp.s4();
  // put everything together
  Complex vertex =
    -0.125*ii*norm*UnitRemoval::InvE*( aspin[3]*dot[3]
				       -aspin[0]*dot[0]
				       -aspin[1]*dot[1]
				       -aspin[2]*dot[2]
				       +4.0*
				       (Psp->mass())*trace*ffbar);
  return vertex;
}

// member function to evaluate an off-shell spinor
SpinorWaveFunction FFTVertex::evaluate(Energy2 q2, int iopt, tcPDPtr out,
				       const SpinorWaveFunction & sp,
				       const TensorWaveFunction & ten,
				       DiracRep dirac) {
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
  complex<Energy> dot[4];
  for(int ix=0;ix<4;++ix)
    {
      // evaluate the products we need
      dot[ix] =(ten(ix,3)+ten(3,ix))*(pout.e()+sp.e());
      dot[ix]-=(ten(ix,0)+ten(0,ix))*(pout.x()+sp.px());
      dot[ix]-=(ten(ix,1)+ten(1,ix))*(pout.y()+sp.py());
      dot[ix]-=(ten(ix,2)+ten(2,ix))*(pout.z()+sp.pz());
    }
  LorentzVector<complex<Energy> > vec(dot[0],dot[1],dot[2],dot[3]);
  vec -= 2.0 * trace * (pout + sp.getMomentum());
  // combinations of the vector
  complex<Energy> a1p2=vec.x()+ii*vec.y();
  complex<Energy> a1m2=vec.x()-ii*vec.y();
  // ensure the correct Dirac representation for the spinor
  LorentzSpinor<double>    spt  =sp  .wave().transformRep(dirac);
  // now compute the first stage of the spinor wavefunction
  // low energy
  if(dirac==HaberDRep)
    {
      complex<Energy> a0=vec.t();
      complex<Energy> a3=vec.z();
      vec.setX( a0*spt.s1()-a3*spt.s3()-a1m2*spt.s4()); 
      vec.setY( a0*spt.s2()+a3*spt.s4()-a1p2*spt.s3());
      vec.setZ(-a0*spt.s3()+a3*spt.s1()+a1m2*spt.s2());
      vec.setT(-a0*spt.s4()-a3*spt.s2()+a1p2*spt.s1());
    }
  else if(dirac==HELASDRep)
    // high energy
    {
      complex<Energy> a0p3=vec.t()+vec.z();
      complex<Energy> a0m3=vec.t()-vec.z();
      vec.setX(a0m3*spt.s3()-a1m2*spt.s4()); 
      vec.setY(a0p3*spt.s4()-a1p2*spt.s3());
      vec.setZ(a0p3*spt.s1()+a1m2*spt.s2());
      vec.setT(a0m3*spt.s2()+a1p2*spt.s1());
    }
  if(mass!=Energy())
    {
      complex<Energy> dot = 4.*mass*trace;
      vec.setX(vec.x() + dot*spt.s1()); 
      vec.setY(vec.y() + dot*spt.s2());
      vec.setZ(vec.z() + dot*spt.s3());
      vec.setT(vec.t() + dot*spt.s4());
    }
  // combinations of the momentum
  complex<Energy> p1p2=pout.x()+ii*pout.y();
  complex<Energy> p1m2=pout.x()-ii*pout.y();
  // finally put everything together as the spinor
  Complex ferm[4];
  // low energy
  if(dirac==HaberDRep)
    {
      Energy p0 = pout.e();
      Energy p3 = pout.z();
      ferm[0] = UnitRemoval::InvE2 * fact*( p0*vec.x()-  p3*vec.z()-p1m2*vec.t());
      ferm[1] = UnitRemoval::InvE2 * fact*( p0*vec.y()-p1p2*vec.z()+  p3*vec.t());
      ferm[2] = UnitRemoval::InvE2 * fact*(-p0*vec.z()+  p3*vec.x()+p1m2*vec.y());
      ferm[3] = UnitRemoval::InvE2 * fact*(-p0*vec.t()+p1p2*vec.x()-  p3*vec.y());
    }
  // high energy
  else if(dirac==HELASDRep)
    {
      complex<Energy> p0p3=pout.e() +   pout.z();
      complex<Energy> p0m3=pout.e() -   pout.z();
      ferm[0] = UnitRemoval::InvE2 * fact*( p0m3*vec.z()-p1m2*vec.t());
      ferm[1] = UnitRemoval::InvE2 * fact*(-p1p2*vec.z()+p0p3*vec.t());
      ferm[2] = UnitRemoval::InvE2 * fact*( p0p3*vec.x()+p1m2*vec.y());
      ferm[3] = UnitRemoval::InvE2 * fact*( p1p2*vec.x()+p0m3*vec.y());
    }
  if(mass!=Energy())
    {
      ferm[0] += UnitRemoval::InvE2 * fact*(mass*vec.x());
      ferm[1] += UnitRemoval::InvE2 * fact*(mass*vec.y());
      ferm[2] += UnitRemoval::InvE2 * fact*(mass*vec.z());
      ferm[3] += UnitRemoval::InvE2 * fact*(mass*vec.t());
    }
  // return the wavefunction
  return SpinorWaveFunction(pout,out,ferm[0],ferm[1],ferm[2],ferm[3],dirac);
}


// member function to evaluate an off-shell spinor bar
SpinorBarWaveFunction FFTVertex::evaluate(Energy2 q2, int iopt, tcPDPtr out,
					  const SpinorBarWaveFunction & sbar,
					  const TensorWaveFunction & ten,
					  DiracRep dirac) {
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
  complex<Energy> dot[4];
  for(int ix=0;ix<4;++ix) {
      // evaluate the products we need
      dot[ix] =-(ten(ix,3)+ten(3,ix))*(pout.e()+sbar.e());
      dot[ix]+= (ten(ix,0)+ten(0,ix))*(pout.x()+sbar.px());
      dot[ix]+= (ten(ix,1)+ten(1,ix))*(pout.y()+sbar.py());
      dot[ix]+= (ten(ix,2)+ten(2,ix))*(pout.z()+sbar.pz());
    }
  LorentzVector<complex<Energy> > vec(dot[0],dot[1],dot[2],dot[3]);
  vec += 2.*trace*(pout+sbar.getMomentum());
  // combinations of the vector
  complex<Energy> a1p2=vec.x()+ii*vec.y();
  complex<Energy> a1m2=vec.x()-ii*vec.y();
  // ensure the correct Dirac representation for the spinor
  LorentzSpinorBar<double> sbart=sbar.wave().transformRep(HELASDRep);
  // now compute the first stage of the spinorbar wavefunction
  // low energy
  if(dirac==HaberDRep)
    {
      complex<Energy> a0=vec.t();
      complex<Energy> a3 = vec.z();
      vec.setX( a0*sbart.s1()+a3*sbart.s3()+a1p2*sbart.s4()); 
      vec.setY( a0*sbart.s2()-a3*sbart.s4()+a1m2*sbart.s3());
      vec.setZ(-a0*sbart.s3()-a3*sbart.s1()-a1p2*sbart.s2());
      vec.setT(-a0*sbart.s4()+a3*sbart.s2()-a1m2*sbart.s1());
    }
  // high energy
  else if(dirac==HELASDRep)
    {
      complex<Energy> a0p3=vec.t()+vec.z();
      complex<Energy> a0m3=vec.t()-vec.z();
      vec.setX(a0p3*sbart.s3()+a1p2*sbart.s4()); 
      vec.setY(a0m3*sbart.s4()+a1m2*sbart.s3());
      vec.setZ(a0m3*sbart.s1()-a1p2*sbart.s2());
      vec.setT(a0p3*sbart.s2()-a1m2*sbart.s1());
    }
  if(mass!=Energy())
    {
      complex<Energy> dot = 4.*mass*trace;
      vec.setX(vec.x() + dot*sbart.s1()); 
      vec.setY(vec.y() + dot*sbart.s2());
      vec.setZ(vec.z() + dot*sbart.s3());
      vec.setT(vec.t() + dot*sbart.s4());
    }
  // combinations of the momentum
  complex<Energy>  p1p2=pout.x()+ii*pout.y();
  complex<Energy>  p1m2=pout.x()-ii*pout.y();
  // finally put everything together as the spinor
  Complex ferm[4];
  // low energy
  if(dirac==HaberDRep)
    {
      Energy p3=pout.z();
      Energy p0=pout.e();
      ferm[0] = UnitRemoval::InvE2 * fact*(-p0*vec.x()-  p3*vec.z()-p1p2*vec.t());
      ferm[1] = UnitRemoval::InvE2 * fact*(-p0*vec.y()-p1m2*vec.z()+  p3*vec.t());
      ferm[2] = UnitRemoval::InvE2 * fact*( p0*vec.z()+  p3*vec.x()+p1p2*vec.y());
      ferm[3] = UnitRemoval::InvE2 * fact*( p0*vec.t()+p1m2*vec.x()-  p3*vec.y());
    }
  // high energy
  else if(dirac==HELASDRep)
    {
      complex<Energy> p0p3=pout.e() +   pout.z();
      complex<Energy> p0m3=pout.e() -   pout.z();
      ferm[0] = UnitRemoval::InvE2 * fact*(-p0p3*vec.z()-p1p2*vec.t());
      ferm[1] = UnitRemoval::InvE2 * fact*(-p1m2*vec.z()-p0m3*vec.t());
      ferm[2] = UnitRemoval::InvE2 * fact*(-p0m3*vec.x()+p1p2*vec.y());
      ferm[3] = UnitRemoval::InvE2 * fact*( p1m2*vec.x()-p0p3*vec.y());
    }
  if(mass!=Energy())
    {
      ferm[0] += UnitRemoval::InvE2 * fact*mass*vec.x();
      ferm[1] += UnitRemoval::InvE2 * fact*mass*vec.y();
      ferm[2] += UnitRemoval::InvE2 * fact*mass*vec.z();
      ferm[3] += UnitRemoval::InvE2 * fact*mass*vec.t();
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
  LorentzSpinorBar<double> sbart=sbar.wave();
  LorentzSpinor<double>    spt  =sp  .wave();
  if(sp.wave().Rep()==HaberDRep&&sbar.wave().Rep()==HaberDRep) {
    aspin[3] = sbart.s1()*spt.s1()+sbart.s2()*spt.s2()
              -sbart.s3()*spt.s3()-sbart.s4()*spt.s4();
  }
  // high energy convention
  else if(sp.wave().Rep()==HELASDRep&&sbar.wave().Rep()==HELASDRep) {
    aspin[3] = sbart.s1()*spt.s3()+sbart.s2()*spt.s4()
              +sbart.s3()*spt.s1()+sbart.s4()*spt.s2();
  }
  else {
    sbart.changeRep(HELASDRep);
    spt.changeRep(HELASDRep);
    aspin[3] = sbart.s1()*spt.s3()+sbart.s2()*spt.s4()
	      +sbart.s3()*spt.s1()+sbart.s4()*spt.s2();
    }
  // spatial components are the same in both conventions
  aspin[0] =     +sbart.s1()*spt.s4()+sbart.s2()*spt.s3()
                 -sbart.s3()*spt.s2()-sbart.s4()*spt.s1();
  aspin[1] = ii*(-sbart.s1()*spt.s4()+sbart.s2()*spt.s3()
		 +sbart.s3()*spt.s2()-sbart.s4()*spt.s1());
  aspin[2] =     +sbart.s1()*spt.s3()-sbart.s2()*spt.s4()
                 -sbart.s3()*spt.s1()+sbart.s4()*spt.s2();
  // mass dependent term
  Complex ffbar;
  if(Psp->mass()!=Energy())
    {
      ffbar = UnitRemoval::InvE * (Psp->mass())*
	(sp.s1()*sbar.s1()+sp.s2()*sbar.s2()+sp.s3()*sbar.s3()+sp.s4()*sbar.s4());
    }
  else
    {
      ffbar = 0.;
    }
  // dot products for the calculation
  complex<Energy> dotka = 
    +aspin[3]*pout.e()-aspin[0]*pout.x()
    -aspin[1]*pout.y()-aspin[2]*pout.z();
  complex<Energy> dot12a = 
    +aspin[3]*(sp.e() -sbar.e() )-aspin[0]*(sp.px()-sbar.px())
    -aspin[1]*(sp.py()-sbar.py())-aspin[2]*(sp.pz()-sbar.pz());
  complex<Energy2> diff=sp.m2()-sbar.m2();
  complex<InvEnergy> dotkam=dotka/mass2;
  Complex diffm =diff/mass2;
  double p2m = p2/mass2;
  // construct the vectors for the first two terms
  Complex veca[4],vecb[4];

  veca[0] = aspin[0]-UnitRemoval::InvE2*dotka*pout.x();
  vecb[0] = UnitRemoval::InvE*(sp.px()-sbar.px()-diffm*pout.x());
  
  veca[1] = aspin[1]-UnitRemoval::InvE2*dotka*pout.y();
  vecb[1] = UnitRemoval::InvE*(sp.py()-sbar.py()-diffm*pout.y());
  
  veca[2] = aspin[2]-UnitRemoval::InvE2*dotka*pout.z();
  vecb[2] = UnitRemoval::InvE*(sp.pz()-sbar.pz()-diffm*pout.z());
  
  veca[3] = aspin[3]-UnitRemoval::InvE2*dotka*pout.e();
  vecb[3] = UnitRemoval::InvE*(sp.e()-sbar.e()-diffm*pout.e());
  
  // coefficients fr hte second two terms
  Complex temp = UnitRemoval::InvE*(p2m*dot12a-dotkam*diff);
  Complex coeff1 = -4./3.*(2.*ffbar*(1.-p2m) + temp);
  temp = UnitRemoval::InvE*(-3.*dot12a+2.*p2m*dot12a+diffm*dotka);
  Complex coeff2 = -4./3./mass2*( 4.*ffbar*(1.-p2m) + temp)*UnitRemoval::E2;
  // construct the tensor
  Complex ten[4][4];
  
  const complex<Energy> pout_tmp[4] 
    = {pout.x(), pout.y(), pout.z(), pout.e()};

  for(int ix=0;ix<4;++ix)
    {
      for(int iy=0;iy<4;++iy)
	{
	  Complex temp = coeff2*pout_tmp[ix]*pout_tmp[iy]*UnitRemoval::InvE2;
	  ten[ix][iy] = 2.*(veca[ix]*vecb[iy]+veca[iy]*vecb[ix]) + temp;
	}
    }
  
  ten[0][0]=ten[0][0]-coeff1;
  ten[1][1]=ten[1][1]-coeff1;
  ten[2][2]=ten[2][2]-coeff1;
  ten[3][3]=ten[3][3]+coeff1;
  // multiply by final prefactor
  for(int ix=0;ix<4;++ix) {
    for(int iy=0;iy<4;++iy) {
      ten[ix][iy] = fact*ten[ix][iy];
    }
  }
  // return the wavefunction
  return TensorWaveFunction(pout,out,
			    ten[0][0],ten[0][1],ten[0][2],ten[0][3],
		 	    ten[1][0],ten[1][1],ten[1][2],ten[1][3],
			    ten[2][0],ten[2][1],ten[2][2],ten[2][3],
			    ten[3][0],ten[3][1],ten[3][2],ten[3][3]);
}
    
}
}

