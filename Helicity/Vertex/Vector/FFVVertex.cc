// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFVVertex class.
//

#include "FFVVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
namespace Helicity{

using namespace ThePEG;
using ThePEG::Helicity::DiracRep;
using ThePEG::Helicity::HELASDRep;
using ThePEG::Helicity::HaberDRep;

FFVVertex::~FFVVertex() {}
    
void FFVVertex::persistentOutput(PersistentOStream & os) const {;}
    
void FFVVertex::persistentInput(PersistentIStream & is, int i) {;}
   
// Definition of the static class description member
ClassDescription<FFVVertex> FFVVertex::initFFVVertex;
    
void FFVVertex::Init() {
      
  static ClassDocumentation<FFVVertex> documentation
    ("The \\classname{FFVVertex} class implements the helicity amplitude"
     "calculations for a fermion-fantifermion gauge boson vertex. Any   "
     "implementation of such a vertex should inherit from in and implement"
     " the virtual setCoupling member to calculate the coupling");
}

// duumy setCoupling member
void FFVVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b,tcPDPtr c){;}

// evalulate the full vertex
Complex FFVVertex::evaluate(Energy2 q2,
			    const SpinorWaveFunction & sp, 
			    const SpinorBarWaveFunction & sbar,
			    const VectorWaveFunction & vec)
{
  // extract the pointers to the particle data objects
  tcPDPtr  Psp=sp.getParticle();
  tcPDPtr  Pvec=vec.getParticle();
  tcPDPtr  Psbar=sbar.getParticle();
  // first calculate the couplings
  setCoupling(q2,Psp,Psbar,Pvec);
  Complex ii(0.,1.);
  Complex vertex(0.);
  // useful combinations of the polarization vector components
  Complex e0p3=vec.t()+vec.z();
  Complex e0m3=vec.t()-vec.z();
  Complex e1p2=vec.x()+ii*vec.y();
  Complex e1m2=vec.x()-ii*vec.y();
  // work out which convention to use 
  // low energy convention
  if(sp.Wave().Rep()==HaberDRep&&sbar.Wave().Rep()==HaberDRep)
    {      
      Complex st[4]={0.,0.,0.,0.};
      // first the left piece as this is virtually always needed
      if(_left!=0.)
	{
	  Complex s1m3=sp.s1()-sp.s3();
	  Complex s2m4=sp.s2()-sp.s4();
	  st[0] = _left*(e0p3*s1m3+e1m2*s2m4);
	  st[1] = _left*(e1p2*s1m3+e0m3*s2m4);
	}
      // then the right piece (often not needed eg W vertex)
      if(_right!=0.)
	{
	  Complex s1p3=sp.s1()+sp.s3();
	  Complex s2p4=sp.s2()+sp.s4();
	  st[2] = _right*( e0m3*s1p3-e1m2*s2p4);
	  st[3] = _right*(-e1p2*s1p3+e0p3*s2p4);
	}
      vertex = ii*0.5*
	(  sbar.s1()*(st[0]+st[2])+sbar.s2()*(st[1]+st[3])
	   +sbar.s3()*(st[0]-st[2])+sbar.s4()*(st[1]-st[3]));
    }
  // high energy convention
  else if(sp.Wave().Rep()==HELASDRep&&sbar.Wave().Rep()==HELASDRep)
    {
      // first the left piece as this is virtually always needed
      if(_left!=0.)
	{
	  vertex = _left*(+sbar.s3()*(sp.s1()*e0p3+sp.s2()*e1m2)
			  +sbar.s4()*(sp.s1()*e1p2+sp.s2()*e0m3)); 
	}
      // then the right piece (often not needed eg W vertex)
      if(_right!=0.)
	{
	  vertex += _right*(+sbar.s1()*(sp.s3()*e0m3-sp.s4()*e1m2)
			    -sbar.s2()*(sp.s3()*e1p2-sp.s4()*e0p3));
	}
      vertex*=ii;
    }
  // mixed conventions
  else
    {
      // get both the spinors in the same representation (using the default choice)
        sp.Wave().changeRep(HELASDRep);
      sbar.Wave().changeRep(HELASDRep);
      // first the left piece as this is virtually always needed
      if(_left!=0.)
	{
	  vertex = _left*(+sbar.s3()*(sp.s1()*e0p3+sp.s2()*e1m2)
			  +sbar.s4()*(sp.s1()*e1p2+sp.s2()*e0m3)); 
	}
      // then the right piece (often not needed eg W vertex)
      if(_right!=0.)
	{
	  vertex += _right*(+sbar.s1()*(sp.s3()*e0m3-sp.s4()*e1m2)
			    -sbar.s2()*(sp.s3()*e1p2-sp.s4()*e0p3));
	}
      vertex*=ii;
    }
  // final factors
  vertex *= getNorm();
  return vertex;
}
// evaluate an off-shell spinor
SpinorWaveFunction FFVVertex::evaluate(Energy2 q2, int iopt,tcPDPtr  out,
				       const SpinorWaveFunction & sp,
				       const VectorWaveFunction &vec,DiracRep dirac)
{
  // extract the pointers to the particle data objects
  tcPDPtr  Psp=sp.getParticle();
  tcPDPtr  Pvec=vec.getParticle();
  // first calculate the couplings
  setCoupling(q2,Psp,out,Pvec);
  Complex ii(0.,1.);
  // work out the momentum of the off-shell particle
  Lorentz5Momentum pout = Lorentz5Momentum(sp.px()+vec.px(),sp.py()+vec.py(),
					   sp.pz()+vec.pz(),sp.e() +vec.e() ); 
  // now evaluate the contribution
  // polarization components
  Complex e0p3=vec.t() +  vec.z();
  Complex e0m3=vec.t() -  vec.z();
  Complex e1p2=vec.x()+ii*vec.y();
  Complex e1m2=vec.x()-ii*vec.y();
  // momentum components
  double mass  = out->mass();
  Complex p1p2=pout.px()+ii*pout.py();
  Complex p1m2=pout.px()-ii*pout.py();
  // complex nos for for the spinor
  Complex s1(0.),s2(0.),s3(0.),s4(0.);
  // overall factor
  double p2 = pout.m2();
  Complex fact=-normPropagator(iopt,p2,out);
  // ensure the spinor is in the correct dirac representation
  sp.Wave().changeRep(dirac);
  // low energy Dirac matrix defn
  if(dirac==HaberDRep)
    {   
      // first step compute the polarization vector combined with the spinor
      Complex st[4]={0.,0.,0.,0.};
      // first the left piece as this is virtually always needed
      if(_left!=0.)
	{
	  st[0] +=_left*( e0p3*(sp.s1()-sp.s3())+e1m2*(sp.s2()-sp.s4()));
	  st[1] +=_left*( e1p2*(sp.s1()-sp.s3())+e0m3*(sp.s2()-sp.s4()));
	  st[2] +=_left*( e0p3*(sp.s1()-sp.s3())+e1m2*(sp.s2()-sp.s4()));
	  st[3] +=_left*( e1p2*(sp.s1()-sp.s3())+e0m3*(sp.s2()-sp.s4()));
	}
      // then the right piece (often not needed eg W vertex)
      if(_right!=0.)
	{
	  st[0] +=_right*( e0m3*(sp.s1()+sp.s3())-e1m2*(sp.s2()+sp.s4()));
	      st[1] +=_right*(-e1p2*(sp.s1()+sp.s3())+e0p3*(sp.s2()+sp.s4()));
	      st[2] +=_right*(-e0m3*(sp.s1()+sp.s3())+e1m2*(sp.s2()+sp.s4()));
	      st[3] +=_right*( e1p2*(sp.s1()+sp.s3())-e0p3*(sp.s2()+sp.s4()));
	}
      // then combine this with pslash+m
      Complex p0pm=pout.e()+mass;
      Complex p0mm=pout.pz()*pout.pz()/p0pm;
      fact *= 0.5;
      s1 =fact*( p0pm*st[0]-pout.pz()*st[2]-     p1m2*st[3]);
      s2 =fact*( p0pm*st[1]-     p1p2*st[2]+pout.pz()*st[3]);
      s3 =fact*(-p0mm*st[2]+pout.pz()*st[0]+     p1m2*st[1]);
      s4 =fact*(-p0mm*st[3]+     p1p2*st[0]-pout.pz()*st[1]);
    }
  // high energy Dirac matrix defn
  else if(dirac==HELASDRep)
    { 
      Complex p0p3=pout.e() +   pout.pz();
      Complex p0m3=pout.e() -   pout.pz();
      // left piece
      if(_left!=0.)
	{
	  Complex a3=_left*fact*( sp.s1()*e0p3+sp.s2()*e1m2);
	  Complex a4=_left*fact*( sp.s1()*e1p2+sp.s2()*e0m3);
	  s1 +=p0m3*a3-p1m2*a4;
	  s2 +=p1p2*a3+p0p3*a4;
	  s3 +=a3*mass;
	  s4 +=a4*mass;
	}
      // right piece
      if(_right!=0.)
	{
	  Complex a1=_right*fact*( sp.s3()*e0m3-sp.s4()*e1m2);
	  Complex a2=_right*fact*(-sp.s3()*e1p2+sp.s4()*e0p3);
	  s1 +=a1*mass;
	  s2 +=a2*mass;
	  s3 +=p0p3*a1+p1m2*a2;
	  s4 +=p1p2*a1+p0m3*a2;
	}
    }
  // return the wavefunction
  return SpinorWaveFunction(pout,out,s1,s2,s3,s4,dirac);
}

// evaluate an off-shell SpinorBar
SpinorBarWaveFunction FFVVertex::evaluate(Energy2 q2,int iopt,tcPDPtr  out,
					  const SpinorBarWaveFunction & sbar,
					  const VectorWaveFunction & vec, DiracRep dirac)
{
  // extract the pointers to the particle data objects
  tcPDPtr  Psbar=sbar.getParticle();
  tcPDPtr  Pvec=vec.getParticle();
  // first calculate the couplings
  setCoupling(q2,out,Psbar,Pvec);
  Complex ii(0.,1.);
  // work out the momentum of the off-shell particle
  Lorentz5Momentum pout = Lorentz5Momentum(sbar.px()+vec.px(),sbar.py()+vec.py(),
					   sbar.pz()+vec.pz(),sbar.e() +vec.e()); 
  // now evaluate the contribution
  // polarization components
  Complex e0p3=vec.t() +  vec.z();
  Complex e0m3=vec.t() -  vec.z();
  Complex e1p2=vec.x()+ii*vec.y();
  Complex e1m2=vec.x()-ii*vec.y();
  // momentum components
  double mass  = out->mass();
  Complex p1p2=pout.px()+ii*pout.py();
  Complex p1m2=pout.px()-ii*pout.py();
  // complex numbers for the spinor
  Complex s1(0.),s2(0.),s3(0.),s4(0.);
  // overall factor
  double p2 = pout.m2();
  Complex fact=-normPropagator(iopt,p2,out);
  // ensure the spinorbar is in the correct dirac representation
  sbar.Wave().changeRep(dirac);
  // low energy convention
  if(dirac==HaberDRep)
    {
      // first step compute the polarization vector combined with the spinorbar
      Complex st[4]={0.,0.,0.,0.};
      // first the left piece as this is virtually always needed
      if(_left!=0.)
	{
	  st[0] += _left*( e0p3*(sbar.s1()+sbar.s3())+e1p2*(sbar.s2()+sbar.s4()));
	  st[1] += _left*( e1m2*(sbar.s1()+sbar.s3())+e0m3*(sbar.s2()+sbar.s4()));
	  st[2] += _left*(-e0p3*(sbar.s1()+sbar.s3())-e1p2*(sbar.s2()+sbar.s4()));
	  st[3] += _left*(-e1m2*(sbar.s1()+sbar.s3())-e0m3*(sbar.s2()+sbar.s4()));
	}
      // then the right piece (often not needed eg W vertex)
      if(_right!=0.)
	{
	  st[0] +=_right*( e0m3*(sbar.s1()-sbar.s3())-e1p2*(sbar.s2()-sbar.s4()));
	  st[1] +=_right*(-e1m2*(sbar.s1()-sbar.s3())+e0p3*(sbar.s2()-sbar.s4()));
	  st[2] +=_right*(-e0m3*(sbar.s1()-sbar.s3())-e1p2*(sbar.s2()-sbar.s4()));
	  st[3] +=_right*( e1m2*(sbar.s1()-sbar.s3())+e0p3*(sbar.s2()-sbar.s4()));
	}
      // then combine this with -pslash+m 
      Complex p0pm=pout.e()+mass;
      Complex p0mm=pout.pz()*pout.pz()/p0pm;
      fact *= 0.5;
      s1 =fact*(-p0mm*st[0]-pout.pz()*st[2]-     p1p2*st[3]);
      s2 =fact*(-p0mm*st[1]-     p1m2*st[2]+pout.pz()*st[3]);
      s3 =fact*( p0pm*st[2]+pout.pz()*st[0]+     p1p2*st[1]);
      s4 =fact*( p0pm*st[3]+     p1m2*st[0]-pout.pz()*st[1]);
    }
  // high energy convention
  else if(dirac==HELASDRep)
    {
      Complex p0p3=pout.e() +   pout.pz();
      Complex p0m3=pout.e() -   pout.pz();
      // left piece
      if(_left!=0.)
	{
	  Complex a1 =  _left*fact*( sbar.s3()*e0p3+sbar.s4()*e1p2);
	  Complex a2 =  _left*fact*( sbar.s3()*e1m2+sbar.s4()*e0m3);
	  s1 += +a1*mass;
	  s2 += +a2*mass;
	  s3 += -p0m3*a1+p1p2*a2;
	  s4 += +p1m2*a1-p0p3*a2;
	}
      // right piece
      if(_right!=0.)
	{
	  Complex a3 = _right*fact*( sbar.s1()*e0m3-sbar.s2()*e1p2);
	  Complex a4 = _right*fact*(-sbar.s1()*e1m2+sbar.s2()*e0p3);
	  s1 += -p0p3*a3-p1p2*a4;
	  s2 += -p1m2*a3-p0m3*a4;
	  s3 += +a3*mass;
	  s4 += +a4*mass;
	}
    }
  return SpinorBarWaveFunction(pout,out,s1,s2,s3,s4,dirac);
}

// off-shell vector
VectorWaveFunction FFVVertex::evaluate(Energy2 q2,int iopt,tcPDPtr  out,
				       const SpinorWaveFunction & sp,
				       const SpinorBarWaveFunction & sbar)
{
  // extract the pointers to the particle data objects
  tcPDPtr  Psbar=sbar.getParticle();
  tcPDPtr  Psp=sp.getParticle();
  // first calculate the couplings
  setCoupling(q2,Psp,Psbar,out);
  Complex ii(0.,1.);
  // work out the momentum of the off-shell particle
  Lorentz5Momentum pout = Lorentz5Momentum(sbar.px()+sp.px(),sbar.py()+sp.py(),
					   sbar.pz()+sp.pz(),sbar.e() +sp.e()); 
  // momentum components
  double mass  = out->mass();
  double mass2=mass*mass;
  // overall factor
  double p2 = pout.m2();
  // the vector for the fermion-antifermion
  Complex vec[4];
  if(sp.Wave().Rep()==HaberDRep&&sbar.Wave().Rep()==HaberDRep)
    {
      // left coupling
      if(_left!=0.)
	{
	  Complex s2m4=sp.s2()-sp.s4();
	  Complex s1m3=sp.s1()-sp.s3();
	  vec[0] =   0.5*_left*(-sbar.s1()*s2m4-sbar.s2()*s1m3
				-sbar.s3()*s2m4-sbar.s4()*s1m3);
	  vec[1] =ii*0.5*_left*(+sbar.s1()*s2m4-sbar.s2()*s1m3
				+sbar.s3()*s2m4-sbar.s4()*s1m3);
	  vec[2] =   0.5*_left*(-sbar.s1()*s1m3+sbar.s2()*s2m4
				-sbar.s3()*s1m3+sbar.s4()*s2m4);
	  vec[3] =   0.5*_left*(+sbar.s1()*s1m3+sbar.s2()*s2m4
				+sbar.s3()*s1m3+sbar.s4()*s2m4);
	}
      // right coupling
      if(_right!=0.)
	{
	  Complex s1p3=sp.s1()+sp.s3();
	  Complex s2p4=sp.s2()+sp.s4();
	  vec[0] +=    +0.5*_right*(+sbar.s1()*s2p4+sbar.s2()*s1p3
				    -sbar.s3()*s2p4-sbar.s4()*s1p3);
	  vec[1] += +ii*0.5*_right*(-sbar.s1()*s2p4+sbar.s2()*s1p3
				    +sbar.s3()*s2p4-sbar.s4()*s1p3);
	  vec[2] +=    +0.5*_right*(+sbar.s1()*s1p3-sbar.s2()*s2p4
				    -sbar.s3()*s1p3+sbar.s4()*s2p4);
	  vec[3] +=    +0.5*_right*(+sbar.s1()*s1p3+sbar.s2()*s2p4
				    -sbar.s3()*s1p3-sbar.s4()*s2p4);
	}
    }
  else if(sp.Wave().Rep()==HELASDRep&&sbar.Wave().Rep()==HELASDRep)
    {
      // left coupling
      if(_left!=0.)
	{
	  vec[0] =   -_left*(sbar.s3()*sp.s2()+sbar.s4()*sp.s1());
	  vec[1] = ii*_left*(sbar.s3()*sp.s2()-sbar.s4()*sp.s1());
	  vec[2] =   -_left*(sbar.s3()*sp.s1()-sbar.s4()*sp.s2());
	  vec[3] =    _left*(sbar.s3()*sp.s1()+sbar.s4()*sp.s2());
	}
      // right coupling
      if(_right!=0.)
	{
	  vec[0] +=    +_right*(sbar.s1()*sp.s4()+sbar.s2()*sp.s3());
	  vec[1] += -ii*_right*(sbar.s1()*sp.s4()-sbar.s2()*sp.s3());
	  vec[2] +=    +_right*(sbar.s1()*sp.s3()-sbar.s2()*sp.s4());
	  vec[3] +=    +_right*(sbar.s1()*sp.s3()+sbar.s2()*sp.s4());
	}
    }
  else
    {
      sbar.Wave().changeRep(HELASDRep);
      sp.Wave().changeRep(HELASDRep);
      // left coupling
      if(_left!=0.)
	{
	  vec[0] =   -_left*(sbar.s3()*sp.s2()+sbar.s4()*sp.s1());
	  vec[1] = ii*_left*(sbar.s3()*sp.s2()-sbar.s4()*sp.s1());
	  vec[2] =   -_left*(sbar.s3()*sp.s1()-sbar.s4()*sp.s2());
	  vec[3] =    _left*(sbar.s3()*sp.s1()+sbar.s4()*sp.s2());
	}
      // right coupling
      if(_right!=0.)
	{
	  vec[0] +=    +_right*(sbar.s1()*sp.s4()+sbar.s2()*sp.s3());
	  vec[1] += -ii*_right*(sbar.s1()*sp.s4()-sbar.s2()*sp.s3());
	  vec[2] +=    +_right*(sbar.s1()*sp.s3()-sbar.s2()*sp.s4());
	  vec[3] +=    +_right*(sbar.s1()*sp.s3()+sbar.s2()*sp.s4());
	}
    }
  // prefactor
  Complex fact = normPropagator(iopt,p2,out);
  // massless boson
  if(mass==0.)
    {
      for(int ix=0;ix<4;++ix){vec[ix]*=fact;}
    }
  // massive boson
  else
    {
      Complex dot = ( pout.e() *vec[3]-pout.px()*vec[0]
			      -pout.py()*vec[1]-pout.pz()*vec[2])/mass2;
      for(int ix=0;ix<4;++ix){vec[ix]=fact*(vec[ix]-dot*pout[ix]);}
    }
  return VectorWaveFunction(pout,out,vec[0],vec[1],vec[2],vec[3]);
}

}
}
