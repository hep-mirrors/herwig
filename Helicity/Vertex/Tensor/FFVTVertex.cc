// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFVTVertex class.
//

#include "FFVTVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;
    
FFVTVertex::~FFVTVertex() {}
    
void FFVTVertex::persistentOutput(PersistentOStream & os) const {}

void FFVTVertex::persistentInput(PersistentIStream & is, int) {}

ClassDescription<FFVTVertex> FFVTVertex::initFFVTVertex;
// Definition of the static class description member.

void FFVTVertex::Init() {
  
  static ClassDocumentation<FFVTVertex> documentation
    ("The \\classname{FFVTVertex} class is the implementation of the"
     "helicity amplitude calculation of the fermion-antifermion-vector-tensor"
     "vertex. All such vertices should inherit from it.");

}
// coupling
void FFVTVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b,tcPDPtr c, tcPDPtr d){;}

// function to evaluate the vertex
Complex FFVTVertex::evaluate(Energy q2, const SpinorWaveFunction & sp,
    				 const SpinorBarWaveFunction & sbar,
    				 const VectorWaveFunction & vec,
    				 const TensorWaveFunction & ten)
{
  // pointers to the particles
  tcPDPtr Psp   = sp.getParticle();
  tcPDPtr Psbar = sbar.getParticle();
  tcPDPtr Pvec  = vec.getParticle();
  tcPDPtr Pten  =  ten.getParticle();
  // set the couplings
  setCoupling(q2,Psp,Psbar,Pvec,Pten);
  Complex norm=getNorm();
  Complex ii(0.,1.);
  // spinor vector
  // low energy convention
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
  // trace of the tensor
  Complex trace=ten.tt()-ten.xx()-ten.yy()-ten.zz();
  // dot product
  Complex dotav = 
    aspin[3]*vec.t()-aspin[0]*vec.x()-aspin[1]*vec.y()-aspin[2]*vec.z();
  // components of the tensor
  Complex tentx = ten.tx()+ten.xt();
  Complex tenty = ten.ty()+ten.yt();
  Complex tentz = ten.tz()+ten.zt();
  Complex tenxy = ten.xy()+ten.yx();
  Complex tenxz = ten.xz()+ten.zx();
  Complex tenyz = ten.yz()+ten.zy();
  // dot product of wavefunctions and momenta with the tensor
  Complex tenav =
    2.*(+ten.tt()*vec.t()*aspin[3]+ten.xx()*vec.x()*aspin[0]
	+ten.yy()*vec.y()*aspin[1]+ten.zz()*vec.z()*aspin[2])
    -tentx*(vec.t()*aspin[0]+vec.x()*aspin[3])
    -tenty*(vec.t()*aspin[1]+vec.y()*aspin[3])
    -tentz*(vec.t()*aspin[2]+vec.z()*aspin[3])
    +tenxy*(vec.x()*aspin[1]+vec.y()*aspin[0])
    +tenxz*(vec.x()*aspin[2]+vec.z()*aspin[0])
    +tenyz*(vec.y()*aspin[2]+vec.z()*aspin[1]);
  // return the vertex
  return ii*0.25*norm*(tenav-2.*trace*dotav);
}

}
}

