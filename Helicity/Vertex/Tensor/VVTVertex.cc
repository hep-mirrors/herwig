// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VVTVertex class.
//

#include "VVTVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
namespace Helicity{
using namespace ThePEG;

AbstractNoPIOClassDescription<VVTVertex> VVTVertex::initVVTVertex;
// Definition of the static class description member.

void VVTVertex::Init() {
  
  
  static ClassDocumentation<VVTVertex> documentation
    ("The VVTVertex class is the implementation of"
     "the vecotr-vector tensor vertices for helicity "
     "amplitude calculations. All such vertices should inherit"
     "from it.");
  
}

// function to evaluate the vertex
Complex VVTVertex::evaluate(Energy q2, const VectorWaveFunction & vec1,
    				const VectorWaveFunction & vec2, 
    				const TensorWaveFunction & ten)
{
  // pointers to the particles
  tcPDPtr Pvec1 = vec1.getParticle();
  tcPDPtr Pvec2 = vec2.getParticle();
  tcPDPtr Pten  =  ten.getParticle();
  // set the couplings
  setCoupling(q2,Pvec1,Pvec2,Pten);
  Complex norm=getNorm();
  Complex ii(0.,1.);
  // mass of the vector
  Energy vmass = Pvec1->mass();
  // mass+k1.k2
  Energy2 mdot =
    + vec1.e()*vec2.e() -vec1.px()*vec2.px()
    -vec1.py()*vec2.py()-vec1.pz()*vec2.pz();
  if(vmass!=0.){mdot=mdot+vmass*vmass;}
  // dot product of wavefunctions and momenta
  Complex dotv1v2  = vec1.t()*vec2.t()- vec1.x()*vec2.x()
    - vec1.y()*vec2.y()- vec1.z()*vec2.z();
  Complex dotk1v2 =  vec1.e()*vec2.t()-vec1.px()*vec2.x()
    -vec1.py()*vec2.y()-vec1.pz()*vec2.z();
  Complex dotk2v1 = vec1.t()*vec2.e() -vec1.x()*vec2.px()
    -vec1.y()*vec2.py()-vec1.z()*vec2.pz();
  // components of the tensor
  Complex tentx = ten.tx()+ten.xt();
  Complex tenty = ten.ty()+ten.yt();
  Complex tentz = ten.tz()+ten.zt();
  Complex tenxy = ten.xy()+ten.yx();
  Complex tenxz = ten.xz()+ten.zx();
  Complex tenyz = ten.yz()+ten.zy();
  // dot product of wavefunctions and momenta with the tensor
  Complex tenv1v2 =
    2.*(+ten.tt()*vec1.t()*vec2.t()+ten.xx()*vec1.x()*vec2.x()
        +ten.yy()*vec1.y()*vec2.y()+ten.zz()*vec1.z()*vec2.z())
    -tentx*(vec1.t()*vec2.x()+vec1.x()*vec2.t())
    -tenty*(vec1.t()*vec2.y()+vec1.y()*vec2.t())
    -tentz*(vec1.t()*vec2.z()+vec1.z()*vec2.t())
    +tenxy*(vec1.x()*vec2.y()+vec1.y()*vec2.x())
    +tenxz*(vec1.x()*vec2.z()+vec1.z()*vec2.x())
    +tenyz*(vec1.y()*vec2.z()+vec1.z()*vec2.y());
  Complex tenk1v2 =
    2.*(+ten.tt()*vec1.e()*vec2.t() +ten.xx()*vec1.px()*vec2.x()
        +ten.yy()*vec1.py()*vec2.y()+ten.zz()*vec1.pz()*vec2.z())
    -tentx*(vec1.e()*vec2.x() +vec1.px()*vec2.t())
    -tenty*(vec1.e()*vec2.y() +vec1.py()*vec2.t())
    -tentz*(vec1.e()*vec2.z() +vec1.pz()*vec2.t())
    +tenxy*(vec1.px()*vec2.y()+vec1.py()*vec2.x())
    +tenxz*(vec1.px()*vec2.z()+vec1.pz()*vec2.x())
    +tenyz*(vec1.py()*vec2.z()+vec1.pz()*vec2.y());
  Complex tenk2v1 =
    2.*(+ten.tt()*vec1.t()*vec2.e() +ten.xx()*vec1.x()*vec2.px()
        +ten.yy()*vec1.y()*vec2.py()+ten.zz()*vec1.z()*vec2.pz())
    -tentx*(vec1.t()*vec2.px()+vec1.x()*vec2.e())
    -tenty*(vec1.t()*vec2.py()+vec1.y()*vec2.e())
    -tentz*(vec1.t()*vec2.pz()+vec1.z()*vec2.e())
    +tenxy*(vec1.x()*vec2.py()+vec1.y()*vec2.px())
    +tenxz*(vec1.x()*vec2.pz()+vec1.z()*vec2.px())
    +tenyz*(vec1.y()*vec2.pz()+vec1.z()*vec2.py());
  Complex tenk1k2 =
    2.*(+ten.tt()*vec1.e()*vec2.e()  +ten.xx()*vec1.px()*vec2.px()
        +ten.yy()*vec1.py()*vec2.py()+ten.zz()*vec1.pz()*vec2.pz())
    -tentx*(vec1.e()*vec2.px() +vec1.px()*vec2.e())
    -tenty*(vec1.e()*vec2.py() +vec1.py()*vec2.e())
    -tentz*(vec1.e()*vec2.pz() +vec1.pz()*vec2.e())
    +tenxy*(vec1.px()*vec2.py()+vec1.py()*vec2.px())
    +tenxz*(vec1.px()*vec2.pz()+vec1.pz()*vec2.px())
    +tenyz*(vec1.py()*vec2.pz()+vec1.pz()*vec2.py());
  // trace of the tensor
  Complex 
    trace = ten.tt()-ten.xx()-ten.yy()-ten.zz();
  // evaluate the vertex
  Complex vertex = -0.5*ii*norm*(trace*(dotk1v2*dotk2v1-dotv1v2*mdot)
    				 +mdot*tenv1v2-dotk2v1*tenk1v2
    				     -dotk1v2*tenk2v1+dotv1v2*tenk1k2);
  return vertex;
}
// evaluate an off-shell vector
VectorWaveFunction VVTVertex::evaluate(Energy2 q2, int iopt, tcPDPtr out,
    				   const VectorWaveFunction & vec,
    				   const TensorWaveFunction & ten)
{
  // pointers to the particle data objects
  tcPDPtr Pvec = vec.getParticle();
  tcPDPtr Pten = ten.getParticle();
  // evaluate the couplings
  setCoupling(q2,Pvec,out,Pten);
  // outgoing momentum
  Lorentz5Momentum pout= Lorentz5Momentum(ten.px()+vec.px(),ten.py()+vec.py(),
    				      ten.pz()+vec.pz(),ten.e() +vec.e());  
  // normalisation factor
  Energy2 mass  = out->mass();
  Energy2 mass2 = mass*mass;
  Energy2 p2=  pout.m2();
  Complex ii(0.,1.);
  Complex fact=-0.5*getNorm()*propagator(iopt,p2,out);
  // dot product of wavefunctions and momenta
  Complex dotk1k2 =+vec.e()*pout.e()  -vec.px()*pout.px()
    -vec.py()*pout.py()-vec.pz()*pout.pz();
  Complex dotk2v1 = vec.t()*pout.e() -vec.x()*pout.px()
    -vec.y()*pout.py()-vec.z()*pout.pz();
  // mass-k1.k2
  Complex mdot = -dotk1k2;
  if(mass!=0.){mdot=mdot+mass2;}
  // components of the tensor
  Complex tentx = ten.tx()+ten.xt();
  Complex tenty = ten.ty()+ten.yt();
  Complex tentz = ten.tz()+ten.zt();
  Complex tenxy = ten.xy()+ten.yx();
  Complex tenxz = ten.xz()+ten.zx();
  Complex tenyz = ten.yz()+ten.zy();
  // dot product of momenta and polarization vectors with the tensor
  Complex tenk2v1 =
    2.*(+ten.tt()*vec.t()*pout.e() +ten.xx()*vec.x()*pout.px()
        +ten.yy()*vec.y()*pout.py()+ten.zz()*vec.z()*pout.pz())
    -tentx*(vec.t()*pout.px()+vec.x()*pout.e())
    -tenty*(vec.t()*pout.py()+vec.y()*pout.e())
    -tentz*(vec.t()*pout.pz()+vec.z()*pout.e())
    +tenxy*(vec.x()*pout.py()+vec.y()*pout.px())
    +tenxz*(vec.x()*pout.pz()+vec.z()*pout.px())
    +tenyz*(vec.y()*pout.pz()+vec.z()*pout.py());
  Complex tenk1k2 =
    2.*(+ten.tt()*vec.e()*pout.e()  +ten.xx()*vec.px()*pout.px()
        +ten.yy()*vec.py()*pout.py()+ten.zz()*vec.pz()*pout.pz())
    -tentx*(vec.e()*pout.px() +vec.px()*pout.e())
    -tenty*(vec.e()*pout.py() +vec.py()*pout.e())
    -tentz*(vec.e()*pout.pz() +vec.pz()*pout.e())
    +tenxy*(vec.px()*pout.py()+vec.py()*pout.px())
    +tenxz*(vec.px()*pout.pz()+vec.pz()*pout.px())
    +tenyz*(vec.py()*pout.pz()+vec.pz()*pout.py());
  // trace of the tensor
  Complex trace = ten.tt()-ten.xx()-ten.yy()-ten.zz();
  // compute the vector
  Complex vec1[3];
  vec1[0] = mdot*(+   tentx*   vec.t()-2.*ten.xx()*vec.x()
    	      -   tenxy*   vec.y()-    tenxz*  vec.z()-trace*vec.x())
    +(tenk2v1-trace*dotk2v1)*vec.px()-tenk1k2*vec.x()
    +dotk2v1*(+   tentx*   vec.e() -2.*ten.xx()*vec.px()
    	  -   tenxy*   vec.py()-    tenxz*  vec.pz());
  vec1[1] = mdot*(+   tenty*   vec.t()-    tenxy*  vec.x()
    	      -2.*ten.yy()*vec.y()-    tenyz*  vec.z()-trace*vec.y())
    +(tenk2v1-trace*dotk2v1)*vec.py()-tenk1k2*vec.y()
    +dotk2v1*(+    tenty*  vec.e() -    tenxy*  vec.px()
    	  -2.*ten.yy()*vec.py()-    tenyz*  vec.pz());
  vec1[2] =  mdot*(+   tentz*   vec.t()-    tenxz*  vec.x()
    	       -   tenyz*   vec.y()-2.*ten.zz()*vec.z()-trace*vec.z())
    +(tenk2v1-trace*dotk2v1)*vec.pz()-tenk1k2*vec.z()
    +dotk2v1*(+   tentz*   vec.e() -    tenxz  *vec.px()
    	  -   tenyz*   vec.py()-2.*ten.zz()*vec.pz());
  vec1[3] = mdot*(+2.*ten.tt()*vec.t()-    tentx*  vec.x()
    	      -   tenty*   vec.y()-    tentz*  vec.z()-trace*vec.t())
    +(tenk2v1-trace*dotk2v1)*vec.e() -tenk1k2*vec.t()
    +dotk2v1*(+2.*ten.tt()*vec.e() -    tentx*  vec.px()
    	  -   tenty*   vec.py()-    tentz*  vec.pz());
  // now add the piece for massive bosons
  if(mass!=0.)
    {
      Complex dot = tenk2v1-dotk1k2*trace;
      vec1[0]-=dot*pout.px();
      vec1[1]-=dot*pout.py();
      vec1[2]-=dot*pout.pz();
      vec1[3]-=dot*pout.e();
    }
  // return the VectorWaveFunction
  for(int ix=0;ix<4;++ix){vec1[ix]=vec1[ix]*fact;}
  return VectorWaveFunction(pout,out,vec1[0],vec1[1],vec1[2],vec1[3]);
}
// offs-shell tensor
TensorWaveFunction VVTVertex::evaluate(Energy2 q2, int iopt,tcPDPtr out,
    				   const VectorWaveFunction & vec1,
    				   const VectorWaveFunction & vec2)
{
  // pointers to the particle data
  tcPDPtr Pvec1=vec1.getParticle();
  tcPDPtr Pvec2=vec2.getParticle();
  // coupling
  setCoupling(q2,Pvec1,Pvec2,out);
  Complex ii(0.,1.);
  // momenta of the outgoing tensor
  // outgoing momentum
  Lorentz5Momentum pout= Lorentz5Momentum(vec1.px()+vec2.px(),vec1.py()+vec2.py(),
    				      vec1.pz()+vec2.pz(),vec1.e() +vec2.e());
  // overall normalisation
  Energy tmass =out->mass();
  Energy2 tmass2=tmass*tmass;   
  Energy vmass   = Pvec1->mass();
  Energy2 vmass2 = vmass*vmass;
  Energy2 p2=pout.m2();
  Complex fact=0.5*getNorm()*propagator(iopt,p2,out);
  // dot products we need to construct the tensor
  Complex dotk1k2=
    +vec1.e() *vec2.e() -vec1.px()*vec2.px()
    -vec1.py()*vec2.py()-vec1.pz()*vec2.pz();
  Complex dotv1k2=
    +vec1.t()*vec2.e() -vec1.x()*vec2.px()
    -vec1.y()*vec2.py()-vec1.z()*vec2.pz();
  Complex dotv1v2=
    +vec1.t()*vec2.t()-vec1.x()*vec2.x()
    -vec1.y()*vec2.y()-vec1.z()*vec2.z();
  Complex dotk1v2=
    +vec1.e() *vec2.t()-vec1.px()*vec2.x()
    -vec1.py()*vec2.y()-vec1.pz()*vec2.z();
  Complex dotkk1=
    +vec1.e() *pout.e() -vec1.px()*pout.px()
    -vec1.py()*pout.py()-vec1.pz()*pout.pz();
  Complex dotkk2=
    +vec2.e() *pout.e() -vec2.px()*pout.px()
    -vec2.py()*pout.py()-vec2.pz()*pout.pz();
  Complex dotkv1=
    +vec1.t()*pout.e() -vec1.x()*pout.px()
    -vec1.y()*pout.py()-vec1.z()*pout.pz();
  Complex dotkv2=
    +vec2.t() *pout.e()-vec2.x()*pout.px()
    -vec2.y()*pout.py()-vec2.z()*pout.pz();
  // dot product ma^2+k1.k2
  Complex mdot=vmass2+dotk1k2;
  // vectors to help construct the tensor
  Complex vecv1[4],vecv2[4],veck1[4],veck2[4];
  for(int ix=0;ix<4;++ix)
    {
      vecv1[ix]=vec1(ix)-pout[ix]*dotkv1/tmass2;
      vecv2[ix]=vec2(ix)-pout[ix]*dotkv2/tmass2;
      veck1[ix]=vec1[ix]-pout[ix]*dotkk1/tmass2;
      veck2[ix]=vec2[ix]-pout[ix]*dotkk2/tmass2;
    }
  // coefficient of g(nu,mu)-k^muk^nu/m^2
  Complex coeff1 =
    +4./3.*mdot*(-2.*dotv1v2+(dotkv1*dotkv2+p2*dotv1v2)/tmass2)
    +4./3.*(+dotv1k2*(dotk1v2-dotkk1*dotkv2/tmass2)
    	+dotk1v2*(dotv1k2-dotkk2*dotkv1/tmass2)
    	-dotv1v2*(dotk1k2-dotkk1*dotkk2/tmass2)
    	+(1.-p2/tmass2)*dotk1v2*dotv1k2);
  // coefficient of g(nu,mu)
  Complex coeff2 =
    +2.*mdot*(1.-p2/tmass2)*dotv1v2
    -2.*(1.-p2/tmass2)*dotk1v2*dotv1k2; 
  // construct the tensor
  Complex ten[4][4];
  for(int ix=0;ix<4;++ix)
    {
      for(int iy=0;iy<4;++iy)
        {
          ten[ix][iy] = 
    	+2.*mdot*   (vecv1[ix]*vecv2[iy]+vecv1[iy]*vecv2[ix])
    	-2.*dotv1k2*(veck1[ix]*vecv2[iy]+veck1[iy]*vecv2[ix])
    	-2.*dotk1v2*(veck2[ix]*vecv1[iy]+veck2[iy]*vecv1[ix])
    	+2.*dotv1v2*(veck1[ix]*veck2[iy]+veck1[iy]*veck2[ix])
    	-coeff1/tmass2*pout[ix]*pout[iy];
        }
    }
  // add the g(mu,nu) term
  ten[0][0] = ten[0][0]-(coeff1+coeff2);
  ten[1][1] = ten[1][1]-(coeff1+coeff2);
  ten[2][2] = ten[2][2]-(coeff1+coeff2);
  ten[3][3] = ten[3][3]+(coeff1+coeff2);
  // overall coefficent
  for(int ix=0;ix<4;++ix)
    {
      for(int iy=0;iy<4;++iy)
        {
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

