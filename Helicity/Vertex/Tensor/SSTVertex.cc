// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSTVertex class.
//

#include "SSTVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;
    
SSTVertex::~SSTVertex() {}

void SSTVertex::persistentOutput(PersistentOStream & os) const { }

void SSTVertex::persistentInput(PersistentIStream & is, int) { }

AbstractClassDescription<SSTVertex> SSTVertex::initSSTVertex;
// Definition of the static class description member.

void SSTVertex::Init() {
  
  static ClassDocumentation<SSTVertex> documentation
("The \\classname{SSTVertex} class is the implementation of the "
 "helicity amplitude calculation for the scalar-scalar-tensor vertex"
 ", all such vertices should inherit from it");
  
}
// set coupling member
void SSTVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b,tcPDPtr c){;}
// evaluate the vertex
Complex SSTVertex::evaluate(Energy2 q2, const ScalarWaveFunction & sca1,
    				const ScalarWaveFunction & sca2,
    				const TensorWaveFunction & ten)
{
  // extract the pointers to the particle data objects
  tcPDPtr Psca1 = sca1.getParticle();
  tcPDPtr Psca2 = sca2.getParticle();
  tcPDPtr Pten = ten.getParticle();
  // obtain the coupling
  setCoupling(q2,Psca1,Psca2,Pten);
  Complex norm=getNorm();
  Complex ii(0.,1.);
  Complex vertex(0.);
  // evaluate the trace of the tensor
  Complex trace = ten.tt()-ten.xx()-ten.yy()-ten.zz();
  // dot product of the two momenta
  Energy2 dot=
    + sca1.e()*sca2.e() -sca1.px()*sca2.px()
    -sca1.py()*sca2.py()-sca1.pz()*sca2.pz();
  Energy mass=Psca1->mass();
  // second term
  Complex second = 
    +2.*ten.tt()*sca1.e()*sca2.e()  +2.*ten.xx()*sca1.px()*sca2.px()
    +2.*ten.yy()*sca1.py()*sca2.py()+2.*ten.zz()*sca1.pz()*sca2.pz()
    -(ten.tx()+ten.xt())*(sca1.e()*sca2.px() +sca1.px()*sca2.e())
    -(ten.ty()+ten.yt())*(sca1.e()*sca2.py() +sca1.py()*sca2.e())
    -(ten.tz()+ten.zt())*(sca1.e()*sca2.pz() +sca1.pz()*sca2.e())
    +(ten.xy()+ten.yx())*(sca1.py()*sca2.px()+sca1.px()*sca2.py())
    +(ten.xz()+ten.zx())*(sca1.pz()*sca2.px()+sca1.px()*sca2.pz())
    +(ten.yz()+ten.zy())*(sca1.pz()*sca2.py()+sca1.py()*sca2.pz());
  // return the answer
  return -0.5*ii*norm*(trace*(mass*mass-dot)+second)*sca1.Wave()*sca2.Wave();
}
// off-shell tensor
TensorWaveFunction SSTVertex::evaluate(Energy q2, int iopt, tcPDPtr out,
    				   const ScalarWaveFunction & sca1,
    				   const ScalarWaveFunction & sca2)
{
  // pointers to the particle data objects
  tcPDPtr Psca1 = sca1.getParticle();
  tcPDPtr Psca2 = sca2.getParticle();
  // obtain the coupling
  setCoupling(q2,Psca1,Psca2,out);
  Complex ii(0.,1.);
  // array for the tensor
  Complex ten[4][4];
  // calculate the outgoing momentum
  Lorentz5Momentum pout = Lorentz5Momentum(sca1.px()+sca2.px(),sca1.py()+sca2.py(),
    				       sca1.pz()+sca2.pz(),sca1.e() +sca2.e() );
  // prefactor
  Energy mass = out->mass();
  Energy2 mass2=mass*mass;
  Energy2 p2 = pout.m2();
  Complex fact=0.5*getNorm()*sca1.Wave()*sca2.Wave()*propagator(iopt,p2,out);
  // dot products we need
  Energy2 dot12 = 
    +sca1.e()*sca2.e()  -sca1.px()*sca2.px()
    -sca1.py()*sca2.py()-sca1.pz()*sca2.pz();
  Energy2 dot1 = 
    +sca1.e()*pout.e()  -sca1.px()*pout.px()
    -sca1.py()*pout.py()-sca1.pz()*pout.pz();
  Energy2 dot2 = 
    +pout.e()*sca2.e()  -pout.px()*sca2.px()
    -pout.py()*sca2.py()-pout.pz()*sca2.pz();
  // the vectors that we need for the tensor
  Energy vec1[4],vec2[4];
  Energy2 a,b;
  Energy2 mphi2=(Psca1->mass())*(Psca1->mass());
  // massive case
  if(mass!=0.)
    {
      double norm1=dot1/mass2;
      double norm2=dot2/mass2;
      a = ((mphi2+dot12)*(2.*p2/mass2-5)
           +4.*(dot12-dot1*dot2/mass2))/3.;
      b = -(-(mphi2+dot12)*(2.+p2/mass2)+4*(dot12-dot1*dot2/mass2))/3./mass2;
      for(int ix=0;ix<4;++ix)
        {
          vec1[ix] =sca1[ix]-norm1*pout[ix];
          vec2[ix] =sca2[ix]-norm2*pout[ix];
        }
    }
  // massless case
  else
    {
      a = (-5.*(mphi2+dot12)+4.*dot12)/3.;
      b = 0.;
      for(int ix=0;ix<4;++ix)
        {
          vec1[ix] =sca1[ix];
          vec2[ix] =sca2[ix];
        }
    }
  // calculate the wavefunction
  for(int ix=0;ix<4;++ix)
    {
      for(int iy=0;iy<4;++iy)
        {
          ten[ix][iy]=-2.*(vec1[ix]*vec2[iy]+vec1[ix]*vec2[iy])-b*pout[ix]*pout[iy];
        }
    }
  ten[3][3]=ten[3][3]-a;
  for(int ix=0;ix<3;++ix){ten[ix][ix]=ten[ix][ix]+a;}
  // prefactor
  for(int ix=0;ix<4;++ix){for(int iy=0;iy<4;++iy){ten[ix][iy]=fact*ten[ix][iy];}}
  // return the wavefunction
  return TensorWaveFunction(pout,out,
    			ten[0][0],ten[0][1],ten[0][2],ten[0][3],
    			ten[1][0],ten[1][1],ten[1][2],ten[1][3],
    			ten[2][0],ten[2][1],ten[2][2],ten[2][3],
    			ten[3][0],ten[3][1],ten[3][2],ten[3][3]);
}
// off-shell scalar
ScalarWaveFunction SSTVertex::evaluate(Energy q2,int iopt, tcPDPtr out,
    				   const ScalarWaveFunction & sca,
    				   const TensorWaveFunction & ten)
{
  // pointers to the particle data objects
  tcPDPtr Psca = sca.getParticle();
  tcPDPtr Pten = ten.getParticle();
  // obtain the coupling
  setCoupling(q2,Psca,out,Pten);
  Complex ii(0.,1.);
  // calculate the outgoing momentum
  Lorentz5Momentum pout = Lorentz5Momentum(sca.px()+ten.px(),sca.py()+ten.py(),
    				       sca.pz()+ten.pz(),sca.e() +ten.e() );
  // prefactors
  Energy mass = out->mass();
  Energy2 mass2=mass*mass;
  Energy2 p2 = pout.m2();
  Complex fact=0.5*getNorm()*sca.Wave()*propagator(iopt,p2,out);
  // trace of the tensor
  Complex trace =ten.tt()-ten.xx()-ten.yy()-ten.zz();
  // dot product of the two momenta
  Energy2 dot = 
    + sca.e()*pout.e() -sca.px()*pout.px()
    -sca.py()*pout.py()-sca.pz()*pout.pz();
  // first term
  trace = trace*(mass2-dot);
  // second term
  Complex second = 
    +2.*ten.tt()*sca.e()*pout.e()  +2.*ten.xx()*sca.px()*pout.px()
    +2.*ten.yy()*sca.py()*pout.py()+2.*ten.zz()*sca.pz()*pout.pz()
    -(ten.tx()+ten.xt())*( sca.e()*pout.px()+sca.px()*pout.e())
    -(ten.ty()+ten.yt())*( sca.e()*pout.py()+sca.py()*pout.e())
    -(ten.tz()+ten.zt())*( sca.e()*pout.pz()+sca.pz()*pout.e())
    +(ten.xy()+ten.yx())*(sca.py()*pout.px()+sca.px()*pout.py())
    +(ten.xz()+ten.zx())*(sca.pz()*pout.px()+sca.px()*pout.pz())
    +(ten.yz()+ten.zy())*(sca.py()*pout.pz()+sca.pz()*pout.py());
  // put it all together
  second  = fact*(trace+second);
  return ScalarWaveFunction(pout,out,second);
}

}
}

