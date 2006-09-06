// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VVVVertex class.
//

#include "VVVVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;
    
AbstractNoPIOClassDescription<VVVVertex> VVVVertex::initVVVVertex;
// Definition of the static class description member.
  
void VVVVertex::Init() {
  
  static ClassDocumentation<VVVVertex> documentation
    ("The VVVVertex class implements the helicity amplitude"
     "calculations for the triple gauge boson vertex. Any   "
     "implementation of such a vertex should inherit from in and implement"
     " the virtual setCoupling member to calculate the coupling");
  
}

// evaluate the vertex
Complex VVVVertex::evaluate(Energy2 q2, const VectorWaveFunction & vec1,
			    const VectorWaveFunction & vec2,
			    const VectorWaveFunction & vec3)
{
  tcPDPtr Pvec1 = vec1.getParticle();
  tcPDPtr Pvec2 = vec2.getParticle();
  tcPDPtr Pvec3 = vec3.getParticle();
  // calculate the coupling
  setCoupling(q2,Pvec1,Pvec2,Pvec3);
  Complex norm=getNorm();
  Complex ii(0.,1.);
  Complex alpha1(0.);
  // decide if we need to use special treatment to avoid gauge cancelations
  // first vector
  if(abs(vec1.t())!=0.)
    {
      if(abs(vec1.t())>0.1*max( max(abs(vec1.x()),abs(vec1.y())),abs(vec1.z())))
	{alpha1=vec1.e()/vec1.t();}
    }
  // second vector
  if(abs(vec2.t())!=0.)
    {
      if(abs(vec2.t())>0.1*max( max(abs(vec2.x()),abs(vec2.y())),abs(vec2.z())))
	{alpha1=vec2.e()/vec2.t();}
    }
  // third vector
  if(abs(vec3.t())!=0.)
    {
      if(abs(vec3.t())>0.1*max( max(abs(vec3.x()),abs(vec3.y())),abs(vec3.z())))
	    {alpha1=vec3.e()/vec3.t();}
    }
  // dot products of the polarization vectors
  Complex dot12 = vec1.t()*vec2.t()-vec1.x()*vec2.x()
    -vec1.y()*vec2.y()-vec1.z()*vec2.z();
  Complex dot13 = vec1.t()*vec3.t()-vec1.x()*vec3.x()
    -vec1.y()*vec3.y()-vec1.z()*vec3.z();
  Complex dot23 = vec3.t()*vec2.t()-vec3.x()*vec2.x()
    -vec3.y()*vec2.y()-vec3.z()*vec2.z();
  // dot products of polarization vectors and momentum
  Complex
    dotp13 = 
    (vec1.e() -alpha1*vec1.t())*vec3.t()-(vec1.px()-alpha1*vec1.x())*vec3.x()
    -(vec1.py()-alpha1*vec1.y())*vec3.y()-(vec1.pz()-alpha1*vec1.z())*vec3.z();
  Complex
    dotp23 =
    (vec2.e() -alpha1*vec2.t())*vec3.t()-(vec2.px()-alpha1*vec2.x())*vec3.x()
    -(vec2.py()-alpha1*vec2.y())*vec3.y()-(vec2.pz()-alpha1*vec2.z())*vec3.z();
  Complex
    dotp21 = 
    (vec2.e() -alpha1*vec2.t())*vec1.t()-(vec2.px()-alpha1*vec2.x())*vec1.x()
    -(vec2.py()-alpha1*vec2.y())*vec1.y()-(vec2.pz()-alpha1*vec2.z())*vec1.z();
  Complex
    dotp31 = 
    (vec3.e() -alpha1*vec3.t())*vec1.t()-(vec3.px()-alpha1*vec3.x())*vec1.x()
    -(vec3.py()-alpha1*vec3.y())*vec1.y()-(vec3.pz()-alpha1*vec3.z())*vec1.z();
  Complex
    dotp32 = 
    (vec3.e() -alpha1*vec3.t())*vec2.t()-(vec3.px()-alpha1*vec3.x())*vec2.x()
    -(vec3.py()-alpha1*vec3.y())*vec2.y()-(vec3.pz()-alpha1*vec3.z())*vec2.z();
  Complex
    dotp12 = 
    (vec1.e() -alpha1*vec1.t())*vec2.t()-(vec1.px()-alpha1*vec1.x())*vec2.x()
    -(vec1.py()-alpha1*vec1.y())*vec2.y()-(vec1.pz()-alpha1*vec1.z())*vec2.z();
  // finally calculate the vertex
  return ii*norm*(dot12*(dotp13-dotp23)+dot23*(dotp21-dotp31)+dot13*(dotp32-dotp12));
}
  
// off-shell vector
VectorWaveFunction VVVVertex::evaluate(Energy2 q2,int iopt, tcPDPtr out,
				       const VectorWaveFunction & vec1,
				       const VectorWaveFunction & vec2)
{
  // pointer to particle data objects
  tcPDPtr Pvec1 = vec1.getParticle();
  tcPDPtr Pvec2 = vec2.getParticle();
  // output momenta
  Lorentz5Momentum pout = Lorentz5Momentum(vec1.px()+vec2.px(),vec1.py()+vec2.py(),
					   vec1.pz()+vec2.pz(),vec1.e() +vec2.e() ); 
  // calculate the coupling
  setCoupling(q2,out,Pvec1,Pvec2);
  // prefactor
  Energy2 p2=pout.m2();
  Complex fact=getNorm()*propagator(iopt,p2,out);
  Energy mass = out->mass();
  Energy2 mass2=mass*mass;
  // compute the polarization vector
  Complex vect[4];
  // dot products we need
  Complex
    dot12 = vec1.t()*vec2.t()-vec1.x()*vec2.x()-vec1.y()*vec2.y()-vec1.z()*vec2.z();
  Complex
	dota = ( pout.e()+vec2.e() )*vec1.t()-(pout.px()+vec2.px())*vec1.x()
    -(pout.py()+vec2.py())*vec1.y()-(pout.pz()+vec2.pz())*vec1.z();
  Complex 
    dotb =  ( pout.e()+vec1.e() )*vec2.t()-(pout.px()+vec1.px())*vec2.x()
    -(pout.py()+vec1.py())*vec2.y()-(pout.pz()+vec1.pz())*vec2.z();
  // construct the vector
  vect[0] = dot12*(vec1.px()-vec2.px())-dotb*vec1.x()+dota*vec2.x();
  vect[1] = dot12*(vec1.py()-vec2.py())-dotb*vec1.y()+dota*vec2.y();
  vect[2] = dot12*(vec1.pz()-vec2.pz())-dotb*vec1.z()+dota*vec2.z();
  vect[3] = dot12*(vec1.e() -vec2.e() )-dotb*vec1.t()+dota*vec2.t();
  // scalar piece for massive case
  if(mass!=0.)
    {
      Complex dot=
	(vect[3]*pout.e()
	 -vect[0]*pout.x()-vect[1]*pout.y()-vect[2]*pout.z())/mass2;
      vect[0] = vect[0] -dot*pout.px();
      vect[1] = vect[1] -dot*pout.py();
      vect[2] = vect[2] -dot*pout.pz();
      vect[3] = vect[3] -dot*pout.e();       
    }
  // multiply by prefactor
  for(int ix=0;ix<4;++ix){vect[ix]=fact*vect[ix];}
  return VectorWaveFunction(pout,out,vect[0],vect[1],vect[2],vect[3]);
}

}
}
