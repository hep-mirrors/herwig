// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VVVTVertex class.
//

#include "VVVTVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;
    
VVVTVertex::~VVVTVertex() {}

void VVVTVertex::persistentOutput(PersistentOStream & os) const { }

void VVVTVertex::persistentInput(PersistentIStream & is, int) { }

AbstractClassDescription<VVVTVertex> VVVTVertex::initVVVTVertex;
// Definition of the static class description member.

void VVVTVertex::Init() {
  
  static ClassDocumentation<VVVTVertex> documentation
    ("The\\classname{VVVTVertex} class is the implementation f the"
     " helicity amplitude calculation of the vector-vector-vector-tensor"
     " vertex. All such vertices should inherit from it.");
}
// coupling
void VVVTVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b,tcPDPtr c, tcPDPtr d){;}
// function to evaluate the vertex
Complex VVVTVertex::evaluate(Energy q2, const VectorWaveFunction & vec1,
    				 const VectorWaveFunction & vec2,
    				 const VectorWaveFunction & vec3,
    				 const TensorWaveFunction & ten)
{
  // pointers to the particles
  tcPDPtr Pvec1 = vec1.getParticle();
  tcPDPtr Pvec2 = vec2.getParticle();
  tcPDPtr Pvec3 = vec3.getParticle();
  tcPDPtr Pten  =  ten.getParticle();
  // set the couplings
  setCoupling(q2,Pvec1,Pvec2,Pvec3,Pten);
  Complex norm=getNorm();
  Complex ii(0.,1.);
  // dot products of the wavefunctions
  Complex dotv1v2 = 
    +vec1.t()*vec2.t()-vec1.x()*vec2.x()
    -vec1.y()*vec2.y()-vec1.z()*vec2.z();
  Complex dotv1v3 = 
    +vec1.t()*vec3.t()-vec1.x()*vec3.x()
    -vec1.y()*vec3.y()-vec1.z()*vec3.z();
  Complex dotv2v3 = 
    +vec2.t()*vec3.t()-vec2.x()*vec3.x()
    -vec2.y()*vec3.y()-vec2.z()*vec3.z();
  // dot product of wavefunctions and momenta
  Complex dotv1k23 = 
    +vec1.t()*(vec2.e() -vec3.e() )-vec1.x()*(vec2.px()-vec3.px())
    -vec1.y()*(vec2.py()-vec3.py())-vec1.z()*(vec2.pz()-vec3.pz());
  Complex dotv2k31 = 
    +vec2.t()*(vec3.e() -vec1.e() )-vec2.x()*(vec3.px()-vec1.px())
    -vec2.y()*(vec3.py()-vec1.py())-vec2.z()*(vec3.pz()-vec1.pz());
  Complex dotv3k12 = 
    +vec3.t()*(vec1.e() -vec2.e() )-vec3.x()*(vec1.px()-vec2.px())
    -vec3.y()*(vec1.py()-vec2.py())-vec3.z()*(vec1.pz()-vec2.pz());
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
  Complex tenv1v3 =
    2.*(+ten.tt()*vec1.t()*vec3.t()+ten.xx()*vec1.x()*vec3.x()
        +ten.yy()*vec1.y()*vec3.y()+ten.zz()*vec1.z()*vec3.z())
    -tentx*(vec1.t()*vec3.x()+vec1.x()*vec3.t())
    -tenty*(vec1.t()*vec3.y()+vec1.y()*vec3.t())
    -tentz*(vec1.t()*vec3.z()+vec1.z()*vec3.t())
    +tenxy*(vec1.x()*vec3.y()+vec1.y()*vec3.x())
    +tenxz*(vec1.x()*vec3.z()+vec1.z()*vec3.x())
    +tenyz*(vec1.y()*vec3.z()+vec1.z()*vec3.y());
  Complex tenv2v3 =
    2.*(+ten.tt()*vec2.t()*vec3.t()+ten.xx()*vec2.x()*vec3.x()
        +ten.yy()*vec2.y()*vec3.y()+ten.zz()*vec2.z()*vec3.z())
    -tentx*(vec2.t()*vec3.x()+vec2.x()*vec3.t())
    -tenty*(vec2.t()*vec3.y()+vec2.y()*vec3.t())
    -tentz*(vec2.t()*vec3.z()+vec2.z()*vec3.t())
    +tenxy*(vec2.x()*vec3.y()+vec2.y()*vec3.x())
    +tenxz*(vec2.x()*vec3.z()+vec2.z()*vec3.x())
    +tenyz*(vec2.y()*vec3.z()+vec2.z()*vec3.y());
  Complex tenv1k23 =
    2.*(+ten.tt()*vec1.t()*(vec2.e() -vec3.e() )
        +ten.xx()*vec1.x()*(vec2.px()-vec3.px())
        +ten.yy()*vec1.y()*(vec2.py()-vec3.py())
        +ten.zz()*vec1.z()*(vec2.pz()-vec3.pz()))
    -tentx*(vec1.t()*(vec2.px()-vec3.px())+vec1.x()*(vec2.e() -vec3.e() ))
    -tenty*(vec1.t()*(vec2.py()-vec3.py())+vec1.y()*(vec2.e() -vec3.e() ))
    -tentz*(vec1.t()*(vec2.pz()-vec3.pz())+vec1.z()*(vec2.e() -vec3.e() ))
    +tenxy*(vec1.x()*(vec2.py()-vec3.py())+vec1.y()*(vec2.px()-vec3.px()))
    +tenxz*(vec1.x()*(vec2.pz()-vec3.pz())+vec1.z()*(vec2.px()-vec3.px()))
    +tenyz*(vec1.y()*(vec2.pz()-vec3.pz())+vec1.z()*(vec2.py()-vec3.py()));
  Complex tenv2k31 =
    2.*(+ten.tt()*vec2.t()*(vec3.e() -vec1.e() )
        +ten.xx()*vec2.x()*(vec3.px()-vec1.px())
        +ten.yy()*vec2.y()*(vec3.py()-vec1.py())
        +ten.zz()*vec2.z()*(vec3.pz()-vec1.pz()))
    -tentx*(vec2.t()*(vec3.px()-vec1.px())+vec2.x()*(vec3.e() -vec1.e() ))
    -tenty*(vec2.t()*(vec3.py()-vec1.py())+vec2.y()*(vec3.e() -vec1.e() ))
    -tentz*(vec2.t()*(vec3.pz()-vec1.pz())+vec2.z()*(vec3.e() -vec1.e() ))
    +tenxy*(vec2.x()*(vec3.py()-vec1.py())+vec2.y()*(vec3.px()-vec1.px()))
    +tenxz*(vec2.x()*(vec3.pz()-vec1.pz())+vec2.z()*(vec3.px()-vec1.px()))
    +tenyz*(vec2.y()*(vec3.pz()-vec1.pz())+vec2.z()*(vec3.py()-vec1.py()));
  Complex tenv3k12 =
    2.*(+ten.tt()*vec3.t()*(vec1.e() -vec2.e() )
        +ten.xx()*vec3.x()*(vec1.px()-vec2.px())
        +ten.yy()*vec3.y()*(vec1.py()-vec2.py())
        +ten.zz()*vec3.z()*(vec1.pz()-vec2.pz()))
    -tentx*(vec3.t()*(vec1.px()-vec2.px())+vec3.x()*(vec1.e() -vec2.e() ))
    -tenty*(vec3.t()*(vec1.py()-vec2.py())+vec3.y()*(vec1.e() -vec2.e() ))
    -tentz*(vec3.t()*(vec1.pz()-vec2.pz())+vec3.z()*(vec1.e() -vec2.e() ))
    +tenxy*(vec3.x()*(vec1.py()-vec2.py())+vec3.y()*(vec1.px()-vec2.px()))
    +tenxz*(vec3.x()*(vec1.pz()-vec2.pz())+vec3.z()*(vec1.px()-vec2.px()))
    +tenyz*(vec3.y()*(vec1.pz()-vec2.pz())+vec3.z()*(vec1.py()-vec2.py()));
  // trace of the tensor
  Complex trace = ten.tt()-ten.xx()-ten.yy()-ten.zz();
  // compute the vertex
  Complex 
    vertex= -0.5*ii*norm*(
    		      +dotv3k12*(tenv1v2-trace*dotv1v2)
    		      +dotv2k31*(tenv1v3-trace*dotv1v3)
    		      +dotv1k23*(tenv2v3-trace*dotv2v3)
    		      +dotv2v3*tenv1k23+dotv1v3*tenv2k31+dotv1v2*tenv3k12);
  // return the answer
  return vertex;
}

}
}

