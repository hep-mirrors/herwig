// -*- C++ -*-
#ifndef HERWIG_VVVVVertex_H
#define HERWIG_VVVVVertex_H
//
// This is the declaration of the <!id>VVVVVertex<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This is the implementation of the four vector vertex. It is based on the VertexBase
// class for the storage of particles which are allowed to interact at the vertex.
// Classes implementation a specific vertex should inherit from this one and implement
// the virtual setCoupling member.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="VertexBase.html">VertexBase.h</a>.
// 

#include "Herwig++/Helicity/Vertex/VertexBase.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"

namespace Herwig {
namespace Helicity{
using namespace ThePEG;

class VVVVVertex: public VertexBase {
  
public:
  
  inline VVVVVertex();
  inline VVVVVertex(const VVVVVertex &);
  virtual ~VVVVVertex();
  // Standard ctors and dtor.
  
public:
  
  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);
      // Standard functions for writing and reading from persistent streams.
  
  static void Init();
  // Standard Init function used to initialize the interfaces.
  
public:
  
  // evaluate the vertex
  Complex evaluate(Energy2, int,
		   const VectorWaveFunction &, const VectorWaveFunction &,
		   const VectorWaveFunction &, const VectorWaveFunction &);
  // calculate the couplings etc
  virtual void setCoupling(Energy2,tcPDPtr,tcPDPtr,tcPDPtr,tcPDPtr);
  
protected:
  
  inline virtual void doupdate() throw(UpdateException);
  inline virtual void doinit() throw(InitException);
  inline virtual void doinitrun();
  inline virtual void dofinish();
  // Standard Interfaced virtual functions.
  
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);
  // Change all pointers to Interfaced objects to corresponding clones.
  
  inline virtual IVector getReferences();
  // Return pointers to all Interfaced objects refered to by this.
  
protected:
  
  // set the order of the particles
  inline void setOrder(int,int,int,int);
  // set the type of the vertex (QCD=1 or electroweak=2)
  inline void setType(int);
  // set the intermediate particles if including s/u/t channel terms
  inline void setIntermediate(tcPDPtr,tcPDPtr,Complex,Complex);
  
private:
  
  static AbstractClassDescription<VVVVVertex> initVVVVVertex;
  // Describe a concrete class with persistent data.
  
  VVVVVertex & operator=(const VVVVVertex &);
  // Private and non-existent assignment operator.
  
private:
  // type of vertex 1=QCD 2=EW
  int _itype;
  // order of the particles
  int _iorder[4];
  tcPDPtr _inter[2];
  Complex _coup[2];
};
}
}
#include "VVVVVertex.icc"

// CLASSDOC OFF

namespace ThePEG {
  
  // The following template specialization informs ThePEG about the
  // base class of VVVVVertex.
  template <>
  struct BaseClassTrait<Herwig::Helicity::VVVVVertex,1> {
    typedef Herwig::Helicity::VertexBase NthBase;
  };
  
  // The following template specialization informs ThePEG about the
  // name of this class and the shared object where it is defined.
  template <>
  struct ClassTraits<Herwig::Helicity::VVVVVertex>
    : public ClassTraitsBase<Herwig::Helicity::VVVVVertex> {
    static string className() { return "/Herwig++/Helicity/VVVVVertex"; }
    // Return the class name.
    static string library() { return "libHwVVertex.so"; }
    // Return the name of the shared library to be loaded to get
    // access to this class and every other class it uses
    // (except the base class).
  };
  
}

#endif /* HERWIG_VVVVVertex_H */
