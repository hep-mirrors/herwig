// -*- C++ -*-
#ifndef HERWIG_RSModelSSGRVertex_H
#define HERWIG_RSModelSSGRVertex_H
//
// This is the declaration of the RSModelSSGRVertex class.

#include "Herwig++/Helicity/Vertex/Tensor/SSTVertex.h"
#include "Herwig++/Models/RSModel/RSModel.h"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;
    
/** \ingroup Helicity
 * 
 *  The RSModelSSGRVertex class is thew implementation of the graviton
 *  coupling to the Higgs in the RSModel. It inherits from the SSTVertex 
 *  and implements the setCoupling member
 *
 *  @see SSTVertex
 *  @see VertexBase
 */
class RSModelSSGRVertex: public SSTVertex {
  
public:
  
  /**
   * Standard ctors and dtor.
   */
  inline RSModelSSGRVertex();
  inline RSModelSSGRVertex(const RSModelSSGRVertex &);
  virtual ~RSModelSSGRVertex();
  
public:
  
  /**
   * Standard functions for writing and reading from persistent streams.
   */
  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);
  
  /**
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();
  
  /**
   * Calculate the couplings.
   */
  void setCoupling(Energy2,tcPDPtr,tcPDPtr,tcPDPtr);

protected:
  
  /**
   * Standard clone methods.
   */
  inline virtual IBPtr clone() const;
  inline virtual IBPtr fullclone() const;
  
protected:
  
  /**
   * Standard Interfaced virtual functions.
   */
  inline virtual void doupdate() throw(UpdateException);
  inline virtual void doinit() throw(InitException);
  inline virtual void doinitrun();
  inline virtual void dofinish();
  
  /**
   * Change all pointers to Interfaced objects to corresponding clones.
   */
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);
  
  /**
   * Return pointers to all Interfaced objects refered to by this.
   */
  inline virtual IVector getReferences();
  
private:
  
  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<RSModelSSGRVertex> initRSModelSSGRVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  RSModelSSGRVertex & operator=(const RSModelSSGRVertex &);
  
  /**
   * Pointer to the Mode.
   */
  SMPtr _theModel;

  /**
   * Coupling.
   */
  double _theKappa;
};

}
}    

#include "RSModelSSGRVertex.icc"

namespace ThePEG {
  
  /**
   * The following template specialization informs ThePEG about the
   * base class of RSModelSSGRVertex.
   */
  template <>
  struct BaseClassTrait<Herwig::Helicity::RSModelSSGRVertex,1> {
    typedef Herwig::Helicity::SSTVertex NthBase;
  };
  
  /**
   * The following template specialization informs ThePEG about the
   * name of this class and the shared object where it is defined.
   */
  template <>
  struct ClassTraits<Herwig::Helicity::RSModelSSGRVertex>
    : public ClassTraitsBase<Herwig::Helicity::RSModelSSGRVertex> {

    /**
     * Return the class name.
     */
    static string className() { return "/Herwig++/Helicity/RSModelSSGRVertex"; }

    /** 
     * Return the name of the shared library to be loaded to get
     * access to this class and every other class it uses
     * (except the base class).
     */
    static string library() { return "libHwRSVertex.so"; }

  };
  
}


#endif /* HERWIG_RSModelSSGRVertex_H */
