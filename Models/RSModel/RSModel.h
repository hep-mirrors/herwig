// -*- C++ -*-
#ifndef HERWIG_RSModel_H
#define HERWIG_RSModel_H
// This is the declaration of the RSModel class.

#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "Herwig++/Helicity/Vertex/Tensor/FFTVertex.h"
#include "Herwig++/Helicity/Vertex/Tensor/VVTVertex.h"
#include "Herwig++/Helicity/Vertex/Tensor/SSTVertex.h"
#include "Herwig++/Helicity/Vertex/Tensor/FFVTVertex.h"
#include "Herwig++/Helicity/Vertex/Tensor/VVVTVertex.h"

namespace Herwig {
using namespace ThePEG;

typedef Ptr<Herwig::Helicity::FFTVertex>::pointer FFTPtr;
typedef Ptr<Herwig::Helicity::FFTVertex>::transient_pointer tFFTPtr;
typedef Ptr<Herwig::Helicity::VVTVertex>::pointer VVTPtr;
typedef Ptr<Herwig::Helicity::VVTVertex>::transient_pointer tVVTPtr;
typedef Ptr<Herwig::Helicity::SSTVertex>::pointer SSTPtr;
typedef Ptr<Herwig::Helicity::SSTVertex>::transient_pointer tSSTPtr;
typedef Ptr<Herwig::Helicity::FFVTVertex>::pointer FFVTPtr;
typedef Ptr<Herwig::Helicity::FFVTVertex>::transient_pointer tFFVTPtr;
typedef Ptr<Herwig::Helicity::VVVTVertex>::pointer VVVTPtr;
typedef Ptr<Herwig::Helicity::VVVTVertex>::transient_pointer tVVVTPtr;

/** \ingroup Models
 *
 *  This is the class to be used instead of the Standard Model class for
 *  the Randell Sundrum model.
 *
 *
 * @see StandardModel
 * @see StandardModelBase
 * 
 */
class RSModel: public StandardModel {

  
public:
  
  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor
   */
  inline RSModel();

  /**
   * Copy constructor
   */
  inline RSModel(const RSModel &);

  /**
   * Destructor
   */
  virtual ~RSModel();
  //@}
  
public:
  
  /**
   * Return the gravition coupling
   */
  inline Energy lambda_pi() const;


  /** @name Vertices */
  //@{
  /**
   * Pointer to the object handling the \f$G\to f\bar{f}\f$ vertex.
   */
  inline tFFTPtr   vertexFFGR() const;

  /**
   * Pointer to the object handling the \f$G\to VV\f$ vertex.
   */
  inline tVVTPtr   vertexVVGR() const;

  /**
   * Pointer to the object handling the \f$G\to SS\f$ vertex.
   */
  inline tSSTPtr   vertexSSGR() const;

  /**
   * Pointer to the object handling the \f$G\to f\bar{f}V\f$ vertex.
   */
  inline tFFVTPtr  vertexFFVGR() const;

  /**
   * Pointer to the object handling the \f$G\to VVV\f\$ vertex.
   */
  inline tVVVTPtr  vertexVVVGR() const;
  
public:
  
  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(PersistentOStream & os) const;

  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void persistentInput(PersistentIStream & is, int version);
  //@}

  /**
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();
  
protected:
  
  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}

protected:
  
  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Check sanity of the object during the setup phase.
   */
  inline virtual void doupdate() throw(UpdateException);

  /**
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  inline virtual void doinit() throw(InitException);

  /**
   * Initialize this object to the begining of the run phase.
   */
  inline virtual void doinitrun();

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  inline virtual void dofinish();

  /**
   * Rebind pointer to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   * @throws RebindException if no cloned object was found for a given pointer.
   */
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);

  /**
   * Return a vector of all pointers to Interfaced objects used in
   * this object.
   * @return a vector of pointers.
   */
  inline virtual IVector getReferences();
  //@}

private:
  
  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<RSModel> initRSModel;
  
    /**
     * Private and non-existent assignment operator.
     */
  RSModel & operator=(const RSModel &);
  
  /**
   * Coupling of the graviton
   */
  Energy _theLambda_pi;

  /**
   * Pointer to the object handling the \f$G\to f\bar{f}\f$ vertex.
   */
  FFTPtr  _theFFGRVertex;

  /**
   * Pointer to the object handling the \f$G\to VV\f$ vertex.
   */
  VVTPtr  _theVVGRVertex;

  /**
   * Pointer to the object handling the \f$G\to SS\f$ vertex.
   */
  SSTPtr  _theSSGRVertex;

  /**
   * Pointer to the object handling the \f$G\to f\bar{f}V\f$ vertex.
   */
  FFVTPtr _theFFVGRVertex;

  /**
   * Pointer to the object handling the \f$G\to VVV\f\$ vertex.
   */
  VVVTPtr _theVVVGRVertex;
  
};
}
#include "RSModel.icc"


namespace ThePEG {
  
  /**
   * The following template specialization informs ThePEG about the
   * base class of RSModel.
   */
  template <>
  struct BaseClassTrait<Herwig::RSModel,1> {
    /** Typedef of the base class of RSModel. */
    typedef Herwig::StandardModel NthBase;
  };
  
  /**
   * The following template specialization informs ThePEG about the
   * name of this class and the shared object where it is defined.
   */
  template <>
  struct ClassTraits<Herwig::RSModel>
    : public ClassTraitsBase<Herwig::RSModel> {
    /** Return the class name.*/
    static string className() { return "Herwig++::RSModel"; }
    /**
     * Return the name of the shared library to be loaded to get
     * access to this class and every other class it uses
     * (except the base class).
     */
    static string library() { return "libHwRSModel.so"; }

  };
  
}


#endif /* HERWIG_RSModel_H */
