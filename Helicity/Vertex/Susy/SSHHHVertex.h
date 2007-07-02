// -*- C++ -*-
#ifndef HERWIG_SSHHHVertex_H
#define HERWIG_SSHHHVertex_H
//
// This is the declaration of the SSHHHVertex class.
//

#include "Herwig++/Helicity/Vertex/Scalar/SSSVertex.h"
#include "Herwig++/Models/Susy/MSSM.h"
#include "SSHHHVertex.fh"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;

/**
 * This is the higgs coupling with other higgs bosons in the MSSM.
 *
 * @see \ref SSHHHVertexInterfaces "The interfaces"
 * defined for SSHHHVertex.
 */
class SSHHHVertex: public SSSVertex {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  SSHHHVertex();

  /**
   * The destructor.
   */
  virtual ~SSHHHVertex();
  //@}

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
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

  /**
   * Calculate the coupling for the vertex
   * @param q2 The scale to at which evaluate the coupling.
   * @param particle1 The first particle in the vertex.
   * @param particle2 The second particle in the vertex.
   * @param particle3 The third particle in the vertex.
   */
  virtual void setCoupling(Energy2 q2, tcPDPtr particle1, tcPDPtr particle2,
			   tcPDPtr particle3);

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
  //@}

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit() throw(InitException);
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<SSHHHVertex> initSSHHHVertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SSHHHVertex & operator=(const SSHHHVertex &);

private:

  /**
   * A pointer to the MSSM object 
   */
  tMSSMPtr theMSSM;
  
  /**
   * The mass of the \f$W\f$.
   */
  Energy theMw;

  /**
   * The factor \f$ \frac{m_Z}{\sin2\theta_W} \f$
   */
  Energy theZfact;

  /**
   * The value of \f$\sin\theta_W\f$
   */
  double theSw;

  /**
   * The value of \f$ \sin(\beta + \alpha) \f$.
   */
  double theSbpa;

  /**
   * The value of \f$ \cos(\beta + \alpha) \f$.
   */
  double theCbpa;

  /**
   * The value of \f$ \sin(\beta - \alpha) \f$.
   */
  double theSbma;

  /**
   * The value of \f$ \cos(\beta - \alpha) \f$.
   */
  double theCbma;
  
  /**
   * The value of \f$ \sin 2\alpha \f$.
   */
  double theS2a;

  /**
   * The value of \f$ \cos 2\alpha \f$.
   */
  double theC2a;

  /**
   * The value of \f$ \sin 2\beta \f$.
   */
  double theS2b;

  /**
   * The value of \f$ \cos 2\beta \f$.
   */
  double theC2b;
  
  /**
   * The value of \f$ \sqrt(4.*\pi*\alpha)\f$ when it was last evaluated.
   */
  double theElast;

  /**
   * The scale at which the coupling was last evaluated.  
   */
  Energy2 theq2last;
  
};

}
}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of SSHHHVertex. */
template <>
struct BaseClassTrait<Herwig::Helicity::SSHHHVertex,1> {
  /** Typedef of the first base class of SSHHHVertex. */
  typedef Herwig::Helicity::SSSVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SSHHHVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::Helicity::SSHHHVertex>
  : public ClassTraitsBase<Herwig::Helicity::SSHHHVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::SSHHHVertex"; }
  /**
   * The name of a file containing the dynamic library where the class
   * SSHHHVertex is implemented. It may also include several, space-separated,
   * libraries if the class SSHHHVertex depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwSusy.so HwSusyVertex.so"; }
};

/** @endcond */

}

#include "SSHHHVertex.icc"

#endif /* HERWIG_SSHHHVertex_H */
