// -*- C++ -*-
#ifndef HERWIG_FFVCurrentDecayer_H
#define HERWIG_FFVCurrentDecayer_H
//
// This is the declaration of the FFVCurrentDecayer class.
//

#include "GeneralCurrentDecayer.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "Herwig++/Helicity/Vertex/Vector/FFVVertex.h"
#include "FFVCurrentDecayer.fh"

namespace Herwig {
using namespace ThePEG;
using Helicity::FFVVertexPtr;

/**
 * Here is the documentation of the FFVCurrentDecayer class.
 *
 * @see \ref FFVCurrentDecayerInterfaces "The interfaces"
 * defined for FFVCurrentDecayer.
 */
class FFVCurrentDecayer: public GeneralCurrentDecayer {

public:

  /**
   * The default constructor.
   */
  inline FFVCurrentDecayer();
  
public:

  /** @name Virtual functions required by the Decayer class. */
  //@{
   /**
   * Return the matrix element squared for a given mode and phase-space channel.  
   * @param vertex Output the information on the vertex for spin correlations
   * @param ichan The channel we are calculating the matrix element for.
   * @param part The decaying Particle.
   * @param decay The particles produced in the decay.
   * @return The matrix element squared for the phase-space configuration.
   */
  virtual double me2(bool vertex, const int ichan, const Particle & part,
		     const ParticleVector & decay) const;
  
  /**
   * Function to return partial Width
   * @param inpart Pointer to incoming particle data object
   * @param outa Pointer to first outgoing particle data object
   * @param currout The outgoing particles from the current
   */
  virtual double partialWidth(tPDPtr inpart, tPDPtr outa,
			      vector<tPDPtr> currout);
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
  inline virtual void doinit() throw(InitException);

  /**
   * Rebind pointer to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   * @throws RebindException if no cloned object was found for a given
   * pointer.
   */
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);

  /**
   * Return a vector of all pointers to Interfaced objects used in this
   * object.
   * @return a vector of pointers.
   */
  inline virtual IVector getReferences();
  //@}


private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<FFVCurrentDecayer> initFFVCurrentDecayer;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  FFVCurrentDecayer & operator=(const FFVCurrentDecayer &);

private:
  
  /**
   * Pointer to FFVVertex
   */
  FFVVertexPtr _theFFVPtr;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of FFVCurrentDecayer. */
template <>
struct BaseClassTrait<Herwig::FFVCurrentDecayer,1> {
  /** Typedef of the first base class of FFVCurrentDecayer. */
  typedef Herwig::GeneralCurrentDecayer NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the FFVCurrentDecayer class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::FFVCurrentDecayer>
  : public ClassTraitsBase<Herwig::FFVCurrentDecayer> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::FFVCurrentDecayer"; }
  /**
   * The name of a file containing the dynamic library where the class
   * FFVCurrentDecayer is implemented. It may also include several, space-separated,
   * libraries if the class FFVCurrentDecayer depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "libHwGeneralDecay.so"; }
};

/** @endcond */

}

#include "FFVCurrentDecayer.icc"

#endif /* HERWIG_FFVCurrentDecayer_H */
