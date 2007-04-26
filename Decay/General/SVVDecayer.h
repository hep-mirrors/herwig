// -*- C++ -*-
#ifndef HERWIG_SVVDecayer_H
#define HERWIG_SVVDecayer_H
//
// This is the declaration of the SVVDecayer class.
//

#include "GeneralTwoBodyDecayer.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "Herwig++/Helicity/Vertex/Scalar/VVSVertex.fh"
#include "SVVDecayer.fh"

namespace Herwig {
using namespace ThePEG;
using Helicity::VVSVertexPtr; 

/** \ingroup Decay
 * This SVVDecayer class implements the decay of a scalar to 
 * 2 vector bosons using the tree level VVSVertex. It inherits from 
 * GeneralTwoBodyDecayer and implements the virtual member functions me2() 
 * and partialWidth(). It also stores a pointer to the VVSVertex.
 *
 * @see GeneralTwoBodyDecayer 
 * 
 */
class SVVDecayer: public GeneralTwoBodyDecayer {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline SVVDecayer();

  /**
   * The destructor.
   */
  virtual ~SVVDecayer();
  //@}

public:

  /** @name Virtual functions required by the Decayer class. */
  //@{
  /**
   * Return the matrix element squared for a given mode and phase-space channel.   * @param vertex Output the information on the vertex for spin correlations
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
   * @param outa Pointer to incoming particle data object
   * @param outb Pointer to incoming particle data object
   */
  virtual double partialWidth(const PDPtr inpart,
			      const PDPtr outa,
			      const PDPtr outb) const;
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
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  inline virtual void doinit() throw(InitException);
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<SVVDecayer> initSVVDecayer;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SVVDecayer & operator=(const SVVDecayer &);

private:
  
  /**
   * Store pointer to VVSVertex that is set in doinit by 
   * typecast from vertex pointer in base class
   */
  VVSVertexPtr _theVVSPtr; 
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/// \if TRAITSPECIALIZATIONS

/** This template specialization informs ThePEG about the
 *  base classes of SVVDecayer. */
template <>
struct BaseClassTrait<Herwig::SVVDecayer,1> {
  /** Typedef of the first base class of SVVDecayer. */
  typedef Herwig::GeneralTwoBodyDecayer NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SVVDecayer class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::SVVDecayer>
  : public ClassTraitsBase<Herwig::SVVDecayer> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::SVVDecayer"; }
  /** Return the name of the shared library be loaded to get
   *  access to the SVVDecayer class and every other class it uses
   *  (except the base class). */
  static string library() { return "libHwGeneralDecay.so"; }
};

/// \endif

}

#include "SVVDecayer.icc"

#endif /* HERWIG_SVVDecayer_H */
