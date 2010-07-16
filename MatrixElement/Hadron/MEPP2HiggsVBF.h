// -*- C++ -*-
#ifndef HERWIG_MEPP2HiggsVBF_H
#define HERWIG_MEPP2HiggsVBF_H
//
// This is the declaration of the MEPP2HiggsVBF class.
//

#include "Herwig++/MatrixElement/MEfftoffH.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The MEPP2HiggsVBF class provides the matrix elements for the
 * production of the Higgs boson via the vector boson fusion mechanism
 * in hadron collisions
 *
 * @see \ref MEPP2HiggsVBFInterfaces "The interfaces"
 * defined for MEPP2HiggsVBF.
 */
class MEPP2HiggsVBF: public MEfftoffH {

public:

  /**
   * The default constructor.
   */
  MEPP2HiggsVBF();

  /** @name Virtual functions required by the MEBase class. */
  //@{

  /**
   * Add all possible diagrams with the add() function.
   */
  virtual void getDiagrams() const;
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

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const { return new_ptr(*this); }

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const { return new_ptr(*this); }
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<MEPP2HiggsVBF> initMEPP2HiggsVBF;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEPP2HiggsVBF & operator=(const MEPP2HiggsVBF &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEPP2HiggsVBF. */
template <>
struct BaseClassTrait<Herwig::MEPP2HiggsVBF,1> {
  /** Typedef of the first base class of MEPP2HiggsVBF. */
  typedef Herwig::MEfftoffH NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEPP2HiggsVBF class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEPP2HiggsVBF>
  : public ClassTraitsBase<Herwig::MEPP2HiggsVBF> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEPP2HiggsVBF"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MEPP2HiggsVBF is implemented. It may also include several, space-separated,
   * libraries if the class MEPP2HiggsVBF depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMEHadron.so"; }
};

/** @endcond */

}

#endif /* HERWIG_MEPP2HiggsVBF_H */
