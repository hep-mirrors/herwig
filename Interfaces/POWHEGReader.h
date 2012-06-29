// -*- C++ -*-
#ifndef HERWIG_POWHEGReader_H
#define HERWIG_POWHEGReader_H
//
// This is the declaration of the POWHEGReader class.
//

#include "ThePEG/LesHouches/LesHouchesFileReader.h"
#include "POWHEGReader.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the POWHEGReader class.
 *
 * @see \ref POWHEGReaderInterfaces "The interfaces"
 * defined for POWHEGReader.
 */
class POWHEGReader: public LesHouchesFileReader {

public:

  /**
   * The default constructor.
   */
  inline POWHEGReader();

  /** @name Virtual functions specified by the LesHouchesReader base class. */
  //@{
  /**
   * Open a file with events. Derived classes should overwrite it and
   * first calling it before reading in the run information into the
   * corresponding protected variables.
   */
  virtual void open();

  /**
   * Close the file from which events have been read.
   */
  virtual void close();

  /**
   * Read the next event from the file or stream into the
   * corresponding protected variables. Return false if there is no
   * more events or if this was not a LHF event file.
   */
  virtual bool doReadEvent();
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

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<POWHEGReader> initPOWHEGReader;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  POWHEGReader & operator=(const POWHEGReader &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of POWHEGReader. */
template <>
struct BaseClassTrait<Herwig::POWHEGReader,1> {
  /** Typedef of the first base class of POWHEGReader. */
  typedef LesHouchesFileReader NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the POWHEGReader class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::POWHEGReader>
  : public ClassTraitsBase<Herwig::POWHEGReader> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::POWHEGReader"; }
  /**
   * The name of a file containing the dynamic library where the class
   * POWHEGReader is implemented. It may also include several, space-separated,
   * libraries if the class POWHEGReader depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwPOWHEG.so"; }
};

/** @endcond */

}

#include "POWHEGReader.icc"

#endif /* HERWIG_POWHEGReader_H */
