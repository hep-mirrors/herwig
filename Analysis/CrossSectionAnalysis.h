// -*- C++ -*-
#ifndef Herwig_CrossSectionAnalysis_H
#define Herwig_CrossSectionAnalysis_H
//
// This is the declaration of the CrossSectionAnalysis class.
//

#include "ThePEG/Handlers/AnalysisHandler.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the CrossSectionAnalysis class.
 *
 * @see \ref CrossSectionAnalysisInterfaces "The interfaces"
 * defined for CrossSectionAnalysis.
 */
class CrossSectionAnalysis: public AnalysisHandler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  CrossSectionAnalysis();

  /**
   * The destructor.
   */
  virtual ~CrossSectionAnalysis();
  //@}

protected:

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  virtual void dofinish();

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
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  CrossSectionAnalysis & operator=(const CrossSectionAnalysis &) = delete;

};

}

#endif /* Herwig_CrossSectionAnalysis_H */
