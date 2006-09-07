// -*- C++ -*-
#ifndef HERWIG_QTildeFinder_H
#define HERWIG_QTildeFinder_H
//
// This is the declaration of the QTildeFinder class.
//

#include "Herwig++/Shower/Base/PartnerFinder.h"
#include "Herwig++/Shower/ShowerConfig.h"
#include "ThePEG/Interface/Interfaced.h"
#include "Herwig++/Shower/ShowerVariables.h"
#include "QTildeFinder.fh"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Shower
 *
 *  The QTildeFinder class is responsible for finding the partners and
 *  setting the initial evolution scales for the shower evolution described
 *  in JHEP 0312:045,2003.
 *
 * @see \ref QTildeFinderInterfaces "The interfaces"
 * defined for QTildeFinder.
 */
class QTildeFinder: public PartnerFinder {

public:

  /**
   * The default constructor.
   */
  inline QTildeFinder();

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

protected:

  /**
   * Given a pair of particles, supposedly partners w.r.t. an interaction,
   * this method returns their initial evolution scales as a pair.
   * If something wrong happens, it returns the null (Energy(),Energy()) pair. 
   * This method is used by the above setXXXInitialEvolutionScales 
   * methods.
   */
  //@{
  /**
   *  Calculate the initial evolution scales for two final-state particles
   */
  virtual pair<Energy,Energy> calculateFinalFinalScales(const ShowerPPair &);

  /**
   *  Calculate the initial evolution scales for two initial-state particles
   */
  virtual pair<Energy,Energy> calculateInitialInitialScales(const ShowerPPair &);

  /**
   *  Calculate the initial evolution scales for one initial 
   *  and one final-state particles
   */
  virtual pair<Energy,Energy> calculateInitialFinalScales(const ShowerPPair &,
							  const bool isDecayCase);
  //@}

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
  static NoPIOClassDescription<QTildeFinder> initQTildeFinder;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  QTildeFinder & operator=(const QTildeFinder &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of QTildeFinder. */
template <>
struct BaseClassTrait<Herwig::QTildeFinder,1> {
  /** Typedef of the first base class of QTildeFinder. */
  typedef Herwig::PartnerFinder NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the QTildeFinder class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::QTildeFinder>
  : public ClassTraitsBase<Herwig::QTildeFinder> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::QTildeFinder"; }
  /**
   * The name of a file containing the dynamic library where the class
   * QTildeFinder is implemented. It may also include several, space-separated,
   * libraries if the class QTildeFinder depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#include "QTildeFinder.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "QTildeFinder.tcc"
#endif

#endif /* HERWIG_QTildeFinder_H */
