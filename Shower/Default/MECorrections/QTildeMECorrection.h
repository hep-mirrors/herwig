// -*- C++ -*-
#ifndef HERWIG_QTildeMECorrection_H
#define HERWIG_QTildeMECorrection_H
//
// This is the declaration of the QTildeMECorrection class.
//

#include "Herwig++/Shower/Base/MECorrectionBase.h"
#include "Herwig++/Shower/Base/ShowerTree.h"
#include "QTildeMECorrection.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the QTildeMECorrection class.
 *
 * @see \ref QTildeMECorrectionInterfaces "The interfaces"
 * defined for QTildeMECorrection.
 */
class QTildeMECorrection: public MECorrectionBase {

public:

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class without persistent data.
   */
  static AbstractNoPIOClassDescription<QTildeMECorrection> initQTildeMECorrection;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  QTildeMECorrection & operator=(const QTildeMECorrection &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of QTildeMECorrection. */
template <>
struct BaseClassTrait<Herwig::QTildeMECorrection,1> {
  /** Typedef of the first base class of QTildeMECorrection. */
  typedef Herwig::MECorrectionBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the QTildeMECorrection class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::QTildeMECorrection>
  : public ClassTraitsBase<Herwig::QTildeMECorrection> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::QTildeMECorrection"; }
  /**
   * The name of a file containing the dynamic library where the class
   * QTildeMECorrection is implemented. It may also include several, space-separated,
   * libraries if the class QTildeMECorrection depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMPI.so HwMPIPDF.so HwRemDecayer.so HwShower.so"; }
};

/** @endcond */

}

#include "QTildeMECorrection.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "QTildeMECorrection.tcc"
#endif

#endif /* HERWIG_QTildeMECorrection_H */
