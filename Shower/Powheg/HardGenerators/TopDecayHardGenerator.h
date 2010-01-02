// -*- C++ -*-
#ifndef HERWIG_TopDecayHardGenerator_H
#define HERWIG_TopDecayHardGenerator_H
//
// This is the declaration of the TopDecayHardGenerator class.
//

#include "Herwig++/Shower/Base/HardestEmissionGenerator.h"
#include "Herwig++/Shower/Couplings/ShowerAlpha.h"
#include "Herwig++/Shower/Base/ShowerProgenitor.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the TopDecayHardGenerator class.
 *
 * @see \ref TopDecayHardGeneratorInterfaces "The interfaces"
 * defined for TopDecayHardGenerator.
 */
class TopDecayHardGenerator: public HardestEmissionGenerator {

public:

  /**
   * The default constructor.
   */
  TopDecayHardGenerator();

  /**
   *  Members which must be overridden in the inheriting classes
   */
  //@{
  /**
   *  Member to generate the hardest emission
   */
  virtual HardTreePtr generateHardest(ShowerTreePtr);

  /**
   *  Member to decide if the inheriting class can handle this process
   */
  virtual bool canHandle(ShowerTreePtr);
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

  /**
   *  Members to calculate the matrix elements
   */
  //@{
  /**
   *  Leading-order matrix element
   */
  double loME(RhoDMatrix rho, pair<tcPDPtr,Lorentz5Momentum> top, 
	      pair<tcPDPtr,Lorentz5Momentum> bottom, 
	      pair<tcPDPtr,Lorentz5Momentum> Wboson,
	      vector<pair<tcPDPtr,Lorentz5Momentum> > leptons);

  /**
   *  Real emission  matrix element
   */
  double realME(RhoDMatrix rho, pair<tcPDPtr,Lorentz5Momentum> top, 
		pair<tcPDPtr,Lorentz5Momentum> bottom, 
		pair<tcPDPtr,Lorentz5Momentum> Wboson, 
		pair<tcPDPtr,Lorentz5Momentum> gluon ,
		vector<pair<tcPDPtr,Lorentz5Momentum> > leptons);
  //@}
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
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<TopDecayHardGenerator> initTopDecayHardGenerator;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  TopDecayHardGenerator & operator=(const TopDecayHardGenerator &);

private:

  /**
   *  Pointer to the object calculating the strong coupling
   */
  ShowerAlphaPtr alphaS_;

  /**
   *  Cut of \f$p_T\f$.
   */
  Energy pTmin_;

  /**
   *  Overestimate for the veto algorithm
   */
  double overEstimate_;

  /**
   *  Power for the sampling
   */
  double power_;

  /**
   *  Vertices for the calculations
   */
  //@{
  /**
   *  W vertex
   */
  AbstractFFVVertexPtr FFWVertex_;

  /**
   *  Gluon vertex
   */
  AbstractFFVVertexPtr FFGVertex_;
  //@}

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of TopDecayHardGenerator. */
template <>
struct BaseClassTrait<Herwig::TopDecayHardGenerator,1> {
  /** Typedef of the first base class of TopDecayHardGenerator. */
  typedef Herwig::HardestEmissionGenerator NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the TopDecayHardGenerator class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::TopDecayHardGenerator>
  : public ClassTraitsBase<Herwig::TopDecayHardGenerator> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::TopDecayHardGenerator"; }
  /**
   * The name of a file containing the dynamic library where the class
   * TopDecayHardGenerator is implemented. It may also include several, space-separated,
   * libraries if the class TopDecayHardGenerator depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "TopDecayHardGenerator.so"; }
};

/** @endcond */

}

#endif /* HERWIG_TopDecayHardGenerator_H */
