// -*- C++ -*-

// (C) 2007-2009 Simon Plaetzer -- sp@particle.uni-karlsruhe.de

#ifndef HERWIG_FourJetCorrelations_H
#define HERWIG_FourJetCorrelations_H
//
// This is the declaration of the FourJetCorrelations class.
//

#include "Analysis2Base.h"
#include "FourJetCorrelations.fh"

namespace Analysis2 {

using namespace ThePEG;

/**\ingroup Analysis2
 * 
 * Four jet correlations at e+e- colliders (preferably LEP)
 *
 * @see \ref FourJetCorrelationsInterfaces "The interfaces"
 * defined for FourJetCorrelations.
 */
class FourJetCorrelations: public Analysis2::Analysis2Base {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline FourJetCorrelations();

  /**
   * The destructor.
   */
  virtual ~FourJetCorrelations();
  //@}

public:

  /**
   * Analyze the given vector of particles. The default version calls
   * analyze(tPPtr) for each of the particles.
   * @param particles the vector of pointers to particles to be analyzed
   */
  virtual void analyze(const tPVector & particles);

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


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).



protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit() throw(InitException);

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  virtual void dofinish();
  //@}

protected:

  /**
   *  Methods to compute the four jet angles, assumes the jets are energy ordered
   */
  //@{
  /**
   *  Compute \f$\cos\chi_{BZ}\f$
   */
  inline double cosChiBZ(const vector<Lorentz5Momentum>&);

  /**
   *  Compute \f$\cos\Phi_{KSW}\f$.
   */ 
  inline double cosPhiKSW(const vector<Lorentz5Momentum>&);
  
  /**
   *  Compute \f$\cos\Theta_{NR}\f$
   */
  inline double cosThetaNR(const vector<Lorentz5Momentum>&); 

  /**
   *  Compute \f$\cos\alpha_{34}\f$
   */
  inline double cosAlpha34(const vector<Lorentz5Momentum>&); 
  //@}

private:

  /**
   * Options string for CosAlpha34
   */
  string _cosAlpha34;

  /**
   * Options string for CosChiBZ
   */
  string _cosChiBZ;

  /**
   * Options string for CosPhiKSW
   */
  string _cosPhiKSW;

  /**
   * Options string for CosThetaNR
   */
  string _cosThetaNR;

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<FourJetCorrelations> initFourJetCorrelations;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  FourJetCorrelations & operator=(const FourJetCorrelations &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of FourJetCorrelations. */
template <>
struct BaseClassTrait<Analysis2::FourJetCorrelations,1> {
  /** Typedef of the first base class of FourJetCorrelations. */
  typedef Analysis2::Analysis2Base NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the FourJetCorrelations class and the shared object where it is defined. */
template <>
struct ClassTraits<Analysis2::FourJetCorrelations>
  : public ClassTraitsBase<Analysis2::FourJetCorrelations> {
  /** Return a platform-independent class name */
  static string className() { return "Analysis2::FourJetCorrelations"; }
  /**
   * The name of a file containing the dynamic library where the class
   * FourJetCorrelations is implemented. It may also include several, space-separated,
   * libraries if the class FourJetCorrelations depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "Analysis2.so"; }
};

/** @endcond */

}

#include "FourJetCorrelations.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "FourJetCorrelations.tcc"
#endif

#endif /* HERWIG_FourJetCorrelations_H */
