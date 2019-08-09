// -*- C++ -*-
#ifndef HERWIG_FFgx2ggxDipoleKernel_H
#define HERWIG_FFgx2ggxDipoleKernel_H
//
// This is the declaration of the FFgx2ggxDipoleKernel class.
//

#include "DipoleSplittingKernel.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup DipoleShower
 * \author Simon Platzer
 *
 * \brief FFgx2ggxDipoleKernel implements the g -> gg
 * splitting off a final-final dipole
 *
 */
class FFgx2ggxDipoleKernel: public DipoleSplittingKernel {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  FFgx2ggxDipoleKernel();

  /**
   * The destructor.
   */
  virtual ~FFgx2ggxDipoleKernel();
  //@}

public:

  /**
   * Return true, if this splitting kernel
   * applies to the given dipole index.
   */
  virtual bool canHandle(const DipoleIndex&) const;

  /**
   * Return true, if this splitting kernel is
   * the same for the given index a, as the given
   * splitting kernel for index b.
   */
  virtual bool canHandleEquivalent(const DipoleIndex& a,
				   const DipoleSplittingKernel& sk,
				   const DipoleIndex& b) const;

  /**
   * Return the emitter data after splitting, given
   * a dipole index.
   */
  virtual tcPDPtr emitter(const DipoleIndex&) const;

  /**
   * Return the emission data after splitting, given
   * a dipole index.
   */
  virtual tcPDPtr emission(const DipoleIndex&) const;

  /**
   * Return the spectator data after splitting, given
   * a dipole index.
   */
  virtual tcPDPtr spectator(const DipoleIndex&) const;

  /**
   * Evaluate this splitting kernel for the given
   * dipole splitting.
   */
  virtual double evaluate(const DipoleSplittingInfo&) const;

  /**
   * Evaluate rho_ii' V_ijk V*_i'jk / equivalent for initial-state splitting,
   * required for generating spin-correlated azimuthal angles.
   **/
  virtual vector< pair<int, Complex> >  generatePhi( const DipoleSplittingInfo& dInfo, const RhoDMatrix& rho) const;
   
  /**
   * Return the completely spin-unaveraged (i.e. spin-indexed) splitting kernel.
   **/
  virtual DecayMEPtr matrixElement(const DipoleSplittingInfo& dInfo) const;
  
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
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<FFgx2ggxDipoleKernel> initFFgx2ggxDipoleKernel;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  FFgx2ggxDipoleKernel & operator=(const FFgx2ggxDipoleKernel &) = delete;

  /**
   * Asymmetry option for final state gluon splittings.
   */
  
  int theAsymmetryOption=1;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of FFgx2ggxDipoleKernel. */
template <>
struct BaseClassTrait<Herwig::FFgx2ggxDipoleKernel,1> {
  /** Typedef of the first base class of FFgx2ggxDipoleKernel. */
  typedef Herwig::DipoleSplittingKernel NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the FFgx2ggxDipoleKernel class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::FFgx2ggxDipoleKernel>
  : public ClassTraitsBase<Herwig::FFgx2ggxDipoleKernel> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::FFgx2ggxDipoleKernel"; }
  /**
   * The name of a file containing the dynamic library where the class
   * FFgx2ggxDipoleKernel is implemented. It may also include several, space-separated,
   * libraries if the class FFgx2ggxDipoleKernel depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwDipoleShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_FFgx2ggxDipoleKernel_H */
