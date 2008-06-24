// -*- C++ -*-
#ifndef HERWIG_LeptonRemnant_H
#define HERWIG_LeptonRemnant_H
//
// This is the declaration of the LeptonRemnant class.
//

#include "ThePEG/PDF/RemnantHandler.h"
#include "LeptonRemnant.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the LeptonRemnant class.
 *
 * @see \ref LeptonRemnantInterfaces "The interfaces"
 * defined for LeptonRemnant.
 */
class LeptonRemnant: public RemnantHandler {

public:

  /**
   * The default constructor.
   */
  inline LeptonRemnant();


  /** @name Virtual functions mandated by the RemnantHandler base class. */
  //@{
  /**
   * Return true if this remnant handler can handle extracting all
   * specified \a partons form the given \a particle.
   */
  virtual bool canHandle(tcPDPtr particle, const cPDVector & partons) const;

  /**
   * Generate momenta. Generates the momenta of the extracted parton
   * in the particle cms (but with the parton \f$x\f$ still the
   * positive light-cone fraction) as given by the last argument, \a
   * p. If the particle is space-like the positive and negative
   * light-cone momenta are \f$\sqrt{-m^2}\f$ and \f$-sqrt{-m^2}\f$
   * respectively. If the \a scale is negative, it means that the \a
   * doScale in the previous call to nDim() was true, otherwise the
   * given scale should be the virtuality of the extracted
   * parton. Generated quantities which are not returned in the
   * momentum may be saved in the PartonBin, \a pb, for later use. In
   * particular, if the nDim() random numbers, \a r, are not enough to
   * generate with weight one, the resulting weight should be stored
   * with the remnantWeight() method of the parton bin.
   */
  virtual Lorentz5Momentum generate(PartonBinInstance & pb, const double * r,
				    Energy2 scale,
				    const LorentzMomentum & p) const;

  /**
   * Generate the momentum of the extracted parton with the \a parent
   * momentum given by the last argument. If the \a scale is negative,
   * it means that the doScale in the previous call to nDim() was
   * true, otherwise the given \a scale should be the virtuality of
   * the extracted parton. \a shat is the total invariant mass squared
   * of the hard sub-system produced by the extracted parton and the
   * primary parton entering from the other side. Generated quantities
   * which are not returned in the momentum may be saved in the
   * PartonBinInstance, \a pb, for later use. In particular, if the
   * nDim() random numbers, \a r, are not enough to generate with
   * weight one, the resulting weight should be stored with the
   * remnantWeight() method of the parton bin.
   */
  virtual Lorentz5Momentum generate(PartonBinInstance & pb, const double * r,
				    Energy2 scale, Energy2 shat,
				    const LorentzMomentum & parent) const;
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

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit() throw(InitException);
  //@}

private:

  /**
   * The minimum energy fraction allowed for a photon remnant.
   */
  double minX;

  /**
   * Easy access to a proton data object.
   */
  tPDPtr photon;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<LeptonRemnant> initLeptonRemnant;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  LeptonRemnant & operator=(const LeptonRemnant &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of LeptonRemnant. */
template <>
struct BaseClassTrait<Herwig::LeptonRemnant,1> {
  /** Typedef of the first base class of LeptonRemnant. */
  typedef RemnantHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the LeptonRemnant class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::LeptonRemnant>
  : public ClassTraitsBase<Herwig::LeptonRemnant> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::LeptonRemnant"; }
  /**
   * The name of a file containing the dynamic library where the class
   * LeptonRemnant is implemented. It may also include several, space-separated,
   * libraries if the class LeptonRemnant depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwLeptonPDF.so"; }
};

/** @endcond */

}

#include "LeptonRemnant.icc"

#endif /* HERWIG_LeptonRemnant_H */
