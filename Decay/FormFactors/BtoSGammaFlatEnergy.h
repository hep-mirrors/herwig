// -*- C++ -*-
#ifndef HERWIG_BtoSGammaFlatEnergy_H
#define HERWIG_BtoSGammaFlatEnergy_H
//
// This is the declaration of the BtoSGammaFlatEnergy class.
//

#include "BtoSGammaHadronicMass.h"
#include "BtoSGammaFlatEnergy.fh"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Decay
 *
 * The BtoSGammaFlatEnergy class is a model of the hadronic mass 
 * is \f$B\to s\gamma\f$ decays
 * which produces a flat photon energy spectrum and as such is only intended for
 * testing purposes.
 *
 * @see BtoSGammaHadronicMass
 */
class BtoSGammaFlatEnergy: public BtoSGammaHadronicMass {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline BtoSGammaFlatEnergy();

  /**
   * The copy constructor.
   */
  inline BtoSGammaFlatEnergy(const BtoSGammaFlatEnergy &);

  /**
   * The destructor.
   */
  virtual ~BtoSGammaFlatEnergy();
  //@}

public:

  /**
   * Return the hadronic mass.
   * @param mb The mass of the decaying B meson.
   * @param mquark The minimum mass of the hadronic system based on the consistuent quark
   * masses.
   * @return The hadronic mass
   */
  virtual Energy hadronicMass(Energy mb, Energy mquark);

  /**
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   * @param create Whether or not to add a statement creating the object
   */
  virtual void dataBaseOutput(ofstream & os,bool header,bool create) const;

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
   * Check sanity of the object during the setup phase.
   */
  inline virtual void doupdate() throw(UpdateException);

  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  inline virtual void doinit() throw(InitException);

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  inline virtual void doinitrun();

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  inline virtual void dofinish();

  /**
   * Rebind pointer to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   * @throws RebindException if no cloned object was found for a given
   * pointer.
   */
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);

  /**
   * Return a vector of all pointers to Interfaced objects used in this
   * object.
   * @return a vector of pointers.
   */
  inline virtual IVector getReferences();
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<BtoSGammaFlatEnergy> initBtoSGammaFlatEnergy;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  BtoSGammaFlatEnergy & operator=(const BtoSGammaFlatEnergy &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** This template specialization informs ThePEG about the
 *  base classes of BtoSGammaFlatEnergy. */
template <>
struct BaseClassTrait<Herwig::BtoSGammaFlatEnergy,1> {
  /** Typedef of the first base class of BtoSGammaFlatEnergy. */
  typedef Herwig::BtoSGammaHadronicMass NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the BtoSGammaFlatEnergy class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::BtoSGammaFlatEnergy>
  : public ClassTraitsBase<Herwig::BtoSGammaFlatEnergy> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::BtoSGammaFlatEnergy"; }
  /** Return the name of the shared library be loaded to get
   *  access to the BtoSGammaFlatEnergy class and every other class it uses
   *  (except the base class). */
  static string library() { return "libHwFormFactor.so"; }
};

}

#include "BtoSGammaFlatEnergy.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "BtoSGammaFlatEnergy.tcc"
#endif

#endif /* HERWIG_BtoSGammaFlatEnergy_H */
