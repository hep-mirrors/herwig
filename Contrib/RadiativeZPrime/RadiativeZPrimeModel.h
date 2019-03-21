// -*- C++ -*-
#ifndef RADIATIVEZPRIME_RadiativeZPrimeModel_H
#define RADIATIVEZPRIME_RadiativeZPrimeModel_H
//
// This is the declaration of the RadiativeZPrimeModel class.
//

#include "Herwig/Models/StandardModel/StandardModel.h"

namespace RadiativeZPrime {

using namespace ThePEG;

/**
 * The RadiativeZPrimeModel class is the main class for the radiative
 * \f$Z'\f$ model of hep-ph/0501154.
 *
 * @see \ref RadiativeZPrimeModelInterfaces "The interfaces"
 * defined for RadiativeZPrimeModel.
 */
class RadiativeZPrimeModel: public Herwig::StandardModel {

public:

  /**
   * The default constructor.
   */
  inline RadiativeZPrimeModel()  :
    _gZprime(1.), _useZcouplings(true),
    _vnu(1.0), _ve(-0.072), _vu(0.3813), _vd(-0.6907),
    _anu(1.0), _ae(-1.0), _au(1.0), _ad(-1.0) {}

  /**
   *  The coupling of the \f$Z'\f$
   */
  inline double gZprime() const {return _gZprime;}

  /**
   * The vector neutrino-\f$Z^0\f$ coupling.
   */
  inline double zPrimevnu() const {return _vnu;}

  /**
   * The vector charged lepton-\f$Z^0\f$ coupling.
   */
  inline double zPrimeve() const {return _ve;}

  /**
   * The vector up-type-\f$Z^0\f$ coupling.
   */
  inline double zPrimevu() const {return _vu;}

  /**
   * The vector down-type-\f$Z^0\f$ coupling.
   */
  inline double zPrimevd() const {return _vd;}

  /**
   * The axial neutrino-\f$Z^0\f$ coupling.
   */
  inline double zPrimeanu() const {return _anu;}

  /**
   * The axial charged lepton-\f$Z^0\f$ coupling.
   */
  inline double zPrimeae() const {return _ae;}

  /**
   * The axial up-type-\f$Z^0\f$ coupling.
   */
  inline double zPrimeau() const {return _au;}

  /**
   * The axial down-type-\f$Z^0\f$ coupling.
   */
  inline double zPrimead() const {return _ad;}

  /**
   *  The \f$f\bar{f}Z'\f$ vertex
   */
  inline AbstractFFVVertexPtr vertexFFZPrime() const {return _ffZPrimeVertex;}

  /**
   *  The \f$\gamma Z'Z\f$ vertex
   */
  inline AbstractVVVVertexPtr vertexGammaZPrimeZ() const {return _gammaZPrimeZVertex;}

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
  inline virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const {return new_ptr(*this);}
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<RadiativeZPrimeModel> initRadiativeZPrimeModel;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  RadiativeZPrimeModel & operator=(const RadiativeZPrimeModel &) = delete;


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
   *  The free coupling \f$g_{Z'}\f$
   */
  double _gZprime;

  /**
   *  Switch to control the setting of the \f$Z'\f$ couplings to fermions
   */
  bool _useZcouplings;

  /**
   * Vector coupling between a fundamental fermion and \f$Z'\f$.
   */
  double _vnu;

  /**
   * Vector coupling between a fundamental fermion and \f$Z'\f$.
   */
  double _ve;

  /**
   * Vector coupling between a fundamental fermion and \f$Z'\f$.
   */
  double _vu;

  /**
   * Vector coupling between a fundamental fermion and \f$Z'\f$.
   */
  double _vd;

  /**
   * Axial coupling between a fundamental fermions and \f$Z'\f$.
   */
  double _anu;

  /**
   * Axial coupling between a fundamental fermions and \f$Z'\f$.
   */
  double _ae;

  /**
   * Axial coupling between a fundamental fermions and \f$Z'\f$.
   */
  double _au;

  /**
   * Axial coupling between a fundamental fermions and \f$Z'\f$.
   */
  double _ad;

  /**
   *  The \f$f\bar{f}\f$ vertex
   */
  AbstractFFVVertexPtr _ffZPrimeVertex;

  /**
   *  The \f$\gamma Z'Z\f$ vertex
   */
  AbstractVVVVertexPtr _gammaZPrimeZVertex;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of RadiativeZPrimeModel. */
template <>
struct BaseClassTrait<RadiativeZPrime::RadiativeZPrimeModel,1> {
  /** Typedef of the first base class of RadiativeZPrimeModel. */
  typedef Herwig::StandardModel NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the RadiativeZPrimeModel class and the shared object where it is defined. */
template <>
struct ClassTraits<RadiativeZPrime::RadiativeZPrimeModel>
  : public ClassTraitsBase<RadiativeZPrime::RadiativeZPrimeModel> {
  /** Return a platform-independent class name */
  static string className() { return "RadiativeZPrime::RadiativeZPrimeModel"; }
  /**
   * The name of a file containing the dynamic library where the class
   * RadiativeZPrimeModel is implemented. It may also include several, space-separated,
   * libraries if the class RadiativeZPrimeModel depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "RadiativeZPrime.so"; }
};

ThePEG_DECLARE_POINTERS(RadiativeZPrime::RadiativeZPrimeModel,RadiativeZPrimeModelPtr);

/** @endcond */
}

#endif /* RADIATIVEZPRIME_RadiativeZPrimeModel_H */
