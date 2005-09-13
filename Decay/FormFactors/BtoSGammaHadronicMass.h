// -*- C++ -*-
#ifndef HERWIG_BtoSGammaHadronicMass_H
#define HERWIG_BtoSGammaHadronicMass_H
//
// This is the declaration of the BtoSGammaHadronicMass class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "BtoSGammaHadronicMass.fh"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Decay
 *
 * The BtoSGammaHadronicMass class is the base class for the implementation
 * of models of the hadronic mass spectrum in \f$B\to s\gamma\f$ decays.
 * Classes inheriting from this class should implement the hadronicMass()
 * member which should return a value of the hadronic mass selected from
 * the distribution.
 *
 * The parameters relating to the minimum and maximum values of the mass are stored
 * in this class.
 *
 */
class BtoSGammaHadronicMass: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline BtoSGammaHadronicMass();

  /**
   * The copy constructor.
   */
  inline BtoSGammaHadronicMass(const BtoSGammaHadronicMass &);

  /**
   * The destructor.
   */
  virtual ~BtoSGammaHadronicMass();
  //@}

public:

  /**
   * Virtual member which must be implemented in classes inheriting from this
   * class to return the hadronic mass.
   * @param mb The mass of the decaying B meson
   * @param mquark The minimum mass of the hadronic system based on the consistuent quark
   * masses.
   * @return The hadronic mass
   */
  virtual Energy hadronicMass(Energy mb,Energy mquark) =0;

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

  /** @name Access to the limits on the mass. */
  //@{
  /**
   *  Minimum mass
   */
  inline Energy minMass() const;

  /**
   *  Maximum mass
   */
  inline Energy maxMass() const;
  //@}

  /** @name Functions for the fermi motion needed in classes inheriting from this */
  //@{
  /**
   * Exponential function of the form, \f$N(1-x)^ae^{-3\bar{\Lambda}^2x/\lambda_1}\f$,
   * where 
   * \f$x=k_+/\bar{\Lambda}\f$ taken from hep-ph/9805303
   * @param scale The energy scale, \f$k_+\f$, at which to evaluate the function.
   * @param lambda The hadronic scale, \f$\bar{\Lambda}\f$
   * @param a The shape parameter, \f$a\f$.
   * @param norm The normalisation, \f$N\f$.
   * @param lambda1 Scale related to kinetic energy of b quark, \f$\lambda_1\f$.
   */
  inline double exponentialFermiFunction(Energy scale,Energy lambda,
					 double a,double norm,Energy2 lambda1 );
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
   * Indicates that this is an abstract class with persistent data.
   */
  static AbstractClassDescription<BtoSGammaHadronicMass> initBtoSGammaHadronicMass;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  BtoSGammaHadronicMass & operator=(const BtoSGammaHadronicMass &);

private:

  /**
   * The minimum value of the hadronic mass
   */
  Energy _minMass;

  /**
   * The maximum value of the hadronic mass
   */
  Energy _maxMass;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** This template specialization informs ThePEG about the
 *  base classes of BtoSGammaHadronicMass. */
template <>
struct BaseClassTrait<Herwig::BtoSGammaHadronicMass,1> {
  /** Typedef of the first base class of BtoSGammaHadronicMass. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the BtoSGammaHadronicMass class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::BtoSGammaHadronicMass>
  : public ClassTraitsBase<Herwig::BtoSGammaHadronicMass> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::BtoSGammaHadronicMass"; }
  /** Return the name of the shared library be loaded to get
   *  access to the BtoSGammaHadronicMass class and every other class it uses
   *  (except the base class). */
  static string library() { return "libHwFormFactor.so"; }
};

}

#include "BtoSGammaHadronicMass.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "BtoSGammaHadronicMass.tcc"
#endif

#endif /* HERWIG_BtoSGammaHadronicMass_H */
