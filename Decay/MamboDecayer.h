// -*- C++ -*-
#ifndef HERWIG_MamboDecayer_H
#define HERWIG_MamboDecayer_H
//
// This is the declaration of the MamboDecayer class.
//

#include "ThePEG/PDT/Decayer.h"
#include "ThePEG/PDT/DecayMode.h"
#include "MamboDecayer.fh"
#include "Herwig++/Utilities/GlobalParameters.h"
namespace Herwig {
  using namespace ThePEG;
  
  /**
   * The MamboDecayer class inherits from the Decayer class in 
   * ThePEG and implements the algorithm of R.Kleiss and 
   * W.J.Stirling for massive multi-particle phase-space
   * decays 
   */
class MamboDecayer: public Decayer {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline MamboDecayer();

  /**
   * The copy constructor.
   */
  inline MamboDecayer(const MamboDecayer &);

  /**
   * The destructor.
   */
  virtual ~MamboDecayer();
  //@}

public:

  /** @name Virtual functions required by the Decayer class. */
  //@{
  /**
   * Check if this decayer can perfom the decay specified by the
   * given decay mode.
   * @param dm the DecayMode describing the decay.
   * @return true if this decayer can handle the given mode, otherwise false.
   */
  virtual bool accept(const DecayMode & dm) const;

  /**
   * Perform a decay for a given DecayMode and a given Particle instance.
   * @param dm the DecayMode describing the decay.
   * @param p the Particle instance to be decayed.
   * @return a ParticleVector containing the decay products.
   */
  virtual ParticleVector decay(const DecayMode & dm, const Particle & p) const;
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

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Check sanity of the object during the setup phase.
   */
  inline virtual void doupdate() throw(UpdateException);

  /**
   * Initialize this object after the setup phase before saving and
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
  static ClassDescription<MamboDecayer> initMamboDecayer;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MamboDecayer & operator=(const MamboDecayer &);

private:

   /**
     *Set array of mometum to particles
     *@param mom Momentum set to be distributed over phase-space
     *@param totalMom Total initial momentum required to know total mass
     *@return wgt The weight of the configuration
     **/
  
  double calculateMomentum(vector<Lorentz5Momentum> & mom,
			   const Energy & comEn) const;
  /**
   *Bessel function of first kind(n=0)
   *@param x Argument of Bessel Function 
   *@return long double
   **/
  inline long double BessI0(const long double & x) const;
  
  /**
   *Bessel function of first kind(n=1)
   *@param x Argument of Bessel Function 
   *@return long double
   **/
  inline long double BessI1(const long double & x) const;
  
  /**
   *Bessel function of second kind(n=0)
   *@param x Argument of Bessel Function 
   *@return long double
   **/
  inline long double BessK0(const long double & x) const;
  
  /**
   *Bessel function of second kind(n=1)
   *@param x Argument of Bessel Function 
   *@return long double
   **/
  inline long double BessK1(const long double & x) const;
  
  /**
   * Derivative of ratio K0/K1
   *@param x Argument of Bessel Function 
   *@return long double
   **/
  inline long double BessKPrime(const long double & x) const;

  /**
   * Maximum weight
   */
  double _maxweight;

  /**
   * Pointer to global paramters object
   */
  Ptr<GlobalParameters>::pointer global;
};

}

// CLASSDOC OFF

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** This template specialization informs ThePEG about the
 *  base classes of MamboDecayer. */
template <>
struct BaseClassTrait<Herwig::MamboDecayer,1> {
  /** Typedef of the first base class of MamboDecayer. */
  typedef Decayer NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MamboDecayer class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MamboDecayer>
  : public ClassTraitsBase<Herwig::MamboDecayer> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::MamboDecayer"; }
  /** Return the name of the shared library be loaded to get
   *  access to the MamboDecayer class and every other class it uses
   *  (except the base class). */
  static string library() { return "MamboDecayer.so"; }
};

}

#include "MamboDecayer.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "MamboDecayer.tcc"
#endif

#endif /* HERWIG_MamboDecayer_H */
