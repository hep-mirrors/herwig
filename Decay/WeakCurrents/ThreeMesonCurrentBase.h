// -*- C++ -*-
#ifndef HERWIG_ThreeMesonCurrentBase_H
#define HERWIG_ThreeMesonCurrentBase_H
// This is the declaration of the ThreeMesonCurrentBase class.

#include "WeakDecayCurrent.h"
#include "ThreeMesonCurrentBase.fh"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 *  This is the base class for the three meson decays of the weak current.
 *  It is designed so that the currents for the following modes can be implemented
 *  in classes inheriting from this. 
 *
 * - \f$    \pi^-  \pi^-    \pi^+ \f$, (imode=0)
 * - \f$    \pi^0  \pi^0    \pi^- \f$, (imode=1)
 * - \f$    K^-   \pi^-    K^+ \f$, (imode=2)
 * - \f$    K^0   \pi^-    \bar{K}^0\f$, (imode=3)
 * - \f$    K^-   \pi^0    K^0 \f$, (imode=4)
 * - \f$    \pi^0  \pi^0    K^- \f$, (imode=5)
 * - \f$    K^-   \pi^-    \pi^+ \f$, (imode=6)
 * - \f$    \pi^-  \bar{K}^0  \pi^0 \f$, (imode=7)
 * - \f$    \pi^-  \pi^0    \eta \f$, (imode=8)
 *
 * obvioulsly there are other modes with three pseudoscalar mesons for the decay
 * of the weak current but this model original came from \f$\tau\f$ decay where
 * these are the only modes.
 *
 *  In this case the current is given by
 *  \f[ J^\mu = \left(g^{\mu\nu}-\frac{q^\mu q^\nu}{q^2}\right)
 *   \left[F_1(p_2-p_3)^\mu +F_2(p_3-p_1)^\mu+F_3(p_1-p_2)^\mu\right]
 *  +q^\mu F_4
 *  +F_5\epsilon^{\mu\alpha\beta\gamma}p_1^\alpha p_2^\beta p_3^\gamma
 *  \f]
 * where
 * - \f$p_{1,2,3}\f$ are the momenta of the mesons in the order given above.
 * - \f$F_1,F_2,F_3,F_4,F_5\f$ are the form factors which must be 
 *  calculated in the calculateFormFactors member which should be implemented
 * in classes inheriting from this.
 *
 * @see WeakDecayCurrent.
 *  
 * \author Peter Richardson
 *
 */
class ThreeMesonCurrentBase: public WeakDecayCurrent {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor
   */
  inline ThreeMesonCurrentBase();

  /**
   * Copy constructor
   */
  inline ThreeMesonCurrentBase(const ThreeMesonCurrentBase &);

  /**
   * Destructor
   */
  virtual ~ThreeMesonCurrentBase();
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
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();

public:


  /**
   * Hadronic current. This version returns the hadronic current described above.
   * @param vertex Construct the information needed for spin correlations
   * @param imode The mode
   * @param ichan The phase-space channel the current is needed for.
   * @param scale The invariant mass of the particles in the current.
   * @param decay The decay products
   * @return The current. 
   */
  virtual vector<LorentzPolarizationVector>  current(bool vertex, const int imode,
						     const int ichan,Energy & scale,
						     const ParticleVector & decay) const;

  /**
   * Accept the decay. Checks the mesons against the list.
   * @param id The id's of the particles in the current.
   * @return Can this current have the external particles specified.
   */
  virtual bool accept(vector<int> id);

  /**
   * Return the decay mode number for a given set of particles in the current. 
   * Checks the mesons against the list.
   * @param id The id's of the particles in the current.
   * @return The number of the mode
   */
  virtual unsigned int decayMode(vector<int> id);

  /**
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   * @param create Whether or not to add a statement creating the object
   */
  virtual void dataBaseOutput(ofstream & os,bool header,bool create) const;

protected:

  /**
   * can a particular decayer handle this type of mode
   * @param imode The mode number as given above
   * @return Whether this mode can be handled.
   */
  virtual bool acceptMode(int imode) const=0;

  /**
   * Calculate the form factor for the current.
   * @param ichan The phase space channel
   * @param imode The mode
   * @param q2 The scale \f$q^2\f$ for the current.
   * @param s1 The invariant mass squared of particles 2 and 3, \f$s_1=m^2_{23}\f$.
   * @param s2 The invariant mass squared of particles 1 and 3, \f$s_2=m^2_{13}\f$.
   * @param s3 The invariant mass squared of particles 1 and 2, \f$s_3=m^2_{12}\f$.
   * @param F1 The form factor \f$F_1\f$.
   * @param F2 The form factor \f$F_2\f$.
   * @param F3 The form factor \f$F_3\f$.
   * @param F4 The form factor \f$F_4\f$.
   * @param F5 The form factor \f$F_5\f$.
   */
  virtual void calculateFormFactors(const int ichan,const int imode,
				    Energy2 q2,Energy2 s1,Energy2 s2,Energy2 s3,
				    Complex&F1,Complex&F2,Complex&F3,
				    Complex&F4,Complex&F5) const=0;


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
   * Initialize this object to the begining of the run phase.
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
   * @throws RebindException if no cloned object was found for a given pointer.
   */
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);

  /**
   * Return a vector of all pointers to Interfaced objects used in
   * this object.
   * @return a vector of pointers.
   */
  inline virtual IVector getReferences();
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class with persistent data.
   */
  static AbstractClassDescription<ThreeMesonCurrentBase> initThreeMesonCurrentBase;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ThreeMesonCurrentBase & operator=(const ThreeMesonCurrentBase &);

};

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of ThreeMesonCurrentBase.
 */
template <>
 struct BaseClassTrait<Herwig::ThreeMesonCurrentBase,1> {
  /** Typedef of the base class of ThreeMesonCurrentBase. */
  typedef Herwig::WeakDecayCurrent NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
 struct ClassTraits<Herwig::ThreeMesonCurrentBase>
  : public ClassTraitsBase<Herwig::ThreeMesonCurrentBase> {
   /** Return the class name.*/
  static string className() { return "Herwig++::ThreeMesonCurrentBase"; }
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwWeakCurrents.so"; }

};

}

#include "ThreeMesonCurrentBase.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ThreeMesonCurrentBase.tcc"
#endif

#endif /* HERWIG_ThreeMesonCurrentBase_H */
