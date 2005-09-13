// -*- C++ -*-
#ifndef HERWIG_KiselevBcFormFactor_H
#define HERWIG_KiselevBcFormFactor_H
//
// This is the declaration of the KiselevBcFormFactor class.
//

#include "ScalarFormFactor.h"
#include "KiselevBcFormFactor.fh"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 * The KiselevBcFormFactor class implements the form factors from hep-ph/0211021
 * for the decays of \f$B_c\f$ mesons.
 *
 * @see ScalarFormFactor
 */
class KiselevBcFormFactor: public ScalarFormFactor {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  KiselevBcFormFactor();

  /**
   * The copy constructor.
   */
  inline KiselevBcFormFactor(const KiselevBcFormFactor &);

  /**
   * The destructor.
   */
  virtual ~KiselevBcFormFactor();
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

public:

  /** @name Form Factors */
  //@{
  /**
   * The form factor for the weak decay of a scalar to a scalar.  
   * This method is virtual and must be implementented in classes
   * inheriting from this which include scalar to scalar form factors.
   * @param q2 The scale \f$q^2\f$.
   * @param iloc The location in the form factor list.
   * @param id0 The PDG code of the incoming meson.
   * @param id1 The PDG code of the outgoing meson.
   * @param m0 The mass of the incoming meson.
   * @param m1 The mass of the outgoing meson.
   * @param f0 The form factor \f$f_0\f$. 
   * @param fp The form factor \f$f_+\f$.
   */
  virtual void ScalarScalarFormFactor(Energy2 q2,unsigned int iloc,int id0,int id1,
				      Energy m0,Energy m1,Complex & f0,
				      Complex & fp) const;

  /**
   * The form factor for the weak decay of a scalar to a vector. This method is virtual
   * and must be implemented in classes inheriting from this which include scalar to
   * vector form factors.
   * @param q2 The scale \f$q^2\f$.
   * @param iloc The location in the form factor list.
   * @param id0 The PDG code of the incoming meson.
   * @param id1 The PDG code of the outgoing meson.
   * @param m0 The mass of the incoming meson.
   * @param m1 The mass of the outgoing meson.
   * @param A0 The form factor \f$A_0\f$
   * @param A1 The form factor \f$A_1\f$
   * @param A2 The form factor \f$A_2\f$
   * @param V  The form factor \f$V\f$
   */
  virtual void ScalarVectorFormFactor(Energy2 q2, unsigned int iloc, int id0, int id1,
				      Energy m0, Energy m1,Complex & A0,
				      Complex & A1,Complex & A2, Complex & V) const;
  //@}

  /**
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   * @param create Whether or not to add a statement creating the object
   */
  virtual void dataBaseOutput(ofstream & os,bool header,bool create) const;

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
  static ClassDescription<KiselevBcFormFactor> initKiselevBcFormFactor;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  KiselevBcFormFactor & operator=(const KiselevBcFormFactor &);

private:

  /**
   *  The value of the \f$f_+\f$ form factor evaluated at \f$q^2=0\f$.
   */
  vector<double> _fp;

  /**
   *  The value of the \f$f_-\f$ form factor evaluated at \f$q^2=0\f$.
   */
  vector<double> _fm;

  /**
   *  The value of the \f$F_V\f$ form factor evaluated at \f$q^2=0\f$.
   */
  vector<InvEnergy> _FV;

  /**
   *  The value of the \f$F_0^A\f$ form factor evaluated at \f$q^2=0\f$.
   */
  vector<Energy> _F0A;

  /**
   *  The value of the \f$F_+^A\f$ form factor evaluated at \f$q^2=0\f$.
   */
  vector<InvEnergy> _FpA;

  /**
   *  The value of the \f$F_-^A\f$ form factor evaluated at \f$q^2=0\f$.
   */
  vector<InvEnergy> _FmA;

  /**
   *  The pole mass for the \f$f_+\f$ form factor
   */
  vector<Energy> _Mfp;

  /**
   *  The pole mass for the \f$f_-\f$ form factor
   */
  vector<Energy> _Mfm;

  /**
   *  The pole mass for the \f$F_V\f$ form factor
   */
  vector<Energy> _MFV;

  /**
   *  The pole mass for the \f$F_0^A\f$ form factor
   */
  vector<Energy> _MF0A;

  /**
   *  The pole mass for the \f$F_+^A\f$ form factor
   */
  vector<Energy> _MFpA;

  /**
   *  The pole mass for the \f$F_-^A\f$ form factor
   */
  vector<Energy> _MFmA;

};

}

// CLASSDOC OFF

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** This template specialization informs ThePEG about the
 *  base classes of KiselevBcFormFactor. */
template <>
 struct BaseClassTrait<Herwig::KiselevBcFormFactor,1> {
  /** Typedef of the first base class of KiselevBcFormFactor. */
   typedef Herwig::ScalarFormFactor NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the KiselevBcFormFactor class and the shared object where it is defined. */
template <>
 struct ClassTraits<Herwig::KiselevBcFormFactor>
  : public ClassTraitsBase<Herwig::KiselevBcFormFactor> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::KiselevBcFormFactor"; }
  /** Return the name of the shared library be loaded to get
   *  access to the KiselevBcFormFactor class and every other class it uses
   *  (except the base class). */
  static string library() { return "libHwFormFactor.so"; }
};

}

#include "KiselevBcFormFactor.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "KiselevBcFormFactor.tcc"
#endif

#endif /* HERWIG_KiselevBcFormFactor_H */
