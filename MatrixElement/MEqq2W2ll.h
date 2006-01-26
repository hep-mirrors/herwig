// -*- C++ -*-
#ifndef HERWIG_MEqq2W2ll_H
#define HERWIG_MEqq2W2ll_H
//
// This is the declaration of the MEqq2W2ll class.
//

#include "ThePEG/MatrixElement/ME2to2Base.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Utilities/Rebinder.h"
#include "MEqq2W2ll.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the MEqq2W2ll class.
 *
 * @see \ref MEqq2W2llInterfaces "The interfaces"
 * defined for MEqq2W2ll.
 */
class MEqq2W2ll: public ME2to2Base {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline MEqq2W2ll();

  /**
   * The copy constructor.
   */
  inline MEqq2W2ll(const MEqq2W2ll &);

  /**
   * The destructor.
   */
  virtual ~MEqq2W2ll();
  //@}

public:

  /** @name Virtual functions required by the MEBase class. */
  //@{
  /**
   * Return the order in \f$\alpha_S\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInAlphaS() const;

  /**
   * Return the order in \f$\alpha_{EW}\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInAlphaEW() const;

  /**
   * The matrix element for the kinematical configuration
   * previously provided by the last call to setKinematics(), suitably
   * scaled by sHat() to give a dimension-less number.
   * @return the matrix element scaled with sHat() to give a
   * dimensionless number.
   */
  virtual double me2() const;

  /**
   * Return the scale associated with the last set phase space point.
   */
  virtual Energy2 scale() const;

  /**
   * Add all possible diagrams with the add() function.
   */
  virtual void getDiagrams() const;

  /**
   * Get diagram selector. With the information previously supplied with the
   * setKinematics method, a derived class may optionally
   * override this method to weight the given diagrams with their
   * (although certainly not physical) relative probabilities.
   * @param dv the diagrams to be weighted.
   * @return a Selector relating the given diagrams to their weights.
   */
  inline virtual Selector<DiagramIndex> diagrams(const DiagramVector & dv) const;

  /**
   * Return a Selector with possible colour geometries for the selected
   * diagram weighted by their relative probabilities.
   * @param diag the diagram chosen.
   * @return the possible colour geometries weighted by their
   * relative probabilities.
   */
  virtual Selector<const ColourLines *>
  colourGeometries(tcDiagPtr diag) const;
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
  static ClassDescription<MEqq2W2ll> initMEqq2W2ll;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEqq2W2ll & operator=(const MEqq2W2ll &);

private:
  // pointers to the SM vertices
  Ptr<Herwig::Helicity::FFVVertex>::transient_pointer _theFFWVertex;
  // pointers to the intermediate particles particle data objects
  tcPDPtr _Wp,_Wm;
  // allowed incoming particles
  int _maxflavour;

  mutable double lastCont;
  mutable double lastBW;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEqq2W2ll. */
template <>
struct BaseClassTrait<Herwig::MEqq2W2ll,1> {
  /** Typedef of the first base class of MEqq2W2ll. */
  typedef ME2to2Base NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEqq2W2ll class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEqq2W2ll>
  : public ClassTraitsBase<Herwig::MEqq2W2ll> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::MEqq2W2ll"; }
  /** Return the name(s) of the shared library (or libraries) be loaded to get
   *  access to the MEqq2W2ll class and any other class on which it depends
   *  (except the base class). */
  static string library() { return "HwME.so"; }
};

/** @endcond */

}

#include "MEqq2W2ll.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "MEqq2W2ll.tcc"
#endif

#endif /* HERWIG_MEqq2W2ll_H */
