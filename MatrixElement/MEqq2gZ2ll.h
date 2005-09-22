// -*- C++ -*-
#ifndef ThePEG_MEqq2gZ2ll_H
#define ThePEG_MEqq2gZ2ll_H
//
// This is the declaration of the MEqq2gZ2ll class.

#include "ThePEG/MatrixElement/ME2to2QCD.h"

using namespace ThePEG;

namespace Herwig {

/** \ingroup MatrixElement
 *  The MEqq2gZ2ll class implements the e+ e- -> gamma/Z0 q qbar 
 *  matrix element. Both the continuum and Z0 pole term as well as 
 *  the interference term is included. Although not a strict QCD 
 *  matrix element the cass inherits from ME2to2Base, mainly to 
 *  inherit the parameter for the number of active quark flavours.
 *
 *  @see ME2to2QCD
 */
class MEqq2gZ2ll: public ME2to2QCD {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  MEqq2gZ2ll();

  /**
   * The copy constructor.
   */
  MEqq2gZ2ll(const MEqq2gZ2ll &);

  /**
   * The destructor.
   */
  virtual ~MEqq2gZ2ll();
  //@}

public:

  /**
   * Return the order in \f$\alpha_S\f$ in which this matrix element
   * is given.
   */
  virtual unsigned int orderInAlphaS() const;

  /**
   * Return the order in \f$\alpha_{EM}\f$ in which this matrix
   * element is given. Returns 0.
   */
  virtual unsigned int orderInAlphaEW() const;

  /**
   * Return the matrix element for the kinematical configuation
   * previously provided by the last call to setKinematics().
   */
  virtual double me2() const;

  /**
   * Add all possible diagrams with the add() function.
   */
  virtual void getDiagrams() const;

  /**
   * Return a Selector with possible colour geometries for the selected
   * diagram weighted by their relative probabilities.
   */
  virtual Selector<const ColourLines *>
  colourGeometries(tcDiagPtr diag) const;

  /**
   * Weight the given diagrams with their relative probabilities.
   */
  virtual Selector<DiagramIndex> diagrams(const DiagramVector &) const;

  /**
   * Return the scale associated with the last set phase space point.
   */
  virtual Energy2 scale() const;


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

  /** 
   * Standard Interfaced virtual functions.
   */
  virtual void doinit() throw(InitException);

protected:

  /**
   * Constants for the different terms set from the StandardModel in
   * the init() function.
   */
  vector<double> coefs;

  /**
   * The mass squared of the \f$Z^0\f$.
   */
  Energy2 mZ2;

  /**
   * The width squared of the \f$Z^0\f$.
   */
  Energy2 GZ2;

  /**
   * The last continuum term to be used to select primary diagram.
   */
  mutable double lastCont;

  /**
   * The last Breit-Wigner term to be used to select primary diagram.
   */
  mutable double lastBW;

private:

  /**
   * Describe a class with persistent data.
   */
  static ClassDescription<MEqq2gZ2ll> initMEqq2gZ2ll;

  /** 
   * Private and non-existent assignment operator.
   */
  MEqq2gZ2ll & operator=(const MEqq2gZ2ll &);

};

}

namespace ThePEG {

/**
 * This template specialization informs ThePEG about the base class of
 * MEBase.
 */
template <>
struct BaseClassTrait<Herwig::MEqq2gZ2ll,1> {
  /** Typedef of the base class of MEqq2gZ2ll. */
  typedef ME2to2QCD NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * MEBase class.
 */
template <>
struct ClassTraits<Herwig::MEqq2gZ2ll>
  : public ClassTraitsBase<Herwig::MEqq2gZ2ll> {
  /** Return the class name. */
  static string className() { return "Herwig++::MEqq2gZ2ll"; }
  /** Return the name of the shared library be loaded to get
   *  access to the WeakPartonicDecayer class and every other class it uses
   *  (except the base class). */
  static string library() { return "HwME.so"; }
};

}

#include "MEqq2gZ2ll.icc"

#endif /* ThePEG_MEqq2gZ2ll_H */
