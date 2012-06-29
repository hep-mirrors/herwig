// -*- C++ -*-
//
// MyEventShapes.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MyEventShapes_H
#define HERWIG_MyEventShapes_H
//
// This is the declaration of the MyEventShapes class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "ThePEG/Handlers/AnalysisHandler.h"
#include "ThePEG/Vectors/Lorentz5Vector.h"
#include "ThePEG/Vectors/ThreeVector.h"
#include "ThePEG/EventRecord/Particle.h"
#include "MyEventShapes.fh"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Analysis
 *
 * The MyEventShapes class is designed so that certain event shapes,
 * such as the thrust are only calculated once per event given the
 * speed of the calculation.
 *
 * @see \ref MyEventShapesInterfaces "The interfaces"
 * defined for MyEventShapes.
 */
class MyEventShapes: public Interfaced {

  //public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  // inline MyEventShapes();

  /**
   * The copy constructor.
   */
  // inline MyEventShapes(const MyEventShapes &);

  /**
   * The destructor.
   */
  //  virtual ~MyEventShapes();
  //@}

public:

  /**
   *  Member to reset the particles to be considered
   */
  inline void reset(const tPVector &part);

public:

  /**
   *  Member functions to return thrust related shapes
   */
  //@{
  /**
   *  The thrust
   */
  inline double thrust();

  /**
   *  The major
   */ 
  inline double thrustMajor();

  /**
   *  The minor
   */ 
  inline double thrustMinor();

  /**
   *  The oblateness
   */ 
  inline double oblateness(); 

  /**
   *  The thrust axis
   */
  inline Axis thrustAxis();

  /**
   *  The major axis
   */ 
  inline Axis majorAxis(); 

  /**
   *  The minor axis
   */
  inline Axis minorAxis(); 
  //@}

  /**
   * Linear momentum tensor related event shapes
   */
  //@{
  /**
   *  The C parameter
   */
  inline double CParameter();

  /**
   *  The D parameter
   */
  inline double DParameter();

  /**
   *  The eigenvalues in descending order
   */
  inline vector<double> linTenEigenValues();

  /**
   *  The eigenvectors in order of descending eigenvalue
   */
  inline vector<Axis> linTenEigenVectors();
  //@}

  /**
   * Quadratic momentum tensor related variables
   */
  //@{
  /**
   *  The sphericity
   */
  inline double sphericity();

  /**
   *  The aplanarity
   */
  inline double aplanarity();

  /**
   *  The planarity
   */
  inline double planarity();

  /**
   *  The sphericity axis
   */
  inline Axis sphericityAxis();

  /**
   *  The sphericity eigenvalues
   */
  inline vector<double> sphericityEigenValues();

  /**
   *  The sphericity eigenvectors
   */
  inline vector<Axis> sphericityEigenVectors();
  //@}

  /**
   * Jet mass related event shapes
   */
  //@{
  /**
   *  The high hemishpere mass squared divided by the visible energy squared
   */
  inline double Mhigh2();

  /**
   *  The low hemishpere mass squared divided by the visible energy squared
   */
  inline double Mlow2();

  /**
   *  The difference between the 
   * hemishpere masses squared divided by the visible energy squared
   */
  inline double Mdiff2();
  //@}

  /**
   * Jet broadening related event shapes
   */
  //@{
  /**
   *  The wide jet broadening
   */
  inline double Bmax();

  /**
   *  The narrow jet broadening
   */
  inline double Bmin();

  /**
   *  The sum of the jet broadenings
   */
  inline double Bsum();

  /**
   *  The difference of the jet broadenings
   */
  inline double Bdiff();
  //@}

  
  /**
   *  The scaled momentum \f$\xi=-\log\left( p/E_{\rm beam}\right)\f$.
   */
  inline double getXi(const Lorentz5Momentum & p, const Energy & Ebeam);

  /**
   *  Transverse momentum with respect to the beam
   */
  inline Energy getPt(const Lorentz5Momentum & p);

  /**
   *  Rapidity with respect to the beam direction
   */
  inline double getRapidity(const Lorentz5Momentum & p);
  //@}

  /**
   * Single particle variables related to one of the shape axis.
   */
  //@{
  /**
   *  Transverse momentum with respect to the thrust axis in the event plane
   */
  inline Energy ptInT(const Lorentz5Momentum & p);

  /**
   *  Transverse momentum with respect to the thrust axis out of the event plane
   */
  inline Energy ptOutT(const Lorentz5Momentum & p);

  /**
   *  Rapidity with respect to the thrust axis
   */
  inline double yT(const Lorentz5Momentum & p);

  /**
   *  Transverse momentum with respect to the sphericity axis in the event plane
   */
  inline Energy ptInS(const Lorentz5Momentum & p);

  /**
   *  Transverse momentum with respect to the sphericity axis out of the event plane
   */
  inline Energy ptOutS(const Lorentz5Momentum & p);

  /**
   *  Rapidity with respect to the sphericity axis
   */
  inline double yS(const Lorentz5Momentum & p);
  //@}


  /**
   * Energy-energy correlation (EEC)
   * @param hi is the histogram and has to be provided externally
   * It is understood that
   * the range of the histogam is -1 < cos(chi) < 1. 
   * hi.front() contains the bin [-1 < cos(chi) < -1+delta] and
   * hi.back() the bin [1-delta < cos(chi) < 1].  delta =
   * 2/hi.size(). We use classical indices to access the vector. 
   */
  void bookEEC(vector<double> & hi);

  /**
   * Before writing the histogram it has to be normalized according to the
   * number of events. 
   */
  inline void normalizeEEC(vector<double> & hi, long evts);

  /**
   * The asymmetry of EEC is calculated from a given \f$\cos\chi\f$ and EEC
   * histogram, which is a vector<double> as described above.
   */
  inline double AEEC(vector<double> & hi, double& coschi);

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
   *  Check whether the initialization of a certain class of event shapes
   *  has been calculated and if not do so
   */
  //@{
  /**
   *  Check if thrust related variables have been calculated and if not do so
   */
  inline void checkThrust();

  /**
   *  Check if the linear tensor related variables have been calculated and if not do so
   */
  inline void checkLinTen();

  /**
   *  Check if the quadratic tensor related variables have been calculated and if not do so
   */
  inline void checkSphericity();

  /**
   *  Check if the hemisphere mass variables and jet broadenings 
   *  have been calculated and if not do so
   */
  inline void checkHemispheres();
  //@}

  /**
   *  Methods that actually calculate the event shapes
   */
  //@{
  /**
   *  Calculate the hemisphere masses and jet broadenings
   */
    inline void calcHemisphereMasses();

  /**
   * Calculate the thrust and related axes
   */
  void calculateThrust();

  /**
   * Diagonalize the tensors
   * @param linear switch between diagonalization of linear/quadratic tensor.
   * @param cmboost tells whether to boost into cm frame of all
   * momenta first, or not (default off, and no interface to this).
   */
  void diagonalizeTensors(bool linear, bool cmboost);

  /**
   * Quite general diagonalization of a symmetric Matrix  T, given as
   * an array of doubles.  The symmetry is not checked explicitly as
   * this is clear in the context.  It uses an explicit generic
   * solution of the eigenvalue problem and no numerical
   * approximation, based on Cardano's formula.
   * @param T Matrix to be diagonalised 
   */
  vector<double> eigenvalues(const double T[3][3]);

  /**
   * The eigenvector of @param T to a given eigenvalue @param lam 
   */
  Axis eigenvector(const double T[3][3], const double &lam);

  /**
   * The eigenvectors of @param T corresponding to the eigenvectors
   * @param lam . The ordering of the vectors corresponds to the
   * ordering of the eigenvalues.
   */
  vector<Axis> eigenvectors(const double T[3][3], const vector<double> &lam);

  /**
   *  Member to calculate the thrust
   * @param p The three vectors
   * @param t The thrust-squared (up to an Energy scale factor)
   * @param taxis The thrust axis
   */
  void calcT(const vector<Momentum3> &p, Energy2 &t, Axis &taxis);

  /**
   *  Member to calculate the major
   * @param p The three vectors
   * @param m The major-squared (up to an Energy scale factor)
   * @param maxis The major axis
   */
  void calcM(const vector<Momentum3> &p, Energy2 &m, Axis &maxis);
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static NoPIOClassDescription<MyEventShapes> initMyEventShapes;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MyEventShapes & operator=(const MyEventShapes &);

private:

  /**
   *  Vector of particle momenta to be analysed
   */
  vector<Lorentz5Momentum> _pv; 
  
  /**
   *  Various event shape axes
   */
  //@{
  /**
   *  The thrust related axes
   */
  vector<Axis> _thrustAxis;

  /**
   *  The sphericity related axes
   */
  vector<Axis> _spherAxis; 

  /**
   *  The linearised tensor axes
   */
  vector<Axis> _linTenAxis; 
  //@}

  /**
   *  Values of axis related event shapes
   */
  //@{
  /**
   *  Values of thrust related variables
   */
  vector<double> _thrust;

  /**
   *  Values of sphericity related variables
   */
  vector<double> _spher;

  /**
   *  Values of linearized tensor related variables
   */
  vector<double> _linTen;
  //@} 

  /**
   *  Whether or not certain event axes have been calculated
   */
  //@{
  /**
   *  Whether or not the thrust is calculated
   */
  bool _thrustDone;

  /**
   *  Whether or not the sphericity is calculated
   */
  bool _spherDone;

  /**
   *  Whether or not the linearizes tensor is calculated 
   */
  bool _linTenDone;

  /**
   *  Whether or not the hemisphere masses have been calculated
   */
  bool _hemDone; 
  //@}

  /**
   *  Whether ot not to boost to the CMS frame for the tensor diagonalizations
   */
  bool _useCmBoost;

  /**
   *  Hemisphere masses
   */
  //@{
  /**
   *  The high hemisphere mass
   */
  double _mPlus;

  /**
   *  The low hemisphere mass
   */
  double _mMinus;
  //@}

  /**
   *  The jet broadenings
   */
  //@{
  /**
   *  The wide jet broadening
   */
  double _bPlus;

  /**
   *  The narrow jet broadening
   */
  double _bMinus; 
  //@}
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MyEventShapes. */
template <>
struct BaseClassTrait<Herwig::MyEventShapes,1> {
  /** Typedef of the first base class of MyEventShapes. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MyEventShapes class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MyEventShapes>
  : public ClassTraitsBase<Herwig::MyEventShapes> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MyEventShapes"; }
  /** Return the name(s) of the shared library (or libraries) be loaded to get
   *  access to the MyEventShapes class and any other class on which it depends
   *  (except the base class). */
   static string library() { return "MyEventShapes.so"; }
  //   static string library() { return "HwAnalysis.so"; }
};

/** @endcond */

}

#include "MyEventShapes.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "MyEventShapes.tcc"
#endif

#endif /* HERWIG_MyEventShapes_H */
