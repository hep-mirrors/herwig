// -*- C++ -*-

// (C) 2007-2009 Simon Plaetzer -- sp@particle.uni-karlsruhe.de
// Copyright (C) 2002-2007 The Herwig Collaboration

#ifndef Analysis2_EventShapes2_H
#define Analysis2_EventShapes2_H
//
// This is the declaration of the EventShapes2 class.
//

#include "Analysis2Base.h"
#include "ThePEG/Vectors/Lorentz5Vector.h"
#include "ThePEG/Vectors/ThreeVector.h"
#include "ThePEG/EventRecord/Particle.h"
#include "EventShapes2.fh"

namespace Analysis2 {

using namespace ThePEG;

/**\ingroup Analysis2
 * 
 * Class to analyse event shapes
 *
 * @see \ref EventShapes2Interfaces "The interfaces"
 * defined for EventShapes2.
 */
class EventShapes2: public Analysis2Base {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline EventShapes2();

  /**
   * The destructor.
   */
  virtual ~EventShapes2();
  //@}

public:

  /**
   * Analyze the given vector of particles. The default version calls
   * analyze(tPPtr) for each of the particles.
   * @param particles the vector of pointers to particles to be analyzed
   */
  virtual void analyze(const tPVector & particles);

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

  /**
   *  Member to reset the particles to be considered
   */
  inline void reset(const vector<Lorentz5Momentum> &part);

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
   *  Single particle variables which do not depend on event shapes axes
   */
  //@{
  /**
   *  Ratio of momentum to beam momentum
   */
  inline double getX(const Lorentz5Momentum & p, const Energy & Ebeam);

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


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


protected:

  /** @name Standard Interfaced functions. */
  //@{

  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit() throw(InitException);

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  virtual void dofinish();

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

private:

  /**
   * Output option strings for \f$1-T\f$ distribution.
   */
  string _omthr;

  /**
   * Output option strings for the major distribution
   */
  string _maj;

  /**
   * Output option strings for the minor distribution
   */
  string _min;

  /**
   *Output option strings for the oblateness distribution
   */
  string _obl;

  /**
   * Output option strings for the sphericity distribution
   */
  string _sph;

  /**
   * Output option strings for the aplanarity distribution
   */
  string _apl;

  /**
   * Output option strings for the planarity distribution
   */
  string _pla;

  /**
   * Output option strings for the C distribution
   */
  string _c;

  /**
   * Output option strings for the D distribution
   */
  string _d;

  /**
   * Output option strings for the \f$M_{\rm high}\f$ distribution
   */
  string _mhi;

  /**
   * Output option strings for the \f$M_{\rm low}\f$ distribution
   */
  string _mlo;

  /**
   * Output option strings for the \f$M_{\rm high}-M_{\rm low}\f$ distribution
   */
  string _mdiff;

  /**
   * Output option strings for the \f$B_{\rm max}\f$ distribution
   */
  string _bmax;

  /**
   * Output option strings for the \f$B_{\rm min}\f$ distribution
   */
  string _bmin;

  /**
   * Output option strings for the \f$B_{\rm max}+B_{\rm min}\f$ distribution
   */
  string _bsum;

  /**
   * Output option strings for the \f$B_{\rm max}-B_{\rm min}\f$ distribution
   */
  string _bdiff;

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<EventShapes2> initEventShapes2;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  EventShapes2 & operator=(const EventShapes2 &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of EventShapes2. */
template <>
struct BaseClassTrait<Analysis2::EventShapes2,1> {
  /** Typedef of the first base class of EventShapes2. */
  typedef Analysis2::Analysis2Base NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the EventShapes2 class and the shared object where it is defined. */
template <>
struct ClassTraits<Analysis2::EventShapes2>
  : public ClassTraitsBase<Analysis2::EventShapes2> {
  /** Return a platform-independent class name */
  static string className() { return "Analysis2::EventShapes2"; }
  /**
   * The name of a file containing the dynamic library where the class
   * EventShapes2 is implemented. It may also include several, space-separated,
   * libraries if the class EventShapes2 depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "Analysis2.so"; }
};

/** @endcond */

}

#include "EventShapes2.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "EventShapes2.tcc"
#endif

#endif /* Analysis2_EventShapes2_H */
