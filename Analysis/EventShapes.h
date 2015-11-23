// -*- C++ -*-
//
// EventShapes.h is a part of Herwig - A multi-purpose Monte Carlo
// event generator Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for
// details.  Please respect the MCnet academic guidelines, see
// GUIDELINES for details.
//
#ifndef HERWIG_EventShapes_H
#define HERWIG_EventShapes_H
//
// This is the declaration of the EventShapes class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "ThePEG/Handlers/AnalysisHandler.h"
#include "ThePEG/Vectors/Lorentz5Vector.h"
#include "ThePEG/Vectors/ThreeVector.h"
#include "ThePEG/EventRecord/Particle.h"
#include "EventShapes.fh"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Analysis
 *
 * The EventShapes class is designed so that certain event shapes, such
 * as the thrust are only calculated once per event given the speed of
 * the calculation.
 *
 * @see \ref EventShapesInterfaces "The interfaces" defined for
 * EventShapes.
 */
class EventShapes: public Interfaced {

public:

  /**
   * The default constructor.
   */
  EventShapes() : _thrustDone(false), _spherDone(false), _linTenDone(false),
		  _hemDone(false), _useCmBoost(false), 
		  _mPlus(), _mMinus(), _bPlus(), _bMinus() 
  {}

  /**
   *  Member to reset the particles to be considered
   */
  void reset(const tPVector &part) {
    _pv.resize(part.size());
    for(unsigned int ix=0;ix<part.size();++ix) _pv[ix]=part[ix]->momentum();
    _thrustDone = false;
    _spherDone  = false;
    _linTenDone = false;
    _hemDone    = false;
    _useCmBoost = false; 
  }


public:

  /**
   *  Member functions to return thrust related shapes
   */
  //@{
  /**
   *  The thrust
   */
  double thrust() {
    checkThrust(); 
    return _thrust[0];
  }

  /**
   *  The major
   */ 
  double thrustMajor() {
    checkThrust(); 
    return _thrust[1];
  }

  /**
   *  The minor
   */ 
  double thrustMinor() {
    checkThrust(); 
    return _thrust[2];
  }

  /**
   *  The oblateness
   */ 
  double oblateness() {
    checkThrust(); 
    return _thrust[1]-_thrust[2];
  }

  /**
   *  The thrust axis
   */
  Axis thrustAxis() {
    checkThrust(); 
    return _thrustAxis[0];
  }

  /**
   *  The major axis
   */ 
  Axis majorAxis() {
    checkThrust(); 
    return _thrustAxis[1];
  }

  /**
   *  The minor axis
   */
  Axis minorAxis() {
    checkThrust(); 
    return _thrustAxis[2];
  }
  //@}

  /**
   * Linear momentum tensor related event shapes
   */
  //@{
  /**
   *  The C parameter
   */
  double CParameter() {
    checkLinTen(); 
    return 3.*(_linTen[0]*_linTen[1]+_linTen[1]*_linTen[2]
	       +_linTen[2]*_linTen[0]); 
  }

  /**
   *  The D parameter
   */
  double DParameter() {
    checkLinTen(); 
    return 27.*(_linTen[0]*_linTen[1]*_linTen[2]); 
  }

  /**
   *  The eigenvalues in descending order
   */
  vector<double> linTenEigenValues() {
    checkLinTen(); 
    return _linTen; 
  }


  /**
   *  The eigenvectors in order of descending eigenvalue
   */
  vector<Axis> linTenEigenVectors() {
    checkLinTen(); 
    return _linTenAxis; 
  }

  //@}

  /**
   * Quadratic momentum tensor related variables
   */
  //@{
  /**
   *  The sphericity
   */
  double sphericity() {
    checkSphericity(); 
    return 3./2.*(_spher[1]+_spher[2]); 
  }

  /**
   *  The aplanarity
   */
  double aplanarity() {
    checkSphericity(); 
    return 3./2.*_spher[2];
  }


  /**
   *  The planarity
   */
  double planarity() {
    checkSphericity(); 
    return _spher[1]-_spher[2]; 
  }

  /**
   *  The sphericity axis
   */
  Axis sphericityAxis() {
    checkSphericity(); 
    return _spherAxis[0]; 
  }


  /**
   *  The sphericity eigenvalues
   */
  vector<double> sphericityEigenValues() {
    checkSphericity(); 
    return _spher; 
  }

  /**
   *  The sphericity eigenvectors
   */
  vector<Axis> sphericityEigenVectors() {
    checkSphericity(); 
    return _spherAxis; 
  }  //@}

  /**
   * Jet mass related event shapes
   */
  //@{
  /**
   *  The high hemishpere mass squared divided by the visible energy
   *  squared
   */
  double Mhigh2() {
    checkHemispheres();
    return _mPlus; 
  } 
  
  /**
   *  The low hemishpere mass squared divided by the visible energy
   *  squared
   */
  double Mlow2() {
    checkHemispheres();
    return _mMinus; 
  } 

  /**
   *  The difference between the 
   * hemishpere masses squared divided by the visible energy squared
   */
  double Mdiff2() {
    checkHemispheres();
    return _mPlus-_mMinus; 
  } 

  //@}

  /**
   * Jet broadening related event shapes
   */
  //@{
  /**
   *  The wide jet broadening
   */
  double Bmax() {
    checkHemispheres(); 
    return _bPlus;
  }

  /**
   *  The narrow jet broadening
   */
  double Bmin() {
    checkHemispheres(); 
    return _bMinus;
  }

  /**
   *  The sum of the jet broadenings
   */
  double Bsum() {
    checkHemispheres(); 
    return _bPlus+_bMinus;
  }


  /**
   *  The difference of the jet broadenings
   */
  double Bdiff() {
    checkHemispheres(); 
    return _bPlus-_bMinus;
  }
  //@}

  /**
   *  Single particle variables which do not depend on event shapes axes
   */
  //@{

  /**
   *  The scaled momentum \f$\xi=-\log\left( p/E_{\rm beam}\right)\f$.
   */
  double getXi(const Lorentz5Momentum & p, 
				   const Energy & Ebeam) {
    return((Ebeam > 0*MeV && p.vect().mag() > 0*MeV) ? 
	   log(Ebeam/p.vect().mag()) : -1.); 
  }

  /**
   *  Transverse momentum with respect to the beam
   */
  Energy getPt(const Lorentz5Momentum & p) {
    return p.perp(); 
  }

  /**
   *  Rapidity with respect to the beam direction
   */
  double getRapidity(const Lorentz5Momentum & p) {
    return (p.t() > p.z() ? p.rapidity() : 1e99); 
  }
  //@}

  /**
   * Single particle variables related to one of the shape axis.
   */
  //@{
  /**
   *  Transverse momentum with respect to the thrust axis in the event plane
   */
  Energy ptInT(const Lorentz5Momentum & p) {
    checkThrust(); 
    return p.vect()*_thrustAxis[1]; 
  }
  
  /**
   *  Transverse momentum with respect to the thrust axis out of the
   *  event plane
   */
  Energy ptOutT(const Lorentz5Momentum & p) {
    checkThrust(); 
    return p.vect()*_thrustAxis[2]; 
  }

  /**
   *  Rapidity with respect to the thrust axis
   */
  double yT(const Lorentz5Momentum & p) {
    checkThrust(); 
    return (p.t() > p.vect()*_thrustAxis[0] ? 
	    p.rapidity(_thrustAxis[0]) : 1e99);
  }

  /**
   *  Transverse momentum with respect to the sphericity axis in the
   *  event plane
   */
  Energy ptInS(const Lorentz5Momentum & p) { 
    checkSphericity(); 
    return p.vect()*_spherAxis[1]; 
  }

  /**
   *  Transverse momentum with respect to the sphericity axis out of the
   *  event plane
   */
  Energy ptOutS(const Lorentz5Momentum & p) {
    checkSphericity(); 
    return p.vect()*_spherAxis[2]; 
  }

  /**
   *  Rapidity with respect to the sphericity axis
   */
  double yS(const Lorentz5Momentum & p) {
    checkSphericity(); 
    return (p.t() > p.vect()*_spherAxis[0] ? 
	    p.rapidity(_spherAxis[0]) : 1e99);
  }
  //@}


  /**
   * Energy-energy correlation (EEC) @param hi is the histogram and has
   * to be provided externally It is understood that the range of the
   * histogam is -1 < cos(chi) < 1.  hi.front() contains the bin [-1 <
   * cos(chi) < -1+delta] and hi.back() the bin [1-delta < cos(chi) <
   * 1].  delta = 2/hi.size(). We use classical indices to access the
   * vector.
   */
  void bookEEC(vector<double> & hi);

  /**
   * Before writing the histogram it has to be normalized according to
   * the number of events.
   */
  void normalizeEEC(vector<double> & hi, long evts) {
    for (unsigned int bin = 0; bin < hi.size(); bin++) bin /= (hi.size()*evts);
  }
  
  /**
   * The asymmetry of EEC is calculated from a given \f$\cos\chi\f$ and
   * EEC histogram, which is a vector<double> as described above.
   */
  double AEEC(vector<double> & hi, double& coschi) {
    if (coschi > 0. && coschi <= 1.) {
      int i = static_cast<int>( floor((-coschi+1.)/2.*hi.size()) ); 
      int j = static_cast<int>( floor(( coschi+1.)/2.*hi.size()) ); 
      return hi[i]-hi[j];
    } else {
      return 1e99;
    }
  }

public:

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or when this class is dynamically
   * loaded.
   */
  static void Init();

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.  @return a pointer to the new
   * object.
   */
  virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.  @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const {return new_ptr(*this);}
  //@}

private:

  /**
   *  Check whether the initialization of a certain class of event shapes
   *  has been calculated and if not do so
   */
  //@{
  /**
   *  Check if thrust related variables have been calculated and if not
   *  do so
   */
  void checkThrust() {
    if (!_thrustDone) {
      _thrustDone = true;
      calculateThrust(); 
    }
  }

  /**
   *  Check if the linear tensor related variables have been calculated
   *  and if not do so
   */
  void checkLinTen() {
    if (!_linTenDone) {
      _linTenDone = true;
      diagonalizeTensors(true, _useCmBoost); 
    }
  }

  /**
   *  Check if the quadratic tensor related variables have been
   *  calculated and if not do so
   */
  void checkSphericity() {
    if (!_spherDone) {
      _spherDone = true;
      diagonalizeTensors(false, _useCmBoost); 
    }
  }

  /**
   *  Check if the hemisphere mass variables and jet broadenings have
   *  been calculated and if not do so
   */
  void checkHemispheres() {
    if (!_hemDone) {
      _hemDone = true;
      calcHemisphereMasses(); 
    }
  }
  //@}

  /**
   *  Methods that actually calculate the event shapes
   */
  //@{
  /**
   *  Calculate the hemisphere masses and jet broadenings
   */
  void calcHemisphereMasses();

  /**
   * Calculate the thrust and related axes
   */
  void calculateThrust();

  /**
   * Diagonalize the tensors @param linear switch between
   * diagonalization of linear/quadratic tensor.  @param cmboost tells
   * whether to boost into cm frame of all momenta first, or not
   * (default off, and no interface to this).
   */
  void diagonalizeTensors(bool linear, bool cmboost);

  /**
   * Quite general diagonalization of a symmetric Matrix T, given as an
   * array of doubles.  The symmetry is not checked explicitly as this
   * is clear in the context.  It uses an explicit generic solution of
   * the eigenvalue problem and no numerical approximation, based on
   * Cardano's formula.  @param T Matrix to be diagonalised
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
  static NoPIOClassDescription<EventShapes> initEventShapes;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  EventShapes & operator=(const EventShapes &);

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
 *  base classes of EventShapes. */
template <>
struct BaseClassTrait<Herwig::EventShapes,1> {
  /** Typedef of the first base class of EventShapes. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the EventShapes class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::EventShapes>
  : public ClassTraitsBase<Herwig::EventShapes> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::EventShapes"; }
  /** Return the name(s) of the shared library (or libraries) be loaded to get
   *  access to the EventShapes class and any other class on which it depends
   *  (except the base class). */
  static string library() { return "HwAnalysis.so"; }
};

/** @endcond */

}

#endif /* HERWIG_EventShapes_H */
