#ifndef _EVT_SHAPES_H
#define _EVT_SHAPES_H

#include "ThePEG/CLHEPWrap/Matrix.h"
#include "ThePEG/CLHEPWrap/ThreeVector.h"
#include "ThePEG/EventRecord/Particle.h"

using namespace ThePEG;

namespace Herwig {

/** \ingroup Interfaces
 *
 *  The purpose of this class is to calculate all typical event shapes.
 *  Thrust is just obtained from a wrapper to EventShape (from hep.).
 *  The linear and quadratic momentum tensors are diagonalized with the
 *  help of CLHEP numerical methods.
 */ 
class EvtShapes {

public:

  /**
   * Constructor takes vector of transient Pythia7 particle pointers 
   * (= tPVector).  Only copies particles. 
   */
  EvtShapes(const tPVector & part);
  ~EvtShapes();

  /**
   * For convenience get just single particle variables here.
   */
  double getX(const Lorentz5Momentum & p, const Energy & Ebeam);
  double getXi(const Lorentz5Momentum & p, const Energy & Ebeam);
  Energy getPt(const Lorentz5Momentum & p);
  Energy getRapidity(const Lorentz5Momentum & p);

  /**
   * Single particle variables related to one of the shape axis.
   */
  Energy ptInT(const Lorentz5Momentum & p);
  Energy ptOutT(const Lorentz5Momentum & p);
  double yT(const Lorentz5Momentum & p);
  Energy ptInS(const Lorentz5Momentum & p);
  Energy ptOutS(const Lorentz5Momentum & p);
  double yS(const Lorentz5Momentum & p);

  /**
   * Get thrust-axis related shapes. 
   */
  double thrust(); 
  double thrustMajor(); 
  double thrustMinor(); 
  double oblateness(); 
  Vector3 thrustAxis(); 
  Vector3 majorAxis(); 
  Vector3 minorAxis(); 

  /** 
   * Linear momentum tensor related
   * get eigenvalues and vectors of the linear momentum tensor.  
   * They are sorted: linTenEigenValue[0] > ...[1] > ...[2].  
   * Yes, eigenvector[i] really does correspond to eigenvalue[i]!
   */
  vector<double> linTenEigenValues();
  vector<Vector3> linTenEigenVectors();
  double CParameter();
  double DParameter();
  
  /**
   * Quadratic momentum tensor related.
   */
  double sphericity();
  double aplanarity();
  double planarity();
  Vector3 sphericityAxis();
  vector<double> sphericityEigenValues();
  vector<Vector3> sphericityEigenVectors();

  /**
   * Jet mass related, high, low hemisphere masses squared divided 
   * by visible energy squared and their difference.
   */
  double Mhigh2();
  double Mlow2();
  double Mdiff2();
  
  /**
   * Jet broadening. 
   */
  double Bmax();
  double Bmin();
  double Bsum();
  double Bdiff();

  /**
   * Energy-energy correlation (EEC)
   * hi is the histogram and has to be provided externally since this
   * class is supposed to die after one event.  It is understood that
   * the range of the histogam is -1 < cos(chi) < 1. 
   * hi.front() contains the bin [-1 < cos(chi) < -1+delta] and
   * hi.back() the bin [1-delta < cos(chi) < 1].  delta =
   * 2/hi.size().  We use classical indices to access the vector. 
   */
  void bookEEC(vector<double> & hi);

  /**
   * Before writing the histogram it has to be normalized acc to the
   * number of events. 
   */
  void normalizeEEC(vector<double> & hi, long evts);

  /**
   * The asymmetry of EEC is calculated from a given cos(chi) and EEC
   * histogram , which is a vector<double> as described above.
   */
  double AEEC(vector<double> & hi, double& coschi);

private: 

  EvtShapes();

  /**
   * Check whether initialization has been done and if not do so.
   */
  inline void checkThrust();
  inline void checkLinTen();
  inline void checkSphericity();
  inline void checkHemispheres();
  inline void checkBroadening();

  /**
   * Actually do something.
   */
  void calcHemisphereMasses();
  void calcBroadening();

  /**
   * 'linear' to switch between diagonalization of linear/quadratic
   * tensor. 'cmboost' tells whether to boost into cm frame of all
   * momenta first, or not (default off, and no interface to this).
   */
  void diagonalizeTensors(bool linear, bool cmboost);

  /**
   * Doing the work.
   */
  void calcT(const vector<Vector3> &p, double &t, Vector3 &taxis);
  void calcM(const vector<Vector3> &p, double &m, Vector3 &maxis);

  /**
   * Fills/gets values and calculates minor.
   */
  void calculateThrust();

private:
 
  tPVector _pv; 
  std::vector<Vector3> _thrustAxis, _spherAxis, _linTenAxis; 
  std::vector<double> _thrust, _spher, _linTen; 
  bool _linTenDone, _spherDone, _thrustDone;
  double _mPlus, _mMinus, _bPlus, _bMinus; 
  bool _hemDone, _broadDone; 
  bool _useCmBoost;

};

}
#endif
