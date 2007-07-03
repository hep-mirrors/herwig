#ifndef EVT_SHAPES_H
#define EVT_SHAPES_H

#include "ThePEG/CLHEPWrap/ThreeVector.h"
#include "ThePEG/EventRecord/Particle.h"

using namespace ThePEG;
using std::vector;

namespace Herwig {

/** \ingroup Analysis
 *
 *  The purpose of this class is to calculate all typical event shapes.
 *  Thrust is just obtained from a wrapper to EventShape (from hep.).
 *  The linear and quadratic momentum tensors are diagonalized with the
 *  help of CLHEP numerical methods.
 */ 
class EvtShapes {

public:

  /**
   * Constructor takes vector of transient particle pointers 
   * (= tPVector).  Only copies particles. 
   */
  EvtShapes(const tPVector & part);

  /**
   *  Destructor
   */
  ~EvtShapes();

  /**
   *  Single particle variables which do not depend on event shapes axes
   */
  //@{
  /**
   *  Ratio of momentum to beam momentum
   */
  double getX(const Lorentz5Momentum & p, const Energy & Ebeam);

  /**
   *  The scaled momentum \f$\xi=-\log\left( p/E_{\rm beam}\right)\f$.
   */
  double getXi(const Lorentz5Momentum & p, const Energy & Ebeam);

  /**
   *  Transverse momentum with respect to the beam
   */
  Energy getPt(const Lorentz5Momentum & p);

  /**
   *  Rapidity with respect to the beam direction
   */
  double getRapidity(const Lorentz5Momentum & p);
  //@}

  /**
   * Single particle variables related to one of the shape axis.
   */
  //@{
  /**
   *  Transverse momentum with respect to the thrust axis in the event plane
   */
  Energy ptInT(const Lorentz5Momentum & p);

  /**
   *  Transverse momentum with respect to the thrust axis out of the event plane
   */
  Energy ptOutT(const Lorentz5Momentum & p);

  /**
   *  Rapidity with respect to the thrust axis
   */
  double yT(const Lorentz5Momentum & p);

  /**
   *  Transverse momentum with respect to the sphericity axis in the event plane
   */
  Energy ptInS(const Lorentz5Momentum & p);

  /**
   *  Transverse momentum with respect to the sphericity axis out of the event plane
   */
  Energy ptOutS(const Lorentz5Momentum & p);

  /**
   *  Rapidity with respect to the sphericity axis
   */
  double yS(const Lorentz5Momentum & p);
  //@}

  /**
   *  Member functions to return thrust related shapes
   */
  //@{
  /**
   *  The thrust
   */
  double thrust(); 

  /**
   *  The major
   */ 
  double thrustMajor(); 

  /**
   *  The minor
   */ 
  double thrustMinor(); 

  /**
   *  The oblateness
   */ 
  double oblateness(); 

  /**
   *  The thrust axis
   */
  Axis thrustAxis(); 

  /**
   *  The major axis
   */ 
  Axis majorAxis(); 

  /**
   *  The minor axis
   */
  Axis minorAxis(); 
  //@}


  /**
   * Linear momentum tensor related event shapes
   */
  //@{
  /**
   *  The C parameter
   */
  double CParameter();

  /**
   *  The D parameter
   */
  double DParameter();

  /**
   *  The eigenvalues in descending order
   */
  vector<double> linTenEigenValues();

  /**
   *  The eigenvectors in order of descending eigenvalue
   */
  vector<Axis> linTenEigenVectors();
  //@}

  /**
   * Quadratic momentum tensor related variables
   */
  //@{
  /**
   *  The sphericity
   */
  double sphericity();

  /**
   *  The aplanarity
   */
  double aplanarity();

  /**
   *  The planarity
   */
  double planarity();

  /**
   *  The sphericity axis
   */
  Axis sphericityAxis();

  /**
   *  The sphericity eigenvalues
   */
  vector<double> sphericityEigenValues();

  /**
   *  The sphericity eigenvectors
   */
  vector<Axis> sphericityEigenVectors();
  //@}

  /**
   * Jet mass related event shapes
   */
  //@{
  /**
   *  The high hemishpere mass squared divided by the visible energy squared
   */
  double Mhigh2();

  /**
   *  The low hemishpere mass squared divided by the visible energy squared
   */
  double Mlow2();

  /**
   *  The difference between the 
   * hemishpere masses squared divided by the visible energy squared
   */
  double Mdiff2();
  //@}

  /**
   * Jet broadening related event shapes
   */
  //@{
  /**
   *  The wide jet broadening
   */
  double Bmax();

  /**
   *  The narrow jet broadening
   */
  double Bmin();

  /**
   *  The sum of the jet broadenings
   */
  double Bsum();

  /**
   *  The difference of the jet broadenings
   */
  double Bdiff();
  //@}

  /**
   *  Energy-Energy correlation related members
   */
  //@{
  /**
   * Energy-energy correlation (EEC)
   * @param hi is the histogram and has to be provided externally since this
   * class is supposed to die after one event.  It is understood that
   * the range of the histogam is -1 < cos(chi) < 1. 
   * hi.front() contains the bin [-1 < cos(chi) < -1+delta] and
   * hi.back() the bin [1-delta < cos(chi) < 1].  delta =
   * 2/hi.size().  We use classical indices to access the vector. 
   */
  void bookEEC(vector<double> & hi);

  /**
   * Before writing the histogram it has to be normalized according to the
   * number of events. 
   */
  void normalizeEEC(vector<double> & hi, long evts);

  /**
   * The asymmetry of EEC is calculated from a given \f$\cos\chi\f$ and EEC
   * histogram, which is a vector<double> as described above.
   */
  double AEEC(vector<double> & hi, double& coschi);
  //@}

private: 

  /**
   *  Default constructor should not be used and is therefore private
   */
  EvtShapes();

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
   *  Check if the quadratic tensor related variables have been calculated
   *  and if not do so
   */
  inline void checkSphericity();

  /**
   *  Check if the hemisphere mass variables have been calculated and if not do so
   */
  inline void checkHemispheres();

  /**
   *  Check if the jet broadenings have been calculated and if not do so
   */
  inline void checkBroadening();
  //@}

  /**
   *  Methods that actually calculate the event shapes
   */
  //@{
  /**
   *  Calculate the hemisphere masses
   */
  void calcHemisphereMasses();

  /**
   *  Calculate the jet broadenings
   */
  void calcBroadening();

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
   *  Member to calculate the thrust
   * @param p The three vectors
   * @param t The thrust
   * @param taxis The thrust axis
   */
  void calcT(const vector<Momentum3 > &p,
	     Energy2 &t, Axis & taxis);
  /**
   *  Member to calculate the major
   * @param p The three vectors
   * @param m The major
   * @param maxis The major axis
   */
  void calcM(const vector<Momentum3 > &p, 
	     Energy2 &m, Axis & maxis);
  //@}


private:
 
  /**
   *  Vector of particles to be analysed
   */
  tPVector _pv; 

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

  /**
   * Whether or not the jet broadenings have been calculated
   */
  bool _broadDone; 
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
#include "EvtShapes.icc"
#endif
