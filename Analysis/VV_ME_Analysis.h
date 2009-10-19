// -*- C++ -*-
#ifndef HERWIG_VV_ME_Analysis_H
#define HERWIG_VV_ME_Analysis_H
//
// This is the declaration of the VV_ME_Analysis class.
//
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Handlers/AnalysisHandler.h"
#include "Herwig++/Utilities/Histogram.h"

#include "fastjet/PseudoJet.hh"

namespace Herwig {
  
using namespace ThePEG;

/**
 * Here is the documentation of the VV_ME_Analysis class.
 *
 * @see \ref VV_ME_AnalysisInterfaces "The interfaces"
 * defined for VV_ME_Analysis.
 */
class VV_ME_Analysis: public AnalysisHandler {
    
public:

  /**
   * The default constructor.
   */
  VV_ME_Analysis();

  /** @name Virtual functions required by the AnalysisHandler class. */
  //@{
  /**
   * Analyze a given Event. Note that a fully generated event
   * may be presented several times, if it has been manipulated in
   * between. The default version of this function will call transform
   * to make a lorentz transformation of the whole event, then extract
   * all final state particles and call analyze(tPVector) of this
   * analysis object and those of all associated analysis objects. The
   * default version will not, however, do anything on events which
   * have not been fully generated, or have been manipulated in any
   * way.
   * @param event pointer to the Event to be analyzed.
   * @param ieve the event number.
   * @param loop the number of times this event has been presented.
   * If negative the event is now fully generated.
   * @param state a number different from zero if the event has been
   * manipulated in some way since it was last presented.
   */
  virtual void analyze(tEventPtr event, long ieve, int loop, int state);
    
  /**
   * Analyze the given vector of particles. The default version calls
   * analyze(tPPtr) for each of the particles.
   * @param particles the vector of pointers to particles to be analyzed
   */
  virtual void analyze(const tPVector & particles);
  
  /**
   * Analyze the given particle.
   * @param particle pointer to the particle to be analyzed.
   */
  virtual void analyze(tPPtr particle);
  
  //@}
    
public:

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

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  virtual void dofinish();

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class without persistent data.
   */
  static NoPIOClassDescription<VV_ME_Analysis> initVV_ME_Analysis;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  VV_ME_Analysis & operator=(const VV_ME_Analysis &);

private: 

  /**
   *   Polar angles of the leptons in their parent vector boson
   *   rest frame.
   */
  Histogram _cos3_h    , _cos4_h    , _cos5_h  , _cos6_h  ;

  /**
   *   Scalar sum of the \f$p_T\f$ of the leptons from the first
   *   vector boson.
   */
  Histogram _HT34_h    ;

  /**
   *   Pseudorapidities and \f$p_T\f$'s of the leptons from the first
   *   vector boson.
   */
  Histogram _eta3_h    , _pt3_h     , _eta4_h  , _pt4_h   ;

  /**
   *   Pseudorapidity, rapidity, \f$p_T\f$ and mass of  the first
   *   vector boson.
   */
  Histogram _eta34_h   , _y34_h     , _pt34_h  , _pt34_mcfm_h  , _m34wz_h ;

  /**
   *   Pseudorapidities and \f$p_T\f$'s of the leptons from the second
   *   vector boson.
   */
  Histogram _eta5_h    , _pt5_h     , _eta6_h  , _pt6_h   ;

  /**
   *   Pseudorapidity, rapidity, \f$p_T\f$ and mass of  the second
   *   vector boson.
   */
  Histogram _eta56_h   , _y56_h     , _pt56_h  , _pt56_mcfm_h  , _m56_h   ;

  /**
   *   Scalar sum of the \f$p_T\f$'s of all of the leptons, the
   *   rapidity Pseudorapidity, rapidity, \f$p_T\f$ and also the
   *   rapidity and mass of the diboson system.
   */
  Histogram _HT3456_h  , _y3456_h   , _m3456_h ;

  /**
   *   The \f$\theta_{1}\f$ born variable and the \f$\theta_{2}\f$
   *   radiative variable.
   */
  Histogram _th1_h     , _th2_h     ;

  /**
   *   \f$p_T\f$ and rapidity of the boson pair system
   */
  Histogram _ptVV_h    , _yVV_h     ;

  /**
   *   \f$p_T\f$ and rapidity of the hardest jet
   */
  Histogram _ptJet_h   ;
  Histogram _yJet_10_h     , _yJet_40_h     , _yJet_80_h    ;

  /**
   *   Rapidity gap between the hardest jet and the first vector boson
   */
  Histogram _yJet_yV1_10_h , _yJet_yV1_40_h , _yJet_yV1_80_h;

  /**
   *   Rapidity gap between the hardest jet and the second vector boson
   */
  Histogram _yJet_yV2_10_h , _yJet_yV2_40_h , _yJet_yV2_80_h;

  /**
   *   Rapidity gap between the hardest jet and the diboson system
   */
  Histogram _yJet_yVV_10_h , _yJet_yVV_40_h , _yJet_yVV_80_h;

  /**
   *   Jet multiplicities subject to \f$p_T\f$ cuts
   */
  Histogram _nJets_10_h    , _nJets_40_h    , _nJets_80_h   ;

public:
    
  /**
   * A pair of the incoming hadrons.
   */
  PPair inbound_;
  
  /**
   * A very small parameter to determine if rotations etc should be applied.
   */
  double epsilon_;
  
protected:
  
  /**
   * Function to return a specific boost to the rest frame of the 
   * leptons to histogram the polar angles in this frame. This is 
   * verbatim c++ translation of mcfm/src/Singletop/boostx.f.
   */
  Lorentz5Momentum boostx(Lorentz5Momentum p_in, Lorentz5Momentum pt,
			  Lorentz5Momentum ptt);
  
  /**
   * This function is stolen from fortran herwig, in order to maintain
   * consistency with the fortran analysis in mcfm. It takes a vector
   * p and returns the rotation matrix which when applied to p aligns
   * it along the z-axis and follows this with an azimuthal rotation 
   * of cp = cos(phi), sp = sin(phi).
   */
  LorentzRotation hwurot(Lorentz5Momentum p, double cp, double sp);
  
  // If needed, insert declarations of virtual function defined in the
  // InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).
  
};
  
}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {
  
/** @cond TRAITSPECIALIZATIONS */
  
/** This template specialization informs ThePEG about the
 *  base classes of VV_ME_Analysis. */
template <>
struct BaseClassTrait<Herwig::VV_ME_Analysis,1> {
  /** Typedef of the first base class of VV_ME_Analysis. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the VV_ME_Analysis class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::VV_ME_Analysis>
  : public ClassTraitsBase<Herwig::VV_ME_Analysis> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::VV_ME_Analysis"; }
  /**
   * The name of a file containing the dynamic library where the class
   * VV_ME_Analysis is implemented. It may also include several, space-separated,
   * libraries if the class VV_ME_Analysis depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "libfastjet.so HwVV_ME_Analysis.so"; }
};
  
/** @endcond */

}

#endif /* HERWIG_VV_ME_Analysis_H */

