// -*- C++ -*-
#ifndef HERWIG_PGSInterface_H
#define HERWIG_PGSInterface_H
//
// This is the declaration of the PGSInterface class.
//

#include "ThePEG/Handlers/AnalysisHandler.h"
#include "ThePEG/Vectors/Lorentz5Vector.h"
#include "HepMC/GenEvent.h"
#include "HepMC/IO_HEPEVT.h"

namespace Herwig {
using namespace ThePEG;

/**
 *  Enumeration for the different types of reconstructed objects
 */
enum ObjectType {
  Photon=0,    /**< Photons               */
  Electron,    /**< Electrons             */
  Muon,        /**< Muons                 */
  Tau,         /**< Taus                  */
  Jet,         /**< Jets                  */
  HeavyCharged /**< Heavy Charged Objects */
};

/**
 *  Enumeration for the type of b-tag
 */
enum bTag {
  None=0, /**< No b-tag    */
  Loose,  /**< Loose b-tag */
  Tight,  /**< Tight b-tag */
  Both ,  /**< Both loose and tight b-tag */
};

/**
 *  Struct to store the properties of the reconstructed objects
 */
struct ReconstructedObject {

  /**
   *  Type of object
   */
  ObjectType type;

  /**
   *  Momentum of the object
   */
  Lorentz5Momentum momentum;

  /**
   *  Charge of the particle
   */
  double charge;

  /**
   *  Electromagnetic energy
   */
  Energy emenergy;

  /**
   * Hadronic energy
   */
  Energy hadenergy;

  /**
   *  Track energy
   */
  Energy trackenergy;

  /**
   *  Number of tracks
   */
  int numtracks;

  /**
   *  PDG Code for the particle
   */
  int PDGcode;

  /**
   *  Whether it is a loose or tight b-tag
   */
  bTag btagging;

  /**
   *  Transverse energy
   */
  Energy ET;

  /**
   *  Transverse energy in the isolation cone
   */
  Energy isolationET;

  /**
   *  Transverse momentum in the isolation cone
   */
  Energy isolationpT;

  /**
   *  Ratio of hadronic to electromagentic energy
   */
  double hadronicem;

  /**
   *  Ratio of electromagnetic energy to track momentum
   */
  double ep;

  /**
   *  track isolation energy
   */
  Energy trkisoEnergy;

  /**
   *  Number of pi0 in cone for tau
   */
  int npi0;

  /**
   *  Sum of pt of pi0 not in cone for tau
   */
  Energy taupTpi0;

  /**
   *  Sum of pt of tracks not in cone for tau
   */
  Energy taupTtracks;

  /**
   *  pT of highest track in tau decays
   */
  Energy ptHightestTrack;

  /**
   *  Cluster width
   */
  Energy clusterWidth;
};

/**
 * Here is the documentation of the PGSInterface class.
 *
 * @see \ref PGSInterfaceInterfaces "The interfaces"
 * defined for PGSInterface.
 */
class PGSInterface: public AnalysisHandler {

public:

  /**
   * The default constructor.
   */
  inline PGSInterface() : _pgs_param_file("lhc.par")
  {}

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
  //@}

  /**
   *  Public access to the reconstructed objects
   */
  inline vector<ReconstructedObject> reconstructedObjects() const {
    return _objects;
  }

  /**
   *  Missing transverse energy in calorimeter and azimuthal angle
   */
  inline pair<Energy,double> missingETCalorimeter() const {
    return _calorimeterMET;
  }
  /**
   *  Missing transverse energy after muon correction and azimuthal angle
   */
  inline pair<Energy,double> missingETCorrected() const {
    return _muonMET;
  }

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


protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  virtual void dofinish();
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<PGSInterface> initPGSInterface;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  PGSInterface & operator=(const PGSInterface &);

private:

  /**
   *  Converter from HepMC to HEPEVT
   */
  HepMC::IO_HEPEVT * _converter;

  /**
   *  file name for the PGS parameters
   */
  string _pgs_param_file;
  
  /**
   *  The reconstructed objects after PGS
   */
  vector<ReconstructedObject> _objects;

  /**
   *  missing ET measured in calorimeter
   */
  pair <Energy,double> _calorimeterMET;

  /**
   *  missing ET aftrer muon correction
   */
  pair <Energy,double> _muonMET;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of PGSInterface. */
template <>
struct BaseClassTrait<Herwig::PGSInterface,1> {
  /** Typedef of the first base class of PGSInterface. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the PGSInterface class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::PGSInterface>
  : public ClassTraitsBase<Herwig::PGSInterface> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::PGSInterface"; }
  /**
   * The name of a file containing the dynamic library where the class
   * PGSInterface is implemented. It may also include several, space-separated,
   * libraries if the class PGSInterface depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwPGSInterface.so"; }
};

/** @endcond */

}

#endif /* HERWIG_PGSInterface_H */
