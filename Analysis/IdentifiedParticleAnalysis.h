// -*- C++ -*-
#ifndef HERWIG_IdentifiedParticleAnalysis_H
#define HERWIG_IdentifiedParticleAnalysis_H
//
// This is the declaration of the IdentifiedParticleAnalysis class.
//

#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Handlers/AnalysisHandler.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "EventShapes.h"
#include "Herwig++/Utilities/Histogram.h"
#include "IdentifiedParticleAnalysis.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the IdentifiedParticleAnalysis class.
 *
 * @see \ref IdentifiedParticleAnalysisInterfaces "The interfaces"
 * defined for IdentifiedParticleAnalysis.
 */
class IdentifiedParticleAnalysis: public AnalysisHandler {

public:

  /** @name Virtual functions required by the AnalysisHandler class. */
  //@{
  /**
   * Analyze the given vector of particles. The default version calls
   * analyze(tPPtr) for each of the particles.
   * @param particles the vector of pointers to particles to be analyzed
   */
  virtual void analyze(const tPVector & particles);
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

  /**
   *  Work out the flavour of the quarks produced
   */
  inline int getFlavour(const tPVector &);

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
  static ClassDescription<IdentifiedParticleAnalysis> initIdentifiedParticleAnalysis;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  IdentifiedParticleAnalysis & operator=(const IdentifiedParticleAnalysis &);

private:

  /**
   *  Single particle spectra
   */
  //@{
  /**
   * Histogram for the \f$\xi\f$ distribution for all particles from all quarks
   */
  HistogramPtr _xpa;

  /**
   * Histogram for the \f$\xi\f$ distribution for all particles from light quarks
   */
  HistogramPtr _xpl;

  /**
   * Histogram for the \f$\xi\f$ distribution for all particles from charm quarks
   */
  HistogramPtr _xpc;

  /**
   * Histogram for the \f$\xi\f$ distribution for all particles from bottom quarks
   */
  HistogramPtr _xpb;

  /**
   * Histogram for the \f$\xi\f$ distribution for charged pions from all quarks
   */
  HistogramPtr _pipma;

  /**
   * Histogram for the \f$\xi\f$ distribution for charged pions from light quarks
   */
  HistogramPtr _pipml;

  /**
   * Histogram for the \f$\xi\f$ distribution for charged pions from charm quarks
   */
  HistogramPtr _pipmc;

  /**
   * Histogram for the \f$\xi\f$ distribution for charged pions from bottom quarks
   */
  HistogramPtr _pipmb;

  /**
   * Histogram for the \f$\xi\f$ distribution for charged pions from OPAL
   */
  HistogramPtr _pipm;

  /**
   * Histogram for the \f$\xi\f$ distribution for charged kaons from all quarks
   */
  HistogramPtr _kpma;

  /**
   * Histogram for the \f$\xi\f$ distribution for charged kaons from light quarks
   */
  HistogramPtr _kpml;

  /**
   * Histogram for the \f$\xi\f$ distribution for charged kaons from charm quarks
   */
  HistogramPtr _kpmc;

  /**
   * Histogram for the \f$\xi\f$ distribution for charged kaons from bottom quarks
   */
  HistogramPtr _kpmb;

  /**
   * Histogram for the \f$\xi\f$ distribution for charged kaons from OPAL
   */
  HistogramPtr _kpm;


  /**
   * Histogram for the \f$\xi\f$ distribution for protons from all quarks
   */
  HistogramPtr _ppma;

  /**
   * Histogram for the \f$\xi\f$ distribution for protons from light quarks
   */
  HistogramPtr _ppml;

  /**
   * Histogram for the \f$\xi\f$ distribution for protons from charm quarks
   */
  HistogramPtr _ppmc;

  /**
   * Histogram for the \f$\xi\f$ distribution for protons from bottom quarks
   */
  HistogramPtr _ppmb;

  /**
   * Histogram for the \f$\xi\f$ distribution for protons from OPAL
   */
  HistogramPtr _ppm;

  /**
   * Histogram for the \f$x\f$ distribution for light quark events (lin)
   */ 
  HistogramPtr _udsxp;

  /**
   * Histogram for the \f$\xi\f$ distribution for light quark events (lin)
   */ 
  HistogramPtr _udsxip;

  /**
   *  Histogram for the \f$\xi\f$ distribution for \f$\Lambda\f$ 
   */
  HistogramPtr _lpm;

  /**
   *  Histogram for the ALEPH \f$K^{*\pm}\f$ \f$x\f$distribution
   */
  HistogramPtr _xpKstarplus;

  /**
   *  Histogram for the OPAL \f$\Xi^-\f$ \f$x\f$ distribution
   */
  HistogramPtr _xpXiminus;

  /**
   *  Histogram for the OPAL \f$\Xi^-\f$ \f$\xi\f$ distribution
   */
  HistogramPtr _xiXiminus;

  /**
   *  Histogram for the OPAL \f$\Sigma^{*+}\f$ \f$x\f$ distribution
   */
  HistogramPtr _xpSigmaplus;

  /**
   *  Histogram for the OPAL \f$\Sigma^{*+}\f$ \f$\xi\f$ distribution
   */
  HistogramPtr _xiSigmaplus;

  /**
   *  Histogram for the OPAL \f$\Sigma^{*-}\f$ \f$x\f$ distribution
   */
  HistogramPtr _xpSigmaminus;

  /**
   *  Histogram for the OPAL \f$\Sigma^{*-}\f$ \f$\xi\f$ distribution
   */
  HistogramPtr _xiSigmaminus;

  /**
   *  Histogram for the OPAL \f$\Xi^{*0}\f$ \f$x\f$ distribution
   */
  HistogramPtr _xpXi0;

  /**
   *  Histogram for the OPAL \f$\Xi^{*0}\f$ \f$\xi\f$ distribution
   */
  HistogramPtr _xiXi0;

  /**
   *  Histogram for \f$\Lambda(1520)\f$ \f$x\f$ distribution
   */
  HistogramPtr _xpLambda1520;

  /**
   *  Histogram for \f$\Lambda(1520)\f$ \f$\xi\f$ distribution
   */
  HistogramPtr _xiLambda1520;

  /**
   *  Histogram for \f$\Delta^{++}\f$ \f$x\f$ distribution
   */
  HistogramPtr _xeDelta;

  /**
   *  Histogram for \f$f_0(980)\f$ \f$x\f$ distribution
   */
  HistogramPtr _xpf980;

  /**
   *  Histogram for \f$\phi\f$ \f$x\f$ distribution
   */
  HistogramPtr _xpphi;

  /**
   *  Histogram for \f$f_2\f$ \f$x\f$ distribution
   */
  HistogramPtr _xpf2;

  /**
   *  Histogram for \f$K^{*0}\f$ \f$x\f$ distribution
   */
  HistogramPtr _xpKstar0;  

  /**
   *  Histogram for \f$K^0\f$ \f$x\f$ distribution
   */
  HistogramPtr _xpK0; 

  /**
   *  Histogram for \f$\rho^0\f$ \f$x\f$ distribution
   */
  HistogramPtr _xerho0;

  /**
   *  Histogram for the OPAL \f$\pi^0\f$ \f$x\f$ distribution
   */
  HistogramPtr _xepi0;

  /**
   *  Histogram for the OPAL \f$\pi^0\f$ \f$\xi\f$ distribution
   */
  HistogramPtr _xipi0;

  /**
   *  Histogram for the OPAL \f$\eta\f$ \f$x\f$ distribution
   */
  HistogramPtr _xeeta;

  /**
   *  Histogram for the OPAL \f$\eta'\f$ \f$\xi\f$ distribution
   */
  HistogramPtr _xieta;

  /**
   *  Histogram for the OPAL \f$\eta\f$ \f$x\f$ distribution
   */
  HistogramPtr _xeetap;

  /**
   *  Histogram for the OPAL \f$\eta'\f$ \f$\xi\f$ distribution
   */
  HistogramPtr _xietap;

  /**
   *  Histogram for the OPAL \f$\omega\f$ \f$x\f$ distribution
   */
  HistogramPtr _xeomega;

  /**
   *  Histogram for the OPAL \f$\omega'\f$ \f$\xi\f$ distribution
   */
  HistogramPtr _xiomega;

  /**
   *  Histogram for the OPAL \f$\rho^+\f$ \f$x\f$ distribution
   */
  HistogramPtr _xerhop;

  /**
   *  Histogram for the OPAL \f$\rho^+\f$ \f$\xi\f$ distribution
   */
  HistogramPtr _xirhop;

  /**
   *  Histogram for the OPAL \f$\a_0^+\f$ \f$x\f$ distribution
   */
  HistogramPtr _xea_0p;

  /**
   *  Histogram for the OPAL \f$\a_0^+\f$ \f$\xi\f$ distribution
   */
  HistogramPtr _xia_0p;

  /**
   *  Histogram for \f$D^0\f$ \f$x\f$ distribution
   */
  HistogramPtr _xeD0; 

  /**
   *  Histogram for \f$D^{*+}\f$ \f$x\f$ distribution
   */
  HistogramPtr _xeDstar;

  /**
   *  Pointer to the event shapes object
   */
  EventShapesPtr _shapes;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of IdentifiedParticleAnalysis. */
template <>
struct BaseClassTrait<Herwig::IdentifiedParticleAnalysis,1> {
  /** Typedef of the first base class of IdentifiedParticleAnalysis. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the IdentifiedParticleAnalysis class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::IdentifiedParticleAnalysis>
  : public ClassTraitsBase<Herwig::IdentifiedParticleAnalysis> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::IdentifiedParticleAnalysis"; }
  /** Return the name(s) of the shared library (or libraries) be loaded to get
   *  access to the IdentifiedParticleAnalysis class and any other class on which it depends
   *  (except the base class). */
  static string library() { return "HwAnalysis.so HwLEPAnalysis.so"; }
};

/** @endcond */

}

#include "IdentifiedParticleAnalysis.icc"

#endif /* HERWIG_IdentifiedParticleAnalysis_H */
