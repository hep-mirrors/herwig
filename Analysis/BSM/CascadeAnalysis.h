// -*- C++ -*-
#ifndef HERWIG_CascadeAnalysis_H
#define HERWIG_CascadeAnalysis_H
//
// This is the declaration of the CascadeAnalysis class.
//

#include "ThePEG/Handlers/AnalysisHandler.h"
#include "Herwig++/Utilities/Histogram.h"
#include "CascadeAnalysis.fh"

namespace Herwig {
using namespace ThePEG;

/**
 * This class is designed to compare the theoretical and numerical
 * invariant mass distributions for the near and far quark-lepton 
 * in either neutral or charged Susy and UED cascade decays.
 *
 * For the neutral decay chain the decay cascade must look like
 * \f[ D \to C\,q \to B\,l_n\,q \to A\,l_f\,l_n\,q \f]. Histograms of
 * the invariant mass spectrum for the quark with the near/far-leptons
 * are booked then the combined quark + lepton/anti-lepton and finally
 * the dilepton.
 *
 * For the charge chain the cascade must have this form:
 * \f[C \to B\,q \to W\,A \to l\,\nu_l\,A \f].
 * Histograms of the quark + lepton invariant mass are then booked.
 * 
 * There is an interface to set the mode of the analysis handler and
 * interfaces to set the PDG codes of the various particles along the chain.
 * The default analysis is the neutral chain in the MSSM.
 *
 * @see \ref CascadeAnalysisInterfaces "The interfaces"
 * defined for CascadeAnalysis.
 */
class CascadeAnalysis: public AnalysisHandler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  CascadeAnalysis();

public:

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
   * Find the particle after it has been showered.
   * @param p The particle to be found
   */
  tPPtr showeredProduct(tPPtr p ) const;
  //@}
  
private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<CascadeAnalysis> initCascadeAnalysis;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  CascadeAnalysis & operator=(const CascadeAnalysis &);

private:
  
  /** @name Cascade analysis functions. */
  //@{
  /**
   * Analyse a neutral decay chain
   * @param inpart The first particle in the cascade
   * @param products The first decay products
   */
  void analyseNeutralChain(tPPtr inpart, ParticleVector & products);

  /**
   * Analyse a charged decay chain
   * @param inpart The first particle in the cascade
   * @param products The first decay products
   */
  void analyseChargedChain(tPPtr inpart, ParticleVector & products);
  //@}
  
  /** @name Calculate analytic results.*/
  //@{
    /**
   * Pick which analytic chain to do
   * @param s Identify the process
   * @param m The rescaled mass value
   */
  double neutralChain(string s, double m);

  /**
   * Plot the charged chain distribution
   * @param s Identify the process
   * @param m The mass value
   */
  double chargedChain(string s, double m);

  /**
   * The analytic functions for the SUSY chain
   * @param s The distribution
   * @param m The value of the rescaled mass
   */
  double neutralChainSUSY(string s, double m);

  /**
   * The analytic functions for the UED chain
   * @param s The distribution
   * @param m The value of the rescaled mass
   */
  double neutralChainUED(string s, double m);
  //@}

private:

  /** @name Switch variables. */
  //@{
  /**
   * The model under study i.e. MSSM or MUED.
   */
  unsigned int theModel;
    
  /**
   * Whether to normalise the invariant mass to their 
   * maximum.
   */
  bool theMassNormalize;
  
  /**
   * The number of bins on each histogram 
   */
  unsigned int theNBins;
  //@}

  /**
   * The PDG codes of the cascade particles D,B,C,A
   */
  vector<long> theResonances;

  /**
   * The histograms indexed by their name.
   */
  map<string, HistogramPtr> theHistograms;
  
  /**
   * The maximum invariant masses. Positions 0,1,2 and 3 give the maxima
   * for the neutral chain and position 4 gives the rescale mass prefactor
   * for the charged chain. 
   */
  vector<Energy> theMassMaxima; 

  /**
   * The values of the dimensionless quantities \f$ \frac{m_i^2}{m_j^2}\f$
   */
  double theX, theY, theZ, theYc, theZc;
    
  /**
   * The number of decay chains indexed with a string identifier.
   */
  vector<long> theNChain;

  /**
   * The fractions of chains initiated by a quark-partner or anti-partner.
   */
  vector<double> theFractions;

  /**
   * The vertex depenedent factor in the charged distributions
   */
  double theAlpha;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of CascadeAnalysis. */
template <>
struct BaseClassTrait<Herwig::CascadeAnalysis,1> {
  /** Typedef of the first base class of CascadeAnalysis. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the CascadeAnalysis class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::CascadeAnalysis>
  : public ClassTraitsBase<Herwig::CascadeAnalysis> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::CascadeAnalysis"; }
  /**
   * The name of a file containing the dynamic library where the class
   * CascadeAnalysis is implemented. It may also include several, space-separated,
   * libraries if the class CascadeAnalysis depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwSusy.so HwUED.so HwBSMAnalysis.so"; }
};

/** @endcond */

}

#include "CascadeAnalysis.icc"

#endif /* HERWIG_CascadeAnalysis_H */
