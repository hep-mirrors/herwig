// -*- C++ -*-
//
// y23.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_y23_H
#define HERWIG_y23_H
//
// This is the declaration of the y23 class.
//

#include "ThePEG/Handlers/AnalysisHandler.h"
#include "Herwig++/Utilities/Histogram.h"
#include "Herwig++/Interfaces/KtJetInterface.h"
#include "Herwig++/Shower/Base/ShowerParticle.h"
#include "KtJet/KtJet.h"
#include "KtJet/KtDistance.h"
#include "KtJet/KtLorentzVector.h"


namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the y23 class.
 *
 * @see \ref y23Interfaces "The interfaces"
 * defined for y23.
 */
class y23: public AnalysisHandler {

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
   * Indicates that this is an concrete class without persistent data.
   */
  static NoPIOClassDescription<y23> inity23;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  y23 & operator=(const y23 &);

private:

  /**
   * Checks to see that the splitting is allowed and finds the
   * Sudakov form factor for the splitting.
   */
  bool getIds(int & qq_pairs, long & emmitter_id,
		    ShowerParticlePtr & part_i, 
		    ShowerParticlePtr & part_j ) ;

  /**
   * Returns the durham jet measure, yij, for the two particles. 
   */
  double getJetMeasure(ShowerParticlePtr part_i, ShowerParticlePtr part_j);

  /**
   * Checks to see that the splitting is allowed.
   */
  bool splittingAllowed( ShowerParticlePtr part_i,
			 ShowerParticlePtr part_j,
			 int qq_pairs);

  /**
   * Clusters the partons and creates a branching history
   * by combining the 2 particles with smallest
   * jet measure out of all allowed pairings until we are left 
   * with \f$q\bar{q}\f$.
   */
  double doClustering( tPVector parts, int clustered_no);

  /**
   *  Histograms of the \f$y\f$ distributions
   */
  //@{
  /**
   *  \f$y_{23}\f$
   */
  HistogramPtr _y23;

  HistogramPtr _y23_luc;

  /**
   *  The interface between Herwig++ and KtJet
   */
  Herwig::KtJetInterface _kint;

  /**
   *  c.o.m energy - used in getjetmeasure
   */
  Energy _s;

  
};

}
namespace KtJet{
 class KtDistanceLuc : public KtDistance {
  public:
    KtDistanceLuc(int collision_type=1);
    virtual ~KtDistanceLuc(){}
    /** Jet Kt */
    KtFloat operator()(const KtLorentzVector &) const;
    /** Pair Kt */
    KtFloat operator()(const KtLorentzVector &, const KtLorentzVector &) const;
    /** Name of scheme */
    std::string name() const;
  private:
    int m_type;
    std::string m_name;
  };
}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of y23. */
template <>
struct BaseClassTrait<Herwig::y23,1> {
  /** Typedef of the first base class of y23. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the y23 class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::y23>
  : public ClassTraitsBase<Herwig::y23> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::y23"; }
  /**
   * The name of a file containing the dynamic library where the class
   * y23 is implemented. It may also include several, space-separated,
   * libraries if the class y23 depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so HwKtJet.so HwLEPJetAnalysis.so"; }
};

/** @endcond */

}

#endif /* HERWIG_y23_H */
