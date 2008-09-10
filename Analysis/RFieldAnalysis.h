// -*- C++ -*-
#ifndef THEPEG_RFieldAnalysis_H
#define THEPEG_RFieldAnalysis_H
//
// This is the declaration of the RFieldAnalysis class.
//
#include "ThePEG/Handlers/AnalysisHandler.h"
#include "ThePEG/Vectors/Lorentz5Vector.h" 
#include "Herwig++/Shower/ShowerHandler.h"
#include "Herwig++/Utilities/Statistic.h"
#include "Herwig++/Utilities/Histogram.h"
#include "ThePEG/Repository/CurrentGenerator.h"

namespace Herwig {
    using namespace ThePEG;

/**
 * Simple struct to read in ascii files as long as histogram(2) is not 
 * capable of doing my sort of analysis
 */
  struct DataContainer{
    /** lower bin border */
    double low;
    /** upper bin border */
    double up;
    /** bin content */
    double content;
    /** bin error */
    double error;

    /** constructor */
    DataContainer(){
      low = 0.0;
      up  = 0.0;
      content = 0.0;
      error = 0.0;
    }
  };

/**
 * Here is the documentation of the RFieldAnalysis class.
 *
 * @see \ref RFieldAnalysisInterfaces "The interfaces"
 * defined for RFieldAnalysis.
 */
class RFieldAnalysis: public AnalysisHandler {

public:

  /**
   * The default constructor.
   */
  RFieldAnalysis() : thelow(30), theup(50), theDir("."), 
                     theRealMult(-0.5, 10.5, 11) {}

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

  /**
   * Transform the event to the desired Lorentz frame and return the
   * corresponding LorentzRotation.
   * @param event a pointer to the Event to be transformed.
   * @return the LorentzRotation used in the transformation.
   */
  virtual LorentzRotation transform(tEventPtr event) const;

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
   * Read in data ascii file.
   */
  vector<DataContainer *> getData(string ascii);

  /**
   * Calculate chi squared compared to data in ascii file.
   */
  void chisq(vector<Statistic> mc, string ascii, pair<double, int> &tot);

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<RFieldAnalysis> initRFieldAnalysis;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  RFieldAnalysis & operator=(const RFieldAnalysis &);

private:

  /** lower border of chi^2 comparison to data */
  int thelow;

  /** upper border of chi^2 comparison to data */
  int theup;

  /** dir where Data files are stored */
  string theDir;

  /**  N_ch in towards region */
  vector<Statistic> theNTow;

  /**  N_ch in transverse region */
  vector<Statistic> theNTrans;

  /**  N_ch in away region */
  vector<Statistic> theNAway;

  /**  p_Tsum in towards region */
  vector<Statistic> thePtsumTow;

  /**  p_Tsum in transverse region */
  vector<Statistic> thePtsumTrans;

  /**  p_Tsum in away region */
  vector<Statistic> thePtsumAway;

  /**
   * Histogram for the real extra scatter multiplicity
   */
  Histogram theRealMult;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of RFieldAnalysis. */
template <>
struct BaseClassTrait<Herwig::RFieldAnalysis,1> {
  /** Typedef of the first base class of RFieldAnalysis. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the RFieldAnalysis class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::RFieldAnalysis>
  : public ClassTraitsBase<Herwig::RFieldAnalysis> {
  /** Return a platform-independent class name */
    static string className() { return "Herwig::RFieldAnalysis"; }
  /** Return the name(s) of the shared library (or libraries) be loaded to get
   *  access to the RFieldAnalysis class and any other class on which it depends
   *  (except the base class). */
  static string library() { return "HwShower.so HwMPIAnalysis.so"; }
};

/** @endcond */

}

#endif /* THEPEG_RFieldAnalysis_H */
