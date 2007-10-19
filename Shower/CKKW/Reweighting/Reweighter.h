// -*- C++ -*-
#ifndef HERWIG_Reweighter_H
#define HERWIG_Reweighter_H
//
// This is the declaration of the Reweighter class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "ThePEG/Utilities/Rebinder.h"
#include "Reweighter.fh"
#include "ThePEG/Repository/EventGenerator.h"
#include "Herwig++/Shower/Couplings/ShowerAlpha.h"
#include "Herwig++/Shower/CKKW/Clustering/CascadeReconstructor.h"
#include "JetMeasure.h"

#ifdef HERWIG_CHECK_CKKW_REWEIGHTING
#include "Herwig++/Utilities/Histogram2/Histogram2.h"
#endif

namespace Herwig {

using namespace ThePEG;

/**\ingroup CKKW
 *
 * Reweighter is the base class for reweighting matrix elements
 * in a ME/PS merging approach.
 *
 *@author Simon Plaetzer
 *
 * @see \ref ReweighterInterfaces "The interfaces"
 * defined for Reweighter.
 */
class Reweighter: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline Reweighter();

  /**
   * The destructor.
   */
  virtual ~Reweighter();
  //@}

public:

  /**
   * Check for unresolved partons. If any found
   * throw a veto.
   */
  void unresolvedCut (PPair in, PVector out);

  /**
   * Perform the reweighting. It calls the Sudakov
   * reweighting and then the alpha_s reweighting.
   */
  double reweight (CascadeHistory, unsigned int, unsigned int);
  
  /**
   * Analyze the given cascade history
   * and set missing scales.
   */
  virtual void analyzeHistory (CascadeHistory) = 0;

  /**
   * Setup this Reweighter. This is to be called
   * during the init phase from the ShowerHandler
   * using it.
   */
  virtual void setup (CascadeReconstructorPtr);

  /**
   * Initialize this reweighter. This is to be called
   * just before the run phase from the ShowerHandler
   * using it.
   */
  virtual void initialize () = 0;

  /**@name Methods used for reweighting. */
  //@{

  /**
   * Perform the coupling reweighting.
   */
  double couplingReweight (CascadeHistory);

  /**
   * Perform the Sudakov reweighting. The default
   * returns 1.
   */
  virtual double sudakovReweight (CascadeHistory, unsigned int, unsigned int);

  //@}

  /**@name Member access. */
  //@{

  /**
   * Return the JetMeasure
   */
  inline tJetMeasurePtr resolution () const;

  /**
   * Return true, if showering off the highest multiplicty is
   * to be vetoed.
   */
  inline bool vetoHighest () const;

  /**
   * Return the fixed running coupling the
   * matrix elements were generated with.
   */
  inline double MEalpha () const;

  /**
   * Set the fixed running coupling the
   * matrix elements were generated with.
   */
  inline void MEalpha (double);

  /**
   * Return the shower qcd running
   * coupling.
   */
  inline tShowerAlphaPtr showerAlpha () const;

  /**
   * Get the associated CascadeReconstructor
   */
  inline CascadeReconstructorPtr reconstructor () const;

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


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).

protected:

  /** @name Standard Interfaced functions. */
  //@{

#ifdef HERWIG_CHECK_CKKW_REWEIGHTING

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  inline virtual void doinitrun();

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  inline virtual void dofinish();

#endif

  /**
   * Rebind pointer to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   * @throws RebindException if no cloned object was found for a given
   * pointer.
   */
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);

  /**
   * Return a vector of all pointers to Interfaced objects used in this
   * object.
   * @return a vector of pointers.
   */
  inline virtual IVector getReferences();
  //@}


private:

  /**
   * The jet resolution to be used.
   */
  JetMeasurePtr _resolution;

  /**
   * The associated CascadeReconstructor
   */
  CascadeReconstructorPtr _reconstructor;

  /**
   * True, if the highest multiplicity is to be vetoed.
   */
  bool _vetoHighest;

  /**
   * The fixed coupling the matrix elements
   * were generated with.
   */
  double _MEalpha;

  /**
   * The QCD running coupling.
   */
  ShowerAlphaPtr _showerAlpha;

#ifdef HERWIG_CHECK_CKKW_REWEIGHTING

  /**
   * Map multiplicities to number of events
   * after cut and average CKKW weight
   */
  map<unsigned int, pair<unsigned long, double> > _stats;

  /**
   * Histogram to store scales
   */
  Histogram2Ptr _clustering_scales;

  /**
   * Histogram to store weights
   */
  Histogram2Ptr _weights;

  /**
   * Map multiplicities to channel names
   */
  map<unsigned int,string> _mult;

  /**
   * Map multiplicities and clustering number
   * to channel names
   */
  map<pair<unsigned int,unsigned int>,string> _mult_cluster;

#endif

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class with persistent data.
   */
  static AbstractClassDescription<Reweighter> initReweighter;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  Reweighter & operator=(const Reweighter &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of Reweighter. */
template <>
struct BaseClassTrait<Herwig::Reweighter,1> {
  /** Typedef of the first base class of Reweighter. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the Reweighter class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::Reweighter>
  : public ClassTraitsBase<Herwig::Reweighter> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::Reweighter"; }
  /**
   * The name of a file containing the dynamic library where the class
   * Reweighter is implemented. It may also include several, space-separated,
   * libraries if the class Reweighter depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#include "Reweighter.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "Reweighter.tcc"
#endif

#endif /* HERWIG_Reweighter_H */
