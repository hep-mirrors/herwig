// -*- C++ -*-
#ifndef HERWIG_DefaultReweighter_H
#define HERWIG_DefaultReweighter_H
//
// This is the declaration of the DefaultReweighter class.
//

#include "ThePEG/Utilities/Rebinder.h"

#include "Herwig++/Shower/ShowerConfig.h"
#include "Herwig++/Shower/SplittingFunctions/SplittingFunction.h"
#include "Reweighter.h"
#include "DefaultReweighter.fh"

#include "DefaultSudakov.h"

namespace Herwig {

using namespace ThePEG;

  struct splittingKey {

    IdList ids;
    bool initial;

    inline bool operator < (const splittingKey& x) const 
    { return ids[0] < x.ids[0] || ids[1] < x.ids[1]; }

  };

  inline PersistentOStream& operator << (PersistentOStream& os, const splittingKey& k) {
    os << k.ids << k.initial;
    return os;
  }

  inline PersistentIStream& operator >> (PersistentIStream& is, splittingKey& k) {
    is >> k.ids >> k.initial;
    return is;
  }


/**\ingroup CKKW
 *
 * DefaultReweighter performs the Sudakov reweighting
 * for standard CKKW approaches.
 *
 *@author Simon Plaetzer
 *
 * @see \ref DefaultReweighterInterfaces "The interfaces"
 * defined for DefaultReweighter.
 */
class DefaultReweighter: public Reweighter {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline DefaultReweighter();

  /**
   * The destructor.
   */
  virtual ~DefaultReweighter();
  //@}

public:

  /**
   * Perform the Sudakov reweighting.
   */
  virtual double sudakovReweight (CascadeHistory, unsigned int, unsigned int);

  /**@name Methods used for setup. */
  //@{
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
  virtual void initialize ();

  /**
   * Insert a branching and splitting function according
   * to which this branching is generated.
   */
  inline void insertSplitting (const IdList&, SplittingFnPtr, bool initial = false);

  //@}


  /**@name Methods related to integration and interpolation */
  //@{

  /**
   * Return the spacing between interpolation points
   */
  inline Energy2 interpolationSpacing () const;

  /**
   * Return the relative accuracy to which branching
   * probabilities should be integrated.
   */
  inline double integrationAccuracy () const;

  /**
   * Return the maximum scale up to which Sudakov
   * from factors are precomputed.
   */
  inline Energy2 sudakovMaxScale () const;

  /**
   * Return the path where Sudakov interpolation data
   * is stored.
   */
  inline string sudakovDataPath () const;

  /**
   * Return true, if mass dependend splitting
   * functions should be used.
   */
  inline bool useMassiveSplittings () const;

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

#ifdef HERWIG_DEBUG_CKKW_CHECK_SUDAKOVS

protected:

  inline virtual void dofinish();

#endif


protected:

  /** @name Standard Interfaced functions. */
  //@{

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
   * A map of branchings to associated splitting functions.
   */
  map<splittingKey,SplittingFnPtr> _splittingMap;

  /**
   * A map of partons (id,initial) to Sudakovs.
   */
  multimap<pair<long,bool>, DefaultSudakovPtr> _sudakovMap;

  /**
   * The spacing between interpolation points
   */
  Energy2 _interpolationSpacing;

  /**
   * The relative integration accuracy
   */
  double _integrationAccuracy;

  /**
   * The maximum scale up to which Sudakov 
   * form factors are precomputed.
   */
  Energy2 _sudakovMaxScale;

  /**
   * The path where Sudakov interpolation data
   * is stored.
   */
  string _sudakovDataPath;

  /**
   * True, if massive splitting functions should be used.
   */
  bool _useMassive;

  /**
   * Wether or not we should deliver the weight
   * divided by the maximum Sudakov weight (i.e.
   * the hard rocess weight).
   */
  bool _sudakovUnweight;

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class with persistent data.
   */
  static AbstractClassDescription<DefaultReweighter> initDefaultReweighter;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DefaultReweighter & operator=(const DefaultReweighter &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of DefaultReweighter. */
template <>
struct BaseClassTrait<Herwig::DefaultReweighter,1> {
  /** Typedef of the first base class of DefaultReweighter. */
  typedef Herwig::Reweighter NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the DefaultReweighter class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::DefaultReweighter>
  : public ClassTraitsBase<Herwig::DefaultReweighter> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::DefaultReweighter"; }
  /**
   * The name of a file containing the dynamic library where the class
   * DefaultReweighter is implemented. It may also include several, space-separated,
   * libraries if the class DefaultReweighter depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#include "DefaultReweighter.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "DefaultReweighter.tcc"
#endif

#endif /* HERWIG_DefaultReweighter_H */
