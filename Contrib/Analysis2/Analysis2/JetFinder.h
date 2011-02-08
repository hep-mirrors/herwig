// -*- C++ -*-

// (C) 2007-2009 Simon Plaetzer -- sp@particle.uni-karlsruhe.de

#ifndef Analysis2_JetFinder_H
#define Analysis2_JetFinder_H
//
// This is the declaration of the JetFinder class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "JetFinder.fh"

#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Vectors/Lorentz5Vector.h"

namespace Analysis2 {

using namespace ThePEG;

/**\ingroup Analysis2
 * 
 * A general interface to a jet finder which
 * may be used by an Analysis2 handler.
 *
 * @author Simon Plaetzer
 *
 * @see \ref JetFinderInterfaces "The interfaces"
 * defined for JetFinder.
 */
class JetFinder: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline JetFinder();

  /**
   * The destructor.
   */
  virtual ~JetFinder();
  //@}

public:

  /**@name Resolution parameters */
  //@{

  /**
   * Get the resolution scale
   */
  inline Energy dCut () const;

  /**
   * Get the dimensionless resolution parameter.
   */
  inline double yCut () const;

  /**
   * Get the optional cone radius.
   */
  inline double R () const;

  //@}

public:

  /**
   * Set the event to be analysed.
   * The implementation of this method may
   * also do some conversion.
   */
  virtual inline void use (const vector<Lorentz5Momentum>&, bool inclusive = false);

  /**
   * Get the event to be analysed
   */
  inline const vector<Lorentz5Momentum>& lastEvent () const;

  /**
   * Return the jets resolved
   */
  inline const vector<Lorentz5Momentum>& jets () const;

public:

  /**@name Jet finding */
  //@{

  /**
   * Do inclusive jet finding
   */
  virtual void findJets () = 0;

  /**
   * Do exclusive jet-finding for nJets jets
   */
  virtual void findJetsN (unsigned int nJets) = 0;

  /**
   * Do jet-finding up to the resolution scale
   */
  virtual void findJetsD () = 0;

  /**
   * Do exclusive jet-finding up to dimensionless resolution
   */
  virtual void findJetsY (double y = -1.) = 0;

  /**
   * Get the number of final state jets
   */
  inline unsigned int getNJets () const;

  /**
   * Resolution scale where n+1 jets merged to n
   */
  virtual Energy getDMerge (unsigned int nJets) const = 0;

  /**
   * Dimensionless scale where n+1 jets merged to n
   */
  virtual double getYMerge (unsigned int nJets) const = 0;

  //@}

public:

  /**@name Jet sorting */
  //@{

  /**
   * Sort the jets in decreasing energy.
   */
  virtual inline void sortEnergy ();

  /**
   * Sort the jets in decreasing transverse energy.
   */
  virtual inline void sortEt ();

  /**
   * Sort the jets in decreasing transverse momentum.
   */
  virtual inline void sortPt ();

  /**
   * Sort the jets in decreasing rapidity.
   */
  virtual inline void sortY ();

  /**
   * Sort the jets in decreasing pseudorapidity.
   */
  virtual inline void sortEta ();

  //@}

protected:

  /**@name Binary predicats for jet sorting. */
  //@{

  /**
   * Sort in energy
   */
  struct GreaterInEnergy {
    /**
     * Return true, if less in energy
     */
    inline bool operator () (const Lorentz5Momentum& a, const Lorentz5Momentum& b) {
      return a.e() > b.e();
    }
  };

  /**
   * Sort in transverse energy
   */
  struct GreaterInEt {
    /**
     * Return true, if less in transverse energy
     */
    inline bool operator () (const Lorentz5Momentum& a, const Lorentz5Momentum& b) {
      return a.et() > b.et();
    }
  };

  /**
   * Sort in pt
   */
  struct GreaterInPt {
    /**
     * Return true, if less in pt
     */
    inline bool operator () (const Lorentz5Momentum& a, const Lorentz5Momentum& b) {
      return a.perp() > b.perp();
    }
  };

  /**
   * Sort in y
   */
  struct GreaterInY {
    /**
     * Return true, if less in y
     */
    inline bool operator () (const Lorentz5Momentum& a, const Lorentz5Momentum& b) {
      return a.rapidity() > b.rapidity();
    }
  };

  /**
   * Sort in eta
   */
  struct GreaterInEta {
    /**
     * Return true, if less in eta
     */
    inline bool operator () (const Lorentz5Momentum& a, const Lorentz5Momentum& b) {
      return a.eta() > b.eta();
    }
  };

  //@}

protected:

  /**
   * Set the result from jet finding
   */
  inline void jets (const vector<Lorentz5Momentum>&);

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


private:

  /**
   * The dimensionful resolution parameter
   */
  Energy _dCut;

  /**
   * The dimensionless resolution parameter
   */
  double _yCut;

  /**
   * The optional cone radius.
   */
  double _R;

  /**
   * The event to be analysed.
   */
  vector<Lorentz5Momentum> _lastEvent;

  /**
   * The jets resolved.
   */
  vector<Lorentz5Momentum> _jets;

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class with persistent data.
   */
  static AbstractClassDescription<JetFinder> initJetFinder;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  JetFinder & operator=(const JetFinder &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of JetFinder. */
template <>
struct BaseClassTrait<Analysis2::JetFinder,1> {
  /** Typedef of the first base class of JetFinder. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the JetFinder class and the shared object where it is defined. */
template <>
struct ClassTraits<Analysis2::JetFinder>
  : public ClassTraitsBase<Analysis2::JetFinder> {
  /** Return a platform-independent class name */
  static string className() { return "Analysis2::JetFinder"; }
  /**
   * The name of a file containing the dynamic library where the class
   * JetFinder is implemented. It may also include several, space-separated,
   * libraries if the class JetFinder depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "Analysis2.so"; }
};

/** @endcond */

}

#include "JetFinder.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "JetFinder.tcc"
#endif

#endif /* Analysis2_JetFinder_H */
