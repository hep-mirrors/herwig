// -*- C++ -*-

// (C) 2007-2009 Simon Plaetzer -- sp@particle.uni-karlsruhe.de

#ifndef Analysis2_FastJetFinder_H
#define Analysis2_FastJetFinder_H
//
// This is the declaration of the FastJetFinder class.
//

#include "JetFinder.h"
#include "FastJetFinder.fh"

#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

#include <memory>

namespace Analysis2 {

using namespace ThePEG;

/**\ingroup Analysis2
 * 
 * A jet finder using the FastJetFinder library.
 * See the <a href="http://www.lpthe.jussieu.fr/~salam/fastjet">FastJet</a>
 * homepage for details.
 *
 * @author Simon Plaetzer
 *
 * @see \ref FastJetFinderInterfaces "The interfaces"
 * defined for FastJetFinder.
 */
class FastJetFinder: public JetFinder {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline FastJetFinder();

  /**
   * The copy constructor
   */
  inline FastJetFinder (const FastJetFinder&);

  /**
   * The destructor.
   */
  virtual ~FastJetFinder();
  //@}

public:

  /**
   * Set the event to be analysed.
   */
  virtual void use (const vector<Lorentz5Momentum>&,bool);

  /**@name Jet finding */
  //@{

  /**
   * Do inclusive jet finding
   */
  virtual inline void findJets ();

  /**
   * Do exclusive jet-finding for nJets jets
   */
  virtual inline  void findJetsN (unsigned int nJets);

  /**
   * Do jet-finding up to the resolution scale
   */
  virtual inline void findJetsD ();

  /**
   * Do exclusive jet-finding up to dimensionless resolution
   */
  virtual inline void findJetsY (double y = -1.);

  /**
   * Resolution scale where n+1 jets merged to n
   */
  virtual inline Energy getDMerge (unsigned int nJets) const;

  /**
   * Dimensionless scale where n+1 jets merged to n
   */
  virtual inline double getYMerge (unsigned int nJets) const;

  //@}

public:

  /**@name FastJet options */
  //@{

  /**
   * Get the jet finder
   */
  inline int jetFinder () const;

  /**
   * Get the recombination scheme
   */
  inline int recombinationScheme () const;

  /**
   * Get the startegy
   */
  inline int strategy () const;

  /**
   * Get the jet definition
   */
  inline const fastjet::JetDefinition& jetDefintion () const;

  //@}

protected:

  /**
   * Convert event to FastJet pseudojet vector
   */
  void convert ();

  /**
   * Convert pseudojet vector to vector of Lorentz5Momenta
   */
  void convert (const vector<fastjet::PseudoJet>&);

  /**
   * Perform clustering
   */
  void cluster ();

  /**
   * Return the vector of pseudojets from
   * given event.
   */
  inline const vector<fastjet::PseudoJet>& lastPseudojets () const;

  /**
   * Return the last cluster sequence obtained.
   */
  inline const fastjet::ClusterSequence& lastClusterSequence () const;

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
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
  //@}


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


protected:

  /** @name Standard Interfaced functions. */
  //@{

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

  //@}

private:

  /**
   * The jet finder
   */
  int _jetFinder;

  /**
   * The strategy used for jet finding
   */
  int _strategy;

  /**
   * The recombination scheme
   */
  int _recombinationScheme;

  /**
   * The last pseudojet vector
   */
  vector<fastjet::PseudoJet> _lastPseudojets;

  /**
   * The last obtained cluster sequence.
   */
  std::auto_ptr<fastjet::ClusterSequence> _lastClusterSequence;

  /**
   * The fastjet::JetDefinition to be used
   */
  fastjet::JetDefinition _jetDefinition;

  /**
   * The visible energy squared as obtained by
   * squaring the sum of final state momenta.
   */
  Energy2 _E2vis;

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<FastJetFinder> initFastJetFinder;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  FastJetFinder & operator=(const FastJetFinder &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of FastJetFinder. */
template <>
struct BaseClassTrait<Analysis2::FastJetFinder,1> {
  /** Typedef of the first base class of FastJetFinder. */
  typedef Analysis2::JetFinder NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the FastJetFinder class and the shared object where it is defined. */
template <>
struct ClassTraits<Analysis2::FastJetFinder>
  : public ClassTraitsBase<Analysis2::FastJetFinder> {
  /** Return a platform-independent class name */
  static string className() { return "Analysis2::FastJetFinder"; }
  /**
   * The name of a file containing the dynamic library where the class
   * FastJetFinder is implemented. It may also include several, space-separated,
   * libraries if the class FastJetFinder depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "FastJetFinder.so Analysis2.so"; }
};

/** @endcond */

}

#include "FastJetFinder.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "FastJetFinder.tcc"
#endif

#endif /* Analysis2_FastJetFinder_H */
