// -*- C++ -*-

// (C) 2007-2009 Simon Plaetzer -- sp@particle.uni-karlsruhe.de

#ifndef Analysis2_KtJetFinder_H
#define Analysis2_KtJetFinder_H
//
// This is the declaration of the KtJetFinder class.
//

#include "JetFinder.h"
#include "KtJetFinder.fh"

#include "KtJet/KtEvent.h"

#include <memory>

namespace Analysis2 {

using namespace ThePEG;

/**\ingroup Analysis2
 *
 * Interface to KtJet
 *
 * @author Simon Plaetzer
 *
 * @see \ref KtJetFinderInterfaces "The interfaces"
 * defined for KtJetFinder.
 */
class KtJetFinder: public JetFinder {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline KtJetFinder();

  /**
   * The copy constructor.
   */
  inline KtJetFinder(const KtJetFinder &);

  /**
   * The destructor.
   */
  virtual ~KtJetFinder();
  //@}

public:

  /**
   * Set the event to be analysed.
   * The implementation of this method may
   * also do some conversion.
   */
  virtual void use (const vector<Lorentz5Momentum>&, bool);

  /**@name Jet finding */
  //@{

  /**
   * Do inclusive jet finding
   */
  virtual inline void findJets ();

  /**
   * Do exclusive jet-finding for nJets jets
   */
  virtual inline void findJetsN (unsigned int nJets);

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

  /**@name KtJet options */
  //@{

  /**
   * Get the collision type
   */
  inline int collisionType () const;

  /**
   * Get the distance scheme
   */
  inline int distanceScheme () const;

  /**
   * Get the recombination scheme
   */
  inline int recombinationScheme () const;

  //@}

protected:

  /**
   * Convert event to KtLorentzVector vector
   */
  void convert ();

  /**
   * Convert KtLorenztVector vector to vector
   * of Lorentz5Momenta
   */
  void convert (const vector<KtJet::KtLorentzVector>&);

  /**
   * Perform clustering in exclusive mode.
   */
  void cluster ();

  /**
   * Perform clustering in inclusive mode.
   */
  void clusterInclusive ();

  /**
   * Return the vector of KtLorenztVector's obtained
   * from the last event.
   */
  inline const vector<KtJet::KtLorentzVector>& lastMomenta () const;

  /**
   * Return the last KtEvent obtained.
   */
  inline const KtJet::KtEvent& lastKtEvent () const;

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


private:

  /**
   * The collision type
   */
  int _collisionType;

  /**
   * The distance scheme
   */
  int _distanceScheme;

  /**
   * The recombination scheme
   */
  int _recombinationScheme;

  /**
   * KtLorentzVector's from last event
   */
  vector<KtJet::KtLorentzVector> _lastMomenta;

  /**
   * Last KtEvent obtained
   */
  std::auto_ptr<KtJet::KtEvent> _lastKtEvent;

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<KtJetFinder> initKtJetFinder;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  KtJetFinder & operator=(const KtJetFinder &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of KtJetFinder. */
template <>
struct BaseClassTrait<Analysis2::KtJetFinder,1> {
  /** Typedef of the first base class of KtJetFinder. */
  typedef Analysis2::JetFinder NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the KtJetFinder class and the shared object where it is defined. */
template <>
struct ClassTraits<Analysis2::KtJetFinder>
  : public ClassTraitsBase<Analysis2::KtJetFinder> {
  /** Return a platform-independent class name */
  static string className() { return "Analysis2::KtJetFinder"; }
  /**
   * The name of a file containing the dynamic library where the class
   * KtJetFinder is implemented. It may also include several, space-separated,
   * libraries if the class KtJetFinder depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "KtJetFinder.so Analysis2.so"; }
};

/** @endcond */

}

#include "KtJetFinder.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "KtJetFinder.tcc"
#endif

#endif /* Analysis2_KtJetFinder_H */
