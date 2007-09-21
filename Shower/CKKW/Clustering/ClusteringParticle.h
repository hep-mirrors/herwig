// -*- C++ -*-
#ifndef HERWIG_ClusteringParticle_H
#define HERWIG_ClusteringParticle_H
//
// This is the declaration of the ClusteringParticle class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "ThePEG/Vectors/Lorentz5Vector.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ClusteringParticle.fh"
#include "Clustering.fh"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {

using namespace ThePEG;

  namespace ClusteringParticleState {
    /**\ingroup CKKW
     *
     * Enumeration to identify the state of a
     * parton.
     *
     *@author Simon Plaetzer
     *
     */
    enum ClusteringParticleState {
      undefined = -2,
      initial,
      intermediate,
      final
    };
  }

  namespace ClusteringInteractionType {
    /**\ingroup CKKW
     *
     * Enumeration to identify interaction types.
     *
     *@author Simon Plaetzer
     *
     */
    enum ClusteringInteractionType {
      undefined = 0,
      QCD,
      QED,
      EWK
    };

  }

  /**\ingroup CKKW
   *
   * Simple struct to hold PDG id and state
   * of a parton.
   *
   *@author Simon Plaetzer
   *
   */
  struct PartonId {

    inline PartonId () : PDGId(0), state(ClusteringParticleState::undefined) { }
    inline explicit PartonId (long id, int st) : PDGId(id), state(st) { }
    inline bool operator < (const PartonId& x) const {
      if (PDGId == x.PDGId) return state < x.state;
      return PDGId < x.PDGId;
    }
    inline bool operator == (const PartonId& x) const {
      return PDGId == x.PDGId && state == x.state;
    }

    long PDGId;
    int state;
  };

  inline PersistentOStream& operator << (PersistentOStream& os, const PartonId&);
  inline PersistentIStream& operator >> (PersistentIStream& is, PartonId&);

  /**\ingroup CKKW
   *
   * ClusteringParticleData is a struct which is used internally
   * for reconstructing a CascadeHistory, storing all relevant
   * data for cascade reconstruction.
   *
   *@author Simon Plaetzer
   *
   */
  struct ClusteringParticleData {

    /**
     * Default constructor
     */
    inline ClusteringParticleData ()
      : partonId (),
	colour(0), antiColour(0) { }

    /**
     * Lexicographic ordering
     * of vectors of ClusteringParticleData
     */
    inline bool operator < (const ClusteringParticleData& x) const {
      if (partonId == x.partonId) {
	if (colour == x.colour) return antiColour < x.antiColour;
	return colour < x.colour;
      }
      return partonId < x.partonId;
    }

    inline bool operator == (const ClusteringParticleData& x) const {
      return colour == x.colour && antiColour == x.antiColour
	&& partonId.PDGId == x.partonId.PDGId && partonId.state == x.partonId.state;
    }

    /**
     * The PDG id of the particle
     */
    PartonId partonId;

    /**
     * The colour-line index this particle
     * is connected to.
     */
    int colour;

    /**
     * The anticolour-line index this particle
     * is connected to.
     */
    int antiColour;

  };


/**\ingroup CKKW
 *
 * ClusteringParticle is a simplified particle class
 * storing all relevant information for particles
 * being handled within a cascade history.
 *
 *@author Simon Plaetzer
 *
 * @see \ref ClusteringParticleInterfaces "The interfaces"
 * defined for ClusteringParticle.
 */
class ClusteringParticle: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline ClusteringParticle();

  /**
   * Constructor giving particle data, momentum and
   * initial state momentum fraction.
   */
  inline ClusteringParticle (const ClusteringParticleData&,
			     const Lorentz5Momentum& mom = Lorentz5Momentum(),
			     double x = 0.);

  /**
   * The destructor.
   */
  virtual ~ClusteringParticle();
  //@}

public:

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

public:
  /**@name Methods used during clustering. */
  //@{

  /**
   * Indicate that this particle has been clustered.
   * The momentum and x are backuped, but no change
   * on indices is done.
   */
  inline void wasClustered (tClusteringPtr);

  /**
   * Indicate that this particle emerged from a clustering
   * giving the new index and the clustering.
   */
  inline void emergedFromClustering (unsigned int, tClusteringPtr);

  /**
   * Indicate that a clustering this particle belongs
   * to has been undone. Momentum and x are restored.
   */
  inline void wasUnclustered ();

  /**
   * Indicate that a clustering has been performed giving
   * the index of this particle after the clustering.
   * The index is pushed to _indexStack and the current
   * momentum and x are backuped.
   */
  inline void performedClustering (unsigned int);

  /**
   * Indicate that a clustering has been undone.
   * Momentum, x and index are restored.
   */
  inline void undoneClustering ();

  /**
   * Return the current index of this particle
   */
  inline unsigned int index () const;

  /**
   * Return the production Clustering of this particle.
   */
  inline tcClusteringPtr production () const;

  /**
   * Return the splitting Clustering of this particle.
   */
  inline tcClusteringPtr splitting () const;

  //@}

public:
  /**@ Access to parents and children */
  //@{

  /**
   * Get the children of this particle.
   */
  vector<tClusteringParticlePtr> children () const;

  /**
   * Get the parents of this particle.
   */
  vector<ClusteringParticlePtr> parents () const;

  //@}

public:
  /**@ Access to data and kinematics */
  //@{

  /**
   * Access the ClusteringParticleData object
   */
  inline ClusteringParticleData& pData ();

  /**
   * Get the momentum
   */
  inline Lorentz5Momentum momentum () const;

  /**
   * Set the momentum
   */
  inline void momentum (const Lorentz5Momentum&);

  /**
   * Get the momentum fraction
   */
  inline double x () const;

  /**
   * Set the momentum fraction
   */
  inline void x (double);

  //@}

public:
  /**@ Access to evolution scales */
  //@{

  /**
   * Get the production scale.
   */
  inline Energy2 productionScale () const;

  /**
   * Set the production scale.
   */
  inline void productionScale (const Energy2&);

  /**
   * Get the splitting scale.
   */
  inline Energy2 splittingScale () const;

  /**
   * Set the splitting scale.
   */
  inline void splittingScale (const Energy2&);

  /**
   * Set the showering scale for this particle.
   */
  inline Energy2 showerScale () const;

  /**
   * Get the showering scale for this particle.
   */
  inline void showerScale (const Energy2&);

  /**
   * Return true, if no weights should be associated
   * with this particle
   */
  inline bool noReweight() const;

  /**
   * Indicate that no weights should be associated
   * with this particle
   */
  inline void setNoReweight();

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


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).

#ifdef HERWIG_DEBUG_CKKW

public:

  void debugDump (ostream&);

#endif

public:

private:

  /**
   * The ClusteringParticleData object
   */
  ClusteringParticleData _data;

  /**
   * The momentum of this particle in the current clustering
   * step.
   */
  Lorentz5Momentum _momentum;

  /**
   * For initial state partons the momentum fraction in the
   * current clustering step.
   */
  double _x;

  /**
   * The evolution scale associated with production
   * of this particle.
   */
  Energy2 _productionScale;

  /**
   * The evolution scale associated with splitting
   * of this particle.
   */
  Energy2 _splittingScale;

  /**
   * The showering scale associated
   * with this particle.
   */
  Energy2 _showerScale;

  /**
   * A pointer to the clustering, where this
   * particle is considered to be produced.
   */
  tClusteringPtr _production;

  /**
   * A pointer to the clustering where this
   * particle is considered to split.
   */
  tClusteringPtr _splitting;

  /**
   * Wether or not to associate weights
   * with this particle.
   */
  bool _noReweight;

  /**
   * The index stack.
   */
  list<unsigned int> _indexStack;

  /**
   * Momentum backups.
   */
  list<Lorentz5Momentum> _momentumBackup;

  /**
   * x backups.
   */
  list<double> _xBackup;

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an concrete class without persistent data.
   */
  static NoPIOClassDescription<ClusteringParticle> initClusteringParticle;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ClusteringParticle & operator=(const ClusteringParticle &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of ClusteringParticle. */
template <>
struct BaseClassTrait<Herwig::ClusteringParticle,1> {
  /** Typedef of the first base class of ClusteringParticle. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the ClusteringParticle class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::ClusteringParticle>
  : public ClassTraitsBase<Herwig::ClusteringParticle> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::ClusteringParticle"; }
  /**
   * The name of a file containing the dynamic library where the class
   * ClusteringParticle is implemented. It may also include several, space-separated,
   * libraries if the class ClusteringParticle depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#include "ClusteringParticle.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ClusteringParticle.tcc"
#endif

#endif /* HERWIG_ClusteringParticle_H */
