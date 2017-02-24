// -*- C++ -*-
#ifndef HERWIG_AlpGenHandlerOL_H
#define HERWIG_AlpGenHandlerOL_H
//
// This is the declaration of the AlpGenHandlerOL class.
//

#include "Herwig/Shower/ShowerHandler.h"
#include "Herwig/Shower/QTilde/QTildeShowerHandler.h"
#include "ThePEG/Config/Pointers.h"
#include "Herwig/Shower/Core/Couplings/ShowerAlpha.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "HiggsPair.h"

namespace Herwig {
  class AlpGenHandlerOL;
}

//declaration of thepeg ptr
namespace ThePEG {
  ThePEG_DECLARE_POINTERS(Herwig::AlpGenHandlerOL,AlpGenHandlerOLPtr);
}

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the AlpGenHandlerOL class.
 *
 * @see \ref AlpGenHandlerOLInterfaces "The interfaces"
 * defined for AlpGenHandlerOL.
 */
class AlpGenHandlerOL: public QTildeShowerHandler {

public:

  /**
   * The default constructor.
   */
  AlpGenHandlerOL();

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
   * Finalize the object
   */
  virtual void dofinish();

  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();
  //@}

public:
  /**
   * Hook to allow vetoing of event after showering hard sub-process
   * as in e.g. MLM merging.
   */
  virtual bool showerHardProcessVeto() const;


protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}

private:
  
  /* 
   * Run MLM jet-parton matching on the 'extra' jets.
   */
  
  bool lightJetPartonVeto();
  
  /* 
   * Function that calculates deltaR between a parton and a jet
   */

  double partonJetDeltaR(ThePEG::tPPtr partonptr, LorentzMomentum jetmom) const;
  
  /**
   * c++ translation of subroutine of same name from alpsho.f.
   * Initialize calorimeter for calsim_m and getjet_m. Note that
   * because initialization is separte calsim_m can be called more
   * than once to simulate pileup of several events.
   */
  void calini_m() const;
  
 
  /**
   * c++ translation of subroutine of same name from alpsho.f.
   * Simple calorimeter simulation - assume uniform Y and phi bins.
   */
  void calsim_m() const;

 /**
   * Find jets using the FastJet package on particlesToCluster_.
   */
  void getFastJets(double rjet, Energy ejcut, double etajcut) const;
  

  /**
   * c++ translation of subroutine of same name from alpsho.f.
   * Simple jet-finding algorithm (similar to UA1). Find highest
   * remaining cell > ETSTOP and sum surrounding cells with --
   * DELTA(Y)**2+DELTA(PHI)**2 < RJET**2 , ET>ECCUT.
   * Keep sets with ET>EJCUT and ABS(ETA)<ETACUT. The UA1
   * parameters are RJET=1.0 and EJCUT=5.0.
   */
  void getjet_m(double rjet, Energy ejcut, double etajcut) const;

  /**
   * Deletes particles from partonsToMatch_ and particlesToCluster_
   * vectors so that these contain only the partons to match to the
   * jets and the particles used to build jets respectively. By and
   * large the candidates for deletion are: vector bosons and their
   * decay products, Higgs bosons, photons as well as _primary_, i.e.
   * present in the lowest multiplicity process, heavy quarks and
   * any related decay products.
   */
  void caldel_m() const;

  /**  
   * c++ translation of subroutine of same name from alpsho.f.
   * Label all particles with status between ISTLO and ISTHI 
   * (until a particle with status ISTOP is found) as final-state,
   * call calsim_m and then put labels back to normal. This 
   * version keeps only all IST=1 particles rejected by caldel as
   * daughters of vetoed heavy-quark mothers: jets complementary
   * to those reconstructed by caldel.
   */
  void caldel_hvq() const;

  /**
   * Get the particles from lastXCombPtr filling the pair
   * preshowerISPs_ and particle pointer vector preshowerFSPs_.
   */
  void getPreshowerParticles() const;

  /**
   * Get the particles from eventHandler()->currentEvent()->...
   * filling the particle pairs showeredISHs_, showeredISPs_,
   * showeredRems_ and the particle pointer vector showeredFSPs_.
   */
  void getShoweredParticles() const;

  /**
   * Allows printing of debug output and sanity checks like
   * total momentum consrvation to be carried out.
   * debugLevel = -1, 0, ...5 
   *            = no debugging, minimal debugging, ... verbose.
   */
  void doSanityChecks(int debugLevel) const;

  /**
   * Given a pointer to a particle this finds all its final state
   * descendents.
   */
  void getDescendents(PPtr theParticle) const;

  /**
   * Accumulates all descendents of tops down to the b and W
   * but not including them.
   */
  void getTopRadiation(PPtr theParticle) const;

  /** 
   * Sorts a given vector of particles by descending pT or ETJET
   */
  
  ParticleVector pTsort(ParticleVector unsortedVec);
  pair< vector<Energy>, vector<Lorentz5Momentum> > ETsort(vector<Energy> unsortedetjet, vector<Lorentz5Momentum> unsortedVec);

  /* 
   * A function that prints a vector of Lorentz5Momenta in a fancy way
   */
  void printMomVec(vector<Lorentz5Momentum> momVec);


  /* 
   * A probability function for varying etclus_ about the mean value
   */
  Energy etclusran_(double petc) const;

private:
  
  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<AlpGenHandlerOL> initAlpGenHandlerOL;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  AlpGenHandlerOL & operator=(const AlpGenHandlerOL &);

private:

  /**
   *  Initial-state incoming partons prior to showering
   *  (i.e. from lastXCombPtr).
   */
  mutable PPair preshowerISPs_;

  /**
   *  Final-state outgoing partICLEs prior to showering
   *  (i.e. from lastXCombPtr).
   */
  mutable ParticleVector preshowerFSPs_;

  /**
   *  Final-state outgoing partICLEs prior to showering _to_be_removed_
   *  from preShowerFSPs_ prior to the light-parton-light-jet matching
   *  step. This same list is the starting point for determining 
   *  partonsToMatch_ for the case of merging in heavy quark production.
   */
  mutable ParticleVector preshowerFSPsToDelete_;

  /**
   *  Initial-state incoming hadrons after shower of hard process
   *  (eventHandler()->currentEvent()->incoming()).
   */
  mutable PPair showeredISHs_;

  /**
   *  Initial-state incoming partons after shower of hard process
   *  (look for partonic children of showeredISHs_).
   */
  mutable PPair showeredISPs_;

  /**
   *  Final-state outgoing partICLEs after shower of hard process
   *  (eventHandler()->currentEvent()->getFinalState()).
   */
  mutable tPVector showeredFSPs_;

  /**
   *  Final-state outgoing partICLEs after shower of hard process
   *  _to_be_removed_ from showeredFSPs_ prior to the
   *  light-parton-light-jet matching step. This same list is the
   *  starting point for determining particlesToCluster_ for the
   *  case of merging in heavy quark production.
   */
  mutable ParticleVector showeredFSPsToDelete_;

  /**
   *  ONLY the final-state partons from preshowerFSPs_ that are
   *  supposed to enter the jet-parton matching.
   */
  mutable ParticleVector partonsToMatch_;
 
  /*
   * The shower progenitors
   */

  mutable PPtr theProgenitor;
  mutable PPtr theLastProgenitor;

  /**
   *  ONLY the final-state particles from showeredFSPs_ (and maybe
   *  also showeredRems_) that are supposed to go for jet clustering.
   */
  mutable tPVector particlesToCluster_;

  /**
   *  Final-state remnants after shower of hard process
   *  (look for remnants initially in showeredFSPs_).
   */
  mutable PPair showeredRems_;

  /**
   *  Pointer to the object calculating the strong coupling
   */
  mutable ShowerAlphaPtr alphaS_;

  /**
   *  Information extracted from the XComb object
   */
  //@{
  /**
   * The fixed factorization scale used in the MEs.
   */
  mutable Energy pdfScale_;

  /**
   *  Centre of mass energy
   */
  mutable Energy2 sHat_;

  /**
   * Constant alphaS used to generate LH events - if not already
   * using CKKW scale (ickkw = 1 in AlpGen for example).
   */
  mutable double alphaSME_;
  //@}

  /*
   * Number of rapidity segments of the calorimeter.
   */
  mutable unsigned int ncy_;

  /*
   * Number of phi segments of the calorimeter.
   */
  mutable unsigned int ncphi_;

  /*
   * Heavy flavour in WQQ,ZQQ,2Q etc (4=c, 5=b, 6=t).
   */
  mutable int ihvy_;

  /*
   * Number of photons in the AlpGen process.
   */
  mutable int nph_;

  /*
   * Number of higgses in the AlpGen process.
   */
  mutable int nh_;

  /*
   * Jet ET cut to apply in jet clustering (in merging).
   */
  mutable Energy etclus_;

  /*
   * Mean Jet ET cut to apply in jet clustering (in merging).
   */
  mutable Energy etclusmean_;

  /*
   * maximum deviation from mean Jet ET cut to apply in jet clustering (in merging).
   */
  mutable Energy epsetclus_;

  /*
   * type of smoothing function to use
   */
  mutable unsigned int smoothingtype_;





  /*
   * Cone size used in jet clustering (in merging).
   */
  mutable double rclus_;

  /*
   * Max |eta| for jets in clustering (in merging).
   */
  double etaclmax_;

  /*
   * Default 1.5 factor used to decide if a jet matches a parton
   * in merging: if DR(parton,jet)<rclusfactor*rclus the parton
   * and jet are said to have been matched.
   */
  double rclusfactor_;

  /*
   * The AlpGen hard process code. Relation to the AlpGen process names:
   * 1: wqq, 2: zqq, 3: wjet, 4: zjet, 5: vbjet, 6: 2Q, 8: QQh, 9: Njet, 
   * 10: wcjet, 11: phjet, 12: hjet, 13: top, 14: wphjet, 15: wphqq, 
   * 16: 2Qph.
   */
  mutable int ihrd_;

  /*
   * The number of light jets in the AlpGen process (i.e. the 'extra' ones).
   */
  mutable int njets_;

  /*
   * Mimimum parton-parton R-sep used for generation (used for hvq merging).
   */
  mutable double drjmin_;

  /*
   * This flags that the highest multiplicity ME-level process is
   * being processed.
   */
  mutable bool highestMultiplicity_;

  /*
   * This is the highest number of jets to be included in the matching.
   * This implies that exclusive rates will be calculated up to highestNjets_-1 and 
   * inclusive for highestNjets_.
   */
  mutable int highestNjets_; 

  /*
   * This flags that the highest NLO multiplicity ME-level process is
   * being processed.
   */
  mutable bool highestNLOMultiplicity_;
  
  /*
   * This flags whether the etclus_ (merging scale) should be fixed or variable according to a prob. distribution around the mean
   */
  mutable bool etclusfixed_;

  /*
   * The forwards rapidity span of the calorimeter.
   */
  mutable double ycmax_;

  /*
   * The backwards rapidity span of the calorimeter.
   */
  mutable double ycmin_;

  /*
   * The jet algorithm used for parton-jet matching in the MLM procedure.
   */
  mutable int jetAlgorithm_;

  /*
   * Allows the vetoing to be turned off completely - just for convenience.
   */
  mutable bool vetoIsTurnedOff_;

  /*
   * Switch between original and OpenLoops implementations of the handler
   * 
   */
  mutable int vetoType_;

  /*
   * Signals that the LH file being read-in is a NLO (Powheg one).
   */
  mutable bool inputIsNLO_;

  /*
   * Cosine of phi values of calorimeter cell centres.
   * Goes phi~=0 to phi~=2*pi          (index = 0 ---> ncphi).
   * ==> Cosine goes from +1 ---> +1   (index = 0 ---> ncphi).
   */
  mutable vector<double> cphcal_;

  /*
   * Sine of phi values of calorimeter cell centres.
   * Goes phi~=0 to phi~=2*pi             (index = 0 ---> ncphi).
   * ==> Sine goes 0 -> 1 -> 0 -> -1 -> 0 (index = 0 ---> ncphi).
   */
  mutable vector<double> sphcal_;

  /*
   * Cosine of theta values of calorimeter cell centres in Y. 
   * Goes bwds th~=pi to fwds th~=0       (index = 0 ---> ncy).
   * ==> Cosine goes from -1 ---> +1      (index = 0 ---> ncy).
   */
  mutable vector<double> cthcal_;

  /*
   * Sine of theta values of calorimeter cell centres in Y. 
   * Goes bwds th~=pi to fwds th~=0       (index = 0 ---> ncy).
   * ==> Sine goes from  0 ---> +1 ---> 0 (index = 0 ---> ncy).
   */
  mutable vector<double> sthcal_;

  /*
   * Transverse energy deposit in a given calorimeter cell.
   * First array index corresponds to rapidity index of cell,
   * second array index corresponds to phi cell index.
   */
  mutable vector<vector<Energy> > et_;

  /*
   * For a given calorimeter cell this holds the index of the jet
   * that the cell was clustered into.
   */
  mutable vector<vector<int> > jetIdx_;

  /*
   * Vector holding the Lorentz 5 momenta of each jet.
   */
  mutable vector<Lorentz5Momentum> pjet_;

  /*
   * Vector holding the list of FS particles resulting from 
   * the particle input to getDescendents.
   */
  mutable ParticleVector tmpList_;

  /*
   * Variables for the C++ translation of the calini_m(), calsim_m(),
   * getjet_m(...) and caldel_m() functions 
   */
  mutable vector<Energy> etjet_;
  mutable double dely_, delphi_;
  
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of AlpGenHandlerOL. */
template <>
struct BaseClassTrait<Herwig::AlpGenHandlerOL,1> {
  /** Typedef of the first base class of AlpGenHandlerOL. */
  typedef Herwig::QTildeShowerHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the AlpGenHandlerOL class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::AlpGenHandlerOL>
  : public ClassTraitsBase<Herwig::AlpGenHandlerOL> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::AlpGenHandlerOL"; }
  /**
   * The name of a file containing the dynamic library where the class
   * AlpGenHandlerOL is implemented. It may also include several, space-separated,
   * libraries if the class AlpGenHandlerOL depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "AlpGenHandlerOL.so"; }
};

/** @endcond */

}

#endif /* HERWIG_AlpGenHandlerOL_H */
