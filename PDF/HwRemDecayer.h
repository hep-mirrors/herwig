// -*- C++ -*-
#ifndef HERWIG_HwRemDecayer_H
#define HERWIG_HwRemDecayer_H
//
// This is the declaration of the HwRemDecayer class.
//

#include "ThePEG/PDT/RemnantDecayer.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/SubProcess.h"
#include "ForcedSplitting.h"
#include "HwRemDecayer.fh"

namespace Herwig {
using namespace ThePEG;
/**
 * Here is the documentation of the HwRemDecayer class.
 *
 * @see \ref HwRemDecayerInterfaces "The interfaces"
 * defined for HwRemDecayer.
 */
class HwRemDecayer: public RemnantDecayer {

public:

  typedef vector<pair<tPPtr, tPPtr> > PartnerMap;

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline HwRemDecayer();

  /**
   * The copy constructor.
   */
  inline HwRemDecayer(const HwRemDecayer &);

  /**
   * The destructor.
   */
  virtual ~HwRemDecayer();
  //@}

public:

  /** @name Virtual functions required by the Decayer class. */
  //@{
  /**
   * Check if this decayer can perfom the decay specified by the
   * given decay mode.
   * @param dm the DecayMode describing the decay.
   * @return true if this decayer can handle the given mode, otherwise false.
   */
  virtual bool accept(const DecayMode & dm) const;

  /**
   * Return true if this decayer can handle the extraction of the \a   
   * extracted parton from the given \a particle.   
   */  
  virtual bool canHandle(tcPDPtr parent, tcPDPtr extracted) const;  

  /**   
   * Return true if this decayed can extract more than one parton from   
   * a particle.   
   */  
  virtual bool multiCapable() const;

  /**
   * Perform a decay for a given DecayMode and a given Particle instance.
   * @param dm the DecayMode describing the decay.
   * @param p the Particle instance to be decayed.
   * @param step the step we are working on.
   * @return a ParticleVector containing the decay products.
   */
  virtual ParticleVector decay(const DecayMode & dm, const Particle & p, Step & step) const;
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

  /**
   * Do several checks and initialization, for remnantdecay inside ShowerHandler.
   */
  void initialize(pair<tRemPPtr, tRemPPtr> rems, Step & step);

  void doSplit(pair<tPPtr, tPPtr> partons, bool first);

  void finalize();

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
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<HwRemDecayer> initHwRemDecayer;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  HwRemDecayer & operator=(const HwRemDecayer &);

private:
  
  /**                                                                           
   * Simple struct to store info about baryon quark and di-quark                
   * constituents.                                                              
   */                                                                           
  struct HadronContent {                                                        
    /**
     * Randomly choose a valence flavour from \a flav.
     */
    int getValence();

    /**
     * manually extract the valence flavour \a id.
     */
    void extract(int id);

    /**
     * Return a proper particle ID assuming that \a id has been removed
     * from the hadron.
     */
    long RemID() const;

    /**
     * Method to determine whether \a parton is a quark from the sea.
     * @return TRUE if \a parton is neither a valence quark nor a gluon.
     */
    bool isSeaQuark(tcPPtr parton) const;

    /**
     * Method to determine whether \a parton is a valence quark.
     */
    bool isValenceQuark(tcPPtr parton) const;

    /** The valence flavours of the corresponding baryon. */                    
    vector<int> flav;                                                           

    /** The array index of the extracted particle. */
    int extracted;

    /** -1 if the particle is an anti-particle. +1 otherwise. */                
    int sign;                                                                   
  }; 

  /**
   * Return a HadronContent struct from a PPtr to a hadron.
   */
  HadronContent getHadronContent(tcPPtr hadron) const;

  /**
   * Do the forced Splitting of the Remnant with respect to the 
   * extracted parton \a parton.
   * @param parton = PPtr to the parton going into the subprocess
   * @param content = HadronContent struct to keep track of flavours.
   * @param used = Momentum vector to keep track of remaining momenta.
   * @param partners = vector of pairs filled with tPPtr to the particles 
   * which should be colour connected.
   */
  void split(tPPtr parton, HadronContent & content, tRemPPtr rem, 
	     Lorentz5Momentum & used, PartnerMap & partners, bool first);

  /**
   * Do all colour connections.
   * @param partners = Object that holds the information which particles to connect.
   * @param anti = flag to indicate, if (anti)colour was extracted as first parton.
   */
  void fixColours(PartnerMap partners, bool anti) const;

  /**
   * Set the momenta of the Remnants properly and boost the decay particles.
   */
  void setRemMasses() const;

  /**
   * This is a poniter to the Herwig::ForcedSplitting object
   */
  ForcedSplittingPtr theForcedSplitter; 

  /**
   * A flag which indicates, whether the extracted valence quark was a 
   * anti particle.
   */
  pair<bool, bool> theanti;

  /**
   * variable to sum up the x values of the extracted particles
   */
  pair<double, double> theX;

  /**Pair of HadronContent structs to know about the quark content of the beams*/
  pair<HadronContent, HadronContent> theContent;

  /**Pair of Lorentz5Momentum to keep track of the forced splitting product momenta*/
  pair<Lorentz5Momentum, Lorentz5Momentum> theUsed;

  /**
   * Pair of PartnerMap's to store the particles, which will be colour
   * connected in the end.
   */
  pair<PartnerMap, PartnerMap> theMaps;

  /**
   * Variable to hold a pointer to the current step. The variable is used to 
   * determine, wether decay(const DecayMode & dm, const Particle & p, Step & step) 
   * has been called in this event or not.
   */
  StepPtr thestep;

  /**
   * Pair of Remnant pointers. This is needed to boost
   * in the Remnant-Remnant CMF after all have been decayed.
   */
  pair<RemPPtr, RemPPtr> theRems;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of HwRemDecayer. */
template <>
struct BaseClassTrait<Herwig::HwRemDecayer,1> {
  /** Typedef of the first base class of HwRemDecayer. */
  typedef RemnantDecayer NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the HwRemDecayer class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::HwRemDecayer>
  : public ClassTraitsBase<Herwig::HwRemDecayer> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::HwRemDecayer"; }
  /**
   * The name of a file containing the dynamic library where the class
   * HwRemDecayer is implemented. It may also include several, space-separated,
   * libraries if the class HwRemDecayer depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwRemDecayer.so"; }
};

/** @endcond */

}

#include "HwRemDecayer.icc"
#ifndef HERWIG_TEMPLATES_IN_CC_FILE
// #include "HwRemDecayer.tcc"
#endif

#endif /* HERWIG_HwRemDecayer_H */
