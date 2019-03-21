// -*- C++ -*-
//
// HwRemDecayer.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_HwRemDecayer_H
#define HERWIG_HwRemDecayer_H
//
// This is the declaration of the HwRemDecayer class.
//

#include "ThePEG/PDT/RemnantDecayer.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/SubProcess.h"
#include "ThePEG/PDF/BeamParticleData.h"
#include "Herwig/Shower/Core/Couplings/ShowerAlpha.h"
#include "Herwig/PDT/StandardMatchers.h"
#include "ThePEG/PDT/StandardMatchers.h"
#include "HwRemDecayer.fh"

namespace Herwig {
using namespace ThePEG;
/**
 * The HwRemDecayer class is responsible for the decay of the remnants. Additional 
 * secondary scatters have to be evolved backwards to a gluon, the
 * first/hard interaction has to be evolved back to a valence quark.
 * This is all generated inside this class,
 * which main methods are then called by the ShowerHandler.
 *
 * A simple forced splitting algorithm is used.
 * This takes the Remnant object produced from the PDF and backward
 * evolution (hadron - parton) and produce partons with the remaining 
 * flavours and with the correct colour connections.
 *
 * The algorithim operates by starting with the parton which enters the hard process.
 * If this is from the sea there is a forced branching to produce the antiparticle
 * from a gluon branching. If the parton entering the hard process was a gluon, or
 * a gluon was produced from the first step of the algorithm, there is then a further
 * branching back to a valence parton. After these partons have been produced a quark or
 * diquark is produced to give the remaining valence content of the incoming hadron.
 *
 * The forced branching are generated using a scale between QSpac and EmissionRange times
 * the minimum scale. The energy fractions are then distributed using
 * \f[\frac{\alpha_S}{2\pi}\frac{P(z)}{z}f(x/z,\tilde{q})\f]
 * with the massless splitting functions.
 *
 * \author Manuel B\"ahr
 *
 * @see \ref HwRemDecayerInterfaces "The interfaces"
 * defined for HwRemDecayer.
 */
class HwRemDecayer: public RemnantDecayer {

public:

  /** Typedef to store information about colour partners */
  typedef vector<pair<tPPtr, tPPtr> > PartnerMap;

public:

  /**
   * The default constructor.
   */
  HwRemDecayer() : allowTop_(false), multiPeriph_(false), quarkPair_(false),
                   ptmin_(-1.*GeV), beta_(ZERO),
		   maxtrySoft_(10), 
		   colourDisrupt_(1.0),
		   ladderbFactor_(0.0),
		   ladderPower_(-0.08),
		   ladderNorm_(1.0),
		   gaussWidth_(0.1),
		   valOfN_(0), 
		   initTotRap_(0),
		   _kinCutoff(0.75*GeV), 
		   _forcedSplitScale(2.5*GeV),
		   _range(1.1), _zbin(0.05),_ybin(0.),
		   _nbinmax(100), DISRemnantOpt_(0),
		   pomeronStructure_(0), mg_(ZERO) {}

  /** @name Virtual functions required by the Decayer class. */
  //@{
  /**
   * Check if this decayer can perfom the decay specified by the
   * given decay mode.
   * @return true if this decayer can handle the given mode, otherwise false.
   */
  virtual bool accept(const DecayMode &) const {
    return true;
  }

  /**
   * Return true if this decayer can handle the extraction of the \a   
   * extracted parton from the given \a particle.   
   */  
  virtual bool canHandle(tcPDPtr particle, tcPDPtr parton) const;
  
  /**   
   * Return true if this decayed can extract more than one parton from   
   * a particle.   
   */  
  virtual bool multiCapable() const {  
    return true;
  }
  
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

  /** 
   * struct that is used to catch exceptions which are thrown
   * due to energy conservation issues of additional soft scatters
   */
  struct ExtraSoftScatterVeto {};

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
  void initialize(pair<tRemPPtr, tRemPPtr> rems, tPPair beam, Step & step,
		  Energy forcedSplitScale);

  /**
   * Initialize the soft scattering machinery.
   * @param ptmin = the pt cutoff used in the UE model
   * @param beta = slope of the soft pt-spectrum
   */
  void initSoftInteractions(Energy ptmin, InvEnergy2 beta);

  /**
   * Perform the acual forced splitting.
   * @param partons is a pair of ThePEG::Particle pointers which store the final 
   * partons on which the shower ends.
   * @param pdfs are pointers to the pdf objects for both beams
   * @param first is a flage wether or not this is the first or a secondary interation
   */
  void doSplit(pair<tPPtr, tPPtr> partons, pair<tcPDFPtr, tcPDFPtr> pdfs, bool first);

  /**
   * Perform the final creation of the diquarks. Set the remnant masses and do 
   * all colour connections.
   * @param colourDisrupt = variable to control how many "hard" scatters
   * are colour isolated
   * @param softInt = parameter for the number of soft scatters
   */
  void finalize(double colourDisrupt=0.0, unsigned int softInt=0);

  /**
   *  Find the children
   */
  void findChildren(tPPtr,vector<PPtr> &) const;

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const {return new_ptr(*this);}
  //@}

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit() {
    Interfaced::doinit();
    _ybin=0.25/_zbin;
    mg_ = getParticleData(ParticleID::g)->constituentMass();
  }
  //@}

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
  HwRemDecayer & operator=(const HwRemDecayer &) = delete;

public:
  
  /**                                                                           
   * Simple struct to store info about baryon quark and di-quark                
   * constituents.                                                              
   */                                                                           
  struct HadronContent {

    /**
     * manually extract the valence flavour \a id.
     */
    inline void extract(int id) {
      for(unsigned int i=0; i<flav.size(); i++) {
	if(id == sign*flav[i]){
	  if(hadron->id() == ParticleID::gamma || 
	     (hadron->id() == ParticleID::pomeron && pomeronStructure==1) ||
	     hadron->id() == ParticleID::reggeon) {
	    flav[0] =  id;
	    flav[1] = -id;
	    extracted = 0;
	    flav.resize(2);
	  }
	  else if (hadron->id() == ParticleID::pomeron && pomeronStructure==0) {
	    extracted = 0;
	  }
	  else {
	    extracted = i;
	  }
	  break;
	}
      }      
    }
    
    /**
     * Return a proper particle ID assuming that \a id has been removed
     * from the hadron.
     */
    long RemID() const;

    /**
     * Method to determine whether \a parton is a quark from the sea.
     * @return TRUE if \a parton is neither a valence quark nor a gluon.
     */
    bool isSeaQuark(tcPPtr parton) const {
      return ((parton->id() != ParticleID::g) && ( !isValenceQuark(parton) ) );
    }

    /**
     * Method to determine whether \a parton is a valence quark.
     */
    bool isValenceQuark(tcPPtr parton) const {
      return isValenceQuark(parton->id());
    }

    /**
     * Method to determine whether \a parton is a quark from the sea.
     * @return TRUE if \a parton is neither a valence quark nor a gluon.
     */
    bool isSeaQuarkData(tcPDPtr partonData) const {
      return ((partonData->id() != ParticleID::g) && ( !isValenceQuarkData(partonData) ) );
    }

    /**
     * Method to determine whether \a parton is a valence quark.
     */
    bool isValenceQuarkData(tcPDPtr partonData) const {
      int id(sign*partonData->id());
      return find(flav.begin(),flav.end(),id) != flav.end();
    }

    /**
     * Method to determine whether \a parton is a valence quark.
     */
    bool isValenceQuark(int id) const {
      return find(flav.begin(),flav.end(),sign*id) != flav.end();
    }

    /** The valence flavours of the corresponding baryon. */                    
    vector<int> flav;                                                           

    /** The array index of the extracted particle. */
    int extracted;

    /** -1 if the particle is an anti-particle. +1 otherwise. */                
    int sign;

    /** The ParticleData objects of the hadron */
    tcPDPtr hadron;

    /** Pomeron treatment */
    unsigned int pomeronStructure;
  }; 

  /**
   * Return the hadron content objects for the incoming particles.
   */
  const pair<HadronContent, HadronContent>& content() const {
    return theContent;
  }

  /**
   * Return a HadronContent struct from a PPtr to a hadron.
   */
  HadronContent getHadronContent(tcPPtr hadron) const;

  /**
   * Set the hadron contents.
   */
  void setHadronContent(tPPair beam) {
    theContent.first  = getHadronContent(beam.first);
    theContent.second = getHadronContent(beam.second);
  }

private:

  /**
   * Do the forced Splitting of the Remnant with respect to the 
   * extracted parton \a parton.
   * @param parton = PPtr to the parton going into the subprocess.
   * @param content = HadronContent struct to keep track of flavours.
   * @param rem = Pointer to the ThePEG::RemnantParticle.
   * @param used = Momentum vector to keep track of remaining momenta.
   * @param partners = Vector of pairs filled with tPPtr to the particles 
   * which should be colour connected.
   * @param pdf pointer to the PDF Object which is used for this particle
   * @param first = Flag for the first interaction.
   */
  void split(tPPtr parton, HadronContent & content, tRemPPtr rem, 
	     Lorentz5Momentum & used, PartnerMap & partners, tcPDFPtr pdf, bool first);

  /**
   * Merge the colour lines of two particles
   * @param p1 = Pointer to particle 1
   * @param p2 = Pointer to particle 2
   * @param anti = flag to indicate, if (anti)colour was extracted as first parton.
   */
  void mergeColour(tPPtr p1, tPPtr p2, bool anti) const;

  /**
   * Set the colour connections.
   * @param partners = Object that holds the information which particles to connect.
   * @param anti = flag to indicate, if (anti)colour was extracted as first parton.
   * @param disrupt parameter for disruption of the colour structure
   */
  void fixColours(PartnerMap partners, bool anti, double disrupt) const;

  /**
   * Set the momenta of the Remnants properly and boost the decay particles.
   */
  void setRemMasses() const;

  /**
   * This creates a parton from the remaining flavours of the hadron. The
   * last parton used was a valance parton, so only 2 (or 1, if meson) flavours
   * remain to be used.
   */
  PPtr finalSplit(const tRemPPtr rem, long remID,
		  Lorentz5Momentum usedMomentum) const {
    // Create the remnant and set its momentum, also reset all of the decay 
    // products from the hadron
    PPtr remnant = new_ptr(Particle(getParticleData(remID)));
    Lorentz5Momentum prem(rem->momentum()-usedMomentum);
    prem.setMass(getParticleData(remID)->constituentMass());
    prem.rescaleEnergy();
    remnant->set5Momentum(prem);
    // Add the remnant to the step, but don't do colour connections
    thestep->addDecayProduct(rem,remnant,false);
    return remnant;
  }
  

  /**
   * This takes the particle and find a splitting for np -> p + child and 
   * creates the correct kinematics and connects for such a split. This
   * Splitting has an upper bound on qtilde given by the energy argument
   * @param rem The Remnant
   * @param child The PDG code for the outgoing particle
   * @param oldQ  The maximum scale for the evolution
   * @param oldx  The fraction of the hadron's momentum carried by the last parton
   * @param pf    The momentum of the last parton at input and after branching at output
   * @param p     The total emitted momentum
   * @param content The content of the hadron
   */
  PPtr forceSplit(const tRemPPtr rem, long child, Energy &oldQ, double &oldx, 
		  Lorentz5Momentum &pf, Lorentz5Momentum &p,
		  HadronContent & content) const;

  /**
   *  Check if a particle is a parton from a hadron or not
   * @param parton The parton to be tested
   */
  bool isPartonic(tPPtr parton) const;

  /** @name Soft interaction methods. */
  //@{

  /**
   * Produce pt values according to dN/dp_T = N p_T exp(-beta_*p_T^2)
   */
  Energy softPt() const;

  /**
   * Get the 2 pairs of 5Momenta for the scattering. Needs calling of
   * initSoftInteractions.
   */
  void softKinematics(Lorentz5Momentum &r1, Lorentz5Momentum &r2, 
		      Lorentz5Momentum &g1, Lorentz5Momentum &g2) const;

  /**
   * Create N soft gluon interactions
   */
  void doSoftInteractions(unsigned int N){
  	if(!multiPeriph_){
  		doSoftInteractions_old(N);}
  	else{
  		doSoftInteractions_multiPeriph(N);
  	}
  }
  
  /**
   * Create N soft gluon interactions (old version)
   */
  void doSoftInteractions_old(unsigned int N);
  
  /**
   * Create N soft gluon interactions - multiperhpheral kinematics
   */
  void doSoftInteractions_multiPeriph(unsigned int N);

  /**
   * Method to add a particle to the step
   * @param parent = pointer to the parent particle
   * @param id = Particle ID of the newly created particle
   * @param p = Lorentz5Momentum of the new particle
   */
  tPPtr addParticle(tcPPtr parent, long id, Lorentz5Momentum p) const;
  //@}

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

  /**
   *  The beam particle data for the current incoming hadron
   */
  mutable tcPPtr theBeam;

  /**
   *  the beam data
   */
  mutable Ptr<BeamParticleData>::const_pointer theBeamData;

  /** 
   *  The PDF for the current initial-state shower 
   */ 
  mutable tcPDFPtr _pdf; 
  
private:

  /**
   *  Switch to control handling of top quarks in proton
   */
  bool allowTop_;
  
  /**
   *  Switch to control using multiperipheral kinemaics
   */
  bool multiPeriph_;
  
  /**
   *  True if kinematics is to be calculated for quarks
   */
  bool quarkPair_;

  /** @name Soft interaction variables. */
  //@{

  /**
   * Pair of soft Remnant pointers, i.e. Diquarks.
   */
  tPPair softRems_;

  /**
   * ptcut of the UE model
   */
  Energy ptmin_;

  /**
   * slope of the soft pt-spectrum: dN/dp_T = N p_T exp(-beta*p_T^2)
   */
  InvEnergy2 beta_;

  /**
   *  Maximum number of attempts for the regeneration of an additional
   *  soft scattering, before the number of scatters is reduced.
   */
  unsigned int maxtrySoft_;

  /**
   * Variable to store the relative number of colour disrupted
   * connections to additional soft subprocesses.
   */
  double colourDisrupt_;
  
  /**
   * Variable to store the additive factor of the 
   multiperipheral ladder multiplicity.
   */
  double ladderbFactor_;
  
  /**
   * Variable of the parameterization of the ladder multiplicity.
   */
  double ladderPower_;

  /**
   * Variable of the parameterization of the ladder multiplicity.
   */
  double ladderNorm_;

  /**
   * Variable to store the gaussian width of the 
   * fluctuation of the longitudinal momentum
   * fraction.
   */
  double gaussWidth_;
  
  /**
   * Variable to store the current total multiplicity 
   of a ladder.
   */
  double valOfN_;
  
  /**
   * Variable to store the initial total rapidity between 
   of the remnants.
   */
  double initTotRap_;

  //@}

  /** @name Forced splitting variables. */
  //@{

  /**
   *  The kinematic cut-off
   */
  Energy _kinCutoff;
  
  /**
   * The PDF freezing scale as set in ShowerHandler
   */
  Energy _forcedSplitScale;

  /**
   *  Range for emission
   */
  double _range;

  /**
   *  Size of the bins in z for the interpolation
   */
  double _zbin;

  /**
   *  Size of the bins in y for the interpolation
   */
  double _ybin;

  /**
   *  Maximum number of bins for the z interpolation
   */
  int _nbinmax;

  /**
   *  Pointer to the object calculating the QCD coupling
   */
  ShowerAlphaPtr _alphaS;

  /**
   *  Pointer to the object calculating the QED coupling
   */
  ShowerAlphaPtr _alphaEM; 

  /**
   *  Option for the DIS remnant
   */
  unsigned int DISRemnantOpt_;

  /**
   *  Option for the treatment of the pomeron structure
   */
  unsigned int pomeronStructure_;
  //@}

  /**
   * The gluon constituent mass.
   */
  Energy mg_;

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
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_HwRemDecayer_H */
