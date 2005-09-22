// -*- C++ -*-
#ifndef HERWIG_SplittingGenerator_H
#define HERWIG_SplittingGenerator_H
//
// This is the declaration of the SplittingGenerator class.

#include "ThePEG/Handlers/HandlerBase.h"
#include "ShowerConfig.h"
#include "Herwig++/Utilities/GlobalParameters.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "ShowerIndex.h"
#include "SudakovFormFactor.h"
#include "ShowerVariables.h"
#include "ThePEG/Utilities/StringUtils.h"
#include "ShowerKinematics.h"

namespace Herwig {

using namespace ThePEG;

class ShowerParticle;          // forward declaration

struct Branching { 
   ShoKinPtr first;
   tSudakovPtr second; 
   IdList third; 
   Branching(ShoKinPtr a, tSudakovPtr b, IdList c) {
      first = a; second = b; third = c; 
   }
   Branching() {}
};

/** \ingroup Shower
 *
 *  This class is responsible for creating, at the beginning of the Run, <BR>
 *  all the splitting function objects and the corresponding Sudakov form <BR>
 *  factors objects, and then of the generation of splittings (radiation <BR>
 *  emissions) during the event. <BR>
 *  Many switches are defined in this class which allowed the user to turn on/off:
 *  <UL>
 *   <LI> each type of interaction (QCD, QED, EWK,...);
 *   <LI> initial and final state radiation for all type of interactions;
 *   <LI> initial and final state radiation for each type of interaction;
 *   <LI> each type of splitting (U->UG, D->DG, ... , G->GG, G->UUbar,...);
 *  </UL>
 *  These switches are useful mainly for debugging, but eventually can <BR>
 *  also be used for a "quick and dirty" estimation of systematic errors. <BR>
 *
 *  ***LOOKHERE*** --- the idea is to keep this class not responsible
 *                     for creating new ShowerParticle objects, and
 *                     independent from ShowerVariables: therefore
 *                     the checking that the chosen candidate branching
 *                     is acceptable according to the vetos in ShowerVariables
 *                     and then, if accepted, the creation of the ShowerParticle
 *                     created from the branching, should all be done in
 *                     in ForwardShowerEvolver and BackwardShowerEvolver.
 *                     The advantages in doing that is that SplittingGenerator
 *                     is kept simpler and easier to manage. <BR>
 *                 --- Using the ShowerParticle object provided in input,
 *                     is should be possible TO IMPLEMENT IN THIS CLASS  
 *                     SplittingGenerator THE (1->2 ONLY) AZIMUTHAL-CORRELATIONS
 *                     FOR SOFT EMISSIONS DUE TO QCD COHERENCE. <BR>
 *                 --- SIMILARLY, HAVING THE RHO-D MATRIX IN THE ShowerParticle
 *                     OBJECT, AND THE SplitFun POINTER IN THE SUDAKOV OBJECT,
 *                     IT SHOULD BE POSSIBLE TO IMPLEMENT THE SPIN-CORRELATION.
 *               
 *  @see ShowerIndex
 *  @see SudakovFormFactor
 *  @see ShowerVariables
 *  @see SplitFun
 */
class SplittingGenerator: public ThePEG::HandlerBase {

public:

  /**
   * Standard ctors and dtor.
   */
  inline SplittingGenerator();
  inline SplittingGenerator(const SplittingGenerator &);
  virtual ~SplittingGenerator();

  /**
   * It chooses a new forward branching for the time-like particle. 
   * The method returns a pair of pointers: <BR>
   * --- a pointer to a ShowerKinematics object, which 
   *     contains the information about the new scale and all other
   *     kinematics variables that need to be generated simultaneously; <BR>
   * --- a pointer to the SudakovFormFactor object associated 
   *     with the chosen emission. <BR>
   * In the case no branching has been generated, both the returned 
   * pointers are null ( ShoKinPtr() , tSudakovFFPtr() ).
   * In the case that reverseAngularOrder is true, 
   * the new scale is greater than the initial one: this is used for 
   * the forward evolution of a on-shell decaying particle.
   * If something goes wrong (but this should never ever happen) 
   * a warning should be sent to the log file: in this case, 
   * exceptions seem unnecessary.
   */
  Branching chooseForwardBranching(tEHPtr, ShowerParticle &,
                                   const bool reverseAngOrd = false) const; 

  /**
   * Similar to the previous method, but for the backward evolution of
   * a space-like input particle. 
   * Notice that the PartialCollisionHandler object, ch, 
   * is necessary to access the PDFs.
   */
  Branching chooseBackwardBranching(tEHPtr, ShowerParticle&) const;

  /**
   * Given the particle, a pointer to the ShowerKinematics
   * object which has been created and at least partially filled during
   * the branching selection, and the pointer to the SudakovFormFactor
   * object associated with such selected branching, the method completes
   * the kinematics of the branching, storing the results in the same
   * ShowerKinematics object specified in input.
   * Notice that, it can be that this method does nothing, because
   * the kinematics of the branching has been already completely
   * determined during chooseForwardBranching or chooseBackwardBranching.
   * As for the latter two methods, if something goes wrong (but
   * this should never ever happen) a warning should be sent to
   * the log file: in this case exceptions seem unnecessary.
   */
  void generateBranchingKinematics(tEHPtr, ShowerParticle&, 
				   Branching&) const;

  /**
   * Access to the Initial and Final State Radiation QCD running alphas.
   */
  inline const tShowerAlphaPtr showerAlphaQCD() const;

  /**
   * Access to the QED coupling (a bit too complicated for my taste...).
   */
  inline const tShowerAlphaPtr showerAlphaQED() const;

  /**
   * Access to the ShowerVariables (maybe is not needed).
   */
  inline const ShowerVarsPtr & showerVariables() const;

  //--- SWITCHES ---

  void setSVtoAlpha(ShowerVarsPtr p);
  void setQED(ShowerAlphaPtr p);
  void setQCD(ShowerAlphaPtr p);

  /**
   * It returns true/false if interaction type specified in input is on/off.
   */
  bool isInteractionON(const ShowerIndex::InteractionType interaction) const;

  /**
   * It returns true/false if the initial or final state radiation is on/off.
   */
  inline bool isISRadiationON() const;  
  inline bool isFSRadiationON() const;  

  /**
   * It returns true/false if the initial or final state radiation for the
   * specified interaction type is on/off. However, they return false, 
   * regardless of the switch, if either the corresponding interaction switch 
   * (see method isInteractionON) is off, or if the global initial or final 
   * state radiation (see overloaded methods above without argument) is off.
   */
  bool isISRadiationON(const ShowerIndex::InteractionType interaction) const;  
  bool isFSRadiationON(const ShowerIndex::InteractionType interaction) const;

  /**
   * This method returns the splitting function associated with the two
   * ids. The first is the particle splitting, the second is the first
   * 'product' of that splitting (see initial vs. final splitting syntax)
   * The final bool value indicates whether this is a search over the final
   * or initial state branchings.
   */
  tSplittingFnPtr getSplittingFunction(long id1, long id2, bool initial=true);

public:

  /**
   * Standard functions for writing and reading from persistent streams.
   */
  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);

  /**
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();

protected:

  /**
   * Standard clone methods.
   */
  inline virtual IBPtr clone() const;
  inline virtual IBPtr fullclone() const;

protected:

  /**
   * Standard Interfaced virtual functions.
   */
  string addFinalSplitting(string);
  string addInitialSplitting(string);
  string addSplitting(string,bool);
  inline virtual void doupdate() throw(UpdateException);
  inline virtual void doinit() throw(InitException);
  inline virtual void doinitrun();
  inline virtual void dofinish();

  /**
   * Change all pointers to Interfaced objects to corresponding clones.
   */
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);

  /**
   * Return pointers to all Interfaced objects refered to by this.
   */
  inline virtual IVector getReferences();

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<SplittingGenerator> initSplittingGenerator;

  /**
   * Private and non-existent assignment operator.
   */
  SplittingGenerator & operator=(const SplittingGenerator &);

  //void initializeRun();
  // It is called once, at the beginning of the run, to create all
  // of the splitting function and Sudakov form factor objects. 
  // The splitting function objects are then kept by the corresponding
  // Sudakov form factor objects. The latter are kept in a multimap
  // with key given by a ShowerIndex object.

  /**
   * Print, in the log file, debugging information.
   */
  void addToMap(IdList &, SudakovPtr &, bool);
  void debuggingInfo();

  int _QCDinteractionMode;
  int _QEDinteractionMode;
  int _EWKinteractionMode;
  int _ISR_Mode;
  int _ISR_QCDMode;
  int _ISR_QEDMode;
  int _ISR_EWKMode;
  int _FSR_Mode;
  int _FSR_QCDMode;
  int _FSR_QEDMode;
  int _FSR_EWKMode;

  ShowerAlphaPtr _showerAlphaQCD;
  ShowerAlphaPtr _showerAlphaQED;
  ShowerVarsPtr _showerVariables;
   
  typedef pair<SudakovPtr,IdList> BranchingElement;
  typedef multimap<long,BranchingElement> BranchingList;
  typedef pair<long, BranchingElement> BranchingInsert; 

  /**  
   * Lists of the branchings and the appropriate Sudakov.
   */
  BranchingList _fbranchings;
  BranchingList _bbranchings;

};

}


namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of SplittingGenerator.
 */
template <>
struct BaseClassTrait<Herwig::SplittingGenerator,1> {
  typedef ThePEG::HandlerBase NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::SplittingGenerator>: public ClassTraitsBase<Herwig::SplittingGenerator> {

  /**  
   * Return the class name.
   */
  static string className() { return "/Herwig++/SplittingGenerator"; }

  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwShower.so"; }
};

}

#include "SplittingGenerator.icc"

#endif /* HERWIG_SplittingGenerator_H */
