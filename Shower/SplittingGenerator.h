// -*- C++ -*-
#ifndef HERWIG_SplittingGenerator_H
#define HERWIG_SplittingGenerator_H
//
// This is the declaration of the <!id>SplittingGenerator<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This class is responsible for creating, at the beginning of the Run, <BR>
// all the splitting function objects and the corresponding Sudakov form <BR>
// factors objects, and then of the generation of splittings (radiation <BR>
// emissions) during the event. <BR>
// Many switches are defined in this class which allowed the user to turn on/off:
// <UL>
//   <LI> each type of interaction (QCD, QED, EWK,...);
//   <LI> initial and final state radiation for all type of interactions;
//   <LI> initial and final state radiation for each type of interaction;
//   <LI> each type of splitting (U->UG, D->DG, ... , G->GG, G->UUbar,...);
// </UL>
// These switches are useful mainly for debugging, but eventually can <BR>
// also be used for a "quick and dirty" estimation of systematic errors. <BR>
//
// ***LOOKHERE*** --- the idea is to keep this class not responsible
//                    for creating new ShowerParticle objects, and
//                    independent from ShowerVariables: therefore
//                    the checking that the chosen candidate branching
//                    is acceptable according to the vetos in ShowerVariables
//                    and then, if accepted, the creation of the ShowerParticle
//                    created from the branching, should all be done in
//                    in ForwardShowerEvolver and BackwardShowerEvolver.
//                    The advantages in doing that is that SplittingGenerator
//                    is kept simpler and easier to manage. <BR>
//                --- Using the ShowerParticle object provided in input,
//                    is should be possible TO IMPLEMENT IN THIS CLASS  
//                    SplittingGenerator THE (1->2 ONLY) AZIMUTHAL-CORRELATIONS
//                    FOR SOFT EMISSIONS DUE TO QCD COHERENCE. <BR>
//                --- SIMILARLY, HAVING THE RHO-D MATRIX IN THE ShowerParticle
//                    OBJECT, AND THE SplitFun POINTER IN THE SUDAKOV OBJECT,
//                    IT SHOULD BE POSSIBLE TO IMPLEMENT THE SPIN-CORRELATION.
//               
// CLASSDOC SUBSECTION See also:
//
// <a href="http:ShowerIndex.html">ShowerIndex.h</a>, <BR>
// <a href="http:SudakovFormFactor.html">SudakovFormFactor.h</a>, <BR>
// <a href="http:ShowerVariables.html">ShowerVariables.h</a>, <BR>
// <a href="http:SplitFun.html">SplitFun.h</a>.
// 

#include "ThePEG/Handlers/HandlerBase.h"
#include "ShowerConfig.h"
#include "Herwig++/Utilities/GlobalParameters.h"
#include "ThePEG/Handlers/PartialCollisionHandler.h"
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

class SplittingGenerator: public ThePEG::HandlerBase {

public:

  inline SplittingGenerator();
  inline SplittingGenerator(const SplittingGenerator &);
  virtual ~SplittingGenerator();
  // Standard ctors and dtor.

  Branching chooseForwardBranching(tPartCollHdlPtr, ShowerParticle &,
                                   const bool reverseAngOrd = false) const; 
  // It chooses a new forward branching for the time-like particle. 
  // The method returns a pair of pointers: <BR>
  // --- a pointer to a <!id>ShowerKinematics<!!id> object, which 
  //     contains the information about the new scale and all other
  //     kinematics variables that need to be generated simultaneously; <BR>
  // --- a pointer to the <!id>SudakovFormFactor<!!id> object associated 
  //     with the chosen emission. <BR>
  // In the case no branching has been generated, both the returned 
  // pointers are null ( ShoKinPtr() , tSudakovFFPtr() ).
  // In the case that <!id>reverseAngularOrder<!!id> is true, 
  // the new scale is greater than the initial one: this is used for 
  // the forward evolution of a on-shell decaying particle.
  // If something goes wrong (but this should never ever happen) 
  // a warning should be sent to the log file: in this case, 
  // exceptions seem unnecessary.

  Branching chooseBackwardBranching(tPartCollHdlPtr, ShowerParticle&) const;
  // Similar to the previous method, but for the backward evolution of
  // a space-like input particle. 
  // Notice that the PartialCollisionHandler object, <!id>ch<!!id>, 
  // is necessary to access the PDFs.

  void generateBranchingKinematics(tPartCollHdlPtr, ShowerParticle&, 
				   Branching&) const;
  // Given the particle, a pointer to the <!id>ShowerKinematics<!!id> 
  // object which has been created and at least partially filled during
  // the branching selection, and the pointer to the <!SudakovFormFactor<!!id>
  // object associated with such selected branching, the method completes
  // the kinematics of the branching, storing the results in the same
  // <!id>ShowerKinematics<!!id> object specified in input.
  // Notice that, it can be that this method does nothing, because
  // the kinematics of the branching has been already completely
  // determined during <!id>chooseForwardBranching<!!id> or
  // <!id>chooseBackwardBranching<!!id>.
  // As for the latter two methods, if something goes wrong (but
  // this should never ever happen) a warning should be sent to
  // the log file: in this case exceptions seem unnecessary.

  inline const tShowerAlphaPtr showerAlphaQCD() const;
  // Access to the Initial and Final State Radiation QCD running alphas.

  inline const tShowerAlphaPtr showerAlphaQED() const;
  // Access to the QED coupling (a bit too complicated for my taste...)

  inline const ShowerVarsPtr & showerVariables() const;
  // Access to the ShowerVariables (maybe is not needed).

  //--- SWITCHES ---

  void setSVtoAlpha(ShowerVarsPtr p);
  void setQED(ShowerAlphaPtr p);
  void setQCD(ShowerAlphaPtr p);

  bool isInteractionON(const ShowerIndex::InteractionType interaction) const;
  // It returns true/false if interaction type specified in input is on/off.

  inline bool isISRadiationON() const;  
  inline bool isFSRadiationON() const;  
  // It returns true/false if the initial or final state radiation is on/off.

  bool isISRadiationON(const ShowerIndex::InteractionType interaction) const;  
  bool isFSRadiationON(const ShowerIndex::InteractionType interaction) const;
  // It returns true/false if the initial or final state radiation for the
  // specified interaction type is on/off. However, they return false, 
  // ragardless of the switch, if either the corresponding interaction switch 
  // (see method <!id>isInteractionON<!!id>) is off, or if the global initial or final 
  // state radiation (see overloaded methods above without argument) is off.

  tSplittingFnPtr getSplittingFunction(long id1, long id2, bool initial=true);
  // This method returns the splitting function associated with the two
  // ids. The first is the particle splitting, the second is the first
  // 'product' of that splitting (see initial vs. final splitting syntax)
  // The final bool value indicates whether this is a search over the final
  // or initial state branchings
public:

  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);
  // Standard functions for writing and reading from persistent streams.

  static void Init();
  // Standard Init function used to initialize the interfaces.

protected:

  inline virtual IBPtr clone() const;
  inline virtual IBPtr fullclone() const;
  // Standard clone methods.

protected:
  string addFinalSplitting(string);
  string addInitialSplitting(string);
  string addSplitting(string,bool);
  inline virtual void doupdate() throw(UpdateException);
  inline virtual void doinit() throw(InitException);
  inline virtual void doinitrun();
  inline virtual void dofinish();
  // Standard Interfaced virtual functions.

  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);
  // Change all pointers to Interfaced objects to corresponding clones.

  inline virtual IVector getReferences();
  // Return pointers to all Interfaced objects refered to by this.

private:

  static ClassDescription<SplittingGenerator> initSplittingGenerator;
  // Describe a concrete class with persistent data.

  SplittingGenerator & operator=(const SplittingGenerator &);
  // Private and non-existent assignment operator.

  //void initializeRun();
  // It is called once, at the beginning of the run, to create all
  // of the splitting function and Sudakov form factor objects. 
  // The splitting function objects are then kept by the corresponding
  // Sudakov form factor objects. The latter are kept in a multimap
  // with key given by a ShowerIndex object.

  void addToMap(IdList &, SudakovPtr &, bool);
  void debuggingInfo();
  // Print, in the log file, debugging information.

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
  // Lists of the branchings and the appropriate sudakov
  BranchingList _fbranchings;
  BranchingList _bbranchings;

};

}

// CLASSDOC OFF

namespace ThePEG {

// The following template specialization informs ThePEG about the
// base class of SplittingGenerator.
template <>
struct BaseClassTrait<Herwig::SplittingGenerator,1> {
  typedef ThePEG::HandlerBase NthBase;
};

// The following template specialization informs ThePEG about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::SplittingGenerator>: public ClassTraitsBase<Herwig::SplittingGenerator> {
  static string className() { return "/Herwig++/SplittingGenerator"; }
  // Return the class name.
  static string library() { return "libHwShower.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#include "SplittingGenerator.icc"

#endif /* HERWIG_SplittingGenerator_H */
