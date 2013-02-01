// -*- C++ -*-
#ifndef Herwig_MatchingHandler_H
#define Herwig_MatchingHandler_H
//
// This is the declaration of the MatchingHandler class.
//

#include "Herwig++/Shower/ShowerHandler.h"
#include "ThePEG/MatrixElement/MEBase.h"
#include "Herwig++/MatrixElement/HwMEBase.h"
#include "ThePEG/MatrixElement/DiagramBase.fh"
#include "ThePEG/PDF/PartonExtractor.h"
#include "ProtoTree.h"
#include "PotentialTree.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the MatchingHandler class.
 *
 * @see \ref MatchingHandlerInterfaces "The interfaces"
 * defined for MatchingHandler.
 */
class MatchingHandler: public ShowerHandler {

public:

  /**
   * The default constructor.
   */
  MatchingHandler(bool reWeight = false);

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
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}

protected:

  /**
   *  Access to the core matrix element
   */
  MEPtr matrixElement() const {return matrixElement_;}

  /**
   *  Access to the core Hw++ matrix element
   */
  HwMEBasePtr HWMatrixElement() const {return HWmatrixElement_;}
  
  /**
   * The PartonExtractor object used to construct remnants.
   */
  PExtrPtr partonExtractor() { return partonExtractor_;}

  /**
   * The pairs of PartonBin objects describing the partons which can
   * be extracted by the PartonExtractor object.
   */
  PartonPairVec partonBins() {return partonBins_;}

  /**
   * The pairs of PartonBin objects describing the partons which can
   * be extracted by the PartonExtractor object.
   */
  void partonBins(PartonPairVec bins) {partonBins_ = bins;}

  /**
   * The Cuts object to be used for this reader.
   */
  CutsPtr & cuts() {return cuts_;}

  /**
   *  Initialise the matching
   */
  void initialiseMatching(int minMult, int maxMult);

protected:

  /**
   *  Pointer to the object calculating the strong coupling
   */
  ShowerAlphaPtr alphaS() {return alphaS_;}

  /**
   * whether to just reject event if there are no angular ordered histories (true)
   * or to choose one of the non angular ordered histories (false)
   */
  bool rejectNonAngularOrdered() {return rejectNonAO_;}

  /**
   * whether to reject events (true) for which there is no parton shower interpretation
   * or just shower them with no reweighting
   */
  bool rejectNoShowerHistory() {return rejectNOHistories_;}

  /**
   *  Access to basic info on the event
   */
  //@{
  /**
   *  Centre of mass energy
   */
  Energy2 sHat() const {return sHat_;}

  /**
   *  Centre of mass energy
   */
  void sHat(Energy2 in) {sHat_ = in;}
  
  /**
   * The fixed factorization scale used in the MEs.
   */
  Energy pdfScale() const {return pdfScale_;}

  /**
   * The fixed factorization scale used in the MEs.
   */
  void pdfScale(Energy in) {pdfScale_ = in;}

  /**
   * The fixed alphaS value that was used to generate parton configurations
   */
  double alphaSMG() const {return alphaSMG_;}

  /**
   * The fixed alphaS value that was used to generate parton configurations
   */
  void alphaSMG(double in) {alphaSMG_ = in;}

  /**
   *  Whether we are treating an event with no shower interpretation
   */
  bool noShowerHistory() {return noShowerHists_;}

  /**
   *  whether the current event is a highest multiplicity event.
   */
  bool highestMult() const {return highestMult_;}

  /**
   *  whether the current event is a lowest multiplicity event
   */
  bool lowestMult() const {return lowestMult_;}
  //@}

  /**
   *  Storage of potential shower interpretations
   */
  //@{
  /**
   *  Access to the potential branchings
   */
  set<ProtoBranchingPtr> & protoBranchings() {return protoBranchings_;}

  /**
   * The ProtoTrees which will become CKKWTrees
   */
  set< ProtoTreePtr > & protoTrees() {return protoTrees_;}

  /**
   *  Access to the select tree
   */
  PotentialTree & hardTree() {return hardTree_;}

  /**
   *  Access to the select tree
   */
  void hardTree(PotentialTree in) {hardTree_ = in;}

  /**
   *  access to the hard tree object
   */
  CKKWTreePtr getCKKWTree() const {return hardTree_.tree();}

  /**
   * The possible shower configurations that are angular-ordered
   */
  vector< pair< PotentialTree, double > > & hardTrees() {
    return hardTrees_;
  }

  /**
   * The possible shower configurations that are not angular-ordered
   */
  vector< pair< PotentialTree, double > > & nonOrderedTrees() {
    return nonOrderedTrees_;
  }
  //@}

  /**
   * Returns the diagram corresponding to the (leading-order) hardTree
   */
  double getDiagram(PotentialTree &);

  /**
   * Recursive function to find all possible clustered trees.
   * Does not produce any repeated trees.
   */
  void fillProtoTrees( ProtoTreePtr );

  /**
   *  update sub process based on the CKKWTree and Diagram
   */
  bool updateSubProcess();

  /**
   * Creates all (ordered) cluster histories and selects one.
   */
  PotentialTree doClustering();
 
  /**
   *  Calculate the Sudakov weight
   */
  virtual double sudakovWeight( CKKWTreePtr ) = 0;

  /**
   *  Select the hard tree
   */
  virtual PotentialTree chooseHardTree(double totalWeight,
				       double nonOrderedWeight) = 0;

  /**
   * Function looks to see if a cluster of the branchings already exists
   * in protoBranchings_ if so returns the pointer to that protoBranching
   * if not creates the hardBranchings, adds it
   * to protoBranchings and returns the pointer
   */
  tProtoBranchingPtr getCluster( tProtoBranchingPtr, tProtoBranchingPtr );

  /**
   * Checks whether a ProtoTree containing the same branchings already 
   * exists in protoTrees_ in which case the current tree is a repeat 
   * and should be removed (and not recursed)
   */
  bool repeatProtoTree( ProtoTreePtr currentProtoTree );

  /**
   * Returns the branching element for an FS-FS clustering
   */
  BranchingElement allowedFinalStateBranching( tProtoBranchingPtr &,
					       tProtoBranchingPtr &);
  
  /**
   * Returns the branching element for an IS-FS clustering
   */
  BranchingElement allowedInitialStateBranching( tProtoBranchingPtr & ,
						 tProtoBranchingPtr &);

  /**
   *  Find decaying particles
   */
  void findDecayingParticles(PPtr parent);

  /**
   *  Check if a decayed particle
   */
  map<PPtr,ParticleVector>::const_iterator parent(PPtr parent); 

  /**
   *  Add the decay products of a particle to the new SubProcess
   */
  void addDecayProducts(SubProPtr subProcess, PPtr parent,
			map<PPtr,ParticleVector>::const_iterator decay,
			const LorentzRotation & boost) const;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MatchingHandler & operator=(const MatchingHandler &);

private:

  /**
   *  Whether or not to reweight
   */
  bool reWeight_;

  /**
   * whether to just reject event if there are no angular ordered histories (true)
   * or to choose one of the non angular ordered histories (false)
   */
  bool rejectNonAO_;

  /**
   * whether to reject events (true) for which there is no parton shower interpretation
   * or just shower them with no reweighting
   */
  bool rejectNOHistories_;

  /**
   *  Include decaying particles
   */
  bool includeDecays_;

  /**
   *  The matrix element for the core process
   */
  MEPtr matrixElement_;

  /**
   *  The cast matrix element
   */
  HwMEBasePtr HWmatrixElement_;
  
  /**
   * The PartonExtractor object used to construct remnants.
   */
  PExtrPtr partonExtractor_;

  /**
   * The pairs of PartonBin objects describing the partons which can
   * be extracted by the PartonExtractor object.
   */
  PartonPairVec partonBins_;

  /**
   * The Cuts object to be used for this reader.
   */
  CutsPtr cuts_;

  /**
   *  Centre of mass energy
   */
  Energy2 sHat_;

  /**
   *  Pointer to the object calculating the strong coupling
   */
  ShowerAlphaPtr alphaS_;
  
  /**
   * The fixed factorization scale used in the MEs.
   */
  Energy pdfScale_;

  /**
   * The fixed alphaS value that was used to generate parton configurations
   */
  double alphaSMG_;

private:

  /**
   *  Storage of the potential shower interpretation
   */
  //@{
  /**
   *  All the Protobranchings used in finding the shower histories
   */
  set<ProtoBranchingPtr> protoBranchings_;

  /**
   * The ProtoTrees which will become CKKWTrees
   */
  set< ProtoTreePtr > protoTrees_;

  /**
   *  The chosen hard tree
   */
  PotentialTree hardTree_;

  /**
   * The possible shower configurations that are angular-ordered
   */
  vector< pair< PotentialTree, double > > hardTrees_;

  /**
   * The possible shower configurations that are not angular-ordered
   */
  vector< pair< PotentialTree, double > > nonOrderedTrees_;

  /**
   *  Any decays
   */
  map<PPtr,ParticleVector> decayingParticles_;

  /**
   *  Whether we are treating an event with no shower interpretation
   */
  bool noShowerHists_;

  /**
   *  Whether we are treating the highest multiplicity contribution
   */
  bool highestMult_;
  
  /**
   *  Whether we are treating the lowest multiplicity contribution
   */
  bool lowestMult_;
  //@}

  /**
   *  Which branchings are allowed?
   */
  //@{
  /**
   *  The allowed final-state branchings
   */
  map<pair<long,long>,pair<SudakovPtr,IdList> > allowedFinal_;

  /**
   *  The allowed initial-state branchings
   */
  multimap<long, pair<SudakovPtr,IdList> > allowedInitial_;
  //@}
};

}

#endif /* Herwig_MatchingHandler_H */
