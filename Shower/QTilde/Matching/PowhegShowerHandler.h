// -*- C++ -*-
//
// PowhegShowerHandler.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_PowhegShowerHandler_H
#define HERWIG_PowhegShowerHandler_H
//
// This is the declaration of the PowhegShowerHandler class.
//

#include "Herwig/Shower/QTilde/QTildeShowerHandler.h"
#include "Herwig/MatrixElement/Matchbox/MatchboxFactory.h"
#include "Herwig/MatrixElement/HwMEBase.h"
#include "Herwig/Shower/Core/Base/HardBranching.h"

#include "Herwig/Shower/QTilde/Matching/CKKWTree.h"
#include "Herwig/Shower/QTilde/Matching/ProtoTree.h"
#include "Herwig/Shower/QTilde/Matching/ProtoBranching.h"
#include "Herwig/Shower/QTilde/Matching/PotentialTree.h"

#include "ThePEG/MatrixElement/DiagramBase.fh"
#include "ThePEG/MatrixElement/MEBase.h"
#include "ThePEG/PDF/PartonExtractor.h"
#include "Herwig/Shower/Core/Base/HardTree.h"

// #include "Herwig/Shower/Core/SplittingFunctions/SplittingGenerator.h"
// #include "Herwig/Shower/QTilde/Base/ShowerModel.h"
// #include "ThePEG/PDF/BeamParticleData.h"
// #include "Herwig/Shower/Core/Base/ShowerTree.h"
// #include "Herwig/Shower/Core/Base/ShowerProgenitor.fh"
// #include "Herwig/Shower/QTilde/QTildeShowerHandler.fh"
// #include "Herwig/Shower/Core/Base/Branching.h"
// #include "Herwig/Shower/QTilde/Base/ShowerVeto.h"
// #include "ThePEG/Handlers/XComb.h"
// #include "Herwig/Decay/HwDecayerBase.h"


namespace Herwig {

using namespace ThePEG;

class PowhegShowerHandler: public QTildeShowerHandler {

public:

  /**
   * The default constructor.
   */
  PowhegShowerHandler() : subtractionIntegral_(false),
			  enforceColourConsistency_(false),
			  forcePartners_(false),
			  decayRadiation_(0)
  {}

public:

  Ptr<MatchboxFactory>::ptr Factory(){return theFactory;}

  Ptr<MatchboxFactory>::ptr Factory() const {return theFactory;}


  /**
   * Return true, if the shower handler can generate a truncated 
   * shower for POWHEG style events generated using Matchbox
   */
  virtual bool canHandleMatchboxTrunc() const { return true; }

protected:

  /**
   *  Generate hard emissions for CKKW etc
   */
  virtual HardTreePtr generateCKKW(ShowerTreePtr tree) const;

protected:

  /**
   *  Access to the core matrix element
   */
  MEPtr matrixElement() const {return matrixElement_;}

 /**
   * Creates all (ordered) cluster histories and selects one.
   */
  PotentialTree doClustering(tSubProPtr sub,ShowerTreePtr showerTree) const;

  /**
   *  Access to the select tree
   */
  PotentialTree & hardTree() {return hardTree_;}
  const PotentialTree & hardTree() const {return hardTree_;}

  /**
   *  Access to the select tree
   */
  void hardTree(PotentialTree in) const {hardTree_ = in;}

  /**
   * Check if two momenta are equal within 1%
   */
  bool fuzzyEqual(Lorentz5Momentum a, Lorentz5Momentum  b) const;

  /**
   *  Access to the potential branchings
   */
  set<ProtoBranchingPtr> & protoBranchings() const {return protoBranchings_;}

  /**
   * The ProtoTrees which will become CKKWTrees
   */
  set< ProtoTreePtr > & protoTrees() const {return protoTrees_;}

  /**
   * Recursive function to find all possible clustered trees.
   * Does not produce any repeated trees.
   */
  void fillProtoTrees( ProtoTreePtr , long id ) const;

 /**
   * Function looks to see if a cluster of the branchings already exists
   * in protoBranchings_ if so returns the pointer to that protoBranching
   * if not creates the hardBranchings, adds it
   * to protoBranchings and returns the pointer
   */
  tProtoBranchingPtr getCluster( tProtoBranchingPtr, tProtoBranchingPtr ) const;

  /**
   * Checks whether a ProtoTree containing the same branchings already 
   * exists in protoTrees_ in which case the current tree is a repeat 
   * and should be removed (and not recursed)
   */
  bool repeatProtoTree( ProtoTreePtr currentProtoTree ) const;

  /**
   * Returns the branching element for an FS-FS clustering
   */
  BranchingElement allowedFinalStateBranching( tProtoBranchingPtr &,
					       tProtoBranchingPtr &) const;
  
  /**
   * Returns the branching element for an IS-FS clustering
   */
  BranchingElement allowedInitialStateBranching( tProtoBranchingPtr & ,
						 tProtoBranchingPtr &) const;

  /**
   * Returns the diagram corresponding to the (leading-order) hardTree
   */
  bool checkDiagram(PotentialTree &,tcDiagPtr) const;

  bool subtractionIntegral() const {return subtractionIntegral_;}

  void setSubtractionIntegral(bool subInt) const { subtractionIntegral_=subInt;}


private:

   /**
   * The factory object to fetch splitting channels from
   */
  Ptr<MatchboxFactory>::ptr theFactory;

 /**
   *  The matrix element for the core process
   */
  mutable MEPtr matrixElement_;

  /**
   *  The chosen hard tree
   */
  mutable PotentialTree hardTree_;

 /**
   *  All the Protobranchings used in finding the shower histories
   */
  mutable set<ProtoBranchingPtr> protoBranchings_;

  /**
   * The ProtoTrees which will become CKKWTrees
   */
  mutable set< ProtoTreePtr > protoTrees_;

 /**
   * The possible shower configurations that are angular-ordered
   */
  mutable vector< pair< PotentialTree, double > > hardTrees_;

  /**
   *  Which branchings are allowed?
   */
  //@{
  /**
   *  The allowed final-state branchings
   */
  mutable map<pair<long,long>,BranchingElement > allowedFinal_;

  /**
   *  The allowed initial-state branchings
   */
  mutable multimap<long, BranchingElement > allowedInitial_;
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

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  PowhegShowerHandler & operator=(const PowhegShowerHandler &);

private:

  /**
   *  Emitter particle from the original generation
   */
  mutable int emitter_;

  /**
   *  Spectator particle from the original generation
   */
  mutable int spectator_;

  /**
   *  Whether or not a subtraction integral
   */
  mutable bool subtractionIntegral_;

  /**
   * Whether or not do enforce consistency of the Born and real colour flows
   */
  bool enforceColourConsistency_;

  /**
   *  Force emitter and spectator partners
   */
  bool forcePartners_;

  /**
   *  Handling of radiation in decays
   */
  unsigned int decayRadiation_;

};

}

#endif /* HERWIG_PowhegShowerHandler_H */
