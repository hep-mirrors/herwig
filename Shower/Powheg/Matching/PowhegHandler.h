// -*- C++ -*-
#ifndef HERWIG_PowhegHandler_H
#define HERWIG_PowhegHandler_H
//
// This is the declaration of the PowhegHandler class.
//

#include "Herwig++/Utilities/Histogram.h"
#include "Herwig++/Shower/ShowerHandler.h"
#include "Herwig++/Shower/Base/Evolver.h"
#include "Herwig++/Utilities/Interpolator.h"
#include "Herwig++/Utilities/Interpolator2d.h"
#include "Herwig++/Utilities/GaussianIntegrator.h"
#include "Herwig++/Shower/Base/ShowerParticle.h"
#include "ThePEG/Config/Pointers.h"
#include "Herwig++/Shower/Base/HardTree.h"
#include "Herwig++/Shower/Couplings/ShowerAlpha.h"
#include "ThePEG/MatrixElement/MEBase.h"
#include "ThePEG/MatrixElement/DiagramBase.fh"
#include "ThePEG/PDF/PartonExtractor.h"
#include "Herwig++/Shower/Base/CKKWVeto.h"


namespace Herwig {
  class PowhegHandler;
  class ProtoTree;
}

//declaration of thepeg ptr
namespace ThePEG {

  ThePEG_DECLARE_POINTERS(Herwig::PowhegHandler,PowhegHandlerPtr);
  ThePEG_DECLARE_POINTERS(Herwig::ProtoTree,ProtoTreePtr);

}

namespace Herwig {

using namespace ThePEG;

  class ProtoTree:public Base {
  public:
    ProtoTree() {}
    ProtoTree(const set< HardBranchingPtr > & newBranchings) : theBranchings( newBranchings ) {}
    
    void addBranching( HardBranchingPtr Branching ){
      theBranchings.insert( Branching );
    }
    bool removeBranching( HardBranchingPtr Branching ){
      if( theBranchings.find( Branching ) != theBranchings.end() ){
	theBranchings.erase( Branching );
	return true;
      }
      else return false;
    }    
    /**
     * Fix the parent-child relations along backward lines in 
     * the hardBranchings 
     */
    bool fixFwdBranchings();

    const set< HardBranchingPtr > & getBranchings() const {
      return  theBranchings;
    }
  private:
    set< HardBranchingPtr > theBranchings;
  };

/**
 * Here is the documentation of the PowhegHandler class.
 *
 * @see \ref PowhegHandlerInterfaces "The interfaces"
 * defined for PowhegHandler.
 */
class PowhegHandler: public ShowerHandler {

public:

  /**
   * The default constructor.
   */
  PowhegHandler() : _npoint(10), _sudopt(0), _sudname("sudakov.data"), _jetMeasureMode(1), _lepton(true), _reweightOpt(0), 
		    _highestMult(false), _lowestMult(false), _testSudakovs(false), _qtildeDist( false ),
		    _yini(0.001), _alphaSMG(0.118), _max_qtilde( 91.2*GeV ), _max_pt_cut( 45.6*GeV ), _min_pt_cut( 0.*GeV ), 
		    _clusterOption( 0 ),  _rejectNonAO( true ), _rejectNoHistories( true ), _pdfScale( 91.2*GeV ),
		    _dynamicSuds( false ), _cutOption( 0 ) {}

  /**
   * Perform CKKW reweighting
   */
  virtual double reweightCKKW(int minMult, int maxMult);

  /**
   *  Main method for the cascade
   */
  virtual void cascade();

public:

  /**
   *  access to the hard tree object
   */
  inline HardTreePtr getHardTree() {
    return _theHardTree;
  }

  /**
   *  access to the merging scale
   */
  inline Energy getMergeScale() {
    return sqrt( _yini * _s );
  }
  /**
   *  whether the current event is a highest multiplicity event.
   */
  inline bool highestMult(){
    return _highestMult;
  }
  /**
   *  whether the current event is a lowest multiplicity event
   */
  inline bool lowestMult(){
    return _lowestMult;
  }
  /**
   *  access to the vetoes
   */
  inline Ptr<CKKWVeto>::pointer theVeto(){
    return _showerVeto;
  }
  /**
   *  access to the evolver
   */
  inline tEvolverPtr getEvolver(){
    return evolver();
  }
  
  /**
   *  access to the jet measure definition being used
   */
  inline unsigned int jetMeasureMode() {
    return _jetMeasureMode;
  }

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

  virtual void dofinish();

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();
  //@}
  
  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

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

  /**
   * All clusters used in finding all shower histories
   */
  map< HardBranchingPtr , pair< ShowerParticlePtr, ShowerParticlePtr > > _all_clusters;

  /**
   * The ProtoTrees which will become HardTrees
   */
  set< ProtoTreePtr > _proto_trees;

  /**
   * The possible shower configurations that are angular-ordered
   */
  vector< pair< HardTreePtr, double > > _hardTrees;
  
  /**
   * The possible shower configurations that are not angular-ordered
   */
  vector< pair< HardTreePtr, double > > _rejectedHardTrees;

  /**
   * Set colour connection in HardTree outwards from the hard process.
   * Currently assumes qq_bar leading order parton configuration (IS or FS) 
   */
  bool simpleColConnections( HardTreePtr currentHardTree );
  
  /**
   * Recursive function to find all possible clustered trees.
   * Does not produce any repeated trees.
   */
  bool fillProtoTrees( map< ShowerParticlePtr, HardBranchingPtr >, 
		       ProtoTreePtr );
  /**
   * Function looks to see if a cluster of the cluster_particles already exists in _all_clusterings
   * if so returns the pointer to that hardBranching if not creates the hardBranchings, adds it
   * to _all_branchings and returns the pointer
   */
   HardBranchingPtr getCluster( pair< ShowerParticlePtr, ShowerParticlePtr >, 
				map< ShowerParticlePtr, HardBranchingPtr >, bool incoming );

  /**
   *checks whether a prototree containing the same branchings exists already in _proto_trees in 
   *which case the current tree is a repeat and should be removed (and not recursed)
   */
  bool repeatProtoTree( ProtoTreePtr currentProtoTree );

  /**
   * Clusters the partons and creates a branching history
   * by combining the 2 particles with smallest
   * jet measure as in traditional CKKW.
   */
  HardTreePtr doClustering( );

  /**
   * Creates all (ordered) cluster histories and selects one.
   */
  HardTreePtr doClusteringOrdered( );

  /**
   * Returns Sudakov factor
   */
  double Sud( Energy scale, long id, Energy pt_cut );

  /**
   * Returns the branching element for an FS-FS clustering
   */
  BranchingElement allowedFinalStateBranching( pair< ShowerParticlePtr, ShowerParticlePtr > & );
  
  /**
   * Returns the branching element for an IS-FS clustering
   */
  BranchingElement allowedInitialStateBranching( pair< ShowerParticlePtr, ShowerParticlePtr > & );

  /**
   * Returns the diagram corresponding to the (leading-order) hardTree
   */
  DiagPtr getDiagram(HardTreePtr);

  /**
   * Returns the hadronic jet measure for two clustered momenta
   */
  Energy2 hadronJetMeasure(const Lorentz5Momentum & p1,
			   const Lorentz5Momentum & p2,
			   bool final=true);
  
  /**
   * General function for assigning a colour flow to a HardTree (not yet implemented
   */
  void createColourFlow(HardTreePtr,DiagPtr);
  
  /**
   * Sets the beam relations in the hard subprocess
   */
  void setBeams(HardTreePtr);
 
  /**
   *  Calculate the Sudakov weight
   */
  double sudakovWeight( HardTreePtr );
 
  /**
   *  Calculate the PDF weight
   */
  double pdfWeight( HardTreePtr );
 
  /**
   *  Calculate the splitting function weight
   */
  double splittingFnWeight( HardTreePtr );

  /**
   * Checks to see that the splitting is allowed.
   */
  bool splittingAllowed( ShowerParticlePtr part_i,
			 ShowerParticlePtr part_j );
  
  /**
   *  Sort's out the colour lines
   */
  void fixColours(tPPtr parent, tPPtr child1, tPPtr child2);
  
  /**
   * Checks to see that the splitting is allowed and finds the
   * Sudakov form factor for the splitting.
   */
  SudakovPtr getSud(long & emmitter_id,
		    ShowerParticlePtr & part_i, 
		    ShowerParticlePtr & part_j ) ;
  
  /**
   * Returns the durham jet measure, yij, for the two particles. 
   */
  double getJetMeasure(ShowerParticlePtr part_i, ShowerParticlePtr part_j);
  
  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<PowhegHandler> initPowhegHandler;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  PowhegHandler & operator=(const PowhegHandler &);

  /**                                                              
   * Outputs some sudakov test histograms 
   */
  void testSudakovs();

  /**
   * Outputs qtilde distributions from the tabulated sudakovs
   * for comparison with those from the shower (single emission)
   */
  void makeQtildeDist();

private:
  /**
   *  update sub process based on _theHardTree
   */
  bool updateSubProcess();

  /**
   *  The chosen hard tree
   */
  HardTreePtr _theHardTree;

  /**
   *  Number of points for the interpolation tables
   */
  unsigned int _npoint;

  /**
   *  Centre of mass energy
   */
  Energy2 _s;

  /**
   *  Map containing the sudakovs for the final-state particles
   */
  multimap< long, pair< Interpolator2d< double, Energy, Energy >::Ptr, Energy > > _fbranchings;

  /**
   *  Pointer to the object calculating the strong coupling
   */
  ShowerAlphaPtr _alphaS;
  
  /*
   * The global alphaS factor to ensure all alphaS weights are less than one
   */
  double _global_alphaS_wgt;
 
  /**
   *  Option for the Sudakov table
   */
  unsigned int _sudopt;

  /**
   *  Filename for the Sudakov table
   */
  string _sudname;
  
  /**
   * The jet measure algorithm we are using.
   */
  unsigned int _jetMeasureMode;

  /**
   *  Whether \f$e^+e^-\f$ or hadron-hadron
   */
  bool _lepton;

  /**
   *  The reweighting option being used
   */
  unsigned int _reweightOpt;

  /**
   *  Whether we are treating the highest multiplicity contribution
   */
  bool _highestMult;
  
  /**
   *  Whether we are treating the lowest multiplicity contribution
   */
  bool _lowestMult;

  /**
   *  Whether we are treating an event with no shower interpretation
   */
  bool _noShowerHists;

  /**                                                                                                   
   *  Whether the sudakovs should be tested                                                                    
   */
  bool _testSudakovs;

  /**
   * Whether to do the qtilde distribution output
   */
  bool _qtildeDist;

  /**
   *  The allowed final-state branchings
   */
  map<pair<long,long>,pair<SudakovPtr,IdList> > _allowedFinal;

  /**
   *  The allowed initial-state branchings
   */
  multimap<long, pair<SudakovPtr,IdList> > _allowedInitial;

  /**
   *  The matrix element for the core process
   */
  MEPtr _matrixElement;

  /**
   * The ckkw merging scale
   */
  double _yini;

  /**
   * The fixed alphaS value that was used to generate mad graph events
   */
  double _alphaSMG;

  /**
   * maximum qtilde scale for sudakov interpolation tables
   */
  Energy _max_qtilde;

  /**
   * minimum qtilde scale for sudakov interpolation tables
   */
  Energy _min_qtilde;

  /**                                                    
   * maximum pt cut for sudakov interpolation tables
   */
  Energy _max_pt_cut;

  /** 
   * minimum pt cut for sudakov interpolation tables
   */
  Energy _min_pt_cut;
  
  /**
   * which clustering scheme to use
   */
  unsigned int _clusterOption; 

  /**
   * whether to just reject event if there are no angular ordered histories (true)
   * or to choose one of the non angular ordered histories (false)
   */
  bool _rejectNonAO;

  /**
   * whether to reject events (true) for which there is no parton shower interpretation
   * or just shower them with no reweighting
   */
  bool _rejectNoHistories;

  /**
   * total number of hard trees created in this run
   */
  unsigned int _trees_created;

  /**
   * total number of events with no ordered histories found
   */
  unsigned int _unorderedEvents;

  /**
   * the number of ordered hard trees
   */
  unsigned int _ordered_trees_created;

  /**
   * the maximum multplcity
   */
  unsigned int _max_mult;
  
  /**
   * the minimum multplcity
   */
  unsigned int _min_mult;
  
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
   * The fixed factorization scale used in the MEs.
   */
  Energy _pdfScale;

  /**
   * The first beam.
   */
  Ptr<BeamParticleData>::pointer _beamA;

  /**
   * The second beam.
   */
  Ptr<BeamParticleData>::pointer _beamB;

  /**
   * The veto applied to the shower.
   */
  Ptr<CKKWVeto>::pointer _showerVeto;

  /**
   * Whether to generate the Sudakov reweighting via event vetoes
   */
  bool _dynamicSuds;

  /**
   * The option for applying dynamic resolution cuts to the hardTree configurations
   * to ensure that all partons are resolved at the merging scale.
   */
  int _cutOption;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of PowhegHandler. */
template <>
struct BaseClassTrait<Herwig::PowhegHandler,1> {
  /** Typedef of the first base class of PowhegHandler. */
  typedef Herwig::ShowerHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the PowhegHandler class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::PowhegHandler>
  : public ClassTraitsBase<Herwig::PowhegHandler> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::PowhegHandler"; }
  /**
   * The name of a file containing the dynamic library where the class
   * PowhegHandler is implemented. It may also include several, space-separated,
   * libraries if the class PowhegHandler depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMEHadron.so HwPowhegME.so HwPowhegShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_PowhegHandler_H */
