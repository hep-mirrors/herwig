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

namespace Herwig {
  class PowhegHandler;
  class PrototypeBranching;
  class ProtoTree;
}

//declaration of thepeg ptr
namespace ThePEG {

  ThePEG_DECLARE_POINTERS(Herwig::PowhegHandler,PowhegHandlerPtr);
  ThePEG_DECLARE_POINTERS(Herwig::PrototypeBranching,PrototypeBranchingPtr);
  ThePEG_DECLARE_POINTERS(Herwig::ProtoTree,ProtoTreePtr);

}

namespace Herwig {

using namespace ThePEG;

  class PrototypeBranching :  public Base {
  public:
    PrototypeBranching(ShowerParticlePtr in) : particle(in) {}
    ShowerParticlePtr particle;
    tSudakovPtr sudakov;
    tPrototypeBranchingPtr parent;
    vector<PrototypeBranchingPtr> children;
    HardBranchingPtr convert();
    PrototypeBranchingPtr reset(PrototypeBranchingPtr,
				map<PrototypeBranchingPtr,PrototypeBranchingPtr> &);
  };

  struct PrototypeTree {
    set<PrototypeBranchingPtr> outgoing;
    set<PrototypeBranchingPtr> incoming;
    set<PrototypeBranchingPtr> currentSpaceLike;
    vector<Energy2> scales;
    HardTreePtr convert();
    map<PrototypeBranchingPtr,PrototypeBranchingPtr> reset();
    DiagPtr diagram;
  };

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
		    _highestMult(false), _testSudakovs(false), _qtildeDist( false ),
		    _yini(0.001), _alphaSMG(0.118), _max_qtilde( 91.2*GeV ), _max_pt_cut( 45.6*GeV ), _min_pt_cut( 0.*GeV ), 
		    _clusterOption( 0 ), _dalitzOn( false ) {}

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

  inline bool highestMult(){
    return _highestMult;
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

  //a list of all clusters used in finding all shower histories
  map< HardBranchingPtr , pair< ShowerParticlePtr, ShowerParticlePtr > > _all_clusters;

  //a set of prototrees (which are a set of HardBranchings)
  //can do find on set to see if it contains a certain hardbranching
  set< ProtoTreePtr > _proto_trees;

  //all possible shower configurations that are angular ordered
  vector< pair< HardTreePtr, double > > _hardTrees;

  //just connect up the progenitors in currentProtoTree
  bool simpleColConnections( ProtoTreePtr currentProtoTree );

  bool simpleColConnections( HardTreePtr currentHardTree );

  
  //recursive method to find all trees
  //does a single clustering on the current tree adding more prototrees if more than one possible branching
  bool fillProtoTrees( map< ShowerParticlePtr, HardBranchingPtr >, 
		       ProtoTreePtr );

  //looks to see if a cluster of the cluster_particles already exists in _all_clusterings
  //if so returns the pointer to that hardBranching if not creates the hardBranchings, adds it
  //to _all_branchings and returns the pointer
  HardBranchingPtr getCluster( pair< ShowerParticlePtr, ShowerParticlePtr >, map< ShowerParticlePtr, HardBranchingPtr > );

  //checks whether a prototree containing the same branchings exists already in _proto_trees in 
  //which case the current tree is a repeat and should be removed (and not recursed)
  bool repeatProtoTree( ProtoTreePtr currentProtoTree );

  /**
   * Clusters the partons and creates a branching history
   * by combining the 2 particles with smallest
   * jet measure out of all allowed pairings until we are left 
   * with \f$q\bar{q}\f$.
   */
  HardTreePtr doClustering( );

  /**
   * Creates all cluster histories and selects one
   * according to its shower probability
   */
  HardTreePtr doClusteringOrdered( );
  
  double Sud( Energy scale, long id, Energy pt_cut );


  HardTreePtr generalClustering();
  BranchingElement allowedFinalStateBranching
  (pair<PrototypeBranchingPtr,PrototypeBranchingPtr> &);
  BranchingElement allowedInitialStateBranching
  (pair<PrototypeBranchingPtr,PrototypeBranchingPtr> &);
  DiagPtr getDiagram(const PrototypeTree &);
  Energy2 hadronJetMeasure(const Lorentz5Momentum & p1,
			   const Lorentz5Momentum & p2,
			   bool final=true);
  void createColourFlow(HardTreePtr,DiagPtr);
  void setBeams(HardTreePtr);
  bool checkTree(HardTreePtr);
  bool checkBranching(HardBranchingPtr);

  /**
   *  Calculate the Sudakov weight
   */
  double sudakovWeight( HardTreePtr );
  
  /**                                                                                                                                                       
   *  Calculate the splitting function weight                                                                                                                          
   */
  double splittingFnWeight( HardTreePtr );


  /**
   * Checks to see that the splitting is allowed.
   */
  bool splittingAllowed( ShowerParticlePtr part_i,
			 ShowerParticlePtr part_j,
			 int qq_pairs);
  
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
   * Fills vectors for dalitz scatter plots (3 jets only)                                                             
   */
  void getDalitz();

  /**
   * Outputs qtilde distributions from the tabulated sudakovs
   * for comparison with those from the shower (single emission)
   */
  void makeQtildeDist();

private:

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
   *  and also the minimum allowed scale of branching type
   */
  multimap< long, pair< Interpolator2d< double, Energy, Energy >::Ptr, Energy > > _fbranchings;

  /**
   *  Pointer to the object calculating the strong coupling
   */
  ShowerAlphaPtr _alphaS;

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
   *  Whether the sudakovs should be tested                                                                    
   */
  bool _testSudakovs;

  /**
   * Whether to do the qtilde distribution output
   */
  bool _qtildeDist;


  vector< pair< double, double > > _dalitz_from_q1;

  vector< pair< double, double > > _dalitz_from_q2;

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
  
  /*
   * The global alphaS factor to ensure all alphaS weights are less than one
   */
  double _global_alphaS_wgt;


  /**`
   * Histogram of sudakov weights
   */
  HistogramPtr _hSud;

  /**
   * Histogram of alphaS weights
   */
  HistogramPtr _halphaS;

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
   * switch for dalitz analysis of hard tree clustering (3 jets only)
   */
  bool _dalitzOn;

  /**
   * total number of hard trees created in this run
   */
  unsigned int _trees_created;

  /**
   * the number of ordered hard trees
   */
  unsigned int _ordered_trees_created;

  /**
   * the maximum multplcity
   */
  unsigned int _max_mult;



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
  static string library() { return "HwPowhegShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_PowhegHandler_H */

