// -*- C++ -*-
#ifndef HERWIG_PowhegHandler_H
#define HERWIG_PowhegHandler_H
//
// This is the declaration of the PowhegHandler class.
//

#include "Herwig++/Shower/ShowerHandler.h"
#include "Herwig++/Shower/Base/Evolver.h"
#include "Herwig++/Utilities/Interpolator.h"
#include "Herwig++/Utilities/GaussianIntegrator.h"
#include "Herwig++/Shower/Base/ShowerParticle.h"
#include "ThePEG/Config/Pointers.h"
#include "Herwig++/Shower/Powheg/HardTree.h"
#include "Herwig++/Shower/Couplings/ShowerAlpha.h"
#include "ThePEG/MatrixElement/MEBase.fh"


namespace Herwig {

  class PowhegHandler;
  class PrototypeBranching;

}

//declaration of thepeg ptr
namespace ThePEG {

  ThePEG_DECLARE_POINTERS(Herwig::PowhegHandler,PowhegHandlerPtr);
  ThePEG_DECLARE_POINTERS(Herwig::PrototypeBranching,PrototypeBranchingPtr);

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
  PowhegHandler() : _npoint(200), _sudopt(0), _sudname("sudakov.data"),
		    _jetMeasureMode(0),_lepton(true)
		  , _yini(0.001), _alphaSMG(0.118) {}

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
  virtual void doinit() throw(InitException);

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
   * Clusters the partons and creates a branching history
   * by combining the 2 particles with smallest
   * jet measure out of all allowed pairings until we are left 
   * with \f$q\bar{q}\f$.
   */
  HardTreePtr doClustering();

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
  double sudakovWeight();
  
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
  SudakovPtr getSud(int & qq_pairs, long & emmitter_id,
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

private:

  /**
   *  The hard tree
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
  multimap< long, pair < Interpolator<double,Energy>::Ptr,
			 Interpolator<Energy,double>::Ptr>  > _fbranchings;
  
  /**
   * Map containing particle and the resolution parameter of 
   * external partons.
   */
  map< ShowerParticlePtr, pair< double, Energy > > _theExternals;

  /**
   * Map containing all nodes with the ingoing partons and their clustered jet 
   * parameter (this is the the intermediates and their ending node).
   */
  map< HardBranchingPtr, double > _theNodes;
  
  /**
   * Map containing all intermediate partons and 
   * their start and end node resolution parameters ( y and qtilde ).
   */
  map< long, pair< pair< double, Energy > , pair< double, Energy > > > 
  _theIntermediates;

  /**
   * Map between the cluster number and the jet parameter.
   */
  map< int, double > _theJetParameters;
  
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


