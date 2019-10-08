  // -*- C++ -*-
#ifndef Herwig_Node_H
#define Herwig_Node_H
  //
  // This is the declaration of the Node class.
  //
#include "Node.fh"
#include "MergingFactory.fh"
#include "Merger.h"



#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/Config/std.h"
#include "ThePEG/Interface/Interfaced.h"
#include "Herwig/MatrixElement/Matchbox/Base/MatchboxMEBase.fh"
#include "Herwig/MatrixElement/Matchbox/Dipoles/SubtractionDipole.fh"
#include "Herwig/Shower/Dipole/Base/DipoleEventRecord.h"
#include "ThePEG/MatrixElement/MEBase.h"

#include <vector>



namespace Herwig {
    
    
    using namespace ThePEG;
  
  /**
   * Here is the documentation of the Node class.
   *
   * @see \ref NodeInterfaces "The interfaces"
   * defined for Node.
   */
  
  
  class Node : public Interfaced {
		public:
    
    /** @name Standard constructors and destructors. */
      //@{
    
    Node (){};
    
      // constructor for first nodes
    Node(MatchboxMEBasePtr nodeME , int cutstage , MergerPtr mh );
      // another constructor for underlying nodes
    Node(NodePtr deephead, 
         NodePtr head, 
         SubtractionDipolePtr dipol, 
         MatchboxMEBasePtr nodeME, 
         int cutstage);
      /// The destructor.
    virtual ~Node();
      //@}
    
  public:
      // get children from vector<MatchboxMEBasePtr>
    void birth(const vector<MatchboxMEBasePtr> & vec);
      /// recursive setXComb. proStage is the number of clusterings
      /// before the projectors get filled.
    void setXComb(tStdXCombPtr xc);
      /// calculate the dipole and ps approximation
    pair<CrossSection, CrossSection> calcDipandPS(Energy scale)const;
      /// calculate the ps approximation
    CrossSection calcPs(Energy scale)const;
      /// calculate the dipole
    CrossSection calcDip(Energy scale)const;
      /// recursive flush caches and clean up XCombs.
    void flushCaches();
      /// recursive clearKinematics
    void clearKinematics();
      /// recursive setKinematics
    void setKinematics();
      /// recursive generateKinematics using tilde kinematics of the dipoles
    bool generateKinematics(const double *r, bool directCut);
      /// generate the kinamatics of the first node
    bool firstgenerateKinematics(const double *r, bool directCut);
      //return the ME
    const MatchboxMEBasePtr nodeME() const;
      //return the node ME
    MatchboxMEBasePtr nodeME();
      //return the parent Node
    NodePtr parent() const {return theparent;}
      /// vector of children nodes created in birth
    vector< NodePtr > children() const {return thechildren;}
      //pick a random child (flat)
    NodePtr  randomChild();
      /// true if all children show scales above pt
    bool allAbove(Energy pt);
      /// return maximum of all child pts.
    Energy maxChildPt();  
      /// true if the node is in the history of other.
    bool isInHistoryOf(NodePtr other);
      /// legsize of the node ME
    int legsize() const;
      /// set the first node (first men). only use in factory
    void deepHead(NodePtr deephead) {theDeepHead = deephead;}
      /// return the first node
    NodePtr deepHead() const {return theDeepHead;}
      /// returns the dipol of the node.
    SubtractionDipolePtr dipole() const;
      /// return the xcomb
    StdXCombPtr xcomb() const;
      /// return the xcomb (if not created, create one from head)
    StdXCombPtr xcomb() ;
      /// return the current running pt
    Energy runningPt() const { return theRunningPt; }
      /// set the current running pt
    void runningPt(Energy x) { theRunningPt=x; }
      /// return the cut stage to cut on merging pt in generate kinematics
    int cutStage() const { return theCutStage; }
      /// get a vector of the next nodes, ordered in pt (and in parton shower phace space)
    vector<NodePtr> getNextOrderedNodes(bool normal=true, double hardscalefactor=1.) const;
      //true if the node is in shower history for a given pt
    bool inShowerPS(Energy hardpt)const;
      //get the history
    NodePtr getHistory(bool normal=true, double hardscalefactor=1.);
      //true if node correspond to a subtracted real.
    bool subtractedReal() const {return theSubtractedReal;}
      /// set if node correspont to a subtracted real.
    void subtractedReal(bool x) { theSubtractedReal = x;}
      //true if node correspond to a virtual contribution.
    bool virtualContribution() const { return theVirtualContribution ;}
      /// set if node correspont to a virtual contribution.
    void virtualContribution(bool x) {theVirtualContribution = x;}
      //pointer to the merging helper
    MergerPtr MH()const{return theMergingHelper;}
      /// set the merging helper
    void MH(MergerPtr a){theMergingHelper=a;}
      ///  pT of the dipole
    Energy pT()const{return dipole()->lastPt();}
      /// get incoming and outgoing particles (TODO: expensive)
    pair<PVector , PVector> getInOut();
    
  private:
      /// the Matrixelement representing this node.
    MatchboxMEBasePtr thenodeMEPtr;
      /// the dipol used to substract
      /// and generate kinematics using tilde kinematics
    SubtractionDipolePtr thedipol;
      /// the parent node
    NodePtr theparent;
      /// The godfather node of whole tree.(Firstnode)
    NodePtr theDeepHead;
      /**
       * The CutStage is number of clusterings which are possible without
       * introducing a merging scale to cut away singularities.
       * -> subtracted MEs have the CutStage 1.
       * -> virtual and normal tree level ME get 0.
       */
    int theCutStage;
      /// tell if node belongs to an ordered history
    // bool isOrdered;
      /// flag to tell if node is subtracted real
    bool theSubtractedReal;
      /// flag to tell if node is virtual contribution
    bool theVirtualContribution;
      /// the merging helper
    MergerPtr theMergingHelper;
      //the xcomb of the node
    StdXCombPtr thexcomb;
      /// vector of the children node
    vector< NodePtr > thechildren;
      /// the current running pt
    Energy theRunningPt;
      /// The nodes of the projection stage.
    NodePtr theProjector;
      /// flag not to enter infinite loop. (There should be a better solution...)
    bool didflush=false;
    
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
    
    
      // If needed, insert declarations of virtual function defined in the
      // InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).
    
    
    
  private:
    /**
     * The assignment operator is private and must never be called.
     * In fact, it should not even be implemented.
     */
    Node & operator=(const Node &) = delete;
    
  };
  
}

#endif /* Herwig_Node_H */





