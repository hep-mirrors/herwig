  // -*- C++ -*-
#ifndef Herwig_Node_H
#define Herwig_Node_H
  //
  // This is the declaration of the Node class.
  //
  //#include "Node.fh"
#include "MFactory.fh"
#include "Merger.h"


#include "ThePEG/Interface/Interfaced.h"
#include "Herwig/MatrixElement/Matchbox/Base/MatchboxMEBase.fh"
#include "Herwig/MatrixElement/Matchbox/Dipoles/SubtractionDipole.fh"
#include "Herwig/Shower/Dipole/Base/DipoleEventRecord.h"
#include "ThePEG/MatrixElement/MEBase.h"
#include "ThePEG/Config/Pointers.h"
#include <vector>

namespace Herwig {
  
  using namespace ThePEG;
  
  /**
   * Here is the documentation of the Node class.
   *
   * @see \ref NodeInterfaces "The interfaces"
   * defined for Node.
   */
  
  
  /**
   * Define the SafeClusterMap type map<pair<pair<emitter,emmision>,spectator >
   *                                    ,pair<first-clustering,second-clustering> >
   */
  typedef map<pair<pair<int,int>,int >,pair<bool,bool> > SafeClusterMap;
  
  
  
 
  
  
  
  
  
  
  
  class Node : public Interfaced {
		public:
    
    /** @name Standard constructors and destructors. */
      //@{
    
    /**
     * The default constructor. Do not use!
     */
    Node();
    
    Node(Ptr<MatchboxMEBase>::ptr nodeME,  int deepprostage,int cutstage);
    
    Node(Ptr<Node>::ptr deephead,
         Ptr<Node>::ptr head,
         Ptr<SubtractionDipole>::ptr dipol,
         Ptr<MatchboxMEBase>::ptr nodeME,
         int deepprostage, int cutstage);
    
    /**
     * The destructor.
     */
    virtual ~Node();
      //@}
    
		public:
    
    /** get children from vector<Ptr<MatchboxMEBase>::ptr>  */
    
    void birth(vector<Ptr<MatchboxMEBase>::ptr> vec);
    
    /** recursive setXComb. proStage is the number of clusterings
     * before the projectors get filled. */
    
    void setXComb(tStdXCombPtr xc, int proStage);
    
    
    bool headContribution(double hardscalefactor);
    
    bool DipolesAboveMergeingScale(Ptr<Node>::ptr& selectedNode,double& sum,Energy& minpt,int& number);
 
    /*   
    bool diffPsDipBelowMergeingScale(Ptr<Node>::ptr& selectedNode,double & sum,Energy& minpt,int& number);
    
    bool psBelowMergeingScale(Ptr<Node>::ptr& selectedNode,double & sum,Energy& minpt,int& number);
    
    bool dipBelowMergeingScale(Ptr<Node>::ptr& selectedNode,double & sum,Energy& minpt,int& number);
    */
    
    pair<CrossSection,CrossSection> calcDipandPS(Energy scale);
    CrossSection calcPs(Energy scale);
    
    /** recursive flush caches and clean up XCombs. */
    
    void flushCaches();
    
    /** recursive clearKinematics */
    
    void clearKinematics();
    
    
    /** recursive setKinematics */
    
    void setKinematics();
    
    /** recursive generateKinematics using tilde kinematics of the dipoles*/
    
    bool generateKinematics(const double *r, int stage,Energy2 shat);
    
    void  firstgenerateKinematics(const double *r, int stage,Energy2 shat);
    
    /** returns the matrix element pointer */
    
    const Ptr<MatchboxMEBase>::ptr nodeME()const;
    
    /** access the matrix element pointer */
    
    Ptr<MatchboxMEBase>::ptr nodeME();
    
    /** the parent node */
    
    Ptr<Node>::ptr parent()const {
      return theparent;
    }
    
    /** map of children nodes*/
    
    vector< Ptr<Node>::ptr > children() {
      return thechildren;
    }
    
    Ptr<Node>::ptr  randomChild();
  
    Energy miniPt()const ;  
    
    bool allAbove(Energy pt);
    
    bool isInHistoryOf(Ptr<Node>::ptr other);
    
    int legsize() const;
    
    
    /** set the first node (godfather). only use in factory*/
    
    void deepHead(Ptr<Node>::ptr deephead) {
      theDeepHead = deephead;
    }
    
    /** return the first node*/
    
    Ptr<Node>::ptr deepHead() const {
      return theDeepHead;
    }
    
    /** insert nodes to projector vector */
    
    void Projector(double a, Ptr<Node>::ptr pro) {
      pair<double,Ptr<Node>::ptr> p;
      p.first  = a;
      p.second = pro;
      theProjectors.push_back(p);
      
        //theProjectors.push_back(make_pair(a,pro));
    }
    
    /** insert nodes to projector vector */
    
    vector< pair <double , Ptr<Node>::ptr > > Projector() {
      return theProjectors;
    }
    /** returns the dipol of the node. */
    
    Ptr<SubtractionDipole>::ptr dipol() const;
    
    /** set the xcomb of the node */
    
    void xcomb(StdXCombPtr xc) {
      thexcomb = xc;
    }
    
    /** return the xcomb */
    
    StdXCombPtr xcomb() const {
      return thexcomb;
    }
    /** set the headxcomb of the node */
    
    void headxcomb(StdXCombPtr xc) {
      thexcomb = xc;
    }
    
    /** return the xcomb */
    
    StdXCombPtr headxcomb() const {
      return theheadxcomb;
    }
    

    
    Energy vetoPt() const {
      return theVetoPt;
    }
    
    void vetoPt(Energy x) {
      theVetoPt=x;
    }
    
    
    
    Energy runningPt(){return theRunningPt;}
    void runningPt(Energy x){theRunningPt=x;}
    
    int cutStage() const {
      return theCutStage;
    }

    unsigned int DeepProStage() const {
      return theDeepProStage;
    }
    
  
    
    SafeClusterMap clusterSafe() const {
      return clustersafer;
    }
    
    bool ordered() {
      return isOrdered;
    }
    
    int orderedSteps() const {
      return theOrderdSteps;
    }
    
    void setOrderedSteps(int x) {
      theOrderdSteps = x;
    }
    
    bool isThisSafe() const {
      return isthissafe;
    }
    
    vector<Ptr<Node>::ptr> getNextOrderedNodes(bool normal=true,double hardscalefactor=1.);
    
    bool inShowerPS(Energy hardpt);
    
    
    vector<Ptr<Node>::ptr> getNextFullOrderedNodes();
    
    bool hasFullHistory();
    
    Ptr<Node>::ptr getHistory(bool normal=true,double hardscalefactor=1.);
    
    bool subtractedReal() {
      return theSubtractedReal;
    }
    void subtractedReal(bool x) {
      theSubtractedReal = x;
    }
    
    bool virtualContribution() {
      return theVirtualContribution ;
    }
    
    
    void virtualContribution(bool x) {
      theVirtualContribution = x;
    }
    
    Ptr<Merger>::ptr MH()const{return theMergingHelper;}
    
     void MH(Ptr<Merger>::ptr a){theMergingHelper=a;}
    

    
		private:
    StdXCombPtr theheadxcomb;
    
    StdXCombPtr thexcomb;
    
    /** the Matrixelement representing this node. */
    
    Ptr<MatchboxMEBase>::ptr thenodeMEPtr;
    
    /** the dipol used to substract
     *  and generate kinematics using tilde kinematics */
    
    Ptr<SubtractionDipole>::ptr thedipol;
    
    /** vector of the children node*/
    
    vector< Ptr<Node>::ptr > thechildren;
    
    /** the parent node*/
    
    Ptr<Node>::ptr theparent;
    
    /** The nodes of the projection stage.    */
    
    vector< pair <double,Ptr<Node>::ptr> > theProjectors;
    
    
    /** The godfather node of whole tree.(Firstnode) */
    
    Ptr<Node>::ptr theDeepHead;
    
    /**
     * The CutStage is number of clusterings which are possible without
     * introducing a merging scale to cut away singularities.
     * -> subtracted MEs have the CutStage 1.
     * -> virtual and normal tree level ME get 0.
     */
    
    int theCutStage;
    
    
    /**
     * all nodes downstream have pt over merged pt.
     */
    
    bool isthissafe;

    
    
    
    /**
     * The DeepProStage is the number of additional partons of the DeepHead
     * compared to the lowest multiplicity.
     */
    
    int theDeepProStage;
    
    
    /**
     * For [[Emitter,Emission],Spectator] the mapped pair gives
     * information if the first and the second cluster is safe.
     */
    
    SafeClusterMap clustersafer;
    
    
    bool isOrdered;
    
    int theOrderdSteps;
    
    bool theNeedFullOrderedHistory;
    
    unsigned int theSudakovSteps;
    
    Energy theVetoPt;
    
    Energy theRunningPt;
    
    bool theSubtractedReal;
    
    bool theVirtualContribution;
    
    
     Ptr<Merger>::ptr theMergingHelper;
    
    
    
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
    Node & operator=(const Node &);
    
  };
  
}

#endif /* Herwig_Node_H */





