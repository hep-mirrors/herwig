  // -*- C++ -*-
#ifndef Herwig_Node_H
#define Herwig_Node_H
  //
  // This is the declaration of the Node class.
  //
  //#include "Node.fh"
#include "MergingFactory.fh"
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
      // another constructor for first nodes
    Node(Ptr<MatchboxMEBase>::ptr nodeME,int cutstage,Ptr<Merger>::ptr mh);
      // another constructor for underlying nodes
    Node(Ptr<Node>::ptr deephead,
         Ptr<Node>::ptr head,
         Ptr<SubtractionDipole>::ptr dipol,
         Ptr<MatchboxMEBase>::ptr nodeME,
         int cutstage);
    
    /**
     * The destructor.
     */
    virtual ~Node();
      //@}
    
  public:
      // get children from vector<Ptr<MatchboxMEBase>::ptr>
    void birth(vector<Ptr<MatchboxMEBase>::ptr> vec);
      // recursive setXComb. proStage is the number of clusterings
      // before the projectors get filled.
    void setXComb(tStdXCombPtr xc, int proStage);
      // calculate the dipole and ps approximation
    pair<CrossSection,CrossSection> calcDipandPS(Energy scale);
      // calculate the ps approximation
    CrossSection calcPs(Energy scale);
      // calculate the dipole
    CrossSection calcDip(Energy scale);
      // recursive flush caches and clean up XCombs.
    void flushCaches();
      // recursive clearKinematics
    void clearKinematics();
      // recursive setKinematics
    void setKinematics();
      // recursive generateKinematics using tilde kinematics of the dipoles
    bool generateKinematics(const double *r, int stage,Energy2 shat);
      // generate the kinamatics of the first node
    void  firstgenerateKinematics(const double *r, int stage);
      //return the ME
    const Ptr<MatchboxMEBase>::ptr nodeME()const;
      //return the node ME
    Ptr<MatchboxMEBase>::ptr nodeME();
      //return the parent Node
    Ptr<Node>::ptr parent()const {return theparent;}
      // vector of children nodes created in birth
    vector< Ptr<Node>::ptr > children()const {return thechildren;}
      //pick a random child (flat)
    Ptr<Node>::ptr  randomChild();
      // true if all children show scales above pt
    bool allAbove(Energy pt);
      // true if the node is in the history of other.
    bool isInHistoryOf(Ptr<Node>::ptr other);
      // legsize of the node ME
    int legsize() const;
      // set the first node (first men). only use in factory
    void deepHead(Ptr<Node>::ptr deephead) {theDeepHead = deephead;}
      // return the first node
    Ptr<Node>::ptr deepHead() const {return theDeepHead;}
      // insert nodes to projector vector
    void Projector(double a, Ptr<Node>::ptr pro) {
      pair<double,Ptr<Node>::ptr> p;
      p.first  = a;
      p.second = pro;
      theProjectors.push_back(p);
    }
      // insert nodes to projector vector
    vector< pair <double , Ptr<Node>::ptr > > Projector() {return theProjectors;}
      // returns the dipol of the node.
    Ptr<SubtractionDipole>::ptr dipole() const;
      // set the xcomb of the node
    void xcomb(StdXCombPtr xc) { thexcomb = xc;}
      // return the xcomb
    StdXCombPtr xcomb() const {return thexcomb;}
      // return the current running pt
    Energy runningPt(){return theRunningPt;}
      // set the current running pt
    void runningPt(Energy x){theRunningPt=x;}
      // return the cut stage to cut on merging pt in generate kinematics
    int cutStage() const {return theCutStage;}
      // get the clustersafe map for this node
    SafeClusterMap clusterSafe() const {return clustersafer;}
      // get a vector of the next nodes, ordered in pt (and in parton shower phace space)
    vector<Ptr<Node>::ptr> getNextOrderedNodes(bool normal=true,double hardscalefactor=1.);
      //true if the node is in shower history for a given pt
    bool inShowerPS(Energy hardpt);
      //get the history
    Ptr<Node>::ptr getHistory(bool normal=true,double hardscalefactor=1.);
      //true if node correspond to a subtracted real.
    bool subtractedReal() {return theSubtractedReal;}
      // set if node correspont to a subtracted real.
    void subtractedReal(bool x) { theSubtractedReal = x;}
      //true if node correspond to a virtual contribution.
    bool virtualContribution() { return theVirtualContribution ;}
      // set if node correspont to a virtual contribution.
    void virtualContribution(bool x) {theVirtualContribution = x;}
      //pointer to the merging helper
    Ptr<Merger>::ptr MH()const{return theMergingHelper;}
      // set the merging helper
    void MH(Ptr<Merger>::ptr a){theMergingHelper=a;}
    
    Energy pT()const{return dipole()->lastPt();}
    
  private:
      //the xcomb of the node
    StdXCombPtr thexcomb;
      // the Matrixelement representing this node.
    Ptr<MatchboxMEBase>::ptr thenodeMEPtr;
      // the dipol used to substract
      // and generate kinematics using tilde kinematics
    Ptr<SubtractionDipole>::ptr thedipol;
      // vector of the children node
    vector< Ptr<Node>::ptr > thechildren;
      // the parent node
    Ptr<Node>::ptr theparent;
      // The nodes of the projection stage.
    vector< pair <double,Ptr<Node>::ptr> > theProjectors;
      // The godfather node of whole tree.(Firstnode)
    Ptr<Node>::ptr theDeepHead;
      // The CutStage is number of clusterings which are possible without
      // introducing a merging scale to cut away singularities.
      // -> subtracted MEs have the CutStage 1.
      // -> virtual and normal tree level ME get 0.
    int theCutStage;
      // For [[Emitter,Emission],Spectator] the mapped pair gives
      // information if the first and the second cluster is safe.
    SafeClusterMap clustersafer;
      // the current running pt
    Energy theRunningPt;
      // tell if node belongs to an ordered history
    bool isOrdered;
      // flag to tell if node is subtracted real
    bool theSubtractedReal;
      // flag to tell if node is virtual contribution
    bool theVirtualContribution;
      // the merging helper
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





