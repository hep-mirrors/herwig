  // -*- C++ -*-
#ifndef Herwig_ClusterNode_H
#define Herwig_ClusterNode_H
  //
  // This is the declaration of the ClusterNode class.
  //
  //#include "ClusterNode.fh"
#include "MergeFactory.fh"

#include "ThePEG/Interface/Interfaced.h"
#include "Herwig/MatrixElement/Matchbox/Base/MatchboxMEBase.fh"
#include "Herwig/MatrixElement/Matchbox/Dipoles/SubtractionDipole.fh"
#include "Herwig/DipoleShower/Base/DipoleEventRecord.h"
#include "ThePEG/MatrixElement/MEBase.h"
#include "ThePEG/Config/Pointers.h"
#include <vector>

namespace Herwig {
  
  using namespace ThePEG;
  
  /**
   * Here is the documentation of the ClusterNode class.
   *
   * @see \ref ClusterNodeInterfaces "The interfaces"
   * defined for ClusterNode.
   */
  
  
  /**
   * Define the SafeClusterMap type map<pair<pair<emitter,emmision>,spectator >
   *                                    ,pair<first-clustering,second-clustering> >
   */
  typedef map<pair<pair<int,int>,int >,pair<bool,bool> > SafeClusterMap;
  
  
  
  
  template <typename T>
  struct ClusterNodecounter
  {
    ClusterNodecounter()
    {
      objects_created++;
      objects_alive++;
    }
    
    virtual ~ClusterNodecounter()
    {
      --objects_alive;
    }
    static int objects_created;
    static int objects_alive;
  };
  template <typename T> int ClusterNodecounter<T>::objects_created( 0 );
  template <typename T> int ClusterNodecounter<T>::objects_alive( 0 );
  
  
  
  
  
  
  
  
  
  class ClusterNode : public Interfaced,ClusterNodecounter<ClusterNode> {
		public:
    
    /** @name Standard constructors and destructors. */
      //@{
    
    /**
     * The default constructor. Do not use!
     */
    ClusterNode();
    
    ClusterNode(Ptr<MatchboxMEBase>::ptr nodeME,  int deepprostage,int cutstage, bool nFOH);
    
    ClusterNode(Ptr<ClusterNode>::ptr deephead, Ptr<ClusterNode>::ptr head, Ptr<SubtractionDipole>::ptr dipol, Ptr<MatchboxMEBase>::ptr nodeME,
                int deepprostage, int cutstage);
    
    /**
     * The destructor.
     */
    virtual ~ClusterNode();
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
    
    
		public:
    
    /** get children from vector<Ptr<MatchboxMEBase>::ptr>  */
    
    void birth(vector<Ptr<MatchboxMEBase>::ptr> vec);
    
    /** recursive setXComb. proStage is the number of clusterings
     * before the projectors get filled. */
    
    void setXComb(tStdXCombPtr xc, int proStage);
    
    
    
    
    double setProjectorStage(bool fast=false);
    
    bool headContribution(double hardscalefactor);
    
    bool DipolesAboveMergeingScale(Ptr<ClusterNode>::ptr& selectedNode,double& sum,Energy& minpt,int& number);
    
    bool diffPsDipBelowMergeingScale(Ptr<ClusterNode>::ptr& selectedNode,double & sum,Energy& minpt,int& number);
    
    bool psBelowMergeingScale(Ptr<ClusterNode>::ptr& selectedNode,double & sum,Energy& minpt,int& number);
    
    bool dipBelowMergeingScale(Ptr<ClusterNode>::ptr& selectedNode,double & sum,Energy& minpt,int& number);
    
    
    
    
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
    
    Ptr<ClusterNode>::ptr parent()const {
      return theparent;
    }
    
    /** map of children nodes*/
    
    vector< Ptr<ClusterNode>::ptr > children() {
      return thechildren;
    }
    
    /** set the first node (godfather). only use in factory*/
    
    void deepHead(Ptr<ClusterNode>::ptr deephead) {
      theDeepHead = deephead;
    }
    
    /** return the first node*/
    
    Ptr<ClusterNode>::ptr deepHead() const {
      return theDeepHead;
    }
    
    /** insert nodes to projector vector */
    
    void Projector(double a, Ptr<ClusterNode>::ptr pro) {
      pair<double,Ptr<ClusterNode>::ptr> p;
      p.first  = a;
      p.second = pro;
      theProjectors.push_back(p);
    }
    
    /** insert nodes to projector vector */
    
    vector< pair <double , Ptr<ClusterNode>::ptr > > Projector() {
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
    
    Energy mergePt() {
      return theMergePt;
    }
    
    void mergePt(Energy x) {
      theMergePt = x;
    }
    
    Energy centralMergePt() {
      return theCentralMergePt;
    }
    
    void centralMergePt(Energy x) {
      theCentralMergePt = x;
    }
    
    double smear();
    
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
    
    void N(unsigned int N) {
      theN = N;
    }
    
    void N0(unsigned int N) {
      theN0 = N;
    }
    
    
    void onlyN( int N) {
      theOnlyN = N;
    }
    int onlyN( ) {
      return theOnlyN ;
    }
    
    
    
    unsigned int N() const {
      return theN;
    }
    
    unsigned int N0() const {
      return theN0;
    }
    
    void M(unsigned int M) {
      theM = M;
    }
    unsigned int M() const {
      return theM;
    }
    
    unsigned int sudakovSteps() const {
      return theSudakovSteps;
    }
    
    unsigned int DeepProStage() const {
      return theDeepProStage;
    }
    
    void irsavePt(Energy x) {
      theIRsafePt = x;
    }
    
    Energy irsavePt() {
      return theIRsafePt;
    }
    
    void irsaveRatio(double x) {
      theIRCsafeRatio = x;
    }
    
    double irsaveRatio() {
      return theIRCsafeRatio;
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
    
    vector<Ptr<ClusterNode>::ptr> getNextOrderedNodes(bool normal=true,double hardscalefactor=1.);
    
    bool inShowerPS(Energy hardpt);
    
    Energy newHardPt();
    
    vector<Ptr<ClusterNode>::ptr> getNextFullOrderedNodes();
    
    bool hasFullHistory();
    
    bool unitarized() const {
      return isUnitarized;
    }
    void unitarized(bool is) {
      isUnitarized = is;
    }
    
    bool NLOunitarized() const {
      return isNLOUnitarized;
    }
    void NLOunitarized(bool is) {
      isNLOUnitarized = is;
    }
    
    bool needFullOrderedHistory() {
      return theNeedFullOrderedHistory;
    }
    
    Ptr<ClusterNode>::ptr getLongestHistory(bool normal=true);
    Ptr<ClusterNode>::ptr getLongestHistory_simple(bool normal=true,double hardscalefactor=1.);
    
    
    tSubProPtr showeredSub(){return theShoweredSub;}
    
    void showeredSub( tSubProPtr x){theShoweredSub=x;}
    
    DipoleEventRecord& eventRec(){return theEventRec;}
    void eventRec(DipoleEventRecord& x){theEventRec=x;}
    
    bool VetoedShower(){return needsVetoedShower;}
    void VetoedShower(bool x){needsVetoedShower = x;}
    bool calculateInNode(){return theCalculateInNode;}
    void calculateInNode(bool x) {
      theCalculateInNode = x;
    }
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
    
    
    bool finiteDipoles() {
      return thefiniteDipoles;
    }
    void finiteDipoles(bool x) {
      thefiniteDipoles = x;
    }
    
    
    
    bool subCorrection() {
      return thesubCorrection;
    }
    void subCorrection(bool x) {
      thesubCorrection = x;
    }
    
    
    
    void renormscale(Energy x) {
      therenormscale = x;
    }
    
    Energy renormscale() {
      return therenormscale;
    }
    
    void chooseHistory(int x){
      theChooseHistory=x;
    }
    
    int chooseHistory(){
      return theChooseHistory;
    }
    Ptr<MergeFactory>::ptr treefactory();
    
			 void treefactory(Ptr<MergeFactory>::ptr x);
			 
			 
			 void addclusternumber(){theclusteredto++;}
    int  clusternumber(){return theclusteredto;}
    
    
		private:
    StdXCombPtr theheadxcomb;
    
    StdXCombPtr thexcomb;
    
    /** the Matrixelement representing this node. */
    
    Ptr<MatchboxMEBase>::ptr thenodeMEPtr;
    
    /** the dipol used to substract
     *  and generate kinematics using tilde kinematics */
    
    Ptr<SubtractionDipole>::ptr thedipol;
    
    /** vector of the children node*/
    
    vector< Ptr<ClusterNode>::ptr > thechildren;
    
    /** the parent node*/
    
    Ptr<ClusterNode>::ptr theparent;
    
    /** The nodes of the projection stage.    */
    
    vector< pair <double,Ptr<ClusterNode>::ptr> > theProjectors;
    
    
    /** The godfather node of whole tree.(Firstnode) */
    
    Ptr<ClusterNode>::ptr theDeepHead;
    
    /**
     * The CutStage is number of clusterings which are possible without
     * introducing a merging scale to cut away singularities.
     * -> subtracted MEs have the CutStage 1.
     * -> virtual and normal tree level ME get 0.
     */
    
    int theCutStage;
    
    
    /**
     * If any clustering below the CutStage has a lower pT than the MergePt
     * the phase space point has to be thrown away.
     */
    
    static Energy theMergePt;
    
    
    static Energy theCentralMergePt;
    
    static double smearing;
    
    /**
     * all nodes downstream have pt over merged pt.
     */
    
    bool isthissafe;
    
    /**
     * If any clustering below the CutStage has a lower pT than the MergePt
     * the phase space point has to be thrown away.
     */
    
    static Energy theIRsafePt;
    
    /**
     * If any clustering below the CutStage has a lower (pT/shat)^2 than this value
     * the phase space point has to be thrown away.
     */
    
    static double theIRCsafeRatio;
    
    
    
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
    
    static unsigned int theN0;
    
    static int  theOnlyN;
    
    static unsigned int theN;
    
    static unsigned int theM;
    
    unsigned int theSudakovSteps;
    
    static bool isUnitarized;
    
    static bool isNLOUnitarized;
    
    Energy theVetoPt;
    
    Energy theRunningPt;
    
    tSubProPtr theShoweredSub;
    
    DipoleEventRecord theEventRec;
    
    bool needsVetoedShower;
    
    bool theCalculateInNode;
    
    
    
    bool theSubtractedReal;
    
    bool theVirtualContribution;
    
    bool thefiniteDipoles;
    
    bool thesubCorrection;
    
    int theclusteredto;
    
    
    Energy therenormscale;
    
    Ptr<MergeFactory>::ptr theTreeFactory;
    
    
      // 1 == always take one of the smallest Pts.
      // 2 == always take one of the highest  Pts.
      // 3 == choose correspondingly to the dipolweight.
    
    static int theChooseHistory;
    
    
    
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
    ClusterNode & operator=(const ClusterNode &);
    
  };
  
}

#endif /* Herwig_ClusterNode_H */





