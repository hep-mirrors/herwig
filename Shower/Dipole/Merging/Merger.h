  // -*- C++ -*-
  //
  // Merger.h is a part of Herwig - A multi-purpose Monte Carlo event generator
  // Copyright (C) 2002-2007 The Herwig Collaboration
  //
  // Herwig is licenced under version 2 of the GPL, see COPYING for details.
  // Please respect the MCnet academic guidelines, see GUIDELINES for details.
  //
#ifndef HERWIG_Merger_H
#define HERWIG_Merger_H
  //
  // This is the declaration of the Merger class.
  //
#include "MFactory.fh"
#include "Node.fh"



#include "ThePEG/Handlers/HandlerBase.h"
#include "ThePEG/Vectors/SpinOneLorentzRotation.h"
#include "Herwig/Shower/Dipole/DipoleShowerHandler.h"
#include "Herwig/Shower/Dipole/Base/DipoleSplittingGenerator.h"
#include "Herwig/MatrixElement/Matchbox/Base/MergerBase.h"

#include "ThePEG/Cuts/JetFinder.h"
#include "ThePEG/Cuts/Cuts.h"



namespace Herwig {
  
  using namespace ThePEG;
  
  typedef Ptr<Node>::ptr NPtr;
  typedef vector<Ptr<Node>::ptr> NPtrVec;
  struct HistoryStep {
    
    HistoryStep(){}
    
    HistoryStep(NPtr cn,double w ,Energy sc){
      node=cn;
      weight=w;
      scale=sc;
    }
    
    NPtr node;
    double weight;
    Energy scale;
  };
  
  typedef vector< HistoryStep > History;
  
  typedef multimap<DipoleIndex,Ptr<DipoleSplittingGenerator>::ptr> GeneratorMap2;
  
  /**
   * \ingroup DipoleShower
   * \author Johannes Bellm
   *
   * \brief Merger handles the Merger ....... //TODO .
   *
   * @see \ref MergerInterfaces "The interfaces"
   * defined for Merger.
   */
  class Merger: public MergerBase {
    
  public:
    
    /** @name Standard constructors and destructors. */
      //@{
    /**
     * The default constructor.
     */
    Merger();
    
    /**
     * The destructor.
     */
    virtual ~Merger();
      //@}
    
  public:
    
    
    double singlesudakov(list<Dipole>::iterator,Energy,Energy,pair<bool,bool>);
    bool   dosudakov(NPtr Born,Energy  running, Energy next, double& sudakov0_n);
  
    //bool   dosudakovold(NPtr, Energy, Energy, double&);
    //bool   sudakov(NPtr Born, Energy  running, Energy next);
    
    
    void   cleanup(NPtr);
    
    bool   projectorStage(NPtr);
    Energy CKKW_StartScale(NPtr);
    void   CKKW_PrepareSudakov(NPtr,Energy);
    CrossSection TreedSigDR(Energy startscale,NPtr,double diffalpha=1.);
    CrossSection LoopdSigDR(Energy startscale,NPtr);
    bool   fillProjector(Energy&);
    void   fillHistory(Energy, NPtr, NPtr );
    
    
    double sumpdfReweightUnlops();
    double sumalphaReweightUnlops();
    double pdfratio(NPtr,Energy,Energy,int);
    double pdfReweight();
    double alphaReweight();
    
    double sumfillHistoryUnlops();
    
    
    
    double alphasUnlops( Energy next,Energy fixedScale);
    double pdfUnlops(tcPDPtr,tcPDPtr,tcPDFPtr,Energy,Energy,double,int,Energy);
    double singleUNLOPS(list<Dipole>::iterator,Energy,Energy,Energy,pair<bool,bool>);
    bool   doUNLOPS(NPtr Born,Energy  running, Energy next,Energy fixedScale, double& UNLOPS);
    
    
    void setME(Ptr<MatchboxMEBase>::ptr me){theCurrentME=me;}
    
    CrossSection MergingDSigDR() ;
    CrossSection MergingDSigDRBornStandard(NPtr Node);
    CrossSection MergingDSigDRBornGamma(NPtr Node);

    CrossSection MergingDSigDRVirtualStandard(NPtr Node);
    
    CrossSection MergingDSigDRRealStandard(NPtr Node);
    
    CrossSection MergingDSigDRRealAllAbove(NPtr Node);
    CrossSection MergingDSigDRRealBelowSubReal(NPtr Node);
    CrossSection MergingDSigDRRealBelowSubInt(NPtr Node);
    
    
    double as(Energy q){return DSH()->as(q);}
    Ptr<DipoleShowerHandler>::ptr DSH(){return theDipoleShowerHandler;}
    
    
    size_t maxLegs() const {return theCurrentMaxLegs;}
    size_t maxLegsLO() const {return N0()+N();}
    size_t maxLegsNLO()const {return N0()+M()+1;}
    
    
    
    Ptr<MFactory>::ptr treefactory();
    
    map<Ptr<MatchboxMEBase>::ptr,NPtr> firstNodeMap();
    
    
    void setXComb(Ptr<MatchboxMEBase>::ptr,tStdXCombPtr,int);
    void setKinematics(Ptr<MatchboxMEBase>::ptr);
    void clearKinematics(Ptr<MatchboxMEBase>::ptr);
    bool generateKinematics(Ptr<MatchboxMEBase>::ptr,const double *);
    bool calculateInNode() const;
    void fillProjectors(Ptr<MatchboxMEBase>::ptr);
    pair<bool,bool> clusterSafe(Ptr<MatchboxMEBase>::ptr,int,int,int);
    bool theCalculateInNode;
    
    
        bool matrixElementRegion(PVector particles,Energy winnerScale=0.*GeV,Energy cutscale=0.*GeV);
    
    
    
    
    
    void renormscale(Energy x) {
      therenormscale = x;
    }
    
    Energy renormscale() {
      return therenormscale;
    }
    
    
    Energy mergingScale()const{
      return theMergePt;
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
    
    
    
    void firstNodeMap(Ptr<MatchboxMEBase>::ptr,NPtr);
   

    
    void smeareMergePt(){theMergePt=centralMergePt()*(1.+0.*(-1.+2.*UseRandom::rnd())*smear());}
    
    
    double gamma()const{return theGamma;}
    
    
    
    
    double smear()const{return theSmearing;}
    
    
    bool MERegionByJetAlg(){return defMERegionByJetAlg;}
    
    
    Ptr<ColourBasis>::ptr largeNBasis(){return theLargeNBasis;}
    
    
    void largeNBasis(Ptr<ColourBasis>::ptr x){theLargeNBasis=x;}
    
    int M()const{return theM;}
    
    int N()const{return theN;}
    
    int N0()const{return theN0;}
    void N0(int n){ theN0=n;}
    
    bool openInitialStateZ(){return theOpenInitialStateZ;}
    
    
  private:
    
    /**
     * Variations
     */
    bool minusL;
    bool Unlopsweights;
    bool theCMWScheme;
    bool projected;
    bool isUnitarized;
     bool isNLOUnitarized;
    bool defMERegionByJetAlg;
    bool theOpenInitialStateZ;
     int theChooseHistory;
    int theN0;
     int  theOnlyN;
    int theN;
    int theM;
    size_t theCurrentMaxLegs;
    double   xiRenME;
    double   xiFacME;
    double   xiRenSh;
    double   xiFacSh;
    double   xiQSh;
    double theNf;
    double weight,weightCB;
    double theGamma;
    double ee_ycut;
    double pp_dcut;
     double theSmearing;
    
    
    Energy therenormscale;
    /**
     * If any clustering below the CutStage has a lower pT than the MergePt
     * the phase space point has to be thrown away.
     */
    Energy theIRSafePT;
     Energy theMergePt;
     Energy theCentralMergePt;

    History history;
    Ptr<JetFinder>::ptr theMergingJetFinder;
    Ptr<ColourBasis>::ptr theLargeNBasis;
    
    Ptr<MatchboxMEBase>::ptr theCurrentME;
    
    

    
    /**
     * The mean of the Gaussian distribution for
     * the intrinsic pt of sea partons.
     */
    
    Ptr<DipoleShowerHandler>::ptr theDipoleShowerHandler;
    
    
    Ptr<MFactory>::ptr theTreeFactory;
    map<Ptr<MatchboxMEBase>::ptr,NPtr> theFirstNodeMap;
    
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
    Merger & operator=(const Merger &);
    
  };
  
}



#endif /* HERWIG_Merger_H */
