  /// -*- C++ -*-
  //
  /// Merger.h is a part of Herwig - A multi-purpose Monte Carlo event generator
  /// Copyright (C) 2002-2019 The Herwig Collaboration
  //
  /// Herwig is licenced under version 3 of the GPL, see COPYING for details.
  /// Please respect the MCnet academic guidelines, see GUIDELINES for details.
  //
#ifndef HERWIG_Merger_H
#define HERWIG_Merger_H
  //
  /// This is the declaration of the Merger class.
  //
#include "MergingFactory.fh"
#include "Node.fh"



#include "ThePEG/Handlers/HandlerBase.h"
#include "Herwig/Shower/Dipole/DipoleShowerHandler.h"
  //#include "Herwig/Shower/Dipole/Base/DipoleSplittingGenerator.h"
#include "Herwig/MatrixElement/Matchbox/Base/MergerBase.h"

#include "Herwig/MatrixElement/Matchbox/Phasespace/FFLightTildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/IFLightTildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/FFMassiveTildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/IFMassiveTildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/FILightTildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/IILightTildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/FIMassiveTildeKinematics.h"

#include "ThePEG/Cuts/JetFinder.h"
#include "ThePEG/Cuts/Cuts.h"



namespace Herwig {
  
  using namespace ThePEG;
  
  
  class Merger;
  
  ThePEG_DECLARE_POINTERS(Merger , MergerPtr );
  
  typedef vector<NodePtr> NodePtrVec;
    //definition of a history step
  struct HistoryStep {
      /// containing the full information
    NodePtr node;
      /// current sudakov weight of the history
    double weight;
      /// current scale of the history
    Energy scale;
  };
  
  typedef vector< HistoryStep > History;
  
  typedef multimap<DipoleIndex, Ptr<DipoleSplittingGenerator>::ptr> GeneratorMap2;
  
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
    
    friend class MergingFactory;
    friend class Node;
    
  public:
    
      // define the ME region for a particle vector.
    bool matrixElementRegion(PVector incoming, 
                             PVector outgoing, 
                             Energy winnerScale = ZERO, 
                             Energy cutscale = ZERO)const;
      /// return the current merging scale, 
      /// gets smeared around the central merging scale in generate kinematics.
    Energy mergingScale()const{return theMergePt;}
      /// return the current merging pt (should be unified with mergingScale)
    Energy mergePt()const {return theMergePt;}
      /// legsize of highest process with NLO corrections
    int M()const;
      /// legsize of the highest LO merged process
    int N()const;
      /// legsize of the production process
    int N0()const{return theN0;}
      /// cross section of as given by the merging
    CrossSection MergingDSigDR();
      /// ***** virtual functions of the base class ****///
      /// set the current xcomb, called from ME
    void setXComb( tStdXCombPtr );
      /// set kinematics, called from ME
    void setKinematics();
      ///  clear kinematics, called from ME
    void clearKinematics();
      /// generate kinematics, called from ME
    bool generateKinematics( const double * );
      /// generate kinematics, called from ME
    void flushCaches();
      /// return the current maximum legs, the shower should veto
    size_t maxLegs() const {return theCurrentMaxLegs;}
      /// set the current ME
    void setME(MatchboxMEBasePtr me){
      theCurrentME=me;
      assert(theFirstNodeMap.count(theCurrentME));
      theCurrentNode=theFirstNodeMap[theCurrentME];
    }
      /// allow emissions with a given probability.in the ME region
    double emissionProbability() const{ return theEmissionProbability; }
      /// set a probability of allowed emission into ME region.
    void setEmissionProbability(double x){theEmissionProbability=x;}


  protected:
      /// the merging factory needs to set the legsize of the production process
    void N0(int n){ theN0=n;}
      /// return the large-N basis (TODO: implement check if born ME works with the choice)
    Ptr<ColourBasis>::ptr largeNBasis()const{return theLargeNBasis;}
      /// smear the merging pt
    void smearMergePt(){
    	const double factor = 1. + (-1. + 2.*UseRandom::rnd() ) * smear();
    	theMergePt = factor * centralMergePt();
    }
      /// true if the phase space for initial emissions should not be restricted in z.
    int openZBoundaries()const{return DSH()->showerPhaseSpaceOption();}
      /// return the current ME
    MatchboxMEBasePtr currentME() const { return theCurrentME; }
      /// return the current Node
    NodePtr currentNode() const { return theCurrentNode; }
      /// the gamma parameter to subtract dipoles above a alpha parameter
      /// and subtract the corresponding IPK operator
    double gamma()const{return theGamma;}
    
  private:
      /// calculate a single sudakov step for a given dipole
    double singlesudakov(Dipole, Energy, Energy, pair<bool, bool>);
      /// calculate the sudakov supression for a clusternode between
      /// the current running scale and next scale
    bool   dosudakov(NodePtr Born, Energy  running, Energy next, double& sudakov0_n);
      /// cleanup
    void   cleanup(NodePtr);
      /// return true if the cluster node has the matching number of
      /// legs to the current projector stage
    bool   isProjectorStage( NodePtr , int )const;
      /** 
       * Calculate the staring scale:
       * if Node is part of the production process, calculate according to the
       * scale choice object in the merging scale objekt, else
       * return max(scale as scalechoice , min(Mass(i, j)))
       */
    Energy CKKW_StartScale(NodePtr) const;
      /// prepare the sudakov calculation
    void   CKKW_PrepareSudakov(NodePtr, Energy);
      /// number of active flavours as given by the shower
    double Nf(Energy scale)const{return DSH()->Nf(scale);}
      /// pointer to the factory
    MergingFactoryPtr treefactory() const;
      /// map from ME to first clusternode
    map<MatchboxMEBasePtr, NodePtr> firstNodeMap() const ;
      /// set the current merging pt, smeared in generate kinematics
    void mergePt(Energy x) {theMergePt = x;}
      /// return the central merging pt
    Energy centralMergePt() const {return theCentralMergePt;}
    
  private:
      /// calculate the history weighted born cross section
    CrossSection MergingDSigDRBornStandard();
      /**
       * calculate the history weighted born cross section
       * add the difference of IPK with and without alpha parameter
       * subtract the dipoles above the alpha parameter
       */ 
    CrossSection MergingDSigDRBornGamma();
      /// calculate the history weighted virtual contribution
    CrossSection MergingDSigDRVirtualStandard();
      /**
       * calculate the history weighted real contribution
       * splitted into 3 differnt contibutions
       */
    CrossSection MergingDSigDRRealStandard();
      /// calculate the history weighted real contribution
      /// all dipoles above:
      /// N*(R rnd(i)-Dip_i) history_i U(\phi^n_i)
    CrossSection MergingDSigDRRealAllAbove();
      /// calculate the history weighted real contribution
      /// not all dipoles above:
      /// (R - sum PS_i) history_rnd U(\phi^n+1)
    CrossSection MergingDSigDRRealBelowSubReal();
      /// calculate the history weighted real contribution
      /// not all dipoles above:
      /// rnd(i)-> N*(PS_i - Dip_i) history_i U(\phi^n_i)
    CrossSection MergingDSigDRRealBelowSubInt();
      /// max legssize the shower should veto for LO
    size_t maxLegsLO() const {return N0()+N();}
      /// Calculate the LO partonic cross section.
      /// if diffalpha != 1, add the difference of IPK(1)-IPK(diffalpha)
    CrossSection TreedSigDR(Energy startscale, double diffalpha=1.);
      /// fill the projecting xcomb
    Energy fillProjector(int);
      /// fill the history, including calculation of sudakov supression
    void   fillHistory(Energy, NodePtr, NodePtr );
      /// calculate the pdf ratio for the given clusternode
    double pdfratio(NodePtr, Energy, Energy, int, bool fromIsME, bool toIsME);
      /// return the pdf-ratio reweight for the history
    double pdfReweight();
      /// return the alpha_s reweight for the history
    double alphaReweight(bool nocmw=false);
      /// max legssize the shower should veto for NLO
    size_t maxLegsNLO()const {return N0()+M();}
      /// calculate the virtual contribution.
    CrossSection LoopdSigDR(Energy startscale );
      /// calculate alpha_s expansion of the pdf-ratios
    double sumPdfReweightExpansion()const;
      /// calculate alpha_s expansion of the alpha_s-ratios, including K_g
    double sumAlphaSReweightExpansion()const;
      /// calculate alpha_s expansion of the sudakov exponents
    double sumFillHistoryExpansion();
      /// calculate alpha_s expansion of the single step alpha_s-ratio, including K_g
    double alphasExpansion( Energy next, Energy fixedScale)const;
      /// calculate alpha_s expansion of the single step pdf-ratio
    double pdfExpansion(NodePtr, int, Energy, Energy, double, int, Energy)const;
      /// calculate alpha_s expansion of the single step sudakov exponent
    bool   doHistExpansion(NodePtr Born, Energy  running, Energy next, Energy fixedScale, double& HistExpansion);
      /// calculate alpha_s expansion of the single dipole sudakov exponent
    double singleHistExpansion(Dipole, Energy, Energy, Energy, pair<bool, bool>);
      //alpha_s as given in the shower
    double as(Energy q)const{return DSH()->as(q);}
      // set the pointer to the Mergingfactory.
    void setFactory(MergingFactoryPtr f){theTreeFactory=f;}
      // set the pointer to the DipoleShower.
    void setDipoleShower(DipoleShowerHandlerPtr dsh){theDipoleShowerHandler=dsh;}
      //return the dipole shower handler
    DipoleShowerHandlerPtr DSH(){return theDipoleShowerHandler;}
      //return the const dipole shower handler
    cDipoleShowerHandlerPtr DSH()const{return theDipoleShowerHandler;}
      /// insert map from ME to first clusternode
    void firstNodeMap(MatchboxMEBasePtr, NodePtr);
      /// history choice: weighted history choice
    int chooseHistory()const {return theChooseHistory;}
      /// the smearing factor for the merging scale
    double smear()const{return theSmearing;}
      /// return the large-N colour basis
    void largeNBasis(Ptr<ColourBasis>::ptr x){theLargeNBasis=x;}
      /// helper function to check the only multi condition.
    bool notOnlyMulti()const;
      /// Calculate the CMW AlphaS
    double cmwAlphaS(Energy q)const;
      /// debug output for virtual
    void debugVirt(double, double, double, double, CrossSection,
                   double, double, double, NodePtr,CrossSection) const;
      /// debug output for reals
    void debugReal( string, double, CrossSection, CrossSection) const;

    
  private:
    
      /// calculate the history expansion
    unsigned int theShowerExpansionWeights = 2;
      /// use CMW scheme
    unsigned int theCMWScheme = 0;
      /// true if current point should be projected
    bool projected = true;
      /// true if LO cross sections should be unitarised
    bool isUnitarized = true;
      /// true if NLO contributions should be unitarised
    bool isNLOUnitarized = true;
      /// history weight choice
    int theChooseHistory = 0;
      /// legsize of production process
    int theN0 = 0;
      /// calculate only the N particle contribution
    int theOnlyN = -1;
      /// the current maxlegs (either LO or NLO maxlegs)
    int theCurrentMaxLegs = -1;
      /// current weight and weight of clustered born
    double weight = 1.0;
    double weightCB = 1.0;
      /// subtract the dipole contribution above a given gamma
    double theGamma = 1.0;
      /// smearing factor for merging scale
    double theSmearing = 0.;

     /** 
      * Allow emissions for the unitarising LO contributions 
      *  according to the prob 1-min(B_n/sum Dip_n,Dip_n/B_n)
      */
    bool emitDipoleMEDiff = false;
      /// The conditional emission probability if emitDipoleMEDiff is true.
    double theEmissionProbability = 0.;
      /// cutoff for real emission contribution
    Energy theIRSafePT = 1_GeV;
      /// current merging scale
    Energy theMergePt = 4_GeV;
      /// central merging scale
    Energy theCentralMergePt = 4_GeV;
      /// below mergingscale/theRealSubtractionRatio the dipoles are used to subtract.
      /// above the shower approximation is in use.
    double theRealSubtractionRatio=3.;
      /// current cluster histoy including sudakov weights
    History history;
      /// pointer to the large-N basis
    Ptr<ColourBasis>::ptr theLargeNBasis;
      /// current Node
    NodePtr theCurrentNode;
      /// current ME
    MatchboxMEBasePtr theCurrentME;
      /// Tilde kinematics pointers, only to use lastPt(emitter, emission, spectator)
    Ptr<FFLightTildeKinematics>::ptr FFLTK = new_ptr( FFLightTildeKinematics() );
    Ptr<FILightTildeKinematics>::ptr FILTK = new_ptr( FILightTildeKinematics() );
    Ptr<IFLightTildeKinematics>::ptr IFLTK = new_ptr( IFLightTildeKinematics() );
    Ptr<IILightTildeKinematics>::ptr IILTK = new_ptr( IILightTildeKinematics() );
    Ptr<FFMassiveTildeKinematics>::ptr FFMTK = new_ptr( FFMassiveTildeKinematics() );
    Ptr<FIMassiveTildeKinematics>::ptr FIMTK = new_ptr( FIMassiveTildeKinematics() );
    Ptr<IFMassiveTildeKinematics>::ptr IFMTK = new_ptr( IFMassiveTildeKinematics() );
      //pointer to the shower handler
    DipoleShowerHandlerPtr theDipoleShowerHandler;
      /// pointer to the MergingFactory
    MergingFactoryPtr theTreeFactory;
      /// map from ME to first Node
    map<MatchboxMEBasePtr, NodePtr> theFirstNodeMap;
      /// map from ME to highest ME weight so far
    map<NodePtr, CrossSection> theHighMeWeightMap;
    
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
    Merger & operator=(const Merger &) = delete;
    
  };
  
}



#endif /* HERWIG_Merger_H */
