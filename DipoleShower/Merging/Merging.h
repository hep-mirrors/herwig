  // -*- C++ -*-
  //
  // Merging.h is a part of Herwig - A multi-purpose Monte Carlo event generator
  // Copyright (C) 2002-2007 The Herwig Collaboration
  //
  // Herwig is licenced under version 2 of the GPL, see COPYING for details.
  // Please respect the MCnet academic guidelines, see GUIDELINES for details.
  //
#ifndef HERWIG_Merging_H
#define HERWIG_Merging_H
  //
  // This is the declaration of the Merging class.
  //

#include "ThePEG/Handlers/HandlerBase.h"
#include "ThePEG/Vectors/SpinOneLorentzRotation.h"
#include "Herwig/DipoleShower/DipoleShowerHandler.h"
#include "Herwig/DipoleShower/Base/DipoleSplittingGenerator.h"
#include "Herwig/MatrixElement/Matchbox/Mergeing/ClusterNode.h"


namespace Herwig {
  
  using namespace ThePEG;
  
  typedef Ptr<ClusterNode>::ptr CNPtr;
  
  struct HistStep {
    CNPtr node;
    double weight;
    Energy scale;
  };
  
  typedef vector< HistStep > Hist;
  
  typedef multimap<DipoleIndex,Ptr<DipoleSplittingGenerator>::ptr> GeneratorMap2;
  
  /**
   * \ingroup DipoleShower
   * \author Johannes Bellm
   *
   * \brief Merging handles the merging ....... //TODO .
   *
   * @see \ref MergingInterfaces "The interfaces"
   * defined for Merging.
   */
  class Merging: public HandlerBase {
    
  public:
    
    /** @name Standard constructors and destructors. */
      //@{
    /**
     * The default constructor.
     */
    Merging();
    
    /**
     * The destructor.
     */
    virtual ~Merging();
      //@}
    
  public:
    
    bool   sudakov(CNPtr Born, Energy & running, Energy next);
    double singlesudakov(list<Dipole>::iterator,Energy,Energy,pair<bool,bool>,bool fast=false);
    bool   dosudakov(Energy & running, Energy next, double& sudakov0_n,bool fast=false);
    bool   dosudakovold(CNPtr, Energy&, Energy, double&);
    
    
    
    void   cleanup(CNPtr);
    
    bool   projectorStage(CNPtr);
    Energy CKKW_StartScale(CNPtr);
    void   CKKW_PrepareSudakov(CNPtr,Energy);
    double matrixElementWeight(Energy startscale,CNPtr);
    bool   fillProjector(Energy&);
    void   fillHistory(Energy&, CNPtr, CNPtr ,bool fast=false);
    
    
    double sumpdfReweightUnlops();
    double sumalphaReweightUnlops();
    double pdfratio(CNPtr,Energy&,Energy,int);
    double pdfReweight();
    double alphaReweight();
    
    double sumfillHistoryUnlops();
    
    
    
    double alphasUnlops( Energy next,Energy fixedScale);
    double pdfUnlops(tcPDPtr,tcPDPtr,tcPDFPtr,Energy,Energy,double,int,Energy);
    double singleUNLOPS(list<Dipole>::iterator,Energy,Energy,Energy,pair<bool,bool>);
    bool   doUNLOPS(Energy  running, Energy next,Energy fixedScale, double& UNLOPS);
    
    
    
    
    bool   reweightCKKWSingle(Ptr<MatchboxXComb>::ptr SX, double & res,bool fast=false) ;
    double reweightCKKWBorn(CNPtr Node,bool fast=false);
    double reweightCKKWBorn2(CNPtr Node,bool fast=false);

    
    double reweightCKKWVirt(CNPtr Node);
    double reweightCKKWReal(CNPtr Node);
    double as(Energy q){return theDipoleShowerHandler->as(q);}
    
    
    Energy mergingScale()const {return MergingScale;}
    
    
    unsigned int maxLegsLO() const {return theMaxLegsLO;}
    unsigned int maxLegsNLO()const {return theMaxLegsNLO;}
    
    
    
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
     * Variations
     */
    double   xiRenME;
    double   xiFacME;
    double   xiRenSh;
    double   xiFacSh;
    double   xiQSh;
    
    
    /**
     *
     **/
    
    double theNf;
    
    bool minusL;
    bool Unlopsweights;
    bool theKImproved;
    Energy MergingScale;
    
    unsigned int theMaxLegsLO;
    unsigned int theMaxLegsNLO;
    double projectorWeight;
    double fastweight;
    
    Hist history;
    CNPtr StartingBorn;
    CNPtr StartingCalcBornBorn;
    CNPtr CalcBorn;
    
    
    /**
     * The mean of the Gaussian distribution for
     * the intrinsic pt of sea partons.
     */
    Ptr<DipoleShowerHandler>::ptr theDipoleShowerHandler;
    
  private:
    
    /**
     * The static object used to initialize the description of this class.
     * Indicates that this is a concrete class with persistent data.
     */
    static ClassDescription<Merging> initMerging;
    
    /**
     * The assignment operator is private and must never be called.
     * In fact, it should not even be implemented.
     */
    Merging & operator=(const Merging &);
    
  };
  
}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {
  
  /** @cond TRAITSPECIALIZATIONS */
  
  /** This template specialization informs ThePEG about the
   *  base classes of Merging. */
  template <>
  struct BaseClassTrait<Herwig::Merging,1> {
    /** Typedef of the first base class of Merging. */
    typedef HandlerBase NthBase;
  };
  
  /** This template specialization informs ThePEG about the name of
   *  the Merging class and the shared object where it is defined. */
  template <>
  struct ClassTraits<Herwig::Merging>
  : public ClassTraitsBase<Herwig::Merging> {
    /** Return a platform-independent class name */
    static string className() { return "Herwig::Merging"; }
    /**
     * The name of a file containing the dynamic library where the class
     * Merging is implemented. It may also include several, space-separated,
     * libraries if the class Merging depends on other classes (base classes
     * excepted). In this case the listed libraries will be dynamically
     * linked in the order they are specified.
     */
    static string library() { return "HwDipoleShower.so"; }
  };
  
  /** @endcond */
  
}

#endif /* HERWIG_Merging_H */
