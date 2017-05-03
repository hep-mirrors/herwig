  /// -*- C++ -*-
  //
  /// MergingFactory.h is a part of Herwig - A multi-purpose Monte Carlo event generator
  /// Copyright (C) 2002-2017 The Herwig Collaboration
  //
  /// Herwig is licenced under version 3 of the GPL, see COPYING for details.
  /// Please respect the MCnet academic guidelines, see GUIDELINES for details.
  //
#ifndef HERWIG_MergingFactory_H
#define HERWIG_MergingFactory_H
  //
  /// This is the declaration of the MergingFactory class.
  //

#include "MergingFactory.fh"
#include "Herwig/MatrixElement/Matchbox/MatchboxFactory.h"
#include "Node.fh"
#include "Merger.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
namespace Herwig {
  
  using namespace ThePEG;
  
  /**
   * \ingroup Matchbox
   * \author Johannes Bellm
   *
   * \brief MergingFactory automatically sets up a NLO
   * QCD merging.
   *
   * @see \ref MergingFactoryInterfaces "The interfaces"
   * defined for MergeboxFactory.
   */
  class MergingFactory : public MatchboxFactory {
  public:
      ///Check consistency and switch to porduction mode.
    void productionMode();
      /// main method to setup the ME vector
    virtual void setup();
      /// fill all amplitudes, stored in pureMEsMap
    void fillMEsMap();
      /// prepare the Born and virtual matrix elements.
    void prepare_BV(int i);
      /// prepare the real emission matrix elements.
    void prepare_R(int i);
      /// push the born contributions to the ME vector.
    void pushB(MatchboxMEBasePtr, int);
      //push the virtual contributions to the ME vector.
    void pushV(MatchboxMEBasePtr, int);
      /// push the real contributions to the ME vector.
    void pushR(MatchboxMEBasePtr, int);
      /// order matrix elements from one loop provider.
    void orderOLPs();
      /// Debugging: push only multiplicities to the ME vector
      /// in range of specified mulltiplicity.
    int onlymulti()const {
       return theonlymulti==-1?-1:(theonlymulti+processMap.find(0)->second.size());
    }
      /// pointer to the merging helper.
    MergerPtr MH() {return theMergingHelper;}
      /// maximal NLO mulitplicity: 0=NLO corrections to the productio process.
    int M() const {return theM-1;}
      /// leg size of highest multiplicity.
    int N() const {return theN;}
     /// Return the Map of matrix elements to be considered
     /// (the Key is the number of additional jets)
    const map<int, vector<MatchboxMEBasePtr> >& pureMEsMap() const {
      return thePureMEsMap;
    }
      /// Access the Map of matrix elements to be considered
      /// (the Key is the number of additional jets)
    map<int, vector<MatchboxMEBasePtr> >& pureMEsMap() {
      return thePureMEsMap;
    }
      //Parse a process description
    virtual vector<string> parseProcess(string);
      // fill the virtuals vector (these are IPK-operators)
    void getVirtuals(MatchboxMEBasePtr nlo, bool clone );
      // In the merged setup we only produce single phase space points. 
    bool subProcessGroups() const { return false;}
      // Cut on non-QCD observables. 
    Ptr<Cuts>::ptr nonQCDCuts(){return theNonQCDCuts;}
 
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
    void persistentInput(PersistentIStream & is, int);
      //@}
    
    static void Init();
    
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
    
  protected:
    
    /** @name Clone Methods. */
      //@{
    /**
     * Make a simple clone of this object.
     * @return a pointer to the new object.
     */
    virtual IBPtr clone() const;
    
    /** 
     * Make a clone of this object, possibly modifying the cloned object
     * to make it sane.
     * @return a pointer to the new object.
     */
    virtual IBPtr fullclone() const;
      //@}
    
  private:
    
      /// unitarise the LO contributions.
    bool unitarized = true;
      /// did run setup.
    bool ransetup = false;
      /// Debugging: push only multiplicities to the ME vector
      /// in range of specified mulltiplicity.
    int theonlymulti = -1;
      /// maximal legsize for NLO corrections.
    int theM = -1;
      /// maximal legsize for LO contributions.
    int theN = -1;
      /// map for processes.
    map< int, vector<string> > processMap;
      //The matrix elements: int = number of additional jets
    map< int, vector<MatchboxMEBasePtr> > thePureMEsMap;
      /// the merging helper
    MergerPtr theMergingHelper;
      /// Cut on non-QCD modified observables. 
    Ptr<Cuts>::ptr theNonQCDCuts;  
    /**
     * The assignment operator is private and must never be called.
     * In fact, it should not even be implemented.
     */
    MergingFactory & operator=(const MergingFactory &) = delete;

    
  };
  
}

#endif /* HERWIG_MergingFactory_H */
