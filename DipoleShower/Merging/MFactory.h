  // -*- C++ -*-
  //
  // MFactory.h is a part of Herwig - A multi-purpose Monte Carlo event generator
  // Copyright (C) 2002-2012 The Herwig Collaboration
  //
  // Herwig is licenced under version 2 of the GPL, see COPYING for details.
  // Please respect the MCnet academic guidelines, see GUIDELINES for details.
  //
#ifndef HERWIG_MFactory_H
#define HERWIG_MFactory_H
  //
  // This is the declaration of the MergeboxFactory class.
  //

#include "Herwig/MatrixElement/Matchbox/MatchboxFactory.h"
#include "Node.fh"
#include "Merger.h"

  //#include "ThePEG/MatrixElement/ReweightConstant.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
namespace Herwig {
  
  using namespace ThePEG;
  
  /**
   * \ingroup Matchbox
   * \author Johannes Bellm
   *
   * \brief MFactory automatically sets up a NLO
   * QCD calculation carried out in dipole subtraction.
   *
   * @see \ref MFactoryInterfaces "The interfaces"
   * defined for MergeboxFactory.
   */
  class MFactory : public MatchboxFactory {
  public:
    
    /** @name Standard constructors and destructors. */
      //@{
    /**
     * The default constructor.
     */
    MFactory();
    
    /**
     * The destructor.
     */
    virtual ~MFactory();
      //@}
    
    virtual void setup();
    
    /**fill all amplitudes, stored in pureMEsMap
     */
    void fill_amplitudes();
    
    /** prepare the Born and virtual matrix elements.
     */
    void prepare_BV(int i);
    
    void prepare_R(int i);
    
    void pushB(Ptr<MatchboxMEBase>::ptr,int);
    
    void pushV(Ptr<MatchboxMEBase>::ptr,int);
    
    void pushProR(Ptr<MatchboxMEBase>::ptr,int);
    
    void orderOLPs();
    
    int onlymulti()const {return theonlymulti==-1?-1:(theonlymulti+processMap.find(0)->second.size());}
    
    bool onlyUnlopsweights() const {return theonlyUnlopsweights;}
    
    Ptr<Merger>::ptr MH(){return theMergingHelper;}
    
    
    /**
     * Return the Map of matrix elements to be considered
     * (the Key is the number of additional jets)
     */
    const map<int, vector<Ptr<MatchboxMEBase>::ptr> >& pureMEsMap() const {
      return thePureMEsMap;
    }
    
    /**
     * Access the Map of matrix elements to be considered
     * (the Key is the number of additional jets)
     */
    map<int, vector<Ptr<MatchboxMEBase>::ptr> >& pureMEsMap() {
      return thePureMEsMap;
    }
    
    
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
    
    bool theonlyNLOParts;
    bool theonlyvirtualNLOParts;
    bool theonlyrealNLOParts;
    bool theonlyUnlopsweights;
    bool theunitarizeNLOParts;
    bool calc_born;
    bool calc_virtual;
    bool calc_real;
    bool unitarized;
    bool NLOunitarized;
    int theonlyk;
    int theonlymulti;
    int theonlyabove;
    int theonlysub;
    int divideSub;
    int divideSubNumber;
    
    
    /**
     * Prefix for subtraction data
     */
    string theSubtractionData;
    
    
    map< int, vector<string> > processMap;
    
    map< int, bool > thevirtualsAreExpandedMap;
    
    
    /**
     * The matrix elements: int = number of additional jets
     */
    
    map< int, vector<Ptr<MatchboxMEBase>::ptr> > thePureMEsMap;
    
    /**
     * The Born matrix elements to be considered
     */
    vector<Ptr<MatchboxMEBase>::ptr> theBornMEs;
    
    /**
     * The virtual corrections to be considered
     */
    map<int, vector<Ptr<MatchboxInsertionOperator>::ptr> > theVirtualsMap;
    
    /**
     * The produced NLO matrix elements
     */
    vector<Ptr<MatchboxMEBase>::ptr> theBornVirtualMEs;
    
    /**
     * The real emission matrix elements to be considered
     */
    
    vector<Ptr<MatchboxMEBase>::ptr> theRealEmissionMEs;
    
    /**
     * A large-N colour basis to be used when reproducing the shower
     * kernels.
     */
    
    Ptr<Merger>::ptr theMergingHelper;
    
    bool ransetup;
    
    /**
     * The assignment operator is private and must never be called.
     * In fact, it should not even be implemented.
     */
    MFactory & operator=(const MFactory &);
    
  };
  
}

#endif /* HERWIG_MFactory_H */
