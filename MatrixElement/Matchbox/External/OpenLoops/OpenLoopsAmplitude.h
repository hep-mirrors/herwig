// -*- C++ -*-
//
// OpenLoopsAmplitude.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_OpenLoopsAmplitude_H
#define Herwig_OpenLoopsAmplitude_H
//
// This is the declaration of the OpenLoopsAmplitude class.
//

#include "Herwig/MatrixElement/Matchbox/Base/MatchboxOLPME.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Johannes Bellm, Simon Platzer
 *
 * \brief Process information for OpenLoops
 */
class OpenLoopsProcInfo{

public:

  /**
   * Default constructor
   */
  OpenLoopsProcInfo() {}

  /**
   * Construct giving data
   */
  OpenLoopsProcInfo(int HID,int GID, string procstr,string typestr)
    : theHOlpId(HID), theGOlpId(GID), theProcstr(procstr), theTypestr(typestr) {}

  /**
   * Document me
   */
  int HID() const { return theHOlpId; }

  /**
   * Document me
   */
  int GID() const { return theGOlpId; }

  /**
   * Document me
   */
  const string& Pstr() const { return theProcstr; }

  /**
   * Document me
   */
  const string& Tstr() const { return theTypestr; }

  /**
   * Document me
   */
  void setGID(int g) { theGOlpId=g; }

  /**
   * Document me
   */
  void setOAs(int i) { orderAlphas=i; }

  /**
   * Document me
   */
  int orderAs() { return orderAlphas; }

private:

  /**
   * Document me
   */
  int theHOlpId;

  /**
   * Document me
   */
  int theGOlpId;

  /**
   * Document me
   */
  string theProcstr;

  /**
   * Document me
   */
  string theTypestr;

  /**
   * Document me
   */
  int orderAlphas;

public:

  /**
   * Write to persistent stream
   */
  void persistentOutput(PersistentOStream & os) const{
    os << theHOlpId << theGOlpId << theProcstr << theTypestr << orderAlphas;
  }

  /**
   * Read from persistent stream
   */
  void persistentInput(PersistentIStream &is) {
    is >> theHOlpId >> theGOlpId >> theProcstr >> theTypestr >> orderAlphas;
  }

};

/**
 * \ingroup Matchbox
 * \author Johannes Bellm, Simon Platzer
 *
 * \brief OpenLoopsAmplitude implements an interface to OpenLoops
 */
class OpenLoopsAmplitude: public MatchboxOLPME {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  OpenLoopsAmplitude();

  /**
   * The destructor.
   */
  virtual ~OpenLoopsAmplitude();
  //@}

public:

  virtual void fillOrderFile(const map<pair<Process,int>,int>& procs);


  virtual bool isCS() const { return false; }
  virtual bool isExpanded() const { return true; }
  virtual bool isBDK() const { return false; }
  //virtual bool isDR() const { return true; }
  /**
   * Start the one loop provider, if appropriate, giving order and
   * contract files
   */

  virtual bool checkOLPContract();

  /**
   * Start the one loop provider, if appropriate
   */
  virtual void startOLP(const string&, int& status);

  /**
   * Return the value of the dimensional regularization
   * parameter. Note that renormalization scale dependence is fully
   * restored in DipoleIOperator.
   */
  // virtual Energy2 mu2() const { return lastSHat(); }

  /**
   * Start the one loop provider, if appropriate. This default
   * implementation writes an BLHA 2.0 order file and starts the OLP
   */
  virtual bool startOLP(const map<pair<Process,int>,int>& procs);


    /**
   * Return true, if this amplitude already includes averaging over
   * incoming parton's quantum numbers.
   */
  virtual bool hasInitialAverage() const { return true; }

  /**
   * Return true, if this amplitude already includes symmetry factors
   * for identical outgoing particles.
   */
  virtual bool hasFinalStateSymmetry() const { return true; }

  /**
   * Call OLP_EvalSubProcess and fill in the results
   */
   void evalSubProcess() const;

  /**
   * Fill in results for the given colour correlator
   */
  virtual void evalColourCorrelator(pair<int,int> ij) const;


  /**
   * Fill in results for the given colour/spin correlator
   */
  virtual void evalSpinColourCorrelator(pair<int,int> ij) const;


  /**
   * Return the colour and spin correlated matrix element.
   */
  virtual double spinColourCorrelatedME2(pair<int,int> ij,
					 const SpinCorrelationTensor& c) const;


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

  virtual void doinitrun();


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  OpenLoopsAmplitude & operator=(const OpenLoopsAmplitude &);

  /**
   * Store colour correlator results
   */
  mutable vector<double> colourCorrelatorResults;

  /**
   * Store spin colour correlator results
   */
  mutable vector<double> spinColourCorrelatorResults;


  /**
   * first is the olp id from herwig, second the answer from openloops
   */
  static vector< int > idpair;

  
  /**
   * Helper map to store information in different procs.
   */
  
  map<int , OpenLoopsProcInfo > processmap;
  
  
  /**
   * Interface for Higgs Effective
   */
  bool theHiggsEff;
  
  /**
   * Complex Mass Scheme.
   */
  
  bool use_cms;
  
 
  /**
   * Use of Collier Lib (arXiv:1604.06792), available since OpenLoops 1.3.0.
   */
  bool theCollierLib=true; 
 
  /**
   * parameter to set Phase space tolerance for massiv particles.
   * Should not be used. Better: set Openloops:Massless 11
   */
  
  int psp_tolerance;

  /**
   *   Location of the OpenLoops libraries
   */
  static string OpenLoopsLibs_;

  /**
   *   Location of the OpenLoops
   */
  static string OpenLoopsPrefix_;
  
  
  /**
   *  Helper functions to make long strings static
   */
  
  void setOpenLoopsLibs(string p);
  string getOpenLoopsLibs() const;
  
  void setOpenLoopsPrefix(string p);
  string getOpenLoopsPrefix() const;

  
  
  
};

}

#endif /* Herwig_OpenLoopsAmplitude_H */
