// -*- C++ -*-
//
// MadGraphAmplitude.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_MadGraphAmplitude_H
#define Herwig_MadGraphAmplitude_H
//
// This is the declaration of the MadGraphAmplitude class.
//

#include "Herwig/MatrixElement/Matchbox/Base/MatchboxAmplitude.h"
#include "Herwig/MatrixElement/Matchbox/Builtin/Amplitudes/MatchboxCurrents.h"
#include "ThePEG/Utilities/DynamicLoader.h"

namespace Herwig {


using namespace ThePEG; 

/**
 * \ingroup Matchbox
 * \author Johannes Bellm, Simon Platzer
 *
 * \brief MadGraphAmplitude implements an interface to MadGraph
 */
class MadGraphAmplitude: 
    public MatchboxAmplitude {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  MadGraphAmplitude();

  /**
   * The destructor.
   */
  virtual ~MadGraphAmplitude();
  //@}

public:

  /**
   * Return true, if this amplitude can handle the given process.
   */
  virtual bool canHandle(const PDVector&,
			 Ptr<MatchboxFactory>::tptr,
			 bool) const;

  /**
   * Return true, if this amplitude already includes symmetry factors
   * for identical outgoing particles.
   */
  virtual bool hasFinalStateSymmetry() const { return false; }

  /**
   * Set the (tree-level) order in \f$g_S\f$ in which this matrix
   * element should be evaluated.
   */
  virtual void orderInGs(unsigned int ogs) { theOrderInGs = ogs; }

  /**
   * Return the (tree-level) order in \f$g_S\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInGs() const { return theOrderInGs; }

  /**
   * Set the (tree-level) order in \f$g_{EM}\f$ in which this matrix
   * element should be evaluated.
   */
  virtual void orderInGem(unsigned int oge) { theOrderInGem = oge; }

  /**
   * Return the (tree-level) order in \f$g_{EM}\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInGem() const { return theOrderInGem; }

  /**
   * Return true, if this amplitude is capable of calculating one-loop
   * (QCD) corrections.
   */
  virtual bool haveOneLoop() const { return true; }

  /**
   * Calculate the tree level amplitudes for the phasespace point
   * stored in lastXComb.
   */
  virtual void prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr);
  
  /**
   * Calculate the one-loop amplitudes for the phasespace point
   * stored in lastXComb, if provided.
   */
  virtual void prepareOneLoopAmplitudes(Ptr<MatchboxMEBase>::tcptr);

  /**
   * Evaluate the amplitude for the given colour tensor id and
   * helicity assignment
   */
  virtual Complex evaluate(size_t, const vector<int>&, Complex&);

   /**
   * Return true, if one-loop contributions will be evaluated at amplitude level.
   */
  virtual bool oneLoopAmplitudes() const { return false; }
  
   /**
   * Return the one-loop/tree interference.
   */
  virtual double oneLoopInterference() const;
  
  void evaloneLoopInterference() const;


  /**
   * Return true, if one loop corrections are given in the conventions
   * of BDK.
   */

  virtual bool isCS() const { return false; }
  virtual bool isExpanded() const { return true; }
  virtual bool isBDK() const { return false; }
  virtual bool isDR() const { return false; }
  virtual bool isDRbar() const {return false;}

  /**
   * Return the value of the dimensional regularization
   * parameter. Note that renormalization scale dependence is fully
   * restored in DipoleIOperator.
   */
  virtual Energy2 mu2() const { return sqr(91.18800000000*GeV);}
  //virtual Energy2 mu2() const { return lastSHat(); }
  
  
  virtual LorentzVector<Complex> plusPolarization(const Lorentz5Momentum& p,
						  const Lorentz5Momentum& n,
						  int id = -1) const;
  /**
   * Return the colour and spin correlated matrix element.
   */
  virtual double spinColourCorrelatedME2(pair<int,int> emitterSpectator,
					 const SpinCorrelationTensor& c) const;

  /**
   * Order process in MadGraph conventions
   */						  
  vector<size_t> gluonsFirst(vector<size_t> i);

  /**
   * Flush all cashes.
   */
  virtual void flushCaches() {
    MatchboxAmplitude::flushCaches();
  }

  /**
   * Return true, if this amplitude needs to initialize an external
   * code.
   */
  virtual bool isExternal() const { return true; }

  /**
   * Initialize this amplitude
   */
  virtual bool initializeExternal();

  /**
   * Return a generic process id for the given process
   */
  virtual int externalId(const cPDVector&);
      
  
  /**
   * Write out BornAmplitudes.dat and VirtualAmplitudes.dat,
   * return true if new files are written.
   */
      
  virtual bool writeAmplitudesDat();
      
  /**
   * CHeck if all amplitudes are in the Libary.
   */
      
  virtual bool checkAmplitudes();
      
  
  string mgProcLibPath();

  /**
   * Return true, if this amplitude is capable of consistently filling
   * the rho matrices for the spin correllations
   */
  virtual bool canFillRhoMatrix() const { return true; }

  /**
   * Return the helicity combination of the physical process in the
   * conventions used by the spin correlation algorithm.
   */
  virtual vector<unsigned int> physicalHelicities(const vector<int>&) const;

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

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();

  /**
   * The (tree-level) order in \f$g_S\f$ in which this matrix
   * element is given.
   */
  unsigned int theOrderInGs;

  /**
   * The (tree-level) order in \f$g_{EM}\f$ in which this matrix
   * element is given.
   */
  unsigned int theOrderInGem;
  
  /**
   * The path for the process libraries.
   */
  static string theProcessPath;
  
  /**
   * The path to generate amplitudes in.
   */
  static string theMGmodel;
  
  bool keepinputtopmass;


  /**
   * Initialize the given process
   */
  void initProcess(const cPDVector&);


  /**
   * Storage for Amplitudes 
   */
  static vector<string>  BornAmplitudes,VirtAmplitudes;
   
    /**
   * Helper for color and crossing handling 
   */ 
  mutable vector<int>  colourindex, crossing;
  
        
  /**
  * Static Variables to handle initialization.
  */
        
  static bool ranMadGraphInitializeExternal;
  static bool initializedMad;
  

  
  //@}

protected:

  /**
   *   Location of the installed executables
   */
  static string bindir_;

  /**
   *   Location of the installed include files
   */
  static string includedir_;

  /**
   *   Location of the data files
   */
  static string pkgdatadir_;

  /**
   *  Location of MADGRAPH
   */
  static string madgraphPrefix_;
      
  /**
   *  Helper functions to make long strings static
   */
      
  void setProcessPath(string );
  string getProcessPath() const;
      
  void setBinDir(string p);
  string getBinDir() const;
      
  void setDataDir(string p);
  string getDataDir() const;
     
  void setModel(string p);
  string getModel() const;
      
  void setMadgraphPrefix(string p);
  string getMadgraphPrefix() const ;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MadGraphAmplitude & operator=(const MadGraphAmplitude &) = delete;


};

}

#endif /* Herwig_MadGraphAmplitude_H */
