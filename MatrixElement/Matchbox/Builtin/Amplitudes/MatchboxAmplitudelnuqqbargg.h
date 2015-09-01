// -*- C++ -*-
//
// MatchboxAmplitudelnuqqbargg.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_MatchboxAmplitudelnuqqbargg_H
#define Herwig_MatchboxAmplitudelnuqqbargg_H
//
// This is the declaration of the MatchboxAmplitudelnuqqbargg class.
//

#include "Herwig/MatrixElement/Matchbox/Base/MatchboxAmplitude.h"
#include "Herwig/MatrixElement/Matchbox/Builtin/Amplitudes/MatchboxCurrents.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Herwig/MatrixElement/Matchbox/Utility/SU2Helper.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Thomas Schuh
 *
 * \brief MatchboxAmplitudelnuqqbargg
 */
class MatchboxAmplitudelnuqqbargg: public MatchboxAmplitude, public MatchboxCurrents {

public:

  /**
   * The default constructor.
   */
  MatchboxAmplitudelnuqqbargg();

  /**
   * The destructor.
   */
  virtual ~MatchboxAmplitudelnuqqbargg();

  /**
   * Return the (tree-level) order in \f$g_S\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInGs() const { return 2; }

  /**
   * Return the (tree-level) order in \f$g_{EM}\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInGem() const { return 2; }

  /**
   * Return true, if this amplitude is capable of calculating one-loop
   * (QCD) corrections.
   */
  virtual bool haveOneLoop() const { return false; }
  
  /**
   * Return true, if this amplitude can handle the given process.
   */
  virtual bool canHandle(const PDVector&) const;
  
  /**
   * Calculate the tree level amplitudes for the phasespace point
   * stored in lastXComb.
   */
  virtual void prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr);
  
  /**
   * Evaluate the amplitude for the given colour tensor id and
   * helicity assignment
   */
  virtual Complex evaluate(size_t, const vector<int>&, Complex&);

  /**
   * Return true, if one loop corrections are given in the conventions
   * of BDK.
   */
  virtual bool isBDK() const { return true; }

  /**
   * Return the value of the dimensional regularization
   * parameter. Note that renormalization scale dependence is fully
   * restored in DipoleIOperator.
   */
  virtual Energy2 mu2() const { return lastSHat(); }

  /**
   * Flush all cashes.
   */
  virtual void flushCaches() {
    MatchboxCurrents::reset();
    MatchboxAmplitude::flushCaches();
  }

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
  
  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();
  
protected:
    
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

private:
  
  /**
   * True, if a diagonal CKM matrix should be assumed. This ignores
   * the CKM object of the StandardModel.
   */
  bool theDiagonal;
  
  /**
   * The ckm.
   */
  vector< vector<Complex> > theCKM;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MatchboxAmplitudelnuqqbargg & operator=(const MatchboxAmplitudelnuqqbargg &);

};

}

#endif /* Herwig_MatchboxAmplitudelnuqqbargg_H */
