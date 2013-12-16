// -*- C++ -*-
//
// VBFNLOAmplitudeBase.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_VBFNLOAmplitudeBase_H
#define HERWIG_VBFNLOAmplitudeBase_H
//
// This is the declaration of the VBFNLOAmplitudeBase class.
//

#include "Herwig++/MatrixElement/Matchbox/Base/MatchboxAmplitude.h"
#include "Diagrams/DiagramContainer.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Ken Arnold
 *
 * \brief VBFNLOAmplitudeBase is the base class for amplitudes
 * that implement interfaces to VBFNLO.
 *
 * @see \ref VBFNLOAmplitudeBaseInterfaces "The interfaces"
 * defined for VBFNLOAmplitudeBase.
 */
class VBFNLOAmplitudeBase: public MatchboxAmplitude {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  VBFNLOAmplitudeBase();

  /**
   * The destructor.
   */
  virtual ~VBFNLOAmplitudeBase();
  //@}

public:

  /** @name Tree-level amplitudes */
  //@{

  /**
   * Return true, if tree-level contributions will be evaluated at amplitude level.
   */
  virtual bool treeAmplitudes() const { return false; }

  //@}

  /** @name Subprocess information */
  //@{
  
  /**
   * Tell whether the outgoing partons should be sorted when determining
   * allowed subprocesses. Otherwise, all permutations are counted as
   * separate subprocesses.
   */
  virtual bool sortOutgoing() { return false; }

  //@}

  /** @name One-loop amplitudes */
  //@{

  /**
   * Return the value of the dimensional regularization
   * parameter. Note that renormalization scale dependence is fully
   * restored in DipoleIOperator.
   */
  virtual Energy2 mu2() const { return lastSHat(); }

    /**
   * Return true, if one-loop contributions will be evaluated at amplitude level.
   */
  virtual bool oneLoopAmplitudes() const { return false; }
  //@}

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

  //@}

public:

  /**
   * Convert momenta to VBFNLO compliant format
   */
  void L5MomToDouble(Lorentz5Momentum L5Mom, double * DoubleMom) const{
    *DoubleMom++ = L5Mom.t()/GeV;
    *DoubleMom++ = L5Mom.x()/GeV;
    *DoubleMom++ = L5Mom.y()/GeV;
    *DoubleMom++ = L5Mom.z()/GeV;
    return;
  };

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

// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).

protected:

  /**
   * Copy the defined coupling constants into the Fortran common blocks.
   */
  void initCouplings();
  
  /**
   * A derieved class may override this method to figure out an integer which
   * may then be used by VBFNLO to calculate a corresponding subprocess.
   */
  virtual int getSubprocessID(int) const = 0;

  /**
   * A derieved class may override this method to initialize values
   * in common blocks once per run.
   */
  virtual void initProcess(const int &) const = 0;

  // /**
  //  * Initialize diagram containers
  //  */
  // virtual void initDiagramContainers() const = 0;

protected:

  /**
   * A vector with all allowed diagrams for all possible 
   * configurations defined by the interfaces.
   */
  mutable vector<Ptr<Tree2toNDiagram>::ptr> allPossibleDiagrams;

  /**
   * A container from which we to the allowed diagrams
   * for a given subprocess.
   */
  mutable RCPtr<DiagramContainer> theDiagramContainer;


private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  VBFNLOAmplitudeBase & operator=(const VBFNLOAmplitudeBase &);

};

}

#endif /* HERWIG_VBFNLOAmplitudeBase_H */
