// -*- C++ -*-
//
// MergerBase.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MergerBase_H
#define HERWIG_MergerBase_H
//
// This is the declaration of the MergerBase class.
//

#include "ThePEG/Handlers/HandlerBase.h"
#include "ThePEG/MatrixElement/MEBase.h"
#include "Herwig/MatrixElement/Matchbox/Utility/LastMatchboxXCombInfo.h"
#include "Herwig/MatrixElement/Matchbox/Utility/MatchboxXComb.h"

#include "MatchboxMEBase.fh"

namespace Herwig {

using namespace ThePEG;
class MergerBase;
  
ThePEG_DECLARE_POINTERS(MergerBase,MergerBasePtr);

/**
 * \ingroup Matchbox
 * \author Johannes Bellm
 *
 * \brief MergerBase is the base class MergingHelpers.
 *
 * @see \ref MergerBaseInterfaces "The interfaces"
 * defined for MergerBase.
 */
class MergerBase: public HandlerBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  MergerBase();

  /**
   * The destructor.
   */
  virtual ~MergerBase();
    //define the ME region for a particle vector.
  virtual bool matrixElementRegion(PVector incoming,
                                   PVector outcoming,
                                   Energy winnerScale = ZERO,
                                   Energy cutscale = ZERO) = 0;
    // cross section of as given by the merging
  virtual CrossSection MergingDSigDR() = 0;
    // set the current xcomb, called from ME
  virtual void setXComb(Ptr<MatchboxMEBase>::ptr, tStdXCombPtr) = 0;
   // set the current ME
  virtual void setME(Ptr<MatchboxMEBase>::ptr) = 0;
    // set kinematics, called from ME
  virtual void setKinematics(Ptr<MatchboxMEBase>::ptr) = 0;
  //  clear kinematics, called from ME
  virtual void clearKinematics(Ptr<MatchboxMEBase>::ptr) = 0;
    // generate kinematics, called from ME
  virtual bool generateKinematics(Ptr<MatchboxMEBase>::ptr, const double *) = 0;
    // fill the projector the subprocess is build from
  virtual void fillProjectors(Ptr<MatchboxMEBase>::ptr) = 0;
    // return true true if a given ME phase space point is TODO
  virtual pair<bool,bool> clusterSafe(Ptr<MatchboxMEBase>::ptr, int, int, int) = 0;
    // return the current maximum legs, the shower should veto
  virtual size_t maxLegs() const = 0;
    // return the current merging scale,
  virtual Energy mergingScale() const = 0;

  //@}

public:

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();



};



}

#endif /* HERWIG_MergerBase_H */
