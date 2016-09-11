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

  virtual void setXComb(Ptr<MatchboxMEBase>::ptr,tStdXCombPtr,int)=0;
  virtual void setME(Ptr<MatchboxMEBase>::ptr)=0;
  
  virtual void setKinematics(Ptr<MatchboxMEBase>::ptr)=0;
  virtual void clearKinematics(Ptr<MatchboxMEBase>::ptr)=0;
  virtual bool generateKinematics(Ptr<MatchboxMEBase>::ptr,const double *)=0;
  virtual bool calculateInNode() const=0;
  virtual void fillProjectors(Ptr<MatchboxMEBase>::ptr)=0;
  virtual pair<bool,bool> clusterSafe(Ptr<MatchboxMEBase>::ptr,int,int,int)=0;
  
  virtual size_t maxLegs() const=0;
  virtual Energy mergingScale() const=0;
  virtual bool matrixElementRegion(PVector particles,Energy winnerScale=0.*GeV,Energy cutscale=0.*GeV)=0;
  
  virtual CrossSection MergingDSigDR() =0;


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



};



}

#endif /* HERWIG_MergerBase_H */
