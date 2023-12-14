// -*- C++ -*-
//
// ClusterFissioner.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_DynamicClusterFissioner_H
#define HERWIG_DynamicClusterFissioner_H

#include <ThePEG/Interface/Interfaced.h>
#include "CluHadConfig.h"
#include "ClusterFissioner.h"
#include "DynamicPartonSplitter.h"
#include "DynamicGluonMassGenerator.h"

namespace Herwig {
using namespace ThePEG;

  
class DynamicClusterFissioner: public ClusterFissioner {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
   DynamicClusterFissioner();

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
   * Standard Init function used to initialize the interfaces.
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

private:

  /**
   * Private and non-existent assignment operator.
   */
  DynamicClusterFissioner & operator=(const DynamicClusterFissioner &) = delete;

  

public:

  //@{
  
  /**
   *  Split three-component cluster
   */
  virtual cutType cutThree(ClusterPtr &, tPVector & finalhadrons, bool softUEisOn){
    assert(false);//implementation missing
  }
  //@}

protected:




  InvEnergy Pqzproposal(double z, Energy qtilde, Energy Qqtilde, Energy mq) const;
  InvEnergy Pqz(double z, Energy Qqtilde, Energy mq) const;


  /**
   * The actual implementation of the dynamic fission
   */
  virtual void dynamicFission(tPPtr & ptrP, tPPtr & ptrPbar, PPtr & ptrQ, PPtr & ptrQbar) const;


  /**
   * call the dynamic fission function here
   */
  virtual void drawNewFlavour(PPtr& newPtr1, PPtr& newPtr2,const ClusterPtr& cluster) const {
    if(!cluster->particle(0)->hasAntiColour()){
      dynamicFission(cluster->particle(0),cluster->particle(1),newPtr2,newPtr1);
    } 
    else{
      dynamicFission(cluster->particle(1),cluster->particle(0),newPtr1,newPtr2);
    }
  }



  bool canSplitMinimally(tcClusterPtr clu, Energy minmass) const;
  bool canSplitMinimally(Energy Mcl, Energy m1, Energy m2, Energy m0) const;
  

  /**
   * set the momenta here
   */
  virtual pair<Energy,Energy> drawNewMasses(Energy Mc, bool soft1, bool soft2,
					    Lorentz5Momentum& pClu1, Lorentz5Momentum& pClu2,
					    tPPtr ptrQ1, Lorentz5Momentum& pQ1, 
					    tPPtr newPtr1, Lorentz5Momentum& pQone,
					    tPPtr newPtr2, Lorentz5Momentum& pQtwo,
					    tPPtr ptrQ2,  Lorentz5Momentum& pQ2) const
  {
  pQ1 = ptrQ1->momentum();
  pQ2 = ptrQ2->momentum();
  pQone = newPtr1->momentum();
  pQtwo = newPtr2->momentum();
  pClu1 = pQ1+pQone;
  pClu2 = pQ2 + pQtwo;

  pClu1.setMass(pClu1.m());
  pClu2.setMass(pClu2.m());

  pQ1.setMass(ptrQ1->data().constituentMass());
  pQ2.setMass(ptrQ2->data().constituentMass());
  pQone.setMass(newPtr1->data().constituentMass());
  pQtwo.setMass(newPtr2->data().constituentMass());
  }

  /**
   * does nothing
   * kinematics already done in drawNewMasses
   */
  virtual void calculateKinematics(const Lorentz5Momentum &pClu,
				   const Lorentz5Momentum &p0Q1,
				   const bool toHadron1, const bool toHadron2,
				   Lorentz5Momentum &pClu1, Lorentz5Momentum &pClu2,
				   Lorentz5Momentum &pQ1, Lorentz5Momentum &pQb,
				   Lorentz5Momentum &pQ2, Lorentz5Momentum &pQ2b) const {}



  /**
   *  pointer to the dynamic parton splitter
   */
  Ptr<DynamicPartonSplitter>::tptr partonSplitter() const {
    Ptr<DynamicPartonSplitter>::tptr partonSplitter_ = dynamic_ptr_cast<Ptr<DynamicPartonSplitter>::tptr>(ClusterHadronizationHandler::currentHandler()->partonSplitter());
    return partonSplitter_;
  }

  /**
   *  pointer to the dynamic gluon mass generator
   */
  Ptr<DynamicGluonMassGenerator>::tptr gluonMassGenerator() const {
    Ptr<DynamicGluonMassGenerator>::tptr gluonMassGenerator_ = dynamic_ptr_cast<Ptr<DynamicGluonMassGenerator>::tptr>(ClusterHadronizationHandler::currentHandler()->gluonMassGenerator());
    return gluonMassGenerator_;
  }


private:

  //scales for quark and gluon splitting in the fission
  Energy _Qqtilde;
  Energy _Qg2tilde;

  //switch for angular ordering in dynamic fission
  int _AngOrdFission;

   //switch for restricting the gluon virtuality in dynamic fission
  int _restrictGluon;
  //end Daniel


  
};
}

#endif /* HERWIG_DynamicClusterFissioner_H */
