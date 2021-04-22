// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FS_QTildeShowerKinematics1to1 class.
//

#include "FS_QTildeShowerKinematics1to1.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "Herwig/Shower/QTilde/Base/ShowerParticle.h"
#include "Herwig/Shower/QTilde/QTildeShowerHandler.h"
#include "Herwig/Shower/QTilde/Kinematics/KinematicsReconstructor.h"

using namespace Herwig;
void FS_QTildeShowerKinematics1to1::
updateChildren(const tShowerParticlePtr parent, 
	       const ShowerParticleVector & children,
	       ShowerPartnerType partnerType) const {
  assert(children.size()==1);
  // update the parameters
  children[0]->showerParameters() = parent->showerParameters();
  // set up the colour connections
  parent->antiColourLine()->join(parent->colourLine());
  // make the product a child of the parent
  parent->addChild(children[0]);
  // set the momenta of the child
  children[0]->showerBasis(parent->showerBasis(),true);
  // // sort out the helicity stuff 
  // if(! dynamic_ptr_cast<tcQTildeShowerHandlerPtr>(ShowerHandler::currentHandler())->correlations()) return;
  // SpinPtr pspin(parent->spinInfo());
  // if(!pspin ||  !dynamic_ptr_cast<tcQTildeShowerHandlerPtr>(ShowerHandler::currentHandler())->spinCorrelations() ) return;
  // Energy2 t = sqr(scale())*z()*(1.-z());
  // IdList ids;
  // ids.push_back(parent->dataPtr());
  // ids.push_back(children[0]->dataPtr());
  // ids.push_back(children[1]->dataPtr());
  // // create the vertex
  // SVertexPtr vertex(new_ptr(ShowerVertex()));
  // // set the matrix element
  // vertex->ME(sudakov1to2_->matrixElement(z(),t,ids,phi(),true));
  // RhoDMatrix mapping;
  // SpinPtr inspin;
  // bool needMapping = parent->getMapping(inspin,mapping);
  // if(needMapping) vertex->incomingBasisTransform(mapping);
  // // set the incoming particle for the vertex
  // parent->spinInfo()->decayVertex(vertex);
  // for(ShowerParticleVector::const_iterator pit=children.begin();
  //     pit!=children.end();++pit) {
  //   // construct the spin info for the children
  //   (**pit).constructSpinInfo(true);
  //   // connect the spinInfo object to the vertex
  //   (*pit)->spinInfo()->productionVertex(vertex);
  // }
}

void FS_QTildeShowerKinematics1to1::updateParent(const tShowerParticlePtr parent, 
						 const ShowerParticleVector & children,
						 unsigned int pTscheme,
						 ShowerPartnerType) const {
  children[0]->virtualMass(children[0]->dataPtr()->mass());
  parent->virtualMass(children[0]->virtualMass());
}

void FS_QTildeShowerKinematics1to1::
resetChildren(const tShowerParticlePtr parent, 
	      const ShowerParticleVector & children) const {
  children[0]->showerParameters() = parent->showerParameters();




  // updateParameters(parent, children[0], children[1], false);
  // for(unsigned int ix=0;ix<children.size();++ix) {
  //   if(children[ix]->children().empty()) continue;
  //   ShowerParticleVector newChildren;
  //   for(unsigned int iy=0;iy<children[ix]->children().size();++iy)
  //     newChildren.push_back(dynamic_ptr_cast<ShowerParticlePtr>
  // 			    (children[ix]->children()[iy]));
  //   children[ix]->showerKinematics()->resetChildren(children[ix],newChildren);
  // }




}

void FS_QTildeShowerKinematics1to1::reconstructLast(const tShowerParticlePtr last,
						    Energy mass) const {
  // set beta component and consequently all missing data from that,
  // using the nominal (i.e. PDT) mass.
  Energy theMass =ZERO;
  if(!(mass > ZERO) && ShowerHandler::currentHandler()->retConstituentMasses())
    theMass = last->data().constituentMass();
  else
    theMass = mass > ZERO ? mass : last->data().mass();
  Lorentz5Momentum pVector = last->showerBasis()->pVector();
  ShowerParticle::Parameters & lastParam = last->showerParameters();
  Energy2 denom = 2. * lastParam.alpha * last->showerBasis()->p_dot_n();
  if(abs(denom)/(sqr(pVector.e())+pVector.rho2())<1e-10) {
    throw KinematicsReconstructionVeto();
  }
  lastParam.beta = ( sqr(theMass) + sqr(lastParam.pt)
  		     - sqr(lastParam.alpha) * pVector.m2() )
    / denom;
  // set that new momentum
  Lorentz5Momentum newMomentum = last->showerBasis()->
    sudakov2Momentum( lastParam.alpha, lastParam.beta,
  		      lastParam.ptx  , lastParam.pty);
  newMomentum.setMass(theMass);
  newMomentum.rescaleEnergy();
  if(last->data().stable()) {
    last->set5Momentum( newMomentum );
  }
  else {
    last->boost(last->momentum().findBoostToCM());
    last->boost(newMomentum.boostVector());
  }
}


void FS_QTildeShowerKinematics1to1::
reconstructParent(const tShowerParticlePtr parent, 
		  const ParticleVector & children ) const {
  ShowerParticlePtr c1 = dynamic_ptr_cast<ShowerParticlePtr>(children[0]);
  parent->showerParameters().beta = c1->showerParameters().beta; 
  Lorentz5Momentum pnew = c1->momentum();
  pnew.setMass(children[0]->mass());
  parent->set5Momentum( pnew );
}
