// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEvv2vs class.
//

#include "MEvv2vs.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using ThePEG::Helicity::ScalarWaveFunction;
using ThePEG::Helicity::VectorWaveFunction;
using ThePEG::Helicity::incoming;
using ThePEG::Helicity::outgoing;

IBPtr MEvv2vs::clone() const {
  return new_ptr(*this);
}

IBPtr MEvv2vs::fullclone() const {
  return new_ptr(*this);
}

void MEvv2vs::persistentOutput(PersistentOStream & os) const {
  os << scalar_ << vector_;
}

void MEvv2vs::persistentInput(PersistentIStream & is, int) {
  is >> scalar_ >> vector_;
  initializeMatrixElements(PDT::Spin1, PDT::Spin1,
			   PDT::Spin1, PDT::Spin0);
}

ClassDescription<MEvv2vs> MEvv2vs::initMEvv2vs;
// Definition of the static class description member.

void MEvv2vs::Init() {

  static ClassDocumentation<MEvv2vs> documentation
    ("The MEvv2vs class implements the general matrix elements"
     " for vector vector -> vector scalar");

}

void MEvv2vs::doinit() {
  GeneralHardME::doinit();
  scalar_.resize(numberOfDiags());
  vector_.resize(numberOfDiags());
  initializeMatrixElements(PDT::Spin1, PDT::Spin1,
			   PDT::Spin1, PDT::Spin0);
  for(size_t i = 0; i < numberOfDiags(); ++i) {
    HPDiagram diag = getProcessInfo()[i];
    tcPDPtr offshell = diag.intermediate;
    assert(offshell);
    if(offshell->iSpin() == PDT::Spin0) {
      AbstractVVSVertexPtr vert1;
      AbstractVSSVertexPtr vert2;
      if(diag.channelType == HPDiagram::sChannel ||
	 (diag.channelType == HPDiagram::tChannel && diag.ordered.second)) {
	vert1 = dynamic_ptr_cast<AbstractVVSVertexPtr>(diag.vertices.first );
	vert2 = dynamic_ptr_cast<AbstractVSSVertexPtr>(diag.vertices.second);
      }
      else {
	vert1 = dynamic_ptr_cast<AbstractVVSVertexPtr>(diag.vertices.second);
	vert2 = dynamic_ptr_cast<AbstractVSSVertexPtr>(diag.vertices.first );
      }
      scalar_[i] = make_pair(vert1, vert2);
    }
    else if(offshell->iSpin() == PDT::Spin1) {
      AbstractVVVVertexPtr vert1;
      AbstractVVSVertexPtr vert2;
      if(diag.channelType == HPDiagram::sChannel ||
	 (diag.channelType == HPDiagram::tChannel && diag.ordered.second)) {
	vert1 = dynamic_ptr_cast<AbstractVVVVertexPtr>(diag.vertices.first );
	vert2 = dynamic_ptr_cast<AbstractVVSVertexPtr>(diag.vertices.second);
      }
      else {
	vert1 = dynamic_ptr_cast<AbstractVVVVertexPtr>(diag.vertices.second);
	vert2 = dynamic_ptr_cast<AbstractVVSVertexPtr>(diag.vertices.first );
      }
      vector_[i] = make_pair(vert1, vert2);
    }
  }
}

double MEvv2vs::me2() const {
  VBVector va(2), vb(2), vc(3);
  for(unsigned int i = 0; i < 2; ++i) {
    va[i] = VectorWaveFunction(rescaledMomenta()[0], mePartonData()[0], 2*i, 
			       incoming);
    vb[i] = VectorWaveFunction(rescaledMomenta()[1], mePartonData()[1], 2*i, 
			       incoming);
  }
  //always 0 and 2 polarisations
  for(unsigned int i = 0; i < 2; ++i) {
    vc[2*i] = VectorWaveFunction(rescaledMomenta()[2], mePartonData()[2], 2*i, 
				 outgoing);
  }
  bool mc  = !(mePartonData()[2]->mass() > ZERO);
  //massive vector, also 1
  if( !mc )
    vc[1] = VectorWaveFunction(rescaledMomenta()[2], mePartonData()[2], 1, 
			       outgoing);
  ScalarWaveFunction sd(rescaledMomenta()[3], mePartonData()[3], outgoing);
  double full_me(0.);
  vv2vsHeME(va, vb, vc, mc, sd, full_me,true);
  return full_me;
}

ProductionMatrixElement 
MEvv2vs::vv2vsHeME(VBVector & vin1, VBVector & vin2, 
		   VBVector & vout1, bool mc, ScalarWaveFunction & sd,
		   double & me2, bool first) const {
  const Energy2 q2(scale());
  const Energy mass = vout1[0].mass();
  // weights for the selection of the diagram
  vector<double> me(numberOfDiags(), 0.);
  // weights for the selection of the colour flow
  vector<double> flow(numberOfFlows(),0.);
  // flow over the helicities and diagrams
  for(unsigned int ihel1 = 0; ihel1 < 2; ++ihel1) { 
    for(unsigned int ihel2 = 0; ihel2 < 2; ++ihel2) {
      for(unsigned int ohel1 = 0; ohel1 < 3; ++ohel1) {
	if(mc && ohel1 == 1) ++ohel1;
	vector<Complex> flows(numberOfFlows(),0.);
	for(HPCount ix = 0; ix < numberOfDiags(); ++ix) {
	  Complex diag(0.);
	  const HPDiagram & current = getProcessInfo()[ix];
	  tcPDPtr offshell = current.intermediate;
	  if(!offshell) continue;
	    if(current.channelType == HPDiagram::sChannel) {
	      if(offshell->iSpin() == PDT::Spin0) {
		ScalarWaveFunction interS = scalar_[ix].first->
		  evaluate(q2, 1, offshell,vin1[ihel1], vin2[ihel2]);
		diag = scalar_[ix].second->
		  evaluate(q2, vout1[ohel1], sd, interS);
	      }
	      else if(offshell->iSpin() == PDT::Spin1) {
		VectorWaveFunction interV = vector_[ix].first->
		  evaluate(q2, 1, offshell, vin1[ihel1], vin2[ihel2]);
		diag = vector_[ix].second->
		  evaluate(q2, vout1[ohel1], interV, sd);
	      }
	      else 
		assert(false);
	    }
	    else if(current.channelType == HPDiagram::tChannel) {
	      if(offshell->iSpin() == PDT::Spin0) {
		if(current.ordered.second) {
		  ScalarWaveFunction interS = scalar_[ix].
		    first->evaluate(q2, 3, offshell, vin1[ihel1],vout1[ohel1], mass);
		  diag = scalar_[ix].second->
		    evaluate(q2, vin2[ihel2], sd, interS);
		}
		else {
		  ScalarWaveFunction interS = scalar_[ix].first->
		    evaluate(q2, 3, offshell, vin2[ihel2],vout1[ohel1], mass);
		  diag = scalar_[ix].second->
		    evaluate(q2, vin1[ihel1], sd, interS);
		}
	      }
	      else if(offshell->iSpin() == PDT::Spin1) {
		if(current.ordered.second) {
		  VectorWaveFunction interV = vector_[ix].
		    first->evaluate(q2, 3, offshell, vin1[ihel1],vout1[ohel1], mass);
		  diag = vector_[ix].second->
		    evaluate(q2, vin2[ihel2], interV, sd);
		}
		else {
		  VectorWaveFunction interV = vector_[ix].first->
		    evaluate(q2, 3, offshell, vin2[ihel2],vout1[ohel1], mass);
		  diag = vector_[ix].second->
		    evaluate(q2, vin1[ihel1], interV, sd);
		}
	      }
	      else 
		assert(false);
	    }
	    else
	      assert(false);
	    me[ix] += norm(diag);
	    diagramME()[ix](2*ihel1, 2*ihel2, ohel1,0) = diag;
	    //Compute flows
	    for(size_t iy = 0; iy < current.colourFlow.size(); ++iy) {
	      assert(current.colourFlow[iy].first<flows.size());
	      flows[current.colourFlow[iy].first] += 
		current.colourFlow[iy].second * diag;
	    }
	}
	// MEs for the different colour flows
	for(unsigned int iy = 0; iy < numberOfFlows(); ++iy) 
	  flowME()[iy](2*ihel1, 2*ihel2, ohel1, 0) = flows[iy];
	//Now add flows to me2 with appropriate colour factors
	for(size_t ii = 0; ii < numberOfFlows(); ++ii)
	  for(size_t ij = 0; ij < numberOfFlows(); ++ij)
	    me2 += getColourFactors()[ii][ij]*(flows[ii]*conj(flows[ij])).real();
	// contribution to the colour flow
	for(unsigned int ii = 0; ii < numberOfFlows(); ++ii) {
	  flow[ii] += getColourFactors()[ii][ii]*norm(flows[ii]);
	}
      }
    }
  }
  // if not computing the cross section return the selected colour flow
  if(!first) return flowME()[colourFlow()];
  me2 = selectColourFlow(flow,me,me2);
  return flowME()[colourFlow()];
}

void MEvv2vs::constructVertex(tSubProPtr sub) {
  ParticleVector ext = hardParticles(sub);
  // set wave functions with real momenta
  VBVector v1, v2, v3;
  VectorWaveFunction(v1, ext[0], incoming, false, true);
  VectorWaveFunction(v2, ext[1], incoming, false, true);
  //function to calculate me2 expects massless incoming vectors
  // and this constructor sets the '1' polarisation at element [2] 
  //in the vector
  bool mc  = !(ext[2]->data().mass() > ZERO);
  VectorWaveFunction(v3, ext[2], outgoing, true, mc);
  ScalarWaveFunction sd(ext[3], outgoing, true);
  // Need to use rescale momenta to calculate matrix element
  setRescaledMomenta(ext);
  // wave functions with rescaled momenta
  VectorWaveFunction vr1(rescaledMomenta()[0],
			 ext[0]->dataPtr(), incoming);
  VectorWaveFunction vr2(rescaledMomenta()[1],
			 ext[1]->dataPtr(), incoming);
  VectorWaveFunction vr3(rescaledMomenta()[2],
			 ext[2]->dataPtr(), outgoing);
  ScalarWaveFunction sr4(rescaledMomenta()[3],
			 ext[3]->dataPtr(), outgoing);
  for( unsigned int ihel = 0; ihel < 2; ++ihel ) {
    vr1.reset(2*ihel);
    v1[ihel] = vr1;
    vr2.reset(2*ihel);
    v2[ihel] = vr2;
    vr3.reset(2*ihel);
    v3[2*ihel] = vr3;
  }
  if( !mc ) {
    vr3.reset(1);
    v3[1] = vr3;
  }
  double dummy(0.);
  ProductionMatrixElement pme = vv2vsHeME(v1, v2, v3, mc, sd, dummy,false);
  createVertex(pme,ext);
}
