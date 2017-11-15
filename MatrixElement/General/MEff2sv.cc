// -*- C++ -*-
//
// MEff2sv.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEff2sv class.
//

#include "MEff2sv.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using ThePEG::Helicity::incoming;
using ThePEG::Helicity::outgoing;

void MEff2sv::doinit() {
  GeneralHardME::doinit();
  scalar_.resize(numberOfDiags());
  vector_.resize(numberOfDiags());
  fermion_.resize(numberOfDiags());
  initializeMatrixElements(PDT::Spin1Half, PDT::Spin1Half,
			   PDT::Spin0,PDT::Spin1);
  for(HPCount i = 0; i < numberOfDiags(); ++i) {
    const HPDiagram & current = getProcessInfo()[i];
    if( current.channelType == HPDiagram::sChannel ) {
      if( current.intermediate->iSpin() == PDT::Spin0 )
	scalar_[i] = 
	  make_pair(dynamic_ptr_cast<AbstractFFSVertexPtr>(current.vertices.first), 
		    dynamic_ptr_cast<AbstractVSSVertexPtr>(current.vertices.second));
      else if( current.intermediate->iSpin() == PDT::Spin1 )
	vector_[i] = 
	  make_pair(dynamic_ptr_cast<AbstractFFVVertexPtr>(current.vertices.first), 
		    dynamic_ptr_cast<AbstractVVSVertexPtr>(current.vertices.second));
    }
    else if( current.channelType == HPDiagram::tChannel ) {
      if(current.intermediate->iSpin() == PDT::Spin1Half) {
	if( current.ordered.second ) 
	  fermion_[i] = 
	    make_pair(dynamic_ptr_cast<AbstractFFSVertexPtr>(current.vertices.first), 
		      dynamic_ptr_cast<AbstractFFVVertexPtr>(current.vertices.second));
	else
	  fermion_[i] = 
	    make_pair(dynamic_ptr_cast<AbstractFFSVertexPtr>(current.vertices.second), 
		      dynamic_ptr_cast<AbstractFFVVertexPtr>(current.vertices.first));
      }
    }
  }
}

void MEff2sv::persistentOutput(PersistentOStream & os) const {
  os << scalar_ << vector_ << fermion_;
}

void MEff2sv::persistentInput(PersistentIStream & is, int) {
  is >> scalar_ >> vector_ >> fermion_;
  initializeMatrixElements(PDT::Spin1Half, PDT::Spin1Half,
			   PDT::Spin0,PDT::Spin1);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEff2sv,GeneralHardME>
describeHerwigMEff2sv("Herwig::MEff2sv", "Herwig.so");

void MEff2sv::Init() {

  static ClassDocumentation<MEff2sv> documentation
    ("MEff2sv implements the ME calculation of the fermion-antifermion "
     "to vector-scalar hard process.");

}

double MEff2sv::me2() const {
  //set up wavefunctions
  SpinorVector ina(2);
  SpinorBarVector inb(2);
  VBVector outa(3);
  ScalarWaveFunction sca(rescaledMomenta()[2], mePartonData()[2], Complex(1.),
			 outgoing);
  for(unsigned int ih = 0; ih < 2; ++ih) {
    ina[ih] = SpinorWaveFunction(rescaledMomenta()[0], mePartonData()[0], 
				 ih, incoming);
    inb[ih] = SpinorBarWaveFunction(rescaledMomenta()[1], mePartonData()[1], 
				    ih, incoming);
    outa[2*ih] = VectorWaveFunction(rescaledMomenta()[3], mePartonData()[3], 
				    2*ih, outgoing);
  }
  if( mePartonData()[2]->mass() > ZERO ) {
    outa[1] = VectorWaveFunction(rescaledMomenta()[3], mePartonData()[3], 
				 1, outgoing);
  }
  double full_me(0.);
  ffb2svHeME(ina, inb, sca, outa, full_me,true);
  return full_me;
}

ProductionMatrixElement
MEff2sv::ffb2svHeME(SpinorVector & sp, SpinorBarVector & spbar,
		    ScalarWaveFunction & sca, VBVector & vec, 
		    double & me2, bool first) const {
  Energy2 m2(scale());
  bool mv = mePartonData()[2]->mass() == ZERO;
  // weights for the selection of the diagram
  vector<double> me(numberOfDiags(), 0.);
  // weights for the selection of the colour flow
  vector<double> flow(numberOfFlows(),0.);
  me2 = 0.;
  for(unsigned int ihel1 = 0; ihel1 < 2; ++ihel1) {
    for(unsigned int ihel2 = 0; ihel2 < 2; ++ihel2) {
      for(unsigned int ovhel = 0; ovhel < 3; ++ovhel) {
	if( mv && ovhel == 1 ) continue;
	vector<Complex> flows(numberOfFlows(),0.);
	for(HPCount ix = 0; ix < numberOfDiags(); ++ix) {
	  Complex diag(0.);
	  const HPDiagram & current = getProcessInfo()[ix];
	  tcPDPtr offshell(current.intermediate);
	  if( current.channelType == HPDiagram::sChannel ) {
	    if( offshell->iSpin() == PDT::Spin0 ) {
	      ScalarWaveFunction interS = scalar_[ix].first->
		evaluate(m2, 1, offshell, sp[ihel1], spbar[ihel2]);
	      diag = scalar_[ix].second->evaluate(m2, vec[ovhel], sca, interS); 
	    }
	    else if( offshell->iSpin() == PDT::Spin1 ) {
	      VectorWaveFunction interV = vector_[ix].first->
		evaluate(m2, 1, offshell, sp[ihel1], spbar[ihel2]);
	      diag = vector_[ix].second->evaluate(m2, vec[ovhel], interV, sca);
	    }
	    else diag = 0.0;
	  }
	  else if( current.channelType == HPDiagram::tChannel ) {
	    if( offshell->iSpin() == PDT::Spin1Half ) {
	      if( current.ordered.second ) {
		SpinorBarWaveFunction interFB = fermion_[ix].second->
		  evaluate(m2, 3, offshell, spbar[ihel2], vec[ovhel]);
		diag = fermion_[ix].first->
		  evaluate(m2, sp[ihel1], interFB, sca);
	      }
	      else {
		SpinorBarWaveFunction interFB = fermion_[ix].first->
		  evaluate(m2, 3, offshell, spbar[ihel2], sca);
		diag = fermion_[ix].second->
		  evaluate(m2, sp[ihel1], interFB, vec[ovhel]);
	      }
	    }
	  }
	  else diag = 0.;
	  me[ix] += norm(diag);
	  diagramME()[ix](ihel1, ihel2, 0, ovhel) = diag;
	  //Compute flows
	  for(size_t iy = 0; iy < current.colourFlow.size(); ++iy) {
	    assert(current.colourFlow[iy].first<flows.size());
	    flows[current.colourFlow[iy].first] += 
	      current.colourFlow[iy].second * diag;
	  }
	}
	// MEs for the different colour flows
	for(unsigned int iy = 0; iy < numberOfFlows(); ++iy) 
	  flowME()[iy](ihel1, ihel2, 0, ovhel) = flows[iy];
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

void MEff2sv::constructVertex(tSubProPtr sub) {
  // Hard proces external particles
  ParticleVector hdp = hardParticles(sub);  
  SpinorVector sp;
  SpinorBarVector spbar;
  VBVector vec;
  bool mv(hdp[3]->dataPtr()->mass() == ZERO);
  SpinorWaveFunction     (sp,    hdp[0], incoming, false);
  SpinorBarWaveFunction  (spbar, hdp[1], incoming, false);
  ScalarWaveFunction sca (       hdp[2], outgoing, true);
  VectorWaveFunction     (vec,   hdp[3], outgoing, true, mv);
  //Need to use rescale momenta to calculate matrix element
  setRescaledMomenta(hdp);
  // wavefunctions with rescaled momenta
  SpinorWaveFunction    spr(rescaledMomenta()[0],
			    hdp[0]->dataPtr(), incoming);
  SpinorBarWaveFunction sbr(rescaledMomenta()[1],
			    hdp[1]->dataPtr(), incoming);
  sca = ScalarWaveFunction (rescaledMomenta()[2],
			    hdp[2]->dataPtr(), outgoing);
  VectorWaveFunction     vr(rescaledMomenta()[3],
			    hdp[3]->dataPtr(), mv, outgoing);
  for( unsigned int ihel = 0; ihel < 2; ++ihel ) {  
    spr.reset(ihel);
    sp[ihel] = spr;
    sbr.reset(ihel);
    spbar[ihel] = sbr;
    vr.reset(2*ihel);
    vec[2*ihel] = vr;
  }
  if( !mv ) {
    vr.reset(1);
    vec[1] = vr;
  }
  double dummy(0.);
  ProductionMatrixElement prodme = ffb2svHeME(sp, spbar, sca, vec, dummy,false);
  createVertex(prodme,hdp);
}


