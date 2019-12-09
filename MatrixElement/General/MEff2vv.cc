// -*- C++ -*-
//
// MEff2vv.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEff2vv class.
//

#include "MEff2vv.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using ThePEG::Helicity::ScalarWaveFunction;
using ThePEG::Helicity::TensorWaveFunction;
using ThePEG::Helicity::incoming;
using ThePEG::Helicity::outgoing;

void MEff2vv::doinit() {
  GeneralHardME::doinit();
  fermion_.resize(numberOfDiags()); 
  vector_.resize(numberOfDiags());
  tensor_.resize(numberOfDiags());
  scalar_.resize(numberOfDiags());
  initializeMatrixElements(PDT::Spin1Half, PDT::Spin1Half,
			   PDT::Spin1, PDT::Spin1);
  for(HPCount i = 0; i < numberOfDiags(); ++i) {
    const HPDiagram & current = getProcessInfo()[i];
    if(current.channelType == HPDiagram::tChannel) {
      if(current.intermediate->iSpin() == PDT::Spin1Half)
	fermion_[i] = 
	  make_pair(dynamic_ptr_cast<AbstractFFVVertexPtr>(current.vertices.first),
		    dynamic_ptr_cast<AbstractFFVVertexPtr>(current.vertices.second));
    }
    else if(current.channelType == HPDiagram::sChannel) {
      if(current.intermediate->iSpin() == PDT::Spin0)
	scalar_[i] = 
	  make_pair(dynamic_ptr_cast<AbstractFFSVertexPtr>(current.vertices.first ),
		    dynamic_ptr_cast<AbstractVVSVertexPtr>(current.vertices.second));
      if(current.intermediate->iSpin() == PDT::Spin1)
	vector_[i] = 
	  make_pair(dynamic_ptr_cast<AbstractFFVVertexPtr>(current.vertices.first),
		    dynamic_ptr_cast<AbstractVVVVertexPtr>(current.vertices.second));
      if(current.intermediate->iSpin() == PDT::Spin2)
	tensor_[i] = 
	  make_pair(dynamic_ptr_cast<AbstractFFTVertexPtr>(current.vertices.first),
		    dynamic_ptr_cast<AbstractVVTVertexPtr>(current.vertices.second));
    }
  }
}

double MEff2vv::me2() const {
  SpinorVector sp(2);
  SpinorBarVector sbar(2);
  // vector 
  bool mc  = !(mePartonData()[2]->mass() > ZERO);
  bool md  = !(mePartonData()[3]->mass() > ZERO);
  VBVector v1(3), v2(3);  
  for( unsigned int i = 0; i < 2; ++i ) {
    sp[i] = SpinorWaveFunction(rescaledMomenta()[0], mePartonData()[0], i, incoming);
    sbar[i] = SpinorBarWaveFunction(rescaledMomenta()[1], mePartonData()[1], i, 
				    incoming);
    v1[2*i] = VectorWaveFunction(rescaledMomenta()[2], mePartonData()[2],2*i , 
    				 outgoing);
    v2[2*i] = VectorWaveFunction(rescaledMomenta()[3], mePartonData()[3], 2*i, 
    				 outgoing);
  }
  if( !mc ) v1[1] = VectorWaveFunction(rescaledMomenta()[2], mePartonData()[2], 1, 
				       outgoing);
  if( !md ) v2[1] = VectorWaveFunction(rescaledMomenta()[3], mePartonData()[3], 1, 
				       outgoing);
  double full_me(0.);
  ff2vvME(sp, sbar, v1, mc, v2, md, full_me,true);

#ifndef NDEBUG
  if(  debugME() ) debug(full_me);
#endif
  return full_me;
}

ProductionMatrixElement 
MEff2vv::ff2vvME(const SpinorVector & sp, const SpinorBarVector sbar, 
		 const VBVector & v1, bool m1, const VBVector & v2, bool m2,
		 double & me2, bool first) const {
  const Energy2 q2 = scale();
  // weights for the selection of the diagram
  vector<double> me(numberOfDiags(), 0.);
  // weights for the selection of the colour flow
  vector<double> flow(numberOfFlows(),0.);
  me2 = 0.;
  //Loop over helicities and calculate diagrams
  for(unsigned int if1 = 0; if1 < 2; ++if1) {
    for(unsigned int if2 = 0; if2 < 2; ++if2) {
      for(unsigned int vh1 = 0; vh1 < 3; ++vh1) {
 	if( vh1 == 1 && m1 ) ++vh1;
	for(unsigned int vh2 = 0; vh2 < 3; ++vh2) {
	  if( vh2 == 1 && m2 ) ++vh2;
	  vector<Complex> flows(numberOfFlows(),0.);
	  //loop and calculate diagrams
	  for(HPCount ix = 0; ix < numberOfDiags(); ++ix) {
	    const HPDiagram & current = getProcessInfo()[ix];
	    tcPDPtr offshell = current.intermediate;
	    Complex diag(0.);
	    if(current.channelType == HPDiagram::tChannel) {
	      unsigned int iopt = abs(offshell->id())==abs(sp[if1].particle()->id()) ? 5 : 3;
	      if(current.intermediate->iSpin() == PDT::Spin1Half) {
		if(current.ordered.second) {
		  SpinorWaveFunction interF = fermion_[ix].first->
		    evaluate(q2, iopt, offshell, sp[if1], v1[vh1]);
		  diag = fermion_[ix].second->
		    evaluate(q2, interF, sbar[if2],v2[vh2]);
		}
		else {
		  SpinorWaveFunction interF = fermion_[ix].first->
		    evaluate(q2, iopt , offshell, sp[if1], v2[vh2]);
		  diag = fermion_[ix].second->
		    evaluate(q2, interF, sbar[if2],v1[vh1]);
		}	      
	      }
	    }
	    else if(current.channelType == HPDiagram::sChannel) {
	      if(current.intermediate->iSpin() == PDT::Spin0) {
		ScalarWaveFunction interS = scalar_[ix].first->
		  evaluate(q2, 1, offshell, sp[if1], sbar[if2]);
		diag = scalar_[ix].second->
		  evaluate(q2, v1[vh1], v2[vh2],interS);
	      }
	      else if(current.intermediate->iSpin() == PDT::Spin1) {
		VectorWaveFunction interV = vector_[ix].first->
		  evaluate(q2, 1, offshell, sp[if1], sbar[if2]);
 		diag = vector_[ix].second->
		  evaluate(q2, v1[vh1],v2[vh2], interV);
	      }
	      else if(current.intermediate->iSpin() == PDT::Spin2) {
		TensorWaveFunction interT = tensor_[ix].first->
		  evaluate(q2, 1, offshell, sp[if1], sbar[if2]);
		diag = tensor_[ix].second->
		  evaluate(q2, v1[vh1], v2[vh2],interT);
	      }
	    }
	    me[ix] += norm(diag);
	    diagramME()[ix](if1, if2, vh1, vh2) = diag;
	    //Compute flows
	    for(size_t iy = 0; iy < current.colourFlow.size(); ++iy) {
	      assert(current.colourFlow[iy].first<flows.size());
	      flows[current.colourFlow[iy].first] += 
		current.colourFlow[iy].second * diag;
	    }
	  }
	  // MEs for the different colour flows
	  for(unsigned int iy = 0; iy < numberOfFlows(); ++iy) 
	    flowME()[iy](if1, if2, vh1, vh2) = flows[iy];
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
  }
  // if not computing the cross section return the selected colour flow
  if(!first) return flowME()[colourFlow()];
  me2 = selectColourFlow(flow,me,me2);
  return flowME()[colourFlow()];
}

void MEff2vv::persistentOutput(PersistentOStream & os) const {
  os << fermion_ << vector_ << tensor_ << scalar_;
}

void MEff2vv::persistentInput(PersistentIStream & is, int) {
  is >> fermion_ >> vector_ >> tensor_ >> scalar_;
  initializeMatrixElements(PDT::Spin1Half, PDT::Spin1Half,
			   PDT::Spin1, PDT::Spin1);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEff2vv,GeneralHardME>
describeHerwigMEff2vv("Herwig::MEff2vv", "Herwig.so");

void MEff2vv::Init() {

  static ClassDocumentation<MEff2vv> documentation
    ("This class implements the matrix element calculation of the 2->2 "
     "process, fermion-antifermion -> vector vector");

}

void MEff2vv::constructVertex(tSubProPtr sub) {
  //get particles
  ParticleVector ext = hardParticles(sub);
  vector<SpinorWaveFunction> sp;
  SpinorWaveFunction(sp, ext[0], incoming, false);
  vector<SpinorBarWaveFunction> sbar;
  SpinorBarWaveFunction(sbar, ext[1], incoming, false);
  vector<VectorWaveFunction> v1, v2;
  bool mc  = !(ext[2]->data().mass() > ZERO);
  bool md  = !(ext[3]->data().mass() > ZERO);
  VectorWaveFunction(v1, ext[2], outgoing, true, mc);
  VectorWaveFunction(v2, ext[3], outgoing, true, md);
  //Need to use rescale momenta to calculate matrix element
  setRescaledMomenta(ext);
  // wavefunctions with rescaled momenta
  SpinorWaveFunction spr   (rescaledMomenta()[0],
			    ext[0]->dataPtr(), incoming);
  SpinorBarWaveFunction sbr(rescaledMomenta()[1],
			    ext[1]->dataPtr(), incoming);
  VectorWaveFunction vr1   (rescaledMomenta()[2],
			    ext[2]->dataPtr(), outgoing);
  VectorWaveFunction vr2   (rescaledMomenta()[3],
			    ext[3]->dataPtr(), outgoing);
  for( unsigned int ihel = 0; ihel < 2; ++ihel ) {
    spr.reset(ihel);
    sp[ihel] = spr;
    sbr.reset(ihel);
    sbar[ihel] = sbr;
    vr1.reset(2*ihel);
    v1[2*ihel] = vr1;
    vr2.reset(2*ihel);
    v2[2*ihel] = vr2;
  }
  if( !mc ) {
    vr1.reset(1);
    v1[1] = vr1;
  }
  if( !md ) {
    vr2.reset(1);
    v2[1] = vr2;
  }
  double dummy(0.);
  ProductionMatrixElement pme = ff2vvME(sp, sbar, v1, mc, v2, md, dummy,false);
  
#ifndef NDEBUG
  if( debugME() ) debug(dummy);
#endif

  createVertex(pme,ext);
}

void MEff2vv::debug(double me2) const {
  if( !generator()->logfile().is_open() ) return;
  if( (mePartonData()[0]->id() != 1 && mePartonData()[0]->id() != 2) ||
      (mePartonData()[1]->id() != -1 && mePartonData()[1]->id() != -2) ||
      mePartonData()[2]->id() != 5100021 || 
      mePartonData()[3]->id() != 5100021 ) return;
  
  tcSMPtr sm = generator()->standardModel();
  double gs4 = sqr( 4.*Constants::pi*sm->alphaS(scale()) );
  Energy2 s(sHat());
  Energy2 mf2 = meMomenta()[2].m2();
  Energy2 t3(tHat() - mf2), u4(uHat() - mf2);
  double analytic = gs4*( mf2*( (57.*s/t3/u4)  - (4.*s*s*s/t3/t3/u4/u4) 
				- (108./s) )  
			  + (20.*s*s/t3/u4) - 93. + (108.*t3*u4/s/s) )/27.;
  
  double diff = abs( analytic - me2 );
  if( diff  > 1e-4 ) {
    generator()->log() 
      << mePartonData()[0]->PDGName() << ","
      << mePartonData()[1]->PDGName() << "->"
      << mePartonData()[2]->PDGName() << ","
      << mePartonData()[3]->PDGName() << "   difference: " 
      << setprecision(10) << diff << "   ratio: " 
      << analytic/me2 << '\n';
  }


}
