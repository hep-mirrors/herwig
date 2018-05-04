// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEfv2rs class.
//

#include "MEfv2rs.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

double MEfv2rs::me2() const {
  //massless vector
  VecWFVector vecIn(2);
  ScalarWaveFunction scaOut(rescaledMomenta()[3], mePartonData()[3],
  			    Complex(1.,0.), outgoing);
  double fullme(0.);
  bool massless = mePartonData()[2]->mass()==ZERO;
  if( mePartonData()[0]->id() > 0 ) {
    SpinorVector spIn(2);
    RSSpinorBarVector spbOut(4);
    for(size_t ih = 0; ih < 2; ++ih) {
      spIn[ih] = SpinorWaveFunction(rescaledMomenta()[0], mePartonData()[0], ih,
   				    incoming);
      vecIn[ih] = VectorWaveFunction(rescaledMomenta()[1], mePartonData()[1], 2*ih,
   				     incoming);
    }
    for(size_t ih = 0; ih < 4; ++ih) {
      if(massless && (ih==1 || ih==2)) continue;
      spbOut[ih] = RSSpinorBarWaveFunction(rescaledMomenta()[2], mePartonData()[2], ih,
					   outgoing);
    }
    fv2rbsHeME(spIn, vecIn, spbOut, scaOut, fullme,true);
  }
  else {
    SpinorBarVector spbIn(2);
    RSSpinorVector spOut(4);
    for(size_t ih = 0; ih < 2; ++ih) {
      spbIn[ih] = SpinorBarWaveFunction(rescaledMomenta()[0], mePartonData()[0], ih,
  					incoming);
      vecIn[ih] = VectorWaveFunction(rescaledMomenta()[1], mePartonData()[1], 2*ih,
  				     incoming);
    }
    for(size_t ih = 0; ih < 4; ++ih) {
      if(massless && (ih==1 || ih==2)) continue;
      spOut[ih] = RSSpinorWaveFunction(rescaledMomenta()[2], mePartonData()[2], ih,
  				     outgoing);
    }
    fbv2rsHeME(spbIn, vecIn, spOut, scaOut, fullme,true);
  }

  return fullme;
}

IBPtr MEfv2rs::clone() const {
  return new_ptr(*this);
}

IBPtr MEfv2rs::fullclone() const {
  return new_ptr(*this);
}

void MEfv2rs::doinit() {
  GeneralHardME::doinit();
  scalar_.resize(numberOfDiags());
  fermion1_.resize(numberOfDiags());
  fermion2_.resize(numberOfDiags());
  vector_.resize(numberOfDiags());
  four_.resize(numberOfDiags());
  initializeMatrixElements(PDT::Spin1Half, PDT::Spin1, 
			   PDT::Spin3Half, PDT::Spin0);
  for(size_t ix = 0; ix < numberOfDiags(); ++ix) {
    HPDiagram curr = getProcessInfo()[ix];
    if(curr.channelType == HPDiagram::tChannel) {
      if( curr.intermediate->iSpin() == PDT::Spin0 ) {
	AbstractRFSVertexPtr rfs = 
	  dynamic_ptr_cast<AbstractRFSVertexPtr>(curr.vertices.first);
	AbstractVSSVertexPtr vss = 
	  dynamic_ptr_cast<AbstractVSSVertexPtr>(curr.vertices.second);
	scalar_[ix] =  make_pair(rfs, vss); 
      }
      else if ( curr.intermediate->iSpin() == PDT::Spin1Half ) {
	AbstractFFSVertexPtr ffs = 
	  dynamic_ptr_cast<AbstractFFSVertexPtr>(curr.vertices.first);
	AbstractRFVVertexPtr rfv = 
	  dynamic_ptr_cast<AbstractRFVVertexPtr>(curr.vertices.second);
	fermion2_[ix] = make_pair(ffs, rfv);
      }
      else if ( curr.intermediate->iSpin() == PDT::Spin1 ) {
	AbstractRFVVertexPtr rfv = 
	  dynamic_ptr_cast<AbstractRFVVertexPtr>(curr.vertices.first);
	AbstractVVSVertexPtr vvs = 
	  dynamic_ptr_cast<AbstractVVSVertexPtr>(curr.vertices.second);
	vector_[ix] = make_pair(rfv, vvs);
      }
      else
	assert(false);
    }
    else if (curr.channelType == HPDiagram::sChannel)  {
      assert(curr.intermediate->iSpin() == PDT::Spin1Half );
      AbstractFFVVertexPtr ffv = 
	dynamic_ptr_cast<AbstractFFVVertexPtr>(curr.vertices.first);
      AbstractRFSVertexPtr rfs = 
	dynamic_ptr_cast<AbstractRFSVertexPtr>(curr.vertices.second);
      fermion1_[ix] = make_pair(ffv,rfs); 
    }
    else if (curr.channelType == HPDiagram::fourPoint)  {
      four_[ix] = dynamic_ptr_cast<AbstractRFVSVertexPtr>(curr.vertices.first);
    }
    else
      assert(false);
  }
}

void MEfv2rs::persistentOutput(PersistentOStream & os) const {
  os << scalar_ << fermion1_ << fermion2_ << vector_ << four_;
}

void MEfv2rs::persistentInput(PersistentIStream & is, int) {
  is >> scalar_ >> fermion1_ >> fermion2_ >> vector_ >> four_;
  initializeMatrixElements(PDT::Spin1Half, PDT::Spin1, 
			   PDT::Spin3Half, PDT::Spin0);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEfv2rs,GeneralHardME>
describeHerwigMEfv2rs("Herwig::MEfv2rs", "Herwig.so");

void MEfv2rs::Init() {

  static ClassDocumentation<MEfv2rs> documentation
    ("The MEfv2rs class implements the general matrix element for "
     "fermion-vector to RS fermion scalar");

}

ProductionMatrixElement 
MEfv2rs::fv2rbsHeME(const SpinorVector & spIn, const VecWFVector & vecIn,
		    const RSSpinorBarVector & spbOut,
		    const ScalarWaveFunction & scaOut,
		    double & me2, bool first) const {
  const Energy2 q2(scale());
  // weights for the selection of the diagram
  vector<double> me(numberOfDiags(), 0.);
  // weights for the selection of the colour flow
  vector<double> flow(numberOfFlows(),0.);
  me2 = 0.;
  bool massless = mePartonData()[2]->mass()==ZERO;
  for(unsigned int ihel1 = 0; ihel1 < 2; ++ihel1) {
    for(unsigned int ihel2 = 0; ihel2 < 2; ++ihel2) {
      for(unsigned int ohel1 = 0; ohel1 < 4; ++ohel1) {
	if(massless && (ohel1==1 || ohel1==2)) continue;
	vector<Complex> flows(numberOfFlows(),0.);
	for(size_t ix = 0; ix < numberOfDiags(); ++ix) {
	  Complex diag(0.);
	  const HPDiagram & current = getProcessInfo()[ix];
	  tcPDPtr offshell = current.intermediate;
	  if( current.channelType == HPDiagram::tChannel ) {
  	    if( offshell->iSpin() == PDT::Spin0 ) {
	      ScalarWaveFunction interS = scalar_[ix].first->
		evaluate(q2, 3, offshell, spIn[ihel1], spbOut[ohel1]);
	      diag = scalar_[ix].second->
		evaluate(q2, vecIn[ihel2], scaOut, interS);
   	    }
   	    else if( offshell->iSpin() == PDT::Spin1Half ) {
	      if(offshell->CC()) offshell = offshell->CC();
	      SpinorBarWaveFunction interFB = fermion2_[ix].second->
		evaluate(q2, 3, offshell,spbOut[ohel1],vecIn[ihel2]);
	      diag = fermion2_[ix].first->
		evaluate(q2, spIn[ihel1], interFB, scaOut);
   	    }
   	    else {
	      VectorWaveFunction interV = vector_[ix].first->
	      	evaluate(q2,3,offshell,spIn[ihel1],spbOut[ohel1]);
	      diag = vector_[ix].second->
	      	evaluate(q2, vecIn[ihel2], interV, scaOut);
  	    }
  	  }
  	  else if( current.channelType == HPDiagram::sChannel ) {
	    // check if take intermediate massless
	    unsigned int propOpt = 
	      abs(offshell->id()) != abs(spIn[ihel1].particle()->id()) ? 1 : 5;
	    SpinorWaveFunction interF = fermion1_[ix].first->
	      evaluate(q2, propOpt, offshell, spIn[ihel1], vecIn[ihel2]);
	    diag = fermion1_[ix].second->
	      evaluate(q2, interF, spbOut[ohel1], scaOut);
  	  }
  	  else if( current.channelType == HPDiagram::fourPoint ) {
	    diag = four_[ix]->
	      evaluate(q2,  spIn[ihel1], spbOut[ohel1], vecIn[ihel2], scaOut);
	  }
  	  else
  	    assert(false);
  	  me[ix] += norm(diag);
  	  diagramME()[ix](ihel1, 2*ihel2, ohel1, 0) = diag;
  	  //Compute flows
  	  for(size_t iy = 0; iy < current.colourFlow.size(); ++iy) {
  	    assert(current.colourFlow[iy].first<flows.size());
  	    flows[current.colourFlow[iy].first] += 
  	      current.colourFlow[iy].second * diag;
  	  }
  	}
  	// MEs for the different colour flows
  	for(unsigned int iy = 0; iy < numberOfFlows(); ++iy)
  	  flowME()[iy](ihel1, 2*ihel2, ohel1, 0) = flows[iy];
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

ProductionMatrixElement 
MEfv2rs::fbv2rsHeME(const SpinorBarVector & spbIn, const VecWFVector & vecIn,
		    const RSSpinorVector & spOut,
		    const ScalarWaveFunction & scaOut,
		    double & me2, bool first) const {
  const Energy2 q2(scale());
  // weights for the selection of the diagram
  vector<double> me(numberOfDiags(), 0.);
  // weights for the selection of the colour flow
  vector<double> flow(numberOfFlows(),0.);
  me2 = 0.;
  bool massless = mePartonData()[2]->mass()==ZERO;
  for(unsigned int ihel1 = 0; ihel1 < 2; ++ihel1) {
    for(unsigned int ihel2 = 0; ihel2 < 2; ++ihel2) {
      for(unsigned int ohel1 = 0; ohel1 < 4; ++ohel1) {
	if(massless && (ohel1==1 || ohel1==2)) continue;
	vector<Complex> flows(numberOfFlows(),0.);
	for(HPCount ix = 0; ix < numberOfDiags(); ++ix) {
	  Complex diag(0.);
	  const HPDiagram & current = getProcessInfo()[ix];
	  tcPDPtr offshell = current.intermediate;
  	  if( current.channelType == HPDiagram::tChannel ) {
   	    if( offshell->iSpin() == PDT::Spin0 ) {
	      ScalarWaveFunction interS = scalar_[ix].first->
		evaluate(q2, 3, offshell, spOut[ohel1], spbIn[ihel1]);
	      diag = scalar_[ix].second->
		evaluate(q2, vecIn[ihel2], interS, scaOut);
  	    }
  	    else if( offshell->iSpin() == PDT::Spin1Half ) {
	      SpinorBarWaveFunction interFB = fermion2_[ix].first->
	      	evaluate(q2, 3, offshell, spbIn[ihel1], scaOut);
	      diag = fermion2_[ix].second->
	      	evaluate(q2, spOut[ohel1], interFB, vecIn[ihel2]);
  	    }
  	    else if( offshell->iSpin() == PDT::Spin1 ) {
	      VectorWaveFunction interV = vector_[ix].first->
		evaluate(q2,3,offshell, spOut[ohel1], spbIn[ihel1]);
	      diag = vector_[ix].second->
		evaluate(q2, vecIn[ihel2], interV, scaOut);
 	    }
  	    else
  	      assert(false);
  	  }
  	  else if( current.channelType == HPDiagram::sChannel ) {
	    // check if take intermediate massless
	    unsigned int propOpt = 
	      abs(offshell->id()) != abs(spbIn[ihel1].particle()->id()) ? 1 : 5;
	    SpinorBarWaveFunction interFB = fermion1_[ix].first->
	      evaluate(q2, propOpt, offshell, spbIn[ihel1], vecIn[ihel2]);
	    diag = fermion1_[ix].second->
	      evaluate(q2, spOut[ohel1], interFB, scaOut);
  	  }
  	  else if( current.channelType == HPDiagram::fourPoint ) {
	    diag = four_[ix]-> evaluate(q2, spOut[ohel1], spbIn[ihel1],
					vecIn[ihel2], scaOut);
	  }
  	  else
  	    assert(false);
  	  me[ix] += norm(diag);
  	  diagramME()[ix](ihel1, 2*ihel2, ohel1, 0) = diag;
  	  //Compute flows
  	  for(size_t iy = 0; iy < current.colourFlow.size(); ++iy) {
  	    assert(current.colourFlow[iy].first<flows.size());
  	    flows[current.colourFlow[iy].first] += 
  	      current.colourFlow[iy].second * diag;
  	  }
  	}
  	// MEs for the different colour flows
  	for(unsigned int iy = 0; iy < numberOfFlows(); ++iy) 
  	  flowME()[iy](ihel1, 2*ihel2, ohel1, 0) = flows[iy];
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

void MEfv2rs::constructVertex(tSubProPtr subp) {
  ParticleVector external = hardParticles(subp);
  //calculate production ME
  VecWFVector vecIn;
  VectorWaveFunction(vecIn, external[1], incoming, false, true);
  ScalarWaveFunction scaOut(external[3], outgoing, true);
  //Need to use rescale momenta to calculate matrix element
  setRescaledMomenta(external);
  bool massless = mePartonData()[2]->mass()==ZERO;
  double dummy(0.);
  if( external[0]->id() > 0 ) {
    SpinorVector spIn;
    RSSpinorBarVector spbOut;
    SpinorWaveFunction   (spIn,   external[0], incoming, false);
    RSSpinorBarWaveFunction(spbOut, external[2], outgoing, true);
    
    SpinorWaveFunction spr     (rescaledMomenta()[0],
    			external[0]->dataPtr(), incoming);
    VectorWaveFunction vr      (rescaledMomenta()[1],
    			external[1]->dataPtr(), incoming);
    RSSpinorBarWaveFunction sbr  (rescaledMomenta()[2],
				  external[2]->dataPtr(), outgoing);
    scaOut = ScalarWaveFunction(rescaledMomenta()[3],
    				external[3]->dataPtr(), outgoing);
    
    for( unsigned int ihel = 0; ihel < 2; ++ihel ) {  
      spr.reset(ihel);
      spIn[ihel] = spr;
      vr.reset(2*ihel);
      vecIn[ihel] = vr;
    }
    for( unsigned int ihel = 0; ihel < 4; ++ihel ) { 
	if(massless && (ihel==1 || ihel==2)) continue;
      sbr.reset(ihel);
      spbOut[ihel] = sbr;
    }
    ProductionMatrixElement prodME = fv2rbsHeME(spIn, vecIn, spbOut, 
						scaOut, dummy,false);
    createVertex(prodME,external);
  }
  else {
    SpinorBarVector spbIn;
    RSSpinorVector spOut;
    SpinorBarWaveFunction(spbIn, external[0], incoming, false);
    RSSpinorWaveFunction   (spOut, external[2], outgoing, true);
    SpinorBarWaveFunction sbr  (rescaledMomenta()[0],
  				external[0]->dataPtr(), incoming);
    VectorWaveFunction     vr  (rescaledMomenta()[1],
  				external[1]->dataPtr(), incoming);
    RSSpinorWaveFunction    spr  (rescaledMomenta()[2],
				  external[2]->dataPtr(), outgoing);
    scaOut = ScalarWaveFunction(rescaledMomenta()[3],
  				external[3]->dataPtr(), outgoing);
    for( unsigned int ihel = 0; ihel < 2; ++ihel ) {  
      sbr.reset(ihel);
      spbIn[ihel] = sbr;
      vr.reset(2*ihel);
      vecIn[ihel] = vr;
    }
    for( unsigned int ihel = 0; ihel < 4; ++ihel ) { 
      if(massless && (ihel==1 || ihel==2)) continue;
      spr.reset(ihel);
      spOut[ihel] = spr;
    }
    ProductionMatrixElement prodME = fbv2rsHeME(spbIn, vecIn, spOut, 
  						scaOut, dummy,false);
    createVertex(prodME,external);
  }
}
