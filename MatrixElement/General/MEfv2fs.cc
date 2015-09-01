// -*- C++ -*-
//
// MEfv2fs.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEfv2fs class.
//

#include "MEfv2fs.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

using ThePEG::Helicity::incoming;
using ThePEG::Helicity::outgoing;

void MEfv2fs::doinit() {
  GeneralHardME::doinit();
  scalar_.resize(numberOfDiags());
  fermion_.resize(numberOfDiags());
  initializeMatrixElements(PDT::Spin1Half, PDT::Spin1, 
			   PDT::Spin1Half, PDT::Spin0);
 for(size_t ix = 0; ix < numberOfDiags(); ++ix) {
   HPDiagram curr = getProcessInfo()[ix];
   if(curr.channelType == HPDiagram::tChannel) {
     AbstractFFSVertexPtr ffs = 
       dynamic_ptr_cast<AbstractFFSVertexPtr>(curr.vertices.first);
     if( curr.intermediate->iSpin() == PDT::Spin0 ) {
       AbstractVSSVertexPtr vss = 
	 dynamic_ptr_cast<AbstractVSSVertexPtr>(curr.vertices.second);
       scalar_[ix] =  make_pair(ffs, vss); 
     }
     else {
       AbstractFFVVertexPtr ffv = 
	 dynamic_ptr_cast<AbstractFFVVertexPtr>(curr.vertices.second);
       fermion_[ix] = make_pair(ffs, ffv);
     }
   }
   else {
     AbstractFFVVertexPtr ffv = 
       dynamic_ptr_cast<AbstractFFVVertexPtr>(curr.vertices.first);
     AbstractFFSVertexPtr ffs = 
       dynamic_ptr_cast<AbstractFFSVertexPtr>(curr.vertices.second);
     fermion_[ix] = make_pair(ffs, ffv); 
   }
 }
}

double MEfv2fs::me2() const {
  //massless vector
  VecWFVector vecIn(2);
  ScalarWaveFunction scaOut(rescaledMomenta()[3], mePartonData()[3],
			    Complex(1.,0.), outgoing);
  double fullme(0.);
  if( mePartonData()[0]->id() > 0 ) {
    SpinorVector spIn(2);
    SpinorBarVector spbOut(2);
    for(size_t ih = 0; ih < 2; ++ih) {
      spIn[ih] = SpinorWaveFunction(rescaledMomenta()[0], mePartonData()[0], ih,
				    incoming);
      spbOut[ih] = SpinorBarWaveFunction(rescaledMomenta()[2], mePartonData()[2], ih,
					 outgoing);
      vecIn[ih] = VectorWaveFunction(rescaledMomenta()[1], mePartonData()[1], 2*ih,
				     incoming);
    }
    fv2fbsHeME(spIn, vecIn, spbOut, scaOut, fullme,true);
  }
  else {
    SpinorBarVector spbIn(2);
    SpinorVector spOut(2);
    for(size_t ih = 0; ih < 2; ++ih) {
      spbIn[ih] = SpinorBarWaveFunction(rescaledMomenta()[0], mePartonData()[0], ih,
					incoming);
      spOut[ih] = SpinorWaveFunction(rescaledMomenta()[2], mePartonData()[2], ih,
				     outgoing);
      vecIn[ih] = VectorWaveFunction(rescaledMomenta()[1], mePartonData()[1], 2*ih,
				     incoming);
    }
    fbv2fsHeME(spbIn, vecIn, spOut, scaOut, fullme,true);
  }

#ifndef NDEBUG
  if( debugME() ) debug(fullme);
#endif

  return fullme;
}

ProductionMatrixElement 
MEfv2fs::fv2fbsHeME(const SpinorVector & spIn, const VecWFVector & vecIn,
		    const SpinorBarVector & spbOut,
		    const ScalarWaveFunction & scaOut,
		    double & me2, bool first) const {
  const Energy2 q2(scale());
  // weights for the selection of the diagram
  vector<double> me(numberOfDiags(), 0.);
  // weights for the selection of the colour flow
  vector<double> flow(numberOfFlows(),0.);
  me2 = 0.;
  for(unsigned int ihel1 = 0; ihel1 < 2; ++ihel1) {
    for(unsigned int ihel2 = 0; ihel2 < 2; ++ihel2) {
      for(unsigned int ohel1 = 0; ohel1 < 2; ++ohel1) {
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
	      SpinorBarWaveFunction interFB = fermion_[ix].second->
		evaluate(q2, 3, offshell,spbOut[ohel1],vecIn[ihel2]);
	      diag = fermion_[ix].first->
		evaluate(q2, spIn[ihel1], interFB, scaOut);
	    }
	    else diag = 0.0;
	  }
	  else if( current.channelType == HPDiagram::sChannel ) {
	    // check if take intermediate massless
	    unsigned int propOpt = 
	      abs(offshell->id()) != abs(spIn[ihel1].particle()->id()) ? 1 : 5;
	    SpinorWaveFunction interF = fermion_[ix].second->
	      evaluate(q2, propOpt, offshell, spIn[ihel1], vecIn[ihel2]);
	    diag = fermion_[ix].first->
	      evaluate(q2, interF, spbOut[ohel1], scaOut);
	  }
	  else diag = 0.0;
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
MEfv2fs::fbv2fsHeME(const SpinorBarVector & spbIn, const VecWFVector & vecIn,
		    const SpinorVector & spOut,
		    const ScalarWaveFunction & scaOut,
		    double & me2, bool first) const {
  const Energy2 q2(scale());
  // weights for the selection of the diagram
  vector<double> me(numberOfDiags(), 0.);
  // weights for the selection of the colour flow
  vector<double> flow(numberOfFlows(),0.);
  me2 = 0.;
  for(unsigned int ihel1 = 0; ihel1 < 2; ++ihel1) {
    for(unsigned int ihel2 = 0; ihel2 < 2; ++ihel2) {
      for(unsigned int ohel1 = 0; ohel1 < 2; ++ohel1) {
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
	      SpinorBarWaveFunction interFB = fermion_[ix].first->
		evaluate(q2, 3, offshell, spbIn[ihel1], scaOut);
	      diag = fermion_[ix].second->
		evaluate(q2, spOut[ohel1], interFB, vecIn[ihel2]);
	    }
	    else diag = 0.0;
	  }
	  else if( current.channelType == HPDiagram::sChannel ) {
	    // check if take intermediate massless
	    unsigned int propOpt = 
	      abs(offshell->id()) != abs(spbIn[ihel1].particle()->id()) ? 1 : 5;
	    SpinorBarWaveFunction interFB = fermion_[ix].second->
	      evaluate(q2, propOpt, offshell, spbIn[ihel1], vecIn[ihel2]);
	    diag = fermion_[ix].first->
	      evaluate(q2, spOut[ohel1], interFB, scaOut);
	  }
	  else diag = 0.0;
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

void MEfv2fs::persistentOutput(PersistentOStream & os) const {
  os << scalar_ << fermion_;
}

void MEfv2fs::persistentInput(PersistentIStream & is, int) {
  is >> scalar_ >> fermion_;
  initializeMatrixElements(PDT::Spin1Half, PDT::Spin1, 
			   PDT::Spin1Half, PDT::Spin0);
}

ClassDescription<MEfv2fs> MEfv2fs::initMEfv2fs;
// Definition of the static class description member.

void MEfv2fs::Init() {

  static ClassDocumentation<MEfv2fs> documentation
    ("This class implements the matrix element for a fermion-vector to "
     "a fermioin-scalar.");

}

void MEfv2fs::constructVertex(tSubProPtr subp) {
  ParticleVector external = hardParticles(subp);
  //calculate production ME
  VecWFVector vecIn;
  VectorWaveFunction(vecIn, external[1], incoming, false, true);
  ScalarWaveFunction scaOut(external[3], outgoing, true);
  //Need to use rescale momenta to calculate matrix element
  setRescaledMomenta(external);
  double dummy(0.);
  if( external[0]->id() > 0 ) {
    SpinorVector spIn;
    SpinorBarVector spbOut;
    SpinorWaveFunction   (spIn,   external[0], incoming, false);
    SpinorBarWaveFunction(spbOut, external[2], outgoing, true);

    SpinorWaveFunction spr     (rescaledMomenta()[0],
				external[0]->dataPtr(), incoming);
    VectorWaveFunction vr      (rescaledMomenta()[1],
				external[1]->dataPtr(), incoming);
    SpinorBarWaveFunction sbr  (rescaledMomenta()[2],
				external[2]->dataPtr(), outgoing);
    scaOut = ScalarWaveFunction(rescaledMomenta()[3],
				external[3]->dataPtr(), outgoing);

    for( unsigned int ihel = 0; ihel < 2; ++ihel ) {  
      spr.reset(ihel);
      spIn[ihel] = spr;
      vr.reset(2*ihel);
      vecIn[ihel] = vr;
      sbr.reset(ihel);
      spbOut[ihel] = sbr;
    }
    ProductionMatrixElement prodME = fv2fbsHeME(spIn, vecIn, spbOut, 
						scaOut, dummy,false);
    createVertex(prodME,external);
  }
  else {
    SpinorBarVector spbIn;
    SpinorVector spOut;
    SpinorBarWaveFunction(spbIn, external[0], incoming, false);
    SpinorWaveFunction   (spOut, external[2], outgoing, true);

    SpinorBarWaveFunction sbr  (rescaledMomenta()[0],
				external[0]->dataPtr(), incoming);
    VectorWaveFunction     vr  (rescaledMomenta()[1],
				external[1]->dataPtr(), incoming);
    SpinorWaveFunction    spr  (rescaledMomenta()[2],
				external[2]->dataPtr(), outgoing);
    scaOut = ScalarWaveFunction(rescaledMomenta()[3],
				external[3]->dataPtr(), outgoing);
    for( unsigned int ihel = 0; ihel < 2; ++ihel ) {  
      sbr.reset(ihel);
      spbIn[ihel] = sbr;
      vr.reset(2*ihel);
      vecIn[ihel] = vr;
      spr.reset(ihel);
      spOut[ihel] = spr;
    }
    ProductionMatrixElement prodME = fbv2fsHeME(spbIn, vecIn, spOut, 
						scaOut, dummy,false);
    createVertex(prodME,external);
  }
  
#ifndef NDEBUG
  if( debugME() ) debug(dummy/96.);
#endif

}

void MEfv2fs::debug(double me2) const {
  if( !generator()->logfile().is_open() ) return;
  long id1 = abs(mePartonData()[0]->id());
  long id4 = abs(mePartonData()[3]->id());
  if( (id1 != 1 && id1 != 2) || mePartonData()[1]->id() != 21 ||
      mePartonData()[2]->id() != 1000021 ||
      (id4 != 1000001 && id4 != 1000002 &&  id4 != 2000001 && 
       id4 != 2000002) ) return;
  tcSMPtr sm = generator()->standardModel();
  double gs4 = sqr( 4.*Constants::pi*sm->alphaS(scale()) );
  int Nc = sm->Nc();
  Energy2 m3s = sqr(mePartonData()[2]->mass());
  Energy2 m4s = sqr(mePartonData()[3]->mass());
  //formula has vf->fs so swap t and u
  Energy2 s(sHat()), t3(uHat() - m3s), u4(tHat() - m4s);
  double analytic = -gs4*( u4 + 2.*(m4s - m3s)*(1. + m3s/t3 + m4s/u4) )*
    ( sqr(u4) + sqr(s) - sqr(t3)/sqr(Nc) )/s/t3/u4/4.;
  double diff = abs( analytic - me2);
  if( diff > 1e-4 ) {
    generator()->log() 
      << mePartonData()[0]->PDGName() << ","
      << mePartonData()[1]->PDGName() << "->"
      << mePartonData()[2]->PDGName() << ","
      << mePartonData()[3]->PDGName() << "   difference: " 
      << setprecision(10) << diff << "  ratio: " << analytic/me2  << '\n';
  }
    
}
