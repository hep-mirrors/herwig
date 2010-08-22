// -*- C++ -*-
//
// MEvv2ss.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEvv2ss class.
//

#include "MEvv2ss.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using ThePEG::Helicity::TensorWaveFunction;
using ThePEG::Helicity::incoming;
using ThePEG::Helicity::outgoing;
using ThePEG::Helicity::SpinfoPtr;

void MEvv2ss::doinit() {
  GeneralHardME::doinit();
  scalar_.resize(numberOfDiags());
  vector_.resize(numberOfDiags()); 
  tensor_.resize(numberOfDiags());
  flowME().resize(numberOfFlows(),
		  ProductionMatrixElement(PDT::Spin1, PDT::Spin1, 
					  PDT::Spin0, PDT::Spin0));
  diagramME().resize(numberOfDiags(),
		     ProductionMatrixElement(PDT::Spin1, PDT::Spin1, 
					     PDT::Spin0, PDT::Spin0));
  for(size_t i = 0; i < numberOfDiags(); ++i ) {
    HPDiagram dg = getProcessInfo()[i];
    if( !dg.intermediate ) {
      contact_ = dynamic_ptr_cast<AbstractVVSSVertexPtr>(dg.vertices.first);
    }
    else if(dg.channelType == HPDiagram::tChannel) {
      AbstractVSSVertexPtr vss1 = 
	dynamic_ptr_cast<AbstractVSSVertexPtr>(dg.vertices.first);
      AbstractVSSVertexPtr vss2 =
	dynamic_ptr_cast<AbstractVSSVertexPtr>(dg.vertices.second);
      scalar_[i] = make_pair(vss1, vss2);
    }
    else {
      if( dg.intermediate->iSpin() == PDT::Spin1 ) {
	AbstractVVVVertexPtr vvv =
	  dynamic_ptr_cast<AbstractVVVVertexPtr>(dg.vertices.first);
	AbstractVSSVertexPtr vss =
	  dynamic_ptr_cast<AbstractVSSVertexPtr>(dg.vertices.second);
	vector_[i] = make_pair(vvv, vss);
      }
      else if( dg.intermediate->iSpin() == PDT::Spin2 ) {
	AbstractVVTVertexPtr vvt =
	  dynamic_ptr_cast<AbstractVVTVertexPtr>(dg.vertices.first);
	AbstractSSTVertexPtr sst = 
	  dynamic_ptr_cast<AbstractSSTVertexPtr>(dg.vertices.second);
	tensor_[i] = make_pair(vvt, sst);
      }
    }
  }
}

void MEvv2ss::doinitrun() {
  GeneralHardME::doinitrun();
  flowME().resize(numberOfFlows(),
		  ProductionMatrixElement(PDT::Spin1, PDT::Spin1, 
					  PDT::Spin0, PDT::Spin0));
  diagramME().resize(numberOfDiags(),
		     ProductionMatrixElement(PDT::Spin1, PDT::Spin1, 
					     PDT::Spin0, PDT::Spin0));
}

double MEvv2ss::me2() const {
  VBVector v1(2), v2(2);
  for( size_t i = 0; i < 2; ++i ) {
    v1[i] = VectorWaveFunction(rescaledMomenta()[0],mePartonData()[0], 2*i,
			       incoming);
    v2[i] = VectorWaveFunction(rescaledMomenta()[1],mePartonData()[1], 2*i,
			       incoming);
  }
  ScalarWaveFunction sca1(rescaledMomenta()[2],mePartonData()[2],
			  Complex(1.,0.),outgoing);
  ScalarWaveFunction sca2(rescaledMomenta()[3],mePartonData()[3],
			  Complex(1.,0.),outgoing);
  double full_me(0.);
  vv2ssME(v1, v2, sca1, sca2, full_me , true);
  
#ifndef NDEBUG
  if( debugME() ) debug(full_me);
#endif

  return full_me;
}

ProductionMatrixElement 
MEvv2ss::vv2ssME(const VBVector & v1, const VBVector & v2,
		 const ScalarWaveFunction & sca1, 
		 const ScalarWaveFunction & sca2,
		 double & me2, bool first) const {
  const Energy2 m2(scale());
  const Energy masst = sca1.mass(), massu = sca2.mass();
  // weights for the selection of the diagram
  vector<double> me(numberOfDiags(), 0.);
  // weights for the selection of the colour flow
  vector<double> flow(numberOfFlows(),0.);
  //loop over vector helicities
  for(unsigned int iv1 = 0; iv1 < 2; ++iv1) {
    for(unsigned int iv2 = 0; iv2 < 2; ++iv2) {
      vector<Complex> flows(numberOfFlows(),0.);
      // loop over diagrams
      for(HPCount ix = 0; ix < numberOfDiags(); ++ix) {
	Complex diag(0.);
	const HPDiagram & current = getProcessInfo()[ix];
	// do four-point diag first
	if(current.channelType == HPDiagram::fourPoint) {
	  diag = contact_->evaluate(m2, v1[iv1], v2[iv2], sca1, sca2);
	}
	else {
	  tcPDPtr offshell = current.intermediate;
	  if(current.channelType == HPDiagram::tChannel) {
	    if(current.ordered.second) {
	      ScalarWaveFunction interS = scalar_[ix].first->
		evaluate(m2, 3, offshell, v1[iv1], sca1, masst);
	      diag = scalar_[ix].second->evaluate(m2, v2[iv2], interS, sca2);
	    }
	    else {
	      ScalarWaveFunction interS = scalar_[ix].first->
		evaluate(m2, 3, offshell, v1[iv1], sca2, massu);
	      diag = scalar_[ix].second->evaluate(m2, v2[iv2], interS, sca1);
	    }
	  }
	  else if(current.channelType == HPDiagram::sChannel) {
	    if(offshell->iSpin() == PDT::Spin1) {
	      VectorWaveFunction interV = vector_[ix].first->
		evaluate(m2, 1, offshell, v1[iv1], v2[iv2]);
	      diag = vector_[ix].second->evaluate(m2, interV, sca1, sca2);
	    }
	    else if(offshell->iSpin() == PDT::Spin2) {
	      TensorWaveFunction interT = tensor_[ix].first->
		evaluate(m2, 1, offshell, v1[iv1], v2[iv2]);
	      diag = tensor_[ix].second->evaluate(m2, sca1, sca2, interT);
	    }
	  }
	  else
	    diag = 0.;
	}
	me[ix] += norm(diag);
	diagramME()[ix](2*iv1, 2*iv2, 0, 0) = diag;
	//Compute flows
	for(size_t iy = 0; iy < current.colourFlow.size(); ++iy) {
	  flows[current.colourFlow[iy].first] += 
	    current.colourFlow[iy].second * diag;
	}
      }
      // MEs for the different colour flows
      for(unsigned int iy = 0; iy < numberOfFlows(); ++iy) 
	flowME()[iy](2*iv1, 2*iv2, 0, 0) = flows[iy];
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
  // if not computing the cross section return the selected colour flow
  if(!first) return flowME()[colourFlow()];
  me2 = selectColourFlow(flow,me,me2);
  return flowME()[colourFlow()];
}

void MEvv2ss::persistentOutput(PersistentOStream & os) const {
  os << scalar_ << vector_ << tensor_ << contact_;
}

void MEvv2ss::persistentInput(PersistentIStream & is, int) {
  is >> scalar_ >> vector_ >> tensor_ >> contact_;
}

ClassDescription<MEvv2ss> MEvv2ss::initMEvv2ss;
// Definition of the static class description member.

void MEvv2ss::Init() {

  static ClassDocumentation<MEvv2ss> documentation
    ("This class implements the ME for the vector-vector to scalar-scalar "
     "hard-process");

}

void MEvv2ss::constructVertex(tSubProPtr sub) {
  ParticleVector ext = hardParticles(sub);
  VBVector v1, v2;
  // set up the wavefunctions with real momenta
  VectorWaveFunction(v1, ext[0], incoming, false, true);
  VectorWaveFunction(v2, ext[1], incoming, false, true);
  ScalarWaveFunction sca1(ext[2], outgoing, true);
  ScalarWaveFunction sca2(ext[3], outgoing, true);
  // calculate rescaled moment
  setRescaledMomenta(ext);
  // wavefunctions with rescaled momenta
  VectorWaveFunction v1r   (rescaledMomenta()[0],
			    ext[0]->dataPtr(), incoming);
  VectorWaveFunction v2r   (rescaledMomenta()[1],
			    ext[1]->dataPtr(), incoming);
  sca1 = ScalarWaveFunction(rescaledMomenta()[2],
			    ext[2]->dataPtr(), outgoing);
  sca2 = ScalarWaveFunction(rescaledMomenta()[3],
			    ext[3]->dataPtr(), outgoing);
  for( unsigned int ihel = 0; ihel < 2; ++ihel ) {  
    v1r.reset(2*ihel);
    v1[ihel] = v1r;
    v2r.reset(2*ihel);
    v2[ihel] = v2r;
  }
  double dummy(0.);
  ProductionMatrixElement pme = vv2ssME(v1, v2, sca1, sca2, dummy , false);

#ifndef NDEBUG
  if( debugME() ) debug(dummy);
#endif

  createVertex(pme,ext);
}

void MEvv2ss::debug(double me2) const {
  if( !generator()->logfile().is_open() ) return;
  //SUSY gg>~q~q
  long id3 = abs(mePartonData()[2]->id());
  long id4 = abs(mePartonData()[3]->id());
  if( mePartonData()[0]->id() != 21 || mePartonData()[1]->id() != 21 ||
      (id3 < 1000001 &&  id3 > 1000006 ) || (id3 < 2000001 && id3 > 2000006 ) ||
      (id4 < 1000001 &&  id4 > 1000006 ) || (id4 < 2000001 &&  id4 > 2000006 ) ) 
    return;
  tcSMPtr sm = generator()->standardModel();
  double gs4 = sqr( 4.*Constants::pi*sm->alphaS(scale()) );
  int Nc = sm->Nc();
  Energy4 s2 = sqr(sHat());
  Energy2 m3s = meMomenta()[2].m2();
  Energy2 m4s = meMomenta()[3].m2();
  Energy4 spt2 = uHat()*tHat() - m3s*m4s;
  Energy4 t3s = sqr(tHat() - m3s);
  Energy4 u4s = sqr(uHat() - m4s);

  double analytic = gs4*Nc*( sqr(spt2) + s2*m3s*m4s ) * 
    ( u4s + t3s - s2/sqr(Nc) )/2./(sqr(Nc) - 1.)/s2/t3s/u4s;
  double diff = abs(analytic - me2);
  if(  diff > 1e-4 ) {
    generator()->log() 
      << mePartonData()[0]->PDGName() << ","
      << mePartonData()[1]->PDGName() << "->"
      << mePartonData()[2]->PDGName() << ","
      << mePartonData()[3]->PDGName() << "   difference: "
      << setprecision(10) << diff << "  ratio: " << analytic/me2 << '\n';
  }
}
