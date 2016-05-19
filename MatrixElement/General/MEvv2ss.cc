// -*- C++ -*-
//
// MEvv2ss.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
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


void MEvv2ss::doinit() {
  GeneralHardME::doinit();
  scalar1_.resize(numberOfDiags());
  scalar2_.resize(numberOfDiags());
  scalar3_.resize(numberOfDiags());
  vector_ .resize(numberOfDiags()); 
  tensor_ .resize(numberOfDiags());
  initializeMatrixElements(PDT::Spin1, PDT::Spin1, 
			   PDT::Spin0, PDT::Spin0);
  for(size_t i = 0; i < numberOfDiags(); ++i ) {
    HPDiagram dg = getProcessInfo()[i];
    if( !dg.intermediate ) {
      contact_ = dynamic_ptr_cast<AbstractVVSSVertexPtr>(dg.vertices.first);
    }
    else if(dg.channelType == HPDiagram::tChannel) {
      if (dg.intermediate->iSpin() == PDT::Spin0 ) {
	AbstractVSSVertexPtr vss1 = 
	  dynamic_ptr_cast<AbstractVSSVertexPtr>(dg.vertices.first);
	AbstractVSSVertexPtr vss2 =
	  dynamic_ptr_cast<AbstractVSSVertexPtr>(dg.vertices.second);
	scalar2_[i] = make_pair(vss1, vss2);
      }
      else if( dg.intermediate->iSpin() == PDT::Spin1 ) {
	AbstractVVSVertexPtr vvs1 = 
	  dynamic_ptr_cast<AbstractVVSVertexPtr>(dg.vertices.first);
	AbstractVVSVertexPtr vvs2 =
	  dynamic_ptr_cast<AbstractVVSVertexPtr>(dg.vertices.second);
	scalar3_[i] = make_pair(vvs1, vvs2);
      }
      else
	assert(false);
    }
    else {
      if( dg.intermediate->iSpin() == PDT::Spin0 ) {
	AbstractVVSVertexPtr vvs =
	  dynamic_ptr_cast<AbstractVVSVertexPtr>(dg.vertices.first);
	AbstractSSSVertexPtr sss =
	  dynamic_ptr_cast<AbstractSSSVertexPtr>(dg.vertices.second);
	scalar1_[i] = make_pair(vvs, sss);
      }
      else if( dg.intermediate->iSpin() == PDT::Spin1 ) {
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
	    if(offshell->iSpin() == PDT::Spin0) {
	      if(current.ordered.second) {
		ScalarWaveFunction interS = scalar2_[ix].first->
		  evaluate(m2, 3, offshell, v1[iv1], sca1, masst);
		diag = scalar2_[ix].second->evaluate(m2, v2[iv2], interS, sca2);
	      }
	      else {
		ScalarWaveFunction interS = scalar2_[ix].first->
		  evaluate(m2, 3, offshell, v1[iv1], sca2, massu);
		diag = scalar2_[ix].second->evaluate(m2, v2[iv2], interS, sca1);
	      }
	    }
	    else {
	      if(current.ordered.second) {
		VectorWaveFunction interV = scalar3_[ix].first->
		  evaluate(m2, 3, offshell, v1[iv1], sca1);
		diag = scalar3_[ix].second->evaluate(m2, v2[iv2], interV, sca2);
	      }
	      else {
		VectorWaveFunction interV = scalar3_[ix].first->
		  evaluate(m2, 3, offshell, v1[iv1], sca2);
		diag = scalar3_[ix].second->evaluate(m2, v2[iv2], interV, sca1);
	      }
	    }
	  }
	  else if(current.channelType == HPDiagram::sChannel) {
	    if(offshell->iSpin() == PDT::Spin0) {
	      ScalarWaveFunction interS = scalar1_[ix].first->
		evaluate(m2, 1, offshell, v1[iv1], v2[iv2]);
	      diag = scalar1_[ix].second->evaluate(m2, interS, sca1, sca2);
	    }
	    else if(offshell->iSpin() == PDT::Spin1) {
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
	  assert(current.colourFlow[iy].first<flows.size());
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
  os << scalar1_ << scalar2_ << scalar3_ << vector_ << tensor_ << contact_;
}

void MEvv2ss::persistentInput(PersistentIStream & is, int) {
  is >> scalar1_ >> scalar2_ >> scalar3_ >> vector_ >> tensor_ >> contact_;
  initializeMatrixElements(PDT::Spin1, PDT::Spin1, 
			   PDT::Spin0, PDT::Spin0);
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
  if( mePartonData()[0]->id() != 21 || mePartonData()[1]->id() != 21) return;
  long id3 = abs(mePartonData()[2]->id());
  long id4 = abs(mePartonData()[3]->id());
  int type = -1;
  //SUSY gg>~q~q
  if( ((id3 >= 1000001 && id3 <= 1000006 ) && (id4 >= 1000001 &&  id4 <= 1000006 ) ) ||
      ((id3 >= 2000001 && id3 <= 2000006 ) && (id4 >= 2000001 &&  id4 <= 2000006 ) ) ) {
    type = 0;
  }
  // Sextet production
  else if(mePartonData()[2]->iColour() == PDT::Colour6 &&
	  mePartonData()[3]->iColour() == PDT::Colour6bar ) {
    type = 1;
  }
  else {
    return;
  }
  double gs4 = sqr( 4.*Constants::pi*SM().alphaS(scale()));
  int Nc = SM().Nc();
  Energy4 s2 = sqr(sHat());
  Energy2 m3s = meMomenta()[2].m2();
  Energy2 m4s = meMomenta()[3].m2();
  Energy4 spt2 = uHat()*tHat() - m3s*m4s;
  Energy2 t3  = tHat()-m3s, u4  = uHat()-m4s;
  Energy4 t3s = sqr(t3)   , u4s = sqr(u4);
  Energy8 pre = gs4*(sqr(spt2) + s2*m3s*m4s);
  // matrix element
  double analytic(0.);
  // triplet scalars
  if(type==0) {
    analytic = pre*Nc*
      ( u4s + t3s - s2/sqr(Nc) )/2./(sqr(Nc) - 1.)/s2/t3s/u4s;
  }
  // sextet scalars
  else if(type==1) {
    analytic = pre*(Nc+2.)/(sqr(Nc)-1.)/Nc*
      ((Nc+2.)*(Nc-1.)/t3s/u4s - sqr(Nc)/t3/u4/s2);
  }
  double diff = abs(analytic - me2)/(analytic+me2);
  if(  diff > 1e-10 ) {
    generator()->log() 
      << mePartonData()[0]->PDGName() << ","
      << mePartonData()[1]->PDGName() << "->"
      << mePartonData()[2]->PDGName() << ","
      << mePartonData()[3]->PDGName() << "   difference: "
      << setprecision(10) << diff << "  ratio: " << analytic/me2 << '\n';
  }
}
