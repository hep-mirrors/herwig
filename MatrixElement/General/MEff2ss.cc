// -*- C++ -*-
//
// MEff2ss.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEff2ss class.
//

#include "MEff2ss.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using ThePEG::Helicity::VectorWaveFunction;
using ThePEG::Helicity::TensorWaveFunction;
using ThePEG::Helicity::incoming;
using ThePEG::Helicity::outgoing;

void MEff2ss::doinit() {
  GeneralHardME::doinit();
  fermion_.resize(numberOfDiags());
  scalar_ .resize(numberOfDiags());
  vector_ .resize(numberOfDiags());
  tensor_ .resize(numberOfDiags());
  initializeMatrixElements(PDT::Spin1Half, PDT::Spin1Half, 
			   PDT::Spin0    , PDT::Spin0    );
  for(HPCount i = 0; i < numberOfDiags(); ++i) {
    const HPDiagram & current = getProcessInfo()[i];
    if(current.channelType == HPDiagram::tChannel) {
      if(current.intermediate->iSpin() == PDT::Spin1Half)
	fermion_[i] = 
	  make_pair(dynamic_ptr_cast<AbstractFFSVertexPtr>(current.vertices.first), 
		    dynamic_ptr_cast<AbstractFFSVertexPtr>(current.vertices.second));
      else
	throw InitException() << "MEFF2ss:doinit() - t-channel"
			      << " intermediate must be a fermion "
			      << Exception::runerror;
    }
    else if(current.channelType == HPDiagram::sChannel) {
      if(current.intermediate->iSpin() == PDT::Spin0)
	scalar_[i] = 
	  make_pair(dynamic_ptr_cast<AbstractFFSVertexPtr>(current.vertices.first), 
		    dynamic_ptr_cast<AbstractSSSVertexPtr>(current.vertices.second));
      else if(current.intermediate->iSpin() == PDT::Spin1)
	vector_[i] = 
	  make_pair(dynamic_ptr_cast<AbstractFFVVertexPtr>(current.vertices.first), 
		    dynamic_ptr_cast<AbstractVSSVertexPtr>(current.vertices.second));
      else if(current.intermediate->iSpin() == PDT::Spin2)
	tensor_[i] = 
	  make_pair(dynamic_ptr_cast<AbstractFFTVertexPtr>(current.vertices.first), 
		    dynamic_ptr_cast<AbstractSSTVertexPtr>(current.vertices.second));
      else
	throw InitException() << "MEFF2ss:doinit() - s-channel"
			      << " intermediate must be a vector or tensor "
			      << Exception::runerror;
    }
    else 
      throw InitException() << "MEFF2ss:doinit() - Cannot find correct "
			    << "channel from diagram. Vertex not cast! "
			    << Exception::runerror;
  }
}

double MEff2ss::me2() const {
  // first setup  wavefunctions for external particles
  SpinorVector sp(2);
  SpinorBarVector sbar(2);
  for( unsigned int i = 0; i < 2; ++i ) {
    sp[i]   = SpinorWaveFunction   (rescaledMomenta()[0], 
				    mePartonData()[0], i, incoming);
    sbar[i] = SpinorBarWaveFunction(rescaledMomenta()[1], 
				    mePartonData()[1], i, incoming);
  }
  ScalarWaveFunction sca1(rescaledMomenta()[2], 
			  mePartonData()[2], 1., outgoing);
  ScalarWaveFunction sca2(rescaledMomenta()[3], 
			  mePartonData()[3], 1., outgoing);
  // calculate the ME
  double full_me(0.);
  ff2ssME(sp, sbar, sca1, sca2, full_me,true);  
  // debugging tests if needed
#ifndef NDEBUG
  if( debugME() ) debug(full_me);
#endif
  // return the answer
  return full_me;
}

ProductionMatrixElement 
MEff2ss::ff2ssME(const SpinorVector & sp, const SpinorBarVector & sbar, 
		 const ScalarWaveFunction & sca1, 
		 const ScalarWaveFunction & sca2,
		 double & me2, bool first) const {
  // scale
  const Energy2 q2(scale());
  // weights for the selection of the diagram
  vector<double> me(numberOfDiags(), 0.);
  // weights for the selection of the colour flow
  vector<double> flow(numberOfFlows(),0.);
  // flow over the helicities and diagrams
  for(unsigned int if1 = 0; if1 < 2; ++if1) {
    for(unsigned int if2 = 0; if2 < 2; ++if2) {
      vector<Complex> flows(numberOfFlows(),0.);
      for(HPCount ix = 0; ix < numberOfDiags(); ++ix) {
	Complex diag(0.);
	const HPDiagram & current = getProcessInfo()[ix];
	tcPDPtr internal(current.intermediate);	
	if(current.channelType == HPDiagram::tChannel &&
	   internal->iSpin() == PDT::Spin1Half) {
	  if(internal->CC()) internal=internal->CC();
	  unsigned int iopt = ( abs(sbar[if2].particle()->id()) == abs(internal->id()) ||
	  			abs(sp[if1]  .particle()->id()) == abs(internal->id())) ? 5 : 3;
	  SpinorBarWaveFunction interFB;
	  if(current.ordered.second) {
	    if(iopt==3) { 
	      interFB = fermion_[ix].second->
		evaluate(q2, iopt, internal, sbar[if2], sca2);
	    }
	    else {
	      interFB = fermion_[ix].second->
		evaluate(q2, iopt, internal, sbar[if2], sca2, 0.*GeV, 0.*GeV);
	    }
	    diag = fermion_[ix].first->evaluate(q2, sp[if1], interFB, sca1);
	  }
	  else {
	    if(iopt==3) { 
	      interFB = fermion_[ix].second->
		evaluate(q2, iopt, internal, sbar[if2], sca1);
	    }
	    else {
	      interFB = fermion_[ix].second->
		evaluate(q2, iopt, internal, sbar[if2], sca1, 0.*GeV, 0.*GeV);
	    }
	    diag = fermion_[ix].first->evaluate(q2, sp[if1], interFB, sca2);
	  }
	}
	else if(current.channelType == HPDiagram::sChannel) {
	  if(internal->iSpin() == PDT::Spin0) {
	    ScalarWaveFunction interS = scalar_[ix].first->
	      evaluate(q2, 1, internal, sp[if1], sbar[if2]);
	    diag = scalar_[ix].second->evaluate(q2, interS, sca2, sca1);
	  }
	  else if(internal->iSpin() == PDT::Spin1) {
	    VectorWaveFunction interV = vector_[ix].first->
	      evaluate(q2, 1, internal, sp[if1], sbar[if2]);
	    diag = vector_[ix].second->evaluate(q2, interV, sca2, sca1);
	  }
	  else if(internal->iSpin() == PDT::Spin2) {
	    TensorWaveFunction interT = tensor_[ix].first->
	      evaluate(q2, 1, internal, sp[if1], sbar[if2]);
	    diag = tensor_[ix].second ->evaluate(q2, sca2, sca1, interT);
	  }
	}
	// diagram
	me[ix] += norm(diag);
	diagramME()[ix](if1,if2,0,0) = diag;
	// contributions to the different colour flows
	for(unsigned int iy = 0; iy < current.colourFlow.size(); ++iy) {
	  assert(current.colourFlow[iy].first<flows.size());
	  flows[current.colourFlow[iy].first] += 
	    current.colourFlow[iy].second * diag;
	}
      }
      // MEs for the different colour flows
      for(unsigned int iy = 0; iy < numberOfFlows(); ++iy) 
	flowME()[iy](if1,if2,0,0) = flows[iy];
      // contribution to the squared matrix element
      for(unsigned int ii = 0; ii < numberOfFlows(); ++ii) 
	for(unsigned int ij = 0; ij < numberOfFlows(); ++ij)
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


void MEff2ss::persistentOutput(PersistentOStream & os) const {
  os << fermion_ << scalar_ << vector_ << tensor_;
}

void MEff2ss::persistentInput(PersistentIStream & is, int) {
  is >> fermion_ >> scalar_ >> vector_ >> tensor_;
  initializeMatrixElements(PDT::Spin1Half, PDT::Spin1Half, 
			   PDT::Spin0    , PDT::Spin0    );
}

ClassDescription<MEff2ss> MEff2ss::initMEff2ss;
// Definition of the static class description member.

void MEff2ss::Init() {

  static ClassDocumentation<MEff2ss> documentation
    ("MEff2ss implements the ME calculation of the fermion-antifermion "
     "to scalar-scalar hard process.");

}

void MEff2ss::constructVertex(tSubProPtr sub) {
  //get particles
  ParticleVector ext = hardParticles(sub);
  //First calculate wave functions with off-shell momenta
  //to calculate correct spin information
  SpinorVector sp;
  SpinorBarVector sbar;
  SpinorWaveFunction     (sp  , ext[0], incoming, false);
  SpinorBarWaveFunction  (sbar, ext[1], incoming, false);
  ScalarWaveFunction sca1(      ext[2], outgoing, true);
  ScalarWaveFunction sca2(      ext[3], outgoing, true);
  // Need to use rescale momenta to calculate matrix element
  setRescaledMomenta(ext);
  // wave functions with rescaled momenta
  SpinorWaveFunction    spr(rescaledMomenta()[0],
			    ext[0]->dataPtr(), incoming);
  SpinorBarWaveFunction sbr(rescaledMomenta()[1],
			    ext[1]->dataPtr(), incoming);
  sca1 = ScalarWaveFunction(rescaledMomenta()[2],
			    ext[2]->dataPtr(), outgoing);
  sca2 = ScalarWaveFunction(rescaledMomenta()[3],
			    ext[3]->dataPtr(), outgoing);
  for( unsigned int ihel = 0; ihel < 2; ++ihel ) {  
    spr.reset(ihel);
    sp[ihel] = spr;
    sbr.reset(ihel);
    sbar[ihel] = sbr;
  }
  double dummy(0.);
  ProductionMatrixElement pme = ff2ssME(sp, sbar, sca1, sca2, dummy,false);
#ifndef NDEBUG
  if( debugME() ) debug(dummy/36.);
#endif
  createVertex(pme,ext);
}

void MEff2ss::debug(double me2) const {
  if( !generator()->logfile().is_open() ) return;
  long id1 = mePartonData()[0]->id();
  long id2 = mePartonData()[1]->id();
  long id3 = mePartonData()[2]->id();
  long id4 = mePartonData()[3]->id();
  if( (abs(id1) != 1 && abs(id1) != 2) || (abs(id2) != 1 && abs(id2) != 2) ||
      ( abs(id3) != 1000001 && abs(id3) != 1000002 && 
        abs(id3) != 2000001 && abs(id3) != 2000002 ) || 
      ( abs(id4) != 1000001 && abs(id4) != 1000002  &&
	abs(id4) != 2000001 && abs(id4) != 2000002 ) ) return;
  tcSMPtr sm = generator()->standardModel();
  double gs4 = sqr( 4.*Constants::pi*sm->alphaS(scale()) );
  int Nc = sm->Nc();
  double Cf = (sqr(Nc) - 1)/2./Nc;
  Energy2 s(sHat());
  Energy2 mgos = sqr( getParticleData(ParticleID::SUSY_g)->mass());
  Energy2 m3s = sqr(mePartonData()[2]->mass());
  Energy2 m4s = sqr(mePartonData()[3]->mass());
  Energy4 spt2 = uHat()*tHat() - m3s*m4s;
  Energy2 tgl(tHat() - mgos), ugl(uHat() - mgos);
  unsigned int alpha = abs(id3)/1000000;
  unsigned int beta = abs(id4)/1000000;
  bool iflav = ( abs(id1) == abs(id2) );
  unsigned int oflav = ( abs(id3) - abs(id1) ) % 10;
  
  double analytic(0.);
  if( alpha != beta ) {
    if( ( id1 > 0 && id2 > 0) ||
	( id1 < 0 && id2 < 0) ) { 
      analytic = spt2/sqr(tgl);
      if( iflav ) analytic += spt2/sqr(ugl);
    }
    else {
      analytic = s*mgos/sqr(tgl);
    }
  }
  else {
    if( oflav != 0 ) {
      analytic = 2.*spt2/sqr(s);
    }
    else if( ( id1 > 0 && id2 > 0) ||
	     ( id1 < 0 && id2 < 0) ) {
      analytic = s*mgos/sqr(tgl);
      if( iflav ) {
	analytic += s*mgos/sqr(ugl) - 2.*s*mgos/Nc/tgl/ugl;
      }
      analytic /= ( iflav ? 2. : 1.);
    }
    else {
      analytic = spt2/sqr(tgl);
      if( iflav ) {
	analytic += 2.*spt2/sqr(s) - 2.*spt2/Nc/s/tgl;
      }
    }
  }
  analytic *= gs4*Cf/2./Nc;
  double diff = abs(analytic - me2);
  if( diff > 1e-4 ) {
    generator()->log() 
      << mePartonData()[0]->PDGName() << "," 	
      << mePartonData()[1]->PDGName() << "->"
      << mePartonData()[2]->PDGName() << ","
      << mePartonData()[3]->PDGName() << "   difference: " 
      << setprecision(10) << diff << "  ratio: " << analytic/me2 
      << '\n';
  }
    
}
