// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEff2ts class.
//

#include "MEff2ts.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IBPtr MEff2ts::clone() const {
  return new_ptr(*this);
}

IBPtr MEff2ts::fullclone() const {
  return new_ptr(*this);
}

void MEff2ts::persistentOutput(PersistentOStream & os) const {
  os << fermion_ << scalar_ << fourPoint_;
}

void MEff2ts::persistentInput(PersistentIStream & is, int) {
  is >> fermion_ >> scalar_ >> fourPoint_;
  initializeMatrixElements(PDT::Spin1Half, PDT::Spin1Half, 
			   PDT::Spin2    , PDT::Spin0);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEff2ts,GeneralHardME>
describeHerwigMEff2ts("Herwig::MEff2ts", "Herwig.so");

void MEff2ts::Init() {

  static ClassDocumentation<MEff2ts> documentation
    ("The MEff2ts class implements the general matrix element for "
     "fermion-antifermion -> tensor scalar");

}

double MEff2ts::me2() const {
  // first setup  wavefunctions for external particles
  SpinorVector sp(2);
  SpinorBarVector sbar(2);
  ScalarWaveFunction sca(rescaledMomenta()[3], mePartonData()[3],outgoing);
  TBVector ten(5);
  bool tMass = meMomenta()[2].mass()!=ZERO;
  for( unsigned int i = 0; i < 5; ++i ) {
    if(i<2) {
      sp[i]   = SpinorWaveFunction   (rescaledMomenta()[0], 
				      mePartonData()[0], i, incoming);
      sbar[i] = SpinorBarWaveFunction(rescaledMomenta()[1], 
				      mePartonData()[1], i, incoming);
    }
    if( tMass || i==0 || i==4) {
      ten[i] = TensorWaveFunction(rescaledMomenta()[2], mePartonData()[2],i , 
				  outgoing);
    }
  }
  // calculate the ME
  double full_me(0.);
  ffb2tsHeME(sp, sbar, ten, sca, full_me,true);  
  // return the answer
  return full_me;
}

void MEff2ts::doinit() {
  GeneralHardME::doinit();
  fermion_   .resize(numberOfDiags());
  scalar_    .resize(numberOfDiags());
  fourPoint_ .resize(numberOfDiags());
  initializeMatrixElements(PDT::Spin1Half, PDT::Spin1Half, 
			   PDT::Spin2    , PDT::Spin0);
  for(HPCount i = 0; i < numberOfDiags(); ++i) {
    const HPDiagram & current = getProcessInfo()[i];
    if(current.channelType == HPDiagram::tChannel) {
      if(current.intermediate->iSpin() != PDT::Spin1Half)
	throw InitException() << "MEff2ts:doinit() - Cannot find correct "
			      << "t-channel from diagram. Vertex not cast! "
			      << Exception::runerror;
      if( current.ordered.second ) 
	fermion_[i] = 
	  make_pair(dynamic_ptr_cast<AbstractFFTVertexPtr>(current.vertices.first), 
		    dynamic_ptr_cast<AbstractFFSVertexPtr>(current.vertices.second));
      else
	fermion_[i] = 
	  make_pair(dynamic_ptr_cast<AbstractFFTVertexPtr>(current.vertices.second), 
		    dynamic_ptr_cast<AbstractFFSVertexPtr>(current.vertices.first));
    }
    else if(current.channelType == HPDiagram::sChannel) {
      if(current.intermediate->iSpin() != PDT::Spin0)
	throw InitException() << "MEff2ts:doinit() - Cannot find correct "
			      << "s-channel from diagram. Vertex not cast! "
			      << Exception::runerror;
      scalar_[i] = 
	make_pair(dynamic_ptr_cast<AbstractFFSVertexPtr>(current.vertices.first), 
		  dynamic_ptr_cast<AbstractSSTVertexPtr>(current.vertices.second));
    }
    else if(current.channelType == HPDiagram::fourPoint) {
      fourPoint_[i] = 
	dynamic_ptr_cast<AbstractFFSTVertexPtr>(current.vertices.first);
    }
  }
}

ProductionMatrixElement MEff2ts::
ffb2tsHeME(SpinorVector & sp, SpinorBarVector & sb,
	   TBVector & ten, ScalarWaveFunction & sca, 
	   double & me2,bool first) const {
  // scale
  const Energy2 q2(scale());
  // weights for the selection of the diagram
  vector<double> me(numberOfDiags(), 0.);
  // weights for the selection of the colour flow
  vector<double> flow(numberOfFlows(),0.);
  bool tMass = meMomenta()[2].mass() != ZERO;
  // flow over the helicities and diagrams
  for(unsigned int if1 = 0; if1 < 2; ++if1) {
    for(unsigned int if2 = 0; if2 < 2; ++if2) {
      for(unsigned int it=0; it<5; ++it) {
	if( (it>0&&it<4) && !tMass ) continue;
	vector<Complex> flows(numberOfFlows(),0.);
	for(HPCount ix = 0; ix < numberOfDiags(); ++ix) {
	  Complex diag(0.);
	  const HPDiagram & current = getProcessInfo()[ix];
	  tcPDPtr internal(current.intermediate);
	  if(current.channelType == HPDiagram::tChannel) {
	    if(current.ordered.second) {
	      if(internal->CC()) internal = internal->CC();
	      SpinorBarWaveFunction interFB = fermion_[ix].second->
		evaluate(q2,5,internal,sb[if2],sca);
	      diag = fermion_[ix].first->
		evaluate(q2,sp[if1],interFB,ten[it]);
	    }
	    else {
	      SpinorWaveFunction interF = fermion_[ix].second->
		evaluate(q2,5,internal,sp[if1],sca);
	      diag = fermion_[ix].first->
		evaluate(q2,interF,sb[if2],ten[it]);
	    }
	  }
	  else if(current.channelType == HPDiagram::sChannel) {
	    ScalarWaveFunction interS = scalar_[ix].first->
	      evaluate(q2, 1, internal, sp[if1], sb[if2],sca.mass());
	    diag = scalar_[ix].second->evaluate(q2, interS, sca,ten[it]);
	  }
	  else if(current.channelType == HPDiagram::fourPoint) {
	    diag = fourPoint_[ix]->
	      evaluate(q2,sp[if1],sb[if2],sca,ten[it]);
	  }
	  // diagram
	  me[ix] += norm(diag);
	  diagramME()[ix](if1,if2,it,0) = diag;
	  // contributions to the different colour flows
	  for(unsigned int iy = 0; iy < current.colourFlow.size(); ++iy) {
	    assert(current.colourFlow[iy].first<flows.size());
	    flows[current.colourFlow[iy].first] += 
	  	current.colourFlow[iy].second * diag;
	  }
	}
	// MEs for the different colour flows
	for(unsigned int iy = 0; iy < numberOfFlows(); ++iy) 
	  flowME()[iy](if1,if2,it,0) = flows[iy];
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
  }
  // if not computing the cross section return the selected colour flow
  if(!first) return flowME()[colourFlow()];
  me2 = selectColourFlow(flow,me,me2);
  return flowME()[colourFlow()];
}

void MEff2ts::constructVertex(tSubProPtr sub) {
  ParticleVector ext = hardParticles(sub);
  vector<SpinorWaveFunction> sp;
  SpinorWaveFunction(sp, ext[0], incoming, false);
  vector<SpinorBarWaveFunction> sbar;
  SpinorBarWaveFunction(sbar, ext[1], incoming, false);
  vector<TensorWaveFunction> t1;
  bool mc  = !(ext[2]->momentum().mass() > ZERO);
  TensorWaveFunction(t1, ext[2], outgoing, true, mc);
  ScalarWaveFunction sca(ext[3], outgoing, true);
  // Need to use rescale momenta to calculate matrix element
  setRescaledMomenta(ext);
  SpinorWaveFunction spr   (rescaledMomenta()[0],
			    ext[0]->dataPtr(), incoming);
  SpinorBarWaveFunction sbr(rescaledMomenta()[1],
			    ext[1]->dataPtr(), incoming);
  TensorWaveFunction tr1   (rescaledMomenta()[2],
			    ext[2]->dataPtr(), outgoing);
  sca = ScalarWaveFunction(rescaledMomenta()[3],
			   ext[3]->dataPtr(), outgoing);
  for( unsigned int ihel = 0; ihel < 2; ++ihel ) {
    spr.reset(ihel);
    sp[ihel] = spr;
    sbr.reset(ihel);
    sbar[ihel] = sbr;
    tr1.reset(4*ihel);
    t1[4*ihel] = tr1;
  }
  if( !mc ) {
    for(unsigned int ihel=1;ihel<4;++ihel) {
      tr1.reset(ihel);
      t1[ihel] = tr1;
    }
  }
  double dummy(0.);
  ProductionMatrixElement pme = ffb2tsHeME(sp, sbar, t1, sca,dummy,false);
  createVertex(pme,ext);
}
