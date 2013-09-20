// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEfv2tf class.
//

#include "MEfv2tf.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IBPtr MEfv2tf::clone() const {
  return new_ptr(*this);
}

IBPtr MEfv2tf::fullclone() const {
  return new_ptr(*this);
}

void MEfv2tf::persistentOutput(PersistentOStream & os) const {
  os << fermion_ << vector_ << fourPoint_;
}

void MEfv2tf::persistentInput(PersistentIStream & is, int) {
  is >> fermion_ >> vector_ >> fourPoint_;
  initializeMatrixElements(PDT::Spin1Half, PDT::Spin1    , 
			   PDT::Spin2    , PDT::Spin1Half);
}

ClassDescription<MEfv2tf> MEfv2tf::initMEfv2tf;
// Definition of the static class description member.

void MEfv2tf::Init() {

  static ClassDocumentation<MEfv2tf> documentation
    ("The MEfv2tf class implements the general matrix element for "
     "fermion-vector -> tensor fermion.");

}

double MEfv2tf::me2() const {
  // check tensor mass
  bool tMass = meMomenta()[2].mass()!=ZERO;
  // first setup  wavefunctions for external particles
  SpinorVector sp(2);
  SpinorBarVector sbar(2);
  VBVector vec(3);
  TBVector ten(5);
  for( unsigned int i = 0; i < 5; ++i ) {
    if(i<2) {
      if(mePartonData()[0]->id()>0) {
	sp[i]   = SpinorWaveFunction   (rescaledMomenta()[0], 
					mePartonData()[0], i, incoming);
	sbar[i] = SpinorBarWaveFunction(rescaledMomenta()[3], 
					mePartonData()[3], i, outgoing);
      }
      else {
	sp[i]   = SpinorWaveFunction   (rescaledMomenta()[3], 
					mePartonData()[3], i, outgoing);
	sbar[i] = SpinorBarWaveFunction(rescaledMomenta()[0], 
					mePartonData()[0], i, incoming);
      }
    }
    if( tMass || i==0 || i==4) {
      ten[i] = TensorWaveFunction(rescaledMomenta()[2], mePartonData()[2],i , 
				  outgoing);
    }
    if(i!=1 &&i<3) {
      vec[i] = VectorWaveFunction(rescaledMomenta()[1], mePartonData()[1],i , 
      				  incoming);
    }
  }
  // calculate the ME
  double full_me(0.);
  if(mePartonData()[0]->id()>0) {
    fv2tfHeME(sp,vec,ten,sbar,full_me,true);
  }
  else {
    fbv2tfbHeME(sbar,vec,ten,sp,full_me,true);
  }
  // debugging tests if needed
#ifndef NDEBUG
  if( debugME() ) debug(full_me);
#endif
  // return the answer
  return full_me;
}

void MEfv2tf::doinit() {
  GeneralHardME::doinit();
  fermion_   .resize(numberOfDiags());
  vector_    .resize(numberOfDiags());
  fourPoint_ .resize(numberOfDiags());
  initializeMatrixElements(PDT::Spin1Half, PDT::Spin1    , 
			   PDT::Spin2    , PDT::Spin1Half);
  for(HPCount i = 0; i < numberOfDiags(); ++i) {
    const HPDiagram & current = getProcessInfo()[i];
    if(current.channelType == HPDiagram::sChannel) {
      if(current.intermediate->iSpin() != PDT::Spin1Half)
	throw InitException() << "MEfv2tf:doinit() - Cannot find correct "
			      << "s-channel from diagram. Vertex not cast! "
			      << Exception::runerror;
      fermion_[i] = 
	make_pair(dynamic_ptr_cast<AbstractFFTVertexPtr>(current.vertices.second), 
		  dynamic_ptr_cast<AbstractFFVVertexPtr>(current.vertices.first));
    }
    else if(current.channelType == HPDiagram::tChannel) {
      if(current.intermediate->iSpin() == PDT::Spin1Half)
	fermion_[i] = 
	  make_pair(dynamic_ptr_cast<AbstractFFTVertexPtr>(current.vertices.first), 
		    dynamic_ptr_cast<AbstractFFVVertexPtr>(current.vertices.second));
      else if(current.intermediate->iSpin() == PDT::Spin1)
	vector_[i] = 
	  make_pair(dynamic_ptr_cast<AbstractFFVVertexPtr>(current.vertices.first), 
		    dynamic_ptr_cast<AbstractVVTVertexPtr>(current.vertices.second));
      else
	throw InitException() << "MEfv2tf:doinit() - Cannot find correct "
			      << "t-channel from diagram. Vertex not cast! "
			      << Exception::runerror;
    }
    else if(current.channelType == HPDiagram::fourPoint) {
      fourPoint_[i] = 
	dynamic_ptr_cast<AbstractFFVTVertexPtr>(current.vertices.first);
    }
  }
}

ProductionMatrixElement MEfv2tf::fv2tfHeME(const SpinorVector & sp, 
					   const VBVector & vec,
					   const TBVector & ten,
					   const SpinorBarVector & sb,
					   double & full_me, bool first) const {
  // scale
  const Energy2 q2(scale());
  // weights for the selection of the diagram
  vector<double> me(numberOfDiags(), 0.);
  // weights for the selection of the colour flow
  vector<double> flow(numberOfFlows(),0.);
  bool tMass = meMomenta()[2].mass() != ZERO;
  // flow over the helicities and diagrams
  for(unsigned int if1 = 0; if1 < 2; ++if1) {
    for(unsigned int iv=0; iv<3;iv+=2) {
      for(unsigned int it=0; it<5; ++it) {
	if( (it>0&&it<4) && !tMass ) continue;
	for(unsigned int if2 = 0; if2 < 2; ++if2) {
	  vector<Complex> flows(numberOfFlows(),0.);
	  for(HPCount ix = 0; ix < numberOfDiags(); ++ix) {
 	    Complex diag(0.);
 	    const HPDiagram & current = getProcessInfo()[ix];
 	    tcPDPtr internal(current.intermediate);
 	    if(current.channelType == HPDiagram::sChannel) {
 	      SpinorWaveFunction interF = fermion_[ix].second->
 		evaluate(q2,5,internal,sp[if1],vec[iv]);
 	      diag = fermion_[ix].first->
 		evaluate(q2,interF,sb[if2],ten[it]);
 	    }
 	    else if(current.channelType == HPDiagram::tChannel) {
	      if(internal->CC()) internal=internal->CC();
 	      if(internal->iSpin()==PDT::Spin1Half) {
		unsigned int iopt = abs(internal->id())==abs(sb[if2].particle()->id()) ? 5 : 3;
 		SpinorBarWaveFunction interFB = fermion_[ix].second->
 		  evaluate(q2,iopt,internal,sb[if2],vec[iv]);
 		diag = fermion_[ix].first->
 		  evaluate(q2,sp[if1],interFB,ten[it]);
 	      }
 	      else {
 		VectorWaveFunction interV = vector_[ix].first->
 		  evaluate(q2, 1, internal, sp[if1], sb[if2]);
 		diag = vector_[ix].second->evaluate(q2, interV, vec[iv],ten[it]);
 	      }
 	    }
 	    else if(current.channelType == HPDiagram::fourPoint) {
 	      diag = fourPoint_[ix]->
 		evaluate(q2,sp[if1],sb[if2],vec[iv],ten[it]);
 	    }
  	    // diagram
  	    me[ix] += norm(diag);
  	    diagramME()[ix](if1,iv,it,if2) = diag;
  	    // contributions to the different colour flows
  	    for(unsigned int iy = 0; iy < current.colourFlow.size(); ++iy) {
	      assert(current.colourFlow[iy].first<flows.size());
  	      flows[current.colourFlow[iy].first] += 
  		current.colourFlow[iy].second * diag;
  	    }
	  }
 	  // MEs for the different colour flows
 	  for(unsigned int iy = 0; iy < numberOfFlows(); ++iy) 
 	    flowME()[iy](if1,iv,it,if2) = flows[iy];
 	  // contribution to the squared matrix element
 	  for(unsigned int ii = 0; ii < numberOfFlows(); ++ii) 
 	    for(unsigned int ij = 0; ij < numberOfFlows(); ++ij)
 	      full_me += getColourFactors()[ii][ij]*(flows[ii]*conj(flows[ij])).real();
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
  full_me = selectColourFlow(flow,me,full_me);
  return flowME()[colourFlow()];
}

ProductionMatrixElement MEfv2tf::fbv2tfbHeME(const SpinorBarVector & sb, 
					     const VBVector & vec,
					     const TBVector & ten,
					     const SpinorVector & sp,
					     double & full_me, bool first) const {
  // scale
  const Energy2 q2(scale());
  // weights for the selection of the diagram
  vector<double> me(numberOfDiags(), 0.);
  // weights for the selection of the colour flow
  vector<double> flow(numberOfFlows(),0.);
  bool tMass = meMomenta()[2].mass() != ZERO;
  // flow over the helicities and diagrams
  for(unsigned int if1 = 0; if1 < 2; ++if1) {
    for(unsigned int iv=0; iv<3;iv+=2) {
      for(unsigned int it=0; it<5; ++it) {
	if( (it>0&&it<4) && !tMass ) continue;
	for(unsigned int if2 = 0; if2 < 2; ++if2) {
	  vector<Complex> flows(numberOfFlows(),0.);
	  for(HPCount ix = 0; ix < numberOfDiags(); ++ix) {
	    Complex diag(0.);
	    const HPDiagram & current = getProcessInfo()[ix];
	    tcPDPtr internal(current.intermediate);	
	    if(current.channelType == HPDiagram::sChannel) {
	      SpinorBarWaveFunction interFB = fermion_[ix].second->
		evaluate(q2,5,internal,sb[if1],vec[iv]);
	      diag = fermion_[ix].first->
		evaluate(q2,sp[if2],interFB,ten[it]);
	    }
	    else if(current.channelType == HPDiagram::tChannel) {
 	      if(internal->iSpin()==PDT::Spin1Half) {
		SpinorWaveFunction interF = fermion_[ix].second->
		  evaluate(q2,5,internal,sp[if2],vec[iv]);
		diag = fermion_[ix].first->
		  evaluate(q2,interF,sb[if1],ten[it]);
	      }
	      else {
		VectorWaveFunction interV = vector_[ix].first->
		  evaluate(q2, 1, internal, sp[if2], sb[if1]);
		diag = vector_[ix].second->evaluate(q2, interV, vec[iv],ten[it]);
	      }
	    }
	    else if(current.channelType == HPDiagram::fourPoint) {
	      diag = fourPoint_[ix]->
		evaluate(q2,sp[if2],sb[if1],vec[iv],ten[it]);
	    }
 	    // diagram
 	    me[ix] += norm(diag);
 	    diagramME()[ix](if1,iv,it,if2) = diag;
 	    // contributions to the different colour flows
 	    for(unsigned int iy = 0; iy < current.colourFlow.size(); ++iy) {
	      assert(current.colourFlow[iy].first<flows.size());
 	      flows[current.colourFlow[iy].first] += 
 		current.colourFlow[iy].second * diag;
 	    }
	  }
 	  // MEs for the different colour flows
 	  for(unsigned int iy = 0; iy < numberOfFlows(); ++iy) 
 	    flowME()[iy](if1,iv,it,if2) = flows[iy];
  	  // contribution to the squared matrix element
  	  for(unsigned int ii = 0; ii < numberOfFlows(); ++ii) 
  	    for(unsigned int ij = 0; ij < numberOfFlows(); ++ij)
  	      full_me += getColourFactors()[ii][ij]*(flows[ii]*conj(flows[ij])).real();
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
  full_me = selectColourFlow(flow,me,full_me);
  return flowME()[colourFlow()];
}

void MEfv2tf::debug(double me2) const {
  if( !generator()->logfile().is_open() ) return;
  long id1 = mePartonData()[0]->id();
  long id2 = mePartonData()[1]->id();
  long id4 = mePartonData()[3]->id();
  if(id1==-id4||abs(id1)>5) return;
  if(id2!=ParticleID::g) return;
  unsigned int iloc(0);
  for(;iloc<vector_.size();++iloc)
    if(vector_[iloc].first) break;
  double gs = abs(vector_[iloc].second->norm());
  InvEnergy kappa = abs(vector_[iloc].first->norm())*UnitRemoval::InvE;
  Energy2 mg2 = sqr(meMomenta()[2].mass());
  double anal2 = -3./8.*sqr(gs)*sqr(kappa)/36.*(4.*sHat()*tHat()+uHat()*mg2)*
    (sqr(tHat()-mg2)+sqr(sHat()-mg2))/sHat()/tHat()/uHat();
  double diff = abs((anal2 - me2)/(anal2+me2));
  if( diff > 1e-4 ) {
    generator()->log() 
      << mePartonData()[0]->PDGName() << "," 	
      << mePartonData()[1]->PDGName() << "->"
      << mePartonData()[2]->PDGName() << ","
      << mePartonData()[3]->PDGName() << "   difference: " 
      << setprecision(10) << diff << "  ratio: " << anal2/me2 << '\n';
  }
}

void MEfv2tf::constructVertex(tSubProPtr sub) {
  ParticleVector ext = hardParticles(sub);
  VBVector v1;
  vector<TensorWaveFunction> t3;
  bool mc  = !(ext[2]->momentum().mass() > ZERO);
  SpinorVector sp;  
  SpinorBarVector sbar;
  VectorWaveFunction(v1, ext[1], incoming, false, true);
  TensorWaveFunction(t3, ext[2], outgoing, true, mc);
  double dummy(0.);
  //Need to use rescale momenta to calculate matrix element
  setRescaledMomenta(ext);
  // wavefunctions with rescaled momenta 
  VectorWaveFunction vec(rescaledMomenta()[1],
			 ext[1]->dataPtr(), incoming);
  TensorWaveFunction ten(rescaledMomenta()[2],
			 ext[2]->dataPtr(), outgoing);
  if( ext[0]->id() > 0 ) {
    SpinorWaveFunction   (sp  , ext[0], incoming, false);
    SpinorBarWaveFunction(sbar, ext[3], outgoing, true);
    SpinorWaveFunction spr   (rescaledMomenta()[0],
			      ext[0]->dataPtr(), incoming);
    SpinorBarWaveFunction sbr(rescaledMomenta()[3],
			      ext[3]->dataPtr(), outgoing);
    for( unsigned int ihel = 0; ihel < 2; ++ihel ) {  
      spr.reset(ihel);
      sp[ihel] = spr;
      vec.reset(2*ihel);
      v1[ihel] = vec;
      ten.reset(4*ihel);
      t3[4*ihel] = ten;
      sbr.reset(ihel);
      sbar[ihel] = sbr;
    }
    if( !mc ) {
      for(unsigned int ihel=1;ihel<4;++ihel) {
	ten.reset(ihel);
	t3[ihel] = ten;
      }
    }
    ProductionMatrixElement pme = fv2tfHeME(sp, v1, t3, sbar, dummy,false);
    createVertex(pme,ext);
  }
  else {
    SpinorBarWaveFunction(sbar, ext[0], incoming, false);
    SpinorWaveFunction(sp, ext[3], outgoing, true);
    SpinorBarWaveFunction sbr(rescaledMomenta()[0],
			      ext[0]->dataPtr(), incoming);
    SpinorWaveFunction spr   (rescaledMomenta()[3],
			      ext[3]->dataPtr(), outgoing);
    for( unsigned int ihel = 0; ihel < 2; ++ihel ) {  
      sbr.reset(ihel);
      sbar[ihel] = sbr;
      vec.reset(2*ihel);
      v1[ihel] = vec;
      ten.reset(4*ihel);
      t3[4*ihel] = ten;
      spr.reset(ihel);
      sp[ihel] = spr;
    }
    if( !mc ) {
      for(unsigned int ihel=1;ihel<4;++ihel) {
	ten.reset(ihel);
	t3[ihel] = ten;
      }
    }
    ProductionMatrixElement pme = fbv2tfbHeME(sbar, v1, t3, sp, dummy,false);
    createVertex(pme,ext);
  }
}
