// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEfv2rv class.
//

#include "MEfv2rv.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using ThePEG::Helicity::incoming;
using ThePEG::Helicity::outgoing;
using ThePEG::Helicity::SpinorWaveFunction;
using ThePEG::Helicity::SpinorBarWaveFunction;
using ThePEG::Helicity::RSSpinorWaveFunction;
using ThePEG::Helicity::RSSpinorBarWaveFunction;
using ThePEG::Helicity::VectorWaveFunction;

void MEfv2rv::doinit() {
  GeneralHardME::doinit();
  fermion_.resize(numberOfDiags());
  vector_.resize(numberOfDiags());
  four_   .resize(numberOfDiags());
  initializeMatrixElements(PDT::Spin1Half, PDT::Spin1, 
			   PDT::Spin3Half, PDT::Spin1);
  for(HPCount ix = 0; ix < numberOfDiags(); ++ix) {
    HPDiagram diagram = getProcessInfo()[ix];
    PDT::Spin offspin = diagram.intermediate->iSpin();
    if ( diagram.channelType == HPDiagram::fourPoint) {
      four_[ix] = dynamic_ptr_cast<AbstractRFVVVertexPtr>(diagram.vertices.first);
    }
    else if(diagram.channelType == HPDiagram::sChannel ||
       ( diagram.channelType == HPDiagram::tChannel 
	 && offspin == PDT::Spin1Half)) {
      AbstractFFVVertexPtr vert1 = dynamic_ptr_cast<AbstractFFVVertexPtr>
		(diagram.vertices.first);
      AbstractRFVVertexPtr vert2 = dynamic_ptr_cast<AbstractRFVVertexPtr>
	(diagram.vertices.second);
      fermion_[ix] = make_pair(vert1, vert2);
    }
    else {
      if(offspin == PDT::Spin1) {
	AbstractRFVVertexPtr vert1 = dynamic_ptr_cast<AbstractRFVVertexPtr>
	  (diagram.vertices.first);
	AbstractVVVVertexPtr vert2 = dynamic_ptr_cast<AbstractVVVVertexPtr>
	  (diagram.vertices.second);
	vector_[ix] = make_pair(vert1, vert2);
      }
    }
  }
}

double MEfv2rv::me2() const {
  double fullme(0.);
  //wavefunctions for the vectors
  VBVector vecIn(2), vecOut(3);
  bool mc = !(mePartonData()[3]->mass() > ZERO);
  bool massless = mePartonData()[2]->mass()!=ZERO;
  for(unsigned int i = 0; i < 2; ++i) {
    vecIn[i] = VectorWaveFunction(rescaledMomenta()[1], mePartonData()[1], 2*i, 
				  incoming);
    vecOut[2*i] = VectorWaveFunction(rescaledMomenta()[3], mePartonData()[3], 2*i, 
				     outgoing);
  }
  if( !mc )
    vecOut[1] = VectorWaveFunction(rescaledMomenta()[3], mePartonData()[3], 1, 
				   outgoing);
  // spinor wavefunctions and me call
  if(mePartonData()[0]->id() > 0) {
    SpinorVector sp(2);
    RSSpinorBarVector spb(4);
    for(unsigned int i = 0; i < 2; ++i) {
      sp[i] = SpinorWaveFunction(rescaledMomenta()[0], mePartonData()[0], i, 
				 incoming);
    }
    for(unsigned int i = 0; i < 4; ++i) {
      if(massless && (i==1 || i==2)) continue;
      spb[i] = RSSpinorBarWaveFunction(rescaledMomenta()[2], mePartonData()[2], i, 
				       outgoing);
    }
    fv2rvHeME(sp, vecIn, spb, vecOut, mc, fullme,true);
  }
  else {
    RSSpinorVector sp(4);
    SpinorBarVector spb(2);
    for(unsigned int i = 0; i < 2; ++i) {
      spb[i] = SpinorBarWaveFunction(rescaledMomenta()[0], mePartonData()[0], i, 
 				     incoming);
    }
    for(unsigned int i = 0; i < 4; ++i) {
      if(massless && (i==1 || i==2)) continue;
      sp[i] = RSSpinorWaveFunction(rescaledMomenta()[2], mePartonData()[2], i, 
				   outgoing);
    }
    fbv2rbvHeME(spb, vecIn, sp, vecOut, mc, fullme,true);
  }
  return fullme;
}

IBPtr MEfv2rv::clone() const {
  return new_ptr(*this);
}

IBPtr MEfv2rv::fullclone() const {
  return new_ptr(*this);
}

void MEfv2rv::persistentOutput(PersistentOStream & os) const {
  os << fermion_ << vector_ << four_;
}

void MEfv2rv::persistentInput(PersistentIStream & is, int) {
  is >> fermion_ >> vector_ >> four_;
  initializeMatrixElements(PDT::Spin1Half, PDT::Spin1, 
			   PDT::Spin3Half, PDT::Spin1);
}


// The following static variable is needed for the type
// description system in ThePEG. 
DescribeClass<MEfv2rv,GeneralHardME>
  describeHerwigMEfv2rv("Herwig::MEfv2rv", "Herwig.so");

void MEfv2rv::Init() {

  static ClassDocumentation<MEfv2rv> documentation
    ("The MEfv2rv class implements the general matrix element for fv -> rv");

}

ProductionMatrixElement
MEfv2rv::fv2rvHeME(const SpinorVector & spIn,  const VBVector & vecIn,
		   const RSSpinorBarVector & spbOut,  
		   const VBVector & vecOut, bool mc, 
		   double & me2, bool first) const {
  const Energy2 q2(scale());
  // weights for the selection of the diagram
  vector<double> me(numberOfDiags(), 0.);
  // weights for the selection of the colour flow
  vector<double> flow(numberOfFlows(),0.);
  bool massless = mePartonData()[2]->mass()!=ZERO;
  me2 = 0.;
  //loop over helicities
  for(unsigned int ifh = 0; ifh < 2; ++ifh) {
    for(unsigned int ivh = 0; ivh < 2; ++ivh) {
      for(unsigned int ofh = 0; ofh < 4; ++ofh) {
	if(massless && (ofh==1 || ofh==2)) continue;
	for(unsigned int ovh = 0; ovh < 3; ++ovh) {
	  if(mc && ovh == 1) ++ovh;
   	  vector<Complex> flows(numberOfFlows(),0.);
   	  for(HPCount ix = 0; ix < numberOfDiags(); ++ix) {
	    Complex diag(0.);
	    const HPDiagram & current = getProcessInfo()[ix];
	    tcPDPtr offshell = current.intermediate;
	    if(current.channelType == HPDiagram::tChannel) {
	      //t-chan spin-1/2
	      if(offshell->iSpin() == PDT::Spin1Half) {
		if(offshell->CC()) offshell = offshell->CC();
		unsigned int iopt = abs(offshell->id())==abs(spIn[ifh].particle()->id()) ? 5 : 3;
		SpinorBarWaveFunction interFB = fermion_[ix].second->
		  evaluate(q2, iopt, offshell, spbOut[ofh], vecIn[ivh]);
		diag = fermion_[ix].first->
		  evaluate(q2, spIn[ifh], interFB, vecOut[ovh]);
	      }
	      else if(offshell->iSpin() == PDT::Spin1) {
		VectorWaveFunction interV = vector_[ix].second->
		  evaluate(q2, 3, offshell, vecIn[ivh], vecOut[ovh]);
		diag = vector_[ix].first->
		  evaluate(q2, spIn[ifh], spbOut[ofh], interV);
	      }
	      else
		assert(false);
	    }
	    else if(current.channelType == HPDiagram::sChannel) {
	      if(offshell->CC()) offshell = offshell->CC();
	      unsigned int iopt = abs(offshell->id())==abs(spIn[ifh].particle()->id()) ? 5 : 1;
	      SpinorBarWaveFunction interFB = fermion_[ix].second->
	    	evaluate(q2, iopt, offshell, spbOut[ofh], vecOut[ovh]);
	      diag = fermion_[ix].first->
	    	evaluate(q2, spIn[ifh], interFB, vecIn[ivh]);
	    }
	    else if(current.channelType == HPDiagram::fourPoint) {
	      diag = four_[ix]->evaluate(q2, spIn[ifh],spbOut[ofh],vecIn[ivh],vecOut[ovh]);
	    }
	    me[ix] += norm(diag);
	    diagramME()[ix](ifh, 2*ivh, ofh, ovh) = diag;
	    //Compute flows
	    for(size_t iy = 0; iy < current.colourFlow.size(); ++iy) {
	      assert(current.colourFlow[iy].first<flows.size());
	      flows[current.colourFlow[iy].first] += 
		current.colourFlow[iy].second * diag;
	    }
   	  }
   	  // MEs for the different colour flows
   	  for(unsigned int iy = 0; iy < numberOfFlows(); ++iy) 
   	    flowME()[iy](ifh, 2*ivh, ovh, ofh) = flows[iy];
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

ProductionMatrixElement
MEfv2rv::fbv2rbvHeME(const SpinorBarVector & spbIn, const VBVector & vecIn, 
		     const RSSpinorVector & spOut,
		     const VBVector & vecOut, bool mc, 
		     double & me2, bool first) const {
  const Energy2 q2(scale());
  // weights for the selection of the diagram
  vector<double> me(numberOfDiags(), 0.);
  // weights for the selection of the colour flow
  vector<double> flow(numberOfFlows(),0.);
  bool massless = mePartonData()[2]->mass()!=ZERO;
  me2 = 0.;
  //loop over helicities
  for(unsigned int ifh = 0; ifh < 2; ++ifh) {
    for(unsigned int ivh = 0; ivh < 2; ++ivh) {
      for(unsigned int ofh = 0; ofh < 4; ++ofh) {
	if(massless && (ofh==1 || ofh==2)) continue;
	for(unsigned int ovh = 0; ovh < 3; ++ovh) {
	  if(mc && ovh == 1) ++ovh;
   	  vector<Complex> flows(numberOfFlows(),0.);
   	  for(HPCount ix = 0; ix < numberOfDiags(); ++ix) {
   	    Complex diag(0.);
   	    const HPDiagram & current = getProcessInfo()[ix];
   	    tcPDPtr offshell = current.intermediate;
  	    if(current.channelType == HPDiagram::tChannel) {
  	      if(offshell->iSpin() == PDT::Spin1Half) {
   		SpinorBarWaveFunction interFB = fermion_[ix].first->
   		  evaluate(q2, 3, offshell, spbIn[ifh], vecOut[ovh]);
   		diag = fermion_[ix].second->
   		  evaluate(q2, spOut[ofh], interFB, vecIn[ivh]);
   	      }
   	      else if(offshell->iSpin() == PDT::Spin1) {
   		VectorWaveFunction interV = vector_[ix].first->
   		  evaluate(q2, 3, offshell, spOut[ofh], spbIn[ifh]);
   		diag = vector_[ix].second->
   		  evaluate(q2, vecIn[ivh], interV, vecOut[ovh]);
  	      }
  	      else
		assert(false);
  	    }
  	    else if(current.channelType == HPDiagram::sChannel) {
  	      if(offshell->iSpin() == PDT::Spin1Half) {
  		SpinorBarWaveFunction interFB = fermion_[ix].first->
  		  evaluate(q2, 1, offshell, spbIn[ifh], vecIn[ivh]);
  		diag = fermion_[ix].second->
  		  evaluate(q2, spOut[ofh], interFB, vecOut[ovh]);
  	      }
	      else
		assert(false);
	    }
	    else if(current.channelType == HPDiagram::fourPoint) {
	      diag = four_[ix]->evaluate(q2, spOut[ofh], spbIn[ifh],vecOut[ovh],vecIn[ivh]);
	    }
   	    me[ix] += norm(diag);
   	    diagramME()[ix](ifh, ivh, ovh, ofh) = diag;
   	    //Compute flows
   	    for(size_t iy = 0; iy < current.colourFlow.size(); ++iy) {
   	      assert(current.colourFlow[iy].first<flows.size());
   	      flows[current.colourFlow[iy].first] += 
   		current.colourFlow[iy].second * diag;
   	    }
  	  }
  	  // MEs for the different colour flows
  	  for(unsigned int iy = 0; iy < numberOfFlows(); ++iy) 
  	    flowME()[iy](ifh, ivh, ovh, ofh) = flows[iy];
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

void MEfv2rv::constructVertex(tSubProPtr sub) {
  ParticleVector ext = hardParticles(sub);
  VBVector v1, v3;
  bool mc = !(ext[2]->data().mass() > ZERO);
  bool md = !(ext[3]->data().mass() > ZERO);
  VectorWaveFunction(v1, ext[1], incoming, false, true);
  VectorWaveFunction(v3, ext[3], outgoing, true, md);
  double dummy(0.);
  // Need to use rescale momenta to calculate matrix element
  setRescaledMomenta(ext);
  // wavefunctions with rescaled momenta
  // vector
  VectorWaveFunction vir(rescaledMomenta()[1],
  			 ext[1]->dataPtr(), incoming);
  VectorWaveFunction vor(rescaledMomenta()[3],
  			 ext[3]->dataPtr(), outgoing);
  for( unsigned int ihel = 0; ihel < 2; ++ihel ) {  
    vir.reset(2*ihel);
    v1[ihel] = vir;
    vor.reset(2*ihel);
    v3[2*ihel] = vor;
  }
  if( !md ) {
    vor.reset(1);
    v3[1] = vor;
  }
  // matrix element and spinor wavefunctions
  if( ext[0]->id() > 0 ) {
    SpinorVector sp;  
    RSSpinorBarVector sbar;
    SpinorWaveFunction   (sp  , ext[0], incoming, false);
    RSSpinorBarWaveFunction(sbar, ext[2], outgoing, true);
    SpinorWaveFunction spr   (rescaledMomenta()[0],
   			      ext[0]->dataPtr(), incoming);
    RSSpinorBarWaveFunction sbr(rescaledMomenta()[2],
				ext[2]->dataPtr(), outgoing);
    for( unsigned int ihel = 0; ihel < 2; ++ihel ) {  
      spr.reset(ihel);
      sp[ihel] = spr;
    }
    for( unsigned int ihel = 0; ihel < 4; ++ihel ) {
      if(mc &&(ihel==1 || ihel==2)) continue;
      sbr.reset(ihel);
      sbar[ihel] = sbr;
    } 
    ProductionMatrixElement pme = fv2rvHeME(sp, v1, sbar, v3, md, dummy,false);
    createVertex(pme,ext);
  }
  else {
    RSSpinorVector sp;  
    SpinorBarVector sbar;
    SpinorBarWaveFunction(sbar, ext[0], incoming, false);
    RSSpinorWaveFunction(sp, ext[2], outgoing, true);
    SpinorBarWaveFunction sbr(rescaledMomenta()[0],
  			      ext[0]->dataPtr(), incoming);
    RSSpinorWaveFunction spr   (rescaledMomenta()[2],
				ext[2]->dataPtr(), outgoing);
    for( unsigned int ihel = 0; ihel < 2; ++ihel ) {  
      sbr.reset(ihel);
      sbar[ihel] = sbr;
    }
    for( unsigned int ihel = 0; ihel < 4; ++ihel ) {
      if(mc &&(ihel==1 || ihel==2)) continue;
      spr.reset(ihel);
      sp[ihel] = spr;
    }
    ProductionMatrixElement pme = fbv2rbvHeME(sbar, v1, sp, v3, md, dummy,false);
    createVertex(pme,ext);
  }
}
