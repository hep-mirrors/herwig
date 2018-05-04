// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEff2rf class.
//

#include "MEff2rf.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;


double MEff2rf::me2() const {
  // calculate the fermion spinors
  for(unsigned int ix=0;ix<3;++ix) {
    spin_[ix].clear();
    sbar_[ix].clear();
    unsigned int iy(ix);
    if(ix==2) ++iy;
    for(unsigned int ih=0;ih<2;++ih) {
      spin_[ix].push_back(SpinorWaveFunction   (rescaledMomenta()[iy],
						mePartonData()[iy],
						ih, ix<2 ? incoming : outgoing));
      sbar_[ix].push_back(SpinorBarWaveFunction(rescaledMomenta()[iy],
						mePartonData()[iy],
						ih, ix<2 ? incoming : outgoing));
    }
  }
  rs_   .clear();
  rsbar_.clear();
  // calculate the rs spinors
  bool massless = mePartonData()[2]->mass()==ZERO;
  for(unsigned int ih=0;ih<4;++ih) {
    if(massless && (ih==1 || ih==2) ) continue;
    rs_   .push_back(RSSpinorWaveFunction(rescaledMomenta()[2],
					  mePartonData()[2],ih, outgoing));
    rsbar_.push_back(RSSpinorBarWaveFunction(rescaledMomenta()[2],
					     mePartonData()[2],ih,outgoing));
  }
  double full_me(0.);
  if( mePartonData()[0]->id() > 0 && mePartonData()[1]->id() < 0) {
    if(mePartonData()[2]->id()>0) 
      ffb2rfbHeME (full_me,true);
    else
      ffb2rbfHeME (full_me,true);
  }
  else if( mePartonData()[0]->id() > 0 && mePartonData()[1]->id() > 0 )
    ff2rfHeME(full_me,true);
  else if( mePartonData()[0]->id() < 0 && mePartonData()[1]->id() < 0 )
    fbfb2rbfbHeME(full_me,true);
  else
    assert(false);
  return full_me;
}

IBPtr MEff2rf::clone() const {
  return new_ptr(*this);
}

IBPtr MEff2rf::fullclone() const {
  return new_ptr(*this);
}

void MEff2rf::doinit() {
  GeneralHardME::doinit();
  scalar_.resize(numberOfDiags());
  vector_.resize(numberOfDiags());
  initializeMatrixElements(PDT::Spin1Half, PDT::Spin1Half, 
			   PDT::Spin3Half, PDT::Spin1Half);
  for(size_t ix = 0;ix < numberOfDiags(); ++ix) {
    const HPDiagram & current = getProcessInfo()[ix];
    tcPDPtr offshell = current.intermediate;
    if(offshell->iSpin() == PDT::Spin0) {
      if(current.channelType == HPDiagram::sChannel ||
	 (current.channelType == HPDiagram::tChannel &&
	  !current.ordered.second)) { 
	AbstractFFSVertexPtr vert1 = dynamic_ptr_cast<AbstractFFSVertexPtr>
	  (current.vertices.first);
	AbstractRFSVertexPtr vert2 = dynamic_ptr_cast<AbstractRFSVertexPtr>
	  (current.vertices.second);
	scalar_[ix] = make_pair(vert1, vert2);
      }
      else {
	AbstractFFSVertexPtr vert1 = dynamic_ptr_cast<AbstractFFSVertexPtr>
	  (current.vertices.second);
	AbstractRFSVertexPtr vert2 = dynamic_ptr_cast<AbstractRFSVertexPtr>
	  (current.vertices.first );
	scalar_[ix] = make_pair(vert1, vert2);
      }
    }
    else if(offshell->iSpin() == PDT::Spin1) {
      if(current.channelType == HPDiagram::sChannel ||
	 (current.channelType == HPDiagram::tChannel &&
	  !current.ordered.second)) { 
	AbstractFFVVertexPtr vert1 = dynamic_ptr_cast<AbstractFFVVertexPtr>
	  (current.vertices.first);
	AbstractRFVVertexPtr vert2 = dynamic_ptr_cast<AbstractRFVVertexPtr>
	  (current.vertices.second);
	vector_[ix] = make_pair(vert1, vert2);
      }
      else {
	AbstractFFVVertexPtr vert1 = dynamic_ptr_cast<AbstractFFVVertexPtr>
	  (current.vertices.second);
	AbstractRFVVertexPtr vert2 = dynamic_ptr_cast<AbstractRFVVertexPtr>
	  (current.vertices.first );
	vector_[ix] = make_pair(vert1, vert2);
      }
    }
  }
}

void MEff2rf::persistentOutput(PersistentOStream & os) const {
  os << scalar_ << vector_;
}

void MEff2rf::persistentInput(PersistentIStream & is, int) {
  is >> scalar_ >> vector_;
  initializeMatrixElements(PDT::Spin1Half, PDT::Spin1Half, 
			   PDT::Spin3Half, PDT::Spin1Half);
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEff2rf,GeneralHardME>
  describeHerwigMEff2rf("Herwig::MEff2rf", "Herwig.so");

void MEff2rf::Init() {

  static ClassDocumentation<MEff2rf> documentation
    ("There is no documentation for the MEff2rf class");

}

ProductionMatrixElement 
MEff2rf::ffb2rfbHeME(double & me2, bool first) const {
  const Energy2 q2(scale());
  // weights for the selection of the diagram
  vector<double> me(numberOfDiags(), 0.);
  // weights for the selection of the colour flow
  vector<double> flow(numberOfFlows(),0.);
  // flow over the helicities and diagrams
  bool massless = mePartonData()[2]->mass()==ZERO;
  for(unsigned int ifhel1 = 0; ifhel1 < 2; ++ifhel1) {
    for(unsigned int ifhel2 = 0; ifhel2 < 2; ++ifhel2) {
      for(unsigned int ofhel1 = 0; ofhel1 < 4; ++ofhel1) {
	if(massless && (ofhel1==1 || ofhel1==2)) continue;
 	for(unsigned int ofhel2 = 0; ofhel2 < 2; ++ofhel2) {
 	  vector<Complex> flows(numberOfFlows(),0.);
 	  for(HPCount ix = 0; ix < numberOfDiags(); ++ix) {
 	    Complex diag(0.);
 	    const HPDiagram & current = getProcessInfo()[ix];
 	    tcPDPtr offshell = current.intermediate;
 	    if(current.channelType == HPDiagram::tChannel) {
 	      if(offshell->iSpin() == PDT::Spin0) {
 		if(current.ordered.second) {
		  ScalarWaveFunction interS = scalar_[ix].first->
		    evaluate(q2, 3, offshell,spin_[2][ofhel2],sbar_[1][ifhel2]);
		  diag = -scalar_[ix].second->
		    evaluate(q2, spin_[0][ifhel1],rsbar_[ofhel1],interS);
 		}
 		else {
		  ScalarWaveFunction interS = scalar_[ix].second->
		    evaluate(q2, 3, offshell,rs_[ofhel1],sbar_[1][ifhel2]);
		  diag = scalar_[ix].first->
		    evaluate(q2, spin_[0][ifhel1],sbar_[2][ofhel2],interS);
 		}
 	      }
 	      else if(offshell->iSpin() == PDT::Spin1) {
 		if(current.ordered.second) {
		  VectorWaveFunction interV = vector_[ix].first->
		    evaluate(q2, 3, offshell,spin_[2][ofhel2],sbar_[1][ifhel2]);
		  diag = -vector_[ix].second->
		    evaluate(q2, spin_[0][ifhel1], rsbar_[ofhel1], interV);
 		}
 		else {
		  VectorWaveFunction interV = vector_[ix].second->
		    evaluate(q2, 3, offshell,rs_[ofhel1],sbar_[1][ifhel2]);
		  diag = vector_[ix].first->
		    evaluate(q2, spin_[0][ifhel1], sbar_[2][ofhel2], interV);
  		}
 	      }
	      else
		assert(false);
 	    }
 	    else if(current.channelType == HPDiagram::sChannel) {
 	      if(offshell->iSpin() == PDT::Spin0) {
		ScalarWaveFunction interS = scalar_[ix].first->
		  evaluate(q2, 1, offshell, spin_[0][ifhel1], sbar_[1][ifhel2]);
		diag = scalar_[ix].second->
		  evaluate(q2,spin_[2][ofhel2],rsbar_[ofhel1], interS);
 	      }
 	      else if(offshell->iSpin() == PDT::Spin1) {
		VectorWaveFunction interV = vector_[ix].first->
		  evaluate(q2, 1, offshell,spin_[0][ifhel1], sbar_[1][ifhel2]);
		diag = vector_[ix].second->
		  evaluate(q2, spin_[2][ofhel2],rsbar_[ofhel1] , interV);
	      }
	      else
		assert(false);
	    }
 	    else
	      assert(false);
 	    me[ix] += norm(diag);
 	    diagramME()[ix](ifhel1, ifhel2, ofhel1, ofhel2) = diag;
 	    //Compute flows
 	    for(size_t iy = 0; iy < current.colourFlow.size(); ++iy) {
 	      assert(current.colourFlow[iy].first<flows.size());
 	      flows[current.colourFlow[iy].first] += 
 		current.colourFlow[iy].second * diag;
 	    }
 	  }
 	  // MEs for the different colour flows
 	  for(unsigned int iy = 0; iy < numberOfFlows(); ++iy) 
 	    flowME()[iy](ifhel1, ifhel2, ofhel1, ofhel2) = flows[iy];
 	  //Now add flows to me2 with appropriate colour factors
 	  for(size_t ii = 0; ii < numberOfFlows(); ++ii) {
 	    for(size_t ij = 0; ij < numberOfFlows(); ++ij) {
 	      me2 += getColourFactors()[ii][ij]*(flows[ii]*conj(flows[ij])).real();
 	    }
 	  }
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
MEff2rf::ffb2rbfHeME(double & me2, bool first) const {
  const Energy2 q2(scale());
  // weights for the selection of the diagram
  vector<double> me(numberOfDiags(), 0.);
  // weights for the selection of the colour flow
  vector<double> flow(numberOfFlows(),0.);
  // flow over the helicities and diagrams
  bool massless = mePartonData()[2]->mass()==ZERO;
  for(unsigned int ifhel1 = 0; ifhel1 < 2; ++ifhel1) {
    for(unsigned int ifhel2 = 0; ifhel2 < 2; ++ifhel2) {
      for(unsigned int ofhel1 = 0; ofhel1 < 4; ++ofhel1) {
	if(massless && (ofhel1==1 || ofhel1==2)) continue;
	for(unsigned int ofhel2 = 0; ofhel2 < 2; ++ofhel2) {
	  vector<Complex> flows(numberOfFlows(),0.);
	  for(HPCount ix = 0; ix < numberOfDiags(); ++ix) {
	    Complex diag(0.);
	    const HPDiagram & current = getProcessInfo()[ix];
	    tcPDPtr offshell = current.intermediate;
	    if(current.channelType == HPDiagram::tChannel) {
	      if(offshell->iSpin() == PDT::Spin0) {
	    	if(current.ordered.second) {
		  ScalarWaveFunction interS = scalar_[ix].first->
		    evaluate(q2, 3, offshell,spin_[2][ofhel2],sbar_[1][ifhel2]);
		  diag = -scalar_[ix].second->
		    evaluate(q2, spin_[0][ifhel1],rsbar_[ofhel1],interS);
		}
		else {
		  ScalarWaveFunction interS = scalar_[ix].second->
		    evaluate(q2, 3, offshell,rs_[ofhel1],sbar_[1][ifhel2]);
		  diag = scalar_[ix].first->
		    evaluate(q2, spin_[0][ifhel1],sbar_[2][ofhel2],interS);
		}
	      }
	      else if(offshell->iSpin() == PDT::Spin1) {
		if(current.ordered.second) {
		  VectorWaveFunction interV = vector_[ix].first->
		    evaluate(q2, 3, offshell,spin_[2][ofhel2],sbar_[1][ifhel2]);
		  diag = -vector_[ix].second->
		    evaluate(q2, spin_[0][ifhel1], rsbar_[ofhel1], interV);
		}
		else {
		  VectorWaveFunction interV = vector_[ix].second->
		    evaluate(q2, 3, offshell,rs_[ofhel1],sbar_[1][ifhel2]);
		  diag = vector_[ix].first->
		    evaluate(q2, spin_[0][ifhel1], sbar_[3][ofhel2], interV);
		}
	      }
	      else
		assert(false);
	    }
	    else if(current.channelType == HPDiagram::sChannel) {
	      if(offshell->iSpin() == PDT::Spin0) {
		ScalarWaveFunction interS = scalar_[ix].second->
		  evaluate(q2, 1, offshell,rs_[ofhel1],sbar_[2][ofhel2]);
		diag = -scalar_[ix].first->
		  evaluate(q2, spin_[0][ifhel1], sbar_[1][ifhel2], interS);
	      }
	      else if(offshell->iSpin() == PDT::Spin1) {
		VectorWaveFunction interV = vector_[ix].second->
		  evaluate(q2, 1, offshell,rs_[ofhel1],sbar_[2][ofhel2]);
		diag = -vector_[ix].first->
		  evaluate(q2, spin_[0][ifhel1], sbar_[1][ifhel2], interV);
	      }
	      else
		assert(false);
	    }
	    else
	      assert(false);
 	    me[ix] += norm(diag);
 	    diagramME()[ix](ifhel1, ifhel2, ofhel1, ofhel2) = diag;
 	    //Compute flows
 	    for(size_t iy = 0; iy < current.colourFlow.size(); ++iy) {
 	      assert(current.colourFlow[iy].first<flows.size());
 	      flows[current.colourFlow[iy].first] += 
 		current.colourFlow[iy].second * diag;
 	    }
 	  }
 	  // MEs for the different colour flows
 	  for(unsigned int iy = 0; iy < numberOfFlows(); ++iy) 
 	    flowME()[iy](ifhel1, ifhel2, ofhel1, ofhel2) = flows[iy];
 	  //Now add flows to me2 with appropriate colour factors
 	  for(size_t ii = 0; ii < numberOfFlows(); ++ii) {
 	    for(size_t ij = 0; ij < numberOfFlows(); ++ij) {
 	      me2 += getColourFactors()[ii][ij]*(flows[ii]*conj(flows[ij])).real();
 	    }
 	  }
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
MEff2rf:: ff2rfHeME(double & me2, bool first) const {
  const Energy2 q2(scale());
  // weights for the selection of the diagram
  vector<double> me(numberOfDiags(), 0.);
  // weights for the selection of the colour flow
  vector<double> flow(numberOfFlows(),0.);
  bool massless = mePartonData()[2]->mass()==ZERO;
  // flow over the helicities and diagrams
  for(unsigned int ifhel1 = 0; ifhel1 < 2; ++ifhel1) {
    for(unsigned int ifhel2 = 0; ifhel2 < 2; ++ifhel2) {
      for(unsigned int ofhel1 = 0; ofhel1 < 4; ++ofhel1) {
	if(massless && (ofhel1==1 || ofhel1==2)) continue;
 	for(unsigned int ofhel2 = 0; ofhel2 < 2; ++ofhel2) {
	  vector<Complex> flows(numberOfFlows(),0.);
	  for(HPCount ix = 0; ix < numberOfDiags(); ++ix) {
	    Complex diag(0.);
	    const HPDiagram & current = getProcessInfo()[ix];
	    tcPDPtr offshell = current.intermediate;
 	    if(current.channelType == HPDiagram::tChannel) {
 	      if(offshell->iSpin() == PDT::Spin0) {
		if(current.ordered.second) {
		  ScalarWaveFunction interS = scalar_[ix].first->
		    evaluate(q2, 3, offshell,spin_[1][ifhel2],sbar_[2][ofhel2]);
		  diag = scalar_[ix].second->
		    evaluate(q2, spin_[0][ifhel1], rsbar_[ofhel1], interS);
		}
		else {
		  ScalarWaveFunction interS = scalar_[ix].second->
		    evaluate(q2, 3, offshell,spin_[1][ifhel2],rsbar_[ofhel1]);
		  diag = -scalar_[ix].first->
		    evaluate(q2, spin_[0][ifhel1], sbar_[2][ofhel2], interS);
		}
 	      }
 	      else if(offshell->iSpin() == PDT::Spin1) {
		if(current.ordered.second) {
		  VectorWaveFunction interV = vector_[ix].first->
		    evaluate(q2, 3, offshell,spin_[1][ifhel2],sbar_[2][ofhel2]);
		  diag = vector_[ix].second->
		    evaluate(q2, spin_[0][ifhel1], rsbar_[ofhel1], interV);
		}
		else {
		  VectorWaveFunction interV = vector_[ix].second->
		    evaluate(q2, 3, offshell,spin_[1][ifhel2],rsbar_[ofhel1]);
		  diag = -vector_[ix].first->
		    evaluate(q2, spin_[0][ifhel1], sbar_[2][ofhel2], interV);
		}
	      }
	      else
		assert(false);
	    }
	    else if(current.channelType == HPDiagram::sChannel) {
 	      if(offshell->iSpin() == PDT::Spin0) {
		ScalarWaveFunction interS =
		  scalar_[ix].first->evaluate(q2, 1, offshell, spin_[0][ifhel1], sbar_[1][ifhel2]);
		diag = scalar_[ix].second->evaluate(q2, spin_[2][ofhel2],rsbar_[ofhel1],interS);
 	      }
 	      else if(offshell->iSpin() == PDT::Spin1) {
		VectorWaveFunction interV = vector_[ix].first->
		  evaluate(q2, 1, offshell,spin_[0][ifhel1], sbar_[1][ifhel2]);
		diag = vector_[ix].second->evaluate(q2,spin_[2][ofhel2],rsbar_[ofhel1],interV);
 	      }
	      else
		assert(false);
	    }
	    else
	      assert(false);
 	    me[ix] += norm(diag);
 	    diagramME()[ix](ifhel1, ifhel2, ofhel1, ofhel2) = diag;
 	    //Compute flows
 	    for(size_t iy = 0; iy < current.colourFlow.size(); ++iy) {
 	      assert(current.colourFlow[iy].first<flows.size());
 	      flows[current.colourFlow[iy].first] += 
 		current.colourFlow[iy].second * diag;
 	    }
	  }
	  // MEs for the different colour flows
	  for(unsigned int iy = 0; iy < numberOfFlows(); ++iy) 
	    flowME()[iy](ifhel1, ifhel2, ofhel1, ofhel2) = flows[iy];
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
MEff2rf::fbfb2rbfbHeME(double & me2, bool first) const {
   const Energy2 q2(scale());
   // weights for the selection of the diagram
   vector<double> me(numberOfDiags(), 0.);
   // weights for the selection of the colour flow
   vector<double> flow(numberOfFlows(),0.);
   // flow over the helicities and diagrams
  bool massless = mePartonData()[2]->mass()==ZERO;
  // flow over the helicities and diagrams
  for(unsigned int ifhel1 = 0; ifhel1 < 2; ++ifhel1) {
    for(unsigned int ifhel2 = 0; ifhel2 < 2; ++ifhel2) {
      for(unsigned int ofhel1 = 0; ofhel1 < 4; ++ofhel1) {
	if(massless && (ofhel1==1 || ofhel1==2)) continue;
 	for(unsigned int ofhel2 = 0; ofhel2 < 2; ++ofhel2) {
 	  vector<Complex> flows(numberOfFlows(),0.);
 	  for(HPCount ix = 0; ix < numberOfDiags(); ++ix) {
 	    Complex diag(0.);
 	    const HPDiagram & current = getProcessInfo()[ix];
 	    tcPDPtr offshell = current.intermediate;
	    if(current.channelType == HPDiagram::tChannel) {
	      if(offshell->iSpin() == PDT::Spin0) {
		if(current.ordered.second) {
		  ScalarWaveFunction interS = scalar_[ix].first->
		    evaluate(q2, 3, offshell,spin_[2][ofhel2],sbar_[1][ifhel2]);
		  diag = scalar_[ix].second->
		    evaluate(q2, rs_[ofhel1], sbar_[0][ifhel1], interS);
		}
		else {
		  ScalarWaveFunction interS = scalar_[ix].second->
		    evaluate(q2, 3, offshell,rs_[ofhel1],sbar_[1][ifhel2]);
		  diag = -scalar_[ix].first->
		    evaluate(q2, spin_[2][ofhel2], sbar_[0][ifhel1], interS);
		}
	      }
	      else if(offshell->iSpin() == PDT::Spin1) {
		if(current.ordered.second) {
		  VectorWaveFunction interV = vector_[ix].first->
		    evaluate(q2, 3, offshell,spin_[2][ofhel2],sbar_[1][ifhel2]);
		  diag = vector_[ix].second->
		    evaluate(q2, rs_[ofhel1], sbar_[0][ifhel1], interV);
		}
		else {
		  VectorWaveFunction interV = vector_[ix].second->
		    evaluate(q2, 3, offshell,rs_[ofhel1],sbar_[1][ifhel2]);
		  diag = -vector_[ix].first->
		    evaluate(q2, spin_[2][ofhel2], sbar_[0][ifhel1], interV);
		}
	      }
	      else
		assert(false);
	    }
	    else if(current.channelType == HPDiagram::sChannel) {
	      if(offshell->CC()) offshell=offshell->CC();
	      if(offshell->iSpin() == PDT::Spin0) {
		ScalarWaveFunction interS = scalar_[ix].second->
		  evaluate(q2, 1, offshell,spin_[2][ofhel2],rsbar_[ofhel1]);
		diag = scalar_[ix].first->
		  evaluate(q2, spin_[0][ifhel1], sbar_[1][ifhel2], interS);
	      }
	      else if(offshell->iSpin() == PDT::Spin1) {
		VectorWaveFunction interV = vector_[ix].second->
		  evaluate(q2, 1, offshell,spin_[2][ofhel2],rsbar_[ofhel1]);
		diag = vector_[ix].first->
		  evaluate(q2, spin_[0][ifhel1], sbar_[1][ifhel2], interV);
	      }
	      else
		assert(false);
	    }
	    else
	      assert(false);
	    me[ix] += norm(diag);
	    diagramME()[ix](ifhel1, ifhel2, ofhel1, ofhel2) = diag;
	    //Compute flows
	    for(size_t iy = 0; iy < current.colourFlow.size(); ++iy) {
	      assert(current.colourFlow[iy].first<flows.size());
	      flows[current.colourFlow[iy].first] += 
		current.colourFlow[iy].second * diag;
	    }
	  }
	  // MEs for the different colour flows
	  for(unsigned int iy = 0; iy < numberOfFlows(); ++iy) 
	    flowME()[iy](ifhel1, ifhel2, ofhel1, ofhel2) = flows[iy];
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

void MEff2rf::constructVertex(tSubProPtr subp) {
  // Hard process external particles
  ParticleVector hardpro = hardParticles(subp);
  //Need to use rescale momenta to calculate matrix element
  setRescaledMomenta(hardpro);
  // calculate the fermion spinors
  for(unsigned int ix=0;ix<3;++ix) {
    spin_[ix].clear();
    sbar_[ix].clear();
    unsigned int iy(ix);
    if(ix==2) ++iy;
    SpinorWaveFunction   (spin_[ix],hardpro[iy],
			  ix<2 ? incoming : outgoing,ix>1);
    SpinorBarWaveFunction(sbar_[ix],hardpro[iy],
			  ix<2 ? incoming : outgoing,ix>1);
  }
  rs_   .clear();
  rsbar_.clear();
  // calculate the rs spinors
  RSSpinorWaveFunction   (rs_   , hardpro[2], outgoing, true);
  RSSpinorBarWaveFunction(rsbar_, hardpro[2], outgoing, true);
  double dummy(0.);
  if( mePartonData()[0]->id() > 0 && mePartonData()[1]->id() < 0) {
    if(mePartonData()[2]->id()>0) {
      ProductionMatrixElement prodME = ffb2rfbHeME (dummy,true);
      createVertex(prodME,hardpro);
    }
    else {
      ProductionMatrixElement prodME = ffb2rbfHeME (dummy,true);
      createVertex(prodME,hardpro);
    }
  }
  else if( mePartonData()[0]->id() > 0 && mePartonData()[1]->id() > 0 ) {
    ProductionMatrixElement prodME = ff2rfHeME(dummy,true);
    createVertex(prodME,hardpro);
  }
  else if( mePartonData()[0]->id() < 0 && mePartonData()[1]->id() < 0 ) {
    ProductionMatrixElement prodME = fbfb2rbfbHeME(dummy,true);
    createVertex(prodME,hardpro);
  }
  else
    assert(false);
}
