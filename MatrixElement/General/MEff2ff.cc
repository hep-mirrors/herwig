// -*- C++ -*-
//
// MEff2ff.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEff2ff class.
//

#include "MEff2ff.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using ThePEG::Helicity::VectorWaveFunction;
using ThePEG::Helicity::ScalarWaveFunction;
using ThePEG::Helicity::TensorWaveFunction;
using ThePEG::Helicity::incoming;
using ThePEG::Helicity::outgoing;


void MEff2ff::doinit() {
  GeneralHardME::doinit();
  scalar_.resize(numberOfDiags());
  vector_.resize(numberOfDiags());
  tensor_.resize(numberOfDiags());
  initializeMatrixElements(PDT::Spin1Half, PDT::Spin1Half, 
			   PDT::Spin1Half, PDT::Spin1Half);
  for(size_t ix = 0;ix < numberOfDiags(); ++ix) {
    const HPDiagram & current = getProcessInfo()[ix];
    tcPDPtr offshell = current.intermediate;
    if(offshell->iSpin() == PDT::Spin0) {
      AbstractFFSVertexPtr vert1 = dynamic_ptr_cast<AbstractFFSVertexPtr>
	(current.vertices.first);
      AbstractFFSVertexPtr vert2 = dynamic_ptr_cast<AbstractFFSVertexPtr>
	(current.vertices.second);
      scalar_[ix] = make_pair(vert1, vert2);
    }
    else if(offshell->iSpin() == PDT::Spin1) {
      AbstractFFVVertexPtr vert1 = dynamic_ptr_cast<AbstractFFVVertexPtr>
	(current.vertices.first);
      AbstractFFVVertexPtr vert2 = dynamic_ptr_cast<AbstractFFVVertexPtr>
	(current.vertices.second);
      vector_[ix] = make_pair(vert1, vert2);
    }
    else if(offshell->iSpin() == PDT::Spin2) {
      AbstractFFTVertexPtr vert1 = dynamic_ptr_cast<AbstractFFTVertexPtr>
	(current.vertices.first);
      AbstractFFTVertexPtr vert2 = dynamic_ptr_cast<AbstractFFTVertexPtr>
	(current.vertices.second);
      tensor_[ix] = make_pair(vert1, vert2);
    }
  }
}

double MEff2ff::me2() const {
  tcPDPtr ina (mePartonData()[0]), inb (mePartonData()[1]);
  tcPDPtr outa(mePartonData()[2]), outb(mePartonData()[3]);
  for(unsigned int ix=0;ix<4;++ix) {
    spin_[ix].clear();
    sbar_[ix].clear();
    for(unsigned int ih=0;ih<2;++ih) {
      spin_[ix].push_back(SpinorWaveFunction   (rescaledMomenta()[ix],
						mePartonData()[ix],
						ih, ix<2 ? incoming : outgoing));
      sbar_[ix].push_back(SpinorBarWaveFunction(rescaledMomenta()[ix],
						mePartonData()[ix],
						ih, ix<2 ? incoming : outgoing));
    }
  }
  double full_me(0.);
  if( ina->id() > 0 && inb->id() < 0) {
    ffb2ffbHeME (full_me,true);
  }
  else if( ina->id() > 0 && inb->id() > 0 )
    ff2ffHeME(full_me,true);
  else if( ina->id() < 0 && inb->id() < 0 )
    fbfb2fbfbHeME(full_me,true);
  else 
    throw MEException() 
      << "MEff2ff::me2() - Cannot find correct function to deal with process " 
      << ina->PDGName() << "," << inb->PDGName() << "->" << outa->PDGName() 
      << "," << outb->PDGName() << "\n";

#ifndef NDEBUG
  if( debugME() ) debug(full_me);
#endif

  return full_me;
}

ProductionMatrixElement 
MEff2ff::ffb2ffbHeME(double & me2, bool first) const {
  const Energy2 q2(scale());
  // weights for the selection of the diagram
  vector<double> me(numberOfDiags(), 0.);
  // weights for the selection of the colour flow
  vector<double> flow(numberOfFlows(),0.);
  // flow over the helicities and diagrams
  for(unsigned int ifhel1 = 0; ifhel1 < 2; ++ifhel1) {
    for(unsigned int ifhel2 = 0; ifhel2 < 2; ++ifhel2) {
      for(unsigned int ofhel1 = 0; ofhel1 < 2; ++ofhel1) {
	for(unsigned int ofhel2 = 0; ofhel2 < 2; ++ofhel2) {
	  vector<Complex> flows(numberOfFlows(),0.);
	  for(HPCount ix = 0; ix < numberOfDiags(); ++ix) {
	    Complex diag(0.);
	    const HPDiagram & current = getProcessInfo()[ix];
	    tcPDPtr offshell = current.intermediate;
	    if(current.channelType == HPDiagram::tChannel) {
	      if(offshell->iSpin() == PDT::Spin0) {
		if(current.ordered.second) {
		  ScalarWaveFunction interS = scalar_[ix].second->
		    evaluate(q2, 3, offshell,spin_[3][ofhel2],sbar_[1][ifhel2]);
		  diag = -scalar_[ix].first->
		    evaluate(q2, spin_[0][ifhel1],sbar_[2][ofhel1],interS);
		}
		else {
		  ScalarWaveFunction interS = scalar_[ix].second->
		    evaluate(q2, 3, offshell,spin_[2][ofhel1],sbar_[1][ifhel2]);
		  diag = scalar_[ix].first->
		    evaluate(q2, spin_[0][ifhel1],sbar_[3][ofhel2],interS);
		}
	      }
	      else if(offshell->iSpin() == PDT::Spin1) {
		if(current.ordered.second) {
		  VectorWaveFunction interV = vector_[ix].second->
		    evaluate(q2, 3, offshell,spin_[3][ofhel2],sbar_[1][ifhel2]);
		  diag = -vector_[ix].first->
		    evaluate(q2, spin_[0][ifhel1], sbar_[2][ofhel1], interV);
		}
		else {
		  VectorWaveFunction interV = vector_[ix].second->
		    evaluate(q2, 3, offshell,spin_[2][ofhel1],sbar_[1][ifhel2]);
		  diag = vector_[ix].first->
		    evaluate(q2, spin_[0][ifhel1], sbar_[3][ofhel2], interV);
 		}
	      }
	      else if(offshell->iSpin() == PDT::Spin2) {
		if(current.ordered.second) {
		  TensorWaveFunction interT = tensor_[ix].second->
		    evaluate(q2, 3, offshell,spin_[3][ofhel2],sbar_[1][ifhel2]);
		  diag = -tensor_[ix].first->
		    evaluate(q2, spin_[0][ifhel1], sbar_[2][ofhel1], interT);
		}
		else {
		  TensorWaveFunction interT = tensor_[ix].second->
		    evaluate(q2, 3, offshell,spin_[2][ofhel1],sbar_[1][ifhel2]);
		  diag = tensor_[ix].first->
		    evaluate(q2, spin_[0][ifhel1], sbar_[3][ofhel2], interT);
 		}
	      }
	    }
	    else if(current.channelType == HPDiagram::sChannel) {
	      if(offshell->iSpin() == PDT::Spin0) {
		ScalarWaveFunction interS = scalar_[ix].second->
		  evaluate(q2, 1, offshell,spin_[3][ofhel2],sbar_[2][ofhel1]);
		diag = scalar_[ix].first->
		  evaluate(q2, spin_[0][ifhel1], sbar_[1][ifhel2], interS);
	      }
	      else if(offshell->iSpin() == PDT::Spin1) {
		VectorWaveFunction interV = vector_[ix].second->
		  evaluate(q2, 1, offshell,spin_[3][ofhel2],sbar_[2][ofhel1]);
		diag = vector_[ix].first->
		  evaluate(q2, spin_[0][ifhel1], sbar_[1][ifhel2], interV);
	      }
	      else if(offshell->iSpin() == PDT::Spin2) {
		TensorWaveFunction interT = tensor_[ix].second->
		  evaluate(q2, 1, offshell,spin_[3][ofhel2],sbar_[2][ofhel1]);
		diag = tensor_[ix].first->
		  evaluate(q2, spin_[0][ifhel1], sbar_[1][ifhel2], interT);
	      }
	    }
	    else assert(false);
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
MEff2ff:: ff2ffHeME(double & me2, bool first) const {
  const Energy2 q2(scale());
  // weights for the selection of the diagram
  vector<double> me(numberOfDiags(), 0.);
  // weights for the selection of the colour flow
  vector<double> flow(numberOfFlows(),0.);
  // flow over the helicities and diagrams
  for(unsigned int ifhel1 = 0; ifhel1 < 2; ++ifhel1) {
    for(unsigned int ifhel2 = 0; ifhel2 < 2; ++ifhel2) {
      for(unsigned int ofhel1 = 0; ofhel1 < 2; ++ofhel1) {
	for(unsigned int ofhel2 = 0; ofhel2 < 2; ++ofhel2) {
	  vector<Complex> flows(numberOfFlows(),0.);
	  for(HPCount ix = 0; ix < numberOfDiags(); ++ix) {
	    Complex diag(0.);
	    const HPDiagram & current = getProcessInfo()[ix];
	    tcPDPtr offshell = current.intermediate;
	    if(current.channelType == HPDiagram::tChannel) {
	      if(offshell->iSpin() == PDT::Spin0) {
		if(current.ordered.second) {
		  ScalarWaveFunction interS = scalar_[ix].second->
		    evaluate(q2, 3, offshell,spin_[1][ifhel2],sbar_[3][ofhel2]);
		  diag = scalar_[ix].first->
		    evaluate(q2, spin_[0][ifhel1], sbar_[2][ofhel1], interS);
		}
		else {
		  ScalarWaveFunction interS = scalar_[ix].second->
		    evaluate(q2, 3, offshell,spin_[1][ifhel2],sbar_[2][ofhel1]);
		  diag = -scalar_[ix].first->
		    evaluate(q2, spin_[0][ifhel1], sbar_[3][ofhel2], interS);
		}
	      }
	      else if(offshell->iSpin() == PDT::Spin1) {
		if(current.ordered.second) {
		  VectorWaveFunction interV = vector_[ix].second->
		    evaluate(q2, 3, offshell,spin_[1][ifhel2],sbar_[3][ofhel2]);
		  diag = vector_[ix].first->
		    evaluate(q2, spin_[0][ifhel1], sbar_[2][ofhel1], interV);
		}
		else {
		  VectorWaveFunction interV = vector_[ix].second->
		    evaluate(q2, 3, offshell,spin_[1][ifhel2],sbar_[2][ofhel1]);
		  diag = -vector_[ix].first->
		    evaluate(q2, spin_[0][ifhel1], sbar_[3][ofhel2], interV);
		}
	      }
	      else if(offshell->iSpin() == PDT::Spin2) {
		if(current.ordered.second) {
		  TensorWaveFunction interT = tensor_[ix].second->
		    evaluate(q2, 3, offshell,spin_[1][ifhel2],sbar_[3][ofhel2]);
		  diag = tensor_[ix].first->
		    evaluate(q2, spin_[0][ifhel1], sbar_[2][ofhel1], interT);
		}
		else {
		  TensorWaveFunction interT = tensor_[ix].second->
		    evaluate(q2, 3, offshell,spin_[1][ifhel2],sbar_[2][ofhel1]);
		  diag = -tensor_[ix].first->
		    evaluate(q2, spin_[0][ifhel1], sbar_[3][ofhel2], interT);
		}
	      }
	    }
	    else if(current.channelType == HPDiagram::sChannel) {
	      if(offshell->CC()) offshell=offshell->CC();
	      if(offshell->iSpin() == PDT::Spin0) {
		ScalarWaveFunction interS = scalar_[ix].second->
		  evaluate(q2, 1, offshell,spin_[3][ofhel2],sbar_[2][ofhel1]);
		diag = scalar_[ix].first->
		  evaluate(q2, spin_[0][ifhel1], sbar_[1][ifhel2], interS);
	      }
	      else if(offshell->iSpin() == PDT::Spin1) {
		VectorWaveFunction interV = vector_[ix].second->
		  evaluate(q2, 1, offshell,spin_[3][ofhel2],sbar_[2][ofhel1]);
		diag = vector_[ix].first->
		  evaluate(q2, spin_[0][ifhel1], sbar_[1][ifhel2], interV);
	      }
	      else if(offshell->iSpin() == PDT::Spin2) {
		TensorWaveFunction interT = tensor_[ix].second->
		  evaluate(q2, 1, offshell,spin_[3][ofhel2],sbar_[2][ofhel1]);
		diag = tensor_[ix].first->
		  evaluate(q2, spin_[0][ifhel1], sbar_[1][ifhel2], interT);
	      }
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
MEff2ff::fbfb2fbfbHeME(double & me2, bool first) const {
  const Energy2 q2(scale());
  // weights for the selection of the diagram
  vector<double> me(numberOfDiags(), 0.);
  // weights for the selection of the colour flow
  vector<double> flow(numberOfFlows(),0.);
  // flow over the helicities and diagrams
  for(unsigned int ifhel1 = 0; ifhel1 < 2; ++ifhel1) {
    for(unsigned int ifhel2 = 0; ifhel2 < 2; ++ifhel2) {
      for(unsigned int ofhel1 = 0; ofhel1 < 2; ++ofhel1) {
	for(unsigned int ofhel2 = 0; ofhel2 < 2; ++ofhel2) {
	  vector<Complex> flows(numberOfFlows(),0.);
	  for(HPCount ix = 0; ix < numberOfDiags(); ++ix) {
	    Complex diag(0.);
	    const HPDiagram & current = getProcessInfo()[ix];
	    tcPDPtr offshell = current.intermediate;
	    if(current.channelType == HPDiagram::tChannel) {
	      if(offshell->iSpin() == PDT::Spin0) {
		if(current.ordered.second) {
		  ScalarWaveFunction interS = scalar_[ix].second->
		    evaluate(q2, 3, offshell,spin_[3][ofhel2],sbar_[1][ifhel2]);
		  diag = scalar_[ix].first->
		    evaluate(q2, spin_[2][ofhel1], sbar_[0][ifhel1], interS);
		}
		else {
		  ScalarWaveFunction interS = scalar_[ix].second->
		    evaluate(q2, 3, offshell,spin_[2][ofhel1],sbar_[1][ifhel2]);
		  diag = -scalar_[ix].first->
		    evaluate(q2, spin_[3][ofhel2], sbar_[0][ifhel1], interS);
		}
	      }
	      else if(offshell->iSpin() == PDT::Spin1) {
		if(current.ordered.second) {
		  VectorWaveFunction interV = vector_[ix].second->
		    evaluate(q2, 3, offshell,spin_[3][ofhel2],sbar_[1][ifhel2]);
		  diag = vector_[ix].first->
		    evaluate(q2, spin_[2][ofhel1], sbar_[0][ifhel1], interV);
		}
		else {
		  VectorWaveFunction interV = vector_[ix].second->
		    evaluate(q2, 3, offshell,spin_[2][ofhel1],sbar_[1][ifhel2]);
		  diag = -vector_[ix].first->
		    evaluate(q2, spin_[3][ofhel2], sbar_[0][ifhel1], interV);
		}
	      }
	      else if(offshell->iSpin() == PDT::Spin2) {
		if(current.ordered.second) {
		  TensorWaveFunction interT = tensor_[ix].second->
		    evaluate(q2, 3, offshell,spin_[3][ofhel2],sbar_[1][ifhel2]);
		  diag = tensor_[ix].first->
		    evaluate(q2, spin_[2][ofhel1], sbar_[0][ifhel1], interT);
		}
		else {
		  TensorWaveFunction interT = tensor_[ix].second->
		    evaluate(q2, 3, offshell,spin_[2][ofhel1],sbar_[1][ifhel2]);
		  diag = -tensor_[ix].first->
		    evaluate(q2, spin_[3][ofhel2], sbar_[0][ifhel1], interT);
		}
	      }
	    }
	    else if(current.channelType == HPDiagram::sChannel) {
	      if(offshell->CC()) offshell=offshell->CC();
	      if(offshell->iSpin() == PDT::Spin0) {
		ScalarWaveFunction interS = scalar_[ix].second->
		  evaluate(q2, 1, offshell,spin_[3][ofhel2],sbar_[2][ofhel1]);
		diag = scalar_[ix].first->
		  evaluate(q2, spin_[0][ifhel1], sbar_[1][ifhel2], interS);
	      }
	      else if(offshell->iSpin() == PDT::Spin1) {
		VectorWaveFunction interV = vector_[ix].second->
		  evaluate(q2, 1, offshell,spin_[3][ofhel2],sbar_[2][ofhel1]);
		diag = vector_[ix].first->
		  evaluate(q2, spin_[0][ifhel1], sbar_[1][ifhel2], interV);
	      }
	      else if(offshell->iSpin() == PDT::Spin2) {
		TensorWaveFunction interT = tensor_[ix].second->
		  evaluate(q2, 1, offshell,spin_[3][ofhel2],sbar_[2][ofhel1]);
		diag = tensor_[ix].first->
		  evaluate(q2, spin_[0][ifhel1], sbar_[1][ifhel2], interT);
	      }
	    }
	    else {
	      assert(false);
	    }
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

void MEff2ff::constructVertex(tSubProPtr subp) {
  // Hard process external particles
  ParticleVector hardpro = hardParticles(subp);
  //Need to use rescale momenta to calculate matrix element
  setRescaledMomenta(hardpro);
  for(unsigned int ix=0;ix<4;++ix) {
    spin_[ix].clear();
    sbar_[ix].clear();
    for(unsigned int ih=0;ih<2;++ih) {
      SpinorWaveFunction   (spin_[ix],hardpro[ix],
			    ix<2 ? incoming : outgoing,ix>1);
      SpinorBarWaveFunction(sbar_[ix],hardpro[ix],
			    ix<2 ? incoming : outgoing,ix>1);
    }
  }
  double dummy(0.);
  //pick which process we are doing
  if( hardpro[0]->id() > 0) {
    //ffbar->ffbar
    if( hardpro[1]->id() < 0 ) {
      ProductionMatrixElement prodME = ffb2ffbHeME(dummy,false);
      createVertex(prodME,hardpro);
    }
    //ff2ff
    else {
      ProductionMatrixElement prodME = ff2ffHeME(dummy,false);
      createVertex(prodME,hardpro);
    }
  } 
  //fbarfbar->fbarfbar
  else {
    ProductionMatrixElement prodME = fbfb2fbfbHeME(dummy,false);
    createVertex(prodME,hardpro);
  }
  
#ifndef NDEBUG
  if( debugME() ) debug(dummy);
#endif

}

void MEff2ff::persistentOutput(PersistentOStream & os) const {
  os << scalar_ << vector_ << tensor_;
}

void MEff2ff::persistentInput(PersistentIStream & is, int) {
  is >> scalar_ >> vector_ >> tensor_;
  initializeMatrixElements(PDT::Spin1Half, PDT::Spin1Half, 
			   PDT::Spin1Half, PDT::Spin1Half);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEff2ff,GeneralHardME>
describeHerwigMEff2ff("Herwig::MEff2ff", "Herwig.so");

void MEff2ff::Init() {

  static ClassDocumentation<MEff2ff> documentation
    ("This is the implementation of the matrix element for fermion-"
     "antifermion -> fermion-antifermion.");

}

void MEff2ff::debug(double me2) const {
  if( !generator()->logfile().is_open() ) return;
  long id1 = mePartonData()[0]->id();
  long id2 = mePartonData()[1]->id();
  long id3 = mePartonData()[2]->id();
  long id4 = mePartonData()[3]->id();
  long aid1 = abs(mePartonData()[0]->id());
  long aid2 = abs(mePartonData()[1]->id());
  long aid3 = abs(mePartonData()[2]->id());
  long aid4 = abs(mePartonData()[3]->id());
  if( (aid1 != 1 && aid1 != 2) || (aid2 != 1 && aid2 != 2) ) return;
  double analytic(0.);
  if( id3 == id4 && id3 == 1000021 ) {
    tcSMPtr sm = generator()->standardModel();
    double gs4 = sqr( 4.*Constants::pi*sm->alphaS(scale()) );
    int Nc = sm->Nc();
    double Cf = (sqr(Nc) - 1.)/2./Nc;
    Energy2 mgo2 = meMomenta()[3].m2();
    long squark = (aid1 == 1) ? 1000001 : 1000002;
    Energy2 muL2 = sqr(getParticleData(squark)->mass());
    Energy2 deltaL = muL2 - mgo2;
    Energy2 muR2 = sqr(getParticleData(squark + 1000000)->mass());
    Energy2 deltaR = muR2 - mgo2;
    Energy2 s(sHat());
    Energy2 m3s = meMomenta()[2].m2();
    Energy2 m4s = meMomenta()[3].m2();
    Energy4 spt2 = uHat()*tHat() - m3s*m4s;
    Energy2 t3(tHat() - m3s), u4(uHat() - m4s);
    
    double Cl = 2.*spt2*( (u4*u4 - deltaL*deltaL) + (t3*t3 - deltaL*deltaL)
			  - (s*s/Nc/Nc) )/s/s/(u4 - deltaL)/(t3 - deltaL);
    Cl += deltaL*deltaL*( (1./sqr(t3 - deltaL)) + (1./sqr(u4 - deltaL))
			  - ( sqr( (1./(t3 - deltaL)) - 
				   (1./(u4 - deltaL)) )/Nc/Nc ) );
    
    double Cr = 2.*spt2*( (u4*u4 - deltaR*deltaR) + (t3*t3 - deltaR*deltaR)
			  - (s*s/Nc/Nc) )/s/s/(u4 - deltaR)/(t3 - deltaR);
    Cr += deltaR*deltaR*( (1./sqr(t3 - deltaR)) + (1./sqr(u4 - deltaR))
			  - ( sqr( (1./(t3 - deltaR)) 
				   - (1./(u4 - deltaR)) )/Nc/Nc ) );
    analytic = gs4*Cf*(Cl + Cr)/4.;
  }
  else if( (aid3 == 5100001 || aid3 == 5100002 ||
	    aid3 == 6100001 || aid3 == 6100002) &&
	   (aid4 == 5100001 || aid4 == 5100002 ||
	    aid4 == 6100001 || aid4 == 6100002) ) {
    tcSMPtr sm = generator()->standardModel();
    double gs4 = sqr( 4.*Constants::pi*sm->alphaS(scale()) );
    Energy2 s(sHat());
    Energy2 mf2 = meMomenta()[2].m2();
    Energy2 t3(tHat() - mf2), u4(uHat() - mf2);
    Energy4 s2(sqr(s)), t3s(sqr(t3)), u4s(sqr(u4));
    
    bool iflav = (aid2 - aid1 == 0);
    int alpha(aid3/1000000), beta(aid4/1000000);
    bool oflav = ((aid3 - aid1) % 10  == 0);
    if( alpha != beta ) {
      if( ( id1 > 0 && id2 > 0) ||
	  ( id1 < 0 && id2 < 0) ) {
	if( iflav )
	  analytic = gs4*( mf2*(2.*s2*s/t3s/u4s - 4.*s/t3/u4) 
			   + 2.*sqr(s2)/t3s/u4s - 8.*s2/t3/u4 + 5. )/9.;
	else
	  analytic = gs4*( -2.*mf2*(1./t3 + u4/t3s) + 0.5 + 2.*u4s/t3s)/9.;
      }
      else
	analytic = gs4*( 2.*mf2*(1./t3 + u4/t3s) + 5./2. + 4.*u4/t3 
			 + 2.*u4s/t3s)/9.;
    }
    else {
      if( ( id1 > 0 && id2 > 0) ||
	  ( id1 < 0 && id2 < 0) ) {
	if( iflav ) {
	  analytic = gs4*( mf2*(6.*t3/u4s + 6.*u4/t3s - s/t3/u4) 
			   + 2.*(3.*t3s/u4s + 3.*u4s/t3s 
				 + 4.*s2/t3/u4 - 5.) )/27.;
	}
	else
	  analytic = 2.*gs4*( -mf2*s/t3s + 0.25 + s2/t3s )/9.;
      }
      else {
	if( iflav ) {
	  if( oflav )	  
	    analytic = gs4*( 2.*mf2*(4./s + s/t3s - 1./t3) + 23./6.+ 2.*s2/t3s
			     + 8.*s/3./t3 + 6.*t3/s + 8.*t3s/s2 )/9.;
	  else
	    analytic = 4.*gs4*( 2.*mf2/s + (t3s + u4s)/s2)/9.;
	}
	else 
	  analytic = gs4*(4.*mf2*s/t3s + 5. + 4.*s2/t3s + 8.*s/t3 )/18.;
      }
    }
    if( id3 == id4 ) analytic /= 2.;

  }
  else return;
  double diff = abs(analytic - me2);
  if( diff  > 1e-4 ) {
    generator()->log() 
      << mePartonData()[0]->PDGName() << ","
      << mePartonData()[1]->PDGName() << "->"
      << mePartonData()[2]->PDGName() << ","
      << mePartonData()[3]->PDGName() << "   difference: " 
      << setprecision(10) << diff  << "  ratio: " << analytic/me2 
      << '\n';
  }
}
