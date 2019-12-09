// -*- C++ -*-
//
// MEvv2vv.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEvv2vv class.
//

#include "MEvv2vv.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Models/StandardModel/StandardModel.h"

using namespace Herwig;
using ThePEG::Helicity::ScalarWaveFunction;
using ThePEG::Helicity::TensorWaveFunction;
using ThePEG::Helicity::incoming;
using ThePEG::Helicity::outgoing;


void MEvv2vv::doinit() {
  GeneralHardME::doinit();
  scalar_.resize(numberOfDiags());
  vector_.resize(numberOfDiags());
  tensor_.resize(numberOfDiags());
  four_  .resize(numberOfDiags());
  initializeMatrixElements(PDT::Spin1, PDT::Spin1,
			   PDT::Spin1, PDT::Spin1);
  for(size_t i = 0; i < numberOfDiags(); ++i) {
    HPDiagram diag = getProcessInfo()[i];
    tcPDPtr offshell = diag.intermediate;
    if(!offshell)
      four_[i] = dynamic_ptr_cast<AbstractVVVVVertexPtr>
	(diag.vertices.first);
    else if(offshell->iSpin() == PDT::Spin0) {
      AbstractVVSVertexPtr vert1 = dynamic_ptr_cast<AbstractVVSVertexPtr>
	(diag.vertices.first);
      AbstractVVSVertexPtr vert2 = dynamic_ptr_cast<AbstractVVSVertexPtr>
	(diag.vertices.second);
      scalar_[i] = make_pair(vert1, vert2);
    }
    else if(offshell->iSpin() == PDT::Spin1) {
      AbstractVVVVertexPtr vert1 = dynamic_ptr_cast<AbstractVVVVertexPtr>
	(diag.vertices.first);
      AbstractVVVVertexPtr vert2 = dynamic_ptr_cast<AbstractVVVVertexPtr>
	(diag.vertices.second);
      vector_[i] = make_pair(vert1, vert2);
    }
    else if(offshell->iSpin() == PDT::Spin2) {
      AbstractVVTVertexPtr vert1 = dynamic_ptr_cast<AbstractVVTVertexPtr>
	(diag.vertices.first);
      AbstractVVTVertexPtr vert2 = dynamic_ptr_cast<AbstractVVTVertexPtr>
	(diag.vertices.second);
      tensor_[i] = make_pair(vert1, vert2);
    }
  }
  if(colour()==Colour88to88||colour()==Colour88to66bar) {
    tcHwSMPtr hwsm= dynamic_ptr_cast<tcHwSMPtr>(standardModel());
    for(size_t i = 0; i < numberOfDiags(); ++i) {
      HPDiagram diag = getProcessInfo()[i];
      if(diag.intermediate) continue;
      vector<int> order;
      for(map<string,pair<unsigned int,int> >::const_iterator it=hwsm->couplings().begin();
	  it!=hwsm->couplings().end();++it) {
	order.push_back(0);
	if(diag.vertices.first ) order.back() += diag.vertices.first ->orderInCoupling(it->second.first);
	if(diag.vertices.second&&diag.vertices.first->getNpoint()==3)
	  order.back() += diag.vertices.second->orderInCoupling(it->second.first);
      }
      vector<unsigned int> matchdiags;
      for(size_t j = 0; j < numberOfDiags(); ++j) {
	HPDiagram diag2 = getProcessInfo()[j];
	if(!diag2.intermediate ||
	   (diag2.intermediate->iColour()==PDT::Colour8 &&
	    diag2.intermediate->iColour()==PDT::Colour6 &&
	    diag2.intermediate->iColour()==PDT::Colour6bar)) continue;
	unsigned int iloc(0);
	bool match=true;
	for(map<string,pair<unsigned int,int> >::const_iterator it=hwsm->couplings().begin();
	    it!=hwsm->couplings().end();++it) {
	  int otemp(0);
	  if(diag2.vertices.first ) otemp += diag2.vertices.first ->orderInCoupling(it->second.first);
	  if(diag2.vertices.second&&diag2.vertices.first->getNpoint()==3)
	    otemp += diag2.vertices.second->orderInCoupling(it->second.first);
	  if(otemp!=order[iloc]) {
	    match = false;
	    break;
	  }
	  iloc+=1;
	}
	if(!match) continue;
	matchdiags.push_back(j);
      }
      double weight = 3./double(matchdiags.size());
      for(unsigned int iy=0;iy<matchdiags.size();++iy)
	if(fourFlow_.find(matchdiags[iy])!=fourFlow_.end())
	  fourFlow_[matchdiags[iy]].push_back(make_pair(i,weight));
	else
	  fourFlow_[matchdiags[iy]] = vector<pair<unsigned int,double> >(1,make_pair(i,weight));
    }
  }
}

double MEvv2vv::me2() const {
  VBVector va(2), vb(2), vc(3), vd(3);  
  for(unsigned int i = 0; i < 2; ++i) {
    va[i] = VectorWaveFunction(rescaledMomenta()[0], mePartonData()[0], 2*i, 
			       incoming);
    vb[i] = VectorWaveFunction(rescaledMomenta()[1], mePartonData()[1], 2*i, 
			       incoming);
  }
  //always 0 and 2 polarisations
  for(unsigned int i = 0; i < 2; ++i) {
    vc[2*i] = VectorWaveFunction(rescaledMomenta()[2], mePartonData()[2], 2*i, 
				 outgoing);
    vd[2*i] = VectorWaveFunction(rescaledMomenta()[3], mePartonData()[3], 2*i, 
				 outgoing);
  }
  bool mc  = !(mePartonData()[2]->mass() > ZERO);
  //massive vector, also 1
  if( !mc )
    vc[1] = VectorWaveFunction(rescaledMomenta()[2], mePartonData()[2], 1, 
			       outgoing);
  bool md  = !(mePartonData()[3]->mass() > ZERO);
  if( !md ) 
    vd[1] = VectorWaveFunction(rescaledMomenta()[3], mePartonData()[3], 1, 
			       outgoing);
  double full_me(0.);
  vv2vvHeME(va, vb, vc, mc, vd, md, full_me,true);

#ifndef NDEBUG
  if( debugME() ) debug(full_me);
#endif

  return full_me;
}

ProductionMatrixElement 
MEvv2vv::vv2vvHeME(VBVector & vin1, VBVector & vin2, 
		   VBVector & vout1, bool mc, VBVector & vout2, bool md,
		   double & me2, bool first) const {
  const Energy2 q2(scale());
  const Energy mass = vout1[0].mass();
  // weights for the selection of the diagram
  vector<double> me(numberOfDiags(), 0.);
  // weights for the selection of the colour flow
  vector<double> flow(numberOfFlows(),0.);
  // flow over the helicities and diagrams
  for(unsigned int ihel1 = 0; ihel1 < 2; ++ihel1) { 
    for(unsigned int ihel2 = 0; ihel2 < 2; ++ihel2) {
      for(unsigned int ohel1 = 0; ohel1 < 3; ++ohel1) {
	if(mc && ohel1 == 1) ++ohel1;
	for(unsigned int ohel2 = 0; ohel2 < 3; ++ohel2) {
	  if(md && ohel2 == 1) ++ohel2;
	  vector<Complex> flows(numberOfFlows(),0.);
	  for(HPCount ix = 0; ix < numberOfDiags(); ++ix) {
	    Complex diag(0.);
	    const HPDiagram & current = getProcessInfo()[ix];
	    tcPDPtr offshell = current.intermediate;
	    if(current.channelType == HPDiagram::sChannel) {
	      if(offshell->iSpin() == PDT::Spin0) {
		ScalarWaveFunction interS = 
		  scalar_[ix].first->evaluate(q2, 1, offshell,
					      vin1[ihel1], vin2[ihel2]);
		diag = scalar_[ix].second->
		  evaluate(q2, vout1[ohel1], vout2[ohel2], interS);
	      }
	      else if(offshell->iSpin() == PDT::Spin1) {
		VectorWaveFunction interV = vector_[ix].first->
		  evaluate(q2, 1, offshell, vin1[ihel1], vin2[ihel2]);
		diag = vector_[ix].second->
		  evaluate(q2, vout1[ohel1], vout2[ohel2], interV);
		if(colour()==Colour88to88)
		  for(unsigned int iy=0;iy<fourFlow_.at(ix).size();++iy) {
		    unsigned int iloc=fourFlow_.at(ix)[iy].first;
		    double wgt = fourFlow_.at(ix)[iy].second; 
		    diag += wgt*four_[iloc]->evaluate(q2, 0, vout1[ohel1], vin2[ihel2], 
						      vout2[ohel2], vin1[ihel1]);
		  }
		else if(colour()==Colour88to66bar)
		  for(unsigned int iy=0;iy<fourFlow_.at(ix).size();++iy) {
		    unsigned int iloc=fourFlow_.at(ix)[iy].first;
		    double wgt = fourFlow_.at(ix)[iy].second; 
		    diag -= wgt*four_[iloc]->evaluate(q2, 0, vout1[ohel1], vin2[ihel2], 
						      vout2[ohel2], vin1[ihel1]);
		  }
	      }
	      else if(offshell->iSpin() == PDT::Spin2) {
		TensorWaveFunction interT = tensor_[ix].first->
		  evaluate(q2, 1, offshell, vin1[ihel1], vin2[ihel2]);
		diag = tensor_[ix].second->
		  evaluate(q2, vout1[ohel1], vout2[ohel2],interT);
	      }
	      else 
		assert(false);
	    }
	    else if(current.channelType == HPDiagram::tChannel) {
	      if(offshell->iSpin() == PDT::Spin0) {
		if(current.ordered.second) {
		  ScalarWaveFunction interS = scalar_[ix].
		    first->evaluate(q2, 3, offshell, vin1[ihel1],vout1[ohel1]);
		  diag = scalar_[ix].second->
		    evaluate(q2, vin2[ihel2], vout2[ohel2], interS);
		}
		else {
		  ScalarWaveFunction interS = scalar_[ix].first->
		    evaluate(q2, 3, offshell, vin2[ihel2],vout1[ohel1]);
		  diag = scalar_[ix].second->
		    evaluate(q2, vin1[ihel1], vout2[ohel2], interS);
		}
	      }
	      else if(offshell->iSpin() == PDT::Spin1) {
		if(current.ordered.second) {
		  VectorWaveFunction interV = vector_[ix].
		    first->evaluate(q2, 3, offshell, vin1[ihel1],vout1[ohel1], mass);
		  diag = vector_[ix].second->
		    evaluate(q2, vin2[ihel2], interV, vout2[ohel2]);
		  if(colour()==Colour88to88 || colour()==Colour88to66bar)
		    for(unsigned int iy=0;iy<fourFlow_.at(ix).size();++iy) {
		      unsigned int iloc=fourFlow_.at(ix)[iy].first;
		      double wgt = fourFlow_.at(ix)[iy].second; 
		      diag += wgt*four_[iloc]->evaluate(q2, 0, vin1[ihel1], vin2[ihel2], 
							vout1[ohel1], vout2[ohel2]);
		    }
		}
		else {
		  if(offshell->CC()) offshell = offshell->CC();
		  VectorWaveFunction interV = vector_[ix].first->
		    evaluate(q2, 3, offshell, vin2[ihel2],vout1[ohel1], mass);
		  diag = vector_[ix].second->
		    evaluate(q2, vin1[ihel1], interV, vout2[ohel2]);
		  if(colour()==Colour88to88 || colour()==Colour88to66bar)
		    for(unsigned int iy=0;iy<fourFlow_.at(ix).size();++iy) {
		      unsigned int iloc=fourFlow_.at(ix)[iy].first;
		      double wgt = fourFlow_.at(ix)[iy].second; 
		      diag += wgt*four_[iloc]->
			evaluate(q2, 0, vin2[ihel2], vin1[ihel1],
				 vout1[ohel1], vout2[ohel2]);
		    }
		}
	      }
	      else if(offshell->iSpin() == PDT::Spin2) {
		if(current.ordered.second) {
		  TensorWaveFunction interT = tensor_[ix].first->
		    evaluate(q2, 3, offshell, vin1[ihel1],vout1[ohel1]);
		 diag = tensor_[ix].second->
		   evaluate(q2, vin2[ihel2], vout2[ohel2], interT);
		}
		else {
		  TensorWaveFunction interT = tensor_[ix].first->
		    evaluate(q2, 3, offshell, vin2[ihel2],vout1[ohel1]);
		  diag = tensor_[ix].second->
		    evaluate(q2, vin1[ihel1], vout2[ohel2], interT);
		}
	      }
	      else 
		assert(false);
	    }
	    else if(current.channelType == HPDiagram::fourPoint) {
	      if(colour()==Colour88to88||colour()==Colour88to66bar)
		diag = 0.;
	      else
		diag = four_[ix]->evaluate(q2, 0, vin1[ihel1], vin2[ihel2],
					   vout1[ohel1], vout2[ohel2]);
	    }
	    else
	      assert(false);
	    me[ix] += norm(diag);
	    diagramME()[ix](2*ihel1, 2*ihel2, ohel1, ohel2) = diag;
	    //Compute flows
	    for(size_t iy = 0; iy < current.colourFlow.size(); ++iy) {
	      assert(current.colourFlow[iy].first<flows.size());
	      flows[current.colourFlow[iy].first] += 
		current.colourFlow[iy].second * diag;
	    }
	  }
	  // MEs for the different colour flows
	  for(unsigned int iy = 0; iy < numberOfFlows(); ++iy) 
	    flowME()[iy](2*ihel1, 2*ihel2, ohel1, ohel2) = flows[iy];
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


void MEvv2vv::persistentOutput(PersistentOStream & os) const {
  os << scalar_ << vector_ << tensor_ << four_ << fourFlow_;
}

void MEvv2vv::persistentInput(PersistentIStream & is, int) {
  is >> scalar_ >> vector_ >> tensor_ >> four_ >> fourFlow_;
  initializeMatrixElements(PDT::Spin1, PDT::Spin1,
			   PDT::Spin1, PDT::Spin1);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEvv2vv,GeneralHardME>
describeHerwigMEvv2vv("Herwig::MEvv2vv", "Herwig.so");

void MEvv2vv::Init() {

  static ClassDocumentation<MEvv2vv> documentation
    ("This is the implementation of the 2 to 2 ME for a pair"
     "of massless vector-bosons to a pair of vector bosons");

}

void MEvv2vv::constructVertex(tSubProPtr sub) {
  ParticleVector ext = hardParticles(sub);
  // set wave functions with real momenta
  VBVector v1, v2, v3, v4;
  VectorWaveFunction(v1, ext[0], incoming, false, true);
  VectorWaveFunction(v2, ext[1], incoming, false, true);
  //function to calculate me2 expects massless incoming vectors
  // and this constructor sets the '1' polarisation at element [2] 
  //in the vector
  bool mc  = !(ext[2]->data().mass() > ZERO);
  bool md  = !(ext[3]->data().mass() > ZERO);
  VectorWaveFunction(v3, ext[2], outgoing, true, mc);
  VectorWaveFunction(v4, ext[3], outgoing, true, md);
  // Need to use rescale momenta to calculate matrix element
  setRescaledMomenta(ext);
  // wave functions with rescaled momenta
  VectorWaveFunction vr1(rescaledMomenta()[0],
			 ext[0]->dataPtr(), incoming);
  VectorWaveFunction vr2(rescaledMomenta()[1],
			 ext[1]->dataPtr(), incoming);
  VectorWaveFunction vr3(rescaledMomenta()[2],
			 ext[2]->dataPtr(), outgoing);
  VectorWaveFunction vr4(rescaledMomenta()[3],
			 ext[3]->dataPtr(), outgoing);
  for( unsigned int ihel = 0; ihel < 2; ++ihel ) {
    vr1.reset(2*ihel);
    v1[ihel] = vr1;
    vr2.reset(2*ihel);
    v2[ihel] = vr2;
    vr3.reset(2*ihel);
    v3[2*ihel] = vr3;
    vr4.reset(2*ihel);
    v4[2*ihel] = vr4;
  }
  if( !mc ) {
    vr3.reset(1);
    v3[1] = vr3;
  }
  if( !md ) {
    vr4.reset(1);
    v4[1] = vr4;
  }
  double dummy(0.);
  ProductionMatrixElement pme = vv2vvHeME(v1, v2, v3, mc, v4, md, dummy,false);

#ifndef NDEBUG
  if( debugME() ) debug(dummy);
#endif

  createVertex(pme,ext);
}

void MEvv2vv::debug(double me2) const {
  if( !generator()->logfile().is_open() ) return;
  if( mePartonData()[0]->id() != 21 || mePartonData()[1]->id() != 21 ||
      mePartonData()[2]->id() != 5100021 || 
      mePartonData()[3]->id() != 5100021 ) return;
  tcSMPtr sm = generator()->standardModel();
  double gs4 = sqr( 4.*Constants::pi*sm->alphaS(scale()) );
  Energy2 s(sHat());
  Energy2 mf2 = meMomenta()[2].m2();
  Energy2 t3(tHat() - mf2), u4(uHat() - mf2);
  Energy4 s2(sqr(s)), t3s(sqr(t3)), u4s(sqr(u4)); 

  Energy4 num = s2 + t3s + u4s;  
  double analytic = 3.*mf2*( mf2*num/t3s/u4s - num/s/t3/u4 ) + 1.
    + sqr(num)*num/4./s2/t3s/u4s - t3*u4/s2;
  analytic *= 9.*gs4/8.;
  
  double diff = abs( analytic - me2 );
  if( diff  > 1e-4 ) {
    generator()->log() 
      << mePartonData()[0]->PDGName() << ","
      << mePartonData()[1]->PDGName() << "->"
      << mePartonData()[2]->PDGName() << ","
      << mePartonData()[3]->PDGName() << "   difference: " 
      << setprecision(10) << diff  << "  ratio: " << analytic/me2  << '\n';
  }
}
