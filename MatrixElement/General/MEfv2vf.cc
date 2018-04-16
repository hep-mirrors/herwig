// -*- C++ -*-
//
// MEfv2vf.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEfv2vf class.
//

#include "MEfv2vf.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using ThePEG::Helicity::incoming;
using ThePEG::Helicity::outgoing;
using ThePEG::Helicity::SpinorWaveFunction;
using ThePEG::Helicity::SpinorBarWaveFunction;
using ThePEG::Helicity::VectorWaveFunction;

void MEfv2vf::doinit() {
  GeneralHardME::doinit();
  fermion_.resize(numberOfDiags());
  vector_.resize(numberOfDiags());
  initializeMatrixElements(PDT::Spin1Half, PDT::Spin1, 
			   PDT::Spin1, PDT::Spin1Half);
  for(HPCount ix = 0; ix < numberOfDiags(); ++ix) {
    HPDiagram diagram = getProcessInfo()[ix];
    PDT::Spin offspin = diagram.intermediate->iSpin();
    if(diagram.channelType == HPDiagram::sChannel ||
       ( diagram.channelType == HPDiagram::tChannel 
	 && offspin == PDT::Spin1Half)) {
      AbstractFFVVertexPtr vert1 = dynamic_ptr_cast<AbstractFFVVertexPtr>
		(diagram.vertices.first);
      AbstractFFVVertexPtr vert2 = dynamic_ptr_cast<AbstractFFVVertexPtr>
	(diagram.vertices.second);
      fermion_[ix] = make_pair(vert1, vert2);
    }
    else {
      if(offspin == PDT::Spin1) {
	AbstractFFVVertexPtr vert1 = dynamic_ptr_cast<AbstractFFVVertexPtr>
	  (diagram.vertices.first);
	AbstractVVVVertexPtr vert2 = dynamic_ptr_cast<AbstractVVVVertexPtr>
	  (diagram.vertices.second);
	vector_[ix] = make_pair(vert1, vert2);
      }
    }
  }
}

double MEfv2vf::me2() const {
  //wavefunctions
  SpinorVector sp(2);
  VBVector vecIn(2), vecOut(3);
  SpinorBarVector spb(2);
  double fullme(0.);
  bool mc = !(mePartonData()[2]->mass() > ZERO);
  if(mePartonData()[0]->id() > 0) {
    for(unsigned int i = 0; i < 2; ++i) {
      sp[i] = SpinorWaveFunction(rescaledMomenta()[0], mePartonData()[0], i, 
				 incoming);
      vecIn[i] = VectorWaveFunction(rescaledMomenta()[1], mePartonData()[1], 2*i, 
				    incoming);
      vecOut[2*i] = VectorWaveFunction(rescaledMomenta()[2], mePartonData()[2], 2*i, 
				     outgoing);
      spb[i] = SpinorBarWaveFunction(rescaledMomenta()[3], mePartonData()[3], i, 
				     outgoing);
    }
    if( !mc )
      vecOut[1] = VectorWaveFunction(rescaledMomenta()[2], mePartonData()[2], 1, 
				     outgoing);
    fv2vfHeME(sp, vecIn, vecOut, mc, spb, fullme,true);
  }
  else {
    for(unsigned int i = 0; i < 2; ++i) {
      spb[i] = SpinorBarWaveFunction(rescaledMomenta()[0], mePartonData()[0], i, 
				     incoming);
      vecIn[i] = VectorWaveFunction(rescaledMomenta()[1], mePartonData()[1], 2*i, 
				    incoming);
      vecOut[2*i] = VectorWaveFunction(rescaledMomenta()[2], mePartonData()[2], 2*i, 
				     outgoing);
      sp[i] = SpinorWaveFunction(rescaledMomenta()[3], mePartonData()[3], i, 
				 outgoing);
    }
    if( !mc )
      vecOut[1] = VectorWaveFunction(rescaledMomenta()[2], mePartonData()[2], 1, 
				     outgoing);
    fbv2vfbHeME(spb, vecIn, vecOut, mc, sp, fullme,true);
  }

#ifndef NDEBUG
  if( debugME() ) debug(fullme);
#endif

  return fullme;
}

ProductionMatrixElement
MEfv2vf::fv2vfHeME(const SpinorVector & spIn,  const VBVector & vecIn, 
		   const VBVector & vecOut, bool mc, 
		   const SpinorBarVector & spbOut, 
		   double & me2, bool first) const {
  const Energy2 q2(scale());
  // weights for the selection of the diagram
  vector<double> me(numberOfDiags(), 0.);
  // weights for the selection of the colour flow
  vector<double> flow(numberOfFlows(),0.);
  me2 = 0.;
  //loop over helicities
  for(unsigned int ifh = 0; ifh < 2; ++ifh) {
    for(unsigned int ivh = 0; ivh < 2; ++ivh) {
      for(unsigned int ovh = 0; ovh < 3; ++ovh) {
	if(mc && ovh == 1) ++ovh;
	for(unsigned int ofh = 0; ofh < 2; ++ofh) {
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
		diag = 0.0;
	    }
	    else if(current.channelType == HPDiagram::sChannel) {
	      if(offshell->CC()) offshell = offshell->CC();
	      unsigned int iopt = abs(offshell->id())==abs(spIn[ifh].particle()->id()) ? 5 : 1;
	      SpinorBarWaveFunction interFB = fermion_[ix].second->
		evaluate(q2, iopt, offshell, spbOut[ofh], vecOut[ovh]);
	      diag = fermion_[ix].first->
		evaluate(q2, spIn[ifh], interFB, vecIn[ivh]);
	    }
	    me[ix] += norm(diag);
	    diagramME()[ix](ifh, 2*ivh, ovh, ofh) = diag;
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
MEfv2vf::fbv2vfbHeME(const SpinorBarVector & spbIn, const VBVector & vecIn, 
		     const VBVector & vecOut, bool mc, 
		     const SpinorVector & spOut,
		     double & me2, bool first) const {
  const Energy2 q2(scale());
  // weights for the selection of the diagram
  vector<double> me(numberOfDiags(), 0.);
  // weights for the selection of the colour flow
  vector<double> flow(numberOfFlows(),0.);
  me2 = 0.;
  //loop over helicities
  for(unsigned int ifh = 0; ifh < 2; ++ifh) {
    for(unsigned int ivh = 0; ivh < 2; ++ivh) {
      for(unsigned int ovh = 0; ovh < 3; ++ovh) {
	if(mc && ovh == 1) ++ovh;
	for(unsigned int ofh = 0; ofh < 2; ++ofh) {
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
	      else diag = 0.0;
	    }
	    else if(current.channelType == HPDiagram::sChannel) {
	      if(offshell->iSpin() == PDT::Spin1Half) {
		SpinorBarWaveFunction interFB = fermion_[ix].first->
		  evaluate(q2, 1, offshell, spbIn[ifh], vecIn[ivh]);
		diag = fermion_[ix].second->
		  evaluate(q2, spOut[ofh], interFB, vecOut[ovh]);
	      }
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

void MEfv2vf::persistentOutput(PersistentOStream & os) const {
  os << fermion_ << vector_;
}

void MEfv2vf::persistentInput(PersistentIStream & is, int) {
  is >> fermion_ >> vector_;
  initializeMatrixElements(PDT::Spin1Half, PDT::Spin1, 
			   PDT::Spin1, PDT::Spin1Half);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEfv2vf,GeneralHardME>
describeHerwigMEfv2vf("Herwig::MEfv2vf", "Herwig.so");

void MEfv2vf::Init() {

  static ClassDocumentation<MEfv2vf> documentation
    ("This is the implementation of the matrix element for a fermion-vector boson"
     "to a vector-fermion.");

}

void MEfv2vf::constructVertex(tSubProPtr sub) {
  ParticleVector ext = hardParticles(sub);
  VBVector v1, v3;
  bool mc = !(ext[2]->data().mass() > ZERO);
  SpinorVector sp;  
  SpinorBarVector sbar;
  VectorWaveFunction(v1, ext[1], incoming, false, true);
  VectorWaveFunction(v3, ext[2], outgoing, true, mc);
  double dummy(0.);
  //Need to use rescale momenta to calculate matrix element
  setRescaledMomenta(ext);
  // wavefunctions with rescaled momenta 
  VectorWaveFunction vir(rescaledMomenta()[1],
			 ext[1]->dataPtr(), incoming);
  VectorWaveFunction vor(rescaledMomenta()[2],
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
      vir.reset(2*ihel);
      v1[ihel] = vir;
      vor.reset(2*ihel);
      v3[2*ihel] = vor;
      sbr.reset(ihel);
      sbar[ihel] = sbr;
    }
    if( !mc ) {
      vor.reset(1);
      v3[1] = vor;
    }
    ProductionMatrixElement pme = fv2vfHeME(sp, v1, v3, mc, sbar, dummy,false);
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
      vir.reset(2*ihel);
      v1[ihel] = vir;
      vor.reset(2*ihel);
      v3[2*ihel] = vor;
      spr.reset(ihel);
      sp[ihel] = spr;
    }
    if( !mc ) {
      vor.reset(1);
      v3[1] = vor;
    }
    ProductionMatrixElement pme = fbv2vfbHeME(sbar, v1, v3, mc, sp, dummy,false);
    createVertex(pme,ext);
  }
  
#ifndef NDEBUG
  if( debugME() ) debug(dummy);
#endif
}

void MEfv2vf::debug(double me2) const {
  if( !generator()->logfile().is_open() ) return;
  long id1 = abs(mePartonData()[0]->id());
  long id4 = abs(mePartonData()[3]->id());
  if( (id1 != 1 && id1 != 2) || mePartonData()[1]->id() != 21 ||
      mePartonData()[2]->id() != 5100021 || 
      (id4 != 5100001 && id4 != 5100002 &&
       id4 != 6100001 && id4 != 6100002) ) return;
  tcSMPtr sm = generator()->standardModel();
  double gs4 = sqr( 4.*Constants::pi*sm->alphaS(scale()) );
  Energy2 s(sHat());
  Energy2 mf2 = meMomenta()[2].m2();
  //  Energy4 spt2 = uHat()*tHat() - sqr(mf2);
  //swap t and u as formula defines process vf->vf
  Energy2 t3(uHat() - mf2), u4(tHat() - mf2);
  Energy4 s2(sqr(s)), t3s(sqr(t3)), u4s(sqr(u4));

  double analytic = -gs4*( 5.*s2/12./t3s + s2*s/t3s/u4 + 11.*s*u4/6./t3s
			   + 5.*u4s/12./t3s + u4s*u4/s/t3s)/3.;
  double diff = abs(analytic - me2);
  if( diff > 1e-4 ) {
    generator()->log() 
      << mePartonData()[0]->PDGName() << ","
      << mePartonData()[1]->PDGName() << "->"
      << mePartonData()[2]->PDGName() << ","
      << mePartonData()[3]->PDGName() << "   difference: " 
      << setprecision(10) << diff << "  ratio: " << analytic/me2 << '\n';
  }
}
