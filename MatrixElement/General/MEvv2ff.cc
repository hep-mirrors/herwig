// -*- C++ -*-
//
// MEvv2ff.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEvv2ff class.
//

#include "MEvv2ff.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using ThePEG::Helicity::TensorWaveFunction;
using ThePEG::Helicity::incoming;
using ThePEG::Helicity::outgoing;

void MEvv2ff::doinit() {
  GeneralHardME::doinit();
  scalar_ .resize(numberOfDiags());
  fermion_.resize(numberOfDiags());
  vector_ .resize(numberOfDiags());
  tensor_ .resize(numberOfDiags());
  initializeMatrixElements(PDT::Spin1    , PDT::Spin1, 
			   PDT::Spin1Half, PDT::Spin1Half);

  for( size_t i = 0; i < numberOfDiags(); ++i ) {
    HPDiagram dg = getProcessInfo()[i];
    if( dg.channelType == HPDiagram::tChannel ) {
      AbstractFFVVertexPtr ffv1 = 
	dynamic_ptr_cast<AbstractFFVVertexPtr>(dg.vertices.first);
      AbstractFFVVertexPtr ffv2 = 
	dynamic_ptr_cast<AbstractFFVVertexPtr>(dg.vertices.second);
      fermion_[i] = make_pair(ffv1, ffv2);
    }
    else if( dg.channelType == HPDiagram::sChannel ) {
      if( dg.intermediate->iSpin() == PDT::Spin0 ) {
	AbstractVVSVertexPtr vvs = 
	  dynamic_ptr_cast<AbstractVVSVertexPtr>(dg.vertices.first );
	AbstractFFSVertexPtr ffs = 
	  dynamic_ptr_cast<AbstractFFSVertexPtr>(dg.vertices.second);
	scalar_[i] = make_pair(vvs,ffs);
      }
      else if( dg.intermediate->iSpin() == PDT::Spin1) {
	AbstractVVVVertexPtr vvv = 
	  dynamic_ptr_cast<AbstractVVVVertexPtr>(dg.vertices.first);
	AbstractFFVVertexPtr ffv = 
	  dynamic_ptr_cast<AbstractFFVVertexPtr>(dg.vertices.second);
	vector_[i] = make_pair(vvv,ffv);
      }
      else if(dg.intermediate->iSpin() == PDT::Spin2) {
	AbstractVVTVertexPtr vvt = 
	  dynamic_ptr_cast<AbstractVVTVertexPtr>(dg.vertices.first);
	AbstractFFTVertexPtr fft = 
	  dynamic_ptr_cast<AbstractFFTVertexPtr>(dg.vertices.second);
	tensor_[i] = make_pair(vvt,fft);
      }
    }
  }
}

double MEvv2ff::me2() const {
  // Set up wavefuctions
  VBVector v1(2), v2(2);
  SpinorVector sp(2); SpinorBarVector sbar(2);
  for( size_t i = 0; i < 2; ++i ) {
    v1[i] = VectorWaveFunction(rescaledMomenta()[0],mePartonData()[0], 2*i,
			       incoming);
    v2[i] = VectorWaveFunction(rescaledMomenta()[1],mePartonData()[1], 2*i,
			       incoming);
    sbar[i] = SpinorBarWaveFunction(rescaledMomenta()[2], mePartonData()[2], i,
				    outgoing);
    sp[i] = SpinorWaveFunction(rescaledMomenta()[3], mePartonData()[3], i,
			       outgoing);
  }
  double full_me(0.);
  vv2ffME(v1, v2, sbar, sp, full_me,true);

#ifndef NDEBUG 
  if( debugME() ) debug(full_me);
#endif

  return full_me;
}

ProductionMatrixElement 
MEvv2ff::vv2ffME(const VBVector & v1, const VBVector & v2,
		 const SpinorBarVector & sbar,const SpinorVector & sp, 
		 double & me2, bool first) const {
  const Energy mass = sp[0].mass();
  const Energy2 q2 = scale();
  // weights for the selection of the diagram
  vector<double> me(numberOfDiags(), 0.);
  // weights for the selection of the colour flow
  vector<double> flow(numberOfFlows(),0.);
  //sum over vector helicities
  for(unsigned int iv1 = 0; iv1 < 2; ++iv1) {
    for(unsigned int iv2 = 0; iv2 < 2; ++iv2) {
      //sum over fermion helicities
      for(unsigned int of1 = 0; of1 < 2; ++of1) {
	for(unsigned int of2 = 0; of2 < 2; ++of2) {
	  vector<Complex> flows(numberOfFlows(),0.);
	  for(HPCount ix = 0; ix < numberOfDiags(); ++ix) {
	    Complex diag(0.);
	    const HPDiagram & current = getProcessInfo()[ix];
	    PDPtr offshell = current.intermediate;
	    if(current.channelType == HPDiagram::tChannel && 
	       offshell->iSpin() == PDT::Spin1Half) {
	      if(current.ordered.second) {
                SpinorBarWaveFunction interF = fermion_[ix].first->
                  evaluate(q2, 3, offshell, sbar[of1], v1[iv1], mass);
                diag = fermion_[ix].second->
                  evaluate(q2, sp[of2], interF, v2[iv2]);
	      }
	      else {
		SpinorWaveFunction interF = fermion_[ix].second->
		  evaluate(q2, 3, offshell, sp[of2], v1[iv1], mass);
		diag = fermion_[ix].first->
		  evaluate(q2, interF, sbar[of1], v2[iv2]);
	      }
	    }
	    else if(current.channelType == HPDiagram::sChannel) {
	      if(offshell->iSpin() == PDT::Spin0) {
		ScalarWaveFunction interS = scalar_[ix].first->
		  evaluate(q2, 1, offshell, v1[iv1], v2[iv2]);
		diag = scalar_[ix].second->
		  evaluate(q2, sp[of2], sbar[of1], interS);
	      }
	      else if(offshell->iSpin() == PDT::Spin1) {
		VectorWaveFunction interV = vector_[ix].first->
		  evaluate(q2, 1, offshell, v1[iv1], v2[iv2]);
		diag = vector_[ix].second->
		  evaluate(q2, sp[of2], sbar[of1], interV);
	      }
	      else if(offshell->iSpin() == PDT::Spin2) {
		TensorWaveFunction interT = tensor_[ix].first->
		  evaluate(q2, 1, offshell, v1[iv1], v2[iv2]);
		diag = tensor_[ix].second->
		  evaluate(q2, sp[of2], sbar[of1], interT);
	      }
	    }
	    else diag = 0.;
	    me[ix] += norm(diag);
	    diagramME()[ix](2*iv1, 2*iv2, of1, of2) = diag;
	    //Compute flows
	    for(size_t iy = 0; iy < current.colourFlow.size(); ++iy) {
	      assert(current.colourFlow[iy].first<flows.size());
	      flows[current.colourFlow[iy].first] += 
		current.colourFlow[iy].second * diag;
	    }
	  }
	  // MEs for the different colour flows
	  for(unsigned int iy = 0; iy < numberOfFlows(); ++iy) 
	    flowME()[iy](2*iv1, 2*iv2, of1, of2) = flows[iy];
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

void MEvv2ff::persistentOutput(PersistentOStream & os) const {
  os << scalar_ << fermion_ << vector_ << tensor_;
}

void MEvv2ff::persistentInput(PersistentIStream & is, int) {
  is >> scalar_ >> fermion_ >> vector_ >> tensor_;
  initializeMatrixElements(PDT::Spin1    , PDT::Spin1, 
			   PDT::Spin1Half, PDT::Spin1Half);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEvv2ff,GeneralHardME>
describeHerwigMEvv2ff("Herwig::MEvv2ff", "Herwig.so");

void MEvv2ff::Init() {

  static ClassDocumentation<MEvv2ff> documentation
    ("The MEvv2ff class handles the ME calculation for the general "
     "spin configuration vector-vector to fermion-antifermion\n.");

}

void MEvv2ff::constructVertex(tSubProPtr sub) {
  ParticleVector ext = hardParticles(sub);
  // wavefunction with real momenta
  VBVector v1, v2;
  VectorWaveFunction(v1, ext[0], incoming, false, true);
  VectorWaveFunction(v2, ext[1], incoming, false, true);
  SpinorBarVector sbar;
  SpinorBarWaveFunction(sbar, ext[2], outgoing, true);
  SpinorVector sp;
  SpinorWaveFunction(sp, ext[3], outgoing, true);
  // rescale momenta
  setRescaledMomenta(ext);
  // wavefuncions with rescaled momenta
  VectorWaveFunction    v1r(rescaledMomenta()[0],
			    ext[0]->dataPtr(), incoming);
  VectorWaveFunction    v2r(rescaledMomenta()[1],
			    ext[1]->dataPtr(), incoming);
  SpinorBarWaveFunction sbr(rescaledMomenta()[2],
			    ext[2]->dataPtr(), outgoing);
  SpinorWaveFunction    spr(rescaledMomenta()[3],
			    ext[3]->dataPtr(), outgoing);
  for( unsigned int ihel = 0; ihel < 2; ++ihel ) {  
    v1r.reset(2*ihel);
    v1[ihel] = v1r;
    v2r.reset(2*ihel);
    v2[ihel] = v2r;
    sbr.reset(ihel);
    sbar[ihel] = sbr;
    spr.reset(ihel);
    sp[ihel] = spr;
  }
  double dummy(0.);
  ProductionMatrixElement pme = vv2ffME(v1, v2, sbar, sp, dummy,false);

#ifndef NDEBUG
  if( debugME() ) debug(dummy);
#endif

  createVertex(pme,ext);
}

void MEvv2ff::debug(double me2) const {
  if( !generator()->logfile().is_open() ) return;
  long id3(abs(mePartonData()[2]->id())), id4(abs(mePartonData()[3]->id()));
  if( mePartonData()[0]->id() != 21 || mePartonData()[1]->id() != 21 ||
      id3 != id4 || (id3 != 1000021 && id3 != 5100002 && id3 != 5100001 &&
		     id3 != 6100002 && id3 != 6100001) )
    return;
  tcSMPtr sm = generator()->standardModel();
  double gs4 = sqr( 4.*Constants::pi*sm->alphaS(scale()) );
  int Nc = sm->Nc();
  Energy2 s(sHat());
  Energy2 mf2 = meMomenta()[2].m2();
  Energy4 spt2 = uHat()*tHat() - sqr(mf2);
  Energy2 t3(tHat() - mf2), u4(uHat() - mf2);
    
  double analytic(0.);
  if( id3 == 1000021 ) {
   analytic = gs4*sqr(Nc)*u4*t3*
     ( sqr(u4) + sqr(t3) + 4.*mf2*s*spt2/u4/t3 ) * 
     ( 1./sqr(s*t3) + 1./sqr(s*u4) + 1./sqr(u4*t3) )/2./(Nc*Nc - 1.);
  }
  else {
    double brac = sqr(s)/6./t3/u4 - 3./8.;
    analytic = gs4*( -4.*sqr(mf2)*brac/t3/u4 + 4.*mf2*brac/s + brac 
		     - 1./3. + 3.*t3*u4/4/s/s);
  }
  double diff = abs(analytic - me2);
  if( diff > 1e-4 ) {
    generator()->log() 
      << mePartonData()[0]->PDGName() << ","
      << mePartonData()[1]->PDGName() << "->"
      << mePartonData()[2]->PDGName() << ","
      << mePartonData()[3]->PDGName() << "   difference: " 
      << setprecision(10) << diff << "  ratio: " << analytic/me2  << '\n';
  }

}
