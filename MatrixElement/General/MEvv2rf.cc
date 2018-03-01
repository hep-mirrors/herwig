// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEvv2rf class.
//

#include "MEvv2rf.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using ThePEG::Helicity::incoming;
using ThePEG::Helicity::outgoing;

IBPtr MEvv2rf::clone() const {
  return new_ptr(*this);
}

IBPtr MEvv2rf::fullclone() const {
  return new_ptr(*this);
}

void MEvv2rf::doinit() {
  GeneralHardME::doinit();
  scalar_ .resize(numberOfDiags());
  fermion_.resize(numberOfDiags());
  vector_ .resize(numberOfDiags());
  four_   .resize(numberOfDiags());
  initializeMatrixElements(PDT::Spin1    , PDT::Spin1, 
			   PDT::Spin3Half, PDT::Spin1Half);

  for( size_t i = 0; i < numberOfDiags(); ++i ) {
    HPDiagram dg = getProcessInfo()[i];
    if( dg.channelType == HPDiagram::tChannel ) {
      AbstractRFVVertexPtr rfv;
      AbstractFFVVertexPtr ffv;
      if(dg.ordered.second) {
	rfv = dynamic_ptr_cast<AbstractRFVVertexPtr>(dg.vertices.first);
	ffv = dynamic_ptr_cast<AbstractFFVVertexPtr>(dg.vertices.second);
      }
      else {
	rfv = dynamic_ptr_cast<AbstractRFVVertexPtr>(dg.vertices.second);
	ffv = dynamic_ptr_cast<AbstractFFVVertexPtr>(dg.vertices.first );
      }
      fermion_[i] = make_pair(rfv, ffv);
    }
    else if( dg.channelType == HPDiagram::sChannel ) {
      if( dg.intermediate->iSpin() == PDT::Spin0 ) {
	AbstractVVSVertexPtr vvs = 
	  dynamic_ptr_cast<AbstractVVSVertexPtr>(dg.vertices.first );
	AbstractRFSVertexPtr rfs = 
	  dynamic_ptr_cast<AbstractRFSVertexPtr>(dg.vertices.second);
	scalar_[i] = make_pair(vvs,rfs);
      }
      else if( dg.intermediate->iSpin() == PDT::Spin1) {
	AbstractVVVVertexPtr vvv = 
	  dynamic_ptr_cast<AbstractVVVVertexPtr>(dg.vertices.first);
	AbstractRFVVertexPtr rfv = 
	  dynamic_ptr_cast<AbstractRFVVertexPtr>(dg.vertices.second);
	vector_[i] = make_pair(vvv,rfv);
      }
    }
    else if ( dg.channelType == HPDiagram::fourPoint) {
      four_[i] = dynamic_ptr_cast<AbstractRFVVVertexPtr>(dg.vertices.first);
    }
  }
}

void MEvv2rf::persistentOutput(PersistentOStream & os) const {
  os << scalar_ << fermion_ << vector_ << four_;
}

void MEvv2rf::persistentInput(PersistentIStream & is, int) {
  is >> scalar_ >> fermion_ >> vector_ >> four_;  
  initializeMatrixElements(PDT::Spin1    , PDT::Spin1, 
			   PDT::Spin3Half, PDT::Spin1Half);
}

//The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEvv2rf,GeneralHardME>
describeHerwigMEvv2rf("Herwig::MEvv2rf", "Herwig.so");

void MEvv2rf::Init() {

  static ClassDocumentation<MEvv2rf> documentation
    ("The MEvv2rf class handes the ME calculation for vv -> rf");

}

double MEvv2rf::me2() const {
  // set up the vector wavefunctions
  VBVector v1(2), v2(2);
  for( size_t i = 0; i < 2; ++i ) {
    v1[i] = VectorWaveFunction(rescaledMomenta()[0],mePartonData()[0], 2*i,
     			       incoming);
    v2[i] = VectorWaveFunction(rescaledMomenta()[1],mePartonData()[1], 2*i,
     			       incoming);
  }
  // setup spinor wavefunctions and decide which case to use
  bool massless = mePartonData()[2]->mass()==ZERO;
  double full_me(0.);
  if(mePartonData()[2]->id()<0) {
    RSSpinorVector sp(4); SpinorBarVector sbar(2);
    for( size_t i = 0; i < 4; ++i ) {
      if(massless && (i==2||i==3)) continue; 
      sp[i] = RSSpinorWaveFunction(rescaledMomenta()[2], mePartonData()[2], i,
				   outgoing);
    }
    for( size_t i = 0; i < 2; ++i ) {
      sbar[i] = SpinorBarWaveFunction(rescaledMomenta()[3], mePartonData()[3], i,
				      outgoing);
    }
    vv2frME(v1, v2, sbar, sp, full_me,true);
  }
  else {
    SpinorVector sp(2); RSSpinorBarVector sbar(4);
    for( size_t i = 0; i < 4; ++i ) {
      if(massless && (i==2||i==3)) continue;
      sbar[i] = RSSpinorBarWaveFunction(rescaledMomenta()[2], mePartonData()[2], i,
					outgoing);
    }
    for( size_t i = 0; i < 2; ++i ) {
      sp[i] = SpinorWaveFunction(rescaledMomenta()[3], mePartonData()[3], i,
				 outgoing);
    }
    vv2rfME(v1, v2, sbar, sp, full_me,true);
  }
  return full_me;
}

ProductionMatrixElement 
MEvv2rf::vv2rfME(const VBVector & v1, const VBVector & v2,
		 const RSSpinorBarVector & sbar,const SpinorVector & sp, 
		 double & me2, bool first) const {
  // scale
  const Energy2 q2 = scale();
  // whether or not rs fermion is massless
  bool massless = mePartonData()[2]->mass()==ZERO;
  // weights for the selection of the diagram
  vector<double> me(numberOfDiags(), 0.);
  // weights for the selection of the colour flow
  vector<double> flow(numberOfFlows(),0.);
  //sum over vector helicities
  for(unsigned int iv1 = 0; iv1 < 2; ++iv1) {
    for(unsigned int iv2 = 0; iv2 < 2; ++iv2) {
      //sum over fermion helicities
      for(unsigned int of1 = 0; of1 < 4; ++of1) {
	if(massless && (of1==1 || of1==2) ) continue;
	for(unsigned int of2 = 0; of2 < 2; ++of2) {
	  vector<Complex> flows(numberOfFlows(),0.);
	  for(HPCount ix = 0; ix < numberOfDiags(); ++ix) {
   	    Complex diag(0.);
	    const HPDiagram & current = getProcessInfo()[ix];
	    tPDPtr offshell = current.intermediate;
	    if(current.channelType == HPDiagram::tChannel && 
	       offshell->iSpin() == PDT::Spin1Half) {
  	      if(current.ordered.second) {
                SpinorBarWaveFunction interF = fermion_[ix].first->
                  evaluate(q2, 3, offshell, sbar[of1], v1[iv1]);
                diag = fermion_[ix].second->
                  evaluate(q2, sp[of2], interF, v2[iv2]);
  	      }
   	      else {
  		SpinorWaveFunction interF = fermion_[ix].second->
  		  evaluate(q2, 3, offshell, sp[of2], v1[iv1]);
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
  	    }
  	    else if(current.channelType == HPDiagram::fourPoint) {
  	      diag = four_[ix]->evaluate(q2, sp[of2], sbar[of1], v1[iv1], v2[iv2]);
  	    }
  	    else
  	      assert(false);
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

ProductionMatrixElement 
MEvv2rf::vv2frME(const VBVector & v1, const VBVector & v2,
		 const SpinorBarVector & sbar,const RSSpinorVector & sp, 
		 double & me2, bool first) const {
  // scale
  const Energy2 q2 = scale();
  // whether or not rs fermion is massless
  bool massless = mePartonData()[2]->mass()==ZERO;
  // weights for the selection of the diagram
  vector<double> me(numberOfDiags(), 0.);
  // weights for the selection of the colour flow
  vector<double> flow(numberOfFlows(),0.);
  //sum over vector helicities
  for(unsigned int iv1 = 0; iv1 < 2; ++iv1) {
    for(unsigned int iv2 = 0; iv2 < 2; ++iv2) {
      //sum over fermion helicities
      for(unsigned int of1 = 0; of1 < 4; ++of1) {
  	for(unsigned int of2 = 0; of2 < 2; ++of2) {
	  if(massless && (of2==1 || of2==2) ) continue;
  	  vector<Complex> flows(numberOfFlows(),0.);
  	  for(HPCount ix = 0; ix < numberOfDiags(); ++ix) {
   	    Complex diag(0.);
  	    const HPDiagram & current = getProcessInfo()[ix];
  	    tPDPtr offshell = current.intermediate;
  	    if(current.channelType == HPDiagram::tChannel && 
  	       offshell->iSpin() == PDT::Spin1Half) {
  	      if(current.ordered.second) {
                SpinorWaveFunction interF = fermion_[ix].first->
                  evaluate(q2, 3, offshell, sp[of2], v1[iv1]);
                diag = fermion_[ix].second->
                  evaluate(q2,  interF, sbar[of1],v2[iv2]);
  	      }
   	      else {
   		SpinorBarWaveFunction interF = fermion_[ix].second->
   		  evaluate(q2, 3, offshell, sbar[of1], v1[iv1]);
  		diag = fermion_[ix].first->
  		  evaluate(q2, sp[of2], interF, v2[iv2]);
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
   	    }
  	    else if(current.channelType == HPDiagram::fourPoint) {
  	      diag = four_[ix]->evaluate(q2, sp[of2], sbar[of1], v1[iv1], v2[iv2]);
  	    }
  	    else
  	      assert(false);
  	    me[ix] += norm(diag);
  	    diagramME()[ix](2*iv1, 2*iv2, of2, of1) = diag;
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

void MEvv2rf::constructVertex(tSubProPtr sub) {
  ParticleVector ext = hardParticles(sub);
  // vector wavefunctions are common do them first
  // wavefunction with real momenta
  VBVector v1, v2;
  VectorWaveFunction(v1, ext[0], incoming, false, true);
  VectorWaveFunction(v2, ext[1], incoming, false, true);
  // rescale momenta
  setRescaledMomenta(ext);
  // wavefunctions with rescaled momenta
  VectorWaveFunction    v1r(rescaledMomenta()[0],
			    ext[0]->dataPtr(), incoming);
  VectorWaveFunction    v2r(rescaledMomenta()[1],
			    ext[1]->dataPtr(), incoming);
  ProductionMatrixElement pme;
  for( unsigned int ihel = 0; ihel < 2; ++ihel ) {  
    v1r.reset(2*ihel);
    v1[ihel] = v1r;
    v2r.reset(2*ihel);
    v2[ihel] = v2r;
  }
  // setup spinor wavefunctions and decide which case to use
  bool massless = mePartonData()[2]->mass()==ZERO;
  double dummy(0.);
  if(ext[2]->id()<0) {
    RSSpinorVector sp;
    RSSpinorWaveFunction(sp, ext[2], outgoing, true);
    SpinorBarVector sbar;
    SpinorBarWaveFunction(sbar, ext[3], outgoing, true);
    // wavefuncions with rescaled momenta
    SpinorBarWaveFunction sbr(rescaledMomenta()[3],
			      ext[3]->dataPtr(), outgoing);
    RSSpinorWaveFunction    spr(rescaledMomenta()[2],
				ext[2]->dataPtr(), outgoing);
    for( unsigned int ihel = 0; ihel < 2; ++ihel ) {  
      sbr.reset(ihel);
      sbar[ihel] = sbr;
    }
    for( unsigned int ihel = 0; ihel < 4; ++ihel ) {
      if(massless && (ihel==1 || ihel==2)) continue;
      spr.reset(ihel);
      sp[ihel] = spr;
    }
    pme = vv2frME(v1, v2, sbar, sp, dummy,false);
  }
  else {
    SpinorVector sp;
    SpinorWaveFunction(sp, ext[3], outgoing, true);
    RSSpinorBarVector sbar;
    RSSpinorBarWaveFunction(sbar, ext[2], outgoing, true);
    // wavefuncions with rescaled momenta
    RSSpinorBarWaveFunction sbr(rescaledMomenta()[2],
			      ext[2]->dataPtr(), outgoing);
    SpinorWaveFunction    spr(rescaledMomenta()[3],
			      ext[3]->dataPtr(), outgoing);
    for( unsigned int ihel = 0; ihel < 4; ++ihel ) {  
      if(massless && (ihel==1 || ihel==2)) continue;
      sbr.reset(ihel);
      sbar[ihel] = sbr;
    }
    for( unsigned int ihel = 0; ihel < 2; ++ihel ) {  
      spr.reset(ihel);
      sp[ihel] = spr;
    }
    pme = vv2rfME(v1, v2, sbar, sp, dummy,false);
  }
  createVertex(pme,ext);
}
