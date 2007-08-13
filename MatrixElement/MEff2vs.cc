// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEff2vs class.
//

#include "MEff2vs.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "HardVertex.h"
#include <numeric>

using namespace Herwig;
using ThePEG::Helicity::incoming;
using ThePEG::Helicity::outgoing;

void MEff2vs::doinit() throw(InitException) {
  GeneralHardME::doinit();
  HPCount ndiags(numberOfDiags());
  theSca.resize(ndiags);
  theVec.resize(ndiags);
  theFerm.resize(ndiags);
  for(HPCount i = 0; i < ndiags; ++i) {
    HPDiagram current = getProcessInfo()[i];
    if( current.channelType == HPDiagram::sChannel ) {
      if( current.intermediate->iSpin() == PDT::Spin0 )
	theSca[i] = 
	  make_pair(dynamic_ptr_cast<FFSVertexPtr>(current.vertices.first), 
		    dynamic_ptr_cast<VSSVertexPtr>(current.vertices.second));
      else if( current.intermediate->iSpin() == PDT::Spin1 )
	theVec[i] = 
	  make_pair(dynamic_ptr_cast<FFVVertexPtr>(current.vertices.first), 
		    dynamic_ptr_cast<VVSVertexPtr>(current.vertices.second));
    }
    else if( current.channelType == HPDiagram::tChannel ) {
      if(current.intermediate->iSpin() == PDT::Spin1Half) {
	if( current.ordered.second ) 
	  theFerm[i] = 
	    make_pair(dynamic_ptr_cast<FFVVertexPtr>(current.vertices.first), 
		      dynamic_ptr_cast<FFSVertexPtr>(current.vertices.second));
	else
	  theFerm[i] = 
	    make_pair(dynamic_ptr_cast<FFVVertexPtr>(current.vertices.second), 
		      dynamic_ptr_cast<FFSVertexPtr>(current.vertices.first));
      }
    }
  }
}


void MEff2vs::persistentOutput(PersistentOStream & os) const {
  os << theSca << theVec << theFerm;
}

void MEff2vs::persistentInput(PersistentIStream & is, int) {
  is >> theSca >> theVec >> theFerm;
}

ClassDescription<MEff2vs> MEff2vs::initMEff2vs;
// Definition of the static class description member.

void MEff2vs::Init() {

  static ClassDocumentation<MEff2vs> documentation
    ("MEff2vs implements the ME calculation of the fermion-antifermion "
     "to vector-scalar hard process.");

}

double MEff2vs::me2() const {
  //set up wavefunctions
  SpinorVector ina(2);
  SpinorBarVector inb(2);
  VBVector outa(3);
  ScalarWaveFunction sca(meMomenta()[3], mePartonData()[3], Complex(1.),
			 outgoing);
  for(unsigned int ih = 0; ih < 2; ++ih) {
    ina[ih] = SpinorWaveFunction(meMomenta()[0], mePartonData()[0], 
				 ih, incoming);
    inb[ih] = SpinorBarWaveFunction(meMomenta()[1], mePartonData()[1], 
				    ih, incoming);
    outa[2*ih] = VectorWaveFunction(meMomenta()[2], mePartonData()[2], 
				    2*ih, outgoing);
  }
  if( mePartonData()[2]->mass() > 0.0*MeV ) {
    outa[1] = VectorWaveFunction(meMomenta()[2], mePartonData()[2], 
				 1, outgoing);
  }
  double full_me(0.);
  ffb2vsHeME(ina, inb, outa, sca, full_me);
  return full_me;
}

ProductionMatrixElement
MEff2vs::ffb2vsHeME(SpinorVector & sp, SpinorBarVector & spbar,
		    VBVector & vec, ScalarWaveFunction & sca, 
		    double & me2) const {
  Energy2 m2(scale());
  const vector<DVector> cfactors = getColourFactors();
  const HPCount ndiags = numberOfDiags();
  const size_t ncf = numberOfFlows();
  vector<double> me(ndiags, 0.);
  vector<Complex> diag(ndiags, Complex(0.)), flows(ncf, Complex(0.));
  me2 = 0.;
  ScalarWaveFunction interS; VectorWaveFunction interV; 
  SpinorBarWaveFunction interFB;
  bool mv(false);
  if( mePartonData()[2]->mass() == 0.0*MeV ) mv = true;
  ProductionMatrixElement prodme(PDT::Spin1Half, PDT::Spin1Half,
				 PDT::Spin1,PDT::Spin0);
  for(unsigned int ihel1 = 0; ihel1 < 2; ++ihel1) {
    for(unsigned int ihel2 = 0; ihel2 < 2; ++ihel2) {
      for(unsigned int ovhel = 0; ovhel < 3; ++ovhel) {
	if( mv ) continue;
	flows = vector<Complex>(ncf, Complex(0.));
	for(HPCount ix = 0; ix < ndiags; ++ix) {
	  HPDiagram current = getProcessInfo()[ix];
	  tcPDPtr offshell(current.intermediate);
	  if( current.channelType == HPDiagram::sChannel ) {
	    if( offshell->iSpin() == PDT::Spin0 ) {
	      interS = theSca[ix].first->evaluate(m2, 1, offshell, sp[ihel1], 
						  spbar[ihel2]);
	      diag[ix] = theSca[ix].second->evaluate(m2, vec[ovhel], sca, interS); 
	    }
	    else if( offshell->iSpin() == PDT::Spin1 ) {
	      interV = theVec[ix].first->evaluate(m2, 3, offshell, sp[ihel1], 
						  spbar[ihel2]);
	      diag[ix] = theVec[ix].second->evaluate(m2, vec[ovhel], interV, sca);
	    }
	    else diag[ix] = 0.0;
	  }
	  else if( current.channelType == HPDiagram::tChannel ) {
	    if( offshell->iSpin() == PDT::Spin1Half ) {
	      if( current.ordered.second ) {
		interFB = theFerm[ix].second->evaluate(m2, 3, offshell, 
						       spbar[ihel2], sca);
		diag[ix] = theFerm[ix].first->evaluate(m2, sp[ihel1], 
						       interFB, vec[ovhel]);
	      }
	      else {
		interFB = theFerm[ix].first->evaluate(m2, 3, offshell, 
						      spbar[ihel2], vec[ovhel]);
		diag[ix] = theFerm[ix].second->evaluate(m2, sp[ihel1], 
							interFB, sca);
	      }
	    }
	  }
	  else diag[ix] = 0.0;
	  
	  me[ix] += norm(diag[ix]);
	  //add to appropriate colourflow(s)
	  for(unsigned int iy = 0; iy < current.colourFlow.size(); ++iy)
	    flows[current.colourFlow[iy].first - 1] += 
	      current.colourFlow[iy].second * diag[ix];
	  
	}//end of diagram loop
	for(unsigned int ii = 0; ii < ncf; ++ii) 
	  for(unsigned int ij = 0; ij < ncf; ++ij)
	    me2 += cfactors[ii][ij]*(flows[ii]*conj(flows[ij])).real();

	prodme(ihel1, ihel2, ovhel, 0) = 
	  std::accumulate(flows.begin(), flows.end(), Complex(0.));
      }
    }
  }
  const double colourAvg = mePartonData()[0]->iColour() == PDT::Colour3 
    ? 1./9. : 1.;
  DVector save(ndiags);
  for(DVector::size_type ix = 0; ix < ndiags; ++ix)
    save[ix] = 0.25*colourAvg*me[ix];
  meInfo(save);
  me2 *= 0.25*colourAvg;
  return prodme;
}

Selector<const ColourLines *>
MEff2vs::colourGeometries(tcDiagPtr diag) const {
  static vector<ColourLines> cf(5);
  //11->11
  cf[0] = ColourLines("");
  //33b->11
  cf[1] = ColourLines("1 -2");
  cf[2] = ColourLines("1 2 -3");
  //33b->81
  cf[3] = ColourLines("1 4, -4 2 -3");
  cf[4] = ColourLines("1 2 4, -4 -3");
  
  HPDiagram current = getProcessInfo()[abs(diag->id()) - 1];
  PDT::Colour ca = getParticleData(current.incoming.first)->iColour();
  PDT::Colour cb = getParticleData(current.incoming.second)->iColour();
  PDT::Colour cc = getParticleData(current.outgoing.first)->iColour();
  PDT::Colour cd = getParticleData(current.outgoing.second)->iColour();
  Selector<const ColourLines *> sel;  
  if( ca == PDT::Colour0 && cb == PDT::Colour0 &&
      cc == PDT::Colour0 && cd == PDT::Colour0 ) 
    sel.insert(1., &cf[0]);
  else if( ca == PDT::Colour3 && cb == PDT::Colour3bar ) {
    if( cc == PDT::Colour0 && cd == PDT::Colour0 ) {
      if( current.channelType == HPDiagram::sChannel )
	sel.insert(1., &cf[1]);
      else
	sel.insert(1., &cf[2]);
    }
    else if( cc == PDT::Colour8 && cd == PDT::Colour0 ) {
      if( current.ordered.second )
	sel.insert(1., &cf[3]);
      else
	sel.insert(1., &cf[4]);
    }
  }
  else
    throw MEException()
      << "MEff2vs::colourGeometries() - Unknown incoming colour configuration "
      << ca << "  " << cb << Exception::runerror;  
  return sel;
}

void MEff2vs::constructVertex(tSubProPtr sub) {
  // Hard proces external particles
  ParticleVector hdp(4);
  hdp[0] = sub->incoming().first; 
  hdp[1] = sub->incoming().second;
  hdp[2] = sub->outgoing()[0]; 
  hdp[3] = sub->outgoing()[1];

  //check ordering
  if( hdp[0]->id() < hdp[1]->id() ) swap(hdp[0], hdp[1]);
  if( hdp[2]->dataPtr()->iSpin() == PDT::Spin0 ) swap(hdp[2], hdp[3]);
  
  SpinorVector sp;
  SpinorWaveFunction(sp, hdp[0], incoming, false, true);
  SpinorBarVector spbar;
  SpinorBarWaveFunction(spbar, hdp[1], incoming, false, true);
  VBVector vec;
  bool mv(hdp[2]->dataPtr()->mass() == 0.0*MeV);
  VectorWaveFunction(vec, hdp[2], outgoing, true, mv,true);
  ScalarWaveFunction sca(hdp[3], outgoing, true, true);
  
  double dummy(0.);
  ProductionMatrixElement prodme = ffb2vsHeME(sp, spbar, vec, sca, dummy);
  HardVertexPtr hardvertex = new_ptr(HardVertex());
  hardvertex->ME(prodme);
  for(ParticleVector::size_type i = 0; i < 4; ++i)
    dynamic_ptr_cast<SpinfoPtr>(hdp[i]->spinInfo())->
      setProductionVertex(hardvertex);
}


