// -*- C++ -*-
//
// MEvv2ss.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEvv2ss class.
//

#include "MEvv2ss.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "ThePEG/Helicity/Vertex/Scalar/SSSVertex.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Helicity/SpinInfo.h"
#include "HardVertex.h"
#include <numeric>

using namespace Herwig;
using ThePEG::Helicity::TensorWaveFunction;
using ThePEG::Helicity::incoming;
using ThePEG::Helicity::outgoing;
using ThePEG::Helicity::SpinfoPtr;

void MEvv2ss::doinit() throw(InitException) {
  GeneralHardME::doinit();
  size_t ndiags = numberOfDiags();
  theSca.resize(ndiags); theVec.resize(ndiags); 
  theTen.resize(ndiags);
  
  for(size_t i = 0; i < ndiags; ++i ) {
    HPDiagram dg = getProcessInfo()[i];
    if( !dg.intermediate ) {
      theContact = dynamic_ptr_cast<AbstractVVSSVertexPtr>(dg.vertices.first);
    }
    else if(dg.channelType == HPDiagram::tChannel) {
      AbstractVSSVertexPtr vss1 = 
	dynamic_ptr_cast<AbstractVSSVertexPtr>(dg.vertices.first);
      AbstractVSSVertexPtr vss2 =
	dynamic_ptr_cast<AbstractVSSVertexPtr>(dg.vertices.second);
      theSca[i] = make_pair(vss1, vss2);
    }
    else {
      if( dg.intermediate->iSpin() == PDT::Spin1 ) {
	AbstractVVVVertexPtr vvv =
	  dynamic_ptr_cast<AbstractVVVVertexPtr>(dg.vertices.first);
	AbstractVSSVertexPtr vss =
	  dynamic_ptr_cast<AbstractVSSVertexPtr>(dg.vertices.second);
	theVec[i] = make_pair(vvv, vss);
      }
      else if( dg.intermediate->iSpin() == PDT::Spin2 ) {
	AbstractVVTVertexPtr vvt =
	  dynamic_ptr_cast<AbstractVVTVertexPtr>(dg.vertices.first);
	AbstractSSTVertexPtr sst = 
	  dynamic_ptr_cast<AbstractSSTVertexPtr>(dg.vertices.second);
	theTen[i] = make_pair(vvt, sst);
      }
    }
  }

}

double MEvv2ss::me2() const {
  VBVector v1(2), v2(2);
  for( size_t i = 0; i < 2; ++i ) {
    v1[i] = VectorWaveFunction(meMomenta()[0],mePartonData()[0], 2*i,
			       incoming);
    v2[i] = VectorWaveFunction(meMomenta()[1],mePartonData()[1], 2*i,
			       incoming);
  }
  ScalarWaveFunction sca1(meMomenta()[2],mePartonData()[2],
			  Complex(1.,0.),outgoing);
  ScalarWaveFunction sca2(meMomenta()[3],mePartonData()[3],
			  Complex(1.,0.),outgoing);
  double full_me(0.);
  vv2ssME(v1, v2, sca1, sca2, full_me);
  
#ifndef NDEBUG
  if( debugME() ) debug(full_me);
#endif

  return full_me;
}

ProductionMatrixElement 
MEvv2ss::vv2ssME(const VBVector & v1, const VBVector & v2,
		 const ScalarWaveFunction & sca1, 
		 const ScalarWaveFunction & sca2, double & me2) const {
  const HPCount ndiags(numberOfDiags());
  const size_t ncf(numberOfFlows());
  const Energy2 m2(scale());
  const vector<vector<double> > cfactors = getColourFactors();
  vector<double> me(ndiags,0.);
  vector<Complex> flows(ncf, Complex(0.)), diag(ndiags, Complex(0.));
  ScalarWaveFunction interS; VectorWaveFunction interV;
  TensorWaveFunction interT;
  //keep location of contact diagram to zero me calc for this diagram
  //so as not to add contact term to diagram selector
  HPCount contact = ndiags;
  //loop over vector helicities
  ProductionMatrixElement pme(PDT::Spin1, PDT::Spin1, PDT::Spin0, 
			      PDT::Spin0);
  for(unsigned int iv1 = 0; iv1 < 2; ++iv1) {
    for(unsigned int iv2 = 0; iv2 < 2; ++iv2) {
      //loop over diagrams
      flows = vector<Complex>(ncf, Complex(0.));
      for(HPCount ix = 0; ix < ndiags; ++ix){
	HPDiagram current = getProcessInfo()[ix];
	//do four-point diag first
	if(current.channelType == HPDiagram::fourPoint) {
	  diag[ix] = theContact->evaluate(m2, v1[iv1], v2[iv2], sca1, sca2);
	  contact = ix;
	}
	else {
	  tcPDPtr offshell = current.intermediate;
	  if(current.channelType == HPDiagram::tChannel) {
	    if(current.ordered.second) {
	      interS = theSca[ix].first->evaluate(m2, 3, offshell, v1[iv1], 
						  sca1);
	      diag[ix] = theSca[ix].second->evaluate(m2, v2[iv2], interS, sca2);
	    }
	    else {
	      interS = theSca[ix].first->evaluate(m2, 3, offshell, v1[iv1], 
						  sca2);
	      diag[ix] = theSca[ix].second->evaluate(m2, v2[iv2], interS, sca1);
	    }
	  }
	  else if(current.channelType == HPDiagram::sChannel) {
	    if(offshell->iSpin() == PDT::Spin1) {
	      interV = theVec[ix].first->evaluate(m2, 1, offshell, v1[iv1], 
						  v2[iv2]);
	      diag[ix] = theVec[ix].second->evaluate(m2, interV, sca1, sca2);
	    }
	    else if(offshell->iSpin() == PDT::Spin2) {
	      interT = theTen[ix].first->evaluate(m2, 1, offshell, v1[iv1], 
						  v2[iv2]);
	      diag[ix] = theTen[ix].second->evaluate(m2, sca1, sca2, interT);
	    }
	  }
	  else
	    diag[ix] = 0.;
	}
	me[ix] += norm(diag[ix]);
	//set appropriate element in ProductionMatrixElement
	pme(iv1, iv2, 0, 0) = 
	  std::accumulate(flows.begin(), flows.end(), Complex(0.0, 0.0));
	
 	//colourflows
	for(size_t iy = 0; iy < current.colourFlow.size(); ++iy)
	  flows[current.colourFlow[iy].first - 1] += 
	    current.colourFlow[iy].second * diag[ix];
	
      }//end of diagram loop
      for(size_t ii = 0; ii < ncf; ++ii)
	for(size_t ij = 0; ij < ncf; ++ij)	  
	  me2 += cfactors[ii][ij]*(flows[ii]*conj(flows[ij])).real();
       
    }//iv2 end
  }//iv1 end
  DVector save(ndiags);
  const double ifact = mePartonData()[2]->id() == mePartonData()[3]->id() ?
    0.5 : 1.;
  // contact must have been set during the loop run
  assert(contact != ndiags);
  me[contact] = 0.;
  for(DVector::size_type ix = 0; ix < ndiags; ++ix)
    save[ix] = ifact*me[ix]/256.;
  meInfo(save);
  me2 *= ifact/256.;
  return pme;
}

Selector<const ColourLines *>
MEvv2ss::colourGeometries(tcDiagPtr diag) const {
  //88->33bar
  static ColourLines cf1("1 4, -1 2 3, -3 -5");
  static ColourLines cf2("-1 -5, 3 4, 1 2 -3");
  static ColourLines cf3("1 3 4, -2 -3 -5, -1 2");
  //88->11
  static ColourLines cf4("1 -2, 2 -1");
  HPDiagram current = getProcessInfo()[abs(diag->id()) - 1];
  Selector<const ColourLines *> sel;
  if(current.channelType == HPDiagram::tChannel ) {
    if(current.ordered.second)
      sel.insert(1., &cf1);
    else
      sel.insert(1., &cf2);
  }
  else if(current.channelType == HPDiagram::sChannel) {
    if(mePartonData()[2]->iColour() == PDT::Colour0)
      sel.insert(1., &cf4);
    else
      sel.insert(1., &cf3);
  }
  else
    throw MEException() << "MEvv2ss::colourGeometries - "
			<< "Cannot find correct colour configuration! \n"
			<< Exception::warning;
  return sel;
}


void MEvv2ss::persistentOutput(PersistentOStream & os) const {
  os << theSca << theVec << theTen << theContact;
}

void MEvv2ss::persistentInput(PersistentIStream & is, int) {
  is >> theSca >> theVec >> theTen >> theContact;
}

ClassDescription<MEvv2ss> MEvv2ss::initMEvv2ss;
// Definition of the static class description member.

void MEvv2ss::Init() {

  static ClassDocumentation<MEvv2ss> documentation
    ("This class implements the ME for the vector-vector to scalar-scalar "
     "hard-process");

}

void MEvv2ss::constructVertex(tSubProPtr sub) {
  ParticleVector ext(4);
  ext[0] = sub->incoming().first;
  ext[1] = sub->incoming().second;
  ext[2] = sub->outgoing()[0];
  ext[3] = sub->outgoing()[1];
  
  VBVector v1, v2;
  VectorWaveFunction(v1, ext[0], incoming, false, true, true);
  VectorWaveFunction(v2, ext[1], incoming, false, true, true);
  //function to calculate me2 expects massless incoming vectors
  // and this constructor sets the '1' polarisation at element [2] 
  //in the vector
  v1[1] = v1[2]; 
  v2[1] = v2[2];
  ScalarWaveFunction sca1(ext[2], outgoing, true, true);
  ScalarWaveFunction sca2(ext[3], outgoing, true, true);

  double dummy(0.);
  ProductionMatrixElement pme = vv2ssME(v1, v2, sca1, sca2, dummy);

#ifndef NDEBUG
  if( debugME() ) debug(dummy);
#endif

  HardVertexPtr hv = new_ptr(HardVertex());
  hv->ME(pme);
  for(unsigned int i = 0; i < 4; ++i )
    dynamic_ptr_cast<SpinfoPtr>(ext[i]->spinInfo())->setProductionVertex(hv);
}

void MEvv2ss::debug(double me2) const {
  if( !generator()->logfile().is_open() ) return;
  //SUSY gg>~q~q
  long id3 = abs(mePartonData()[2]->id());
  long id4 = abs(mePartonData()[3]->id());
  if( mePartonData()[0]->id() != 21 || mePartonData()[1]->id() != 21 ||
      (id3 < 1000001 &&  id3 > 1000006 ) || (id3 < 2000001 && id3 > 2000006 ) ||
      (id4 < 1000001 &&  id4 > 1000006 ) || (id4 < 2000001 &&  id4 > 2000006 ) ) 
    return;
  tcSMPtr sm = generator()->standardModel();
  double gs4 = sqr( 4.*Constants::pi*sm->alphaS(scale()) );
  int Nc = sm->Nc();
  Energy4 s2 = sqr(sHat());
  Energy2 m3s = meMomenta()[2].m2();
  Energy2 m4s = meMomenta()[3].m2();
  Energy4 spt2 = uHat()*tHat() - m3s*m4s;
  Energy4 t3s = sqr(tHat() - m3s);
  Energy4 u4s = sqr(uHat() - m4s);

  double analytic = gs4*Nc*( sqr(spt2) + s2*m3s*m4s ) * 
    ( u4s + t3s - s2/sqr(Nc) )/2./(sqr(Nc) - 1.)/s2/t3s/u4s;
  double diff = abs(analytic - me2);
  if(  diff > 1e-4 ) {
    generator()->log() 
      << mePartonData()[0]->PDGName() << ","
      << mePartonData()[1]->PDGName() << "->"
      << mePartonData()[2]->PDGName() << ","
      << mePartonData()[3]->PDGName() << "   difference: "
      << setprecision(10) << diff << "  ratio: " << analytic/me2 << '\n';
  }
}
