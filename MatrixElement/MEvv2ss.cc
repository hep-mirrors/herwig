// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEvv2ss class.
//

#include "MEvv2ss.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/TensorWaveFunction.h"
#include "Herwig++/Helicity/Vertex/Scalar/VVSVertex.h"
#include "Herwig++/Helicity/Vertex/Scalar/VSSVertex.h"
#include "Herwig++/Helicity/Vertex/Vector/VVVVertex.h"
#include "Herwig++/Helicity/Vertex/Scalar/SSSVertex.h"
#include "Herwig++/Helicity/Vertex/Tensor/VVTVertex.h"
#include "Herwig++/Helicity/Vertex/Tensor/SSTVertex.h"
#include "Herwig++/Helicity/Vertex/Scalar/VVSSVertex.h"
#include <cassert>

using namespace Herwig;
using Herwig::Helicity::VectorWaveFunction;
using Herwig::Helicity::ScalarWaveFunction;
using Herwig::Helicity::TensorWaveFunction;
using Herwig::Helicity::VVSVertexPtr;
using Herwig::Helicity::VSSVertex;
using Herwig::Helicity::VVVVertexPtr;
using Herwig::Helicity::SSSVertex;
using Herwig::Helicity::VVTVertexPtr;
using Herwig::Helicity::SSTVertexPtr;
using Herwig::Helicity::VVSSVertex;
using Herwig::Helicity::incoming;
using Herwig::Helicity::outgoing;

typedef Ptr<VSSVertex>::pointer VSSVertexPtr;
typedef Ptr<SSSVertex>::pointer SSSVertexPtr;
typedef Ptr<VVSSVertex>::pointer VVSSVertexPtr;

double MEvv2ss::me2() const {
  //setup wavefunctions
  VectorWaveFunction vec1(meMomenta()[0],mePartonData()[0],incoming);
  VectorWaveFunction vec2(meMomenta()[1],mePartonData()[1],incoming);
  //2 is positve id, 3 is negative id
  ScalarWaveFunction sca1(meMomenta()[2],mePartonData()[2],
			  Complex(1.,0.),outgoing);
  ScalarWaveFunction sca2(meMomenta()[3],mePartonData()[3],
			  Complex(1.,0.),outgoing);
  const HPCount ndiags(numberOfDiags());
  const size_t ncf(numberOfFlows());
  const Energy2 m2(scale());
  const vector<vector<double> > cfactors = getColourFactors();
  vector<double> me(ndiags,0.);
  vector<Complex> flows(ncf, Complex(0.)), diag(ndiags, Complex(0.));
  double full_me(0.);
  ScalarWaveFunction interS; VectorWaveFunction interV;
  TensorWaveFunction interT;
  //keep location of contact diagram to zero me calc for this diagram
  //so as not to add contact term to diagram selector
  HPCount contact = ndiags;
  //loop over vector helicities
  for(unsigned int vhel1 = 0; vhel1 < 3; vhel1 += 2) {
    vec1.reset(vhel1);
    for(unsigned int vhel2 = 0; vhel2 < 3;vhel2 += 2) {
      vec2.reset(vhel2);
      //loop over diagrams
      flows = vector<Complex>(ncf, Complex(0.));
      for(HPCount ix = 0; ix < ndiags; ++ix){
	HPDiagram current = getProcessInfo()[ix];
	//do four-point diag first
	if(current.channelType == HPDiagram::fourPoint) {
	  VVSSVertexPtr vert = dynamic_ptr_cast<VVSSVertexPtr>
	    (current.vertices.first);
	  diag[ix] = vert->evaluate(m2, vec1, vec2, sca1, sca2);
	  contact = ix;
	}
	else {
	  tcPDPtr offshell = current.intermediate;
	  if(current.channelType == HPDiagram::tChannel) {
	    VSSVertexPtr vert = dynamic_ptr_cast<VSSVertexPtr>
	      (current.vertices.first);
	    if(current.ordered.second) {
	      interS = vert->evaluate(m2, 3, offshell, vec1, sca1);
	      diag[ix] = vert->evaluate(m2, vec2, interS, sca2);	      
	    }
	    else {
	      interS = vert->evaluate(m2, 3, offshell, vec1, sca2);
	      diag[ix] = vert->evaluate(m2, vec2, interS, sca1);
	    }
	  }
	  else if(current.channelType == HPDiagram::sChannel) {
	    if(offshell->iSpin() == PDT::Spin1) {
	      VVVVertexPtr vert1 = dynamic_ptr_cast<VVVVertexPtr>
		(current.vertices.first);
	      VSSVertexPtr vert2 = dynamic_ptr_cast<VSSVertexPtr>
		(current.vertices.second);
	      interV = vert1->evaluate(m2, 1, offshell, vec1, vec2);
	      diag[ix] = vert2->evaluate(m2, interV, sca1, sca2);
	    }
	    else if(offshell->iSpin() == PDT::Spin2) {
	      VVTVertexPtr vert1 = dynamic_ptr_cast<VVTVertexPtr>
		(current.vertices.first);
	      SSTVertexPtr vert2 = dynamic_ptr_cast<SSTVertexPtr>
		(current.vertices.second);
	      interT = vert1->evaluate(m2, 1, offshell, vec1, vec2);
	      diag[ix] = vert2->evaluate(m2, sca1, sca2, interT);
	    }
	  }
	  else
	    diag[ix] = 0.;
	}
	me[ix] += norm(diag[ix]);
 	//colourflows
	for(size_t iy = 0; iy < current.colourFlow.size(); ++iy)
	  flows[current.colourFlow[iy].first - 1] += 
	    current.colourFlow[iy].second * diag[ix];
	
      }//end of diagram loop
      for(size_t ii = 0; ii < ncf; ++ii)
	for(size_t ij = 0; ij < ncf; ++ij)	  
	  full_me += cfactors[ii][ij]*(flows[ii]*conj(flows[ij])).real();
       
    }//vhel2 end
  }//vhel1 end
  DVector save(ndiags);
  const double ifact = mePartonData()[2]->id() == mePartonData()[3]->id() ?
    0.5 : 1.;
  // contact must have been set during the loop run
  assert(contact != ndiags);
  me[contact] = 0.;
  for(DVector::size_type ix = 0; ix < ndiags; ++ix)
    save[ix] = ifact*me[ix]/256.;
  meInfo(save);
  return ifact*full_me/256.;
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


void MEvv2ss::persistentOutput(PersistentOStream &) const {
}

void MEvv2ss::persistentInput(PersistentIStream &, int) {
}

ClassDescription<MEvv2ss> MEvv2ss::initMEvv2ss;
// Definition of the static class description member.

void MEvv2ss::Init() {

  static ClassDocumentation<MEvv2ss> documentation
    ("This class implements the ME for the vector-vector to scalar-scalar "
     "hard-process");

}

