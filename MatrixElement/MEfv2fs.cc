// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEfv2fs class.
//

#include "MEfv2fs.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig++/Helicity/Vertex/Vector/FFVVertex.h"
#include "Herwig++/Helicity/Vertex/Scalar/FFSVertex.h"
#include "Herwig++/Helicity/Vertex/Scalar/VSSVertex.h"
#include "ThePEG/Handlers/StandardXComb.h"

using namespace Herwig;
using Herwig::Helicity::SpinorWaveFunction;
using Herwig::Helicity::SpinorBarWaveFunction;
using Herwig::Helicity::ScalarWaveFunction;
using Herwig::Helicity::VectorWaveFunction;
using Herwig::Helicity::FFVVertexPtr;
using Herwig::Helicity::FFSVertexPtr;
using Herwig::Helicity::VSSVertex;
using Herwig::Helicity::incoming;
using Herwig::Helicity::outgoing;

typedef Ptr<VSSVertex>::pointer VSSVertexPtr;

double MEfv2fs::me2() const {
  bool cc = mePartonData()[0]->id() < 0;
  SpinorWaveFunction sp;
  SpinorBarWaveFunction spbar;
  if(cc) {
    spbar = SpinorBarWaveFunction(meMomenta()[0], mePartonData()[0],
				  incoming);
    sp = SpinorWaveFunction(meMomenta()[2], mePartonData()[2],
			    outgoing);
  }
  else {
    sp = SpinorWaveFunction(meMomenta()[0], mePartonData()[0],
			    incoming);
    spbar = SpinorBarWaveFunction(meMomenta()[2], mePartonData()[2],
				  outgoing);
  }
  VectorWaveFunction vecIn(meMomenta()[1], mePartonData()[1], incoming);  
  ScalarWaveFunction scaOut(meMomenta()[3], mePartonData()[3],
			    Complex(1.,0.), outgoing);
  double clrAvg(1./3.);
  if(mePartonData()[0]->iColour() == PDT::Colour0) clrAvg = 1.;
  const Energy2 m2(scale());
  const HPCount ndiags(numberOfDiags());
  const size_t ncf(numberOfFlows());
  const vector<vector<double> > cfactors = getColourFactors();
  vector<double> me(ndiags, 0.);
  vector<Complex> diag(ndiags, Complex(0.)), flows(ncf, Complex(0.));
  double full_me(0.);
  //intermediate wavefunctions
  SpinorWaveFunction interF; ScalarWaveFunction interS;
  SpinorBarWaveFunction interFB;
  FFVVertexPtr ffv; FFSVertexPtr ffs;
  VSSVertexPtr vss;
  for(unsigned int fhel1 = 0; fhel1 < 2; ++fhel1) {
    sp.reset(fhel1);
    for(unsigned int vhel1 = 0; vhel1 < 3; vhel1 += 2) {
      vecIn.reset(vhel1);
      for(unsigned int fhel2 = 0; fhel2 < 2; ++fhel2) {
	spbar.reset(fhel2);
	flows = vector<Complex>(ncf, Complex(0.));
	for(HPCount ix = 0; ix < ndiags; ++ix) {
	  HPDiagram current = getProcessInfo()[ix];
	  tcPDPtr offshell = current.intermediate;
	  if(current.channelType == HPDiagram::tChannel &&
	     offshell->iSpin() == PDT::Spin1Half) {  
	    if(dynamic_ptr_cast<FFVVertexPtr>(current.vertices.first)) {
	      ffv = dynamic_ptr_cast<FFVVertexPtr>(current.vertices.first);
	      ffs = dynamic_ptr_cast<FFSVertexPtr>
		(current.vertices.second);
	    }
	    else {
	      ffs = dynamic_ptr_cast<FFSVertexPtr>(current.vertices.first);
	      ffv = dynamic_ptr_cast<FFVVertexPtr>
		(current.vertices.second);
	    }
	    if(cc) {
	      interF = ffv->evaluate(m2, 3, offshell, sp, vecIn);
	      diag[ix] = ffs->evaluate(m2, interF, spbar, scaOut);
	    }
	    else {
	      interFB = ffv->evaluate(m2, 3, offshell, spbar, vecIn);
	      diag[ix] = ffs->evaluate(m2, sp, interFB, scaOut);
	    }
	  }
	  else if(current.channelType == HPDiagram::tChannel && 
		  offshell->iSpin() == PDT::Spin0) {
	    if(dynamic_ptr_cast<FFSVertexPtr>(current.vertices.first)) {
	      ffs = dynamic_ptr_cast<FFSVertexPtr>(current.vertices.first);
	      vss = dynamic_ptr_cast<VSSVertexPtr>
		(current.vertices.second);
	    }
	    else {
	      vss = dynamic_ptr_cast<VSSVertexPtr>(current.vertices.first);
	      ffs = dynamic_ptr_cast<FFSVertexPtr>
		(current.vertices.second);
	    }
	    interS = ffs->evaluate(m2, 3, offshell, sp, spbar);
	    if(cc)
	      diag[ix] = vss->evaluate(m2, vecIn, interS, scaOut);
	    else
	      diag[ix] = vss->evaluate(m2, vecIn, scaOut, interS);
	  }
	  else {
	    if(dynamic_ptr_cast<FFVVertexPtr>(current.vertices.first)) {
	      ffv = dynamic_ptr_cast<FFVVertexPtr>(current.vertices.first);
	      ffs = dynamic_ptr_cast<FFSVertexPtr>
		(current.vertices.second);
	    }
	    else {
	      ffs = dynamic_ptr_cast<FFSVertexPtr>(current.vertices.first);
	      ffv = dynamic_ptr_cast<FFVVertexPtr>
		(current.vertices.second);
	    }
	    if(cc) {
	      interFB = ffv->evaluate(m2, 5, offshell, spbar, vecIn);
	      diag[ix] = ffs->evaluate(m2, sp, interFB, scaOut);
	    }
	    else {
	      interF = ffv->evaluate(m2, 5, offshell, sp, vecIn);
	      diag[ix] = ffs->evaluate(m2, interF, spbar, scaOut);
	    }
	  }
	  me[ix] += norm(diag[ix]);

	  for(size_t iy = 0; iy < current.colourFlow.size(); ++iy)
	    flows[current.colourFlow[iy].first - 1] += 
	      current.colourFlow[iy].second * diag[ix];
	  
	}//end of diag loop
	for(size_t ii = 0; ii < ncf; ++ii)
	  for(size_t ij = 0; ij < ncf; ++ij) 
	    full_me += cfactors[ii][ij]*(flows[ii]*conj(flows[ij])).real();

      }
    }
  }
  DVector save(ndiags);
  for(DVector::size_type ix = 0; ix < ndiags; ++ix)
    save[ix] = clrAvg*me[ix]/32.;
  meInfo(save);
  full_me *= clrAvg/32.;
  return full_me;
}

Selector<const ColourLines *>
MEfv2fs::colourGeometries(tcDiagPtr diag) const {
  static vector<ColourLines> cl(10);
  //38->83
  cl[0] = ColourLines("1 4, -4 2 -3, 3 5");
  cl[1] = ColourLines("1 -2, 2 3 4, 5 -4 ");
  cl[2] = ColourLines("1 2 -3, -4 -2  5, 3 4");
  //3b8->83b
  cl[3] = ColourLines("-1 -4, 3 2 4, -3 -5");
  cl[4] = ColourLines("-1 2, -4 -3 -2, 4 -5");
  cl[5] = ColourLines("-1 -2 3, -5 2 4, -3 -4");
  //38->13
  cl[6] = ColourLines("1 2 -3, 3 5");
  cl[7] = ColourLines("1-2, 2 3 5");
  //3b8->13b
  cl[8] = ColourLines("-1 2 3, -3 -5");
  cl[9] = ColourLines("-1 2, -5 -3 -2");
  vector<ColourLines>::size_type offset;
  if(mePartonData()[0]->id() > 0 && 
     mePartonData()[2]->iColour() == PDT::Colour8 ) offset = 0;
  else if(mePartonData()[0]->id() < 0 && 
	  mePartonData()[2]->iColour() == PDT::Colour8 ) offset = 3;
  else if(mePartonData()[0]->id() > 0 && 
	  mePartonData()[2]->iColour() == PDT::Colour0 ) offset = 6;
  else offset = 8;
  HPDiagram current = getProcessInfo().at(abs(diag->id()) - 1); 
  Selector<const ColourLines *> sel;
  if(current.channelType == HPDiagram::tChannel && 
     (current.intermediate->iColour() == PDT::Colour3 ||
      current.intermediate->iColour() == PDT::Colour3bar))
    sel.insert(1., &cl[offset]);
  else if(current.channelType == HPDiagram::sChannel)
    sel.insert(1., &cl[offset + 1]);
  else
    sel.insert(1., &cl[offset + 2]);
  return sel;

}

void MEfv2fs::persistentOutput(PersistentOStream &) const {
}

void MEfv2fs::persistentInput(PersistentIStream &, int) {
}

ClassDescription<MEfv2fs> MEfv2fs::initMEfv2fs;
// Definition of the static class description member.

void MEfv2fs::Init() {

  static ClassDocumentation<MEfv2fs> documentation
    ("This class implements the matrix element for a fermion-vector to "
     "a fermioin-scalar.");

}

