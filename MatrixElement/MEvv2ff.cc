// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEvv2ff class.
//

#include "MEvv2ff.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"
#include "ThePEG/Helicity/Vertex/Vector/VVVVertex.h"
#include "ThePEG/Helicity/Vertex/Scalar/VVSVertex.h"
#include "ThePEG/Helicity/Vertex/Scalar/FFSVertex.h"
#include "ThePEG/Helicity/Vertex/Tensor/VVTVertex.h"
#include "ThePEG/Helicity/Vertex/Tensor/FFTVertex.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

using namespace Herwig;
using ThePEG::Helicity::SpinorWaveFunction;
using ThePEG::Helicity::SpinorBarWaveFunction;
using ThePEG::Helicity::VectorWaveFunction;
using ThePEG::Helicity::ScalarWaveFunction;
using ThePEG::Helicity::TensorWaveFunction;
using ThePEG::Helicity::FFVVertexPtr;
using ThePEG::Helicity::VVVVertexPtr;
using ThePEG::Helicity::VVSVertexPtr;
using ThePEG::Helicity::FFSVertexPtr;
using ThePEG::Helicity::VVTVertexPtr;
using ThePEG::Helicity::FFTVertexPtr;
using ThePEG::Helicity::incoming;
using ThePEG::Helicity::outgoing;

double MEvv2ff::me2() const {
  // Set up wavefuctions
  VectorWaveFunction vec1(meMomenta()[0],mePartonData()[0],incoming);
  VectorWaveFunction vec2(meMomenta()[1],mePartonData()[1],incoming);
  SpinorBarWaveFunction spbar(meMomenta()[2], mePartonData()[2],
			      outgoing);
  SpinorWaveFunction sp(meMomenta()[3], mePartonData()[3], outgoing);
  //Define constants

  const HPCount ndiags(numberOfDiags());
  const size_t ncf(numberOfFlows());
  const vector<vector<double> > cfactors = getColourFactors();
  //vectors to store results
  vector<double> me(ndiags, 0.), sumflow(ncf, 0.);
  vector<Complex> diag(ndiags, Complex(0.)), flows(ncf, Complex(0.));
  double full_me(0.);
  //offshell wavefunctions
  SpinorWaveFunction interF; VectorWaveFunction interV;
  ScalarWaveFunction interS; TensorWaveFunction interT;
  const Energy2 m2 = scale();
  //sum over vector helicities
  for(unsigned int vhel1 = 0; vhel1 < 3; vhel1 += 2) {
    vec1.reset(vhel1);
    for(unsigned int vhel2 = 0; vhel2 < 3; vhel2 += 2) {
      vec2.reset(vhel2);
      //sum over fermion helicities
      for(unsigned int fhel1 = 0; fhel1 < 2; ++fhel1) {
	spbar.reset(fhel1);
	for(unsigned int fhel2 = 0; fhel2 < 2; ++ fhel2) {
	  sp.reset(fhel2); 
	  //for each helicity calculate diagram 
	  flows = vector<Complex>(ncf, Complex(0,0.));
	  for(HPCount ix = 0; ix < ndiags; ++ix) {
	    HPDiagram current = getProcessInfo()[ix];
	    PDPtr offshell = current.intermediate;
	    if(current.channelType == HPDiagram::tChannel && 
	       offshell->iSpin() == PDT::Spin1Half) {
	      FFVVertexPtr vert1 = dynamic_ptr_cast<FFVVertexPtr>
		(current.vertices.first);
	      FFVVertexPtr vert2 = dynamic_ptr_cast<FFVVertexPtr>
		(current.vertices.second);
	      if(current.ordered.second) {
		interF = vert2->evaluate(m2, 3, offshell, sp, vec2);
		diag[ix] = vert1->evaluate(m2, interF, spbar, vec1);
	      }
	      else {
		interF = vert2->evaluate(m2, 3, offshell, sp, vec1);
		diag[ix] = vert1->evaluate(m2, interF, spbar, vec2);
	      }
	    }
	    else if(current.channelType == HPDiagram::sChannel && 
		    offshell->iSpin() == PDT::Spin1) {
	      VVVVertexPtr vert1 = dynamic_ptr_cast<VVVVertexPtr>
		(current.vertices.first);
	      FFVVertexPtr vert2 = dynamic_ptr_cast<FFVVertexPtr>
		(current.vertices.second);
	      interV = vert1->evaluate(m2, 1, offshell, vec1, vec2);
	      diag[ix] = vert2->evaluate(m2, sp, spbar, interV);
	    }
	    else if(current.channelType == HPDiagram::sChannel && 
		    offshell->iSpin() == PDT::Spin2) {
	      VVTVertexPtr vert1 = dynamic_ptr_cast<VVTVertexPtr>
		(current.vertices.first);
	      FFTVertexPtr vert2 = dynamic_ptr_cast<FFTVertexPtr>
		(current.vertices.second);
	      interT = vert1->evaluate(m2, 1, offshell, vec1, vec2);
	      diag[ix] = vert2->evaluate(m2, sp, spbar, interT);
	    }
	    else diag[ix] = 0.;
	    me[ix] += norm(diag[ix]);
	    for(size_t iy = 0; iy < current.colourFlow.size(); 
		++iy) {
	      flows[current.colourFlow[iy].first - 1] += 
		current.colourFlow[iy].second * diag[ix];
	      sumflow[iy] += norm(flows[iy]);
	    }
	  }//end of diag loop
	  for(size_t ii = 0; ii < ncf; ++ii)
	    for(size_t ij = 0; ij < ncf; ++ij)
	      full_me += cfactors[ii][ij]*(flows[ii]*conj(flows[ij])).real();
	
	}
      }
    }
  }
  const double identFact = mePartonData()[2]->id() == mePartonData()[3]->id()
    ? 0.5 : 1.;
  DVector save(ndiags);
  for(DVector::size_type ix = 0; ix < ndiags; ++ix)
    save[ix] = identFact*me[ix]/256.;
  meInfo(save);
  full_me *= identFact/256;

#ifndef NDEBUG 
  if( debugME() ) debug(full_me);
#endif

  return full_me;
}

Selector<const ColourLines *>
MEvv2ff::colourGeometries(tcDiagPtr diag) const {
  static vector<ColourLines> cf(17);
  //88->33bar
  cf[0] = ColourLines("1 4, -3 -5, 3 2 -1");
  cf[1] = ColourLines("1 -5, 1 2 -3, 3 4 ");
  cf[2] = ColourLines("2 -1, 1 3 4, -2 -3 -5");
  cf[3] = ColourLines("1 -2, -1 -3 -5, 2 3 4");
  cf[4] = ColourLines("1 -2, 2 -1, 4 -5");
  //88->88
  cf[5] = ColourLines("1 4, -3 -5, 3 -2 -1, -4 2 5");
  cf[6] = ColourLines("-1 -4, 3 5, -3 2 1, 4 -2 -5");
  cf[7] = ColourLines("1 4, 3 5, -3 2 -4, -1 -2 -5");
  cf[8] = ColourLines("-1 -4, -3 -5, 3 -2 4, 1 2 5");
  cf[9] = ColourLines("1 5, -3 -4, 3 -2 -1, -5 2 4");
  cf[10] = ColourLines("-1 -5, 3 4, -3 2 1, 5 -2 -4");
  cf[11] = ColourLines("1 5, 3 4, -3 2 -5, -1 -2 -4");
  cf[12] = ColourLines("-1 -5, -3 -4, 3 -2 5, 1 2 4");
  cf[13] = ColourLines("1 -2, 2 3 5, -1 -3 -4, -5 4");
  cf[14] = ColourLines("-1 2, -2 -3 -5, 1 3 4, 5 -4");
  cf[15] = ColourLines("-1 2, -2 -3 -4, 1 3 5, -5 4");
  cf[16] = ColourLines("1 -2, 2 3 4, -1 -3 -5, 5 -4");

  HPDiagram current = getProcessInfo()[abs(diag->id()) - 1];
  vector<ColourLines>::size_type offset(0);
  if(current.channelType == HPDiagram::tChannel && current.ordered.first)
    offset = 5;
  else if(current.channelType == HPDiagram::tChannel && !current.ordered.first)
    offset = 9;
  else 
    offset = 13;  
  Selector<const ColourLines *> sel;
  if(current.intermediate->iColour() == PDT::Colour0)
    sel.insert(1., &cf[4]);
  else if(getParticleData(current.outgoing.first)->iColour() == PDT::Colour8) {
    sel.insert(0.25, &cf[offset]);
    sel.insert(0.25, &cf[offset + 1]);
    sel.insert(0.25, &cf[offset + 2]);
    sel.insert(0.25, &cf[offset + 3]);
  }
  else {
    if(offset == 5 || offset == 9)
      sel.insert(1., &cf[(int)(offset/4) - 1]);
    else {
      sel.insert(0.5, &cf[2]);
      sel.insert(0.5, &cf[3]);
    }
  }    
  return sel;
}

void MEvv2ff::persistentOutput(PersistentOStream &) const {
}

void MEvv2ff::persistentInput(PersistentIStream &, int) {
}

ClassDescription<MEvv2ff> MEvv2ff::initMEvv2ff;
// Definition of the static class description member.

void MEvv2ff::Init() {

  static ClassDocumentation<MEvv2ff> documentation
    ("The MEvv2ff class handles the ME calculation for the general "
     "spin configuration vector-vector to fermion-antifermion\n.");

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
