// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEvv2vv class.
//

#include "MEvv2vv.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "ThePEG/Helicity/Vertex/Tensor/VVTVertex.h"
#include <numeric>

using namespace Herwig;
using ThePEG::Helicity::ScalarWaveFunction;
using ThePEG::Helicity::TensorWaveFunction;
using ThePEG::Helicity::VVTVertexPtr;
using ThePEG::Helicity::incoming;
using ThePEG::Helicity::outgoing;

double MEvv2vv::me2() const {
  VBVector va(2), vb(2), vc(3), vd(3);  
  for(unsigned int i = 0; i < 2; ++i) {
    va[i] = VectorWaveFunction(meMomenta()[0], mePartonData()[0], 2*i, incoming);
    vb[i] = VectorWaveFunction(meMomenta()[1], mePartonData()[1], 2*i, incoming);
  }
  //always 0 and 2 polarisations
  for(unsigned int i = 0; i < 2; ++i) {
    vc[2*i] = VectorWaveFunction(meMomenta()[2], mePartonData()[2], 2*i, outgoing);
    vd[2*i] = VectorWaveFunction(meMomenta()[3], mePartonData()[3], 2*i, outgoing);
  }
  //massive vector, also 1
  if(mePartonData()[2]->mass() > 0.0*MeV)
    vc[1] = VectorWaveFunction(meMomenta()[2], mePartonData()[2], 1, outgoing);
  if(mePartonData()[3]->mass() > 0.0*MeV)
    vd[1] = VectorWaveFunction(meMomenta()[3], mePartonData()[3], 1, outgoing);
  double full_me(0.);
  vv2vvHeME(va, vb, vc, vd, full_me);
  return full_me;
}

ProductionMatrixElement 
MEvv2vv::vv2vvHeME(VBVector & vin1, VBVector & vin2, 
		   VBVector & vout1, VBVector & vout2,
		   double & mesq) const {
  const HPCount ndiags = numberOfDiags();
  const size_t ncf(numberOfFlows());
  const vector<vector<double> > cfactors(getColourFactors());
  const Energy2 q2(scale());
  vector<Complex> diag(ndiags, Complex(0.));
  vector<double> me(ndiags, 0.);
  bool masslessC = (mePartonData()[2]->mass() == 0.0*MeV);
  bool masslessD = (mePartonData()[3]->mass() == 0.0*MeV);
  ScalarWaveFunction interS; VectorWaveFunction interV;
  TensorWaveFunction interT;
  ProductionMatrixElement prodME(PDT::Spin1, PDT::Spin1, PDT::Spin1, PDT::Spin1);
  for(unsigned int ihel1 = 0; ihel1 < 2; ++ihel1) { 
    for(unsigned int ihel2 = 0; ihel2 < 2; ++ihel2) {
      for(unsigned int ohel1 = 0; ohel1 < 3; ++ohel1) {
	if(masslessC && ohel1 == 1) ++ohel1;
	for(unsigned int ohel2 = 0; ohel2 < 3; ++ohel2) {
	  if(masslessD && ohel2 == 1) ++ohel2;
	  vector<Complex> cflows(ncf, Complex(0.0));
	  for(HPCount ix = 0; ix < ndiags; ++ix) {
	    HPDiagram current = getProcessInfo()[ix];
	    tcPDPtr offshell = current.intermediate;
	    if(!offshell) continue;
	    if(current.channelType == HPDiagram::sChannel) {
	      if(offshell->iSpin() == PDT::Spin1) {
		interV = theVecV[ix].first->evaluate(q2, 3, offshell,
						     vin1[ihel1], vin2[ihel2]);
		diag[ix] = theVecV[ix].second->evaluate(q2, vout1[ohel1],
							vout2[ohel2], interV);
		diag[ix] += theFPVertex->evaluate(q2, 0, vout1[ohel1], vin2[ihel2], 
 						  vout2[ohel2], vin1[ihel1]);
	      }
	      else if(offshell->iSpin() == PDT::Spin2) {
		interT = theTenV[ix].first->evaluate(q2, 3, offshell,
						     vin1[ihel1], vin2[ihel2]);
		diag[ix] = theTenV[ix].second->evaluate(q2, vout1[ohel1], 
							vout2[ohel2],interT);
	      }
	    }
	    else if(current.channelType == HPDiagram::tChannel) {
	      if(offshell->iSpin() == PDT::Spin1) {
		if(current.ordered.second) {
		  interV = theVecV[ix].first->evaluate(q2, 3, offshell, vin1[ihel1],
						       vout1[ohel1]);
		  diag[ix] = theVecV[ix].second->evaluate(q2, vin2[ihel2], interV, 
							  vout2[ohel2]);
		  diag[ix] += theFPVertex->evaluate(q2, 0, vin1[ihel1], vin2[ihel2], 
						    vout1[ohel1], vout2[ohel2]);
		}
		else {
		  interV = theVecV[ix].first->evaluate(q2, 3, offshell, vin2[ihel2],
						       vout1[ohel1]);
		  diag[ix] = theVecV[ix].second->evaluate(q2, vin1[ihel1], interV, 
							  vout2[ohel2]);
		  diag[ix] += theFPVertex->evaluate(q2, 0, vin2[ihel2], vin1[ihel1],
						    vout1[ohel1], vout2[ohel2]);
		}
	      }
	      else if(offshell->iSpin() == PDT::Spin2) {
		if(current.ordered.second) {
		  interT = theTenV[ix].first->evaluate(q2, 3, offshell, vin1[ihel1],
						       vout1[ohel1]);
		  diag[ix] = theTenV[ix].second->evaluate(q2, vin2[ihel2], 
							  vout2[ohel2], interT);
		}
		else {
		  interT = theTenV[ix].first->evaluate(q2, 3, offshell, vin2[ihel2],
						       vout1[ohel1]);
		  diag[ix] = theTenV[ix].second->evaluate(q2, vin1[ihel1], 
							  vout2[ohel2], interT);
		}
	      }
	    }
	    else 
	      diag[ix] = 0.0;
	    
	    me[ix] += norm(diag[ix]);
	    //Compute flows

	    for(size_t iy = 0; iy < current.colourFlow.size(); ++iy)
	      cflows[current.colourFlow[iy].first - 1] += 
		current.colourFlow[iy].second * diag[ix];
	    
	  } //end of diagram loop

	  for(size_t ii = 0; ii < ncf; ++ii)
	    for(size_t ij = 0; ij < ncf; ++ij)
	      mesq += cfactors[ii][ij]*(cflows[ii]*conj(cflows[ij])).real();
	
	}
      }
    }
  }
  const double identfact = mePartonData()[2]->id() == mePartonData()[3]->id() 
    ? 0.5 : 1.;
  const double colfact = mePartonData()[0]->iColour() == PDT::Colour8 ? 
    1./64. : 1;
  DVector save(ndiags);
  for(DVector::size_type ix = 0; ix < ndiags; ++ix)
    save[ix] = 0.25*identfact*colfact*me[ix];
  meInfo(save);
  mesq *= 0.25*identfact*colfact;
  return prodME;
}

Selector<const ColourLines *>
MEvv2vv::colourGeometries(tcDiagPtr diag) const {
  static vector<ColourLines> colourflows(15);
  //88->8->88
  colourflows[0] = ColourLines("1 -2, -1 -3 -4, 4 -5, 2 3 5");
  colourflows[1] = ColourLines("-1 2, 1 3 4, -4 5, -2 -3 -5");
  colourflows[2] = ColourLines("1 -2, -1 -3 -5, 5 -4, 2 3 4");
  colourflows[3] = ColourLines("-1 2, 1 3 5, -5 4, -2 -3 -4");
  colourflows[4] = ColourLines("1 4, -1 -2 3, -3 -5, -4 2 5");
  colourflows[5] = ColourLines("-1 -4, 1 2 -3, 3 5, 4 -2 -5");
  colourflows[6] = ColourLines("1 4, -1 -2 -5, 3 5, -3 2 -4");
  colourflows[7] = ColourLines("-1 -4, 1 2 5, -3 -5, 3 -2 4");
  colourflows[8] = ColourLines("1 5, -1 -2 3, -3 -4, -5 2 4");
  colourflows[9] = ColourLines("-1 -5, 1 2 -3, 3 4, 5 -2 -4");
  colourflows[10] = ColourLines("1 5, -1 -2 -4, 3 4, -3 2 -5");
  colourflows[11] = ColourLines("-1 -5, 1 2 4, -3 -4, 3 -2 5");
  //88->0->88
  colourflows[12] = ColourLines("1 -2, 2 -1, 4 -5, 5 -4");
  colourflows[13] = ColourLines("1 4, -1 -4, 3 5, -5 -3");
  colourflows[14] = ColourLines("1 5, -1 -5, 3 4, -3 -4");
  
  HPDiagram current = getProcessInfo()[abs(diag->id()) - 1];
  Selector<const ColourLines *> select;
  if(current.channelType == HPDiagram::sChannel) {
    if(current.intermediate->iColour() == PDT::Colour8) {
      select.insert(0.25, &colourflows[0]);
      select.insert(0.25, &colourflows[1]);
      select.insert(0.25, &colourflows[2]);
      select.insert(0.25, &colourflows[3]);
    }
    else
      select.insert(1., &colourflows[12]);
  }
  else if(current.channelType == HPDiagram::tChannel) {
    if(current.ordered.second) {
      if(current.intermediate->iColour() == PDT::Colour8) {
	select.insert(0.25, &colourflows[4]);
	select.insert(0.25, &colourflows[5]);
	select.insert(0.25, &colourflows[6]);
	select.insert(0.25, &colourflows[7]);
      }
      else 
	select.insert(1., &colourflows[13]);
    }
    else {
      if(current.intermediate->iColour() == PDT::Colour8) {
      select.insert(0.25, &colourflows[8]);
      select.insert(0.25, &colourflows[9]);
      select.insert(0.25, &colourflows[10]);
      select.insert(0.25, &colourflows[11]);
      }
      else
	select.insert(1., &colourflows[14]);
    }      
  }
  else 
    throw MEException() << "MEvv2vv::colourGeometries - Trying to set ColourLines "
			<< "for an unknown diagram type. " 
			<< Exception::warning;
  return select;
}


void MEvv2vv::persistentOutput(PersistentOStream & os) const {
  os << theScaV << theVecV << theTenV << theFPVertex;
}

void MEvv2vv::persistentInput(PersistentIStream & is, int) {
  is >> theScaV >> theVecV >> theTenV >> theFPVertex;
}

ClassDescription<MEvv2vv> MEvv2vv::initMEvv2vv;
// Definition of the static class description member.

void MEvv2vv::Init() {

  static ClassDocumentation<MEvv2vv> documentation
    ("This is the implementation of the 2 to 2 ME for a pair"
     "of massless vector-bosons to a pair of vector bosons");

}

