// -*- C++ -*-
//
// MEvv2vv.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
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
#include "ThePEG/Helicity/SpinInfo.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "Herwig++/MatrixElement/HardVertex.h"
#include <numeric>

using namespace Herwig;
using ThePEG::Helicity::ScalarWaveFunction;
using ThePEG::Helicity::TensorWaveFunction;
using ThePEG::Helicity::incoming;
using ThePEG::Helicity::outgoing;
using ThePEG::Helicity::SpinfoPtr;

void MEvv2vv::doinit() throw(InitException) {
  GeneralHardME::doinit();
  size_t ndiags = numberOfDiags();
  theScaV.resize(ndiags);
  theVecV.resize(ndiags);
  theTenV.resize(ndiags);
  for(size_t i = 0; i < ndiags; ++i) {
    HPDiagram diag = getProcessInfo()[i];
    tcPDPtr offshell = diag.intermediate;
    if(!offshell)
      theFPVertex = dynamic_ptr_cast<AbstractVVVVVertexPtr>
	(diag.vertices.first);
    else if(offshell->iSpin() == PDT::Spin0) {
      AbstractVVSVertexPtr vert1 = dynamic_ptr_cast<AbstractVVSVertexPtr>
	(diag.vertices.first);
      AbstractVVSVertexPtr vert2 = dynamic_ptr_cast<AbstractVVSVertexPtr>
	(diag.vertices.second);
      theScaV[i] = make_pair(vert1, vert2);
    }
    else if(offshell->iSpin() == PDT::Spin1) {
      AbstractVVVVertexPtr vert1 = dynamic_ptr_cast<AbstractVVVVertexPtr>
	(diag.vertices.first);
      AbstractVVVVertexPtr vert2 = dynamic_ptr_cast<AbstractVVVVertexPtr>
	(diag.vertices.second);
      theVecV[i] = make_pair(vert1, vert2);
    }
    else if(offshell->iSpin() == PDT::Spin2) {
      AbstractVVTVertexPtr vert1 = dynamic_ptr_cast<AbstractVVTVertexPtr>
	(diag.vertices.first);
      AbstractVVTVertexPtr vert2 = dynamic_ptr_cast<AbstractVVTVertexPtr>
	(diag.vertices.second);
      theTenV[i] = make_pair(vert1, vert2);
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
  bool mc  = !(mePartonData()[2]->mass() > 0.*MeV);
  //massive vector, also 1
  if( !mc )
    vc[1] = VectorWaveFunction(rescaledMomenta()[2], mePartonData()[2], 1, 
			       outgoing);
  bool md  = !(mePartonData()[3]->mass() > 0.*MeV);
  if( !md ) 
    vd[1] = VectorWaveFunction(rescaledMomenta()[3], mePartonData()[3], 1, 
			       outgoing);
  double full_me(0.);
  vv2vvHeME(va, vb, vc, mc, vd, md, full_me);

#ifndef NDEBUG
  if( debugME() ) debug(full_me);
#endif

  return full_me;
}

ProductionMatrixElement 
MEvv2vv::vv2vvHeME(VBVector & vin1, VBVector & vin2, 
		   VBVector & vout1, bool mc, VBVector & vout2, bool md,
		   double & mesq) const {
  const HPCount ndiags = numberOfDiags();
  const size_t ncf(numberOfFlows());
  const vector<vector<double> > cfactors(getColourFactors());
  const Energy2 q2(scale());
  const Energy mass = vout1[0].mass();
  vector<Complex> diag(ndiags, Complex(0.));
  vector<double> me(ndiags, 0.);
  ScalarWaveFunction interS; VectorWaveFunction interV;
  TensorWaveFunction interT;
  ProductionMatrixElement prodME(PDT::Spin1, PDT::Spin1, PDT::Spin1, PDT::Spin1);
  for(unsigned int ihel1 = 0; ihel1 < 2; ++ihel1) { 
    for(unsigned int ihel2 = 0; ihel2 < 2; ++ihel2) {
      for(unsigned int ohel1 = 0; ohel1 < 3; ++ohel1) {
	if(mc && ohel1 == 1) ++ohel1;
	for(unsigned int ohel2 = 0; ohel2 < 3; ++ohel2) {
	  if(md && ohel2 == 1) ++ohel2;
	  vector<Complex> cflows(ncf, Complex(0.0));
	  for(HPCount ix = 0; ix < ndiags; ++ix) {
	    HPDiagram current = getProcessInfo()[ix];
	    tcPDPtr offshell = current.intermediate;
	    if(!offshell) continue;
	    if(current.channelType == HPDiagram::sChannel) {
	      if(offshell->iSpin() == PDT::Spin1) {
		interV = theVecV[ix].first->evaluate(q2, 1, offshell,
						     vin1[ihel1], vin2[ihel2]);
		diag[ix] = theVecV[ix].second->evaluate(q2, vout1[ohel1],
							vout2[ohel2], interV);
		diag[ix] += theFPVertex->evaluate(q2, 0, vout1[ohel1], vin2[ihel2], 
 						  vout2[ohel2], vin1[ihel1]);
	      }
	      else if(offshell->iSpin() == PDT::Spin2) {
		interT = theTenV[ix].first->evaluate(q2, 1, offshell,
						     vin1[ihel1], vin2[ihel2]);
		diag[ix] = theTenV[ix].second->evaluate(q2, vout1[ohel1], 
							vout2[ohel2],interT);
	      }
	    }
	    else if(current.channelType == HPDiagram::tChannel) {
	      if(offshell->iSpin() == PDT::Spin1) {
		if(current.ordered.second) {
		  interV = theVecV[ix].first->evaluate(q2, 3, offshell, vin1[ihel1],
						       vout1[ohel1], mass);
		  diag[ix] = theVecV[ix].second->evaluate(q2, vin2[ihel2], interV, 
							  vout2[ohel2]);
		  diag[ix] += theFPVertex->evaluate(q2, 0, vin1[ihel1], vin2[ihel2], 
						    vout1[ohel1], vout2[ohel2]);
		}
		else {
		  interV = theVecV[ix].first->evaluate(q2, 3, offshell, vin2[ihel2],
						       vout1[ohel1], mass);
		  diag[ix] = theVecV[ix].second->evaluate(q2, vin1[ihel1], interV, 
							  vout2[ohel2]);
		  diag[ix] += theFPVertex->evaluate(q2, 0, vin2[ihel2], vin1[ihel1],
						    vout1[ohel1], vout2[ohel2]);
		}
	      }
	      else if(offshell->iSpin() == PDT::Spin2) {
		if(current.ordered.second) {
		  interT = theTenV[ix].first->evaluate(q2, 3, offshell, vin1[ihel1],
						       vout1[ohel1], mass);
		  diag[ix] = theTenV[ix].second->evaluate(q2, vin2[ihel2], 
							  vout2[ohel2], interT);
		}
		else {
		  interT = theTenV[ix].first->evaluate(q2, 3, offshell, vin2[ihel2],
						       vout1[ohel1], mass);
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
	  //set the appropriate element in ProductionMatrixElement
	  prodME(ihel1, ihel2, ohel1, ohel2) = 
	    std::accumulate(cflows.begin(), cflows.end(), Complex(0.0, 0.0));
	  
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
  static vector<ColourLines> colourflows(16);
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
  //88->0->00
  colourflows[15] = ColourLines("1 -2,2 -1");

  HPDiagram current = getProcessInfo()[abs(diag->id()) - 1];
  Selector<const ColourLines *> select;
  if(current.channelType == HPDiagram::sChannel) {
    if(current.intermediate->iColour() == PDT::Colour8) {
      select.insert(0.25, &colourflows[0]);
      select.insert(0.25, &colourflows[1]);
      select.insert(0.25, &colourflows[2]);
      select.insert(0.25, &colourflows[3]);
    }
    else {
      if(getParticleData(current.outgoing.first)->iColour() == PDT::Colour0)
	select.insert(1., &colourflows[15]);
      else
	select.insert(1., &colourflows[12]);
    }
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

void MEvv2vv::constructVertex(tSubProPtr sub) {
  ParticleVector ext(4);
  ext[0] = sub->incoming().first;
  ext[1] = sub->incoming().second;
  ext[2] = sub->outgoing()[0];
  ext[3] = sub->outgoing()[1];

  //Ensure particles are in the same order as specified in the diagrams
  if( ext[0]->id() != getIncoming().first ) swap(ext[0], ext[1]);
  if( ext[2]->id() != getOutgoing().first ) swap(ext[2], ext[3]);

  VBVector v1, v2, v3, v4;
  VectorWaveFunction(v1, ext[0], incoming, false, true, true);
  VectorWaveFunction(v2, ext[1], incoming, false, true, true);
  //function to calculate me2 expects massless incoming vectors
  // and this constructor sets the '1' polarisation at element [2] 
  //in the vector
  bool mc  = !(ext[2]->data().mass() > 0.*MeV);
  bool md  = !(ext[3]->data().mass() > 0.*MeV);
  VectorWaveFunction(v3, ext[2], outgoing, true, mc, true);
  VectorWaveFunction(v4, ext[3], outgoing, true, md, true);
    
  double dummy(0.);

  //Need to use rescale momenta to calculate matrix element
  cPDVector data(4);
  vector<Lorentz5Momentum> momenta(4);
  for( size_t i = 0; i < 4; ++i ) {
    data[i] = ext[i]->dataPtr();
    momenta[i] = ext[i]->momentum();
  }
  rescaleMomenta(momenta, data);

  VectorWaveFunction vr1(rescaledMomenta()[0], data[0], incoming);
  VectorWaveFunction vr2(rescaledMomenta()[1], data[1], incoming);
  VectorWaveFunction vr3(rescaledMomenta()[2], data[2], outgoing);
  VectorWaveFunction vr4(rescaledMomenta()[3], data[3], outgoing);

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
  ProductionMatrixElement pme = vv2vvHeME(v1, v2, v3, mc, v4, md, dummy);

#ifndef NDEBUG
  if( debugME() ) debug(dummy);
#endif

  HardVertexPtr hv = new_ptr(HardVertex());
  hv->ME(pme);
  for(unsigned int i = 0; i < 4; ++i) 
    dynamic_ptr_cast<SpinfoPtr>(ext[i]->spinInfo())->setProductionVertex(hv);
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
