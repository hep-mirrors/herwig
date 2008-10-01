// -*- C++ -*-
//
// MEfv2fs.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEfv2fs class.
//

#include "MEfv2fs.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/MatrixElement/HardVertex.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/Helicity/SpinInfo.h"
#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"
#include "ThePEG/Helicity/Vertex/Scalar/FFSVertex.h"
#include "ThePEG/Helicity/Vertex/Scalar/VSSVertex.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include <numeric>

using namespace Herwig;
using ThePEG::Helicity::SpinfoPtr;
using ThePEG::Helicity::incoming;
using ThePEG::Helicity::outgoing;

void MEfv2fs::doinit() throw(InitException) {
 GeneralHardME::doinit();
 size_t ndiags(numberOfDiags());
 theScaV.resize(ndiags);
 theFermV.resize(ndiags);
 for(size_t ix = 0; ix < ndiags; ++ix) {
   HPDiagram curr = getProcessInfo()[ix];
   if(curr.channelType == HPDiagram::tChannel) {
     AbstractFFSVertexPtr ffs = 
       dynamic_ptr_cast<AbstractFFSVertexPtr>(curr.vertices.first);
     if( curr.intermediate->iSpin() == PDT::Spin0 ) {
       AbstractVSSVertexPtr vss = 
	 dynamic_ptr_cast<AbstractVSSVertexPtr>(curr.vertices.second);
       theScaV[ix] =  make_pair(ffs, vss); 
     }
     else {
       AbstractFFVVertexPtr ffv = 
	 dynamic_ptr_cast<AbstractFFVVertexPtr>(curr.vertices.second);
       theFermV[ix] = make_pair(ffs, ffv);
     }
   }
   else {
     AbstractFFVVertexPtr ffv = 
       dynamic_ptr_cast<AbstractFFVVertexPtr>(curr.vertices.first);
     AbstractFFSVertexPtr ffs = 
       dynamic_ptr_cast<AbstractFFSVertexPtr>(curr.vertices.second);
     theFermV[ix] = make_pair(ffs, ffv); 
   }
 }

}


double MEfv2fs::me2() const {
  //massless vector
  VecWFVector vecIn(2);
  ScalarWaveFunction scaOut(rescaledMomenta()[3], mePartonData()[3],
			    Complex(1.,0.), outgoing);
  double fullme(0.);
  if( mePartonData()[0]->id() > 0 ) {
    SpinorVector spIn(2);
    SpinorBarVector spbOut(2);
    for(size_t ih = 0; ih < 2; ++ih) {
      spIn[ih] = SpinorWaveFunction(rescaledMomenta()[0], mePartonData()[0], ih,
				    incoming);
      spbOut[ih] = SpinorBarWaveFunction(rescaledMomenta()[2], mePartonData()[2], ih,
					 outgoing);
      vecIn[ih] = VectorWaveFunction(rescaledMomenta()[1], mePartonData()[1], 2*ih,
				     incoming);
    }
    fv2fbsHeME(spIn, vecIn, spbOut, scaOut, fullme);
  }
  else {
    SpinorBarVector spbIn(2);
    SpinorVector spOut(2);
    for(size_t ih = 0; ih < 2; ++ih) {
      spbIn[ih] = SpinorBarWaveFunction(rescaledMomenta()[0], mePartonData()[0], ih,
					incoming);
      spOut[ih] = SpinorWaveFunction(rescaledMomenta()[2], mePartonData()[2], ih,
				     outgoing);
      vecIn[ih] = VectorWaveFunction(rescaledMomenta()[1], mePartonData()[1], 2*ih,
				     incoming);
    }
    fbv2fsHeME(spbIn, vecIn, spOut, scaOut, fullme);
  }

#ifndef NDEBUG
  if( debugME() ) debug(fullme);
#endif

  return fullme;
}

ProductionMatrixElement 
MEfv2fs::fv2fbsHeME(const SpinorVector & spIn, const VecWFVector & vecIn,
		    const SpinorBarVector & spbOut,
		    const ScalarWaveFunction & scaOut,
		    double & fullme) const {
  const Energy2 q2(scale());
  const HPCount ndiags(numberOfDiags());
  const size_t ncf(numberOfFlows());
  const vector<vector<double> > cfactors = getColourFactors();
  vector<double> me(ndiags, 0.);
  vector<Complex> diag(ndiags, Complex(0.));
  fullme = 0.;
  //intermediate wave functions
  SpinorWaveFunction interF; ScalarWaveFunction interS;
  SpinorBarWaveFunction interFB;
  //vertex pointers
  FFSVertexPtr ffs; FFVVertexPtr ffv;
  VSSVertexPtr vss;
  ProductionMatrixElement prodME(PDT::Spin1Half, PDT::Spin1, 
				 PDT::Spin1Half, PDT::Spin0);
  for(unsigned int ihel1 = 0; ihel1 < 2; ++ihel1) {
    for(unsigned int ihel2 = 0; ihel2 < 2; ++ihel2) {
      for(unsigned int ohel1 = 0; ohel1 < 2; ++ohel1) {
	vector<Complex> flows = vector<Complex>(ncf, Complex(0.));
	for(size_t ix = 0; ix < ndiags; ++ix) {
	  HPDiagram current = getProcessInfo()[ix];
	  tcPDPtr offshell = current.intermediate;
	  if( current.channelType == HPDiagram::tChannel ) {
	    if( offshell->iSpin() == PDT::Spin0 ) {
	      interS = theScaV[ix].first->evaluate(q2, 3, offshell, spIn[ihel1], 
						   spbOut[ohel1]);
	      diag[ix] = theScaV[ix].second->evaluate(q2, vecIn[ihel2], 
						      scaOut, interS);
	    }
	    else if( offshell->iSpin() == PDT::Spin1Half ) {
	      interFB = theFermV[ix].second->evaluate(q2, 3, offshell, 
						      spbOut[ohel1], 
						      vecIn[ihel2]);
	      diag[ix] = theFermV[ix].first->evaluate(q2, spIn[ihel1], 
						      interFB, scaOut);
	    }
	    else diag[ix] = 0.0;
	  }
	  else if( current.channelType == HPDiagram::sChannel ) {
	    interF = theFermV[ix].second->evaluate(q2, 1, offshell, spIn[ihel1], 
						   vecIn[ihel2]);
	    diag[ix] = theFermV[ix].first->evaluate(q2, interF, spbOut[ohel1], 
						    scaOut);
	  }
	  else diag[ix] = 0.0;
	  // save 
	  me[ix] += norm(diag[ix]);
	  //add to correct flow
	  for(size_t iy = 0; iy < current.colourFlow.size(); ++iy)
	    flows[current.colourFlow[iy].first - 1] += 
	      current.colourFlow[iy].second * diag[ix];
	  
	}//end of diag loop
	for(size_t ii = 0; ii < ncf; ++ii)
	  for(size_t ij = 0; ij < ncf; ++ij) 
	    fullme += cfactors[ii][ij]*(flows[ii]*conj(flows[ij])).real();

	prodME(ihel1, 2*ihel2, ohel1, 0) = 
	  std::accumulate(flows.begin(), flows.end(), Complex(0., 0.));
	
      }
    }
  }
  double clrAvg = (mePartonData()[0]->iColour() == PDT::Colour3) ? 1./3. : 1.; 
  DVector save(ndiags);
  for(DVector::size_type ix = 0; ix < ndiags; ++ix)
    save[ix] = clrAvg*me[ix]/32.;
  meInfo(save);
  fullme *= clrAvg/32.;
  return prodME;
}

ProductionMatrixElement 
MEfv2fs::fbv2fsHeME(const SpinorBarVector & spbIn, const VecWFVector & vecIn,
		    const SpinorVector & spOut,
		    const ScalarWaveFunction & scaOut,
		    double & fullme) const {
  const Energy2 q2(scale());
  const HPCount ndiags(numberOfDiags());
  const size_t ncf(numberOfFlows());
  const vector<vector<double> > cfactors = getColourFactors();
  vector<double> me(ndiags, 0.);
  vector<Complex> diag(ndiags, Complex(0.));

  fullme = 0.;
  //intermediate wave functions
  ScalarWaveFunction interS;
  SpinorBarWaveFunction interFB;
  //vertex pointers
  ProductionMatrixElement prodME(PDT::Spin1Half, PDT::Spin1, PDT::Spin1Half,
				 PDT::Spin0);
  for(unsigned int ihel1 = 0; ihel1 < 2; ++ihel1) {
    for(unsigned int ihel2 = 0; ihel2 < 2; ++ihel2) {
      for(unsigned int ohel1 = 0; ohel1 < 2; ++ohel1) {
	vector<Complex> flows = vector<Complex>(ncf, Complex(0.));
	for(HPCount ix = 0; ix < ndiags; ++ix) {
	  HPDiagram current = getProcessInfo()[ix];
	  tcPDPtr offshell = current.intermediate;
	  if( current.channelType == HPDiagram::tChannel ) {
	    if( offshell->iSpin() == PDT::Spin0 ) {
	      interS = theScaV[ix].first->evaluate(q2, 3, offshell, 
						   spOut[ohel1],
						   spbIn[ihel1]);
	      diag[ix] = theScaV[ix].second->evaluate(q2, vecIn[ihel2], 
						     interS, scaOut);
	    }
	    else if( offshell->iSpin() == PDT::Spin1Half ) {
	      interFB = theFermV[ix].first->evaluate(q2, 3, offshell, 
						     spbIn[ihel1], scaOut);
	      diag[ix] = theFermV[ix].second->evaluate(q2, spOut[ohel1], 
						       interFB, vecIn[ihel2]);
	    }
	    else diag[ix] = 0.0;
	  }
	  else if( current.channelType == HPDiagram::sChannel ) {
	    interFB = theFermV[ix].second->evaluate(q2, 1, offshell, 
						    spbIn[ihel1], 
						    vecIn[ihel2]);
	    diag[ix] = theFermV[ix].first->evaluate(q2, spOut[ohel1], 
						    interFB, scaOut);
	  }
	  else diag[ix] = 0.0;
	  // save 
	  me[ix] += norm(diag[ix]);
	  //add to correct flow
	  for(size_t iy = 0; iy < current.colourFlow.size(); ++iy)
	    flows[current.colourFlow[iy].first - 1] += 
	      current.colourFlow[iy].second * diag[ix];
	  
	}//end of diag loop
	for(size_t ii = 0; ii < ncf; ++ii)
	  for(size_t ij = 0; ij < ncf; ++ij) 
	    fullme += cfactors[ii][ij]*(flows[ii]*conj(flows[ij])).real();

	prodME(ihel1, 2*ihel2, ohel1, 0) = 
	  std::accumulate(flows.begin(), flows.end(), Complex(0., 0.));
	
      }
    }
  }
  double clrAvg = (mePartonData()[0]->iColour() == PDT::Colour3bar) ? 1./3. : 1.; 
  DVector save(ndiags);
  for(DVector::size_type ix = 0; ix < ndiags; ++ix)
    save[ix] = clrAvg*me[ix]/32.;
  meInfo(save);
  fullme *= clrAvg/32.;
  return prodME;
}

Selector<const ColourLines *>
MEfv2fs::colourGeometries(tcDiagPtr diag) const {
  static vector<ColourLines> cl(14);
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
  //38->31
  cl[10] = ColourLines("1 2 -3, 3 4");
  cl[11] = ColourLines("1-2, 2 3 4");
  //3b8->3b1
  cl[12] = ColourLines("-1 2 3, -3 -4");
  cl[13] = ColourLines("-1 2, -4 -3 -2");
  vector<ColourLines>::size_type offset = 999;
  if(mePartonData()[0]->id() > 0 && 
     mePartonData()[2]->iColour() == PDT::Colour8 ) offset = 0;
  else if(mePartonData()[0]->id() < 0 && 
	  mePartonData()[2]->iColour() == PDT::Colour8 ) offset = 3;
  else if(mePartonData()[0]->id() > 0 && 
	  mePartonData()[2]->iColour() == PDT::Colour0 ) offset = 6;
  else if(mePartonData()[0]->id() < 0 && 
	  mePartonData()[2]->iColour() == PDT::Colour0 ) offset = 8;
  else if(mePartonData()[0]->id() > 0 && 
	  mePartonData()[3]->iColour() == PDT::Colour0 ) offset = 10;
  else if(mePartonData()[0]->id() < 0 && 
	  mePartonData()[3]->iColour() == PDT::Colour0 ) offset = 12;
  assert(offset != 999);
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

void MEfv2fs::persistentOutput(PersistentOStream & os) const {
  os << theScaV << theFermV;
}

void MEfv2fs::persistentInput(PersistentIStream & is, int) {
  is >> theScaV >> theFermV;
}

ClassDescription<MEfv2fs> MEfv2fs::initMEfv2fs;
// Definition of the static class description member.

void MEfv2fs::Init() {

  static ClassDocumentation<MEfv2fs> documentation
    ("This class implements the matrix element for a fermion-vector to "
     "a fermioin-scalar.");

}

void MEfv2fs::constructVertex(tSubProPtr subp) {
  //get external particles
  ParticleVector external(4);
  external[0] = subp->incoming().first; 
  external[1] = subp->incoming().second;
  external[2] = subp->outgoing()[0]; 
  external[3] = subp->outgoing()[1];

  //make sure the order is correct
  if( external[0]->dataPtr()->iSpin() > external[1]->dataPtr()->iSpin() ) 
    swap(external[0], external[1]);

  if( external[2]->dataPtr()->iSpin() < external[3]->dataPtr()->iSpin() ) 
    swap(external[2], external[3]);

  //calculate production ME
  VecWFVector vecIn;
  VectorWaveFunction(vecIn, external[1], incoming, false, true);
  //function to calculate me2 expects massless incoming vectors
  // and this constructor sets the '1' polarisation at element [2] 
  //in the vector
  vecIn[1] = vecIn[2];
  ScalarWaveFunction scaOut(external[3], outgoing, true);

  //Need to use rescale momenta to calculate matrix element
  cPDVector data(4);
  vector<Lorentz5Momentum> momenta(4);
  for( size_t i = 0; i < 4; ++i ) {
    data[i] = external[i]->dataPtr();
    momenta[i] = external[i]->momentum();
  }
  rescaleMomenta(momenta, data);
  double dummy(0.);
  if( external[0]->id() > 0 ) {
    SpinorVector spIn;
    SpinorWaveFunction(spIn, external[0], incoming, false);
    SpinorBarVector spbOut;
    SpinorBarWaveFunction(spbOut, external[2], outgoing, true);

    SpinorWaveFunction spr(rescaledMomenta()[0], data[0], incoming);
    VectorWaveFunction vr(rescaledMomenta()[1], data[1], incoming);
    SpinorBarWaveFunction sbr(rescaledMomenta()[2], data[2], outgoing);

    for( unsigned int ihel = 0; ihel < 2; ++ihel ) {  
      spr.reset(ihel);
      spIn[ihel] = spr;
      vr.reset(2*ihel);
      vecIn[ihel] = vr;
      sbr.reset(ihel);
      spbOut[ihel] = sbr;
    }
    scaOut = ScalarWaveFunction(rescaledMomenta()[3], data[3], outgoing);
    ProductionMatrixElement prodME = fv2fbsHeME(spIn, vecIn, spbOut, 
						scaOut, dummy);
    HardVertexPtr hardvertex = new_ptr(HardVertex());
    hardvertex->ME(prodME);
    for(ParticleVector::size_type i = 0; i < 4; ++i)
      dynamic_ptr_cast<SpinfoPtr>(external[i]->spinInfo())->
	setProductionVertex(hardvertex);
  }
  else {
    SpinorBarVector spbIn;
    SpinorBarWaveFunction(spbIn, external[0], incoming, false);
    SpinorVector spOut;
    SpinorWaveFunction(spOut, external[2], outgoing, true);

    SpinorBarWaveFunction sbr(rescaledMomenta()[0], data[0], incoming);
    VectorWaveFunction vr(rescaledMomenta()[1], data[1], incoming);
    SpinorWaveFunction spr(rescaledMomenta()[2], data[2], outgoing);
    for( unsigned int ihel = 0; ihel < 2; ++ihel ) {  
      sbr.reset(ihel);
      spbIn[ihel] = sbr;
      vr.reset(2*ihel);
      vecIn[ihel] = vr;
      spr.reset(ihel);
      spOut[ihel] = spr;
    }
    scaOut = ScalarWaveFunction(rescaledMomenta()[3], data[3], outgoing);
    ProductionMatrixElement prodME = fbv2fsHeME(spbIn, vecIn, spOut, 
						scaOut, dummy);
    HardVertexPtr hardvertex = new_ptr(HardVertex());
    hardvertex->ME(prodME);
    for(ParticleVector::size_type i = 0; i < 4; ++i)
      dynamic_ptr_cast<SpinfoPtr>(external[i]->spinInfo())->
	setProductionVertex(hardvertex);
  }
  
#ifndef NDEBUG
  if( debugME() ) debug(dummy);
#endif

}

void MEfv2fs::debug(double me2) const {
  if( !generator()->logfile().is_open() ) return;
  long id1 = abs(mePartonData()[0]->id());
  long id4 = abs(mePartonData()[3]->id());
  if( (id1 != 1 && id1 != 2) || mePartonData()[1]->id() != 21 ||
      mePartonData()[2]->id() != 1000021 ||
      (id4 != 1000001 && id4 != 1000002 &&  id4 != 2000001 && 
       id4 != 2000002) ) return;
  tcSMPtr sm = generator()->standardModel();
  double gs4 = sqr( 4.*Constants::pi*sm->alphaS(scale()) );
  int Nc = sm->Nc();
  Energy2 m3s = sqr(mePartonData()[2]->mass());
  Energy2 m4s = sqr(mePartonData()[3]->mass());
  //formula has vf->fs so swap t and u
  Energy2 s(sHat()), t3(uHat() - m3s), u4(tHat() - m4s);
  double analytic = -gs4*( u4 + 2.*(m4s - m3s)*(1. + m3s/t3 + m4s/u4) )*
    ( sqr(u4) + sqr(s) - sqr(t3)/sqr(Nc) )/s/t3/u4/4.;
  double diff = abs( analytic - me2);
  if( diff > 1e-4 ) {
    generator()->log() 
      << mePartonData()[0]->PDGName() << ","
      << mePartonData()[1]->PDGName() << "->"
      << mePartonData()[2]->PDGName() << ","
      << mePartonData()[3]->PDGName() << "   difference: " 
      << setprecision(10) << diff << "  ratio: " << analytic/me2  << '\n';
  }
    
}
