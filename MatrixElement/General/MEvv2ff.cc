// -*- C++ -*-
//
// MEvv2ff.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEvv2ff class.
//

#include "MEvv2ff.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "HardVertex.h"
#include <numeric>

using namespace Herwig;
using ThePEG::Helicity::TensorWaveFunction;
using ThePEG::Helicity::incoming;
using ThePEG::Helicity::outgoing;

void MEvv2ff::doinit() throw(InitException) {
  GeneralHardME::doinit();
  size_t ndiags = numberOfDiags();
  theFerm.resize(ndiags);  theVec.resize(ndiags);
  theTen.resize(ndiags);
  for( size_t i = 0; i < ndiags; ++i ) {
    HPDiagram dg = getProcessInfo()[i];
    if( dg.channelType == HPDiagram::tChannel ) {
      AbstractFFVVertexPtr ffv1 = 
	dynamic_ptr_cast<AbstractFFVVertexPtr>(dg.vertices.first);
      AbstractFFVVertexPtr ffv2 = 
	dynamic_ptr_cast<AbstractFFVVertexPtr>(dg.vertices.second);
      theFerm[i] = make_pair(ffv1, ffv2);
    }
    else if( dg.channelType == HPDiagram::sChannel ) {
      if( dg.intermediate->iSpin() == PDT::Spin1) {
	AbstractVVVVertexPtr vvv = 
	  dynamic_ptr_cast<AbstractVVVVertexPtr>(dg.vertices.first);
	AbstractFFVVertexPtr ffv = 
	  dynamic_ptr_cast<AbstractFFVVertexPtr>(dg.vertices.second);
	theVec[i] = make_pair(vvv,ffv);
      }
      else if(dg.intermediate->iSpin() == PDT::Spin2) {
	AbstractVVTVertexPtr vvt = 
	  dynamic_ptr_cast<AbstractVVTVertexPtr>(dg.vertices.first);
	AbstractFFTVertexPtr fft = 
	  dynamic_ptr_cast<AbstractFFTVertexPtr>(dg.vertices.second);
	theTen[i] = make_pair(vvt,fft);
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
  vv2ffME(v1, v2, sbar, sp, full_me);

#ifndef NDEBUG 
  if( debugME() ) debug(full_me);
#endif

  return full_me;
}

ProductionMatrixElement 
MEvv2ff::vv2ffME(const VBVector & v1, const VBVector & v2,
		 const SpinorBarVector & sbar,const SpinorVector & sp, 
		 double & me2) const {
  const HPCount ndiags(numberOfDiags());
  const size_t ncf(numberOfFlows());
  const vector<vector<double> > cfactors = getColourFactors();
  const Energy mass = sp[0].mass();
  //vectors to store results
  vector<double> me(ndiags, 0.);
  vector<Complex> diag(ndiags, Complex(0.));
  //offshell wavefunctions
  SpinorWaveFunction interF; VectorWaveFunction interV;
  TensorWaveFunction interT;
  const Energy2 q2 = scale();
  ProductionMatrixElement pme(PDT::Spin1, PDT::Spin1,
			      PDT::Spin1Half, PDT::Spin1Half);
  //sum over vector helicities
  for(unsigned int iv1 = 0; iv1 < 2; ++iv1) {
    for(unsigned int iv2 = 0; iv2 < 2; ++iv2) {
      //sum over fermion helicities
      for(unsigned int of1 = 0; of1 < 2; ++of1) {
	for(unsigned int of2 = 0; of2 < 2; ++of2) {
	  //for each helicity calculate diagram 
	  vector<Complex> flows = vector<Complex>(ncf, Complex(0,0.));
	  for(HPCount ix = 0; ix < ndiags; ++ix) {
	    HPDiagram current = getProcessInfo()[ix];
	    PDPtr offshell = current.intermediate;
	    if(current.channelType == HPDiagram::tChannel && 
	       offshell->iSpin() == PDT::Spin1Half) {
	      if(current.ordered.second) {
		interF = theFerm[ix].second->evaluate(q2, 3, offshell, sp[of2], 
						      v2[iv2], mass);
		diag[ix] = theFerm[ix].first->evaluate(q2, interF, sbar[of1], 
						       v1[iv1]);
	      }
	      else {
		interF = theFerm[ix].second->evaluate(q2, 3, offshell, sp[of2], 
						      v1[iv1], mass);
		diag[ix] = theFerm[ix].first->evaluate(q2, interF, sbar[of1], 
						       v2[iv2]);
	      }
	    }
	    else if(current.channelType == HPDiagram::sChannel) {
	      if(offshell->iSpin() == PDT::Spin1) {
		interV = theVec[ix].first->evaluate(q2, 1, offshell, v1[iv1], 
						    v2[iv2]);
		diag[ix] = theVec[ix].second->evaluate(q2, sp[of2], sbar[of1], 
						       interV);
	      }
	      else if(offshell->iSpin() == PDT::Spin2) {
		interT = theTen[ix].first->evaluate(q2, 1, offshell, v1[iv1], 
						    v2[iv2]);
		diag[ix] = theTen[ix].second->evaluate(q2, sp[of2], sbar[of1], 
						       interT);
	      }
	    }
	    else diag[ix] = 0.;
	    me[ix] += norm(diag[ix]);
	    for(size_t iy = 0; iy < current.colourFlow.size(); ++iy) {
	      flows[current.colourFlow[iy].first - 1] += 
		current.colourFlow[iy].second * diag[ix];
	    }
	  }//end of diag loop
	  //set appropriate element of ProductionMatrixElement
	  pme(iv1, iv2, of1, of2) = 
	    std::accumulate(flows.begin(), flows.end(), Complex(0.0, 0.0));
	  
	  for(size_t ii = 0; ii < ncf; ++ii)
	    for(size_t ij = 0; ij < ncf; ++ij)
	      me2 += cfactors[ii][ij]*(flows[ii]*conj(flows[ij])).real();
	
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
  me2 *= identFact/256;
  return pme;
}

Selector<const ColourLines *>
MEvv2ff::colourGeometries(tcDiagPtr diag) const {
  static vector<ColourLines> cf(18);
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
  //88->00
  cf[17] = ColourLines("1 -2,2 -1");

  HPDiagram current = getProcessInfo()[abs(diag->id()) - 1];
  vector<ColourLines>::size_type offset(0);
  if(current.channelType == HPDiagram::tChannel && current.ordered.first)
    offset = 5;
  else if(current.channelType == HPDiagram::tChannel && !current.ordered.first)
    offset = 9;
  else 
    offset = 13;  
  Selector<const ColourLines *> sel;
  if(current.intermediate->iColour() == PDT::Colour0) {
    if(getParticleData(current.outgoing.first)->iColour() == PDT::Colour0)
      sel.insert(1., &cf[17]);
    else
      sel.insert(1., &cf[4]);
  }
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

void MEvv2ff::persistentOutput(PersistentOStream & os) const {
  os << theFerm << theVec << theTen;
}

void MEvv2ff::persistentInput(PersistentIStream & is, int) {
  is >> theFerm >> theVec >> theTen;
}

ClassDescription<MEvv2ff> MEvv2ff::initMEvv2ff;
// Definition of the static class description member.

void MEvv2ff::Init() {

  static ClassDocumentation<MEvv2ff> documentation
    ("The MEvv2ff class handles the ME calculation for the general "
     "spin configuration vector-vector to fermion-antifermion\n.");

}

void MEvv2ff::constructVertex(tSubProPtr sub) {
  ParticleVector ext(4);
  ext[0] = sub->incoming().first;
  ext[1] = sub->incoming().second;
  ext[2] = sub->outgoing()[0];
  ext[3] = sub->outgoing()[1];

  if( ext[2]->id() < ext[3]->id() ) swap(ext[2], ext[3]);
  
  VBVector v1, v2;
  VectorWaveFunction(v1, ext[0], incoming, false, true, true);
  VectorWaveFunction(v2, ext[1], incoming, false, true, true);
  //function to calculate me2 expects massless incoming vectors
  // and this constructor sets the '1' polarisation at element [2] 
  //in the vector
  v1[1] = v1[2]; 
  v2[1] = v2[2];
  SpinorBarVector sbar;
  SpinorBarWaveFunction(sbar, ext[2], outgoing, true, true);
  SpinorVector sp;
  SpinorWaveFunction(sp, ext[3], outgoing, true, true);

  //Need to use rescale momenta to calculate matrix element
  cPDVector data(4);
  vector<Lorentz5Momentum> momenta(4);
  for( size_t i = 0; i < 4; ++i ) {
    data[i] = ext[i]->dataPtr();
    momenta[i] = ext[i]->momentum();
  }
  rescaleMomenta(momenta, data);

  VectorWaveFunction v1r(rescaledMomenta()[0], data[0], incoming),
    v2r(rescaledMomenta()[1], data[1], incoming);
  SpinorBarWaveFunction sbr(rescaledMomenta()[2], data[2], outgoing);
  SpinorWaveFunction spr(rescaledMomenta()[3], data[3], outgoing);
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
  ProductionMatrixElement pme = vv2ffME(v1, v2, sbar, sp, dummy);

#ifndef NDEBUG
  if( debugME() ) debug(dummy);
#endif

  HardVertexPtr hv = new_ptr(HardVertex());
  hv->ME(pme);
  for(unsigned int i = 0; i < 4; ++i) 
    dynamic_ptr_cast<SpinfoPtr>(ext[i]->spinInfo())->setProductionVertex(hv);
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
