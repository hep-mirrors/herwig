// -*- C++ -*-
//
// MEff2ss.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEff2ss class.
//

#include "MEff2ss.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "HardVertex.h"
#include <numeric>

using namespace Herwig;
using ThePEG::Helicity::VectorWaveFunction;
using ThePEG::Helicity::TensorWaveFunction;
using ThePEG::Helicity::incoming;
using ThePEG::Helicity::outgoing;

void MEff2ss::doinit() throw(InitException) {
  GeneralHardME::doinit();
  HPCount ndiags(numberOfDiags());
  theFerm.resize(ndiags);
  theVec.resize(ndiags);
  theTen.resize(ndiags);
  for(HPCount i = 0; i < ndiags; ++i) {
    HPDiagram current = getProcessInfo()[i];
    if(current.channelType == HPDiagram::tChannel) {
      if(current.intermediate->iSpin() == PDT::Spin1Half)
	theFerm[i] = 
	  make_pair(dynamic_ptr_cast<AbstractFFSVertexPtr>(current.vertices.first), 
		    dynamic_ptr_cast<AbstractFFSVertexPtr>(current.vertices.second));
    }
    else if(current.channelType == HPDiagram::sChannel) {
      if(current.intermediate->iSpin() == PDT::Spin1)
	theVec[i] = 
	  make_pair(dynamic_ptr_cast<AbstractFFVVertexPtr>(current.vertices.first), 
		    dynamic_ptr_cast<AbstractVSSVertexPtr>(current.vertices.second));

      if(current.intermediate->iSpin() == PDT::Spin2)
	theTen[i] = 
	  make_pair(dynamic_ptr_cast<AbstractFFTVertexPtr>(current.vertices.first), 
		    dynamic_ptr_cast<AbstractSSTVertexPtr>(current.vertices.second));
    }
    else 
      throw InitException() << "MEFF2ss:doinit() - Cannot find correct "
			    << "channel from diagram. Vertex not cast! "
			    << Exception::runerror;
  }
}

double MEff2ss::me2() const {
  //first setup  wavefunctions for external particles
  SpinorVector sp(2);
  SpinorBarVector sbar(2);
  for( unsigned int i = 0; i < 2; ++i ) {
    sp[i] = SpinorWaveFunction(rescaledMomenta()[0], mePartonData()[0], i,
			       incoming);
    sbar[i] = SpinorBarWaveFunction(rescaledMomenta()[1], mePartonData()[1], i,
				    incoming);
  }
  ScalarWaveFunction sca1(rescaledMomenta()[2], mePartonData()[2],
			  Complex(1.), outgoing);
  ScalarWaveFunction sca2(rescaledMomenta()[3], mePartonData()[3],
			  Complex(1.), outgoing);
  double full_me(0.);
  ff2ssME(sp, sbar, sca1, sca2, full_me);  
 
#ifndef NDEBUG
  if( debugME() ) debug(full_me);
#endif
  
  return full_me;
}

ProductionMatrixElement 
MEff2ss::ff2ssME(const SpinorVector & sp, const SpinorBarVector & sbar, 
		 const ScalarWaveFunction & sca1, 
		 const ScalarWaveFunction & sca2, double & me2) const {
  //Define factors
  const Energy2 q2(scale());
  const vector<vector<double> > cfactors = getColourFactors();
  const HPCount ndiags = numberOfDiags();
  const size_t ncf = numberOfFlows();
  vector<double> me(ndiags, 0.);
  vector<Complex> diag(ndiags, Complex(0.)), flows(ncf, Complex(0.));
  ScalarWaveFunction interS; VectorWaveFunction interV; 
  SpinorBarWaveFunction interFB; TensorWaveFunction interT;
  ProductionMatrixElement pme(PDT::Spin1Half, PDT::Spin1Half, 
			      PDT::Spin0, PDT::Spin0);
  for(unsigned int if1 = 0; if1 < 2; ++if1) {
    for(unsigned int if2 = 0; if2 < 2; ++if2) {
      flows = vector<Complex>(ncf, Complex(0.));
      for(HPCount ix = 0; ix < ndiags; ++ix) {
	HPDiagram current = getProcessInfo()[ix];
	tcPDPtr internal(current.intermediate);	
	if(current.channelType == HPDiagram::tChannel &&
	   internal->iSpin() == PDT::Spin1Half) {
	  if(current.ordered.second) {
	    interFB = theFerm[ix].second->evaluate(q2, 3, internal, sbar[if2], 
						   sca2);
	    diag[ix] = theFerm[ix].first->evaluate(q2, sp[if1], interFB, sca1);
	  }
	  else {
	    interFB = theFerm[ix].second->evaluate(q2, 3, internal, sbar[if2], 
						   sca1);
	    diag[ix] = theFerm[ix].first->evaluate(q2, sp[if1], interFB, sca2);
	  }
	}
	else if(current.channelType == HPDiagram::sChannel) {
	  if(internal->iSpin() == PDT::Spin1) {
	    interV = theVec[ix].first->evaluate(q2, 1, internal, sp[if1], 
						sbar[if2]);
	    diag[ix] = theVec[ix].second->evaluate(q2, interV, sca2, sca1);
	  }
	  else if(internal->iSpin() == PDT::Spin2) {
	    interT = theTen[ix].first->evaluate(q2, 1, internal, sp[if1], 
						sbar[if2]);
	    diag[ix] = theTen[ix].second ->evaluate(q2, sca2, sca1, interT);
	  }
	  else 
	    diag[ix] = 0.;
	}
	me[ix] += norm(diag[ix]);
	//colourflows
	for(unsigned int iy = 0; iy < current.colourFlow.size(); ++iy)
	  flows[current.colourFlow[iy].first - 1] += 
	    current.colourFlow[iy].second * diag[ix];

      }//end of diag loop
      //set the appropriate element of the ProductionMatrixElement
      pme(if1, if2, 0, 0) = 
	std::accumulate(flows.begin(), flows.end(), Complex(0.0,0.0));

      for(unsigned int ii = 0; ii < ncf; ++ii) 
	for(unsigned int ij = 0; ij < ncf; ++ij)
	  me2 += cfactors[ii][ij]*(flows[ii]*conj(flows[ij])).real();
    }
  }
  const double identFact = mePartonData()[2]->id() == mePartonData()[3]->id() 
    ? 0.5 : 1;
  int cola = mePartonData()[0]->iColour();
  const double colourAvg = ( abs(cola) == 3 ) ? 1./9. : 1.;
  DVector save(ndiags);
  for(DVector::size_type ix = 0; ix < ndiags; ++ix)
    save[ix] = 0.25*identFact*colourAvg*me[ix];
  meInfo(save);
  me2 *= 0.25*identFact*colourAvg;
  return pme;
}

Selector<const ColourLines *>
MEff2ss::colourGeometries(tcDiagPtr diag) const {
  static vector<ColourLines> cf(16);
  //33b->33b
  cf[0] = ColourLines("1 2 -3, 4 -2 -5");
  cf[1] = ColourLines("1 3 4, -2 -3 -5");
  cf[2] = ColourLines("1 4, -3 -5");	 
  cf[3] = ColourLines("1 -2, 4 -5");
  //33->33
  cf[4] = ColourLines("1 2 5, 3 -2 4");  
  cf[5] = ColourLines("1 2 4, 3 -2 5");
  cf[6] = ColourLines("1 4, 3 5");	
  cf[7] = ColourLines("1 5, 3 4");
  //3b3b->3b3b
  cf[8] = ColourLines("-1 -2 -5, -3 2 -4");
  cf[9] = ColourLines("-1 -2 -4, -3 2 -5");
  cf[10] = ColourLines("-1 -4, -3 -5");
  cf[11] = ColourLines("-1 -5, -3 -4");  
  //33b->11
  cf[12] = ColourLines("1 2 -3");
  cf[13] = ColourLines("1 -2");
  //11->11
  cf[14] = ColourLines("");
  //11->33b
  cf[15] = ColourLines("4 -5");
  HPDiagram current = getProcessInfo()[abs(diag->id()) - 1];
  vector<ColourLines>::size_type cl(0);
  PDT::Colour inac(mePartonData()[0]->iColour());
  PDT::Colour inbc(mePartonData()[1]->iColour());
  PDT::Colour outac(mePartonData()[2]->iColour());
  if(inac == PDT::Colour0 && inbc == PDT::Colour0) {
    cl = outac == PDT::Colour0 ?  14 : 15;
  }
  else if(inac == PDT::Colour3 && inbc == PDT::Colour3) {
    if(current.intermediate->iColour() == PDT::Colour8)
      cl = current.ordered.second ? 4 : 5;
    else
      cl = current.ordered.second ? 6 : 7;
  }
  else if(inac == PDT::Colour3bar && inbc == PDT::Colour3bar) {
    if(current.intermediate->iColour() == PDT::Colour8)
      cl = current.ordered.second ? 8 : 9;
    else
      cl = current.ordered.second ? 10 : 11 ;
  }
  else {
    if(outac == PDT::Colour3) {
      if(current.intermediate->iColour() == PDT::Colour8)
	cl = current.channelType == HPDiagram::tChannel ? 0 : 1;
      else
	cl = current.channelType == HPDiagram::tChannel ? 2 : 3;
    }
    else
      cl = current.channelType == HPDiagram::tChannel ? 12 : 13;
  }
  Selector<const ColourLines *> sel;
  sel.insert(1., &cf[cl]);
  return sel;
}

void MEff2ss::persistentOutput(PersistentOStream & os) const {
  os << theFerm << theVec << theTen;
}

void MEff2ss::persistentInput(PersistentIStream & is, int) {
  is >> theFerm >> theVec >> theTen;
}

ClassDescription<MEff2ss> MEff2ss::initMEff2ss;
// Definition of the static class description member.

void MEff2ss::Init() {

  static ClassDocumentation<MEff2ss> documentation
    ("MEff2ss implements the ME calculation of the fermion-antifermion "
     "to scalar-scalar hard process.");

}

void MEff2ss::constructVertex(tSubProPtr sub) {
  //get particles
  ParticleVector ext(4);
  ext[0] = sub->incoming().first;
  ext[1] = sub->incoming().second;
  ext[2] = sub->outgoing()[0];
  ext[3] = sub->outgoing()[1];

  if( ext[0]->id() != mePartonData()[0]->id() ) swap(ext[0], ext[1]);
  if( ext[2]->id() != mePartonData()[2]->id() ) swap(ext[2], ext[3]);

  //First calculate wave functions with off-shell momenta
  //to calculate correct spin information
  SpinorVector sp;
  SpinorWaveFunction(sp, ext[0], incoming, false, true);
  SpinorBarVector sbar;
  SpinorBarWaveFunction(sbar, ext[1], incoming, false, true);
  ScalarWaveFunction sca1(ext[2], outgoing, true, true);
  ScalarWaveFunction sca2(ext[3], outgoing, true, true);
  //Need to use rescale momenta to calculate matrix element
  cPDVector data(4);
  vector<Lorentz5Momentum> momenta(4);
  for( size_t i = 0; i < 4; ++i ) {
    data[i] = ext[i]->dataPtr();
    momenta[i] = ext[i]->momentum();
  }
  rescaleMomenta(momenta, data);
  SpinorWaveFunction spr(rescaledMomenta()[0], data[0], incoming);
  SpinorBarWaveFunction sbr(rescaledMomenta()[1], data[1],incoming);
  for( unsigned int ihel = 0; ihel < 2; ++ihel ) {  
    spr.reset(ihel);
    sp[ihel] = spr;

    sbr.reset(ihel);
    sbar[ihel] = sbr;
  }
  sca1 = ScalarWaveFunction(rescaledMomenta()[2], data[2], outgoing);
  sca2 = ScalarWaveFunction(rescaledMomenta()[3], data[3], outgoing);
  double dummy(0.);
  ProductionMatrixElement pme = ff2ssME(sp, sbar, sca1, sca2, dummy);
#ifndef NDEBUG
  if( debugME() ) debug(dummy);
#endif
  HardVertexPtr hv = new_ptr(HardVertex());
  hv->ME(pme);
  for(unsigned int i = 0; i < 4; ++i )
    dynamic_ptr_cast<SpinfoPtr>(ext[i]->spinInfo())->setProductionVertex(hv);
  
}

void MEff2ss::debug(double me2) const {
  if( !generator()->logfile().is_open() ) return;
  long id1 = mePartonData()[0]->id();
  long id2 = mePartonData()[1]->id();
  long id3 = mePartonData()[2]->id();
  long id4 = mePartonData()[3]->id();
  if( (abs(id1) != 1 && abs(id1) != 2) || (abs(id2) != 1 && abs(id2) != 2) ||
      ( abs(id3) != 1000001 && abs(id3) != 1000002 && 
        abs(id3) != 2000001 && abs(id3) != 2000002 ) || 
      ( abs(id4) != 1000001 && abs(id4) != 1000002  &&
	abs(id4) != 2000001 && abs(id4) != 2000002 ) ) return;
  tcSMPtr sm = generator()->standardModel();
  double gs4 = sqr( 4.*Constants::pi*sm->alphaS(scale()) );
  int Nc = sm->Nc();
  double Cf = (sqr(Nc) - 1)/2./Nc;
  Energy2 s(sHat());
  Energy2 mgos = sqr( getParticleData(ParticleID::SUSY_g)->mass());
  Energy2 m3s = sqr(mePartonData()[2]->mass());
  Energy2 m4s = sqr(mePartonData()[3]->mass());
  Energy4 spt2 = uHat()*tHat() - m3s*m4s;
  Energy2 tgl(tHat() - mgos), ugl(uHat() - mgos);
  unsigned int alpha = abs(id3)/1000000;
  unsigned int beta = abs(id4)/1000000;
  bool iflav = ( abs(id1) == abs(id2) );
  unsigned int oflav = ( abs(id3) - abs(id1) ) % 10;
  
  double analytic(0.);
  if( alpha != beta ) {
    if( ( id1 > 0 && id2 > 0) ||
	( id1 < 0 && id2 < 0) ) { 
      analytic = spt2/sqr(tgl);
      if( iflav ) analytic += spt2/sqr(ugl);
    }
    else {
      analytic = s*mgos/sqr(tgl);
    }
  }
  else {
    if( oflav != 0 ) {
      analytic = 2.*spt2/sqr(s);
    }
    else if( ( id1 > 0 && id2 > 0) ||
	     ( id1 < 0 && id2 < 0) ) {
      analytic = s*mgos/sqr(tgl);
      if( iflav ) {
	analytic += s*mgos/sqr(ugl) - 2.*s*mgos/Nc/tgl/ugl;
      }
      analytic /= ( iflav ? 2. : 1.);
    }
    else {
      analytic = spt2/sqr(tgl);
      if( iflav ) {
	analytic += 2.*spt2/sqr(s) - 2.*spt2/Nc/s/tgl;
      }
    }
  }
  analytic *= gs4*Cf/2./Nc;
  double diff = abs(analytic - me2);
  if( diff > 1e-4 ) {
    generator()->log() 
      << mePartonData()[0]->PDGName() << "," 	
      << mePartonData()[1]->PDGName() << "->"
      << mePartonData()[2]->PDGName() << ","
      << mePartonData()[3]->PDGName() << "   difference: " 
      << setprecision(10) << diff << "  ratio: " << analytic/me2 
      << '\n';
  }
    
}
