// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEff2ss class.
//

#include "MEff2ss.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace Herwig;
using ThePEG::Helicity::SpinorWaveFunction;
using ThePEG::Helicity::SpinorBarWaveFunction;
using ThePEG::Helicity::VectorWaveFunction;
using ThePEG::Helicity::ScalarWaveFunction;
using ThePEG::Helicity::TensorWaveFunction;
using ThePEG::Helicity::incoming;
using ThePEG::Helicity::outgoing;

double MEff2ss::me2() const {
  //first setup  wavefunctions for external particles
  SpinorWaveFunction sp(meMomenta()[0], mePartonData()[0], incoming);
  SpinorBarWaveFunction spbar(meMomenta()[1], mePartonData()[1], incoming);
  ScalarWaveFunction sca1(meMomenta()[2], mePartonData()[2],
			  Complex(1.), outgoing);
  ScalarWaveFunction sca2(meMomenta()[3], mePartonData()[3],
			  Complex(1.), outgoing);
  //Define factors
  const Energy2 m2(scale());
  const vector<vector<double> > cfactors = getColourFactors();
  const HPCount ndiags = numberOfDiags();
  const size_t ncf = numberOfFlows();
  vector<double> me(ndiags, 0.);
  vector<Complex> diag(ndiags, Complex(0.)), flows(ncf, Complex(0.));
  double full_me(0.);
  ScalarWaveFunction interS; VectorWaveFunction interV; 
  SpinorBarWaveFunction interFB; TensorWaveFunction interT;
  for(unsigned int ihel1 = 0; ihel1 < 2; ++ihel1) {
    sp.reset(ihel1);
    for(unsigned int ihel2 = 0; ihel2 < 2; ++ihel2) {
      spbar.reset(ihel2);
      flows = vector<Complex>(ncf, Complex(0.));
      for(HPCount ix = 0; ix < ndiags; ++ix) {
	HPDiagram current = getProcessInfo()[ix];
	tcPDPtr internal(current.intermediate);	
	if(current.channelType == HPDiagram::tChannel &&
	   internal->iSpin() == PDT::Spin1Half) {
	  if(current.ordered.second) {
	    interFB = theFerm[ix].second->evaluate(m2, 3, internal, spbar, sca2);
	    diag[ix] = theFerm[ix].first->evaluate(m2, sp, interFB, sca1);
	  }
	  else {
	    interFB = theFerm[ix].second->evaluate(m2, 3, internal, spbar, sca1);
	    diag[ix] = theFerm[ix].first->evaluate(m2, sp, interFB, sca2);
	  }
	}
	else if(current.channelType == HPDiagram::sChannel) {
	  if(internal->iSpin() == PDT::Spin1) {
	    interV = theVec[ix].first->evaluate(m2, 1, internal, sp, spbar);
	    diag[ix] = theVec[ix].second->evaluate(m2, interV, sca2, sca1);
	  }
	  else if(internal->iSpin() == PDT::Spin2) {
	    interT = theTen[ix].first->evaluate(m2, 1, internal, sp, spbar);
	    diag[ix] = theTen[ix].second ->evaluate(m2, sca2, sca1, interT);
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
      for(unsigned int ii = 0; ii < ncf; ++ii) 
	for(unsigned int ij = 0; ij < ncf; ++ij)
	  full_me += cfactors[ii][ij]*(flows[ii]*conj(flows[ij])).real();
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
  full_me *= 0.25*identFact*colourAvg;
 
#ifndef NDEBUG
  if( debugME() ) debug(full_me);
#endif
  
  return full_me;
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
  Energy4 spt2 = uHat()*tHat() - meMomenta()[2].m2()*meMomenta()[3].m2();
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
  if( diff > 1e-8 ) {
    generator()->log() 
      << mePartonData()[0]->PDGName() << "," 	
      << mePartonData()[1]->PDGName() << "->"
      << mePartonData()[2]->PDGName() << ","
      << mePartonData()[3]->PDGName() << "   difference: " 
      << setprecision(10) << diff << '\n';
  }
    
}
