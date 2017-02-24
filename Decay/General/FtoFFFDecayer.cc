// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FtoFFFDecayer class.
//

#include "FtoFFFDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Decay/DecayPhaseSpaceMode.h"
#include "Herwig/PDT/ThreeBodyAllOnCalculator.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"
#include <numeric>

using namespace Herwig;

IBPtr FtoFFFDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr FtoFFFDecayer::fullclone() const {
  return new_ptr(*this);
}

void FtoFFFDecayer::persistentOutput(PersistentOStream & os) const {
  os << _sca << _vec << _ten;
}

void FtoFFFDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _sca >> _vec >> _ten;
}

ClassDescription<FtoFFFDecayer> FtoFFFDecayer::initFtoFFFDecayer;
// Definition of the static class description member.

void FtoFFFDecayer::Init() {

  static ClassDocumentation<FtoFFFDecayer> documentation
    ("The FtoFFFDecayer class implements the general decay of a fermion to "
     "three fermions.");

}

void FtoFFFDecayer::doinit() {
  GeneralThreeBodyDecayer::doinit();
  unsigned int ndiags = getProcessInfo().size();
  _sca.resize(ndiags);
  _vec.resize(ndiags);
  _ten.resize(ndiags);
  for(unsigned int ix = 0;ix < ndiags; ++ix) {
    TBDiagram current = getProcessInfo()[ix];
    tcPDPtr offshell = current.intermediate;
    if( offshell->CC() ) offshell = offshell->CC();
    if(offshell->iSpin() == PDT::Spin0) {
      AbstractFFSVertexPtr vert1 = dynamic_ptr_cast<AbstractFFSVertexPtr>
	(current.vertices.first);
      AbstractFFSVertexPtr vert2 = dynamic_ptr_cast<AbstractFFSVertexPtr>
	(current.vertices.second);
      if(!vert1||!vert2) throw Exception() 
	<< "Invalid vertices for a scalar diagram in FtoFFFDecayer::doinit()"
	<< Exception::runerror;
      _sca[ix] = make_pair(vert1, vert2);
    }
    else if(offshell->iSpin() == PDT::Spin1) {
      AbstractFFVVertexPtr vert1 = dynamic_ptr_cast<AbstractFFVVertexPtr>
	(current.vertices.first);
      AbstractFFVVertexPtr vert2 = dynamic_ptr_cast<AbstractFFVVertexPtr>
	(current.vertices.second);
      if(!vert1||!vert2) throw Exception() 
	<< "Invalid vertices for a vector diagram in FtoFFFDecayer::doinit()"
	<< Exception::runerror;
      _vec[ix] = make_pair(vert1, vert2);
    }
    else if(offshell->iSpin() == PDT::Spin2) {
      AbstractFFTVertexPtr vert1 = dynamic_ptr_cast<AbstractFFTVertexPtr>
	(current.vertices.first);
      AbstractFFTVertexPtr vert2 = dynamic_ptr_cast<AbstractFFTVertexPtr>
	(current.vertices.second);
      if(!vert1||!vert2) throw Exception() 
	<< "Invalid vertices for a tensor diagram in FtoFFFDecayer::doinit()"
	<< Exception::runerror;
      _ten[ix] = make_pair(vert1, vert2);
    }
  }
}

double  FtoFFFDecayer::me2(const int ichan, const Particle & inpart,
			   const ParticleVector & decay,
			   MEOption meopt) const {
  // particle or CC of particle
  bool cc = (*getProcessInfo().begin()).incoming != inpart.id();
  // special handling or first/last call
  const vector<vector<double> > cfactors(getColourFactors());
  const vector<vector<double> > nfactors(getLargeNcColourFactors());
  const size_t ncf(numberOfFlows());
  Energy2 scale(sqr(inpart.mass()));
  if(meopt==Initialize) {
    SpinorWaveFunction::
      calculateWaveFunctions(_inwave.first,_rho,const_ptr_cast<tPPtr>(&inpart),
			    Helicity::incoming);
    _inwave.second.resize(2);
    if(_inwave.first[0].wave().Type() == SpinorType::u) {
      for(unsigned int ix = 0; ix < 2; ++ix) {
	_inwave.second[ix] = _inwave.first[ix].bar();
	_inwave.second[ix].conjugate();
      }
    }
    else {
      for(unsigned int ix = 0; ix < 2; ++ix) {
	_inwave.second[ix] = _inwave.first[ix].bar();
	_inwave.first[ix].conjugate();
      }
    }
  }
  // setup spin info when needed
  if(meopt==Terminate) {
    // for the decaying particle
    if(inpart.id()<0) 
      SpinorWaveFunction::constructSpinInfo(_inwave.first,
					    const_ptr_cast<tPPtr>(&inpart),
					    Helicity::incoming,true);
    else
      SpinorBarWaveFunction::constructSpinInfo(_inwave.second,
					       const_ptr_cast<tPPtr>(&inpart),
					       Helicity::incoming,true);
    // outgoing particles
    for(unsigned int ix = 0; ix < 3; ++ix) {
      SpinorWaveFunction::
      constructSpinInfo(_outwave[ix].first,decay[ix],Helicity::outgoing,true);
    }
  }
  // outgoing particles
  for(unsigned int ix = 0; ix < 3; ++ix) {
    SpinorWaveFunction::
      calculateWaveFunctions(_outwave[ix].first,decay[ix],Helicity::outgoing);
    _outwave[ix].second.resize(2);
    if(_outwave[ix].first[0].wave().Type() == SpinorType::u) {
      for(unsigned int iy = 0; iy < 2; ++iy) {
	_outwave[ix].second[iy] = _outwave[ix].first[iy].bar();
	_outwave[ix].first[iy].conjugate();
      }
    }
    else {
      for(unsigned int iy = 0; iy < 2; ++iy) {
	_outwave[ix].second[iy] = _outwave[ix].first[iy].bar();
	_outwave[ix].second[iy].conjugate();
      }
    }
  }
  bool ferm = inpart.id()>0; 
  vector<Complex> flows(ncf, Complex(0.)),largeflows(ncf, Complex(0.)); 
  static const unsigned int out2[3]={1,0,0},out3[3]={2,2,1};
  vector<GeneralDecayMEPtr> mes(ncf,new_ptr(GeneralDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,
								      PDT::Spin1Half,PDT::Spin1Half)));
  vector<GeneralDecayMEPtr> mel(ncf,new_ptr(GeneralDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,
								      PDT::Spin1Half,PDT::Spin1Half)));
  unsigned int ihel[4];
  for(ihel[0] = 0; ihel[0] < 2; ++ihel[0]) {
    for(ihel[1] = 0; ihel[1] < 2; ++ihel[1]) {
      for(ihel[2] = 0; ihel[2] < 2; ++ihel[2]) {
	for(ihel[3] = 0; ihel[3] < 2; ++ihel[3]) {
	  flows = vector<Complex>(ncf, Complex(0.));
	  largeflows = vector<Complex>(ncf, Complex(0.));
	  unsigned int idiag=0;
	  for(vector<TBDiagram>::const_iterator dit=getProcessInfo().begin();
	      dit!=getProcessInfo().end();++dit) {
	    if(ichan>=0&&diagramMap()[ichan]!=idiag) {
	      ++idiag;
	      continue;
	    }
	    // the sign from normal ordering
	    double sign = ferm ? 1. : -1;
	    // outgoing wavefunction and NO sign
	    if     (dit->channelType==TBDiagram::channel23) sign *= -1.;
	    else if(dit->channelType==TBDiagram::channel13) sign *=  1.;
	    else if(dit->channelType==TBDiagram::channel12) sign *= -1.;
	    else throw Exception()
	      << "Unknown diagram type in FtoFFFDecayer::me2()" << Exception::runerror;
	    // wavefunctions
	    SpinorWaveFunction    w0,w3;
	    SpinorBarWaveFunction w1,w2;
	    // incoming wavefunction
	    if(ferm) {
	      w0 = _inwave.first [ihel[0]];
	      w1 = _outwave[dit->channelType].second[ihel[dit->channelType+1]];
	    }
	    else {
	      w0 = _outwave[dit->channelType].first [ihel[dit->channelType+1]];
	      w1 = _inwave.second[ihel[0]];
	    }
	    if(decay[out2[dit->channelType]]->id()<0&&
	       decay[out3[dit->channelType]]->id()>0) {
	      w2 = _outwave[out3[dit->channelType]].second[ihel[out3[dit->channelType]+1]];
	      w3 = _outwave[out2[dit->channelType]].first [ihel[out2[dit->channelType]+1]];
	      sign *= -1.;
	    }
	    else {
	      w2 = _outwave[out2[dit->channelType]].second[ihel[out2[dit->channelType]+1]];
	      w3 = _outwave[out3[dit->channelType]].first [ihel[out3[dit->channelType]+1]];
	    }
	    tcPDPtr offshell = dit->intermediate;
	    if(cc&&offshell->CC()) offshell=offshell->CC();
	    Complex diag(0.);
	    // intermediate scalar
	    if     (offshell->iSpin() == PDT::Spin0) { 
	      ScalarWaveFunction inters = _sca[idiag].first->
		evaluate(scale, widthOption(), offshell, w0, w1);
	      diag = _sca[idiag].second->evaluate(scale,w3,w2,inters);
	    }
	    // intermediate vector
	    else if(offshell->iSpin() == PDT::Spin1) {
	      VectorWaveFunction interv = _vec[idiag].first->
		evaluate(scale, widthOption(), offshell, w0, w1);
	      diag = _vec[idiag].second->evaluate(scale,w3,w2,interv);
	    }
	    // intermediate tensor
	    else if(offshell->iSpin() == PDT::Spin2) {
	      TensorWaveFunction intert = _ten[idiag].first->
		evaluate(scale, widthOption(), offshell, w0, w1);
	      diag = _ten[idiag].second->evaluate(scale,w3,w2,intert);
	    }
	    // unknown
	    else throw Exception()
	      << "Unknown intermediate in FtoFFFDecayer::me2()" 
	      << Exception::runerror;
	    // apply NO sign
	    diag *= sign;
	    // matrix element for the different colour flows
	    if(ichan<0) {
	      for(unsigned iy = 0; iy < dit->colourFlow.size(); ++iy) {
		flows[dit->colourFlow[iy].first - 1] += 
		  dit->colourFlow[iy].second * diag;
	      }
	      for(unsigned iy = 0; iy < dit->largeNcColourFlow.size(); ++iy) {
		largeflows[dit->largeNcColourFlow[iy].first - 1] += 
		  dit->largeNcColourFlow[iy].second * diag;
	      }
	    }
	    else {
	      for(unsigned iy = 0; iy < dit->colourFlow.size(); ++iy) {
		if(dit->colourFlow[iy].first - 1!=colourFlow()) continue;
		flows[dit->colourFlow[iy].first - 1] += 
		  dit->colourFlow[iy].second * diag;
	      }
	      for(unsigned iy = 0; iy < dit->largeNcColourFlow.size(); ++iy) {
		if(dit->colourFlow[iy].first - 1!=colourFlow()) continue;
		largeflows[dit->largeNcColourFlow[iy].first - 1] += 
		  dit->largeNcColourFlow[iy].second * diag;
	      }
	    }
	    ++idiag;
	  }
	  // now add the flows to the me2 with appropriate colour factors
	  for(unsigned int ix = 0; ix < ncf; ++ix) {
	    (*mes[ix])(ihel[0],ihel[1],ihel[2],ihel[3]) =      flows[ix];
	    (*mel[ix])(ihel[0],ihel[1],ihel[2],ihel[3]) = largeflows[ix];
	  }
	}
      }
    }
  }
  double me2(0.);
  if(ichan<0) {
    vector<double> pflows(ncf,0.);
    for(unsigned int ix = 0; ix < ncf; ++ix) {
      for(unsigned int iy = 0; iy < ncf; ++ iy) {
	double con = cfactors[ix][iy]*(mes[ix]->contract(*mes[iy],_rho)).real();
	me2 += con;
	if(ix==iy) {
	  con = nfactors[ix][iy]*(mel[ix]->contract(*mel[iy],_rho)).real();
	  pflows[ix] += con;
	}
      }
    }
    double ptotal(std::accumulate(pflows.begin(),pflows.end(),0.));
    ptotal *=UseRandom::rnd();
    for(unsigned int ix=0;ix<pflows.size();++ix) {
      if(ptotal<=pflows[ix]) {
	colourFlow(ix);
	ME(mes[ix]);
	break;
      }
      ptotal-=pflows[ix];
    }
  }
  else {
    unsigned int iflow = colourFlow();
    me2 = nfactors[iflow][iflow]*(mel[iflow]->contract(*mel[iflow],_rho)).real();
  }
  // return the matrix element squared
  return me2;
}

WidthCalculatorBasePtr FtoFFFDecayer::
threeBodyMEIntegrator(const DecayMode & ) const {
  vector<int> intype;
  vector<Energy> inmass,inwidth;
  vector<double> inpow,inweights;
  constructIntegratorChannels(intype,inmass,inwidth,inpow,inweights);
  return new_ptr(ThreeBodyAllOnCalculator<FtoFFFDecayer>
		 (inweights,intype,inmass,inwidth,inpow,*this,0,
		  outgoing()[0]->mass(),outgoing()[1]->mass(),outgoing()[2]->mass(),
		 relativeError()));
}
