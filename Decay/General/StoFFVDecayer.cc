// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the StoFFVDecayer class.
//

#include "StoFFVDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/PDT/ThreeBodyAllOnCalculator.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"
#include <numeric>

using namespace Herwig;
using namespace ThePEG;
using namespace ThePEG::Helicity;

IBPtr StoFFVDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr StoFFVDecayer::fullclone() const {
  return new_ptr(*this);
}

void StoFFVDecayer::persistentOutput(PersistentOStream & os) const {
  os << _sca << _fer << _vec;
}

void StoFFVDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _sca >> _fer >> _vec;
}

ClassDescription<StoFFVDecayer> StoFFVDecayer::initStoFFVDecayer;
// Definition of the static class description member.

void StoFFVDecayer::Init() {

  static ClassDocumentation<StoFFVDecayer> documentation
    ("The StoFFVDecayer class implements the general decay of a scalar to "
     "a two fermions and a vector.");

}

WidthCalculatorBasePtr StoFFVDecayer::
threeBodyMEIntegrator(const DecayMode & ) const {
  vector<int> intype;
  vector<Energy> inmass,inwidth;
  vector<double> inpow,inweights;
  constructIntegratorChannels(intype,inmass,inwidth,inpow,inweights);
  return new_ptr(ThreeBodyAllOnCalculator<StoFFVDecayer>
		 (inweights,intype,inmass,inwidth,inpow,*this,0,
		  outgoing()[0]->mass(),outgoing()[1]->mass(),
		  outgoing()[2]->mass(),relativeError()));
}

void StoFFVDecayer::doinit() {
  GeneralThreeBodyDecayer::doinit(); 
  unsigned int ndiags = getProcessInfo().size();
  _sca.resize(ndiags);
  _fer.resize(ndiags);
  _vec.resize(ndiags);
  for(unsigned int ix = 0;ix < ndiags; ++ix) {
    TBDiagram current = getProcessInfo()[ix];
    tcPDPtr offshell = current.intermediate;
    if( offshell->CC() ) offshell = offshell->CC();
    if(offshell->iSpin() == PDT::Spin0) {
      AbstractVSSVertexPtr vert1 = dynamic_ptr_cast<AbstractVSSVertexPtr>
	(current.vertices.first);
      AbstractFFSVertexPtr vert2 = dynamic_ptr_cast<AbstractFFSVertexPtr>
	(current.vertices.second);
      if(!vert1||!vert2) throw Exception() 
	<< "Invalid vertices for a scalar diagram in StoFFVDecayer::doinit()"
	<< Exception::runerror;
      _sca[ix] = make_pair(vert1, vert2);
    }
    else if(offshell->iSpin() == PDT::Spin1Half) {
      AbstractFFSVertexPtr vert1 = dynamic_ptr_cast<AbstractFFSVertexPtr>
	(current.vertices.first);
      AbstractFFVVertexPtr vert2 = dynamic_ptr_cast<AbstractFFVVertexPtr>
	(current.vertices.second);
      if(!vert1||!vert2) throw Exception() 
	<< "Invalid vertices for a fermion diagram in StoFFVDecayer::doinit()"
	<< Exception::runerror;
      _fer[ix] = make_pair(vert1, vert2);
    }
    else if(offshell->iSpin() == PDT::Spin1) {
      AbstractVVSVertexPtr vert1 = dynamic_ptr_cast<AbstractVVSVertexPtr>
	(current.vertices.first);
      AbstractFFVVertexPtr vert2 = dynamic_ptr_cast<AbstractFFVVertexPtr>
	(current.vertices.second);
      if(!vert1||!vert2) throw Exception() 
	<< "Invalid vertices for a vector diagram in StoFFVDecayer::doinit()"
	<< Exception::runerror;
      _vec[ix] = make_pair(vert1, vert2);
    }
  }
}

double StoFFVDecayer::me2(const int ichan, const Particle & inpart,
			  const ParticleVector & decay, MEOption meopt) const {
  // particle or CC of particle
  bool cc = (*getProcessInfo().begin()).incoming != inpart.id();
  // special handling or first/last call
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(_rho,const_ptr_cast<tPPtr>(&inpart),
			     Helicity::incoming);
    _swave = ScalarWaveFunction(inpart.momentum(),inpart.dataPtr(),
				Helicity::incoming);
  }
  if(meopt==Terminate) {
    ScalarWaveFunction::
      constructSpinInfo(const_ptr_cast<tPPtr>(&inpart),
			Helicity::incoming,true);
    for(unsigned int ix=0;ix<decay.size();++ix) {
      if(decay[ix]->dataPtr()->iSpin()==PDT::Spin1) {
	VectorWaveFunction::constructSpinInfo(_outVector,decay[ix],
					      Helicity::outgoing,true,false);
      }
      else {
	SpinorWaveFunction::
	  constructSpinInfo(_outspin[ix].first,decay[ix],Helicity::outgoing,true);
      }
    }
  }
  unsigned int ivec(0);
  bool massless(false);
  for(unsigned int ix = 0; ix < decay.size();++ix) {
    if(decay[ix]->dataPtr()->iSpin() == PDT::Spin1) {
      ivec = ix;
      massless = decay[ivec]->mass()==ZERO;
      VectorWaveFunction::
	calculateWaveFunctions(_outVector, decay[ix], Helicity::outgoing,massless);
    }
    else {
      SpinorWaveFunction::
	calculateWaveFunctions(_outspin[ix].first,decay[ix],Helicity::outgoing);
      _outspin[ix].second.resize(2);
      // Need a ubar and a v spinor
      if(_outspin[ix].first[0].wave().Type() == u_spinortype) {
	for(unsigned int iy = 0; iy < 2; ++iy) {
	  _outspin[ix].second[iy] = _outspin[ix].first[iy].bar();
	  _outspin[ix].first[iy].conjugate();
	}
      }
      else {
	for(unsigned int iy = 0; iy < 2; ++iy) {
	  _outspin[ix].second[iy] = _outspin[ix].first[iy].bar();
	  _outspin[ix].second[iy].conjugate();
	}
      }
    }
  }
  const vector<vector<double> > cfactors(getColourFactors());
  const vector<vector<double> > nfactors(getLargeNcColourFactors());
  Energy2 scale(sqr(inpart.mass()));
  const size_t ncf(numberOfFlows());
  vector<Complex> flows(ncf, Complex(0.)), largeflows(ncf, Complex(0.));
  // setup the DecayMatrixElement
  vector<GeneralDecayMEPtr> 
    mes(ncf,new_ptr(GeneralDecayMatrixElement(PDT::Spin0,
					      ivec == 0 ? PDT::Spin1 : PDT::Spin1Half,
					      ivec == 1 ? PDT::Spin1 : PDT::Spin1Half,
					      ivec == 2 ? PDT::Spin1 : PDT::Spin1Half)));
  vector<GeneralDecayMEPtr> 
    mel(ncf,new_ptr(GeneralDecayMatrixElement(PDT::Spin0,
					      ivec == 0 ? PDT::Spin1 : PDT::Spin1Half,
					      ivec == 1 ? PDT::Spin1 : PDT::Spin1Half,
					      ivec == 2 ? PDT::Spin1 : PDT::Spin1Half)));
  //the channel possiblities
  static const unsigned int out2[3] = {1,0,0}, out3[3] = {2,2,1};
  for(unsigned int s1 = 0; s1 < 2; ++s1) {
    for(unsigned int s2 = 0; s2 < 2; ++s2) {
      for(unsigned int v1 = 0; v1 < 3; ++v1) {
	if(massless&&v1==1) ++v1;
	flows = vector<Complex>(ncf, Complex(0.));
	largeflows = vector<Complex>(ncf, Complex(0.));
	unsigned int idiag(0);
	for(vector<TBDiagram>::const_iterator dit=getProcessInfo().begin();
	    dit!=getProcessInfo().end();++dit) {
	  // channels if selecting
	  if( ichan >= 0 && diagramMap()[ichan] != idiag ) {
	    ++idiag;
	    continue;
	  }
	  tcPDPtr offshell = dit->intermediate;
	  if(cc&&offshell->CC()) offshell=offshell->CC();
	  Complex diag;
	  unsigned int o2(out2[dit->channelType]), o3(out3[dit->channelType]);
	  double sign = (o3 < o2) ? 1. : -1.;
	  // intermediate scalar
	  if(offshell->iSpin() == PDT::Spin0) {
	    ScalarWaveFunction inters = _sca[idiag].first->
	      evaluate(scale, widthOption(), offshell, _outVector[v1], _swave);
	    unsigned int h1(s1),h2(s2);
	    if(o2 > o3) swap(h1, h2);
	    if(decay[o2]->id() < 0 &&  decay[o3]->id() > 0) {
	      diag = -sign*_sca[idiag].second->
		evaluate(scale,_outspin[o2].first[h1],
			 _outspin[o3].second[h2],inters);
	    }
	    else {
	      diag = sign*_sca[idiag].second->
		evaluate(scale, _outspin[o3].first [h2],
			 _outspin[o2].second[h1],inters);
	    }
	  }
	  // intermediate fermion
	  else if(offshell->iSpin() == PDT::Spin1Half) {
	    int iferm = (decay[o2]->dataPtr()->iSpin() == PDT::Spin1Half) 
	      ? o2 : o3;
	    unsigned int h1(s1),h2(s2);
	    if(dit->channelType > iferm) swap(h1, h2);
	    sign = iferm < dit->channelType ? 1. : -1.;
	    if((decay[dit->channelType]->id() < 0 && decay[iferm]->id() > 0 ) ||
	       (decay[dit->channelType]->id()*offshell->id()>0)) {
	      SpinorWaveFunction inters = _fer[idiag].first->
		evaluate(scale,widthOption(),offshell,
			 _outspin[dit->channelType].first[h1], _swave);
	      diag = -sign*_fer[idiag].second->
		evaluate(scale,inters,_outspin[iferm].second[h2], _outVector[v1]);
	    }
	    else {
	      SpinorBarWaveFunction inters = _fer[idiag].first->
		evaluate(scale,widthOption(),offshell,
			 _outspin[dit->channelType].second[h1],_swave);
	      diag =  sign*_fer[idiag].second->
		evaluate(scale,_outspin[iferm].first [h2],inters, _outVector[v1]);
	    }
	  }
	  // intermediate vector
	  else if(offshell->iSpin() == PDT::Spin1) {
	    VectorWaveFunction interv = _vec[idiag].first->
	      evaluate(scale, widthOption(), offshell, _outVector[v1], _swave);
	    unsigned int h1(s1),h2(s2);
	    if(o2 > o3) swap(h1,h2);
	    if(decay[o2]->id() < 0 && decay[o3]->id() > 0) {
	      diag =-sign*_vec[idiag].second->
		evaluate(scale, _outspin[o2].first[h1],
			 _outspin[o3].second[h2], interv);
	    }
	    else {
	      diag = sign*_vec[idiag].second->
		evaluate(scale, _outspin[o3].first[h2],
			 _outspin[o2].second[h1], interv);
	    }
	  }
	  // unknown
	  else throw Exception()
	    << "Unknown intermediate in StoFFVDecayer::me2()" 
	    << Exception::runerror;

	  // matrix element for the different colour flows
	  if(ichan < 0) {
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
	      if(dit->colourFlow[iy].first - 1 != colourFlow()) continue;
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
	} //end of diagrams
	// now add the flows to the me2 with appropriate colour factors
	for(unsigned int ix = 0; ix < ncf; ++ix) {
	  if ( ivec == 0 ) { 
	    (*mes[ix])(0, v1, s1, s2) = flows[ix];
	    (*mel[ix])(0, v1, s1, s2) = largeflows[ix];
	  }
	  else if( ivec == 1 ) {
	    (*mes[ix])(0, s1, v1, s2) = flows[ix];
	    (*mel[ix])(0, s1, v1, s2) = largeflows[ix];
	  }
	  else if( ivec == 2 ) {
	    (*mes[ix])(0, s1, s2, v1) = flows[ix];
	    (*mel[ix])(0, s1, s2, v1) = largeflows[ix];
	  }
	}
      }
    }
  }
  double me2(0.);
  if(ichan < 0) {
    vector<double> pflows(ncf,0.);
    for(unsigned int ix = 0; ix < ncf; ++ix) {
      for(unsigned int iy = 0; iy < ncf; ++ iy) {
	double con = cfactors[ix][iy]*(mes[ix]->contract(*mes[iy],_rho)).real();
	me2 += con;
	if(ix == iy) {
	  con = nfactors[ix][iy]*(mel[ix]->contract(*mel[iy],_rho)).real();
	  pflows[ix] += con;
	}
      }
    }
    double ptotal(std::accumulate(pflows.begin(),pflows.end(),0.));
    ptotal *= UseRandom::rnd();
    for(unsigned int ix = 0;ix < pflows.size(); ++ix) {
      if(ptotal <= pflows[ix]) {
	colourFlow(ix);
	ME(mes[ix]);
	break;
      }
      ptotal -= pflows[ix];
    }
  }
  else {
    unsigned int iflow = colourFlow();
    me2 = nfactors[iflow][iflow]*(mel[iflow]->contract(*mel[iflow],_rho)).real();
  }
  // return the matrix element squared
  return me2;
}
