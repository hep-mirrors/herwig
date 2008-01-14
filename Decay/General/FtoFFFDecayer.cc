// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FtoFFFDecayer class.
//

#include "FtoFFFDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Decay/DecayPhaseSpaceMode.h"
#include "Herwig++/PDT/ThreeBodyAllOnCalculator.h"
#include <numeric>

using namespace Herwig;

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

void FtoFFFDecayer::doinit() throw(InitException) {
  GeneralThreeBodyDecayer::doinit();
  unsigned int ndiags = getProcessInfo().size();
  _sca.resize(ndiags);
  _vec.resize(ndiags);
  _ten.resize(ndiags);
  for(unsigned int ix = 0;ix < ndiags; ++ix) {
    TBDiagram current = getProcessInfo()[ix];
    tcPDPtr offshell = current.intermediate;
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

double  FtoFFFDecayer::me2(bool vertex, const int ichan, const Particle & inpart,
			   const ParticleVector & decay) const {
  const vector<vector<double> > cfactors(getColourFactors());
  Energy2 scale(sqr(inpart.mass()));
  // spin density matrix
  RhoDMatrix rhoin(PDT::Spin1Half);
  rhoin.average();
  // get the wavefunctions for all the particles
  pair<vector<SpinorWaveFunction>,vector<SpinorBarWaveFunction> > inwave;
  pair<vector<SpinorWaveFunction>,vector<SpinorBarWaveFunction> > outwave[3];
  // incoming particle
  SpinorWaveFunction(inwave.first,rhoin,const_ptr_cast<tPPtr>(&inpart),
		     Helicity::incoming,true,vertex);
  if(inwave.first[0].wave().Type() == u_spinortype) {
    for(unsigned int ix = 0; ix < 2; ++ix) {
      inwave.second.push_back(inwave.first[ix].bar());
      inwave.second[ix].conjugate();
    }
  }
  else {
    for(unsigned int ix = 0; ix < 2; ++ix) {
      inwave.second.push_back(inwave.first[ix].bar());
      inwave.first[ix].conjugate();
    }
  }
  // outgoing particles
  for(unsigned int ix = 0; ix < 3; ++ix) {
    SpinorWaveFunction(outwave[ix].first,decay[ix],Helicity::outgoing,true,vertex);
    if(outwave[ix].first[0].wave().Type() == u_spinortype) {
      for(unsigned int iy = 0; iy < 2; ++iy) {
	outwave[ix].second.push_back(outwave[ix].first[iy].bar());
	outwave[ix].first[iy].conjugate();
      }
    }
    else {
      for(unsigned int iy = 0; iy < 2; ++iy) {
	outwave[ix].second.push_back(outwave[ix].first[iy].bar());
	outwave[ix].second[iy].conjugate();
      }
    }
  }
  bool ferm = inpart.id()>0; 
  const size_t ncf(numberOfFlows());
  vector<Complex> flows(ncf, Complex(0.)); 
  static const unsigned int out2[3]={1,0,0},out3[3]={2,2,1};
  vector<DecayMatrixElement> mes(ncf,DecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,
							PDT::Spin1Half,PDT::Spin1Half));
  unsigned int ihel[4];
  for(ihel[0] = 0; ihel[0] < 2; ++ihel[0]) {
    for(ihel[1] = 0; ihel[1] < 2; ++ihel[1]) {
      for(ihel[2] = 0; ihel[2] < 2; ++ihel[2]) {
	for(ihel[3] = 0; ihel[3] < 2; ++ihel[3]) {
	  flows = vector<Complex>(ncf, Complex(0.));
	  int nchan=-1;
	  unsigned int idiag=0;
	  for(vector<TBDiagram>::const_iterator dit=getProcessInfo().begin();
	      dit!=getProcessInfo().end();++dit) {
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
	      w0 = inwave.first [ihel[0]];
	      w1 = outwave[dit->channelType].second[ihel[dit->channelType+1]];
	    }
	    else {
	      w0 = outwave[dit->channelType].first [ihel[dit->channelType+1]];
	      w1 = inwave.second[ihel[0]];
	    }
	    if(decay[out2[dit->channelType]]->id()<0&&
	       decay[out3[dit->channelType]]->id()>0) {
	      w2 = outwave[out3[dit->channelType]].second[ihel[out3[dit->channelType]+1]];
	      w3 = outwave[out2[dit->channelType]].first [ihel[out2[dit->channelType]+1]];
	      sign *= -1.;
	    }
	    else {
	      w2 = outwave[out2[dit->channelType]].second[ihel[out2[dit->channelType]+1]];
	      w3 = outwave[out3[dit->channelType]].first [ihel[out3[dit->channelType]+1]];
	    }
	    // channels if selecting
	    ++nchan;
	    if(ichan>0&&ichan!=nchan) continue;
	    tcPDPtr offshell = dit->intermediate;
	    Complex diag;
	    // intermediate scalar
	    if     (offshell->iSpin() == PDT::Spin0) { 
	      ScalarWaveFunction inters = _sca[idiag].first->
		evaluate(scale, 1, offshell, w0, w1);
	      diag = _sca[idiag].second->evaluate(scale,w3,w2,inters);
	    }
	    // intermediate vector
	    else if(offshell->iSpin() == PDT::Spin1) {
	      VectorWaveFunction interv = _vec[idiag].first->
		evaluate(scale, 1, offshell, w0, w1);
	      diag = -_vec[idiag].second->evaluate(scale,w3,w2,interv);
	    }
	    // intermediate tensor
	    else if(offshell->iSpin() == PDT::Spin2) {
	      TensorWaveFunction intert = _ten[idiag].first->
		evaluate(scale, 1, offshell, w0, w1);
	      diag = _ten[idiag].second->evaluate(scale,w3,w2,intert);
	    }
	    // unknown
	    else throw Exception()
	      << "Unknown intermediate in FtoFFFDecayer::me2()" 
	      << Exception::runerror;
	    // apply NO sign
	    diag *= sign;
	    // matrix element for the different colour flows
	    for(unsigned iy = 0; iy < dit->colourFlow.size(); ++iy) {
	      flows[dit->colourFlow[iy].first - 1] += 
		dit->colourFlow[iy].second * diag;
	    }
	    ++idiag;
	  }
	  // now add the flows to the me2 with appropriate colour factors
	  for(unsigned int ix = 0; ix < ncf; ++ix) {
	    mes[ix](ihel[0],ihel[1],ihel[2],ihel[3]) = flows[ix];
	  }
	}
      }
    }
  }
  vector<double> pflows(ncf,0.);
  double me2(0.);
  for(unsigned int ix = 0; ix < ncf; ++ix) {
    for(unsigned int iy = 0; iy < ncf; ++ iy) {
      double con = cfactors[ix][iy]*(mes[ix].contract(mes[iy],rhoin)).real();
      if(ix==iy) pflows[ix] += con;
      me2 += con;
    }
  }
  // select the matrix element according to the colour flow
  if(vertex) {
    double ptotal(std::accumulate(pflows.begin(),pflows.end(),0.));
    ptotal *=UseRandom::rnd();
    for(unsigned int ix=0;ix<pflows.size();++ix) {
      if(ptotal<=pflows[ix]) {
	ME(mes[ix]);
	break;
      }
      ptotal-=pflows[ix];
    }
  }
  // make the colour connections
  colourConnections(inpart, decay);
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
		  outgoing()[0]->mass(),outgoing()[1]->mass(),outgoing()[2]->mass()));
}
