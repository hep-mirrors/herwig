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
    ("There is no documentation for the FtoFFFDecayer class");

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
      FFSVertexPtr vert1 = dynamic_ptr_cast<FFSVertexPtr>
	(current.vertices.first);
      FFSVertexPtr vert2 = dynamic_ptr_cast<FFSVertexPtr>
	(current.vertices.second);
      _sca[ix] = make_pair(vert1, vert2);
    }
    else if(offshell->iSpin() == PDT::Spin1) {
      FFVVertexPtr vert1 = dynamic_ptr_cast<FFVVertexPtr>
	(current.vertices.first);
      FFVVertexPtr vert2 = dynamic_ptr_cast<FFVVertexPtr>
	(current.vertices.second);
      _vec[ix] = make_pair(vert1, vert2);
    }
    else if(offshell->iSpin() == PDT::Spin2) {
      FFTVertexPtr vert1 = dynamic_ptr_cast<FFTVertexPtr>
	(current.vertices.first);
      FFTVertexPtr vert2 = dynamic_ptr_cast<FFTVertexPtr>
	(current.vertices.second);
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
  static unsigned int out2[3]={1,0,0},out3[3]={2,2,1};
  vector<DecayMatrixElement> mes(ncf,DecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,
							PDT::Spin1Half,PDT::Spin1Half));
  vector<double> pflows(ncf,0.);
  double me2(0.);
  for(unsigned int s0 = 0; s0 < 2; ++s0) {
    for(unsigned int s1 = 0; s1 < 2; ++s1) {
      for(unsigned int s2 = 0; s2 < 2; ++s2) {
	for(unsigned int s3 = 0; s3 < 2; ++s3) {
	  flows = vector<Complex>(ncf, Complex(0.));
	  int nchan=-1;
	  unsigned int idiag=0;
	  for(vector<TBDiagram>::const_iterator dit=getProcessInfo().begin();
	      dit!=getProcessInfo().end();++dit) {
	    // the sign from normal ordering
	    double sign = ferm ? 1. : -1;
	    // outgoing wavefunction and NO sign
	    if     (dit->channelType==TBDiagram::channel23) sign *= -1.;
	    else if(dit->channelType==TBDiagram::channel13) sign *= -1.;
	    else if(dit->channelType==TBDiagram::channel12) sign *=  1.;
	    else throw Exception()
	      << "Unknown diagram type in FtoFFFDecayer::me2()" << Exception::runerror;
	    // wavefunctions
	    SpinorWaveFunction    w0,w3;
	    SpinorBarWaveFunction w1,w2;
	    // incoming wavefunction
	    if(ferm) {
	      w0 = inwave.first [s0];
	      w1 = outwave[dit->channelType].second[s1];
	    }
	    else {
	      w0 = outwave[dit->channelType].first [s1];
	      w1 = inwave.second[s0];
	    }
	    if(decay[out2[dit->channelType]]->id()<0&&
	       decay[out3[dit->channelType]]->id()>0) {
	      w2 = outwave[out3[dit->channelType]].second[s3];
	      w3 = outwave[out2[dit->channelType]].first [s2];
	    }
	    else {
	      w2 = outwave[out2[dit->channelType]].second[s2];
	      w3 = outwave[out3[dit->channelType]].first [s3];
	    }
	    // channels if selecting
	    ++nchan;
	    if(ichan>0&&ichan!=nchan) continue;
	    tcPDPtr offshell = dit->intermediate;
	    // intermediate scalar
	    Complex diag;
	    if     (offshell->iSpin() == PDT::Spin0) {
	      ScalarWaveFunction inters = _sca[idiag].first->
		evaluate(scale, 1, offshell, w0, w1);
	      diag = _sca[idiag].second->evaluate(scale,w3,w2,inters);
	    }
	    // intermediate vector
	    else if(offshell->iSpin() == PDT::Spin1) {
	      VectorWaveFunction interv = _vec[idiag].first->
		evaluate(scale, 1, offshell, w0, w1);
	      diag = _vec[idiag].second->evaluate(scale,w3,w2,interv);
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
	    // matrix element for the different colour flows
	    for(unsigned iy = 0; iy < dit->colourFlow.size(); ++iy) {
	      flows[dit->colourFlow[iy].first - 1] += 
		dit->colourFlow[iy].second * diag;
	    }
	    ++idiag;
	  }
	  // now add the flows to the me2 with appropriate colour factors
	  for(unsigned int ix = 0; ix < ncf; ++ix) {
	    mes[ix](s0,s1,s2,s3) = cfactors[ix][ix]*flows[ix];
	    pflows[ix] += cfactors[ix][ix]*norm(flows[ix]);
	    for(unsigned int iy = 0; iy < ncf; ++iy) {
	      me2 += cfactors[ix][iy]*(flows[ix]*conj(flows[iy])).real();
	    }
	  }
	}
      }
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
  // return the matrix element squared
  return me2;
}

WidthCalculatorBasePtr FtoFFFDecayer::
threeBodyMEIntegrator(const DecayMode & ) const {
  vector<int> intype;
  vector<Energy> inmass,inwidth;
  vector<double> inpow;
  int nchannel(0);
  for(unsigned int ix=0;ix<getProcessInfo().size();++ix) {
    if(getProcessInfo()[ix].channelType==TBDiagram::fourPoint) continue;
    else if(getProcessInfo()[ix].channelType==TBDiagram::channel23) intype.push_back(3);
    else if(getProcessInfo()[ix].channelType==TBDiagram::channel13) intype.push_back(2);
    else if(getProcessInfo()[ix].channelType==TBDiagram::channel12) intype.push_back(1);
    if(getProcessInfo()[ix].intermediate->id()!=ParticleID::gamma) {
      inpow.push_back(0.);
      inmass.push_back(getProcessInfo()[ix].intermediate->mass());
      inwidth.push_back(getProcessInfo()[ix].intermediate->width());
    }
    else {
      inpow.push_back(-2.);
      inmass.push_back(-1.*GeV);
      inwidth.push_back(-1.*GeV);
    }
    ++nchannel;
  }
  vector<double> inweights(nchannel,1./double(nchannel));
  return new_ptr(ThreeBodyAllOnCalculator<FtoFFFDecayer>
		 (inweights,intype,inmass,inwidth,inpow,*this,0,
		  outgoing()[0]->mass(),outgoing()[1]->mass(),outgoing()[2]->mass()));
}
