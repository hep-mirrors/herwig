// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VVVVVertex class.
//

#include "VVVVVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;

VVVVVertex::~VVVVVertex() {}
    
void VVVVVertex::persistentOutput(PersistentOStream & os) const {}

void VVVVVertex::persistentInput(PersistentIStream & is, int) {}
    
AbstractClassDescription<VVVVVertex> VVVVVertex::initVVVVVertex;
// Definition of the static class description member.
    
void VVVVVertex::Init() {
      
static ClassDocumentation<VVVVVertex> documentation
  ("The VVVVVertex class is the implementation of the 4-vector vertex");
 
}

// calculate the vertex
Complex VVVVVertex::evaluate(Energy2 q2 , int iopt, 
			     const VectorWaveFunction & vec1,
			     const VectorWaveFunction & vec2,
			     const VectorWaveFunction & vec3,
			     const VectorWaveFunction & vec4)
{
  // get the pointers to the particles
  tcPDPtr Pvec[4]={vec1.getParticle(),vec2.getParticle(),
		   vec3.getParticle(),vec4.getParticle()};
  // workout the coupling
  setCoupling(q2,Pvec[0],Pvec[1],Pvec[2],Pvec[3]);
  Complex norm=getNorm();
  Complex vertex,ii(0.,1.);
  // calculate the vertex
  // QCD type vertex
  if(_itype==1)
    {
      // dot products we need
      Complex dotv1v2 =
	+vec1.t()*vec2.t()-vec1.x()*vec2.x()
	-vec1.y()*vec2.y()-vec1.z()*vec2.z();
      Complex dotv3v4 =
	+vec3.t()*vec4.t()-vec3.x()*vec4.x()
	-vec3.y()*vec4.y()-vec3.z()*vec4.z();
      Complex dotv1v4 =
	+vec1.t()*vec4.t()-vec1.x()*vec4.x()
	-vec1.y()*vec4.y()-vec1.z()*vec4.z();
      Complex dotv2v3 =
	+vec3.t()*vec2.t()-vec3.x()*vec2.x()
	-vec3.y()*vec2.y()-vec3.z()*vec2.z();
      // first the 4-point part of the vertex
      vertex = dotv1v2*dotv3v4-dotv1v4*dotv2v3;
      // now the virtual gluon exchange if needed
      if(iopt!=0)
	{
	  // dot products
	  Complex dotv1v3 =
	    +vec1.t()*vec3.t()-vec1.x()*vec3.x()
	    -vec1.y()*vec3.y()-vec1.z()*vec3.z();
	  Complex dotv2v4 =
	    +vec4.t()*vec2.t()-vec4.x()*vec2.x()
	    -vec4.y()*vec2.y()-vec4.z()*vec2.z();
	  Complex dotv1p13 =
	    +vec1.t()*(2.*vec3.e() +vec1.e() )-vec1.x()*(2.*vec3.px()+vec1.px())
	    -vec1.y()*(2.*vec3.py()+vec1.py())-vec1.z()*(2.*vec3.pz()+vec1.pz()); 
	  Complex dotv2p24 =
	    +vec2.t()*(2.*vec4.e() +vec2.e() )-vec2.x()*(2.*vec4.px()+vec2.px())
	    -vec2.y()*(2.*vec4.py()+vec2.py())-vec2.z()*(2.*vec4.pz()+vec2.pz());
	  Complex dotv3p13 =
	    +vec3.t()*(2.*vec1.e() +vec3.e() )-vec3.x()*(2.*vec1.px()+vec3.px())
	    -vec3.y()*(2.*vec1.py()+vec3.py())-vec3.z()*(2.*vec1.pz()+vec3.pz());
	  Complex dotv4p24 =
	    +vec4.t()*(2.*vec2.e() +vec4.e() )-vec4.x()*(2.*vec2.px()+vec4.px())
	    -vec4.y()*(2.*vec2.py()+vec4.py())-vec4.z()*(2.*vec2.pz()+vec4.pz());
	  // construct the vectors
	  Complex veca[4],vecb[4];
	  for(int ix=0;ix<4;++ix)
	    {
	      veca[ix] = dotv3p13*vec1(ix)
		-dotv1p13*vec3(ix)+dotv1v3*(vec3[ix]-vec1[ix]);
	      vecb[ix] = dotv4p24*vec2(ix)
		-dotv2p24*vec4(ix)+dotv2v4*(vec4[ix]-vec2[ix]);
	    }
	  Energy2 numerator = 1./( 
				  +(vec1.e() +vec3.e() )*(vec1.e() +vec3.e() )
				  -(vec1.px()+vec3.px())*(vec1.px()+vec3.px())
				  -(vec1.py()+vec3.py())*(vec1.py()+vec3.py())
				  -(vec1.pz()+vec3.pz())*(vec1.pz()+vec3.pz()));
	  vertex = vertex+numerator*(+veca[3]*vecb[3]-veca[0]*vecb[0]
				     -veca[1]*vecb[1]-veca[2]*vecb[2]);
	}
      // final coupling factors
      vertex = -ii*norm*vertex;
    }
  // EW type vertex
  else if(_itype==2)
    {
      Complex dotv1v2 = 
	+vec1.t()*vec2.t()-vec1.x()*vec2.x()
	-vec1.y()*vec2.y()-vec1.z()*vec2.z();
      Complex dotv1v3 = 
	+vec1.t()*vec3.t()-vec1.x()*vec3.x()
	-vec1.y()*vec3.y()-vec1.z()*vec3.z();
      Complex dotv1v4 = 
	+vec1.t()*vec4.t()-vec1.x()*vec4.x()
	-vec1.y()*vec4.y()-vec1.z()*vec4.z();
      Complex dotv2v3 = 
	+vec2.t()*vec3.t()-vec2.x()*vec3.x()
	-vec2.y()*vec3.y()-vec2.z()*vec3.z();
      Complex dotv2v4 = 
	+vec2.t()*vec4.t()-vec2.x()*vec4.x()
	-vec2.y()*vec4.y()-vec2.z()*vec4.z();
      Complex dotv3v4 = 
	+vec3.t()*vec4.t()-vec3.x()*vec4.x()
	-vec3.y()*vec4.y()-vec3.z()*vec4.z();
      // evaluate the vertex
      // need to sort the order out here
      if(( _iorder[0]==0 && _iorder[1]==1 && _iorder[2]==2 && _iorder[3]==3)||
	 ( _iorder[0]==1 && _iorder[1]==0 && _iorder[2]==2 && _iorder[3]==3)||
	 ( _iorder[0]==0 && _iorder[1]==1 && _iorder[2]==3 && _iorder[3]==2)||
	 ( _iorder[0]==1 && _iorder[1]==0 && _iorder[2]==3 && _iorder[3]==2)||
	 ( _iorder[0]==2 && _iorder[1]==3 && _iorder[2]==0 && _iorder[3]==1)||
	 ( _iorder[0]==2 && _iorder[1]==3 && _iorder[2]==1 && _iorder[3]==0)||
	 ( _iorder[0]==3 && _iorder[1]==2 && _iorder[2]==0 && _iorder[3]==1)||
	 ( _iorder[0]==3 && _iorder[1]==2 && _iorder[2]==1 && _iorder[3]==0))
	{
	  // contact term
	  vertex = 2.*dotv1v2*dotv3v4-dotv1v3*dotv2v4-dotv1v4*dotv2v3;
	  // now for the u- and t-channel terms if needed
	  if(iopt!=0)
	    {
	      // dot products of momenta and wavefunction
	      Complex dotv1p13 =
		+vec1.t()*(vec1.e() +2.*vec3.e() )-vec1.x()*(vec1.px()+2.*vec3.px())
		-vec1.y()*(vec1.py()+2.*vec3.py())-vec1.z()*(vec1.pz()+2.*vec3.pz());
	      Complex dotv1p14 =
		+vec1.t()*(vec1.e() +2.*vec4.e() )-vec1.x()*(vec1.px()+2.*vec4.px())
		-vec1.y()*(vec1.py()+2.*vec4.py())-vec1.z()*(vec1.pz()+2.*vec4.pz());
	      Complex dotv2p23 =
		+vec2.t()*(vec2.e() +2.*vec3.e() )-vec2.x()*(vec2.px()+2.*vec3.px())
		-vec2.y()*(vec2.py()+2.*vec3.py())-vec2.z()*(vec2.pz()+2.*vec3.pz());
	      Complex dotv2p24 =
		+vec2.t()*(vec2.e() +2.*vec4.e() )-vec2.x()*(vec2.px()+2.*vec4.px())
		-vec2.y()*(vec2.py()+2.*vec4.py())-vec2.z()*(vec2.pz()+2.*vec4.pz());
	      Complex dotv3p31 = 
		+vec3.t()*(vec3.e() +2.*vec1.e() )-vec3.x()*(vec3.px()+2.*vec1.px())
		-vec3.y()*(vec3.py()+2.*vec1.py())-vec3.z()*(vec3.pz()+2.*vec1.pz());
	      Complex dotv3p32 = 
		+vec3.t()*(vec3.e() +2.*vec2.e() )-vec3.x()*(vec3.px()+2.*vec2.px())
		-vec3.y()*(vec3.py()+2.*vec2.py())-vec3.z()*(vec3.pz()+2.*vec2.pz());
	      Complex dotv4p41 = 
		+vec4.t()*(vec4.e() +2.*vec1.e() )-vec4.x()*(vec4.px()+2.*vec1.px())
		-vec4.y()*(vec4.py()+2.*vec1.py())-vec4.z()*(vec4.pz()+2.*vec1.pz());
	      Complex dotv4p42 = 
		+vec4.t()*(vec4.e() +2.*vec2.e() )-vec4.x()*(vec4.px()+2.*vec2.px())
		-vec4.y()*(vec4.py()+2.*vec2.py())-vec4.z()*(vec4.pz()+2.*vec2.pz());
	      Complex ja[4],jb[4],jc[4],jd[4];
	      // vectors
	      for(int ix=0;ix<4;++ix)
		{
		  ja[ix] =(vec3[ix]-vec1[ix])*dotv1v3
		    +dotv3p31*vec1(ix)-dotv1p13*vec3(ix);
		  jc[ix] =(vec4[ix]-vec1[ix])*dotv1v4
		    +dotv4p41*vec1(ix)-dotv1p14*vec4(ix);
		  jb[ix] =(vec4[ix]-vec2[ix])*dotv2v4
		    +dotv4p42*vec2(ix)-dotv2p24*vec4(ix);
		  jd[ix] =(vec3[ix]-vec2[ix])*dotv2v3
		    +dotv3p32*vec2(ix)-dotv2p23*vec3(ix);
		}
	      // dot products of these vectors
	      Complex dotjajb = ja[3]*jb[3]
		-ja[0]*jb[0]-ja[1]*jb[1]-ja[2]*jb[2];
	      Complex dotjcjd = jc[3]*jd[3]
		-jc[0]*jd[0]-jc[1]*jd[1]-jc[2]*jd[2];
	      Complex dotjaq = 
		+ja[3]*(vec1.e() +vec3.e() )-ja[0]*(vec1.px()+vec3.px())
		-ja[1]*(vec1.py()+vec3.py())-ja[2]*(vec1.pz()+vec3.pz());
	      Complex dotjbq = 
		+jb[3]*(vec1.e() +vec3.e() )-jb[0]*(vec1.px()+vec3.px())
		-jb[1]*(vec1.py()+vec3.py())-jb[2]*(vec1.pz()+vec3.pz());
	      Complex dotjck = 
		+jc[3]*(vec1.e() +vec4.e() )-jc[0]*(vec1.px()+vec4.px())
		-jc[1]*(vec1.py()+vec4.py())-jc[2]*(vec1.pz()+vec4.pz());
	      Complex dotjdk = 
		+jd[3]*(vec1.e() +vec4.e() )-jd[0]*(vec1.px()+vec4.px())
		-jd[1]*(vec1.py()+vec4.py())-jd[2]*(vec1.pz()+vec4.pz());
	      Energy2 q2 = 
		+(vec1.e() +vec3.e() )*(vec1.e() +vec3.e() )
		-(vec1.px()+vec3.px())*(vec1.px()+vec3.px())
		-(vec1.py()+vec3.py())*(vec1.py()+vec3.py())
		-(vec1.pz()+vec3.pz())*(vec1.pz()+vec3.pz());
	      Energy2 k2 = 
		+(vec1.e() +vec4.e() )*(vec1.e() +vec4.e() )
		-(vec1.px()+vec4.px())*(vec1.px()+vec4.px())
		-(vec1.py()+vec4.py())*(vec1.py()+vec4.py())
		-(vec1.pz()+vec4.pz())*(vec1.pz()+vec4.pz());
	      // compute the term we need
	      Energy mass2;
	      for(int ix=0;ix<2;++ix)
		{
		  if(_inter[ix])
		    {
		      mass2 = (_inter[ix]->mass())*(_inter[ix]->mass());
		      if(mass2!=0.)
			{
			  vertex=vertex+_coup[ix]*propagator(iopt,q2,_inter[ix])*
			    (dotjajb-dotjaq*dotjbq/mass2);
			  vertex=vertex+_coup[ix]*propagator(iopt,k2,_inter[ix])*
			    (dotjcjd-dotjck*dotjdk/mass2);
			}
		      else
			{
			  vertex+=_coup[ix]*propagator(iopt,q2,_inter[ix])*dotjajb;
			  vertex+=_coup[ix]*propagator(iopt,k2,_inter[ix])*dotjcjd;
			}
		    }
		}
	    }
	}
      else if(( _iorder[0]==0 && _iorder[1]==2 && _iorder[2]==1 && _iorder[3]==3)||
	      ( _iorder[0]==2 && _iorder[1]==0 && _iorder[2]==1 && _iorder[3]==3)||
	      ( _iorder[0]==0 && _iorder[1]==2 && _iorder[2]==3 && _iorder[3]==1)||
	      ( _iorder[0]==2 && _iorder[1]==0 && _iorder[2]==3 && _iorder[3]==1)||
	      ( _iorder[0]==1 && _iorder[1]==3 && _iorder[2]==0 && _iorder[3]==2)||
	      ( _iorder[0]==1 && _iorder[1]==3 && _iorder[2]==2 && _iorder[3]==0)||
	      ( _iorder[0]==3 && _iorder[1]==1 && _iorder[2]==0 && _iorder[3]==2)||
	      ( _iorder[0]==3 && _iorder[1]==1 && _iorder[2]==2 && _iorder[3]==0))
	{
	  // contact term
	  vertex = 2.*dotv1v3*dotv2v4-dotv1v2*dotv3v4-dotv1v4*dotv2v3;
	  // now for the u- and t-channel terms if needed
	  if(iopt!=0)
	    {
	      // dot products of momenta and wavefunction
	      Complex dotv1p12 =
		+vec1.t()*(vec1.e() +2.*vec2.e() )-vec1.x()*(vec1.px()+2.*vec2.px())
		-vec1.y()*(vec1.py()+2.*vec2.py())-vec1.z()*(vec1.pz()+2.*vec2.pz());
	      Complex dotv1p14 =
		+vec1.t()*(vec1.e() +2.*vec4.e() )-vec1.x()*(vec1.px()+2.*vec4.px())
		-vec1.y()*(vec1.py()+2.*vec4.py())-vec1.z()*(vec1.pz()+2.*vec4.pz());
	      Complex dotv3p32 =
		+vec3.t()*(vec3.e() +2.*vec2.e() )-vec3.x()*(vec3.px()+2.*vec2.px())
		-vec3.y()*(vec3.py()+2.*vec2.py())-vec3.z()*(vec3.pz()+2.*vec2.pz());
	      Complex dotv3p34 =
		+vec3.t()*(vec3.e() +2.*vec4.e() )-vec3.x()*(vec3.px()+2.*vec4.px())
		-vec3.y()*(vec3.py()+2.*vec4.py())-vec3.z()*(vec3.pz()+2.*vec4.pz());
	      Complex dotv2p21 = 
		+vec2.t()*(vec2.e() +2.*vec1.e() )-vec2.x()*(vec2.px()+2.*vec1.px())
		-vec2.y()*(vec2.py()+2.*vec1.py())-vec2.z()*(vec2.pz()+2.*vec1.pz());
	      Complex dotv2p23 = 
		+vec2.t()*(vec2.e() +2.*vec3.e() )-vec2.x()*(vec2.px()+2.*vec3.px())
		-vec2.y()*(vec2.py()+2.*vec3.py())-vec2.z()*(vec2.pz()+2.*vec3.pz());
	      Complex dotv4p41 = 
		+vec4.t()*(vec4.e() +2.*vec1.e() )-vec4.x()*(vec4.px()+2.*vec1.px())
		-vec4.y()*(vec4.py()+2.*vec1.py())-vec4.z()*(vec4.pz()+2.*vec1.pz());
	      Complex dotv4p43 = 
		+vec4.t()*(vec4.e() +2.*vec3.e() )-vec4.x()*(vec4.px()+2.*vec3.px())
		-vec4.y()*(vec4.py()+2.*vec3.py())-vec4.z()*(vec4.pz()+2.*vec3.pz());
	      Complex ja[4],jb[4],jc[4],jd[4];
	      // vectors
	      for(int ix=0;ix<4;++ix)
		{
		  ja[ix] =(vec2[ix]-vec1[ix])*dotv1v2
		    +dotv2p21*vec1(ix)-dotv1p12*vec2(ix);
		  jc[ix] =(vec4[ix]-vec1[ix])*dotv1v4
		    +dotv4p41*vec1(ix)-dotv1p14*vec4(ix);
		  jb[ix] =(vec4[ix]-vec3[ix])*dotv3v4
		    +dotv4p43*vec3(ix)-dotv3p34*vec4(ix);
		  jd[ix] =(vec2[ix]-vec3[ix])*dotv2v3
		    +dotv2p23*vec3(ix)-dotv3p32*vec2(ix);
		}
	      // dot products of these vectors
	      Complex dotjajb = ja[3]*jb[3]
		-ja[0]*jb[0]-ja[1]*jb[1]-ja[2]*jb[2];
	      Complex dotjcjd = jc[3]*jd[3]
		-jc[0]*jd[0]-jc[1]*jd[1]-jc[2]*jd[2];
	      Complex dotjaq = 
		+ja[3]*(vec1.e() +vec2.e() )-ja[0]*(vec1.px()+vec2.px())
		-ja[1]*(vec1.py()+vec2.py())-ja[2]*(vec1.pz()+vec2.pz());
	      Complex dotjbq = 
		+jb[3]*(vec1.e() +vec2.e() )-jb[0]*(vec1.px()+vec2.px())
		-jb[1]*(vec1.py()+vec2.py())-jb[2]*(vec1.pz()+vec2.pz());
	      Complex dotjck = 
		+jc[3]*(vec1.e() +vec4.e() )-jc[0]*(vec1.px()+vec4.px())
		-jc[1]*(vec1.py()+vec4.py())-jc[2]*(vec1.pz()+vec4.pz());
	      Complex dotjdk = 
		+jd[3]*(vec1.e() +vec4.e() )-jd[0]*(vec1.px()+vec4.px())
		-jd[1]*(vec1.py()+vec4.py())-jd[2]*(vec1.pz()+vec4.pz());
	      Energy2 q2 = 
		+(vec1.e() +vec2.e() )*(vec1.e() +vec2.e() )
		-(vec1.px()+vec2.px())*(vec1.px()+vec2.px())
		-(vec1.py()+vec2.py())*(vec1.py()+vec2.py())
		-(vec1.pz()+vec2.pz())*(vec1.pz()+vec2.pz());
	      Energy2 k2 = 
		+(vec1.e() +vec4.e() )*(vec1.e() +vec4.e() )
		-(vec1.px()+vec4.px())*(vec1.px()+vec4.px())
		-(vec1.py()+vec4.py())*(vec1.py()+vec4.py())
		-(vec1.pz()+vec4.pz())*(vec1.pz()+vec4.pz());
	      // compute the term we need
	      Energy2 mass2;
	      for(int ix=0;ix<2;++ix)
		{
		  if(_inter[ix])
		    {
		      mass2 = (_inter[ix]->mass())*(_inter[ix]->mass());
		      if(mass2!=0.)
			{
			  vertex+=_coup[ix]*propagator(iopt,q2,_inter[ix])*
			    (dotjajb-dotjaq*dotjbq/mass2);
			  vertex+=_coup[ix]*propagator(iopt,k2,_inter[ix])*
			    (dotjcjd-dotjck*dotjdk/mass2);
			}
		      else
			{
			  vertex+=_coup[ix]*propagator(iopt,q2,_inter[ix])*dotjajb;
			  vertex+=_coup[ix]*propagator(iopt,k2,_inter[ix])*dotjcjd;
			}
		    }
		}
	    }
	}
      else if(( _iorder[0]==0 && _iorder[1]==3 && _iorder[2]==1 && _iorder[3]==2)||
	      ( _iorder[0]==3 && _iorder[1]==0 && _iorder[2]==1 && _iorder[3]==2)||
	      ( _iorder[0]==0 && _iorder[1]==3 && _iorder[2]==2 && _iorder[3]==1)||
	      ( _iorder[0]==3 && _iorder[1]==0 && _iorder[2]==2 && _iorder[3]==1)||
	      ( _iorder[0]==1 && _iorder[1]==2 && _iorder[2]==0 && _iorder[3]==3)||
	      ( _iorder[0]==2 && _iorder[1]==1 && _iorder[2]==0 && _iorder[3]==3)||
	      ( _iorder[0]==1 && _iorder[1]==2 && _iorder[2]==3 && _iorder[3]==0)||
	      ( _iorder[0]==2 && _iorder[1]==1 && _iorder[2]==3 && _iorder[3]==0))
	{
	  // contact term
	  vertex = 2.*dotv1v4*dotv2v3-dotv1v3*dotv2v4-dotv1v2*dotv3v4;
	  // now for the u- and t-channel terms if needed
	  if(iopt!=0)
	    {
	      // dot products of momenta and wavefunction
	      Complex dotv1p12 =
		+vec1.t()*(vec1.e() +2.*vec2.e() )-vec1.x()*(vec1.px()+2.*vec2.px())
		-vec1.y()*(vec1.py()+2.*vec2.py())-vec1.z()*(vec1.pz()+2.*vec2.pz());
	      Complex dotv1p13 =
		+vec1.t()*(vec1.e() +2.*vec3.e() )-vec1.x()*(vec1.px()+2.*vec3.px())
		-vec1.y()*(vec1.py()+2.*vec3.py())-vec1.z()*(vec1.pz()+2.*vec3.pz());
	      Complex dotv2p24 = 
		+vec2.t()*(vec2.e() +2.*vec4.e() )-vec2.x()*(vec2.px()+2.*vec4.px())
		-vec2.y()*(vec2.py()+2.*vec4.py())-vec2.z()*(vec2.pz()+2.*vec4.pz());
	      Complex dotv2p21 = 
		+vec2.t()*(vec2.e() +2.*vec1.e() )-vec2.x()*(vec2.px()+2.*vec1.px())
		-vec2.y()*(vec2.py()+2.*vec1.py())-vec2.z()*(vec2.pz()+2.*vec1.pz());
	      Complex dotv3p31 = 
		+vec3.t()*(vec3.e() +2.*vec1.e() )-vec3.x()*(vec3.px()+2.*vec1.px())
		-vec3.y()*(vec3.py()+2.*vec1.py())-vec3.z()*(vec3.pz()+2.*vec1.pz());
	      Complex dotv3p34 =
		+vec3.t()*(vec3.e() +2.*vec4.e() )-vec3.x()*(vec3.px()+2.*vec4.px())
		-vec3.y()*(vec3.py()+2.*vec4.py())-vec3.z()*(vec3.pz()+2.*vec4.pz());
	      Complex dotv4p43 = 
		+vec4.t()*(vec4.e() +2.*vec3.e() )-vec4.x()*(vec4.px()+2.*vec3.px())
		-vec4.y()*(vec4.py()+2.*vec3.py())-vec4.z()*(vec4.pz()+2.*vec3.pz());
	      Complex dotv4p42 =
		+vec4.t()*(vec4.e() +2.*vec2.e() )-vec4.x()*(vec4.px()+2.*vec2.px())
		-vec4.y()*(vec4.py()+2.*vec2.py())-vec4.z()*(vec4.pz()+2.*vec2.pz());
	      Complex ja[4],jb[4],jc[4],jd[4];
	      // vectors
	      for(int ix=0;ix<4;++ix)
		{
		  ja[ix] =(vec2[ix]-vec1[ix])*dotv1v2
		    +dotv2p21*vec1(ix)-dotv1p12*vec2(ix);
		  jc[ix] =(vec3[ix]-vec1[ix])*dotv1v3
		    +dotv3p31*vec1(ix)-dotv1p13*vec3(ix);
		  jb[ix] =(vec3[ix]-vec4[ix])*dotv3v4
		    +dotv3p34*vec4(ix)-dotv4p43*vec3(ix);
		  jd[ix] =(vec2[ix]-vec4[ix])*dotv2v4
		    +dotv2p24*vec4(ix)-dotv4p42*vec2(ix);
		}
	      // dot products of these vectors
	      Complex dotjajb = ja[3]*jb[3]
		-ja[0]*jb[0]-ja[1]*jb[1]-ja[2]*jb[2];
	      Complex dotjcjd = jc[3]*jd[3]
		-jc[0]*jd[0]-jc[1]*jd[1]-jc[2]*jd[2];
	      Complex dotjaq = 
		+ja[3]*(vec1.e() +vec2.e() )-ja[0]*(vec1.px()+vec2.px())
		-ja[1]*(vec1.py()+vec2.py())-ja[2]*(vec1.pz()+vec2.pz());
	      Complex dotjbq = 
		+jb[3]*(vec1.e() +vec2.e() )-jb[0]*(vec1.px()+vec2.px())
		-jb[1]*(vec1.py()+vec2.py())-jb[2]*(vec1.pz()+vec2.pz());
	      Complex dotjck = 
		+jc[3]*(vec1.e() +vec3.e() )-jc[0]*(vec1.px()+vec3.px())
		-jc[1]*(vec1.py()+vec3.py())-jc[2]*(vec1.pz()+vec3.pz());
	      Complex dotjdk = 
		+jd[3]*(vec1.e() +vec3.e() )-jd[0]*(vec1.px()+vec3.px())
		-jd[1]*(vec1.py()+vec3.py())-jd[2]*(vec1.pz()+vec3.pz());
	      Energy2 q2 = 
		+(vec1.e() +vec2.e() )*(vec1.e() +vec2.e() )
		-(vec1.px()+vec2.px())*(vec1.px()+vec2.px())
		-(vec1.py()+vec2.py())*(vec1.py()+vec2.py())
		-(vec1.pz()+vec2.pz())*(vec1.pz()+vec2.pz());
	      Energy2 k2 = 
		+(vec1.e() +vec3.e() )*(vec1.e() +vec3.e() )
		-(vec1.px()+vec3.px())*(vec1.px()+vec3.px())
		-(vec1.py()+vec3.py())*(vec1.py()+vec3.py())
		-(vec1.pz()+vec3.pz())*(vec1.pz()+vec3.pz());
	      // compute the term we need
	      Energy2 mass2;
	      for(int ix=0;ix<2;++ix)
		{
		  if(_inter[ix])
		    {
		      mass2 =(_inter[ix]->mass())*(_inter[ix]->mass());
		      if(mass2!=0.)
			{
			  vertex=vertex+_coup[ix]*propagator(iopt,q2,_inter[ix])*
			    (dotjajb-dotjaq*dotjbq/mass2);
			  vertex=vertex+_coup[ix]*propagator(iopt,k2,_inter[ix])*
			    (dotjcjd-dotjck*dotjdk/mass2);
			}
		      else
			{
			  vertex=vertex+_coup[ix]*propagator(iopt,q2,_inter[ix])*dotjajb;
			  vertex=vertex+_coup[ix]*propagator(iopt,k2,_inter[ix])*dotjcjd;
			}
		    }
		}
	    }
	}
      else
	{std::cerr << "Unknown order of particles in VVVV Vertex"<< std::endl;}
      vertex=-ii*norm*vertex;
    }
  else
    {
      std::cerr << "unknown type of VVVV Vertex" << std::endl;
      vertex=0.;
    }
  // return the answer
  return vertex;
}

}
}

