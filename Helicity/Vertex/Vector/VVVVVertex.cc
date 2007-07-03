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
    
AbstractNoPIOClassDescription<VVVVVertex> VVVVVertex::initVVVVVertex;
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
	  complex<Energy> dotv1p13 =
	    +vec1.t()*(2.*vec3.e() +vec1.e() )
	    -vec1.x()*(2.*vec3.px()+vec1.px())
	    -vec1.y()*(2.*vec3.py()+vec1.py())
	    -vec1.z()*(2.*vec3.pz()+vec1.pz()); 
	  complex<Energy> dotv2p24 =
	    +vec2.t()*(2.*vec4.e() +vec2.e() )
	    -vec2.x()*(2.*vec4.px()+vec2.px())
	    -vec2.y()*(2.*vec4.py()+vec2.py())
	    -vec2.z()*(2.*vec4.pz()+vec2.pz());
	  complex<Energy> dotv3p13 =
	    +vec3.t()*(2.*vec1.e() +vec3.e() )
	    -vec3.x()*(2.*vec1.px()+vec3.px())
	    -vec3.y()*(2.*vec1.py()+vec3.py())
	    -vec3.z()*(2.*vec1.pz()+vec3.pz());
	  complex<Energy> dotv4p24 =
	    +vec4.t()*(2.*vec2.e() +vec4.e() )
	    -vec4.x()*(2.*vec2.px()+vec4.px())
	    -vec4.y()*(2.*vec2.py()+vec4.py())
	    -vec4.z()*(2.*vec2.pz()+vec4.pz());
	  // construct the vectors
	  complex<Energy> veca[4],vecb[4];

	  veca[0] = dotv3p13*vec1.x()
	    -dotv1p13*vec3.x()+dotv1v3*(vec3.px()-vec1.px());

	  veca[1] = dotv3p13*vec1.y()
	    -dotv1p13*vec3.y()+dotv1v3*(vec3.py()-vec1.py());

	  veca[2] = dotv3p13*vec1.z()
	    -dotv1p13*vec3.z()+dotv1v3*(vec3.pz()-vec1.pz());

	  veca[3] = dotv3p13*vec1.t()
	    -dotv1p13*vec3.t()+dotv1v3*(vec3.e()-vec1.e());


	  vecb[0] = dotv4p24*vec2.x()
	    -dotv2p24*vec4.x()+dotv2v4*(vec4.px()-vec2.px());

	  vecb[1] = dotv4p24*vec2.y()
	    -dotv2p24*vec4.y()+dotv2v4*(vec4.py()-vec2.py());

	  vecb[2] = dotv4p24*vec2.z()
	    -dotv2p24*vec4.z()+dotv2v4*(vec4.pz()-vec2.pz());

	  vecb[3] = dotv4p24*vec2.t()
	    -dotv2p24*vec4.t()+dotv2v4*(vec4.e()-vec2.e());

	  InvEnergy2 numerator = 1./( 
				   sqr(vec1.e() +vec3.e() )
				  -sqr(vec1.px()+vec3.px())
				  -sqr(vec1.py()+vec3.py())
				  -sqr(vec1.pz()+vec3.pz())
				  );

	  vertex += numerator*(veca[3]*vecb[3]-veca[0]*vecb[0]
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
	      complex<Energy> dotv1p13 =
		+vec1.t()*(vec1.e() +2.*vec3.e() )
		-vec1.x()*(vec1.px()+2.*vec3.px())
		-vec1.y()*(vec1.py()+2.*vec3.py())
		-vec1.z()*(vec1.pz()+2.*vec3.pz());
	      complex<Energy> dotv1p14 =
		+vec1.t()*(vec1.e() +2.*vec4.e() )
		-vec1.x()*(vec1.px()+2.*vec4.px())
		-vec1.y()*(vec1.py()+2.*vec4.py())
		-vec1.z()*(vec1.pz()+2.*vec4.pz());
	      complex<Energy> dotv2p23 =
		+vec2.t()*(vec2.e() +2.*vec3.e() )
		-vec2.x()*(vec2.px()+2.*vec3.px())
		-vec2.y()*(vec2.py()+2.*vec3.py())
		-vec2.z()*(vec2.pz()+2.*vec3.pz());
	      complex<Energy> dotv2p24 =
		+vec2.t()*(vec2.e() +2.*vec4.e() )
		-vec2.x()*(vec2.px()+2.*vec4.px())
		-vec2.y()*(vec2.py()+2.*vec4.py())
		-vec2.z()*(vec2.pz()+2.*vec4.pz());
	      complex<Energy> dotv3p31 = 
		+vec3.t()*(vec3.e() +2.*vec1.e() )
		-vec3.x()*(vec3.px()+2.*vec1.px())
		-vec3.y()*(vec3.py()+2.*vec1.py())
		-vec3.z()*(vec3.pz()+2.*vec1.pz());
	      complex<Energy> dotv3p32 = 
		+vec3.t()*(vec3.e() +2.*vec2.e() )
		-vec3.x()*(vec3.px()+2.*vec2.px())
		-vec3.y()*(vec3.py()+2.*vec2.py())
		-vec3.z()*(vec3.pz()+2.*vec2.pz());
	      complex<Energy> dotv4p41 = 
		+vec4.t()*(vec4.e() +2.*vec1.e() )
		-vec4.x()*(vec4.px()+2.*vec1.px())
		-vec4.y()*(vec4.py()+2.*vec1.py())
		-vec4.z()*(vec4.pz()+2.*vec1.pz());
	      complex<Energy> dotv4p42 = 
		+vec4.t()*(vec4.e() +2.*vec2.e() )
		-vec4.x()*(vec4.px()+2.*vec2.px())
		-vec4.y()*(vec4.py()+2.*vec2.py())
		-vec4.z()*(vec4.pz()+2.*vec2.pz());
	      complex<Energy> ja[4],jb[4],jc[4],jd[4];
	      // vectors
	      ja[0] =(vec3.px()-vec1.px())*dotv1v3
		+dotv3p31*vec1.x()-dotv1p13*vec3.x();
	      jc[0] =(vec4.px()-vec1.px())*dotv1v4
		+dotv4p41*vec1.x()-dotv1p14*vec4.x();
	      jb[0] =(vec4.px()-vec2.px())*dotv2v4
		+dotv4p42*vec2.x()-dotv2p24*vec4.x();
	      jd[0] =(vec3.px()-vec2.px())*dotv2v3
		+dotv3p32*vec2.x()-dotv2p23*vec3.x();

	      ja[1] =(vec3.py()-vec1.py())*dotv1v3
		+dotv3p31*vec1.y()-dotv1p13*vec3.y();
	      jc[1] =(vec4.py()-vec1.py())*dotv1v4
		+dotv4p41*vec1.y()-dotv1p14*vec4.y();
	      jb[1] =(vec4.py()-vec2.py())*dotv2v4
		+dotv4p42*vec2.y()-dotv2p24*vec4.y();
	      jd[1] =(vec3.py()-vec2.py())*dotv2v3
		+dotv3p32*vec2.y()-dotv2p23*vec3.y();

	      ja[2] =(vec3.pz()-vec1.pz())*dotv1v3
		+dotv3p31*vec1.z()-dotv1p13*vec3.z();
	      jc[2] =(vec4.pz()-vec1.pz())*dotv1v4
		+dotv4p41*vec1.z()-dotv1p14*vec4.z();
	      jb[2] =(vec4.pz()-vec2.pz())*dotv2v4
		+dotv4p42*vec2.z()-dotv2p24*vec4.z();
	      jd[2] =(vec3.pz()-vec2.pz())*dotv2v3
		+dotv3p32*vec2.z()-dotv2p23*vec3.z();

	      ja[3] =(vec3.e()-vec1.e())*dotv1v3
		+dotv3p31*vec1.t()-dotv1p13*vec3.t();
	      jc[3] =(vec4.e()-vec1.e())*dotv1v4
		+dotv4p41*vec1.t()-dotv1p14*vec4.t();
	      jb[3] =(vec4.e()-vec2.e())*dotv2v4
		+dotv4p42*vec2.t()-dotv2p24*vec4.t();
	      jd[3] =(vec3.e()-vec2.e())*dotv2v3
		+dotv3p32*vec2.t()-dotv2p23*vec3.t();
		
	      // dot products of these vectors
	      complex<Energy2> dotjajb = ja[3]*jb[3]
		-ja[0]*jb[0]-ja[1]*jb[1]-ja[2]*jb[2];
	      complex<Energy2> dotjcjd = jc[3]*jd[3]
		-jc[0]*jd[0]-jc[1]*jd[1]-jc[2]*jd[2];
	      complex<Energy2> dotjaq = 
		+ja[3]*(vec1.e() +vec3.e() )-ja[0]*(vec1.px()+vec3.px())
		-ja[1]*(vec1.py()+vec3.py())-ja[2]*(vec1.pz()+vec3.pz());
	      complex<Energy2> dotjbq = 
		+jb[3]*(vec1.e() +vec3.e() )-jb[0]*(vec1.px()+vec3.px())
		-jb[1]*(vec1.py()+vec3.py())-jb[2]*(vec1.pz()+vec3.pz());
	      complex<Energy2> dotjck = 
		+jc[3]*(vec1.e() +vec4.e() )-jc[0]*(vec1.px()+vec4.px())
		-jc[1]*(vec1.py()+vec4.py())-jc[2]*(vec1.pz()+vec4.pz());
	      complex<Energy2> dotjdk = 
		+jd[3]*(vec1.e() +vec4.e() )-jd[0]*(vec1.px()+vec4.px())
		-jd[1]*(vec1.py()+vec4.py())-jd[2]*(vec1.pz()+vec4.pz());
	      Energy2 q2 = 
		 (vec1.e() +vec3.e() )*(vec1.e() +vec3.e() )
		-(vec1.px()+vec3.px())*(vec1.px()+vec3.px())
		-(vec1.py()+vec3.py())*(vec1.py()+vec3.py())
		-(vec1.pz()+vec3.pz())*(vec1.pz()+vec3.pz());
	      Energy2 k2 = 
		 (vec1.e() +vec4.e() )*(vec1.e() +vec4.e() )
		-(vec1.px()+vec4.px())*(vec1.px()+vec4.px())
		-(vec1.py()+vec4.py())*(vec1.py()+vec4.py())
		-(vec1.pz()+vec4.pz())*(vec1.pz()+vec4.pz());
	      // compute the term we need
	      Energy2 mass2;
	      for(int ix=0;ix<2;++ix)
		{
		  if(_inter[ix])
		    {
		      mass2 = sqr(_inter[ix]->mass());
		      if(mass2!=Energy2())
			{
			 
			  vertex += UnitRemoval::InvE2 *
			    _coup[ix]*propagator(iopt,q2,_inter[ix])*
			    (dotjajb-dotjaq*dotjbq/mass2);
			  
			  vertex += UnitRemoval::InvE2 *
			    _coup[ix]*propagator(iopt,k2,_inter[ix])*
			    (dotjcjd-dotjck*dotjdk/mass2);
			    
			}
		      else
			{
			  vertex+=UnitRemoval::InvE2 *_coup[ix]*propagator(iopt,q2,_inter[ix])*dotjajb;
			  vertex+=UnitRemoval::InvE2 *_coup[ix]*propagator(iopt,k2,_inter[ix])*dotjcjd;
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
	      complex<Energy> dotv1p12 =
		+vec1.t()*(vec1.e() +2.*vec2.e() )
		-vec1.x()*(vec1.px()+2.*vec2.px())
		-vec1.y()*(vec1.py()+2.*vec2.py())
		-vec1.z()*(vec1.pz()+2.*vec2.pz());
	      complex<Energy> dotv1p14 =
		+vec1.t()*(vec1.e() +2.*vec4.e() )
		-vec1.x()*(vec1.px()+2.*vec4.px())
		-vec1.y()*(vec1.py()+2.*vec4.py())
		-vec1.z()*(vec1.pz()+2.*vec4.pz());
	      complex<Energy> dotv3p32 =
		+vec3.t()*(vec3.e() +2.*vec2.e() )
		-vec3.x()*(vec3.px()+2.*vec2.px())
		-vec3.y()*(vec3.py()+2.*vec2.py())
		-vec3.z()*(vec3.pz()+2.*vec2.pz());
	      complex<Energy> dotv3p34 =
		+vec3.t()*(vec3.e() +2.*vec4.e() )
		-vec3.x()*(vec3.px()+2.*vec4.px())
		-vec3.y()*(vec3.py()+2.*vec4.py())
		-vec3.z()*(vec3.pz()+2.*vec4.pz());
	      complex<Energy> dotv2p21 = 
		+vec2.t()*(vec2.e() +2.*vec1.e() )
		-vec2.x()*(vec2.px()+2.*vec1.px())
		-vec2.y()*(vec2.py()+2.*vec1.py())
		-vec2.z()*(vec2.pz()+2.*vec1.pz());
	      complex<Energy> dotv2p23 = 
		+vec2.t()*(vec2.e() +2.*vec3.e() )
		-vec2.x()*(vec2.px()+2.*vec3.px())
		-vec2.y()*(vec2.py()+2.*vec3.py())
		-vec2.z()*(vec2.pz()+2.*vec3.pz());
	      complex<Energy> dotv4p41 = 
		+vec4.t()*(vec4.e() +2.*vec1.e() )
		-vec4.x()*(vec4.px()+2.*vec1.px())
		-vec4.y()*(vec4.py()+2.*vec1.py())
		-vec4.z()*(vec4.pz()+2.*vec1.pz());
	      complex<Energy> dotv4p43 = 
		+vec4.t()*(vec4.e() +2.*vec3.e() )
		-vec4.x()*(vec4.px()+2.*vec3.px())
		-vec4.y()*(vec4.py()+2.*vec3.py())
		-vec4.z()*(vec4.pz()+2.*vec3.pz());
	      complex<Energy> ja[4],jb[4],jc[4],jd[4];
	      // vectors
	      ja[0] =(vec2.px()-vec1.px())*dotv1v2
		+dotv2p21*vec1.x()-dotv1p12*vec2.x();
	      jc[0] =(vec4.px()-vec1.px())*dotv1v4
		+dotv4p41*vec1.x()-dotv1p14*vec4.x();
	      jb[0] =(vec4.px()-vec3.px())*dotv3v4
		+dotv4p43*vec3.x()-dotv3p34*vec4.x();
	      jd[0] =(vec2.px()-vec3.px())*dotv2v3
		+dotv2p23*vec3.x()-dotv3p32*vec2.x();
		
	      ja[1] =(vec2.py()-vec1.py())*dotv1v2
		+dotv2p21*vec1.y()-dotv1p12*vec2.y();
	      jc[1] =(vec4.py()-vec1.py())*dotv1v4
		+dotv4p41*vec1.y()-dotv1p14*vec4.y();
	      jb[1] =(vec4.py()-vec3.py())*dotv3v4
		+dotv4p43*vec3.y()-dotv3p34*vec4.y();
	      jd[1] =(vec2.py()-vec3.py())*dotv2v3
		+dotv2p23*vec3.y()-dotv3p32*vec2.y();
		
	      ja[2] =(vec2.pz()-vec1.pz())*dotv1v2
		+dotv2p21*vec1.z()-dotv1p12*vec2.z();
	      jc[2] =(vec4.pz()-vec1.pz())*dotv1v4
		+dotv4p41*vec1.z()-dotv1p14*vec4.z();
	      jb[2] =(vec4.pz()-vec3.pz())*dotv3v4
		+dotv4p43*vec3.z()-dotv3p34*vec4.z();
	      jd[2] =(vec2.pz()-vec3.pz())*dotv2v3
		+dotv2p23*vec3.z()-dotv3p32*vec2.z();
		
	      ja[3] =(vec2.e()-vec1.e())*dotv1v2
		+dotv2p21*vec1.t()-dotv1p12*vec2.t();
	      jc[3] =(vec4.e()-vec1.e())*dotv1v4
		+dotv4p41*vec1.t()-dotv1p14*vec4.t();
	      jb[3] =(vec4.e()-vec3.e())*dotv3v4
		+dotv4p43*vec3.t()-dotv3p34*vec4.t();
	      jd[3] =(vec2.e()-vec3.e())*dotv2v3
		+dotv2p23*vec3.t()-dotv3p32*vec2.t();
		
	      // dot products of these vectors
	      complex<Energy2> dotjajb = ja[3]*jb[3]
		-ja[0]*jb[0]-ja[1]*jb[1]-ja[2]*jb[2];
	      complex<Energy2> dotjcjd = jc[3]*jd[3]
		-jc[0]*jd[0]-jc[1]*jd[1]-jc[2]*jd[2];
	      complex<Energy2> dotjaq = 
		+ja[3]*(vec1.e() +vec2.e() )-ja[0]*(vec1.px()+vec2.px())
		-ja[1]*(vec1.py()+vec2.py())-ja[2]*(vec1.pz()+vec2.pz());
	      complex<Energy2> dotjbq = 
		+jb[3]*(vec1.e() +vec2.e() )-jb[0]*(vec1.px()+vec2.px())
		-jb[1]*(vec1.py()+vec2.py())-jb[2]*(vec1.pz()+vec2.pz());
	      complex<Energy2> dotjck = 
		+jc[3]*(vec1.e() +vec4.e() )-jc[0]*(vec1.px()+vec4.px())
		-jc[1]*(vec1.py()+vec4.py())-jc[2]*(vec1.pz()+vec4.pz());
	      complex<Energy2> dotjdk = 
		+jd[3]*(vec1.e() +vec4.e() )-jd[0]*(vec1.px()+vec4.px())
		-jd[1]*(vec1.py()+vec4.py())-jd[2]*(vec1.pz()+vec4.pz());
	      Energy2 q2 = 
		 (vec1.e() +vec2.e() )*(vec1.e() +vec2.e() )
		-(vec1.px()+vec2.px())*(vec1.px()+vec2.px())
		-(vec1.py()+vec2.py())*(vec1.py()+vec2.py())
		-(vec1.pz()+vec2.pz())*(vec1.pz()+vec2.pz());
	      Energy2 k2 = 
		 (vec1.e() +vec4.e() )*(vec1.e() +vec4.e() )
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
		      if(mass2!=Energy2())
			{
			  vertex+=UnitRemoval::InvE2 *
			    _coup[ix]*propagator(iopt,q2,_inter[ix])*
			    (dotjajb-dotjaq*dotjbq/mass2);
			  vertex+=UnitRemoval::InvE2 *
			    _coup[ix]*propagator(iopt,k2,_inter[ix])*
			    (dotjcjd-dotjck*dotjdk/mass2);
			}
		      else
			{
			  vertex+=UnitRemoval::InvE2 *_coup[ix]*propagator(iopt,q2,_inter[ix])*dotjajb;
			  vertex+=UnitRemoval::InvE2 *_coup[ix]*propagator(iopt,k2,_inter[ix])*dotjcjd;
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
	      complex<Energy> dotv1p12 =
		+vec1.t()*(vec1.e() +2.*vec2.e() )
		-vec1.x()*(vec1.px()+2.*vec2.px())
		-vec1.y()*(vec1.py()+2.*vec2.py())
		-vec1.z()*(vec1.pz()+2.*vec2.pz());
	      complex<Energy> dotv1p13 =
		+vec1.t()*(vec1.e() +2.*vec3.e() )
		-vec1.x()*(vec1.px()+2.*vec3.px())
		-vec1.y()*(vec1.py()+2.*vec3.py())
		-vec1.z()*(vec1.pz()+2.*vec3.pz());
	      complex<Energy> dotv2p24 = 
		+vec2.t()*(vec2.e() +2.*vec4.e() )
		-vec2.x()*(vec2.px()+2.*vec4.px())
		-vec2.y()*(vec2.py()+2.*vec4.py())
		-vec2.z()*(vec2.pz()+2.*vec4.pz());
	      complex<Energy> dotv2p21 = 
		+vec2.t()*(vec2.e() +2.*vec1.e() )
		-vec2.x()*(vec2.px()+2.*vec1.px())
		-vec2.y()*(vec2.py()+2.*vec1.py())
		-vec2.z()*(vec2.pz()+2.*vec1.pz());
	      complex<Energy> dotv3p31 = 
		+vec3.t()*(vec3.e() +2.*vec1.e() )
		-vec3.x()*(vec3.px()+2.*vec1.px())
		-vec3.y()*(vec3.py()+2.*vec1.py())
		-vec3.z()*(vec3.pz()+2.*vec1.pz());
	      complex<Energy> dotv3p34 =
		+vec3.t()*(vec3.e() +2.*vec4.e() )
		-vec3.x()*(vec3.px()+2.*vec4.px())
		-vec3.y()*(vec3.py()+2.*vec4.py())
		-vec3.z()*(vec3.pz()+2.*vec4.pz());
	      complex<Energy> dotv4p43 = 
		+vec4.t()*(vec4.e() +2.*vec3.e() )
		-vec4.x()*(vec4.px()+2.*vec3.px())
		-vec4.y()*(vec4.py()+2.*vec3.py())
		-vec4.z()*(vec4.pz()+2.*vec3.pz());
	      complex<Energy> dotv4p42 =
		+vec4.t()*(vec4.e() +2.*vec2.e() )
		-vec4.x()*(vec4.px()+2.*vec2.px())
		-vec4.y()*(vec4.py()+2.*vec2.py())
		-vec4.z()*(vec4.pz()+2.*vec2.pz());
	      complex<Energy> ja[4],jb[4],jc[4],jd[4];
	      // vectors
	      ja[0] =(vec2.px()-vec1.px())*dotv1v2
		+dotv2p21*vec1.x()-dotv1p12*vec2.x();
	      jc[0] =(vec3.px()-vec1.px())*dotv1v3
		+dotv3p31*vec1.x()-dotv1p13*vec3.x();
	      jb[0] =(vec3.px()-vec4.px())*dotv3v4
		+dotv3p34*vec4.x()-dotv4p43*vec3.x();
	      jd[0] =(vec2.px()-vec4.px())*dotv2v4
		+dotv2p24*vec4.x()-dotv4p42*vec2.x();
		
	      ja[1] =(vec2.py()-vec1.py())*dotv1v2
		+dotv2p21*vec1.y()-dotv1p12*vec2.y();
	      jc[1] =(vec3.py()-vec1.py())*dotv1v3
		+dotv3p31*vec1.y()-dotv1p13*vec3.y();
	      jb[1] =(vec3.py()-vec4.py())*dotv3v4
		+dotv3p34*vec4.y()-dotv4p43*vec3.y();
	      jd[1] =(vec2.py()-vec4.py())*dotv2v4
		+dotv2p24*vec4.y()-dotv4p42*vec2.y();
		
	      ja[2] =(vec2.pz()-vec1.pz())*dotv1v2
		+dotv2p21*vec1.z()-dotv1p12*vec2.z();
	      jc[2] =(vec3.pz()-vec1.pz())*dotv1v3
		+dotv3p31*vec1.z()-dotv1p13*vec3.z();
	      jb[2] =(vec3.pz()-vec4.pz())*dotv3v4
		+dotv3p34*vec4.z()-dotv4p43*vec3.z();
	      jd[2] =(vec2.pz()-vec4.pz())*dotv2v4
		+dotv2p24*vec4.z()-dotv4p42*vec2.z();
		
	      ja[3] =(vec2.e()-vec1.e())*dotv1v2
		+dotv2p21*vec1.t()-dotv1p12*vec2.t();
	      jc[3] =(vec3.e()-vec1.e())*dotv1v3
		+dotv3p31*vec1.t()-dotv1p13*vec3.t();
	      jb[3] =(vec3.e()-vec4.e())*dotv3v4
		+dotv3p34*vec4.t()-dotv4p43*vec3.t();
	      jd[3] =(vec2.e()-vec4.e())*dotv2v4
		+dotv2p24*vec4.t()-dotv4p42*vec2.t();
		
	      // dot products of these vectors
	      complex<Energy2> dotjajb = ja[3]*jb[3]
		-ja[0]*jb[0]-ja[1]*jb[1]-ja[2]*jb[2];
	      complex<Energy2> dotjcjd = jc[3]*jd[3]
		-jc[0]*jd[0]-jc[1]*jd[1]-jc[2]*jd[2];
	      complex<Energy2> dotjaq = 
		+ja[3]*(vec1.e() +vec2.e() )-ja[0]*(vec1.px()+vec2.px())
		-ja[1]*(vec1.py()+vec2.py())-ja[2]*(vec1.pz()+vec2.pz());
	      complex<Energy2> dotjbq = 
		+jb[3]*(vec1.e() +vec2.e() )-jb[0]*(vec1.px()+vec2.px())
		-jb[1]*(vec1.py()+vec2.py())-jb[2]*(vec1.pz()+vec2.pz());
	      complex<Energy2> dotjck = 
		+jc[3]*(vec1.e() +vec3.e() )-jc[0]*(vec1.px()+vec3.px())
		-jc[1]*(vec1.py()+vec3.py())-jc[2]*(vec1.pz()+vec3.pz());
	      complex<Energy2> dotjdk = 
		+jd[3]*(vec1.e() +vec3.e() )-jd[0]*(vec1.px()+vec3.px())
		-jd[1]*(vec1.py()+vec3.py())-jd[2]*(vec1.pz()+vec3.pz());
	      Energy2 q2 = 
		 (vec1.e() +vec2.e() )*(vec1.e() +vec2.e() )
		-(vec1.px()+vec2.px())*(vec1.px()+vec2.px())
		-(vec1.py()+vec2.py())*(vec1.py()+vec2.py())
		-(vec1.pz()+vec2.pz())*(vec1.pz()+vec2.pz());
	      Energy2 k2 = 
		 (vec1.e() +vec3.e() )*(vec1.e() +vec3.e() )
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
		      if(mass2!=Energy2())
			{
			  vertex+=UnitRemoval::InvE2 *
			    _coup[ix]*propagator(iopt,q2,_inter[ix])*
			    (dotjajb-dotjaq*dotjbq/mass2);
			  vertex+=UnitRemoval::InvE2 *
			    _coup[ix]*propagator(iopt,k2,_inter[ix])*
			    (dotjcjd-dotjck*dotjdk/mass2);
			}
		      else
			{
			  vertex+=UnitRemoval::InvE2 *_coup[ix]*propagator(iopt,q2,_inter[ix])*dotjajb;
			  vertex+=UnitRemoval::InvE2 *_coup[ix]*propagator(iopt,k2,_inter[ix])*dotjcjd;
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

