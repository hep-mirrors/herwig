// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NMSSMFFHVertex class.
//

#include "NMSSMFFHVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "NMSSM.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

NMSSMFFHVertex::NMSSMFFHVertex() : _mw(0.*MeV), _sinb(0.), _cosb(0.), 
				   _tanb(0.), _idlast(make_pair(0,0)),
				   _q2last(0.*MeV2), 
				   _masslast(make_pair(0.*MeV,0*MeV)),
				   _couplast(0.) {
  // the quarks and neutral higgs
  int in[5]={25,35,45,36,46};
  for(unsigned int iy=0;iy<5;++iy)
    for(int ix=1;ix<7;++ix)
      addToList( -ix, ix, in[iy] );

  // leptons and neutral higgs
  for(unsigned int iy=0;iy<5;++iy)
    for(int ix=11;ix<17;ix+=2)
      addToList( -ix, ix, in[iy] );

  // the quarks  and the charged higgs
  //H-
  for(int ix=0;ix<3;++ix) 
    addToList(2*ix+2, -2*ix-1, -37);

  //H+
  for(int ix=0;ix<3;++ix)
    addToList(-(2*ix+2), 2*ix+1, 37);

  // the leptons and the charged higgs
  //H-
  for(int ix=0;ix<3;++ix)
    addToList( 2*ix+12, -2*ix-11, -37 );

  //H+
  for(int ix=0;ix<3;++ix)
    addToList( -(2*ix+12), 2*ix+11, 37 );
}

void NMSSMFFHVertex::persistentOutput(PersistentOStream & os) const {
  os << _mixS << _mixP << ounit(_mw,GeV)
     << _sinb << _cosb << _tanb << _sw << _theSM;
}

void NMSSMFFHVertex::persistentInput(PersistentIStream & is, int) {
  is >> _mixS >> _mixP >> iunit(_mw,GeV)
     >> _sinb >> _cosb >> _tanb >> _sw >> _theSM;
}

void NMSSMFFHVertex::doinit() {
  // cast to NMSSM model
  tcNMSSMPtr model=dynamic_ptr_cast<tcNMSSMPtr>(generator()->standardModel());
  if(!model) 
    throw InitException() << "Must have the NMSSM Model in NMSSMFFHVertex::doinit()"
			  << Exception::runerror;
  _theSM = model;
  // sin theta_W
  double sw2=sin2ThetaW();
  _sw = sqrt(sw2);
  // get the mixing matrices
  _mixS=model->CPevenHiggsMix();
  if(!_mixS) throw InitException() << "Mixing matrix for CP-even neutral Higgs"
				   << " bosons is not set in NMSSMFFHVertex::doinit()" 
				   << Exception::runerror;
  _mixP=model->CPoddHiggsMix();
  if(!_mixP) throw InitException() << "Mixing matrix for CP-odd neutral Higgs"
				   << " bosons is not set in NMSSMFFHVertex::doinit()" 
				   << Exception::runerror;
  // Mass of the W boson
  _mw=getParticleData(ParticleID::Wplus)->mass();
  // sin and cos beta
  _tanb = model->tanBeta();
  double beta = atan(_tanb);
  _sinb=sin(beta);
  _cosb=cos(beta);
  // order in couplings
  orderInGem(1);
  orderInGs(0);
  // base class
  FFSVertex::doinit();
}

ClassDescription<NMSSMFFHVertex> NMSSMFFHVertex::initNMSSMFFHVertex;
// Definition of the static class description member.

void NMSSMFFHVertex::Init() {

  static ClassDocumentation<NMSSMFFHVertex> documentation
    ("The NMSSMFFHVertex class implements the vertex for the couplings"
     " of the Higgs bosons of the NMSSM to Standard Model fermions");

}
//calulate the couplings
void NMSSMFFHVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b, tcPDPtr c) {
  int ihiggs=c->id();
  int id(abs(a->id()));
  Complex output(1.);
  double _pi= Constants::pi ;
  Energy _mh = getParticleData(ihiggs)->mass();
  double _alphas =_theSM->alphaS(q2);
  // neutral Higgs
  if(ihiggs==25||ihiggs==35||ihiggs==45||ihiggs==36||ihiggs==46) {
    if(_idlast.first!=id||q2!=_q2last) {
      _idlast.first=id;
      _masslast.first = _theSM->mass(q2,a);
    }
	
	Energy _mp = getParticleData(id)->mass();
    Energy _mtp = getParticleData(6)->mass();

    double temp = sqr(_masslast.first/_mh);
	double temp2 = log(temp);
	double ratio = sqr(_mp/_mh);
	double beta= sqrt(1. - 4*ratio);
	double higtop = sqr(_mh/_mtp);
    double rat = 2.*_mp/_mh;
	Complex _cu, _cd, ratcoup, HQCD,QCDH, QCD0, HQCDM, TQCDH;

	QCD0 = (1.0 + sqr(beta))*(4.0*SPfunction((1.0 - beta)/(1.0 + beta))
     	 + 2.0*SPfunction((beta - 1.0)/(beta + 1.0))
     	 - 3.0*log((1.0 + beta)/(1.0 - beta))*log(2.0/(1.0 + beta))
     	 - 2.*log((1.0 + beta)/(1.0 - beta))*log(beta))
     	 - 3.0*beta*log(4.0/(1.0 - sqr(beta))) - 4.0*beta*log(beta);

    // CP-even 
    if(ihiggs==25||ihiggs==35||ihiggs==45) {
      int iloc = (ihiggs-25)/10;
     _cu = (*_mixS)(iloc,1)/_sinb;
	 _cd = (*_mixS)(iloc,0)/_cosb;
	 if ( id == 4){  //c
	  ratcoup = 1.0;
	 }
	 else if (id == 6){  //t
	 ratcoup = 0.0;
	 }
	 else {  
	 ratcoup = _cu/_cd;
	 }
	 
	HQCDM= QCD0/beta + (3.0 + 34.0*sqr(beta) 
		- 13.0*sqr(beta)*sqr(beta))/(16.0*beta*sqr(beta))
     	 *log((1.0 + beta)/(1.0 - beta)) 
		 + 3.0/(8.0*sqr(beta))*(7.0*sqr(beta) - 1.0);
	TQCDH= 1.0 + 4.0/3.0*HQCDM*_alphas/_pi;
	HQCD= 5.67*(_alphas/_pi)
     	 + (29.14 + ratcoup*(1.57 - 2.0/3.0*log(higtop)
     	 + sqr(temp2)/9.))*sqr(_alphas/_pi)
     	 + (164.14 - 25.77*5.0 + 0.259*25.0)*(_alphas/_pi)*sqr(_alphas/_pi);
   QCDH = (1.0 + HQCD);
 
   output *= (id%2==0) ? (*_mixS)(iloc,1)/_sinb : (*_mixS)(iloc,0)/_cosb;

   if ( id == 6 && _mh <= 2.*_masslast.first) {
     output *=   rat*sqrt(TQCDH)*_mp/_mw;
   }
   else if  ( real(sqr(output*QCDH)) < 0.0){
       output *=   rat*sqrt(TQCDH)*_mp/_mw;
   }
   else if ( id < 7) {
   output *= sqrt(sqr(rat)*TQCDH*sqr(_mp)/sqr(_mw) + (1.0 - sqr(rat))*QCDH*sqr(_masslast.first)/sqr(_mw));
   }	 
   else if ( id > 7 ){
   output *= _masslast.first/_mw;
   }
   
   left(1.); right(1.);
} 
    // CP-odd
    else {
      int iloc = (ihiggs-36)/10;
	  _cu = (*_mixP)(iloc,1)/_sinb;
	 _cd  = (*_mixP)(iloc,0)/_cosb;
	 
	if ( id == 4){  //c
	  ratcoup = 1.0;
	 }
	 else if (id == 6){  //t
	 ratcoup = 0.00000;
	 }
	else{  
	 ratcoup = _cu/_cd;
	 }
	  
	 HQCD= 5.67*(_alphas/_pi)
     	 + (29.14 + ratcoup*(3.83 - log(higtop)
		 + sqr(temp2)/6.))*sqr(_alphas/_pi)
     	 + (164.14 - 25.77*5.0 + 0.259*25)*(_alphas/_pi)*sqr(_alphas/_pi);
	 QCDH = 1.0 + HQCD;
	 HQCDM= QCD0/beta + (19.0 + 2.0*sqr(beta)
		        + 3.0*sqr(beta)*sqr(beta))/16.0/beta
               *log((1.0 + beta)/(1.0 - beta)) + 3.0/8.0*(7.0 - sqr(beta));
	 TQCDH= 1.0 + 4.0/3.0*HQCDM*_alphas/_pi;
	 
      output *= (id%2==0) ? (*_mixP)(iloc,1)/_sinb : (*_mixP)(iloc,0)/_cosb;
	  
    if ( id == 6 && _mh <= 2.*_masslast.first) {
     output *=   rat*sqrt(TQCDH)*_mp/_mw;
     }
     else if ( id < 7) {
	 if ( real(sqr(output*QCDH)) < 0.0){
       output *=   rat*sqrt(TQCDH)*_mp/_mw;
	 }
	else{
	output *=sqrt(sqr(rat)*TQCDH*sqr(_mp)/sqr(_mw) + (1.0 - sqr(rat))*QCDH*sqr(_masslast.first)/sqr(_mw));
	} }
	else if (id > 7){
	output*=_masslast.first/_mw;
	}
	  
      left(1.); right(-1.);
      output *= Complex(0.,-1.);
	  
    }
  }
  // Charged higgs
  else if(abs(ihiggs)==37) {
    output*=-sqrt(2.);
    int id2=abs(b->id());
    if(id2<id) swap(id,id2);
    if(_idlast.first!=id||_idlast.second!=id2||q2!=_q2last) {
      _idlast.first =id ;
      _idlast.second=id2;
      _masslast.first  = _theSM->mass(q2,a);
      _masslast.second = _theSM->mass(q2,b);
    }
    double rgt=_masslast.first *_tanb/_mw;
    double lft =_masslast.second/_tanb/_mw;
    if(ihiggs<0) swap(lft,rgt);
    right(rgt);
    left (lft);
  }
  else {
    throw Exception() << "Unknown Higgs boson, PDG code = " << ihiggs 
		      << "in NMSSMFFHVertex::setCoupling()"
		      << Exception::runerror;
  }
  // prefactor
  if(q2!=_q2last) {
    _couplast = 0.5*weakCoupling(q2);
    _q2last=q2;
  }
  norm(-_couplast*output);
}


  //Below taken fron NMHDecay
  //*   Spence function and auxiliary functions as in HDECAY
double NMSSMFFHVertex::SPfunction( Complex X) {

//*  REAL dilogarithm (Spence-function)


//	DOUBLE COMPLEX FUNCTION LI2(X)


    Complex Y;
	 double _pi= Constants::pi;
    double ZETA2 = sqr(_pi)/6.0;
	double ZERO=1.0000000000000000;
	double XR = real(X);
	Complex XI = imag(X);
	Complex R2=XR*XR+XI*XI;
	Complex LI2=0;
	if (real(R2) <= ZERO) {
	  LI2=real(X);
     }
	   Complex RR=XR/R2;
	if(real(R2) == 1.0 && XI==0.0){
	if(XR == 1.0){
	       LI2=ZETA2;
		   }
	else  {
	    LI2 = -ZETA2/2.0;
	   }
	   }
	else if (real(R2) > 1.0 && real(RR) > 0.5){
	   Y=(real(X) - 1.0)/real(X);
	  LI2=CLI2function(Y) + ZETA2 - log(real(X))*log(1.0 - real(X))+0.5*sqr(log(real(X)));
	  }
	else if (real(R2) > 1.0 && real(RR) <= 0.5){
	  Y = 1.0/real(X);
	  LI2 = -CLI2function(Y) - ZETA2 - 0.5*sqr(log(-real(X)));
	  }
	else if (real(R2) <= 1.0 && XR > 0.5){
	  Y = 1.0 - real(X);
	  LI2 = -CLI2function(Y) + ZETA2 - log(real(X))*log(1.0 - real(X));
	 }
	else if (real(R2) <= 1.0 && XR <= 0.5){
	  Y = real(X);
	  LI2 = CLI2function(Y);
	  }
   
 //*  REAL dilogarithm (Spence-function)
    double SP = real(LI2);
    return  SP;
  }		


//*  COMPLEX dilogarithm (Spence-function)
Complex NMSSMFFHVertex::CLI2function( Complex Y) {
double B[18], B2[18];
  	B[0]=-1.0/2.0;
  	B[1]=1.0/6.0;
	B[2]=0.0;
 	B[3]=-1.0/30.0;
	B[4]=0.0;
	B[5]=1.0/42.0;
	B[6]=0.0;
  	B[7]=-1.0/30.0;
	B[8]=0.0;
   	B[9]=5.0/66.0;
 	B[10]=0.0;
    B[11]=-691.0/2730.0;
    B[12]=0.0;
    B[13]=7.0/6.0;
 	B[14]=0.0;
   	B[15]=-3617.0/510.0;
	B[16]=0.0;
  	B[17]=43867.0/798.0;
//	DOUBLE COMPLEX FUNCTION CLI2(X)
//  Taylor-expansion for complex dilogarithm (Spence-function)
    Complex	Z=-log (1.0 - real(Y));
    Complex	CLI2 = B[17]/factorial(18);
for (int p = 16 ; p = 0; --p){
          B2[p]=B[p]/factorial(p + 1);
		
		
		  CLI2=Z*CLI2+B2[p];
			 }
	CLI2 = sqr(Z)*CLI2 + Z;
return CLI2;

}


int NMSSMFFHVertex::factorial(int num){ 
 int result=1;
 for (int i=1; i<=num; ++i)
    result=result*=i;
 return result;
}

double NMSSMFFHVertex::BIJ(double X, double Y){ 
	double LAMB= sqrt(sqr((1.0- X - Y)) - 4.0*X*Y);	
	double		XIX = 2.0*X/(1.0 - X - Y + LAMB);
	double		XIY = 2.0*Y/(1.0 - X - Y + LAMB);
	
		 
	double	BIJ= (1.0 - Y - X)/LAMB
     	 * (4.0*SPfunction(XIX*XIY)
     	 - 2.0*SPfunction(-XIX) - 2.0*SPfunction(-XIY)
     	 + 2.0*log(XIX*XIY)*log(1.0 - XIX*XIY)
     	 - log(XIX)*log(1.0 + XIX)
     	 - log(XIY)*log(1.0 + XIY))
     	 - 4.0*(log(1.0 - XIX*XIY)
     	 + XIX*XIY/(1.0 - XIX*XIY)*log(XIX*XIY))
     	 + (LAMB + X -Y)/LAMB*(log(1.0 + XIX)
     	 - XIX/(1.0 + XIX)*log(abs(XIX)))
     	 + (LAMB- X + Y)/LAMB*(log(1.0 + XIY)
     	 - XIY/(1.0 + XIY)*log(abs(XIY)));	
		 
		 return BIJ;
		 
		 }
