// -*- C++ -*-
//
// DifractivePDF.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#include "PomeronPDF.h"
#include <ThePEG/PDT/ParticleData.h>
#include <ThePEG/PDT/EnumParticles.h>
#include <ThePEG/Persistency/PersistentOStream.h>
#include <ThePEG/Persistency/PersistentIStream.h>
#include <ThePEG/Repository/EventGenerator.h>
#include <ThePEG/Interface/ClassDocumentation.h>
#include <ThePEG/Interface/Parameter.h>
#include <ThePEG/Interface/Switch.h>
#include <istream>
#include <iostream>
#include <string>

using namespace std;
using namespace ThePEG;
using namespace Herwig;

const int  PomeronPDF::nPDFFlavour_ = 3;

PomeronPDF::PomeronPDF() :
  pdfTable_(nPDFFlavour_, vector<vector<double> >() ),
  lxGrid_(nPDFFlavour_, vector<double>() ),
  lqqGrid_ (nPDFFlavour_, vector<double>() ),
  fileName_(3),
  nxPoints_(100),
  nqPoints_(88),
  PDFFit_(0),
  boundary_(0)
{}

bool PomeronPDF::canHandleParticle(tcPDPtr particle) const {
  return ( abs(particle->id()) == ParticleID::pomeron );
}

cPDVector PomeronPDF::partons(tcPDPtr p) const {
  // Return the parton types which are described by these parton
  // densities.
  cPDVector ret;
  if ( canHandleParticle(p) ) {
    ret.push_back(getParticleData( ParticleID::g));
    if(PDFFit_==0) {
      ret.push_back(getParticleData( ParticleID::c));
      ret.push_back(getParticleData( ParticleID::cbar));
    }
    ret.push_back(getParticleData( ParticleID::d));
    ret.push_back(getParticleData( ParticleID::dbar));
    ret.push_back(getParticleData( ParticleID::u));
    ret.push_back(getParticleData( ParticleID::ubar));
    ret.push_back(getParticleData( ParticleID::s));
    ret.push_back(getParticleData( ParticleID::sbar));
  }
  return ret;
}

#ifdef NDEBUG
double PomeronPDF::xfx(tcPDPtr, tcPDPtr parton, Energy2 qq,
		       double x, double, Energy2) const {
#else
double PomeronPDF::xfx(tcPDPtr particle, tcPDPtr parton, Energy2 qq,
		       double x, double, Energy2) const {
#endif
  // assert particle is pomeron
  assert(particle->id()==ParticleID::pomeron);
  switch(parton->id()) {
  case ParticleID::g:
    return getPDFValue(gluon  , x, qq);
  case ParticleID::u: case ParticleID::ubar:
  case ParticleID::d: case ParticleID::dbar:
  case ParticleID::s: case ParticleID::sbar: 
    return PDFFit_==0 ? getPDFValue(singlet, x, qq)/6. : getPDFValue(singlet, x, qq);
  case ParticleID::c: case ParticleID::cbar:
    return PDFFit_==0 ? getPDFValue(charm  , x, qq)*9./8. : 0.; 
  case ParticleID::b: case ParticleID::bbar:
    return 0;
  default:
    return 0.;
  }
}

#ifdef NDEBUG
double PomeronPDF::xfvx(tcPDPtr, tcPDPtr parton, Energy2 qq,
			    double x, double, Energy2) const {
#else
double PomeronPDF::xfvx(tcPDPtr particle, tcPDPtr parton, Energy2 qq,
			    double x, double, Energy2) const {
#endif
  // assert particle is pomeron
  assert(particle->id()==ParticleID::pomeron);
  // valence parton is just gluon
  if(parton->id() ==  ParticleID::g){ 
    return getPDFValue(gluon  , x, qq);
  }
  else
    return 0.;
}

double PomeronPDF::getPDFValue(PDFFlavour flPDF, double x, Energy2 qq) const {
  // compute log(x) and log(q2)
  // Fit 2007 contains l= logx. Fit 2006 contains l = logl/log10  
  double l, lmax, lmin, lqqmin, lqqmax; 
  if(!PDFFit_) {
    l = log(x);
    lmin = lxGrid_.at(flPDF).at(0);
    lmax = lxGrid_.at(flPDF).at(nxPoints_- 1);
    lqqmin = lqqGrid_.at(flPDF).at(0);
    lqqmax = lqqGrid_.at(flPDF).at(nqPoints_ -1);
  } 
  else {
    l = log(x)/log(10);
    lmin = lxGrid_.at(flPDF).at(0);
    lmax = lxGrid_.at(flPDF).at(nxPoints_- 1);
    lqqmin = lqqGrid_.at(flPDF).at(0);
    lqqmax = lqqGrid_.at(flPDF).at(nqPoints_ -1);
  }
  double lqq = log(qq/GeV2);
  // find the x location in the table
  
  if(l <= lmin ) {
    if( lqq <= lqqmin) {
      
      double pdf00 = pdfTable_.at(flPDF).at(0).at(0);
      double pdf10 = pdfTable_.at(flPDF).at(1).at(0);
      
      if(boundary_==0){
	
	double lqq1 = lqqGrid_.at(flPDF).at(1);
	double lqq0 = lqqGrid_.at(flPDF).at(0);
	
	double dq  = lqq1 - lqq0;
	double qdq = (lqq - lqq0)/dq;
	
	return pdf00 + qdq*(pdf10 - pdf00);
      }
      else
	return pdf00;
    }
    else if( lqq >= lqqmax) {
      
      double pdf10  = pdfTable_.at(flPDF).at(nqPoints_-1).at(0);
      double pdf00 = pdfTable_.at(flPDF).at(nqPoints_-2).at(0);
      
      if(boundary_==0){
	
	double lqq0 = lqqGrid_.at(flPDF).at(nqPoints_-2);
	double lqq1 = lqqGrid_.at(flPDF).at(nqPoints_-1);
	double dq  = lqq1 - lqq0;
	double qdq = (lqq - lqq1)/dq;
	
	return pdf10 + qdq*(pdf10 - pdf00);
	
      }
      else return pdf10;
    } // (x< xmin)
    
    unsigned int qIndx = 0;
    for(; qIndx < lqqGrid_[flPDF].size(); ++qIndx )
      if(lqq < lqqGrid_[flPDF][qIndx] ) break;
    
    assert(qIndx > 0);
    
    double pdf10 = pdfTable_.at(flPDF).at(qIndx).at(0);
    double pdf00 = pdfTable_.at(flPDF).at(qIndx-1).at(0);
    
    double lqq0 = lqqGrid_.at(flPDF).at(qIndx-1);
    double lqq1 = lqqGrid_.at(flPDF).at(qIndx);
    double dq  =  lqq1 - lqq0;
    double qdq = (lqq - lqq0)/dq;
    
    return pdf00 + qdq*(pdf10 - pdf00);
  } 
  else if (l >= lmax ) {
    if( lqq <= lqqmin){
      
      double pdf01 = pdfTable_.at(flPDF).at(0).at(nxPoints_-1);
      double pdf11 = pdfTable_.at(flPDF).at(1).at(nxPoints_-1);
      
      if(boundary_==0){
	
	double lqq1 = lqqGrid_.at(flPDF).at(1);
	double lqq0 = lqqGrid_.at(flPDF).at(0);
	
	double dq  = lqq1 - lqq0;
	double qdq = (lqq - lqq0)/dq;
	
	return pdf01 + qdq*(pdf11 - pdf01);
	
      }
      else return pdf01;
    }
    else if( lqq >= lqqmax){
      
      double pdf11  = pdfTable_.at(flPDF).at(nqPoints_-1).at(nxPoints_-1);
      double pdf01 = pdfTable_.at(flPDF).at(nqPoints_-2).at(nxPoints_-1);
      
      if(boundary_==0){
	
	double lqq1  = lqqGrid_.at(flPDF).at(nqPoints_-1);
	double lqq0 = lqqGrid_.at(flPDF).at(nqPoints_-2);
	
	double dq  =  lqq1 - lqq0;
	double qdq = ( lqq - lqq1 )/dq;
	
	return pdf11 + qdq*( pdf11 - pdf01);
	
      } else return pdf11;
    }
    
    unsigned int qIndx = 0;
    for(; qIndx < lqqGrid_.at(flPDF).size(); ++qIndx )
      if(lqq < lqqGrid_.at(flPDF).at(qIndx) ) break;
    
    assert(qIndx > 0);
    
    double pdf01 = pdfTable_.at(flPDF).at(qIndx-1).at(nxPoints_-1);
    double pdf11 = pdfTable_.at(flPDF).at(qIndx  ).at(nxPoints_-1); 
    
    double lqq0 = lqqGrid_.at(flPDF).at(qIndx-1);
    double lqq1 = lqqGrid_.at(flPDF).at(qIndx);
    
    double dq  =  lqq1 - lqq0;
    double qdq = (lqq - lqq0)/dq;
    
    return  pdf01 + qdq*(pdf11 - pdf01);
  } 
  else { 
    unsigned int xIndx = 0;
    for(; xIndx < lxGrid_.at(flPDF).size(); ++xIndx ) 
      if (l < lxGrid_.at(flPDF).at(xIndx) ) break;
    
    assert(xIndx > 0);
    
    double lx0 = lxGrid_.at(flPDF).at(xIndx-1);
    double lx1 = lxGrid_.at(flPDF).at(xIndx);
    
    double dl  = lx1 - lx0;
    double ldl = (l - lx0)/dl;
    
    if( lqq <= lqqmin) {
      double pdf00 = pdfTable_.at(flPDF).at(0).at(xIndx-1);
      double pdf10 = pdfTable_.at(flPDF).at(1).at(xIndx-1);
      double pdf01 = pdfTable_.at(flPDF).at(0).at(xIndx); 
      double pdf11 = pdfTable_.at(flPDF).at(1).at(xIndx);
      
      if(boundary_==0){ 
	
	double lqq0 = lqqGrid_.at(flPDF).at(0);
	double lqq1 = lqqGrid_.at(flPDF).at(1);
	double dq   = lqq1 - lqq0; 
	double qdq  = (lqq - lqq0)/dq;
	
	return  pdf00 + ldl*(pdf01 - pdf00) + qdq*(pdf10 - pdf00)
	  + qdq*ldl*(pdf11 - pdf01 + pdf00 - pdf10);
	
      }
      else return pdf00 + ldl*(pdf01 - pdf00);
    } 
    else if(lqq >= lqqmax){
      double pdf00 = pdfTable_.at(flPDF).at(nqPoints_-2).at(xIndx-1);
      double pdf10 = pdfTable_.at(flPDF).at(nqPoints_-1).at(xIndx-1);
      double pdf01 = pdfTable_.at(flPDF).at(nqPoints_-2).at(xIndx);
      double pdf11 = pdfTable_.at(flPDF).at(nqPoints_-1).at(xIndx);
      
      if(boundary_==0){ 
	
	double lqq0 = lqqGrid_.at(flPDF).at(nqPoints_-2);
	double lqq1 = lqqGrid_.at(flPDF).at(nqPoints_-1);
	
	double dq  = lqq1 - lqq0 ;
	double qdq = ( lqq - lqq0)/dq;
	
	return pdf10 + ldl*(pdf11 - pdf10) + qdq*(pdf10 - pdf00)
	  + qdq*ldl*(pdf11 - pdf01 + pdf00 - pdf10);
	
      }
      else return pdf10 + ldl*(pdf11 - pdf10);
    } 
    else {
      unsigned int qIndx = 0;
      for(; qIndx < lqqGrid_.at(flPDF).size(); ++qIndx )
	if(lqq < lqqGrid_.at(flPDF).at(qIndx) ) break;
      
      double lqq0 = lqqGrid_.at(flPDF).at(qIndx-1);
      double lqq1 = lqqGrid_.at(flPDF).at(qIndx);
      
      double dq  = lqq1 - lqq0; 
      double qdq = (lqq - lqq0)/dq;
      
      double pdf00 = pdfTable_.at(flPDF).at(qIndx-1).at(xIndx-1);
      double pdf01 = pdfTable_.at(flPDF).at(qIndx-1).at(xIndx);
      double pdf10 = pdfTable_.at(flPDF).at(qIndx).at(xIndx-1);
      double pdf11 = pdfTable_.at(flPDF).at(qIndx).at(xIndx);
      
      return pdf00 + ldl*(pdf01 - pdf00) + qdq*(pdf10 - pdf00) 
	+ ldl*qdq*(pdf11 - pdf01 + pdf00 - pdf10);
    }
  }
  throw Exception() << "Problem in PomeronPDF::getPDFValue() "
		    << Exception::runerror;
  return 0; 
}

void PomeronPDF::loadTables() const {
  double dtemp; // helper variable for reading
  ifstream datafile;
  for(int pdfIndx = 0; pdfIndx < nPDFFlavour_; ++pdfIndx ) {
    // Set correct size for vectors
    pdfTable_.at(pdfIndx).resize(nqPoints_);
    lqqGrid_.at(pdfIndx).resize(nqPoints_);
    lxGrid_.at(pdfIndx).resize(nxPoints_);
    string name = rootName_ + fileName_.at(pdfIndx);
    datafile.open(name.c_str(), ios::in );
    if(!datafile.is_open()) 
      throw Exception() << "Could not open file '" << name
			<< "' in PomeronPDF::loadTables()"
			<< Exception::runerror;
    // read lx points
    for (unsigned int xIndx = 0; xIndx < lxGrid_ .at(pdfIndx).size(); ++xIndx) 
      datafile >> lxGrid_.at(pdfIndx).at(xIndx);
    // read qq points                                                   
    for (unsigned int qIndx = 0; qIndx < lqqGrid_.at(pdfIndx).size(); ++qIndx) {
      datafile >> dtemp;
      lqqGrid_.at(pdfIndx).at(qIndx) = log(dtemp);
    }
    // read xf(lx, qq) table 
    for (unsigned int qIndx = 0; qIndx < pdfTable_.at(pdfIndx).size(); ++qIndx) {
       pdfTable_.at(pdfIndx).at(qIndx).resize(nxPoints_);
      for (unsigned int xIndx = 0; xIndx < pdfTable_.at(pdfIndx).at(qIndx).size();
	   ++xIndx) {
        datafile >> pdfTable_.at(pdfIndx).at(qIndx).at(xIndx); 
	if(datafile.eof()) 
	  throw Exception() << "Error while reading " << fileName_.at(pdfIndx) 
			    << " too few data points in file" 
			    << "in PomeronPDF::loadTables()"
			    << Exception::runerror;
      }   
    }
    datafile >> dtemp;
    if(!datafile.eof()) 
      throw Exception() << "Error reading end of " << fileName_.at(pdfIndx)  
			<< " too many data points in file" 
			<< "in PomeronPDF::loadTables()"
			<< Exception::runerror;
    datafile.close();
  }
}

void PomeronPDF::doinit() {
  PDFBase::doinit();  
  switch(PDFFit_){
  case 0:
    fileName_.at(0) = "2007/h12007jetsdpdf_charm.data";
    fileName_.at(1) = "2007/h12007jetsdpdf_gluon.data";
    fileName_.at(2) = "2007/h12007jetsdpdf_singlet.data";
    break;
  case 1:
    fileName_.at(0) = "2006/h12006jetspdf_singlet_fitA.data";
    fileName_.at(1) = "2006/h12006jetspdf_gluon_fitA.data";
    fileName_.at(2) = "2006/h12006jetspdf_singlet_fitA.data";
    nxPoints_ = 100;
    nqPoints_ = 30;
    break;
  case 2:
    fileName_.at(0) = "2006/h12006jetspdf_singlet_fitB.data";
    fileName_.at(1) = "2006/h12006jetspdf_gluon_fitB.data";
    fileName_.at(2) = "2006/h12006jetspdf_singlet_fitB.data";
    nxPoints_ = 100;
    nqPoints_ = 30;
    break;
  default:
    assert(false);
  }
  loadTables();
}

ClassDescription<PomeronPDF> PomeronPDF::initPomeronPDF; 
// Definition of the static class description member.

void PomeronPDF::Init(){

static ClassDocumentation<PomeronPDF> documentation
    ("Implementation of the Diffractive PDFs");
 
 static Switch<PomeronPDF,int> interfacePDFFit
   ("PDFFit",
    "Switch between different PDF fits.",
    &PomeronPDF::PDFFit_, 0, true, false);
 static SwitchOption interfacePDFFit2007
   (interfacePDFFit,
    "2007",
    "Hera fit 2007",
    0);
 static SwitchOption interfacePDFFit2006A
   (interfacePDFFit,
    "2006A",
    "Hera fit A 2006",
    1);
 static SwitchOption interfacePDFFit2006B
   (interfacePDFFit,
    "2006B",
    "Hera fit B 2006",
    2);

  static Parameter<PomeronPDF,string> interfaceRootName
    ("RootName",
     "Root name for the input files",
     &PomeronPDF::rootName_, "",
     false, false);

  static Switch<PomeronPDF,int> interfaceBoundary
    ("Boundary",
     "Switch between different beheviours at out of PDF range.",
     &PomeronPDF::boundary_, 0, true, false);
  static SwitchOption interfaceBounderyExtrapolate
    (interfaceBoundary,
     "Extrapolate",
     "Extrapolate the PDF when it is out of range.",
     0);
  static SwitchOption interfaceBounderyFreeze
    (interfaceBoundary,
     "Freeze",
     "Freeze the value of PDF at the boundary.",
     1);
}

IBPtr PomeronPDF::clone() const {
  return new_ptr(*this);
}

IBPtr PomeronPDF::fullclone() const {
  return new_ptr(*this);
}

void PomeronPDF::persistentOutput(PersistentOStream & os) const {
  os << lxGrid_ << lqqGrid_ << pdfTable_ 
     << PDFFit_ << nxPoints_ << nqPoints_
     << rootName_ << boundary_;
}

void PomeronPDF::persistentInput(PersistentIStream & is, int) {
  is >> lxGrid_ >> lqqGrid_ >> pdfTable_ 
     >> PDFFit_ >> nxPoints_>> nqPoints_
     >> rootName_ >> boundary_;
}
