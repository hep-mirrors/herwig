// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ReweightEW class.
//

#include "ReweightEW.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/Parameter.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

ReweightEW::ReweightEW()
  : filename() {}

ReweightEW::~ReweightEW() {}

IBPtr ReweightEW::clone() const {
  return new_ptr(*this);
}

IBPtr ReweightEW::fullclone() const {
  return new_ptr(*this);
}

double ReweightEW::weight() const{
  int N = subProcess()->outgoing().size();
  if (N != 2) {
    std::cerr << "< 2 particles in FS!\n"; 
    return 1.0; 
  } else {
    Energy2 s = (subProcess()->outgoing()[0]->momentum()
		 +subProcess()->outgoing()[1]->momentum()).m2();
    Energy2 t = (subProcess()->incoming().first->momentum()
		 -subProcess()->outgoing()[0]->momentum()).m2();
    double sg = s/(1.0*GeV2), tg = t/(1.0*GeV2);
    int f1, f2; 
    f1 = subProcess()->incoming().first->id(); 
    f2 = subProcess()->incoming().second->id(); 

    // Flavour of incoming quarks: there's only one correction factor
    // for all flavours in WZ production, so call with u quark.  ZZ, WW,
    // AA: always same and opposite flavours.  Value of f1 matters.
    if (f1 != -f2) {
      f1 = 2; 
    }

    double val = EWKFac(abs(f1), sg, tg); 

    thelasts = sg;
    thelastt = tg; 
    thelastk = val;
    
    return val;
  }
}

void ReweightEW::setSTK(const double &s, const double &t, const double &K) {
  thelasts = s;
  thelastt = t; 
  thelastk = K; 
}

void ReweightEW::inittable() {
  int n = 200; 
  //  int n = 500; 
  FILE *fp;

  // initialize table, n^2 rows in file
  fp = fopen(filename.c_str(),"r");
  if (!fp) {
    cerr << "ReweightEW: Can't open EW file "+filename+". \n"
	 << "You can download EW grid files from\n"
	 << "http://www.hepforge.org/archive/herwig/ewgrids\n";
    // put files to login.hepforge.org:herwig/downloads/ewgrids/
  } else {
    for (int i = 0; i < n*n; ++i) {
      int out = fscanf(fp,"%lf %lf %lf %lf %lf",&tab[i][1],&tab[i][2],&tab[i][3],&tab[i][4],&tab[i][5]);
      if (!out) cerr << "Problems reading EW file " << filename << ".\n";
    }
    fclose(fp);
  }
}

void ReweightEW::doinitrun() {
  inittable();
}

void ReweightEW::doinit() {
  inittable();
}

double ReweightEW::EWKFac(unsigned int f, double s, double t) const {
  int kmin, i, imin, j, jmin, n = 200; //n = 500;
  double tmin, smin, delsmin, deltmin; 
  
  // initialization 
  delsmin=1E30;
  deltmin=1E30;
  imin=0;
  jmin=0;
  smin=0E0;
  tmin=0E0;

  // Find best values for s,t
  for (i=1; i<=n; i++) {
    if (fabs(s-pow(tab[n*i-1][1],2)) < delsmin) {
      smin = pow(tab[n*i-1][1],2);
      delsmin = fabs(s-smin);
      imin = i;
    }
  }
  for (j=1; j<=n; j++) {
    if (fabs(fabs(t)-fabs(tab[n*(imin-1)+j-1][2])) < deltmin) {
      tmin = tab[n*(imin-1)+j-1][2];
      deltmin = fabs(fabs(t)-fabs(tmin));
      jmin = j;
    }
  }
  // Compute correct column
  kmin=n*(imin-1)+jmin-1;

  // return K-factor = 1 + (value delta from table). 
  double val = tab[kmin][5];
  if (f > 0 && f < 6) {
    if (f == 2 || f == 4) val = tab[kmin][3]; // u, c
    if (f == 1 || f == 3) val = tab[kmin][4]; // d, s
    if (f == 5) val = tab[kmin][5]; // b
  } else val = 0.0; // don't reweight. 
  
  return 1.+val;
}




// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void ReweightEW::persistentOutput(PersistentOStream & os) const {
  os << filename;
}

void ReweightEW::persistentInput(PersistentIStream & is, int) {
  is >> filename;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<ReweightEW,ReweightBase>
  describeHerwigReweightEW("Herwig::ReweightEW", "Herwig.so");

void ReweightEW::Init() {

  static ClassDocumentation<ReweightEW> documentation
    ("There is no documentation for the ReweightEW class");

  static Parameter<ReweightEW,string> interfaceFilename
    ("Filename", "Name of the EW K factor file",
     &ReweightEW::filename, "");

}

