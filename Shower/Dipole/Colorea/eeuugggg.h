  //==========================================================================
  // This file has been automatically generated for C++ Standalone by
  // MadGraph5_aMC@NLO v. 2.5.4, 2017-03-28
  // By the MadGraph5_aMC@NLO Development Team
  // Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
  //==========================================================================
  // and was then modified by J. Bellm.
#ifndef MG5_Sigma_sm_epem_uuxgggg_H
#define MG5_Sigma_sm_epem_uuxgggg_H
#include <complex> 
#include <vector>
#include "Parameters_sm.h"
using namespace std;
class eeuugggg
{
  public:
  eeuugggg(){ 
    initProc("param_card.dat");
  }
    // Retrun the prefered gluon order
  vector<int> producePermutation(double r,vector < double * > & momenta);
  void setMass(int i,double M){mME[i]=M;}
  private:
    ///////// Process Dependent
  static const int ninitial = 2;
  static const int nexternal = 8;
  static const int nprocesses = 1;
  static const int nwavefuncs = 205;
  static const int namplitudes = 1032;
  std::complex<double> w[nwavefuncs][18];
  double matrix_1_epem_uuxgggg();
    ///////// Generic
  double matrix_element[nprocesses];
  std::complex<double> amp[namplitudes];
  double * jamp2[nprocesses];
  Parameters_sm * pars;
  vector<double> mME;
  vector < double * > p;
  int id1, id2;
  void initProc(string param_card_name);
  void sigmaKin();
  const vector<double> & getMasses() const {return mME;}
  vector < double * > getMomenta(){return p;}
  void setMomenta(vector < double * > & momenta){p = momenta;}
  void setInitial(int inid1, int inid2){id1 = inid1; id2 = inid2;}
  void calculate_wavefunctions(const int perm[], const int hel[]);
    // New function
  double get_jamp2(int i);
  int colorstring(int i, int j);
  int NCol();
}; 
#endif  // MG5_Sigma_sm_epem_uuxgggg_H
