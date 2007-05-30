// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DtoPiPiPiFOCUS class.
//

#include "DtoPiPiPiFOCUS.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

DtoPiPiPiFOCUS::DtoPiPiPiFOCUS() {
  // parameters for the K-matrix
  _malpha.resize(5); _gKK .resize(5); _getaeta .resize(5);
  _gpipi .resize(5); _g4pi.resize(5); _getaetap.resize(5);
  _malpha[0]  = 0.65100*GeV; _gpipi[0]    = 0.24844*GeV;
  _malpha[1]  = 1.20720*GeV; _gpipi[1]    = 0.91779*GeV;
  _malpha[2]  = 1.56122*GeV; _gpipi[2]    = 0.37024*GeV;
  _malpha[3]  = 1.21257*GeV; _gpipi[3]    = 0.34501*GeV;
  _malpha[4]  = 1.81746*GeV; _gpipi[4]    = 0.15770*GeV;
  _gKK[0]     =-0.52523*GeV; _g4pi[0]     = 0.00000*GeV;
  _gKK[1]     = 0.55427*GeV; _g4pi[1]     = 0.00000*GeV;
  _gKK[2]     = 0.23591*GeV; _g4pi[2]     = 0.62605*GeV;
  _gKK[3]     = 0.39642*GeV; _g4pi[3]     = 0.97644*GeV;
  _gKK[4]     =-0.17915*GeV; _g4pi[4]     =-0.90100*GeV;
  _getaeta[0] =-0.38878*GeV; _getaetap[0] =-0.36397*GeV;
  _getaeta[1] = 0.38705*GeV; _getaetap[1] = 0.29448*GeV;
  _getaeta[2] = 0.18409*GeV; _getaetap[2] = 0.18923*GeV;
  _getaeta[3] = 0.19746*GeV; _getaetap[3] = 0.00357*GeV;
  _getaeta[4] =-0.00931*GeV; _getaetap[4] = 0.20689*GeV;

//   /**
//    *  \f$s^{\rm scatt}_0\f$
//    */
//   Energy2 _s0scatt;

//   /**
//    *  \f$s_{A_0}\f$
//    */
//   Energy2 _sA0;

//   /**
//    * \f$s_A\f$
//    */
//   double _sA;

//   /**
//    *  Pion mass
//    */
//   Energy _mpi;

//   /**
//    *  Kaon mass
//    */
//   Energy _mK;

//   /**
//    * \f$\eta\f$ mass
//    */
//   Energy _meta;

//   /**
//    *  \f$\eta'\f$ mass
//    */
//   Energy _metap;

//   /**
//    *  The magnitude of the \f$f_{\rm prod}\f$ couplings for \f$D\f$
//    */
//   vector<double> _fprodmagD;

//   /**
//    *  The phase of the \f$f_{\rm prod}\f$ couplings for \f$D\f$
//    */
//   vector<double> _fprodphaseD;

//   /**
//    * The \f$f_{\rm prod}\f$ couplings for \f$D\f$
//    */
//   vector<Complex> _fprodD;

//   /**
//    *  The magnitude of the \f$f_{\rm prod}\f$ couplings for \f$D_s\f$
//    */
//   vector<double> _fprodmagDs;

//   /**
//    *  The phase of the \f$f_{\rm prod}\f$ couplings for \f$D_s\f$
//    */
//   vector<double> _fprodphaseDs;

//   /**
//    * The \f$f_{\rm prod}\f$ couplings for \f$D_s\f$
//    */
//   vector<Complex> _fprodDs;

//   /**
//    *  The magnitude of the \f$\beta\f$ couplings for \f$D\f$
//    */
//   vector<Energy> _betamagD;

//   /**
//    *  The phase of the \f$\beta\f$ couplings for \f$D\f$
//    */
//   vector<double> _betaphaseD;

//   /**
//    * The \f$\beta\f$ couplings for \f$D\f$
//    */
//   vector<complex<Energy> > _betaD;

//   /**
//    *  The magnitude of the \f$\beta\f$ couplings for \f$D_s\f$
//    */
//   vector<Energy> _betamagDs;

//   /**
//    *  The phase of the \f$\beta\f$ couplings for \f$D_s\f$
//    */
//   vector<double> _betaphaseDs;

//   /**
//    * The \f$\beta\f$ couplings for \f$D_s\f$
//    */
//   vector<complex<Energy> > _betaDs;
//   //@}

//   /**
//    *  Amplitudes and phases for the vector and tensor mesons
//    */
//   //@{
//   /**
//    *  Magnitude for the \f$f_2(1270)\f$ component for \f$D_s^+\to\pi^+\pi^-\pi^+\f$
//    */
//   InvEnergy2 _aDsf2;

//   /**
//    *  Phase for the \f$f_2(1270)\f$ component for \f$D_s^+\to\pi^+\pi^-\pi^+\f$
//    */
//   double _phiDsf2;

//   /**
//    *  Magnitude for the \f$\rho^0(1450)\f$ component for \f$D_s^+\to\pi^+\pi^-\pi^+\f$
//    */
//   double _aDsrho1450;

//   /**
//    *  Phase for the \f$\rho^0(1450)\f$ component for \f$D_s^+\to\pi^+\pi^-\pi^+\f$
//    */
//   double _phiDsrho1450;

//   /**
//    *  Amplitude for the \f$f_2(1270)\f$ component for \f$D_s^+\to\pi^+\pi^-\pi^+\f$
//    */
//   complex<InvEnergy2> _cDsf2;

//   /**
//    *  Amplitude for the \f$\rho^0(1450)\f$ component for \f$D_s^+\to\pi^+\pi^-\pi^+\f$
//    */
//   Complex _cDsrho1450;

//   /**
//    *  Magnitude for the \f$f_2(1270)\f$ component for \f$D^+\to\pi^+\pi^-\pi^+\f$
//    */
//   InvEnergy2 _aDf2;

//   /**
//    *  Phase for the \f$f_2(1270)\f$ component for \f$D^+\to\pi^+\pi^-\pi^+\f$
//    */
//   double _phiDf2;

//   /**
//    *  Magnitude for the \f$\rho^0(770)\f$ component for \f$D^+\to\pi^+\pi^-\pi^+\f$
//    */
//   double _aDrho770;

//   /**
//    *  Phase for the \f$\rho^0(770)\f$ component for \f$D^+\to\pi^+\pi^-\pi^+\f$
//    */
//   double _phiDrho770;

//   /**
//    *  Amplitude for the \f$f_2(1270)\f$ component for \f$D^+\to\pi^+\pi^-\pi^+\f$
//    */
//   complex<InvEnergy2> _cDf2;

//   /**
//    *  Amplitude for the \f$\rho^0(770)\f$ component for \f$D^+\to\pi^+\pi^-\pi^+\f$
//    */
//   Complex _cDrho770;
//   //@}


//   /**
//    *  Masses and Widths for the vector and tensor mesons
//    */
//   //@{
//   /**
//    *  Mass of the \f$\rho(770)\f$
//    */
//   Energy _mrho770;

//   /**
//    *  Mass of the \f$\rho(1450)\f$
//    */
//   Energy _mrho1450;

//   /**
//    *  Mass of the \f$f_2(1270)\f$
//    */
//   Energy _mf2;

//   /**
//    *  Width of the \f$\rho(770)\f$
//    */
//   Energy _wrho770;

//   /**
//    *  Width of the \f$\rho(1450)\f$
//    */
//   Energy _wrho1450;

//   /**
//    *  Width of the \f$f_2(1270)\f$
//    */
//   Energy _wf2;
//   //@}

//   /**
//    * Parameters for the phase-space integration
//    */
//   //@{
//   /**
//    *  Maximum weights for the different modes
//    */
//   vector<double> _maxwgt;

//   /**
//    *  Weights for the different phase-space channels
//    */
//   vector<double> _weights;
//   //@}

}

void DtoPiPiPiFOCUS::persistentOutput(PersistentOStream & os) const {

//   /**
//    *  Parameters for the K-matrix
//    */
//   //@{
//   /**
//    *  Masses of the poles for the K-matrix
//    */
//   vector<Energy> _malpha;

//   /**
//    *  The \f$g_{\pi\pi}\f$ coupling
//    */
//   vector<Energy> _gpipi;

//   /**
//    *  The \f$g_{K\bar{K}}\f$ coupling
//    */
//   vector<Energy> _gKK;

//   /**
//    *  The \f$g_{4\pi}\f$ coupling
//    */
//   vector<Energy> _g4pi;

//   /*
//    *  The \f$g_{\eta\eta}\f$ coupling
//    */
//   vector<Energy> _getaeta;

//   /**
//    *  The \f$g_{\eta\eta'}\f$ coupling
//    */
//   vector<Energy> _getaetap;

//   /**
//    *  The g couplings for easy access
//    */
//   vector<vector<Energy> > _gcoup;

//   /**
//    *  \f$s^{\rm scatt}_0\f$
//    */
//   Energy2 _s0scatt;

//   /**
//    *  \f$s_{A_0}\f$
//    */
//   Energy2 _sA0;

//   /**
//    * \f$s_A\f$
//    */
//   double _sA;

//   /**
//    *  Pion mass
//    */
//   Energy _mpi;

//   /**
//    *  Kaon mass
//    */
//   Energy _mK;

//   /**
//    * \f$\eta\f$ mass
//    */
//   Energy _meta;

//   /**
//    *  \f$\eta'\f$ mass
//    */
//   Energy _metap;

//   /**
//    *  The magnitude of the \f$f_{\rm prod}\f$ couplings for \f$D\f$
//    */
//   vector<double> _fprodmagD;

//   /**
//    *  The phase of the \f$f_{\rm prod}\f$ couplings for \f$D\f$
//    */
//   vector<double> _fprodphaseD;

//   /**
//    * The \f$f_{\rm prod}\f$ couplings for \f$D\f$
//    */
//   vector<Complex> _fprodD;

//   /**
//    *  The magnitude of the \f$f_{\rm prod}\f$ couplings for \f$D_s\f$
//    */
//   vector<double> _fprodmagDs;

//   /**
//    *  The phase of the \f$f_{\rm prod}\f$ couplings for \f$D_s\f$
//    */
//   vector<double> _fprodphaseDs;

//   /**
//    * The \f$f_{\rm prod}\f$ couplings for \f$D_s\f$
//    */
//   vector<Complex> _fprodDs;

//   /**
//    *  The magnitude of the \f$\beta\f$ couplings for \f$D\f$
//    */
//   vector<Energy> _betamagD;

//   /**
//    *  The phase of the \f$\beta\f$ couplings for \f$D\f$
//    */
//   vector<double> _betaphaseD;

//   /**
//    * The \f$\beta\f$ couplings for \f$D\f$
//    */
//   vector<complex<Energy> > _betaD;

//   /**
//    *  The magnitude of the \f$\beta\f$ couplings for \f$D_s\f$
//    */
//   vector<Energy> _betamagDs;

//   /**
//    *  The phase of the \f$\beta\f$ couplings for \f$D_s\f$
//    */
//   vector<double> _betaphaseDs;

//   /**
//    * The \f$\beta\f$ couplings for \f$D_s\f$
//    */
//   vector<complex<Energy> > _betaDs;
//   //@}

//   /**
//    *  Amplitudes and phases for the vector and tensor mesons
//    */
//   //@{
//   /**
//    *  Magnitude for the \f$f_2(1270)\f$ component for \f$D_s^+\to\pi^+\pi^-\pi^+\f$
//    */
//   InvEnergy2 _aDsf2;

//   /**
//    *  Phase for the \f$f_2(1270)\f$ component for \f$D_s^+\to\pi^+\pi^-\pi^+\f$
//    */
//   double _phiDsf2;

//   /**
//    *  Magnitude for the \f$\rho^0(1450)\f$ component for \f$D_s^+\to\pi^+\pi^-\pi^+\f$
//    */
//   double _aDsrho1450;

//   /**
//    *  Phase for the \f$\rho^0(1450)\f$ component for \f$D_s^+\to\pi^+\pi^-\pi^+\f$
//    */
//   double _phiDsrho1450;

//   /**
//    *  Amplitude for the \f$f_2(1270)\f$ component for \f$D_s^+\to\pi^+\pi^-\pi^+\f$
//    */
//   complex<InvEnergy2> _cDsf2;

//   /**
//    *  Amplitude for the \f$\rho^0(1450)\f$ component for \f$D_s^+\to\pi^+\pi^-\pi^+\f$
//    */
//   Complex _cDsrho1450;

//   /**
//    *  Magnitude for the \f$f_2(1270)\f$ component for \f$D^+\to\pi^+\pi^-\pi^+\f$
//    */
//   InvEnergy2 _aDf2;

//   /**
//    *  Phase for the \f$f_2(1270)\f$ component for \f$D^+\to\pi^+\pi^-\pi^+\f$
//    */
//   double _phiDf2;

//   /**
//    *  Magnitude for the \f$\rho^0(770)\f$ component for \f$D^+\to\pi^+\pi^-\pi^+\f$
//    */
//   double _aDrho770;

//   /**
//    *  Phase for the \f$\rho^0(770)\f$ component for \f$D^+\to\pi^+\pi^-\pi^+\f$
//    */
//   double _phiDrho770;

//   /**
//    *  Amplitude for the \f$f_2(1270)\f$ component for \f$D^+\to\pi^+\pi^-\pi^+\f$
//    */
//   complex<InvEnergy2> _cDf2;

//   /**
//    *  Amplitude for the \f$\rho^0(770)\f$ component for \f$D^+\to\pi^+\pi^-\pi^+\f$
//    */
//   Complex _cDrho770;
//   //@}


//   /**
//    *  Masses and Widths for the vector and tensor mesons
//    */
//   //@{
//   /**
//    *  Mass of the \f$\rho(770)\f$
//    */
//   Energy _mrho770;

//   /**
//    *  Mass of the \f$\rho(1450)\f$
//    */
//   Energy _mrho1450;

//   /**
//    *  Mass of the \f$f_2(1270)\f$
//    */
//   Energy _mf2;

//   /**
//    *  Width of the \f$\rho(770)\f$
//    */
//   Energy _wrho770;

//   /**
//    *  Width of the \f$\rho(1450)\f$
//    */
//   Energy _wrho1450;

//   /**
//    *  Width of the \f$f_2(1270)\f$
//    */
//   Energy _wf2;
//   //@}

//   /**
//    * Parameters for the phase-space integration
//    */
//   //@{
//   /**
//    *  Maximum weights for the different modes
//    */
//   vector<double> _maxwgt;

//   /**
//    *  Weights for the different phase-space channels
//    */
//   vector<double> _weights;
//   //@}

}

void DtoPiPiPiFOCUS::persistentInput(PersistentIStream & is, int) {
}

ClassDescription<DtoPiPiPiFOCUS> DtoPiPiPiFOCUS::initDtoPiPiPiFOCUS;
// Definition of the static class description member.

void DtoPiPiPiFOCUS::Init() {

  static ClassDocumentation<DtoPiPiPiFOCUS> documentation
    ("There is no documentation for the DtoPiPiPiFOCUS class");

//   /**
//    *  Parameters for the K-matrix
//    */
//   //@{
//   /**
//    *  Masses of the poles for the K-matrix
//    */
//   vector<Energy> _malpha;

//   /**
//    *  The \f$g_{\pi\pi}\f$ coupling
//    */
//   vector<Energy> _gpipi;

//   /**
//    *  The \f$g_{K\bar{K}}\f$ coupling
//    */
//   vector<Energy> _gKK;

//   /**
//    *  The \f$g_{4\pi}\f$ coupling
//    */
//   vector<Energy> _g4pi;

//   /*
//    *  The \f$g_{\eta\eta}\f$ coupling
//    */
//   vector<Energy> _getaeta;

//   /**
//    *  The \f$g_{\eta\eta'}\f$ coupling
//    */
//   vector<Energy> _getaetap;

//   /**
//    *  The g couplings for easy access
//    */
//   vector<vector<Energy> > _gcoup;

//   /**
//    *  \f$s^{\rm scatt}_0\f$
//    */
//   Energy2 _s0scatt;

//   /**
//    *  \f$s_{A_0}\f$
//    */
//   Energy2 _sA0;

//   /**
//    * \f$s_A\f$
//    */
//   double _sA;

//   /**
//    *  Pion mass
//    */
//   Energy _mpi;

//   /**
//    *  Kaon mass
//    */
//   Energy _mK;

//   /**
//    * \f$\eta\f$ mass
//    */
//   Energy _meta;

//   /**
//    *  \f$\eta'\f$ mass
//    */
//   Energy _metap;

//   /**
//    *  The magnitude of the \f$f_{\rm prod}\f$ couplings for \f$D\f$
//    */
//   vector<double> _fprodmagD;

//   /**
//    *  The phase of the \f$f_{\rm prod}\f$ couplings for \f$D\f$
//    */
//   vector<double> _fprodphaseD;

//   /**
//    * The \f$f_{\rm prod}\f$ couplings for \f$D\f$
//    */
//   vector<Complex> _fprodD;

//   /**
//    *  The magnitude of the \f$f_{\rm prod}\f$ couplings for \f$D_s\f$
//    */
//   vector<double> _fprodmagDs;

//   /**
//    *  The phase of the \f$f_{\rm prod}\f$ couplings for \f$D_s\f$
//    */
//   vector<double> _fprodphaseDs;

//   /**
//    * The \f$f_{\rm prod}\f$ couplings for \f$D_s\f$
//    */
//   vector<Complex> _fprodDs;

//   /**
//    *  The magnitude of the \f$\beta\f$ couplings for \f$D\f$
//    */
//   vector<Energy> _betamagD;

//   /**
//    *  The phase of the \f$\beta\f$ couplings for \f$D\f$
//    */
//   vector<double> _betaphaseD;

//   /**
//    * The \f$\beta\f$ couplings for \f$D\f$
//    */
//   vector<complex<Energy> > _betaD;

//   /**
//    *  The magnitude of the \f$\beta\f$ couplings for \f$D_s\f$
//    */
//   vector<Energy> _betamagDs;

//   /**
//    *  The phase of the \f$\beta\f$ couplings for \f$D_s\f$
//    */
//   vector<double> _betaphaseDs;

//   /**
//    * The \f$\beta\f$ couplings for \f$D_s\f$
//    */
//   vector<complex<Energy> > _betaDs;
//   //@}

//   /**
//    *  Amplitudes and phases for the vector and tensor mesons
//    */
//   //@{
//   /**
//    *  Magnitude for the \f$f_2(1270)\f$ component for \f$D_s^+\to\pi^+\pi^-\pi^+\f$
//    */
//   InvEnergy2 _aDsf2;

//   /**
//    *  Phase for the \f$f_2(1270)\f$ component for \f$D_s^+\to\pi^+\pi^-\pi^+\f$
//    */
//   double _phiDsf2;

//   /**
//    *  Magnitude for the \f$\rho^0(1450)\f$ component for \f$D_s^+\to\pi^+\pi^-\pi^+\f$
//    */
//   double _aDsrho1450;

//   /**
//    *  Phase for the \f$\rho^0(1450)\f$ component for \f$D_s^+\to\pi^+\pi^-\pi^+\f$
//    */
//   double _phiDsrho1450;

//   /**
//    *  Amplitude for the \f$f_2(1270)\f$ component for \f$D_s^+\to\pi^+\pi^-\pi^+\f$
//    */
//   complex<InvEnergy2> _cDsf2;

//   /**
//    *  Amplitude for the \f$\rho^0(1450)\f$ component for \f$D_s^+\to\pi^+\pi^-\pi^+\f$
//    */
//   Complex _cDsrho1450;

//   /**
//    *  Magnitude for the \f$f_2(1270)\f$ component for \f$D^+\to\pi^+\pi^-\pi^+\f$
//    */
//   InvEnergy2 _aDf2;

//   /**
//    *  Phase for the \f$f_2(1270)\f$ component for \f$D^+\to\pi^+\pi^-\pi^+\f$
//    */
//   double _phiDf2;

//   /**
//    *  Magnitude for the \f$\rho^0(770)\f$ component for \f$D^+\to\pi^+\pi^-\pi^+\f$
//    */
//   double _aDrho770;

//   /**
//    *  Phase for the \f$\rho^0(770)\f$ component for \f$D^+\to\pi^+\pi^-\pi^+\f$
//    */
//   double _phiDrho770;

//   /**
//    *  Amplitude for the \f$f_2(1270)\f$ component for \f$D^+\to\pi^+\pi^-\pi^+\f$
//    */
//   complex<InvEnergy2> _cDf2;

//   /**
//    *  Amplitude for the \f$\rho^0(770)\f$ component for \f$D^+\to\pi^+\pi^-\pi^+\f$
//    */
//   Complex _cDrho770;
//   //@}


//   /**
//    *  Masses and Widths for the vector and tensor mesons
//    */
//   //@{
//   /**
//    *  Mass of the \f$\rho(770)\f$
//    */
//   Energy _mrho770;

//   /**
//    *  Mass of the \f$\rho(1450)\f$
//    */
//   Energy _mrho1450;

//   /**
//    *  Mass of the \f$f_2(1270)\f$
//    */
//   Energy _mf2;

//   /**
//    *  Width of the \f$\rho(770)\f$
//    */
//   Energy _wrho770;

//   /**
//    *  Width of the \f$\rho(1450)\f$
//    */
//   Energy _wrho1450;

//   /**
//    *  Width of the \f$f_2(1270)\f$
//    */
//   Energy _wf2;
//   //@}

//   /**
//    * Parameters for the phase-space integration
//    */
//   //@{
//   /**
//    *  Maximum weights for the different modes
//    */
//   vector<double> _maxwgt;

//   /**
//    *  Weights for the different phase-space channels
//    */
//   vector<double> _weights;
//   //@}


}

void DtoPiPiPiFOCUS::doinit() throw(InitException) {
  DecayIntegrator::doinit();

//   /**
//    *  Parameters for the K-matrix
//    */
//   //@{
//   /**
//    *  Masses of the poles for the K-matrix
//    */
//   vector<Energy> _malpha;

//   /**
//    *  The \f$g_{\pi\pi}\f$ coupling
//    */
//   vector<Energy> _gpipi;

//   /**
//    *  The \f$g_{K\bar{K}}\f$ coupling
//    */
//   vector<Energy> _gKK;

//   /**
//    *  The \f$g_{4\pi}\f$ coupling
//    */
//   vector<Energy> _g4pi;

//   /*
//    *  The \f$g_{\eta\eta}\f$ coupling
//    */
//   vector<Energy> _getaeta;

//   /**
//    *  The \f$g_{\eta\eta'}\f$ coupling
//    */
//   vector<Energy> _getaetap;

//   /**
//    *  The g couplings for easy access
//    */
//   vector<vector<Energy> > _gcoup;

//   /**
//    *  \f$s^{\rm scatt}_0\f$
//    */
//   Energy2 _s0scatt;

//   /**
//    *  \f$s_{A_0}\f$
//    */
//   Energy2 _sA0;

//   /**
//    * \f$s_A\f$
//    */
//   double _sA;

//   /**
//    *  Pion mass
//    */
//   Energy _mpi;

//   /**
//    *  Kaon mass
//    */
//   Energy _mK;

//   /**
//    * \f$\eta\f$ mass
//    */
//   Energy _meta;

//   /**
//    *  \f$\eta'\f$ mass
//    */
//   Energy _metap;

//   /**
//    *  The magnitude of the \f$f_{\rm prod}\f$ couplings for \f$D\f$
//    */
//   vector<double> _fprodmagD;

//   /**
//    *  The phase of the \f$f_{\rm prod}\f$ couplings for \f$D\f$
//    */
//   vector<double> _fprodphaseD;

//   /**
//    * The \f$f_{\rm prod}\f$ couplings for \f$D\f$
//    */
//   vector<Complex> _fprodD;

//   /**
//    *  The magnitude of the \f$f_{\rm prod}\f$ couplings for \f$D_s\f$
//    */
//   vector<double> _fprodmagDs;

//   /**
//    *  The phase of the \f$f_{\rm prod}\f$ couplings for \f$D_s\f$
//    */
//   vector<double> _fprodphaseDs;

//   /**
//    * The \f$f_{\rm prod}\f$ couplings for \f$D_s\f$
//    */
//   vector<Complex> _fprodDs;

//   /**
//    *  The magnitude of the \f$\beta\f$ couplings for \f$D\f$
//    */
//   vector<Energy> _betamagD;

//   /**
//    *  The phase of the \f$\beta\f$ couplings for \f$D\f$
//    */
//   vector<double> _betaphaseD;

//   /**
//    * The \f$\beta\f$ couplings for \f$D\f$
//    */
//   vector<complex<Energy> > _betaD;

//   /**
//    *  The magnitude of the \f$\beta\f$ couplings for \f$D_s\f$
//    */
//   vector<Energy> _betamagDs;

//   /**
//    *  The phase of the \f$\beta\f$ couplings for \f$D_s\f$
//    */
//   vector<double> _betaphaseDs;

//   /**
//    * The \f$\beta\f$ couplings for \f$D_s\f$
//    */
//   vector<complex<Energy> > _betaDs;
//   //@}

//   /**
//    *  Amplitudes and phases for the vector and tensor mesons
//    */
//   //@{
//   /**
//    *  Magnitude for the \f$f_2(1270)\f$ component for \f$D_s^+\to\pi^+\pi^-\pi^+\f$
//    */
//   InvEnergy2 _aDsf2;

//   /**
//    *  Phase for the \f$f_2(1270)\f$ component for \f$D_s^+\to\pi^+\pi^-\pi^+\f$
//    */
//   double _phiDsf2;

//   /**
//    *  Magnitude for the \f$\rho^0(1450)\f$ component for \f$D_s^+\to\pi^+\pi^-\pi^+\f$
//    */
//   double _aDsrho1450;

//   /**
//    *  Phase for the \f$\rho^0(1450)\f$ component for \f$D_s^+\to\pi^+\pi^-\pi^+\f$
//    */
//   double _phiDsrho1450;

//   /**
//    *  Amplitude for the \f$f_2(1270)\f$ component for \f$D_s^+\to\pi^+\pi^-\pi^+\f$
//    */
//   complex<InvEnergy2> _cDsf2;

//   /**
//    *  Amplitude for the \f$\rho^0(1450)\f$ component for \f$D_s^+\to\pi^+\pi^-\pi^+\f$
//    */
//   Complex _cDsrho1450;

//   /**
//    *  Magnitude for the \f$f_2(1270)\f$ component for \f$D^+\to\pi^+\pi^-\pi^+\f$
//    */
//   InvEnergy2 _aDf2;

//   /**
//    *  Phase for the \f$f_2(1270)\f$ component for \f$D^+\to\pi^+\pi^-\pi^+\f$
//    */
//   double _phiDf2;

//   /**
//    *  Magnitude for the \f$\rho^0(770)\f$ component for \f$D^+\to\pi^+\pi^-\pi^+\f$
//    */
//   double _aDrho770;

//   /**
//    *  Phase for the \f$\rho^0(770)\f$ component for \f$D^+\to\pi^+\pi^-\pi^+\f$
//    */
//   double _phiDrho770;

//   /**
//    *  Amplitude for the \f$f_2(1270)\f$ component for \f$D^+\to\pi^+\pi^-\pi^+\f$
//    */
//   complex<InvEnergy2> _cDf2;

//   /**
//    *  Amplitude for the \f$\rho^0(770)\f$ component for \f$D^+\to\pi^+\pi^-\pi^+\f$
//    */
//   Complex _cDrho770;
//   //@}


//   /**
//    *  Masses and Widths for the vector and tensor mesons
//    */
//   //@{
//   /**
//    *  Mass of the \f$\rho(770)\f$
//    */
//   Energy _mrho770;

//   /**
//    *  Mass of the \f$\rho(1450)\f$
//    */
//   Energy _mrho1450;

//   /**
//    *  Mass of the \f$f_2(1270)\f$
//    */
//   Energy _mf2;

//   /**
//    *  Width of the \f$\rho(770)\f$
//    */
//   Energy _wrho770;

//   /**
//    *  Width of the \f$\rho(1450)\f$
//    */
//   Energy _wrho1450;

//   /**
//    *  Width of the \f$f_2(1270)\f$
//    */
//   Energy _wf2;
//   //@}

//   /**
//    * Parameters for the phase-space integration
//    */
//   //@{
//   /**
//    *  Maximum weights for the different modes
//    */
//   vector<double> _maxwgt;

//   /**
//    *  Weights for the different phase-space channels
//    */
//   vector<double> _weights;
//   //@}

}
int DtoPiPiPiFOCUS::modeNumber(bool & cc,const DecayMode & dm) const {
}

double DtoPiPiPiFOCUS::me2(bool vertex, const int ichan,const Particle & part,
			   const ParticleVector & decay) const {
}

void DtoPiPiPiFOCUS::dataBaseOutput(ofstream & os,bool header) const {
}

