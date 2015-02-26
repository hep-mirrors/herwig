//calculate the LO matrix element squared as given by hpair, constants are taken care of by me2()
//both top and bottom loops are considered.
double MEHiggsPair::MATRIX(double S, double T,double U, double M1, double M2) const {
  //C--LO MATRIX ELEMENT FOR GG -> HH
  using Constants::pi;

  Complex A1,A2,H1,H2,Z1,Z2,AA1,AA2,HH1,HH2,AH1,AH2;
  Complex F1,F2,PROH, PROHEAVY;
  
  //bottom and top masses
  double AMT = _topmass/GeV;
  double AMB = _bottommass/GeV;
  
  //W mass
  double MW = _Wmass/GeV;

  //higgs width and mass
  double FACH=1; 
  double GAMH= _wh/GeV; 
  double AMH = _mh/GeV;

  //vacuum expectation value, self-coupling, top and bottom Yukawas
  double V =  1./sqrt(sqrt(2) * SM().fermiConstant()*GeV*GeV);
  double GHHH =  _selfcoupling * 3.*sqr(AMH) / V;
  double GHT=AMT/V;
  double GHB=AMB/V;

  double ALPS=SM().alphaS(scale());

  double eftscale = _EFTScale/GeV;

  //Higgs propagator
  PROH = Complex(S-sqr(AMH),AMH*GAMH*FACH);

  //EFT-modified self-coupling
  double GHHH_EFT = 0;
  if(_process == 4) { 
    GHHH_EFT = (3.*sqr(AMH) / V) * ( - (3./2.) * _cH + _c6 ) ;  
  }
  
  //EFT-modified top/bottom Yukawas
  double GHT_EFT = 0, GHB_EFT = 0;
  if(_process == 4) { 
    GHT_EFT = AMT/V * ( _ct1 - 0.5 * _cH );
    GHB_EFT = AMB/V * ( _cb1 - 0.5 * _cH );
  }
  //Heavy Higgs
  double AMHEAVY = _heavyHmass/GeV;
  double GAMHEAVY = _heavyHwidth/GeV;			      
  double FACHEAVY = 1.;
  PROHEAVY = Complex(S-sqr(AMHEAVY),AMHEAVY*GAMHEAVY*FACHEAVY);
  double GHhh = _hhHcoupling * 3.*sqr(AMHEAVY) / V;  // Hhh coupling (for _process == 3)
  
  //calculation of scalar integrals C_ij and D_ijk
  vector<Complex> ScalFacs;
  ScalFacs = iniscal(AMT, S, T, U, M1, M2);
  double amq[1] = {AMT};
  double m1[1] = {M1};
  double m2[1] = {M2};
  double ss[1] = {S};
  double tt[1] = {T};
  double uu[1] = {U};

  Complex C0AB[1] = {ScalFacs[0]};
  Complex C0AC[1] = {ScalFacs[1]};
  Complex C0AD[1] = {ScalFacs[2]};
  Complex C0BC[1] = {ScalFacs[3]};
  Complex C0BD[1] = {ScalFacs[4]};
  Complex C0CD[1] = {ScalFacs[5]};
  Complex D0ABC[1] = {ScalFacs[6]};
  Complex D0BAC[1] = {ScalFacs[7]};
  Complex D0ACB[1] = {ScalFacs[8]};
    //form factors for top quark contributions
  formfac_(amq, ss, tt, uu, m1, m2, C0AB, C0AC, C0AD, C0BC, C0BD, C0CD, D0ABC, D0BAC, D0ACB);

  // this factor needs to be divided out to go from SM EFT to DIM-6 EFT
  double CH = ALPS / (8. * pi);
  
  F1 = Complex(0.,0.);
  F2 = Complex(0.,0.);
  //triangle
  if(_process==0||_process==1||_process==3||_process==4) { 
    F1 = F1 + (form_.H1*(GHT*GHHH/PROH)); //(1)
  }
  //triangle: EFT
  if(_process==4) { 
    F1 = F1 + (form_.H1*(GHT*GHHH_EFT+GHT_EFT*GHHH)/PROH);
  }


  //box
  if(_process==0||_process==2||_process==3||_process==4) { 
    F1 = F1 + (form_.HH1*GHT*GHT); //(2)
    F2 = F2 + (form_.HH2*GHT*GHT);
  }
  //box: EFT
  if(_process==4) { 
    F1 = F1 + 2. * (form_.HH1*GHT_EFT*GHT);
    F2 = F2 + 2. * (form_.HH2*GHT_EFT*GHT);
  }

  //additional EFT from gg -> h -> hh diagram:
  if(_process==4) { 
    F1 = F1 + _cg1 * (2.*S/V)*2*(GHHH/PROH); 
    /* divide (1) by form_.H1*AMT/(2*S) to remove F_triangle
     * multiply by 2
     * and multiply by cg1
     */
  }
  //additional EFT from gg -> hh diagram:
  if(_process==4) { 
    F1 = F1 + _cg2 * (2.*S/sqr(V))*2.;
     /* divide (2) by form_.HH1*sqr(AMT)/(2*S) to remove F_box
      * multiply by 2 
      * and multiply by cg2
      */
  }
  //EFT: new t-tbar-h-h diagram
  if(_process==4) { 
    double GttHH = (3.*_ct2 - _cH) *AMT/(2.*sqr(V)); // we have already replaced yf/sqrt(2) = mf/V
    F1 = F1 + (form_.H1*2.*GttHH); //factor of 2 from tthh feynman rule
  }
  //TESTING BELOW

  /* cout << "form_.H1 = " << form_.H1 << " form_.HH1 = " << form_.HH1 << endl; 
  cout << "Ftri = " << real(form_.H1)*AMT/(2*S) << endl; 
  cout << "Fbox = " << real(form_.HH1)*sqr(AMT)/(2*S) << endl; */

  //form factors for bottom quark contributions
  ScalFacs = iniscal(AMB, S, T, U, M1, M2);
  amq[0] = AMB;
  C0AB[0] = ScalFacs[0];
  C0AC[0] = ScalFacs[1];
  C0AD[0] = ScalFacs[2];
  C0BC[0] = ScalFacs[3];
  C0BD[0] = ScalFacs[4];
  C0CD[0] = ScalFacs[5];
  D0ABC[0] = ScalFacs[6];
  D0BAC[0] = ScalFacs[7];
  D0ACB[0] = ScalFacs[8];
  formfac_(amq, ss, tt, uu, m1, m2, C0AB, C0AC, C0AD, C0BC, C0BD, C0CD, D0ABC, D0BAC, D0ACB);
  //triangle 
  if(_process==0||_process==1||_process==3||_process==4) { 
    F1 = F1 + (form_.H1*(GHB*GHHH/PROH));
  }
  //triangle: EFT
  if(_process==4) { 
    F1 = F1 + (form_.H1*(GHB*GHHH_EFT+GHB_EFT*GHHH)/PROH);
  }

  //box
  if(_process==0||_process==2||_process==3||_process==4) { 
    F1 = F1 + (form_.HH1*GHB*GHB);
    F2 = F2 + (form_.HH2*GHB*GHB);
  }
  //box: EFT
  if(_process==4) { 
    F1 = F1 + 2. * (form_.HH1*GHB_EFT*GHB);
    F2 = F2 + 2. * (form_.HH2*GHB_EFT*GHB);
  }

  //EFT: new b-bbar-h-h diagram
  if(_process==4) { 
    double GbbHH = (3.*_cb2 - _cH) *AMB/(2.*sqr(V));// we have already replaced yf/sqrt(2) = mf/V
    F1 = F1 + (form_.H1*2.*GbbHH);
  }
  
  if(_process==3) { 
    m1[0] = AMH;
    m2[0] = AMHEAVY;
    ScalFacs = iniscal(AMT, S, T, U, M1, M2);
    amq[0] = AMT;
    C0AB[0] = ScalFacs[0];
    C0AC[0] = ScalFacs[1];
    C0AD[0] = ScalFacs[2];
    C0BC[0] = ScalFacs[3];
    C0BD[0] = ScalFacs[4];
    C0CD[0] = ScalFacs[5];
    D0ABC[0] = ScalFacs[6];
    D0BAC[0] = ScalFacs[7];
    D0ACB[0] = ScalFacs[8];
    formfac_(amq, ss, tt, uu, m1, m2, C0AB, C0AC, C0AD, C0BC, C0BD, C0CD, D0ABC, D0BAC, D0ACB);
    F1 = F1 + (form_.H1*(GHT*GHhh/PROHEAVY));
    ScalFacs = iniscal(AMB, S, T, U, M1, M2);
    amq[0] = AMB;
    C0AB[0] = ScalFacs[0];
    C0AC[0] = ScalFacs[1];
    C0AD[0] = ScalFacs[2];
    C0BC[0] = ScalFacs[3];
    C0BD[0] = ScalFacs[4];
    C0CD[0] = ScalFacs[5];
    D0ABC[0] = ScalFacs[6];
    D0BAC[0] = ScalFacs[7];
    D0ACB[0] = ScalFacs[8];
    formfac_(amq, ss, tt, uu, m1, m2, C0AB, C0AC, C0AD, C0BC, C0BD, C0CD, D0ABC, D0BAC, D0ACB);
    F1 = F1 + (form_.H1*(GHB*GHhh/PROHEAVY));
  }
  
  //square
  double DMAT = 2. * (norm(F1) + norm(F2));  
  return DMAT;
}


 
