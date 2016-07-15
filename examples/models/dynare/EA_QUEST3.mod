//**************************************************************************
// A New Comparative Approach to Macroeconomic Modeling and Policy Analysis
//
// Volker Wieland, Tobias Cwik, Gernot J. Mueller, Sebastian Schmidt and
// Maik Wolters
//
// Working Paper, 2009
//**************************************************************************

// Model: EA_QUEST3

// Further references:
// Ratto, M., W. Roeger, and J. in �t Veld. 2009. "QUEST III: An Estimated Open-Economy DSGE
// Model of the Euro Area with Fiscal and Monetary Policy." Economic Modelling 26(1) pp. 222-233.

// Last edited: 10/08/26 by S. Schmidt


var
E_BGYN E_BWRY E_CLCSN E_DBGYN  E_LER E_ETA E_GC E_GCL E_GCLC E_GCNLC E_GE E_GEX E_GEXL E_GG E_GGL E_GI E_GIG E_GIL
E_GIM E_GIML E_GK E_GKG E_GL E_GSN E_GTAX E_GTFP E_GTFPUCAP E_GTR E_GUC E_GUCAP E_GWRY E_GY E_GYL E_GYPOT E_GYW
E_INOM E_INOMW E_LL E_LL0 E_LBGYN E_LCSN E_LCLCSN  E_LCNLCSN  E_LEXYN E_LGSN  E_LIGSN E_LIMYN  E_LIK E_LIKG
E_LISN E_LOL E_LPCP E_LPMP E_LPXP E_LTRYN E_LUCYN E_LUCLCYN E_LYGAP E_LYKPPI E_LYWR E_LYWY E_MRY E_PHI
E_PHIC E_PHIPI E_PHIM E_PHIML E_PHIW E_PHIX E_PHIXL E_Q E_R E_TAXYN E_TBYN E_TRTAXYN E_TRW E_TRYN E_TW E_UCAP
E_UCAP0 E_VL E_VLLC E_WPHI E_WRPHI E_WS  E_WSW E_ZEPS_C E_ZEPS_ETA E_ZEPS_ETAM E_ZEPS_ETAX E_ZEPS_EX E_ZEPS_G
E_ZEPS_IG E_ZEPS_L E_ZEPS_M E_ZEPS_PPI E_ZEPS_RPREME E_ZEPS_RPREMK	E_ZEPS_TR   E_ZEPS_W  E_ZPHIT E_LCY E_LGY E_LWS E_EPS_G


//**************************************************************************
// Modelbase Variables                                                   //*
    interest inflation inflationq outputgap fispol;                      //*
//**************************************************************************


varexo
E_EPS_C  E_EPS_ETA  E_EPS_ETAM E_EPS_ETAX  E_EPS_EX  E_EPS_IG E_EPS_INOMW E_EPS_L E_EPS_LOL E_EPS_M  E_EPS_PPI
E_EPS_PW    E_EPS_RPREME	E_EPS_RPREMK	E_EPS_TR    E_EPS_W   E_EPS_Y  E_EPS_YW


//**************************************************************************
// Modelbase Shocks                                                      //*
    interest_ fiscal_  ;                                                 //*
//**************************************************************************


parameters
//**************************************************************************
// Modelbase Parameters                                                  //*
                                                                         //*
        cofintintb1 cofintintb2 cofintintb3 cofintintb4                  //*
        cofintinf0 cofintinfb1 cofintinfb2 cofintinfb3 cofintinfb4       //*
        cofintinff1 cofintinff2 cofintinff3 cofintinff4                  //*
        cofintout cofintoutb1 cofintoutb2 cofintoutb3 cofintoutb4        //*
        cofintoutf1 cofintoutf2 cofintoutf3 cofintoutf4                  //*
        cofintoutp cofintoutpb1 cofintoutpb2 cofintoutpb3 cofintoutpb4   //*
        cofintoutpf1 cofintoutpf2 cofintoutpf3 cofintoutpf4              //*
        std_r_ std_r_quart coffispol                                     //*
//**************************************************************************
A1E A2E  ALPHAX  ALPHAE  ALPHAGE  BETAE BGADJ1 BGADJ2 BGTAR DELTAE DELTAGE DGEX DGIM DGPM  DGPX DDYN  E_EX_INOMW
E_EX_R   E_EX_RW G1E  GAMI2E  GAMIE  GAMLE GAMPE  GAMPME GAMPXE  GAMWE GP0 GPCPI0 GPOP0  GPW0 GSLAG  GVECM
GFLAG  GLAGFLAG IGLAGFLAG IGFLAG TRFLAG TWFLAG  GEXOFLAG IGEXOFLAG TREXOFLAG GAMIFLAG SLCFLAG GSN  GTFP0 GY0
GYW0 HABE HABLE IGSLAG IGVECM ILAGE INFLAGE IG1E IGSN ISN KAPPAE L0 LOL LYWY0 OMEGE RHOCE RHOETA RHOETAM  RHOETAX
RHOEXE RHOGE RHOIG RHOL0 RHOLE  RHOLOL RHOPPI1 RHOPPI2 RHOPPI3 RHOPPI4  RHOPCPM  RHOPWPX  RHORPE  RHORPK  RHOUCAP0
RII  RIP  RIX RPI RPP RPX RXI RXP RXX RXY  RPREME RPREMK SE SFPE SFPME SFPXE SFWE SIGC SIGEXE SIGIME  SLC SSC
TAUE TP  THETAE TINFE TR1E TRSN RHOTR TYE1 TYE2 TVAT TW0 TW1 UCAP0 WRLAG ZETE interestq_exog inflationannual_exog;

//estimated parameters (mean posterior distribution)
A2E        =   0.0453;
G1E        =   -0.0754;
GAMIE      =   76.0366;
GAMI2E     =  1.1216;
GAMLE      =   58.2083;
GAMPE      =   61.4415;
GAMPME     =   1.6782;
GAMPXE     =   26.1294;
GAMWE      =   1.2919;
GSLAG      =   -0.4227;
GVECM =        -0.1567;
HABE       =   0.5634;
HABLE       =  0.8089;
IG1E     =   0.1497;
IGSLAG   = 0.4475;
IGVECM = -0.1222;
ILAGE      =   0.9009;
KAPPAE     =   1.9224;
RHOCE      =   0.9144;
RHOETA    =   0.1095 ;
RHOETAM     =   0.9557;
RHOETAX     =   0.8109;
RHOGE     =   0.2983;
RHOIG   = 0.8530;
RHOLE     =   0.9750;
RHOL0     =   0.9334;
RHOPCPM  =   0.6652;
RHOPWPX =  0.2159;
RHORPE     =   0.9842;
RHORPK     =   0.9148;
RHOUCAP0    =   0.9517;
RPREME     =   0.0200;
RPREMK     =   0.0245;
SE         =   0.8588;
SFPE       =   0.8714;
SFPME      =   0.7361;
SFPXE      =   0.9180;
SFWE       =   0.7736;
SIGC       =  4.0962;
SIGEXE     =  2.5358;
SIGIME     =  1.1724;
SLC        =  0.3507;
TINFE      =   1.9590;
TR1E       =   0.9183;
RHOTR     =   0.8636;
TYE1        =   0.4274;
TYE2       =   0.0783;
WRLAG      =   0.2653;


// calibrated parameters;
ALPHAX        =   0.5;
BETAE      =   0.996;
ALPHAE     =   0.52;
ALPHAGE    =   0.9;
BGADJ2      =  0.004;
BGADJ1      =  0.001*BGADJ2;
BGTAR      =   2.4;
DELTAE     =   0.025;
DELTAGE     =   0.0125;
DGIM = 0.00738619107021;  // log of dgim in ddd
DGEX = 0.00738619107021;
DGPM = -0.00396650612294;
DGPX = -0.00396650612294;
INFLAGE    =   0.98;
LOL        =   0;
RHOEXE     =   0.975;
RHOLOL     =   0.99;
SSC        =   0.2;
TAUE       =   0.1000;
THETAE     =   1.6;
TP         =   0.2;
TRSN       =   0.36;
TVAT       =   0.2;
TW0        =   0.2;
TW1        =   0.8;
ZETE       =   0.4000;

// SWITCHES for various simulations
DDYN = 1;
GFLAG=1;
GLAGFLAG=1;
IGLAGFLAG=1;
IGFLAG=1;
TRFLAG=1;
TWFLAG=1;
GEXOFLAG=1;
IGEXOFLAG=1;
TREXOFLAG=1;
SLCFLAG=1;
GAMIFLAG=1;

// steady state parameters
A1E        =    0.0669; % this is actually a function of the steady state and estimated parameters to get UCAP0=1;
OMEGE      =   1.4836; % this is actually a function of the steady state and estimated parameters to get L0=0.65;
GSN        =   0.203;
IGSN      =    0.025;
//ISN        =   0.17; % used only when DDYN=0, and RPREMK would be a function of ISN;
GPCPI0     =   0;
GP0        =   0.005;
GPW0       =   GP0;
GPOP0      =   0.00113377677398;
GY0        =   0.003;
GTFP0      =   (ALPHAE+ALPHAGE-1)/ALPHAE*GY0-(2-ALPHAE-ALPHAGE)/ALPHAE*GPCPI0;
GYW0       =   GY0;
UCAP0      =   1;
L0         =  0.65;
E_EX_R     =   1/BETAE-1;    %   -log(BETAE);
E_EX_INOMW =  E_EX_R+GP0;
E_EX_RW = E_EX_INOMW - GPW0;
LYWY0 = 0;

// VAR model for exogenous world variables, not jointly estimated
RXY = -0.0001;
RII = 0.887131978334279;
RIP = 0.147455589872832;
RIX = 0.120095599681076;
RPI = 0.112067224979767;
RPP = 0.502758194258928;
RPX = 0.082535400836409;
RXI = 0.073730131176521;
RXP = -0.302655015002645;
RXX = 0.495000246553550;

// AR(4) for exogenous PPI process, not jointly estimated
RHOPPI1 =   0.24797097628284;
RHOPPI2 =  0.13739098460472;
RHOPPI3 =  0.10483962746747;
RHOPPI4 =  0.09282876044442;

interestq_exog=1.00901606;
inflationannual_exog=1.02;

//**************************************************************************
// Specification of Modelbase Parameters                                 //*
                                                                         //*
// Load Modelbase Monetary Policy Parameters                             //*
//thispath = cd;
//cd('..');
//load policy_param.mat;
//for i=1:33
//    deep_parameter_name = M_.param_names(i,:);
//    eval(['M_.params(i)  = ' deep_parameter_name ' ;'])
//end
//cd(thispath);
//std_r_ = 100;                                                            //*
                                                                         //*
// Definition of Discretionary Fiscal Policy Parameter                   //*
coffispol = 1/GSN;                                                       //*
//**************************************************************************

model;

//**************************************************************************
// Definition of Modelbase Variables in Terms of Original Model Variables //*

interest   = 100*((E_INOM+1)^4-interestq_exog^4)/interestq_exog^4;                //*
inflation = (1/2)*(inflationq+inflationq(-1));      //*
inflationq = 100*(4*E_PHIC+1-inflationannual_exog)/inflationannual_exog;          //*
outputgap  = 100*E_LYGAP;                                                         //*
fispol = E_EPS_G;                                                                 //*
 //**************************************************************************


//**************************************************************************
// Policy Rule                                                           //*
                                                                         //*
// Monetary Policy                                                       //*
                                                                         //*
interest = cofintintb1*interest(-1)                                      //*
          //  + cofintintb2*interest(-2)                                    //*
          //  + cofintintb3*interest(-3)                                    //*
          //  + cofintintb4*interest(-4)                                    //*
           + cofintinf0*inflationq                                       //*
           + cofintinfb1*inflationq(-1)                                  //*
          //  + cofintinfb2*inflationq(-2)                                  //*
          //  + cofintinfb3*inflationq(-3)                                  //*
          //  + cofintinfb4*inflationq(-4)                                  //*
           + cofintinff1*inflationq(+1)                                  //*
          //  + cofintinff2*inflationq(+2)                                  //*
          //  + cofintinff3*inflationq(+3)                                  //*
          //  + cofintinff4*inflationq(+4)                                  //*
           + cofintout*outputgap                                         //*
           + cofintoutb1*outputgap(-1)                                   //*
          //  + cofintoutb2*outputgap(-2)                                   //*
          //  + cofintoutb3*outputgap(-3)                                   //*
          //  + cofintoutb4*outputgap(-4)                                   //*
           + cofintoutf1*outputgap(+1)                                   //*
          //  + cofintoutf2*outputgap(+2)                                   //*
          //  + cofintoutf3*outputgap(+3)                                   //*
          //  + cofintoutf4*outputgap(+4)                                   //*
           + std_r_ * interest_;                                         //*
                                                                         //*
// Discretionary Government Spending                                     //*
                                                                         //*
fispol = coffispol*fiscal_;                                              //*
//**************************************************************************

//old Monetary policy rule
//E_INOM = ILAGE*E_INOM(-1)+(1-ILAGE)*(E_EX_R + E_ZPHIT + TINFE*(E_PHIC-E_ZPHIT) + TYE1*E_LYGAP(-1) )
  //            + TYE2*(E_LYGAP-E_LYGAP(-1)) + E_ZEPS_M;

//Derivative of utility function of unconstrained hh w.r.t consumption
exp(E_LUCYN) = exp(E_ZEPS_C)*(exp(E_LCNLCSN)*(1-HABE/(1+E_GCNLC-GY0)))^(-SIGC)*
    (1-OMEGE*exp(E_ZEPS_L)*(exp(E_LL)-HABLE*exp(E_LL(-1)))^KAPPAE)^(1-SIGC);                       //   1: EUC

//Derivative of utility function of constrained hh w.r.t consumption
exp(E_LUCLCYN) = exp(E_LCLCSN)^(-SIGC)*
    (1-OMEGE*exp(E_ZEPS_L)*(exp(E_LL)-HABLE*exp(E_LL(-1)))^KAPPAE)^(1-SIGC);                       //   2: EUCLC

//Derivative of utility function of unconstrained hh w.r.t leisure
E_VL = exp(E_ZEPS_C)*(exp(E_LCNLCSN)*(1-HABE/(1+E_GCNLC-GY0)))^(1-SIGC)
    *(1-OMEGE*exp(E_ZEPS_L)*(exp(E_LL)-HABLE*exp(E_LL(-1)))^KAPPAE)^(-SIGC)
    *KAPPAE*OMEGE*(exp(E_LL)-HABLE*exp(E_LL(-1)))^(KAPPAE-1)*exp(E_ZEPS_L);                        //   3: EVL

//Derivative of utility function of constrained hh w.r.t leisure
E_VLLC = exp(E_LCLCSN)^(1-SIGC)
    *(1-OMEGE*exp(E_ZEPS_L)*(exp(E_LL)-HABLE*exp(E_LL(-1)))^KAPPAE)^(-SIGC)
    *KAPPAE*OMEGE*(exp(E_LL)-HABLE*exp(E_LL(-1)))^(KAPPAE-1)*exp(E_ZEPS_L);                        //   4: EVL


1/BETAE-1 = E_GUC(1) + E_INOM - E_PHIC(1);                                                         //   5: EC

E_LUCYN-E_LUCYN(-1) = E_GUC + SIGC*(E_GY-GY0 - E_PHIC + E_PHI);                                    //   6: ELAM

//Non-Ricardian hh consume after-tax income
(1+TVAT)*exp(E_LCLCSN) = (1-E_TW-SSC)*E_WS + E_TRW*E_WS - E_TAXYN;                                 //   7: CLC

E_WS = exp(E_LL-E_LYWR);                                                                           //   8: WS

//aggregate consumption
exp(E_LCSN)=SLCFLAG*SLC*exp(E_LCLCSN)+(1-SLCFLAG*SLC)*exp(E_LCNLCSN);                              //   9: EC

//Wage adjustment with wage markup included
//(1+TVAT)*((1-SLC)*E_VL+SLC*E_VLLC )/((1-SLC)*exp(E_LUCYN)+SLC*exp(E_LUCLCYN)) =
(1+TVAT)*(((1-SLC)*E_VL+SLC*E_VLLC )/((1-SLC)*exp(E_LUCYN)+SLC*exp(E_LUCLCYN)))^(1-WRLAG)*((1-E_TW-SSC)/(1+TVAT)*(THETAE-1)/THETAE/exp(E_LYWR(-1))/(1+E_GY-GY0))^WRLAG =
   ((THETAE-1)/THETAE/exp(E_LYWR)*(1-E_TW-SSC)+GAMIFLAG*GAMWE/THETAE/exp(E_LYWR)*((E_WPHI-GP0-GY0)-(1-SFWE)*(E_PHI(-1)-GP0)))-
   BETAE*GAMIFLAG*GAMWE/THETAE/exp(E_LYWR)*((E_WPHI(1)-GP0-GY0)-(1-SFWE)*(E_PHI-GP0));             //  10: EW

(1+E_ZEPS_W)/(exp(E_LYWR)) = E_ETA*ALPHAE/exp(E_LL)*(1+E_LOL) - 1/(exp(E_LYWR))*GAMLE*(E_LL-E_LL(-1)) +
           1/(exp(E_LYWR(+1)))*(1+E_GY(+1)-GY0)*GAMLE/(1+E_R)*(E_LL(+1)-E_LL);                     //  11: EL

//Derivative of non-ricardian hh w.r.t. physical investment
GAMIE*(exp(E_LIK)-(DELTAE +GPOP0+GY0+GPCPI0)) + GAMI2E*(E_GI-GY0-GPCPI0)
        - GAMI2E/(1+E_INOM)*(E_GI(+1)-GY0-GPCPI0) = E_Q-1;                                         //  12: ENU

//present discounted value of the rental of return from investing in real assets (Q_t)
E_ETA*(1-TP)*(1-ALPHAE)*(exp(E_LYKPPI)) =
    E_Q-(1-E_R-DELTAE-RPREMK-E_ZEPS_RPREMK-E_PHIPI(+1)+GPCPI0)*E_Q(1) +
    (A1E*(E_UCAP-UCAP0)+A2E*(E_UCAP-UCAP0)^2)*(1-TP);                                              //  13: ERK

(E_ETA*(1-ALPHAE)*(exp(E_LYKPPI) ) ) =
        (A1E+2*A2E*(E_UCAP-UCAP0))*E_UCAP;                                                         //  14: EK

//aggregate price markup
E_ETA = 1 - (TAUE+E_ZEPS_ETA) -
     GAMIFLAG*GAMPE*(BETAE*(SFPE*E_PHI(1)+(1-SFPE)*E_PHI(-1)-GP0)- (E_PHI-GP0));                   //  15: EETA

E_MRY = (1+E_INOM)^(-ZETE);                                                                        //  16: EM

E_GK-(GY0+GPCPI0) = exp(E_LIK) -(DELTAE +GPOP0+GY0+GPCPI0);                                        //  17: EI

E_GKG-(GY0+GPCPI0)= exp(E_LIKG)-(DELTAGE+GPOP0+GY0+GPCPI0);                                        //  18: EIG

E_LISN = -E_LYKPPI+E_LIK+GY0+GPCPI0-E_GK;                                                          //  19: ISN

//production function
E_GY = (1-ALPHAE)*(E_GK+E_GUCAP) +  ALPHAE*(E_GTFP + E_GL*(1+LOL)) + E_GKG*(1-ALPHAGE);            //  20: EY

//Fisher Gleichung
E_R = E_INOM-E_PHI(1);                                                                             //  21: EINOM

//output gap
E_LYGAP = (1-ALPHAE)*(log(E_UCAP)-log(E_UCAP0))+ALPHAE*(E_LL-E_LL0);                               //  22: EYGAP

//SS employment
E_LL0 = RHOL0*E_LL0(-1) + (1-RHOL0)*E_LL;                                                          //  23: L0

//SS capacity utilisation
E_UCAP0 = RHOUCAP0*E_UCAP0(-1) + (1-RHOUCAP0)*E_UCAP;                                              //  24: UCAP0

                                            //  25: ETAYLOR

//relative consumer price index
exp(E_LPCP) = (SE+(1-SE)*(exp(E_LPMP))^(1-SIGIME))^(1/(1-SIGIME));                                 //  26: EPCP

//import prices including markup
exp(E_LPMP) =
  (1 + E_ZEPS_ETAM + GAMIFLAG*GAMPME*( BETAE*(SFPME*E_PHIM(1)+(1-SFPME)*E_PHIM(-1)-GP0) -E_PHIM+GP0))*
   exp(E_LER)^ALPHAX;                                                                              //  27: EPX

//export prices including markup
exp(E_LPXP) =
  (1 + E_ZEPS_ETAX + GAMIFLAG*GAMPXE*( BETAE*(SFPXE*E_PHIX(1)+(1-SFPXE)*E_PHIX(-1)-GP0) -E_PHIX+GP0));   //  28: EPX

//GDP normalised to 1?
1 = exp(E_LCSN)+exp(E_LISN)+exp(E_LIGSN)+exp(E_LGSN)+E_TBYN;                                       //  29: EQ

//trade balance (share in GDP)
E_TBYN = exp(E_LEXYN)-exp(E_LIMYN) + E_ZEPS_EX;                                                    //  30: TB

//exp(E_LIMYN) = (1-SE)*(exp(E_LPCP-E_LPMP))^(SIGIME-1)*(exp(E_LCSN)+exp(E_LISN)+exp(E_LIGSN)+exp(E_LGSN));  //  31: IM

//exp(E_LEXYN) = (1-SE)*exp(E_LER*(ALPHAX*SE)-E_LPXP)^SIGEXE*exp(E_LPXP)*exp(E_LYWY)^ALPHAX;               //  32: EX

//aggregate imports
exp(E_LIMYN) = (1-SE)*(exp(RHOPCPM*(E_LPCP(-1)-E_LPMP(-1)) + (1-RHOPCPM)*(E_LPCP-E_LPMP)))^SIGIME*
     exp(E_LPMP-E_LPCP)*(exp(E_LCSN)+exp(E_LISN)+exp(E_LIGSN)+exp(E_LGSN));                        //  31: IM

//aggregate exports
exp(E_LEXYN) = (1-SE)*
    exp(RHOPWPX*(E_LER(-1)*(ALPHAX*SE)-E_LPXP(-1))+(1-RHOPWPX)*(E_LER*(ALPHAX*SE)-E_LPXP))^SIGEXE*
    exp(E_LPXP)*exp(E_LYWY)^ALPHAX;                                                                //  32: EX

//UIP condition
E_INOM = E_INOMW + E_GE(+1) - RPREME*E_BWRY + E_ZEPS_RPREME;                                       //  33: EE

//Evolution of net foreign assets (domestic currency)
E_BWRY = (1+E_INOM-E_PHI(+1)-E_GY-GPOP0)*E_BWRY(-1) + E_TBYN;                                      //  34: ECA

//Government debt
exp(E_LBGYN) = (1+E_R-E_GY-GPOP0)*exp(E_LBGYN(-1))+exp(E_LGSN) + exp(E_LIGSN) + E_TRW*exp(E_LL-E_LYWR)
     -(E_TW+SSC)*E_WS -TP*(1-E_WS) -TVAT*exp(E_LCSN)-E_TAXYN;                                     //  35: BG

//Fiscal policy rule (government consumption)
E_GG-GY0 = GSLAG*(E_GG(-1)-GY0) + GFLAG*GVECM*(E_LGSN(-1)-log(GSN)) +
    GFLAG*G1E*(E_LYGAP-E_LYGAP(-1)) + GFLAG*GEXOFLAG*E_ZEPS_G + (1-GFLAG)*(E_ZEPS_G-E_ZEPS_G(-1)); //  36: G

//Fiscal policy rule (government investment)
E_GIG-GY0-GPCPI0 = IGSLAG*(E_GIG(-1)-GY0-GPCPI0) + IGFLAG*IGVECM*(E_LIGSN(-1)-log(IGSN)) +
    IGFLAG*IG1E*(E_LYGAP-E_LYGAP(-1)) + IGFLAG*IGEXOFLAG*E_ZEPS_IG + (1-IGFLAG)*(E_ZEPS_IG-E_ZEPS_IG(-1)); // 37: IG

//growth rate of government investment minus growth rate of investment
E_GIG-E_GI = E_LIGSN-E_LISN-E_LIGSN(-1)+E_LISN(-1);                                                //  38: IGSN

//Linear transfer rule
E_TRW = TRSN + TRFLAG*(TR1E*(1-exp(E_LL)-(1-L0)))  + E_ZEPS_TR;                                    //  39: ETRS

//Lump-sum tax
E_TAXYN-E_TAXYN(-1) = BGADJ1*(exp(E_LBGYN(-1))-BGTAR)+BGADJ2*(exp(E_LBGYN)-exp(E_LBGYN(-1)));      //  40: TAX

//Labour income tax
E_TW=TW0*(1+TW1*TWFLAG*E_LYGAP);                                                                   //  41: TW

E_TRYN = E_TRW*exp(E_LL-E_LYWR) ;                                                                  //  42: TRYN

E_TRTAXYN = E_TRW*exp(E_LL-E_LYWR) - E_TAXYN;                                                      //  43: TRTAXYN

E_WSW = (1-E_TW-SSC)*E_WS;                                                                        //  44: WSW

//world interest rate
E_INOMW = (1-RII)*E_EX_INOMW + RII*E_INOMW(-1)
        + RIP*(E_PHIW(-1)-GPW0) + RIX*(E_GYW(-1)-GYW0) + E_EPS_INOMW;                              //  45 INOMW (VAR)

//world inflation
E_PHIW-GPW0  = RPI*(E_INOMW(-1) - E_EX_INOMW)
        + RPP*(E_PHIW(-1)-GPW0) + RPX*(E_GYW(-1)-GYW0) + E_EPS_PW;                                 //  46: PHIW (VAR)

//world output
E_GYW-GYW0    = RXI*(E_INOMW(-1) - E_EX_INOMW)
        + RXP*(E_PHIW(-1)-GPW0) + RXX*(E_GYW(-1)-GYW0) + RXY*(E_LYWY(-1)-LYWY0) + E_EPS_YW;        //  47: GYW (VECM)

E_LYWY-E_LYWY(-1)=E_GYW-E_GY;                                                                      //  48: YWY



// EXOGENOUS

(E_GTFP -GTFP0) = E_EPS_Y;                                                                         // 49: ETFP

//Share of overhead labour in total employment
E_LOL-LOL = RHOLOL*(E_LOL(-1)-LOL)+E_EPS_LOL;                                                      // 49b: E_LOL

E_PHIPI = GPCPI0 + E_ZEPS_PPI;                                                                     // 50: PHIPI

//consumption preference shock
E_ZEPS_C = RHOCE*E_ZEPS_C(-1)+E_EPS_C;                                                             // 51: E_ZEPS_C

//price markup shock
E_ZEPS_ETA = RHOETA*E_ZEPS_ETA(-1)+E_EPS_ETA;                                                      // 52: E_ZEPS_ETA

//import price markup shock
E_ZEPS_ETAM = RHOETAM*E_ZEPS_ETAM(-1)+E_EPS_ETAM;                                                  // 53: E_ZEPS_ETAM

//export price markup shock
E_ZEPS_ETAX = RHOETAX*E_ZEPS_ETAX(-1)+E_EPS_ETAX;                                                 // 54: E_ZEPS_ETAX

E_ZEPS_EX = RHOEXE*E_ZEPS_EX(-1)+E_EPS_EX;                                                         // 55: E_ZEPS_EX

//Fiscal Policy shock (government consumption)
E_ZEPS_G = RHOGE*E_ZEPS_G(-1) + E_EPS_G;                                                           // 56: ZEPS_G

//Fiscal Policy shock (government investment)
E_ZEPS_IG = RHOIG*E_ZEPS_IG(-1) + IGEXOFLAG*E_EPS_IG;                                              // 57: ZEPS_IG

//labour supply shock
E_ZEPS_L = RHOLE*E_ZEPS_L(-1)+E_EPS_L;                                                             // 58: E_ZEPS_L

//Monetary policy shock
E_ZEPS_M = E_EPS_M;                                                                                // 59: E_ZEPS_M

// E_ZEPS_PPI = RHOPPI1*E_ZEPS_PPI(-1)+RHOPPI2*E_ZEPS_PPI(-2)+RHOPPI3*E_ZEPS_PPI(-3)+RHOPPI4*E_ZEPS_PPI(-4)+E_EPS_PPI;  //  60: ZEPS_PPI

E_ZEPS_PPI = 0;

E_ZEPS_RPREME = RHORPE*E_ZEPS_RPREME(-1)+E_EPS_RPREME;                                             //  61: E_ZEPS_RPREME

E_ZEPS_RPREMK = RHORPK*E_ZEPS_RPREMK(-1)+E_EPS_RPREMK;                                             //  62: E_ZEPS_RPREMK

E_ZEPS_W = E_EPS_W;                                                                                //  63: E_ZEPS_W

//Shock to transfers
E_ZEPS_TR = RHOTR*E_ZEPS_TR(-1)+TREXOFLAG*E_EPS_TR;                                                //  64: ZEPS_TR

//Inflation target shock (2 percent annually)
E_ZPHIT-GP0 = 0; //INFLAGE*(E_ZPHIT(-1)-GP0); // + E_EPS_INF;                                      //  65: ZPHIT



// IDENTITIES

E_GC-E_GY+E_PHIC-E_PHI = E_LCSN-E_LCSN(-1);                                                        //  66: GC

E_GCLC-E_GY+E_PHIC-E_PHI = E_LCLCSN-E_LCLCSN(-1);                                                  //  67: GC

E_GCNLC-E_GY+E_PHIC-E_PHI = E_LCNLCSN-E_LCNLCSN(-1);                                               //  68: GC

E_GE + E_PHIW-E_PHI = E_LER-E_LER(-1);                                                             //  69: GE

E_GEX-E_GY+E_PHIX-E_PHI = E_LEXYN-E_LEXYN(-1);                                                     //  70: GEX

//growth rate of government consumption
E_GG-E_GY+E_PHIC-E_PHI = E_LGSN-E_LGSN(-1);                                                        //  71: GG

//Growth rate of investment
E_GI - E_GK(-1) = E_LIK - E_LIK(-1);                                                               //  72: GI

//growth rate of government investment
E_GIG-E_GKG(-1) = E_LIKG-E_LIKG(-1);                                                               //  73: GIG

E_GIM-E_GY+E_PHIM-E_PHI = E_LIMYN-E_LIMYN(-1);                                                     //  74: GIM

(E_GY - E_GK + E_PHIPI) = E_LYKPPI - E_LYKPPI(-1);                                                 //  75: GK

E_GL = E_LL - E_LL(-1);                                                                            //  76: GL

E_GTAX-E_GY-E_PHI = log(E_TAXYN)-log(E_TAXYN(-1));                                                 //  77: GTAX

E_GTFPUCAP = (1-ALPHAE)*E_GUCAP+ALPHAE*E_GTFP;                                                     //  78: E_GTFPUCAP

E_GTR-E_GL-E_WRPHI=log(E_TRW)-log(E_TRW(-1));                                                      //  79: GTR

E_GUCAP = log(E_UCAP)-log(E_UCAP(-1));                                                             //  80: GUCAP

E_GWRY=-E_LYWR+E_LYWR(-1);                                                                         //  81: GWRY

E_GY - E_GYPOT = E_LYGAP-E_LYGAP(-1);                                                              //  82: GYPOT

E_BGYN = exp(E_LBGYN);                                                                             //  83: BGYN

E_DBGYN = E_BGYN-E_BGYN(-1);                                                                       //  84: DBGYN

E_CLCSN = exp(E_LCLCSN);                                                                           //  85: CLCSN

E_GSN = exp(E_LGSN);                                                                               //  86: GSN


E_LTRYN = log(E_TRYN);                                                                             //  87: log(TRYN)

//Definition of CPI_Inflation
E_PHIC-E_PHI = E_LPCP-E_LPCP(-1);                                                                  //  88: PHIC

//Definition of Import price- Inflation
E_PHIM-E_PHI = E_LPMP-E_LPMP(-1);                                                                  //  89: PHIM

//Definition of Export price - inflation
E_PHIX-E_PHI = E_LPXP-E_LPXP(-1);                                                                  //  90: PHIX

-E_WPHI + E_GY + E_PHI = E_LYWR - E_LYWR(-1);                                                      //  91: YWR

E_WPHI = E_WRPHI+E_PHI;                                                                            //  92: EWPHI

E_GYL = E_GY+GPOP0;                                                                                //  93: GYL

E_GCL = E_GC+GPOP0;                                                                                //  94: GCL

E_GIL = E_GI+GPOP0;                                                                                //  95: GIL

E_GGL = E_GG+GPOP0;                                                                                //  96: GGL

E_GEXL = E_GEX+GPOP0+DGEX;                                                                         //  97: GEXL

E_GIML = E_GIM+GPOP0+DGIM;                                                                         //  98: GIML

E_PHIML = E_PHIM + DGPM;                                                                           //  99: PHIML

E_PHIXL = E_PHIX + DGPX;                                                                           //  100: PHIXL

E_LCY = E_LCSN - E_LPCP;

E_LGY = E_LGSN - E_LPCP;

E_LWS = E_LL-E_LYWR;

end;


initval;
E_BGYN    =    2.4000;
E_BWRY    =    0;
E_CLCSN   =      0.3862;
E_DBGYN   =   0;
E_LER     =   0;
E_ETA     =   0.9000;
E_GC      =     0.0030;
E_GCL     =     0.0041;
E_GCLC    =      0.0030;
E_GCNLC   =      0.0030;
E_GE      =     0;
E_GEX     =     0.0030;
E_GEXL    =      0.0115;
E_GG      =      0.0030;
E_GGL      =      0.0041;
E_GI      =       0.0030;
E_GIG     =       0.0030;
E_GIL     =     0.0041;
E_GIM      =    0.0030;
E_GIML     =   0.0115;
E_GK       =      0.0030;
E_GKG       =     0.0030;
E_GL        =  0;
E_GSN       = 0.2030;
E_GTAX      =    0.0080;
E_GTFP       =0.0024;
E_GTFPUCAP   = 0.0013;
E_GTR       =   0.0030;
E_GUC      =   0;
E_GUCAP     = 0;
E_GWRY     =  0;
E_GY       =  0.0030;
E_GYL     =       0.0041;
E_GYPOT   =    0.0030;
E_GYW     =   0.0030;
E_INOM    =         0.0090;
E_INOMW   =    0.0090;
E_LL      =   -0.4308;
E_LL0     =    -0.4308;
E_LBGYN   =     0.8755;
E_LCSN     =     -0.5381;
E_LCLCSN   =   -0.9514;
E_LCNLCSN  =     -0.3702;
E_LEXYN    =  -1.9576;
E_LGSN    =      -1.5945;
E_LIGSN   =      -3.6889;
E_LIMYN   =      -1.9576;
E_LIK      =    -3.5359;
E_LIKG    =      -4.0963;
E_LISN     =   -1.6706;
E_LOL     =0;
E_LPCP     =0  ;
E_LPMP       =0;
E_LPXP      = 0;
E_LTRYN    =  -1.7809;
E_LUCYN     =   4.9955;
E_LUCLCYN   =   3.9815;
E_LYGAP    =  0;
E_LYKPPI    =   -1.8653;
E_LYWR     =     0.3285;
E_LYWY     =  0;
E_MRY      =   0.9964;
E_PHI      =  0.0050;
E_PHIC     =   0.0050;
E_PHIPI    =  0;
E_PHIM     =   0.0050;
E_PHIML    =         0.0010;
E_PHIW     =   0.0050;
E_PHIX     =    0.0050;
E_PHIXL    =      0.0010;
E_Q        =      1.0000;
E_R        =      0.0040;
E_TAXYN    =  -0.0142;
E_TBYN     =    0;
E_TRTAXYN  =    0.1827;
E_TRW       =  0.3600;
E_TRYN     =      0.1685;
E_TW        =   0.2000;
E_UCAP      =   1.0000;
E_UCAP0    =   1.0000;
E_VL        = 19.0679;
E_VLLC      =  8.8600;
E_WPHI      =  0.0080;
E_WRPHI      = 0.0030;
E_WS         = 0.4680;
E_WSW        = 0.2808;
E_ZEPS_C     =0;
E_ZEPS_ETA   =0;
E_ZEPS_ETAM  =0;
E_ZEPS_ETAX  =0;
E_ZEPS_EX    =0;
E_ZEPS_G     =0;
E_ZEPS_IG    =0;
E_ZEPS_L     =0;
E_ZEPS_M     =0;
E_ZEPS_PPI   =0;
E_ZEPS_RPREME=0;
E_ZEPS_RPREMK=0;
E_ZEPS_TR    =0;
E_ZEPS_W     =0;
E_ZPHIT      =0.0050;
E_LCY        =-0.5381;
E_LGY        =-1.5945;
E_LWS        =-0.7593;
interest     = 0;
inflation    = 0;
outputgap    =0;
end;


// set exogenous shocks NOT jointly estimated with the DSGE model
shocks;
var E_EPS_INOMW;
stderr  0.0023;
var E_EPS_PW;
stderr  0.0029;
var E_EPS_YW;
stderr 0.0044;
var E_EPS_PPI;
stderr 0.00312216772065;

var E_EPS_C;
stderr  0.0597;
var  E_EPS_ETA;
stderr  0.1500;
var E_EPS_ETAM;
stderr 0.0202;
var E_EPS_ETAX;
stderr 0.0648;
var E_EPS_EX;
stderr  0.0044;
//var  E_EPS_G;
//stderr  0.0048;
var E_EPS_IG;
stderr 0.0056;
var E_EPS_L;
stderr 0.0283;
var  E_EPS_LOL;
stderr  0.0048;
var E_EPS_M;
stderr 0.0013;
var E_EPS_RPREME;
stderr 0.0017;
var E_EPS_RPREMK;
stderr  0.0070;
var  E_EPS_TR;
stderr  0.0022;
var E_EPS_W;
stderr 0.0437;
var E_EPS_Y;
stderr 0.0121;

var fiscal_;    // taken from E_EPS_G shock
stderr  0.0048;

end;