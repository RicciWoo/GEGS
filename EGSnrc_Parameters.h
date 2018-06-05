

#ifndef MC_PARA
#define MCPARA          //avoide double definition
#define MXMED   6       //Maximum number of media (exclude Vac.)
#define MXREG   20      //MAXIMUM NO. OF REGIONS ALLOCATED
#define MXSTACK 40      //STACK SIZE
#define MXVRT1  1000    //NO. OF REPRESENTATIVE ANGLES IN VERT1
#define FMXVRT1 1000.   //FLOATING $MXVRT1
//#define MXPWR2I 50      //SIZE OF TABLE OF INVERSE POWERS OF TWO
#define MXJREFF 200     //SIZE OF MULTIPLE SCATTERING JREFF MAP
#define MSSTEPS 16      //NO. OF MULTIPLE SCATTERING STEP SIZES
#define MXVRT2  100     //DISTRBTN OF NONOVERLAPPING PARTS OF VERT

//Following define max. no. of intervals for fit to functions
#define MXSGE   400     //GAMMA SMALL ENERGY INTERVALS
#define MXGE    200     //GAMMA MAPPED ENERGY INTERVALS
#define MXSEKE  300     //ELECTRON SMALL ENERGY INTERVALS
#define MXEKE   150     //ELECTRON MAPPED ENERGY INTERVALS
#define MXEL    10      //MAXIMUM # OF ELEMENTS IN A MEDIUM
#define MXLEKE  100     //ELECTRON ENERGY INTERVALS BELOW EKELIM
#define MXCMFP  100     //CUMULATIVE ELECTRON MEAN_FREE_PATH
#define MXRANGE 100     //ELECTRON RANGE                                                                         
#define MXBLC   20      //MOLIERE'S LOWER CASE B
#define MXRNTH  20      //RANDOM NUMBER FOR SUBMEDIAN ANGLES
#define MXRNTHI 20      //RANDOM NUMBER FOR SUPERMEDIAN ANGLES
#define MXRAYFF 100     //RAYLEIGH ATOMIC FORM FACTOR
#define RAYCDFSIZE 100  //CDF from RAYLEIGH FORM FACTOR SQUARED
#define MXSINC  1002    //ANGLE INTERVALS IN (0,5*PI/2) FOR SINE

//parameters for AUSGAB calls
//#define MXAUS 29           //CHANGE IF MORE AUSGAB CALLS ARE , may not be needed here.
//#define MXAUSM5 (MXAUS-5)  //SET THIS TO $MXAUS VALUE LESS 5"

//BREMPR--Bremsstrahlung and Pair Production data
//#define MXBREN  57
//#define MXBRXX  54
//#define MXBREL  100
//#define MXGAUSS 64
//#define MXBRES  100
//#define MXBRXS  50
//#define MXBRXS_PLUS_1 (MXBRXS+1)  //due to a fortran array defined as (0:MXBRXS)

//#define NIST_ENERGY_SCALE 1.0
//#define  NRC_PAIR_NXX  65
//#define  NRC_PAIR_NEE  84
//#define  NRC_PAIR_NX_1 (NRC_PAIR_NXX-1)
//#define  NRC_PAIR_NE_1 (NRC_PAIR_NEE-1)
//#define  CDUM_SIZE  (4*(NRC-PAIR-NXX-4)-1) //REPLACE {$cdum_size} WITH {{COMPUTE 4*($NRC-PAIR-NXX-4)-1}}; which equals 243.

//Triplet data
//#define MAX_TRIPLET 250

//Compton data--Incoherent Scattering data
#define MXTOTSH 1538  //Total number of shells for Z=1..100
#define MXMDSH  50    //Max. number of shells per medium

//EDGE data
#define MXELEMENT 100  //Number of elements
#define MXSHELL   6    //Number of shells treated
#define MXINTER   5    //$MXSHELL_1
#define MXTRANS   39   //Number of possible transitions
#define MXEDGE    16   //max. number of edges above 1 keV

//EII data
#define MAX_EII_SHELLS 40   //Maximum number of shells participating
#define N_EII_BINS     250  //Number of bins for EII x_section
#define MAX_EII_BINS   (N_EII_BINS*MAX_EII_SHELLS)

//Screened Rutherford MS data
#define MAXL_MS  63
#define MAXQ_MS  7
#define MAXU_MS  31
#define MAXL_MS_PLUS_1  (MAXL_MS+1)   //due to that some fortran array was defined as AAA(0:MAXL_MS)
#define MAXQ_MS_PLUS_1  (MAXQ_MS+1)
#define MAXU_MS_PLUS_1  (MAXU_MS+1)
#define LAMBMIN_MS  1.
#define LAMBMAX_MS  1e5
#define QMIN_MS     1e-3
#define QMAX_MS     0.5
#define MAXE_SPIN   15
#define MAXE_SPI1   (2*MAXE_SPIN+1)
#define MAXQ_SPIN   15
#define MAXU_SPIN   31
#define MAXE_SPI1_PLUS_1  (MAXE_SPI1+1)
#define MAXQ_SPIN_PLUS_1  (MAXQ_SPIN+1)
#define MAXU_SPIN_PLUS_1  (MAXU_SPIN+1)

//the commons used in each subprogram
#define ENEPS   0.0001  //difference between Ecut and end point energy for range calculation
#define EPSEMFP 1.e-5   //SMALLEST ELECTRON MFP VALUE
#define EPSGMFP 1.e-5   //SMALLEST GAMMA MFP VALUE

//transport algorithm related stuff

//the various transport algorithms, these numbers just have to be distinct
#define PRESTA_II  0
#define PRESTA__I  1
#define VMC        2    //here does not include the VMC option

//#define SKIN_DEPTH_FOR_BCA  3    //already used one sdfb in mscati!
#define EXACT_BCA_XIMAX     0.5
#define INEXACT_BCA_XIMAX   0.5  //this is not realy neccessary
#define	MAX_ELOSS    0.25
#define MAX_SMAX     1e10
#define GLOBAL_ECUT  0.
#define GLOBAL_PCUT  0.

//the NRC auxilliary get_inputs subroutine which is part of the standard NRC user-codes
//#define NMAX         100
//#define NVALUE       100
//#define STRING80     80
//#define STRING32     32
//#define STRING40     40
//#define STRING256    256
//#define MXALINP      5
//#define MXALINP_PLUS_1  MXALINP+1
//#define MAX_EII_BINS 10000

//#define RADC_NE 128         //Number of initialized energies"            
//#define RADC_NE_PLUS_1   (RADC_NE+1)
//#define RADC_NU 32         //Number of angular bins"                    
//#define RADC_NU_PLUS_1   (RADC_NU+1)
//#define RADC_NBOX 13917    //Number of cells (boxes) for all energies"  
//#define RADC_NX 8929       //Number of cell boundaries for all energies"

#define RELAX_CUTOFF 0.001    //the minimal relaxiation cut off energy
#define MXVAC 10
#define IBRDST_SIMPLE 0
#define IBRDST_KM 1
#define IPRDST_SIMPLE 0
#define IPRDST_KM 1

#define E_MAX_RR 0.8

//#define CTUnitNumber 45  //assign a uint number for the CT data
#define CTIMAX 512
#define CTJMAX 512
#define CTKMAX 270
#define   IMAX 128  //Maximum number of x cells
#define   JMAX 128  //Maximum number of y cells
#define   KMAX 128  //Maximum number of z cells
#define CT_REG CTIMAX*CTJMAX*CTKMAX
#define XYZREG IMAX*JMAX*KMAX    //"MAXIMUM NO. OF REGIONS ALLOCATED"
//The above define the largest CT data set we can read in"
//You can make the code require much less space by reducing to your"
//              local maximum needs"
//IMAX, JMAX and KMAX defined in dosxyznrc_user_macros.mortran"
//              define the number of the phantom/calculational voxels"

#endif  // for #ifndef MC_PARA