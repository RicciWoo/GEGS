

#include <stdlib.h>
#include <stdio.h>
#include <cstdarg>
#include <string>
#include <time.h>
#include <math.h>
#include <cuda_runtime.h>
#include <curand.h>
#include <curand_kernel.h>
#include <Windows.h>
#include <iostream>
#include <fstream>
#include "Polygon.h"    //added by Wu 20130124
#include "read_egsinp.h"
#include "leaf_end.h"
#include "read_egsphsp.h"


using namespace std;
//Defination related to MLC information , added by Wu 2013
#include "read_mlcinp.h"

#define ROTATE_SecJaws_XY
//#define ROTATE_MLC_XY
#define INVERSE_MLC_Bank_AB
#define ROTATE_Block_XY
#define ROTATE_Wedge_XY

//#define VERIFICATION_RUN

//EGS Simulation Parameters
#define EPSTEP 0.0005
//#define BCA_EXACT   //Boundary Crossing Algorithm : Exact
                    //  else                      : PRESTA-I
#define SKIN_DEPTH_FOR_BCA  0
//#define ESA_PII     //Electron Step Algorithm : PRESTA-II
                    //  else                  : PRESTA-I
//#define SPIN_FLAG   //Spin Effects: On
#define BAS_KM      //Brems Angular Sampling : KM
                    //  else                 : Simple
#define IPRDST 1    //Pair Angular Sampling : Simple(default)
                    //  0 => Off   , iprdst=0, EGS4 default angle selection
                    //  1 => Simple, iprdst=1, lowest order angular distribution
                    //  2 => KM    , iprdst=2, Motz, Olsen and Koch (1969) eq. 3D-2003
#define EII_FLAG    //Electron Impact Ionization: On
#define RELAX_FLAG  //Atomic Relaxations: On
#define RAYLE_FLAG  //Rayleigh Scattering: On
#define BCOMP_FLAG  //Bound Compton Scattering: On
#define PHTER_FLAG  //Photoelectron Angular Sampling: On


//#define AdjustGeomPASSMLC  // adjust the particle geometry right after it pass the MLC to take into account for the gantry roation and beam isocenters.

/*****************
 * CONFIGURATION *
 *****************/

//if you only want to simulate MLC, enable this macro
//#define SKIP_PHANTOM

//if you don't want to simualte MLC, enable this macro
//#define SKIP_MLC

//if particles start from phase space, enable this macro, now this is necessary
#define USE_PHASE_SPACE

//if you want to use energy scoring, enable this macro, now this is necessary
#define USE_ENERGY_SCORE

//if you want to calculate dose uncertainty, enable this macro
//#define GET_UNCERTAINTY

//if you want to use multiple GPUs, enable this macro
#define USE_MULTIPLE_GPU

//if you want to print particle number in stacks, enable this macro
//#define PRINT_PART_NUM

//if you want to do step count, enable this macro
//#define DO_STEP_COUNT

#define ECUT_MLC 0.7
#define PCUT_MLC 0.01

#define ECUT_PHANTOM 0.7
#define PCUT_PHANTOM 0.01

#define INFI_DIST 1e5
#define OUTMLC -99
#define PASSMLC -98

#ifndef SKIP_MLC
#define INIT_MLC_STAT OUTMLC
#else
#define INIT_MLC_STAT PASSMLC
#endif

#define PHANTOM_STARTED -97
//#define SourceToIsocenter 100

float h_SourceToIsocenter = 100;
__device__ float SourceToIsocenter;

// Use the sample program "Device Query" included with the CUDA SDK to find the warp size and 
// number of multiprocessors of the GPU you are going to use. Set the following two constants
// to those values.

// warp size
#define WARP_SIZE 32

// number of multiprocessors
#define NUM_MULTIPROC 15

// The following is the number of warps in each block. This should be large enough to ensure
// a good occupancy, but it is limited by the available registers and shared memory.
#define SIMULATION_WARPS_PER_BLOCK 4

// The following is the number of blocks that are launched for each multiprocessor. There is
// probably no reason for this to be much larger than 1. More blocks require more global memory.
#define SIMULATION_BLOCKS_PER_MULTIPROC 5

// the number of blocks used to run the simulation kernel
#define SIMULATION_NUM_BLOCKS (SIMULATION_BLOCKS_PER_MULTIPROC * NUM_MULTIPROC)
#define SIMULATION_NUM_THREADS (SIMULATION_NUM_BLOCKS*SIMULATION_WARPS_PER_BLOCK*WARP_SIZE)

// The following is the number of iterations of the outer loop per simulation kernel. A larger
// number will increase the performance of the simulation because fewer kernels launches will
// be necessary, which all have an overhead cost. However, a larger number will also increase 
// the accumulative effect of single precision rounding errors, thus potentially decreasing
// the accuracy of the simulation.
#define SIMULATION_ITERATIONS 2000
// The following is the number of particles that are launched for each thread.
#define SIMULATION_PARTICLES_PER_THREAD 50
//The total size of the MAIN STACK
#define MAX_STACK_SIZE 10000000
#define PHSP_STACK_SIZE SIMULATION_PARTICLES_PER_THREAD*SIMULATION_NUM_THREADS



/*************
 * CONSTANTS *
 *************/

#define PI                          3.1415926535F
#define ELECTRON_REST_MASS_FLOAT    0.5110034F          // MeV * c^(-2)
#define ELECTRON_REST_MASS_DOUBLE   0.5110034           // MeV * c^(-2)
#define HC_INVERSE                  80.65506856998F     // (hc)^(-1) in (Angstrom * MeV)^(-1)
#define TWICE_HC2                   0.000307444456F     // 2*(hc)^2 in (Angstrom * Mev)^2



/**********************************
 * MISCELLANEOUS TYPE DEFINITIONS *
 **********************************/

typedef unsigned char           uchar;
typedef unsigned short          ushort;
typedef unsigned int            uint;
typedef unsigned long long int  ulong;

// the different indices of a thread
typedef struct indices {
    uint b;     // index of the block in the grid
    uchar w;    // index of the warp in the block
    uchar t;    // index of the thread in the warp
    uint p;     // index of the particle on the stack
} indices;



// category names
// p primary (never scattered)
// c compton (Compton scattered once)
// r rayleigh (Rayleigh scattered once)
// m multiple (scattered more than once)
// t total (all photons)
const char categories[] = "pcrmt";

// buffer to read or write strings
#define CHARLEN	1024
char charBuffer[CHARLEN];

const char *input_file;
const char *egsphant_file;
const char *pegs_file;
const char *egsinp_file;    //added by Wu 20130124
const char *secjaws_file;
const char *phsp_file;  //added by Tong Xu 20130313
const char *mlcinfo_file;   //added by Wu 20130228
const char *block_file;
const char *wedge_file;
const char *output_factor_file;
const char *meas_data_file;
const char *ct_ramps_file;
const char *CTFileName;
FILE *logfile;



//FILE *phsp;
//phsp_header_t phsp_header;
int nPhsp;
uint *npphsp;
int npphsp_total;
//string *phsp_file_name;

//#define MAX_GPU 4
int GPUNo = 1;
float z_phsp;

struct SecJaws_t {
	float rmax;
	float zminy;
	float zmaxy;
	//float yfp;
	//float ybp;
	//float yfn;
	//float ybn;
	float zminx;
	float zmaxx;
	//float xfp;
	//float xbp;
	//float xfn;
	//float xbn;
} SecJaws;

struct VarianMLC_t {
	//bool switch_on;
	//int N_full;
    //int N_half;
    float rmax;
	float zmin;
    float zthick;
    //float x_start;
    //float x_end;
    //float x_mid_1;
    //float x_mid_2;
    //float air_gap;
    //float rad_end;
} VarianMLC;

struct Wedge_t {
	float rmax;
	float zminy;
	float zmaxy;
	//float yfp;
	//float ybp;
	//float yfn;
	//float ybn;
	float zminx;
	float zmaxx;
	//float xfp;
	//float xbp;
	//float xfn;
	//float xbn;
} Wedge;

//bool Wedge_on;
bool h_Change_xyz;
__constant__ bool Change_xyz;

//Information related to the MLC field and isocenter location within the patient
//added by Wu 20130124
#define MAXFIELDS 20
float3 h_isocenter_location;
Field_t field[MAXFIELDS];
int  N_Field;


//ulong histories_per_MU;
//ulong histories_per_cGy;
int total_MU;


//Gantry angle information , by Tong Xu 20130301
__constant__ float gantry_angle;
__constant__ float cosAngle;
__constant__ float sinAngle;
//__device__ float leaf_pos[120];
__constant__ float3 	isocenter_location;

point2D **d_endnode;
__constant__ point2D *endnode[TOTAL_MLC_LEAF/2];


/**************************
 * SIMULATION DEFINITIONS *
 **************************/


// all data of one particle
typedef struct particle_t {
    uchar   status;     // the current (or next) simulation step for this particle
	uchar   nRepeat;    // the remaining times that a primary particle has to repeat. 
    char    charge;     // charge of the particle 
    bool    process;    // bool indicating whether the particle needs to perform the current 
                        // simulation step
    uint    region;     // current region
	char    leafIndex;
    char    bank;       //the bank of MLC  (Bank A : 1, Bank B: 2, in the air: 0)
	uchar   module;
	uchar   latch;
    //ushort  latch;      // variable for tracking scatter events

    float   wt;         // statistical weight
	float   e;          // energy

    // position
    float   x;
    float   y;
    float   z;

    // direction
    float   u;
    float   v;
    float   w;
	// Added by Tong, the leafIndex of MLC that corresponding to the current location. equals OUTMLC if it is outside of MLC complex
	float dnear;

} particle_t;

typedef particle_t *stack_t;

/*
// we split up the data for each particle into 16-byte (128-bit) blocks (one uint4) so that we get
// coalesced global memory accesses
typedef struct stack_t {
    
    // 1st block
    uint4   *a;
    // consists of
    //uchar   status;     // 1 byte
    //uchar   nRepeat;    // 1 byte
    //char    charge;     // 1 byte
    //bool    process;    // 1 byte
    //float   e;          // 4 bytes
    //float   wt;         // 4 bytes
    //uint    region;     // 4 bytes

    // 2nd block
    uint4   *b;
    // consists of
    //char    leafIndex;  // 1 byte
	//char    bank;       // 1 byte
	//ushort  latch;      // 2 byte2
    //float   x;          // 4 bytes
    //float   y;          // 4 bytes
    //float   z;          // 4 bytes

    // 3rd block
    uint4   *c;
    // consists of
    //float   u;          // 4 bytes
    //float   v;          // 4 bytes
    //float   w;          // 4 bytes
    //float   dnear;      // 4 bytes . Modifyed by Wu, it is now dnear
} stack_t;
*/

#define NUM_SIMU_CAT 10

stack_t *d_stack;
__constant__ stack_t stack[NUM_SIMU_CAT];
//int h_np[NUM_SIMU_CAT];
__device__ int np[NUM_SIMU_CAT], np_max[NUM_SIMU_CAT];

__device__ bool simu_stack_is_empty;

#define max_oth_med 5
#define default_med 0
#define secjaws_med 1
#define MLC_medium  2
#define block_med   3
#define wedge_med   4

#define vac_and_oth 6
#define max_oth_reg 5
#define vacuum_reg  0
#define default_reg 1
#define secjaws_reg 2
#define MLC_region  3
#define block_reg   4
#define wedge_reg   5

//#define m_FirstMd 0
//#define m_SecJawY 1
//#define m_SecJawX 2
//#define m_VarMLCs 3
//#define m_WedgeMd 4
//#define m_BlockMd 5
//#define m_Phantom 6


enum module_status {
	m_FirstMd = 0x00,
	m_SecJawY = 0x01,
	m_SecJawX = 0x02,
	m_VarMLCs = 0x03,
	m_WedgeMd = 0x04,
	m_BlockMd = 0x05,
	m_Phantom = 0x06
};


#define MAX_MODULE_PLUS_1 7
char h_module_order[MAX_MODULE_PLUS_1];
__device__ char module_order[MAX_MODULE_PLUS_1];
bool SecJaws_on;
bool Wedge_on;
bool Block_on;
bool VarMLC_on;


enum particle_status {
    p_photon_step       = 0x00,
    p_compton           = 0x01,
	p_photo             = 0x02,
    p_rayleigh          = 0x03,
    p_pair              = 0x04,

    e_electron_step     = 0x05,
	e_moller            = 0x06,
    e_brems             = 0x07,
	e_bhabha            = 0x08,
	e_annih             = 0x09,

	p_new_particle      = 0x0A,
	e_new_particle      = 0x0B,

	p_empty             = 0x0C,
    e_empty             = 0x0D,

    p_cutoff_discard    = 0x0E,
    p_user_discard      = 0x0F,
    e_cutoff_discard    = 0x10,
    e_user_discard      = 0x11
};

// number of different particle statuses
#define NUM_CAT 12
// number of different detector categories (primary, compton, rayleigh, multiple)
#define NUM_DETECTOR_CAT 4

// number of blocks of the summing kernel
#define SUM_DETECTOR_NUM_BLOCKS (2 * NUM_DETECTOR_CAT)
// warps per block of the summing kernel
#define SUM_DETECTOR_WARPS_PER_BLOCK 32

/*
__shared__ uint step_counters_shared[SIMULATION_WARPS_PER_BLOCK][NUM_CAT];
typedef ulong total_step_counts_t[SIMULATION_NUM_BLOCKS][NUM_CAT];
total_step_counts_t **d_total_step_counts, *h_total_step_counts;
__constant__ total_step_counts_t *total_step_counts;
*/

float *voxel_mass;

typedef double score_t;
score_t **d_eng_score;
__constant__ score_t *eng_score;

//#ifdef GET_UNCERTAINTY
score_t *eng_dep;
score_t *eng_dep2;
//#endif


/************************
 * GEOMETRY DEFINITIONS *
 ************************/

/*
// source
typedef struct source_t {
	float   energy;
    float3  source_point;
    float   rectangle_z;
    float2  rectangle_min;
    float2  rectangle_max;
    float2  rectangle_size;
    float   rectangle_area;
} source_t; 

source_t    h_source;
__constant__    source_t    source;
*/

/*
typedef struct detector_t {
    float3  center;
    float2  d;
    uint2   N;
} detector_t;

detector_t  h_detector;
__constant__    detector_t  detector;
*/

bool create_on;
float3 ctsubmin,ctsubmax;
float3 xyz_thickness;

typedef struct phantom_t {
    uint3   N;
    float   *x_bounds;
    float   *y_bounds;
    float   *z_bounds;
} phantom_t;

//phantom_t h_phantom;
//__device__ phantom_t phantom;
__constant__ uint3 phantom_N;
__constant__ float *phantom_x_bounds;
__constant__ float *phantom_y_bounds;
__constant__ float *phantom_z_bounds;

uint3 h_phantom_N;
float *x_bounds, *y_bounds, *z_bounds;
float **d_phantom_x_bounds;
float **d_phantom_y_bounds;
float **d_phantom_z_bounds;

uint h_size_phantom;
__constant__ uint size_phantom;


typedef struct interp_t {
	float3 start;
	int3  dimen;
	float3 sizes;
} interp_t;

bool  calibr_on;
float dose_unit=1.6E-8;  //the multiple coefficient to change dose unit from MeV/g to cGy

bool  outfac_on;

bool  interp_on;
interp_t interp;

bool  AbsDos_on;
float3 dose_pos;
float dose_meas;

bool  RelDos_on;
interp_t RelDos;
char array_type;
float GmIdx_len;
float GmIdx_dos;
char  norm_type;


/*******************************
 * SIMULATION DATA DEFINITIONS *
 *******************************/

#define BOUND_COMPTON_MASK 0x000EU
#define VACUUM 0xFFFFU
#define VACUUM_STEP 1E8F
//#define EPSGMFP 1E-5F
#define SMALL_POLAR_ANGLE_THRESHOLD 1E-20F

enum region_flags {
    f_rayleigh              = 0x0001U,                  // 0000 0000 0000 0001
    f_bound_compton         = 0x0002U,                  // 0000 0000 0000 0010
    f_bound_compton_2       = 0x0006U,                  // 0000 0000 0000 0110
    f_bound_compton_3       = 0x000AU,                  // 0000 0000 0000 1010
    f_bound_compton_4       = 0x000EU,                  // 0000 0000 0000 1110
    f_atomic_relaxation     = 0x0010U,                  // 0000 0000 0001 0000
    f_photo_electron_angular_distribution = 0x0020U,    // 0000 0000 0010 0000
    f_range_rejection       = 0x0040U                   // 0000 0000 0100 0000
};

//typedef struct __align__(16) region_data_t {
typedef struct __align__(16) region_data_t {
    ushort  med;
    ushort  flags;
    float   rhof;
	//float   rho;
    float   pcut;
    float   ecut;
} region_data_t;

region_data_t *h_region_data;
region_data_t **d_region_data;
__constant__ region_data_t *region_data;

uint total_number_regions;

/*
//we now try to use texture menory to store the region_data
cudaExtent volumeSize;
cudaMemcpy3DParms copyParams = {0};
int2 *med_flags_h;
cudaArray *med_flags;
texture<int2, 3, cudaReadModeElementType> med_flags_tex;

float4 * rhof_rho_pcut_ecut_h;
cudaArray *rhof_rho_pcut_ecut;
texture<float4, 3, cudaReadModeElementType> rhof_rho_pcut_ecut_tex;
*/

#define MAXREG 5   // the maximum number of regions, can be adjusted by user. Tong Xu Jan 2013

// __constant__ region_data_t region_data[MAXREG];




/*****************************    //added by Wu 20130124
 * GLOBAL Geometry Variables *
 *****************************/

__constant__  float MLC_zmin, MLC_zthick, MLC_x_start, MLC_x_end, MLC_x_mid_1,
	 MLC_x_mid_2, MLC_air_gap, MLC_wl_full, MLC_wl_target, MLC_wl_isocenter, MLC_wt,
	 MLC_rmax, MLC_rad_end;

//__constant__ point2D boundaryXZ[4];  //define the pointer to first address of boundaryXZ on device constant memory
//__constant__ point2D boundaryYZ[4];  //define the pointer to first address of boundaryYZ on device constant memory
//__constant__ point2D crossXZ[960];
//__constant__ point2D crossYZ[34];

__constant__ point2D crossXZ[960];
__constant__ point2D boundaryXZ[4];  //define the pointer to first address of boundaryXZ on device constant memory
__constant__ point2D boundaryYZ[4];  //define the pointer to first address of boundaryYZ on device constant memory
//__device__ point2D crossYZ[34];  //this is for function leaf_position
point2D h_crossYZ[34];
point2D **d_crossYZ;

__constant__ point2D boundMinY[4];
__constant__ point2D boundMinX[4];
__constant__ point2D boundMaxY[4];
__constant__ point2D boundMaxX[4];
__constant__ point2D jawsGeomY[8];
__constant__ point2D jawsGeomX[8];

__device__ char nVer_Wed;
__constant__ point2D boundWed[4];
__constant__ point2D goem_Wed[6];
__constant__ point2D blockXZ[4];
__constant__ point2D blockYZ[4];

/**************************
 * MEDIA DATA DEFINITIONS *
 **************************/

const char *data_dir;
const char *output_dir;  //added by Wu, 20151201
//const char *dicom_dir;   //added by Wu, 20151201  //for relative path in dicom file list
const char *photon_xsections;
const char *atomic_ff_file;

//clock_t tic, tic2,tac;

//For CUDA Random number generator
__constant__ curandState *devStates;
curandState **d_devStates;

/*
#define DUMP_INT(P1,P2)    iDump=fopen("dump_"#P1"\\"#P2".txt","w");fprintf(iDump,"%8d\n",P2);fprintf(iSum,"%8d  is the value for "#P2"\n",P2);fclose(iDump);
#define DUMP_REAL(P1,P2)    iDump=fopen("dump_"#P1"\\"#P2".txt","w");fprintf(iDump,"%15e\n",P2);;fprintf(iSum,"%15e  is the value for "#P2"\n",P2);fclose(iDump);
#define DUMP_1D_INT(P1,P2,P3)    iDump=fopen("dump_"#P1"\\"#P2".txt","w");sum=0.0;for(int i=0;i<P3;i++){fprintf(iDump,"%8d\n",P2[i]);sum+=P2[i];}fprintf(iDump,"sum=  %15e\n",sum);fprintf(iSum,"%15e  is the sum for "#P2"\n",sum);fclose(iDump);
#define DUMP_1D_REAL(P1,P2,P3)    iDump=fopen("dump_"#P1"\\"#P2".txt","w");sum=0.0;for(int i=0;i<P3;i++){fprintf(iDump,"%15e\n",P2[i]);sum+=P2[i];}fprintf(iDump,"sum=  %15e\n",sum);fprintf(iSum,"%15e  is the sum for "#P2"\n",sum);fclose(iDump);
#define DUMP_1D_FLOAT2(P1,P2,P3)    iDump=fopen("dump_"#P1"\\"#P2"_x.txt","w");sum=0.0;for(int i=0;i<P3;i++){fprintf(iDump,"%15e\n",P2[i].x);sum+=P2[i].x;}fprintf(iDump,"sum=  %15e\n",sum);fprintf(iSum,"%15e  is the sum for "#P2".x\n",sum);fclose(iDump);\
	iDump=fopen("dump_"#P1"\\"#P2"_y.txt","w");sum=0.0;for(int i=0;i<P3;i++){fprintf(iDump,"%15e\n",P2[i].y);sum+=P2[i].y;}fprintf(iDump,"sum=  %15e\n",sum);fprintf(iSum,"%15e  is the sum for "#P2".y\n",sum);fclose(iDump);
*/