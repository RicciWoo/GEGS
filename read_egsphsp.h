
#ifndef READ_EGSPHSP_H
#define READ_EGSPHSP_H

//#define WITH_ZLAST

struct  phsp_header_t {
    char MODE_RW[8];
	unsigned int NPPHSP;
	unsigned int NPHOTPHSP;
	float EKMAXPHSP;
	float EKMINPHSP;
	unsigned int NINCPHSP;
};

// total length of the array
//#define LEN SIMULATION_PARTICLES_PER_THREAD*WARP_SIZE*SIMULATION_WARPS_PER_BLOCK*SIMULATION_BLOCKS_PER_MULTIPROC*NUM_MULTIPROC

struct phsp_particle_t{
	unsigned int LATCH;
	float E,X,Y,U,V,WT;
#ifdef WITH_ZLAST
	float ZLAST;
#endif
};

//int read_egsphsp_data(FILE *phsp, phsp_header_t &phsp_header, phsp_particle_t * data, unsigned int istart, unsigned len);
//FILE * read_egsphsp_header(const char *egsphsp_file, phsp_header_t &phsp_header);

#endif