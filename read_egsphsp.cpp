
#include <stdio.h>
#include <string>
#include "read_egsphsp.h"

using namespace std;


int SwapBytesInt(int swap){
	int swap_temp,swaped=0;
	for(int i=0;i<4;i++){
		swap_temp=swap&255;
		swap=swap>>8;
		swaped=swaped<<8;
		swaped=swaped&(-255);
		swaped+=swap_temp;
	}
	return swaped;
}

float SwapBytesFloat(float swapf){
	int swap;
	memcpy(&swap,&swapf,sizeof(float));
	int swap_temp,swaped=0;
	for(int i=0;i<4;i++){
		swap_temp=swap&255;
		swap=swap>>8;
		swaped=swaped<<8;
		swaped=swaped&(-255);
		swaped+=swap_temp;
	}
	float swapedf;
	memcpy(&swapedf,&swaped,sizeof(int));

	return swapedf;
}


//Open the phase space file and read the header information
// input:   the file name of the phase space file
// return:  number of total particles in the phase space file
// also assign the FILE pointer phsp to the phase space file and put the header information in the phsp_header structure. 
//FILE * read_egsphsp_header(const char *egsphsp_file, phsp_header_t &phsp_header){
void read_egsphsp_header(const char *phsp_file){

	FILE *phsp;
    char MODE_RW[8];
	npphsp = (uint*)malloc(nPhsp*sizeof(uint));
	npphsp_total = 0;
	//unsigned int NPPHSP;
	//unsigned int NPHOTPHSP;
	//float EKMAXPHSP;
	//float EKMINPHSP;
	//unsigned int NINCPHSP;

	if(nPhsp>1){

		for(int iPhsp=0; iPhsp<nPhsp; iPhsp++){

			char buf[20];
			string phsp_file_name = string(phsp_file)+"_p"+itoa(iPhsp,buf,10)+".egsphsp1";
			logstr("  Phase space file #%d . . . . %s\n", iPhsp, phsp_file_name.c_str());
			fopen_s(&phsp, phsp_file_name.c_str(), "rb");
			if (!phsp) 
			{
				logstr("Could not open the egsphsp file \"%s\".\n", phsp_file_name.c_str());
				exit(-1);
			}

			// read MODE_RW of the phase space file
			//memset(MODE_RW, 0, 8);
			fread(MODE_RW, sizeof(char), 5, phsp);
			//printf("MODE_RW is: %s\n", phsp_header.MODE_RW);
			string mode_rw(MODE_RW);
#ifdef WITH_ZLAST
			if (mode_rw != "MODE2")
			{
				logstr("ERROR! The code is compiled with option of WITH_ZLAST (MODE2), however the input phase space file does not have ZLAST");
				exit(-1);
				//return -1;
			}
#endif

			// read NPPHSP of the phase space file
			fread(&npphsp[iPhsp], sizeof(unsigned int), 1, phsp);
			npphsp_total += npphsp[iPhsp];
			//logstr("    number of particles . . . %d\n", npphsp[iPhsp]);
			//printf("NPPHSP is: %d\n", phsp_header.NPPHSP);

			// read NPHOTPHSP of the phase space file
			//fread(&phsp_header.NPHOTPHSP, sizeof(unsigned int), 1, phsp);
			//printf("NPHOTPHSP is: %d\n", phsp_header.NPHOTPHSP);

			// read EKMAXPHSP of the phase space file
			//fread(&phsp_header.EKMAXPHSP, sizeof(float), 1, phsp);
			//printf("EKMAXPHSP is: %f\n", phsp_header.EKMAXPHSP);

			// read EKMINPHSP of the phase space file
			//fread(&phsp_header.EKMINPHSP, sizeof(float), 1, phsp);
			//printf("EKMINPHSP is: %f\n",phsp_header.EKMINPHSP);

			// read NINCPHSP of the phase space file
			//fread(&phsp_header.NINCPHSP, sizeof(unsigned int), 1, phsp);
			//printf("NINCPHSP is: %d\n", phsp_header.NINCPHSP);

			fclose(phsp);

		}

	}
	else{

		string phsp_file_name = string(phsp_file)+".egsphsp1";
		logstr("  Phase space file  . . . . . %s\n", phsp_file_name.c_str());
		fopen_s(&phsp, phsp_file_name.c_str(), "rb");
		if (!phsp) 
		{
			logstr("Could not open the egsphsp file \"%s\".\n", phsp_file_name.c_str());
			exit(-1);
		}

		// read MODE_RW of the phase space file
		//memset(MODE_RW, 0, 8);
		fread(MODE_RW, sizeof(char), 5, phsp);
		//printf("MODE_RW is: %s\n", phsp_header.MODE_RW);
		string mode_rw(MODE_RW);
#ifdef WITH_ZLAST
		if (mode_rw != "MODE2")
		{
			logstr("ERROR! The code is compiled with option of WITH_ZLAST (MODE2), however the input phase space file does not have ZLAST");
			exit(-1);
			//return -1;
		}
#endif

		// read NPPHSP of the phase space file
		fread(&npphsp[0], sizeof(unsigned int), 1, phsp);
		npphsp_total += npphsp[0];
		//logstr("    number of particles . . . %d\n", npphsp[0]);
		//printf("NPPHSP is: %d\n", phsp_header.NPPHSP);

		// read NPHOTPHSP of the phase space file
		//fread(&phsp_header.NPHOTPHSP, sizeof(unsigned int), 1, phsp);
		//printf("NPHOTPHSP is: %d\n", phsp_header.NPHOTPHSP);

		// read EKMAXPHSP of the phase space file
		//fread(&phsp_header.EKMAXPHSP, sizeof(float), 1, phsp);
		//printf("EKMAXPHSP is: %f\n", phsp_header.EKMAXPHSP);

		// read EKMINPHSP of the phase space file
		//fread(&phsp_header.EKMINPHSP, sizeof(float), 1, phsp);
		//printf("EKMINPHSP is: %f\n",phsp_header.EKMINPHSP);

		// read NINCPHSP of the phase space file
		//fread(&phsp_header.NINCPHSP, sizeof(unsigned int), 1, phsp);
		//printf("NINCPHSP is: %d\n", phsp_header.NINCPHSP);

		fclose(phsp);

	}

	logstr("  total particles . . . . . . %d\n", npphsp_total);

	//return phsp;
}

//read the phase space data
//  from the event "istart" and read "len" number of events
int read_egsphsp_data(phsp_particle_t * data, unsigned int istart, unsigned len){

	FILE *phsp;

	// the offset due to the header + 3 bytes (because the MODE_RW only have 5 bytes, so 3 bytes was added ) .
	int data_offset=5+4+4+4+4+4+3; //total 28 bytes

	int npread;
	if(istart>=npphsp_total){
		printf("The number of particles in phase space is smaller than %d.\n", istart);
		npread=-1;
	}
	else if(npphsp_total-istart<len){
		printf("The remaining particles are fewer than %d.\n", len);
		//printf("Read %d particles, start from %d.\n", phsp_header.NPPHSP-istart, istart);
		//fseek(phsp, istart*sizeof(phsp_particle_t) + data_offset, SEEK_SET); //SEEK_SET mease seek from the beginning of the file
		//fread(data, sizeof(phsp_particle_t), phsp_header.NPPHSP-istart, phsp);
		//npread=phsp_header.NPPHSP-istart;
		npread=-1;
	}
	else{
		// read LEN particles of the phase space file
		int iPhsp;
		for(iPhsp=0; iPhsp<nPhsp; iPhsp++){
			if(istart<npphsp[iPhsp]) break;
			istart -= npphsp[iPhsp];
		}
		if(iPhsp>=nPhsp) return -1;
		string phsp_file_name;
		char buf[20];
		int npremain = npphsp[iPhsp]-istart;
		if(npremain<len){
			//logstr("Read %d particles from phspfile #%d, start from %d.\n", npremain, iPhsp, istart);
			phsp_file_name = string(phsp_file)+"_p"+itoa(iPhsp,buf,10)+".egsphsp1";
			fopen_s(&phsp,phsp_file_name.c_str(),"rb");
			fseek(phsp, istart*sizeof(phsp_particle_t) + data_offset, SEEK_SET);
			fread(data, sizeof(phsp_particle_t), npremain, phsp);

			//logstr("Read %d particles from phspfile #%d, start from %d.\n", len-npremain, iPhsp+1, 0);
			phsp_file_name = string(phsp_file)+"_p"+itoa(iPhsp,buf,10)+".egsphsp1";
			fopen_s(&phsp,phsp_file_name.c_str(),"rb");
			fseek(phsp, data_offset, SEEK_SET);
			fread(data+npremain, sizeof(phsp_particle_t), len-npremain, phsp);
		}
		else{
			//logstr("Read %d particles from phspfile #%d, start from %d.\n",len,iPhsp,istart);
			phsp_file_name = string(phsp_file)+"_p"+itoa(iPhsp,buf,10)+".egsphsp1";
			fopen_s(&phsp,phsp_file_name.c_str(),"rb");
			fseek(phsp, istart*sizeof(phsp_particle_t) + data_offset, SEEK_SET);
			fread(data, sizeof(phsp_particle_t), len, phsp);
		}
		npread=len;
	}
	/*
	if(npread>0){
		for(int i=0;i<npread;i=i+100){
			printf("Particle %d: LATCH=%d, E=%f, X=%f, Y=%f, U=%f, V=%f, WT=%f,\n",
				i, particle[i].LATCH, particle[i].E, particle[i].X, particle[i].Y,
				particle[i].U, particle[i].V, particle[i].WT);
		}
	}
	*/
    fclose(phsp);
    return npread;
}