

#define CUDA_EGS

#include "CUDA_EGS.h"
#include "EGSnrc_Parameters.h"
#include "commonBlocks host.h"
#include "commonBlocks dev.h"
#include "output.c"
#include "media.cu"
#include "PhysicsInit.cu"
#include "copy.cu"
//#include "dump_phys.cu"
#include "read_egsphsp.cpp"
#include "ReadCT_DICOM.c"
#include "ctcreate.c"
#include "init.cu"
#include "geom.cu"
#include "utility.cu"
#include "scoring.cu"
#include "annihatrest.cu"
#include "pair.cu"
#include "relax.cu"
#include "eii_sample.cu"
#include "moller.cu"
#include "compton.cu"
#include "photo.cu"
#include "rayleigh.cu"
#include "photon.cu"
//#include "PrepareProcessForVerfication.h"
#include "brems.cu"
//#include "EndOfVerification.h"
#include "bhabha.cu"
#include "annih.cu"
#include "electron.cu"
//#include "read_verfication.cu"

//#include "Geometry.cu"    //added by Wu 20120124

int load_stack_from_phspfile(uint istart, int stack_size, uchar nRepeat)
{
	phsp_particle_t *particles = (phsp_particle_t *)malloc(stack_size * sizeof(phsp_particle_t));

	//stack_size = read_egsphsp_data(phsp, phsp_header, particles,istart,stack_size);
	stack_size = read_egsphsp_data(particles,istart,stack_size);
	if(stack_size <0) return -1;
	
	//stack_t h_stack;
	//h_stack.a=(uint4*)malloc(stack_size * sizeof(uint4));
	//h_stack.b=(uint4*)malloc(stack_size * sizeof(uint4)); 
	//h_stack.c=(uint4*)malloc(stack_size * sizeof(uint4));
	particle_t *h_stack =(particle_t*)malloc(stack_size*sizeof(particle_t));

	//float output_factor = field[iField].output_factor;

	//uint4 tmp;
	char  q;
	//char  leafIndex;
	float w; //The Z-direction cosine
	// write particle to stack
	for(int i = 0; i< stack_size; i++)
	{
		if( particles[i].LATCH & 0x40000000 )  q = -1;   //latch bit 29/30 indicate the charge of the particle
		else if( particles[i].LATCH & 0x20000000 ) q = 1; 
		else q = 0;
		w=sqrt(1-particles[i].U*particles[i].U - particles[i].V*particles[i].V);
		if(particles[i].WT<0) {w=-w, particles[i].WT *= -1;}  //the weight WT also carries the sign of the SZ-direction cosine

		/*
		uchar status;
		if(q==0) status = p_new_particle;  //start as p_new_paritcle to distinguish from secondary particle
		else status = e_new_particle;

		((uchar*)&tmp.x)[0] = status;
		// was called "reserved". Now, used to store the MLC Bank that the particle is in.
		((uchar*)&tmp.x)[1] = nRepeat;  //p.nRepeat
		((char*)&tmp.x)[2] = q; // p.charge
		((uchar*)&tmp.x)[3] = 1;// p.process;
		if(particles[i].E<0) particles[i].E=-particles[i].E;
		tmp.y = *(uint*)&particles[i].E;
		tmp.z = *(uint*)&particles[i].WT;
		tmp.w = 0;// p.region
		h_stack.a[i] = tmp;

		leafIndex = INIT_MLC_STAT;
		((char*)&tmp.x)[0] = leafIndex;  // p.leafIndex;  modified by Wu, 20130324, note that leafIndex can be negative , but it will still be safely transfered. 
		((char*)&tmp.x)[1] = -1;   // p.bank;  start with -1, means to be determined
		((ushort*)&tmp.x)[1] = 0;  // reset p.latch
		tmp.y = *(uint*)&particles[i].X;
		tmp.z = *(uint*)&particles[i].Y;
		tmp.w = *(uint*)&h_source.source_point.z; //source_point is z location of the source
		tmp.w = *(uint*)&z_phsp;
		h_stack.b[i] = tmp;

		tmp.x = *(uint*)&particles[i].U;
		tmp.y = *(uint*)&particles[i].V;
		tmp.z = *(uint*)&w;
		float dnear = INFI_DIST;
		tmp.w = *(uint*)&dnear; //p.dnear set to infinity to start with.
		h_stack.c[i] = tmp;
		*/

		//uchar status;
		if(q==0) h_stack[i].status = p_new_particle;
		else     h_stack[i].status = e_new_particle;

		// was called "reserved". Now, used to store the MLC Bank that the particle is in.
		h_stack[i].nRepeat = nRepeat;  //p.nRepeat
		h_stack[i].charge = q; // p.charge
		h_stack[i].process = 1;// p.process;
		if(particles[i].E<0) particles[i].E = -particles[i].E;
		h_stack[i].e = particles[i].E;
		h_stack[i].wt = particles[i].WT;
		//h_stack[i].wt = particles[i].WT*output_factor;
		h_stack[i].region = 0;// p.region

		//leafIndex = INIT_MLC_STAT;
		h_stack[i].leafIndex = INIT_MLC_STAT;  // p.leafIndex;  modified by Wu, 20130324, note that leafIndex can be negative , but it will still be safely transfered. 
		h_stack[i].bank = -1;   // p.bank;  start with -1, means to be determined
		h_stack[i].module = h_module_order[m_FirstMd];
		//h_stack[i].module = m_SecJawY;
		//h_stack[i].module = m_VarMLCs;
		//h_stack[i].module = m_Phantom;
		h_stack[i].latch = 0;  // reset p.latch
		h_stack[i].x = particles[i].X;
		h_stack[i].y = particles[i].Y;
		//h_stack[i].z = h_source.source_point.z; //source_point is z location of the source
		h_stack[i].z = z_phsp;

		h_stack[i].u = particles[i].U;
		h_stack[i].v = particles[i].V;
		h_stack[i].w = w;
		//float dnear = INFI_DIST;
		h_stack[i].dnear = INFI_DIST; //p.dnear set to infinity to start with.
	}

	//Put all the primary particles in the photon stack. Of cause there will be some electrons in the phase space file
	// we are going to move them into electron stack during first running of photon simulation kernel.
	cudaMemcpy(d_stack[p_photon_step],h_stack,stack_size*sizeof(particle_t),cudaMemcpyHostToDevice); ce(4021);
    //cudaMemcpy(d_stack[p_photon_step].a,h_stack.a,stack_size*sizeof(uint4),cudaMemcpyHostToDevice); ce(4021);
    //cudaMemcpy(d_stack[p_photon_step].b,h_stack.b,stack_size*sizeof(uint4),cudaMemcpyHostToDevice); ce(4022);
    //cudaMemcpy(d_stack[p_photon_step].c,h_stack.c,stack_size*sizeof(uint4),cudaMemcpyHostToDevice); ce(4023);

	//np_phsp=stack_size-1;
	//h_np[p_photon_step] = stack_size-1;
	int h_np_phsp = stack_size-1;
	cudaMemcpyToSymbol(np,&h_np_phsp,sizeof(int)); ce(40005);
	//cudaMemcpyToSymbol(np,h_np,NUM_SIMU_CAT*sizeof(int)); ce(40005);
	
	free(particles);
	free(h_stack);
	//free(h_stack.a);
	//free(h_stack.b);
	//free(h_stack.c);

	return stack_size;

}

/*
void launch_simu_kernel(int N_simu, int N_photon, int N_electr, float &time_sim, float &time_sum){
	if(N_simu<0||N_simu>=4) return;
	if(N_photon<0||N_photon>=4) return;
	if(N_electr<0||N_electr>=4) return;
	if(np_stack[N_simu]<0) {
		logstr("stack %d is empty, skip this kernel.\n", N_simu);
		return;
	}

	cudaEvent_t start, stop; 
    cudaEventCreate(&start); ce(11001);
    cudaEventCreate(&stop); ce(11002);

	float elapsed;

	cudaMemcpyToSymbol(stack_simu, &d_stack[N_simu], sizeof(stack_t)); ce(4024);
	cudaMemcpyToSymbol(np_simu, &np_stack[N_simu], sizeof(int)); ce(4025);
	cudaMemcpyToSymbol(stack_photon, &d_stack[N_photon], sizeof(stack_t)); ce(4026);
	cudaMemcpyToSymbol(np_photon, &np_stack[N_photon], sizeof(int)); ce(4027);
	cudaMemcpyToSymbol(stack_electr, &d_stack[N_electr], sizeof(stack_t)); ce(4028);
	cudaMemcpyToSymbol(np_electr, &np_stack[N_electr], sizeof(int)); ce(4029);
	int NP_Simu = np_stack[N_simu]+1;

	cudaEventRecord(start); ce(11003);

	if(N_simu%2==0)
		photon_kernel<<<dim3(SIMULATION_BLOCKS_PER_MULTIPROC,NUM_MULTIPROC), SIMULATION_WARPS_PER_BLOCK*WARP_SIZE>>>();
	else
		electr_kernel<<<dim3(SIMULATION_BLOCKS_PER_MULTIPROC,NUM_MULTIPROC), SIMULATION_WARPS_PER_BLOCK*WARP_SIZE>>>();
	
	cudaEventRecord(stop); ce(11004);
    cudaEventSynchronize(stop); ce(11005);
    cudaEventElapsedTime(&elapsed, start, stop); ce(11006);
    time_sim += elapsed;

	cudaMemcpyFromSymbol(&np_stack[N_simu],np_simu,sizeof(int)); ce(4030);
	cudaMemcpyFromSymbol(&np_stack[N_photon],np_photon,sizeof(int)); ce(4031);
	cudaMemcpyFromSymbol(&np_stack[N_electr],np_electr,sizeof(int)); ce(4032);

	logstr("kernel(%d,%d,%d),np=%d,%d,%d,%d, time=%.2f, Hs/ms: %.2f\n",N_simu,N_photon,N_electr,np_stack[0],np_stack[1],np_stack[2],np_stack[3],elapsed,NP_Simu/elapsed);

    cudaEventDestroy(start); ce(11013);
    cudaEventDestroy(stop); ce(11014);
}
*/

/*
uint read_step_counts(ulong *this_total, ulong *grand_total, int GPUId) {
    clock_t start = clock();
    cudaMemcpy(h_total_step_counts, d_total_step_counts[GPUId], sizeof(total_step_counts_t), cudaMemcpyDeviceToHost); ce(12001);
    
    for (uchar i = 0; i < NUM_CAT; i++) {
        this_total[i] = 0;
        for (uchar j = 0; j < SIMULATION_NUM_BLOCKS; j++)
            this_total[i] += (*h_total_step_counts)[j][i];
            
        grand_total[i] += this_total[i];
    }

    clock_t stop = clock();

    return stop - start;
}
*/

bool launch_kernels(uchar status, int *h_np, float *time_sim){

	cudaEvent_t  *start =(cudaEvent_t*) malloc(GPUNo*sizeof(cudaEvent_t));
	cudaEvent_t  *stop  =(cudaEvent_t*) malloc(GPUNo*sizeof(cudaEvent_t));
	cudaStream_t *stream=(cudaStream_t*)malloc(GPUNo*sizeof(cudaStream_t));
	for(int GPUId=0; GPUId<GPUNo; GPUId++){
#ifdef USE_MULTIPLE_GPU
		cudaSetDevice(GPUId); ce(58020);
#endif
		cudaEventCreate (&start[GPUId] ); ce(11001);
		cudaEventCreate (&stop[GPUId]  ); ce(11002);
		cudaStreamCreate(&stream[GPUId]); ce(11003);
	}

	int *Number=(int*)malloc(GPUNo*sizeof(int));
	float *elapsed=(float*)malloc(GPUNo*sizeof(float));
	bool all_stack_empty = true;

	for(int GPUId=0; GPUId<GPUNo; GPUId++){
#ifdef USE_MULTIPLE_GPU
		cudaSetDevice(GPUId); ce(58022);
#endif
		Number[GPUId] = h_np[GPUId*NUM_SIMU_CAT+p_photon_step]+1;
		cudaEventRecord(start[GPUId],stream[GPUId]); ce(11004);
		switch(status){
		case p_photon_step:
			photon_step<<<dim3(SIMULATION_BLOCKS_PER_MULTIPROC,NUM_MULTIPROC), SIMULATION_WARPS_PER_BLOCK*WARP_SIZE, 0, stream[GPUId]>>>();
			break;
		case p_compton:
			compton<<<dim3(SIMULATION_BLOCKS_PER_MULTIPROC,NUM_MULTIPROC), SIMULATION_WARPS_PER_BLOCK*WARP_SIZE, 0, stream[GPUId]>>>();
			break;
		case p_photo:
			photo<<<dim3(SIMULATION_BLOCKS_PER_MULTIPROC,NUM_MULTIPROC), SIMULATION_WARPS_PER_BLOCK*WARP_SIZE, 0, stream[GPUId]>>>();
			break;
		case p_rayleigh:
			rayleigh<<<dim3(SIMULATION_BLOCKS_PER_MULTIPROC,NUM_MULTIPROC), SIMULATION_WARPS_PER_BLOCK*WARP_SIZE, 0, stream[GPUId]>>>();
			break;
		case p_pair:
			pair_production<<<dim3(SIMULATION_BLOCKS_PER_MULTIPROC,NUM_MULTIPROC), SIMULATION_WARPS_PER_BLOCK*WARP_SIZE, 0, stream[GPUId]>>>();
			break;
		case e_electron_step:
			electr_step<<<dim3(SIMULATION_BLOCKS_PER_MULTIPROC,NUM_MULTIPROC), SIMULATION_WARPS_PER_BLOCK*WARP_SIZE, 0, stream[GPUId]>>>();
			break;
		case e_moller:
			moller<<<dim3(SIMULATION_BLOCKS_PER_MULTIPROC,NUM_MULTIPROC), SIMULATION_WARPS_PER_BLOCK*WARP_SIZE, 0, stream[GPUId]>>>();
			break;
		case e_brems:
			brems<<<dim3(SIMULATION_BLOCKS_PER_MULTIPROC,NUM_MULTIPROC), SIMULATION_WARPS_PER_BLOCK*WARP_SIZE, 0, stream[GPUId]>>>();
			break;
		case e_bhabha:
			bhabha<<<dim3(SIMULATION_BLOCKS_PER_MULTIPROC,NUM_MULTIPROC), SIMULATION_WARPS_PER_BLOCK*WARP_SIZE, 0, stream[GPUId]>>>();
			break;
		case e_annih:
			annih<<<dim3(SIMULATION_BLOCKS_PER_MULTIPROC,NUM_MULTIPROC), SIMULATION_WARPS_PER_BLOCK*WARP_SIZE, 0, stream[GPUId]>>>();
			break;
		default:
			logstr("Unknown kernel number %d, no kernel is launched!!",status);
			break;
		}
		cudaEventRecord(stop[GPUId]); ce(11005);
	}
	for(int GPUId=0; GPUId<GPUNo; GPUId++){
#ifdef USE_MULTIPLE_GPU
		cudaSetDevice(GPUId); ce(58023);
#endif
		//cudaEventSynchronize(stop); ce(11006);
		cudaStreamSynchronize(stream[GPUId]); ce(11006);
		cudaEventElapsedTime(&elapsed[GPUId], start[GPUId], stop[GPUId]); ce(11007);
		time_sim[GPUId] += elapsed[GPUId];
		//time_copy += read_step_counts(this_total, grand_total, GPUId);
		//all_stack_empty = true;
#ifdef PRINT_PART_NUM
		logstr("ID%d:",GPUId);
		switch(status){
		case p_photon_step:
			logstr("photon:");
			break;
		case p_compton:
			logstr("compt :");
			break;
		case p_photo:
			logstr("photo :");
			break;
		case p_rayleigh:
			logstr("raylei:");
			break;
		case p_pair:
			logstr("pair  :");
			break;
		case e_electron_step:
			logstr("electr:");
			break;
		case e_moller:
			logstr("moller:");
			break;
		case e_brems:
			logstr("brems :");
			break;
		case e_bhabha:
			logstr("bhabha:");
			break;
		case e_annih:
			logstr("annih :");
			break;
		default:
			logstr("Unknown kernel number %d, no kernel is launched!!",status);
			exit(0);
		}
#endif
		cudaMemcpyFromSymbol(&h_np[GPUId*NUM_SIMU_CAT],np,NUM_SIMU_CAT*sizeof(int)); ce(11008);
		for(uchar i=0;i<NUM_SIMU_CAT;i++){
#ifdef PRINT_PART_NUM
			if(i==p_photon_step||i==p_compton||i==e_electron_step)
				logstr("%8d,",h_np[GPUId*NUM_SIMU_CAT+i]);
			else if(i==p_photo||i==p_rayleigh||i==e_moller)
				logstr("%7d,",h_np[GPUId*NUM_SIMU_CAT+i]);
			else
				logstr("%6d,",h_np[GPUId*NUM_SIMU_CAT+i]);
#endif
			all_stack_empty &= (h_np[GPUId*NUM_SIMU_CAT+i]<0);
		}
#ifdef PRINT_PART_NUM
		logstr("%9.2f\n",Number[GPUId]/elapsed[GPUId]);
#endif
	}

	for(int GPUId=0; GPUId<GPUNo; GPUId++){
#ifdef USE_MULTIPLE_GPU
		cudaSetDevice(GPUId); ce(58048);
#endif
		cudaEventDestroy ( start[GPUId]); ce(11069);
		cudaEventDestroy (  stop[GPUId]); ce(11070);
		cudaStreamDestroy(stream[GPUId]); ce(11071);
	}
	free(start);
	free(stop);
	free(stream);

	return all_stack_empty;
}

void runSimulation(ulong histories_per_MU, uint iField)
{
	
	// record start time
    clock_t run_tic = clock();
	//clock_t tic2 = clock();

	getTime();
	logstr("Field #%d started simulation at %s\n", iField,charBuffer);

	ulong num_histories = (ulong)(field[iField].MU*histories_per_MU);
	uint batch_size = PHSP_STACK_SIZE;

	//"repeat" means "in addition to the first time the particle being simulated"
	//int NRepeat = (num_histories-1)/phsp_header.NPPHSP;  //how many times each phase space particle need to repeat to get the total number of histories.
	int NRepeat = (num_histories-1)/npphsp_total;
	if(NRepeat<0||NRepeat>255){
		logstr("The repeat number is not acceptable, try another number of hitories or phase space file\n");
		return;
	}

	uchar nRepeat = NRepeat;
	//logstr("The %d phase space particles will be repeated by %d times to get %I64u histories\n", npphsp_total, nRepeat, num_histories);
	logstr("  History number  . . . . . . %I64u\n", num_histories);
	logstr("  Repeat times  . . . . . . . %d\n", nRepeat+1);
	//the phase_space particles will be simualted by batches

	batch_size = init_stack(nRepeat);
	// the batch_size is initially set to PHSP_STACK_SIZE, but init_stack can adjust it if biggest stack size if bigger than MAX_STACK_SIZE

	//int   NBatch = (num_histories-1)/(PHSP_STACK_SIZE*(nRepeat+1));  //here should be PHSP_STACK_SIZE
	int   nBatch = (num_histories-1)/(batch_size*(nRepeat+1));
	/*
	if(NBatch<0||NBatch>255){
		logstr("The number of batch is not acceptable, try another number of hitories or phase space file\n");
		return;
	}
	*/
	//uchar nBatch = NBatch;
	//uint  lastBatch = (uint)(num_histories-nBatch*(PHSP_STACK_SIZE*(nRepeat+1)))/(nRepeat+1)+1;
	uint  lastBatch = (uint)(num_histories-nBatch*(batch_size*(nRepeat+1)))/(nRepeat+1)+1;
    //logstr("These phase space particles will be devided into %d batches,\n",nBatch+1);
    //logstr("%d batches of %d particles and last batch of %d particles.\n",nBatch,batch_size,lastBatch);
	logstr("  Batch number  . . . . . . . %d\n", nBatch+1);
	logstr("  Batch size  . . . . . . . . %d\n", batch_size);

	//transfer the MLC leaf position information and gantry angle to GPU.
	//cudaMemcpyToSymbol(leaf_pos, &field[iField].pos, TOTAL_MLC_LEAF* sizeof(float)); ce(20001);
	setup_secjaws(iField);
	init_leaf_pos(iField);

	logstr("  Monitor unit  . . . . . . . %d\n", field[iField].MU);
	logstr("  Gantry angle  . . . . . . . %.1f\n", field[iField].angle);
	field[iField].angle = field[iField].angle*3.1415926/180;
	float h_cosAngle = cos(field[iField].angle);
	float h_sinAngle = sin(field[iField].angle);
	if(outfac_on)
		logstr("  Output factor . . . . . . . %f\n", field[iField].output_factor);

	for(int GPUId=0; GPUId<GPUNo; GPUId++) {
#ifdef USE_MULTIPLE_GPU
		cudaSetDevice(GPUId); ce(58010);
#endif
		cudaMemcpyToSymbol(gantry_angle, &field[iField].angle, sizeof(float)); ce(20002);
		cudaMemcpyToSymbol(cosAngle, &h_cosAngle, sizeof(float)); ce(20003);
		cudaMemcpyToSymbol(sinAngle, &h_sinAngle, sizeof(float)); ce(20004);

	}

	//int h_np[NUM_SIMU_CAT];
	int *h_np=(int*)malloc(GPUNo*NUM_SIMU_CAT*sizeof(int));
	memset(h_np,-1,GPUNo*NUM_SIMU_CAT*sizeof(int));

    float *time_sim=(float*)malloc(GPUNo*sizeof(float));
	memset(time_sim,0,GPUNo*sizeof(float));

	bool all_stack_empty = false;

#ifdef PRINT_PART_NUM
	logstr("GPU:kernel:  photon, compton,  photo, raylei,  pair,electron, moller, brems,bhabha, annih,    Hs/ms\n");
#else
	float complete;
#endif

        bool load_stack=true;
		int  stack_size=0;
        uint iBatch=0;
        uint istart=0;

		//bool h_simu_stack_is_empty = false;

#ifdef GET_UNCERTAINTY
		for( ; load_stack; ){
#endif
			
			//load a batch of particles when GET_UNCERTAINTY or else only load first batch of particles
			for(int GPUId=0; GPUId<GPUNo; GPUId++){
#ifdef USE_MULTIPLE_GPU
				cudaSetDevice(GPUId); ce(58021);
#endif
				istart = iBatch*batch_size;
				if(iBatch<nBatch){
					stack_size = load_stack_from_phspfile(istart,batch_size,nRepeat);
					//logstr("Bat%3d:%8d\n",iBatch,h_np[p_photon_step]);
				}
				else if(iBatch==nBatch){
					stack_size = load_stack_from_phspfile(istart,lastBatch,nRepeat);
					load_stack=false;
					//logstr("Bat%3d:%8d\n",iBatch,h_np[p_photon_step]);
				}
				h_np[GPUId*NUM_SIMU_CAT+p_photon_step] = stack_size-1;
#ifdef PRINT_PART_NUM
				logstr("ID%d:Bat%3d:%8d\n",GPUId,iBatch,h_np[GPUId*NUM_SIMU_CAT+p_photon_step]);
#endif
				iBatch++;
			}
#ifndef PRINT_PART_NUM
			complete = (float)istart*(nRepeat+1)/(float)num_histories;
			printf("\rField #%d, %.0f %% completed",iField,complete*100.0F);
#endif
			
			// clear momory to store engery deposition
			for(int GPUId=0; GPUId<GPUNo; GPUId++){
#ifdef USE_MULTIPLE_GPU
				cudaSetDevice(GPUId); ce(58056);
#endif
				cudaMemset(d_eng_score[GPUId], 0, h_size_phantom*sizeof(score_t)); ce(7005);
			}

			bool dosimu=false;
			bool done = false;
			for(int j=0;true;j++){
			
				done = false;
				for(char k=0;k<9;k++){

					//photon_step
					dosimu=false;
					for(int GPUId=0; GPUId<GPUNo; GPUId++){
						dosimu|=(h_np[GPUId*NUM_SIMU_CAT+p_photon_step]>=0);
					}
					//if(h_np[p_photon_step]>=0){
					if(dosimu){
						all_stack_empty = launch_kernels(p_photon_step, h_np, time_sim);
#ifndef GET_UNCERTAINTY
						if(all_stack_empty&&!load_stack){
							done = true;
							break;
						}
#else
						if(all_stack_empty){
							done = true;
							break;
						}
#endif
					}
					
					//compton
					dosimu=false;
					for(int GPUId=0; GPUId<GPUNo; GPUId++){
						dosimu|=(h_np[GPUId*NUM_SIMU_CAT+p_compton]>=0);
					}
					//if(h_np[p_compton]>=0){
					if(dosimu){
						all_stack_empty = launch_kernels(p_compton, h_np, time_sim);
#ifndef GET_UNCERTAINTY
						if(all_stack_empty&&!load_stack){
							done = true;
							break;
						}
#else
						if(all_stack_empty){
							done = true;
							break;
						}
#endif
					}

					//photo
					dosimu=false;
					for(int GPUId=0; GPUId<GPUNo; GPUId++){
						dosimu|=(h_np[GPUId*NUM_SIMU_CAT+p_photo]>=0);
					}
					//if(h_np[p_photo]>=0){
					if(dosimu){
						all_stack_empty = launch_kernels(p_photo, h_np, time_sim);
#ifndef GET_UNCERTAINTY
						if(all_stack_empty&&!load_stack){
							done = true;
							break;
						}
#else
						if(all_stack_empty){
							done = true;
							break;
						}
#endif
					}

					//electron_step
					dosimu=false;
					for(int GPUId=0; GPUId<GPUNo; GPUId++){
						dosimu|=(h_np[GPUId*NUM_SIMU_CAT+e_electron_step]>=0);
					}
					//if(h_np[e_electron_step]>=0){
					if(dosimu){
						all_stack_empty = launch_kernels(e_electron_step, h_np, time_sim);
#ifndef GET_UNCERTAINTY
						if(all_stack_empty&&!load_stack){
							done = true;
							break;
						}
#else
						if(all_stack_empty){
							done = true;
							break;
						}
#endif
					}


					//moller
					dosimu=false;
					for(int GPUId=0; GPUId<GPUNo; GPUId++){
						dosimu|=(h_np[GPUId*NUM_SIMU_CAT+e_moller]>=0);
					}
					//if(h_np[e_moller]>=0){
					if(dosimu){
						all_stack_empty = launch_kernels(e_moller, h_np, time_sim);
#ifndef GET_UNCERTAINTY
						if(all_stack_empty&&!load_stack){
							done = true;
							break;
						}
#else
						if(all_stack_empty){
							done = true;
							break;
						}
#endif
					}

				}
				if(done) break;

				//photon_step
				dosimu=false;
				for(int GPUId=0; GPUId<GPUNo; GPUId++){
					dosimu|=(h_np[GPUId*NUM_SIMU_CAT+p_photon_step]>=0);
				}
				//if(h_np[p_photon_step]>=0){
				if(dosimu){
					all_stack_empty = launch_kernels(p_photon_step, h_np, time_sim);
#ifndef GET_UNCERTAINTY
					if(all_stack_empty&&!load_stack) break;
#else
					if(all_stack_empty) break;
#endif
				}

#ifndef GET_UNCERTAINTY
				//load_from_phsp
				dosimu=true;
				for(int GPUId=0; GPUId<GPUNo; GPUId++){
					dosimu&=(h_np[GPUId*NUM_SIMU_CAT+p_photon_step]<0);
				}
				//if(h_np[p_photon_step]<0){
				if(dosimu){
					for(int GPUId=0; GPUId<GPUNo; GPUId++){
#ifdef USE_MULTIPLE_GPU
						cudaSetDevice(GPUId); ce(58049);
#endif
						if(load_stack){
							istart=iBatch*batch_size;
							if(iBatch<nBatch){
								stack_size = load_stack_from_phspfile(istart,batch_size,nRepeat);
								//logstr("Batch%2d: np=%7d\n",iBatch,h_np[p_photon_step]);
							}
							else if(iBatch==nBatch){
								stack_size = load_stack_from_phspfile(istart,lastBatch,nRepeat);
								load_stack=false;
								//logstr("Batch%2d: np=%7d\n",iBatch,h_np[p_photon_step]);
							}
							h_np[GPUId*NUM_SIMU_CAT+p_photon_step] = stack_size-1;
#ifdef PRINT_PART_NUM
							logstr("ID%d:Bat%3d:%8d\n",GPUId,iBatch,h_np[GPUId*NUM_SIMU_CAT+p_photon_step]);
#endif
							iBatch++;
						}
					}
#ifndef PRINT_PART_NUM
					complete = (float)istart*(nRepeat+1)/(float)num_histories;
					printf("\rField #%d, %.0f %% completed",iField,complete*100.0F);
#endif
				}
#endif

				//compton
				dosimu=false;
				for(int GPUId=0; GPUId<GPUNo; GPUId++){
					dosimu|=(h_np[GPUId*NUM_SIMU_CAT+p_compton]>=0);
				}
				//if(h_np[p_compton]>=0){
				if(dosimu){
					all_stack_empty = launch_kernels(p_compton, h_np, time_sim);
#ifndef GET_UNCERTAINTY
					if(all_stack_empty&&!load_stack) break;
#else
					if(all_stack_empty) break;
#endif
				}

				//photo
				dosimu=false;
				for(int GPUId=0; GPUId<GPUNo; GPUId++){
					dosimu|=(h_np[GPUId*NUM_SIMU_CAT+p_photo]>=0);
				}
				//if(h_np[p_photo]>=0){
				if(dosimu){
					all_stack_empty = launch_kernels(p_photo, h_np, time_sim);
#ifndef GET_UNCERTAINTY
					if(all_stack_empty&&!load_stack) break;
#else
					if(all_stack_empty) break;
#endif
				}

				//rayleigh
				dosimu=false;
				for(int GPUId=0; GPUId<GPUNo; GPUId++){
					dosimu|=(h_np[GPUId*NUM_SIMU_CAT+p_rayleigh]>=0);
				}
				//if(h_np[p_rayleigh]>=0){
				if(dosimu){
					all_stack_empty = launch_kernels(p_rayleigh, h_np, time_sim);
#ifndef GET_UNCERTAINTY
					if(all_stack_empty&&!load_stack) break;
#else
					if(all_stack_empty) break;
#endif
				}

				//pair
				dosimu=false;
				for(int GPUId=0; GPUId<GPUNo; GPUId++){
					dosimu|=(h_np[GPUId*NUM_SIMU_CAT+p_pair]>=0);
				}
				//if(h_np[p_pair]>=0){
				if(dosimu){
					all_stack_empty = launch_kernels(p_pair, h_np, time_sim);
#ifndef GET_UNCERTAINTY
					if(all_stack_empty&&!load_stack) break;
#else
					if(all_stack_empty) break;
#endif
				}

				//electron_step
				dosimu=false;
				for(int GPUId=0; GPUId<GPUNo; GPUId++){
					dosimu|=(h_np[GPUId*NUM_SIMU_CAT+e_electron_step]>=0);
				}
				//if(h_np[e_electron_step]>=0){
				if(dosimu){
					all_stack_empty = launch_kernels(e_electron_step, h_np, time_sim);
#ifndef GET_UNCERTAINTY
					if(all_stack_empty&&!load_stack) break;
#else
					if(all_stack_empty) break;
#endif
				}

				//moller
				dosimu=false;
				for(int GPUId=0; GPUId<GPUNo; GPUId++){
					dosimu|=(h_np[GPUId*NUM_SIMU_CAT+e_moller]>=0);
				}
				//if(h_np[e_moller]>=0){
				if(dosimu){
					all_stack_empty = launch_kernels(e_moller, h_np, time_sim);
#ifndef GET_UNCERTAINTY
					if(all_stack_empty&&!load_stack) break;
#else
					if(all_stack_empty) break;
#endif
				}

				//brems
				dosimu=false;
				for(int GPUId=0; GPUId<GPUNo; GPUId++){
					dosimu|=(h_np[GPUId*NUM_SIMU_CAT+e_brems]>=0);
				}
				//if(h_np[e_brems]>=0){
				if(dosimu){
					all_stack_empty = launch_kernels(e_brems, h_np, time_sim);
#ifndef GET_UNCERTAINTY
					if(all_stack_empty&&!load_stack) break;
#else
					if(all_stack_empty) break;
#endif
				}

				//bhabha
				dosimu=false;
				for(int GPUId=0; GPUId<GPUNo; GPUId++){
					dosimu|=(h_np[GPUId*NUM_SIMU_CAT+e_bhabha]>=0);
				}
				//if(h_np[e_bhabha]>=0){
				if(dosimu){
					all_stack_empty = launch_kernels(e_bhabha, h_np, time_sim);
#ifndef GET_UNCERTAINTY
					if(all_stack_empty&&!load_stack) break;
#else
					if(all_stack_empty) break;
#endif
				}

				//annih
				dosimu=false;
				for(int GPUId=0; GPUId<GPUNo; GPUId++){
					dosimu|=(h_np[GPUId*NUM_SIMU_CAT+e_annih]>=0);
				}
				//if(h_np[e_annih]>=0){
				if(dosimu){
					all_stack_empty = launch_kernels(e_annih, h_np, time_sim);
#ifndef GET_UNCERTAINTY
					if(all_stack_empty&&!load_stack) break;
#else
					if(all_stack_empty) break;
#endif
				}

			}

			/*
#ifndef GET_UNCERTAINTY
			for(int GPUId=0; GPUId<GPUNo; GPUId++) {
#ifdef USE_MULTIPLE_GPU
				cudaSetDevice(GPUId); ce(58015);
#endif
				//cudaMemcpyFromSymbol(&d_eng_score,eng_score,sizeof(score_t*)); ce(110181);
				cudaMemcpy(eng_dep2, d_eng_score[GPUId], h_size_phantom*sizeof(score_t), cudaMemcpyDeviceToHost); ce(11018);
				for(int i=0; i<h_size_phantom; i++){
					eng_dep[i] += eng_dep2[i]*output_factor;
					//eng_dep[i] += eng_dep2[i]/voxel_mass[i]*dose_unit*output_factor;
				}
			}
#else
			*/

#ifdef GET_UNCERTAINTY
			score_t *h_eng_score = (score_t *)malloc(h_size_phantom*sizeof(score_t));

			if(load_stack){  //here we don't take the last launch into account
				for(int GPUId=0; GPUId<GPUNo; GPUId++){
#ifdef USE_MULTIPLE_GPU
					cudaSetDevice(GPUId); ce(58061);
#endif
					//copy the total dose from GPU to Host for analysis and sum 
					cudaMemcpy(h_eng_score, d_eng_score[GPUId], h_size_phantom*sizeof(score_t), cudaMemcpyDeviceToHost); ce(11018); 

					// calculate energy deposition uncertianty
					for(int iReg=0;iReg<h_size_phantom;iReg++){
						//devided by the mass to get dose, then change dose unit from MeV/g to cGy
						//h_eng_score[iReg] = h_eng_score[iReg]*output_factor;
						h_eng_score[iReg] = h_eng_score[iReg];
						//h_eng_score[iReg] = h_eng_score[iReg]/voxel_mass[iReg]*dose_unit*output_factor;
						eng_dep[iReg] += h_eng_score[iReg];
						eng_dep2[iReg] += h_eng_score[iReg]*h_eng_score[iReg];
					}
				}
			}
			else{  //take care of last launch
				int BatNo = (iBatch-1)%GPUNo;
				int GPUId = 0;
				for(; GPUId<BatNo; GPUId++){
#ifdef USE_MULTIPLE_GPU
					cudaSetDevice(GPUId); ce(58062);
#endif
					//copy the total dose from GPU to Host for analysis and sum 
					cudaMemcpy(h_eng_score, d_eng_score[GPUId], h_size_phantom*sizeof(score_t), cudaMemcpyDeviceToHost); ce(11018); 

					// calculate energy deposition uncertianty
					for(int iReg=0;iReg<h_size_phantom;iReg++){
						//devided by the mass to get dose, then change dose unit from MeV/g to cGy
						//h_eng_score[iReg] = h_eng_score[iReg]*output_factor;
						h_eng_score[iReg] = h_eng_score[iReg];
						//h_eng_score[iReg] = h_eng_score[iReg]/voxel_mass[iReg]*dose_unit*output_factor;
						eng_dep[iReg] += h_eng_score[iReg];
						eng_dep2[iReg] += h_eng_score[iReg]*h_eng_score[iReg];
					}
				}

				//now this is the last batch
#ifdef USE_MULTIPLE_GPU
				cudaSetDevice(GPUId); ce(58063);
#endif
				//copy the total dose from GPU to Host for analysis and sum 
				cudaMemcpy(h_eng_score, d_eng_score[GPUId], h_size_phantom*sizeof(score_t), cudaMemcpyDeviceToHost); ce(11018); 

				for(int iReg=0;iReg<h_size_phantom;iReg++){  //nBatch is the batch number without last batch
					eng_dep2[iReg] = sqrt((eng_dep2[iReg]/nBatch - (eng_dep[iReg]/nBatch)*(eng_dep[iReg]/nBatch))/(nBatch-1));
					eng_dep2[iReg] = eng_dep2[iReg]/eng_dep[iReg]*(float)nBatch;  //change the uncertainty to percentage
					//devided by the mass to get dose, then change dose unit from MeV/g to cGy
					//h_eng_score[iReg] = h_eng_score[iReg]*output_factor;
					h_eng_score[iReg] = h_eng_score[iReg];
					//h_eng_score[iReg] = h_eng_score[iReg]/voxel_mass[iReg]*dose_unit*output_factor;
					eng_dep[iReg] += h_eng_score[iReg];  //then add the last batch dose (batch # = nbatch+1)
				}
			}

			free(h_eng_score);

		}
#endif

    printf("\r");

    // log time
	getTime();
	logstr("Field #%d Finished simulation at %s\n", iField, charBuffer);

	clock_t run_tac = clock();

    float total_time = (float)(run_tac - run_tic) / (float)CLOCKS_PER_SEC * 1000.0F;
	logstr("  Simulation time . . . . . . %.2f min\n\n", total_time/6.0E+4F);

	free(h_np);
	free(time_sim);

	free_stack();
	free_leaf_pos(iField);
	//end of simulation of ONE MLC field
}

int main(int argc, char **argv) {
    
	clock_t tic = clock();
	
	// get input file
	int idxI = -1;
	for (int i = 0; i < argc; i++) {
		if (string(argv[i]).compare("-i") == 0) {
			idxI = i;
			break;
		}
	}

	if ((idxI == -1) || (idxI == argc - 1)) {
		printf("ERROR (1001): No input file was specified. Use the default input file, CUDA_EGS_example.ini\n");
		input_file = "CUDA_EGS_example.ini";
	}
	else {
		input_file = argv[idxI + 1];
	}

	FILE *input_f;
	if (fopen_s(&input_f, input_file, "r")) {
		printf("ERROR (1002): The input file \"%s\" could not be opened.\n", input_file);
		exit(1002);
	}
	fclose(input_f);

    // get the full path of the input file, otherwise GetPrivateProfileString etc. might not 
    // be able to read it
    GetFullPathName(input_file, CHARLEN, charBuffer, NULL);
    string input_file_input = string(charBuffer);
    input_file = input_file_input.c_str();

	bool defaultOutputPrefixUsed = false;
	string output_prefix_input = "";
	GetPrivateProfileString("General", "output prefix", "<<missing>>", charBuffer, CHARLEN, input_file);
	if (string(charBuffer) == "<<missing>>") {
		defaultOutputPrefixUsed = true;
		output_prefix_input = "GEGS";
	}
	else {
		output_prefix_input = string(charBuffer);
		//checkOutputPrefix(2001, output_prefix_input);
	}

	bool defaultOutputDirectoryUsed = false;
	string output_directory_input = "";
	GetPrivateProfileString("General", "output directory", "<<missing>>", charBuffer, CHARLEN, input_file);
	if (string(charBuffer) == "<<missing>>" || string(charBuffer) == "") {
		defaultOutputDirectoryUsed = true;
		// get current working directory
		GetCurrentDirectory(CHARLEN, charBuffer);
		output_directory_input = string(charBuffer);
		output_directory_input += "\\";
	}
	else {
		output_directory_input = string(charBuffer);
		output_directory_input += "\\";
		checkOutputPrefix(2001, output_directory_input);
		// get the absolute path of the output directory
		GetFullPathName(output_directory_input.c_str(), CHARLEN, charBuffer, NULL);
		output_directory_input = string(charBuffer);
	}
	output_dir = output_directory_input.c_str();
	string rootName = output_directory_input + output_prefix_input;

	// create log file
	fopen_s(&logfile, (rootName + "_log.txt").c_str(), "w");
	
	// log time
	getTime();
	logstr("\nGEGS with MLC version 1.0.0 started at %s\n", charBuffer);

	// log info
    logstr("\nSimulation\n");
	logstr("  Input file  . . . . . . . . %s\n", input_file);
	if (defaultOutputPrefixUsed)
		logstr("WARNING: No output prefix was specified. Default value is used.\n");
	logstr("  Output prefix . . . . . . . %s\n", output_prefix_input.c_str());
	if (defaultOutputDirectoryUsed)
		logstr("WARNING: No output directory was specified. Current Working Directory is used.\n");
	logstr("  Output directory  . . . . . %s\n", output_dir);

#ifdef USE_MULTIPLE_GPU
    // get GPU number
    bool defaultGPUNo = false;
    // first check whether GPU number is given
    GetPrivateProfileString("General", "GPU number", "<<missing>>", charBuffer, CHARLEN, input_file);
	if (string(charBuffer) == "<<missing>>") {
            defaultGPUNo = true;
    }
	else {
		int givenNo;
		if (sscanf(charBuffer, "%d", &givenNo) != 1)
			error(2004, "The given GPU number \"%s\" could not be parsed as a 4-byte integer.\n", charBuffer);
        GPUNo = givenNo;
	}

    if (defaultGPUNo)
		logstr("WARNING: No GPU number was given, default to be 1.\n");

#else

	// get GPU id
	bool defaultGPUUsed = false;
	int GPU_Id = 0;
	// first check whether a GPU id is given
	GetPrivateProfileString("General", "GPU id", "<<missing>>", charBuffer, CHARLEN, input_file);
	if (string(charBuffer) == "<<missing>>") {
		// no GPU id is given, check whether a GPU name is given
			defaultGPUUsed = true;
			GPU_Id = 0;
	}
	else {
		int givenId;
		if (sscanf(charBuffer, "%d", &givenId) != 1)
			error(2004, "The given GPU id \"%s\" could not be parsed as a 4-byte integer.\n", charBuffer);
		GPU_Id = givenId;
	}

	if (defaultGPUUsed)
		logstr("WARNING: No GPU was specified. Default GPU is used.\n");
#endif

    // check whether GPU has CUDA capability 2.0
    int deviceCount = 0;
    cudaGetDeviceCount(&deviceCount); ce(2005);
    if (deviceCount <= 0)
        error(2006, "No CUDA capabale GPUs found.\n");
    if (GPUNo > deviceCount)
        error(2007, "The number of CUDA capable GPU is smaller than GPU number %d.\n", GPUNo);
    
    cudaDeviceProp deviceProp;
#ifdef USE_MULTIPLE_GPU
	for(int GPUId=0; GPUId<GPUNo; GPUId++) {
		cudaGetDeviceProperties(&deviceProp, GPUId); ce(2021);
		logstr("  GPU . . . . . . . . . . . . %s (id %d)\n", deviceProp.name, GPUId);

		if (deviceProp.major < 2)
			error(2008, "CUDA capability 2.0 or higher is required, but the specified GPU only supports CUDA capability %d.%d.\n", deviceProp.major, deviceProp.minor);
		if (deviceProp.warpSize != WARP_SIZE)
			error(2009, "The code was compiled for a warp size of %d, but the GPU uses a warp size of %d. Please change the value of the constant WARP_SIZE in CUDA_EGS.h and recompile.\n", WARP_SIZE, deviceProp.warpSize);

		if (deviceProp.multiProcessorCount != NUM_MULTIPROC)
			error(2010, "The code was compiled for %d multiprocessors, but the GPU has %d multiprocessors. Please change the value of the constant NUM_MULTIPROC in CUDA_EGS.h and recompile.\n", NUM_MULTIPROC, deviceProp.multiProcessorCount);

		logstr("    Warp size . . . . . . . . %d\n", deviceProp.warpSize);
		logstr("    Multiprocessors . . . . . %d\n", deviceProp.multiProcessorCount);

	}
#else
	cudaGetDeviceProperties(&deviceProp, GPU_Id); ce(2021);
	logstr("  GPU . . . . . . . . . . . . %s (id %d)\n", deviceProp.name, GPU_Id);

	if (deviceProp.major < 2)
		error(2008, "CUDA capability 2.0 or higher is required, but the specified GPU only supports CUDA capability %d.%d.\n", deviceProp.major, deviceProp.minor);
	if (deviceProp.warpSize != WARP_SIZE)
		error(2009, "The code was compiled for a warp size of %d, but the GPU uses a warp size of %d. Please change the value of the constant WARP_SIZE in CUDA_EGS.h and recompile.\n", WARP_SIZE, deviceProp.warpSize);

	if (deviceProp.multiProcessorCount != NUM_MULTIPROC)
		error(2010, "The code was compiled for %d multiprocessors, but the GPU has %d multiprocessors. Please change the value of the constant NUM_MULTIPROC in CUDA_EGS.h and recompile.\n", NUM_MULTIPROC, deviceProp.multiProcessorCount);

	logstr("    Warp size . . . . . . . . %d\n", deviceProp.warpSize);
	logstr("    Multiprocessors . . . . . %d\n", deviceProp.multiProcessorCount);

    cudaSetDevice(GPU_Id); ce(58001);
#endif

    // get histories
    GetPrivateProfileString("General", "histories", "<<missing>>", charBuffer, CHARLEN, input_file);
	if (string(charBuffer) == "<<missing>>")
		error(2013, "Property \"histories\" is not speciefied in the section \"General\" in the input file \"%s\".\n", input_file);
	float histories_per_MU_input;
	if (sscanf(charBuffer, "%f", &histories_per_MU_input) != 1)
		error(2014, "Could not parse the number of histories \"%s\" in the input file \"%s\".\n", charBuffer, input_file);
    if (histories_per_MU_input <= 0)
		error(2015, "The number of histories must be positive, but %f was specified in the input file \"%s\".\n", histories_per_MU_input, input_file);
	//logstr("  Total Histories . . . . . . %.3E\n", histories_per_MU_input);
	logstr("  Histories per MU  . . . . . %.3E\n", histories_per_MU_input);
	ulong histories_per_MU = (ulong)histories_per_MU_input;

    // get data specifications
    GetPrivateProfileString("Data", "data directory", "<<missing>>", charBuffer, CHARLEN, input_file);
	string data_dir_input = string(charBuffer);
    if (data_dir_input == "<<missing>>")
		error(2017, "Property \"data directory\" is not speciefied in the section \"Data\" in the input file \"%s\".\n", input_file);
    data_dir_input += "\\";

    // get the absolute path of the data directory
    GetFullPathName(data_dir_input.c_str(), CHARLEN, charBuffer, NULL);
    data_dir_input = string(charBuffer);
    data_dir = data_dir_input.c_str();
    logstr("  Data directory  . . . . . . %s\n", data_dir);

    // get egsinp file
    GetPrivateProfileString("Data", "egsinp file", "<<missing>>", charBuffer, CHARLEN, input_file);
	string egsinp_file_input = string(charBuffer);
	if (egsinp_file_input == "<<missing>>")
		error(2011, "Property \"egsinp file\" is not speciefied in the section \"Data\" in the input file \"%s\".\n", input_file);
	/*
    // get the absolute path of the egsinp file
    GetFullPathName(egsinp_file_input.c_str(), CHARLEN, charBuffer, NULL);
    egsinp_file_input = string(charBuffer);
    egsinp_file = egsinp_file_input.c_str();
    logstr("  egsinp file . . . . . . . . %s\n", egsinp_file);
	*/
    egsinp_file_input = data_dir_input + egsinp_file_input;
    egsinp_file = egsinp_file_input.c_str();
    logstr("  egsinp file . . . . . . . . %s\n", egsinp_file);

    // get pegs file
    GetPrivateProfileString("Data", "pegs file", "<<missing>>", charBuffer, CHARLEN, input_file);
	string pegs_file_input = string(charBuffer);
	if (pegs_file_input == "<<missing>>")
		error(2012, "Property \"pegs file\" is not speciefied in the section \"Data\" in the input file \"%s\".\n", input_file);
    pegs_file_input = data_dir_input + pegs_file_input;
    pegs_file = pegs_file_input.c_str();
    logstr("  PEGS4 file  . . . . . . . . %s\n", pegs_file);
	
    GetPrivateProfileString("Data", "photon xsections", "<<missing>>", charBuffer, CHARLEN, input_file);
	string photon_xsections_input = string(charBuffer);
    if (photon_xsections_input == "<<missing>>")
		error(2019, "Property \"photon xsections\" is not speciefied in the section \"Data\" in the input file \"%s\".\n", input_file);
    photon_xsections = photon_xsections_input.c_str();
    logstr("  Photon xsections  . . . . . %s\n", photon_xsections);

    GetPrivateProfileString("Data", "atomic ff file", "<<missing>>", charBuffer, CHARLEN, input_file);
	string atomic_ff_file_input = string(charBuffer);
    if (atomic_ff_file_input == "<<missing>>")
		error(2020, "Property \"atomic ff file\" is not speciefied in the section \"Data\" in the input file \"%s\".\n", input_file);
    atomic_ff_file = atomic_ff_file_input.c_str();
    logstr("  Atomic ff file  . . . . . . %s\n", atomic_ff_file);

	/*
    // get output factor file
    GetPrivateProfileString("Data", "output factor", "<<missing>>", charBuffer, CHARLEN, input_file);
	string output_factor_file_input = string(charBuffer);
	if (output_factor_file_input == "<<missing>>")
		error(2012, "Property \"output factor\" is not speciefied in the section \"Data\" in the input file \"%s\".\n", input_file);
    output_factor_file_input = data_dir_input + output_factor_file_input;
    output_factor_file = output_factor_file_input.c_str();
    logstr("  Output factor file  . . . . %s\n", output_factor_file);
	*/

    // read default medium
    GetPrivateProfileString("Data", "default medium", "<<missing>>", charBuffer, CHARLEN, input_file);
	string default_medium = string(charBuffer);
	if (default_medium == "<<missing>>")
		error(2013, "Property \"default medium\" is not speciefied in the section \"Data\" in the input file \"%s\".\n", input_file);
	logstr("  default medium  . . . . . . %s\n", default_medium.c_str());

	// read phase space file number from the input file
	GetPrivateProfileString("Phase space", "phsp file number", "<<missing>>", charBuffer, CHARLEN, input_file);
	if (string(charBuffer) == "<<missing>>")
		error(5013, "Property \"phsp file number\" is not speciefied in the section \"Phase space\" in the input file \"%s\".\n", input_file);
	if (sscanf(charBuffer, "%d", &nPhsp) != 1)
		error(5014, "Could not parse phsp file number \"%s\" in the input file \"%s\".\n", charBuffer, input_file);
	logstr("  number of phsp file . . . . %d\n", nPhsp);

    // get Phase space file
    GetPrivateProfileString("Phase space", "phase space file", "<<missing>>", charBuffer, CHARLEN, input_file);
	string phsp_file_input = string(charBuffer);
	if (phsp_file_input == "<<missing>>")
		error(2012, "Property \"phase space file\" is not speciefied in the section \"Phase space\" in the input file \"%s\".\n", input_file);
    phsp_file_input = data_dir_input + phsp_file_input;
    phsp_file = phsp_file_input.c_str();

	read_egsphsp_header(phsp_file);

	// read z position of phase space from the input file
	GetPrivateProfileString("Phase space", "z position", "<<missing>>", charBuffer, CHARLEN, input_file);
	if (string(charBuffer) == "<<missing>>")
		error(5013, "Property \"position\" is not speciefied in the section \"Phase space\" in the input file \"%s\".\n", input_file);
	float z;
	if (sscanf(charBuffer, "%f", &z) != 1)
		error(5014, "Could not parse the z position \"%s\" in the input file \"%s\".\n", charBuffer, input_file);
    z_phsp = z;
	logstr("  z position of phsp  . . . . %.2f\n", z_phsp);
	
	// read flag of module SecJaws from the input file
	SecJaws_on = false;
	GetPrivateProfileString("SecJaws", "module SecJaws", "<<missing>>", charBuffer, CHARLEN, input_file);
	if (string(charBuffer) == "<<missing>>")
		logstr("Property \"module SecJaws\" is not speciefied in the section \"SecJaws\", default to be Off.\n");
	if (string(charBuffer) == "On")
		SecJaws_on = true;
	else if(string(charBuffer) != "Off")
		logstr("WARNING: Default is not to use ctcreate.\n");
	if (SecJaws_on)
		logstr("  module SecJaws  . . . . . . On\n");
	else
		logstr("  module SecJaws  . . . . . . Off\n");

	// get SecJaws file
	string secjaws_file_input;
	//if(SecJaws_on){
		GetPrivateProfileString("SecJaws", "secjaws file", "<<missing>>", charBuffer, CHARLEN, input_file);
		secjaws_file_input = string(charBuffer);
		if (secjaws_file_input == "<<missing>>")
			error(2011, "Property \"secjaws file\" is not speciefied in the section \"SecJaws\" in the input file \"%s\".\n", input_file);
		secjaws_file_input = data_dir_input + secjaws_file_input;
		logstr("  secjaws file  . . . . . . . %s\n", secjaws_file_input.c_str());
	//}
	secjaws_file = secjaws_file_input.c_str();

	// read flag of module MLC from the input file
	VarMLC_on = false;
	GetPrivateProfileString("MLC", "module MLC", "<<missing>>", charBuffer, CHARLEN, input_file);
	if (string(charBuffer) == "<<missing>>")
		logstr("Property \"module MLC\" is not speciefied in the section \"MLC\", default to be Off.\n");
	if (string(charBuffer) == "On")
		VarMLC_on = true;
	else if(string(charBuffer) != "Off")
		logstr("WARNING: Default is not to use ctcreate.\n");
	if (VarMLC_on)
		logstr("  module MLC  . . . . . . . . On\n");
	else
		logstr("  module MLC  . . . . . . . . Off\n");

	// get mlcinfo file
	GetPrivateProfileString("MLC", "mlcinfo file", "<<missing>>", charBuffer, CHARLEN, input_file);
	string mlcinfo_file_input = string(charBuffer);
	if (mlcinfo_file_input == "<<missing>>")
		error(2011, "Property \"mlcinfo file\" is not speciefied in the section \"MLC\" in the input file \"%s\".\n", mlcinfo_file);
	mlcinfo_file_input = data_dir_input + mlcinfo_file_input;
	mlcinfo_file = mlcinfo_file_input.c_str();
	logstr("  mlcinfo file  . . . . . . . %s\n", mlcinfo_file);

	// read flag of module Block from the input file
	Block_on = false;
	GetPrivateProfileString("Block", "module Block", "<<missing>>", charBuffer, CHARLEN, input_file);
	if (string(charBuffer) == "<<missing>>")
		logstr("Property \"module Block\" is not speciefied in the section \"Block\", default to be Off.\n");
	if (string(charBuffer) == "On")
		Block_on = true;
	if (Wedge_on)
		logstr("  module Block  . . . . . . . On\n");
	else
		logstr("  module Block  . . . . . . . Off\n");
	
	// get block file
	string block_file_input;
	if(Block_on){
		GetPrivateProfileString("Block", "block file", "<<missing>>", charBuffer, CHARLEN, input_file);
		block_file_input = string(charBuffer);
		if (block_file_input == "<<missing>>")
			error(2011, "Property \"block file\" is not speciefied in the section \"Block\" in the input file \"%s\".\n", input_file);
		block_file_input = data_dir_input + block_file_input;
		logstr("  block file  . . . . . . . . %s\n", block_file_input.c_str());
	}
	block_file = block_file_input.c_str();

	// read flag of module Wedge from the input file
	Wedge_on = false;
	GetPrivateProfileString("Wedge", "module Wedge", "<<missing>>", charBuffer, CHARLEN, input_file);
	if (string(charBuffer) == "<<missing>>")
		logstr("Property \"module Wedge\" is not speciefied in the section \"Wedge\", default to be Off.\n");
	if (string(charBuffer) == "On")
		Wedge_on = true;
	if (Wedge_on)
		logstr("  module Wedge  . . . . . . . On\n");
	else
		logstr("  module Wedge  . . . . . . . Off\n");
	
	// get wedge file
	string wedge_file_input;
	if(Wedge_on){
		GetPrivateProfileString("Wedge", "wedge file", "<<missing>>", charBuffer, CHARLEN, input_file);
		wedge_file_input = string(charBuffer);
		if (wedge_file_input == "<<missing>>")
			error(2011, "Property \"wedge file\" is not speciefied in the section \"Wedge\" in the input file \"%s\".\n", input_file);
		wedge_file_input = data_dir_input + wedge_file_input;
		logstr("  wedge file  . . . . . . . . %s\n", wedge_file_input.c_str());
	}
	wedge_file = wedge_file_input.c_str();

	//set the first module
	char next_module = m_Phantom;
	if(Wedge_on){
		h_module_order[m_WedgeMd] = next_module;
		next_module = m_WedgeMd;
	}
	if(Block_on){
		h_module_order[m_BlockMd] = next_module;
		next_module = m_BlockMd;
	}
	if(VarMLC_on){
		h_module_order[m_VarMLCs] = next_module;
		next_module = m_VarMLCs;
	}
	if(SecJaws_on){
		h_module_order[m_SecJawX] = next_module;
		h_module_order[m_SecJawY] = m_SecJawX;
		next_module = m_SecJawY;
	}
	h_module_order[m_FirstMd] = next_module;

	for(int GPUId=0; GPUId<GPUNo; GPUId++) {
#ifdef USE_MULTIPLE_GPU
		cudaSetDevice(GPUId); ce(58060);
#endif
		cudaMemcpyToSymbol(module_order,h_module_order,MAX_MODULE_PLUS_1*sizeof(char)); ce(1919);
	}

	// read flag of use ctcreate from the input file
	create_on = false;
	GetPrivateProfileString("CTcreate", "use ctcreate", "<<missing>>", charBuffer, CHARLEN, input_file);
	if (string(charBuffer) == "<<missing>>")
		logstr("Property \"use ctcreate\" is not speciefied in the section \"CTcreate\", default to be No.\n");
	if (string(charBuffer) == "Yes")
		create_on = true;
	else if(string(charBuffer) != "No")
		logstr("WARNING: Default is not to use ctcreate.\n");
	if (create_on)
		logstr("  use ctcreate  . . . . . . . Yes\n");
	else
		logstr("  use ctcreate  . . . . . . . No\n");

	string ct_ramps_file_input = "";
	//string dicom_directory_input = "";
	string CTFileName_input = "";
	if(create_on){
		// get CT ramps file
		GetPrivateProfileString("CTcreate", "CT ramps file", "<<missing>>", charBuffer, CHARLEN, input_file);
		ct_ramps_file_input = string(charBuffer);
		if (ct_ramps_file_input == "<<missing>>")
			error(2012, "Property \"CT ramps file\" is not speciefied in the section \"CTcreate\" in the input file \"%s\".\n", input_file);
		ct_ramps_file_input = data_dir_input + ct_ramps_file_input;
		logstr("  CT ramps file . . . . . . . %s\n", ct_ramps_file_input.c_str());

		/*
		//for relative path in dicom file list
		// get dicom directory
		GetPrivateProfileString("CTcreate", "dicom directory", "<<missing>>", charBuffer, CHARLEN, input_file);
		dicom_directory_input = string(charBuffer);
		if (dicom_directory_input == "<<missing>>")
			error(2012, "Property \"dicom directory\" is not speciefied in the section \"CTcreate\" in the input file \"%s\".\n", input_file);
		dicom_directory_input += "\\";
		// get the absolute path of dicom file list
		GetFullPathName(dicom_directory_input.c_str(), CHARLEN, charBuffer, NULL);
		dicom_directory_input = string(charBuffer);
		logstr("  dicom directory . . . . . . %s\n", dicom_directory_input.c_str());
		*/

		// get dicom file list
		GetPrivateProfileString("CTcreate", "dicom file list", "<<missing>>", charBuffer, CHARLEN, input_file);
		CTFileName_input = string(charBuffer);
		if (CTFileName_input == "<<missing>>")
			error(2012, "Property \"dicom file list\" is not speciefied in the section \"CTcreate\" in the input file \"%s\".\n", input_file);
		// get the absolute path of dicom file list
		GetFullPathName(CTFileName_input.c_str(), CHARLEN, charBuffer, NULL);
		CTFileName_input = string(charBuffer);
		//CTFileName_input = dicom_directory_input + CTFileName_input;  //for relative path in dicom file list
		logstr("  dicom file list . . . . . . %s\n", CTFileName_input.c_str());
	}
	ct_ramps_file = ct_ramps_file_input.c_str();
	//dicom_dir = dicom_directory_input.c_str();  //for relative path in dicom file list
	CTFileName = CTFileName_input.c_str();

	if(create_on){
		// read boundaries of the subset of CT from the input file
		GetPrivateProfileString("CTcreate", "xyz boundaries", "<<missing>>", charBuffer, CHARLEN, input_file);
		if (string(charBuffer) == "<<missing>>")
			error(2011, "Property \"xyz boundaries\" is not speciefied in the section \"CTcreate\" in the input file \"%s\".\n", input_file);
		float x1,x2,y1,y2,z1,z2;
		if (sscanf(charBuffer, "%f %f %f %f %f %f", &x1, &x2, &y1, &y2, &z1, &z2) != 6)
			error(2011, "Could not parse boundaries of the subset of CT \"%s\" in the input file \"%s\".\n", charBuffer, input_file);
		ctsubmin = make_float3(x1,y1,z1);
		ctsubmax = make_float3(x2,y2,z2);
		logstr("  phantom boundaries  . . . . %.2f, %.2f, %.2f, %.2f, %.2f, %.2f\n", ctsubmin.x, ctsubmin.y, ctsubmin.z, ctsubmax.x, ctsubmax.y, ctsubmax.z);

		// read voxel dimensions from the input file
		GetPrivateProfileString("CTcreate", "3d voxel sizes", "<<missing>>", charBuffer, CHARLEN, input_file);
		if (string(charBuffer) == "<<missing>>")
			error(2011, "Property \"3d voxel sizes\" is not speciefied in the section \"CTcreate\" in the input file \"%s\".\n", input_file);
		if (sscanf(charBuffer, "%f %f %f", &x1, &y1, &z1) != 3)
			error(2011, "Could not parse the voxel dimensions of CT \"%s\" in the input file \"%s\".\n", charBuffer, input_file);
		xyz_thickness = make_float3(x1,y1,z1);
		logstr("  phantom voxel sizes . . . . %.2f, %.2f, %.2f\n", xyz_thickness.x, xyz_thickness.y, xyz_thickness.z);

		ctcreate(rootName);
		
	}

	string egsphant_file_input;
	if(create_on){
		egsphant_file_input = rootName + ".egsphant";
		logstr("  Phantom from ctcreate . . . %s\n", egsphant_file_input.c_str());
	}
	else{
		// get egsphant file
		GetPrivateProfileString("Phantom", "egsphant file", "<<missing>>", charBuffer, CHARLEN, input_file);
		egsphant_file_input = string(charBuffer);
		if (egsphant_file_input == "<<missing>>")
			error(2011, "Property \"egsphant file\" is not speciefied in the section \"Phantom\" in the input file \"%s\".\n", input_file);
		egsphant_file_input = data_dir_input + egsphant_file_input;
		logstr("  Phantom from input  . . . . %s\n", egsphant_file_input.c_str());
	}
    egsphant_file = egsphant_file_input.c_str();

	// read flag of change coordinate from the input file
	h_Change_xyz = true;
	GetPrivateProfileString("Phantom", "change coordinate", "<<missing>>", charBuffer, CHARLEN, input_file);
	if (string(charBuffer) == "<<missing>>")
		logstr("Property \"change coordinate\" is not speciefied in the section \"Phantom\", default to be Yes.\n");
	if (string(charBuffer) == "No")
		h_Change_xyz = false;
	else if(string(charBuffer) != "Yes")
		logstr("WARNING: Default is to change coordinate system in phantom.\n");
	if (h_Change_xyz)
		logstr("  change xyz in phantom . . . Yes\n");
	else
		logstr("  change xyz in phantom . . . No\n");
	for(int GPUId=0; GPUId<GPUNo; GPUId++){
#ifdef USE_MULTIPLE_GPU
		cudaSetDevice(GPUId); ce(58057);
#endif
		cudaMemcpyToSymbol(Change_xyz,&h_Change_xyz,sizeof(bool)); ce(30003);
	}

	// read isocenter location from the input file
	GetPrivateProfileString("Phantom", "isocenter location", "<<missing>>", charBuffer, CHARLEN, input_file);
	if (string(charBuffer) == "<<missing>>")
		error(2011, "Property \"isocenter location\" is not speciefied in the section \"Phantom\" in the input file \"%s\".\n", input_file);
	float p1, p2, p3;
	if (sscanf(charBuffer, "%f %f %f", &p1, &p2, &p3) != 3)
		error(2011, "Could not parse the isocenter location \"%s\" in the input file \"%s\".\n", charBuffer, input_file);
    h_isocenter_location = make_float3(p1, p2, p3);
	logstr("  isocenter location  . . . . %.2f, %.2f, %.2f\n", h_isocenter_location.x, h_isocenter_location.y, h_isocenter_location.z);
	for(int GPUId=0; GPUId<GPUNo; GPUId++){
#ifdef USE_MULTIPLE_GPU
		cudaSetDevice(GPUId); ce(58058);
#endif
		cudaMemcpyToSymbol(isocenter_location, &h_isocenter_location,  sizeof(float3)); ce(30002);
	}

	// read distance from source to isocenter
	GetPrivateProfileString("Phantom", "source to isocenter", "<<missing>>", charBuffer, CHARLEN, input_file);
	if (string(charBuffer) == "<<missing>>")
		error(2011, "Property \"source to isocenter\" is not speciefied in the section \"Phantom\" in the input file \"%s\".\n", input_file);
	float SID;
	if (sscanf(charBuffer, "%f", &SID) != 1)
		error(2011, "Could not parse the isocenter location \"%s\" in the input file \"%s\".\n", charBuffer, input_file);
    h_SourceToIsocenter = SID;
	logstr("  source to isocenter . . . . %.2f\n", h_SourceToIsocenter);
	for(int GPUId=0; GPUId<GPUNo; GPUId++){
#ifdef USE_MULTIPLE_GPU
		cudaSetDevice(GPUId); ce(58059);
#endif
		cudaMemcpyToSymbol(SourceToIsocenter, &h_SourceToIsocenter,  sizeof(float)); ce(30002);
	}

	logstr("\nDose Calculation\n");
	// read flag of dose calibration from the input file
	calibr_on = false;
	GetPrivateProfileString("Calibration", "dose calibration", "<<missing>>", charBuffer, CHARLEN, input_file);
	if (string(charBuffer) == "<<missing>>")
		logstr("Property \"dose calibration\" is not speciefied in the section \"Calibration\", default to be No.\n");
	if (string(charBuffer) == "Yes")
		calibr_on = true;
	else if(string(charBuffer) != "No")
		logstr("WARNING: Default is not to do dose calibration.\n");
	if (calibr_on)
		logstr("  do dose calibration . . . . Yes\n");
	else
		logstr("  do dose calibration . . . . No\n");

	if(calibr_on){
		// read histories per MU from the input file
		GetPrivateProfileString("Calibration", "histories per cGy", "<<missing>>", charBuffer, CHARLEN, input_file);
		if (string(charBuffer) == "<<missing>>")
			error(5013, "Property \"histories per cGy\" is not speciefied in the section \"Calibration\" in the input file \"%s\".\n", input_file);
		float histories_per_cGy_input;
		if (sscanf(charBuffer, "%f", &histories_per_cGy_input) != 1)
			error(5014, "Could not parse histories per cGy \"%s\" in the input file \"%s\".\n", charBuffer, input_file);
		dose_unit = dose_unit*histories_per_cGy_input/histories_per_MU_input;
		logstr("  histories per cGy . . . . . %.3E\n", histories_per_cGy_input);
	}
	else
		logstr("  dose unit factor  . . . . . %.3E\n", dose_unit);

	// read flag of use output factor from the input file
	outfac_on = false;
	GetPrivateProfileString("Output factor", "use output factor", "<<missing>>", charBuffer, CHARLEN, input_file);
	if (string(charBuffer) == "<<missing>>")
		logstr("Property \"use output factor\" is not speciefied in the section \"Output factor\", default to be No.\n");
	if (string(charBuffer) == "Yes")
		outfac_on = true;
	else if(string(charBuffer) != "No")
		logstr("WARNING: Default is not to use output factor.\n");
	if (outfac_on)
		logstr("  use output factor . . . . . Yes\n");
	else
		logstr("  use output factor . . . . . No\n");

	// get output factor file
	string output_factor_file_input;
	if(outfac_on){
		GetPrivateProfileString("Output factor", "output factor file", "<<missing>>", charBuffer, CHARLEN, input_file);
		output_factor_file_input = string(charBuffer);
		if (output_factor_file_input == "<<missing>>")
			error(2012, "Property \"output factor file\" is not speciefied in the section \"Output factor\" in the input file \"%s\".\n", input_file);
		output_factor_file_input = data_dir_input + output_factor_file_input;
		logstr("  Output factor file  . . . . %s\n", output_factor_file_input.c_str());
	}
	output_factor_file = output_factor_file_input.c_str();

	// read flag of do interpolate from the input file
	interp_on = false;
	GetPrivateProfileString("Interpolation", "do interpolate", "<<missing>>", charBuffer, CHARLEN, input_file);
	if (string(charBuffer) == "<<missing>>")
		logstr("Property \"do interpolate\" is not speciefied in the section \"Interpolation\", default to be No.\n");
	if (string(charBuffer) == "Yes")
		interp_on = true;
	else if(string(charBuffer) != "No")
		logstr("WARNING: Default is not to do results interpolation.\n");
	if (interp_on)
		logstr("  results interpolate . . . . Yes\n");
	else
		logstr("  results interpolate . . . . No\n");
	

	if(interp_on){
	//if(!create_on && interp_on){
		// read interpolation starting point from the input file
		GetPrivateProfileString("Interpolation", "starting point", "<<missing>>", charBuffer, CHARLEN, input_file);
		if (string(charBuffer) == "<<missing>>")
			error(2011, "Property \"starting point\" is not speciefied in the section \"Interpolation\" in the input file \"%s\".\n", input_file);
		if (sscanf(charBuffer, "%f %f %f", &p1, &p2, &p3) != 3)
			error(2011, "Could not parse the starting point of interpolation \"%s\" in the input file \"%s\".\n", charBuffer, input_file);
		interp.start = make_float3(p1, p2, p3);
		logstr("  CT starting point . . . . . %g, %g, %g\n", interp.start.x, interp.start.y, interp.start.z);

		// read interpolation dimensions from the input file
		GetPrivateProfileString("Interpolation", "xyz dimensions", "<<missing>>", charBuffer, CHARLEN, input_file);
		if (string(charBuffer) == "<<missing>>")
			error(2011, "Property \"xyz dimensions\" is not speciefied in the section \"Interpolation\" in the input file \"%s\".\n", input_file);
		int d1,d2,d3;
		if (sscanf(charBuffer, "%d %d %d", &d1, &d2, &d3) != 3)
			error(2011, "Could not parse the xyz dimensions of interpolation \"%s\" in the input file \"%s\".\n", charBuffer, input_file);
		interp.dimen = make_int3(d1, d2, d3);
		logstr("  CT xyz dimensions . . . . . %d, %d, %d\n", interp.dimen.x, interp.dimen.y, interp.dimen.z);

		// read interpolation 3d voxel sizes from the input file
		GetPrivateProfileString("Interpolation", "3d voxel sizes", "<<missing>>", charBuffer, CHARLEN, input_file);
		if (string(charBuffer) == "<<missing>>")
			error(2011, "Property \"3d voxel sizes\" is not speciefied in the section \"Interpolation\" in the input file \"%s\".\n", input_file);
		if (sscanf(charBuffer, "%f %f %f", &p1, &p2, &p3) != 3)
			error(2011, "Could not parse the 3d voxel sizes of interpolation \"%s\" in the input file \"%s\".\n", charBuffer, input_file);
		interp.sizes = make_float3(p1, p2, p3);
		logstr("  CT 3d voxel sizes . . . . . %g, %g, %g\n", interp.sizes.x, interp.sizes.y, interp.sizes.z);
		interp.start.x += interp.sizes.x/2;
		interp.start.y += interp.sizes.y/2;
		interp.start.z += interp.sizes.z/2;
	}
	

	// read flag of output absolute dose from the input file
	AbsDos_on = false;
	GetPrivateProfileString("Dose output", "absolute dose", "<<missing>>", charBuffer, CHARLEN, input_file);
	if (string(charBuffer) == "<<missing>>")
		logstr("Property \"absolute dose\" is not speciefied in the section \"Dose output\", default to be No.\n");
	if (string(charBuffer) == "Yes")
		AbsDos_on = true;
	else if(string(charBuffer) != "No")
		logstr("WARNING: Default is not to output absolute dose.\n");
	if (AbsDos_on)
		logstr("  output absolute dose  . . . Yes\n");
	else
		logstr("  output absolute dose  . . . No\n");

	if(AbsDos_on){
		// read dose position of absolute dose from the input file
		GetPrivateProfileString("Dose output", "dose position", "<<missing>>", charBuffer, CHARLEN, input_file);
		if (string(charBuffer) == "<<missing>>")
			error(2011, "Property \"dose position\" is not speciefied in the section \"Dose output\" in the input file \"%s\".\n", input_file);
		if (sscanf(charBuffer, "%f %f %f", &p1, &p2, &p3) != 3)
			error(2011, "Could not parse the dose position of absolute dose \"%s\" in the input file \"%s\".\n", charBuffer, input_file);
		dose_pos = make_float3(p1, p2, p3);
		logstr("    dose position . . . . . . %.2f, %.2f, %.2f\n", dose_pos.x, dose_pos.y, dose_pos.z);

		// read dose measured of absolute dose from the input file
		GetPrivateProfileString("Dose output", "dose measured", "<<missing>>", charBuffer, CHARLEN, input_file);
		if (string(charBuffer) == "<<missing>>")
			error(2011, "Property \"dose measured\" is not speciefied in the section \"Dose output\" in the input file \"%s\".\n", input_file);
		if (sscanf(charBuffer, "%f", &p1) != 1)
			error(2011, "Could not parse the dose measured of absolute dose \"%s\" in the input file \"%s\".\n", charBuffer, input_file);
		dose_meas = p1;
		logstr("    dose measured . . . . . . %.1f\n", dose_meas);
	}

	// read flag of output relative dose from the input file
	RelDos_on = false;
	GetPrivateProfileString("Dose output", "relative dose", "<<missing>>", charBuffer, CHARLEN, input_file);
	if (string(charBuffer) == "<<missing>>")
		logstr("Property \"relative dose\" is not speciefied in the section \"Dose output\", default to be No.\n");
	if (string(charBuffer) == "Yes")
		RelDos_on = true;
	else if(string(charBuffer) != "No")
		logstr("WARNING: Default is not to output relative dose.\n");
	if (RelDos_on)
		logstr("  output relative dose  . . . Yes\n");
	else
		logstr("  output relative dose  . . . No\n");

	string meas_data_file_input;
	if(RelDos_on){
		// read starting point of relative dose from the input file
		GetPrivateProfileString("Dose output", "starting point", "<<missing>>", charBuffer, CHARLEN, input_file);
		if (string(charBuffer) == "<<missing>>")
			error(2011, "Property \"starting point\" is not speciefied in the section \"Dose output\" in the input file \"%s\".\n", input_file);
		if (sscanf(charBuffer, "%f %f %f", &p1, &p2, &p3) != 3)
			error(2011, "Could not parse the starting point of relative dose \"%s\" in the input file \"%s\".\n", charBuffer, input_file);
		RelDos.start = make_float3(p1, p2, p3);
		RelDos.start.x += h_isocenter_location.x;
		RelDos.start.y += h_isocenter_location.y;
		RelDos.start.z += h_isocenter_location.z;
		logstr("    starting point  . . . . . %.2f, %.2f, %.2f\n", RelDos.start.x, RelDos.start.y, RelDos.start.z);

		// read dimensions of relative dose from the input file
		GetPrivateProfileString("Dose output", "xyz dimensions", "<<missing>>", charBuffer, CHARLEN, input_file);
		if (string(charBuffer) == "<<missing>>")
			error(2011, "Property \"xyz dimensions\" is not speciefied in the section \"Dose output\" in the input file \"%s\".\n", input_file);
		int d1,d2,d3;
		if (sscanf(charBuffer, "%d %d %d", &d1, &d2, &d3) != 3)
			error(2011, "Could not parse the xyz dimensions of relative dose \"%s\" in the input file \"%s\".\n", charBuffer, input_file);
		RelDos.dimen = make_int3(d1, d2, d3);
		logstr("    xyz dimensions  . . . . . %d, %d, %d\n", RelDos.dimen.x, RelDos.dimen.y, RelDos.dimen.z);

		// read voxel sizes of relative dose from the input file
		GetPrivateProfileString("Dose output", "3d voxel sizes", "<<missing>>", charBuffer, CHARLEN, input_file);
		if (string(charBuffer) == "<<missing>>")
			error(2011, "Property \"3d voxel sizes\" is not speciefied in the section \"Dose output\" in the input file \"%s\".\n", input_file);
		if (sscanf(charBuffer, "%f %f %f", &p1, &p2, &p3) != 3)
			error(2011, "Could not parse the 3d voxel sizes of relative dose \"%s\" in the input file \"%s\".\n", charBuffer, input_file);
		RelDos.sizes = make_float3(p1, p2, p3);
		logstr("    3d voxel sizes  . . . . . %.2f, %.2f, %.2f\n", RelDos.sizes.x, RelDos.sizes.y, RelDos.sizes.z);

		// read ion chamber array type for relative dose measurement from the input file
		array_type = 0;
		GetPrivateProfileString("Dose output", "PTW729/PTW1000", "<<missing>>", charBuffer, CHARLEN, input_file);
		if (string(charBuffer) == "<<missing>>")
			logstr("Property \"PTW729/PTW1000\" is not speciefied in the section \"Dose output\".\n");
		if (string(charBuffer) == "PTW729")
			array_type = 1;
		else if(string(charBuffer) == "PTW1000")
			array_type = 2;
		else
			logstr("WARNING: Property \"PTW729/PTW1000\" is not properly speciefied in the section \"Dose output\".\n");
		if (array_type == 1)
			logstr("  ion chamber array . . . . . PTW729\n");
		else if(array_type == 2)
			logstr("  ion chamber array . . . . . PTW1000\n");

		// read normalization type for relative dose measurement from the input file
		norm_type = 0;
		GetPrivateProfileString("Dose output", "norm to Max/Center", "<<missing>>", charBuffer, CHARLEN, input_file);
		if (string(charBuffer) == "<<missing>>")
			logstr("Property \"norm to Max/Center\" is not speciefied in the section \"Dose output\".\n");
		if (string(charBuffer) == "Max")
			norm_type = 1;
		else if(string(charBuffer) == "Center")
			norm_type = 2;
		else
			logstr("  normalization . . . . . . . No\n");
		if (norm_type == 1)
			logstr("  normalization to  . . . . . Maximum dose\n");
		else if(norm_type == 2)
			logstr("  normalization to  . . . . . Center point\n");

		// read gamma index parameters from the input file
		GetPrivateProfileString("Dose output", "gamma index (mm,%)", "<<missing>>", charBuffer, CHARLEN, input_file);
		if (string(charBuffer) == "<<missing>>")
			error(2011, "Property \"gamma index (mm,%)\" is not speciefied in the section \"Dose output\" in the input file \"%s\".\n", input_file);
		if (sscanf(charBuffer, "%f %f", &p2, &p3) != 2)
			error(2011, "Could not parse the gamma index paremeters \"%s\" in the input file \"%s\".\n", charBuffer, input_file);
		GmIdx_len = p2;
		GmIdx_dos = p3;
		logstr("  gamma index paremters . . . %d mm, %d %%\n", (int)p2, (int)p3);
		GmIdx_len *= 0.1F;   //change from mm to cm
		GmIdx_dos *= 0.01F;  //change from % to float

		// get measured data file
		GetPrivateProfileString("Dose output", "measured data", "<<missing>>", charBuffer, CHARLEN, input_file);
		meas_data_file_input = string(charBuffer);
		if (meas_data_file_input == "<<missing>>")
			error(2012, "Property \"measured data\" is not speciefied in the section \"Dose output\" in the input file \"%s\".\n", input_file);
		meas_data_file_input = data_dir_input + meas_data_file_input;
		logstr("  measured data . . . . . . . %s\n", meas_data_file_input.c_str());

	}
	meas_data_file = meas_data_file_input.c_str();

    // write settings
    logstr("\nSettings\n");
    logstr("  Warps per block . . . . . . %d\n", SIMULATION_WARPS_PER_BLOCK);
	logstr("  Blocks per processor  . . . %d\n", SIMULATION_BLOCKS_PER_MULTIPROC);

	// perform initialization
    init(default_medium);

	logstr("  Total number of MU  . . . . %d\n", total_MU);
	logstr("  Histories per MU  . . . . . %I64u\n\n", histories_per_MU);

#ifndef GET_UNCERTAINTY
	//output the total dose
	score_t *h_eng_score = (score_t*)malloc(h_size_phantom*sizeof(score_t));
	memset(h_eng_score, 0, h_size_phantom*sizeof(score_t));
	score_t *h_score_tmp = (score_t*)malloc(h_size_phantom*sizeof(score_t));
#endif

	//now run simulation for each fields
	for(int iField = 0; iField<N_Field; iField++)
	{
		runSimulation(histories_per_MU,iField);

#ifndef GET_UNCERTAINTY
		float output_factor;
		if(outfac_on) output_factor = field[iField].output_factor;

		for(int GPUId=0; GPUId<GPUNo; GPUId++) {
#ifdef USE_MULTIPLE_GPU
			cudaSetDevice(GPUId); ce(58015);
#endif
			//cudaMemcpyFromSymbol(&d_eng_score,eng_score,sizeof(score_t*)); ce(110181);
			cudaMemcpy(h_score_tmp, d_eng_score[GPUId], h_size_phantom*sizeof(score_t), cudaMemcpyDeviceToHost); ce(11018);
			if(outfac_on)
				for(int i=0; i<h_size_phantom; i++)
					h_eng_score[i] += h_score_tmp[i]*output_factor;
			else
				for(int i=0; i<h_size_phantom; i++)
					h_eng_score[i] += h_score_tmp[i];
		}
#else
		if(iField==0){
			if(N_Field>1) logstr("The uncertainty of only the first field is output.\n");
			write_uncertainty(rootName,eng_dep2);
		}
#endif
	}

#ifndef GET_UNCERTAINTY
	write_dose(rootName,h_eng_score);
	free(h_eng_score);
	free(h_score_tmp);
#else
	write_dose(rootName,eng_dep);
#endif

	free_all();

	for(int GPUId=0; GPUId<GPUNo; GPUId++) {
#ifdef USE_MULTIPLE_GPU
		cudaSetDevice(GPUId); ce(58016);
#endif
		cudaDeviceReset();
	}

	clock_t tac = clock();
	float total_cpu_time = (float)(tac - tic) / (float)CLOCKS_PER_SEC * 1000.0F;
	logstr("\nTotal CPU/GPU . . . . . . . . %.2f ms (%.2f min)\n", total_cpu_time, total_cpu_time/6.0E+4F);

    // log time
	getTime();
	logstr("\n\nCUDA EGS with MLC finished at %s\n", charBuffer);

    fclose(logfile);

    return 0;
}

