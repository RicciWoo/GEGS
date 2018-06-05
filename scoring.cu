
/*__device__ void escore(uint region, float energy)
{
	score_t tmp = energy;
	if(region <3) return; //if the energy deposit other than the phantom regions, return.
	float density = region_data[region].rho;
	//if(density>0.9) tmp=tmp/density; //check here !! 
	tmp=tmp/density;
	atomicAdd(&eng_score[region-3], tmp);
	
}
*/

__device__ void escore(uint region, float energy)
{
	score_t tmp = energy;
	if(region < vac_and_oth) return; //if the energy deposit other than the phantom regions, return.
	//if(region < 5) return; //if the energy deposit other than the phantom regions, return.
	//if(region < 5)
		//atomicAdd(&eng_score[region], tmp);
	//else
		//atomicAdd(&eng_score[region-5], tmp);
		atomicAdd(&eng_score[region-vac_and_oth], tmp);
	//if(region==278359) printf("GPUId=%d, energy=%f\n",GPUId,tmp);
	//atomicAdd(&eng_score[region-5], 0.1);
}