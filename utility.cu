/* * * * * * * * * * * * * * *
 * General Helper Functions  *
 * * * * * * * * * * * * * * */

// check whether two floats are almost equal
// taken from http://www.cygnus-software.com/papers/comparingfloats/comparingfloats.htm

// the maximum number of ULPs that two floats may differ and still be considered "almost equal"
#define MAX_ULPS 20

__device__ bool almostEqual(float A, float B) {
    // Make sure maxUlps is non-negative and small enough that the
    // default NAN won't compare as equal to anything.

    // Make aInt lexicographically ordered as a twos-complement int
    int aInt = *(int*)&A;
    if (aInt < 0)
        aInt = 0x80000000 - aInt;

    // Make bInt lexicographically ordered as a twos-complement int
    int bInt = *(int*)&B;
    if (bInt < 0)
        bInt = 0x80000000 - bInt;

    int intDiff = abs(aInt - bInt);
    if (intDiff <= MAX_ULPS)
        return true;

    return false;
}

// calculate indices of this thread
__device__ indices get_indices() {
    indices idx;
    // index of the block in the grid
    idx.b = blockIdx.y * gridDim.x + blockIdx.x;
    
    // index of the particle on the stack
    idx.p = idx.b * blockDim.x + threadIdx.x;
    
    // index of the warp in the block
    idx.w = threadIdx.x / WARP_SIZE;
    
    // index of the thread in the warp
    idx.t = threadIdx.x % WARP_SIZE;
	
    return idx;
}


//this function seams to have problem, need to be checked. Tong
__device__ double atomicAdd(double* address, double val)
{
	unsigned long long int* address_as_ull = (unsigned long long int*)address;
	unsigned long long int old = *address_as_ull, assumed;
	do {
			assumed = old;
			old = atomicCAS(address_as_ull, assumed,	__double_as_longlong(val +	__longlong_as_double(assumed)));
	} while (assumed != old);
	return __longlong_as_double(old);
}



// This is the subroutine UPHI(IENTRY,LVL) with IENTRY = 2 and LVL = 1 in the file 
// egsnrc.mortran (v 1.72 2011/05/05) of the EGSnrc code system.
// However, note that we are not using the box method implemented in the macro 
// $SELECT-AZIMUTHAL-ANGLE, because that involves a sampling loop and then all
// threads would have to wait until the last thread has finished the loop. Calculating
// sin and cos with __sincosf is not that expensive, so we use that instead.
// Note that this is based on assumption and was not experimentally verified.
// This is the subroutine UPHI(IENTRY,LVL) with IENTRY = 2 and LVL = 1 in the file 
// egsnrc.mortran (v 1.72 2011/05/05) of the EGSnrc code system.
// However, note that we are not using the box method implemented in the macro 
// $SELECT-AZIMUTHAL-ANGLE, because that involves a sampling loop and then all
// threads would have to wait until the last thread has finished the loop. Calculating
// sin and cos with __sincosf is not that expensive, so we use that instead.
// Note that this is based on assumption and was not experimentally verified.

//Additions done by Tong Xu:
// The following function is equivalent to UPHI(2,1)
// It is called to set the direction vector for the FIRST secondary particle.
// incident particle --> two secondary particles.
// The direction of the secondary particle was sampled first in the frame at which the incident paricle 
// always go along z direction.The angles to be sampled are theta and PHI. Then we rotate the direction
// to world frame
// We now add sinphi and cosphi as two returned parameters to get the PHI angle result back.
// This is because, in many case, there are two secondary particles, and the second one should have the 
// same PHI as the first. The information of PHI was originally passed by COMMON/UPHIOT/ 
// however, in GPU we can't use them as __device__ globels because we have many threads running at the same time.
// So we have to pass the sinphi and cosphi as parameter back to the calling function. Then the uphi32 will be called 
// with known sinphi and cosphi.

__device__ void uphi21(float costhe, float sinthe, float &cosphi, float &sinphi, particle_t &p, curandState *randStat_ptr) {

    //if (p.process) {
		float r1 = curand_uniform(randStat_ptr);
        float phi = 2.0F * 3.14159265F * r1;
        //float cosphi, sinphi;
        __sincosf(phi, &sinphi, &cosphi);

		float A, B, C;
		A = p.u; B = p.v; C = p.w;
        float sinps2 = A*A + B*B;
        // small polar angle
        if (sinps2 < SMALL_POLAR_ANGLE_THRESHOLD) {
            p.u = sinthe * cosphi;
            p.v = sinthe * sinphi;
            p.w = p.w * costhe;
        }
        else {
            float sinpsi_i = rsqrtf(sinps2);
            //float sinpsi = sqrtf(sinps2);
            float us = sinthe * cosphi;
            float vs = sinthe * sinphi;
            //float sindel = B / sinpsi;
            //float cosdel = A / sinpsi;
            float sindel = B * sinpsi_i;
            float cosdel = A * sinpsi_i;
            
            p.u = C * cosdel * us - sindel * vs + A * costhe; 
            p.v = C * sindel * us + cosdel * vs + B * costhe;
            //p.w = -sinpsi * us + C * costhe;
			p.w = - us/sinpsi_i + C * costhe;
        }
    //}
}

/*
 The following function is equivalent to UPHI(3,2)
 It is called to set the direction vector for the SECOND secondary particle.
 incident particle --> two secondary particles.
 The direction of the secondary particle was sampled first in the frame at which the incident paricle 
 always go along z direction.The angles to be sampled are theta and PHI. Then we rotate the direction
 to world frame.

 The scatter theta angle of the secondary particles has already been calculated before calling UPHI.  The theta is
 And the azimuthal scattering angle PHI has already sampled when calling UPHI(2,1) for the first secondary particle.
 So we expect to have a sinphi and cosphi as inputs rather than sample PHI.
 Also we asume that the direction of the original particle is passed in  through the "particle_t &p".
 We just need to rotate the secondary particle direction to the world frame.
 */

__device__ void uphi32( float costhe, float sinthe, float cosphi, float sinphi, particle_t &p) {

    //if (p.process) {
		
		//we already have PHI, so just do the rotation.
		//p carries the incident particle directions
		float A, B, C;
		A = p.u; B = p.v; C = p.w;
        float sinps2 = A*A + B*B;
        // small polar angle
        if (sinps2 < SMALL_POLAR_ANGLE_THRESHOLD) {
            p.u = sinthe * cosphi;
            p.v = sinthe * sinphi;
            p.w = p.w * costhe;
        }
        else {
            float sinpsi_i = rsqrtf(sinps2);
            //float sinpsi = sqrtf(sinps2);
            float us = sinthe * cosphi;
            float vs = sinthe * sinphi;
            //float sindel = B / sinpsi;
            //float cosdel = A / sinpsi;
            float sindel = B * sinpsi_i;
            float cosdel = A * sinpsi_i;
            
            p.u = C * cosdel * us - sindel * vs + A * costhe; 
            p.v = C * sindel * us + cosdel * vs + B * costhe;
            //p.w = -sinpsi * us + C * costhe;
			p.w = - us/sinpsi_i + C * costhe;
        }
    //}
	//when return, the p will carry the secondary particle direction vector
}

// The following function is equivalent to UPHI(1,1)
// it assume that the theta (rather than sinthe and costhe) is already sampled before calling.
// The only different of uphi(1,1) as compared to uphi(2,1) is that it calcualte sinthe and costhe here. 
__device__ void uphi11(float theta, float &costhe, float &sinthe, float &cosphi, float &sinphi, particle_t &p, curandState *randStat_ptr)
{
	//if (p.process) {
		float r1 =curand_uniform(randStat_ptr);
        __sincosf(theta, &sinthe, &costhe);
		//sinthe=sin(theta);
		//float cthet=pi5d2-theta;
		//costhe=sin(cthet);
        float phi = 2.0F * 3.14159265F * r1;
        //float cosphi, sinphi;
        __sincosf(phi, &sinphi, &cosphi);

		float A, B, C;
		A = p.u; B = p.v; C = p.w;
        float sinps2 = A*A + B*B;
        // small polar angle
        if (sinps2 < SMALL_POLAR_ANGLE_THRESHOLD) {
            p.u = sinthe * cosphi;
            p.v = sinthe * sinphi;
            p.w = p.w * costhe;
        }
        else {
            float sinpsi_i = rsqrtf(sinps2);
            //float sinpsi = sqrtf(sinps2);
            float us = sinthe * cosphi;
            float vs = sinthe * sinphi;
            //float sindel = B / sinpsi;
            //float cosdel = A / sinpsi;
            float sindel = B * sinpsi_i;
            float cosdel = A * sinpsi_i;
            
            p.u = C * cosdel * us - sindel * vs + A * costhe; 
            p.v = C * sindel * us + cosdel * vs + B * costhe;
            //p.w = -sinpsi * us + C * costhe;
			p.w = - us/sinpsi_i + C * costhe;
        }
    //}
}

//load a particle from the simulation stack. if there is no more, set simu_stack_is_empty to true.
//return true if successfuly load a particle.
__device__ bool pop_stack(particle_t &p, uchar status)
{
	int np_local = atomicSub(&np[status],1);
	if(np_local<0)
	{
		np[status] = -1;
		simu_stack_is_empty = true;  //means simulation stack is empty
		return false;
	}

	/*
    //read particle from stack
    uint4 tmp = stack[status].a[np_local];
	p.status = ((uchar*)&tmp.x)[0];
    p.nRepeat = ((uchar*)&tmp.x)[1];
    p.charge = ((char*)&tmp.x)[2];
    p.process = ((uchar*)&tmp.x)[3];
    p.e = *(float*)&tmp.y;
    p.wt = *(float*)&tmp.z;
    p.region = tmp.w;
    
    tmp = stack[status].b[np_local];
	p.leafIndex = ((char*)&tmp.x)[0];  //modified by Wu, 20130324
	p.bank = ((char*)&tmp.x)[1];
    p.latch = ((ushort*)&tmp.x)[1];
    p.x = *(float*)&tmp.y;
    p.y = *(float*)&tmp.z;
    p.z = *(float*)&tmp.w;

    tmp = stack[status].c[np_local];
    p.u = *(float*)&tmp.x;
    p.v = *(float*)&tmp.y;
    p.w = *(float*)&tmp.z;
	p.dnear= *(float*)&tmp.w;  
	*/

    //read particle from stack
	p = stack[status][np_local];

	return true;
}

//save a particle to the stack whose index is status.
//return true if successfuly save a particle.
__device__ bool  push_stack(particle_t &p, uchar status)
{
    int np_local = atomicAdd(&np[status],1);
	np_local++;   //We add 1 here, because the atomicAdd only return the old value before atomicAdd
	if(np_local>=np_max[status])
	{
		np[status]=np_max[status]-1;
		printf("Error: stack[%d] overflow!!",status);
		return false;
	}

	/*
	uint4 tmp;
	((uchar*)&tmp.x)[0] = p.status;
    ((uchar*)&tmp.x)[1] = p.nRepeat;
    ((char*)&tmp.x)[2] = p.charge;
    ((uchar*)&tmp.x)[3] = p.process;
    tmp.y = *(uint*)&p.e;
    tmp.z = *(uint*)&p.wt;
    tmp.w = p.region;
    stack[status].a[np_local]=tmp;
   
	((char*)&tmp.x)[0] = p.leafIndex;  //modified by Wu, 20130323
	((char*)&tmp.x)[1] = p.bank;  //modified by Wu, 20130323
    ((ushort*)&tmp.x)[1] = p.latch;
    tmp.y = *(uint*)&p.x;
    tmp.z = *(uint*)&p.y;
    tmp.w = *(uint*)&p.z;
    stack[status].b[np_local]=tmp;

    tmp.x = *(uint*)&p.u;
    tmp.y = *(uint*)&p.v;
    tmp.z = *(uint*)&p.w;
	tmp.w = *(uint*)&p.dnear;
    stack[status].c[np_local]=tmp;
	*/

	stack[status][np_local]=p;

	return true;
}
