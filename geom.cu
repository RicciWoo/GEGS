//declaration of geometry functions, added by Wu 20130125
//__device__ bool pointInPolygon(point2D p, point2D *polygon,int nVertex);
//__device__ uchar pointInLeaf(point2D p, point2D *endnode);
//__device__ char  indexOfLeaf(point2D p);
//__device__ void leaf_position(point2D *endnode, char index);
//__device__ float rayIntersectSegment(point2D p0, point2D v, point2D *segment);
//__device__ float rayIntersectPolygon(point2D p0, point2D v, point2D *polygon, int nVertex, bool & insdie);
//__device__ uint  inwhere_MLC(particle_t &q);
//__device__ uint  howfar_MLC(particle_t &p, float &t);
//__device__ uint  howfar_phantom(particle_t &p, float &t,int &ix, int &iy, int &iz);

/* * * * * * * * * * * *
 * Geometry Functions  *
 * * * * * * * * * * * */

// Determine the index i such that the point p lies between bounds[i] and bounds[i+1].
// Code was taken from the function isWhere of the class EGS_PlanesT in the file
// egs_planes.h (v 1.17 2009/07/06) and the function findRegion of the class EGS_BaseGeometry
// in the file egs_base_geometry.h (v 1.26 2008/09/22) of the EGSnrc C++ Class Library.
__device__ int isWhere(float p, uint nreg, float *bounds) {
    if ((p < bounds[0]) || (p > bounds[nreg]))
        return -1;
    if (nreg == 1)
        return 0;

    int ml = 0;
    int mu = nreg;
    while (mu - ml > 1) {
        int mav = (ml + mu) / 2;
        if (p <= bounds[mav])
            mu = mav; 
        else 
            ml = mav;
    }
    return  mu - 1;
}

/*
__device__ uint howfar(particle_t &p, float &t)
{
	if(p.leafIndex == PASSMLC)
		return howfar_phantom(p, t);
	else
		return howfar_MLC(p, t);
}
*/
// Determine the distance t to the next voxel boundary for the particle p and return
// the region index that the particle will enter.
// Code was taken from the function howfar of the class EGS_XYZGeometry in the file
// egs_nd_geometry.h (v 1.26 2009/07/06) of the EGSnrc C++ Class Library.
// Our region indices are shifted by 1 because our outside region is 0, while
// the outside region in the EGSnrc C++ Class Library is -1.
//__device__ uint howfar_phantom(particle_t &p0, float &t, int &ix, int &iy, int &iz) {
__device__ uint howfar_phantom(particle_t &p0, float &t) {
  			  particle_t p;
#ifndef 	AdjustGeomPASSMLC   //if the particle geometry is not adjust right after PASSMLC, then adjust it every time when doing geometry in Phantom
				//p.bank=p0.bank;
				//p.leafIndex=p0.leafIndex;
				if(Change_xyz){
					/*
					p.x = p0.x;
					//p.y =p0.z;
					p.y =p0.z - SourceToIsocenter;
					p.z =-p0.y; //switch the z and y
					p.u=p0.u;
					p.v= p0.w; //switch the v and w
					p.w=-p0.v;
					*/
					p.x =   p0.y;
					p.y =   (p0.z - SourceToIsocenter);
					p.z =   p0.x;
					p.u =   p0.v;
					p.v =   p0.w;
					p.w =   p0.u;
				}
				else{
					p.x = p0.x;
					p.y = p0.y;
					p.z = p0.z - SourceToIsocenter;
					p.u = p0.u;
					p.v = p0.v;
					p.w = p0.w;
				}
				p.region = p0.region;


//now rotate the incident particles by gantry angle
				float3 tmp;

				//tmp.x=   p.x*cosAngle + p.y*sinAngle;
				//tmp.y= - p.x*sinAngle + p.y*cosAngle;
				tmp.x = p.x*cosAngle - p.y*sinAngle;
				tmp.y = p.x*sinAngle + p.y*cosAngle;

				p.x=tmp.x;
				p.y=tmp.y;

				//tmp.x=   p.u*cosAngle + p.v*sinAngle;
				//tmp.y= - p.u*sinAngle + p.v*cosAngle;
				tmp.x = p.u*cosAngle - p.v*sinAngle;
				tmp.y = p.u*sinAngle + p.v*cosAngle;

				p.u=tmp.x;
				p.v=tmp.y;

				//now shift
				p.x = p.x + isocenter_location.x;
				p.y = p.y + isocenter_location.y;
				p.z = p.z + isocenter_location.z;
#else
			  p=p0;  //this seams able to copy the contents of p0 to p, no need to do it individually, by Tong Xu May 2013
#endif
/*
		if(p.leafIndex ==PASSMLC)
		{
			  //first time outMLC, check if it already inside Phantom Geometry
				int ix = isWhere(p.x, phantom_N.x, phantom_x_bounds);
				int iy = isWhere(p.y, phantom_N.y, phantom_y_bounds);
				int iz = isWhere(p.z, phantom_N.z, phantom_z_bounds);
				if(ix>=0 && iy>=0 && iz>=0)
				{
					 p.region=ix + iy * phantom_N.x + iz * phantom_N.x * phantom_N.y + 1 + 2;
				}
				p.leafIndex = PHANTOM_STARTED;
		}
	*/	

	//if (p.region > 4) {
	if (p.region >= vac_and_oth) {    //the particle is currently in the phantom, commented by Tong Xu, Jan 2013, remember the first Voxel in the phantom has region number =2. Because the first 2 regions are used for MLC.
        // because of the above mentioned shift, we have to substract 1
		// ir is the actual voxel index
		//uint ir = p.region - 1 - 2 - 2;  // " - 2" because the region index is shift by 2, used by the MLC regions.
		uint ir = p.region - vac_and_oth;  // " - 2" because the region index is shift by 2, used by the MLC regions.


        int iz = ir / (phantom_N.x * phantom_N.y); 
        ir -= iz * phantom_N.x * phantom_N.y; 
        int iy = ir / phantom_N.x;
        int ix = ir - iy * phantom_N.x;

		uint inew = p.region ;

        if (p.u > 0.0F) {
            float d = (phantom_x_bounds[ix + 1] - p.x) / p.u;
            if (d <= t) { 
                t = d; 
                if (++ix  < phantom_N.x) 
                    inew = p.region + 1; 
                else 
					//inew = 0;
                    inew = vacuum_reg;  //going out the phantom, commented by Tong Xu, Jan 2013
            }
        }
        else if (p.u < 0.0F) {
            float d = (phantom_x_bounds[ix] - p.x) / p.u;
            if (d <= t) { 
                t = d; 
                if (--ix >= 0) 
                    inew = p.region - 1; 
                else 
                    //inew = 0;
					inew = vacuum_reg;
            }
        }

        if (p.v > 0.0F) {
            float d = (phantom_y_bounds[iy + 1] - p.y) / p.v;
            if (d <= t) { 
                t = d; 
                if (++iy < phantom_N.y) 
                    inew = p.region + phantom_N.x; 
                else 
                    //inew = 0;
					inew = vacuum_reg;
            }
        }
        else if (p.v < 0.0F) {
            float d = (phantom_y_bounds[iy] - p.y) / p.v;
            if (d <= t) { 
                t = d; 
                if (--iy  >= 0) 
                    inew = p.region - phantom_N.x; 
                else 
                    //inew = 0;
					inew = vacuum_reg;
            }
        }

        if (p.w > 0.0F) {
            float d = (phantom_z_bounds[iz + 1] - p.z) / p.w;
            if (d <= t) { 
                t = d; 
                if (++iz < phantom_N.z) 
                    inew = p.region + phantom_N.x * phantom_N.y; 
                else 
                    //inew = 0;
					inew = vacuum_reg;
            }
        }
        else if (p.w < 0.0F) {
            float d = (phantom_z_bounds[iz] - p.z) / p.w;
            if (d <= t) { 
                t = d; 
                if (--iz  >= 0) 
                    inew = p.region - phantom_N.x * phantom_N.y; 
                else 
                    //inew = 0;
					inew = vacuum_reg;
            }
        }

        if(t<0)
        {
                t=EPSTEP;  //added by Tong Xu, 20131120
        }

        return inew ;
    }
    // this part corresponds to the function howfarFromOut of the class EGS_XYZGeometry 
    // in the file egs_nd_geometry.h (v 1.26 2009/07/06) of the EGSnrc C++ Class Library
	// if p.region =0, means currently in in VACUUM,  i.e. outside of phantom. check if the particle will go into phantom
    else {
        int ix, iy, iz;
        float t1;

        ix = -1;
        if ((p.x <= phantom_x_bounds[0]) && (p.u > 0.0F)) {     //on the -x side of the phantom, but  going to the +x direction, commented by Tong Xu, 2013 Jan
            t1 = (phantom_x_bounds[0] - p.x) / p.u; 
            ix = 0;
        }
        else if ((p.x >= phantom_x_bounds[phantom_N.x]) && (p.u < 0.0F)) {  //on the +x side of the phantom, but  going to the -x direction. commented by Tong Xu, 2013 Jan
            t1 = (phantom_x_bounds[phantom_N.x] - p.x) / p.u; 
            ix = phantom_N.x - 1;
        }
        
        if ((ix >= 0) && (t1 <= t)) {
            float y1 = p.y + p.v * t1;
            iy = isWhere(y1, phantom_N.y, phantom_y_bounds);     // now find if the paticle will intersect with x boundaries, where y, and z it will land on, and is it within the phantom's y and z range? commented by Tong Xu, 2013 Jan
            
            if (iy >= 0) {
                float z1 = p.z + p.w * t1;
                iz = isWhere(z1, phantom_N.z, phantom_z_bounds);
                
                if (iz >= 0) {
                    t = t1; 
					//return ix + iy * phantom_N.x + iz * phantom_N.x * phantom_N.y + 1 + 2 + 2;
                    return ix + iy * phantom_N.x + iz * phantom_N.x * phantom_N.y + vac_and_oth;             // if it will land on a voxel, step right to the edge and return the region number of that voxel. commented by Tong Xu, 2013 Jan
                }																																	//"+2" because the shift of region index caused by the region number used by MLC
            }
        }

        iy = -1;					//simillary, check the y boundaries , commented by Tong Xu, 2013 Jan
        if ((p.y <= phantom_y_bounds[0]) && (p.v > 0.0F)) {
            t1 = (phantom_y_bounds[0] - p.y) / p.v; 
            iy = 0;
        }
        else if ((p.y >= phantom_y_bounds[phantom_N.y]) && (p.v < 0.0F)) {
            t1 = (phantom_y_bounds[phantom_N.y] - p.y) / p.v; 
            iy = phantom_N.y - 1;
        }
        
        if ((iy >= 0) && (t1 <= t)) {
            float x1 = p.x + p.u * t1;
            ix = isWhere(x1, phantom_N.x, phantom_x_bounds);
            
            if (ix >= 0) {
                float z1 = p.z + p.w * t1;
                iz = isWhere(z1, phantom_N.z, phantom_z_bounds);
                
                if (iz >= 0) {
                    t = t1; 
					//return ix + iy * phantom_N.x + iz * phantom_N.x * phantom_N.y + 1 + 2 + 2;
                    return ix + iy * phantom_N.x + iz * phantom_N.x * phantom_N.y + vac_and_oth;
                }
            }
        }

        iz = -1;			//simillary, check the y boundaries , commented by Tong Xu, 2013 Jan
        if ((p.z <= phantom_z_bounds[0]) && (p.w > 0.0F)) {
            t1 = (phantom_z_bounds[0] - p.z) / p.w; 
            iz = 0;
        }
        else if ((p.z >= phantom_z_bounds[phantom_N.z]) && (p.w < 0.0F)) {
            t1 = (phantom_z_bounds[phantom_N.z] - p.z) / p.w; 
            iz = phantom_N.z - 1;
        }
        
        if ((iz >= 0) && (t1 <= t)) {
            float x1 = p.x + p.u * t1;
            ix = isWhere(x1, phantom_N.x, phantom_x_bounds);
            
            if (ix >= 0) {
                float y1 = p.y + p.v * t1;
                iy = isWhere(y1, phantom_N.y, phantom_y_bounds);
                
                if (iy >= 0) {
                    t = t1; 
					//return ix + iy * phantom_N.x + iz * phantom_N.x * phantom_N.y + 1 + 2 + 2;
                    return ix + iy * phantom_N.x + iz * phantom_N.x * phantom_N.y + vac_and_oth;
                }
            }
        }
        
		//return 0;
        return vacuum_reg;
    }
}

//definition of geometry functions, added by Wu 20130125
__device__ bool pointInPolygon(point2D p, point2D *polygon,int nVertex)
{

  int polySides = nVertex;
  int      i, j =polySides-1;
  bool  oddNodes=false      ;

  for (i=0; i<polySides; i++)
  {
    if (((polygon[i].y< p.y && polygon[j].y>=p.y)
    ||  (polygon[j].y< p.y && polygon[i].y>=p.y))
    &&  (polygon[i].x<=p.x || polygon[j].x<=p.x)) {
      oddNodes^=(polygon[i].x+(p.y-polygon[i].y)/(polygon[j].y-polygon[i].y)*(polygon[j].x-polygon[i].x)<p.x); }
    j=i;
  }

  return oddNodes;

}

__device__  uchar pointInLeaf(point2D p, point2D *endnode_loc)
{
    uchar bank=0; // by defalt in the field/air

    if(pointInPolygon(p,endnode_loc,17))
    {
        bank=1;
    }
    else if(pointInPolygon(p,(endnode_loc+17),17))
    {
        bank=2;
    }

    return bank;
}

__device__ char  indexOfLeaf(point2D p)
{
    char i=OUTMLC;
    char index=OUTMLC;    //by default set the index to OUTMLC (-99) to mean that the point is outside the zone of interest

    //float x=p.x/(p.y+MLC_zmin) *MLC_zmin;
	float x=p.x/p.y*MLC_zmin;
    //float z=p.y;
	float z=p.y-MLC_zmin;

    if(z<0 ||z> MLC_zthick )
        return index;
    if(x<MLC_x_start||x>MLC_x_end)
        return index;

    if(x<MLC_x_start+MLC_wt)
    {
        i=0;
    }
    else if(x<MLC_x_mid_1+MLC_wt)
    {
        i=int((x-MLC_x_start-MLC_wt)/(MLC_wl_full+MLC_air_gap));
    }
    else if(x<MLC_x_mid_2+MLC_wt)
    {
        i=int((x-MLC_x_mid_1-MLC_wt)/(MLC_wl_target+MLC_wl_isocenter+MLC_air_gap*2));
        if(x-MLC_x_mid_1-MLC_wt-(MLC_wl_target+MLC_wl_isocenter+MLC_air_gap*2)*i<MLC_wl_target+MLC_air_gap)
        {
            i=10+i*2;
        }
        else
        {
            i=11+i*2;
        }
    }
    else if(x<MLC_x_end-MLC_wl_full)
    {
        i=int((x-MLC_x_mid_2-MLC_wt)/(MLC_wl_full+MLC_air_gap))+50;
    }
    else if(x<=MLC_x_end)
    {
        i=58;
    }
    if(i>=0)
    {
        if(pointInPolygon(p,(crossXZ+i*16),14))
        {
            index=i;
        }
        else if (i<59)
		{
			if(pointInPolygon(p,(crossXZ+(i+1)*16),14))
			{
				index=i+1;
			}
			else
			{
				index=-(i+1);
			}
		}
		else index=OUTMLC;
    }

    return index;
}

__device__ float rayIntersectSegment(point2D p0, point2D v, point2D *segment)
{
    float Dist = 1e9;
    float x1,x2,y1,y2,x0,y0;
    float u;
    x0 = p0.x;
    y0 = p0.y;

    //the line segment can be expressed by x = x1 +u * (x2-x1), y=y1+u*(y2-y1)
    // u must be between [0,1]
    float tmp1, dist;
    x1 = segment[0].x;
    y1 = segment[0].y;
    x2 = segment[1].x;
    y2 = segment[1].y;
    //checking if the segment is parralle with ray
    tmp1 = v.x*(y2-y1)-v.y*(x2-x1);
    if(tmp1 != 0) // if not paralle
    {
        dist = ((x2-x1)*(y0-y1)-(y2-y1)*(x0-x1))/tmp1;
        u = ((x1-x0)*v.y-(y1-y0)*v.x)/tmp1;
        if(dist>=0 && (u>=0)&& (u<=1))
        {
            Dist=dist;
        }
    }

    return Dist;
}

__device__ float rayIntersectPolygon(point2D p0, point2D v, point2D *polygon, int nVertex, bool &inside )
{
    int i,j;
    float minDist = INFI_DIST+1;
    float x1,x2,y1,y2,x0,y0;
    float u;
    x0 = p0.x;
    y0 = p0.y;

    //the line segment can be expressed by x = x1 +u * (x2-x1), y=y1+u*(y2-y1)
    // u must be between [0,1]
    float tmp1, dist;
    int countIntersect =0;

    for (i=0;i<nVertex;i++)
    {
        x1 = polygon[i].x;
        y1 = polygon[i].y;
        j=i+1;
        if(j>=nVertex) j=0;
        x2 = polygon[j].x;
        y2 = polygon[j].y;
        //checking if the segment is parralle with ray
        tmp1 = v.x*(y2-y1)-v.y*(x2-x1);
        if(tmp1 != 0) // if not paralle
        {
            dist = ((x2-x1)*(y0-y1)-(y2-y1)*(x0-x1))/tmp1;
            u = ((x1-x0)*v.y-(y1-y0)*v.x)/tmp1;
            if(dist>=0 && (u>=0)&& (u<1)) // Changed "u<=1" to "u<1" to avoid double intersection when the ray passes the vertex.
            {
                countIntersect++;
                if(dist<minDist)
                    minDist=dist;
            }
        }
    }

	//if countIntersect is odd number, the point is inside the polygon
 inside = countIntersect & 1;  //check the lowest bit.
 return minDist;
//    if(countIntersect>0) return minDist;
//       else return INFI_DIST+1;
}

__device__ uint  inwhere_MLC(float3 &q, char & index, char &Bank)
{
 //   int index;
//	uint inew = 0;  //default in VACUUM
	//uint inew = 1;  //default in air
	uint inew = default_reg;  //default in air
	//point2D endnode[34];
    point2D p;
    p.x=q.x; p.y=q.z;
    if(index==OUTMLC) index=indexOfLeaf(p);  //now the Howfar_MLC will try to figure out the index by monitoring how the particle crossing the MLC polygons, 
	                                         //so we don't have to do indexOfLeaf every time.
    p.x=q.y; //p.y=q.z;
    if(index>=0)
    {
        //leaf_position(endnode,index);
		Bank = pointInLeaf(p,endnode[index]);
        if(Bank >0)   //Bank =1 : Bank A, Bank = 2, BANK B
		{
			//inew = 2;   //in Tungsten
			inew = MLC_region;   //in Tungsten
		}
		else
		{
			//inew = 1;    //in Air
			inew = default_reg;    //in Air
		}
    }
	
    return inew;
}


__device__ uint  howfar_MLC(particle_t &p0, float &t)
{

				// float zmin=crossXZ[15].x;            //cross[15].x=zmin;

				particle_t p;
				p=p0; 
#ifdef ROTATE_MLC_XY
				p.x = p0.y;
				p.y =p0.x;
#endif
				p.z =p0.z;
				//p.z =p0.z -(crossXZ[15].x -SourceToIsocenter); //shift the z, because the MLC geometry code is based on from surface of MLC is at z=0;
#ifdef ROTATE_MLC_XY
				p.u=p0.v;
				p.v=p0.u;
#endif

		//after the shift, the following code stays almost the same as the GPUMLC code ,expect the p0.index and p0.bank assignment  before return.

				uint inew = p.region;
				char bank = p.bank;
				char index = p.leafIndex;
				char indexLimit1=0, indexLimit2=60;

				float3 q;  //use to be particle_t. 
				if(p.y<-MLC_rmax||p.y>MLC_rmax)
					//return 0;  //outside area is VACUUM
					return vacuum_reg;  //outside area is VACUUM

/*
if(inew>=phantom_N.x*phantom_N.y*phantom_N.z)
	{
		inew=0;
	    p0.status = p_user_discard;
		return inew;
	}

*/
				float minDist=INFI_DIST+1; //set the minDist to very large
				float minBoundaryDist=INFI_DIST+1; //set the minDist to very large
				float dis_xz=0;
				float dis_yz=0;
				float temp1=INFI_DIST+1;
//				float temp2=INFI_DIST+1;
//				point2D endnode[34];
				point2D p1,p2,p3,dir1,dir2;

				bool inside_XZ,inside_YZ;

				p1.x=p.x; p1.y=p.z;
				p2.x=p.y; p2.y=p.z;
				dir1.x=p.u; dir1.y=p.w;    //we do not normalize the 2D direction vectors
				dir2.x=p.v; dir2.y=p.w;

				dis_xz=rayIntersectPolygon(p1,dir1,boundaryXZ,4,inside_XZ);
				dis_yz=rayIntersectPolygon(p2,dir2,boundaryYZ,4,inside_YZ);

				if(dis_xz >= INFI_DIST || dis_yz >= INFI_DIST)  //either INFI_DIST,  means the particle will not have intersection with the 3D boundary of the MLC comples,  
				{
					//inew = 0;
					inew = vacuum_reg;
					return inew;
				}

				if (dis_xz < dis_yz) minBoundaryDist = dis_xz;
				else minBoundaryDist = dis_yz;

				if(inside_XZ &&          //check if the particle is inside the MLC complex in the XZ plane
					          inside_YZ )   //check if the particle is inside the MLC complex in the YZ plane
				{
	  				p3.x = p1.x + dis_xz * dir1.x;  //the x coordinate of the point that the ray intersect with the XZ boundary
						                                                       // we use this to limit the search range of leaves for intersection 
					p3.y=p1.y+(dis_xz-0.001)*dir1.y;  //subtract a little to make sure it stays inside the boundary;
					indexLimit1=indexOfLeaf(p3); //find the limit of  search range for intersecting leaf
					if(indexLimit1 <0) indexLimit1 = -indexLimit1;
					
				} else                    // the particle is not inside the MLC complex boundaries
				{
					//inew = 0;  //defualt going into VACUUM
					inew = vacuum_reg;  //defualt going into VACUUM

					if (minBoundaryDist >= INFI_DIST)  //no intersection with the boundary
						return inew;

					if(t>minBoundaryDist)   
					{
						t = minBoundaryDist + EPSTEP ;   //step the particle in to the boundary
						goto NewRegion;
					/*
						q.x = p.x + t * p.u ;
						q.y = p.y + t * p.v ;
						q.z = p.z  + t * p.w ;
						inew=inwhere_MLC(q, index, bank);  //find where it is will be when first step inside the boundary
						p0.leafIndex = index;
						p0.bank = bank;
						*/
					}
					return inew;
				}

// now the partile is should be inside the MLC boundary
// And should have the leafIndex from previous step.
// However, just to be safe , check first.
				if(p.leafIndex > -60 ) index =p.leafIndex;  //corrected by Tong Xu, 20131208
				else 	index=indexOfLeaf(p1);

				if(index>=0)
				{
					dis_xz=rayIntersectPolygon(p1,dir1,(crossXZ+index*16),14, inside_XZ);
					//dir1 is the 2D projection vector on xz plane of the 3D direction v,
					//the second parameter of function rayIntersectPolygon() should be the nornalized vector,
					//if we use a non-normalized vector, then the result returned is length/norm of vector.
					//so in xz plane the value returned from the function is the length in xz plane/the norm of vector dir1 (in xz plane),
					//and because 3D vector v is a unit vector, the result of length in xz plane/norm of vector in xz plane comes back to length in 3D,
					//that is to say, if dir1 is a projection of a 3D unit vector, the result is the length in 3D. so is in yz plane.
					//point2D  endnode[34] ;
					//leaf_position(endnode,index);
				    //leafend is array of positions nodes in yz cross section when leaf is closed.
					//function leaf_position is to calculate the opening positions of nodes of leaf specified by the index
					//opening position at zmin+zthick/2 is saved in array cross at 15th vertex in each leaf, cross[14+i*16]
				   point2D * YZ_Polygon = endnode[index]; //default: Bank A
				   int YZ_nVertex = 17;

					switch ( p.bank)
					{
						case 0 :  YZ_Polygon =  endnode[index]+2 ; YZ_nVertex = 30; break; // in the field, between the MLC pairs 
						case 2 :  YZ_Polygon =  endnode[index]+17; break;  //bank B
					}
						dis_yz=rayIntersectPolygon(p2,dir2,YZ_Polygon,YZ_nVertex,inside_YZ) ;  
										if(dis_yz>dis_xz)  //if cross the side (xz cross section) of leaf first, then it going to the airgap.
					{
						if(dir1.x <0&&index>0) index=-index;
						else if(index<59)index= -(index+1);
					}
				}
				else if(index>-60) // the particle is in the air gaps.  60 is the number of  MLC leaves in each bank
				{
							char sx = 1;
							char sl = 16;
							if (dir1.x <0){sx =-1;sl=-16; if(indexLimit1>0) indexLimit1-=1;  indexLimit2=-index+1; }//indexLimit2 may not need to +1 ??
							else {  indexLimit2 = indexLimit1; indexLimit1 = -index-1;}
							char i= -index -sx;
							if(indexLimit1<0) indexLimit1 =0;
							if(indexLimit2>59) indexLimit2 =59;
							if(i<0) i=0;
							if(i>59) i=59;
//							bx -= crossXZ[79].x;  //this is the width of full leaf. 
							point2D *leafPolygon = (crossXZ+i*16);
							do
							{
									temp1=rayIntersectPolygon(p1,dir1,leafPolygon,14,inside_XZ);
									i+= sx;
									leafPolygon +=sl;
							} while(i>=indexLimit1&&i<=indexLimit2&&temp1>=INFI_DIST );
							if(temp1 < dis_xz) 
							{
									dis_xz = temp1;
									index = i-sx;  // the index of leaf it is going to.
//									index1 = i;
							}else{
								index=PASSMLC; //if in the air gap and does not intercept with a leaf before it go out the boundary.
								p0.module = module_order[m_VarMLCs];
								//p0.module = m_Phantom;
								//p0.module = m_Discard;
							}
	/*

							if(dir1.x<0)
							{
								i=-index+1;    //check here  later! 
								bx -= crossXZ[79].x;  //this is the width of full leaf. 
								do
								{
									temp1=rayIntersectPolygon(p1,dir1,(crossXZ+i*16),14,inside_XZ);
									i--;
								}
								while(i>=0&&temp1>=INFI_DIST && (crossXZ+i*16)->x  >= bx);
								if(temp1 < dis_xz) dis_xz = temp1;
							}
							else
							{
								i=-index-1;  //check here later! 
								do
								{
									temp2=rayIntersectPolygon(p1,dir1,(crossXZ+i*16),14,inside_XZ);
									i++;
								}
								while(i<60&&temp2>=INFI_DIST&&  (crossXZ+i*16)->x   <= bx);
								if(temp2 < dis_xz) dis_xz = temp2;
							}
					*/
					}

//finaly give the minDist the shortest intersectin distance
	  		  if (dis_xz < dis_yz) minDist = dis_xz;
		  			else minDist = dis_yz;

			 bool steppingOut = false;
				if( (minDist + EPSTEP) >= minBoundaryDist) 
				{
					minDist = minBoundaryDist;
					steppingOut=true;
				}
				if(t>minDist)
				{
						t = minDist + EPSTEP ; //step just across the interface
						if(steppingOut)
						{
							//inew = 1;  //Not VACUUM, Here we step out to air gap. ready to enter the phantom
							inew = default_reg;  //Not VACUUM, Here we step out to air gap. ready to enter the phantom
#ifdef SKIP_PHANTOM
							//inew = 0;  //set to VACUUM, this will cause the program to skip the phantom.
							inew = vacuum_reg;  //set to VACUUM, this will cause the program to skip the phantom.
#endif

#ifdef AdjustGeomPASSMLC
							// now adjust the particle coordinate and direction to account for gantry rotation and phantom isocenter.
							// however, the particle will not be able to hit the detector correctly.
									float3 tmp;
									tmp.x = p0.x;
									tmp.y=p0.y;
									tmp.z=p0.z;
									p0.x =tmp.x;  //due to a bug in the MLC field defination in HWTPS
									p0.y =tmp.z;
									p0.z =-tmp.y; //switch the z and y

									tmp.x = p0.u;
									tmp.y= p0.v;
									tmp.z= p0.w;
									p0.u= tmp.x;
									p0.v= tmp.z; //switch the v and w
									p0.w=-tmp.y;  

					//now rotate the incident particles by gantry angle

									tmp.x=   p0.x*cosAngle + p0.y*sinAngle;
									tmp.y= -p0.x*sinAngle  + p0.y*cosAngle;

									p0.x=tmp.x;
									p0.y=tmp.y;

									tmp.x= p0.u*cosAngle + p0.v*sinAngle;
									tmp.y= -p0.u*sinAngle + p0.v*cosAngle;

									p0.u=tmp.x;
									p0.v=tmp.y;

									//now shift , still need to check, should we shift before rotation? or should we subtract isocenter location?
									p0.x += isocenter_location.x;
									p0.y += isocenter_location.y;
									p0.z += isocenter_location.z;
#endif

							p0.leafIndex = PASSMLC;
							p0.module = module_order[m_VarMLCs];
							//p0.module = m_Phantom;
							return inew;
						}  //else means that particle cross some leaf but not steppingOut,
						   //this will go to NewRegion to get new inew, index and bank
				} else return inew;

NewRegion:
				q.x = p.x + t * p.u ;
				q.y = p.y + t * p.v ;
				q.z = p.z + t * p.w ;
				inew=inwhere_MLC(q, index,bank); // find out which region we will step in
				p0.leafIndex = index;   //crossing boundary cause the leafIndex to change
				p0.bank = bank;  //record its bank information

				return inew;

}


__device__ uint  howfar_SecJawsY(particle_t &p, float &t)
{
	//uint inew = 0;
	uint inew = vacuum_reg;
	point2D p0;
	//float minDist = INFI_DIST;

	point2D p1,p2,dir1,dir2;
#ifndef ROTATE_SecJaws_XY
	p1.x = p.x;  p1.y = p.z;
	p2.x = p.y;  p2.y = p.z;
	dir1.x = p.u;  dir1.y = p.w;
	dir2.x = p.v;  dir2.y = p.w;
#else
	p1.x = p.y;  p1.y = p.z;
	p2.x = p.x;  p2.y = p.z;
	dir1.x = p.v;  dir1.y = p.w;
	dir2.x = p.u;  dir2.y = p.w;
#endif

	bool inside_XZ,inside_YZ;
	float dis_xz,dis_yz;
	dis_xz=rayIntersectPolygon(p1,dir1,boundMaxY,4,inside_XZ);
	dis_yz=rayIntersectPolygon(p2,dir2,boundMaxY,4,inside_YZ);

	float minBoundaryDist = INFI_DIST;

	if(dis_xz >= INFI_DIST || dis_yz >= INFI_DIST)  //either INFI_DIST,  means the particle will not have intersection with the 3D boundary of the MLC comples,  
	{
		//inew = 0;
		return inew;
	}

	if (dis_xz < dis_yz) minBoundaryDist = dis_xz;
	else minBoundaryDist = dis_yz;

	if(inside_XZ &&              //check if the particle is inside the MLC complex in the XZ plane
		            inside_YZ )  //check if the particle is inside the MLC complex in the YZ plane
	{
		if(pointInPolygon(p2,jawsGeomY+2,4)){  //in field opening
			//inew = 3;
			inew = default_reg;
			dis_xz = rayIntersectPolygon(p2,dir2,jawsGeomY,4,inside_YZ);
			dis_yz = rayIntersectPolygon(p2,dir2,jawsGeomY+4,4,inside_YZ);
			if(dis_xz<INFI_DIST){
				if(t>dis_xz){
					//inew = 4;  //step into W of SecJaws
					inew = secjaws_reg;  //step into W of SecJaws
					t = dis_xz + EPSTEP;
				}
			}
			else if(dis_yz<INFI_DIST){
				if(t>dis_yz){
					//inew = 4;  //step into W of SecJaws
					inew = secjaws_reg;  //step into W of SecJaws
					t = dis_yz + EPSTEP;
				}
			}
			else{
				if(t>minBoundaryDist){
					//inew = 0;  //step outside from field opening
					t = minBoundaryDist + EPSTEP;
					p.module = module_order[m_SecJawY];
					//p.module = m_SecJawsX;
				}
			}
		}
		else{                                  //in W of SecJaws
			//inew = 4;
			inew = secjaws_reg;
			dis_yz = rayIntersectPolygon(p2,dir2,jawsGeomY+2,4,inside_YZ);
			if(dis_yz<INFI_DIST){
				if(t>dis_yz){
					//inew = 3;  //scatering back into field opening
					inew = default_reg;  //scatering back into field opening
					t = dis_yz + EPSTEP;
				}
			}
			else{
				if(t>minBoundaryDist){
					//inew = 0;  //step outside from W of SecJaws
					inew = vacuum_reg;  //step outside from W of SecJaws
					t = minBoundaryDist + EPSTEP;
				}
			}
		}
	}
	else                   // the particle is not inside the SecJaws complex boundaries
	{
		//inew = 3;
		inew = default_reg;
		if(t>minBoundaryDist)
		{
			t = minBoundaryDist + EPSTEP;   //step the particle into the SecJaws
#ifndef ROTATE_SecJaws_XY
			p0.x = p.x + t * p.u;
			p0.y = p.y + t * p.v;
#else
			p0.x = p.y + t * p.v;
			p0.y = p.x + t * p.u;
#endif
			if(!pointInPolygon(p0,boundMinY,4))
				//inew = 0;  //if outside the interesting area, discard!
				inew = vacuum_reg;  //if outside the interesting area, discard!
			else{
				p0.x = p0.y;
				p0.y = p.z + t * p.w;
				if(pointInPolygon(p0,jawsGeomY+2,4))
					//inew = 3;  //step into field opening from up front surface
					inew = default_reg;  //step into field opening from up front surface
				else
					//inew = 4;  //step into W from up front surface
					inew = secjaws_reg;  //step into W from up front surface
			}
		}
	}

	return inew;
}

__device__ uint  howfar_SecJawsX(particle_t &p, float &t)
{
	//uint inew = 0;
	uint inew = vacuum_reg;
	point2D p0;
	//float minDist = INFI_DIST;

	point2D p1,p2,dir1,dir2;
#ifndef ROTATE_SecJaws_XY
	p1.x = p.x;  p1.y = p.z;
	p2.x = p.y;  p2.y = p.z;
	dir1.x = p.u;  dir1.y = p.w;
	dir2.x = p.v;  dir2.y = p.w;
#else
	p1.x = p.y;  p1.y = p.z;
	p2.x = p.x;  p2.y = p.z;
	dir1.x = p.v;  dir1.y = p.w;
	dir2.x = p.u;  dir2.y = p.w;
#endif

	bool inside_XZ,inside_YZ;
	float dis_xz,dis_yz;
	dis_xz=rayIntersectPolygon(p1,dir1,boundMaxX,4,inside_XZ);
	dis_yz=rayIntersectPolygon(p2,dir2,boundMaxX,4,inside_YZ);

	float minBoundaryDist = INFI_DIST;

	if(dis_xz >= INFI_DIST || dis_yz >= INFI_DIST)  //either INFI_DIST,  means the particle will not have intersection with the 3D boundary of the MLC comples,  
	{
		//inew = 0;
		return inew;
	}

	if (dis_xz < dis_yz) minBoundaryDist = dis_xz;
	else minBoundaryDist = dis_yz;

	if(inside_XZ &&              //check if the particle is inside the MLC complex in the XZ plane
		            inside_YZ )  //check if the particle is inside the MLC complex in the YZ plane
	{
		if(pointInPolygon(p1,jawsGeomX+2,4)){  //in field opening
			//inew = 3;
			inew = default_reg;
			dis_xz = rayIntersectPolygon(p1,dir1,jawsGeomX,4,inside_XZ);
			dis_yz = rayIntersectPolygon(p1,dir1,jawsGeomX+4,4,inside_XZ);
			if(dis_xz<INFI_DIST){
				if(t>dis_xz){
					//inew = 4;  //step into W of SecJaws
					inew = secjaws_reg;  //step into W of SecJaws
					t = dis_xz + EPSTEP;
				}
			}
			else if(dis_yz<INFI_DIST){
				if(t>dis_yz){
					//inew = 4;  //step into W of SecJaws
					inew = secjaws_reg;  //step into W of SecJaws
					t = dis_yz + EPSTEP;
				}
			}
			else{
				if(t>minBoundaryDist){
					//inew = 0;  //step outside from field opening
					t = minBoundaryDist + EPSTEP;
#ifndef SKIP_MLC
					p.module = module_order[m_SecJawX];
					//p.module = m_VarMLCs;
#else
					p.module = m_Phantom;
#endif
				}
			}
		}
		else{                                  //in W of SecJaws
			//inew = 4;
			inew = secjaws_reg;
			dis_xz = rayIntersectPolygon(p1,dir1,jawsGeomX+2,4,inside_XZ);
			if(dis_xz<INFI_DIST){
				if(t>dis_xz){
					//inew = 3;  //scatering back into field opening
					inew = default_reg;  //scatering back into field opening
					t = dis_yz + EPSTEP;
				}
			}
			else{
				if(t>minBoundaryDist){
					//inew = 0;  //step outside from W of SecJaws
					inew = vacuum_reg;  //step outside from W of SecJaws
					t = minBoundaryDist + EPSTEP;
				}
			}
		}
	}
	else                   // the particle is not inside the SecJaws complex boundaries
	{
		//inew = 3;
		inew = default_reg;
		if(t>minBoundaryDist)
		{
			t = minBoundaryDist + EPSTEP;   //step the particle into the SecJaws
#ifndef ROTATE_SecJaws_XY
			p0.x = p.x + t * p.u;
			p0.y = p.y + t * p.v;
#else
			p0.x = p.y + t * p.v;
			p0.y = p.x + t * p.u;
#endif
			if(!pointInPolygon(p0,boundMinX,4))
				//inew = 0;  //if outside the interesting area, discard!
				inew = vacuum_reg;  //if outside the interesting area, discard!
			else{
				//p0.x = p0.x;
				p0.y = p.z + t * p.w;
				if(pointInPolygon(p0,jawsGeomX+2,4))
					//inew = 3;  //step into field opening from up front surface
					inew = default_reg;  //step into field opening from up front surface
				else
					//inew = 4;  //step into W from up front surface
					inew = secjaws_reg;  //step into W from up front surface
			}
		}
	}

	return inew;
}

__device__ uint  howfar_wedge(particle_t &p, float &t)
{
	uint inew = vacuum_reg;

	point2D p1,p2,dir1,dir2;
#ifndef ROTATE_Wedge_XY
	p1.x = p.x;  p1.y = p.z;
	p2.x = p.y;  p2.y = p.z;
	dir1.x = p.u;  dir1.y = p.w;
	dir2.x = p.v;  dir2.y = p.w;
#else
	p1.x = p.y;  p1.y = p.z;
	p2.x = p.x;  p2.y = p.z;
	dir1.x = p.v;  dir1.y = p.w;
	dir2.x = p.u;  dir2.y = p.w;
#endif

	bool inside_XZ,inside_YZ;
	float dis_xz,dis_yz;
	dis_xz=rayIntersectPolygon(p1,dir1,goem_Wed,nVer_Wed,inside_XZ);
	dis_yz=rayIntersectPolygon(p2,dir2,boundWed,4,inside_YZ);

	float minBoundaryDist = INFI_DIST;

	if(dis_xz >= INFI_DIST || dis_yz >= INFI_DIST)  //either INFI_DIST,  means the particle will not have intersection with the wedge comples,  
	{
		inew = default_reg;
		p.module = module_order[m_WedgeMd];
		return inew;
	}

	if (dis_xz < dis_yz) minBoundaryDist = dis_xz;
	else minBoundaryDist = dis_yz;

	if(inside_XZ &&              //check if the particle is inside the wedge complex in the XZ plane
		            inside_YZ )  //check if the particle is inside the wedge complex in the YZ plane
	{
		inew = wedge_reg;
		if(t>minBoundaryDist){
			inew = default_reg;
			t = minBoundaryDist + EPSTEP;
			p.module = module_order[m_WedgeMd];
		}
	}
	else                   // the particle is not inside the wedge complex
	{
		inew = default_reg;
		if(t>minBoundaryDist)
		{
			inew = wedge_reg;
			t = minBoundaryDist + EPSTEP;   //step the particle into the wedge
		}
	}

	return inew;
}

__device__ uint  howfar_block(particle_t &p, float &t)
{
	uint inew = vacuum_reg;

	point2D p1,p2,dir1,dir2;
#ifndef ROTATE_Block_XY
	p1.x = p.x;  p1.y = p.z;
	p2.x = p.y;  p2.y = p.z;
	dir1.x = p.u;  dir1.y = p.w;
	dir2.x = p.v;  dir2.y = p.w;
#else
	p1.x = p.y;  p1.y = p.z;
	p2.x = p.x;  p2.y = p.z;
	dir1.x = p.v;  dir1.y = p.w;
	dir2.x = p.u;  dir2.y = p.w;
#endif

	bool inside_XZ,inside_YZ;
	float dis_xz,dis_yz;
	dis_xz=rayIntersectPolygon(p1,dir1,blockXZ,4,inside_XZ);
	dis_yz=rayIntersectPolygon(p2,dir2,blockYZ,4,inside_YZ);

	float minBoundaryDist = INFI_DIST;

	if(dis_xz >= INFI_DIST || dis_yz >= INFI_DIST)  //either INFI_DIST,  means the particle will not have intersection with the block comples,  
	{
		inew = default_reg;
		p.module = module_order[m_BlockMd];
		return inew;
	}

	if (dis_xz < dis_yz) minBoundaryDist = dis_xz;
	else minBoundaryDist = dis_yz;

	if(inside_XZ &&              //check if the particle is inside the block complex in the XZ plane
		            inside_YZ )  //check if the particle is inside the block complex in the YZ plane
	{
		inew = block_reg;
		if(t>minBoundaryDist){
			inew = default_reg;
			t = minBoundaryDist + EPSTEP;
			p.module = module_order[m_BlockMd];
		}
	}
	else                   // the particle is not inside the block complex
	{
		inew = default_reg;
		if(t>minBoundaryDist)
		{
			inew = block_reg;
			t = minBoundaryDist + EPSTEP;   //step the particle into the block
		}
	}

	return inew;
}

//definition of hownear_MLC for electron simulation
__device__ float dist2(point2D v, point2D w)
{
    return (v.x - w.x)*(v.x - w.x) + (v.y - w.y)*(v.y - w.y) ;
}

__device__ float dist(point2D v, point2D w)
{
    return sqrtf((v.x - w.x)*(v.x - w.x) + (v.y - w.y)*(v.y - w.y)) ;
}

__device__ float pointDistToSegment(point2D p, point2D v, point2D w)
{
  float l2 = dist2(v, w);
  if (l2 == 0) return dist(p, v);
  float t = ((p.x - v.x) * (w.x - v.x) + (p.y - v.y) * (w.y - v.y)) / l2;
  if (t < 0) return dist(p, v);
  if (t > 1) return dist(p, w);
  point2D anchor;
  anchor.x = v.x + t * (w.x - v.x);
  anchor.y = v.y + t * (w.y - v.y);
  return dist(p, anchor);
}

__device__ float pointMinDistToPolygon(point2D p0, point2D * polygon,int nVertex) {

  int i,j;
  float minDist = 1e9;

  float dist;

  for (i=0;i<nVertex;i++)
  {
      j=i+1;
      if (j>=nVertex) j=0;
      dist = pointDistToSegment(p0,polygon[i],polygon[j]);
      if(dist<minDist)
                minDist=dist;
  }
  return minDist;
}

__device__ void hownear_MLC(particle_t &p, float &tperp)
{

				// float zmin=crossXZ[15].x;            //cross[15].x=zmin;
				point2D p1,p2;
			
#ifdef ROTATE_MLC_XY
				p1.x=p.y; p1.y=p.z;
				p2.x=p.x; p2.y=p.z;
#else
				p1.x=p.x; p1.y=p.z;
				p2.x=p.y; p2.y=p.z;
#endif
				//p1.y=p.z;
				//p1.y=p.z -(crossXZ[15].x -SourceToIsocenter); //shift the z, because the MLC geometry code is based on from surface of MLC is at z=0;
				//p2.y=p1.y;

//after the shift, the following code stays almost the same as the GPUMLC code ,expect the p0.index and p0.bank assignment  before return.

				short index ;

				float tperp_xz=INFI_DIST+1;
				float tperp_yz=INFI_DIST+1;

				float bound_xz=INFI_DIST+1;
				float bound_yz=INFI_DIST+1;
				float bound = INFI_DIST+1;

				//point2D endnode[34];
				bool inside_XZ,inside_YZ;

				if(p.leafIndex<-60)
				{
					inside_XZ=pointInPolygon(p1, boundaryXZ, 4);
					inside_YZ=pointInPolygon(p2, boundaryYZ, 4);  //bug corrected by Wu, 20131210
				}else
				{
					inside_XZ=true;
					inside_YZ=true;
				}

				bound_xz=pointMinDistToPolygon(p1,boundaryXZ,4);
				bound_yz=pointMinDistToPolygon(p2,boundaryYZ,4);
				if(bound_xz<bound_yz) bound = bound_xz;
				else bound = bound_yz;

				if(inside_XZ && inside_YZ )
				{
					if(p.leafIndex>PASSMLC) index = p.leafIndex;
					else 
						{ 
							index=indexOfLeaf(p1);
							p.leafIndex=index;
						}
					
					if(index>=0)
					{
						tperp_xz = pointMinDistToPolygon(p1,(crossXZ+index*16),14);
						//leaf_position(endnode,index);
						//leaf_position(leafend,endnode,crossXZ,index);
						if(p.bank<0) p.bank = pointInLeaf(p2,endnode[index]);

						point2D * YZ_Polygon = endnode[index]; //default: Bank A
						int YZ_nVertex = 17;
						switch ( p.bank)
						{
						    case 0 :  YZ_Polygon =  endnode[index]+2 ; YZ_nVertex = 30; break; // in the field, between the MLC pairs 
						    case 2 :  YZ_Polygon =  endnode[index]+17; break;  //bank B
						}
						tperp_yz=pointMinDistToPolygon(p2,YZ_Polygon,YZ_nVertex) ;  

					}
					else if(index>-60) // the particle is in the air gaps.  60 is the number of  MLC leaves in each bank
					{
						float temp1_xz = INFI_DIST+1;
						int i=-index;  //index=(-1,-59), i=(1,59)
						point2D *leafPolygon = (crossXZ+i*16);
						temp1_xz = pointMinDistToPolygon(p1,leafPolygon,14);
						i--;
						leafPolygon = (crossXZ+i*16);
						tperp_xz = pointMinDistToPolygon(p1,leafPolygon,14);
						if(temp1_xz<tperp_xz) tperp_xz = temp1_xz;
					}

					if(tperp_xz<tperp_yz) tperp = tperp_xz;
					else tperp = tperp_yz;
				}
				else
				{
					tperp = bound;
				}
			
				return;
}

__device__ void hownear_SecJawsY(particle_t &p, float &tperp)
{
	float bound_xz=INFI_DIST+1;
	float bound_yz=INFI_DIST+1;
	float dist_yz =INFI_DIST+1;

	point2D p1,p2;
#ifndef ROTATE_SecJaws_XY
	p1.x = p.x;  p1.y = p.z;
	p2.x = p.y;  p2.y = p.z;
#else
	p1.x = p.y;  p1.y = p.z;
	p2.x = p.x;  p2.y = p.z;
#endif

	bool inside_XZ,inside_YZ;
	inside_XZ=pointInPolygon(p1, boundMaxY, 4);
	inside_YZ=pointInPolygon(p2, boundMaxY, 4);

	bound_xz=pointMinDistToPolygon(p1,boundMaxY,4);
	bound_yz=pointMinDistToPolygon(p2,boundMaxY,4);

	if(inside_XZ && inside_YZ )
	{
		dist_yz = pointMinDistToPolygon(p2,jawsGeomY+2,4);
		if(pointInPolygon(p2,jawsGeomY+2,4)){
			if(dist_yz<bound_xz) tperp = dist_yz;
			else tperp = bound_xz;
		}
		else{
			if(dist_yz>bound_yz) dist_yz = bound_yz;
			if(dist_yz<bound_xz) tperp = dist_yz;
			else tperp = bound_xz;
		}
	}
	else
	{
		if(bound_xz<bound_yz) tperp = bound_xz;
		else tperp = bound_yz;
	}

	return;
}

__device__ void hownear_SecJawsX(particle_t &p, float &tperp)
{
	float bound_xz=INFI_DIST+1;
	float bound_yz=INFI_DIST+1;
	float dist_xz =INFI_DIST+1;

	point2D p1,p2;
#ifndef ROTATE_SecJaws_XY
	p1.x = p.x;  p1.y = p.z;
	p2.x = p.y;  p2.y = p.z;
#else
	p1.x = p.y;  p1.y = p.z;
	p2.x = p.x;  p2.y = p.z;
#endif

	bool inside_XZ,inside_YZ;
	inside_XZ=pointInPolygon(p1, boundMaxX, 4);
	inside_YZ=pointInPolygon(p2, boundMaxX, 4);

	bound_xz=pointMinDistToPolygon(p1,boundMaxX,4);
	bound_yz=pointMinDistToPolygon(p2,boundMaxX,4);

	if(inside_XZ && inside_YZ )
	{
		dist_xz = pointMinDistToPolygon(p1,jawsGeomX+2,4);
		if(pointInPolygon(p1,jawsGeomX+2,4)){
			if(dist_xz<bound_yz) tperp = dist_xz;
			else tperp = bound_yz;
		}
		else{
			if(dist_xz>bound_xz) dist_xz = bound_xz;
			if(dist_xz<bound_yz) tperp = dist_xz;
			else tperp = bound_yz;
		}
	}
	else
	{
		if(bound_xz<bound_yz) tperp = bound_xz;
		else tperp = bound_yz;
	}

	return;
}

__device__ void hownear_wedge(particle_t &p, float &tperp)
{
	point2D p1,p2;
#ifndef ROTATE_Wedge_XY
	p1.x = p.x;  p1.y = p.z;
	p2.x = p.y;  p2.y = p.z;
#else
	p1.x = p.y;  p1.y = p.z;
	p2.x = p.x;  p2.y = p.z;
#endif

	float bound_xz=pointMinDistToPolygon(p1,goem_Wed,nVer_Wed);
	float bound_yz=pointMinDistToPolygon(p2,boundWed,4);

	if(bound_xz<bound_yz) tperp = bound_xz;
	else tperp = bound_yz;

	return;
}

__device__ void hownear_block(particle_t &p, float &tperp)
{
	point2D p1,p2;
#ifndef ROTATE_Block_XY
	p1.x = p.x;  p1.y = p.z;
	p2.x = p.y;  p2.y = p.z;
#else
	p1.x = p.y;  p1.y = p.z;
	p2.x = p.x;  p2.y = p.z;
#endif

	float bound_xz=pointMinDistToPolygon(p1,blockXZ,4);
	float bound_yz=pointMinDistToPolygon(p2,blockYZ,4);

	if(bound_xz<bound_yz) tperp = bound_xz;
	else tperp = bound_yz;

	return;
}

//definition of hownear_phantom for electron simulation
__device__ float dist3D(point3D v, point3D w)
{
	return sqrtf((v.x - w.x)*(v.x - w.x) + (v.y - w.y)*(v.y - w.y) + (v.z - w.z)*(v.z - w.z)) ;
}

__device__ void hownear_phantom(particle_t &p0, float &tperp) {

	particle_t p;
#ifndef 	AdjustGeomPASSMLC   //if the particle geometry is not adjust right after PASSMLC, then adjust it every time when doing geometry in Phantom
	p.bank=p0.bank;
	p.leafIndex=p0.leafIndex;
	if(Change_xyz){
		/*
		p.x = p0.x;
		//p.y =p0.z;
		p.y =p0.z - SourceToIsocenter;
		p.z =-p0.y; //switch the z and y
		p.u=p0.u;
		p.v= p0.w; //switch the v and w
		p.w=-p0.v; 
		*/
		p.x =   p0.y;
		p.y =   (p0.z - SourceToIsocenter);
		p.z =   p0.x;
		p.u =   p0.v;
		p.v =   p0.w;
		p.w =   p0.u;
	}
	else{
		p.x = p0.x;
		p.y = p0.y;
		p.z = p0.z - SourceToIsocenter;
		p.u = p0.u;
		p.v = p0.v;
		p.w = p0.w;
	}
	p.region = p0.region;

//now rotate the incident particles by gantry angle
	float3 tmp;

	//tmp.x=   p.x*cosAngle + p.y*sinAngle;
	//tmp.y= - p.x*sinAngle + p.y*cosAngle;
	tmp.x = p.x*cosAngle - p.y*sinAngle;
	tmp.y = p.x*sinAngle + p.y*cosAngle;

	p.x=tmp.x;
	p.y=tmp.y;

	//tmp.x=   p.u*cosAngle + p.v*sinAngle;
	//tmp.y= - p.u*sinAngle + p.v*cosAngle;
	tmp.x = p.u*cosAngle - p.v*sinAngle;
	tmp.y = p.u*sinAngle + p.v*cosAngle;

	p.u=tmp.x;
	p.v=tmp.y;

	//now shift
	p.x = p.x + isocenter_location.x;
	p.y = p.y + isocenter_location.y;
	p.z = p.z + isocenter_location.z;
#else
	p=p0;  //this seams able to copy the contents of p0 to p, no need to do it individually, by Tong Xu May 2013
#endif

	tperp = INFI_DIST;  //default is a very large distance
	float tperp_x = INFI_DIST;
	float tperp_y = INFI_DIST;
	float tperp_z = INFI_DIST;

	if (p.region > 4) {    //the particle is currently in the phantom, commented by Tong Xu, Jan 2013, remember the first Voxel in the phantom has region number =2. Because the first 2 regions are used for MLC.
        // because of the above mentioned shift, we have to substract 1
		// ir is the actual voxel index
		uint ir = p.region - 1 - 2 - 2;  // " - 2" because the region index is shift by 2, used by the MLC regions.
				
        int iz = ir / (phantom_N.x * phantom_N.y); 
        ir -= iz * phantom_N.x * phantom_N.y; 
        int iy = ir / phantom_N.x;
        int ix = ir - iy * phantom_N.x;

		float temp1_x = phantom_x_bounds[ix + 1] - p.x;
		tperp_x = p.x - phantom_x_bounds[ix];
		if(temp1_x<tperp_x) tperp_x = temp1_x;

		float temp1_y = phantom_y_bounds[iy + 1] - p.y;
		tperp_y = p.y - phantom_y_bounds[iy];
		if(temp1_y<tperp_y) tperp_y = temp1_y;

		float temp1_z = phantom_z_bounds[iz + 1] - p.z;
		tperp_z = p.z - phantom_z_bounds[iz];
		if(temp1_z<tperp_z) tperp_z = temp1_z;

		if(tperp_x < tperp_y) tperp = tperp_x;
		else tperp = tperp_y;
		if(tperp > tperp_z) tperp = tperp_z;
    }
    // this part corresponds to the function howfarFromOut of the class EGS_XYZGeometry 
    // in the file egs_nd_geometry.h (v 1.26 2009/07/06) of the EGSnrc C++ Class Library
	// if p.region =0, means currently in in VACUUM,  i.e. outside of phantom. check if the particle will go into phantom
    else {
        //int zone_x=0,zone_y=0,zone_z=0;
		//point3D p0;
		point3D p1,p2;

		p1.x = phantom_x_bounds[0];
		p1.y = phantom_y_bounds[0];
		p1.z = phantom_z_bounds[0];
		p2.x = phantom_x_bounds[phantom_N.x];
		p2.y = phantom_y_bounds[phantom_N.y];
		p2.z = phantom_z_bounds[phantom_N.z];

		if(p.x-p1.x>0 && p.x-p2.x>0) p2.x = p2.x;
		else if(p.x-p1.x<0 && p.x-p2.x<0) p2.x = p1.x;
		else p2.x = p.x;
		if(p.y-p1.y>0 && p.y-p2.y>0) p2.y = p2.y;
		else if(p.y-p1.y<0 && p.y-p2.y<0) p2.y = p1.y;
		else p2.y = p.y;
		if(p.z-p1.z>0 && p.z-p2.z>0) p2.z = p2.z;
		else if(p.z-p1.z<0 && p.z-p2.z<0) p2.z = p1.z;
		else p2.z = p.z;

		p1.x = p.x;
		p1.y = p.y;
		p1.z = p.z;
		tperp = dist3D(p1,p2);

    }
}