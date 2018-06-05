//******************************************************************************
//         ****************************************
//         *                                      *
//         *    ctcreate.c                        *
//         *    rewrite from ctcreate.mortran     *
//         *                                      *
//         ****************************************
//
//
// ctcreate reads in CT data supplied by the user (Currently, ctcreate can only
// handle one CT data formats: DICOM--see description 
// of input and manual for more details) and outputs a CT phantom that can be 
// used as a direct input to dosxyznrc.  ctcreate outputs the CT phantom in ASCII
// format to the file (prefix of CT input file).egsphant.  In the .egsphant
// file, the media numbers in each voxel are output in such a way 
// as to allow the user to see rough slice-by-slice views of the original CT 
// data.
//
//******************************************************************************

/*
//#define CTUnitNumber 45  //assign a uint number for the CT data
#define CTIMAX 512
#define CTJMAX 512
#define CTKMAX 270
#define   IMAX 128  //Maximum number of x cells
#define   JMAX 128  //Maximum number of y cells
#define   KMAX 128  //Maximum number of z cells
//The above define the largest CT data set we can read in"
//You can make the code require much less space by reducing to your"
//              local maximum needs"
//IMAX, JMAX and KMAX defined in dosxyznrc_user_macros.mortran"
//              define the number of the phantom/calculational voxels"
*/

//    "ctcreate.mortran - start of subroutine SubsizeCT
//******************************************************************************
//=CT-Data-Processing-Subroutines============================================

void SubsizeCT(int *asize, short *ctdata, float *vsize, float *offset){

	//This subroutine allows the user to specify a subset of the original CT
	//volume.  The user enters the minimum and maximum X,Y and Z coordinates
	//of the planes defining the subvolume, and this subroutine determines
	//the matrix of voxel boundaries and the new size of the CT data.

		  //int           asize(3);   //The size of the array for CT data set.                              "
		  //float         vsize(3);   //The size of the voxels in the CT volume.
		  //float         offset(3);  //The offset distance of the voxels in the CT
									  // data set. This identifies the position in
									  // space of the lowest bound of the CT voxels.
		  //short         ctdata(CTIMAX,CTJMAX,CTKMAX);  //The CT data set.
		  int           i,j,k,ii,jj,kk;  //Misc indices.
		  int           imax,jmax,kmax;  //actual max i,j,k indices for CT data
		  //float         ctsubmin.x,ctsubmax.x,ctsubmin.y,ctsubmax.y,ctsubmin.z,ctsubmax.z;
		  int           ictsubmin,
						ictsubmax,
						jctsubmin,
						jctsubmax,
						kctsubmin,
						kctsubmax;
		  float xbounds[CTIMAX+1],
			    ybounds[CTJMAX+1],
			    zbounds[CTKMAX+1];

		  // Read in the dimension limits the user wants to use. If zeroes then  "
		  //  use entire volume.                                                 "
		  //OUTPUT61; ('--------------------------------------------');
		  ////printf("CT Volume subset selection.\n");
		  ////printf("Please enter the positions of limiting planes (cm):\n");
		  //OUTPUT61; ('CT Volume subset selection.');
		  //OUTPUT61; ('Please enter the positions of limiting ');
		  //OUTPUT61; ('planes (cm):');
		  //OUTPUT61; 
		  //('  ctsubmin.x,ctsubmax.x,ctsubmin.y,ctsubmax.y,ctsubmin.z,ctsubmax.z');
		  //read(5,*) ctsubmin.x,ctsubmax.x,
			//				 ctsubmin.y,ctsubmax.y,
			//				 ctsubmin.z,ctsubmax.z;
		  //OUTPUT61 ctsubmin.x,ctsubmax.x,ctsubmin.y,ctsubmax.y,ctsubmin.z,ctsubmax.z;
			//	   (' : ',6F10.4);
		  imax=asize[0];
		  jmax=asize[1];
		  kmax=asize[2];
		  //Crank out the bounds."
		  xbounds[0]=offset[0]; for(i=1;i<=imax;i++) xbounds[i]=xbounds[i-1]+vsize[0];
		  ybounds[0]=offset[1]; for(i=1;i<=jmax;i++) ybounds[i]=ybounds[i-1]+vsize[1];
		  zbounds[0]=offset[2]; for(i=1;i<=kmax;i++) zbounds[i]=zbounds[i-1]+vsize[2];
		  if (ctsubmin.x==0.0 && ctsubmax.x==0.0 && ctsubmin.y==0.0 && ctsubmax.y==0.0 &&
			  ctsubmin.z==0.0 && ctsubmax.z==0.0) {
				  ////printf(" No subset requested, will use entire CT volume.\n");
			//OUTPUT61; (' No subset requested, will use entire CT volume.');
			//OUTPUT61; ('--------------------------------------------');
				  ctsubmin.x=xbounds[0]; ctsubmax.x=xbounds[imax];
				  ctsubmin.y=ybounds[0]; ctsubmax.y=ybounds[jmax];
				  ctsubmin.z=zbounds[0]; ctsubmax.z=zbounds[kmax];
		  }
		  if ((ctsubmin.x<xbounds[0] && ctsubmax.x<xbounds[0]) ||
			  (ctsubmax.x>xbounds[imax] && ctsubmin.x>xbounds[imax])){
				  ////printf(" ***WARNING: X range does not intersect the original\n");
				  ////printf(" CT data.  Will use entire X range of original data.\n");
			//OUTPUT61; (' ***WARNING: X range does not intersect the original');
			//OUTPUT61; (' CT data.  Will use entire X range of original data.');
			//OUTPUT61; ('----------------------------------------------');
			ctsubmin.x=xbounds[0];
			ctsubmax.x=xbounds[imax];
		  } 
		  if ((ctsubmin.y<ybounds[0] && ctsubmax.y<ybounds[0]) ||
			  (ctsubmax.y>ybounds[jmax] && ctsubmin.y>ybounds[jmax])){
				  ////printf(" ***WARNING: Y range does not intersect the original\n");
				  ////printf(" CT data.  Will use entire Y range of original data.\n");
			//OUTPUT61; (' ***WARNING: Y range does not intersect the original');
			//OUTPUT61; (' CT data.  Will use entire Y range of original data.');
			//OUTPUT61; ('----------------------------------------------');
			ctsubmin.y=ybounds[0];
			ctsubmax.y=ybounds[jmax];
		  }
		  if ((ctsubmin.z<zbounds[0] && ctsubmax.z<zbounds[0]) ||
			  (ctsubmax.z>zbounds[kmax] && ctsubmin.z>zbounds[kmax])){
				  ////printf(" ***WARNING: Z range does not intersect the original\n");
				  ////printf(" CT data.  Will use entire Z range of original data.\n");
			//OUTPUT61; (' ***WARNING: Z range does not intersect the original');
			//OUTPUT61; (' CT data.  Will use entire Z range of original data.');
			//OUTPUT61; ('----------------------------------------------');
			ctsubmin.z=zbounds[0];
			ctsubmax.z=zbounds[kmax];
		  }
		  if (ctsubmax.x <= ctsubmin.x){
			  ////printf("***WARNING: X upper limit is <=  X lower limit\n");
			  ////printf(" will use entire X range.                     \n");
			//OUTPUT61; ('***WARNING: X upper limit is <=  X lower limit');
			//OUTPUT61; (' will use entire X range.                     ');
			//OUTPUT61; ('----------------------------------------------');
			ctsubmin.x=xbounds[0];
			ctsubmax.x=xbounds[imax];
		  }
		  if (ctsubmax.y <= ctsubmin.y){
			  ////printf("***WARNING: Y upper limit is <=  Y lower limit\n");
			  ////printf(" will use entire Y range.                     \n");
			//OUTPUT61; ('***WARNING: Y upper limit is <=  Y lower limit');
			//OUTPUT61; (' will use entire Y range.                     ');
			//OUTPUT61; ('----------------------------------------------');
			ctsubmin.y=ybounds[0];
			ctsubmax.y=ybounds[jmax];
		  }
		  if (ctsubmax.z <= ctsubmin.z){
			  ////printf("***WARNING: Z upper limit is <=  Z lower limit\n");
			  ////printf(" will use entire Z range.                     \n");
			//OUTPUT61; ('***WARNING: Z upper limit is <=  Z lower limit');
			//OUTPUT61; (' will use entire Z range.                     ');
			//OUTPUT61; ('----------------------------------------------');
			ctsubmin.z=zbounds[0];
			ctsubmax.z=zbounds[kmax];
		  }
		  //Bounds checking.
		  if (ctsubmin.x<xbounds[0])      {
			  ////printf("X lower limit out of bounds, will be set to \n");
			  ////printf(" lower bound.                               \n");
			//OUTPUT61; ('X lower limit out of bounds, will be set to ');
			//OUTPUT61; (' lower bound.                               ');
			//OUTPUT61; ('--------------------------------------------');
			ctsubmin.x=xbounds[0];
		  }
		  if (ctsubmax.x>xbounds[imax])   {
			  ////printf("X upper limit out of bounds, will be set to \n");
			  ////printf(" upper bound.                               \n");
			//OUTPUT61; ('X upper limit out of bounds, will be set to ');
			//OUTPUT61; (' upper bound.                               ');
			//OUTPUT61; ('--------------------------------------------');
			ctsubmax.x=xbounds[imax];
		  }
		  if (ctsubmin.y<ybounds[0])      {
			  ////printf("Y lower limit out of bounds, will be set to \n");
			  ////printf(" lower bound.                               \n");
			//OUTPUT61; ('Y lower limit out of bounds, will be set to ');
			//OUTPUT61; (' lower bound.                               ');
			//OUTPUT61; ('--------------------------------------------');
			ctsubmin.y=ybounds[0];
		  }
		  if (ctsubmax.y>ybounds[jmax])   {
			  ////printf("Y upper limit out of bounds, will be set to \n");
			  ////printf(" upper bound.                               \n");
			//OUTPUT61; ('Y upper limit out of bounds, will be set to ');
			//OUTPUT61; (' upper bound.                               ');
			//OUTPUT61; ('--------------------------------------------');
			ctsubmax.y=ybounds[jmax];
		  }
		  if (ctsubmin.z<zbounds[0])      {
			  ////printf("Z lower limit out of bounds, will be set to \n");
			  ////printf(" lower bound.                               \n");
			//OUTPUT61; ('Z lower limit out of bounds, will be set to ');
			//OUTPUT61; (' lower bound.                               ');
			//OUTPUT61; ('--------------------------------------------');
			ctsubmin.z=zbounds[0];
		  }
		  if (ctsubmax.z>zbounds[kmax]) {
			  ////printf("Z upper limit out of bounds, will be set to \n");
			  ////printf(" upper bound.                               \n");
			//OUTPUT61; ('Z upper limit out of bounds, will be set to ');
			//OUTPUT61; (' upper bound.                               ');
			//OUTPUT61; ('--------------------------------------------');
			ctsubmax.z=zbounds[kmax];
		  }
		  //Do the shift of the arrays
		  //Determine the indices that will include the desired dimensions.
		  ////printf("The voxel index limits are as follows:\n");
		  //OUTPUT61;('The voxel index limits are as follows:');
		  for(i=0;i<imax;i++) {
			  if(xbounds[i]<=ctsubmin.x && xbounds[i+1]>ctsubmin.x)  ictsubmin=i;
			  if(xbounds[i]<=ctsubmax.x && xbounds[i+1]>=ctsubmax.x) ictsubmax=i;
		  }
		  ////printf("I Limits -> i=%6d to i=%6d\n",ictsubmin,ictsubmax);
		  //OUTPUT61 ictsubmin,ictsubmax;('I Limits -> i=',I6,' to i=',I6);
		  for(j=0;j<jmax;j++) {
			  if(ybounds[j]<=ctsubmin.y && ybounds[j+1]>ctsubmin.y)  jctsubmin=j;
			  if(ybounds[j]<=ctsubmax.y && ybounds[j+1]>=ctsubmax.y) jctsubmax=j;
		  }
		  ////printf("J Limits -> j=%6d to j=%6d\n",jctsubmin,jctsubmax);
		  //OUTPUT61 jctsubmin,jctsubmax;('J Limits -> j=',I6,' to j=',I6);
		  for(k=0;k<kmax;k++) {
			  if(zbounds[k]<=ctsubmin.z && zbounds[k+1]>ctsubmin.z)  kctsubmin=k;
			  if(zbounds[k]<=ctsubmax.z && zbounds[k+1]>=ctsubmax.z) kctsubmax=k;
		  }
		  if (zbounds[kmax-1]==ctsubmax.z ) kctsubmax=kmax-1;  //check here later!!!!!!!!!!!!!!!
		  ////printf("K Limits -> k=%6d to k=%6d\n",kctsubmin,kctsubmax);
		  //OUTPUT61 kctsubmin,kctsubmax;('K Limits -> k=',I6,' to k=',I6);
		  ctsubmin.x=xbounds[ictsubmin]; ctsubmax.x=xbounds[ictsubmax+1]; 
		  ctsubmin.y=ybounds[jctsubmin]; ctsubmax.y=ybounds[jctsubmax+1]; 
		  ctsubmin.z=zbounds[kctsubmin]; ctsubmax.z=zbounds[kctsubmax+1]; 
		  ////printf(" xctsubmin,xctsubmax,yctsubmin,yctsubmax,zctsubmin,zctsubmax (cm)\n");
		  ////printf(" after adjustment to fit integer no. of voxels \n");
		  ////printf(" :%10.4f,%10.4f,%10.4f,%10.4f,%10.4f,%10.4f\n",ctsubmin.x,ctsubmax.x,ctsubmin.y,ctsubmax.y,ctsubmin.z,ctsubmax.z);
		  //OUTPUT61 ctsubmin.x,ctsubmax.x,ctsubmin.y,ctsubmax.y,ctsubmin.z,ctsubmax.z;
			//  (' ctsubmin.x,ctsubmax.x,ctsubmin.y,ctsubmax.y,ctsubmin.z,ctsubmax.z (cm)'/
			//   ' after adjustment to fit integer no. of voxels '/
			//   ' : ',6F10.4);
		  //Copy(reset?) new CT volume into a temporary array.
		  ii=0;
		  jj=0;
		  kk=0;
		  for(k=kctsubmin;k<=kctsubmax;k++) {
			  for(j=jctsubmin;j<=jctsubmax;j++) {
				for(i=ictsubmin;i<=ictsubmax;i++) {
				  ctdata[ii+jj*CTIMAX+kk*CTJMAX*CTIMAX]=ctdata[i+j*CTIMAX+k*CTJMAX*CTIMAX];
				   /*  This check is now done in the ReadCT routine
				   for the format of the data being read in and/or
				   the subroutine CTToMedium 
				   if (ctdata(ii,jj,kk)<0 |
					  ii<1 | ii>CTIMAX |
					  jj<1 | jj>CTJMAX |
					  kk<1 | kk>CTKMAX 
					  ) [
					OUTPUT61 ii,jj,kk,i,j,k,ctdata(ii,jj,kk);
							 ('CT Number out of bounds at',6I5,I10);
				  ];
				  */
				  ii=ii+1;
				}
				ii=0;
				jj=jj+1;
			  }
			  jj=0;
			  kk=kk+1;
		  }
			//Reset offsets
		  offset[0]=offset[0]+vsize[0]*ictsubmin;
		  offset[1]=offset[1]+vsize[1]*jctsubmin;
		  offset[2]=offset[2]+vsize[2]*kctsubmin;
		  asize[0]=ictsubmax-ictsubmin+1;
		  asize[1]=jctsubmax-jctsubmin+1;
		  asize[2]=kctsubmax-kctsubmin+1;
            
		  ////printf("============================================\n");
		  //OUTPUT61; ('============================================');
}

//-----Subroutine-ReSampleCT-------------------------------------------------"
void ResampleCT(int ct_imax, int ct_jmax, int ct_kmax,
				float ct_xthickness, float ct_ythickness, float ct_zthickness,
				short *ct_data,
				int &xyz_imax, int &xyz_jmax, int &xyz_kmax,
				float *xyz_xbounds, float *xyz_ybounds, float *xyz_zbounds,
				short *xyz_ct,
				float *CTOffset){
	//---------------------------------------------------------------------------"
	//This subroutine should allow up- and down-sizing of CT data
	//limited by IMAX,JMAX,KMAX--the max number of voxels allowed in
	//the dosxyznrc phantom
	//
	//Finding the xyz boundaries for each ct voxel, determining the weight to 
	//apply to the ct voxel, and calculating the resized xyz data by multiplying 
	//the ct data by its appropriate weights is done in a single set of nested loops.
	//This allows the weight array to be 1-D and
	//also allows the upper and lower bounds in the x,y and z directions
	//to be single variables
	//
	//Also, the user inputs the dimensions of the dosxyznrc voxels.
	//This means that this subroutine actually calculates imax, jmax, kmax
	//and xbound(i),ybound(i),zbound(i) to be used by dosxyznrc from now on.  
	//Of course, imax, jmax, kmax are limited by IMAX, JMAX, KMAX.
	//This method assumes that the dimensions of either the total CT
	//data or of each CT voxel (or both) are known and can be input to the 
	//subroutine. 
	//---------------------------------------------------------------------------

	    //int ct_imax;     //the max # of ct voxels in x-direction"
		//int ct_jmax;     //the max # of ct voxels in y-direction"
		//int ct_kmax;     //the max # of ct voxels in z-direction"
		//int xyz_imax;    //max # of xyz voxels in x-direction"
		//int xyz_jmax;    //max # of xyz voxels in y-direction"
		//int xyz_kmax;    //max # of xyz voxels in z-direction"
		int	i_lower_xyz; //index for lower x bounds of ct voxels
		int	i_upper_xyz; //index for upper x bounds of ct voxels
		int	j_lower_xyz; //index for lower y bounds of ct voxels
		int	j_upper_xyz; //index for upper y bounds of ct voxels
		int	k_lower_xyz; //index for lower z bounds of ct voxels
		int	k_upper_xyz; //index for upper z bounds of ct voxels
		int	param;       //parameter for the 1-D weight array
		int	i_ct,j_ct,k_ct,i_xyz,j_xyz,k_xyz,i,j,k; //loop indices

	//short   ct_data[CTIMAX*CTJMAX*CTKMAX]; //Original ct data."
	//short   xyz_ct[IMAX*JMAX*KMAX];        //Resampled ct data."

	//float   xyz_xbounds(IMAX+1); //xyz voxel boundaries in x direction"
	//float   xyz_ybounds(JMAX+1); //xyz voxel boundaries in y direction"
	//float   xyz_zbounds(KMAX+1); //xyz voxel boundaries in z direction"

	//float   ct_xthickness;  //thickness of each ct voxel in x-direction"
	//float   ct_ythickness;  //thickness of each ct voxel in y-direction"   
	//float   ct_zthickness;  //thickness of each ct voxel in z-direction"
	                          //Current incarnation assumes that these are"
	                          //available and can be passed to the subroutine."
	                          //If not, they can be calculated from the total"
	                          //dimensions of the CT data"

	//float	xyz_thickness.x;   //thickness of each xyz voxel in x-direction"
	//float xyz_thickness.y;   //thickness of each xyz voxel in y-direction" 
	//float	xyz_thickness.z;   //thickness of each xyz voxel in z-direction"
	float	weight_xyz[CTIMAX+CTJMAX+CTKMAX]; //fraction of xyz"
	                          //voxel in ct voxel as calculated separately in x,y,"
	                          //and z-directions"
	//float   CTOffset(3);    //The offset of the CT data from the origin."
	float	*realxyz_ct = (float*)malloc(XYZREG*sizeof(float));
	//float	realxyz_ct[IMAX*JMAX*KMAX];

	//char iorjork; //stores character i or j or k for use in FIND_WEIGHT macro" 
	//---------------------------------------------------------------------------"
	/*
	REPLACE {$FIND_UPPER_LOWER_#_#_BNDS;} WITH {
	;
	DO i_xyz=1,xyz_{P1}max[
	   "find the xyz_voxel that the ct voxel lower bound lies in"
	   if((xyz_{P2}bounds(i_xyz) <=
		  xyz_{P2}bounds(1)+ct_{P2}thickness*({P1}_ct-1)) &
		  (xyz_{P2}bounds(i_xyz+1) >=
		   xyz_{P2}bounds(1)+ct_{P2}thickness*({P1}_ct-1)))[
			  {P1}_lower_xyz=i_xyz;
	   ]
	   "find the xyz_voxel that the ct voxel upper bound lies in"
	   if((xyz_{P2}bounds(i_xyz) <= 
		   xyz_{P2}bounds(1)+ct_{P2}thickness*FLOAT({P1}_ct)) &
		  (xyz_{P2}bounds(i_xyz+1) >= 
		   xyz_{P2}bounds(1)+ct_{P2}thickness*FLOAT({P1}_ct))
		  )[
			 {P1}_upper_xyz=i_xyz;
			 EXIT;
	   ]
	]
	;
	}
	;
	*/
	/*
	REPLACE {$FIND_WEIGHT_#_#;} WITH {
	;
	"param below allows the use of a 1-D weight array"
	iorjork='{P1}'; 
	if(iorjork='i')[ param=0; ]
	ELSEIF(iorjork='j')[ param=xyz_imax; ]
	ELSEIF(iorjork='k')[ param=xyz_imax+xyz_jmax; ]

	if({P1}_lower_xyz={P1}_upper_xyz)["lower and upper bounds are in"
											 "the same xyz voxel"
	   weight_xyz(param+{P1}_lower_xyz)=ct_{P2}thickness/
										(xyz_{P2}bounds({P1}_lower_xyz+1)-
										 xyz_{P2}bounds({P1}_lower_xyz));
	]
	ELSE["lower and upper bounds NOT in the same xyz voxel"
	  DO i_xyz={P1}_lower_xyz,{P1}_upper_xyz[
		if((xyz_{P2}bounds(i_xyz)>=
		   (xyz_{P2}bounds(1)+ct_{P2}thickness*({P1}_ct-1)))&
		   (xyz_{P2}bounds(i_xyz+1)<=
		   (xyz_{P2}bounds(1)+ct_{P2}thickness*({P1}_ct))))[
		   "the xyz voxel is entirely in the ct voxel"
			  weight_xyz(param+i_xyz)=1.00;
		]
		ELSEIF((xyz_{P2}bounds(i_xyz)<=
			   (xyz_{P2}bounds(1)+ct_{P2}thickness*({P1}_ct-1)))&
			   (xyz_{P2}bounds(i_xyz+1)<=
			   (xyz_{P2}bounds(1)+ct_{P2}thickness*({P1}_ct))))[
			  "the ct voxel straddles the upper bound of the xyz voxel"
			  weight_xyz(param+i_xyz)=
							(xyz_{P2}bounds(i_xyz+1)-(xyz_{P2}bounds(1)+
							 ct_{P2}thickness*({P1}_ct-1)))/
							 (xyz_{P2}bounds(i_xyz+1)-xyz_{P2}bounds(i_xyz));
		]
		ELSE[ "the ct voxel straddles the lower bound of the xyz voxel"
			  weight_xyz(param+i_xyz)=
							((xyz_{P2}bounds(1)+ct_{P2}thickness*({P1}_ct))-
									  xyz_{P2}bounds(i_xyz))/
							 (xyz_{P2}bounds(i_xyz+1)-xyz_{P2}bounds(i_xyz));
		]
	  ]
	]
	;
	}
	;
	*/


	/*
	for(i_xyz=0;i_xyz<xyz_{P1}max;i_xyz++){
	   //find the xyz_voxel that the ct voxel lower bound lies in"
	   if((xyz_{P2}bounds[i_xyz] <=
		  xyz_{P2}bounds[0]+ct_{P2}thickness*{P1}_ct) &&
		  (xyz_{P2}bounds[i_xyz+1] >=
		   xyz_{P2}bounds[0]+ct_{P2}thickness*{P1}_ct)){
			  {P1}_lower_xyz=i_xyz;
	   }
	   //find the xyz_voxel that the ct voxel upper bound lies in"
	   if((xyz_{P2}bounds[i_xyz] <= 
		   xyz_{P2}bounds[0]+ct_{P2}thickness*float({P1}_ct+1)) &&
		  (xyz_{P2}bounds[i_xyz+1] >= 
		   xyz_{P2}bounds[0]+ct_{P2}thickness*float({P1}_ct+1))
		  ){
			 {P1}_upper_xyz=i_xyz;
			 break;
	   }
	}

	//param below allows the use of a 1-D weight array"
	iorjork='{P1}'; 
	if(iorjork=='i')      param=0;
	else if(iorjork=='j') param=xyz_imax;
	else if(iorjork=='k') param=xyz_imax+xyz_jmax;

	if({P1}_lower_xyz=={P1}_upper_xyz){//lower and upper bounds are in the same xyz voxel"
	   weight_xyz[param+{P1}_lower_xyz]=ct_{P2}thickness/
										(xyz_{P2}bounds[{P1}_lower_xyz+1]-
										 xyz_{P2}bounds[{P1}_lower_xyz]);
	}
	else{//lower and upper bounds NOT in the same xyz voxel"
	  for(i_xyz={P1}_lower_xyz;i_xyz<={P1}_upper_xyz;i_xyz++){
		if((xyz_{P2}bounds[i_xyz]>=
		   (xyz_{P2}bounds[0]+ct_{P2}thickness*{P1}_ct))&&
		   (xyz_{P2}bounds[i_xyz+1]<=
		   (xyz_{P2}bounds[0]+ct_{P2}thickness*({P1}_ct+1)))){
		   //the xyz voxel is entirely in the ct voxel"
			  weight_xyz[param+i_xyz]=1.00;
		}
		else if((xyz_{P2}bounds[i_xyz]<=
			   (xyz_{P2}bounds[0]+ct_{P2}thickness*{P1}_ct))&&
			   (xyz_{P2}bounds[i_xyz+1]<=
			   (xyz_{P2}bounds[0]+ct_{P2}thickness*({P1}_ct+1)))){
			  //the ct voxel straddles the upper bound of the xyz voxel"
			  weight_xyz[param+i_xyz]=
							(xyz_{P2}bounds[i_xyz+1]-(xyz_{P2}bounds[0]+
							 ct_{P2}thickness*{P1}_ct))/
							 (xyz_{P2}bounds[i_xyz+1]-xyz_{P2}bounds[i_xyz]);
		}
		else{ //the ct voxel straddles the lower bound of the xyz voxel"
			  weight_xyz[param+i_xyz]=
							((xyz_{P2}bounds[0]+ct_{P2}thickness*({P1}_ct+1))-
									  xyz_{P2}bounds[i_xyz])/
							 (xyz_{P2}bounds[i_xyz+1]-xyz_{P2}bounds[i_xyz]);
		}
	  }
	}
	*/


	//---------------------------------------------------------------------------"
	////printf(" Resample CT data for dosxyznrc \n");
	//OUTPUT61;(/' Resample CT data for dosxyznrc '/
	//		   ' --------------------------- ');

	//Have the user enter the x,y,z dimensions that he/she wants.  Previously"
	//these dimensions were passed as arguments to this subroutine"

	//:INPUT_DIMENSIONS:
	////printf(" Input the x,y,z dimensions (cm) of the dosxyznrc voxels on one line\n");
	////printf(" (min= %12.5fx%12.5fx%12.5f cm):\n",ct_imax*ct_xthickness/IMAX,ct_jmax*ct_ythickness/JMAX,ct_kmax*ct_zthickness/KMAX);
	//OUTPUT61 ct_imax*ct_xthickness/IMAX,ct_jmax*ct_ythickness/JMAX,
	//	   ct_kmax*ct_zthickness/KMAX; 
	//	 (/' Input the x,y,z dimensions (cm) of the dosxyznrc voxels on one line'/
	//		' (min= ',F12.5,' x',F12.5,' x',F12.5,' cm)'/
	//		' :',$);
	//read(5,*) xyz_thickness.x,xyz_thickness.y,xyz_thickness.z;
	//read in from .ini file
	////printf("%12.5f, %12.5f, %12.5f\n",xyz_thickness.x,xyz_thickness.y,xyz_thickness.z);
	//OUTPUT61 xyz_thickness.x,xyz_thickness.y,xyz_thickness.z; (3F12.5);

	if(xyz_thickness.x <=0.0 || xyz_thickness.y <=0.0 || xyz_thickness.z <=0.0){
		error(5002," Dimensions must all be positive.  Try again.\n");
	   //OUTPUT61; (' Dimensions must all be positive.  Try again.');
	   //STOP; 
	}
	else if(xyz_thickness.x > ct_imax*ct_xthickness ||
		    xyz_thickness.y > ct_jmax*ct_ythickness ||
		    xyz_thickness.z > ct_kmax*ct_zthickness){
				error(5003," Dimension in a direction cannot be greater than total size of\n CT data in that direction.  Try again.\n");
	   //OUTPUT61;(' Dimension in a direction cannot be greater than total size of'/
		//	   ' CT data in that direction.  Try again.');
	   //STOP; 
	}
	else if(xyz_thickness.x < ct_imax*ct_xthickness/IMAX ||
		    xyz_thickness.y < ct_jmax*ct_ythickness/JMAX ||
		    xyz_thickness.z < ct_kmax*ct_zthickness/KMAX){
				error(5004," Dimensions in at least one direction < min allowed. Either increase\n dimension(s) or increase IMAX, JMAX and/or KMAX.\n");
	   //OUTPUT61;
	   //(' Dimensions in at least one direction < min allowed. Either increase'/
		//' dimension(s) or go into dosxyznrc_user_macros.mortran and increase IMAX,'/
		//' JMAX and/or KMAX');
	   //STOP; 
	}

	xyz_imax=int(ct_imax*ct_xthickness/xyz_thickness.x);
	xyz_jmax=int(ct_jmax*ct_ythickness/xyz_thickness.y);
	xyz_kmax=int(ct_kmax*ct_zthickness/xyz_thickness.z);

	//adjust x,y,z dimensions so that the voxels fit exactly on CT data"

	xyz_thickness.x=float(ct_imax)*ct_xthickness/xyz_imax;
	xyz_thickness.y=float(ct_jmax)*ct_ythickness/xyz_jmax;
	xyz_thickness.z=float(ct_kmax)*ct_zthickness/xyz_kmax;
	//printf("New X voxel thickness -> %10.2f\n",xyz_thickness.x);
	//printf("New Y voxel thickness -> %10.2f\n",xyz_thickness.y);
	//printf("New Z voxel thickness -> %10.2f\n",xyz_thickness.z);
	//printf("New number X voxels   -> %10d\n",xyz_imax);
	//printf("New number Y voxels   -> %10d\n",xyz_jmax);
	//printf("New number Z voxels   -> %10d\n",xyz_kmax);
	//OUTPUT61 xyz_thickness.x; ('New X voxel thickness -> ',F10.2);
	//OUTPUT61 xyz_thickness.y; ('New Y voxel thickness -> ',F10.2);
	//OUTPUT61 xyz_thickness.z; ('New Z voxel thickness -> ',F10.2);
	//OUTPUT61 xyz_imax;        ('New number X voxels   -> ',I10);
	//OUTPUT61 xyz_jmax;        ('New number Y voxels   -> ',I10);
	//OUTPUT61 xyz_kmax;        ('New number Z voxels   -> ',I10);

	//printf(" Final x,y,z dimensions of dosxyznrc voxels in cm (adjusted so that an\n");
	//printf(" integer number fit exactly on the CT data):\n");
	//printf(" %12.5f, %12.5f, %12.5f\n",xyz_thickness.x,xyz_thickness.y,xyz_thickness.z);
	//OUTPUT61 xyz_thickness.x,xyz_thickness.y,xyz_thickness.z; 
	//(' Final x,y,z dimensions of dosxyznrc voxels in cm (adjusted so that an', 
	// ' integer'/
	// ' number fit exactly on the CT data):',3F12.5);

	//calculate the bounds of the xyz voxels"
	xyz_xbounds[0]=CTOffset[0];
	xyz_ybounds[0]=CTOffset[1];
	xyz_zbounds[0]=CTOffset[2];
	for(i_xyz=1;i_xyz<=xyz_imax;i_xyz++){
	   xyz_xbounds[i_xyz]=xyz_xbounds[i_xyz-1]+xyz_thickness.x;
	}
	for(i_xyz=1;i_xyz<=xyz_jmax;i_xyz++){
	   xyz_ybounds[i_xyz]=xyz_ybounds[i_xyz-1]+xyz_thickness.y;
	}
	for(i_xyz=1;i_xyz<=xyz_kmax;i_xyz++){
	   xyz_zbounds[i_xyz]=xyz_zbounds[i_xyz-1]+xyz_thickness.z;
	}

	//zero out the weights and the resized ct data
	for(i_xyz=0;i_xyz<xyz_imax+xyz_jmax+xyz_kmax;i_xyz++){
	   weight_xyz[i_xyz]=0.0;
	}

	//zero out the resized ct data
	for(i_xyz=0;i_xyz<xyz_imax;i_xyz++){
	  for(j_xyz=0;j_xyz<xyz_jmax;j_xyz++){
		for(k_xyz=0;k_xyz<xyz_kmax;k_xyz++){
		  realxyz_ct[i_xyz+j_xyz*IMAX+k_xyz*JMAX*IMAX]=0.0;
		}
	  }
	}

	//printf("Calculating bounds and new CT values\n");
	//OUTPUT61; ('Calculating bounds and new CT values');
	for(i_ct=0;i_ct<ct_imax;i_ct++){
		//$FIND_UPPER_LOWER_i_x_BNDS;
		//$FIND_WEIGHT_i_x;
		for(i_xyz=0;i_xyz<xyz_imax;i_xyz++){
		   //find the xyz_voxel that the ct voxel lower bound lies in"
		   if((xyz_xbounds[i_xyz] <=
			  xyz_xbounds[0]+ct_xthickness*i_ct) &&
			  (xyz_xbounds[i_xyz+1] >=
			   xyz_xbounds[0]+ct_xthickness*i_ct)){
				  i_lower_xyz=i_xyz;
		   }
		   //find the xyz_voxel that the ct voxel upper bound lies in"
		   if((xyz_xbounds[i_xyz] <= 
			   xyz_xbounds[0]+ct_xthickness*float(i_ct+1)) &&
			  (xyz_xbounds[i_xyz+1] >= 
			   xyz_xbounds[0]+ct_xthickness*float(i_ct+1))
			  ){
				 i_upper_xyz=i_xyz;
				 break;
		   }
		}

		//param below allows the use of a 1-D weight array"
		//iorjork='i'; 
		//if(iorjork=='i')      param=0;
		//else if(iorjork=='j') param=xyz_imax;
		//else if(iorjork=='k') param=xyz_imax+xyz_jmax;
		param=0;

		if(i_lower_xyz==i_upper_xyz){//lower and upper bounds are in the same xyz voxel"
		   weight_xyz[param+i_lower_xyz]=ct_xthickness/
											(xyz_xbounds[i_lower_xyz+1]-
											 xyz_xbounds[i_lower_xyz]);
		}
		else{//lower and upper bounds NOT in the same xyz voxel"
		  for(i_xyz=i_lower_xyz;i_xyz<=i_upper_xyz;i_xyz++){
			if((xyz_xbounds[i_xyz]>=
			   (xyz_xbounds[0]+ct_xthickness*i_ct))&&
			   (xyz_xbounds[i_xyz+1]<=
			   (xyz_xbounds[0]+ct_xthickness*(i_ct+1)))){
			   //the xyz voxel is entirely in the ct voxel"
				  weight_xyz[param+i_xyz]=1.00;
			}
			else if((xyz_xbounds[i_xyz]<=
				   (xyz_xbounds[0]+ct_xthickness*i_ct))&&
				   (xyz_xbounds[i_xyz+1]<=
				   (xyz_xbounds[0]+ct_xthickness*(i_ct+1)))){
				  //the ct voxel straddles the upper bound of the xyz voxel"
				  weight_xyz[param+i_xyz]=
								(xyz_xbounds[i_xyz+1]-(xyz_xbounds[0]+
								 ct_xthickness*i_ct))/
								 (xyz_xbounds[i_xyz+1]-xyz_xbounds[i_xyz]);
			}
			else{ //the ct voxel straddles the lower bound of the xyz voxel"
				  weight_xyz[param+i_xyz]=
								((xyz_xbounds[0]+ct_xthickness*(i_ct+1))-
										  xyz_xbounds[i_xyz])/
								 (xyz_xbounds[i_xyz+1]-xyz_xbounds[i_xyz]);
			}
		  }
		}

		for(j_ct=0;j_ct<ct_jmax;j_ct++){
			//$FIND_UPPER_LOWER_j_y_BNDS;
			//$FIND_WEIGHT_j_y;
			for(i_xyz=0;i_xyz<xyz_jmax;i_xyz++){
			   //find the xyz_voxel that the ct voxel lower bound lies in"
			   if((xyz_ybounds[i_xyz] <=
				  xyz_ybounds[0]+ct_ythickness*j_ct) &&
				  (xyz_ybounds[i_xyz+1] >=
				   xyz_ybounds[0]+ct_ythickness*j_ct)){
					  j_lower_xyz=i_xyz;
			   }
			   //find the xyz_voxel that the ct voxel upper bound lies in"
			   if((xyz_ybounds[i_xyz] <= 
				   xyz_ybounds[0]+ct_ythickness*float(j_ct+1)) &&
				  (xyz_ybounds[i_xyz+1] >= 
				   xyz_ybounds[0]+ct_ythickness*float(j_ct+1))
				  ){
					 j_upper_xyz=i_xyz;
					 break;
			   }
			}

			//param below allows the use of a 1-D weight array"
			//iorjork='j'; 
			//if(iorjork=='i')      param=0;
			//else if(iorjork=='j') param=xyz_imax;
			//else if(iorjork=='k') param=xyz_imax+xyz_jmax;
			param=xyz_imax;

			if(j_lower_xyz==j_upper_xyz){//lower and upper bounds are in the same xyz voxel"
			   weight_xyz[param+j_lower_xyz]=ct_ythickness/
												(xyz_ybounds[j_lower_xyz+1]-
												 xyz_ybounds[j_lower_xyz]);
			}
			else{//lower and upper bounds NOT in the same xyz voxel"
			  for(i_xyz=j_lower_xyz;i_xyz<=j_upper_xyz;i_xyz++){
				if((xyz_ybounds[i_xyz]>=
				   (xyz_ybounds[0]+ct_ythickness*j_ct))&&
				   (xyz_ybounds[i_xyz+1]<=
				   (xyz_ybounds[0]+ct_ythickness*(j_ct+1)))){
				   //the xyz voxel is entirely in the ct voxel"
					  weight_xyz[param+i_xyz]=1.00;
				}
				else if((xyz_ybounds[i_xyz]<=
					   (xyz_ybounds[0]+ct_ythickness*j_ct))&&
					   (xyz_ybounds[i_xyz+1]<=
					   (xyz_ybounds[0]+ct_ythickness*(j_ct+1)))){
					  //the ct voxel straddles the upper bound of the xyz voxel"
					  weight_xyz[param+i_xyz]=
									(xyz_ybounds[i_xyz+1]-(xyz_ybounds[0]+
									 ct_ythickness*j_ct))/
									 (xyz_ybounds[i_xyz+1]-xyz_ybounds[i_xyz]);
				}
				else{ //the ct voxel straddles the lower bound of the xyz voxel"
					  weight_xyz[param+i_xyz]=
									((xyz_ybounds[0]+ct_ythickness*(j_ct+1))-
											  xyz_ybounds[i_xyz])/
									 (xyz_ybounds[i_xyz+1]-xyz_ybounds[i_xyz]);
				}
			  }
			}

			for(k_ct=0;k_ct<ct_kmax;k_ct++){
				//$FIND_UPPER_LOWER_k_z_BNDS;
				//$FIND_WEIGHT_k_z;
				for(i_xyz=0;i_xyz<xyz_kmax;i_xyz++){
				   //find the xyz_voxel that the ct voxel lower bound lies in"
				   if((xyz_zbounds[i_xyz] <=
					  xyz_zbounds[0]+ct_zthickness*k_ct) &&
					  (xyz_zbounds[i_xyz+1] >=
					   xyz_zbounds[0]+ct_zthickness*k_ct)){
						  k_lower_xyz=i_xyz;
				   }
				   //find the xyz_voxel that the ct voxel upper bound lies in"
				   if((xyz_zbounds[i_xyz] <= 
					   xyz_zbounds[0]+ct_zthickness*float(k_ct+1)) &&
					  (xyz_zbounds[i_xyz+1] >= 
					   xyz_zbounds[0]+ct_zthickness*float(k_ct+1))
					  ){
						 k_upper_xyz=i_xyz;
						 break;
				   }
				}

				//param below allows the use of a 1-D weight array"
				//iorjork='k'; 
				//if(iorjork=='i')      param=0;
				//else if(iorjork=='j') param=xyz_imax;
				//else if(iorjork=='k') param=xyz_imax+xyz_jmax;
				param=xyz_imax+xyz_jmax;

				if(k_lower_xyz==k_upper_xyz){//lower and upper bounds are in the same xyz voxel"
				   weight_xyz[param+k_lower_xyz]=ct_zthickness/
													(xyz_zbounds[k_lower_xyz+1]-
													 xyz_zbounds[k_lower_xyz]);
				}
				else{//lower and upper bounds NOT in the same xyz voxel"
				  for(i_xyz=k_lower_xyz;i_xyz<=k_upper_xyz;i_xyz++){
					if((xyz_zbounds[i_xyz]>=
					   (xyz_zbounds[0]+ct_zthickness*k_ct))&&
					   (xyz_zbounds[i_xyz+1]<=
					   (xyz_zbounds[0]+ct_zthickness*(k_ct+1)))){
					   //the xyz voxel is entirely in the ct voxel"
						  weight_xyz[param+i_xyz]=1.00;
					}
					else if((xyz_zbounds[i_xyz]<=
						   (xyz_zbounds[0]+ct_zthickness*k_ct))&&
						   (xyz_zbounds[i_xyz+1]<=
						   (xyz_zbounds[0]+ct_zthickness*(k_ct+1)))){
						  //the ct voxel straddles the upper bound of the xyz voxel"
						  weight_xyz[param+i_xyz]=
										(xyz_zbounds[i_xyz+1]-(xyz_zbounds[0]+
										 ct_zthickness*k_ct))/
										 (xyz_zbounds[i_xyz+1]-xyz_zbounds[i_xyz]);
					}
					else{ //the ct voxel straddles the lower bound of the xyz voxel"
						  weight_xyz[param+i_xyz]=
										((xyz_zbounds[0]+ct_zthickness*(k_ct+1))-
												  xyz_zbounds[i_xyz])/
										 (xyz_zbounds[i_xyz+1]-xyz_zbounds[i_xyz]);
					}
				  }
				}

				for(i_xyz=i_lower_xyz;i_xyz<=i_upper_xyz;i_xyz++){
					for(j_xyz=j_lower_xyz;j_xyz<=j_upper_xyz;j_xyz++){
						for(k_xyz=k_lower_xyz;k_xyz<=k_upper_xyz;k_xyz++){
						   realxyz_ct[i_xyz+j_xyz*IMAX+k_xyz*JMAX*IMAX]+=
												ct_data[i_ct+j_ct*CTIMAX+k_ct*CTJMAX*CTIMAX]* 
												weight_xyz[i_xyz]*
												weight_xyz[xyz_imax+j_xyz]*
												weight_xyz[xyz_imax+xyz_jmax+k_xyz];
							//                   OUTPUT61 i_xyz,j_xyz,k_xyz,realxyz_ct(i_xyz,j_xyz,k_xyz),
							//                            ct_data(i_ct,j_ct,k_ct),
							//                            weight_xyz(i_xyz);
							//                   ('In Loop ->',3I4,F12.4,I8,F12.6,F12.6);

						}
					}
				}
			}
		}
	}

	for(i=0;i<xyz_imax;i++){
	  for(j=0;j<xyz_jmax;j++){
		for(k=0;k<xyz_kmax;k++){
		  xyz_ct[i+j*IMAX+k*JMAX*IMAX]=int(realxyz_ct[i+j*IMAX+k*JMAX*IMAX]);
		}
	  }
	}

	/*
	"a debugging loop
	"j_ct=0;
	"LOOP[
	"  j_ct=j_ct+1;
	"  i_ct=0;
	"  i_xyz=0;
	"  LOOP[
	"    i_ct=i_ct+2;
	"    i_xyz=i_xyz+1;
	"    OUTPUT ct_data(i_ct-1,j_ct,1),ct_data(i_ct,j_ct,1),
	"           xyz_ct(i_xyz,j_ct,1);
	"     (' ct_1,ct_2,xyz ',3I8);
	"  ]UNTIL (i_ct=256);
	"]UNTIL (j_ct=256);
    */

	free(realxyz_ct);

}

//    "ctcreate.mortran - start of subroutine CTToMedium"
//******************************************************************************
//
void CTToMedium(int &new_x_dim, int &new_y_dim, int &new_z_dim,
						  short *New_CT_Data,
						  int &num_material,
						  int *material_region,
						  float *density_region,
						  string *material_name){
	//---------------------------------------------------------------------------"
	//A subroutine for converting CT Hounsfield numbers into density data.
	//The user specifies the number of media (<= $MXMED), the type     
	//				 "potential CT offset value        "
	//of media, and the CT-density ramp for each medium.
	//The total ramp is assumed continuous in CT number, but may 
	//have discontinuities in density where one medium ends and the next begins. 
	//If the user enters 0 for the number of media, a default, hard-wired 
	//CT-density ramp is used.  The default ramp uses 4 media.  Note that in 
	//previous versions of ctcreate, this ramp was optimized for Pinnacle format 
	//data, which has CT no. > 0, but has since been shifted down by 1024 to 
	//the CT range of typical DICOM data
	//
	//  medium                  (CTmax-CTmin)/(max density - min density)
	//  ------                  --------------------------------------
	// AIR700ICRU                       (-974 - -1024)/(0.044-0.001)
	// LUNG700ICRU                      (-724 - -974)/(0.302-0.044)
	// ICRUTISSUE700ICRU                (101 - -724)/(1.101-0.302)
	// ICRPBONE700ICRU                  (1976 - 101)/(2.088-1.101)
	//
	//short New_CT_Data(IMAX,JMAX,KMAX);    "resized CT data                  "
	//int   new_x_dim;                      "resized number of x voxels       "
	//int   new_y_dim;                      "resized number of y voxels       "
	//int   new_z_dim;                      "resized number of z voxels       "
	//int   material_region($MXREG);        "matrix of ints identifying media "
	//int   num_material;                   "the number of materials in the   "
	//                                      " ramp.                           "
	int	  material_ct_upper_bound[MXMED]; //"ct upper bounds of ramps         "
	int	  material_ct_lower_bound;        //"min. ct number of ramp           "
	int	  i_material,I,J,K;               //"indices                          "
	int	  ct_low,ct_high;                 //"for outputting warnings          "
	//float density_region[MXREG];         //"matrix of material densities     "
	float material_density_lower_bound[MXMED]; //"lower bounds of ramps      "
	float material_density_upper_bound[MXMED]; //"upper bounds of ramps      "
	float dummy_input;
	char  ch;
	//CHARACTER*4 material_name(24,$MXMED);     "names of the media               " 
	//---------------------------------------------------------------------------"
	//REPLACE {$IRCTM(#,#,#)} WITH {(1 + {P1} + ({P2}-1)*new_x_dim + 
	//						   ({P3}-1)*new_x_dim*new_y_dim)};
	//---------------------------------------------------------------------------"
	//zero out the arrays first"
	for(I=0;I<new_x_dim;I++){
	   for(J=0;J<new_y_dim;J++){
		  for(K=0;K<new_z_dim;K++){
			 material_region[I+J*IMAX+K*JMAX*IMAX]=0;
			 density_region[I+J*IMAX+K*JMAX*IMAX]=0.0;
		  }
	   }
	}

	//Now, get the ramp info either interactively or using the hard-wired 
	//ramp function (user enters 0 for number of materials)"
	//This will eventually go into the subroutine GetCTConvRamp"
	//printf(" The CT-Density Ramp\n");
	//OUTPUT61;(/' The CT-Density Ramp'/
	//		' -------------------'/); 
	//:GETRAMPS:
	//printf(" Number of media (max %4d), min. CT number of ramp\n",MXMED);
	//printf(" (0,0 if you want to use the hard-wired ramp function): \n");
	//OUTPUT61 $MXMED;
	//(' Number of media (max ',I4,'), min. CT number of ramp'/
	// ' (0,0 if you want to use the hard-wired ramp function): ',$);
	FILE *f = fopen(ct_ramps_file, "r");
	//READ(5,*) num_material,material_ct_lower_bound;
	fscanf(f,"%d,%d\n",&num_material,&material_ct_lower_bound);
	//do{ ch = fgetc(f); }while(ch!='\n');
	//printf(" %5d, %5d\n",num_material,material_ct_lower_bound);
	//OUTPUT61 num_material,material_ct_lower_bound;(2I10);
	if(num_material<0 || num_material > MXMED){//error
		error(5005," Number of materials out of range. Try again. ");
		//OUTPUT61;(' Number of materials out of range. Try again. ');
		//STOP; 
	}
	else if(num_material==0){//use the hard-wired ramp"
		//printf(" Using the following default CT ramp.\n");
		//printf(" Note: This is optimized for display and example \n");
		//printf(" calculations.  It is recommended that you enter \n");
		//printf(" the CT conversion ramp for your own imager.\n");
		//OUTPUT61; (/' Using the following default CT ramp.'/
		//			' Note: This is optimized for display and example '/
		//			' calculations.  It is recommended that you enter '/
		//			' the CT conversion ramp for your own imager.'/);
		//DO i=1,4 [ DO j=1,24 [ material_name(j,i)=' '; ]; ];
		for(I=0;I<4;I++) material_name[I] = "";
		num_material=4;
		material_ct_lower_bound=-1024;

		//CT ramp information for Air"
		material_name[0] = "AIR700ICRU              ";
		//material_name(1,1)='A'; material_name(2,1)='I'; material_name(3,1)='R';
		//material_name(4,1)='7'; material_name(5,1)='0'; material_name(6,1)='0';
		//material_name(7,1)='I'; material_name(8,1)='C'; material_name(9,1)='R';
		//material_name(10,1)='U';
		material_ct_upper_bound[0] = - 974;
		material_density_lower_bound[0] = 0.001;
		material_density_upper_bound[0] = 0.044;
		//commented out density upper and lower bounds are for a direct output"
		//of Hounsfield numbers"
		//material_density_lower_bound(1)=-1024;"
		//material_density_upper_bound(1)=-974;"

		//CT ramp information for Lung"
		material_name[1] = "LUNG700ICRU             ";
		//material_name(1,2)='L'; material_name(2,2)='U'; material_name(3,2)='N';
		//material_name(4,2)='G'; material_name(5,2)='7'; material_name(6,2)='0';
		//material_name(7,2)='0'; material_name(8,2)='I'; material_name(9,2)='C';
		//material_name(10,2)='R'; material_name(11,2)='U'; 
		material_ct_upper_bound[1] = - 724;
		material_density_lower_bound[1] = 0.044;
		material_density_upper_bound[1] = 0.302;
		//commented out density upper and lower bounds are for a direct output"
		//of Hounsfield numbers"
		//material_density_lower_bound(2)=-974;"
		//material_density_upper_bound(2)=-724;"

		//CT ramp information for Tissue"
		material_name[2] = "ICRUTISSUE700ICRU       ";
		//material_name(1,3)='I'; material_name(2,3)='C'; material_name(3,3)='R';
		//material_name(4,3)='U'; material_name(5,3)='T'; material_name(6,3)='I';
		//material_name(7,3)='S'; material_name(8,3)='S'; material_name(9,3)='U';
		//material_name(10,3)='E'; material_name(11,3)='7'; material_name(12,3)='0';
		//material_name(13,3)='0'; material_name(14,3)='I'; material_name(15,3)='C';
		//material_name(16,3)='R'; material_name(17,3)='U';
		material_ct_upper_bound[2] = 101;
		material_density_lower_bound[2] = 0.302;
		material_density_upper_bound[2] = 1.101;
		//material_density_lower_bound(3)=-724;"
		//material_density_upper_bound(3)=101;"

		//CT ramp information for Bone"
		material_name[3] = "ICRPBONE700ICRU         ";
		//material_name(1,4)='I'; material_name(2,4)='C'; material_name(3,4)='R';
		//material_name(4,4)='P'; material_name(5,4)='B'; material_name(6,4)='O';
		//material_name(7,4)='N'; material_name(8,4)='E'; material_name(9,4)='7';
		//material_name(10,4)='0'; material_name(11,4)='0'; material_name(12,4)='I';
		//material_name(13,4)='C'; material_name(14,4)='R'; material_name(15,4)='U';
		material_ct_upper_bound[3] = 1976;
		material_density_lower_bound[3] = 1.101;
		material_density_upper_bound[3] = 2.088;
		//material_density_lower_bound(4)=101;"
		//material_density_upper_bound(4)=1976;"
		//printf("' CT no. lower bound of ramp = %5d\n",material_ct_lower_bound);
		//OUTPUT61 material_ct_lower_bound; 
		//(' CT no. lower bound of ramp = ',I5);
		for(i_material=0;i_material<num_material;i_material++){
			//printf(" Medium %4d: %s",i_material,material_name[i_material].c_str());
		   //OUTPUT61 i_material,(material_name(j,i_material),j=1,24);
			//	  (/' Medium ',I4,' : ',24a1);
			//printf(" CT no. upper bound, density lower bound (g/cm^3),\n");
			//printf(" density upper bound (g/cm^3)--all on one line:\n");
			/*
			printf(" %5d, %12.5f, %12.5f\n",material_ct_upper_bound[i_material],
				  material_density_lower_bound[i_material],
				  material_density_upper_bound[i_material]);
			*/
		   //OUTPUT61 material_ct_upper_bound(i_material),
			//	  material_density_lower_bound(i_material),
			//	  material_density_upper_bound(i_material);
		  //(' CT no. upper bound, density lower bound (g/cm^3),'/
		  // ' density upper bound (g/cm^3)--all on one line'/
		  // ' : ',I5,3F12.5);
		}
	}
	else{//user enters ramp"
		for(i_material=0;i_material<num_material;i_material++){
			ch = fgetc(f);
			for(I=0;I<24;I++){
				if(ch!='\n'){
					material_name[i_material] += ch;
					ch = fgetc(f);
				}
				else
					material_name[i_material] += ' ';
			}
			/*
			ch = fgetc(f);
			while(ch!='\n'){
				material_name[i_material] += ch;
				ch = fgetc(f);
			}
			*/
			fscanf(f,"%d,%f,%f,%f\n",&material_ct_upper_bound[i_material],
				  &material_density_lower_bound[i_material],
				  &material_density_upper_bound[i_material],
				  &dummy_input);
			//printf(" Medium %4d: %s\n",i_material,material_name[i_material].c_str());
			//printf(" CT no. upper bound, density lower bound (g/cm^3),\n");
			//printf(" density upper bound (g/cm^3)--all on one line:\n");
			/*
			printf(" %5d, %12.5f, %12.5f\n",material_ct_upper_bound[i_material],
				  material_density_lower_bound[i_material],
				  material_density_upper_bound[i_material]);
			*/
		   //OUTPUT61 i_material;(/' Medium ',I4,' : ',$);
		   //INPUT (material_name(j,i_material),j=1,24);(24a1);
		   //OUTPUT61 (material_name(j,i_material),j=1,24);(24a1);
		   //OUTPUT61;
		   //(' CT no. upper bound, density lower bound (g/cm^3),'/
			//' density upper bound (g/cm^3)--all on one line'/
			//' : ',$);
		   //READ(5,*) material_ct_upper_bound(i_material),
			//	  material_density_lower_bound(i_material),
			//	  material_density_upper_bound(i_material);
		   //OUTPUT61 material_ct_upper_bound(i_material),
			//	  material_density_lower_bound(i_material),
			//	  material_density_upper_bound(i_material);
			//	  (I5,2F12.5);
		   if(material_ct_upper_bound[i_material] < material_ct_lower_bound ||
			   material_density_lower_bound[i_material] < 0 ||
			   material_density_upper_bound[i_material] <
			   material_density_lower_bound[i_material]){
				   error(5006," CT no. or density out of range, or density upper bound <\n  density lower bound.  Try again.");
			//OUTPUT61;
			// (' CT no. or density out of range, or density upper bound <'/
			//  ' density lower bound.  Try again.');
			//	  STOP; 
		   }

		}//finished entering ramp values"
	}

	//Convert CT data into density data and input it into the correct regions

	ct_low=0;
	ct_high=0;

	for(K=0;K<new_z_dim;K++){
	   for(J=0;J<new_y_dim;J++){
		  for(I=0;I<new_x_dim;I++){
			 if(New_CT_Data[I+J*IMAX+K*JMAX*IMAX]<material_ct_lower_bound) {
				  //assume this will be vacuum"
				  New_CT_Data[I+J*IMAX+K*JMAX*IMAX]=0;
				  material_region[I+J*new_x_dim+K*new_y_dim*new_x_dim]=0;
				  density_region[I+J*new_x_dim+K*new_y_dim*new_x_dim]=0.0;
				  ct_low=ct_low+1;
			 }
			 else if(New_CT_Data[I+J*IMAX+K*JMAX*IMAX]>material_ct_upper_bound[num_material-1]) {
				  New_CT_Data[I+J*IMAX+K*JMAX*IMAX]=material_ct_upper_bound[num_material-1];
				  material_region[I+J*new_x_dim+K*new_y_dim*new_x_dim]=num_material;
				  density_region[I+J*new_x_dim+K*new_y_dim*new_x_dim]=
					  material_density_upper_bound[num_material-1];
				  ct_high=ct_high+1;
			 } 
			 else{
			   for(i_material=0;i_material<num_material;i_material++){
				  if(New_CT_Data[I+J*IMAX+K*JMAX*IMAX]<=material_ct_upper_bound[i_material]){
					material_region[I+J*new_x_dim+K*new_y_dim*new_x_dim]=i_material+1;
					if (i_material==0){ //on the first ramp"
					   density_region[I+J*new_x_dim+K*new_y_dim*new_x_dim]=
						   material_density_lower_bound[i_material]+
						   (New_CT_Data[I+J*IMAX+K*JMAX*IMAX]-material_ct_lower_bound)*
						   (material_density_upper_bound[i_material]-
							material_density_lower_bound[i_material])/
							(material_ct_upper_bound[i_material]-
							 material_ct_lower_bound);
					   material_region[I+J*new_x_dim+K*new_y_dim*new_x_dim]=1;
					}
					else{ //on a higher ramp" 
					   density_region[I+J*new_x_dim+K*new_y_dim*new_x_dim]=
						   material_density_lower_bound[i_material]+
					(New_CT_Data[I+J*IMAX+K*JMAX*IMAX]-material_ct_upper_bound[i_material-1])*
					  (material_density_upper_bound[i_material]-
					   material_density_lower_bound[i_material])/
					  (material_ct_upper_bound[i_material]-
					   material_ct_upper_bound[i_material-1]);
					   material_region[I+J*new_x_dim+K*new_y_dim*new_x_dim]=i_material+1;
					} 
					break;
				  }
			   }
			 } 
		  } 
	   } 
	}
	if(ct_low>0){
		//printf(" Warning: CT number in %10d voxels is < min. CT number of \n",ct_low);
		//printf(" ramp (%10d).  Medium in these voxels is set to 0 (vacuum).\n",material_ct_lower_bound);
	  //OUTPUT61 ct_low,material_ct_lower_bound;
	  //(/' Warning: CT number in ',I10,' voxels is < min. CT number of '/
		//' ramp (',I10,').  Medium in these voxels is set to 0 (vacuum).'/);
	}
	if(ct_high>0){
		//printf(" Warning: CT number in %10d voxels is > max. CT number of \n",ct_high);
		//printf(" ramp (%10d).  Medium in these voxels is set to medium no. %4d\n",material_ct_upper_bound[num_material-1],num_material);
	  //OUTPUT61 ct_high,material_ct_upper_bound(num_material),num_material;
	  //(/' Warning: CT number in ',I10,' voxels is > max. CT number of '/
		//' ramp (',I10,').  Medium in these voxels is set to medium no. ',I4,/);
	}
	//warning below is no longer necessary"
	/*
	DO I=1,new_x_dim[
	   DO J=1,new_y_dim[
		  DO K=1,new_z_dim[
			 if (New_CT_Data(i,j,k)=0 | 
				 material_region($IRCTM(I,J,K))=0 |
				 density_region($IRCTM(I,J,K))=0) [
			   OUTPUT61;(/' CT data and/or material no. and/or density is zero');
			   OUTPUT61 i,j,k,New_CT_Data(i,j,k);
					('i, j, k, CT no.      ',3I5,I8);
			   OUTPUT61 i,j,k,material_region($IRCTM(I,J,K));
					('i, j, k, material no.',3I5,I8);
			   OUTPUT61 i,j,k,density_region($IRCTM(I,J,K));
					('i, j, k, density     ',3I5,I8/);
			   STOP;
			 ];
		  ]
	   ]
	]
	*/
}

/*
//    "ctcreate.mortran - start of subroutine write_material"
//******************************************************************************
//---------------------------------------------------------------------------"
void write_material(fname,iimax,jjmax,kkmax,xbnd,ybnd,zbnd,
							  rho,med,CT,CTformat){
	"---------------------------------------------------------------------------"
	"This writes a file in .3ddose format with arrays of density, medium number "
	"and CT number for eventual display using PAW."
	;IMPLICIT NONE;
	Character*256 fname;               "name of header file for CT data"
	Character*80 ddataname;            "name of output file"
	Character*60 CTformat;             "format of original binary CT data"
	int iimax,                     "max number of x cells"
			jjmax,                     "max number of y cells"
			kkmax,                     "max number of z cells"
			lnblnk1,
			rindex,
			ii,jj,kk;                  "indices"
	float    xbnd(IMAX+1),             "voxel x boundaries"
			ybnd(JMAX+1),             "voxel y boundaries"
			zbnd(KMAX+1),             "voxel z boundaries"
			rho($MXREG);               "linear array of doses"
	int med($MXREG);               "linear array of dose uncertainties"
	short CT(IMAX,JMAX,KMAX);   "linear array of dose uncertainties"

	REPLACE {$PLOTOUT#;} WITH {write(15,*){P1};}
	REPLACE {$IRDWD(#,#,#)} WITH {({P1}+({P2}-1)*iimax+({P3}-1)*iimax*jjmax)}
	;
	"check format of CT data"
	if(CTformat='CADPLAN')[
	   ddataname=fname(rindex(fname,'/')+1:lnblnk1(fname)) // '.CTforPAW';
	]
	ELSEIF(CTformat='Pinnacle')[
		ddataname=fname(rindex(fname,'/')+1:index(fname,'header')-1) // 'CTforPAW';
	]
	ELSEIF(CTformat='DICOM')[
	   ddataname=fname(rindex(fname,'/')+1:lnblnk1(fname)-1) // '.CTforPAW';
	]
	OUTPUT61 ddataname(:lnblnk1(ddataname));(/' Writing CT phantom data into ',
	A,' for display.'/);

	Open (15,file=ddataname,Status='new',ERR=:CTforPAWOPENERROR:);

	$PLOTOUT iimax,jjmax,kkmax;
	$PLOTOUT (xbnd(ii),ii=1,iimax+1);
	$PLOTOUT (ybnd(jj),jj=1,jjmax+1);
	$PLOTOUT (zbnd(kk),kk=1,kkmax+1);
	$PLOTOUT (((rho($IRDWD(ii,jj,kk)),ii=1,iimax),jj=1,jjmax),kk=1,kkmax);
	$PLOTOUT (((CT(ii,jj,kk),ii=1,iimax),jj=1,jjmax),kk=1,kkmax);
	close(15);
	RETURN;
	:CTforPAWOPENERROR:
	OUTPUT61 ddataname(:lnblnk1(ddataname));
		   (//' ***ERROR: '/
		' Cannot write to ',A,'.  File already exists.'//);
}
*/

//    "ctcreate.mortran - start of subroutine write_phantom"
//******************************************************************************
//
void write_phantom(string baseName, int nmed, string *media,
					   int iimax, int jjmax, int kkmax, float *xbnd, float *ybnd, float *zbnd,
					   float *rho, int *med){
	//
	//This subroutine writes out CT phantom data for reading by dosxyznrc"
	//******************************************************************************
	//;IMPLICIT NONE;
	//CHARACTER*256 fname;               //file name containing CT data"
	char phantname[256];                 //output file"
	//Character*4  media(24,$MXMED);     //Media names"
	//char media_name[24];
	//Character*60 CTformat;             //format of original binary CT data"
	//int nmed;                      //number of media"
	//int iimax;                     //max number of x cells"
	//int jjmax;                     //max number of y cells"
	//int kkmax;                     //max number of z cells"
	//int	med[MXREG];               //linear array of media numbers"
	//int	lnblnk1;
	//int	rindex;
	int	ii,jj,kk;                  //indices"
	//float xbnd[IMAX+1];            //voxel x boundaries"
	//float ybnd[JMAX+1];            //voxel y boundaries"
	//float zbnd[KMAX+1];            //voxel z boundaries"
	//float rho[MXREG];              //linear array of densities"
	float estepe[MXMED];          //linear array of estepe values"


	//REPLACE {$IR(#,#,#)} WITH {(1 + {P1} + ({P2}-1)*iimax + ({P3}-1)*iimax*jjmax)}

	//;
	//open up the .egsphant file for outputting in the current directory"
	//find out where the last / is"
	/*
	if(CTformat='CADPLAN')[
	   phantname=fname(rindex(fname,'/')+1:lnblnk1(fname)) // '.egsphant';
	]
	ELSEIF(CTformat='Pinnacle')[
	   phantname=fname(rindex(fname,'/')+1:index(fname,'header')-1) // 'egsphant';
	]
	ELSEIF(CTformat='DICOM')[
	   phantname=fname(rindex(fname,'/')+1:lnblnk1(fname)-1) // '.egsphant';
	]
	OUTPUT61 phantname(:lnblnk1(phantname));(/' Writing CT phantom data into ',
	A,' to be read by dosxyznrc.'/);
	*/
	sprintf(phantname, "%s.egsphant",baseName.c_str());
	//printf("%s\n",phantname);
	
	//set estepe to the default value.  This is a dummy line in the .egsphant"
	//file, since this is controlled using the EGSnrc inputs." 
	for(ii=0;ii<nmed;ii++) estepe[ii]=1.0F;
 
	FILE *out = fopen(phantname, "w");
	//Open (15,file=phantname,Status='new',ERR=:egsphantOPENERROR:);
	fprintf(out,"%2d\n",nmed);
	//WRITE(15,'(i2)') nmed;
	for(ii=0;ii<nmed;ii++){
		fprintf(out,"%s\n",media[ii].c_str());
	   //Write(15,'(24a1)') (media(jj,ii),jj=1,24);
	}
	for(ii=0;ii<nmed;ii++){
		fprintf(out," %3.1f",estepe[ii]);
	}
	fprintf(out,"\n");
	//WRITE(15,*) (estepe(ii),ii=1,nmed);
	fprintf(out," %4d %4d %4d",iimax,jjmax,kkmax);
	fprintf(out,"\n");
	//WRITE(15,'(3i5)') iimax,jjmax,kkmax;
	for(ii=0;ii<=iimax;ii++){
		fprintf(out," %g",xbnd[ii]);
		if(ii%8==7) fprintf(out,"\n");
	}
	if(iimax%8!=7) fprintf(out,"\n");
	for(jj=0;jj<=jjmax;jj++){
		fprintf(out," %g",ybnd[jj]);
		if(jj%8==7) fprintf(out,"\n");
	}
	if(jjmax%8!=7) fprintf(out,"\n");
	for(kk=0;kk<=kkmax;kk++){
		fprintf(out," %g",zbnd[kk]);
		if(kk%8==7) fprintf(out,"\n");
	}
	if(kkmax%8!=7) fprintf(out,"\n");
	//WRITE(15,*) (xbnd(ii),ii=1,iimax+1);
	//WRITE(15,*) (ybnd(jj),jj=1,jjmax+1);
	//WRITE(15,*) (zbnd(kk),kk=1,kkmax+1);
	for(kk=0;kk<kkmax;kk++){
		for(jj=0;jj<jjmax;jj++){
			for(ii=0;ii<iimax;ii++){
				fprintf(out,"%1d",med[ii+jj*iimax+kk*jjmax*iimax]);
			}
			fprintf(out,"\n");
			//WRITE(15,1399) (med($IR(ii,jj,kk)),ii=1,iimax);
		}
		fprintf(out,"\n");
		//WRITE(15,*);
	}
	//1399 FORMAT(IMAXi1);
	for(kk=0;kk<kkmax;kk++){
		for(jj=0;jj<jjmax;jj++){
			for(ii=0;ii<iimax;ii++){
				fprintf(out,"  %g",rho[ii+jj*iimax+kk*jjmax*iimax]);
				if((ii+jj*iimax)%8==7) fprintf(out,"\n");
			}
			//fprintf(out,"\n");
			//WRITE(15,*) (rho($IR(ii,jj,kk)),ii=1,iimax);
		}
		if((iimax*jjmax)%8!=0) fprintf(out,"\n");
		fprintf(out,"\n");
		//WRITE(15,*);
	}
	fclose(out);
	//Close(15);
	//RETURN;
	/*
	:egsphantOPENERROR:
	OUTPUT61 phantname(:lnblnk1(phantname));
		   (//' ***ERROR: '/
		' Cannot write to ',A,'.  File already exists.'//);
	end;
	;
	*/
	/*
	FUNCTION rindex(c,a);
	CHARACTER c*(*);
	CHARACTER a*1;
	int j,rindex;
	DO j=LEN(c),1,-1[
	  if (c(j:j) = a) [
		 rindex=j;
		 RETURN;
	  ] 
	]
	rindex=0;
	RETURN;
	*/
}

void ctcreate(string rootName){

	//Character*60 ctformat;                 "format of CT data"
	//Character*40 machine;                  "type of machine being run on"

	//char    CTFileName[256];            //main CT file name"

	//int     nmed;                       //Number of media
	//string  *media = new string[MXMED];
	//char    media[MXMED][24];           //Media names
	//float   rhor[MXREG];                //Density distribution
	//int     med[MXREG];                 //Media distribution
	//float   xbound[IMAX+1];             //Voxel X bounds for dosxyznrc
	//float   ybound[JMAX+1];             //Voxel Y bounds for dosxyznrc
	//float   zbound[KMAX+1];             //Voxel Z bounds for dosxyznrc

	//CT Data Variables.
	//int     CTArraySize[3];                //Size of CT array
	//float   CTVoxelSize[3];                //Original CT voxel size
	//float   CTOffset[3];                   //The posn of the (1,1,1) voxel's corner
	//short   CTHounsData[CTIMAX*CTJMAX*CTKMAX];  //The original CT data
	//short   CTResizeData[IMAX*JMAX*KMAX];       //The resampled CT data
	//int     CTErrorCode;                   //Error code for CT"
	//int     imax,jmax,kmax;                //Max indices for resized CT data"

	//Integer lnblnk1,l;                     "have to declare this for implicit none"
	//integer iargc,narg;
	//character input_file*256;
	//logical unknown_format,file_exists;
	//bool    file_exists;


	/*
	$HAVE_C_COMPILER(#);

	narg = iargc();
	%F
	#if defined HAVE_C_COMPILER
	#define WITH_DICOM
	#endif
	%M

	IF( narg > 0 ) [
		call getarg(1,input_file);
		inquire(file=input_file,exist=file_exists);
		IF( ~file_exists ) [
			write(6,'(/a,a,a/)') 'File ',input_file(:lnblnk1(input_file)),
			  ' does not exist';
			$CALL_EXIT(1);
		]
		open(5,file=input_file);
		open(1,file=input_file(:lnblnk1(input_file))//'.ctlst');
	]
	ELSE [
		write(6,'(/a/)') 'No input file -> all input to come from the terminal';
		open(1,file='interactiv.ctlst');
	]
	*/
	//OUTPUT61; ('=============================================================');
	////printf(" Running ctcreate\n");
	//OUTPUT61; (' Running ctcreate');
	//OUTPUT61; ('=============================================================');

	//determine what format the CT data is in"
	/*
	OUTPUT61; (//' ************************************************************'/
				 ' '/
				 '           CT formats currently supported:'/
				 '           ------------------------------ '/
				 ' '/
				 '           1. Pinnacle                    '/
				 '           2. CADPLAN                     ');
	#ifdef WITH_DICOM;     
	OUTPUT61;  ( '           3. DICOM                       ');
	#endif;
	OUTPUT61; (//' ************************************************************');  
	OUTPUT61; (//' Input the format of your CT data'/ ' : ',$);
	read(5,'(A60)') ctformat; 
	OUTPUT61 ctformat(:lnblnk1(ctformat));(A60);
	unknown_format = .true.;
	IF(ctformat='pinnacle'|ctformat='PINNACLE'|ctformat='Pinnacle')[
	   ctformat='Pinnacle'; unknown_format = .false.;
	]
	ELSEIF(ctformat='cadplan'|ctformat='Cadplan'|ctformat='CadPlan'|
		   ctformat='CADPLAN')[
		   ctformat='CADPLAN'; unknown_format = .false.;
	]
	ELSEIF(ctformat='aapm'| ctformat='AAPM')[
		 OUTPUT61; (//' Convert CT data from AAPM format to Pinnacle format using'/
		 ' $OMEGA_HOME/dosxyznrc/CT/aapm2pinnacle.'//);
		 $CALL_EXIT(1);
	]

	#ifdef WITH_DICOM;     
	IF(ctformat='dicom'|ctformat='Dicom'|ctformat='DICOM'|ctformat='3')[
	  ctformat='DICOM'; unknown_format = .false.;
	]
	#else;
	IF(ctformat='dicom'|ctformat='Dicom'|ctformat='DICOM'|ctformat='3')[
	  OUTPUT; 
		(//' The DICOM routines are written in C. You need a working'/
		   ' C compiler to use them'//);
	  $CALL_EXIT(1);
	]
	#endif;

	IF( unknown_format ) [
	   OUTPUT ctformat;(//' CT data format ',a,' not currently handled.'//);
	   $CALL_EXIT(1);
	] 
	*/

	//Now filename is read in from .ini file
	//Read in the header filename from the egsinp file.
	/*
	IF(ctformat='Pinnacle')[
		OUTPUT61; (//' Input the full name of the header file for the CT data'/
			   ' : ',$);
	]
	ELSEIF(ctformat='CADPLAN')[
		OUTPUT61; (//' Input the full name of the file of CT data file names'/
			   ' : ',$);
	]
	ELSEIF(ctformat='DICOM')[
		OUTPUT61; (//' Input the full name of the file of DICOM file names'/
			   ' : ',$);
	]
	read(5,'(A256)') CTFilename;
	l = lnblnk1(CTFilename);
	OUTPUT61 CTFilename(:l); (A);
	*/

	//Call ReadCT to get original ct volume.              "

	//==========================================================="
	//  Input the CT data using various format types"
	//==========================================================="

	/*
	IF(ctformat='Pinnacle')[
	   Call ReadCT_Pinnacle(CTFileName, CTArraySize, CTHounsData, 
							CTOffset, CTVoxelSize,CTErrorCode);
	]
	"==========================================================="

	ELSEIF(ctformat='CADPLAN')[
	   Call ReadCT_CADPLAN(CTFileName, CTArraySize, CTHounsData, 
							CTOffset, CTVoxelSize,CTErrorCode);
	]
	"==========================================================="
	#ifdef WITH_DICOM;
	IF(ctformat='DICOM')[
	   CTFilename(l+1:l+1) = char(0);
	   Call ReadCT_DICOM(CTFileName, CTArraySize, CTHounsData,
						 CTOffset, CTVoxelSize,CTErrorCode);
	]
	#endif;
	*/

	//CT Data Variables.
	int     CTArraySize[3];                //Size of CT array
	float   CTVoxelSize[3];                //Original CT voxel size
	float   CTOffset[3];                   //The posn of the (1,1,1) voxel's corner
	int     CTErrorCode;                   //Error code for CT"
	short   *CTHounsData = (short*)malloc(CT_REG*sizeof(short));  //The original CT data

	//Now we only support dicom format
	readct_dicom(CTArraySize, CTHounsData, CTOffset, CTVoxelSize, &CTErrorCode);

	/*
	//pass out CT information for interpolation
	interp.start = make_float3(CTOffset[0], CTOffset[1], CTOffset[2]);
	logstr("  CT starting point . . . . . %g, %g, %g\n", interp.start.x, interp.start.y, interp.start.z);

	interp.dimen = make_int3(CTArraySize[0], CTArraySize[1], CTArraySize[2]);
	logstr("  CT xyz dimensions . . . . . %d, %d, %d\n", interp.dimen.x, interp.dimen.y, interp.dimen.z);

	interp.sizes = make_float3(CTVoxelSize[1], CTVoxelSize[1], CTVoxelSize[2]);
	logstr("  CT 3d voxel sizes . . . . . %g, %g, %g\n", interp.sizes.x, interp.sizes.y, interp.sizes.z);
	interp.start.x += interp.sizes.x/2;
	interp.start.y += interp.sizes.y/2;
	interp.start.z += interp.sizes.z/2;
	*/

	//check to see if CT array is too big to handle
	if(CTArraySize[0]>CTIMAX || CTArraySize[1]>CTJMAX || CTArraySize[2]>CTKMAX){
		  error(5001, "CT data array is %4d, %4d, %4d\n Max. ctcreate can deal with is %4d, %4d, %4d\n",
			  CTArraySize[0],CTArraySize[1],CTArraySize[2],CTIMAX,CTJMAX,CTKMAX);
	}
	/*
	IF(CTArraySize(1)>$CTIMAX | CTArraySize(2)>$CTJMAX | 
	   CTArraySize(3)>$CTKMAX)[
	   OUTPUT CTArraySize(1),CTArraySize(2),CTArraySize(3),
			  $CTIMAX,$CTJMAX,$CTKMAX;
	   (//' ***ERROR:'/
		  ' CT data array is ',I4,'x',I4,'x',I4/
		  ' Max. ctcreate can deal with is ',I4,'x',I4,'x',I4/
		  ' Go into ctcreate and change $CTIMAX,$CTJMAX,$CTKMAX to deal'/
		  ' with the size of the CT data array, recompile and try again.'//);
	   STOP;
	]
	*/

	//===========================================================
	//  add another ELSEIF block for your format
	//  and then send us the routine so we can make it available
	//===========================================================

	// Call subsizeCT to get subset of CT volume.
	SubsizeCT(CTArraySize,CTHounsData,CTVoxelSize,CTOffset);

	//====================================================================

	int imax=IMAX;  // Temporary, this will be reset in ResampleCT"
	int jmax=JMAX;
	int kmax=KMAX;
	float   xbound[IMAX+1];             //Voxel X bounds for dosxyznrc
	float   ybound[JMAX+1];             //Voxel Y bounds for dosxyznrc
	float   zbound[KMAX+1];             //Voxel Z bounds for dosxyznrc
	short *CTResizeData = (short*)malloc(XYZREG*sizeof(short));

	// Resample data based on user input for phantom voxels
	//               note there is user input within the routine
	ResampleCT(CTArraySize[0],CTArraySize[1],CTArraySize[2],
						 CTVoxelSize[0],CTVoxelSize[1],CTVoxelSize[2],
						 CTHounsData,
						 imax,jmax,kmax,
						 xbound,ybound,zbound,
						 CTResizeData,
						 CTOffset);
	free(CTHounsData);

	//====================================================================

	int    nmed;
	float  *rhor  = (float*)malloc(XYZREG*sizeof(float));
	int    *med   =   (int*)malloc(XYZREG*sizeof(int));
	string *media = new string[MXMED];
	//Convert the CT data to arrays of MED and RHOR values for each voxel
	CTToMedium(imax,jmax,kmax,              //input array dimensions
						CTResizeData,       //input CT data
						nmed,               //output number of media
						med,                //output media array
						rhor,               //output densities
						media);             //output names of media


	//====================================================================

	//output CT data that can be read by PAW for display
	//comment this out if you do not use PAW
	//the output file for PAW is CTFileName.CTforPAW

	//call write_material(CTFileName,imax,jmax,kmax,
	//                    xbound,ybound,zbound,
	//                    rhor,med,CTResizeData,ctformat);

	//====================================================================

	//now output the CT phantom in a form that can be read by dosxyznrc
	//the phantom will be in file, CTFileName.egsphant
	//Note: estepe is no longer passed 
	write_phantom(rootName,nmed,media,
					   imax,jmax,kmax,xbound,ybound,zbound,
					   rhor,med);

	//====================================================================

	free(CTResizeData);
	free(rhor);
	free(med);

}

//******************r************************************************************
//     end of ctcreate.c
