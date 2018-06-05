



#ifdef CUDA_EGS

#include "random.cu"
//int mscati(void);

//int MLC_med_IDs[2];  //the medium IDs for the 2 media of MLC


//void free_all();
//void read_egsinp(string *media_names);

// remove all control characters (possibly including spaces) from a string
string trim(string str, bool allowSpace) {
    string trimmed;
    for (uint i = 0; i < str.length(); i++) {
        char c = str[i];
        if ((c > 32) || (allowSpace && (c == 32)))
            trimmed += c;
    }
    return trimmed;
}


// allocate memory for the stack and associated counters
uint init_stack(uchar nRepeat) {
	
	float relative_stack_size[NUM_SIMU_CAT] = {1.5, 0.9, 1.2, 0.3, 0.1, 1.6, 0.5, 0.6, 0.1, 0.1};
	//the relative stack size for the 10 stacks: photon_step,compton,photo,rayleigh,pair,electr_step,moller,brems,bhabha,annih

	uint actual_batch_size = PHSP_STACK_SIZE;
	nRepeat += 1;  //include the first use of that particle
	if(actual_batch_size*nRepeat*relative_stack_size[p_photon_step]>=MAX_STACK_SIZE){
		actual_batch_size = (uint)(MAX_STACK_SIZE/(nRepeat*relative_stack_size[p_photon_step]));
		//logstr("actual batch size is adjusted to %d.\n",actual_batch_size);
	}
	//else
		//logstr("actual batch size is %d.\n",actual_batch_size);

	uint stack_size;
	int h_np_max[NUM_SIMU_CAT];
	int h_np[NUM_SIMU_CAT];

	d_stack = (stack_t*)malloc(GPUNo*NUM_SIMU_CAT*sizeof(stack_t));

    for(uchar i=0;i<NUM_SIMU_CAT;i++){
		/*
		if(i==p_photon_step)
			stack_size = (uint)((actual_batch_size*(nRepeat+1))*relative_stack_size[i]);
		else
			stack_size = (uint)((actual_batch_size*nRepeat)*relative_stack_size[i]);
		*/
		stack_size = (uint)((actual_batch_size*nRepeat)*relative_stack_size[i]);
		//logstr("stack_size[%d] = %d.\n",i,stack_size);
		/*
		if(i==p_photon_step||i==p_compton)
			stack_size = HUGE_STACK_SIZE;
		else if(i==e_electron_step)
			stack_size = MAIN_STACK_SIZE;
		else if(i==p_photo||i==p_rayleigh||i==e_moller||e_brems)
			stack_size = HALF_STACK_SIZE;
		else stack_size = PHSP_STACK_SIZE;
		*/
		for(int GPUId=0; GPUId<GPUNo; GPUId++) {
#ifdef USE_MULTIPLE_GPU
			cudaSetDevice(GPUId); ce(58017);
#endif
			cudaMalloc(&d_stack[i+GPUId*NUM_SIMU_CAT],stack_size*sizeof(particle_t)); ce(4051);
			//cudaMalloc(&d_stack[i].a,stack_size*sizeof(uint4)); ce(4051);
			//cudaMalloc(&d_stack[i].b,stack_size*sizeof(uint4)); ce(4052);
			//cudaMalloc(&d_stack[i].c,stack_size*sizeof(uint4)); ce(4053);
		}
		h_np_max[i] = stack_size;
		h_np[i] = -1;
    }
	for(int GPUId=0; GPUId<GPUNo; GPUId++) {
#ifdef USE_MULTIPLE_GPU
		cudaSetDevice(GPUId); ce(58018);
#endif
		cudaMemcpyToSymbol(stack,&d_stack[GPUId*NUM_SIMU_CAT],NUM_SIMU_CAT*sizeof(stack_t)); ce(4054);
		cudaMemcpyToSymbol(np_max,h_np_max,NUM_SIMU_CAT*sizeof(int)); ce(4055);
		cudaMemcpyToSymbol(np,h_np,NUM_SIMU_CAT*sizeof(int)); ce(4056);
	}

	/*
    // allocate step counts on host
    h_total_step_counts = (total_step_counts_t*)malloc(sizeof(total_step_counts_t));
	d_total_step_counts = (total_step_counts_t**)malloc(GPUNo*sizeof(total_step_counts_t*));

    // init step counts on the device
	for(int GPUId=0; GPUId<GPUNo; GPUId++) {
#ifdef USE_MULTIPLE_GPU
		cudaSetDevice(GPUId); ce(58019);
#endif
		cudaMalloc(&d_total_step_counts[GPUId], sizeof(total_step_counts_t)); ce(4005);
		cudaMemset(d_total_step_counts[GPUId], 0, sizeof(total_step_counts_t)); ce(4006);
		cudaMemcpyToSymbol(total_step_counts, &d_total_step_counts[GPUId], sizeof(total_step_counts_t*)); ce(4007);
	}
	*/

	return actual_batch_size; // the maximum number of primary particles allowed to load at each batch
}

void free_stack(){

	for(int GPUId=0; GPUId<GPUNo; GPUId++) {

#ifdef USE_MULTIPLE_GPU
		cudaSetDevice(GPUId); ce(58014);
#endif

		for(uchar i=0;i<NUM_SIMU_CAT;i++){
			cudaFree(d_stack[i+GPUId*NUM_SIMU_CAT]); ce(9051);
			//cudaFree(d_stack[i].a); ce(9051);
			//cudaFree(d_stack[i].b); ce(9052);
			//cudaFree(d_stack[i].c); ce(9053);
		}
	}
}

/*
// read the source parameters from the input file and copy the relevant data to the device
void init_source() {

    // get Phase space file
    GetPrivateProfileString("Source", "phase space file", "<<missing>>", charBuffer, CHARLEN, input_file);
	string phsp_file_input = string(charBuffer);
	if (phsp_file_input == "<<missing>>")
		error(2012, "Property \"pegs file\" is not speciefied in the section \"General\" in the input file \"%s\".\n", input_file);
    
    // get the absolute path of the pegs file
    GetFullPathName(phsp_file_input.c_str(), CHARLEN, charBuffer, NULL);
    phsp_file_input = string(charBuffer);
    phsp_file = phsp_file_input.c_str();
    logstr("  Phase space file  . . . . . %s\n", phsp_file);
	phsp = read_egsphsp_header(phsp_file, phsp_header);
	if(!phsp) 
	{
		logstr("Read phase space header error! \n");
		exit(-1);
	}

	
	// read source energy from the input file
	GetPrivateProfileString("Source", "energy", "<<missing>>", charBuffer, CHARLEN, input_file);
	if (string(charBuffer) == "<<missing>>")
		error(5010, "Property \"energy\" is not speciefied in the section \"Source\" in the input file \"%s\".\n", input_file);
	float e;
	if (sscanf(charBuffer, "%f", &e) != 1)
		error(5011, "Could not parse the source energy \"%s\" in the input file \"%s\".\n", charBuffer, input_file);
	if (e <= 0.0F)
		error(5012, "The source energy must be positive, but %f was specified in the input file \"%s\".\n", e, input_file);
	h_source.energy = e;

	// read source position from the input file
	GetPrivateProfileString("Source", "position", "<<missing>>", charBuffer, CHARLEN, input_file);
	if (string(charBuffer) == "<<missing>>")
		error(5013, "Property \"position\" is not speciefied in the section \"Source\" in the input file \"%s\".\n", input_file);
	float p1, p2, p3;
	if (sscanf(charBuffer, "%f %f %f", &p1, &p2, &p3) != 3)
		error(5014, "Could not parse the source position \"%s\" in the input file \"%s\".\n", charBuffer, input_file);
    h_source.source_point = make_float3(p1, p2, p3);

	// read collimator rectangle from the input file
	GetPrivateProfileString("Source", "collimator rectangle", "<<missing>>", charBuffer, CHARLEN, input_file);
	if (string(charBuffer) == "<<missing>>")
		error(5015, "Property \"collimator rectangle\" is not speciefied in the section \"Source\" in the input file \"%s\".\n", input_file);
	float r1, r2, r3, r4;
	if (sscanf(charBuffer, "%f %f %f %f", &r1, &r2, &r3, &r4) != 4)
		error(5016, "Could not parse the source collimator rectangle \"%s\" in the input file \"%s\".\n", charBuffer, input_file);
    if (r3 <= r1)
		error(5017, "The maximum x coordinate %f of the source collimator must be greater than the minimum x coordinate %f in the input file \"%s\".\n", r3, r1, input_file);
	if (r4 <= r2)
		error(5018, "The maximum y coordinate %f of the source collimator must be greater than the minimum y coordinate %f in the input file \"%s\".\n", r4, r2, input_file);
	h_source.rectangle_min = make_float2(r1, r2);
    h_source.rectangle_max = make_float2(r3, r4);

	// read source collimator z coordinate from the input file
	GetPrivateProfileString("Source", "collimator z", "<<missing>>", charBuffer, CHARLEN, input_file);
	if (string(charBuffer) == "<<missing>>")
		error(5019, "Property \"collimator z\" is not speciefied in the section \"Source\" in the input file \"%s\".\n", input_file);
	float z;
	if (sscanf(charBuffer, "%f", &z) != 1)
		error(5020, "Could not parse the source collimator z coordinate \"%s\" in the input file \"%s\".\n", charBuffer, input_file);
    h_source.rectangle_z = z;
    
    h_source.rectangle_size = make_float2(h_source.rectangle_max.x - h_source.rectangle_min.x, 
                                          h_source.rectangle_max.y - h_source.rectangle_min.y);
    h_source.rectangle_area = h_source.rectangle_size.x * h_source.rectangle_size.y;

    cudaMemcpyToSymbol(source, &h_source, sizeof(source)); ce(5027);

    logstr("\nSource\n");
    logstr("  Energy  . . . . . . . . . . %f\n", h_source.energy);
    logstr("  Position  . . . . . . . . . x = %f, y = %f, z = %f\n", h_source.source_point.x, h_source.source_point.y, h_source.source_point.z);
    logstr("  Collimator\n");
    logstr("    x . . . . . . . . . . . . min = %f, max = %f\n", h_source.rectangle_min.x, h_source.rectangle_max.x);
    logstr("    y . . . . . . . . . . . . min = %f, max = %f\n", h_source.rectangle_min.y, h_source.rectangle_max.y);
    logstr("    z . . . . . . . . . . . . %f\n", h_source.rectangle_z);

}
*/

/*
// read the detector parameters from the input file and copy the relevant data to the device
void init_detector() {
    // read detector position from the input file
	GetPrivateProfileString("Detector", "position", "<<missing>>", charBuffer, CHARLEN, input_file);
	if (string(charBuffer) == "<<missing>>")
		error(6001, "Property \"position\" is not speciefied in the section \"Detector\" in the input file \"%s\".\n", input_file);
	float p1, p2, p3;
	if (sscanf(charBuffer, "%f %f %f", &p1, &p2, &p3) != 3)
		error(6002, "Could not parse the detector position \"%s\" in the input file \"%s\".\n", charBuffer, input_file);
    h_detector.center = make_float3(p1, p2, p3);
	
	// read detector size from the input file
	GetPrivateProfileString("Detector", "size", "<<missing>>", charBuffer, CHARLEN, input_file);
	if (string(charBuffer) == "<<missing>>")
		error(6003, "Property \"size\" is not speciefied in the section \"Detector\" in the input file \"%s\".\n", input_file);
	int s1, s2;
	if (sscanf(charBuffer, "%d %d", &s1, &s2) != 2)
		error(6004, "Could not parse the detector size \"%s\" in the input file \"%s\".\n", charBuffer, input_file);
	if (s1 <= 0)
		error(6005, "The number of pixels in the x direction of the detector must be positive, but %d was specified in the input file \"%s\".\n", s1, input_file);
    if (s2 <= 0)
		error(6006, "The number of pixels in the y direction of the detector must be positive, but %d was specified in the input file \"%s\".\n", s2, input_file);
	h_detector.N = make_uint2(s1, s2);

	// read detector pixel size from the input file
	GetPrivateProfileString("Detector", "pixel size", "<<missing>>", charBuffer, CHARLEN, input_file);
	if (string(charBuffer) == "<<missing>>")
		error(6007, "Property \"pixel size\" is not speciefied in the section \"Detector\" in the input file \"%s\".\n", input_file);
	float ps1, ps2;
	if (sscanf(charBuffer, "%f %f", &ps1, &ps2) != 2)
		error(6008, "Could not parse the detector pixel size \"%s\" in the input file \"%s\".\n", charBuffer, input_file);
	if (s1 <= 0)
		error(6009, "The detector pixels size in the x direction must be positive, but %f was specified in the input file \"%s\".\n", ps1, input_file);
    if (s2 <= 0)
		error(6010, "The detector pixels size in the y direction must be positive, but %f was specified in the input file \"%s\".\n", ps2, input_file);
	h_detector.d = make_float2(ps1, ps2);

	// copy detector parameters to the device
    cudaMemcpyToSymbol(detector, &h_detector, sizeof(detector)); ce(6011);

    logstr("\nDetector\n");
    logstr("  Position  . . . . . . . . . x = %f, y = %f, z = %f\n", h_detector.center.x, h_detector.center.y, h_detector.center.z);
    logstr("  Size (pixels) . . . . . . . x = %d, y = %d\n", h_detector.N.x, h_detector.N.y);
    logstr("  Pixel size  . . . . . . . . x = %f, y = %f\n", h_detector.d.x, h_detector.d.y);
}
*/

// initialize the leaf end shape
void leaf_end()
{
	//float zmin = h_crossXZ[15].x;
	//float zthick=h_crossXZ[15].y;            //h_crossXZ[15].y=h_zthick;
    //float zcenter= zthick/2;
    //float rmax=h_crossXZ[111].x;             //h_crossXZ[111].x=h_rmax;
	float zmin = VarianMLC.zmin;
	float zthick = VarianMLC.zthick;
	float zcenter = zthick/2;
	float rmax = VarianMLC.rmax;

	point2D n[30];
	//leaf in bank A, with position NEG
	n[0].y=-zcenter;  n[0].x= -0.86708;  //n[0].y=-3.38095;  n[0].x= -0.86708;
	n[1].y=-2.92273;  n[1].x= -0.64016;
	n[2].y=-2.45316;  n[2].x= -0.44640;
	n[3].y=-1.97406;  n[3].x= -0.28667;
	n[4].y=-1.48729;  n[4].x= -0.16167;
	n[5].y=-0.99475;  n[5].x= -0.07199;
	n[6].y=-0.49834;  n[6].x= -0.01802;
	n[7].y= 0.00000;  n[7].x=  0.00000;
	n[8].y= 0.49834;  n[8].x= -0.01802;
	n[9].y= 0.99475;  n[9].x= -0.07199;
	n[10].y=1.48729;  n[10].x=-0.16167;
	n[11].y=1.97406;  n[11].x=-0.28667;
	n[12].y=2.45316;  n[12].x=-0.44640;
	n[13].y=2.92273;  n[13].x=-0.64016;
	n[14].y=zcenter;  n[14].x=-0.86708;  //n[14].y=3.38095;  n[14].x=-0.86708;

	//leaf in bank B, with position POS
	n[29].y=-zcenter;  n[29].x=0.86708;  //n[29].y=-3.38095;  n[29].x=0.86708;
	n[28].y=-2.92273;  n[28].x=0.64016;
	n[27].y=-2.45316;  n[27].x=0.44640;
	n[26].y=-1.97406;  n[26].x=0.28667;
	n[25].y=-1.48729;  n[25].x=0.16167;
	n[24].y=-0.99475;  n[24].x=0.07199;
	n[23].y=-0.49834;  n[23].x=0.01802;
	n[22].y= 0.00000;  n[22].x=0.00000;
	n[21].y= 0.49834;  n[21].x=0.01802;
	n[20].y= 0.99475;  n[20].x=0.07199;
	n[19].y= 1.48729;  n[19].x=0.16167;
	n[18].y= 1.97406;  n[18].x=0.28667;
	n[17].y= 2.45316;  n[17].x=0.44640;
	n[16].y= 2.92273;  n[16].x=0.64016;
	n[15].y= zcenter;  n[15].x=0.86708;  //n[15].y= 3.38095;  n[15].x=0.86708;

	//point2D *h_crossYZ=(point2D*)malloc(34*sizeof(point2D));

    int i;
    h_crossYZ[0].x=-rmax;
    h_crossYZ[0].y= zthick+zmin;
    h_crossYZ[1].x=-rmax;
    h_crossYZ[1].y= zmin;
    for(i=0;i<15;i++)
    {
        h_crossYZ[2+i].x=n[i].x;
        h_crossYZ[2+i].y=n[i].y+zcenter+zmin;
    }
    //h_crossYZ[2].y=0;
    //h_crossYZ[16].y= zthick;

    for(i=0;i<15;i++)
    {
        h_crossYZ[17+i].x=n[15+i].x;
        h_crossYZ[17+i].y=n[15+i].y+zcenter+zmin;
    }
    //h_crossYZ[17].y= zthick;
    //h_crossYZ[31].y=0;
    h_crossYZ[32].x=rmax;
    h_crossYZ[32].y=zmin;
    h_crossYZ[33].x=rmax;
    h_crossYZ[33].y=zthick+zmin;

	d_crossYZ = (point2D**)malloc(GPUNo*sizeof(point2D*));

	for(int GPUId=0; GPUId<GPUNo; GPUId++) {
	
#ifdef USE_MULTIPLE_GPU
		cudaSetDevice(GPUId); ce(58004);
#endif

		cudaMalloc(&d_crossYZ[GPUId],34*sizeof(point2D));
		cudaMemcpy(d_crossYZ[GPUId],h_crossYZ,34*sizeof(point2D),cudaMemcpyHostToDevice);
	
	}

}

// read egsinp file to get crossXZ array and other information
void read_egsinp(string &media_names)
{

	int i,j,N;
    float X;
    char char_temp[100];
    char ch;
    char CM_title[32]="*********** start of CM DYNVMLC";  //31 characters
    int N_full;
    int N_half;
    float h_rad_end;
	float h_rmax,h_zmin, h_zthick, h_x_start, h_x_end, h_x_mid_1, h_x_mid_2, h_air_gap;
	//float wl_full, wl_target, wl_isocenter, wt;

	Leaf full,target,isocenter;
	//int nmed=2;

	point2D *h_crossXZ=(point2D*)malloc(960*sizeof(point2D));

	point2D h_boundaryXZ[4];
	point2D h_boundaryYZ[4];


	ifstream infile;
    infile.open(egsinp_file);
    if(infile.is_open())//check if file is opened
    {
        while(infile.good()&&!infile.eof())//while file is in good condition and not end of file
        {
            //memset(char_temp,0,100);  //should include <cstring>
            i=0;
            do
            {
                ch=infile.get();
                char_temp[i]=ch;
                i++;
            }
            while(ch!='\n'&&!infile.eof());//while in each line and not end of file
            i=0;
            while(char_temp[i]==CM_title[i]&&i<31)
            {
                i++;
            }
            if(i==31)
            {
                infile>>h_rmax>>ch;   infile.ignore(100,'\n');
                for(i=0;i<2;i++)
                {
                    infile.ignore(100,'\n');
                }
                infile>>h_zmin>>ch;   infile.ignore(100,'\n');
                infile>>h_zthick>>ch; infile.ignore(100,'\n');
                infile>>full.wl>>ch>>full.wt>>ch>>full.wg>>ch>>full.wtip>>ch>>full.wts>>ch>>full.wbs>>ch;
                infile>>full.ztip>>ch>>full.zl>>ch>>full.zt>>ch>>full.zg>>ch>>full.zth>>ch>>full.zbh>>ch;
                infile>>full.holepos>>ch;
                infile>>full.zts>>ch>>full.zbs>>ch;
                infile>>target.wl>>ch>>target.wt>>ch>>target.wg>>ch>>target.wtip>>ch>>target.wts>>ch>>target.wbs>>ch;
                infile>>target.zts>>ch>>target.zbs>>ch>>target.zth>>ch>>target.zbh>>ch;
                infile>>target.holepos>>ch;
                infile>>target.zt>>ch>>target.zg>>ch>>target.zl>>ch>>target.ztip>>ch;
                infile>>isocenter.wl>>ch>>isocenter.wt>>ch>>isocenter.wg>>ch>>isocenter.wtip>>ch>>isocenter.wts>>ch>>isocenter.wbs>>ch;
                infile>>isocenter.ztip>>ch>>isocenter.zl>>ch>>isocenter.zt>>ch>>isocenter.zg>>ch>>isocenter.zth>>ch>>isocenter.zbh>>ch;
                infile>>isocenter.holepos>>ch;
                infile>>isocenter.zts>>ch>>isocenter.zbs>>ch;

                infile>>N_full>>ch;   infile.ignore(100,'\n'); //number of full leaves on one side(both left and right are the same)
                infile>>N_half>>ch;   infile.ignore(100,'\n'); //total number of half leaves in the middle
                                      infile.ignore(100,'\n'); //skip number of full leaves on the right
                infile>>h_x_start>>ch;  infile.ignore(100,'\n');
                infile>>h_air_gap>>ch;  infile.ignore(100,'\n');
                                      infile.ignore(100,'\n'); //skip ENDTYPE. suppose end type is round
                infile>>h_rad_end>>ch;  infile.ignore(100,'\n');
                                      infile.ignore(100,'\n'); //skip ZFOCUS of leaf sides. suppose z_focus is 0

				link_field *list,*head,*before; //define a list to store the leaf bank. pointer *head is point to the first
                list=new link_field;
                infile>>list->NEG>>ch>>list->POS>>ch>>list->NUM;
                list->next=NULL;
                head=list;
                before=NULL;
                i=list->NUM;
                while(i<60)
                {
                    before=list;
                    list=new link_field;
                    before->next=list;
                    infile>>list->NEG>>ch>>list->POS>>ch>>list->NUM;
                    list->next=NULL;
                    i=i+list->NUM;
                }
				infile.ignore(100,'\n');  //skip the '\n' of the last line of leaf opening definition
				//string name;
				infile.ignore(100,'\n');
				infile.ignore(100,'\n');
				//for(i=0;i<nmed;i++)
				//{
					infile.ignore(100,'\n');
					j=0;
					ch=infile.get();
					while(ch!='\n'&&!infile.eof())//while in each line and not end of file
					{
						media_names+=ch;
						ch=infile.get();
						j++;
					}
				//}

                full.zmin=target.zmin=isocenter.zmin=h_zmin;
                full.zthick=target.zthick=isocenter.zthick=h_zthick;
                full.VerCal_full();
                target.VerCal_target();
                isocenter.VerCal_isocenter();

                //there are 60 leaves, 14 vertice for each leaf, we use an array of 16*60 vertice to save these positions
                //in order to adapt the process of data transportation of GPU, each time it transports a data block of 64 bytes,
                //we extend from 14 vertice to 16 for each leaf, so the data for each leaf is saved within one specific data block
                //here for each leaf, we save the position data in the first 14 vertice in array for each leaf
                X=h_x_start;    //h_x_start is the start position in x coordinate of leaves
                N=N_full;
                for(i=0;i<N;i++)
                {
                    for(j=0;j<14;j++)
                    {
                        h_crossXZ[i*16+j].x=(X+full.v[j].x)/h_zmin*full.v[j].z;
						h_crossXZ[i*16+j].y=full.v[j].z;
                        //h_crossXZ[i*16+j].y=full.v[j].z -h_zmin;         //we no longer substract the h_zmin, every thing is now relative to isocenter
                    }
                    X=X+full.wl+h_air_gap;
                }
                h_x_mid_1=X;    //h_x_mid_1 is the start position of half leaves
                N=N+N_half;
                for(;i<N;i++)
                {
                    for(j=0;j<14;j++)
                    {
                        h_crossXZ[i*16+j].x=(X+target.v[j].x)/h_zmin*target.v[j].z;
						h_crossXZ[i*16+j].y=target.v[j].z;
                        //h_crossXZ[i*16+j].y=target.v[j].z -h_zmin;
                    }
                    X=X+target.wl+h_air_gap;
                    i++;
                    for(j=0;j<14;j++)
                    {
                        h_crossXZ[i*16+j].x=(X+isocenter.v[j].x)/h_zmin*isocenter.v[j].z;
                        h_crossXZ[i*16+j].y=isocenter.v[j].z;
						//h_crossXZ[i*16+j].y=isocenter.v[j].z -h_zmin;
                    }
                    X=X+isocenter.wl+h_air_gap;
                }
                h_x_mid_2=X;    //h_x_mid_2 is the start position of second part of full leaves
                N=N+N_full;
                for(;i<N;i++)
                {
                    for(j=0;j<14;j++)
                    {
                        h_crossXZ[i*16+j].x=(X+full.v[j].x)/h_zmin*full.v[j].z;
                        h_crossXZ[i*16+j].y=full.v[j].z;
						//h_crossXZ[i*16+j].y=full.v[j].z -h_zmin;
                    }
                    X=X+full.wl+h_air_gap;
                }
                h_x_end=X-h_air_gap+full.wg;    //h_x_end is the end position of leaves

                //we use the fifteenth vertex of each leaf to save the open position of the leaf
                i=0;N=0;
                list=head;
                do
                {
                    N=N+list->NUM;
                    for(;i<N;i++)
                    {
                        h_crossXZ[14+i*16].x=list->NEG;
                        h_crossXZ[14+i*16].y=list->POS;
                    }
                    before=list;
                    list=list->next;
                    delete before;
                }
                while(list!=NULL);

                //we save some relevant data in the remaining vertice of the array
                h_crossXZ[15].x=h_zmin ;
                h_crossXZ[15].y=h_zthick;
                h_crossXZ[31].x=h_x_start;
                h_crossXZ[31].y=h_x_end;
                h_crossXZ[47].x=h_x_mid_1;
                h_crossXZ[47].y=h_x_mid_2;
                //h_crossXZ[63].x=h_rad_end;
                h_crossXZ[63].y=h_air_gap;
                h_crossXZ[79].x=full.wl;
                h_crossXZ[79].y=target.wl;
                h_crossXZ[95].x=isocenter.wl;
                h_crossXZ[95].y=full.wt;
                h_crossXZ[111].x=h_rmax;

			    VarianMLC.zmin=h_zmin ;
				VarianMLC.zthick=h_zthick;
                //VarianMLC.x_start=h_x_start;
				//VarianMLC.x_mid_1=h_x_mid_1;
				//VarianMLC.x_mid_2=h_x_mid_2;
				//VarianMLC.rad_end=h_rad_end;
				//VarianMLC.air_gap=h_air_gap;
				VarianMLC.rmax=h_rmax;

				break;
            }
        }
    }
    infile.close();//close the file

	// define leaf end shape in yz cross section
	leaf_end();

	//define the boundary in xz plane
    h_boundaryXZ[0].x=h_x_start/h_zmin*(h_zmin+h_zthick);
	h_boundaryXZ[0].y=h_zthick+h_zmin;  // add h_zmin
	h_boundaryXZ[1].x=h_x_start;
	h_boundaryXZ[1].y=h_zmin;
	h_boundaryXZ[2].x=h_x_end;
	h_boundaryXZ[2].y=h_zmin;
	h_boundaryXZ[3].x=h_x_end/h_zmin*(h_zmin+h_zthick);
	h_boundaryXZ[3].y=h_zthick+h_zmin;

	//define the boundary in xz plane
    h_boundaryYZ[0].x=-h_rmax;
	h_boundaryYZ[0].y=h_zthick+h_zmin;  // add h_zmin
	h_boundaryYZ[1].x=-h_rmax;
	h_boundaryYZ[1].y=h_zmin;
	h_boundaryYZ[2].x=h_rmax;
	h_boundaryYZ[2].y=h_zmin;
	h_boundaryYZ[3].x=h_rmax;
	h_boundaryYZ[3].y=h_zthick+h_zmin;

	for(int GPUId=0; GPUId<GPUNo; GPUId++) {

#ifdef USE_MULTIPLE_GPU
		cudaSetDevice(GPUId); ce(58005);
#endif

		cudaMemcpyToSymbol(MLC_rmax, &h_rmax, sizeof(float)); ce(17006);
		cudaMemcpyToSymbol(MLC_zmin, &h_zmin, sizeof(float)); ce(17006);
		cudaMemcpyToSymbol(MLC_zthick, &h_zthick, sizeof(float)); ce(17006);
		cudaMemcpyToSymbol(MLC_x_start, &h_x_start, sizeof(float)); ce(17006);
		cudaMemcpyToSymbol(MLC_x_end, &h_x_end, sizeof(float)); ce(17006);
		cudaMemcpyToSymbol(MLC_x_mid_1, &h_x_mid_1, sizeof(float)); ce(17006);
		cudaMemcpyToSymbol(MLC_x_mid_2, &h_x_mid_2, sizeof(float)); ce(17006);
		//cudaMemcpyToSymbol(MLC_rad_end, &h_rad_end, sizeof(float)); ce(17006);
		cudaMemcpyToSymbol(MLC_air_gap, &h_air_gap, sizeof(float)); ce(17006);
		cudaMemcpyToSymbol(MLC_wl_full, &full.wl, sizeof(float)); ce(17006);
		cudaMemcpyToSymbol(MLC_wl_target, &target.wl, sizeof(float)); ce(17006);
		cudaMemcpyToSymbol(MLC_wl_isocenter, &isocenter.wl, sizeof(float)); ce(17006);
		cudaMemcpyToSymbol(MLC_wt, &full.wt, sizeof(float)); ce(17006);

		cudaMemcpyToSymbol(crossXZ,h_crossXZ, 960*sizeof(point2D));

		cudaMemcpyToSymbol(boundaryXZ,h_boundaryXZ, 4*sizeof(point2D));

		cudaMemcpyToSymbol(boundaryYZ,h_boundaryYZ, 4*sizeof(point2D));
	
	}

	free(h_crossXZ);

    return;
}

// read egsinp file for SecJaws to get geometric data
void read_egsinp_secjaws(string &media_names)
{

	int i;
    char char_temp[100];
    char ch;
    char CM_title[53]="*********** start of CM JAWS with identifier SECJAWS";  //52 characters
    float h_rmax;
	float h_zminy, h_zmaxy;  // h_yfp, h_ybp, h_yfn, h_ybn;
	float h_zminx, h_zmaxx;  // h_xfp, h_xbp, h_xfn, h_xbn;

	ifstream infile;
    infile.open(secjaws_file);
    if(infile.is_open())//check if file is opened
    {
        while(infile.good()&&!infile.eof())//while file is in good condition and not end of file
        {
            //memset(char_temp,0,100);  //should include <cstring>
            i=0;
            do
            {
                ch=infile.get();
                char_temp[i]=ch;
                i++;
            }
            while(ch!='\n'&&!infile.eof());//while in each line and not end of file
            i=0;
            while(char_temp[i]==CM_title[i]&&i<52)
            {
                i++;
            }
            if(i==52)
            {
                infile>>h_rmax>>ch;   infile.ignore(100,'\n');
                for(i=0;i<3;i++)
                {
                    infile.ignore(100,'\n');
                }
                infile>>h_zminy>>ch>>h_zmaxy>>ch; infile.ignore(100,'\n');
                infile.ignore(100,'\n');
				infile>>h_zminx>>ch>>h_zmaxx>>ch; infile.ignore(100,'\n');
				infile.ignore(100,'\n');
				//for(i=0;i<nmed;i++)
				//{
					infile.ignore(100,'\n');
					int j=0;
					ch=infile.get();
					while(ch!='\n'&&!infile.eof())//while in each line and not end of file
					{
						media_names+=ch;
						ch=infile.get();
						j++;
					}
				//}

                //save variables
				SecJaws.rmax  = h_rmax;
				SecJaws.zminy = h_zminy;
				SecJaws.zmaxy = h_zmaxy;
				//SecJaws.yfp   = h_yfp;
				//SecJaws.ybp   = h_ybp;
				//SecJaws.yfn   = h_yfn;
				//SecJaws.ybn   = h_ybn;
				SecJaws.zminx = h_zminx;
				SecJaws.zmaxx = h_zmaxx;
				//SecJaws.xfp   = h_xfp;
				//SecJaws.xbp   = h_xbp;
				//SecJaws.xfn   = h_xfn;
				//SecJaws.xbn   = h_xbn;

				break;
            }
        }
    }
    infile.close();//close the file

    return;
}

// read egsinp file for Wedge to get geometric data
void read_egsinp_wedge(string &media_names)
{
	
	int i;
    char char_temp[100];
    char ch;
    char CM_title[55]="*********** start of CM PYRAMIDS with identifier WEDGE";  //55 characters
    float h_rmax;
	int h_layers;
	float *coord;
	point2D h_boundWed[4];
	point2D h_goem_Wed[6];

	ifstream infile;
    infile.open(wedge_file);
    if(infile.is_open())//check if file is opened
    {
        while(infile.good()&&!infile.eof())//while file is in good condition and not end of file
        {
            memset(char_temp,0,100);  //should include <cstring>
            i=0;
            do
            {
                ch=infile.get();
                char_temp[i]=ch;
                i++;
            }
            while(ch!='\n'&&!infile.eof());//while in each line and not end of file
            i=0;
            while(char_temp[i]==CM_title[i]&&i<54)
            {
                i++;
            }
            if(i==54)
            {
                infile>>h_rmax>>ch;   infile.ignore(100,'\n');
				infile.ignore(100,'\n');
				infile>>h_layers>>ch;   infile.ignore(100,'\n');
				if(h_layers>4){
					error(9005, "The number of layers %d is larger than 4!!\n", h_layers);
				}

				coord = (float*)malloc(6*h_layers*sizeof(float));
				for(int j=0; j<h_layers; j++){
					for(int k=0; k<6; k++){
						infile>>coord[j*6+k]>>ch;
					}
					infile.ignore(100,'\n');
				}

				infile.ignore(100,'\n');
				//for(i=0;i<nmed;i++)
				//{
					infile.ignore(100,'\n');
					int j=0;
					ch=infile.get();
					while(ch!='\n'&&!infile.eof())//while in each line and not end of file
					{
						media_names+=ch;
						ch=infile.get();
						j++;
					}
				//}
				
				int start,step,offset;
				if(coord[2]-coord[3]==0){
					if(coord[4]-coord[5]>0){ start=4; step=6; offset=4; }
					else{ start=6*h_layers-1; step=-6; offset=4; }
				}
				else if(coord[4]-coord[5]==0){
					if(coord[2]-coord[3]<0){ start=2; step=6; offset=2; }
					else{ start=6*h_layers-3; step=-6; offset=2; }
				}

				h_goem_Wed[0].x = coord[start];
				
				for(j=1;start>=0&&start<6*h_layers;){
					h_goem_Wed[j].x = coord[start];
					h_goem_Wed[j].y = coord[start-offset];
					start+=step;
					j++;
				}
				start=start-step+(step/6);
				h_goem_Wed[j].x = coord[start];
				h_goem_Wed[j].y = coord[start-offset];
				h_goem_Wed[0].y = coord[start-offset];

				char h_nVer_Wed = h_layers + 2;

				h_boundWed[0].x = - h_rmax;
				h_boundWed[0].y = coord[0];
				h_boundWed[1].x =   h_rmax;
				h_boundWed[1].y = coord[0];
				h_boundWed[2].x =   h_rmax;
				h_boundWed[2].y = coord[6*h_layers-5];
				h_boundWed[3].x = - h_rmax;
				h_boundWed[3].y = coord[6*h_layers-5];

				for(int GPUId=0; GPUId<GPUNo; GPUId++) {
#ifdef USE_MULTIPLE_GPU
					cudaSetDevice(GPUId); ce(58064);
#endif
					cudaMemcpyToSymbol(nVer_Wed, &h_nVer_Wed, sizeof(char)); ce(1916);
					cudaMemcpyToSymbol(goem_Wed, h_goem_Wed, h_nVer_Wed*sizeof(point2D)); ce(1917);
					cudaMemcpyToSymbol(boundWed, h_boundWed, 4*sizeof(point2D)); ce(1918);
				}

				free(coord);

				break;
            }
        }
    }
    infile.close();//close the file
	

    return;
}

// read egsinp file for Block to get geometric data
void read_egsinp_block(string &media_names)
{
	
	int i;
    char char_temp[100];
    char ch;
    char CM_title[55]="*********** start of CM PYRAMIDS with identifier BLOCK";  //55 characters
    float h_rmax;
	float *coord;
	point2D h_blockXZ[4];
	point2D h_blockYZ[4];

	ifstream infile;
    infile.open(block_file);
    if(infile.is_open())//check if file is opened
    {
        while(infile.good()&&!infile.eof())//while file is in good condition and not end of file
        {
            memset(char_temp,0,100);  //should include <cstring>
            i=0;
            do
            {
                ch=infile.get();
                char_temp[i]=ch;
                i++;
            }
            while(ch!='\n'&&!infile.eof());//while in each line and not end of file
            i=0;
            while(char_temp[i]==CM_title[i]&&i<54)
            {
                i++;
            }
            if(i==54)
            {
                infile>>h_rmax>>ch;   infile.ignore(100,'\n');
				infile.ignore(100,'\n');
				infile.ignore(100,'\n');

				coord = (float*)malloc(10*sizeof(float));
				for(int j=0; j<10; j++){
					infile>>coord[j]>>ch;
				}
				infile.ignore(100,'\n');

				infile.ignore(100,'\n');
				//for(i=0;i<nmed;i++)
				//{
					infile.ignore(100,'\n');
					int j=0;
					ch=infile.get();
					while(ch!='\n'&&!infile.eof())//while in each line and not end of file
					{
						media_names+=ch;
						ch=infile.get();
						j++;
					}
				//}

				h_blockXZ[0].x = coord[2];
				h_blockXZ[0].y = coord[0];
				h_blockXZ[1].x = coord[3];
				h_blockXZ[1].y = coord[1];
				h_blockXZ[2].x = coord[5];
				h_blockXZ[2].y = coord[1];
				h_blockXZ[3].x = coord[4];
				h_blockXZ[3].y = coord[0];

				h_blockYZ[0].x = coord[6];
				h_blockYZ[0].y = coord[0];
				h_blockYZ[1].x = coord[7];
				h_blockYZ[1].y = coord[1];
				h_blockYZ[2].x = coord[9];
				h_blockYZ[2].y = coord[1];
				h_blockYZ[3].x = coord[8];
				h_blockYZ[3].y = coord[0];

				for(int GPUId=0; GPUId<GPUNo; GPUId++) {
#ifdef USE_MULTIPLE_GPU
					cudaSetDevice(GPUId); ce(58065);
#endif
					cudaMemcpyToSymbol(blockXZ, h_blockXZ, 4*sizeof(point2D)); ce(1914);
					cudaMemcpyToSymbol(blockYZ, h_blockYZ, 4*sizeof(point2D)); ce(1915);
				}

				free(coord);

				break;
            }
        }
    }
    infile.close();//close the file
	

    return;
}

// populate the region data with the values read from the egsphant file and copy it to the device
void init_regions(uint nreg, medium_t **media, int *media_indices, float *densities, int *other_med_IDs) {
    
	//nreg is actually equals the number of Voxels in the phantom
	//total_number_regions= nreg+1+2+2;// "+2" is for the two addition regions of MLC, "+1" is for the VACUUM, "+2" for region of SecJaws
	total_number_regions = nreg + vac_and_oth;// "+2" is for the two addition regions of MLC, "+1" is for the VACUUM, "+2" for region of SecJaws
	uint size = total_number_regions * sizeof(region_data_t);   
    h_region_data = (region_data_t*)malloc(size);
	
    h_region_data[0].med = VACUUM;
    h_region_data[0].flags = 0;
    h_region_data[0].rhof = 0.0F;
    h_region_data[0].pcut = 0.0F;
    h_region_data[0].ecut = 0.0F;
	ushort med;

	//for(int i=1;i<=2;i++)
	for(int i=1;i<=max_oth_reg;i++)   //region 1 and 2 corresponding to the air and W of MLC, their med number is recored in MLC_med_IDs
	{
		med = other_med_IDs[i-1];// should be i-1 , because there is no media for VACUUM by Tong 20130230
		if(i==MLC_region)h_i_do_rr[med]=1;                 //we only do range rejection in W of MLC
		h_e_max_rr[med]=E_MAX_RR;
		h_region_data[i].med = med;
		h_region_data[i].flags = 0x0000U;
#ifdef RAYLE_FLAG
	    h_region_data[i].flags |= f_rayleigh;
#endif
#ifdef RELAX_FLAG
		h_region_data[i].flags |= f_atomic_relaxation;
#endif
#ifdef PHTER_FLAG
	    h_region_data[i].flags |= f_photo_electron_angular_distribution;
#endif
#ifdef BCOMP_FLAG
	    h_region_data[i].flags |= f_bound_compton;
#endif
        h_region_data[i].rhof = 1.0; //media[i-1]->rho;  //note that the indices in h_region_data & media is different 
        h_region_data[i].pcut = media[med]->ap;   //by Wu 20130126
        h_region_data[i].ecut = media[med]->ae;
        //h_region_data[i].rho = media[med]->rho;
	}
	/*
	for(int i=3;i<=4;i++){  //region 3 and 4 corresponding to the air and W of SecJaws
		// i==3 for region of SecJaws
		med = MLC_med_IDs[i-3];// the medium index of SecJaws is the same of MLC
		h_region_data[i].med = med;
		h_region_data[i].flags = 0x0000U;
#ifdef RAYLE_FLAG
		h_region_data[i].flags |= f_rayleigh;
#endif
#ifdef RELAX_FLAG
		h_region_data[i].flags |= f_atomic_relaxation;
#endif
#ifdef PHTER_FLAG
		h_region_data[i].flags |= f_photo_electron_angular_distribution;
#endif
#ifdef BCOMP_FLAG
		h_region_data[i].flags |= f_bound_compton;
#endif
		h_region_data[i].rhof = 1.0; //media[i-1]->rho;  //note that the indices in h_region_data & media is different 
		h_region_data[i].pcut = media[med]->ap;   //by Wu 20130126
		h_region_data[i].ecut = media[med]->ae;
		//h_region_data[3].rho  = media[med]->rho;
	}
	*/
	//for (uint i = 1 + 2 + 2; i < nreg + 1 + 2 + 2; i++) {
    for (uint i = vac_and_oth; i < nreg + vac_and_oth; i++) {
		//int voxelIndex = i - 1 - 2 - 2;
		int voxelIndex = i - vac_and_oth;
        med =  media_indices[voxelIndex] - 1;  // "-1" because in the egsphantom the media number start at 1, but the media array index start at 0
        h_region_data[i].med = med;
        if (med == VACUUM) {
            h_region_data[i].flags = 0;
            h_region_data[i].rhof = 0.0F;
            h_region_data[i].pcut = 0.0F;
            h_region_data[i].ecut = 0.0F;
            //h_region_data[i].rho = 0.0F;
        } 
        else {
			h_region_data[i].flags = 0x0000U;
#ifdef RAYLE_FLAG
			h_region_data[i].flags |= f_rayleigh;
#endif
#ifdef RELAX_FLAG
			h_region_data[i].flags |= f_atomic_relaxation;
#endif
#ifdef PHTER_FLAG
			h_region_data[i].flags |= f_photo_electron_angular_distribution;
#endif
#ifdef BCOMP_FLAG
			h_region_data[i].flags |= f_bound_compton;
#endif
			if (densities[voxelIndex] == 0.0F)
                h_region_data[i].rhof = 1.0F;
            else
                h_region_data[i].rhof = densities[voxelIndex] / (float)media[med]->rho; //density relative to the standard form of the material
                //h_region_data[i].rhof = (float)media[med]->rho; //just use the normal density specified by the pegs4data
            h_region_data[i].pcut = (float)media[med]->ap;
            h_region_data[i].ecut = (float)media[med]->ae;
            //h_region_data[i].rho =  (float) densities[voxelIndex]; //the original material density
        }
    }

	d_region_data = (region_data_t**)malloc(GPUNo*sizeof(region_data_t*));
	d_eng_score = (score_t**)malloc(GPUNo*sizeof(score_t*));
	//score_t *add_e_score;

	for(int GPUId=0; GPUId<GPUNo; GPUId++) {

#ifdef USE_MULTIPLE_GPU
		cudaSetDevice(GPUId); ce(58008);
#endif

		//region_data_t * d_region_data;

		// allocate device memory region data
		cudaMalloc(&d_region_data[GPUId], size); ce(7001);

		// copy region data to device
		cudaMemcpy(d_region_data[GPUId], h_region_data, size, cudaMemcpyHostToDevice); ce(7002);
		cudaMemcpyToSymbol(region_data, &d_region_data[GPUId], sizeof(region_data_t*)); ce(7003);
		//cudaMemcpyToSymbol(region_data, h_region_data, (nreg+1 + 2)*sizeof(region_data_t)); ce(7003);
	
		//now allocate momory to store engery deposition
		cudaMalloc(&d_eng_score[GPUId], nreg*sizeof(score_t)); ce(7004);
		//cudaMemset(d_eng_score, 0, nreg*sizeof(score_t)); ce(7005);
		cudaMemcpyToSymbol(eng_score, &d_eng_score[GPUId], sizeof(score_t*)); ce(7006);
		//cudaGetSymbolAddress((void**)&add_e_score,eng_score); ce(7007);
		//cudaMemcpy(add_e_score, &d_eng_score[GPUId], sizeof(score_t*),cudaMemcpyHostToDevice); ce(7008);

	}

#ifdef GET_UNCERTAINTY
	//allocate memory on CPU to store sum and square of energy deposition
	eng_dep = (score_t*)malloc(nreg*sizeof(score_t));
	memset(eng_dep, 0, nreg*sizeof(score_t));
	eng_dep2 = (score_t*)malloc(nreg*sizeof(score_t));
	memset(eng_dep2, 0, nreg*sizeof(score_t));
#endif

    // free host memory
    //free(h_region_data);   // we not going to free it right now, because the some of the Physics Initialziation still need the data

}

/*
void init_MLC_region(int nreg, medium_t **media){

	    region_data_t  h_region_data[MAXREG] ;

		if( MAXREG < (nreg+1)) 
	{
		printf (" number of regions required is larger than the MAXREG allocated, please adjust MAXREG. \n STOP here");
		free_all();
		exit( 1);
	}

    //uint size = (nreg + 1) * sizeof(region_data_t);
    //region_data_t *h_region_data = (region_data_t*)malloc(size);


    // region 0 is outside of simulation geometry
    h_region_data[0].med = VACUUM;
    h_region_data[0].flags = 0;
    h_region_data[0].rhof = 0.0F;
    h_region_data[0].pcut = 0.0F;
    h_region_data[0].ecut = 0.0F;
	for(int i=1;i<=nreg;i++)
	{
		h_region_data[i].med = i-1; // med should be i-1 , by Tong 20130230
	    h_region_data[i].flags = f_rayleigh;
        h_region_data[i].rhof = media[i-1]->rho;  //note that the indices in h_region_data & media is different 
        h_region_data[i].pcut = media[i-1]->ap;   //by Wu 20130126
        h_region_data[i].ecut = media[i-1]->ae;
	}

	//cudaMalloc(&d_region_data, size); ce(7001);
    // copy region data to device
    //cudaMemcpy(d_region_data, h_region_data, size, cudaMemcpyHostToDevice); ce(7002);
	//cudaMemcpyToSymbol(region_data, &d_region_data, sizeof(region_data_t*)); ce(7003);

	cudaMemcpyToSymbol(region_data, h_region_data, (nreg+1)*sizeof(region_data_t)); ce(7003);

    // free host memory
    //free(h_region_data);
}
*/

/*
void init_phantom_texture( medium_t **media, int *media_indices, float *densities)
{
	int nx, ny, nz;
	nx=h_phantom_N.x;
	ny=h_phantom_N.y;
	nz=h_phantom_N.z;

	med_flags_h = (int2*) malloc(nz*ny*nx*sizeof(int2));
	rhof_rho_pcut_ecut_h = (float4*)malloc(nz*ny*nx*sizeof(float4));

	for (int z = 0; z < nz; z++) {
        for (int y = 0; y < ny; y++) {
            for (int x = 0; x < nx; x++) {
				int voxelIndex = x + y*nx + z*nx*ny;
				int med =media_indices[voxelIndex] -1;  
				med_flags_h[voxelIndex].x=med;
				med_flags_h[voxelIndex].y=f_rayleigh;
				rhof_rho_pcut_ecut_h[voxelIndex].x = densities[voxelIndex] / (float)media[med]->rho;
				rhof_rho_pcut_ecut_h[voxelIndex].y = densities[voxelIndex] ;
				rhof_rho_pcut_ecut_h[voxelIndex].z = media[med]->ap ;
				rhof_rho_pcut_ecut_h[voxelIndex].w = media[med]->ae;
			}
		}
	}
				

	// in the original CUDAEGS, the phantom voxel index is cacluated by ix + iy*Nx + iz *Nx*Ny.
	// this conversion actually corresponding to the same memory arrangement of a 3D array  voxel[Nz][Ny][Nx]
	// however, the CUDA array has the reverse memory arrangement on the dimensions 
	// so to be consistant, we define the CUDA array following the same order of dimesions.
	volumeSize = make_cudaExtent(nx,ny,nz);

	cudaMalloc3DArray(&med_flags, &med_flags_tex.channelDesc, volumeSize) ;
	cudaMalloc3DArray(&rhof_rho_pcut_ecut, &rhof_rho_pcut_ecut_tex.channelDesc, volumeSize);
	// create a 3d array on device

	copyParams.srcPtr   = make_cudaPitchedPtr((void*)med_flags_h, volumeSize.width*sizeof(int2), volumeSize.width, volumeSize.height);
	copyParams.dstArray = med_flags;
	copyParams.extent   = volumeSize;
	copyParams.kind     = cudaMemcpyHostToDevice;
	cudaMemcpy3D(&copyParams) ;
	// copy data from host to device
	
	med_flags_tex.normalized = false;
	med_flags_tex.filterMode = cudaFilterModePoint;
	cudaBindTextureToArray(med_flags_tex, med_flags, med_flags_tex.channelDesc);
	// bind to texture memory

	copyParams.srcPtr   = make_cudaPitchedPtr((void*)rhof_rho_pcut_ecut_h, volumeSize.width*sizeof(float4), volumeSize.width, volumeSize.height);
	copyParams.dstArray = rhof_rho_pcut_ecut;
	copyParams.extent   = volumeSize;
	copyParams.kind     = cudaMemcpyHostToDevice;
	cudaMemcpy3D(&copyParams) ;
	// copy data from host to device
	
	med_flags_tex.normalized = false;
	med_flags_tex.filterMode = cudaFilterModePoint;
	cudaBindTextureToArray(rhof_rho_pcut_ecut_tex, rhof_rho_pcut_ecut, rhof_rho_pcut_ecut_tex.channelDesc);
	// bind to texture memory
	
	free(med_flags_h);
	free(rhof_rho_pcut_ecut_h);

}
*/

// read the egsphant file and copy the phantom data to the device
void init_phantom_and_MLC(string default_medium) {

    // read egsphant file
    FILE *f;
    if (fopen_s(&f, egsphant_file, "r"))
		error(8001, "Could not open the phantom file \"%s\".\n", egsphant_file);

    int nmed = 0; // number of media
    if (fscanf(f, "%d\n", &nmed) != 1)
		error(8002, "Could not read the number of media form the phantom file \"%s\".\n", egsphant_file);

    // read media names
    string *media_names_ALL = new string[nmed+5];  // " nmed +2 " , since there will be 2 media for the MLC,
	//string *media_names = media_names_ALL +2;         //reserve the first two  Media for the MLC, this way, the following code can stay unchanged.
	string *media_names = media_names_ALL;         //reserve the first two  Media for the MLC, this way, the following code can stay unchanged.
    for (int i = 0; i < nmed; i++) {
        char mBuf[1024];
        if (!fgets(mBuf, 1024, f))
            error(8003, "Could not read medium name (i = %d) in the phantom file \"%s\".\n", i , egsphant_file);
        media_names[i] = trim(string(mBuf), false);  //changed from "allow space" = true to false. Tong Xu Feb 2013.
        if (media_names[i].compare("") == 0)
            error(8004, "Medium name cannot be empty (i = %d) in the phantom file \"%s\".\n", i , egsphant_file);
    }

    // skip next line (contains dummy input)
    char dummyBuf[1024];
    fgets(dummyBuf, 1024, f);

    // read voxel numbers
    int nx, ny, nz;
    if (fscanf(f, "%d %d %d\n", &nx, &ny, &nz) != 3)
        error(8006, "Could not read the voxel numbers in the phantom file \"%s\".\n", egsphant_file);

    // read voxel boundaries
    //float *x_bounds = (float*)malloc((nx + 1) * sizeof(float));
	x_bounds = (float*)malloc((nx + 1) * sizeof(float));
    for (int i = 0; i <= nx; i++) {
        if (fscanf(f, "%f", x_bounds + i) != 1)
            error(8007, "Could not read x coordinate of voxel boundary (i = %d) in the phantom file \"%s\".\n", i , egsphant_file);
    }
    //float *y_bounds = (float*)malloc((ny + 1) * sizeof(float));
	y_bounds = (float*)malloc((ny + 1) * sizeof(float));
    for (int i = 0; i <= ny; i++) {
        if (fscanf(f, "%f", y_bounds + i) != 1)
            error(8008, "Could not read y coordinate of voxel boundary (i = %d) in the phantom file \"%s\".\n", i , egsphant_file);
    }
    //float *z_bounds = (float*)malloc((nz + 1) * sizeof(float));
	z_bounds = (float*)malloc((nz + 1) * sizeof(float));
    for (int i = 0; i <= nz; i++) {
        if (fscanf(f, "%f", z_bounds + i) != 1)
            error(8009, "Could not read z coordinate of voxel boundary (i = %d) in the phantom file \"%s\".\n", i , egsphant_file);
    }

    // read media indices of voxels in Phantom
    int *media_indices = (int*)malloc(nx * ny * nz * sizeof(int));
    for (int z = 0; z < nz; z++) {
        for (int y = 0; y < ny; y++) {
            for (int x = 0; x < nx; x++) {
                if (fscanf(f, "%1d", media_indices + x + y * nx + z * nx * ny) != 1)
                    error(8010, "Could not read media index (x = %d, y = %d, z = %d) in the phantom file \"%s\".\n", x, y, z , egsphant_file);
            }
        }

        // skip blank line
        char dummyBuf[1024];
        fgets(dummyBuf, 1024, f);
        if (trim(string(dummyBuf), false).compare("") != 0)
            error(8011, "Expected empty line but found \"%s\" (z = %d) in the phantom file \"%s\".\n", dummyBuf, z, egsphant_file);
    }
    
    // read densities
    float *densities = (float*)malloc(nx * ny * nz * sizeof(float));
    for (int z = 0; z < nz; z++) {
        for (int y = 0; y < ny; y++) {
            for (int x = 0; x < nx; x++) {
                if (fscanf(f, "%f", densities + x + y * nx + z * nx * ny) != 1)
                    error(8012, "Could not read density (x = %d, y = %d, z = %d) in the phantom file \"%s\".\n", x, y, z , egsphant_file);
            }
        }

        // skip blank line
        char dummyBuf[1024];
        fgets(dummyBuf, 1024, f);
        if (trim(string(dummyBuf), false).compare("") != 0)
            error(8013, "Expected empty line but found \"%s\" (z = %d) in the phantom file \"%s\".\n", dummyBuf, z, egsphant_file);
    }

    fclose(f);


	//now calculate the mass of each voxel, will be used to convert energy (KeV) to dose (Gy) when write_dose().
	//we don't need to copy voxel 
	voxel_mass = (float*)malloc(nx * ny * nz * sizeof(float));
	float voxel_volum = 0;
	float mass;
    for (int z = 0; z < nz; z++) {
        for (int y = 0; y < ny; y++) {
            for (int x = 0; x < nx; x++) {
				voxel_volum = (x_bounds[x+1]-x_bounds[x])*(y_bounds[y+1]-y_bounds[y])*(z_bounds[z+1]-z_bounds[z]);
				if(voxel_volum<0) voxel_volum = - voxel_volum;  //incase the bounds are not alway increaseing with index.
				mass = voxel_volum*densities[x + y * nx + z * nx * ny];
				//if(densities[x + nx*y + nx*ny*z]<0.02) mass = 1e9;  //if density is very low, it is air, we are not interested
				                                                    //effectively make the dose very small to avoid display issues
				if(mass==0) mass = 1;  //to avoid divided by zero
				voxel_mass[x + y * nx + z * nx * ny] = mass;
			}
		}
	}


    // copy phantom data to device
    h_phantom_N.x = nx;
    h_phantom_N.y = ny;
    h_phantom_N.z = nz;

	h_size_phantom = h_phantom_N.x*h_phantom_N.y*h_phantom_N.z;
	
	d_phantom_x_bounds = (float**)malloc(GPUNo*sizeof(float*));
	d_phantom_y_bounds = (float**)malloc(GPUNo*sizeof(float*));
	d_phantom_z_bounds = (float**)malloc(GPUNo*sizeof(float*));

	for(int GPUId=0; GPUId<GPUNo; GPUId++) {

#ifdef USE_MULTIPLE_GPU
		cudaSetDevice(GPUId); ce(58003);
#endif

		cudaMalloc(&d_phantom_x_bounds[GPUId], (nx + 1) * sizeof(float)); ce(8014);
		cudaMalloc(&d_phantom_y_bounds[GPUId], (ny + 1) * sizeof(float)); ce(8015);
		cudaMalloc(&d_phantom_z_bounds[GPUId], (nz + 1) * sizeof(float)); ce(8016);

		cudaMemcpy(d_phantom_x_bounds[GPUId], x_bounds, (nx + 1) * sizeof(float), cudaMemcpyHostToDevice); ce(8017);
		cudaMemcpy(d_phantom_y_bounds[GPUId], y_bounds, (ny + 1) * sizeof(float), cudaMemcpyHostToDevice); ce(8018);
		cudaMemcpy(d_phantom_z_bounds[GPUId], z_bounds, (nz + 1) * sizeof(float), cudaMemcpyHostToDevice); ce(8019);

		//cudaMemcpyToSymbol(phantom, &h_phantom, sizeof(phantom_t)); ce(8020);
		cudaMemcpyToSymbol(phantom_N, &h_phantom_N, sizeof(uint3)); ce(80201);
		cudaMemcpyToSymbol(phantom_x_bounds, &d_phantom_x_bounds[GPUId], sizeof(float*)); ce(80202);
		cudaMemcpyToSymbol(phantom_y_bounds, &d_phantom_y_bounds[GPUId], sizeof(float*)); ce(80203);
		cudaMemcpyToSymbol(phantom_z_bounds, &d_phantom_z_bounds[GPUId], sizeof(float*)); ce(80204);

		cudaMemcpyToSymbol(size_phantom, &h_size_phantom, sizeof(uint)); ce(8021);

	}

    for (int x = 0; x < nx; x++)
		x_bounds[x] = (x_bounds[x]+x_bounds[x+1])/2;
    for (int y = 0; y < ny; y++)
		y_bounds[y] = (y_bounds[y]+y_bounds[y+1])/2;
    for (int z = 0; z < nz; z++)
		z_bounds[z] = (z_bounds[z]+z_bounds[z+1])/2;

    //free(x_bounds);
    //free(y_bounds);
    //free(z_bounds);

	//h_size_phantom = h_phantom_N.x*h_phantom_N.y*h_phantom_N.z;
    //cudaMemcpyToSymbol(size_phantom, &h_size_phantom, sizeof(uint)); ce(8021);

	string media_names_other[max_oth_med];  //2 for MLC, 0th is air in MLC, 1st is W in MLC
	media_names_other[default_med]=default_medium;
	read_egsinp_secjaws(media_names_other[secjaws_med]);
	read_egsinp(media_names_other[MLC_medium]); //this will read the egsinp file and the two medias names of the MLC (air and W) 
	if(Block_on) read_egsinp_block(media_names_other[block_med]);
	else media_names_other[block_med]=default_medium;
	if(Wedge_on) read_egsinp_wedge(media_names_other[wedge_med]);
	else media_names_other[wedge_med]=default_medium;

	//now check if any of the media_names_MLC already exist in and media_names from phantom
	int other_med_IDs[max_oth_med];
	int k=nmed;
	for(int i=0; i<max_oth_med; i++)
	{
		bool exist=false;
		for(int j=0; j<k;j++)
		{
			if(media_names_other[i]==media_names_ALL[j])
			{
				other_med_IDs[i]=j;  
				exist = true;
			}
		}
		if(!exist) //IF NOT FOUND, ADD THE NAME TO THE FINAL NAME LIST
		{
			media_names_ALL[k]=media_names_other[i];
			other_med_IDs[i] = k;   //assign the med number for this MLC material.
			k++;
		}
	}
	int nmed_all = k;
	h_nmed=nmed_all;

	// init media
	medium_t **media = init_media(nmed_all, media_names_ALL);  // "+2" is for the First 2 media that used by MLC
	for(int i=0; i<h_nmed;i++)  //range rejection medium by medium
	{
		h_i_do_rr[i]=0;
		h_e_max_rr[i]=E_MAX_RR;
	}

    // init phantom regions (now, the MLC regions are also setup in "init_regions".
    init_regions(nx * ny * nz, media, media_indices, densities, other_med_IDs);
	//init_phantom_texture(media, media_indices, densities);

	//write_density("output\\Debug", densities);
	//write_media_indices("output\\Debug", media_indices);

    free(media);
    free(media_indices);
    free(densities);

}

void read_mlcinp()
{
	int i,j,temp;
    char ch;

	float *field_size;
	float *out_factor;
	int Num_field;
	float actual_field_size;

	ifstream infile;
	
	if(outfac_on){
		infile.open(output_factor_file);
		if(infile.is_open())//check if file is opened
		{
			//while(infile.good()&&!infile.eof())//while file is in good condition and not end of file
			if(infile.good())
			{
				infile.ignore(100,'\n');
				infile>>Num_field;        infile.ignore(100,'\n');//ignore the '\n' of this line
				field_size=(float*)malloc(Num_field*sizeof(float));
				out_factor=(float*)malloc(Num_field*sizeof(float));
				for(i=0;i<Num_field&&!infile.eof();i++)
				{
					infile>>field_size[i]>>out_factor[i];
					infile.ignore(100,'\n');//ignore the '\n' of this line
					//field_size[i]*=field_size[i];
				}
			}
		}
		else
		{ 
			logstr("Error: output factor file reading error!\n");
		}
		infile.close();//close the file
	}

	total_MU = 0;

    infile.open(mlcinfo_file);
    if(infile.is_open())//check if file is opened
    {
		//while(infile.good()&&!infile.eof())//while file is in good condition and not end of file
		if(infile.good())
		{
			infile.ignore(100,'\n');
			infile>>N_Field;        infile.ignore(100,'\n');//ignore the '/n' of this line
			if(N_Field >MAXFIELDS){
				N_Field=MAXFIELDS;
				printf("Warning! number of MLC fields is larger than the MAXFIELDS=%d\n",MAXFIELDS);
			}
			for(i=0;i<N_Field&&!infile.eof();i++)
			{
				infile>>temp>>ch>>field[i].angle>>ch>>field[i].x1>>ch>>field[i].x2>>ch>>field[i].y1>>ch>>field[i].y2>>ch>>field[i].MU>>ch;
				//field[i].angle = field[i].angle *3.1415926/180;
				total_MU += field[i].MU;
				string mlcinp_file_input;
                j=0;
				ch=infile.get();
				while(ch!='\n'&&!infile.eof())//while in each line and not end of file
                {
					mlcinp_file_input+=ch;
					ch=infile.get();
                    j++;
                }
				mlcinp_file_input=trim(mlcinp_file_input,false);
				//GetFullPathName(mlcinp_file_input.c_str(), CHARLEN, charBuffer, NULL);
				//mlcinp_file_input = string(charBuffer);
				mlcinp_file_input = string(data_dir) + mlcinp_file_input;
				field[i].mlcinp_file = mlcinp_file_input;

				if(outfac_on){
					actual_field_size = sqrt((field[i].x2 - field[i].x1) * (field[i].y2 - field[i].y1));
					for(j=0;j<Num_field-1;j++){
						if(actual_field_size>field_size[j]&&actual_field_size<=field_size[j+1]){
							field[i].output_factor = out_factor[j] + (actual_field_size - field_size[j])
								* (out_factor[j+1] - out_factor[j]) / (field_size[j+1] - field_size[j]);
							break;
						}
					}
					/*
					if(j==Num_field-1){
						field[i].output_factor = 1.0;
						logstr("Warning: size of field #%d is not between %d and %d\n",i,(int)sqrt(field_size[0]),(int)sqrt(field_size[Num_field-1]));
					}
					*/
					//logstr("output factor of field #%d is %.5f\n",i,field[i].output_factor);
				}
			}
			N_Field =i;
		}
	}
	else
	{ 
		logstr("Error: MLC info  file reading error! No simulation will be done!\n");
		N_Field= 0;
	}
    infile.close();//close the file

	for(i=0;i<N_Field;i++)
	{
		string mlcinp_file=field[i].mlcinp_file;
		infile.open(mlcinp_file.c_str());
		if(infile.is_open())//check if file is opened
		{
            for(j=0;j<14;j++)
				infile.ignore(100,'\n');    //ignore 14 lines
			for(j=0;j<120;j++)
			{
				infile.ignore(10,'\n');     //ignore 10 characters before the number we need,  should be carefullly check!!!!
				infile>>field[i].pos[j];
				infile.ignore(100,'\n');    //ignore the '\n' of this line
				//field[i].pos[j] =field[i].pos[j] * (VarianMLC.zmin + VarianMLC.zthick/2.0F)/SourceToIsocenter;
			}
		}
		infile.close();
	}

	if(outfac_on){
		free(field_size);
		free(out_factor);
	}

    return;
}

void setup_secjaws(uint iField){
	point2D h_boundMinY[4],h_boundMinX[4];
	point2D h_boundMaxY[4],h_boundMaxX[4];
	point2D h_jawsGeomY[8],h_jawsGeomX[8];

	h_boundMinY[0].x = field[iField].x1*SecJaws.zminy/h_SourceToIsocenter - 0.5;
	h_boundMinY[0].y = field[iField].y1*SecJaws.zminy/h_SourceToIsocenter - 0.5;
	h_boundMinY[1].x = field[iField].x2*SecJaws.zminy/h_SourceToIsocenter + 0.5;
	h_boundMinY[1].y = field[iField].y1*SecJaws.zminy/h_SourceToIsocenter - 0.5;
	h_boundMinY[2].x = field[iField].x2*SecJaws.zminy/h_SourceToIsocenter + 0.5;
	h_boundMinY[2].y = field[iField].y2*SecJaws.zminy/h_SourceToIsocenter + 0.5;
	h_boundMinY[3].x = field[iField].x1*SecJaws.zminy/h_SourceToIsocenter - 0.5;
	h_boundMinY[3].y = field[iField].y2*SecJaws.zminy/h_SourceToIsocenter + 0.5;

	h_boundMinX[0].x = field[iField].x1*SecJaws.zminx/h_SourceToIsocenter - 0.5;
	h_boundMinX[0].y = field[iField].y1*SecJaws.zminx/h_SourceToIsocenter - 0.5;
	h_boundMinX[1].x = field[iField].x2*SecJaws.zminx/h_SourceToIsocenter + 0.5;
	h_boundMinX[1].y = field[iField].y1*SecJaws.zminx/h_SourceToIsocenter - 0.5;
	h_boundMinX[2].x = field[iField].x2*SecJaws.zminx/h_SourceToIsocenter + 0.5;
	h_boundMinX[2].y = field[iField].y2*SecJaws.zminx/h_SourceToIsocenter + 0.5;
	h_boundMinX[3].x = field[iField].x1*SecJaws.zminx/h_SourceToIsocenter - 0.5;
	h_boundMinX[3].y = field[iField].y2*SecJaws.zminx/h_SourceToIsocenter + 0.5;

	h_boundMaxY[0].x = - SecJaws.rmax;
	h_boundMaxY[0].y =   SecJaws.zminy;
	h_boundMaxY[1].x =   SecJaws.rmax;
	h_boundMaxY[1].y =   SecJaws.zminy;
	h_boundMaxY[2].x =   SecJaws.rmax;
	h_boundMaxY[2].y =   SecJaws.zmaxy;
	h_boundMaxY[3].x = - SecJaws.rmax;
	h_boundMaxY[3].y =   SecJaws.zmaxy;

	h_boundMaxX[0].x = - SecJaws.rmax;
	h_boundMaxX[0].y =   SecJaws.zminx;
	h_boundMaxX[1].x =   SecJaws.rmax;
	h_boundMaxX[1].y =   SecJaws.zminx;
	h_boundMaxX[2].x =   SecJaws.rmax;
	h_boundMaxX[2].y =   SecJaws.zmaxx;
	h_boundMaxX[3].x = - SecJaws.rmax;
	h_boundMaxX[3].y =   SecJaws.zmaxx;

	h_jawsGeomY[0].x = - SecJaws.rmax;
	h_jawsGeomY[0].y =   SecJaws.zmaxy;
	h_jawsGeomY[1].x = - SecJaws.rmax;
	h_jawsGeomY[1].y =   SecJaws.zminy;
	h_jawsGeomY[2].x =   field[iField].y1*SecJaws.zminy/h_SourceToIsocenter;
	h_jawsGeomY[2].y =   SecJaws.zminy;
	h_jawsGeomY[3].x =   field[iField].y1*SecJaws.zmaxy/h_SourceToIsocenter;
	h_jawsGeomY[3].y =   SecJaws.zmaxy;
	h_jawsGeomY[4].x =   field[iField].y2*SecJaws.zmaxy/h_SourceToIsocenter;
	h_jawsGeomY[4].y =   SecJaws.zmaxy;
	h_jawsGeomY[5].x =   field[iField].y2*SecJaws.zminy/h_SourceToIsocenter;
	h_jawsGeomY[5].y =   SecJaws.zminy;
	h_jawsGeomY[6].x =   SecJaws.rmax;
	h_jawsGeomY[6].y =   SecJaws.zminy;
	h_jawsGeomY[7].x =   SecJaws.rmax;
	h_jawsGeomY[7].y =   SecJaws.zmaxy;

	h_jawsGeomX[0].x = - SecJaws.rmax;
	h_jawsGeomX[0].y =   SecJaws.zmaxx;
	h_jawsGeomX[1].x = - SecJaws.rmax;
	h_jawsGeomX[1].y =   SecJaws.zminx;
	h_jawsGeomX[2].x =   field[iField].x1*SecJaws.zminx/h_SourceToIsocenter;
	h_jawsGeomX[2].y =   SecJaws.zminx;
	h_jawsGeomX[3].x =   field[iField].x1*SecJaws.zmaxx/h_SourceToIsocenter;
	h_jawsGeomX[3].y =   SecJaws.zmaxx;
	h_jawsGeomX[4].x =   field[iField].x2*SecJaws.zmaxx/h_SourceToIsocenter;
	h_jawsGeomX[4].y =   SecJaws.zmaxx;
	h_jawsGeomX[5].x =   field[iField].x2*SecJaws.zminx/h_SourceToIsocenter;
	h_jawsGeomX[5].y =   SecJaws.zminx;
	h_jawsGeomX[6].x =   SecJaws.rmax;
	h_jawsGeomX[6].y =   SecJaws.zminx;
	h_jawsGeomX[7].x =   SecJaws.rmax;
	h_jawsGeomX[7].y =   SecJaws.zmaxx;

	for(int GPUId=0; GPUId<GPUNo; GPUId++) {

#ifdef USE_MULTIPLE_GPU
		cudaSetDevice(GPUId); ce(58054);
#endif
		cudaMemcpyToSymbol(boundMinY, h_boundMinY, 4*sizeof(point2D));
		cudaMemcpyToSymbol(boundMinX, h_boundMinX, 4*sizeof(point2D));
		cudaMemcpyToSymbol(boundMaxY, h_boundMaxY, 4*sizeof(point2D));
		cudaMemcpyToSymbol(boundMaxX, h_boundMaxX, 4*sizeof(point2D));
		cudaMemcpyToSymbol(jawsGeomY, h_jawsGeomY, 8*sizeof(point2D));
		cudaMemcpyToSymbol(jawsGeomX, h_jawsGeomX, 8*sizeof(point2D));
	}

}

void init_leaf_pos(uint iField)
{
	float zmin  = VarianMLC.zmin;
	float zthick = VarianMLC.zthick;
	float zcenter = zmin+zthick/2;

	d_endnode = (point2D**)malloc(GPUNo*TOTAL_MLC_LEAF/2*sizeof(point2D*));

	for(uchar j=0;j<TOTAL_MLC_LEAF/2;j++){  //corrected by Wu, 20130727
		//if(field[iField].pos[j]>0||field[iField].pos[j+TOTAL_MLC_LEAF/2]>0){
		if(field[iField].pos[j]+field[iField].pos[j+TOTAL_MLC_LEAF/2]>0){
			//float NEG=crossXZ[14+index*16].x;    //leaf position from .egsinp file
			//float POS=crossXZ[14+index*16].y;
#ifndef INVERSE_MLC_Bank_AB
			float NEG = - field[iField].pos[j];    //leaf position from .mlcinfo file
			float POS = field[iField].pos[j+TOTAL_MLC_LEAF/2];
#else
			float POS = field[iField].pos[j];    //leaf position from .mlcinfo file
			float NEG = - field[iField].pos[j+TOTAL_MLC_LEAF/2];
#endif
			char i;
			point2D h_endnode[34];
			h_endnode[0].x=h_crossYZ[0].x;
			h_endnode[0].y=h_crossYZ[0].y;
			h_endnode[1].x=h_crossYZ[1].x;
			h_endnode[1].y=h_crossYZ[1].y;
			for(i=2;i<17;i++)
			{
				//h_endnode[i].x=h_crossYZ[i].x+NEG;
				h_endnode[i].x=h_crossYZ[i].x+NEG*zcenter/h_SourceToIsocenter;
				h_endnode[i].y=h_crossYZ[i].y;
			}
			for(i=17;i<32;i++)
			{
				//h_endnode[i].x=h_crossYZ[i].x+POS;
				h_endnode[i].x=h_crossYZ[i].x+POS*zcenter/h_SourceToIsocenter;
				h_endnode[i].y=h_crossYZ[i].y;
			}
			h_endnode[32].x=h_crossYZ[32].x;
			h_endnode[32].y=h_crossYZ[32].y;
			h_endnode[33].x=h_crossYZ[33].x;
			h_endnode[33].y=h_crossYZ[33].y;

			for(int GPUId=0; GPUId<GPUNo; GPUId++) {
#ifdef USE_MULTIPLE_GPU
				cudaSetDevice(GPUId); ce(58011);
#endif
				cudaMalloc(&d_endnode[j+GPUId*TOTAL_MLC_LEAF/2],34*sizeof(point2D)); ce(20011);
				cudaMemcpy(d_endnode[j+GPUId*TOTAL_MLC_LEAF/2],h_endnode,34*sizeof(point2D),cudaMemcpyHostToDevice); ce(20012);
			}
		}
		else
		{
			for(int GPUId=0; GPUId<GPUNo; GPUId++) {
#ifdef USE_MULTIPLE_GPU
				cudaSetDevice(GPUId); ce(58055);
#endif
				d_endnode[j+GPUId*TOTAL_MLC_LEAF/2]=d_crossYZ[GPUId];
			}
		}
	}

	for(int GPUId=0; GPUId<GPUNo; GPUId++) {
#ifdef USE_MULTIPLE_GPU
		cudaSetDevice(GPUId); ce(58012);
#endif
		cudaMemcpyToSymbol(endnode,&d_endnode[GPUId*TOTAL_MLC_LEAF/2],TOTAL_MLC_LEAF/2*sizeof(point2D*)); ce(20013);
	}
}

void free_leaf_pos(uint iField)
{
	for(uchar i=0;i<TOTAL_MLC_LEAF/2;i++){
		//if(field[iField].pos[i]>0||field[iField].pos[i+TOTAL_MLC_LEAF/2]>0){
		if(field[iField].pos[i]+field[iField].pos[i+TOTAL_MLC_LEAF/2]>0){
			for(int GPUId=0; GPUId<GPUNo; GPUId++) {
#ifdef USE_MULTIPLE_GPU
				cudaSetDevice(GPUId); ce(58013);
#endif
				cudaFree(d_endnode[i+GPUId*TOTAL_MLC_LEAF/2]); ce(20013);
			}
		}
	}

	free(d_endnode);
}

// call the above functions to perfrom the initialization
void init(string default_medium) {

	launch_setup_kernel();
    // init_MTs(seed);  //should remove this line when cuRand works.
	egs_set_defaults();
    //init_stack();
    //init_source();
    //init_detector();
    init_phantom_and_MLC(default_medium);  // MLC initialization is now part of init_phantom
	mscati();
	edgset(1,1);
	init_compton();
	fix_brems();
	eii_init();
	//dump_phys();
	//copyToDevice();
	read_mlcinp();

	for(int GPUId=0; GPUId<GPUNo; GPUId++) {

#ifdef USE_MULTIPLE_GPU
		cudaSetDevice(GPUId); ce(58009);
#endif

		copyToDevice();

	}

}

// free all allocated memory on the host and the device
void free_all() {

    // free host memory
#ifdef GET_UNCERTAINTY
	free(eng_dep);
	free(eng_dep2);
#endif

	free(npphsp);
    free(x_bounds);
    free(y_bounds);
    free(z_bounds);
    free(h_region_data);
	//free(h_total_step_counts);

	for(int GPUId=0; GPUId<GPUNo; GPUId++) {

#ifdef USE_MULTIPLE_GPU
		cudaSetDevice(GPUId); ce(58014);
#endif

		/*
		cudaUnbindTexture(rhof_rho_pcut_ecut_tex);ce(9101);
		cudaFreeArray(rhof_rho_pcut_ecut);ce(9102);

		cudaUnbindTexture(med_flags_tex);
		cudaFreeArray(med_flags);
		*/

		/*
		for(uchar i=0;i<NUM_SIMU_CAT;i++){
			cudaFree(d_stack[i]); ce(9051);
			//cudaFree(d_stack[i].a); ce(9051);
			//cudaFree(d_stack[i].b); ce(9052);
			//cudaFree(d_stack[i].c); ce(9053);
		}
		*/
		
		cudaFree(d_devStates[GPUId]);ce(9103);
		
		cudaFree(d_crossYZ[GPUId]); ce(9031);

		//cudaFree(d_total_step_counts[GPUId]); ce(9007);

		cudaFree(d_phantom_x_bounds[GPUId]); ce(9016);
		cudaFree(d_phantom_y_bounds[GPUId]); ce(9017);
		cudaFree(d_phantom_z_bounds[GPUId]); ce(9018);

		cudaFree(d_region_data[GPUId]); ce(9019);

		cudaFree(d_eng_score[GPUId]);ce(9030);

		cudaFree(d_ge[GPUId]); ce(9020);
		cudaFree(d_gmfp[GPUId]); ce(9021);
		cudaFree(d_gbr1[GPUId]); ce(9022);
		cudaFree(d_gbr2[GPUId]); ce(9023);
		cudaFree(d_cohe[GPUId]); ce(9024);
		
		cudaFree(d_rayleigh_data[GPUId]); ce(9025);
		cudaFree(d_i_array[GPUId]); ce(9026);
		cudaFree(d_pmax[GPUId]); ce(9027);

	}

	free(d_devStates);
	free(d_crossYZ);
	//free(d_total_step_counts);
	free(d_phantom_x_bounds);
	free(d_phantom_y_bounds);
	free(d_phantom_z_bounds);
	free(d_region_data);
	free(d_eng_score);

	free(d_ge);
	free(d_gmfp);
	free(d_gbr1);
	free(d_gbr2);
	free(d_cohe);

	free(d_rayleigh_data);
	free(d_i_array);
	free(d_pmax);
}


#endif
