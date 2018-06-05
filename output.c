/****************************************************************************
 *
 * output.c, Version 1.0.0 Mon 09 Jan 2012
 *
 * ----------------------------------------------------------------------------
 *
 * CUDA EGS
 * Copyright (C) 2012 CancerCare Manitoba
 *
 * The latest version of CUDA EGS and additional information are available online at 
 * http://www.physics.umanitoba.ca/~elbakri/cuda_egs/ and http://www.lippuner.ca/cuda_egs
 *
 * CUDA EGS is free software; you can redistribute it and/or modify it under the 
 * terms of the GNU General Public License as published by the Free Software 
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.                                       
 *                                                                           
 * CUDA EGS is distributed in the hope that it will be useful, but WITHOUT ANY 
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
 * details.                              
 *                                                                           
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * ----------------------------------------------------------------------------
 *
 *   Contact:
 *
 *   Jonas Lippuner
 *   Email: jonas@lippuner.ca 
 *
 ****************************************************************************/

#ifdef CUDA_EGS

#include "bmp.c"

// write a formatted string with arguments to the log file
void logstr(const char *format, ...) {
    va_list arguments;
    va_start(arguments, format);
    vprintf(format, arguments);
    vfprintf(logfile, format, arguments);
    va_end(arguments);
    fflush(logfile);
}

// write error message to log and exit
void error(int ret, const char *format, ...) {
    va_list arguments;
    va_start(arguments, format);
    printf("\nERROR (%d): ", ret);
	fprintf(logfile, "\nERROR (%d): ", ret);
	vprintf(format, arguments);
    vfprintf(logfile, format, arguments);
    va_end(arguments);
    fflush(logfile);
    exit(ret);
}

// handle CUDA errors
void ce(int ret) {
    cudaError_t err = cudaGetLastError();
    if (err == cudaSuccess)
        return;
    
    printf("\nCUDA ERROR (%d): %s (%d)\n", ret, cudaGetErrorString(err), (int)err);
	fprintf(logfile, "\nCUDA ERROR (%d): %s (%d)\n", ret, cudaGetErrorString(err), (int)err);
	fflush(logfile);
    exit(ret);
}

void ce(int ret, const char *format, ...) {
    cudaError_t err = cudaGetLastError();
    if (err == cudaSuccess)
        return;

    va_list arguments;
    va_start(arguments, format);
    printf("\nCUDA ERROR (%d): %s (%d)\nError message: ", ret, cudaGetErrorString(err), (int)err);
	fprintf(logfile, "\nCUDA ERROR (%d): %s (%d)\nError message: ", ret, cudaGetErrorString(err), (int)err);
    vprintf(format, arguments);
    vfprintf(logfile, format, arguments);
    va_end(arguments);
	fflush(logfile);
    exit(ret);
}

// get time
void getTime() {
    time_t rawtime;
	struct tm *timeinfo;
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	strftime(charBuffer, CHARLEN, "%I:%M:%S %p on %A, %B %d, %Y", timeinfo);
}

// check output prefix
void checkOutputPrefix(int ret, string output_prefix) {
	
    // get the absolute path of the output prefix
    GetFullPathName(output_prefix.c_str(), CHARLEN, charBuffer, NULL);
    output_prefix = string(charBuffer);
    
    string testFile = output_prefix + "test";
	uint i = 0;
	char buf[10];
	FILE *t;

    while (true) {
		sprintf(buf, "%d", i);
		fopen_s(&t, (testFile + string(buf)).c_str(), "r");
		// the test file already exists, don't change it
		if (t) {
			fclose(t);
			t = NULL;
			i++;
		}
		// the test file does not exist
		else {
			bool success = true;
			
			// try creating the file
			fopen_s(&t, (testFile + string(buf)).c_str(), "w");

			// file successfully created
			if (t) {
				// try writing to the file
				if (fprintf(t, "test") != 4)
					success = false;
				if (fclose(t))
					success = false;
				if (remove((testFile + string(buf)).c_str()))
					success = false;
				
			}
			// file could not be created
			else
				success = false;

			if (!success) {
				printf("ERROR (%d): Could not create output files with output prefix \"%s\". Please make sure the specified path exists and is writable.\n", ret, output_prefix.c_str());
				exit(ret);
			}
			else
				break;
		}
	}
}

/*
// create a bitmap file of the detector data with a value of 0 being black and the largest data value being white
void saveBitmap(const char *filePath, float *data, int num) {
    // find highest photon count
    double highestCount = data[0];
    int i;
    for (i = 0; i < num; i++) {
        if (data[i] > highestCount)
            highestCount = data[i];
    }

    // calculate rgb values
    char *rgb = (char*)malloc(3 * num);
    for (i = 0; i < num; i++) {
        char val = char (data[i] / highestCount * 255);
        rgb[3 * i] = val;
        rgb[3 * i + 1] = val;
        rgb[3 * i + 2] = val;
    }

    // write bitmap
    write_bmp(filePath, h_detector.N.x, h_detector.N.y, rgb);

    free(rgb);
}
*/

/*
// write the detector data as into a binary file
void write_output(string baseName, string type, double *data[NUM_DETECTOR_CAT + 1]) {

	uchar istart = NUM_DETECTOR_CAT; //only output the "total" image
    for (uchar i = istart; i <= NUM_DETECTOR_CAT; i++) {
        string fname = baseName + type + "_" + categories[i];
        
        FILE *out = fopen((fname + ".bin").c_str(), "wb");
        fwrite(&h_detector.N, sizeof(uint2), 1, out);
        
        // convert to float
        uint num_pixels = h_detector.N.x * h_detector.N.y;
        float *data_f = (float*)malloc(num_pixels * sizeof(float));
        for (uint j = 0; j < num_pixels; j++)
            data_f[j] = (float)data[i][j];
        
		// write data output
        fwrite(data_f, sizeof(float), num_pixels, out);
        fclose(out);

		// write bitmap output
        saveBitmap((fname + ".bmp").c_str(), data_f, num_pixels);

        free(data_f);
    }
}
*/

/*void write_dose(string baseName,  score_t *data) {

        //string fname ;
		char fname[256];
		sprintf(fname, "%s_dose_%03dx%03dx%03d.bin",baseName.c_str(),h_phantom_N.x,h_phantom_N.y,h_phantom_N.z);
		printf("%s\n",fname);
        FILE *out = fopen(fname, "wb");
		fwrite(data, sizeof(score_t), h_size_phantom,out);

		fclose(out);
}*/

void interp_dose(string baseName, float *out_data, interp_t interp){
	char fname[256];
	sprintf(fname, "%s_Interpolate_%03dx%03dx%03d.bin",baseName.c_str(),interp.dimen.x,interp.dimen.y,interp.dimen.z);
	//sprintf(fname,"%sdose.dos",output_dir);
	logstr("%s\n",fname);
	FILE *out = fopen(fname, "wb");
	uint interp_size=interp.dimen.x*interp.dimen.y*interp.dimen.z;
	float *dose_array = (float*)malloc(interp_size*sizeof(float));
	float x,y,z;
	//interpolation of z direction, results lie in lower x and lower y of each interpolated z plane
	uint offset_xy = interp.dimen.x*(interp.dimen.y-h_phantom_N.y)+interp.dimen.x-h_phantom_N.x;
	uint interp_xy = interp.dimen.x*interp.dimen.y;
	uint phant_xy = h_phantom_N.x*h_phantom_N.y;
	for(uint i=0; i<h_phantom_N.x; i++){
		for(uint j=0; j<h_phantom_N.y; j++){
			for(uint iz=0; iz<interp.dimen.z; iz++){
				z = interp.start.z + interp.sizes.z*iz;
				if(z<=z_bounds[0]||z>=z_bounds[h_phantom_N.z-1]){
					dose_array[offset_xy+i + interp.dimen.x*j + interp_xy*iz] = 0;
				}
				else{
					for(uint k=0; k<h_phantom_N.z-1; k++){
						if(z>z_bounds[k]&&z<=z_bounds[k+1]){
							dose_array[offset_xy+i + interp.dimen.x*j + interp_xy*iz]
								= out_data[i+h_phantom_N.x*j+phant_xy*k]
								+ (z-z_bounds[k])
									* (out_data[i+h_phantom_N.x*j+phant_xy*(k+1)]-out_data[i+h_phantom_N.x*j+phant_xy*k])
									/ (z_bounds[k+1]-z_bounds[k]);
							break;
						}
					}
				}
			}
		}
	}
	//interpolation of y direction, results lie in lower x of each interpolated x line(fix y and fix z)
	uint offset_x = interp.dimen.x-h_phantom_N.x;
	for(uint iz=0; iz<interp.dimen.z; iz++){
		for(uint i=0; i<h_phantom_N.x; i++){
			for(uint iy=0; iy<interp.dimen.y; iy++){
				y = interp.start.y + interp.sizes.y*iy;
				if(y<=y_bounds[0]||y>=y_bounds[h_phantom_N.y-1]){
					dose_array[offset_x+i + interp.dimen.x*iy + interp_xy*iz] = 0;
				}
				else{
					for(uint j=0; j<h_phantom_N.y-1; j++){
						if(y>y_bounds[j]&&y<=y_bounds[j+1]){
							dose_array[offset_x+i + interp.dimen.x*iy + interp_xy*iz]
								= dose_array[offset_xy+i + interp.dimen.x*j + interp_xy*iz]
								+ (y-y_bounds[j])
									*(dose_array[offset_xy+i + interp.dimen.x*(j+1) + interp_xy*iz]
										- dose_array[offset_xy+i + interp.dimen.x*j + interp_xy*iz])
									/ (y_bounds[j+1]-y_bounds[j]);
							break;
						}
					}
				}
			}
		}
	}
	//interpolation of x direction
	for(uint iz=0; iz<interp.dimen.z; iz++){
		for(uint iy=0; iy<interp.dimen.y; iy++){
			for(uint ix=0; ix<interp.dimen.x; ix++){
				x = interp.start.x + interp.sizes.x*ix;
				if(x<=x_bounds[0]||x>=x_bounds[h_phantom_N.x-1]){
					dose_array[ix + interp.dimen.x*iy + interp_xy*iz] = 0;
				}
				else{
					for(uint i=0; i<h_phantom_N.x-1; i++){
						if(x>x_bounds[i]&&x<=x_bounds[i+1]){
							dose_array[ix + interp.dimen.x*iy + interp_xy*iz]
								= dose_array[offset_x+i + interp.dimen.x*iy + interp_xy*iz]
								+ (x-x_bounds[i])
									* (dose_array[offset_x+(i+1) + interp.dimen.x*iy + interp_xy*iz]
										- dose_array[offset_x+i + interp.dimen.x*iy + interp_xy*iz])
									/ (x_bounds[i+1]-x_bounds[i]);
							break;
						}
					}
				}
			}
		}
	}
	int dimen[3];
	dimen[0]=interp.dimen.x;
	dimen[1]=interp.dimen.y;
	dimen[2]=interp.dimen.z;
	fwrite(dimen,sizeof(int),3,out);
	fwrite(dose_array, sizeof(float), interp_size, out);
	free(dose_array);

	fclose(out);
}


void output_RelDos(string baseName, float *out_data, interp_t dosout){
	float3 dosout_end;
	dosout_end.x = dosout.start.x+dosout.sizes.x*(dosout.dimen.x-1);
	dosout_end.y = dosout.start.y+dosout.sizes.y*(dosout.dimen.y-1);
	dosout_end.z = dosout.start.z+dosout.sizes.z*(dosout.dimen.z-1);
	bool no_output = false;
	if(dosout.start.x<=x_bounds[0]||dosout_end.x>=x_bounds[h_phantom_N.x-1]){
		logstr("Warning: some dose point in x direction outside phantom:\n");
		logstr("points start at %f, end at %f\n",dosout.start.x,dosout_end.x);
		logstr("phantom boundary is from %f to %f\n",x_bounds[0],x_bounds[h_phantom_N.x-1]);
		no_output = true;
	}
	if(dosout.start.y<=y_bounds[0]||dosout_end.y>=y_bounds[h_phantom_N.y-1]){
		logstr("Warning: some dose point in y direction outside phantom:\n");
		logstr("points start at %f, end at %f\n",dosout.start.y,dosout_end.y);
		logstr("phantom boundary is from %f to %f\n",y_bounds[0],y_bounds[h_phantom_N.y-1]);
		no_output = true;
	}
	if(dosout.start.z<=z_bounds[0]||dosout_end.z>=z_bounds[h_phantom_N.z-1]){
		logstr("Warning: some dose point in z direction outside phantom:\n");
		logstr("points start at %f, end at %f\n",dosout.start.z,dosout_end.z);
		logstr("phantom boundary is from %f to %f\n",z_bounds[0],z_bounds[h_phantom_N.z-1]);
		no_output = true;
	}
	if(no_output) return;

	char fname[256];
	sprintf(fname, "%s_RelativeDose_%03dx%03dx%03d.txt",baseName.c_str(),dosout.dimen.x,dosout.dimen.y,dosout.dimen.z);
	logstr("%s\n",fname);
	FILE *out = fopen(fname, "w");
	uint start_x,dimen_x;
	uint start_y,dimen_y;
	//uint start_z,dimen_z;
	//find out x and y indices of phantom which dose points stay in
	for(start_x=0; start_x<h_phantom_N.x-1; start_x++)
		if(dosout.start.x>x_bounds[start_x]&&dosout.start.x<=x_bounds[start_x+1]) break;
	if(dosout.dimen.x==1) dimen_x=1;
	else{
		for(dimen_x=start_x; dimen_x<h_phantom_N.x-1; dimen_x++)
			if(dosout_end.x>x_bounds[dimen_x]&&dosout_end.x<=x_bounds[dimen_x+1]) break;
		dimen_x = dimen_x-start_x+1;
	}
	for(start_y=0; start_y<h_phantom_N.y-1; start_y++)
		if(dosout.start.y>y_bounds[start_y]&&dosout.start.y<=y_bounds[start_y+1]) break;
	if(dosout.dimen.y==1) dimen_y=1;
	else{
		for(dimen_y=start_y; dimen_y<h_phantom_N.y-1; dimen_y++)
			if(dosout_end.y>y_bounds[dimen_y]&&dosout_end.y<=y_bounds[dimen_y+1]) break;
		dimen_y = dimen_y-start_y+1;
	}
	/*
	for(start_z=0; start_z<h_phantom_N.z-1; start_z++)
		if(dosout.start.z>z_bounds[start_z]&&dosout.start.z<=z_bounds[start_z+1]) break;
	if(dosout.dimen.z==1) dimen_z=1;
	else{
		for(dimen_z=start_z; dimen_z<h_phantom_N.z-1; dimen_z++)
			if(dosout_end.z>z_bounds[dimen_z]&&dosout_end.z<=z_bounds[dimen_z+1]) break;
		dimen_z = dimen_z-start_z+1;
	}
	*/
	//interpolation in z direction, results lie in interp_z
	uint dosout_size = (dimen_x+1)*(dimen_y+1)*dosout.dimen.z;
	float *interp_z = (float*)malloc(dosout_size*sizeof(float));
	float x,y,z;
	uint dimen_xy = (dimen_x+1)*(dimen_y+1);
	uint phant_xy = h_phantom_N.x*h_phantom_N.y;
	for(uint i=0; i<=dimen_x; i++){
		for(uint j=0; j<=dimen_y; j++){
			uint k=0;
			for(uint iz=0; iz<dosout.dimen.z; iz++){
				z = dosout.start.z + dosout.sizes.z*iz;
				for(;k<h_phantom_N.z-1;k++)
					if(z>z_bounds[k]&&z<=z_bounds[k+1]) break;
				interp_z[i + (dimen_x+1)*j + dimen_xy*iz]
					= out_data[(start_x+i)+h_phantom_N.x*(start_y+j)+phant_xy*k]
					+ (z-z_bounds[k])
						* (out_data[(start_x+i)+h_phantom_N.x*(start_y+j)+phant_xy*(k+1)]
							- out_data[(start_x+i)+h_phantom_N.x*(start_y+j)+phant_xy*k])
						/ (z_bounds[k+1]-z_bounds[k]);
			}
		}
	}
	//interpolation in y direction, results lie in interp_y
	dosout_size = (dimen_x+1)*dosout.dimen.y*dosout.dimen.z;
	uint dimen_xj = (dimen_x+1)*dosout.dimen.y;
	float *interp_y = (float*)malloc(dosout_size*sizeof(float));
	for(uint iz=0; iz<dosout.dimen.z; iz++){
		for(uint i=0; i<=dimen_x; i++){
			uint j=0;
			for(uint iy=0; iy<dosout.dimen.y; iy++){
				y = dosout.start.y + dosout.sizes.y*iy;
				for(;j<h_phantom_N.y-1;j++)
					if(y>y_bounds[j]&&y<=y_bounds[j+1]) break;
				interp_y[i + (dimen_x+1)*iy + dimen_xj*iz]
					= interp_z[i + (dimen_x+1)*(j-start_y) + dimen_xy*iz]
					+ (y-y_bounds[j])
						* (interp_z[i + (dimen_x+1)*(j+1-start_y) + dimen_xy*iz]
							- interp_z[i + (dimen_x+1)*(j-start_y) + dimen_xy*iz])
						/ (y_bounds[j+1]-y_bounds[j]);
			}
		}
	}
	free(interp_z);
	//interpolation in x direction, results lie interp_x
	dosout_size = dosout.dimen.x*dosout.dimen.y*dosout.dimen.z;
	uint dimen_ij = dosout.dimen.x*dosout.dimen.y;
	float *interp_x = (float*)malloc(dosout_size*sizeof(float));
	for(uint iz=0; iz<dosout.dimen.z; iz++){
		for(uint iy=0; iy<dosout.dimen.y; iy++){
			uint i=0;
			for(uint ix=0; ix<dosout.dimen.x; ix++){
				x = dosout.start.x + dosout.sizes.x*ix;
				for(;i<h_phantom_N.x-1;i++)
					if(x>x_bounds[i]&&x<=x_bounds[i+1]) break;
				interp_x[ix + dosout.dimen.x*iy + dimen_ij*iz]
					= interp_y[(i-start_x) + (dimen_x+1)*iy + dimen_xj*iz]
					+ (x-x_bounds[i])
						* (interp_y[(i+1-start_x) + (dimen_x+1)*iy + dimen_xj*iz]
							- interp_y[(i-start_x) + (dimen_x+1)*iy + dimen_xj*iz])
						/ (x_bounds[i+1]-x_bounds[i]);
			}
		}
	}
	free(interp_y);

	float *meas_data;
	float dose_meas_max = 0.0F;
	//float dose_calc_max = 0.0F;
	bool ok = true;
	ifstream infile;

	infile.open(meas_data_file);
    if(infile.is_open())//check if file is opened
    {
		//while(infile.good()&&!infile.eof())//while file is in good condition and not end of file
		if(infile.good())
		{
			meas_data = (float*)malloc(dosout_size*sizeof(float));
			//memset(meas_data,0,dosout_size*sizeof(float));
			if(norm_type==1){  //normalization to maximum dose
				float dose_calc_max = 0.0F;
				for (uint i = 0; i < dosout_size; i++) {
					infile>>meas_data[i];
					infile.ignore(100,'\n');//ignore the '\n' of this line
					if(dose_meas_max<meas_data[i]) dose_meas_max = meas_data[i];  //find out the maximum of measured dose
					if(dose_calc_max<interp_x[i])  dose_calc_max = interp_x[i];   //find out the maximum of calculated dose
				}
				for (uint i = 0; i < dosout_size; i++) {
					meas_data[i] /= dose_meas_max;  //normal to maximum dose
					interp_x[i]  /= dose_calc_max;  //normal to maximum dose
				}
				dose_meas_max = 1.0F;
			}
			else if(norm_type==2){  //normalization to center point
				for (uint i = 0; i < dosout_size; i++) {
					infile>>meas_data[i];
					infile.ignore(100,'\n');//ignore the '\n' of this line
					if(dose_meas_max<meas_data[i]) dose_meas_max = meas_data[i];  //find out the maximum of measured dose
				}
				//uint center = dosout_size/2;
				uint center = 0;
				if(dosout.sizes.x!=0)
					center += int((h_isocenter_location.x-dosout.start.x)/dosout.sizes.x+0.1);
				if(dosout.sizes.y!=0)
					center += int((h_isocenter_location.y-dosout.start.y)/dosout.sizes.y+0.1)*dosout.dimen.x;
				if(dosout.sizes.z!=0)
					center += int((h_isocenter_location.z-dosout.start.z)/dosout.sizes.z+0.1)*dosout.dimen.x*dosout.dimen.y;
				printf("center=%d",center);
				float dose_meas_center = meas_data[center];
				float dose_calc_center = interp_x[center];
				for (uint i = 0; i < dosout_size; i++) {
					meas_data[i] /= dose_meas_center;  //normal to center point
					interp_x[i]  /= dose_calc_center;  //normal to center point
				}
				dose_meas_max /= dose_meas_center;
			}
			else{  //output directly
				for (uint i = 0; i < dosout_size; i++) {
					infile>>meas_data[i];
					infile.ignore(100,'\n');//ignore the '\n' of this line
					meas_data[i] *= 100.0F;           //change the measured dose from cGy to Gy
					if(dose_meas_max<meas_data[i]) dose_meas_max = meas_data[i];  //find out the maximum of measured dose
				}
			}
		}
	}
	else
	{
		ok = false;
		logstr("Error: dose data file reading error!\n");
	}
    infile.close();//close the file

	if(ok){
		int pass;
		float percent;
		uint total = 0;
		uint count = 0;
		dose_meas_max *= 0.1F;
		if(array_type==1){
			for (uint i = 0; i < dosout_size; i++) {
				pass = 0;
				if(meas_data[i]>dose_meas_max){
					total ++;
					if(abs(interp_x[i]-meas_data[i])/meas_data[i]<=GmIdx_dos){
						pass = 1;
						count ++;
					}
					else
						pass = -1;
				}
				fprintf(out,"%f\t %f\t %d\n",meas_data[i],interp_x[i],pass);
			}
			percent = (float)count/total;
		}
		else{
			//int iz,ix;
			for (uint i = 0; i < dosout_size; i++) {
				pass = 0;
				float value;
				if(meas_data[i]>dose_meas_max){
					total ++;
					if(abs(interp_x[i]-meas_data[i])/meas_data[i]<=GmIdx_dos){
						pass ++;
					}
					int iz = i / dosout.dimen.x;
					int ix = i - iz * dosout.dimen.x;
					if(ix>0&&pass<=0){
						value = abs(interp_x[i-1]-meas_data[i])/meas_data[i]/GmIdx_dos;
						value *= value;
						value = sqrt(value+0.25*0.25/GmIdx_len/GmIdx_len);
						if(value<=1)
							pass ++;
					}
					if(ix+1<dosout.dimen.x&&pass<=0){
						value = abs(interp_x[i+1]-meas_data[i])/meas_data[i]/GmIdx_dos;
						value *= value;
						value = sqrt(value+0.25*0.25/GmIdx_len/GmIdx_len);
						if(value<=1)
							pass ++;
					}
					if(iz>0){
						value = abs(interp_x[i-dosout.dimen.x]-meas_data[i])/meas_data[i]/GmIdx_dos;
						value *= value;
						value = sqrt(value+0.25*0.25/GmIdx_len/GmIdx_len);
						if(value<=1)
							pass ++;
					}
					if(iz+1<dosout.dimen.z&&pass<=0){
						value = abs(interp_x[i+dosout.dimen.x]-meas_data[i])/meas_data[i]/GmIdx_dos;
						value *= value;
						value = sqrt(value+0.25*0.25/GmIdx_len/GmIdx_len);
						if(value<=1)
							pass ++;
					}
					if(pass>0)
						count ++;
					else
						pass = -1;
				}
				fprintf(out,"%f\t %f\t %d\n",meas_data[i],interp_x[i],pass);
			}
			percent = (float)count/total;
		}
		fprintf(out,"Totally %d points compared, %d passed, percentage = %.0f %%\n",total,count,percent*100);
		logstr("Totally %d points compared, %d passed, percentage = %.0f %%\n",total,count,percent*100);
	}
	else{
		for(uint i=0;i<dosout_size;i++){
			fprintf(out,"%f\n",interp_x[i]);
		}
	}

	if(ok) free(meas_data);
	free(interp_x);
	fclose(out);
}


void output_RelDos_II(string baseName, float *out_data, interp_t dosout){

	char fname[256];
	sprintf(fname, "%s_RelativeDose_%03dx%03dx%03d.txt",baseName.c_str(),dosout.dimen.x,dosout.dimen.y,dosout.dimen.z);
	logstr("%s\n",fname);
	FILE *out = fopen(fname, "w");

	float pos_x,pos_y,pos_z;
	int   idx_x,idx_y,idx_z;
	uint  idx_phan,idx_intp;
	uint dosout_size = dosout.dimen.x*dosout.dimen.y*dosout.dimen.z;
	float *interp_x = (float*)malloc(dosout_size*sizeof(float));
	for(int k=0; k<dosout.dimen.z; k++){
		for(int j=0; j<dosout.dimen.y; j++){
			for(int i=0; i<dosout.dimen.x; i++){
				pos_x = dosout.start.x + dosout.sizes.x*i;
				pos_y = dosout.start.y + dosout.sizes.y*j;
				pos_z = dosout.start.z + dosout.sizes.z*k;
				idx_x = int((pos_x-x_bounds[0])/(x_bounds[1]-x_bounds[0])+0.1);
				idx_y = int((pos_y-y_bounds[0])/(y_bounds[1]-y_bounds[0])+0.1);
				idx_z = int((pos_z-z_bounds[0])/(z_bounds[1]-z_bounds[0])+0.1);
				idx_phan = h_phantom_N.x*h_phantom_N.y*idx_z + h_phantom_N.x*idx_y + idx_x;
				idx_intp = dosout.dimen.x*dosout.dimen.y*k + dosout.dimen.x*j + i;
				interp_x[idx_intp] = out_data[idx_phan];
				//fprintf(out,"%f\t %f\t %f\n",x_bounds[idx_x],y_bounds[idx_y],z_bounds[idx_z]);
				//fprintf(out,"%d\t %d\t %d\n",idx_x,idx_y,idx_z);
			}
		}
	}


	float *meas_data;
	float dose_meas_max = 0.0F;
	//float dose_calc_max = 0.0F;
	bool ok = true;
	ifstream infile;

	infile.open(meas_data_file);
    if(infile.is_open())//check if file is opened
    {
		//while(infile.good()&&!infile.eof())//while file is in good condition and not end of file
		if(infile.good())
		{
			meas_data = (float*)malloc(dosout_size*sizeof(float));
			//memset(meas_data,0,dosout_size*sizeof(float));
			if(norm_type==1){  //normalization to maximum dose
				float dose_calc_max = 0.0F;
				for (uint i = 0; i < dosout_size; i++) {
					infile>>meas_data[i];
					infile.ignore(100,'\n');//ignore the '\n' of this line
					if(dose_meas_max<meas_data[i]) dose_meas_max = meas_data[i];  //find out the maximum of measured dose
					if(dose_calc_max<interp_x[i])  dose_calc_max = interp_x[i];   //find out the maximum of calculated dose
				}
				for (uint i = 0; i < dosout_size; i++) {
					meas_data[i] /= dose_meas_max;  //normal to maximum dose
					interp_x[i]  /= dose_calc_max;  //normal to maximum dose
				}
				dose_meas_max = 1.0F;
			}
			else if(norm_type==2){  //normalization to center point
				for (uint i = 0; i < dosout_size; i++) {
					infile>>meas_data[i];
					infile.ignore(100,'\n');//ignore the '\n' of this line
					if(dose_meas_max<meas_data[i]) dose_meas_max = meas_data[i];  //find out the maximum of measured dose
				}
				//uint center = dosout_size/2;
				uint center = 0;
				if(dosout.sizes.x!=0)
					center += int((h_isocenter_location.x-dosout.start.x)/dosout.sizes.x+0.1);
				if(dosout.sizes.y!=0)
					center += int((h_isocenter_location.y-dosout.start.y)/dosout.sizes.y+0.1)*dosout.dimen.x;
				if(dosout.sizes.z!=0)
					center += int((h_isocenter_location.z-dosout.start.z)/dosout.sizes.z+0.1)*dosout.dimen.x*dosout.dimen.y;
				//printf("center=%d",center);
				float dose_meas_center = meas_data[center];
				float dose_calc_center = interp_x[center];
				for (uint i = 0; i < dosout_size; i++) {
					meas_data[i] /= dose_meas_center;  //normal to center point
					interp_x[i]  /= dose_calc_center;  //normal to center point
				}
				dose_meas_max /= dose_meas_center;
			}
			else{  //output directly
				for (uint i = 0; i < dosout_size; i++) {
					infile>>meas_data[i];
					infile.ignore(100,'\n');//ignore the '\n' of this line
					meas_data[i] *= 100.0F;           //change the measured dose from cGy to Gy
					if(dose_meas_max<meas_data[i]) dose_meas_max = meas_data[i];  //find out the maximum of measured dose
				}
			}
		}
	}
	else
	{
		ok = false;
		logstr("Error: dose data file reading error!\n");
	}
    infile.close();//close the file

	if(ok){
		int pass;
		float percent;
		uint total = 0;
		uint count = 0;
		dose_meas_max *= 0.1F;
		if(array_type==1){
			for (uint i = 0; i < dosout_size; i++) {
				pass = 0;
				if(meas_data[i]>dose_meas_max){
					total ++;
					if(abs(interp_x[i]-meas_data[i])/meas_data[i]<=GmIdx_dos){
						pass = 1;
						count ++;
					}
					else
						pass = -1;
				}
				fprintf(out,"%f\t %f\t %d\n",meas_data[i],interp_x[i],pass);
			}
			percent = (float)count/total;
		}
		else{
			//int iz,ix;
			for (uint i = 0; i < dosout_size; i++) {
				pass = 0;
				float value;
				if(meas_data[i]>dose_meas_max){
					total ++;
					if(abs(interp_x[i]-meas_data[i])/meas_data[i]<=GmIdx_dos){
						pass ++;
					}
					int iz = i / dosout.dimen.x;
					int ix = i - iz * dosout.dimen.x;
					if(ix>0&&pass<=0){
						value = abs(interp_x[i-1]-meas_data[i])/meas_data[i]/GmIdx_dos;
						value *= value;
						value = sqrt(value+0.25*0.25/GmIdx_len/GmIdx_len);
						if(value<=1)
							pass ++;
					}
					if(ix+1<dosout.dimen.x&&pass<=0){
						value = abs(interp_x[i+1]-meas_data[i])/meas_data[i]/GmIdx_dos;
						value *= value;
						value = sqrt(value+0.25*0.25/GmIdx_len/GmIdx_len);
						if(value<=1)
							pass ++;
					}
					if(iz>0){
						value = abs(interp_x[i-dosout.dimen.x]-meas_data[i])/meas_data[i]/GmIdx_dos;
						value *= value;
						value = sqrt(value+0.25*0.25/GmIdx_len/GmIdx_len);
						if(value<=1)
							pass ++;
					}
					if(iz+1<dosout.dimen.z&&pass<=0){
						value = abs(interp_x[i+dosout.dimen.x]-meas_data[i])/meas_data[i]/GmIdx_dos;
						value *= value;
						value = sqrt(value+0.25*0.25/GmIdx_len/GmIdx_len);
						if(value<=1)
							pass ++;
					}
					if(pass>0)
						count ++;
					else
						pass = -1;
				}
				fprintf(out,"%f\t %f\t %d\n",meas_data[i],interp_x[i],pass);
			}
			percent = (float)count/total;
		}
		fprintf(out,"Totally %d points compared, %d passed, percentage = %.0f %%\n",total,count,percent*100);
		logstr("Totally %d points compared, %d passed, percentage = %.0f %%\n",total,count,percent*100);
	}
	else{
		for(uint i=0;i<dosout_size;i++){
			fprintf(out,"%f\n",interp_x[i]);
		}
	}

	if(ok) free(meas_data);
	free(interp_x);
	fclose(out);
}

void output_AbsDos(float *out_data, float3 dose_pos){
	uint i,j,k;
	for(i=1; i<h_phantom_N.x-2; i++){
		if(dose_pos.x>x_bounds[i]&&dose_pos.x<=x_bounds[i+1]) break;
	}
	if(i==h_phantom_N.x-2){
		if(dose_pos.x<=x_bounds[1]) i=0;
		else i=h_phantom_N.x-2;
	}
	for(j=1; j<h_phantom_N.y-2; j++){
		if(dose_pos.y>y_bounds[j]&&dose_pos.y<=y_bounds[j+1]) break;
	}
	if(j==h_phantom_N.y-2){
		if(dose_pos.y<=y_bounds[1]) j=0;
		else j=h_phantom_N.y-2;
	}
	for(k=1; k<h_phantom_N.z-2; k++){
		if(dose_pos.z>z_bounds[k]&&dose_pos.z<=z_bounds[k+1]) break;
	}
	if(k==h_phantom_N.z-2){
		if(dose_pos.z<=z_bounds[1]) k=0;
		else k=h_phantom_N.z-2;
	}
	//logstr("i=%d, j=%d, k=%d\n",i,j,k);
	float Abs_Dose[4];
	uint phant_xy = h_phantom_N.x*h_phantom_N.y;
	//logstr("x_bounds[%d]=%f, x_bounds[%d]=%f\n",i,x_bounds[i],i+1,x_bounds[i+1]);
	//logstr("y_bounds[%d]=%f, y_bounds[%d]=%f\n",j,y_bounds[j],j+1,y_bounds[j+1]);
	//logstr("z_bounds[%d]=%f, z_bounds[%d]=%f\n",k,z_bounds[k],k+1,z_bounds[k+1]);
	//logstr("data[%d,%d,%d]=%.1f, data[%d,%d,%d]=%.1f\n",i,  j,  k,out_data[i  +h_phantom_N.x*j    +phant_xy*k],i,  j,  k+1,out_data[i  +h_phantom_N.x*j    +phant_xy*(k+1)]);
	//logstr("data[%d,%d,%d]=%.1f, data[%d,%d,%d]=%.1f\n",i+1,j,  k,out_data[i+1+h_phantom_N.x*j    +phant_xy*k],i+1,j,  k+1,out_data[i+1+h_phantom_N.x*j    +phant_xy*(k+1)]);
	//logstr("data[%d,%d,%d]=%.1f, data[%d,%d,%d]=%.1f\n",i,  j+1,k,out_data[i  +h_phantom_N.x*(j+1)+phant_xy*k],i,  j+1,k+1,out_data[i  +h_phantom_N.x*(j+1)+phant_xy*(k+1)]);
	//logstr("data[%d,%d,%d]=%.1f, data[%d,%d,%d]=%.1f\n",i+1,j+1,k,out_data[i+1+h_phantom_N.x*(j+1)+phant_xy*k],i+1,j+1,k+1,out_data[i+1+h_phantom_N.x*(j+1)+phant_xy*(k+1)]);
	Abs_Dose[0] = out_data[i+h_phantom_N.x*j+phant_xy*k] + (dose_pos.z-z_bounds[k])
		* (out_data[i+h_phantom_N.x*j+phant_xy*(k+1)] - out_data[i+h_phantom_N.x*j+phant_xy*k])
		/ (z_bounds[k+1]-z_bounds[k]);
	Abs_Dose[1] = out_data[i+1+h_phantom_N.x*j+phant_xy*k] + (dose_pos.z-z_bounds[k])
		* (out_data[i+1+h_phantom_N.x*j+phant_xy*(k+1)] - out_data[i+1+h_phantom_N.x*j+phant_xy*k])
		/ (z_bounds[k+1]-z_bounds[k]);
	Abs_Dose[2] = out_data[i+h_phantom_N.x*(j+1)+phant_xy*k] + (dose_pos.z-z_bounds[k])
		* (out_data[i+h_phantom_N.x*(j+1)+phant_xy*(k+1)] - out_data[i+h_phantom_N.x*(j+1)+phant_xy*k])
		/ (z_bounds[k+1]-z_bounds[k]);
	Abs_Dose[3] = out_data[i+1+h_phantom_N.x*(j+1)+phant_xy*k] + (dose_pos.z-z_bounds[k])
		* (out_data[i+1+h_phantom_N.x*(j+1)+phant_xy*(k+1)] - out_data[i+1+h_phantom_N.x*(j+1)+phant_xy*k])
		/ (z_bounds[k+1]-z_bounds[k]);
	Abs_Dose[0] = Abs_Dose[0] + (dose_pos.y-y_bounds[j]) * (Abs_Dose[2]-Abs_Dose[0]) / (y_bounds[j+1]-y_bounds[j]);
	Abs_Dose[1] = Abs_Dose[1] + (dose_pos.y-y_bounds[j]) * (Abs_Dose[3]-Abs_Dose[1]) / (y_bounds[j+1]-y_bounds[j]);
	Abs_Dose[0] = Abs_Dose[0] + (dose_pos.x-x_bounds[i]) * (Abs_Dose[1]-Abs_Dose[0]) / (x_bounds[i+1]-x_bounds[i]);
	//Abs_Dose[0] = Abs_Dose[0]/100.0F;

	//bollow the last variable Abs_Dose[3]
	Abs_Dose[3] = (x_bounds[i]+x_bounds[i+1])/2;
	if(dose_pos.x>Abs_Dose[3]) i++;
	Abs_Dose[3] = (y_bounds[j]+y_bounds[j+1])/2;
	if(dose_pos.y>Abs_Dose[3]) j++;
	Abs_Dose[3] = (z_bounds[k]+x_bounds[k+1])/2;
	if(dose_pos.z>Abs_Dose[3]) k++;
	Abs_Dose[2] = out_data[i+h_phantom_N.x*j+phant_xy*k];

	if(calibr_on){
		if(dose_meas==0.0F){
			logstr("Absolute dose at (%.2f, %.2f, %.2f) = %.1f cGy\n",dose_pos.x,dose_pos.y,dose_pos.z,Abs_Dose[0]);
			logstr("   voxel dose of (i=%2d, j=%2d, k=%2d) = %.1f cGy\n",i,j,k,Abs_Dose[2]);
		}
		else{

			//dose_meas *= 100.0F;
			Abs_Dose[1] = abs(Abs_Dose[0]-dose_meas)/dose_meas;
			Abs_Dose[3] = abs(Abs_Dose[2]-dose_meas)/dose_meas;

			logstr("Absolute dose at (%.2f, %.2f, %.2f):\n",dose_pos.x,dose_pos.y,dose_pos.z);
			logstr("  measured = %.1f cGy, calculated = %.1f cGy, difference = %.1f %%\n",dose_meas,Abs_Dose[0],Abs_Dose[1]*100.0F);
			logstr("  voxel dose of index (%2d, %2d, %2d) = %.1f cGy, difference = %.1f %%\n",i,j,k,Abs_Dose[2],Abs_Dose[3]*100.0F);

			/*
			Abs_Dose[2] = Abs_Dose[0];
			Abs_Dose[3] = Abs_Dose[1];
			while(Abs_Dose[3]>=0.02) Abs_Dose[3] *= 0.95;
			if(Abs_Dose[1]>=0.02F){
				if(Abs_Dose[2]>dose_meas) Abs_Dose[2] = dose_meas + dose_meas*Abs_Dose[3];
				else                      Abs_Dose[2] = dose_meas - dose_meas*Abs_Dose[3];
			}

			logstr("  measured = %.1f cGy, calculated = %.1f cGy, difference = %.1f %%\n",dose_meas,Abs_Dose[2],Abs_Dose[3]*100.0F);
			*/
		}
	}
	else{
		logstr("Absolute dose at (%.2f, %.2f, %.2f) = %.3E cGy\n",dose_pos.x,dose_pos.y,dose_pos.z,Abs_Dose[0]);
		logstr("   voxel dose of (i=%2d, j=%2d, k=%2d) = %.3E cGy\n",i,j,k,Abs_Dose[2]);
	}
}

void write_dose(string baseName,  score_t *data) {

	//string fname ;
	char fname[256];
	sprintf(fname, "%s_DoseResults_%03dx%03dx%03d.bin",baseName.c_str(),h_phantom_N.x,h_phantom_N.y,h_phantom_N.z);
	logstr("%s\n",fname);
	FILE *out = fopen(fname, "wb");

	float *out_data;
	if(sizeof(score_t)!=sizeof(float))  //if data is in double, convert them to float
		out_data = (float*)malloc(h_size_phantom*sizeof(float));
	else
		out_data = (float*)data;

	//float dose_calc_max = 0.0F;
	for(uint i=0; i<h_size_phantom; i++){
		//out_data[i] = (float)data[i];
		out_data[i] = (float)data[i]/voxel_mass[i]*dose_unit;
		//if(dose_calc_max<out_data[i]) dose_calc_max=out_data[i];
	}
	//logstr("Maximum calculated dose = %f\n",dose_calc_max);

	/*
	dose_calc_max *= 0.1F;
	for(uint i=0; i<h_size_phantom; i++){
		if(out_data[i]<=dose_calc_max) out_data[i] = 0.0F;
	}
	*/

	fwrite(out_data, sizeof(float), h_size_phantom, out);
	fclose(out);
	
	if(interp_on) interp_dose(baseName,out_data,interp);

	if(AbsDos_on)
		output_AbsDos(out_data,dose_pos);

	if(RelDos_on){
		if(RelDos.dimen.x==1&&RelDos.dimen.y==1&&RelDos.dimen.z==1)
			output_AbsDos(out_data,RelDos.start);
		else
			output_RelDos_II(baseName,out_data,RelDos);
	}

	if(sizeof(score_t)!= sizeof(float)) free(out_data);  //if data is in double, convert them to float
}

void write_uncertainty(string baseName,  score_t *data) {

        //string fname ;
		char fname[256];
		sprintf(fname, "%s_uncertianty_%03dx%03dx%03d.bin",baseName.c_str(),h_phantom_N.x,h_phantom_N.y,h_phantom_N.z);
		logstr("%s\n\n",fname);
        FILE *out = fopen(fname, "wb");
		float *out_data;
		if(sizeof(score_t)!=sizeof(float))  //if data is in double, convert them to float 
		{
			out_data = (float*)malloc(h_size_phantom*sizeof(float));
			for(int i=0; i<h_size_phantom; i++)
				out_data[i] = (float)data[i];
		}
		else
		{
			out_data = (float*)data;
			for(int i=0; i<h_size_phantom; i++)
				out_data[i] = (float)data[i];
		}

		fwrite(out_data, sizeof(float), h_size_phantom, out);
		if(sizeof(score_t)!= sizeof(float)) free(out_data);  //if data is in double, convert them to float 

		fclose(out);
}

void write_density(string baseName,  float *data) {

//        string fname ;
		char fname[256];
		sprintf(fname, "%s_Density_%03dx%03dx%03d.bin",baseName.c_str(),h_phantom_N.x,h_phantom_N.y,h_phantom_N.z);
		printf("%s\n",fname);
        FILE *out = fopen(fname, "wb");
		fwrite(data, sizeof(float), h_size_phantom,out);

		fclose(out);
}
void write_media_indices(string baseName,  int *data) {

//        string fname ;
		char fname[256];
		sprintf(fname, "%s_MediaIndex_%03dx%03dx%03d.bin",baseName.c_str(),h_phantom_N.x,h_phantom_N.y,h_phantom_N.z);
		printf("%s\n",fname);
        FILE *out = fopen(fname, "wb");
		fwrite(data, sizeof(int), h_size_phantom,out);

		fclose(out);
}

#endif
