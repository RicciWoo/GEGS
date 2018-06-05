/* This file basically implements the subroutine HATCH of the original EGSnrc code.
 * HATCH is a pretty long subroutine and calls various other subroutines which are
 * also implemented here. HATCH reads the PEGS4 file and initializes all the media
 * data. Here, this is done with the function init_media at the end of the file.
 * Most of the code was taken from the file egsnrc.mortran (v 1.72 2011/05/05) of 
 * the EGSnrc code system. Additional code was added to copy data to the device and
 * read certain files.
 */

#ifdef CUDA_EGS

#define MXGE        200
#define MXRAYFF     100
#define RAYCDFSIZE  100
#define MXELEMENT   100
#define MXSHELL     6

typedef struct element_t {
    char    symbol[2];
    double  z;
    double  wa;
    double  pz;
    double  rhoz;
} element_t;

/*
//orignal defination
typedef struct medium_t {
    string  *name;
    double  rho;
    double  rlc;
    double  ae, ap;
    double  ue, up;
    uint    ne;
    element_t *elements;
    int     iunrst, epstfl, iaprim;
} medium_t;
*/

typedef struct medium_t {
    string  *name;
    double  rho; //mass density of a given medium
    double  rlc; //radiation length in centimeters for a given medium
    double  ae,ap; //ae electron creation threshold energy   ap photon creation threshold energy
    double  ue,up; // ue upper electron energy in PEGS4 data set  up upper photon energy in PEGS4 data set
    uint    ne;  //number of energy bins
    element_t *elements;
    //int iunrst,epstfl,iaprim;
	int msge,mge,mseke,meke,mleke,mcmfp,mrange,irayl;   //inurst flag for type of stopping power  epstfl "flag for ICRU37 collision stopping powers iaprim flag for ICRU37 radiative stopping powers
	double te, thmoll;
	double dl1[6],dl2[6],dl3[6],dl4[6],dl5[6],dl6[6];//alphi[2],bpar[2],delpos[2];
	double alphi[2],bpar[2],delpos[2];
	double delcm;
	double xr0,teff0,blcc,xcc;
	double eke0,eke1;
    double esig0[150],esig1[150],psig0[150],psig1[150],ededx0[150],ededx1[150],pdedx0[150],pdedx1[150],ebr10[150],ebr11[150],pbr10[150],pbr11[150],pbr20[150],pbr21[150],tmxs0[150],tmxs1[150];
	double ebinda,ge0,ge1;
	double gmfp0[200],gmfp1[200],gbr10[200],gbr11[200],gbr20[200],gbr21[200];
	int ngr;
	double rco0,rco1;
	double rsct0[100],rsct1[100];
	double cohe0[200],cohe1[200];
	double rldu;
	double zbrang;
	double pz[MXEL],zelem[MXEL];
	int nmed;
} medium_t;

typedef struct __align__(16) rayleigh_data_t {
    float   xgrid;
    float   fcum;
    float   b_array;
    float   c_array;
} rayleigh_data_t;

double binding_energies2D[MXSHELL][MXELEMENT];
/*
common /PHOTOIN/


EBINDA, "energy of the K-edge for a given medium"
GE0,GE1, "used for indexing in logarithmic interpolations"
GMFP0,GMFP1, "used for gamma MFP interpolation"
GBR10,GBR11, "used for branching into pair interpolation"
GBR20,GBR21, "used for branching into Compton interpolation"
RCO0,RCO1, "used for indexing in momentum trans. sampling in Rayleigh"
RSCT0,RSCT1, "used for interpolation of momentum trans. func. in R"
COHE0,COHE1, "used for Rayleigh modification interpolation"
DPMFP; "number of MFP's to go to the next interaction"
$INTEGER
MPGEM, "??? "
NGR; "array size for Rayleigh scattering data"
*/
float2          **d_ge, **d_gmfp, **d_gbr1, **d_gbr2, **d_cohe, **d_pmax;
rayleigh_data_t **d_rayleigh_data;
uint            **d_i_array;

__device__ float2             *ge, *gmfp, *gbr1, *gbr2, *cohe, *pmax;
__device__ rayleigh_data_t    *rayleigh_data;
__device__ uint               *i_array;
/*
__constant__ float2             *ge, *gmfp, *gbr1, *gbr2, *cohe, *pmax;
__constant__ rayleigh_data_t    *rayleigh_data;
__constant__ uint               *i_array;
*/

ulong abs_diff(ulong a, ulong b) {
    if (a > b)
        return a - b;
    else
        return b - a;
}

void read_xsec_data(const char *file, uint ndat[MXELEMENT], double2 *data[MXELEMENT]) {
    FILE *f = fopen(file, "r");
    bool ok = f > 0;
    
    if (ok) {
        for (uint i = 0; i < MXELEMENT; i++) {
            uint n;
            if (fscanf(f, "%u\n", &n) != 1) {
                ok = false;
                break;
            }
            ndat[i] = n;      
            data[i] = (double2*)malloc(n * sizeof(double2));

            for (uint j = 0; j < n; j++) {
                double2 dat;
                if (fscanf(f, "%lf %lf", &dat.x, &dat.y) != 2) {
                    ok = false;
                    break;
                }
                data[i][j] = dat;
            }

            if (!ok)
                break;  
        }
    }

    if (f)
        fclose(f);

    if (!ok)
        error(10001, "Could not read the data file \"%s\".\n", file);
}

void read_ff_data(const char *file, double xval[MXRAYFF], double aff[MXELEMENT][MXRAYFF]) {
    // read atomic form factors from file
    string ff_file = string(data_dir) + string(file);
    FILE *f = fopen(ff_file.c_str(), "r");

    bool ok = f > 0;
    
    if (ok) {
        // read momentum transfer values
        for (uint i = 0; i < MXRAYFF; i++) {
            if (fscanf(f, "%lf", &xval[i]) != 1) {
                ok = false;
                break;
            }
        }
    }

    if (ok) {
        // read elemental form factors
        for (uint i = 0; i < MXELEMENT; i++) {
            for (uint j = 0; j < MXRAYFF; j++) {
                if (fscanf(f, "%lf", &aff[i][j]) != 1) {
                    ok = false;
                    break;
                }
            }
            if (!ok)
                break;
        }
    }

    if (f)
        fclose(f);

    if (!ok) 
        error(10002, "Could not read the atomic form factors file \"%s\".\n", ff_file.c_str());
}

void heap_sort(uint n, double *values, uint *indices) {
    for (uint i = 0; i < n; i++)
        indices[i] = i + 1;

    if (n < 2)
        return;

    uint l = n / 2 + 1;
    uint idx = n;

    uint i, j;
    double last_value;
    uint last_idx;

    do {
        if (l > 1) {
            l--;
            last_value = values[l - 1];
            last_idx = l;
        }
        else {
            last_value = values[idx - 1];
            last_idx = indices[idx - 1];
            values[idx - 1] = values[0];
            indices[idx - 1] = indices[0];
            idx--;
            if (idx == 0) {
                values[0] = last_value;
                indices[0] = last_idx;
                return;
            }
        }

        i = l;
        j = 2 * l;

        do {
            if (j > idx)
                break;
            if (j < idx) {
                if (values[j - 1] < values[j])
                    j++;
            }
            if (last_value < values[j - 1]) {
                values[i - 1] = values[j - 1];
                indices[i - 1] = indices[j - 1];
                i = j;
                j = 2 * j;
            }
            else
                j = idx + 1;
        } while (true);

        values[i - 1] = last_value;
        indices[i - 1] = last_idx;
    } while (true);
}

double *get_data(uint flag, uint ne, uint ndat[MXELEMENT], double2 *data[MXELEMENT], double *z_sorted, double *pz_sorted, double2 ge) {
    double *res = (double*)malloc(MXGE * sizeof(double));
    
    for (uint i = 0; i < MXGE; i++)
        res[i] = 0.0F;

    for (uint i = 0; i < ne; i++) {
        uint z = (uint)(z_sorted[i] + 0.5F) - 1;
        uint n = ndat[z];
        double2 *in_dat;
        double eth;

        if (flag == 0) {
            in_dat = (double2*)malloc(n * sizeof(double2));
            for (uint j = 0; j < n; j++)
                in_dat[j] = data[z][j];
        }
        else {
            in_dat = (double2*)malloc((n + 1) * sizeof(double2));
            for (uint j = 0; j < n; j++)
                in_dat[j + 1] = data[z][j];

            if (flag == 1)
                eth = 2.0 * ((double)ELECTRON_REST_MASS_FLOAT);
            else 
                eth = 4.0 * ((double)ELECTRON_REST_MASS_FLOAT);

            n++;
            
            for (uint j = 1; j < n; j++)
                in_dat[j].y -= 3.0 * log(1.0 - eth / exp(in_dat[j].x));

            in_dat[0] = make_double2(log(eth), in_dat[1].y);
        }

        for (uint j = 0; j < MXGE; j++) {
            double gle = ((double)j - ge.x) / ge.y;
            double e = exp(gle);
            double sig;
            
            if ((gle < in_dat[0].x) || (gle >= in_dat[n - 1].x)) {
                if (flag == 0) 
                    error(10003, "Energy %f is outside the available data range of %f to %f.\n", e, exp(in_dat[0].x), exp(in_dat[n - 1].x));
                else {
                    if (gle < in_dat[0].x)
                        sig = 0.0F;
                    else
                        sig = exp(in_dat[n - 1].y);
                }
            }
            else {
                uint k;
                for (k = 0; k < n - 1; k++) {
                    if ((gle >= in_dat[k].x) && (gle < in_dat[k + 1].x))
                        break;
                }
                double p = (gle - in_dat[k].x) / (in_dat[k + 1].x - in_dat[k].x);
                sig = exp(p * in_dat[k + 1].y + (1.0F - p) * in_dat[k].y);
            }
            if ((flag != 0) && (e > eth))
                sig *= (1.0F - eth / e) * (1.0F - eth / e) * (1.0F - eth / e);

            res[j] += pz_sorted[i] * sig;
        }

        free(in_dat);
    }

    return res;
}

double kn_sigma0(double e) {
    float con = 0.1274783851F;

    double ko = e / ((double)ELECTRON_REST_MASS_FLOAT);
    if (ko < 0.01)
        return 8.0 * con / 3.0 * (1.0 - ko * (2.0 - ko * (5.2 -13.3 * ko))) / ((double)ELECTRON_REST_MASS_FLOAT);

    double c1 = 1.0 / (ko * ko);
    double c2 = 1.0 - 2.0 * (1.0 + ko) * c1;
    double c3 = (1.0 + 2.0 * ko) * c1;
    double eps2 = 1.0;
    double eps1 = 1.0 / (1.0 + 2.0 * ko);

    return (c1 * (1.0 / eps1 - 1.0 / eps2) + c2 * log(eps2 / eps1) + eps2 * (c3 + 0.5 * eps2) - eps1 * (c3 + 0.5 * eps1)) / e * ((double)con);
}

uint read_pegs_file(const char *pegs_file, uint nmed, string *media_names, medium_t **media, bool *found) {
    FILE *pegs = fopen(pegs_file, "r");

    if (!pegs) 
        error(10004, "Could not open the PEGS file \"%s\".\n", pegs_file);

    uint media_found = 0;
	int med_idx;

    do {  //medium header search loop
        // read line from pegs file
        char buffer[80];
        fgets(buffer, 80, pegs);
        string line(buffer);

        // here starts a medium definition
        if (line.find(" MEDIUM=") == 0) {
            string name_with_spaces = line.substr(8, 24);
            string name = "";
            // read name up to first space
            for (uint i = 0; i < 24; i++) {
                if (name_with_spaces[i] != ' ')
                    name += name_with_spaces[i];
                else
                    break;
            }

            // see whether this is required medium
            bool required = false;

			for (uint i = 0; i < nmed; i++) {
                if (name == media_names[i]) {
					med_idx = i;
                    required = true;
                    break;
                }
            }
            
            if (!required)
                continue;

            // we have found the i'th required medium
            medium_t *medium = (medium_t*)malloc(sizeof(medium_t));
            medium->name = new string(name);
            medium->ne = 0;

            // read the next line containing the density, number of elements and flags
            fgets(buffer, 80, pegs);
            line = string(buffer);
            uint idx = 0;
            bool ok = true;

            do {
                uint pos_comma = line.find(',', idx);
                if (pos_comma == string::npos)
                    pos_comma = line.length();

                string entry = line.substr(idx, pos_comma - idx);
                idx = pos_comma + 1;

                uint pos_equal = entry.find('=');
                if (pos_equal == string::npos)
                    continue;

                string name_with_spaces = entry.substr(0, pos_equal);
                string name = "";
                for (uint i = 0; i < name_with_spaces.length(); i++) {
                    if (name_with_spaces[i] != ' ')
                        name += name_with_spaces[i];
                }
                string value = entry.substr(pos_equal + 1, entry.length() - pos_equal - 1);

                if (name == "RHO") {
                    double d;
                    if (sscanf(value.c_str(), "%lf", &d) != 1) {
                        ok = false;
                        break;
                    }
                    medium->rho = d;
					h_rho[med_idx]=medium->rho;

                }
                else if (name == "NE") {
                    uint u;
                    if (sscanf(value.c_str(), "%u", &u) != 1) {
                        ok = false;
                        break;
                    }
                    medium->ne = u;
					h_nne[med_idx]=medium->ne;
                }
                else if (name == "IUNRST") {
                    int i;
                    if (sscanf(value.c_str(), "%d", &i) != 1) {
                        ok = false;
                        break;
                    }
                    //medium->iunrst = i;
					//h_iunrst[med_idx]=medium->iunrst;
                }
                else if (name == "EPSTFL") {
                    int i;
                    if (sscanf(value.c_str(), "%d", &i) != 1) {
                        ok = false;
                        break;
                    }
                    //medium->epstfl = i;
					//h_epstfl[med_idx]=medium->epstfl;
                }
                else if (name == "IAPRIM") {
                    int i;
                    if (sscanf(value.c_str(), "%d", &i) != 1) {
                        ok = false;
                        break;
                    }
                    //medium->iaprim = i;
					//h_iaprim[med_idx]=medium->iaprim;
                }

            } while (idx < line.length());

            if (!ok)
                continue;

            // read elements
            medium->elements = (element_t*)malloc(medium->ne * sizeof(element_t));
            for (uint ie = 0; ie < medium->ne; ie++) {
                element_t element;
                
                fgets(buffer, 80, pegs);
                line = string(buffer);
                idx = 0;

                do {
                    uint pos_comma = line.find(',', idx);
                    if (pos_comma == string::npos)
                        pos_comma = line.length();

                    string entry = line.substr(idx, pos_comma - idx);
                    idx = pos_comma + 1;

                    uint pos_equal = entry.find('=');
                    if (pos_equal == string::npos)
                        continue;

                    string name_with_spaces = entry.substr(0, pos_equal);
                    string name = "";
                    for (uint i = 0; i < name_with_spaces.length(); i++) {
                        if (name_with_spaces[i] != ' ')
                            name += name_with_spaces[i];
                    }
                    string value = entry.substr(pos_equal + 1, entry.length() - pos_equal - 1);

                    if (name == "ASYM") {
                        if (value.length() < 2) {
                            ok = false;
                            break;
                        }
                        element.symbol[0] = value[0];
                        element.symbol[1] = value[1];
                    }
                    else if (name == "Z") {
                        double d;
                        if (sscanf(value.c_str(), "%lf", &d) != 1) {
                            ok = false;
                            break;
                        }
                        element.z = d;
						h_zelem[MXMED*ie+med_idx]=element.z;
                    }
                    else if (name == "A") {
                        double d;
                        if (sscanf(value.c_str(), "%lf", &d) != 1) {
                            ok = false;
                            break;
                        }
                        element.wa = d;
						h_wa[MXMED*ie+med_idx]=element.wa;
                    }
                    else if (name == "PZ") {
                        double d;
                        if (sscanf(value.c_str(), "%lf", &d) != 1) {
                            ok = false;
                            break;
                        }
                        element.pz = d;
						h_pz[MXMED*ie+med_idx]=element.pz;
                    }
                    else if (name == "RHOZ") {
                        double d;
                        if (sscanf(value.c_str(), "%lf", &d) != 1) {
                            ok = false;
                            break;
                        }
                        element.rhoz = d;
						h_rhoz[MXMED*ie+med_idx]=element.rhoz;
                    }

                } while (idx < line.length());

                if (!ok)
                    break;
                
                medium->elements[ie] = element;
            }

            if (!ok)
                continue;

			// media and thresh
            // read next line that contines rlc, ae, ap, ue, up
			fscanf(pegs, "%lf %lf %lf %lf %lf", &medium->rlc, &medium->ae, &medium->ap, &medium->ue, &medium->up);
			h_rlc[med_idx]=medium->rlc;
			h_ae[med_idx]=medium->ae;
			h_ap[med_idx]=medium->ap;
			h_ue[med_idx]=medium->ue;
			h_up[med_idx]=medium->up;
            medium->te=medium->ae-((double)ELECTRON_REST_MASS_FLOAT);
			medium->thmoll=medium->te*2+((double)ELECTRON_REST_MASS_FLOAT);
			h_te[med_idx]=h_ae[med_idx]-((double)ELECTRON_REST_MASS_FLOAT);
			h_thmoll[med_idx]=h_te[med_idx]*2+((double)ELECTRON_REST_MASS_FLOAT);
			
			// actual array sizes from pegs
			int dummy,nge,neke;
			fscanf(pegs,"%d %d %d %d %d %d %d %d %d",&medium->msge,
				&medium->mge,&medium->mseke,&medium->meke,
				&medium->mleke,&medium->mcmfp,
				&medium->mrange,
				&medium->irayl, &dummy);
			//h_msge[med_idx]=medium->msge;
			//h_mge[med_idx]=medium->mge;
			//h_mseke[med_idx]=medium->mseke;
			h_meke[med_idx]=medium->meke;
			//h_mleke[med_idx]=medium->mleke;
			//h_mcmfp[med_idx]=medium->mcmfp;
			//h_mrange[med_idx]=medium->mrange;
			//irayl=medium->irayl;  
			
			//nsge=medium->msge;
			nge=medium->mge;
			//nseke=medium->mseke;
			neke=medium->meke;
			//nleke=medium->mleke;
			//ncmfp=medium->mcmfp;
			//nrange=medium->mrange;
					
			// brempr
			for (int ia=0;ia<6;ia++)
			{
				fscanf (pegs,"%lf %lf %lf %lf %lf %lf",&(medium->dl1[ia]),&medium->dl2[ia],&medium->dl3[ia],&medium->dl4[ia],&medium->dl5[ia],&medium->dl6[ia]);
				h_dl1[8*med_idx+ia]=medium->dl1[ia];
				h_dl2[8*med_idx+ia]=medium->dl2[ia];
				h_dl3[8*med_idx+ia]=medium->dl3[ia];
				h_dl4[8*med_idx+ia]=medium->dl4[ia];
				h_dl5[8*med_idx+ia]=medium->dl5[ia];
				h_dl6[8*med_idx+ia]=medium->dl6[ia];
			}
			fscanf(pegs,"%lf",&medium->delcm);
			h_delcm[med_idx]=medium->delcm;
			for(int ib=0;ib<2;ib++)
			{
				fscanf (pegs,"%lf %lf %lf",&medium->alphi[ib],&medium->bpar[ib],&medium->delpos[ib]);
				//h_alphi[2*med_idx+ib]=medium->alphi[ib];
				h_bpar[2*med_idx+ib]=medium->bpar[ib];
				//h_delpos[2*med_idx+ib]=medium->delpos[ib];
			}

			// elecin
			fscanf (pegs,"%lf %lf %lf %lf",&medium->xr0,&medium->teff0,&medium->blcc,&medium->xcc);
			//h_xr0[med_idx]=medium->xr0;
			//h_teff0[med_idx]=medium->teff0;
			h_blcc[med_idx]=medium->blcc;
			h_xcc[med_idx]=medium->xcc;
			fscanf (pegs,"%lf %lf", &medium->eke0,&medium->eke1); 
			h_eke01[med_idx].x=medium->eke0;
			h_eke01[med_idx].y=medium->eke1;
			for (int ic=0;ic<neke;ic++)
			{
				fscanf (pegs,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&medium->esig0[ic],&medium->esig1[ic],&medium->psig0[ic],&medium->psig1[ic],&medium->ededx0[ic],&medium->ededx1[ic],&medium->pdedx0[ic],&medium->pdedx1[ic],&medium->ebr10[ic],&medium->ebr11[ic],&medium->pbr10[ic],&medium->pbr11[ic],&medium->pbr20[ic],&medium->pbr21[ic],&medium->tmxs0[ic],&medium->tmxs1[ic]);
				//h_esig[MXEKE*med_idx+ic].x=medium->esig0[ic];
				h_sig[MXEKE*med_idx+ic].x=medium->esig0[ic];
				//h_esig[MXEKE*med_idx+ic].y=medium->esig1[ic];
				h_sig[MXEKE*med_idx+ic].y=medium->esig1[ic];
				//h_psig[MXEKE*med_idx+ic].x=medium->psig0[ic];
				h_sig[MXEKE*MXMED+MXEKE*med_idx+ic].x=medium->psig0[ic];
				//h_psig[MXEKE*med_idx+ic].y=medium->psig1[ic];
				h_sig[MXEKE*MXMED+MXEKE*med_idx+ic].y=medium->psig1[ic];
				//h_ededx[MXEKE*med_idx+ic].x=medium->ededx0[ic];
				h_dedx[MXEKE*med_idx+ic].x=medium->ededx0[ic];
				//h_ededx[MXEKE*med_idx+ic].y=medium->ededx1[ic];
				h_dedx[MXEKE*med_idx+ic].y=medium->ededx1[ic];
				//h_pdedx[MXEKE*med_idx+ic].x=medium->pdedx0[ic];
				h_dedx[MXEKE*MXMED+MXEKE*med_idx+ic].x=medium->pdedx0[ic];
				//h_pdedx[MXEKE*med_idx+ic].y=medium->pdedx1[ic];
				h_dedx[MXEKE*MXMED+MXEKE*med_idx+ic].y=medium->pdedx1[ic];
				h_ebr1[MXEKE*med_idx+ic].x=medium->ebr10[ic];
				h_ebr1[MXEKE*med_idx+ic].y=medium->ebr11[ic];
				h_pbr1[MXEKE*med_idx+ic].x=medium->pbr10[ic];
				h_pbr1[MXEKE*med_idx+ic].y=medium->pbr11[ic];
				h_pbr2[MXEKE*med_idx+ic].x=medium->pbr20[ic];
				h_pbr2[MXEKE*med_idx+ic].y=medium->pbr21[ic];
				h_tmxs[MXEKE*med_idx+ic].x=medium->tmxs0[ic];
				h_tmxs[MXEKE*med_idx+ic].y=medium->tmxs1[ic];
			}

			// photin
			fscanf (pegs,"%lf %lf %lf",&medium->ebinda,&medium->ge0,&medium->ge1);
			//h_ebinda[med_idx]=medium->ebinda;
			//h_ge[med_idx].x=medium->ge0;
			//h_ge[med_idx].y=medium->ge1;
			for(int i0=0;i0<nge;i0++)
			{
				fscanf(pegs,"%lf %lf %lf %lf %lf %lf",&medium->gmfp0[i0],&medium->gmfp1[i0],&medium->gbr10[i0],&medium->gbr11[i0],&medium->gbr20[i0],&medium->gbr21[i0]);
				//h_gmfp[MXGE*med_idx+i0].x=medium->gmfp0[i0];
				//h_gmfp[MXGE*med_idx+i0].y=medium->gmfp1[i0];
				//h_gbr1[MXGE*med_idx+i0].x=medium->gbr10[i0];
				//h_gbr1[MXGE*med_idx+i0].y=medium->gbr11[i0];
				//h_gbr2[MXGE*med_idx+i0].x=medium->gbr20[i0];
				//h_gbr2[MXGE*med_idx+i0].y=medium->gbr21[i0];
			}

			// photin (continued) -- optional rayleigh scattering input
			int ngrim;
			if (medium->irayl==1)    //Rayleigh data available for this medium  in PEGS4 data set.'
			{
				fscanf (pegs,"%d",&medium->ngr);
				//h_ngr[med_idx]=medium->ngr;
				ngrim=medium->ngr;
				fscanf (pegs,"%lf %lf",&medium->rco0,&medium->rco1);
				//h_rco[med_idx].x=medium->rco0;  
				//h_rco[med_idx].y=medium->rco1;
				for (int i1=0;i1<ngrim;i1++)
				{
					fscanf (pegs,"%lf %lf",&medium->rsct0[i1],&medium->rsct1[i1]);
					//h_rsct[MXRAYFF*med_idx+i1].x=medium->rsct0[i1];
					//h_rsct[MXRAYFF*med_idx+i1].y=medium->rsct1[i1];
				}
				for (int i2=0;i2<nge;i2++)
				{
					fscanf (pegs,"%lf %lf",&medium->cohe0[i2],&medium->cohe1[i2]);
					//h_cohe[MXGE*med_idx+i2].x=medium->cohe0[i2];
					//h_cohe[MXGE*med_idx+i2].y=medium->cohe1[i2];
				}
			}

            // save the medium and mark it found
            found[med_idx] = true;
            media[med_idx] = medium;
            media_found++;
        
		} //if (line.find(" MEDIUM=") == 0)
    } while ((media_found < nmed) && (!feof(pegs)));

	// We now have data for all media requested. Now do distance unit change.
	// Data from pegs is in units of radiation lengths.
	// EGS is run in uints of 'dunit' centimeters if dunit>0, or in uints of rlc[-dunit] centimeters if dunit<0.
	// That is, a negative dunit means unit is to be the radiation length of the medium whose index is -dunit.
	double pznorm;
	double dunitr;
	//int id;
	//int c1,c2;
	dunitr=h_dunit;
	if(h_dunit<0.0f)
	{
		/*
		if(int(-h_dunit)<=10)
			c1=int(-h_dunit);
		else c1=10;
	    if(c1>=1)
			c2=c1;
		else c2=1;
		id=c2;
		*/
		int id=max(1,min(MXMED,int(-h_dunit)));  //corrected by Wu, 20140703
		h_dunit=h_rlc[id-1];
		if(h_dunit!=1.0f)
		{
			logstr("DUNIT REQUESTED&USED ARE: %f, %f\n" ,dunitr, h_dunit);
		}
	}
	//int im,i5,i6;  // Sisen 2013 3 21
    double dfact,dfacti;
	for (int im=0;im<nmed;im++)
	{
		//dfact=media[im]->rlc/dunit;
		dfact=h_rlc[im]/h_dunit;	
		dfacti=1.0/dfact;
		h_smaxir[im]=MAX_SMAX;  //corrected by Wu, 20130610
		for (int i5=0;i5<media[im]->meke;i5++)
		{
			/*
			media[im]->esig0[i5]=media[im]->esig0[i5]*dfacti;
			media[im]->esig1[i5]=media[im]->esig1[i5]*dfacti;
			media[im]->psig0[i5]=media[im]->psig0[i5]*dfacti;
			media[im]->psig1[i5]=media[im]->psig1[i5]*dfacti;
			media[im]->ededx0[i5]=media[im]->ededx0[i5]*dfacti;
			media[im]->ededx1[i5]=media[im]->ededx1[i5]*dfacti;
			media[im]->pdedx0[i5]=media[im]->pdedx0[i5]*dfacti;
			media[im]->pdedx1[i5]=media[im]->pdedx1[i5]*dfacti;
			media[im]->tmxs0[i5]=media[im]->tmxs0[i5]*dfacti;
			media[im]->tmxs1[i5]=media[im]->tmxs1[i5]*dfacti;
			*/
			//h_esig[MXEKE*im+i5].x*=dfacti;
			h_sig[MXEKE*im+i5].x*=dfacti;
			//h_esig[MXEKE*im+i5].y*=dfacti;
			h_sig[MXEKE*im+i5].y*=dfacti;
			//h_psig[MXEKE*im+i5].x*=dfacti;
			h_sig[MXEKE*MXMED+MXEKE*im+i5].x*=dfacti;
			//h_psig[MXEKE*im+i5].y*=dfacti;
			h_sig[MXEKE*MXMED+MXEKE*im+i5].y*=dfacti;
			//h_ededx[MXEKE*im+i5].x*=dfacti;
			h_dedx[MXEKE*im+i5].x*=dfacti;
			//h_ededx[MXEKE*im+i5].y*=dfacti;
			h_dedx[MXEKE*im+i5].y*=dfacti;
			//h_pdedx[MXEKE*im+i5].x*=dfacti;
			h_dedx[MXEKE*MXMED+MXEKE*im+i5].x*=dfacti;
			//h_pdedx[MXEKE*im+i5].y*=dfacti;
			h_dedx[MXEKE*MXMED+MXEKE*im+i5].y*=dfacti;
			h_tmxs[MXEKE*im+i5].x*=dfacti;
			h_tmxs[MXEKE*im+i5].y*=dfacti;
		}
		/*
		media[im]->teff0=media[im]->teff0*dfact;
		media[im]->blcc=media[im]->blcc*dfacti;
		media[im]->xcc=media[im]->xcc*sqrtf(dfacti);
		media[im]->rldu=media[im]->rlc/DUNIT;
		*/
		//h_teff0[im]*=dfact;
		h_blcc[im]*=dfacti;
		h_xcc[im]*=sqrtf(dfacti);
		//h_rldu[im]=h_rldu[im]/h_dunit;
		//h_rldu[im]=h_rlc[im]/h_dunit;  //corrected by Wu, 20140703
		/*
		for (int i6=0;i6<media[im]->mge;i6++)
		{
			//media[im]->gmfp0[i6]=media[im]->gmfp0[i6]*dfact;
			//media[im]->gmfp1[i6]=media[im]->gmfp1[i6]*dfact;
			h_gmfp[MXGE*im+i6].x*=dfact;
			h_gmfp[MXGE*im+i6].y*=dfact;
		} // the line 2599 to 2598 is not translated£¡£¡
		*/
	}// im
	// scale vacdst. undo previous scale, then do new.
	//float dunito=1.0;  //corrected by Wu, 20130610
	//h_vacdst=h_vacdst*dunito/h_dunit;
	//dunito=h_dunit;
	/*
	//the following is to correct the region data to be consistant with physics data, need to do it in init_region.
	for(int i5=0;i5<MXREG;i5++)
	{
		int md=h_med[i5];
		if((md>=1)&&(md<=nmed))
		{
			if(h_ecut[i5]>=h_ae[md-1]) h_ecut[i5]=h_ecut[i5];
			else h_ecut[i5]=h_ae[md-1];
			if(h_pcut[i5]<=h_ap[md-1]) h_pcut[i5]=h_ap[md-1];
			else h_pcut[i5]=h_pcut[i5];
			if(h_rhor[i5]==0.0f)
				h_rhor[i5]=h_rho[md-1];
		}
	}
	*/
	if(h_ibrdst==1||h_iprdst>0)
	{
		for(int im=0;im<nmed;im++)
		{
			h_zbrang[im]=0.0f;
			pznorm=0.0f;
			for (int ie=0;ie<h_nne[im];ie++)
			{
				h_zbrang[im]+=h_pz[MXMED*ie+im]*h_zelem[MXMED*ie+im]*(h_zelem[MXMED*ie+im]+1.0f);
				pznorm+=h_pz[MXMED*ie+im];		
			}
			h_zbrang[im]=(8.116224e-05)*pow((h_zbrang[im]/pznorm),1.0/3.0);
			h_lzbrang[im]=-log(h_zbrang[im]);
		}
	}
	/*
	if(h_iprdst>0)  // this is merged into above, by Wu, 20140703
	{
		for(int im=0;im<nmed;im++)
		{
			h_zbrang[im]=0.0f;
			pznorm=0.0f;
			for(int ie=0;ie<h_nne[im];ie++)
			{
				h_zbrang[im]+=h_pz[MXMED*ie+im]*h_zelem[MXMED*ie+im]*(h_zelem[MXMED*ie+im]+1.0f);
				pznorm+=h_pz[MXMED*ie+im];
			}
			h_zbrang[im]=(8.116224e-05)*pow((h_zbrang[im]/pznorm),1.0/3.0);
			h_lzbrang[im]=-log(h_zbrang[im]);
		}
	}
	*/

    fclose(pegs);

    return media_found;
}

void init_user_photon(uint nmed, medium_t **media, const char *photon_xsections, const char *comp_xsections) {
    // read photon cross sections data
    uint photo_ndat[MXELEMENT];
    uint rayleigh_ndat[MXELEMENT];
    uint pair_ndat[MXELEMENT];
    uint triplet_ndat[MXELEMENT];

    double2 *photo_xsec_data[MXELEMENT];
    double2 *rayleigh_xsec_data[MXELEMENT];
    double2 *pair_xsec_data[MXELEMENT];
    double2 *triplet_xsec_data[MXELEMENT];

    string prefix = string(data_dir) + string(photon_xsections);
    read_xsec_data((prefix + "_photo.data").c_str(), photo_ndat, photo_xsec_data);
    read_xsec_data((prefix + "_rayleigh.data").c_str(), rayleigh_ndat, rayleigh_xsec_data);
    read_xsec_data((prefix + "_pair.data").c_str(), pair_ndat, pair_xsec_data);
    read_xsec_data((prefix + "_triplet.data").c_str(), triplet_ndat, triplet_xsec_data);

    for (uint i = 0; i < MXELEMENT; i++) {
        uint n = photo_ndat[i];
        uint k = 0;
        for (uint j = n - 1; j > 1; j--) {
            if (photo_xsec_data[i][j].x - photo_xsec_data[i][j - 1].x < 1.0E-5F)
                binding_energies2D[k++][i] = exp(photo_xsec_data[i][j].x);
            if (k >= 3)
                break;
        }
    }

    // allocate memory
    double2 *h_ge = (double2*)malloc(nmed * sizeof(double2));
    double2 *h_gmfp = (double2*)malloc(nmed * MXGE * sizeof(double2));
    double2 *h_gbr1 = (double2*)malloc(nmed * MXGE * sizeof(double2));
    double2 *h_gbr2 = (double2*)malloc(nmed * MXGE * sizeof(double2));
    double2 *h_cohe = (double2*)malloc(nmed * MXGE * sizeof(double2));

    for (uint i = 0; i < nmed; i++) {
        medium_t m = *media[i];
        h_ge[i].y = (double)(MXGE-1 ) / log(m.up / m.ap);

        // indexing starts at 0 and not at 1 as in FORTRAN, i.e. subtract 1
        h_ge[i].x = -h_ge[i].y * log(m.ap);
        
        double sumA = 0.0F;
        double sumZ = 0.0F;

        double *z_sorted = (double*)malloc(m.ne * sizeof(double));
        for (uint j = 0; j < m.ne; j++) {
            z_sorted[j] = m.elements[j].z;
            sumA += m.elements[j].pz * m.elements[j].wa;
            sumZ += m.elements[j].pz * m.elements[j].z;
        }

        //double con1 = sumZ * m.rho / (sumA * 1.6605655F); // apparently not used
        double con2 = m.rho / (sumA * ((double)1.6605655F));

        uint *sorted = (uint*)malloc(m.ne * sizeof(uint));
        heap_sort(m.ne, z_sorted, sorted);
        // indexing starts at 0 and not at 1 as in FORTRAN, i.e. subtract 1
        for (uint j = 0; j < m.ne; j++)
            sorted[j] -= 1;

        double *pz_sorted = (double*)malloc(m.ne * sizeof(double));
        for (uint j = 0; j < m.ne; j++)
            pz_sorted[j] = m.elements[sorted[j]].pz;

        double *sig_photo = get_data(0, m.ne, photo_ndat, photo_xsec_data, z_sorted, pz_sorted, h_ge[i]);
        double *sig_rayleigh = get_data(0, m.ne, rayleigh_ndat, rayleigh_xsec_data, z_sorted, pz_sorted, h_ge[i]);
        double *sig_pair = get_data(1, m.ne, pair_ndat, pair_xsec_data, z_sorted, pz_sorted, h_ge[i]);
        double *sig_triplet = get_data(2, m.ne, triplet_ndat, triplet_xsec_data, z_sorted, pz_sorted, h_ge[i]);

        // do bound compton here

        double gle;
        double gmfp;
        double gbr1;
        double gbr2;
        double cohe;

        double gmfp_old = 0.0F;
        double gbr1_old = 0.0F;
        double gbr2_old = 0.0F;
        double cohe_old = 0.0F;

        for (uint j = 0; j < MXGE; j++) {
            gle = ((double)j - h_ge[i].x) / h_ge[i].y;
            double e = exp(gle);
            double sig_kn = sumZ * kn_sigma0(e);

            // do bound compton here

            double sig_p = sig_pair[j] + sig_triplet[j];
            double sigma = sig_kn + sig_p + sig_photo[j];
            gmfp = 1.0 / (sigma * con2);
            gbr1 = sig_p / sigma;
            gbr2 = gbr1 + sig_kn / sigma;
            cohe = sigma / (sig_rayleigh[j] + sigma);

            if (j > 0) {
                uint idx = i * MXGE + j - 1;
                h_gmfp[idx].y = (gmfp - gmfp_old) * h_ge[i].y;
                h_gmfp[idx].x = gmfp - h_gmfp[idx].y * gle;
                h_gbr1[idx].y = (gbr1 - gbr1_old) * h_ge[i].y;
                h_gbr1[idx].x = gbr1 - h_gbr1[idx].y * gle;
                h_gbr2[idx].y = (gbr2 - gbr2_old) * h_ge[i].y;
                h_gbr2[idx].x = gbr2 - h_gbr2[idx].y * gle;
                h_cohe[idx].y = (cohe - cohe_old) * h_ge[i].y;
                h_cohe[idx].x = cohe - h_cohe[idx].y * gle;
            }

            gmfp_old = gmfp;
            gbr1_old = gbr1;
            gbr2_old = gbr2;
            cohe_old = cohe;
        }

        uint idx = i * MXGE + MXGE - 1;

        h_gmfp[idx].y = h_gmfp[idx - 1].y;
        h_gmfp[idx].x = gmfp - h_gmfp[idx].y * gle;
        h_gbr1[idx].y = h_gbr1[idx - 1].y;
        h_gbr1[idx].x = gbr1 - h_gbr1[idx].y * gle;
        h_gbr2[idx].y = h_gbr2[idx - 1].y;
        h_gbr2[idx].x = gbr2 - h_gbr2[idx].y * gle;
        h_cohe[idx].y = h_cohe[idx - 1].y;
        h_cohe[idx].x = cohe - h_cohe[idx].y * gle;

        free(z_sorted);
        free(sorted);
        free(pz_sorted);

        free(sig_photo);
        free(sig_rayleigh);
        free(sig_pair);
        free(sig_triplet);
    }

    for (uint i = 0; i < MXELEMENT; i++) {
        free(photo_xsec_data[i]);
        free(rayleigh_xsec_data[i]);
        free(pair_xsec_data[i]);
        free(triplet_xsec_data[i]);
    }

    // convert to floats
    float2 *h_ge_f = (float2*)malloc(nmed * sizeof(float2));
    float2 *h_gmfp_f = (float2*)malloc(nmed * MXGE * sizeof(float2));
    float2 *h_gbr1_f = (float2*)malloc(nmed * MXGE * sizeof(float2));
    float2 *h_gbr2_f = (float2*)malloc(nmed * MXGE * sizeof(float2));
    float2 *h_cohe_f = (float2*)malloc(nmed * MXGE * sizeof(float2));

    for (uint i = 0; i < nmed; i++)
        h_ge_f[i] = make_float2((float)h_ge[i].x, (float)h_ge[i].y);

    for (uint i = 0; i < nmed * MXGE; i++) {
        h_gmfp_f[i] = make_float2((float)h_gmfp[i].x, (float)h_gmfp[i].y);
        h_gbr1_f[i] = make_float2((float)h_gbr1[i].x, (float)h_gbr1[i].y);
        h_gbr2_f[i] = make_float2((float)h_gbr2[i].x, (float)h_gbr2[i].y);
        h_cohe_f[i] = make_float2((float)h_cohe[i].x, (float)h_cohe[i].y);
    }

	/*
	FILE *iSum=fopen("sum_for_init_user_photon.txt","wt");
	FILE *iDump;
	double sum;
	fprintf(iSum,"\n%s\n","now the photin...");
	DUMP_1D_REAL(photin,h_ebinda,MXMED); 
	DUMP_1D_FLOAT2(photin,h_ge,MXMED); 
	DUMP_1D_FLOAT2(photin,h_gmfp,MXGE*MXMED);
	DUMP_1D_FLOAT2(photin,h_gbr1,MXGE*MXMED);
	DUMP_1D_FLOAT2(photin,h_gbr2,MXGE*MXMED);
	DUMP_1D_FLOAT2(photin,h_rco,MXMED); 
	DUMP_1D_FLOAT2(photin,h_rsct,MXRAYFF*MXMED);
	DUMP_1D_FLOAT2(photin,h_cohe,MXGE*MXMED);
	DUMP_1D_INT(photin, h_mpgem,MXSGE*MXMED);
	DUMP_1D_INT(photin, h_ngr,MXMED);
	fclose(iSum);
	*/

    // free host double memory
    free(h_ge);
    free(h_gmfp);
    free(h_gbr1);
    free(h_gbr2);
    free(h_cohe);

	d_ge   = (float2**)malloc(GPUNo*sizeof(float2*));
	d_gmfp = (float2**)malloc(GPUNo*sizeof(float2*));
	d_gbr1 = (float2**)malloc(GPUNo*sizeof(float2*));
	d_gbr2 = (float2**)malloc(GPUNo*sizeof(float2*));
	d_cohe = (float2**)malloc(GPUNo*sizeof(float2*));

	for(int GPUId=0; GPUId<GPUNo; GPUId++) {

#ifdef USE_MULTIPLE_GPU
		cudaSetDevice(GPUId); ce(58006);
#endif

		// allocate device memory
		cudaMalloc(&d_ge[GPUId],   nmed * sizeof(float2)); ce(10005);
		cudaMalloc(&d_gmfp[GPUId], nmed * MXGE * sizeof(float2)); ce(10006);
		cudaMalloc(&d_gbr1[GPUId], nmed * MXGE * sizeof(float2)); ce(10007);
		cudaMalloc(&d_gbr2[GPUId], nmed * MXGE * sizeof(float2)); ce(10008);
		cudaMalloc(&d_cohe[GPUId], nmed * MXGE * sizeof(float2)); ce(10001);

		// copy data to device
		cudaMemcpy(d_ge[GPUId],   h_ge_f,   nmed * sizeof(float2), cudaMemcpyHostToDevice); ce(10010);
		cudaMemcpy(d_gmfp[GPUId], h_gmfp_f, nmed * MXGE * sizeof(float2), cudaMemcpyHostToDevice); ce(10011);
		cudaMemcpy(d_gbr1[GPUId], h_gbr1_f, nmed * MXGE * sizeof(float2), cudaMemcpyHostToDevice); ce(10012);
		cudaMemcpy(d_gbr2[GPUId], h_gbr2_f, nmed * MXGE * sizeof(float2), cudaMemcpyHostToDevice); ce(10013);
		cudaMemcpy(d_cohe[GPUId], h_cohe_f, nmed * MXGE * sizeof(float2), cudaMemcpyHostToDevice); ce(10014);

		cudaMemcpyToSymbol(ge,   &d_ge[GPUId],   sizeof(float2*)); ce(10015);
		cudaMemcpyToSymbol(gmfp, &d_gmfp[GPUId], sizeof(float2*)); ce(10016);
		cudaMemcpyToSymbol(gbr1, &d_gbr1[GPUId], sizeof(float2*)); ce(10017);
		cudaMemcpyToSymbol(gbr2, &d_gbr2[GPUId], sizeof(float2*)); ce(10018);
		cudaMemcpyToSymbol(cohe, &d_cohe[GPUId], sizeof(float2*)); ce(10019);

	}

    // free host single memory
    free(h_ge_f);
    free(h_gmfp_f);
    free(h_gbr1_f);
    free(h_gbr2_f);
    free(h_cohe_f);
}

void init_rayleigh_data(uint nmed, medium_t **media) {
    double xval[MXRAYFF];
    double aff[MXELEMENT][MXRAYFF];

    read_ff_data(atomic_ff_file, xval, aff);

    // allocate memory
    double *h_ff = (double*)malloc(MXRAYFF * nmed * sizeof(double));
    double *h_xgrid = (double*)malloc(MXRAYFF * nmed * sizeof(double));
    double *h_fcum = (double*)malloc(MXRAYFF * nmed * sizeof(double)); 
    double *h_b_array = (double*)malloc(MXRAYFF * nmed * sizeof(double)); 
    double *h_c_array = (double*)malloc(MXRAYFF * nmed * sizeof(double)); 
    uint *h_i_array = (uint*)malloc(RAYCDFSIZE * nmed * sizeof(uint)); 
    double *h_pe_array = (double*)malloc(MXGE * nmed * sizeof(double)); 
    double2 *h_pmax = (double2*)malloc(MXGE * nmed * sizeof(double2));

    for (uint i = 0; i < nmed; i++) {
        // calculate form factor using independent atom model
        for (uint j = 0; j < MXRAYFF; j++) {
            double ff_val = 0.0;
            h_xgrid[i * MXRAYFF + j] = xval[j];

            for (uint k = 0; k < media[i]->ne; k ++) {
                uint z = (uint)media[i]->elements[k].z - 1;
                ff_val += media[i]->elements[k].pz * aff[z][j] * aff[z][j];
            }

            h_ff[i * MXRAYFF + j] = sqrtf(ff_val);
        }

        // to avoid log(0)
        /*if (*((ulong*)&h_xgrid[i * MXRAYFF]) == 0) {
            ulong zero = 1;
            h_xgrid[i * MXRAYFF] = *((double*)&zero);
        }*/
        if (h_xgrid[i * MXRAYFF] < 1E-6)
            h_xgrid[i * MXRAYFF] = ((double)0.0001F);

        // calculate rayleigh data (subroutine prepare_rayleigh_data)

        double2 ge;
        ge.y = (double)(MXGE - 1) / log(media[i]->up / media[i]->ap);

        // indexing starts at 0 and not at 1 as in FORTRAN, i.e. subtract 1
        ge.x = -ge.y * log(media[i]->ap);

        double emin = exp(-ge.x / ge.y);
        double emax = exp(((double)MXGE - 1.0 - ge.x) / ge.y);

        // to avoid log (0)
        for (uint j = 0; j < MXRAYFF; j++) {
            if (*((ulong*)&h_ff[i * MXRAYFF + j]) == 0) {
                ulong zero = 1;
                h_ff[i * MXRAYFF + j] = *((double*)&zero);
            }
        }

        /**********************************************************
         * Calculate the cumulative distribution
         *********************************************************/
        double sum0 = 0.0;
        h_fcum[i * MXRAYFF] = 0.0;
        
        for (uint j = 0; j < MXRAYFF - 1; j++) {
            double b = log(h_ff[i * MXRAYFF + j + 1] / h_ff[i * MXRAYFF + j]) / log(h_xgrid[i * MXRAYFF + j + 1] / h_xgrid[i * MXRAYFF + j]);
            h_b_array[i * MXRAYFF + j] = b;
            double x1 = h_xgrid[i * MXRAYFF + j];
            double x2 = h_xgrid[i * MXRAYFF + j + 1];
            double pow_x1 = pow(x1, 2.0 * b);
            double pow_x2 = pow(x2, 2.0 * b);
            sum0 += h_ff[i * MXRAYFF + j] * h_ff[i * MXRAYFF + j] * (x2 * x2 * pow_x2 - x1 * x1 * pow_x1) / ((1.0 + b) * pow_x1);
            h_fcum[i * MXRAYFF + j + 1] = sum0;
        }

        h_b_array[i * MXRAYFF + MXRAYFF - 1] = 0.0;

        /*************************************************************
         * Now the maximum cumulative probability as a function of
         * incident photon energy. We have xmax = 2*E*20.60744/m, so
         * pe_array(E) = fcum(xmax)
         **************************************************************/
        //double dle = log(media[i]->up / media[i]->ap) / ((double)MXGE - 1.0);
        double dle = log(emax / emin) / ((double)MXGE - 1.0);
        uint idx = 1;

        for (uint j = 1; j <= MXGE; j++) {
            //double e = media[i]->ap * exp(dle * ((double)j - 1.0));
            double e = emin * exp(dle * ((double)j - 1.0));
            double xmax = 20.607544 * 2.0 * e / ELECTRON_REST_MASS_DOUBLE;
            uint k = 1;
            for (k = 1; k <= MXRAYFF - 1; k++) {
                if ((xmax >= h_xgrid[i * MXRAYFF + k - 1]) && (xmax < h_xgrid[i * MXRAYFF + k]))
                    break;
            }
            idx = k;
            double b = h_b_array[i * MXRAYFF + idx - 1];
            double x1 = h_xgrid[i * MXRAYFF + idx - 1];
            double x2 = xmax;
            double pow_x1 = pow(x1, 2.0 * b);
            double pow_x2 = pow(x2, 2.0 * b);
            h_pe_array[i * MXGE + j - 1] = h_fcum[i * MXRAYFF + idx - 1] + h_ff[i * MXRAYFF + idx - 1] * h_ff[i * MXRAYFF + idx - 1] * 
                (x2 * x2 * pow_x2 - x1 * x1 * pow_x1) / ((1.0 + b) * pow_x1);
        }

        h_i_array[i * RAYCDFSIZE + RAYCDFSIZE - 1] = idx;

        /***********************************************************************
         * Now renormalize data so that pe_array(emax)=1
         * Note that we make pe_array(j) slightly larger so that fcum(xmax) is
         * never underestimated when interpolating
         ***********************************************************************/
        double anorm = 1.0 / sqrtf(h_pe_array[i * MXGE + MXGE - 1]);
        double anorm1 = 1.005 / h_pe_array[i * MXGE + MXGE - 1];
        double anorm2 = 1.0 / h_pe_array[i * MXGE + MXGE - 1];

        for (uint j = 0; j < MXGE; j++) {
            h_pe_array[i * MXGE + j] *= anorm1;
            if (h_pe_array[i * MXGE + j] > 1.0)
                h_pe_array[i * MXGE + j] = 1.0;
        }
        
        for (uint j = 0; j < MXRAYFF; j++) {
            h_ff[i * MXRAYFF + j] *= anorm;
            h_fcum[i * MXRAYFF + j] *= anorm2;
            h_c_array[i * MXRAYFF + j] = (1.0 + h_b_array[i * MXRAYFF + j]) / 
                ((h_xgrid[i * MXRAYFF + j] * h_ff[i * MXRAYFF + j]) * (h_xgrid[i * MXRAYFF + j] * h_ff[i * MXRAYFF + j]));
        }

        /***********************************************************************
         * Now prepare uniform cumulative bins
         ***********************************************************************/
        double dw = 1.0 / ((double)RAYCDFSIZE - 1.0);
        double xold = h_xgrid[i * MXRAYFF + 0];
        uint ibin = 1;
        double b = h_b_array[i * MXRAYFF + 0];
        double pow_x1 = pow(h_xgrid[i * MXRAYFF + 0], 2.0 * b);
        h_i_array[i * MXRAYFF + 0] = 1;

        for (uint j = 2; j <= RAYCDFSIZE - 1; j++) {
            double w = dw;
            do {
                double x1 = xold;
                double x2 = h_xgrid[i * MXRAYFF + ibin];
                double t = x1 * x1 * pow(x1, 2.0 * b);
                double pow_x2 = pow(x2, 2.0 * b);
                double aux = h_ff[i * MXRAYFF + ibin - 1] * h_ff[i * MXRAYFF + ibin - 1] * (x2 * x2 * pow_x2 - t) / ((1.0 + b) * pow_x1);
                if (aux > w) {
                    xold = exp(log(t + w * (1.0 + b) * pow_x1 / (h_ff[i * MXRAYFF + ibin - 1] * h_ff[i * MXRAYFF + ibin - 1])) / (2.0 + 2.0 * b));
                    h_i_array[i * RAYCDFSIZE + j - 1] = ibin;
                    break;
                }
                w -= aux;
                xold = x2;
                ibin++;
                b = h_b_array[i * MXRAYFF + ibin - 1];
                pow_x1 = pow(xold, 2.0 * b);
            } while (true);
        }

        /*************************************************************************
         * Change definition of b_array because that's what is needed at run time
         **************************************************************************/
        for (uint j = 0; j < MXRAYFF; j++)
            h_b_array[i * MXRAYFF + j] = 0.5 / (1.0 + h_b_array[i * MXRAYFF + j]);

        // prepare coefficients for pmax interpolation
        //dle = log(media[i]->up / media[i]->ap) / ((double)MXGE - 1.0);
        //double dlei = 1.0 / dle;
        
        for (uint j = 0; j < MXGE - 1; j++) {
            double gle = ((double)j - ge.x) / ge.y;
            h_pmax[i * MXGE + j].y = (h_pe_array[i * MXGE + j + 1] - h_pe_array[i * MXGE + j]) * ge.y;
            h_pmax[i * MXGE + j].x = h_pe_array[i * MXGE + j] - h_pmax[i * MXGE + j].y * gle;
        }

        h_pmax[i * MXGE + MXGE - 1] = h_pmax[i * MXGE + MXGE - 2];

    }

    // convert to floats
    rayleigh_data_t *h_rayleigh_data = (rayleigh_data_t*)malloc(nmed * MXRAYFF * sizeof(rayleigh_data_t));
    float2 *h_pmax_f = (float2*)malloc(nmed * MXGE * sizeof(float2));

    for (uint i = 0; i < nmed * MXRAYFF; i++) {
        h_rayleigh_data[i].xgrid = (float)h_xgrid[i];
        h_rayleigh_data[i].fcum = (float)h_fcum[i];
        h_rayleigh_data[i].b_array = (float)h_b_array[i];
        h_rayleigh_data[i].c_array = (float)h_c_array[i];
    }

    for (uint i = 0; i < nmed * MXGE; i++)
        h_pmax_f[i] = make_float2((float)h_pmax[i].x, (float)h_pmax[i].y);

    // free host double memory
    free(h_ff);
    free(h_xgrid);
    free(h_fcum);
    free(h_b_array);
    free(h_c_array);
    free(h_pe_array);
    free(h_pmax);

	d_rayleigh_data = (rayleigh_data_t**)malloc(GPUNo*sizeof(rayleigh_data_t*));
	d_i_array = (uint**)malloc(GPUNo*sizeof(uint*));
	d_pmax = (float2**)malloc(GPUNo*sizeof(float2*));

	for(int GPUId=0; GPUId<GPUNo; GPUId++) {

#ifdef USE_MULTIPLE_GPU
		cudaSetDevice(GPUId); ce(58007);
#endif

		// allocate device memory
		cudaMalloc(&d_rayleigh_data[GPUId], nmed * MXRAYFF * sizeof(rayleigh_data_t)); ce(10020);
		cudaMalloc(&d_i_array[GPUId], nmed * RAYCDFSIZE * sizeof(uint)); ce(10021);
		cudaMalloc(&d_pmax[GPUId], nmed * MXGE * sizeof(float2)); ce(10022);

		// copy data to device
		cudaMemcpy(d_rayleigh_data[GPUId], h_rayleigh_data, nmed * MXRAYFF * sizeof(rayleigh_data_t), cudaMemcpyHostToDevice); ce(10023);
		cudaMemcpy(d_i_array[GPUId], h_i_array, nmed * RAYCDFSIZE * sizeof(uint), cudaMemcpyHostToDevice); ce(10024);
		cudaMemcpy(d_pmax[GPUId], h_pmax_f, nmed * MXGE * sizeof(float2), cudaMemcpyHostToDevice); ce(10025);

		cudaMemcpyToSymbol(rayleigh_data, &d_rayleigh_data[GPUId], sizeof(rayleigh_data_t*)); ce(10026);
		cudaMemcpyToSymbol(i_array, &d_i_array[GPUId], sizeof(int*)); ce(10027);
		cudaMemcpyToSymbol(pmax, &d_pmax[GPUId], sizeof(float2*)); ce(10028);
	
	}

    // free host single memory
    free(h_rayleigh_data);
    free(h_i_array);
    free(h_pmax_f);
}

medium_t** init_media(uint nmed, string *media_names) {
    // read the data of the required media from the pegs file    
    medium_t **media = (medium_t**)malloc(nmed * sizeof(medium_t*));
    
    bool *found = (bool*)malloc(nmed * sizeof(bool));
    for (uint i = 0; i < nmed; i++)
        found[i] = false;

    uint media_found = read_pegs_file(pegs_file, nmed, media_names, media, found);

    // did not find all media
    if (media_found < nmed) {
        if (nmed - media_found > 1)
            logstr("\nERROR (10029): The following media were not found or could not be read from the PEGS file:");
        else
            logstr("\nERROR (10029): The following mediun was not found or could not be read from the PEGS file:");

        for (uint i = 0; i < nmed; i++) {
            if (!found[i])
                printf(" %s", media_names[i].c_str());
        }

        free(found);
        
        fflush(logfile);
        exit(10029);
    }

    free(found);

    // at this point we have found and read all required media

    // init the photon data using the specified cross sections files
    init_user_photon(nmed, media, photon_xsections, "");

    // init rayleigh data
    init_rayleigh_data(nmed, media);

    return media;
}

#endif
