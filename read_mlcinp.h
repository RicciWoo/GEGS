#define TOTAL_MLC_LEAF 120
struct Field_t
{
    float angle;
	float x1;
	float x2;
	float y1;
	float y2;
	int   MU;
	float output_factor;
	string mlcinp_file;
	float pos[TOTAL_MLC_LEAF];
};
