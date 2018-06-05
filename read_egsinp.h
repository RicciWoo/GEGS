
struct vertex
{
    float x;
    float z;
};

struct link_field
{
    float NEG;
    float POS;
    int NUM;
    link_field *next;
};

struct Leaf       //basic class to store the geometric data of each leaf
{
  public:
    float zmin;
    float zthick;
    float wl;    //leaf (excl. tongue)
    float wt;    //tongue
    float wg;    //groove
    float wtip;  //tip at top/bottom of leaf
    float wts;   //top support rail
    float wbs;   //bottom support rail
    float ztip;  //top/bottom of tip
    float zl;    //top/bottom of leaf
    float zt;    //bottom/top of tongue
    float zg;    //bottom/top of groove
    float zth;   //top of driving screw hole
    float zbh;   //bottom of driving screw hole
    float holepos;  //distance of hole from leaf end
    float zts;   //top of support rail
    float zbs;   //bottom of support rail
	vertex v[14];
    void VerCal_full()       //vertex position calculation function for full leaf
    {
        v[0].z=ztip;   v[0].x=0;
        v[1].z=ztip;   v[1].x=wtip;
        v[2].z=zl;     v[2].x=wtip;
        v[3].z=zl;     v[3].x=wt+wl-wg;
        v[4].z=zg;     v[4].x=wt+wl-wg;
        v[5].z=zg;     v[5].x=wt+wl;
        v[6].z=zmin+zthick; v[6].x=wt+wl;
        v[7].z=zmin+zthick; v[7].x=wt+wts-wbs;
        v[8].z=zbs;    v[8].x=wt+wts-wbs;
        v[9].z=zbs;    v[9].x=wt+wts;
        v[10].z=zts;   v[10].x=wt+wts;
        v[11].z=zts;   v[11].x=wt;
        v[12].z=zt;    v[12].x=wt;
        v[13].z=zt;    v[13].x=0;
    }
    void VerCal_target()         //vertex position calculation function for target leaf
    {
        v[0].z=zmin;   v[0].x=0;
        v[1].z=zmin;   v[1].x=wt+wl-wbs+wts;
        v[2].z=zts;    v[2].x=wt+wl-wbs+wts;
        v[3].z=zts;    v[3].x=wt+wl-wbs;
        v[4].z=zbs;    v[4].x=wt+wl-wbs;
        v[5].z=zbs;    v[5].x=wt+wl;
        v[6].z=zg;     v[6].x=wt+wl;
        v[7].z=zg;     v[7].x=wt+wl-wg;
        v[8].z=ztip;   v[8].x=wt+wl-wg;
        v[9].z=ztip;   v[9].x=wt+wl-wg-wtip;
        v[10].z=zl;    v[10].x=wt+wl-wg-wtip;
        v[11].z=zl;    v[11].x=wt;
        v[12].z=zt;    v[12].x=wt;
        v[13].z=zt;    v[13].x=0;
    }
    void VerCal_isocenter()            //vertex position calculation function for isocenter leaf
    {
        v[0].z=zt;     v[0].x=0;
        v[1].z=zt;     v[1].x=wt;
        v[2].z=ztip;   v[2].x=wt;
        v[3].z=ztip;   v[3].x=wt+wtip;
        v[4].z=zl;     v[4].x=wt+wtip;
        v[5].z=zl;     v[5].x=wt+wl-wg;
        v[6].z=zg;     v[6].x=wt+wl-wg;
        v[7].z=zg;     v[7].x=wt+wl;
        v[8].z=zmin+zthick; v[8].x=wt+wl;
        v[9].z=zmin+zthick; v[9].x=wts-wbs;
        v[10].z=zbs;   v[10].x=wts-wbs;
        v[11].z=zbs;   v[11].x=wts;
        v[12].z=zts;   v[12].x=wts;
        v[13].z=zts;   v[13].x=0;
    }
};
