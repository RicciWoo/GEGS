#include <fstream>
//#include <cstring>
using namespace std;
#include "Polygon.h"
#include "read_egsinp.h"

void read_egsinp(point2D *cross)
{
    int i,j,N;
    float X;
    char char_temp[100];
    char ch;
    char CM_title[32]="*********** start of CM DYNVMLC";  //31 characters
    int N_full;
    int N_half;
    float rmax;
	float zmin;
    float zthick;
    float x_start;
    float x_end;
    float x_mid_1;
    float x_mid_2;
    float air_gap;
    float rad_end;
	Leaf full,target,isocenter;

	ifstream infile;
    infile.open("VarianMLC.egsinp");
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
                infile>>rmax>>ch;   infile.ignore(100,'\n');
                for(i=0;i<2;i++)
                {
                    infile.ignore(100,'\n');
                }
                infile>>zmin>>ch;   infile.ignore(100,'\n');
                infile>>zthick>>ch; infile.ignore(100,'\n');
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
                infile>>x_start>>ch;  infile.ignore(100,'\n');
                infile>>air_gap>>ch;  infile.ignore(100,'\n');
                                      infile.ignore(100,'\n'); //skip ENDTYPE. suppose end type is round
                infile>>rad_end>>ch;  infile.ignore(100,'\n');
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

                full.zmin=target.zmin=isocenter.zmin=zmin;
                full.zthick=target.zthick=isocenter.zthick=zthick;
                full.VerCal_full();
                target.VerCal_target();
                isocenter.VerCal_isocenter();

                //there are 60 leaves, 14 vertice for each leaf, we use an array of 16*60 vertice to save these positions
                //in order to adapt the process of data transportation of GPU, each time it transports a data block of 64 bytes,
                //we extend from 14 vertice to 16 for each leaf, so the data for each leaf is saved within one specific data block
                //here for each leaf, we save the position data in the first 14 vertice in array for each leaf
                X=x_start;    //x_start is the start position in x coordinate of leaves
                N=N_full;
                for(i=0;i<N;i++)
                {
                    for(j=0;j<14;j++)
                    {
                        cross[i*16+j].x=(X+full.v[j].x)/zmin*full.v[j].z;
                        cross[i*16+j].y=full.v[j].z-zmin;
                    }
                    X=X+full.wl+air_gap;
                }
                x_mid_1=X;    //x_mid_1 is the start position of half leaves
                N=N+N_half;
                for(;i<N;i++)
                {
                    for(j=0;j<14;j++)
                    {
                        cross[i*16+j].x=(X+target.v[j].x)/zmin*target.v[j].z;
                        cross[i*16+j].y=target.v[j].z-zmin;
                    }
                    X=X+target.wl+air_gap;
                    i++;
                    for(j=0;j<14;j++)
                    {
                        cross[i*16+j].x=(X+isocenter.v[j].x)/zmin*isocenter.v[j].z;
                        cross[i*16+j].y=isocenter.v[j].z-zmin;
                    }
                    X=X+isocenter.wl+air_gap;
                }
                x_mid_2=X;    //x_mid_2 is the start position of second part of full leaves
                N=N+N_full;
                for(;i<N;i++)
                {
                    for(j=0;j<14;j++)
                    {
                        cross[i*16+j].x=(X+full.v[j].x)/zmin*full.v[j].z;
                        cross[i*16+j].y=full.v[j].z-zmin;
                    }
                    X=X+full.wl+air_gap;
                }
                x_end=X-air_gap+full.wg;    //x_end is the end position of leaves

                //we use the fifteenth vertex of each leaf to save the open position of the leaf
                i=0;N=0;
                list=head;
                do
                {
                    N=N+list->NUM;
                    for(;i<N;i++)
                    {
                        cross[14+i*16].x=list->NEG;
                        cross[14+i*16].y=list->POS;
                    }
                    before=list;
                    list=list->next;
                    delete before;
                }
                while(list!=NULL);
                //we save some relevant data in the remaining vertice of the array
                cross[15].x=zmin;
                cross[15].y=zthick;
                cross[31].x=x_start;
                cross[31].y=x_end;
                cross[47].x=x_mid_1;
                cross[47].y=x_mid_2;
                cross[63].x=rad_end;
                cross[63].y=air_gap;
                cross[79].x=full.wl;
                cross[79].y=target.wl;
                cross[95].x=isocenter.wl;
                cross[95].y=full.wt;
                cross[111].x=rmax;

                infile.close();
                return;
            }
        }
    }
    infile.close();
    return;
}

