#include "MDFunctions.h"

void linearmom(const int &natom,vector<double> &vel)
{
    int i,j;
    double vcm[3];

    vcm[0]=0.0;vcm[1]=0.0;vcm[2]=0.0;
    for(i=1;i<=natom;i++)
    {
        vcm[0]+=vel[(i-1)*3+1-1]/natom;
        vcm[1]+=vel[(i-1)*3+2-1]/natom;
        vcm[2]+=vel[(i-1)*3+3-1]/natom;
    }

    for(i=1;i<=natom;i++)
    {
        vel[(i-1)*3+0]-=vcm[0];
        vel[(i-1)*3+1]-=vcm[1];
        vel[(i-1)*3+2]-=vcm[2];
    }
}