#include "MDFunctions.h"

void get_temp(const int &natom,vector<double> &vel,double &avtemp,const double &mass,const double &kb)
{
    int i,j;
    double ke;

    ke=0.0;
    for(i=1;i<=natom;i++)
    {
        ke+=vel[(i-1)*3+1-1]*vel[(i-1)*3+1-1]
           +vel[(i-1)*3+2-1]*vel[(i-1)*3+2-1]
           +vel[(i-1)*3+3-1]*vel[(i-1)*3+3-1];
    }

    avtemp=mass*ke/(3.0*kb*(natom-1));
}