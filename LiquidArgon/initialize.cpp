#include "MDFunctions.h"

void initialize(int natom,int npartdim,double alat,vector<double> &coord_t0,vector<double> &vel_t0)
{
    double rcell[3][4];// for FCC unit cell

    rcell[1-1][1-1]=0.0;
    rcell[1-1][2-1]=0.5*alat;
    rcell[1-1][3-1]=0.0;
    rcell[1-1][4-1]=0.5*alat;

    rcell[2-1][1-1]=0.0;
    rcell[2-1][2-1]=0.5*alat;
    rcell[2-1][3-1]=0.5*alat;
    rcell[2-1][4-1]=0.0;

    rcell[3-1][1-1]=0.0;
    rcell[3-1][2-1]=0.0;
    rcell[3-1][3-1]=0.5*alat;
    rcell[3-1][4-1]=0.5*alat;

    int n,i,j,k,l;

    // For a FCC crystal structure()
    // TODO: find the literature for the definition of FCC
    n=1;
    for(i=1;i<=npartdim;i++)
    {
        for(j=1;j<=npartdim;j++)
        {
            for(k=1;k<=npartdim;k++)
            {
                for(l=1;l<=4;l++)
                {
                    coord_t0[(n-1)*3+1-1]=alat*(i-1)+rcell[1-1][l-1];
                    coord_t0[(n-1)*3+2-1]=alat*(j-1)+rcell[2-1][l-1];
                    coord_t0[(n-1)*3+3-1]=alat*(k-1)+rcell[3-1][l-1];
                    n+=1;
                }
            }
        }
    }

    // initializing the random velocity
    for(i=1;i<=natom;i++)
    {
        vel_t0[(i-1)*3+1-1]=gasdev();
        vel_t0[(i-1)*3+2-1]=gasdev();
        vel_t0[(i-1)*3+3-1]=gasdev();
    }
    
}