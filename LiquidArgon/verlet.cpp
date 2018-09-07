#include "MDFunctions.h"

void verlet(vector<double> &coord,const vector<double> &coor_t0,
            vector<double> &vel,const vector<double> &vel_t0,
            vector<double> &acc,const vector<double> &acc_t0,
            vector<double> &force,double &pener,
            const int &natom,const double &mass,
            const double &dt,const double (&boxl)[3])
{
    int i,j,k;
    double r[3],f[3];
    double rr,r2,r6;

    // initializing
    for(i=1;i<=natom;i++)
    {
        coord[(i-1)*3+1-1]=0.0;
        coord[(i-1)*3+2-1]=0.0;
        coord[(i-1)*3+3-1]=0.0;

        vel[(i-1)*3+1-1]=0.0;
        vel[(i-1)*3+2-1]=0.0;
        vel[(i-1)*3+3-1]=0.0;

        acc[(i-1)*3+1-1]=0.0;
        acc[(i-1)*3+2-1]=0.0;
        acc[(i-1)*3+3-1]=0.0;

        force[(i-1)*3+1-1]=0.0;
        force[(i-1)*3+2-1]=0.0;
        force[(i-1)*3+3-1]=0.0;
    }
    pener=0.0;

    for(i=1;i<=natom;i++)
    {
        coord[(i-1)*3+1-1]=coor_t0[(i-1)*3+1-1]+vel_t0[(i-1)*3+1-1]*dt
                          +0.5*acc_t0[(i-1)*3+1-1]*dt*dt;
        coord[(i-1)*3+2-1]=coor_t0[(i-1)*3+2-1]+vel_t0[(i-1)*3+2-1]*dt
                          +0.5*acc_t0[(i-1)*3+2-1]*dt*dt;
        coord[(i-1)*3+3-1]=coor_t0[(i-1)*3+3-1]+vel_t0[(i-1)*3+3-1]*dt
                          +0.5*acc_t0[(i-1)*3+3-1]*dt*dt;

        // Apply PBC to coordinates
        for(j=1;j<=3;j++)
        {
            if(coord[(i-1)*3+j-1]>boxl[j-1])
            {
                coord[(i-1)*3+j-1]-=boxl[j-1];
            }
            else if(coord[(i-1)*3+j-1]<0.0)
            {
                coord[(i-1)*3+j-1]+=boxl[j-1];
            }
        }
    }

    //*************************************************************
    // Get force at new atom positions
    // Using Lennard Jones Potential
    // Hint: you might want to also seperate the potential and force calculation into a separate subroutine
    // this will be useful if you want to use other potentials
    //*************************************************************
    for(i=1;i<=natom-1;i++)
    {
        for(j=i+1;j<=natom;j++)
        {
            r[1-1]=coord[(i-1)*3+1-1]-coord[(j-1)*3+1-1];
            r[2-1]=coord[(i-1)*3+2-1]-coord[(j-1)*3+2-1];
            r[3-1]=coord[(i-1)*3+3-1]-coord[(j-1)*3+3-1];

            //minimum image criterion
            // fortran version is: r = r - nint( r / boxl ) * boxl
            r[0]=r[0]-nint(r[0]/boxl[0])*boxl[0];
            r[1]=r[1]-nint(r[1]/boxl[1])*boxl[1];
            r[2]=r[2]-nint(r[2]/boxl[2])*boxl[2];

            rr=r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
            r2=1.0/rr;
            r6=r2*r2*r2;
            //Lennard Jones Potential
            // V = 4 * epsilon * [ (sigma/r)**12 - (sigma/r)**6 ]
            //   = 4 * epsilon * (sigma/r)**6 * [ (sigma/r)**2 - 1 ]
            //   = 4 * r**(-6) * [ r**(-2) - 1 ] for epsilon=sigma=1
            // F_i = 48 * epsilon * (sigma/r)**6 * (1/r**2) * [ ( sigma/r)** 6 - 0.5 ] * i where i = x,y,z
            //     = 48 * r**(-8) * [ r**(-6) - 0.5 ] * i  for epsilon=sigma=1

            pener=pener+4.0*r6*(r6-1.0);

            f[0]=48.0*r2*r6*(r6-0.5)*r[0];
            f[1]=48.0*r2*r6*(r6-0.5)*r[1];
            f[2]=48.0*r2*r6*(r6-0.5)*r[2];

            force[(i-1)*3+1-1]+=f[1-1];
            force[(i-1)*3+2-1]+=f[2-1];
            force[(i-1)*3+3-1]+=f[3-1];

            force[(j-1)*3+1-1]-=f[1-1];
            force[(j-1)*3+2-1]-=f[2-1];
            force[(j-1)*3+3-1]-=f[3-1];
        }
    }

    //*********************************************************
    // Calculate Acceleration and Velocity  at current time step
    //*********************************************************
    for(i=1;i<=natom;i++)
    {
        acc[(i-1)*3+1-1]=force[(i-1)*3+1-1]/mass;
        acc[(i-1)*3+2-1]=force[(i-1)*3+2-1]/mass;
        acc[(i-1)*3+3-1]=force[(i-1)*3+3-1]/mass;

        vel[(i-1)*3+1-1]=vel_t0[(i-1)*3+1-1]+0.5*(acc[(i-1)*3+1-1]+acc_t0[(i-1)*3+1-1])*dt;
        vel[(i-1)*3+2-1]=vel_t0[(i-1)*3+2-1]+0.5*(acc[(i-1)*3+2-1]+acc_t0[(i-1)*3+2-1])*dt;
        vel[(i-1)*3+3-1]=vel_t0[(i-1)*3+3-1]+0.5*(acc[(i-1)*3+3-1]+acc_t0[(i-1)*3+3-1])*dt;
    }
}