#include <iostream>
#include <iomanip>
#include <ctime>
#include <cmath>
#include <fstream>
#include <string>


#include "MDFunctions.h"

using namespace std;

int main()
{
    cout<<"Hi, I'm liquid argon simulator!"<<endl;
    srand(time(0));// start the random number seed at the beginning
    
    // MD Parameters
    int npartdim,natom;
    int nstep;
    double tempK,dt;
    int istep,n,i,j,k,l;
    double boxl[3],alat;
    vector<double> coord_t0,coord;
    vector<double> vel_t0,vel;
    vector<double> acc_t0,acc;
    vector<double> force;
    double pener,mass;
    double vcm[3],r[3],rr,r2,r6,f[3];
    double avtemp,ke,kb,epsilon,sigma,scale;
    //For timing analysis
    //double init_time,start_time,end_time;
    //int value;

    cout<<"*** Please input the npartdim=";
    cin>>npartdim;
    cout<<"*** Please input the time step=";
    cin>>nstep;

    natom=npartdim*npartdim*npartdim*4;
    // initializing array
    coord.resize(3*natom,0.0);coord_t0.resize(3*natom,0.0);
    vel.resize(3*natom,0.0);vel_t0.resize(3*natom,0.0);
    acc.resize(3*natom,0.0);acc_t0.resize(3*natom,0.0);
    force.resize(3*natom,0.0);

    // set temperature and delta t
    tempK=10.0;dt=1.0e-3;

    alat=pow(2.0,2.0/3.0);
    boxl[0]=npartdim*alat;boxl[1]=npartdim*alat;boxl[2]=npartdim*alat;
    kb=1.0;mass=1.0;
    epsilon=1.0;sigma=1.0; // what's this two parameter means?

    // start MD initializing
    cout<<"*** Start md initilaizing...     ***"<<endl;
    initialize(natom,npartdim,alat,coord_t0,vel_t0); 
    cout<<"*** MD initilaizing finished!    ***"<<endl;

    // write out the initial position of atoms
    const string filename="atom.xyz";
    ofstream out;
    out.open(filename,ios::out);
    out<<natom<<endl<<endl;
    for(i=1;i<=natom;i++)
    {
        out<<"Ar ";
        out<<scientific<<setprecision(6)
           <<coord_t0[(i-1)*3+0]<<" "
           <<coord_t0[(i-1)*3+1]<<" "
           <<coord_t0[(i-1)*3+2]<<"\n";
    }
    out<<endl;

    // set Linear Momentum to zero
    linearmom(natom,vel_t0);

    // get current temperature
    get_temp(natom,vel_t0,avtemp,mass,kb);
    cout<<"*** Initial average temperature is:"<<scientific<<setprecision(6)<<avtemp<<endl;

    // scale initial velocity do desired temperature
    scale=sqrt(tempK/avtemp);
    for(i=1;i<=natom;i++)
    {
        vel_t0[(i-1)*3+0]*=scale;
        vel_t0[(i-1)*3+1]*=scale;
        vel_t0[(i-1)*3+2]*=scale;
    }
    get_temp(natom,vel_t0,avtemp,mass,kb);
    cout<<"*** Initial scaled average temperature is:"<<scientific<<setprecision(6)<<avtemp<<endl;
    // it seems scaled temperature is always 1.0 ?

    //**********************************************
    //*** Now we start MD simulation
    //**********************************************
    for(istep=1;istep<=nstep;istep++)
    {
        // Get new atom positions from Velocity Verlet Algorithm
        verlet(coord,coord_t0,vel,vel_t0,acc,acc_t0,force,pener,natom,mass,dt,boxl);

        // Set Linear Momentum to zero
        linearmom(natom,vel);// this is not zero, just kind of shift with average velocity?

        // compute average temperature
        get_temp(natom,vel,avtemp,mass,kb);
        cout<<"*** step="<<istep<<":";
        cout<<" averate temp="<<scientific<<setprecision(6)<<avtemp;
        cout<<", pener="<<scientific<<setprecision(6)<<pener<<endl;

        scale=sqrt(tempK/avtemp);
        //*************************
        //*** update u,v,a
        for(i=1;i<=natom;i++)
        {
            coord_t0[(i-1)*3+0]=coord[(i-1)*3+0];
            coord_t0[(i-1)*3+1]=coord[(i-1)*3+1];
            coord_t0[(i-1)*3+2]=coord[(i-1)*3+2];

            vel_t0[(i-1)*3+0]=vel[(i-1)*3+0]*scale;
            vel_t0[(i-1)*3+1]=vel[(i-1)*3+1]*scale;
            vel_t0[(i-1)*3+2]=vel[(i-1)*3+2]*scale;

            acc_t0[(i-1)*3+0]=acc[(i-1)*3+0];
            acc_t0[(i-1)*3+1]=acc[(i-1)*3+1];
            acc_t0[(i-1)*3+2]=acc[(i-1)*3+2];
        }

        // output
        out<<natom<<endl<<endl;
        for(i=1;i<=natom;i++)
        {
            out<<"Ar ";
            out<<scientific<<setprecision(6)
            <<coord_t0[(i-1)*3+0]<<" "
            <<coord_t0[(i-1)*3+1]<<" "
            <<coord_t0[(i-1)*3+2]<<"\n";
        }
        out<<endl;

    }

    out.close();

    return 0;
}