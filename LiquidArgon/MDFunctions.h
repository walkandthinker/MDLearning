#ifndef LIQUIDARGON_MDFUNCTIONS_H
#define LIQUIDARGON_MDFUNCTIONS_H

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

int nint(double val);

double gasdev();
void initialize(int natom,int npartdim,double alat,vector<double> &coord_t0,vector<double> &vel_t0);

void verlet(vector<double> &coord,const vector<double> &coor_t0,
            vector<double> &vel,const vector<double> &vel_t0,
            vector<double> &acc,const vector<double> &acc_t0,
            vector<double> &force,double &pener,
            const int &natom,const double &mass,
            const double &dt,const double (&boxl)[3]);

void linearmom(const int &natom,vector<double> &vel);

void get_temp(const int &natom,vector<double> &vel,double &avtemp,const double &mass,const double &kb);

#endif