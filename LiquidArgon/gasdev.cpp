#include "MDFunctions.h"

double gasdev()
{
    double result,v1,v2,fac,rsq;
    static double gset;
    static bool available=false;

    if(available)
    {
        result=gset;
        available=false;
    }
    else
    {
        while(1)
        {
            v1=1.0*random()/RAND_MAX;
            v2=1.0*random()/RAND_MAX;

            v1=2.0*v1-1.0;
            v2=2.0*v2-1.0;
            rsq=v1*v1+v2*v2;
            if(rsq>0.0 && rsq<1.0) break;
        }

        fac=sqrt(-2.0*log(rsq)/rsq);
        result=v1*fac;
        gset=v2*fac;
        available=true;
    }
    return result;
}