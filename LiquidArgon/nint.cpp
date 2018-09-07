#include "MDFunctions.h"

int nint(double val)
{
    // it seems c++ should include nint function
    // as a standard one
    if(val<0.0)
    {
        return int(val-0.5);
    }
    else
    {
        return int(val+0.5);
    }
}