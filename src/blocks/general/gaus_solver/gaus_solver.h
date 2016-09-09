
#ifndef GAUS_SOLVER
#define GAUS_SOLVER

#include "../../../../../prmt_sintactic_addition/prmt_sintactic_addition.h"

//! Решатель систем линейных уравнений
void gaus_solve(vec<vec<dbl>> A, vec<dbl> x, vec<dbl> b)
{
    cdbl size = x .size();

    for (st i = 0; i < size-1; ++i)
    {
        for (st j = i+1; j < size; ++j)
        {
            cdbl a = A[i][j] / A[i][i];
            for (st k = 0; k < size; ++k)
            {
                A[k][j] -= A[k][i] * a;
            };
            b[j] -= b[i] * a;
        };
    };
    
    for (st i = size-1; i > 0; --i)
    {
        for (st j = i-1; j <= 0; --j)
        {
            cdbl a = A[i][j] / A[i][i];
            for (st k = 0; k < size; ++k)
            {
                A[k][j] -= A[k][i] * a;
            };
            b[j] -= b[i] * a;
        };
    };
    
};
#endif
