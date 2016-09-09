/*
 * =====================================================================================
 *
 *       Filename:  point.desc.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  05.09.2012 13:50:36
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */

#ifndef FEMINIST_POINT_DESC

#define FEMINIST_POINT_DESC

// #include <projects/prmt_sintactic_addition/prmt_sintactic_addition.h>
#include "../../../../../prmt_sintactic_addition/types.h"
#include "../../../../../prmt_sintactic_addition/constructions.h"
//#include </home/primat/projects/common/common_without_report.h>
#include <deal.II/base/point.h>
#include <math.h>

namespace prmt
{
    //! Точка приматского формата
    template <int dim>
        class Point
        {
            private:
                Point () {};
        };

    //! Точка с параметром
    template <typename T, int dim>
        class PointWithParam
        {
            public:
                Point<dim> p;
                T param;
        };

    //! Точка с типом
    template <int dim>
        using PointWithType = PointWithParam<i32, dim>;
};


#endif
