/*
 * =====================================================================================
 *
 *       Filename:  point.impl.h
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

#ifndef FEMENIST_POINT_IMPL

#define FEMENIST_POINT_IMPL

#include "./point.desc.h"

namespace prmt
{
    template <>
        class Point<1>
        {
            public:
                Point () : _x(0.0) {};
                Point (dbl x) { _x = x; };

                dbl& x() { return _x; };

                dbl x() const { return _x; };
            private:
                dbl _x;
        };

    template <>
        class Point<2>
        {
            public:
                Point () : _x(0.0), _y(0.0) {};
                Point (cdbl x, cdbl y) : _x(x), _y(y) {};
                Point (const dealii::Point<2, dbl> &p) : _x(p(0)), _y(p(1)) {};

                dbl& x() { return _x; };
                dbl& y() { return _y; };

                dbl x() const { return _x; };
                dbl y() const { return _y; };

                operator dealii::Point<2, dbl> () const
                {
                    return dealii::Point<2, dbl>(_x, _y);
                };

                void operator= (const dealii::Point<2, dbl> &p)
                {
                    _x = p(0);
                    _y = p(1);
                };
                
                dbl distance (const Point<2> &p) const
                {
                    return sqrt(pow(_x - p.x(), 2.0) + pow(_y - p.y(), 2.0));
                };
                        
                bool operator==(const Point &b) const
                {
                    if (sqrt(pow(_x - b.x(), 2.0) + pow(_y - b.y(), 2.0)) < imprecision)
                        return true;
                    else
                        return false;
                };
                                
                bool operator==(const dealii::Point<2, dbl> &b) const
                {
                    if (sqrt(pow(_x - b(0), 2.0) + pow(_y - b(1), 2.0)) < imprecision)
                        return true;
                    else
                        return false;
                };
                
            private:
                dbl _x;
                dbl _y;
                dbl imprecision = 1e-10;
        };

    template <>
        class Point<3>
        {
            public:
                Point () : _x(0.0), _y(0.0), _z(0.0) {};
                Point (dbl x, dbl y, dbl z) { _x = x; _y = y; _z = z;};

                dbl& x() { return _x; };
                dbl& y() { return _y; };
                dbl& z() { return _z; };

                dbl x() const { return _x; };
                dbl y() const { return _y; };
                dbl z() const { return _z; };
            private:
                dbl _x;
                dbl _y;
                dbl _z;
        };
};

#endif
