#ifndef points_function_def
#define points_function_def 1

#include <cmath>
#include "../../../../../prmt_sintactic_addition/prmt_sintactic_addition.h"
#include "../point/point.h"

//! Скалярная функция, построенныя по 4 известным точкам
template <u8 dim>
class Scalar4PointsFunc
{
};

//! Скалярная функция, построенныя по 4 известным точкам, двумерный случай
/*! Глобальная система координат (x, y) выражается через локальную (s, r), в
 * которой 4 точки образуют единичный квадрат.
 * По формулам:
 * \f[
 * \mathrm{x(s, r) = a_x + b_x * s + c_x * r + d_x * s * r} \\
 * \mathrm{y(s, r) = a_y + b_y * s + c_y * r + d_y * s * r}
 * \f]
 * \image html isoparametric_coordinate_transformation.img
 * Функця строющиеся по 4 точкам в локальной системе координат является
 * билинейной:
 * \f[
 * \mathrm{f(s, r) =  a_f + b_f * s + c_f * r + d_f * s * r}
 * \f]
 * где
 * \f[
\begin{equation}
    \begin{array}{l}
        \mathrm{a_x = x_1} \\
        \mathrm{b_x = x_2 - x_1} \\
        \mathrm{c_x = x_4 - x_1} \\
        \mathrm{d_x = (x_1 + x_3) - (x_4 + x_2)} \\

        \mathrm{a_y = y_1} \\
        \mathrm{b_y = y_2 - y_1} \\
        \mathrm{c_y = y_4 - y_1} \\
        \mathrm{d_y = (y_1 + y_3) - (y_4 + y_2)} \\

        \mathrm{a_f = f_1} \\
        \mathrm{b_f = f_2 - f_1} \\
        \mathrm{c_f = f_4 - f_1} \\
        \mathrm{d_f = (f_1 + f_3) - (f_4 + f_2)} \\
    \end{array} 
\end{equation}
 * \f]
 * Выражаем локальные координаты через глобальные.
  \f[
 * \mathrm{s(x, y)} = \left\{ \begin{eqnarray} \mathrm{\frac{B_s(x,y)-\sqrt{B_s(x,y)^2-4A_sC_s(x,y)}}{2A_s}; A_s \neq 0, A_r \neq 0} \\ \mathrm{\frac{C_s(x, y)}{B_s(x, y)}; A_s = 0, A_r = 0}  \end{eqnarray} \right.  \\
 * \f]
  \f[
 * \mathrm{r(x, y)} = \left\{ \begin{eqnarray} \mathrm{\frac{B_r(x,y)-\sqrt{B_r(x,y)^2-4A_rC_r(x,y)}}{2A_r}; A_s \neq 0, A_r \neq 0} \\ \mathrm{\frac{C_r(x, y)}{B_r(x, y)}; A_s = 0, A_r = 0}  \end{eqnarray} \right.  \\
 * \f]
  \f[
\begin{equation}
    \begin{array}{l}
        \mathrm{A_s = b_yd_x - b_xd_y} \\
        \mathrm{A_r = c_yd_x - c_xd_y} \\

        \mathrm{B_s(x,y) = a_xd_y - a_yd_x - b_yc_x + b_xc_y - d_yx + d_xy} \\
        \mathrm{B_r(x,y) = a_xd_y - a_yd_x + b_yc_x - b_xc_y - d_yx + d_xy} \\

        \mathrm{C_s(x,y) = a_yc_x - a_xc_y + c_yx - c_xy} \\
        \mathrm{C_r(x,y) = a_yb_x - a_xb_y + b_yx - b_xy}
    \end{array} 
\end{equation}
 * \f]
 * Получаем функцию в глобальных координатах:
 * \f[
 * \mathrm{f(x,y) = f(s(x,y),r(x,y))}
 * \f]
 * Градиент функции в локальных координатах будет иметь вид:
 * \f[
        \mathrm{\frac{\partial f}{\partial x} = b_f\frac{\partial s}{\partial x}+c_f\frac{\partial r}{\partial x}+d_f(\frac{\partial s}{\partial x}r + s\frac{\partial r}{\partial x})} \\ 
        \mathrm{\frac{\partial f}{\partial y} = b_f\frac{\partial s}{\partial y}+c_f\frac{\partial r}{\partial y}+d_f(\frac{\partial s}{\partial y}r + s\frac{\partial r}{\partial y})}
 * \f]
  \f[
 * \mathrm{\frac{\partial s}{\partial x}} = \left\{\begin{eqnarray} \mathrm{\frac{\partial B_s}{\partial x}\frac{1}{2A} - \frac{B\frac{\partial B_s}{\partial x}-2A\frac{\partial C_s}{\partial x}}{\sqrt{B_s^2-4A_sC_s}}; A_s \neq 0, A_r \neq 0} \\ \mathrm{\frac{\frac{\partial C_s}{\partial x}B_s-C_s\frac{\partial B_s}{\partial x}}{B_s^2}; A_s = 0, A_r = -1}  \end{eqnarray} \right.  \\
 * \f]
  \f[
 * \mathrm{\frac{\partial r}{\partial x}} = \left\{\begin{eqnarray} \mathrm{\frac{\partial B_r}{\partial x}\frac{1}{2A} - \frac{B\frac{\partial B_r}{\partial x}-2A\frac{\partial C_r}{\partial x}}{\sqrt{B_r^2-4A_rC_r}}; A_s \neq 0, A_r \neq 0} \\ \mathrm{\frac{\frac{\partial C_r}{\partial x}B_r-C_r\frac{\partial B_r}{\partial x}}{B_r^2}; A_s = 0, A_r = 0}  \end{eqnarray} \right.  \\
 * \f]
  \f[
 * \mathrm{\frac{\partial s}{\partial y}} = \left\{\begin{eqnarray} \mathrm{\frac{\partial B_s}{\partial y}\frac{1}{2A} - \frac{B\frac{\partial B_s}{\partial y}-2A\frac{\partial C_s}{\partial y}}{\sqrt{B_s^2-4A_sC_s}}; A_s \neq 0, A_r \neq 0} \\ \mathrm{\frac{\frac{\partial C_s}{\partial y}B_s-C_s\frac{\partial B_s}{\partial y}}{B_s^2}; A_s = 0, A_r = 0}  \end{eqnarray} \right.  \\
 * \f]
  \f[
 * \mathrm{\frac{\partial r}{\partial y}} = \left\{\begin{eqnarray} \mathrm{\frac{\partial B_r}{\partial y}\frac{1}{2A} - \frac{B\frac{\partial B_r}{\partial y}-2A\frac{\partial C_r}{\partial y}}{\sqrt{B_r^2-4A_rC_r}}; A_s \neq 0, A_r \neq 0} \\ \mathrm{\frac{\frac{\partial C_r}{\partial y}B_r-C_r\frac{\partial B_r}{\partial y}}{B_r^2}; A_s = 0, A_r = 0}  \end{eqnarray} \right.  \\
 * \f]
*/
template <>
class Scalar4PointsFunc<2>
{
    public:
    class LinFunc
    {
            LinFunc () = delete;
        public:
            LinFunc (cdbl _a, cdbl _b, cdbl _c);

            dbl operator() (cdbl x, cdbl y) const;
            dbl dx (cdbl x, cdbl y) const;
            dbl dy (cdbl x, cdbl y) const;

            cdbl a, b, c;
        private:
    };

    class LocVar
    {
            LocVar () = delete;
        public:
            LocVar (cdbl _A, const LinFunc _B, const LinFunc _C);

            arr<dbl,2> operator() (cdbl x, cdbl y) const;
            arr<dbl,2> dx (cdbl x, cdbl y) const;
            arr<dbl,2> dy (cdbl x, cdbl y) const;

            cdbl A;
            const LinFunc B, C;
        private:
            arr<dbl,2> quad_roots (cdbl x, cdbl y) const;

    };

    class LocFunc
    {
            LocFunc () = delete;
        public:
            LocFunc (arr<dbl, 4> glob_val);

        dbl operator() (cdbl s, cdbl r) const;
        dbl dx (cdbl s, cdbl r, cdbl s_dx, cdbl r_dx) const;
        dbl dy (cdbl s, cdbl r, cdbl s_dy, cdbl r_dy) const;

        private:
            cdbl a, b, c, d;
    };

        Scalar4PointsFunc () = delete;
    public:
        Scalar4PointsFunc (arr<prmt::Point<2>, 4> nodes,
                           arr<dbl, 4> f_values);

        dbl operator() (cdbl x, cdbl y) const;
        dbl dx (cdbl x, cdbl y) const;
        dbl dy (cdbl x, cdbl y) const;

        dbl operator() (const prmt::Point<2> &p) const;
        dbl dx (const prmt::Point<2> &p) const;
        dbl dy (const prmt::Point<2> &p) const;

        cdbl a_x, b_x, c_x, d_x, a_y, b_y, c_y, d_y;
        const LocVar s;
        const LocVar r;
        const LocFunc f_;
    private:

};

/////////////////////////////Scalar4PointsFunc
Scalar4PointsFunc<2>::Scalar4PointsFunc (arr<prmt::Point<2>, 4> nodes,
                                         arr<dbl, 4> f_values) :
    a_x (nodes[0].x()),
    b_x (nodes[1].x() - nodes[0].x()),
    c_x (nodes[3].x() - nodes[0].x()),
    d_x ((nodes[0].x() + nodes[2].x()) - (nodes[1].x() + nodes[3].x())),
    a_y (nodes[0].y()),
    b_y (nodes[1].y() - nodes[0].y()),
    c_y (nodes[3].y() - nodes[0].y()),
    d_y ((nodes[0].y() + nodes[2].y()) - (nodes[1].y() + nodes[3].y())),

    s (
            b_y * d_x - b_x * d_y, 
            LinFunc(a_x * d_y - a_y * d_x - b_y * c_x + b_x * c_y, -d_y,  d_x), 
            LinFunc(a_y * c_x - a_x * c_y,                          c_y, -c_x)
      ),

    r (
            c_y * d_x - c_x * d_y, 
            LinFunc(a_x * d_y - a_y * d_x + b_y * c_x - b_x * c_y, -d_y,  d_x), 
            LinFunc(a_y * b_x - a_x * b_y,                          b_y, -b_x)
      ),

    f_(f_values)
{
    // printf("%f %f %f %f %f %f %f %f %f %f %f\n", a_x, b_x, c_x, d_x, a_y, b_y, c_y, d_y,
    //         a_y * b_x - a_x * b_y, b_y, -b_x);
};

dbl Scalar4PointsFunc<2>::operator() (cdbl x, cdbl y) const
{
    /*
     * одной точке (x, y) кооответствуют две точки (s, r),
     * одна в правильном, единичном, квадрате, вторая в другом прямоугольнике
     * вроде бы если точки ячейки идут по кругу, правилиная точка (s[0], r[1]),
     * но на всякий случай проверяю
    */
    const auto t_s = s(x, y);
    const auto t_r = r(x, y);
    dbl valid_s;
    dbl valid_r;
    // if (
    //         (t_s[0] > -1e-10) and (t_s[0] < (1.0 + 1e-10)) and
    //         (t_r[1] > -1e-10) and (t_r[1] < (1.0 + 1e-10))
    //    )
    // {
        valid_s = t_s[0];
        valid_r = t_r[1];
    // // printf("s0 %f r1 %f s1 %f r0 %f s %f r %f\n", 
    // //         valid_s, valid_r, t_s[0], t_r[1], t_s[1], t_r[0]);
    // }
    // else
    // {
    //     valid_s = t_s[1];
    //     valid_r = t_r[0];
    // };
    // printf("s0 %f r1 %f s1 %f r0 %f s %f r %f\n", 
    //        t_s[0], t_r[1], t_s[1], t_r[0], valid_s, valid_r);

    return f_(valid_s, valid_r);
    // return f_(s(x, y)[0], r(x, y)[1]);
};

dbl Scalar4PointsFunc<2>::dx (cdbl x, cdbl y) const
{
    const auto t_s = s(x, y);
    const auto t_r = r(x, y);
    dbl valid_s;
    dbl valid_r;
    dbl valid_s_dx;
    dbl valid_r_dx;
    if (
            (t_s[0] > -1e-10) and (t_s[0] < (1.0 + 1e-10)) and
            (t_r[1] > -1e-10) and (t_r[1] < (1.0 + 1e-10))
       )
    {
        valid_s = t_s[0];
        valid_r = t_r[1];
        valid_s_dx = s.dx(x, y)[0];
        valid_r_dx = r.dx(x, y)[1];
    }
    else
    {
        valid_s = t_s[1];
        valid_r = t_r[0];
        valid_s_dx = s.dx(x, y)[1];
        valid_r_dx = r.dx(x, y)[0];
    };

    return f_.dx(valid_s, valid_r, valid_s_dx, valid_r_dx);
    // return f_.dx(s(x, y)[0], r(x, y)[1], s.dx(x, y)[0], r.dx(x, y)[1]);
};

dbl Scalar4PointsFunc<2>::dy (cdbl x, cdbl y) const
{
    const auto t_s = s(x, y);
    const auto t_r = r(x, y);
    dbl valid_s;
    dbl valid_r;
    dbl valid_s_dy;
    dbl valid_r_dy;
    if (
            (t_s[0] > -1e-10) and (t_s[0] < (1.0 + 1e-10)) and
            (t_r[1] > -1e-10) and (t_r[1] < (1.0 + 1e-10))
       )
    {
        valid_s = t_s[0];
        valid_r = t_r[1];
        valid_s_dy = s.dy(x, y)[0];
        valid_r_dy = r.dy(x, y)[1];
    }
    else
    {
        valid_s = t_s[1];
        valid_r = t_r[0];
        valid_s_dy = s.dy(x, y)[1];
        valid_r_dy = r.dy(x, y)[0];
    };

    return f_.dy(valid_s, valid_r, valid_s_dy, valid_r_dy);
    // return f_.dy(s(x, y)[0], r(x, y)[1], s.dy(x, y)[0], r.dy(x, y)[1]);
};

dbl Scalar4PointsFunc<2>::operator() (const prmt::Point<2> &p) const
{
    return operator()(p.x(), p.y());
};

dbl Scalar4PointsFunc<2>::dx (const prmt::Point<2> &p) const
{
    return dx(p.x(), p.y());
};

dbl Scalar4PointsFunc<2>::dy (const prmt::Point<2> &p) const
{
    return dy(p.x(), p.y());
};
///////////////////////////

//////////////////////////LocFunc
Scalar4PointsFunc<2>::LocFunc::LocFunc (arr<dbl, 4> glob_val) :
    a (glob_val[0]),
    b (glob_val[1] - glob_val[0]),
    c (glob_val[3] - glob_val[0]),
    d ((glob_val[0] + glob_val[2]) - (glob_val[1] + glob_val[3]))
{};

dbl Scalar4PointsFunc<2>::LocFunc::operator() (cdbl s, cdbl r) const
{
    return a + b * s + c * r + d * s * r;
};

dbl Scalar4PointsFunc<2>::LocFunc::dx (cdbl s, cdbl r, cdbl s_dx, cdbl r_dx) const
{
    return b * s_dx + c * r_dx + d * (s_dx * r + s * r_dx);
};

dbl Scalar4PointsFunc<2>::LocFunc::dy (cdbl s, cdbl r, cdbl s_dy, cdbl r_dy) const
{
    return b * s_dy + c * r_dy + d * (s_dy * r + s * r_dy);
};
/////////////////////////

///////////////////////////LocVar
Scalar4PointsFunc<2>::LocVar::LocVar (cdbl _A, const LinFunc _B, const LinFunc _C) :
    A (_A),
    B (_B),
    C (_C)
{};

arr<dbl,2> Scalar4PointsFunc<2>::LocVar::operator() (cdbl x, cdbl y) const
{
    if (std::abs(A) > 1e-10)
    {
        return quad_roots(x, y);
    }
    else
    {
    // printf("C=%f B=%f res=%f\n",
    //         C(x,y), B(x, y), C(x, y) / B(x, y));
        cdbl res = C(x, y) / B(x, y);
        return arr<dbl, 2>{res, res};
    };
};

arr<dbl,2> Scalar4PointsFunc<2>::LocVar::dx (cdbl x, cdbl y) const
{
    if (std::abs(A) > 1e-10)
    {
        cdbl m1 = B.dx(x, y) / (2.0 * A);
        cdbl numerator = B(x, y) * B.dx(x, y) - 2.0 * A * C.dx(x, y);
        cdbl denominator = std::sqrt(B(x, y) * B(x, y) - 4.0 * A * C(x, y)) * 2.0 * A;
        cdbl m2 = numerator / denominator;
        return arr<dbl, 2>{m1 - m2, m1 + m2};
    }
    else
    {
        cdbl res = (C.dx(x, y) * B(x, y) -  C(x, y) * B.dx(x, y)) / (B(x, y) * B(x, y));
        return arr<dbl, 2>{res, res};
    };
};

arr<dbl,2> Scalar4PointsFunc<2>::LocVar::dy (cdbl x, cdbl y) const
{
    if (std::abs(A) > 1e-10)
    {
        cdbl m1 = B.dy(x, y) / (2.0 * A);
        cdbl numerator = B(x, y) * B.dy(x, y) - 2.0 * A * C.dy(x, y);
        cdbl denominator = std::sqrt(B(x, y) * B(x, y) - 4.0 * A * C(x, y)) * 2.0 * A;
        cdbl m2 = numerator / denominator;
        return arr<dbl, 2>{m1 - m2, m1 + m2};
    }
    else
    {
        dbl res = (C.dy(x, y) * B(x, y) -  C(x, y) * B.dy(x, y)) / (B(x, y) * B(x, y));
        return arr<dbl, 2>{res, res};
    };
};

arr<dbl,2> Scalar4PointsFunc<2>::LocVar::quad_roots (cdbl x, cdbl y) const
{
    // puts("!!!");
    cdbl sq_root = std::sqrt(B(x, y) * B(x, y) - 4.0 * A * C(x, y));
    // printf("C=%f B=%f sq_root=%f under_root=%f res=%f\n",
    //         C(x,y), B(x, y), sq_root, 
    //         B(x, y) * B(x, y) - 4.0 * A * C(x, y),
    //         (B(x, y) + sq_root) / (2.0 * A));
    
    return arr<dbl, 2>{
        (B(x, y) - sq_root) / (2.0 * A),
        (B(x, y) + sq_root) / (2.0 * A)};
};
//////////////////////////

//////////////////////////LinFunc
Scalar4PointsFunc<2>::LinFunc::LinFunc (cdbl _a, cdbl _b, cdbl _c) :
    a (_a),
    b (_b),
    c (_c)
{
    // printf("%f %f %f\n", _a, _b, _c);
};

dbl Scalar4PointsFunc<2>::LinFunc::operator() (cdbl x, cdbl y) const
{
    return a + b * x + c * y;
};

dbl Scalar4PointsFunc<2>::LinFunc::dx (cdbl x, cdbl y) const
{
    return b;
};

dbl Scalar4PointsFunc<2>::LinFunc::dy (cdbl x, cdbl y) const
{
    return c;
};
/////////////////////////

#endif
