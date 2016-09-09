#ifndef geometric_tools_def
#define geometric_tools_def 1

#include "../../../../../prmt_sintactic_addition/prmt_sintactic_addition.h"
#include "../point/point.h"

//! Геометрический инструментарий
namespace GTools
{

    //! Отрезок
    /*! Формирует кусочнолинейную кривую curve в виде отрезка. \n
     *  Число вершин полученной кривой num_points + 1. \n
     *  first и second - начальная и конечная вершины отрезка (сама конечная
     *  вершина не попадает в состав кривой)
     */
    void give_line_without_end_point(
            vec<prmt::Point<2>> &curve,
            cst num_points,
            prmt::Point<2> first,
            prmt::Point<2> second)
    {
        dbl dx = (second.x() - first.x()) / num_points;
        dbl dy = (second.y() - first.y()) / num_points;
        dbl x = first.x();
        dbl y = first.y();
        FOR_I(0, num_points - 0)
        {
            curve .push_back (prmt::Point<2>(x, y)); 
            x += dx;
            y += dy;
        };
    };

    //! Прямоугольник
    /*! Формирует кусочнолинейную кривую curve в виде прямоугольника. \n
     *  Число вершин полученной кривой (num_points_on_edge + 1) * 4. \n
     *  first и second - две противоположных вершины, первая в точке (x0, y0),
     *  вторая в точке (x0+a, y0+b); где a - ширина прямоугольника, b - его
     *  высота.
     */
    void give_rectangle(
            vec<prmt::Point<2>> &curve,
            cst num_points_on_edge,
            prmt::Point<2> first,
            prmt::Point<2> second)
    {
        give_line_without_end_point(curve, num_points_on_edge,
                first,
                prmt::Point<2>(first.x(), second.y()));

        give_line_without_end_point(curve, num_points_on_edge,
                prmt::Point<2>(first.x(), second.y()),
                second);

        give_line_without_end_point(curve, num_points_on_edge,
                second,
                prmt::Point<2>(second.x(), first.y()));

        give_line_without_end_point(curve, num_points_on_edge,
                prmt::Point<2>(second.x(), first.y()),
                first);

    };

    //! Прямоугольник с пронумерованными границами (для задания на них граничных условий)
    /*! Формирует кусочнолинейную кривую curve в виде прямоугольника. \n
     *  Число вершин полученной кривой (num_points_on_edge + 1) * 4. \n
     *  first и second - две противоположных вершины, первая в точке (x0, y0),
     *  вторая в точке (x0+a, y0+b); где a - ширина прямоугольника, b - его
     *  высота. \n
     *  piace_id - выходной параметр; нумерация кусков. \n
     *  side_id - входной параметр; номера четырёх границ прямоугольника.
     */
    void give_rectangle_with_border_condition(
            vec<prmt::Point<2>> &curve,
            vec<st> &piace_id,
            const arr<st, 4> side_id,
            cst num_points_on_edge,
            const prmt::Point<2> first,
            const prmt::Point<2> second)
    {
        give_line_without_end_point(curve, num_points_on_edge,
                first,
                prmt::Point<2>(first.x(), second.y()));

        give_line_without_end_point(curve, num_points_on_edge,
                prmt::Point<2>(first.x(), second.y()),
                second);

        give_line_without_end_point(curve, num_points_on_edge,
                second,
                prmt::Point<2>(second.x(), first.y()));

        give_line_without_end_point(curve, num_points_on_edge,
                prmt::Point<2>(second.x(), first.y()),
                first);

        cst n_edge_on_border = curve.size() / 4;
        // printf("type %d\n", n_edge_on_border);
        piace_id.resize(curve.size());

        FOR(i, 0, 4)
            FOR(j, 0 + n_edge_on_border * i, n_edge_on_border + n_edge_on_border * i)
            piace_id[j] = side_id[i];
    };

};

#endif
