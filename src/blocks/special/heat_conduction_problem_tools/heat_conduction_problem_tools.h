#ifndef heat_conduction_problem_def
#define heat_conduction_problem_def 1

#include <fstream>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/numerics/data_out.h>

#include "../../general/4_points_function/4_points_function.h"

//! Инструменты для решения задачи теплопроводности
namespace HCPTools
{
    //! Задать тензор теплопроводности 
    template<u8 dim>
    void set_thermal_conductivity(
            arr<arr<vec<dbl>, dim>, dim> &tensor, const arr<arr<vec<dbl>, dim>, dim> &coef)
    {
        for (st i = 0; i < dim; ++i)
        {
            for (st j = 0; j < dim; ++j)
            {
                tensor[i][j] .resize (coef[i][j].size());
                for (st k = 0; k < coef[i][j].size(); ++k)
                {
                    tensor[i][j][k] = coef[i][j][k];
                };
            };
        };
    };

    //! Распечатать в файл температуру
    template<u8 dim>
    void print_temperature (const dealii::Vector<dbl> &temperature, 
                            const dealii::DoFHandler<dim> &dof_handler,
                            const str file_name,
                            const dealii::DataOutBase::OutputFormat output_format = dealii::DataOutBase::gnuplot)
    {
        dealii::DataOut<dim> data_out;
        data_out.attach_dof_handler (dof_handler);
        data_out.add_data_vector (temperature, "temperature");
        data_out.build_patches ();

        auto name = file_name;

        std::ofstream output (name);
        data_out.write (output, output_format);
    };

    //! Распечатать в файл 2d срез температуры для 3d случая
    void print_temperature_slice (const dealii::Vector<dbl> &temperature, 
                            const dealii::DoFHandler<3> &dof_handler,
                            const str file_name,
                            cst ort,
                            cdbl slice_coor)
    {
        // struct LocPoint
        // {
        //     bool flag = false;
        //     dealii::Point<2> coor;
        //     dbl val;
        // };
        //     
        // cst num_point_on_slice = (2_pow(dof_handler.get_tria().n_levels() - 1) + 1)_pow(2);
        // vec<LocPoint> slice (num_point_on_slice);

        FILE* f_out;
        f_out = fopen (file_name.c_str(), "w");
        for (auto cell = dof_handler.begin_active (); cell != dof_handler.end (); ++cell)
        {
            for (st i = 0; i < dealii::GeometryInfo<3>::vertices_per_cell; ++i)
            {
                if (std::abs(cell->vertex(i)(ort) - slice_coor) < 1e-10)
                {
                    // dealii::Point<2> point;
                    // for (st j = 0; j < 3; ++j)
                    // {
                    //     if (j != ort)
                    //     {
                    //         point(j) = cell->vetrex(i)(j);
                    //     };
                    // };
                    // cdbl val = temperature(cell->vertex_dof_index(i, 0));
                    // fprintf(f_out, "%f %f %f\n", point(0), point(1), val);
                    fprintf(f_out, "%f %f %f %f\n",
                            cell->vertex(i)(0),
                            cell->vertex(i)(1),
                            cell->vertex(i)(2),
                            temperature(cell->vertex_dof_index(i, 0)));
                };
            };
        };
        fclose(f_out);
    };

    //! Распечатать в файл 1d срез температуры для 2d случая
    void print_temperature_slice (const dealii::Vector<dbl> &temperature, 
                            const dealii::DoFHandler<2> &dof_handler,
                            const str file_name,
                            cst ort,
                            cdbl slice_coor)
    {
        FILE* f_out;
        f_out = fopen (file_name.c_str(), "w");
        for (auto cell = dof_handler.begin_active (); cell != dof_handler.end (); ++cell)
        {
            for (st i = 0; i < dealii::GeometryInfo<3>::vertices_per_cell; ++i)
            {
                if (std::abs(cell->vertex(i)(ort) - slice_coor) < 1e-10)
                {
                    fprintf(f_out, "%f %f %f\n",
                            cell->vertex(i)(0),
                            cell->vertex(i)(1),
                            temperature(cell->vertex_dof_index(i, 0)));
                };
            };
        };
        fclose(f_out);
    };

    //! Распечатать в файл тепловые потоки
    template<u8 dim>
    void print_heat_conductions (const dealii::Vector<dbl> &temperature, 
                                 vec<arr<arr<dbl, 3>, 3>> &coef,
                                 const Domain<2> &domain,
                                 const str file_name)
    {
    };

    template<>
    void print_heat_conductions<2> (const dealii::Vector<dbl> &temperature, 
                                    vec<arr<arr<dbl, 3>, 3>> &coef,
                                    const Domain<2> &domain,
                                    const str file_name)
    {
        cu32 dofs_per_cell = domain.dof_handler.get_fe().dofs_per_cell;

        arr<dealii::Vector<dbl>, 2> grad;
        arr<dealii::Vector<dbl>, 2> hc;
        
        grad[0] .reinit (domain.dof_handler.n_dofs());
        grad[1] .reinit (domain.dof_handler.n_dofs());

        hc[0] .reinit (domain.dof_handler.n_dofs());
        hc[1] .reinit (domain.dof_handler.n_dofs());

        vec<st> N(domain.dof_handler.n_dofs());

        for (
                auto cell = domain.dof_handler.begin_active(); 
                cell     != domain.dof_handler.end(); 
                ++cell
            )
        {
            /* Точки 3 и 2 переставленны местами, потому что в диле у них
             * порядок зигзагом, а мне надо по кругу
            */
            arr<prmt::Point<2>, 4> points = {
                prmt::Point<2>(cell->vertex(0)(0), cell->vertex(0)(1)),
                prmt::Point<2>(cell->vertex(1)(0), cell->vertex(1)(1)),
                prmt::Point<2>(cell->vertex(3)(0), cell->vertex(3)(1)),
                prmt::Point<2>(cell->vertex(2)(0), cell->vertex(2)(1))};
            arr<dbl, 4> values = {
                temperature(cell->vertex_dof_index (0, 0)),
                temperature(cell->vertex_dof_index (1, 0)),
                temperature(cell->vertex_dof_index (3, 0)),
                temperature(cell->vertex_dof_index (2, 0))};

            Scalar4PointsFunc<2> function_on_cell(points, values);

            for (st i = 0; i < dofs_per_cell; ++i)
            {
                auto indx = cell->vertex_dof_index(i, 0);

                grad[0][indx] += function_on_cell.dx(cell->vertex(i));
                grad[1][indx] += function_on_cell.dy(cell->vertex(i));
                ++(N[indx]); 

                hc[0][indx] = 
                    coef[cell->material_id()][0][0] * grad[0][indx] +
                    coef[cell->material_id()][0][1] * grad[1][indx];
                hc[1][indx] = 
                    coef[cell->material_id()][1][0] * grad[0][indx] +
                    coef[cell->material_id()][1][1] * grad[1][indx];

                // hc[0][indx] = grad[0][indx];
                //     // coef[cell->material_id()][0][0] * grad[0][indx] +
                //     // coef[cell->material_id()][0][1] * grad[1][indx];
                // hc[1][indx] = grad[1][indx];
                //     // coef[cell->material_id()][1][0] * grad[0][indx] +
                //     // coef[cell->material_id()][1][1] * grad[1][indx];
            };
        };
        for (st i = 0; i < N.size(); ++i)
        {
            hc[0][i] /= N[i];
            hc[1][i] /= N[i];
        };

        {
            dealii::DataOut<2> data_out;
            data_out.attach_dof_handler (domain.dof_handler);
            data_out.add_data_vector (hc[0], "hc_dx");
            data_out.add_data_vector (hc[1], "hc_dy");
            // data_out.add_data_vector (grad[0], "grad_dx");
            // data_out.add_data_vector (grad[1], "drad_dy");
            data_out.build_patches ();

            auto name = file_name;
            // name += ".gpd";

            std::ofstream output (name);
            data_out.write_gnuplot (output);
        };
    }

    //! Распечатать в файл тепловые градиенты
    template<u8 dim>
    void print_heat_gradient (const dealii::Vector<dbl> &temperature, 
                                 const Domain<2> &domain,
                                 const str file_name)
    {
    };

    template<>
    void print_heat_gradient<2> (const dealii::Vector<dbl> &temperature, 
                                    const Domain<2> &domain,
                                    const str file_name)
    {
        cu32 dofs_per_cell = domain.dof_handler.get_fe().dofs_per_cell;

        arr<dealii::Vector<dbl>, 2> grad;
        
        grad[0] .reinit (domain.dof_handler.n_dofs());
        grad[1] .reinit (domain.dof_handler.n_dofs());

        // dealii::Vector<dbl> T(domain.dof_handler.n_dofs());
        // vec<bool> flg(domain.dof_handler.n_dofs());
        // for (st i = 0; i < flg.size(); ++i)
        // {
        //     flg[i] = false;
        // };

        // vec<dbl> gradx;
        // vec<prmt::Point<2>> ps;

        vec<st> N(domain.dof_handler.n_dofs());

        // FILE *F;
        // F = fopen("test_iso_3.gpd", "w");
        for (
                auto cell = domain.dof_handler.begin_active(); 
                cell     != domain.dof_handler.end(); 
                ++cell
            )
        {
            /* Точки 3 и 2 переставленны местами, потому что в диле у них
             * порядок зигзагом, а мне надо по кругу
            */
            arr<prmt::Point<2>, 4> points = {
                prmt::Point<2>(cell->vertex(0)(0), cell->vertex(0)(1)),
                prmt::Point<2>(cell->vertex(1)(0), cell->vertex(1)(1)),
                prmt::Point<2>(cell->vertex(3)(0), cell->vertex(3)(1)),
                prmt::Point<2>(cell->vertex(2)(0), cell->vertex(2)(1))};
            arr<dbl, 4> values = {
                temperature(cell->vertex_dof_index (0, 0)),
                temperature(cell->vertex_dof_index (1, 0)),
                temperature(cell->vertex_dof_index (3, 0)),
                temperature(cell->vertex_dof_index (2, 0))};

            Scalar4PointsFunc<2> function_on_cell(points, values);

            // dbl midl_x = (cell->vertex(0)(0) + cell->vertex(1)(0) + cell->vertex(2)(0) + cell->vertex(3)(0)) / 4.0;
            // dbl midl_y = (cell->vertex(0)(1) + cell->vertex(1)(1) + cell->vertex(2)(1) + cell->vertex(3)(1)) / 4.0;
            //
            // gradx .push_back (function_on_cell.dy(midl_x, midl_y));
            // ps .push_back (prmt::Point<2>(midl_x, midl_y));
            // st n = 5;
            // for (st i = 0; i < n+1; ++i)
            // {
            //     cdbl s = 1.0 / n * i;
            //     for (st j = 0; j < n+1; ++j)
            //     {
            //         cdbl r = 1.0 / n * j;
            //         cdbl x = function_on_cell.a_x +
            //                  function_on_cell.b_x*s +
            //                  function_on_cell.c_x*r +
            //                  function_on_cell.d_x*s*r;
            //         cdbl y = function_on_cell.a_y +
            //                  function_on_cell.b_y*s +
            //                  function_on_cell.c_y*r +
            //                  function_on_cell.d_y*s*r;
            //         cdbl f = function_on_cell.f_(s, r);
            //         cdbl dx = function_on_cell.f_.dx(
            //                   x,
            //                   y,
            //                   function_on_cell.s.dx(x, y)[0],
            //                   function_on_cell.r.dx(x, y)[1]);
            //         cdbl dy = function_on_cell.f_.dy(
            //                   x,
            //                   y,
            //                   function_on_cell.s.dy(x, y)[0],
            //                   function_on_cell.r.dy(x, y)[1]);
            //         fprintf(F, "%f %f %f %f %f\n", x, y, f, dx, dy);
            //     };
            // };
            // if (std::abs(function_on_cell.dy(midl_x, midl_y)) > 1e-5)
            //     fprintf(F, "%f %f %f %f %f %f %f %f %f\n",
            //             cell->vertex(0)(0), cell->vertex(0)(1),
            //             cell->vertex(1)(0), cell->vertex(1)(1),
            //             cell->vertex(3)(0), cell->vertex(3)(1),
            //             cell->vertex(2)(0), cell->vertex(2)(1),
            //             function_on_cell.dy(midl_x, midl_y));

                // auto indx = cell->vertex_dof_index(0, 0);

                // grad[0][indx] = function_on_cell.dx(cell->vertex(0));
                // grad[1][indx] = function_on_cell.dy(cell->vertex(0));
            for (st i = 0; i < dofs_per_cell; ++i)
            {
                auto indx = cell->vertex_dof_index(i, 0);

                grad[0][indx] += function_on_cell.dx(cell->vertex(i));
                grad[1][indx] += function_on_cell.dy(cell->vertex(i));
                ++(N[indx]); 

                // if (!flg[indx])
                // {
                // T[indx] = 
                //     // temperature(cell->vertex_dof_index (i, 0));
                //     // function_on_cell(cell->vertex(i));
                //     // function_on_cell(prmt::Point<2>(cell->vertex(i)(0), cell->vertex(i)(1)));
                //     // cell->vertex(i)(0) * cell->vertex(i)(0);
                //     function_on_cell(prmt::Point<2>(cell->vertex(i)(0), cell->vertex(i)(1)))-
                // temperature(cell->vertex_dof_index (i, 0));
                // // if (std::abs(T[indx]) > 1e-5)
                // //     fprintf(F, "%f %f %f %f %f %f %f %f %f %f %f %f %f %d\n",
                // //             cell->vertex(0)(0), cell->vertex(0)(1),
                // //             cell->vertex(1)(0), cell->vertex(1)(1),
                // //             cell->vertex(3)(0), cell->vertex(3)(1),
                // //             cell->vertex(2)(0), cell->vertex(2)(1),
                // //             temperature(cell->vertex_dof_index (0, 0)),
                // //             temperature(cell->vertex_dof_index (1, 0)),
                // //             temperature(cell->vertex_dof_index (3, 0)),
                // //             temperature(cell->vertex_dof_index (2, 0)),
                // //             T[indx],
                // //             i
                // //             );
                // flg[indx] = true;
                // };
            };
        };
        // fclose(F);
        for (st i = 0; i < N.size(); ++i)
        {
            // printf("%d\n", N[i]);
            grad[0][i] /= N[i];
            grad[1][i] /= N[i];
            // T[i] /= N[i];
        };

        {
            dealii::DataOut<2> data_out;
            data_out.attach_dof_handler (domain.dof_handler);
            data_out.add_data_vector (grad[0], "grad_dx");
            data_out.add_data_vector (grad[1], "drad_dy");
            // data_out.add_data_vector (T, "T");
            data_out.build_patches ();

            auto name = file_name;
            // name += ".gpd";

            std::ofstream output (name);
            data_out.write_gnuplot (output);
        };
        // FILE *F1;
        // F1 = fopen("gradx.gpd", "w");
        // for (st i = 0; i < ps.size(); ++i)
        // {
        //     fprintf(F1, "%lf %lf %lf\n", ps[i].x(), ps[i].y(), gradx[i]);
        // };
        // fclose(F1);
    }
};

#endif
