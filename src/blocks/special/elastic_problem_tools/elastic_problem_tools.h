#ifndef ELASTIC_PROBLEM_TOOLS
#define ELASTIC_PROBLEM_TOOLS 1

#include <fstream>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/numerics/data_out.h>
// #include "../problem_on_cell/calculate_meta_coefficients/calculate_meta_coefficients.h"

// #include "../../general/4_points_function/4_points_function.h"

//! Инструменты для решения задачи упругости
namespace EPTools
{
    //! Задать изотропный тэнзор упругости 

    struct set_isotropic_elascity
    {
        cdbl yung;
        cdbl puasson;

        void operator() (ATools::FourthOrderTensor &coef)
        {
            enum {x, y, z};
            for (st i = 0; i < 3; ++i)
            {
                for (st j = 0; j < 3; ++j)
                {
                    for (st k = 0; k < 3; ++k)
                    {
                        for (st l = 0; l < 3; ++l)
                        {
                            coef[i][j][k][l] = 0.0;
                        };
                    };
                };
            };

            const double lambda = 
                (puasson * yung) / ((1 + puasson) * (1 - 2 * puasson));
            const double mu     = 
                yung / (2 * (1 + puasson));

            coef[x][x][x][x] = lambda + 2 * mu;
            coef[y][y][y][y] = lambda + 2 * mu;
            coef[z][z][z][z] = lambda + 2 * mu;

            coef[x][x][y][y] = lambda;
            coef[y][y][x][x] = lambda;

            coef[x][y][x][y] = mu;
            coef[y][x][y][x] = mu;
            coef[x][y][y][x] = mu;
            coef[y][x][x][y] = mu;

            coef[x][x][z][z] = lambda;
            coef[z][z][x][x] = lambda;

            coef[x][z][x][z] = mu;
            coef[z][x][z][x] = mu;
            coef[x][z][z][x] = mu;
            coef[z][x][x][z] = mu;

            coef[z][z][y][y] = lambda;
            coef[y][y][z][z] = lambda;

            coef[z][y][z][y] = mu;
            coef[y][z][y][z] = mu;
            coef[z][y][y][z] = mu;
            coef[y][z][z][y] = mu;
        };
    };

    //! Задать ортотропный тэнзор упругости 

    void set_ortotropic_elascity (
            const arr<arr<dbl, 3>, 3> E, 
            const arr<dbl, 3> G, 
            ATools::FourthOrderTensor &C)
    {
        enum {x, y, z};
        for (st i = 0; i < 3; ++i)
        {
            for (st j = 0; j < 3; ++j)
            {
                for (st k = 0; k < 3; ++k)
                {
                    for (st l = 0; l < 3; ++l)
                    {
                        C[i][j][k][l] = 0.0;
                    };
                };
            };
        };

        cdbl c = E[x][y] * E[y][x] + E[z][y] * E[y][z] + E[x][z] * E[z][x];
        cdbl d = E[x][y] * E[y][z] * E[z][x] + E[x][z] * E[z][y] * E[y][x];
        cdbl sm = 1.0 / (1 - c - d);

        auto a = [&E] (cst i){
            cst m = (i + 1) % 3;
            cst n = (i + 2) % 3;
            return (1.0 - E[m][n] * E[n][m]);
        };
        auto b = [&E] (cst i, cst j){
            for (st k = 0; k < 3; ++k)
            {
                if ((k != i) and (k != j))
                    return E[i][j] + E[i][k] * E[k][j];
            };
        };

        for (st i = 0; i < 3; ++i)
        {
            for (st j = 0; j < 3; ++j)
            {
                if (i == j)
                {
                    C[i][i][j][j] = E[i][i] * a(i) * sm;
                }
                else
                {
                    C[i][i][j][j] = E[i][i] * b(i, j) * sm;
                    C[i][j][i][j] = G[i + j - 1];
                    C[i][j][j][i] = G[i + j - 1];
                };
                // C[j][j][i][i] = C[i][i][j][j];
                // C[j][i][j][i] = C[i][j][i][j];

            };
        };
    };

    //! Распечатать в файл перемещения
    template<u8 dim>
        void print_move (const dealii::Vector<dbl> &move, 
                const dealii::DoFHandler<dim> &dof_handler,
                const str file_name,
                const dealii::DataOutBase::OutputFormat output_format = dealii::DataOutBase::gnuplot)
        {
            dealii::DataOut<dim> data_out;
            data_out.attach_dof_handler (dof_handler);
            vec<str> solution_names;
            switch (dim)
            {
                case 1:
                    solution_names.push_back ("move_x");
                    break;
                case 2:
                    solution_names.push_back ("move_x");
                    solution_names.push_back ("move_y");
                    break;
                case 3:
                    solution_names.push_back ("move_x");
                    solution_names.push_back ("move_y");
                    solution_names.push_back ("move_z");
                    break;
            };
            data_out.add_data_vector (move, solution_names);
            data_out.build_patches ();

            auto name = file_name;
            // name += ".gpd";

            std::ofstream output (name);
            data_out.write (output, output_format);
        }

    //! Распечатать в файл 2d срез перемещений для 3d случая
    void print_move_slice (const dealii::Vector<dbl> &move, 
            const dealii::DoFHandler<3> &dof_handler,
            const str file_name,
            cst ort,
            cdbl slice_coor)
    {
        // dbl int_x = 0.0;
        // dbl int_y = 0.0;
        // dbl integ = 0.0;
        FILE* f_out;
        // FILE* f_out_x;
        // FILE* f_out_y;
        f_out = fopen (file_name.c_str(), "w");
        // f_out_x = fopen ("approx_2_x", "w");
        // f_out_y = fopen ("approx_2_y", "w");
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
                    fprintf(f_out, "%f %f %f %f %f %f\n",
                            cell->vertex(i)(0),
                            cell->vertex(i)(1),
                            cell->vertex(i)(2),
                            move(cell->vertex_dof_index(i, 0)),
                            move(cell->vertex_dof_index(i, 1)),
                            move(cell->vertex_dof_index(i, 2))
                           );
                    // if (std::abs(cell->vertex(i)(1) - 1.0) < 1e-10)
                    // {
                    //     // fprintf(f_out_x, "%f %f\n",
                    //     //         cell->vertex(i)(0), move(cell->vertex_dof_index(i, 0)));
                    //     int_x += move(cell->vertex_dof_index(i, 0));
                    // }
                    // if (std::abs(cell->vertex(i)(0) - 1.0) < 1e-10)
                    // {
                    //     // fprintf(f_out_y, "%f %f\n",
                    //     //         cell->vertex(i)(1), move(cell->vertex_dof_index(i, 0)));
                    //     int_y += move(cell->vertex_dof_index(i, 0));
                    // };
                    // integ += move(cell->vertex_dof_index(i, 0));
                };
            };
        };
        fclose(f_out);
        // fclose(f_out_x);
        // fclose(f_out_y);
        // printf("INTEGRALS %f %f %f\n", int_x, int_y, integ);
    };
    //
    // void print_cell_stress_slice ( 
    //         const arr<i32, 3> approximation,
    //         cst nu,
    //         const OnCell::ArrayWithAccessToVector<arr<dealii::Vector<dbl>, 3>> &cell_func
    //         const dealii::DoFHandler<3> &dof_handler,
    //         const str file_name,
    //         cst ort,
    //         cdbl slice_coor)
    // {
    //     dealii::Vector<dbl> stress;
    //     
    //     grad .reinit (dof_handler.n_dofs());
    //
    //     vec<st> N(dof_handler.n_dofs());
    //
    //     for (
    //             auto cell = dof_handler.begin_active(); 
    //             cell     != dof_handler.end(); 
    //             ++cell
    //         )
    //     {
    //         cst material_id = cell->material_id();
    //         cell->get_dof_indices (global_dof_indices);
    //         /* Точки 3 и 2 переставленны местами, потому что в диле у них
    //          * порядок зигзагом, а мне надо по кругу
    //          */
    //         arr<prmt::Point<2>, 4> points = {
    //             prmt::Point<2>(cell->vertex(0)(0), cell->vertex(0)(1)),
    //             prmt::Point<2>(cell->vertex(1)(0), cell->vertex(1)(1)),
    //             prmt::Point<2>(cell->vertex(3)(0), cell->vertex(3)(1)),
    //             prmt::Point<2>(cell->vertex(2)(0), cell->vertex(2)(1))};
    //
    //         for (st alpha = 0; alpha < 3; ++alpha)
    //         {
    //             for (st beta = 0; beta < 3; ++beta)
    //                 // for (st phi = 0; phi < 3; ++phi)
    //                 // cst beta = nu;
    //             {
    //                 for (st m = 0; m < dofs_per_cell; ++m)
    //                 {
    //                     // cst beta = m % 3;
    //                     cst phi = m % 3;
    //                     for (st psi = 0; psi < 3; ++psi)
    //                     {
    //                         // arr<i32, 3> k_psi = {k[x], k[y], k[z]}; // k - э_psi
    //                         // bool cell_k_psi_exist = true;
    //                         // if (k[psi] == 0)
    //                         //     cell_k_psi_exist = false;
    //                         // else
    //                         //     k_psi[psi]--;
    //
    //                         arr<i32, 3> k_beta = {k[x], k[y], k[z]}; // k - э_beta
    //                         bool cell_k_beta_exist = true;
    //                         if (k[beta] == 0)
    //                             cell_k_beta_exist = false;
    //                         else
    //                             k_beta[beta]--;
    //
    //                         arr<i32, 3> k_psi_beta = {k_beta[x], k_beta[y], k_beta[z]}; // k - э_psi - э_beta
    //                         bool cell_k_psi_beta_exist = true;
    //                         if (cell_k_beta_exist)
    //                         {
    //                             if (k_beta[psi] == 0)
    //                                 cell_k_psi_beta_exist = false;
    //                             else
    //                                 k_psi_beta[psi]--;
    //                         }
    //                         else
    //                             cell_k_psi_beta_exist = false;
    //
    //                         dbl cell_k_psi      = 0.0;
    //                         dbl cell_k_beta     = 0.0;
    //                         dbl cell_k_psi_beta = 0.0;
    //
    //                         // if (cell_k_psi_exist)
    //                         //     cell_k_psi = cell_func[k_psi][nu][global_dof_indices[m]];
    //                         // else
    //                         //     cell_k_psi = 0.0;
    //                         //
    //                         // printf("%d %d %d  %d %d %d  %d %d  %d %d\n", 
    //                         //         k_beta[0], k_beta[1], k_beta[2],
    //                         //         k_psi_beta[0], k_psi_beta[1], k_psi_beta[2],
    //                         //         cell_k_beta_exist, cell_k_psi_beta_exist,
    //                         //         beta, psi
    //                         //         );
    //                         if (cell_k_beta_exist)
    //                             cell_k_beta = cell_func[k_beta][nu][global_dof_indices[m]];
    //                         else
    //                             cell_k_beta = 0.0;
    //
    //                         if (cell_k_psi_beta_exist)
    //                             cell_k_psi_beta = cell_func[k_psi_beta][nu][global_dof_indices[m]];
    //                         else
    //                             cell_k_psi_beta = 0.0;
    //                         // // printf("%f\n", cell_k_psi_beta);
    //
    //                         res[alpha] += 
    //                             (
    //                              this->coef[material_id][alpha][beta][phi][psi] *
    //                              (
    //                               // 0.0 *
    //                               cell_k_beta *
    //                               this->fe_values.shape_grad (m, q_point)[psi] +
    //                               // 0.0 *
    //                               cell_k_psi_beta *
    //                               this->fe_values.shape_value (m, q_point)
    //                              )
    //                             ) *
    //                             this->fe_values.JxW(q_point);
    //                     };
    //                 };
    //             };
    //         };
    //     };
    //         for (st component = 0; component < 2; ++component)
    //         {
    //             arr<dbl, 4> values = {
    //                 move(cell->vertex_dof_index (0, component)),
    //                 move(cell->vertex_dof_index (1, component)),
    //                 move(cell->vertex_dof_index (3, component)),
    //                 move(cell->vertex_dof_index (2, component))};
    //
    //             Scalar4PointsFunc<2> function_on_cell(points, values);
    //
    //             for (st i = 0; i < 4; ++i)
    //             {
    //                 auto indx = cell->vertex_dof_index(i, component);
    //
    //                 if ((cell->material_id() == 0) or (cell->material_id() == 2))
    //                     // if (cell->material_id() == 0)
    //                 {
    //                     grad[0][indx] += function_on_cell.dx(cell->vertex(i));
    //                     grad[1][indx] += function_on_cell.dy(cell->vertex(i));
    //                     // grad[0][indx] = 1.0;
    //                     // grad[1][indx] = 1.0;
    //                 }
    //                 else
    //                 {
    //                     grad[0][indx] = 0.0;
    //                     grad[1][indx] = 0.0;
    //                 };
    //                 ++(N[indx]); 
    //             };
    //         };
    //     };
    //
    //     arr<dealii::Vector<dbl>, 2> stress;
    //     
    //     stress[0] .reinit (dof_handler.n_dofs());
    //     stress[1] .reinit (dof_handler.n_dofs());
    //
    //     for (
    //             auto cell = dof_handler.begin_active(); 
    //             cell     != dof_handler.end(); 
    //             ++cell
    //         )
    //     {
    //         for (st i = 0; i < 4; ++i)
    //         {
    //                 auto indx_1 = cell->vertex_dof_index(i, 0);
    //                 auto indx_2 = cell->vertex_dof_index(i, 1);
    //                 stress[0][indx_1] = 
    //                     C[0][0][0][0] * grad[0][indx_1] +
    //                     C[0][0][0][1] * grad[0][indx_2] +
    //                     C[0][0][1][0] * grad[1][indx_1] +
    //                     C[0][0][1][1] * grad[1][indx_2];
    //                 stress[0][indx_2] = 
    //                     C[0][1][0][0] * grad[0][indx_1] +
    //                     C[0][1][0][1] * grad[0][indx_2] +
    //                     C[0][1][1][0] * grad[1][indx_1] +
    //                     C[0][1][1][1] * grad[1][indx_2];
    //                 stress[1][indx_1] = 
    //                     C[1][0][0][0] * grad[0][indx_1] +
    //                     C[1][0][0][1] * grad[0][indx_2] +
    //                     C[1][0][1][0] * grad[1][indx_1] +
    //                     C[1][0][1][1] * grad[1][indx_2];
    //                 stress[1][indx_2] = 
    //                     C[1][1][0][0] * grad[0][indx_1] +
    //                     C[1][1][0][1] * grad[0][indx_2] +
    //                     C[1][1][1][0] * grad[1][indx_1] +
    //                     C[1][1][1][1] * grad[1][indx_2];
    //         };
    //     };
    //
    //     // {
    //     //     dealii::DataOut<2> data_out;
    //     //     data_out.attach_dof_handler (dof_handler);
    //     //     data_out.add_data_vector (stress[0], "grad_dx");
    //     //     data_out.add_data_vector (stress[1], "drad_dy");
    //     //     // data_out.add_data_vector (T, "T");
    //     //     data_out.build_patches ();
    //     //
    //     //     auto name = file_name;
    //     //     // name += ".gpd";
    //     //
    //     //     std::ofstream output (name);
    //     //     data_out.write_gnuplot (output);
    //     // };
    //     //
    //     // dbl Int = 0.0;
    //     // {
    //     //     cu32 dofs_per_cell = dof_handler.get_fe().dofs_per_cell;
    //     //
    //     //     dealii::QGauss<2>  quadrature_formula(2);
    //     //
    //     //     dealii::FEValues<2> fe_values (dof_handler.get_fe(), quadrature_formula,
    //     //             dealii::update_values | dealii::update_gradients | dealii::update_quadrature_points |
    //     //             dealii::update_JxW_values);
    //     //
    //     //     cst num_quad_points = quadrature_formula.size();
    //     //     for (
    //     //             auto cell = dof_handler.begin_active(); 
    //     //             cell     != dof_handler.end(); 
    //     //             ++cell
    //     //         )
    //     //     {
    //     //         fe_values .reinit (cell);
    //     //     for (st i = 0; i < 4; ++i)
    //     //     {
    //     //             auto indx_1 = cell->vertex_dof_index(i, 0);
    //     //             for (st q_point = 0; q_point < num_quad_points; ++q_point)
    //     //             {
    //     //             Int += 
    //     //                 stress[0](indx_1) *
    //     //                 // move(indx_1) *
    //     //                 // 1.0 *
    //     //                 // fe_values.shape_grad (i, q_point)[0] *
    //     //                 fe_values.shape_value (i, q_point) *
    //     //                 fe_values.JxW(q_point);
    //     //             };
    //     //     };
    //     //     };
    //     //
    //     // };
    // };

    // //! Распечатать в файл тепловые потоки
    // template<u8 dim>
    // void print_heat_conductions (const dealii::Vector<dbl> &temperature, 
    //                              arr<arr<vec<dbl>, 2>, 2> &coef,
    //                              const Domain<2> &domain,
    //                              const str file_name)
    // {
    // };

    // template<>
    // void print_heat_conductions<2> (const dealii::Vector<dbl> &temperature, 
    //                                 arr<arr<vec<dbl>, 2>, 2> &coef,
    //                                 const Domain<2> &domain,
    //                                 const str file_name)
    // {
    //     cu32 dofs_per_cell = domain.dof_handler.get_fe().dofs_per_cell;

    //     arr<dealii::Vector<dbl>, 2> grad;
    //     arr<dealii::Vector<dbl>, 2> hc;
    //     
    //     grad[0] .reinit (domain.dof_handler.n_dofs());
    //     grad[1] .reinit (domain.dof_handler.n_dofs());

    //     hc[0] .reinit (domain.dof_handler.n_dofs());
    //     hc[1] .reinit (domain.dof_handler.n_dofs());

    //     for (
    //             auto cell = domain.dof_handler.begin_active(); 
    //             cell     != domain.dof_handler.end(); 
    //             ++cell
    //         )
    //     {
    //         arr<prmt::Point<2>, 4> points = {
    //             cell->vertex(0), cell->vertex(1), 
    //             cell->vertex(2), cell->vertex(3)};
    //         arr<dbl, 4> values = {
    //             temperature(cell->vertex_dof_index (0, 0)),
    //             temperature(cell->vertex_dof_index (1, 0)),
    //             temperature(cell->vertex_dof_index (2, 0)),
    //             temperature(cell->vertex_dof_index (3, 0))};

    //         Scalar4PointsFunc<2> function_on_cell(points, values);

    //         for (st i = 0; i < dofs_per_cell; ++i)
    //         {
    //             auto indx = cell->vertex_dof_index(i, 0);

    //             grad[0][indx] = function_on_cell.dx(cell->vertex(i));
    //             grad[1][indx] = function_on_cell.dy(cell->vertex(i));

    //             hc[0][indx] = 
    //                 coef[0][0][cell->material_id()] * grad[0][indx] +
    //                 coef[0][1][cell->material_id()] * grad[1][indx];
    //             hc[1][indx] = 
    //                 coef[1][0][cell->material_id()] * grad[0][indx] +
    //                 coef[1][1][cell->material_id()] * grad[1][indx];
    //         };
    //     };

    //     {
    //         dealii::DataOut<2> data_out;
    //         data_out.attach_dof_handler (domain.dof_handler);
    //         data_out.add_data_vector (hc[0], "hc_dx");
    //         data_out.add_data_vector (hc[1], "hc_dy");
    //         data_out.build_patches ();

    //         auto name = file_name;
    //         name += ".gpd";

    //         std::ofstream output (name);
    //         data_out.write_gnuplot (output);
    //     };
    // }

    // //! Распечатать в файл тепловые градиенты
    // template<u8 dim>
    // void print_heat_gradient (const dealii::Vector<dbl> &temperature, 
    //                              arr<arr<vec<dbl>, 2>, 2> &coef,
    //                              const Domain<2> &domain,
    //                              const str file_name)
    // {
    // };

    // template<>
    // void print_heat_gradient<2> (const dealii::Vector<dbl> &temperature, 
    //                                 arr<arr<vec<dbl>, 2>, 2> &coef,
    //                                 const Domain<2> &domain,
    //                                 const str file_name)
    // {
    //     cu32 dofs_per_cell = domain.dof_handler.get_fe().dofs_per_cell;

    //     arr<dealii::Vector<dbl>, 2> grad;
    //     
    //     grad[0] .reinit (domain.dof_handler.n_dofs());
    //     grad[1] .reinit (domain.dof_handler.n_dofs());

    //     dealii::Vector<dbl> T(domain.dof_handler.n_dofs());
    //     vec<bool> flg(domain.dof_handler.n_dofs());
    //     for (st i = 0; i < flg.size(); ++i)
    //     {
    //         flg[i] = false;
    //     };

    //     vec<dbl> gradx;
    //     vec<prmt::Point<2>> ps;

    //     vec<st> N(domain.dof_handler.n_dofs());

    //     FILE *F;
    //     F = fopen("test_iso_3.gpd", "w");
    //     for (
    //             auto cell = domain.dof_handler.begin_active(); 
    //             cell     != domain.dof_handler.end(); 
    //             ++cell
    //         )
    //     {
    //         /* Точки 3 и 2 переставленны местами, потому что в диле у них
    //          * порядок зигзагом, а мне надо по кругу
    //         */
    //         arr<prmt::Point<2>, 4> points = {
    //             prmt::Point<2>(cell->vertex(0)(0), cell->vertex(0)(1)),
    //             prmt::Point<2>(cell->vertex(1)(0), cell->vertex(1)(1)),
    //             prmt::Point<2>(cell->vertex(3)(0), cell->vertex(3)(1)),
    //             prmt::Point<2>(cell->vertex(2)(0), cell->vertex(2)(1))};
    //         arr<dbl, 4> values = {
    //             temperature(cell->vertex_dof_index (0, 0)),
    //             temperature(cell->vertex_dof_index (1, 0)),
    //             temperature(cell->vertex_dof_index (3, 0)),
    //             temperature(cell->vertex_dof_index (2, 0))};

    //         Scalar4PointsFunc<2> function_on_cell(points, values);

    //         dbl midl_x = (cell->vertex(0)(0) + cell->vertex(1)(0) + cell->vertex(2)(0) + cell->vertex(3)(0)) / 4.0;
    //         dbl midl_y = (cell->vertex(0)(1) + cell->vertex(1)(1) + cell->vertex(2)(1) + cell->vertex(3)(1)) / 4.0;

    //         gradx .push_back (function_on_cell.dy(midl_x, midl_y));
    //         ps .push_back (prmt::Point<2>(midl_x, midl_y));
    //         st n = 5;
    //         for (st i = 0; i < n+1; ++i)
    //         {
    //             cdbl s = 1.0 / n * i;
    //             for (st j = 0; j < n+1; ++j)
    //             {
    //                 cdbl r = 1.0 / n * j;
    //                 cdbl x = function_on_cell.a_x +
    //                          function_on_cell.b_x*s +
    //                          function_on_cell.c_x*r +
    //                          function_on_cell.d_x*s*r;
    //                 cdbl y = function_on_cell.a_y +
    //                          function_on_cell.b_y*s +
    //                          function_on_cell.c_y*r +
    //                          function_on_cell.d_y*s*r;
    //                 cdbl f = function_on_cell.f_(s, r);
    //                 cdbl dx = function_on_cell.f_.dx(
    //                           x,
    //                           y,
    //                           function_on_cell.s.dx(x, y)[0],
    //                           function_on_cell.r.dx(x, y)[1]);
    //                 cdbl dy = function_on_cell.f_.dy(
    //                           x,
    //                           y,
    //                           function_on_cell.s.dy(x, y)[0],
    //                           function_on_cell.r.dy(x, y)[1]);
    //                 fprintf(F, "%f %f %f %f %f\n", x, y, f, dx, dy);
    //             };
    //         };
    //         // if (std::abs(function_on_cell.dy(midl_x, midl_y)) > 1e-5)
    //         //     fprintf(F, "%f %f %f %f %f %f %f %f %f\n",
    //         //             cell->vertex(0)(0), cell->vertex(0)(1),
    //         //             cell->vertex(1)(0), cell->vertex(1)(1),
    //         //             cell->vertex(3)(0), cell->vertex(3)(1),
    //         //             cell->vertex(2)(0), cell->vertex(2)(1),
    //         //             function_on_cell.dy(midl_x, midl_y));

    //             // auto indx = cell->vertex_dof_index(0, 0);

    //             // grad[0][indx] = function_on_cell.dx(cell->vertex(0));
    //             // grad[1][indx] = function_on_cell.dy(cell->vertex(0));
    //         for (st i = 0; i < dofs_per_cell; ++i)
    //         {
    //             auto indx = cell->vertex_dof_index(i, 0);

    //             grad[0][indx] += function_on_cell.dx(cell->vertex(i));
    //             grad[1][indx] += function_on_cell.dy(cell->vertex(i));
    //             ++(N[indx]); 

    //             if (!flg[indx])
    //             {
    //             T[indx] = 
    //                 // temperature(cell->vertex_dof_index (i, 0));
    //                 // function_on_cell(cell->vertex(i));
    //                 // function_on_cell(prmt::Point<2>(cell->vertex(i)(0), cell->vertex(i)(1)));
    //                 // cell->vertex(i)(0) * cell->vertex(i)(0);
    //                 function_on_cell(prmt::Point<2>(cell->vertex(i)(0), cell->vertex(i)(1)))-
    //             temperature(cell->vertex_dof_index (i, 0));
    //             // if (std::abs(T[indx]) > 1e-5)
    //             //     fprintf(F, "%f %f %f %f %f %f %f %f %f %f %f %f %f %d\n",
    //             //             cell->vertex(0)(0), cell->vertex(0)(1),
    //             //             cell->vertex(1)(0), cell->vertex(1)(1),
    //             //             cell->vertex(3)(0), cell->vertex(3)(1),
    //             //             cell->vertex(2)(0), cell->vertex(2)(1),
    //             //             temperature(cell->vertex_dof_index (0, 0)),
    //             //             temperature(cell->vertex_dof_index (1, 0)),
    //             //             temperature(cell->vertex_dof_index (3, 0)),
    //             //             temperature(cell->vertex_dof_index (2, 0)),
    //             //             T[indx],
    //             //             i
    //             //             );
    //             flg[indx] = true;
    //             };
    //         };
    //     };
    //     fclose(F);
    //     for (st i = 0; i < N.size(); ++i)
    //     {
    //         // printf("%d\n", N[i]);
    //         grad[0][i] /= N[i];
    //         grad[1][i] /= N[i];
    //         // T[i] /= N[i];
    //     };

    //     {
    //         dealii::DataOut<2> data_out;
    //         data_out.attach_dof_handler (domain.dof_handler);
    //         data_out.add_data_vector (grad[0], "grad_dx");
    //         data_out.add_data_vector (grad[1], "drad_dy");
    //         data_out.add_data_vector (T, "T");
    //         data_out.build_patches ();

    //         auto name = file_name;
    //         name += ".gpd";

    //         std::ofstream output (name);
    //         data_out.write_gnuplot (output);
    //     };
    //     FILE *F1;
    //     F1 = fopen("gradx.gpd", "w");
    //     for (st i = 0; i < ps.size(); ++i)
    //     {
    //         fprintf(F1, "%lf %lf %lf\n", ps[i].x(), ps[i].y(), gradx[i]);
    //     };
    //     fclose(F1);
    // }

void print_elastic_deformation (const dealii::Vector<dbl> &move, 
        const dealii::DoFHandler<2> &dof_handler,
        const str file_name)
{
    cu32 dofs_per_cell = dof_handler.get_fe().dofs_per_cell;

    arr<dealii::Vector<dbl>, 2> grad;

    grad[0] .reinit (dof_handler.n_dofs());
    grad[1] .reinit (dof_handler.n_dofs());

    vec<st> N(dof_handler.n_dofs());

    for (
            auto cell = dof_handler.begin_active(); 
            cell     != dof_handler.end(); 
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
        for (st component = 0; component < 2; ++component)
        {
            arr<dbl, 4> values = {
                move(cell->vertex_dof_index (0, component)),
                move(cell->vertex_dof_index (1, component)),
                move(cell->vertex_dof_index (3, component)),
                move(cell->vertex_dof_index (2, component))};

            Scalar4PointsFunc<2> function_on_cell(points, values);

            for (st i = 0; i < 4; ++i)
            {
                auto indx = cell->vertex_dof_index(i, component);

                if ((cell->material_id() == 0) or (cell->material_id() == 2))
                    // if (cell->material_id() == 0)
                {
                    grad[0][indx] += function_on_cell.dx(cell->vertex(i));
                    grad[1][indx] += function_on_cell.dy(cell->vertex(i));
                    // grad[0][indx] = 1.0;
                    // grad[1][indx] = 1.0;
                }
                else
                {
                    grad[0][indx] = 0.0;
                    grad[1][indx] = 0.0;
                };
                ++(N[indx]); 
            };
        };
    };
    for (st i = 0; i < N.size(); ++i)
    {
        grad[0][i] /= N[i];
        grad[1][i] /= N[i];
    };

    {
        dealii::DataOut<2> data_out;
        data_out.attach_dof_handler (dof_handler);
        data_out.add_data_vector (grad[0], "grad_dx");
        data_out.add_data_vector (grad[1], "drad_dy");
        // data_out.add_data_vector (T, "T");
        data_out.build_patches ();

        auto name = file_name;
        // name += ".gpd";

        std::ofstream output (name);
        data_out.write_gnuplot (output);
    };

    {
        std::ofstream out ("hole/deform_hole_x.bin", std::ios::out | std::ios::binary);
        for (st i = 0; i < grad[0].size(); ++i)
        {
            out.write ((char *) &(grad[0][i]), sizeof grad[0][0]);
        };
        out.close ();
    };

    {
        std::ofstream out ("hole/deform_hole_y.bin", std::ios::out | std::ios::binary);
        for (st i = 0; i < grad[0].size(); ++i)
        {
            out.write ((char *) &(grad[1][i]), sizeof grad[0][0]);
        };
        out.close ();
    };
};

void print_elastic_deformation (const dealii::Vector<dbl> &move, 
        const dealii::DoFHandler<2> &dof_handler,
        const str file_name, const str file_name_bin)
{
    cu32 dofs_per_cell = dof_handler.get_fe().dofs_per_cell;

    arr<dealii::Vector<dbl>, 2> grad;

    grad[0] .reinit (dof_handler.n_dofs());
    grad[1] .reinit (dof_handler.n_dofs());

    vec<st> N(dof_handler.n_dofs());

    for (
            auto cell = dof_handler.begin_active(); 
            cell     != dof_handler.end(); 
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
        for (st component = 0; component < 2; ++component)
        {
            arr<dbl, 4> values = {
                move(cell->vertex_dof_index (0, component)),
                move(cell->vertex_dof_index (1, component)),
                move(cell->vertex_dof_index (3, component)),
                move(cell->vertex_dof_index (2, component))};

            Scalar4PointsFunc<2> function_on_cell(points, values);

            for (st i = 0; i < 4; ++i)
            {
                auto indx = cell->vertex_dof_index(i, component);

                // if ((cell->material_id() == 0) or (cell->material_id() == 2))
                if (cell->material_id() == 0)
                {
                    grad[0][indx] += function_on_cell.dx(cell->vertex(i));
                    grad[1][indx] += function_on_cell.dy(cell->vertex(i));
                    // grad[0][indx] = 1.0;
                    // grad[1][indx] = 1.0;
                }
                else
                {
                    grad[0][indx] = 0.0;
                    grad[1][indx] = 0.0;
                };
                ++(N[indx]); 
            };
        };
    };
    for (st i = 0; i < N.size(); ++i)
    {
        grad[0][i] /= N[i];
        grad[1][i] /= N[i];
    };

    {
        dealii::DataOut<2> data_out;
        data_out.attach_dof_handler (dof_handler);
        data_out.add_data_vector (grad[0], "grad_dx");
        data_out.add_data_vector (grad[1], "drad_dy");
        // data_out.add_data_vector (T, "T");
        data_out.build_patches ();

        auto name = file_name;
        // name += ".gpd";

        std::ofstream output (name);
        data_out.write_gnuplot (output);
    };

    {
        std::ofstream out ("hole/"+file_name_bin+"/deform_hole_x.bin", std::ios::out | std::ios::binary);
        for (st i = 0; i < grad[0].size(); ++i)
        {
            out.write ((char *) &(grad[0][i]), sizeof grad[0][0]);
        };
        out.close ();
    };

    {
        std::ofstream out ("hole/"+file_name_bin+"/deform_hole_y.bin", std::ios::out | std::ios::binary);
        for (st i = 0; i < grad[0].size(); ++i)
        {
            out.write ((char *) &(grad[1][i]), sizeof grad[0][0]);
        };
        out.close ();
    };
};

void print_elastic_stress (const dealii::Vector<dbl> &move, 
        const dealii::DoFHandler<2> &dof_handler,
        const ATools::FourthOrderTensor &C,
        const str file_name)
{
    cu32 dofs_per_cell = dof_handler.get_fe().dofs_per_cell;

    arr<dealii::Vector<dbl>, 2> grad;

    grad[0] .reinit (dof_handler.n_dofs());
    grad[1] .reinit (dof_handler.n_dofs());

    vec<st> N(dof_handler.n_dofs());

    for (
            auto cell = dof_handler.begin_active(); 
            cell     != dof_handler.end(); 
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
        for (st component = 0; component < 2; ++component)
        {
            arr<dbl, 4> values = {
                move(cell->vertex_dof_index (0, component)),
                move(cell->vertex_dof_index (1, component)),
                move(cell->vertex_dof_index (3, component)),
                move(cell->vertex_dof_index (2, component))};

            Scalar4PointsFunc<2> function_on_cell(points, values);

            for (st i = 0; i < 4; ++i)
            {
                auto indx = cell->vertex_dof_index(i, component);

                // if ((cell->material_id() == 0) or (cell->material_id() == 2))
                //     // if (cell->material_id() == 0)
                // {
                grad[0][indx] += function_on_cell.dx(cell->vertex(i));
                grad[1][indx] += function_on_cell.dy(cell->vertex(i));
                // grad[0][indx] = 1.0;
                // grad[1][indx] = 1.0;
                // }
                // else
                // {
                //     grad[0][indx] = 0.0;
                //     grad[1][indx] = 0.0;
                // };
                ++(N[indx]); 
            };
        };
    };
    for (st i = 0; i < N.size(); ++i)
    {
        grad[0][i] /= N[i];
        grad[1][i] /= N[i];
    };

    arr<dealii::Vector<dbl>, 2> stress;

    stress[0] .reinit (dof_handler.n_dofs());
    stress[1] .reinit (dof_handler.n_dofs());

    for (
            auto cell = dof_handler.begin_active(); 
            cell     != dof_handler.end(); 
            ++cell
        )
    {
        for (st i = 0; i < 4; ++i)
        {
            auto indx_1 = cell->vertex_dof_index(i, 0);
            auto indx_2 = cell->vertex_dof_index(i, 1);
            stress[0][indx_1] = 
                C[0][0][0][0] * grad[0][indx_1] +
                C[0][0][0][1] * grad[0][indx_2] +
                C[0][0][1][0] * grad[1][indx_1] +
                C[0][0][1][1] * grad[1][indx_2];
            stress[0][indx_2] = 
                C[0][1][0][0] * grad[0][indx_1] +
                C[0][1][0][1] * grad[0][indx_2] +
                C[0][1][1][0] * grad[1][indx_1] +
                C[0][1][1][1] * grad[1][indx_2];
            stress[1][indx_1] = 
                C[1][0][0][0] * grad[0][indx_1] +
                C[1][0][0][1] * grad[0][indx_2] +
                C[1][0][1][0] * grad[1][indx_1] +
                C[1][0][1][1] * grad[1][indx_2];
            stress[1][indx_2] = 
                C[1][1][0][0] * grad[0][indx_1] +
                C[1][1][0][1] * grad[0][indx_2] +
                C[1][1][1][0] * grad[1][indx_1] +
                C[1][1][1][1] * grad[1][indx_2];
        };
    };

    {
        dealii::DataOut<2> data_out;
        data_out.attach_dof_handler (dof_handler);
        data_out.add_data_vector (stress[0], "grad_dx");
        data_out.add_data_vector (stress[1], "drad_dy");
        // data_out.add_data_vector (T, "T");
        data_out.build_patches ();

        auto name = file_name;
        // name += ".gpd";

        std::ofstream output (name);
        data_out.write_gnuplot (output);
    };

    {
        std::ofstream out ("hole/stress_hole_x.bin", std::ios::out | std::ios::binary);
        for (st i = 0; i < stress[0].size(); ++i)
        {
            out.write ((char *) &(stress[0][i]), sizeof stress[0][0]);
        };
        out.close ();
    };

    {
        std::ofstream out ("hole/stress_hole_y.bin", std::ios::out | std::ios::binary);
        for (st i = 0; i < stress[0].size(); ++i)
        {
            out.write ((char *) &(stress[1][i]), sizeof stress[0][0]);
        };
        out.close ();
    };

    dbl Int = 0.0;
    {
        cu32 dofs_per_cell = dof_handler.get_fe().dofs_per_cell;

        dealii::QGauss<2>  quadrature_formula(2);

        dealii::FEValues<2> fe_values (dof_handler.get_fe(), quadrature_formula,
                dealii::update_values | dealii::update_gradients | dealii::update_quadrature_points |
                dealii::update_JxW_values);

        cst num_quad_points = quadrature_formula.size();
        for (
                auto cell = dof_handler.begin_active(); 
                cell     != dof_handler.end(); 
                ++cell
            )
        {
            fe_values .reinit (cell);
            for (st i = 0; i < 4; ++i)
            {
                auto indx_1 = cell->vertex_dof_index(i, 0);
                for (st q_point = 0; q_point < num_quad_points; ++q_point)
                {
                    Int += 
                        stress[1](indx_1) *
                        // move(indx_1) *
                        // 1.0 *
                        // fe_values.shape_grad (i, q_point)[0] *
                        fe_values.shape_value (i, q_point) *
                        fe_values.JxW(q_point);
                };
            };
        };

    };

    // for (st i = 0; i < stress[0].size(); ++i)
    // {
    //     Int += stress[0][i];
    // };
    printf("Integral %f\n", Int);
}

void print_elastic_stress (const dealii::Vector<dbl> &move, 
        const dealii::DoFHandler<2> &dof_handler,
        const ATools::FourthOrderTensor &C,
        const str file_name, const str file_name_bin)
{
    cu32 dofs_per_cell = dof_handler.get_fe().dofs_per_cell;

    arr<dealii::Vector<dbl>, 2> grad;

    grad[0] .reinit (dof_handler.n_dofs());
    grad[1] .reinit (dof_handler.n_dofs());

    vec<st> N(dof_handler.n_dofs());

    for (
            auto cell = dof_handler.begin_active(); 
            cell     != dof_handler.end(); 
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
        for (st component = 0; component < 2; ++component)
        {
            arr<dbl, 4> values = {
                move(cell->vertex_dof_index (0, component)),
                move(cell->vertex_dof_index (1, component)),
                move(cell->vertex_dof_index (3, component)),
                move(cell->vertex_dof_index (2, component))};

            Scalar4PointsFunc<2> function_on_cell(points, values);

            for (st i = 0; i < 4; ++i)
            {
                auto indx = cell->vertex_dof_index(i, component);

                // if ((cell->material_id() == 0) or (cell->material_id() == 2))
                //     // if (cell->material_id() == 0)
                // {
                grad[0][indx] += function_on_cell.dx(cell->vertex(i));
                grad[1][indx] += function_on_cell.dy(cell->vertex(i));
                // grad[0][indx] = 1.0;
                // grad[1][indx] = 1.0;
                // }
                // else
                // {
                //     grad[0][indx] = 0.0;
                //     grad[1][indx] = 0.0;
                // };
                ++(N[indx]); 
            };
        };
    };
    for (st i = 0; i < N.size(); ++i)
    {
        grad[0][i] /= N[i];
        grad[1][i] /= N[i];
    };

    arr<dealii::Vector<dbl>, 2> stress;

    stress[0] .reinit (dof_handler.n_dofs());
    stress[1] .reinit (dof_handler.n_dofs());

    for (
            auto cell = dof_handler.begin_active(); 
            cell     != dof_handler.end(); 
            ++cell
        )
    {
        for (st i = 0; i < 4; ++i)
        {
            auto indx_1 = cell->vertex_dof_index(i, 0);
            auto indx_2 = cell->vertex_dof_index(i, 1);
            stress[0][indx_1] = 
                C[0][0][0][0] * grad[0][indx_1] +
                C[0][0][0][1] * grad[0][indx_2] +
                C[0][0][1][0] * grad[1][indx_1] +
                C[0][0][1][1] * grad[1][indx_2];
            stress[0][indx_2] = 
                C[0][1][0][0] * grad[0][indx_1] +
                C[0][1][0][1] * grad[0][indx_2] +
                C[0][1][1][0] * grad[1][indx_1] +
                C[0][1][1][1] * grad[1][indx_2];
            stress[1][indx_1] = 
                C[1][0][0][0] * grad[0][indx_1] +
                C[1][0][0][1] * grad[0][indx_2] +
                C[1][0][1][0] * grad[1][indx_1] +
                C[1][0][1][1] * grad[1][indx_2];
            stress[1][indx_2] = 
                C[1][1][0][0] * grad[0][indx_1] +
                C[1][1][0][1] * grad[0][indx_2] +
                C[1][1][1][0] * grad[1][indx_1] +
                C[1][1][1][1] * grad[1][indx_2];
        };
    };

    {
        dealii::DataOut<2> data_out;
        data_out.attach_dof_handler (dof_handler);
        data_out.add_data_vector (stress[0], "grad_dx");
        data_out.add_data_vector (stress[1], "drad_dy");
        // data_out.add_data_vector (T, "T");
        data_out.build_patches ();

        auto name = file_name;
        // name += ".gpd";

        std::ofstream output (name);
        data_out.write_gnuplot (output);
    };

    {
        std::ofstream out ("hole/"+file_name_bin+"/stress_hole_x.bin", std::ios::out | std::ios::binary);
        for (st i = 0; i < stress[0].size(); ++i)
        {
            out.write ((char *) &(stress[0][i]), sizeof stress[0][0]);
        };
        out.close ();
    };

    {
        std::ofstream out ("hole/"+file_name_bin+"/stress_hole_y.bin", std::ios::out | std::ios::binary);
        for (st i = 0; i < stress[0].size(); ++i)
        {
            out.write ((char *) &(stress[1][i]), sizeof stress[0][0]);
        };
        out.close ();
    };

    dbl Int = 0.0;
    {
        cu32 dofs_per_cell = dof_handler.get_fe().dofs_per_cell;

        dealii::QGauss<2>  quadrature_formula(2);

        dealii::FEValues<2> fe_values (dof_handler.get_fe(), quadrature_formula,
                dealii::update_values | dealii::update_gradients | dealii::update_quadrature_points |
                dealii::update_JxW_values);

        cst num_quad_points = quadrature_formula.size();
        for (
                auto cell = dof_handler.begin_active(); 
                cell     != dof_handler.end(); 
                ++cell
            )
        {
            fe_values .reinit (cell);
            for (st i = 0; i < 4; ++i)
            {
                auto indx_1 = cell->vertex_dof_index(i, 0);
                for (st q_point = 0; q_point < num_quad_points; ++q_point)
                {
                    Int += 
                        stress[1](indx_1) *
                        // move(indx_1) *
                        // 1.0 *
                        // fe_values.shape_grad (i, q_point)[0] *
                        fe_values.shape_value (i, q_point) *
                        fe_values.JxW(q_point);
                };
            };
        };

    };

    // for (st i = 0; i < stress[0].size(); ++i)
    // {
    //     Int += stress[0][i];
    // };
    printf("Integral %f\n", Int);
}

void print_elastic_deformation_2 (const dealii::Vector<dbl> &move, 
        const dealii::DoFHandler<2> &dof_handler,
        const str file_name)
{
    cu32 dofs_per_cell = dof_handler.get_fe().dofs_per_cell;

    arr<dealii::Vector<dbl>, 2> grad;

    grad[0] .reinit (dof_handler.n_dofs());
    grad[1] .reinit (dof_handler.n_dofs());

    vec<st> N(dof_handler.n_dofs());

    for (
            auto cell = dof_handler.begin_active(); 
            cell     != dof_handler.end(); 
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
        for (st component = 0; component < 2; ++component)
        {
            arr<dbl, 4> values = {
                move(cell->vertex_dof_index (0, component)),
                move(cell->vertex_dof_index (1, component)),
                move(cell->vertex_dof_index (3, component)),
                move(cell->vertex_dof_index (2, component))};

            Scalar4PointsFunc<2> function_on_cell(points, values);

            for (st i = 0; i < 4; ++i)
            {
                auto indx = cell->vertex_dof_index(i, component);

                // if (cell->material_id() == 0)
                if ((cell->material_id() == 0) or (cell->material_id() == 2))
                {
                    grad[0][indx] += function_on_cell.dx(cell->vertex(i));
                    grad[1][indx] += function_on_cell.dy(cell->vertex(i));
                }
                else
                {
                    grad[0][indx] = 0.0;
                    grad[1][indx] = 0.0;
                };
                ++(N[indx]); 
            };
        };
    };
    for (st i = 0; i < N.size(); ++i)
    {
        grad[0][i] /= N[i];
        grad[1][i] /= N[i];
    };

    {
        arr<arr<dealii::Vector<dbl>, 2>, 2> grad_2;

        grad_2[0][0] .reinit (dof_handler.n_dofs());
        grad_2[0][1] .reinit (dof_handler.n_dofs());
        grad_2[1][0] .reinit (dof_handler.n_dofs());
        grad_2[1][1] .reinit (dof_handler.n_dofs());

        for (
                auto cell = dof_handler.begin_active(); 
                cell     != dof_handler.end(); 
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
            for (st g_comp = 0; g_comp < 2; ++g_comp)
            {
                for (st component = 0; component < 2; ++component)
                {
                    arr<dbl, 4> values = {
                        grad[g_comp](cell->vertex_dof_index (0, component)),
                        grad[g_comp](cell->vertex_dof_index (1, component)),
                        grad[g_comp](cell->vertex_dof_index (3, component)),
                        grad[g_comp](cell->vertex_dof_index (2, component))};

                    Scalar4PointsFunc<2> function_on_cell(points, values);

                    for (st i = 0; i < 4; ++i)
                    {
                        auto indx = cell->vertex_dof_index(i, component);

                        if ((cell->material_id() == 0) or (std::abs(grad[0][indx]) < 1e-5))
                            // if (cell->material_id() == 0)
                        {
                            grad_2[g_comp][0][indx] += function_on_cell.dx(cell->vertex(i));
                            grad_2[g_comp][1][indx] += function_on_cell.dy(cell->vertex(i));
                        }
                        else
                        {
                            grad_2[g_comp][0][indx] = 0.0;
                            grad_2[g_comp][1][indx] = 0.0;
                        };
                        ++(N[indx]); 
                    };
                };
            };
        };
        for (st i = 0; i < N.size(); ++i)
        {
            grad_2[0][0][i] /= N[i];
            grad_2[0][1][i] /= N[i];
            grad_2[1][0][i] /= N[i];
            grad_2[1][1][i] /= N[i];
        };

        {
            dealii::DataOut<2> data_out;
            data_out.attach_dof_handler (dof_handler);
            data_out.add_data_vector (grad_2[0][0], "grad_dx_dx");
            data_out.add_data_vector (grad_2[0][1], "grad_dx_dy");
            data_out.add_data_vector (grad_2[1][0], "grad_dy_dx");
            data_out.add_data_vector (grad_2[1][1], "grad_dy_dy");
            // data_out.add_data_vector (T, "T");
            data_out.build_patches ();

            auto name = file_name;
            // name += ".gpd";

            std::ofstream output (name);
            data_out.write_gnuplot (output);
        };

        {
            std::ofstream out ("hole/deform_hole_xx.bin", std::ios::out | std::ios::binary);
            for (st i = 0; i < grad_2[0][0].size(); ++i)
            {
                out.write ((char *) &(grad_2[0][0][i]), sizeof grad_2[0][0][0]);
            };
            out.close ();
        };

        {
            std::ofstream out ("hole/deform_hole_xy.bin", std::ios::out | std::ios::binary);
            for (st i = 0; i < grad_2[0][0].size(); ++i)
            {
                out.write ((char *) &(grad_2[0][1][i]), sizeof grad_2[0][1][0]);
            };
            out.close ();
        };

        {
            std::ofstream out ("hole/deform_hole_yx.bin", std::ios::out | std::ios::binary);
            for (st i = 0; i < grad_2[0][0].size(); ++i)
            {
                out.write ((char *) &(grad_2[1][0][i]), sizeof grad_2[1][0][0]);
            };
            out.close ();
        };

        {
            std::ofstream out ("hole/deform_hole_yy.bin", std::ios::out | std::ios::binary);
            for (st i = 0; i < grad_2[0][0].size(); ++i)
            {
                out.write ((char *) &(grad_2[1][1][i]), sizeof grad_2[1][1][0]);
            };
            out.close ();
        };
    };
    {
        cu32 dofs_per_cell = dof_handler.get_fe().dofs_per_cell;

        dealii::QGauss<2>  quadrature_formula(2);

        dealii::FEValues<2> fe_values (dof_handler.get_fe(), quadrature_formula,
                dealii::update_gradients | dealii::update_quadrature_points | dealii::update_JxW_values);

        cst num_quad_points = quadrature_formula.size();

        arr<dealii::Vector<dbl>, 2> grad_2;

        grad_2[0] .reinit (dof_handler.get_tria().n_active_cells());
        grad_2[1] .reinit (dof_handler.get_tria().n_active_cells());

        vec<prmt::Point<2>> midle_point (dof_handler.get_tria().n_active_cells());
        vec<dealii::types::global_dof_index> local_dof_indices (dofs_per_cell);

        st number_of_cell = 0;
        for (
                auto cell = dof_handler.begin_active(); 
                cell     != dof_handler.end(); 
                ++cell
            )
        {
            fe_values .reinit (cell);

            if (cell->material_id() == 0)
            {
                dbl area = 0.0;
                for (st q_point = 0; q_point < num_quad_points; ++q_point)
                {
                    area +=
                        fe_values.JxW(q_point);
                };

                cell->get_dof_indices (local_dof_indices);

                for (st component = 0; component < 2; ++component)
                {
                    arr<dbl, 8> tmp;
                    for (st i = 0; i < dofs_per_cell; ++i)
                    {
                        tmp[i] = 0.0;
                        for (st q_point = 0; q_point < num_quad_points; ++q_point)
                        {
                            tmp[i] += grad[0](local_dof_indices[i]) *
                                fe_values.shape_grad (i, q_point)[component] *
                                fe_values.JxW(q_point);
                        };
                    };
                    grad_2[component][number_of_cell] = (tmp[0] + tmp[2] + tmp[4] + tmp[6]) / (area);
                };
            }
            else
            {
                grad[0][number_of_cell] = 0.0;
                grad[1][number_of_cell] = 0.0;
            };

            midle_point[number_of_cell].x() = 
                (cell->vertex(0)(0) +
                 cell->vertex(1)(0) +
                 cell->vertex(2)(0) +
                 cell->vertex(3)(0)) / 4.0;
            midle_point[number_of_cell].y() = 
                (cell->vertex(0)(1) +
                 cell->vertex(1)(1) +
                 cell->vertex(2)(1) +
                 cell->vertex(3)(1)) / 4.0;

            ++number_of_cell;
        };
        FILE *F;
        F = fopen("deform_mean_2.gpd", "w");
        for (st i = 0; i < grad_2[0].size(); ++i)
        {
            fprintf(F, "%f %f %f %f\n", 
                    midle_point[i].x(), midle_point[i].y(), grad_2[0][i], grad_2[1][i]);
        };
        fclose(F);

    };
}

void print_elastic_deformation_2 (const dealii::Vector<dbl> &move, 
        const dealii::DoFHandler<2> &dof_handler,
        const str file_name, const str file_name_bin)
{
    cu32 dofs_per_cell = dof_handler.get_fe().dofs_per_cell;

    arr<dealii::Vector<dbl>, 2> grad;

    grad[0] .reinit (dof_handler.n_dofs());
    grad[1] .reinit (dof_handler.n_dofs());

    vec<st> N(dof_handler.n_dofs());

    for (
            auto cell = dof_handler.begin_active(); 
            cell     != dof_handler.end(); 
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
        for (st component = 0; component < 2; ++component)
        {
            arr<dbl, 4> values = {
                move(cell->vertex_dof_index (0, component)),
                move(cell->vertex_dof_index (1, component)),
                move(cell->vertex_dof_index (3, component)),
                move(cell->vertex_dof_index (2, component))};

            Scalar4PointsFunc<2> function_on_cell(points, values);

            for (st i = 0; i < 4; ++i)
            {
                auto indx = cell->vertex_dof_index(i, component);

                // if (cell->material_id() == 0)
                if ((cell->material_id() == 0) or (cell->material_id() == 2))
                {
                    grad[0][indx] += function_on_cell.dx(cell->vertex(i));
                    grad[1][indx] += function_on_cell.dy(cell->vertex(i));
                }
                else
                {
                    grad[0][indx] = 0.0;
                    grad[1][indx] = 0.0;
                };
                ++(N[indx]); 
            };
        };
    };
    for (st i = 0; i < N.size(); ++i)
    {
        grad[0][i] /= N[i];
        grad[1][i] /= N[i];
    };

    {
        arr<arr<dealii::Vector<dbl>, 2>, 2> grad_2;

        grad_2[0][0] .reinit (dof_handler.n_dofs());
        grad_2[0][1] .reinit (dof_handler.n_dofs());
        grad_2[1][0] .reinit (dof_handler.n_dofs());
        grad_2[1][1] .reinit (dof_handler.n_dofs());

        for (
                auto cell = dof_handler.begin_active(); 
                cell     != dof_handler.end(); 
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
            for (st g_comp = 0; g_comp < 2; ++g_comp)
            {
                for (st component = 0; component < 2; ++component)
                {
                    arr<dbl, 4> values = {
                        grad[g_comp](cell->vertex_dof_index (0, component)),
                        grad[g_comp](cell->vertex_dof_index (1, component)),
                        grad[g_comp](cell->vertex_dof_index (3, component)),
                        grad[g_comp](cell->vertex_dof_index (2, component))};

                    Scalar4PointsFunc<2> function_on_cell(points, values);

                    for (st i = 0; i < 4; ++i)
                    {
                        auto indx = cell->vertex_dof_index(i, component);

                        if ((cell->material_id() == 0) or (std::abs(grad[0][indx]) < 1e-5))
                            // if (cell->material_id() == 0)
                        {
                            grad_2[g_comp][0][indx] += function_on_cell.dx(cell->vertex(i));
                            grad_2[g_comp][1][indx] += function_on_cell.dy(cell->vertex(i));
                        }
                        else
                        {
                            grad_2[g_comp][0][indx] = 0.0;
                            grad_2[g_comp][1][indx] = 0.0;
                        };
                        ++(N[indx]); 
                    };
                };
            };
        };
        for (st i = 0; i < N.size(); ++i)
        {
            grad_2[0][0][i] /= N[i];
            grad_2[0][1][i] /= N[i];
            grad_2[1][0][i] /= N[i];
            grad_2[1][1][i] /= N[i];
        };

        {
            dealii::DataOut<2> data_out;
            data_out.attach_dof_handler (dof_handler);
            data_out.add_data_vector (grad_2[0][0], "grad_dx_dx");
            data_out.add_data_vector (grad_2[0][1], "grad_dx_dy");
            data_out.add_data_vector (grad_2[1][0], "grad_dy_dx");
            data_out.add_data_vector (grad_2[1][1], "grad_dy_dy");
            // data_out.add_data_vector (T, "T");
            data_out.build_patches ();

            auto name = file_name;
            // name += ".gpd";

            std::ofstream output (name);
            data_out.write_gnuplot (output);
        };

        {
            std::ofstream out ("hole/"+file_name_bin+"/deform_hole_xx.bin", std::ios::out | std::ios::binary);
            for (st i = 0; i < grad_2[0][0].size(); ++i)
            {
                out.write ((char *) &(grad_2[0][0][i]), sizeof grad_2[0][0][0]);
            };
            out.close ();
        };

        {
            std::ofstream out ("hole/"+file_name_bin+"/deform_hole_xy.bin", std::ios::out | std::ios::binary);
            for (st i = 0; i < grad_2[0][0].size(); ++i)
            {
                out.write ((char *) &(grad_2[0][1][i]), sizeof grad_2[0][1][0]);
            };
            out.close ();
        };

        {
            std::ofstream out ("hole/"+file_name_bin+"/deform_hole_yx.bin", std::ios::out | std::ios::binary);
            for (st i = 0; i < grad_2[0][0].size(); ++i)
            {
                out.write ((char *) &(grad_2[1][0][i]), sizeof grad_2[1][0][0]);
            };
            out.close ();
        };

        {
            std::ofstream out ("hole/"+file_name_bin+"/deform_hole_yy.bin", std::ios::out | std::ios::binary);
            for (st i = 0; i < grad_2[0][0].size(); ++i)
            {
                out.write ((char *) &(grad_2[1][1][i]), sizeof grad_2[1][1][0]);
            };
            out.close ();
        };
    };
    {
        cu32 dofs_per_cell = dof_handler.get_fe().dofs_per_cell;

        dealii::QGauss<2>  quadrature_formula(2);

        dealii::FEValues<2> fe_values (dof_handler.get_fe(), quadrature_formula,
                dealii::update_gradients | dealii::update_quadrature_points | dealii::update_JxW_values);

        cst num_quad_points = quadrature_formula.size();

        arr<dealii::Vector<dbl>, 2> grad_2;

        grad_2[0] .reinit (dof_handler.get_tria().n_active_cells());
        grad_2[1] .reinit (dof_handler.get_tria().n_active_cells());

        vec<prmt::Point<2>> midle_point (dof_handler.get_tria().n_active_cells());
        vec<dealii::types::global_dof_index> local_dof_indices (dofs_per_cell);

        st number_of_cell = 0;
        for (
                auto cell = dof_handler.begin_active(); 
                cell     != dof_handler.end(); 
                ++cell
            )
        {
            fe_values .reinit (cell);

            if (cell->material_id() == 0)
            {
                dbl area = 0.0;
                for (st q_point = 0; q_point < num_quad_points; ++q_point)
                {
                    area +=
                        fe_values.JxW(q_point);
                };

                cell->get_dof_indices (local_dof_indices);

                for (st component = 0; component < 2; ++component)
                {
                    arr<dbl, 8> tmp;
                    for (st i = 0; i < dofs_per_cell; ++i)
                    {
                        tmp[i] = 0.0;
                        for (st q_point = 0; q_point < num_quad_points; ++q_point)
                        {
                            tmp[i] += grad[0](local_dof_indices[i]) *
                                fe_values.shape_grad (i, q_point)[component] *
                                fe_values.JxW(q_point);
                        };
                    };
                    grad_2[component][number_of_cell] = (tmp[0] + tmp[2] + tmp[4] + tmp[6]) / (area);
                };
            }
            else
            {
                grad[0][number_of_cell] = 0.0;
                grad[1][number_of_cell] = 0.0;
            };

            midle_point[number_of_cell].x() = 
                (cell->vertex(0)(0) +
                 cell->vertex(1)(0) +
                 cell->vertex(2)(0) +
                 cell->vertex(3)(0)) / 4.0;
            midle_point[number_of_cell].y() = 
                (cell->vertex(0)(1) +
                 cell->vertex(1)(1) +
                 cell->vertex(2)(1) +
                 cell->vertex(3)(1)) / 4.0;

            ++number_of_cell;
        };
        FILE *F;
        F = fopen("deform_mean_2.gpd", "w");
        for (st i = 0; i < grad_2[0].size(); ++i)
        {
            fprintf(F, "%f %f %f %f\n", 
                    midle_point[i].x(), midle_point[i].y(), grad_2[0][i], grad_2[1][i]);
        };
        fclose(F);

    };
}

void print_elastic_mean_micro_stress (const dealii::Vector<dbl> &move, 
        const dealii::DoFHandler<2> &dof_handler,
        const ATools::FourthOrderTensor &C,
        const str file_name,
        cdbl epsilon)
{
    cu32 dofs_per_cell = dof_handler.get_fe().dofs_per_cell;

    arr<dealii::Vector<dbl>, 2> grad;

    grad[0] .reinit (dof_handler.n_dofs());
    grad[1] .reinit (dof_handler.n_dofs());

    vec<st> N(dof_handler.n_dofs());

    for (
            auto cell = dof_handler.begin_active(); 
            cell     != dof_handler.end(); 
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
        for (st component = 0; component < 2; ++component)
        {
            arr<dbl, 4> values = {
                move(cell->vertex_dof_index (0, component)),
                move(cell->vertex_dof_index (1, component)),
                move(cell->vertex_dof_index (3, component)),
                move(cell->vertex_dof_index (2, component))};

            Scalar4PointsFunc<2> function_on_cell(points, values);

            for (st i = 0; i < 4; ++i)
            {
                auto indx = cell->vertex_dof_index(i, component);

                // if (cell->material_id() == 0)
                if ((cell->material_id() == 0) or (cell->material_id() == 2))
                {
                    grad[0][indx] += function_on_cell.dx(cell->vertex(i));
                    grad[1][indx] += function_on_cell.dy(cell->vertex(i));
                }
                else
                {
                    grad[0][indx] = 0.0;
                    grad[1][indx] = 0.0;
                };
                ++(N[indx]); 
            };
        };
    };
    for (st i = 0; i < N.size(); ++i)
    {
        grad[0][i] /= N[i];
        grad[1][i] /= N[i];
    };

    arr<dealii::Vector<dbl>, 2> stress;

    stress[0] .reinit (dof_handler.n_dofs());
    stress[1] .reinit (dof_handler.n_dofs());

    for (
            auto cell = dof_handler.begin_active(); 
            cell     != dof_handler.end(); 
            ++cell
        )
    {
        for (st i = 0; i < 4; ++i)
        {
            auto indx_1 = cell->vertex_dof_index(i, 0);
            auto indx_2 = cell->vertex_dof_index(i, 1);
            stress[0][indx_1] = 
                C[0][0][0][0] * grad[0][indx_1] +
                C[0][0][0][1] * grad[0][indx_2] +
                C[0][0][1][0] * grad[1][indx_1] +
                C[0][0][1][1] * grad[1][indx_2];
            stress[0][indx_2] = 
                C[0][1][0][0] * grad[0][indx_1] +
                C[0][1][0][1] * grad[0][indx_2] +
                C[0][1][1][0] * grad[1][indx_1] +
                C[0][1][1][1] * grad[1][indx_2];
            stress[1][indx_1] = 
                C[1][0][0][0] * grad[0][indx_1] +
                C[1][0][0][1] * grad[0][indx_2] +
                C[1][0][1][0] * grad[1][indx_1] +
                C[1][0][1][1] * grad[1][indx_2];
            stress[1][indx_2] = 
                C[1][1][0][0] * grad[0][indx_1] +
                C[1][1][0][1] * grad[0][indx_2] +
                C[1][1][1][0] * grad[1][indx_1] +
                C[1][1][1][1] * grad[1][indx_2];
        };
    };

    {
        arr<arr<dealii::Vector<dbl>, 2>, 2> stress_2;

        stress_2[0][0] .reinit (dof_handler.n_dofs());
        stress_2[0][1] .reinit (dof_handler.n_dofs());
        stress_2[1][0] .reinit (dof_handler.n_dofs());
        stress_2[1][1] .reinit (dof_handler.n_dofs());

        for (
                auto cell = dof_handler.begin_active(); 
                cell     != dof_handler.end(); 
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
            for (st g_comp = 0; g_comp < 2; ++g_comp)
            {
                for (st component = 0; component < 2; ++component)
                {
                    arr<dbl, 4> values = {
                        stress[g_comp](cell->vertex_dof_index (0, component)),
                        stress[g_comp](cell->vertex_dof_index (1, component)),
                        stress[g_comp](cell->vertex_dof_index (3, component)),
                        stress[g_comp](cell->vertex_dof_index (2, component))};

                    Scalar4PointsFunc<2> function_on_cell(points, values);

                    for (st i = 0; i < 4; ++i)
                    {
                        auto indx = cell->vertex_dof_index(i, component);

                        // if ((cell->material_id() == 0) or (std::abs(grad[0][indx]) < 1e-5))
                        if (cell->material_id() == 0)
                        {
                            stress_2[g_comp][0][indx] += function_on_cell.dx(cell->vertex(i));
                            stress_2[g_comp][1][indx] += function_on_cell.dy(cell->vertex(i));
                        }
                        else
                        {
                            stress_2[g_comp][0][indx] = 0.0;
                            stress_2[g_comp][1][indx] = 0.0;
                        };
                        ++(N[indx]); 
                    };
                };
            };
        };
        for (st i = 0; i < N.size(); ++i)
        {
            stress_2[0][0][i] /= N[i];
            stress_2[0][1][i] /= N[i];
            stress_2[1][0][i] /= N[i];
            stress_2[1][1][i] /= N[i];
        };
        // for (st i = 0; i < N.size(); ++i)
        // {
        //     stress_2[0][0][i] /= epsilon;
        //     stress_2[0][1][i] /= epsilon;
        //     stress_2[1][0][i] /= epsilon;
        //     stress_2[1][1][i] /= epsilon;
        // };

        arr<dealii::Vector<dbl>, 2> mean_micro_stress;

        mean_micro_stress[0] .reinit (dof_handler.n_dofs());
        mean_micro_stress[1] .reinit (dof_handler.n_dofs());

        for (st i = 0; i < N.size(); ++i)
        {
            mean_micro_stress[0][i] = stress[0][i] + (stress_2[0][0][i] + stress_2[0][1][i]) * epsilon;
            mean_micro_stress[1][i] = stress[1][i] + (stress_2[1][0][i] + stress_2[1][1][i]) * epsilon;
        };

        {
            dealii::DataOut<2> data_out;
            data_out.attach_dof_handler (dof_handler);
            data_out.add_data_vector (mean_micro_stress[0], "dx");
            data_out.add_data_vector (mean_micro_stress[1], "dy");
            // data_out.add_data_vector (T, "T");
            data_out.build_patches ();

            auto name = file_name;
            // name += ".gpd";

            std::ofstream output (name);
            data_out.write_gnuplot (output);
        };
        {
            dealii::DataOut<2> data_out;
            data_out.attach_dof_handler (dof_handler);
            data_out.add_data_vector (stress_2[0][0], "dx_dx");
            data_out.add_data_vector (stress_2[0][1], "dx_dy");
            data_out.add_data_vector (stress_2[1][0], "dy_dx");
            data_out.add_data_vector (stress_2[1][1], "dy_dy");
            // data_out.add_data_vector (T, "T");
            data_out.build_patches ();

            auto name = str("stress_2.gpd");
            // name += ".gpd";

            std::ofstream output (name);
            data_out.write_gnuplot (output);
        };
    };

    // {
    // cu32 dofs_per_cell = dof_handler.get_fe().dofs_per_cell;
    //
    // dealii::QGauss<2>  quadrature_formula(2);
    //
    // dealii::FEValues<2> fe_values (dof_handler.get_fe(), quadrature_formula,
    //         dealii::update_gradients | dealii::update_quadrature_points | dealii::update_JxW_values);
    //
    // cst num_quad_points = quadrature_formula.size();
    //
    // arr<dealii::Vector<dbl>, 2> grad_2;
    //
    // grad_2[0] .reinit (dof_handler.get_tria().n_active_cells());
    // grad_2[1] .reinit (dof_handler.get_tria().n_active_cells());
    //
    // vec<prmt::Point<2>> midle_point (dof_handler.get_tria().n_active_cells());
    // vec<dealii::types::global_dof_index> local_dof_indices (dofs_per_cell);
    //
    // st number_of_cell = 0;
    // for (
    //         auto cell = dof_handler.begin_active(); 
    //         cell     != dof_handler.end(); 
    //         ++cell
    //     )
    // {
    //     fe_values .reinit (cell);
    //
    //     if (cell->material_id() == 0)
    //     {
    //     dbl area = 0.0;
    //     for (st q_point = 0; q_point < num_quad_points; ++q_point)
    //     {
    //         area +=
    //             fe_values.JxW(q_point);
    //     };
    //
    //     cell->get_dof_indices (local_dof_indices);
    //
    //     for (st component = 0; component < 2; ++component)
    //     {
    //         arr<dbl, 8> tmp;
    //         for (st i = 0; i < dofs_per_cell; ++i)
    //         {
    //             tmp[i] = 0.0;
    //             for (st q_point = 0; q_point < num_quad_points; ++q_point)
    //             {
    //                 tmp[i] += grad[0](local_dof_indices[i]) *
    //                     fe_values.shape_grad (i, q_point)[component] *
    //                     fe_values.JxW(q_point);
    //             };
    //         };
    //         grad_2[component][number_of_cell] = (tmp[0] + tmp[2] + tmp[4] + tmp[6]) / (area);
    //     };
    //     }
    //     else
    //     {
    //         grad[0][number_of_cell] = 0.0;
    //         grad[1][number_of_cell] = 0.0;
    //     };
    //
    //     midle_point[number_of_cell].x() = 
    //         (cell->vertex(0)(0) +
    //          cell->vertex(1)(0) +
    //          cell->vertex(2)(0) +
    //          cell->vertex(3)(0)) / 4.0;
    //     midle_point[number_of_cell].y() = 
    //         (cell->vertex(0)(1) +
    //          cell->vertex(1)(1) +
    //          cell->vertex(2)(1) +
    //          cell->vertex(3)(1)) / 4.0;
    //
    //     ++number_of_cell;
    // };
    // FILE *F;
    // F = fopen("deform_mean_2.gpd", "w");
    // for (st i = 0; i < grad_2[0].size(); ++i)
    // {
    //     fprintf(F, "%f %f %f %f\n", 
    //             midle_point[i].x(), midle_point[i].y(), grad_2[0][i], grad_2[1][i]);
    // };
    // fclose(F);
    //
    // };
}

void print_elastic_deformation_mean (const dealii::Vector<dbl> &move, 
        const dealii::DoFHandler<2> &dof_handler,
        const str file_name)
{
    cu32 dofs_per_cell = dof_handler.get_fe().dofs_per_cell;

    dealii::QGauss<2>  quadrature_formula(2);

    dealii::FEValues<2> fe_values (dof_handler.get_fe(), quadrature_formula,
            dealii::update_gradients | dealii::update_quadrature_points | dealii::update_JxW_values);

    cst num_quad_points = quadrature_formula.size();

    arr<dealii::Vector<dbl>, 2> grad;

    grad[0] .reinit (dof_handler.get_tria().n_active_cells());
    grad[1] .reinit (dof_handler.get_tria().n_active_cells());

    vec<prmt::Point<2>> midle_point (dof_handler.get_tria().n_active_cells());
    printf("!!!!!!!!!!!!!! %d\n", dofs_per_cell);
    // auto cell = dof_handler.begin_active(); 
    // printf("%d %d %d %d %d\n", 
    //         cell->vertex_dof_index (0, 0),
    //         cell->vertex_dof_index (1, 0),
    //         cell->vertex_dof_index (4, 0),
    //         cell->vertex_dof_index (4, 1),
    //         cell->vertex_dof_index (0, 1)
    //         );
    vec<dealii::types::global_dof_index> local_dof_indices (dofs_per_cell);


    st number_of_cell = 0;
    for (
            auto cell = dof_handler.begin_active(); 
            cell     != dof_handler.end(); 
            ++cell
        )
    {
        fe_values .reinit (cell);

        if (cell->material_id() == 0)
        {
            dbl area = 0.0;
            for (st q_point = 0; q_point < num_quad_points; ++q_point)
            {
                area +=
                    fe_values.JxW(q_point);
            };

            cell->get_dof_indices (local_dof_indices);

            for (st component = 0; component < 2; ++component)
            {
                arr<dbl, 8> tmp;
                for (st i = 0; i < dofs_per_cell; ++i)
                {
                    tmp[i] = 0.0;
                    for (st q_point = 0; q_point < num_quad_points; ++q_point)
                    {
                        tmp[i] += move(local_dof_indices[i]) *
                            fe_values.shape_grad (i, q_point)[component] *
                            fe_values.JxW(q_point);
                    };
                };
                grad[component][number_of_cell] = (tmp[0] + tmp[2] + tmp[4] + tmp[6]) / (area);
            };
        }
        else
        {
            grad[0][number_of_cell] = 0.0;
            grad[1][number_of_cell] = 0.0;
        };

        midle_point[number_of_cell].x() = 
            (cell->vertex(0)(0) +
             cell->vertex(1)(0) +
             cell->vertex(2)(0) +
             cell->vertex(3)(0)) / 4.0;
        midle_point[number_of_cell].y() = 
            (cell->vertex(0)(1) +
             cell->vertex(1)(1) +
             cell->vertex(2)(1) +
             cell->vertex(3)(1)) / 4.0;

        ++number_of_cell;
    };
    FILE *F;
    F = fopen("deform_mean.gpd", "w");
    for (st i = 0; i < grad[0].size(); ++i)
    {
        fprintf(F, "%f %f %f %f\n", midle_point[i].x(), midle_point[i].y(), grad[0][i], grad[1][i]);
    };
    fclose(F);

    dbl Int = 0.0;
    for (st i = 0; i < grad[0].size(); ++i)
    {
        Int += grad[0][i];
    };
    printf("Integral %f\n", Int/grad[0].size());
    //
    // {
    //     dealii::DataOut<2> data_out;
    //     data_out.attach_dof_handler (dof_handler);
    //     data_out.add_data_vector (grad[0], "grad_dx");
    //     data_out.add_data_vector (grad[1], "drad_dy");
    //     // data_out.add_data_vector (T, "T");
    //     data_out.build_patches ();
    //
    //     auto name = file_name;
    //     // name += ".gpd";
    //
    //     std::ofstream output (name);
    //     data_out.write_gnuplot (output);
    // };
}


//
// void print_elastic_deformation_mean_along_line (const dealii::Vector<dbl> &move, 
//         const dealii::DoFHandler<2> &dof_handler,
//         const str file_name,
//         Line &line)
// {
//     cu32 dofs_per_cell = dof_handler.get_fe().dofs_per_cell;
//
//     dealii::QGauss<2>  quadrature_formula(2);
//
//     dealii::FEValues<2> fe_values (dof_handler.get_fe(), quadrature_formula,
//             dealii::update_gradients | dealii::update_quadrature_points | dealii::update_JxW_values);
//
//     cst num_quad_points = quadrature_formula.size();
//
//     arr<dealii::Vector<dbl>, 2> grad;
//
//     grad[0] .reinit (dof_handler.get_tria().n_active_cells());
//     grad[1] .reinit (dof_handler.get_tria().n_active_cells());
//
//     vec<prmt::Point<2>> midle_point (dof_handler.get_tria().n_active_cells());
//     printf("!!!!!!!!!!!!!! %d\n", dofs_per_cell);
//     // auto cell = dof_handler.begin_active(); 
//     // printf("%d %d %d %d %d\n", 
//     //         cell->vertex_dof_index (0, 0),
//     //         cell->vertex_dof_index (1, 0),
//     //         cell->vertex_dof_index (4, 0),
//     //         cell->vertex_dof_index (4, 1),
//     //         cell->vertex_dof_index (0, 1)
//     //         );
//     vec<dealii::types::global_dof_index> local_dof_indices (dofs_per_cell);
//
//
//     st number_of_cell = 0;
//     for (
//             auto cell = dof_handler.begin_active(); 
//             cell     != dof_handler.end(); 
//             ++cell
//         )
//     {
//         fe_values .reinit (cell);
//
//         if (cell->material_id() == 0)
//         {
//             dbl area = 0.0;
//             for (st q_point = 0; q_point < num_quad_points; ++q_point)
//             {
//                 area +=
//                     fe_values.JxW(q_point);
//             };
//
//             cell->get_dof_indices (local_dof_indices);
//
//             for (st component = 0; component < 2; ++component)
//             {
//                 arr<dbl, 8> tmp;
//                 for (st i = 0; i < dofs_per_cell; ++i)
//                 {
//                     tmp[i] = 0.0;
//                     for (st q_point = 0; q_point < num_quad_points; ++q_point)
//                     {
//                         tmp[i] += move(local_dof_indices[i]) *
//                             fe_values.shape_grad (i, q_point)[component] *
//                             fe_values.JxW(q_point);
//                     };
//                 };
//                 grad[component][number_of_cell] = (tmp[0] + tmp[2] + tmp[4] + tmp[6]) / (area);
//             };
//         }
//         else
//         {
//             grad[0][number_of_cell] = 0.0;
//             grad[1][number_of_cell] = 0.0;
//         };
//
//         midle_point[number_of_cell].x() = 
//             (cell->vertex(0)(0) +
//              cell->vertex(1)(0) +
//              cell->vertex(2)(0) +
//              cell->vertex(3)(0)) / 4.0;
//         midle_point[number_of_cell].y() = 
//             (cell->vertex(0)(1) +
//              cell->vertex(1)(1) +
//              cell->vertex(2)(1) +
//              cell->vertex(3)(1)) / 4.0;
//
//         ++number_of_cell;
//     };
//     FILE *F;
//     F = fopen("deform_mean.gpd", "w");
//     for (st i = 0; i < grad[0].size(); ++i)
//     {
//         fprintf(F, "%f %f %f %f\n", midle_point[i].x(), midle_point[i].y(), grad[0][i], grad[1][i]);
//     };
//     fclose(F);
//
//     dbl Int = 0.0;
//     for (st i = 0; i < grad[0].size(); ++i)
//     {
//         Int += grad[0][i];
//     };
//     printf("Integral %f\n", Int/grad[0].size());
//     //
//     // {
//     //     dealii::DataOut<2> data_out;
//     //     data_out.attach_dof_handler (dof_handler);
//     //     data_out.add_data_vector (grad[0], "grad_dx");
//     //     data_out.add_data_vector (grad[1], "drad_dy");
//     //     // data_out.add_data_vector (T, "T");
//     //     data_out.build_patches ();
//     //
//     //     auto name = file_name;
//     //     // name += ".gpd";
//     //
//     //     std::ofstream output (name);
//     //     data_out.write_gnuplot (output);
//     // };
// }

void print_elastic_deformation_mean_other (const dealii::Vector<dbl> &move, 
        const dealii::DoFHandler<2> &dof_handler,
        const str file_name)
{
    cu32 dofs_per_cell = dof_handler.get_fe().dofs_per_cell;

    arr<dealii::Vector<dbl>, 2> grad;

    grad[0] .reinit (dof_handler.n_dofs());
    grad[1] .reinit (dof_handler.n_dofs());

    vec<dbl> area(dof_handler.n_dofs());

    vec<st> N(dof_handler.n_dofs());

    dealii::QGauss<2> quadrature_formula(2);

    dealii::FEValues<2> fe_values (dof_handler.get_fe(), quadrature_formula,
            dealii::update_gradients | dealii::update_values | dealii::update_quadrature_points | dealii::update_JxW_values);

    cst num_quad_points = quadrature_formula.size();

    // vec<prmt::Point<2>> midle_point (dof_handler.get_tria().n_active_cells());
    // printf("!!!!!!!!!!!!!! %d\n", dofs_per_cell);
    // auto cell = dof_handler.begin_active(); 
    // printf("%d %d %d %d %d\n", 
    //         cell->vertex_dof_index (0, 0),
    //         cell->vertex_dof_index (1, 0),
    //         cell->vertex_dof_index (4, 0),
    //         cell->vertex_dof_index (4, 1),
    //         cell->vertex_dof_index (0, 1)
    //         );
    vec<dealii::types::global_dof_index> local_dof_indices (dofs_per_cell);


    st number_of_cell = 0;
    for (
            auto cell = dof_handler.begin_active(); 
            cell     != dof_handler.end(); 
            ++cell
        )
    {
        fe_values .reinit (cell);

        if (cell->material_id() == 0)
        {
            dbl cell_area = 0.0;
            for (st q_point = 0; q_point < num_quad_points; ++q_point)
            {
                cell_area +=
                    fe_values.JxW(q_point);
            };

            cell->get_dof_indices (local_dof_indices);

            arr<dbl, 8> tmp;
            for (st i = 0; i < dofs_per_cell; ++i)
            {
                tmp[i] = 0.0;
                for (st q_point = 0; q_point < num_quad_points; ++q_point)
                {
                    // grad[0][local_dof_indices[i]] += 
                    //     // cell_area;
                    //     move(local_dof_indices[i]) *
                    //     // fe_values.shape_grad (i, q_point)[0] *
                    //     fe_values.shape_value (i, q_point) *
                    //     fe_values.JxW(q_point);
                    //     // 1.0;
                    // grad[1][local_dof_indices[i]] += move(local_dof_indices[i]) *
                    //     fe_values.shape_grad (i, q_point)[1] *
                    //     fe_values.JxW(q_point);
                    tmp[i] += 
                        // cell_area;
                        move(local_dof_indices[i]) *
                        fe_values.shape_grad (i, q_point)[0] *
                        // fe_values.shape_value (i, q_point) *
                        fe_values.JxW(q_point);
                    // 1.0;
                    // tmp[i] += move(local_dof_indices[i]) *
                    //     fe_values.shape_grad (i, q_point)[1] *
                    //     fe_values.JxW(q_point);
                };
                // area[local_dof_indices[i]] += cell_area;
                ++(N[local_dof_indices[i]]); 
            };
            grad[0][local_dof_indices[0]] += (tmp[0] + tmp[2] + tmp[4] + tmp[6]) / (cell_area);
            grad[0][local_dof_indices[2]] += (tmp[0] + tmp[2] + tmp[4] + tmp[6]) / (cell_area);
            grad[0][local_dof_indices[4]] += (tmp[0] + tmp[2] + tmp[4] + tmp[6]) / (cell_area);
            grad[0][local_dof_indices[6]] += (tmp[0] + tmp[2] + tmp[4] + tmp[6]) / (cell_area);
        }
        else
        {
            grad[0][number_of_cell] = 0.0;
            grad[1][number_of_cell] = 0.0;
        };
    };
    for (st i = 0; i < area.size(); ++i)
    {
        // area[i] /= N[i];
        // area[i] /= N[i];
        // grad[0][i] /= area[i];
        // grad[1][i] /= area[i];
        grad[0][i] /= N[i];
        grad[1][i] /= N[i];
    };
    // FILE *F;
    // F = fopen("deform_mean.gpd", "w");
    // for (st i = 0; i < grad[0].size(); ++i)
    // {
    //     fprintf(F, "%f %f %f %f\n", midle_point[i].x(), midle_point[i].y(), grad[0][i], grad[1][i]);
    // };
    // fclose(F);
    //
    // dbl Int = 0.0;
    // for (st i = 0; i < grad[0].size(); ++i)
    // {
    //     Int += grad[0][i];
    // };
    // printf("Integral %f\n", Int/grad[0].size());
    //
    {
        dealii::DataOut<2> data_out;
        data_out.attach_dof_handler (dof_handler);
        data_out.add_data_vector (grad[0], "grad_dx");
        data_out.add_data_vector (grad[1], "drad_dy");
        // data_out.add_data_vector (T, "T");
        data_out.build_patches ();

        auto name = file_name;
        // name += ".gpd";

        std::ofstream output (name);
        data_out.write_gnuplot (output);
    };
}

void print_elastic_deformation_mean_other_3d_binary (const dealii::Vector<dbl> &move, 
        const dealii::DoFHandler<3> &dof_handler,
        const str file_name)
{
    cu32 dofs_per_cell = dof_handler.get_fe().dofs_per_cell;

    arr<dealii::Vector<dbl>, 3> grad;

    grad[0] .reinit (dof_handler.n_dofs());
    grad[1] .reinit (dof_handler.n_dofs());
    grad[2] .reinit (dof_handler.n_dofs());

    // vec<dbl> area(dof_handler.n_dofs());

    vec<st> N(dof_handler.n_dofs());

    dealii::QGauss<3> quadrature_formula(2);

    dealii::FEValues<3> fe_values (dof_handler.get_fe(), quadrature_formula,
            dealii::update_gradients | dealii::update_values | dealii::update_quadrature_points | dealii::update_JxW_values);

    cst num_quad_points = quadrature_formula.size();

    // vec<prmt::Point<2>> midle_point (dof_handler.get_tria().n_active_cells());
    // printf("!!!!!!!!!!!!!! %d\n", dofs_per_cell);
    // auto cell = dof_handler.begin_active(); 
    // printf("%d %d %d %d %d\n", 
    //         cell->vertex_dof_index (0, 0),
    //         cell->vertex_dof_index (1, 0),
    //         cell->vertex_dof_index (4, 0),
    //         cell->vertex_dof_index (4, 1),
    //         cell->vertex_dof_index (0, 1)
    //         );
    vec<dealii::types::global_dof_index> global_dof_indices (dofs_per_cell);


    st number_of_cell = 0;
    for (
            auto cell = dof_handler.begin_active(); 
            cell     != dof_handler.end(); 
            ++cell
        )
    {
        fe_values .reinit (cell);

        if (cell->material_id() == 0)
        {
            dbl cell_area = 0.0;
            for (st q_point = 0; q_point < num_quad_points; ++q_point)
            {
                cell_area +=
                    fe_values.JxW(q_point);
            };

            cell->get_dof_indices (global_dof_indices);

            arr<arr<dbl, 3>, 3> tmp;
            for (st i = 0; i < 3; ++i)
            {
                for (st j = 0; j < 3; ++j)
                {
                    tmp[i][j] = 0.0;
                };
            };
            for (st i = 0; i < 3; ++i)
            {
                for (st j = 0; j < 24; ++j)
                {
                    for (st q_point = 0; q_point < num_quad_points; ++q_point)
                    {
                        tmp[i][j % 3] += 
                            move(global_dof_indices[j]) *
                            fe_values.shape_grad (j, q_point)[i] *
                            // fe_values.shape_value (i, q_point) *
                            fe_values.JxW(q_point);
                        // 1.0;
                        // tmp[i] += move(local_dof_indices[i]) *
                        //     fe_values.shape_grad (i, q_point)[1] *
                        //     fe_values.JxW(q_point);
                    };
                    // area[local_dof_indices[i]] += cell_area;
                    ++(N[global_dof_indices[j]]); 
                };
            };
            for (st i = 0; i < 3; ++i)
            {
                for (st j = 0; j < 3; ++j)
                {
                    tmp[i][j] /= cell_area;
                };
            };
            // grad[0][local_dof_indices[0]] += (tmp[0] + tmp[2] + tmp[4] + tmp[6]) / (cell_area);
            // grad[0][local_dof_indices[2]] += (tmp[0] + tmp[2] + tmp[4] + tmp[6]) / (cell_area);
            // grad[0][local_dof_indices[4]] += (tmp[0] + tmp[2] + tmp[4] + tmp[6]) / (cell_area);
            // grad[0][local_dof_indices[6]] += (tmp[0] + tmp[2] + tmp[4] + tmp[6]) / (cell_area);
            for (st i = 0; i < 3; ++i)
            {
                for (st j = 0; j < 8; ++j)
                {
                    for (st k = 0; k < 3; ++k)
                    {
                        grad[i][global_dof_indices[j*3+k]] += tmp[i][k];
                    };
                };
            };
        }
        else
        {
            for (st i = 0; i < dofs_per_cell; ++i)
            {
                grad[0][global_dof_indices[i]] = 0.0;
            };
        };
        cell->get_dof_indices (global_dof_indices);
        for (st i = 0; i < 24; ++i)
        {
            grad[0][global_dof_indices[i]] =  cell->material_id();
        };
    };
    for (st i = 0; i < dof_handler.n_dofs(); ++i)
    {
        // area[i] /= N[i];
        // area[i] /= N[i];
        // grad[0][i] /= area[i];
        // grad[1][i] /= area[i];
        // grad[0][i] /= (N[i] / 3);
        // grad[1][i] /= (N[i] / 3);
        // grad[2][i] /= (N[i] / 3);
    };
    // FILE *F;
    // F = fopen("deform_mean.gpd", "w");
    // for (st i = 0; i < grad[0].size(); ++i)
    // {
    //     fprintf(F, "%f %f %f %f\n", midle_point[i].x(), midle_point[i].y(), grad[0][i], grad[1][i]);
    // };
    // fclose(F);
    //
    // dbl Int = 0.0;
    // for (st i = 0; i < grad[0].size(); ++i)
    // {
    //     Int += grad[0][i];
    // };
    // printf("Integral %f\n", Int/grad[0].size());
    //
    print_move_slice (grad[0], dof_handler, file_name, 2, 0.5);
    // {
    //     dealii::DataOut<3> data_out;
    //     data_out.attach_dof_handler (dof_handler);
    //     data_out.add_data_vector (grad[0], "grad_dx");
    //     // data_out.add_data_vector (grad[1], "drad_dy");
    //     // data_out.add_data_vector (grad[2], "drad_dz");
    //     // data_out.add_data_vector (T, "T");
    //     data_out.build_patches ();
    //
    //     auto name = file_name;
    //     // name += ".gpd";
    //
    //     std::ofstream output (name);
    //     data_out.write_gnuplot (output);
    // };
}

void print_elastic_deformation_quad (const dealii::Vector<dbl> &move, 
        const dealii::DoFHandler<2> &dof_handler,
        const str file_name)
{
    cu32 dofs_per_cell = dof_handler.get_fe().dofs_per_cell;

    arr<dealii::Vector<dbl>, 2> grad;

    grad[0] .reinit (dof_handler.n_dofs());
    grad[1] .reinit (dof_handler.n_dofs());

    vec<st> N(dof_handler.n_dofs());

    for (
            auto cell = dof_handler.begin_active(); 
            cell     != dof_handler.end(); 
            ++cell
        )
    {
        double x1 = cell->vertex(0)(0);
        double x2 = cell->vertex(1)(0);
        double x3 = cell->vertex(2)(0);
        double x4 = cell->vertex(3)(0);

        double y1 = cell->vertex(0)(1);
        double y2 = cell->vertex(1)(1);
        double y3 = cell->vertex(2)(1);
        double y4 = cell->vertex(3)(1);

        for (st component = 0; component < 2; ++component)
        {
            double f1 = move(cell->vertex_dof_index (0, component));
            double f2 = move(cell->vertex_dof_index (1, component));
            double f3 = move(cell->vertex_dof_index (2, component));
            double f4 = move(cell->vertex_dof_index (3, component));
            for (st i = 0; i < 4; ++i)
            {
                auto indx = cell->vertex_dof_index(i, component);


                if ((cell->material_id() == 0) or (cell->material_id() == 2))
                    // if (cell->material_id() == 0)
                {
                    dbl b=-(x1*y1*y2*f3-x1*y1*y2*f4-x1*y1*f3*y4+x1*y1*y4*f2-x1*y1*f2*y3+x1*y1*y3*f4+y3*x3*y2*f4-y2*x3*y3*f1+y3*x4*y4*f2-y3*x4*y4*f1+x3*y3*f1*y4-x3*y3*f2*y4+f3*x2*y2*y4-x2*y2*f1*y4+x2*y2*f1*y3-y3*x2*y2*f4-y1*x4*y4*f2-y1*y3*x3*f4+y1*x2*y2*f4+y1*x3*y3*f2-y1*x2*y2*f3-f3*y2*x4*y4+y1*f3*y4*x4+f1*y2*x4*y4)/(-y1*y3*x3*x2-x4*x2*y2*y1-x3*y3*x4*y2+x2*y3*x3*y4-x2*y2*x3*y4+x2*y1*x4*y4-x2*y3*x4*y4+y3*x4*x2*y2+y1*x3*x2*y2+y2*x3*x4*y4-y1*x3*x4*y4+x1*y1*x3*y4+y2*x1*y3*x3-x1*y4*x2*y1+x1*x2*y2*y4-y2*x3*x1*y1+y3*y1*x3*x4+x1*y3*x4*y4-x1*y3*x2*y2+x1*y1*y3*x2-x1*x3*y3*y4-x1*y3*x4*y1-x1*y2*x4*y4+y2*x4*x1*y1);
                    dbl c=(x1*x2*y2*f4-x1*f3*x2*y2+x3*x1*y1*f4-x1*y1*x2*f4+x1*y1*x2*f3-x1*x4*y4*f2+x4*x1*y1*f2+x1*f3*y4*x4-x4*x1*y1*f3+x1*x3*y3*f2-x3*x1*y1*f2-x1*y3*x3*f4-x3*x2*y2*f4+x3*x2*y2*f1-x4*x2*y2*f1+x4*x2*y2*f3-f3*y4*x4*x2+x2*x4*y4*f1+y3*x3*x2*f4-x4*x3*y3*f2-x2*x3*y3*f1+x3*x4*y4*f2-x3*x4*y4*f1+x4*x3*y3*f1)/(-y1*y3*x3*x2-x4*x2*y2*y1-x3*y3*x4*y2+x2*y3*x3*y4-x2*y2*x3*y4+x2*y1*x4*y4-x2*y3*x4*y4+y3*x4*x2*y2+y1*x3*x2*y2+y2*x3*x4*y4-y1*x3*x4*y4+x1*y1*x3*y4+y2*x1*y3*x3-x1*y4*x2*y1+x1*x2*y2*y4-y2*x3*x1*y1+y3*y1*x3*x4+x1*y3*x4*y4-x1*y3*x2*y2+x1*y1*y3*x2-x1*x3*y3*y4-x1*y3*x4*y1-x1*y2*x4*y4+y2*x4*x1*y1);
                    dbl d=(-x3*y1*f4+x3*y1*f2+y1*f3*x4-x4*y1*f2-x1*y2*f4+x3*y2*f4-x3*f1*y2+y2*x1*f3+f1*y2*x4-f3*y2*x4+x4*y3*f2-x1*y3*f2-x4*y3*f1+x1*y4*f2-x3*y4*f2+x3*y4*f1-x2*f1*y4+f3*x2*y4+x2*f1*y3-y3*x2*f4+x1*y3*f4-x1*f3*y4+y1*x2*f4-y1*x2*f3)/(-y1*y3*x3*x2-x4*x2*y2*y1-x3*y3*x4*y2+x2*y3*x3*y4-x2*y2*x3*y4+x2*y1*x4*y4-x2*y3*x4*y4+y3*x4*x2*y2+y1*x3*x2*y2+y2*x3*x4*y4-y1*x3*x4*y4+x1*y1*x3*y4+y2*x1*y3*x3-x1*y4*x2*y1+x1*x2*y2*y4-y2*x3*x1*y1+y3*y1*x3*x4+x1*y3*x4*y4-x1*y3*x2*y2+x1*y1*y3*x2-x1*x3*y3*y4-x1*y3*x4*y1-x1*y2*x4*y4+y2*x4*x1*y1);
                    grad[0][indx] += b + d * cell->vertex(i)(1);
                    grad[1][indx] += c + d * cell->vertex(i)(0);
                }
                else
                {
                    grad[0][indx] = 0.0;
                    grad[1][indx] = 0.0;
                };
                ++(N[indx]); 
            };
        };
    };
    for (st i = 0; i < N.size(); ++i)
    {
        grad[0][i] /= N[i];
        grad[1][i] /= N[i];
    };

    {
        dealii::DataOut<2> data_out;
        data_out.attach_dof_handler (dof_handler);
        data_out.add_data_vector (grad[0], "grad_dx");
        data_out.add_data_vector (grad[1], "drad_dy");
        // data_out.add_data_vector (T, "T");
        data_out.build_patches ();

        auto name = file_name;
        // name += ".gpd";

        std::ofstream output (name);
        data_out.write_gnuplot (output);
    };
};
void print_elastic_stress_quad (const dealii::Vector<dbl> &move, 
        const dealii::DoFHandler<2> &dof_handler,
        const ATools::FourthOrderTensor &C,
        const str file_name)
{
    cu32 dofs_per_cell = dof_handler.get_fe().dofs_per_cell;

    arr<dealii::Vector<dbl>, 2> grad;

    grad[0] .reinit (dof_handler.n_dofs());
    grad[1] .reinit (dof_handler.n_dofs());

    vec<st> N(dof_handler.n_dofs());

    for (
            auto cell = dof_handler.begin_active(); 
            cell     != dof_handler.end(); 
            ++cell
        )
    {
        double x1 = cell->vertex(0)(0);
        double x2 = cell->vertex(1)(0);
        double x3 = cell->vertex(2)(0);
        double x4 = cell->vertex(3)(0);

        double y1 = cell->vertex(0)(1);
        double y2 = cell->vertex(1)(1);
        double y3 = cell->vertex(2)(1);
        double y4 = cell->vertex(3)(1);

        for (st component = 0; component < 2; ++component)
        {
            double f1 = move(cell->vertex_dof_index (0, component));
            double f2 = move(cell->vertex_dof_index (1, component));
            double f3 = move(cell->vertex_dof_index (2, component));
            double f4 = move(cell->vertex_dof_index (3, component));
            for (st i = 0; i < 4; ++i)
            {
                auto indx = cell->vertex_dof_index(i, component);


                if ((cell->material_id() == 0) or (cell->material_id() == 2))
                    // if (cell->material_id() == 0)
                {
                    dbl b=-(x1*y1*y2*f3-x1*y1*y2*f4-x1*y1*f3*y4+x1*y1*y4*f2-x1*y1*f2*y3+x1*y1*y3*f4+y3*x3*y2*f4-y2*x3*y3*f1+y3*x4*y4*f2-y3*x4*y4*f1+x3*y3*f1*y4-x3*y3*f2*y4+f3*x2*y2*y4-x2*y2*f1*y4+x2*y2*f1*y3-y3*x2*y2*f4-y1*x4*y4*f2-y1*y3*x3*f4+y1*x2*y2*f4+y1*x3*y3*f2-y1*x2*y2*f3-f3*y2*x4*y4+y1*f3*y4*x4+f1*y2*x4*y4)/(-y1*y3*x3*x2-x4*x2*y2*y1-x3*y3*x4*y2+x2*y3*x3*y4-x2*y2*x3*y4+x2*y1*x4*y4-x2*y3*x4*y4+y3*x4*x2*y2+y1*x3*x2*y2+y2*x3*x4*y4-y1*x3*x4*y4+x1*y1*x3*y4+y2*x1*y3*x3-x1*y4*x2*y1+x1*x2*y2*y4-y2*x3*x1*y1+y3*y1*x3*x4+x1*y3*x4*y4-x1*y3*x2*y2+x1*y1*y3*x2-x1*x3*y3*y4-x1*y3*x4*y1-x1*y2*x4*y4+y2*x4*x1*y1);
                    dbl c=(x1*x2*y2*f4-x1*f3*x2*y2+x3*x1*y1*f4-x1*y1*x2*f4+x1*y1*x2*f3-x1*x4*y4*f2+x4*x1*y1*f2+x1*f3*y4*x4-x4*x1*y1*f3+x1*x3*y3*f2-x3*x1*y1*f2-x1*y3*x3*f4-x3*x2*y2*f4+x3*x2*y2*f1-x4*x2*y2*f1+x4*x2*y2*f3-f3*y4*x4*x2+x2*x4*y4*f1+y3*x3*x2*f4-x4*x3*y3*f2-x2*x3*y3*f1+x3*x4*y4*f2-x3*x4*y4*f1+x4*x3*y3*f1)/(-y1*y3*x3*x2-x4*x2*y2*y1-x3*y3*x4*y2+x2*y3*x3*y4-x2*y2*x3*y4+x2*y1*x4*y4-x2*y3*x4*y4+y3*x4*x2*y2+y1*x3*x2*y2+y2*x3*x4*y4-y1*x3*x4*y4+x1*y1*x3*y4+y2*x1*y3*x3-x1*y4*x2*y1+x1*x2*y2*y4-y2*x3*x1*y1+y3*y1*x3*x4+x1*y3*x4*y4-x1*y3*x2*y2+x1*y1*y3*x2-x1*x3*y3*y4-x1*y3*x4*y1-x1*y2*x4*y4+y2*x4*x1*y1);
                    dbl d=(-x3*y1*f4+x3*y1*f2+y1*f3*x4-x4*y1*f2-x1*y2*f4+x3*y2*f4-x3*f1*y2+y2*x1*f3+f1*y2*x4-f3*y2*x4+x4*y3*f2-x1*y3*f2-x4*y3*f1+x1*y4*f2-x3*y4*f2+x3*y4*f1-x2*f1*y4+f3*x2*y4+x2*f1*y3-y3*x2*f4+x1*y3*f4-x1*f3*y4+y1*x2*f4-y1*x2*f3)/(-y1*y3*x3*x2-x4*x2*y2*y1-x3*y3*x4*y2+x2*y3*x3*y4-x2*y2*x3*y4+x2*y1*x4*y4-x2*y3*x4*y4+y3*x4*x2*y2+y1*x3*x2*y2+y2*x3*x4*y4-y1*x3*x4*y4+x1*y1*x3*y4+y2*x1*y3*x3-x1*y4*x2*y1+x1*x2*y2*y4-y2*x3*x1*y1+y3*y1*x3*x4+x1*y3*x4*y4-x1*y3*x2*y2+x1*y1*y3*x2-x1*x3*y3*y4-x1*y3*x4*y1-x1*y2*x4*y4+y2*x4*x1*y1);
                    grad[0][indx] += b + d * cell->vertex(i)(1);
                    grad[1][indx] += c + d * cell->vertex(i)(0);
                }
                else
                {
                    grad[0][indx] = 0.0;
                    grad[1][indx] = 0.0;
                };
                ++(N[indx]); 
            };
        };
    };
    for (st i = 0; i < N.size(); ++i)
    {
        grad[0][i] /= N[i];
        grad[1][i] /= N[i];
    };

    arr<dealii::Vector<dbl>, 2> stress;

    stress[0] .reinit (dof_handler.n_dofs());
    stress[1] .reinit (dof_handler.n_dofs());

    for (
            auto cell = dof_handler.begin_active(); 
            cell     != dof_handler.end(); 
            ++cell
        )
    {
        for (st i = 0; i < 4; ++i)
        {
            auto indx_1 = cell->vertex_dof_index(i, 0);
            auto indx_2 = cell->vertex_dof_index(i, 1);
            stress[0][indx_1] = 
                C[0][0][0][0] * grad[0][indx_1] +
                C[0][0][0][1] * grad[0][indx_2] +
                C[0][0][1][0] * grad[1][indx_1] +
                C[0][0][1][1] * grad[1][indx_2];
            stress[0][indx_2] = 
                C[0][1][0][0] * grad[0][indx_1] +
                C[0][1][0][1] * grad[0][indx_2] +
                C[0][1][1][0] * grad[1][indx_1] +
                C[0][1][1][1] * grad[1][indx_2];
            stress[1][indx_1] = 
                C[1][0][0][0] * grad[0][indx_1] +
                C[1][0][0][1] * grad[0][indx_2] +
                C[1][0][1][0] * grad[1][indx_1] +
                C[1][0][1][1] * grad[1][indx_2];
            stress[1][indx_2] = 
                C[1][1][0][0] * grad[0][indx_1] +
                C[1][1][0][1] * grad[0][indx_2] +
                C[1][1][1][0] * grad[1][indx_1] +
                C[1][1][1][1] * grad[1][indx_2];
        };
    };

    // {
    //     dealii::DataOut<2> data_out;
    //     data_out.attach_dof_handler (dof_handler);
    //     data_out.add_data_vector (grad[0], "grad_dx");
    //     data_out.add_data_vector (grad[1], "drad_dy");
    //     // data_out.add_data_vector (T, "T");
    //     data_out.build_patches ();
    //
    //     auto name = file_name;
    //     // name += ".gpd";
    //
    //     std::ofstream output (name);
    //     data_out.write_gnuplot (output);
    // };

    {
        dealii::DataOut<2> data_out;
        data_out.attach_dof_handler (dof_handler);
        data_out.add_data_vector (stress[0], "grad_dx");
        data_out.add_data_vector (stress[1], "drad_dy");
        // data_out.add_data_vector (T, "T");
        data_out.build_patches ();

        auto name = file_name;
        // name += ".gpd";

        std::ofstream output (name);
        data_out.write_gnuplot (output);
    };
};

    template <st dim>
void print_coor_bin (const dealii::DoFHandler<dim> &dof_handler, const str file_name)
{
    vec<dealii::Point<dim>> coor(dof_handler.n_dofs());
    {
        cu8 dofs_per_cell =  dof_handler.get_fe().dofs_per_cell;

        std::vector<u32> local_dof_indices (dofs_per_cell);

        auto cell = dof_handler.begin_active();
        auto endc = dof_handler.end();
        for (; cell != endc; ++cell)
        {
            cell ->get_dof_indices (local_dof_indices);

            // FOR (i, 0, dofs_per_cell)
            //     indexes(local_dof_indices[i]) = cell ->vertex_dof_index (i, 0);
            FOR (i, 0, dofs_per_cell / dim)
            {
                for (st j = 0; j < dim; ++j)
                {
                    coor[local_dof_indices[i*dim+j]] = cell ->vertex (i);
                };
            };
        };
    };

    {
        std::ofstream out (file_name, std::ios::out | std::ios::binary);
        for (st i = 0; i < coor.size(); ++i)
        {
            for (st j = 0; j < dim; ++j)
            {
                out.write ((char *) &(coor[i](j)), sizeof(dbl));
            };
        };
        out.close ();
    };
};

};

#endif
