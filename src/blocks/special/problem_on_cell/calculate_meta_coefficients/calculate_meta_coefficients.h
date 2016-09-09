#ifndef CALCULATION_META_COEFFICIENTS_ON_CELL
#define CALCULATION_META_COEFFICIENTS_ON_CELL

#include "../../../../../../prmt_sintactic_addition/prmt_sintactic_addition.h"
#include "../../../general/additional_tools/types/types.h"
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/lac/vector.h>

namespace OnCell
{
    //! Расчет среднего коэффициентов
    template<u8 dim>
    dbl calculate_area_of_domain(
            const dealii::DoFHandler<dim> &dof_handler)
    {
        dbl area_of_domain = 0.0;

        {
            dealii::QGauss<dim>  quadrature_formula(2);

            dealii::FEValues<dim> fe_values (dof_handler.get_fe(), quadrature_formula,
                    dealii::update_quadrature_points | dealii::update_JxW_values);

            cst n_q_points = quadrature_formula.size();

            auto cell = dof_handler.begin_active();
            auto endc = dof_handler.end();
            for (; cell != endc; ++cell)
            {
                fe_values .reinit (cell);

                for (st q_point = 0; q_point < n_q_points; ++q_point)
                    area_of_domain += fe_values.JxW(q_point);
            };
        };
        return  area_of_domain;
    };

    //! Расчет среднего коэффициентов
    template<u8 dim>
    ATools::SecondOrderTensor calculate_mean_coefficients(
            const dealii::DoFHandler<dim> &dof_handler,
            const vec<ATools::SecondOrderTensor> &coef)
    {
        ATools::SecondOrderTensor mean_coefficient;
        for (st i = 0; i < dim; ++i)
            for (st j = 0; j < dim; ++j)
                mean_coefficient[i][j] = 0.0;
        dbl area_of_domain = 0.0;

        {
            dealii::QGauss<dim>  quadrature_formula(2);

            dealii::FEValues<dim> fe_values (dof_handler.get_fe(), quadrature_formula,
                    dealii::update_quadrature_points | dealii::update_JxW_values);

            cst n_q_points = quadrature_formula.size();

            auto cell = dof_handler.begin_active();
            auto endc = dof_handler.end();
            for (; cell != endc; ++cell)
            {
                fe_values .reinit (cell);

                for (st q_point = 0; q_point < n_q_points; ++q_point)
                    for (st i = 0; i < dim; ++i)
                        for (st j = 0; j < dim; ++j)
                        mean_coefficient[i][j] += 
                            coef[cell->material_id()][i][j] *
                            fe_values.JxW(q_point);

                for (st q_point = 0; q_point < n_q_points; ++q_point)
                    area_of_domain += fe_values.JxW(q_point);
            };

            for (st i = 0; i < dim; ++i)
                for (st j = 0; j < dim; ++j)
                    mean_coefficient[i][j] /= area_of_domain;
        };
        return  mean_coefficient;
    };

    //! Расчет среднего коэффициентов
    template<u8 dim>
    ATools::FourthOrderTensor calculate_mean_coefficients(
            const dealii::DoFHandler<dim> &dof_handler,
            const vec<ATools::FourthOrderTensor> &coef)
    {
        ATools::FourthOrderTensor mean_coefficient;
        for (st i = 0; i < dim; ++i)
            for (st j = 0; j < dim; ++j)
                for (st k = 0; k < dim; ++k)
                    for (st l = 0; l < dim; ++l)
                        mean_coefficient[i][j][k][l] = 0.0;
        dbl area_of_domain = 0.0;

        {
            dealii::QGauss<dim>  quadrature_formula(2);

            dealii::FEValues<dim> fe_values (dof_handler.get_fe(), quadrature_formula,
                    dealii::update_quadrature_points | dealii::update_JxW_values);

            cst n_q_points = quadrature_formula.size();

            auto cell = dof_handler.begin_active();
            auto endc = dof_handler.end();
            for (; cell != endc; ++cell)
            {
                fe_values .reinit (cell);

                for (st q_point = 0; q_point < n_q_points; ++q_point)
                    for (st i = 0; i < dim; ++i)
                        for (st j = 0; j < dim; ++j)
                            for (st k = 0; k < dim; ++k)
                                for (st l = 0; l < dim; ++l)
                                    mean_coefficient[i][j][k][l] += 
                                        coef[cell->material_id()][i][j][k][l] *
                                        fe_values.JxW(q_point);

                for (st q_point = 0; q_point < n_q_points; ++q_point)
                    area_of_domain += fe_values.JxW(q_point);
            };

            for (st i = 0; i < dim; ++i)
                for (st j = 0; j < dim; ++j)
                    for (st k = 0; k < dim; ++k)
                        for (st l = 0; l < dim; ++l)
                            mean_coefficient[i][j][k][l] /= area_of_domain;
        };
        return  mean_coefficient;
    };

    //! Расчет метакоэффициентов, скалярный случай
    template<u8 dim>
    ATools::SecondOrderTensor calculate_meta_coefficients_scalar(
            const dealii::DoFHandler<dim> &dof_handler,
            const arr<dealii::Vector<dbl>, dim> &solution,
            const arr<dealii::Vector<dbl>, dim> &heat_flow,
            const vec<ATools::SecondOrderTensor> &coef)
    {
        ATools::SecondOrderTensor mean_coefficient;
        for (st i = 0; i < dim; ++i)
            for (st j = 0; j < dim; ++j)
                mean_coefficient[i][j] = 0.0;
        dbl area_of_domain = 0.0;

        {
            dealii::QGauss<dim>  quadrature_formula(2);

            dealii::FEValues<dim> fe_values (dof_handler.get_fe(), quadrature_formula,
                    dealii::update_quadrature_points | dealii::update_JxW_values);

            cst n_q_points = quadrature_formula.size();


            auto cell = dof_handler.begin_active();
            auto endc = dof_handler.end();
            for (; cell != endc; ++cell)
            {
                fe_values .reinit (cell);

                for (st q_point = 0; q_point < n_q_points; ++q_point)
                    for (st i = 0; i < dim; ++i)
                        for (st j = 0; j < dim; ++j)
                        mean_coefficient[i][j] += 
                            coef[cell->material_id()][i][j] *
                            fe_values.JxW(q_point);

                for (st q_point = 0; q_point < n_q_points; ++q_point)
                    area_of_domain += fe_values.JxW(q_point);
            };

            for (st i = 0; i < dim; ++i)
                for (st j = 0; j < dim; ++j)
                    mean_coefficient[i][j] /= area_of_domain;
        };

        st len_vector_solution = dof_handler.n_dofs();
        ATools::SecondOrderTensor mean_heat_flow;
        ATools::SecondOrderTensor meta_coefficient;

        for (st i = 0; i < dim; ++i)
            for (st j = 0; j < dim; ++j)
            {
                mean_heat_flow[i][j] = 0.0;

                for (st k = 0; k < len_vector_solution; ++k)
                    mean_heat_flow[i][j] += solution[i](k) * (-heat_flow[j](k));

                mean_heat_flow[i][j] /= area_of_domain;

                meta_coefficient[i][j] = mean_coefficient[i][j] + mean_heat_flow[i][j];
            };

        // mean_heat_flow[1][1] = mean_coefficient[1][1];
        // mean_heat_flow[0][1] = meta_coefficient[1][1];
        // return  mean_heat_flow; 
        return  meta_coefficient;
    };

    // //! Расчет теплового потока
    // template<u8 dim>
    // void calculate_heat_flow(
    //         const arr<dealii::Vector<dbl>, dim> &solution,
    //         const arr<dealii::Vector<dbl>, dim> &rhs,
    //         const vec<ATools::SecondOrderTensor> &coef,
    //         arr<arr<dealii::Vector<dbl>, dim>, dim> &heat_flow)
    // {
    //     for (st i = 0; i < dim; ++i)
    //         for (st j = 0; j < dim; ++j)
    //         {
    //             heat_flow[i][j].reinit(solution[i].size());
    //
    //             for (st k = 0; k < solution[i].size(); ++k)
    //                 heat_flow[i][j](k) = solution[i](k) * (-rhs[j](k));// + coef[i][j];
    //         };
    // };

    //! Расчет метакоэффициентов, плоская задача упругости
    template<u8 dim>
    ATools::FourthOrderTensor calculate_meta_coefficients_2d_elastic(
            const dealii::DoFHandler<dim> &dof_handler,
            const SystemsLinearAlgebraicEquations<4> &elastic_slae,
            const SystemsLinearAlgebraicEquations<2> &heat_slae,
            const vec<ATools::FourthOrderTensor> &coef)
    {
        ATools::FourthOrderTensor mean_coefficient;
        for (st i = 0; i < 3; ++i)
            for (st j = 0; j < 3; ++j)
                for (st k = 0; k < 3; ++k)
                    for (st l = 0; l < 3; ++l)
                        mean_coefficient[i][j][k][l] = 0.0;
        dbl area_of_domain = 0.0;

        {
            dealii::QGauss<dim>  quadrature_formula(2);

            dealii::FEValues<dim> fe_values (dof_handler.get_fe(), quadrature_formula,
                    dealii::update_quadrature_points | dealii::update_JxW_values);

            cst n_q_points = quadrature_formula.size();


            auto cell = dof_handler.begin_active();
            auto endc = dof_handler.end();
            for (; cell != endc; ++cell)
            {
                fe_values .reinit (cell);

                for (st q_point = 0; q_point < n_q_points; ++q_point)
                    for (st i = 0; i < 3; ++i)
                        for (st j = 0; j < 3; ++j)
                            for (st k = 0; k < 3; ++k)
                                for (st l = 0; l < 3; ++l)
                                    mean_coefficient[i][j][k][l] += 
                                        coef[cell->material_id()][i][j][k][l] *
                                        fe_values.JxW(q_point);

                for (st q_point = 0; q_point < n_q_points; ++q_point)
                    area_of_domain += fe_values.JxW(q_point);
            };

            for (st i = 0; i < 3; ++i)
                for (st j = 0; j < 3; ++j)
                    for (st k = 0; k < 3; ++k)
                        for (st l = 0; l < 3; ++l)
                            mean_coefficient[i][j][k][l] /= area_of_domain;
        };

        st len_vector_solution = dof_handler.n_dofs() + heat_slae.solution[0].size();
        ATools::FourthOrderTensor mean_stress;
        ATools::FourthOrderTensor meta_coefficient;

        auto solution = [&elastic_slae, &heat_slae] (cst i, cst j, cst n) {
            enum {x, y, z};
            cu8 indx[3][3] = {
                {0, 3, x},
                {3, 1, y},
                {x, y, 2}
            };

            cst elastic_size = elastic_slae.solution[0].size();
            cst heat_size    = heat_slae.solution[0].size();

            if ((i != j) and ((i == z) or (j == z)))
            {
                cst n_var = n - elastic_size;
                return n < elastic_size ? 0.0 : heat_slae.solution[indx[i][j]](n_var);
            }
            else
            {
                return n < elastic_size ? elastic_slae.solution[indx[i][j]](n) : 0.0;
            };
        };

        auto stress = [&elastic_slae, &heat_slae] (cst i, cst j, cst n) {
            enum {x, y, z};
            cu8 indx[3][3] = {
                {0, 3, x},
                {3, 1, y},
                {x, y, 2}
            };

            cst elastic_size = elastic_slae.solution[0].size();
            cst heat_size    = heat_slae.solution[0].size();

            if ((i != j) and ((i == z) or (j == z)))
            {
                cst n_var = n - elastic_size;
                return n < elastic_size ? 0.0 : heat_slae.rhsv[indx[i][j]](n_var);
            }
            else
            {
                return n < elastic_size ? elastic_slae.rhsv[indx[i][j]](n) : 0.0;
            };
        };

    for (auto i : {x, y, z}) 
        for (auto j : {x, y, z}) 
            for (auto k : {x, y, z}) 
                for (auto l : {x, y, z}) 
            {
                mean_stress[i][j][k][l] = 0.0;

                for (st m = 0; m < len_vector_solution; ++m)
                    mean_stress[i][j][k][l] += solution(i, j, m) * (-stress(k, l, m));

                mean_stress[i][j][k][l] /= area_of_domain;

                // printf("%ld %ld %ld %ld %f %f %f\n", i, j, k, l,
                //         mean_stress[i][j][k][l], mean_coefficient[i][j][k][l],
                //         mean_coefficient[i][j][k][l] + mean_stress[i][j][k][l]);

                meta_coefficient[i][j][k][l] = 
                    mean_coefficient[i][j][k][l] + mean_stress[i][j][k][l];
            };
    
    // for (st m = 0; m < len_vector_solution; ++m)
    //     printf("%f %f %f\n", 
    //             solution(1, 2, m), -stress(1,2,m), solution(1, 2, m) * (-stress(1, 2, m)));

        return  meta_coefficient;
    };

    //! Расчет метакоэффициентов, обьёмная задача упругости
    template<u8 dim>
    ATools::FourthOrderTensor calculate_meta_coefficients_3d_elastic(
            const dealii::DoFHandler<dim> &dof_handler,
            const SystemsLinearAlgebraicEquations<6> &slae,
            const vec<ATools::FourthOrderTensor> &coef)
    {
        ATools::FourthOrderTensor mean_coefficient;
        for (st i = 0; i < 3; ++i)
            for (st j = 0; j < 3; ++j)
                for (st k = 0; k < 3; ++k)
                    for (st l = 0; l < 3; ++l)
                        mean_coefficient[i][j][k][l] = 0.0;
        dbl area_of_domain = 0.0;

        {
            dealii::QGauss<dim>  quadrature_formula(2);

            dealii::FEValues<dim> fe_values (dof_handler.get_fe(), quadrature_formula,
                    dealii::update_quadrature_points | dealii::update_JxW_values);

            cst n_q_points = quadrature_formula.size();

            auto cell = dof_handler.begin_active();
            auto endc = dof_handler.end();
            for (; cell != endc; ++cell)
            {
                fe_values .reinit (cell);

                for (st q_point = 0; q_point < n_q_points; ++q_point)
                    for (st i = 0; i < 3; ++i)
                        for (st j = 0; j < 3; ++j)
                            for (st k = 0; k < 3; ++k)
                                for (st l = 0; l < 3; ++l)
                                    mean_coefficient[i][j][k][l] += 
                                        coef[cell->material_id()][i][j][k][l] *
                                        fe_values.JxW(q_point);

                for (st q_point = 0; q_point < n_q_points; ++q_point)
                    area_of_domain += fe_values.JxW(q_point);
            };

            for (st i = 0; i < 3; ++i)
                for (st j = 0; j < 3; ++j)
                    for (st k = 0; k < 3; ++k)
                        for (st l = 0; l < 3; ++l)
                            mean_coefficient[i][j][k][l] /= area_of_domain;
        };

        st len_vector_solution = dof_handler.n_dofs();
        ATools::FourthOrderTensor mean_stress;
        ATools::FourthOrderTensor meta_coefficient;

        cu8 indx[3][3] = {
            {0, 3, 4},
            {3, 1, 5},
            {4, 5, 2}
        };

        for (auto i : {x, y, z}) 
            for (auto j : {x, y, z}) 
                for (auto k : {x, y, z}) 
                    for (auto l : {x, y, z}) 
                    {
                        mean_stress[i][j][k][l] = 0.0;

                        for (st m = 0; m < len_vector_solution; ++m)
                            mean_stress[i][j][k][l] += slae.solution[indx[i][j]](m) * (-slae.rhsv[indx[k][l]](m));

                        mean_stress[i][j][k][l] /= area_of_domain;

                        // printf("%ld %ld %ld %ld %f %f %f\n", i, j, k, l,
                        //         mean_stress[i][j][k][l], mean_coefficient[i][j][k][l],
                        //         mean_coefficient[i][j][k][l] + mean_stress[i][j][k][l]);

                        meta_coefficient[i][j][k][l] = 
                            mean_coefficient[i][j][k][l]
                            + 
                            mean_stress[i][j][k][l];
                    };

        // for (st m = 0; m < len_vector_solution; ++m)
        //     printf("%f %f %f\n", 
        //             solution(1, 2, m), -stress(1,2,m), solution(1, 2, m) * (-stress(1, 2, m)));

        return  meta_coefficient;
    };

    dbl calculate_meta_coefficients_3d_elastic_from_stress(
            const dealii::DoFHandler<3> &dof_handler,
            const dealii::Vector<dbl> &stress,
            cst component)
    {
        dbl mean_coefficient = 0.0;
        dbl area_of_domain = 0.0;

        {
            dealii::QGauss<3>  quadrature_formula(2);

            dealii::FEValues<3> fe_values (dof_handler.get_fe(), quadrature_formula,
                    dealii::update_values | dealii::update_quadrature_points |
                    dealii::update_JxW_values);

            cst n_q_points = quadrature_formula.size();

            auto cell = dof_handler.begin_active();
            auto endc = dof_handler.end();
            for (; cell != endc; ++cell)
            {
                fe_values .reinit (cell);

                for (st q_point = 0; q_point < n_q_points; ++q_point)
                    for (st i = 0; i < 8; ++i)
                        mean_coefficient += 
                            stress[cell->vertex_dof_index(i, component)] *
                            fe_values.shape_value (i, q_point) *
                            fe_values.JxW(q_point);

                for (st q_point = 0; q_point < n_q_points; ++q_point)
                    area_of_domain += fe_values.JxW(q_point);
            };

            mean_coefficient /= area_of_domain;
        };
        return  mean_coefficient;
    };

    //! Расчет метакоэффициентов, обьёмная задача упругости, произвольное приближение
    class MetaCoefficientElasticCalculator 
    {
        public:
            MetaCoefficientElasticCalculator (const dealii::FiniteElement<3> &fe);

            MetaCoefficientElasticCalculator(
                    const dealii::DoFHandler<3> &dof_handler,
                    const vec<ATools::FourthOrderTensor> &coefficient,
                    const dealii::FiniteElement<3> &fe);

            arr<dbl, 3> calculate(
                    const arr<i32, 3> approximation,
                    cst nu,
                    const dealii::DoFHandler<3> &dof_handler,
                    const OnCell::ArrayWithAccessToVector<arr<dealii::Vector<dbl>, 3>> &cell_func);

            static cst x = 0;
            static cst y = 1;
            static cst z = 2;

            vec<ATools::FourthOrderTensor> coef; //!< Коэффициенты материла, например теплопроводность.
            dealii::QGauss<3>          quadrature_formula; //!< Формула интегрирования в квадратурах.
            dealii::FEValues<3, 3>     fe_values; //!< Тип функций формы.
            vec<u32>                   global_dof_indices; //!< Список глобальных значений индексов с ячейки.
            cu8                        dofs_per_cell; //!< Количество узлов в ячейке (зависит от типа функций формы).
            cu8                        num_quad_points; //!< Количество точек по которым считается квадратура.
            dbl                        area_of_domain = 0.0;
    };

        MetaCoefficientElasticCalculator::MetaCoefficientElasticCalculator (
                const dealii::FiniteElement<3> &fe) :
            quadrature_formula (2),
            fe_values (fe, quadrature_formula,
                    dealii::update_values | dealii::update_gradients | dealii::update_quadrature_points | 
                    dealii::update_JxW_values),
            dofs_per_cell (fe.dofs_per_cell),
            global_dof_indices (fe.dofs_per_cell),
            num_quad_points (quadrature_formula.size())
    {};

    MetaCoefficientElasticCalculator::MetaCoefficientElasticCalculator(
            const dealii::DoFHandler<3> &dof_handler,
            const vec<ATools::FourthOrderTensor> &coefficient,
            const dealii::FiniteElement<3> &fe) :
        MetaCoefficientElasticCalculator(fe)
    {
        auto cell = dof_handler.begin_active();
        auto endc = dof_handler.end();
        for (; cell != endc; ++cell)
        {
            fe_values .reinit (cell);

            for (st q_point = 0; q_point < num_quad_points; ++q_point)
                area_of_domain += fe_values.JxW(q_point);
        };

        coef .resize (coefficient.size());
        for (st n = 0; n < coefficient.size(); ++n)
            for (st i = 0; i < 3; ++i)
                for (st j = 0; j < 3; ++j)
                    for (st k = 0; k < 3; ++k)
                        for (st l = 0; l < 3; ++l)
                            coef[n][i][j][k][l] = coefficient[n][i][j][k][l];
    };

    arr<dbl, 3> MetaCoefficientElasticCalculator::calculate(
            const arr<i32, 3> approximation,
            cst nu,
            const dealii::DoFHandler<3> &dof_handler,
            const OnCell::ArrayWithAccessToVector<arr<dealii::Vector<dbl>, 3>> &cell_func)
    {
        arr<dbl, 3> res = {0.0, 0.0, 0.0};
        arr<i32, 3> k = {approximation[0], approximation[1], approximation[2]};

        auto cell = dof_handler.begin_active();
        auto endc = dof_handler.end();
        for (; cell != endc; ++cell)
        {
            fe_values .reinit (cell);
            cst material_id = cell->material_id();
            cell->get_dof_indices (global_dof_indices);

            for (st alpha = 0; alpha < 3; ++alpha)
            {
                for (st q_point = 0; q_point < num_quad_points; ++q_point)
                {
                    for (st beta = 0; beta < 3; ++beta)
                    // for (st phi = 0; phi < 3; ++phi)
                    // cst beta = nu;
                    {
                        for (st m = 0; m < dofs_per_cell; ++m)
                        {
                            // cst beta = m % 3;
                            cst phi = m % 3;
                            for (st psi = 0; psi < 3; ++psi)
                            {
                                // arr<i32, 3> k_psi = {k[x], k[y], k[z]}; // k - э_psi
                                // bool cell_k_psi_exist = true;
                                // if (k[psi] == 0)
                                //     cell_k_psi_exist = false;
                                // else
                                //     k_psi[psi]--;

                                arr<i32, 3> k_beta = {k[x], k[y], k[z]}; // k - э_beta
                                bool cell_k_beta_exist = true;
                                if (k[beta] == 0)
                                    cell_k_beta_exist = false;
                                else
                                    k_beta[beta]--;

                                arr<i32, 3> k_psi_beta = {k_beta[x], k_beta[y], k_beta[z]}; // k - э_psi - э_beta
                                bool cell_k_psi_beta_exist = true;
                                if (cell_k_beta_exist)
                                {
                                    if (k_beta[psi] == 0)
                                        cell_k_psi_beta_exist = false;
                                    else
                                        k_psi_beta[psi]--;
                                }
                                else
                                    cell_k_psi_beta_exist = false;

                                dbl cell_k_psi      = 0.0;
                                dbl cell_k_beta     = 0.0;
                                dbl cell_k_psi_beta = 0.0;

                                // if (cell_k_psi_exist)
                                //     cell_k_psi = cell_func[k_psi][nu][global_dof_indices[m]];
                                // else
                                //     cell_k_psi = 0.0;
                                //
                                // printf("%d %d %d  %d %d %d  %d %d  %d %d\n", 
                                //         k_beta[0], k_beta[1], k_beta[2],
                                //         k_psi_beta[0], k_psi_beta[1], k_psi_beta[2],
                                //         cell_k_beta_exist, cell_k_psi_beta_exist,
                                //         beta, psi
                                //         );
                                if (cell_k_beta_exist)
                                    cell_k_beta = cell_func[k_beta][nu][global_dof_indices[m]];
                                else
                                    cell_k_beta = 0.0;

                                if (cell_k_psi_beta_exist)
                                    cell_k_psi_beta = cell_func[k_psi_beta][nu][global_dof_indices[m]];
                                else
                                    cell_k_psi_beta = 0.0;
                                // // printf("%f\n", cell_k_psi_beta);

                                res[alpha] += 
                                    (
                                    this->coef[material_id][alpha][beta][phi][psi] *
                                     (
                                      // 0.0 *
                                      cell_k_beta *
                                      this->fe_values.shape_grad (m, q_point)[psi] +
                                      // 0.0 *
                                      cell_k_psi_beta *
                                      this->fe_values.shape_value (m, q_point)
                                     )
                                    ) *
                                    this->fe_values.JxW(q_point);
                            };
                        };
                    };
                };
            };
        };

        for (st i = 0; i < 3; ++i)
            res[i] /= area_of_domain;

        return  res;
    };

    //! Расчет напряжений, обьёмная задача упругости, произвольное приближение
    class StressCalculator 
    {
        public:
            StressCalculator (const dealii::FiniteElement<3> &fe);

            StressCalculator(
                    const dealii::DoFHandler<3> &dof_handler,
                    const vec<ATools::FourthOrderTensor> &coefficient,
                    const dealii::FiniteElement<3> &fe);

            void calculate(
                    const arr<i32, 3> approximation,
                    cst nu,
                    cst alpha,
                    cst beta,
                    const dealii::DoFHandler<3> &dof_handler,
                    const OnCell::ArrayWithAccessToVector<arr<dealii::Vector<dbl>, 3>> &cell_func,
                    dealii::Vector<dbl> &stress);

            static cst x = 0;
            static cst y = 1;
            static cst z = 2;

            vec<ATools::FourthOrderTensor> coef; //!< Коэффициенты материла, например теплопроводность.
            dealii::QGauss<3>          quadrature_formula; //!< Формула интегрирования в квадратурах.
            dealii::FEValues<3, 3>     fe_values; //!< Тип функций формы.
            vec<u32>                   global_dof_indices; //!< Список глобальных значений индексов с ячейки.
            cu8                        dofs_per_cell; //!< Количество узлов в ячейке (зависит от типа функций формы).
            cu8                        num_quad_points; //!< Количество точек по которым считается квадратура.
            dbl                        area_of_domain = 0.0;
    };

        StressCalculator::StressCalculator (
                const dealii::FiniteElement<3> &fe) :
            quadrature_formula (2),
            fe_values (fe, quadrature_formula,
                    dealii::update_values | dealii::update_gradients | dealii::update_quadrature_points | 
                    dealii::update_JxW_values),
            dofs_per_cell (fe.dofs_per_cell),
            global_dof_indices (fe.dofs_per_cell),
            num_quad_points (quadrature_formula.size())
    {};

    StressCalculator::StressCalculator(
            const dealii::DoFHandler<3> &dof_handler,
            const vec<ATools::FourthOrderTensor> &coefficient,
            const dealii::FiniteElement<3> &fe) :
        StressCalculator(fe)
    {
        auto cell = dof_handler.begin_active();
        auto endc = dof_handler.end();
        for (; cell != endc; ++cell)
        {
            fe_values .reinit (cell);

            for (st q_point = 0; q_point < num_quad_points; ++q_point)
                area_of_domain += fe_values.JxW(q_point);
        };

        coef .resize (coefficient.size());
        for (st n = 0; n < coefficient.size(); ++n)
            for (st i = 0; i < 3; ++i)
                for (st j = 0; j < 3; ++j)
                    for (st k = 0; k < 3; ++k)
                        for (st l = 0; l < 3; ++l)
                            coef[n][i][j][k][l] = coefficient[n][i][j][k][l];
    };

    void StressCalculator::calculate(
            const arr<i32, 3> approximation,
            cst nu,
            cst alpha,
            cst beta,
            const dealii::DoFHandler<3> &dof_handler,
            const OnCell::ArrayWithAccessToVector<arr<dealii::Vector<dbl>, 3>> &cell_func,
            dealii::Vector<dbl> &stress)
    {
        arr<i32, 3> k = {approximation[0], approximation[1], approximation[2]};

        vec<st> N(dof_handler.n_dofs());

        auto cell = dof_handler.begin_active();
        auto endc = dof_handler.end();
        for (; cell != endc; ++cell)
        {
            fe_values .reinit (cell);
            cst material_id = cell->material_id();
            cell->get_dof_indices (global_dof_indices);

            dbl cell_area = 0.0;
            for (st q_point = 0; q_point < num_quad_points; ++q_point)
            {
                cell_area += fe_values.JxW(q_point);
            };

            dbl tmp = 0.0;

            for (st m = 0; m < dofs_per_cell; ++m)
            {
                for (st q_point = 0; q_point < num_quad_points; ++q_point)
                {
                    cst phi = m % 3;
                    for (st psi = 0; psi < 3; ++psi)
                    {
                        arr<i32, 3> k_psi = {k[x], k[y], k[z]}; // k - э_beta
                        bool cell_k_psi_exist = true;
                        if (k[psi] == 0)
                            cell_k_psi_exist = false;
                        else
                            k_psi[psi]--;

                        cdbl cell_k = cell_func[k][nu][global_dof_indices[m]];

                        dbl cell_k_psi = 0.0;
                        if (cell_k_psi_exist)
                            cell_k_psi = cell_func[k_psi][nu][global_dof_indices[m]];
                        else
                            cell_k_psi = 0.0;

                        // tmp[m] +=
                        tmp +=
                            (
                             this->coef[material_id][alpha][beta][phi][psi] *
                             (
                              // 0.0 *
                              cell_k *
                              this->fe_values.shape_grad (m, q_point)[psi] +
                              // 0.0 *
                              cell_k_psi *
                              this->fe_values.shape_value (m, q_point)
                             )
                            ) *
                            this->fe_values.JxW(q_point);
                    };
                };
                ++(N[global_dof_indices[m]]); 
            };
            // stress[0] = tmp[0];
            // global_dof_indices[0] = 1;
            // stress[global_dof_indices[0]] = tmp[0];
            // puts("f1");
            // stress[global_dof_indices[0]] = tmp[2];
            // puts("f2");
            // stress[global_dof_indices[0]] = tmp[4];
            // puts("f3");
            // stress[global_dof_indices[0]] = tmp[6];
            // puts("f4");
            // stress[global_dof_indices[0]] = 1.0;
            // printf("tmp %f %d\n", tmp[0], dofs_per_cell);
            // tmp[0] = 1.0;
            // arr<dbl, 3> component =
            // arr<dbl, 3>{
            //     tmp[0] + tmp[3] + tmp[6] + tmp[9] + tmp[12] + tmp[15] + tmp[18] + tmp[21],
            //     tmp[1] + tmp[4] + tmp[7] + tmp[10] + tmp[14] + tmp[16] + tmp[19] + tmp[22],
            //     tmp[2] + tmp[5] + tmp[8] + tmp[11] + tmp[14] + tmp[17] + tmp[20] + tmp[23]
            // };
            // arr<dbl, 3> component;
            // component[0] = (tmp[0] + tmp[3] + tmp[6] + tmp[9] + tmp[12] + tmp[15] + tmp[18] + tmp[21]) / (cell_area);
            // component[1] = (tmp[1] + tmp[4] + tmp[7] + tmp[10] + tmp[14] + tmp[16] + tmp[19] + tmp[22]) / (cell_area);
            // component[2] = (tmp[2] + tmp[5] + tmp[8] + tmp[11] + tmp[14] + tmp[17] + tmp[20] + tmp[23]) / (cell_area);
            // printf("%f %f %f %f %f %f %f %f %f\n",
            //          tmp[0] , tmp[3] , tmp[6] , tmp[9] , tmp[12] , tmp[15] , tmp[18] , tmp[21],
            //          component[0]);
            // printf("comp %f %f %f\n", component[0], component[1], component[2]);
            for (st i = 0; i < 8; ++i)
            {
                // stress[global_dof_indices[i*3]] += component[0];
                // stress[global_dof_indices[i*3+1]] += component[1];
                // stress[global_dof_indices[i*3+2]] += component[2];
                stress[global_dof_indices[i*3+beta]] += (tmp / cell_area);// / N[i];
            };

            // stress[global_dof_indices[0]] += (tmp[0] + tmp[2] + tmp[4] + tmp[6])/ (cell_area);
            // stress[global_dof_indices[2]] += (tmp[0] + tmp[2] + tmp[4] + tmp[6]) / (cell_area);
            // stress[global_dof_indices[4]] += (tmp[0] + tmp[2] + tmp[4] + tmp[6]) / (cell_area);
            // stress[global_dof_indices[6]] += (tmp[0] + tmp[2] + tmp[4] + tmp[6]) / (cell_area);
            //
            // stress[global_dof_indices[1]] += (tmp[1] + tmp[3] + tmp[5] + tmp[7]) / (cell_area);
            // stress[global_dof_indices[3]] += (tmp[1] + tmp[3] + tmp[5] + tmp[7]) / (cell_area);
            // stress[global_dof_indices[5]] += (tmp[1] + tmp[3] + tmp[5] + tmp[7]) / (cell_area);
            // stress[global_dof_indices[7]] += (tmp[1] + tmp[3] + tmp[5] + tmp[7]) / (cell_area);
        };
        // puts("333");

        for (st i = 0; i < stress.size(); ++i)
        {
            if ((i%3) == beta)
            stress[i] /= N[i];
        };

    };








    //! Расчет деформаций, обьёмная задача упругости, произвольное приближение
    class DeformCalculator 
    {
        public:
            DeformCalculator (const dealii::FiniteElement<3> &fe);

            DeformCalculator(
                    const dealii::DoFHandler<3> &dof_handler,
                    const dealii::FiniteElement<3> &fe);

            void calculate(
                    const arr<i32, 3> approximation,
                    cst nu,
                    cst beta,
                    const dealii::DoFHandler<3> &dof_handler,
                    const OnCell::ArrayWithAccessToVector<arr<dealii::Vector<dbl>, 3>> &cell_func,
                    dealii::Vector<dbl> &deform);

            static cst x = 0;
            static cst y = 1;
            static cst z = 2;

            dealii::QGauss<3>          quadrature_formula; //!< Формула интегрирования в квадратурах.
            dealii::FEValues<3, 3>     fe_values; //!< Тип функций формы.
            vec<u32>                   global_dof_indices; //!< Список глобальных значений индексов с ячейки.
            cu8                        dofs_per_cell; //!< Количество узлов в ячейке (зависит от типа функций формы).
            cu8                        num_quad_points; //!< Количество точек по которым считается квадратура.
            dbl                        area_of_domain = 0.0;
    };

        DeformCalculator::DeformCalculator (
                const dealii::FiniteElement<3> &fe) :
            quadrature_formula (2),
            fe_values (fe, quadrature_formula,
                    dealii::update_values | dealii::update_gradients | dealii::update_quadrature_points | 
                    dealii::update_JxW_values),
            dofs_per_cell (fe.dofs_per_cell),
            global_dof_indices (fe.dofs_per_cell),
            num_quad_points (quadrature_formula.size())
    {};

    DeformCalculator::DeformCalculator(
            const dealii::DoFHandler<3> &dof_handler,
            const dealii::FiniteElement<3> &fe) :
        DeformCalculator(fe)
    {
        auto cell = dof_handler.begin_active();
        auto endc = dof_handler.end();
        for (; cell != endc; ++cell)
        {
            fe_values .reinit (cell);

            for (st q_point = 0; q_point < num_quad_points; ++q_point)
                area_of_domain += fe_values.JxW(q_point);
        };
    };

    void DeformCalculator::calculate(
            const arr<i32, 3> approximation,
            cst nu,
            cst beta,
            const dealii::DoFHandler<3> &dof_handler,
            const OnCell::ArrayWithAccessToVector<arr<dealii::Vector<dbl>, 3>> &cell_func,
            dealii::Vector<dbl> &deform)
    {
        arr<i32, 3> k = {approximation[0], approximation[1], approximation[2]};

        vec<st> N(dof_handler.n_dofs());

        auto cell = dof_handler.begin_active();
        auto endc = dof_handler.end();
        // puts("2222");
        for (; cell != endc; ++cell)
        {
            fe_values .reinit (cell);
            cst material_id = cell->material_id();
            cell->get_dof_indices (global_dof_indices);

            dbl cell_area = 0.0;
            for (st q_point = 0; q_point < num_quad_points; ++q_point)
            {
                cell_area += fe_values.JxW(q_point);
            };

            arr<dbl,3> tmp = {0.0, 0.0, 0.0};

            for (st m = 0; m < dofs_per_cell; ++m)
            {
                cst alpha = m % 3;
                for (st q_point = 0; q_point < num_quad_points; ++q_point)
                {
                        arr<i32, 3> k_beta = {k[x], k[y], k[z]}; // k - э_beta
                        bool cell_k_beta_exist = true;
                        if (k[beta] == 0)
                            cell_k_beta_exist = false;
                        else
                            k_beta[beta]--;

                        cdbl cell_k = cell_func[k][nu][global_dof_indices[m]];

                        dbl cell_k_beta = 0.0;
                        if (cell_k_beta_exist)
                            cell_k_beta = cell_func[k_beta][nu][global_dof_indices[m]];
                        else
                            cell_k_beta = 0.0;

                        tmp[alpha] +=
                            (
                              cell_k *
                              this->fe_values.shape_grad (m, q_point)[beta] +
                              cell_k_beta *
                              this->fe_values.shape_value (m, q_point)
                            ) *
                            this->fe_values.JxW(q_point);
                };
                ++(N[global_dof_indices[m]]); 
            };
            for (st i = 0; i < 8; ++i)
            {
                for (st alpha = 0; alpha < 3; ++alpha)
                {
                    deform[global_dof_indices[i*3+alpha]] += (tmp[alpha] / cell_area);
                };
            };
        };

        for (st i = 0; i < deform.size(); ++i)
        {
            if ((i%3) == beta)
            deform[i] /= N[i];
        };

    };
};

#endif
