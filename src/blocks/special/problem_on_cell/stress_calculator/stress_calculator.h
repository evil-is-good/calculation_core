#ifndef STRESS_CALCULATOR_BD_ON_CELL
#define STRESS_CALCULATOR_BD_ON_CELL 1

#include "../../../../../../prmt_sintactic_addition/prmt_sintactic_addition.h"
#include "../../../general/additional_tools/types/types.h"
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/lac/vector.h>

namespace OnCell
{
    //! Расчет напряжений, обьёмная задача упругости, произвольное приближение
    class StressCalculatorBD 
    {
        public:
            StressCalculatorBD (
                    const dealii::DoFHandler<3> &dof_handler_p,
                    const dealii::FiniteElement<3> &fe);

            StressCalculatorBD(
                    const dealii::DoFHandler<3> &dof_handler_p,
                    const vec<ATools::FourthOrderTensor> &coefficient,
                    const dealii::FiniteElement<3> &fe);

            void calculate(
                    const arr<i32, 3> approximation,
                    cst nu,
                    cst alpha,
                    cst beta,
                    const OnCell::ArrayWithAccessToVector<arr<dealii::Vector<dbl>, 3>> &move,
                    dealii::Vector<dbl> &stress);

            static cst x = 0;
            static cst y = 1;
            static cst z = 2;

            const dealii::DoFHandler<3>& dof_handler;
            vec<ATools::FourthOrderTensor> coef; //!< Коэффициенты материла, например теплопроводность.
            const dealii::QGauss<3>          quadrature_formula; //!< Формула интегрирования в квадратурах.
            dealii::FEValues<3, 3>     fe_values; //!< Тип функций формы.
            vec<u32>                   global_dof_indices; //!< Список глобальных значений индексов с ячейки.
            cu8                        dofs_per_cell; //!< Количество узлов в ячейке (зависит от типа функций формы).
            cu8                        num_quad_points; //!< Количество точек по которым считается квадратура.
            dbl                        area_of_domain = 0.0;
    };

        StressCalculatorBD::StressCalculatorBD (
                const dealii::DoFHandler<3> &dof_handler_p,
                const dealii::FiniteElement<3> &fe) :
            dof_handler (dof_handler_p),
            quadrature_formula (2),
            fe_values (fe, quadrature_formula,
                    dealii::update_values | dealii::update_gradients | dealii::update_quadrature_points | 
                    dealii::update_JxW_values),
            dofs_per_cell (fe.dofs_per_cell),
            global_dof_indices (fe.dofs_per_cell),
            num_quad_points (quadrature_formula.size())
    {};

    StressCalculatorBD::StressCalculatorBD(
            const dealii::DoFHandler<3> &dof_handler_p,
            const vec<ATools::FourthOrderTensor> &coefficient,
            const dealii::FiniteElement<3> &fe) :
        StressCalculatorBD(dof_handler_p, fe)
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

    void StressCalculatorBD::calculate(
            const arr<i32, 3> approximation,
            cst nu,
            cst alpha,
            cst beta,
            const OnCell::ArrayWithAccessToVector<arr<dealii::Vector<dbl>, 3>> &move,
            dealii::Vector<dbl> &stress)
    {
        arr<i32, 3> k = {approximation[0], approximation[1], approximation[2]};

        vec<st> N(stress.size());

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

                        cdbl cell_k = move[k][nu][global_dof_indices[m]];

                        dbl cell_k_psi = 0.0;
                        if (cell_k_psi_exist)
                            cell_k_psi = move[k_psi][nu][global_dof_indices[m]];
                        else
                            cell_k_psi = 0.0;

                        tmp +=
                            (
                             this->coef[material_id][alpha][beta][phi][psi] *
                             (
                              cell_k *
                              this->fe_values.shape_grad (m, q_point)[psi] +
                              cell_k_psi *
                              this->fe_values.shape_value (m, q_point)
                             )
                            ) *
                            this->fe_values.JxW(q_point);
                    };
                };
                ++(N[global_dof_indices[m]]); 
            };

            for (st i = 0; i < 8; ++i)
            {
                // stress[global_dof_indices[i*3] / 3] += (tmp / cell_area);// / N[i];
                // ++(N[global_dof_indices[i*3] / 3]); 
                stress[global_dof_indices[i*3+beta]] += (tmp / cell_area);// / N[i];
            };
        };

        for (st i = 0; i < stress.size(); ++i)
        {
            if ((i%3) == beta)
                stress[i] /= N[i];
        };

    };
}
#endif
