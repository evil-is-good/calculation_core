#ifndef DEFORM_CALCULATOR_BD_ON_CELL
#define DEFORM_CALCULATOR_BD_ON_CELL 1

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
    class DeformCalculatorBD 
    {
        public:
            DeformCalculatorBD (
                    const dealii::DoFHandler<3> &dof_handler_p,
                    const dealii::FiniteElement<3> &fe);

            // DeformCalculatorBD(
            //         const dealii::DoFHandler<3> &dof_handler_p,
            //         const dealii::FiniteElement<3> &fe);

            void calculate(
                    const arr<i32, 3> approximation,
                    cst nu,
                    cst alpha,
                    cst beta,
                    const OnCell::ArrayWithAccessToVector<arr<dealii::Vector<dbl>, 3>> &move,
                    dealii::Vector<dbl> &deform);

            static cst x = 0;
            static cst y = 1;
            static cst z = 2;

            const dealii::DoFHandler<3>& dof_handler;
            const dealii::QGauss<3>          quadrature_formula; //!< Формула интегрирования в квадратурах.
            dealii::FEValues<3, 3>     fe_values; //!< Тип функций формы.
            vec<u32>                   global_dof_indices; //!< Список глобальных значений индексов с ячейки.
            cu8                        dofs_per_cell; //!< Количество узлов в ячейке (зависит от типа функций формы).
            cu8                        num_quad_points; //!< Количество точек по которым считается квадратура.
            dbl                        area_of_domain = 0.0;
    };

    //     DeformCalculatorBD::DeformCalculatorBD (
    //             const dealii::DoFHandler<3> &dof_handler_p,
    //             const dealii::FiniteElement<3> &fe) :
    //         dof_handler (dof_handler_p),
    //         quadrature_formula (2),
    //         fe_values (fe, quadrature_formula,
    //                 dealii::update_values | dealii::update_gradients | dealii::update_quadrature_points | 
    //                 dealii::update_JxW_values),
    //         dofs_per_cell (fe.dofs_per_cell),
    //         global_dof_indices (fe.dofs_per_cell),
    //         num_quad_points (quadrature_formula.size())
    // {};

    DeformCalculatorBD::DeformCalculatorBD(
            const dealii::DoFHandler<3> &dof_handler_p,
            const dealii::FiniteElement<3> &fe) :
        // DeformCalculatorBD(dof_handler_p, fe)
            dof_handler (dof_handler_p),
            quadrature_formula (2),
            fe_values (fe, quadrature_formula,
                    dealii::update_values | dealii::update_gradients | dealii::update_quadrature_points | 
                    dealii::update_JxW_values),
            dofs_per_cell (fe.dofs_per_cell),
            global_dof_indices (fe.dofs_per_cell),
            num_quad_points (quadrature_formula.size())
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

    void DeformCalculatorBD::calculate(
            const arr<i32, 3> approximation,
            cst nu,
            cst alpha,
            cst beta,
            const OnCell::ArrayWithAccessToVector<arr<dealii::Vector<dbl>, 3>> &move,
            dealii::Vector<dbl> &deform)
    {
        arr<i32, 3> k = {approximation[0], approximation[1], approximation[2]};

        vec<st> N(deform.size());

        auto cell = dof_handler.begin_active();
        auto endc = dof_handler.end();
        for (; cell != endc; ++cell)
        {
            fe_values .reinit (cell);
            cell->get_dof_indices (global_dof_indices);

            dbl cell_area = 0.0;
            for (st q_point = 0; q_point < num_quad_points; ++q_point)
            {
                cell_area += fe_values.JxW(q_point);
            };

            dbl tmp = 0.0;

            for (st m = 0; m < dofs_per_cell; ++m)
            {
                if (alpha == (m % 3))
                {
                    for (st q_point = 0; q_point < num_quad_points; ++q_point)
                    {
                        arr<i32, 3> k_beta = {k[x], k[y], k[z]}; // k - э_beta
                        bool cell_k_beta_exist = true;
                        if (k[beta] == 0)
                            cell_k_beta_exist = false;
                        else
                            k_beta[beta]--;

                        cdbl cell_k = move[k][nu][global_dof_indices[m]];

                        dbl cell_k_beta = 0.0;
                        if (cell_k_beta_exist)
                            cell_k_beta = move[k_beta][nu][global_dof_indices[m]];
                        else
                            cell_k_beta = 0.0;

                        tmp +=
                            (
                              cell_k *
                              this->fe_values.shape_grad (m, q_point)[beta] +
                              cell_k_beta *
                              this->fe_values.shape_value (m, q_point)
                            ) *
                            this->fe_values.JxW(q_point);
                    };
                };
                ++(N[global_dof_indices[m]]); 
            };

            for (st i = 0; i < 8; ++i)
            {
                // deform[global_dof_indices[i*3] / 3] += (tmp / cell_area);// / N[i];
                // ++(N[global_dof_indices[i*3] / 3]); 
                deform[global_dof_indices[i*3+beta]] += (tmp / cell_area);// / N[i];
            };
        };

        for (st i = 0; i < deform.size(); ++i)
        {
            if ((i%3) == beta)
                deform[i] /= N[i];
        };

    };
}
#endif
