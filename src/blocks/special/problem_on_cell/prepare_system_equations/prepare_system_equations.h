#ifndef TRIVIAL_PREPARE_SYSTEM_EQUATIONS_ON_CELL
#define TRIVIAL_PREPARE_SYSTEM_EQUATIONS_ON_CELL

#include <deal.II/lac/vector.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>

#include "../../../general/domain/domain.h"
#include "../../../general/boundary_value/boundary_value.h"
#include "../system_linear_algebraic_equations/system_linear_algebraic_equations.h"
#include "../domain_looper/domain_looper.h"
// #include "../../../../../../deal/main/domain_looper/sources/domain_looper.h"
// #include "../black_on_white_substituter/black_on_white_substituter.h"

//! Задача на ячейке
/*!
 * Сюда входят дополнительные средства необходимые для решения задачи на
 * определение ячейковых функций. 
 */
namespace OnCell
{
    //! Формирование СЛАУ (se) по расчетной области (domain), в случае задачи на ячейке
    template<bool type_space, u8 dim, u8 num_tasks>
    void prepare_system_equations (
            ::OnCell::SystemsLinearAlgebraicEquations<num_tasks> &se,
            ::OnCell::BlackOnWhiteSubstituter &bows,
            const Domain<dim> &domain)
    {
        dealii::SparsityPattern sp;
        {
            dealii::DynamicSparsityPattern dsp (domain.dof_handler.n_dofs());
            dealii::DoFTools::make_sparsity_pattern (domain.dof_handler, dsp);

            // {
            // std::ofstream output ("csp1.gpd");
            // c_sparsity .print_gnuplot (output);
            // };

            ::OnCell::DomainLooper<dim, type_space> dl;
            // DomainLooper<dim, 0> dl;
            dl .loop_domain(
                    domain.dof_handler,
                    bows,
                    dsp);

            sp .copy_from (dsp);
        };

        // {
        // std::ofstream output ("csp_old.gpd");
        // c_sparsity .print_gnuplot (output);
        // };

        sp .compress ();

        se .matrix_reinit (sp);

        for (st i = 0; i < num_tasks; ++i)
        {
            se.solution[i] .reinit (domain.dof_handler .n_dofs());
            se.rhsv[i]     .reinit (domain.dof_handler .n_dofs());
        };
    };
};

#endif
