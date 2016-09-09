#ifndef ALTERNATE_PREPARE_SYSTEM_EQUATIONS_ON_CELL
#define ALTERNATE_PREPARE_SYSTEM_EQUATIONS_ON_CELL 1

#include <deal.II/lac/vector.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>

#include "../../../general/domain/domain.h"
#include "../../../general/boundary_value/boundary_value.h"
#include "../system_linear_algebraic_equations/system_linear_algebraic_equations.h"
#include "../domain_looper_alternate/domain_looper_alternate.h"

//! Задача на ячейке
/*!
 * Сюда входят дополнительные средства необходимые для решения задачи на
 * определение ячейковых функций. 
 */
namespace OnCell
{
    //! Формирование СЛАУ (se) по расчетной области (domain), в случае задачи на ячейке новый способ
    template<u8 dim, u8 spec_dim, u8 num_tasks>
    void prepare_system_equations_alternate (
            ::OnCell::SystemsLinearAlgebraicEquations<num_tasks> &se,
            ::OnCell::BlackOnWhiteSubstituter &bows,
            const Domain<dim> &domain)
    {
        dealii::SparsityPattern sp;
        {
            dealii::DynamicSparsityPattern dsp (domain.dof_handler.n_dofs());
            dealii::DoFTools::make_sparsity_pattern (domain.dof_handler, dsp);

            loop_domain_alternate<dim, spec_dim>(
                    domain.dof_handler,
                    bows,
                    dsp);

            sp .copy_from (dsp);
        };

        {
            std::ofstream output ("csp_alternate.gpd");
            sp .print_gnuplot (output);
        };

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
