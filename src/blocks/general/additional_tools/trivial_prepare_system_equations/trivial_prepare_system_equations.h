#ifndef TRIVIAL_PREPARE_SYSTEM_EQUATIONS
#define TRIVIAL_PREPARE_SYSTEM_EQUATIONS

#include <deal.II/lac/vector.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>

#include "../../domain/domain.h"
#include "../../system_linear_algebraic_equations/system_linear_algebraic_equations.h"
#include "../../boundary_value/boundary_value.h"

//! Дополнительный инструментарий
/*!
 * Additional tools \n
 * Сюда входят инструменты (функции) не вошедшие в другие, специализированные
 * пространства.
 */
namespace ATools
{
    //! Формирование СЛАУ (se) по расчетной области (domain), тривиальный случай
    template<u8 dim>
    void trivial_prepare_system_equations (
            SystemsLinearAlgebraicEquations &se,
            const Domain<dim> &domain)
    {
        dealii::SparsityPattern sp;
        {
            dealii::DynamicSparsityPattern dsp (domain.dof_handler.n_dofs());
            dealii::DoFTools::make_sparsity_pattern (domain.dof_handler, dsp);

            sp .copy_from (dsp);
        };
        se .matrix_reinit (sp);

        se.solution .reinit (domain.dof_handler .n_dofs());
        se.rhsv     .reinit (domain.dof_handler .n_dofs());
    };
};

#endif
