#ifndef SYSTEM_LINEAR_ALGEBRAIC_EQUATION_ON_CELL

#define SYSTEM_LINEAR_ALGEBRAIC_EQUATION_ON_CELL

#include "../../../../../../prmt_sintactic_addition/prmt_sintactic_addition.h"
#include "../../../general/system_linear_algebraic_equations/system_linear_algebraic_equations.h"
#include <deal.II/lac/vector.h>
#include <deal.II/lac/sparse_matrix.h>

namespace OnCell
{
//! Система линейных алгебраических уравнений на ячейке
/*!
 * Отличается от обыкновенной системы наличием num_rhsv векторов решения и правой части.
 * Потому как в задаче на ячейке нужно решать несколько систем, отличающихся
 * друг от друга тольуо правыми частями.
 */
template <u8 num_rhsv>
struct  SystemsLinearAlgebraicEquations : public ::SystemsLinearAlgebraicEquations
{
    arr<dealii::Vector<dbl>, num_rhsv> solution; //!< вектора решения
    arr<dealii::Vector<dbl>, num_rhsv> rhsv; //!< вектора правой части
};
};

#endif
