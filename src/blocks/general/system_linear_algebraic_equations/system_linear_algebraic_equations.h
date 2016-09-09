#ifndef SYSTEM_LINEAR_ALGEBRAIC_EQUATION

#define SYSTEM_LINEAR_ALGEBRAIC_EQUATION

#include "../../../../../prmt_sintactic_addition/prmt_sintactic_addition.h"
#include <deal.II/lac/vector.h>
#include <deal.II/lac/sparse_matrix.h>

//! Система линейных алгебраических уравнений (ваш кэп)
struct  SystemsLinearAlgebraicEquations 
{
    //! Формирование матрици из паттерна разреженности
    void matrix_reinit (dealii::SparsityPattern &csp)
    {
        this->sparsity_pattern .copy_from (csp);
        this->matrix .reinit (this->sparsity_pattern);
    };

    dealii::SparsityPattern   sparsity_pattern; //!< паттерн разреженности
    dealii::SparseMatrix<dbl> matrix; //!< матрица
    dealii::Vector<dbl>       solution; //!< вектор решения
    dealii::Vector<dbl>       rhsv; //!< вектор правой части
};

#endif
