#ifndef ASSEMBLER_MATRIX_AND_VECTOR
#define ASSEMBLER_MATRIX_AND_VECTOR

#include <deal.II/lac/full_matrix.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/lac/sparse_matrix.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/lac/vector.h>

#include "../laplacian/interface/laplacian_interface.h"
#include "../source/interface/source_interface.h"

//! Сборка матриц и векторов
namespace Assembler
{
    template <u8 dim>
    prmt::Report assemble_matrix (
            dealii::SparseMatrix<dbl> &matrix,
            LaplacianInterface<dim> &laplacian,
            const dealii::DoFHandler<dim> &dof_handler)
    {
        cu8 dofs_per_cell = laplacian .get_dofs_per_cell ();

        dealii::FullMatrix<dbl> cell_matrix (dofs_per_cell, dofs_per_cell);
        std::vector<u32> local_dof_indices (dofs_per_cell);

        auto cell = dof_handler.begin_active();
        auto endc = dof_handler.end();
        for (; cell != endc; ++cell)
        {
            cell_matrix = 0;
            laplacian .update_on_cell (cell);

            FOR (i, 0, dofs_per_cell)
                FOR (j, 0, dofs_per_cell)
                    cell_matrix(i,j) += laplacian(i,j);

            cell ->get_dof_indices (local_dof_indices);

            FOR (i, 0, dofs_per_cell)
                FOR (j, 0, dofs_per_cell)
                    matrix .add (local_dof_indices[i],
                            local_dof_indices[j],
                            cell_matrix(i,j));
        };

        REPORT_USE( 
                prmt::Report report;
                report.result = true;
                _return (report););
    };

    template <u8 dim>
    prmt::Report assemble_rhsv (
            dealii::Vector<dbl> &rhsv,
            SourceInterface<dim> &source,
            const dealii::DoFHandler<dim> &dof_handler)
    {
        cu8 dofs_per_cell = source .get_dofs_per_cell ();

        std::vector<u32> local_dof_indices (dofs_per_cell);

        auto cell = dof_handler.begin_active();
        auto endc = dof_handler.end();
        for (; cell != endc; ++cell)
        {
            source .update_on_cell (cell);

            cell ->get_dof_indices (local_dof_indices);

            FOR (i, 0, dofs_per_cell)
                    rhsv(local_dof_indices[i]) += source(i);
        // puts("111");
        };
        // puts("222");

        REPORT_USE( 
                prmt::Report report;
                report.result = true;
                _return (report););
    };

};

#endif
