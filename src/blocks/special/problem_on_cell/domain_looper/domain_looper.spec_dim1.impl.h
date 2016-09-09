/*
 * =====================================================================================
 *
 *       Filename:  domain_looper.spec_dim1.desc.h
 *
 *    Description:  for dim = 1
 *
 *        Version:  1.0
 *        Created:  31.08.2012 11:53:41
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */
#ifndef DOMAIN_LOOPER_SPEC_DIM1_IMPL

#define DOMAIN_LOOPER_SPEC_DIM1_IMPL

//#include "./domain_looper.spec_dim1.desc.h"

using namespace DOMAIN_LOOPER_TOOLS;

    template<bool type_space>
DomainLooper<1, type_space>::DomainLooper()
{

};

template<bool type_space>
void DomainLooper<1, type_space>::loop_domain (
        const dealii::DoFHandler<1> &dof_h, 
        OnCell::BlackOnWhiteSubstituter &bows,
        dealii::CompressedSparsityPattern &csp)
{
    // if ((dof_h.n_boundary_dofs() > 0) && (not csp.empty()))
    // {
        Border<dim> border;
        border = get_borders<dim>(dof_h, csp);

        bows .set_size_of_data (1);

        size_t neighbor;

        for (
                typename dealii::DoFHandler<dim>::active_cell_iterator cell = dof_h.begin_active ();
                cell != dof_h.end ();
                ++cell
            )
        {
            std::cout << cell->vertex(0)[x] << " " << cell ->vertex_dof_index(0,0)  << std::endl;
            std::cout << cell->vertex(1)[x] << " " << cell ->vertex_dof_index(1,0)  << std::endl;
            if (cell->at_boundary())
            {
                if (std::abs(cell->vertex(0)[x] - border[x][wht]) 
                        < MIN_DISTANCE) then
                {
                    bows.white[0] = cell->vertex_dof_index(0,0);
                    continue;
                };

                if (std::abs(cell->vertex(1)[x] == border[x][wht])
                        < MIN_DISTANCE) then
                {
                    bows.white[0] = cell->vertex_dof_index(1,0);
                    continue;
                };

                if (std::abs(cell->vertex(0)[x] == border[x][blc])
                        < MIN_DISTANCE) then
                {
                    bows.black[0] = cell->vertex_dof_index(0,0);
                    neighbor = cell->vertex_dof_index(1,0);
                    continue;
                };

                if (std::abs(cell->vertex(1)[x] == border[x][blc])
                        < MIN_DISTANCE) then
                {
                    bows.black[0] = cell->vertex_dof_index(1,0);
                    neighbor = cell->vertex_dof_index(0,0);
                    continue;
                };
            };
        };

        csp .add (bows.white[0], neighbor);
        csp .add (neighbor, bows.white[0]);

        std::cout << bows.white[0] << std::endl;
        std::cout << bows.black[0] << std::endl;
        std::cout << neighbor << std::endl;

    //     REPORT_USE(
    //             prmt::Report report;
    //             report.result = true;
    //             _return(report);
    //             );
    // }
    // else
    // {
    //     REPORT_USE(
    //             prmt::Report report;
    //             report.result = false;
    //             report.ther_is_a_message = true;
    //             report.message = "";
    //             if (dof_h.n_boundary_dofs() == 0)
    //             report.message += " dof_h is empty";
    //             if (csp.empty())
    //             report.message += " csp is empty";
    //             _return(report);
    //             );
    // };
};

    template<bool type_space>
DomainLooper<1, type_space>::~DomainLooper()
{
};

#endif
