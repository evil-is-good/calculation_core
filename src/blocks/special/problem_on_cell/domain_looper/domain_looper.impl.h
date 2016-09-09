#ifndef DOMAIN_LOOPER_IMPL

#define DOMAIN_LOOPER_IMPL

namespace OnCell
{
namespace DOMAIN_LOOPER_TOOLS
{
    template<uint8_t dim>
    Border<dim> get_borders (const dealii::DoFHandler<dim> &dof_h,
            dealii::DynamicSparsityPattern &csp)
    {
        Border<dim> border;

        const uint8_t x = 0;
        const uint8_t y = 1;
        const uint8_t z = 2;
        const uint8_t wht = 0;
        const uint8_t blc = 1;


        for (uint8_t i = 0; i < dim; ++i)
        {
            border.coor[i][wht] =  0xffffffff * 1.0;
            border.coor[i][blc] = -0xffffffff * 1.0;
        };

        
		for (
				typename dealii::DoFHandler<dim>::active_cell_iterator cell = dof_h.begin_active ();
				cell != dof_h.end ();
				++cell
			)
        {
            if (cell->at_boundary())
			    for(uint8_t i = 0; i < dealii::GeometryInfo<dim>::vertices_per_cell; ++i)
                    for (uint8_t j = 0; j < dim; ++j)
                    {
                        if (cell->vertex(i)[j] < border.coor[j][wht])
                            border.coor[j][wht] = cell->vertex(i)[j];
                        if (cell->vertex(i)[j] > border.coor[j][blc])
                            border.coor[j][blc] = cell->vertex(i)[j];
                    };
        };

        return border;
    };
};


template<uint8_t dim, bool type_space>
DomainLooper<dim, type_space>::DomainLooper()
{

};

template<uint8_t dim, bool type_space>
void DomainLooper<dim, type_space>::loop_domain (
        const dealii::DoFHandler<dim> &dof_h, 
        OnCell::BlackOnWhiteSubstituter &bows,
        dealii::DynamicSparsityPattern &csp)
{
//    std::cout << "ssdf" << std::endl;
//    if ((size == 0) && (dof_h.n_boundary_dofs() > 0) && (not csp.empty()))
//    {
//        
//        DomainLoopsAdditionalTools::IndexAndCoor<dim, type_space>
//            temp_container_black_dofs;
//
//        DomainLoopsAdditionalTools::IndexAndCoor<dim, type_space>
//            temp_container_white_dofs;
//
//        DomainLoopsAdditionalTools::Neighbor<dim, type_space>
//            neighbors_for_bleck_dofs;
//
//
//
//        REPORT DomainLoopsAdditionalTools 
//            ::get_list_black_and_white_and_blecks_neighbors <dim, type_space> (
//                dof_h,
//                and_assigned_to
//                temp_container_black_dofs,
//                temp_container_white_dofs,
//                neighbors_for_bleck_dofs);
//
//
//
//
//
//
//
//
//        REPORT DomainLoopsAdditionalTools 
//            ::get_blecks_neighbors <dim, type_space> (
//                dof_h,
//                temp_container_black_dofs,
//                and_assigned_to
//                neighbors_for_bleck_dofs);
//
//
//
//
//        size_t n = DomainLoopsAdditionalTools 
//            ::num_black_or_white_dofs<dim, type_space> (); 
//
//        black = new size_t [n];
//        white = new size_t [n];
//
//        
//
//
//        REPORT DomainLoopsAdditionalTools 
//            ::get_ratio_black_to_white_and_modify_csp <dim,type_space> (
//                    temp_container_black_dofs,
//                    temp_container_white_dofs,
//                    neighbors_for_bleck_dofs,
//                    and_assigned_to
//                    white,
//                    black,
//                    csp);
//
//
//
//
//        REPORT_USE(
//                prmt::Report report;
//                report.result = _report.result;
//                _return(report);
//                );
//    }
//    else
//    {
//       REPORT_USE(
//               prmt::Report report;
//               report.result = false;
//               report.ther_is_a_message = true;
//               report.message = "";
//               if (size > 0)
//                   report.message = "second call";
//               if (dof_h.n_boundary_dofs() == 0)
//                   report.message += " dof_h is empty";
//               if (csp.empty())
//                   report.message += " csp is empty";
//               _return(report);
//               );
//    };
       REPORT_USE(
               prmt::Report report;
               report.result = false;
               report.ther_is_a_message = true;
               report.message = "Dimension is not supported";
               _return(report);
               )
};

template<uint8_t dim, bool type_space>
DomainLooper<dim, type_space>::~DomainLooper()
{
};
};


#endif

