#ifndef DOMAIN_LOOPER_DESC

#define DOMAIN_LOOPER_DESC

#include "../../../../../../prmt_sintactic_addition/prmt_sintactic_addition.h"
#include "../black_on_white_substituter/black_on_white_substituter.h"

#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_levels.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/base/geometry_info.h>

namespace OnCell
{
namespace DOMAIN_LOOPER_TOOLS
{

    const uint8_t x = 0;
    const uint8_t y = 1;
    const uint8_t z = 2;
    const uint8_t wht = 0;
    const uint8_t blc = 1;

    const double MIN_DISTANCE = 1e-10; // расстояние между "одинаковыми" точками

    template<uint8_t dim>
    struct Border
    {
        double coor[dim][2]; // 0 - black, 1 - white;

        inline double* operator[] (int i)
        {
            return coor[i];
        };
    };

    template<uint8_t dim>
    Border<dim> get_borders (const dealii::DoFHandler<dim> &dof_h,
            dealii::DynamicSparsityPattern &csp);
    
    template<uint8_t dim>
    bool at_boundary (
            const typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
            const Border<dim> &border)
    {
        FOR_I (0, 4)//dealii::GeometryInfo<dim>::vertices_per_cell)
        {
//                if (cell->vertex_dof_index(i,0) == 6)
//                    printf("%ld\n", i);
            FOR_J (0, dim)
            {
//                if (cell->vertex_dof_index(i,0) == 6)
//                    printf("%ld %ld %f %f %f %d\n", 
//                            i, j, fabs(cell->vertex(i)[j] - border.coor[j][0]),
//                            dealii::GeometryInfo<dim>::vertices_per_cell);
                if (
                        (fabs(cell->vertex(i)[j] - border.coor[j][0]) < MIN_DISTANCE) or
                        (fabs(cell->vertex(i)[j] - border.coor[j][1]) < MIN_DISTANCE)
                   )
                    return true;
            };
        };

//            if (cell->at_boundary(i))
//                return true;
        return false;
    };

};

template<uint8_t dim, bool type_space>
class DomainLooper
{
    public:
        DomainLooper ();
        ~DomainLooper ();

    //METHODS
    public:
        void loop_domain (const dealii::DoFHandler<dim> &dof_h,
                            OnCell::BlackOnWhiteSubstituter &bows,
                            dealii::DynamicSparsityPattern &csp);

    OPERATOR_REPORT;
};
};
#endif
