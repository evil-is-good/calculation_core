/*
 * =====================================================================================
 *
 *       Filename:  domain_looper.spec_dim1.desc.h
 *
 *    Description:  for dim = 2
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
#ifndef DOMAIN_LOOPER_SPEC_DIM2_DESC

#define DOMAIN_LOOPER_SPEC_DIM2_DESC

//#include "./domain_looper.desc.h"

namespace OnCell
{
template<bool type_space>
class DomainLooper<2, type_space>
{
    public:
        DomainLooper ();
        ~DomainLooper ();

        const static uint8_t dim = 2;

    //METHODS
    public:
        void loop_domain (const dealii::DoFHandler<dim> &dof_h,
                            OnCell::BlackOnWhiteSubstituter &bows,
                            dealii::DynamicSparsityPattern &csp);
    private:
//        size_t get_num_dofs_on_edge (const dealii::DoFHandler<dim> &dof_h) const;
        
//        void set_size_of_data ();
        
        void add_white (const typename 
                dealii::DoFHandler<dim>::active_cell_iterator &cell, 
                const uint8_t index_of_vertex, const uint8_t side);
        
        void add_black (const typename 
                dealii::DoFHandler<dim>::active_cell_iterator &cell, 
                const uint8_t index_of_vertex, const uint8_t side);;

        void add_blacks_neigbors (const typename 
                dealii::DoFHandler<dim>::active_cell_iterator &cell, 
                const uint8_t index_of_vertex, const uint8_t side);

        void add_white_angle (const typename 
                dealii::DoFHandler<dim>::active_cell_iterator &cell, 
                const uint8_t index_of_vertex);

        void add_black_angle (const typename 
                dealii::DoFHandler<dim>::active_cell_iterator &cell, 
                const uint8_t index_of_vertex);

        void add_black_angle_neigbors (const typename 
                dealii::DoFHandler<dim>::active_cell_iterator &cell, 
                const uint8_t index_of_vertex);

        void set_ratio_black_to_white (OnCell::BlackOnWhiteSubstituter &bows);

        void add_nodes_in_csp (const OnCell::BlackOnWhiteSubstituter &bows,
                dealii::DynamicSparsityPattern &csp);



        const static uint8_t x_is_white = 1;
        const static uint8_t y_is_white = 2;
        const static uint8_t x_is_black = 4;
        const static uint8_t y_is_black = 8;

        const static uint8_t is_white_angle  = x_is_white | y_is_white;
        const static uint8_t is_black_angle  = x_is_black | y_is_black;
        const static uint8_t is_x_gray_angle = x_is_black | y_is_white;
        const static uint8_t is_y_gray_angle = x_is_white | y_is_black;


    //FIELDS
    private:

        const static uint8_t num_angles = 4;
        size_t count_ratios;
        bool this_black_is_exists;
//        size_t num_dofs_on_edge;
//
//        size_t* list_black_dofs[dim][type_space + 1];
//        size_t* list_white_dofs[dim][type_space + 1];
//
//        double* list_black_coor[dim];
//        double* list_white_coor[dim];

        std::vector<size_t> list_black_dofs[dim][type_space + 1];
        std::vector<size_t> list_white_dofs[dim][type_space + 1];

        std::vector<double> list_black_coor[dim];
        std::vector<double> list_white_coor[dim];

        size_t list_black_angular_dofs[type_space + 1][num_angles];
        size_t white_angular_dof[type_space + 1];

//        std::vector<size_t>* 
//            list_blacks_neighbors[dim][type_space + 1];
        std::vector<std::vector<size_t> >
            list_blacks_neighbors[dim][type_space + 1];
        std::vector<size_t> 
            list_black_angular_neighbors[type_space + 1][num_angles];

//        size_t white_count[dim];
//        size_t black_count[dim];
        size_t black_angular_count;

    OPERATOR_REPORT;
};
};
#endif
