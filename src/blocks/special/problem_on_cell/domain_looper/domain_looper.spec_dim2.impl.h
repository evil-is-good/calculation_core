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
#ifndef DOMAIN_LOOPER_SPEC_DIM2_IMPL

#define DOMAIN_LOOPER_SPEC_DIM2_IMPL

//#include "./domain_looper.spec_dim2.desc.h"

//namespace DOMAIN_LOOPER_TOOLS
//{
namespace OnCell
{
using namespace DOMAIN_LOOPER_TOOLS;

    template<bool type_space>
DomainLooper<2, type_space>::DomainLooper()
//    :
//        num_dofs_on_edge(0)
{
    black_angular_count = 0;
    for(size_t i = 0; i < num_angles; ++i)
        list_black_angular_dofs[0][i] = 1000000000000;

};


//template<bool type_space>
//size_t DomainLooper<2, type_space>::
//get_num_dofs_on_edge (const dealii::DoFHandler<dim> &dof_h) const
//{
//    return ((dof_h.n_boundary_dofs() / (type_space + 1)) - num_angles) / (dim * 2);
//};

//template<bool type_space>
//void DomainLooper<2, type_space>::set_size_of_data ()
//{
//    for (uint8_t i = 0; i < type_space + 1; ++i) // for vector
//    {
//        for (uint8_t j = 0; j < dim; ++j) // sides
//        {
//            list_black_dofs[j][i] = new size_t [num_dofs_on_edge];
//            list_white_dofs[j][i] = new size_t [num_dofs_on_edge];
//            list_blacks_neighbors[j][i] = 
//                new std::vector<size_t> [num_dofs_on_edge]; 
//
//            for (size_t k = 0; k < num_dofs_on_edge; ++k)
//            {
//                list_white_dofs[j][i][k] = 0;
//                list_black_dofs[j][i][k] = 0;
//            };
//        };
//    };
//
//    for (uint8_t i = 0; i < dim; ++i)
//    {
//        list_black_coor[i] = new double [num_dofs_on_edge];
//        list_white_coor[i] = new double [num_dofs_on_edge];
//
//        for (size_t j = 0; j < num_dofs_on_edge; ++j)
//        {
//            list_white_coor[i][j] = 0.0;
//            list_black_coor[i][j] = 0.0;
//        };
//
//        white_count[i] = 0;
//        black_count[i] = 0;
//    };
//
//    black_angular_count = 0;
//};


template<bool type_space>
void DomainLooper<2, type_space>::
add_white (const typename 
        dealii::DoFHandler<dim>::active_cell_iterator &cell, 
        const uint8_t index_of_vertex, const uint8_t side)
{
    bool this_white_is_exists = false;
//    printf("111\n");

    for (size_t i = 0; i < list_white_dofs[side][0].size(); ++i)
    {
        if (list_white_dofs[side][0][i] == 
                cell ->vertex_dof_index (index_of_vertex, 0))
            this_white_is_exists = true;
    };
//    printf("111\n");

    if (not this_white_is_exists)
    {
        for (st i = 0; i < type_space + 1; ++i) 
        {
            list_white_dofs[side][i] .push_back (
                cell ->vertex_dof_index (index_of_vertex, i));
        };

        list_white_coor[side] .push_back (
            cell ->vertex (index_of_vertex)[not side]);
        // так как это сторона, например, х, то записать нужно координату у
    };
//    printf("111\n");
};

template<bool type_space>
void DomainLooper<2, type_space>::
add_black (const typename 
        dealii::DoFHandler<dim>::active_cell_iterator &cell, 
        const uint8_t index_of_vertex, const uint8_t side)
{
    this_black_is_exists = false;

    for (size_t i = 0; i < list_black_dofs[side][0].size(); ++i) //black_count[side]
    {
        if (list_black_dofs[side][0][i] == 
                cell ->vertex_dof_index (index_of_vertex, 0))
            this_black_is_exists = true;
    };

    if (not this_black_is_exists)
    {
        for (st i = 0; i < type_space + 1; ++i) 
        {
//            list_black_dofs[side][i][black_count[side]] =
//                cell ->vertex_dof_index (index_of_vertex, i);
            list_black_dofs[side][i] .push_back (
                cell ->vertex_dof_index (index_of_vertex, i));

            list_blacks_neighbors[side][i] .push_back (
                    std::vector<size_t>(0));
        };

//        list_black_coor[side][black_count[side]] =
//            cell ->vertex (index_of_vertex)[ not side];
        list_black_coor[side] .push_back (
            cell ->vertex (index_of_vertex)[not side]);
        // так как это сторона, например, х, то записать нужно координату у


//        ++black_count[side];

    };
};

template<bool type_space>
void DomainLooper<2, type_space>::
add_blacks_neigbors (const typename 
        dealii::DoFHandler<dim>::active_cell_iterator &cell, 
        const uint8_t index_of_vertex, const uint8_t side)
{
    size_t n = 0;

    if (not this_black_is_exists)
        n = list_black_dofs[side][0].size() - 1;
    else
    {
        for (size_t i = 0; i < list_black_dofs[side][0].size(); ++i) 
        {
            if (list_black_dofs[side][0][i] == 
                    cell ->vertex_dof_index (index_of_vertex, 0))
                n = i;
        };
    };

//    if (not this_black_is_exists)
        for (size_t i = 0; i < type_space + 1; ++i) 
            for (size_t j = 0; 
                    j < dealii::GeometryInfo<dim>::vertices_per_cell; ++j)
                if (j != index_of_vertex)
                {
//                    printf("%ld %ld\n", j, k);
                    // list_blacks_neighbors[side][k][i] .push_back (
                    // cell ->vertex_dof_index (j, k));
                    for (size_t k = 0; k < type_space + 1; ++k) 
                        list_blacks_neighbors[side][i][n] .push_back (
                                cell ->vertex_dof_index (j, k));
//                    printf("|||||||\n");
                };
//    printf("bl %ld %d\n", black_count[side], type_space);
//    for (size_t i = 0; i < list_black_dofs[side][0].size() + 1; ++i)
//    {
//        printf("%ld\n", i);
//        if (list_black_dofs[side][0][i] == 
//                cell ->vertex_dof_index (index_of_vertex, 0))
//        {
//            for (size_t j = 0; j < type_space + 1; ++j) 
//            {
//                list_blacks_neighbors[side][j] .push_back (
//                        std::vector<size_t>(0));
//
//                for (size_t k = 0; 
//                        k < dealii::GeometryInfo<dim>::vertices_per_cell; ++k)
//                    if (k != index_of_vertex)
//                    {
//                        printf("%ld %ld\n", j, k);
//                        // list_blacks_neighbors[side][k][i] .push_back (
//                        // cell ->vertex_dof_index (j, k));
//                        list_blacks_neighbors[size][j][i] .push_back (
//                                cell ->vertex_dof_index (j, k));
//                        printf("|||||||\n");
//                    };
//            };
//
//            break;
//        };
//    };
};

template<bool type_space>
void DomainLooper<2, type_space>::
add_white_angle (const typename 
        dealii::DoFHandler<dim>::active_cell_iterator &cell, 
        const uint8_t index_of_vertex)
{
    for (st i = 0; i < type_space + 1; ++i) 
        white_angular_dof[i] = cell ->vertex_dof_index (index_of_vertex, i);
};

template<bool type_space>
void DomainLooper<2, type_space>::
add_black_angle (const typename 
        dealii::DoFHandler<dim>::active_cell_iterator &cell, 
        const uint8_t index_of_vertex)
{
    this_black_is_exists = false;

    for (size_t i = 0; i < num_angles; ++i) 
    {
//        printf("%u %ld\n", list_black_angular_dofs[0][i],
//                cell ->vertex_dof_index (index_of_vertex, 0));
        if (list_black_angular_dofs[0][i] == 
                cell ->vertex_dof_index (index_of_vertex, 0))
            this_black_is_exists = true;
    };
//    printf("//////////////////////////////////////olololo_1\n");
    if (not this_black_is_exists)
        for (st i = 0; i < type_space + 1; ++i) 
            list_black_angular_dofs[i][black_angular_count] =
                cell-> vertex_dof_index (index_of_vertex, i);
};

template<bool type_space>
void DomainLooper<2, type_space>::
add_black_angle_neigbors (const typename 
        dealii::DoFHandler<dim>::active_cell_iterator &cell, 
        const uint8_t index_of_vertex)
{
//    printf("////////////////////////////olololo_2 %ld\n", black_angular_count);
    if (not this_black_is_exists)
    {
        for (st i = 0; i < type_space + 1; ++i) 
            for (st j = 0; j < dealii::GeometryInfo<dim>::vertices_per_cell; ++j)
                if (j != index_of_vertex)
                    for (st k = 0; k < type_space + 1; ++k) 
                        list_black_angular_neighbors[i][black_angular_count] .push_back (
                                cell -> vertex_dof_index (j, k));

        ++black_angular_count;
    };
};

template<bool type_space>
void DomainLooper<2, type_space>::
set_ratio_black_to_white (OnCell::BlackOnWhiteSubstituter &bows)
{
    //            std::cout << "white_count_x " << white_count[x] << std::endl;
    //            std::cout << "white_count_y " << white_count[y] << std::endl;
    //            std::cout << "black_count_x " << black_count[x] << std::endl;
    //            std::cout << "black_count_y " << black_count[y] << std::endl;
    //            std::cout << "angle_count " << black_angular_count << std::endl;
    //            std::cout << "num_dofs_on_edge " << num_dofs_on_edge << std::endl;
    //
    //            for (uint8_t i = 0; i < dim; ++i)
    //                for (uint8_t j = 0; j < type_space + 1; ++j)
    //                    for (uint8_t k = 0; k < num_dofs_on_edge ; ++k)
    //                    {
    //                        std::cout 
    //                            << "white_index " 
    //                            << list_white_dofs[i][j][k] 
    //                            << " white_coor "
    //                            << list_white_coor[i][k]
    //                            << std::endl;
    //                        std::cout 
    //                            << "black_index " 
    //                            << list_black_dofs[i][j][k] 
    //                            << " black_coor "
    //                            << list_black_coor[i][k]
    //                            << std::endl;
    //                        printf("Neighbors: ");
    //                        for (uint8_t l = 0; l < list_blacks_neighbors[i][j][k].size() ; ++l)
    //                            printf("%d ", list_blacks_neighbors[i][j][k][l]);
    //                        printf("\n");
    //                    };
    //            printf("Anular\n");
    //            for (uint8_t i = 0; i < type_space + 1; ++i)
    //                for (uint8_t j = 0; j < (num_angles - 1); ++j) 
    //                {
    //                    printf("%d: ", j);
    //                    for (uint8_t k = 0; k < list_black_angular_neighbors[i][j].size(); ++k)
    //                        printf("%d ", list_black_angular_neighbors[i][j][k]);
    //                    printf("\n");
    //                }
    //            printf("\n");
//////////////
    count_ratios = 0;

    for (st i = 0; i < type_space + 1; ++i)
        for (st j = 0; j < num_angles - 1; ++j)
        {
            // bows.white[count_ratios] = white_angular_dof[i];
            // bows.black[count_ratios] = list_black_angular_dofs[i][j];
            // bows.white.push_back(white_angular_dof[i]);
            // bows.black.push_back(list_black_angular_dofs[i][j]);
            bows.add_white_and_black(white_angular_dof[i], list_black_angular_dofs[i][j]);

            ++count_ratios;
        };

    for (st i = 0; i < dim; ++i)
        for (size_t j = 0; j < list_white_coor[i].size(); ++j)
        {
            size_t n = 0;
            for (size_t k = 0; k < list_black_coor[i].size(); ++k)
                if (std::abs(list_white_coor[i][j] - list_black_coor[i][k]) < 
                        MIN_DISTANCE)
                    n = k;

            for (st k = 0; k < type_space + 1; ++k)
            {
                // bows.white[count_ratios] = list_white_dofs[i][k][j];
                // bows.black[count_ratios] = list_black_dofs[i][k][n];
                // bows.white.push_back(list_white_dofs[i][k][j]);
                // bows.black.push_back(list_black_dofs[i][k][n]);
                bows.add_white_and_black(list_white_dofs[i][k][j], list_black_dofs[i][k][n]);

                ++count_ratios;
//                printf("C_RATIOS %ld\n", count_ratios);
            };
        };
///////////////////
//            for (size_t i = 0; i < bows.size; ++i)
//                printf("w=%ld b=%ld\n", bows.white[i], bows.black[i]);

};

template<bool type_space>
void DomainLooper<2, type_space>::
add_nodes_in_csp (const OnCell::BlackOnWhiteSubstituter &bows, 
        dealii::DynamicSparsityPattern &csp)
{
    for (st i = 0; i < type_space + 1; ++i)
        for (st j = 0; j < (num_angles - 1); ++j) 
            for (st k = 0; k < list_black_angular_neighbors[i][j].size(); ++k)
                list_black_angular_neighbors[i][j][k] = 
                    bows .subst (list_black_angular_neighbors[i][j][k]);

    for (st i = 0; i < dim; ++i)
        for (st j = 0; j < type_space + 1; ++j)
            for (st k = 0; k < list_black_dofs[i][j].size() ; ++k)
                for (st l = 0; l < list_blacks_neighbors[i][j][k].size(); ++l)
                    list_blacks_neighbors[i][j][k][l] = 
                        bows .subst (list_blacks_neighbors[i][j][k][l]);

    //            printf("Anular\n");
    //            for (uint8_t i = 0; i < type_space + 1; ++i)
    //                for (uint8_t j = 0; j < (num_angles - 1); ++j) 
    //                {
    //                    printf("%d: ", j);
    //                    for (uint8_t k = 0; k < list_black_angular_neighbors[i][j].size(); ++k)
    //                        printf("%d ", list_black_angular_neighbors[i][j][k]);
    //                    printf("\n");
    //                }
    //            printf("\n");
    //
    //            for (uint8_t i = 0; i < dim; ++i)
    //                for (uint8_t j = 0; j < type_space + 1; ++j)
    //                    for (uint8_t k = 0; k < num_dofs_on_edge ; ++k)
    //                    {
    //                        printf("Neighbors to %d: ", list_black_dofs[i][j][k]);
    //                        for (uint8_t l = 0; l < list_blacks_neighbors[i][j][k].size() ; ++l)
    //                            printf("%d ", list_blacks_neighbors[i][j][k][l]);
    //                        printf("\n");
    //                    };

//    size_t count_edge_dof = 0;
//
//    for (uint8_t i = 0; i < (num_angles - 1); ++i)
//        for (uint8_t j = 0; j < type_space + 1; ++j)
//        {
//            for (uint8_t k = 0; k < type_space + 1; ++k)
//                for (uint8_t l = 0; l < list_black_angular_neighbors[k][i].size(); ++l)
//                {
//                    csp .add (bows.white[count_edge_dof], 
//                            list_black_angular_neighbors[k][i][l]);
//                    csp .add (list_black_angular_neighbors[k][i][l], 
//                            bows.white[count_edge_dof]);
//
//                };
//
//            ++count_edge_dof;
//        };
//
//    for (uint8_t i = 0; i < dim; ++i)
//        for (size_t j = 0; j < list_white_coor[i].size(); ++j)
//        {
//            size_t n = 0;
//            for (size_t k = 0; k < list_white_coor[i].size(); ++k)
//                if (bows.black[count_edge_dof] == list_black_dofs[i][0][k])
//                    n = k;
//
//            for (uint8_t k = 0; k < type_space + 1; ++k)
//            {
//                for (uint8_t l = 0; l < type_space + 1; ++l)
//                    for (uint8_t m = 0; m < list_blacks_neighbors[i][l][n].size() ; ++m)
//                    {
//                        csp .add (bows.white[count_edge_dof],
//                                list_blacks_neighbors[i][l][n][m]);
//                        csp .add (list_blacks_neighbors[i][l][n][m],
//                                bows.white[count_edge_dof]);
//                    };
//
//                ++count_edge_dof;
//            };
//        };

//////////////////
    size_t n = 0;
    for (st i = 0; i < type_space + 1; ++i)
        for (st j = 0; j < num_angles - 1; ++j)
        {
            for (st k = 0; k < list_black_angular_neighbors[i][j].size();
                    ++k)
            {
//                printf("w=%ld n=%ld\n", bows.white[n],list_black_angular_neighbors[i][j][k]);
                csp .add (bows.white[n], 
                        list_black_angular_neighbors[i][j][k]);

                csp .add (list_black_angular_neighbors[i][j][k], 
                        bows.white[n]);
            };

            ++n;
        };
//    printf("N=%ld\n", n);

    for (st i = 0; i < dim; ++i)
        for (st j = 0; j < type_space + 1; ++j)
            for (st k = 0; k < count_ratios; ++k)
                for (st l = 0; l < list_black_dofs[i][j].size(); ++l)
                    if (bows.black[k] == list_black_dofs[i][j][l])
                    {
                        for (st m = 0; m < 
                                list_blacks_neighbors[i][j][l].size() ; ++m)
                        {
//                            printf("w=%ld n=%ld\n", 
//                                    bows.white[k],
//                                    list_blacks_neighbors[i][j][l][m]);
                            csp .add (bows.white[k],
                                    list_blacks_neighbors[i][j][l][m]);
                            csp .add (list_blacks_neighbors[i][j][l][m],
                                    bows.white[k]);
                        };

                        break;
                    };

};


template<bool type_space>
void DomainLooper<2, type_space>::loop_domain (
        const dealii::DoFHandler<2> &dof_h, 
        OnCell::BlackOnWhiteSubstituter &bows,
        dealii::DynamicSparsityPattern &csp)
{
    // if ((dof_h.n_boundary_dofs() > 0) && (not csp.empty()))
    // {
        Border<dim> border;

        border = get_borders<dim> (dof_h, csp);

//        num_dofs_on_edge = get_num_dofs_on_edge (dof_h);

//        set_size_of_data ();

        FILE *F;
        F = fopen ("deb-out","w");
//        printf("2.1\n");
        int count = 0;

        for (
                typename dealii::DoFHandler<dim>::active_cell_iterator cell = dof_h.begin_active ();
                cell != dof_h.end ();
                ++cell
            )
        {
            ++count;
            // printf("pre ");
            // for (uint8_t i = 0; 
            //         i < dealii::GeometryInfo<dim>::vertices_per_cell; ++i)
            //     printf("%d ", cell ->vertex_dof_index (i,0));
            // printf("\n");
//            printf("count %d\n", count);
            if (at_boundary<dim>(cell, border))
//            if (cell->has_boundary_lines())
//            if (cell->at_boundary())
            {
            // printf("post ");
            //     for (uint8_t i = 0; 
            //             i < dealii::GeometryInfo<dim>::vertices_per_cell; ++i)
            //         printf("%d ", cell ->vertex_dof_index (i,0));
            // printf("\n");
                for (st i = 0; 
                        i < dealii::GeometryInfo<dim>::vertices_per_cell; ++i)
                {
//                    printf("////////////////////////\n");
//                    printf("i=%d coorx=%f coory=%f\n", i, cell -> vertex (i)[0], cell -> vertex (i)[1]);
//                    printf("index_1=%d\n", cell ->vertex_dof_index (i,0));
//                    fprintf(:bF,"////////////////////////\n");
//                    fprintf(F,"i=%d j=%d coor=%f\n", i, 0, cell -> vertex (i)[0]);
//                    fprintf(F,"i=%d j=%d coor=%f\n", i, 1, cell -> vertex (i)[1]);
//                    fprintf(F,"index_1=%d\n", cell ->vertex_dof_index (i,0));
//                    printf("index_2=%d\n", cell ->vertex_dof_index (i,1));
                    uint8_t x_check_on_white = 
                        (uint8_t)(std::abs(cell->vertex(i)[x] - border[x][wht])
                         < MIN_DISTANCE);
                    uint8_t y_check_on_white = 
                        (uint8_t)(std::abs(cell->vertex(i)[y] - border[y][wht])
                         < MIN_DISTANCE);
                    uint8_t x_check_on_black = 
                        (uint8_t)(std::abs(cell->vertex(i)[x] - border[x][blc])
                         < MIN_DISTANCE);
                    uint8_t y_check_on_black = 
                        (uint8_t)(std::abs(cell->vertex(i)[y] - border[y][blc])
                         < MIN_DISTANCE);

                    x_check_on_white <<= 0;
                    y_check_on_white <<= 1;
                    x_check_on_black <<= 2;
                    y_check_on_black <<= 3;
//                    printf("%d %d %d %d \n", x_check_on_white,
//                                             y_check_on_white,
//                                             x_check_on_black,
//                                             y_check_on_black);
//            printf("count2 %d\n", count);


//                    if (count == 77)
//        printf("2.1.1 %d\n",
//                    x_check_on_white + y_check_on_white +
//                            x_check_on_black + y_check_on_black);
                    switch (x_check_on_white + y_check_on_white +
                            x_check_on_black + y_check_on_black)
                    {
                        case x_is_white:      add_white (cell, i, x);
//                                              printf("white_x\n");
                                              //                                                      printf("x_is_white\n");
                                              break;
                        case y_is_white:      add_white (cell, i, y);
//                                              printf("white_y\n");
                                              //                                                      printf("y_is_white\n");
                                              break;
                        case x_is_black:      add_black (cell, i, x);
                                              add_blacks_neigbors (cell, i, x);
//                                              printf("black_x\n");
                                              //                                                      printf("x_is_black\n");
                                              break;
                        case y_is_black:      
//                                              printf("d1\n");
                                              add_black (cell, i, y);
//                                              printf("d2\n");
                                              add_blacks_neigbors (cell, i, y);
//                                              printf("black_y\n");
                                              //                                                      printf("y_is_black\n");
                                              break;
                        case is_white_angle:  add_white_angle (cell, i);
//                                              printf("e\n");
                                              //                                                      printf("is_white_angle\n");
                                              break;
                        case is_black_angle:  add_black_angle (cell, i);
                                              add_black_angle_neigbors (cell, i);
//                                              printf("f\n");
                                              //                                                      printf("is_black_angle\n");
                                              break;
                        case is_x_gray_angle: add_black_angle (cell, i);
                                              add_black_angle_neigbors (cell, i);
//                                              printf("g\n");
                                              //                                                      printf("is_x_gray_angle\n");
                                              break;
                        case is_y_gray_angle: add_black_angle (cell, i);
                                              add_black_angle_neigbors (cell, i);
//                                              printf("l\n");
                                              //                                                      printf("is_y_gray_angle\n");
                                              break;
                    };
//        printf("2.1.2\n");
                };
            };
        };
//    printf("2.2\n");

    fclose (F);

        size_t num_ratios = 0;
        for (st i = 0; i < dim; ++i)
            for (st j = 0; j < type_space + 1; ++j)
                num_ratios += list_white_dofs[i][j].size();
        if (type_space)
            num_ratios += (num_angles - 1) * dim;
        else
            num_ratios +=  num_angles - 1;

        // bows .set_size_of_data (num_ratios);//((2 * num_dofs_on_edge + (num_angles - 1)) * (type_space + 1));
//    printf("2.3 %ld %d\n", num_ratios, int(type_space));
 
        set_ratio_black_to_white (assigned_to bows);
//    printf("2.4\n");

        add_nodes_in_csp (bows, assigned_to csp);
//    printf("2.5\n");


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
DomainLooper<2, type_space>::~DomainLooper()
{
//    if (num_dofs_on_edge > 0)
//    {
//        for (uint8_t i = 0; i < dim; ++i)
//        {
//            for (uint8_t j = 0; j < type_space + 1; ++j)
//            {
//                delete [] list_black_dofs[i][j];
//
//                delete [] list_white_dofs[i][j];
//
//                for (uint8_t k = 0; k < num_dofs_on_edge; ++k)
//                    list_blacks_neighbors[i][j][k].clear();
//                delete [] list_blacks_neighbors[i][j];
//            };
//
//            delete [] list_black_coor[i];
//
//            delete [] list_white_coor[i];
//        };
//    };
};
};
#endif
