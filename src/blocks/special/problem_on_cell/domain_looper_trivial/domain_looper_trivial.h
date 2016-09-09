#ifndef DOMAIN_LOOPER_TRIVIAL_H_4C8LSHXA
#define DOMAIN_LOOPER_TRIVIAL_H_4C8LSHXA

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
    template<u8 dim, u8 type_space>
        class CubeOfNumbers
        {
            private:
                CubeOfNumbers() {};
        };

    template<u8 type_space>
        class CubeOfNumbers<2, type_space>
        {
            public:
                CubeOfNumbers(cst level) : cube_level(level), size(2_pow(level) + 1) 
            {
                number .resize (size);
                for (auto && line : number)
                {
                    line .resize (size);
                    for (auto && point : line)
                    {
                        FOR(i, 0, type_space)
                        {
                            point[i] = 0;
                        };
                    };
                };
                numbered_pice(level, 0, 2_pow(level), 0, 2_pow(level));
            };

                void numbered_pice (cst lvl, cst b_x, cst e_x, cst b_y, cst e_y)
                {
                    if (lvl > 0)
                    {
                        // |----->x  Coor
                        // |
                        // |
                        // \/y
                        //
                        // 1 0 | 0 1  Order
                        // 0 0 | 0 0
                        // ---------
                        // 0 0 | 0 0
                        // 1 0 | 0 1

                        cst m_x = (e_x - b_x) / 2 + b_x;
                        cst m_y = (e_y - b_y) / 2 + b_y;

                        numbered_pice(lvl-1, b_x, m_x, b_y, m_y); numbered_pice(lvl-1, m_x, e_x, b_y, m_y);
                        numbered_pice(lvl-1, b_x, m_x, m_y, e_y); numbered_pice(lvl-1, m_x, e_x, m_y, e_y);
                    }
                    else
                    {
                        numbered_point(b_x, b_y); numbered_point (e_x, b_y);
                        numbered_point(b_x, e_y); numbered_point (e_x, e_y);
                        //if (number[b_x][b_y] == 0) {FOR(i,0,ts){number[b_x][b_y][i] = n++;};}; if (number[e_x][b_y] == 0) {number[e_x][b_y] = n++;};
                        //if (number[b_x][e_y] == 0) {number[b_x][e_y] = n++;}; if (number[e_x][e_y] == 0) {number[e_x][e_y] = n++;};
                    };
                };

                void numbered_point (cst nx, cst ny)
                {
                    if (number[nx][ny][0] == 0) 
                    {
                        FOR(i, 0, type_space)
                        {
                            number[nx][ny][i] = n++;
                        };
                    };
                };

                vec<vec<arr<st, type_space>>> number;
                st n = 0;
                cst cube_level;
                cst size;

            private:
                CubeOfNumbers() : cube_level(0) {};
        };

    template<u8 type_space>
        class CubeOfNumbers<3, type_space>
        {
            public:
                CubeOfNumbers(cst level) : cube_level(level), size(2_pow(level) + 1) 
            {
                number .resize (size);
                for (auto && layer : number)
                {
                    layer .resize (size);
                    for (auto && line : layer)
                    {
                        line .resize (size);
                        for (auto && point : line)
                        {
                            // FOR(i, 0, type_space)
                            //{
                            point[0] = 0;
                            // };
                        };
                    };
                };
                numbered_pice(level, 0, 2_pow(level), 0, 2_pow(level), 0, 2_pow(level));
            };

                void numbered_pice (cst lvl, cst b_x, cst e_x, cst b_y, cst e_y, cst b_z, cst e_z)
                {
                    if (lvl > 0)
                    {
                        // |----->x  Coor
                        // |
                        // |
                        // \/y
                        //
                        // 1 0 | 0 1  Order
                        // 0 0 | 0 0
                        // ---------
                        // 0 0 | 0 0
                        // 1 0 | 0 1

                        cst m_x = (e_x - b_x) / 2 + b_x;
                        cst m_y = (e_y - b_y) / 2 + b_y;
                        cst m_z = (e_z - b_z) / 2 + b_z;

                        // bottom layer
                        numbered_pice(lvl-1, b_x, m_x, b_y, m_y, b_z, m_z); numbered_pice(lvl-1, m_x, e_x, b_y, m_y, b_z, m_z);
                        numbered_pice(lvl-1, b_x, m_x, m_y, e_y, b_z, m_z); numbered_pice(lvl-1, m_x, e_x, m_y, e_y, b_z, m_z);

                        // top layer
                        numbered_pice(lvl-1, b_x, m_x, b_y, m_y, m_z, e_z); numbered_pice(lvl-1, m_x, e_x, b_y, m_y, m_z, e_z);
                        numbered_pice(lvl-1, b_x, m_x, m_y, e_y, m_z, e_z); numbered_pice(lvl-1, m_x, e_x, m_y, e_y, m_z, e_z);
                    }
                    else
                    {
                        numbered_point(b_x, b_y, b_z); numbered_point (e_x, b_y, b_z);
                        numbered_point(b_x, e_y, b_z); numbered_point (e_x, e_y, b_z);

                        numbered_point(b_x, b_y, e_z); numbered_point (e_x, b_y, e_z);
                        numbered_point(b_x, e_y, e_z); numbered_point (e_x, e_y, e_z);
                    };
                };

                void numbered_point (cst nx, cst ny, cst nz)
                {
                    if (number[nx][ny][nz][0] == 0) 
                    {
                        FOR(i, 0, type_space)
                        {
                            number[nx][ny][nz][i] = n++;
                        };
                    };
                };

                vec<vec<vec<arr<st, type_space>>>> number;
                st n = 0;
                cst cube_level;
                cst size;

            private:
                CubeOfNumbers() : cube_level(0) {};
        };

    template<u8 dim, u8 type_space>
        class DomainLooperTrivial
        {
            private:
                DomainLooperTrivial () {};
        };

    template<u8 type_space>
        class DomainLooperTrivial<2, type_space>
        {
            public:
                DomainLooperTrivial () {};

                void loop_domain (const dealii::DoFHandler<2> &dof_h,
                        OnCell::BlackOnWhiteSubstituter &bows,
                        dealii::DynamicSparsityPattern &csp)
                {
                    CubeOfNumbers<dim, type_space> cube(dof_h.get_tria().n_global_levels()-1);

                    cst lst = cube.size - 1;

                    // Добавляем в bows белые и соответствующие им черные точки
                    for (st i = 0; i < type_space; ++i)
                    {
                        // точки на углах
                        bows .add_white_and_black (cube.number[0][0][i], cube.number[0][lst][i]);
                        bows .add_white_and_black (cube.number[0][0][i], cube.number[lst][0][i]);
                        bows .add_white_and_black (cube.number[0][0][i], cube.number[lst][lst][i]);

                        // точки на рёбрах
                        for (st j = 1; j < lst; ++j)
                        {
                            bows .add_white_and_black (cube.number[0][j][i], cube.number[lst][j][i]);
                            bows .add_white_and_black (cube.number[j][0][i], cube.number[j][lst][i]);
                        };
                    };

                    // Добавляем связей в csp между белыми точками и соседями черных

                    auto csp_add = [&csp] (cst white_point, cst black_point_neighbor){
                        csp .add (white_point, black_point_neighbor);
                        csp .add (black_point_neighbor, white_point);
                    };

                    // Пробегаем ячейку в которой ноходится черная точка и добавляем её соседей как соседей белой точки
                    auto tie_w_point_with_b_cell = [&csp, &cube, &bows, &csp_add] (cst wnx, cst wny, cst bnx, cst bny){
                        // (bnx, bny) - левый нижний угол ячейки
                        // пробегаем все точки ячейки
                        for (st i = bnx; i < bnx + 2; ++i)
                            for (st j = bny; j < bny + 2; ++j)
                                if (bows .subst(cube.number[i][j][0]) != cube.number[wnx][wny][0]) // если это не та же самая точка
                                        for (st k = 0; k < type_space; ++k)
                                            for (st l = 0; l < type_space; ++l)
                                                csp_add (cube.number[wnx][wny][k], bows .subst(cube.number[i][j][l]));
                    };

                    // точки на углах
                    tie_w_point_with_b_cell(0, 0, 0, lst - 1);
                    tie_w_point_with_b_cell(0, 0, lst - 1, 0);
                    tie_w_point_with_b_cell(0, 0, lst - 1, lst - 1);

                    // точки на рёбрах
                    for (st i = 1; i < lst; ++i)
                    {
                        // не забываем, что у каждой точки на границе две соседние ячейки
                        tie_w_point_with_b_cell(0, i, lst - 1, i - 1);
                        tie_w_point_with_b_cell(0, i, lst - 1, i);
                        tie_w_point_with_b_cell(i, 0, i - 1, lst - 1);
                        tie_w_point_with_b_cell(i, 0, i, lst - 1);
                    };

                };

            private:
                static const u8 dim = 2;

        };

    template<u8 type_space>
        class DomainLooperTrivial<3, type_space>
        {
            public:
                DomainLooperTrivial () {};

                void loop_domain (const dealii::DoFHandler<3> &dof_h,
                        OnCell::BlackOnWhiteSubstituter &bows,
                        dealii::DynamicSparsityPattern &csp)
                {
                    puts("i1");
                    CubeOfNumbers<dim, type_space> cube(dof_h.get_tria().n_global_levels()-1);
                    puts("i2");

                    cst lst = cube.size - 1;

                    // Добавляем в bows белые и соответствующие им черные точки
                    for (st i = 0; i < type_space; ++i)
                    {
                        // точки на углах
                        for (st j = 1; j < 8; ++j)
                        {
                            st indx_1 = ((j & 4) >> 2) * lst;
                            st indx_2 = ((j & 2) >> 1) * lst;
                            st indx_3 = (j & 1) * lst;
                            bows .add_white_and_black (cube.number[0][0][0][i], cube.number[indx_1][indx_2][indx_3][i]);
                        };

                        // bows .add_white_and_black (cube.number[0][0][0][i], cube.number[0][0][lst][i]);
                        // bows .add_white_and_black (cube.number[0][0][0][i], cube.number[0][lst][0][i]);
                        // bows .add_white_and_black (cube.number[0][0][0][i], cube.number[0][lst][lst][i]);
                        // bows .add_white_and_black (cube.number[0][0][0][i], cube.number[lst][0][0][i]);
                        // bows .add_white_and_black (cube.number[0][0][0][i], cube.number[lst][0][lst][i]);
                        // bows .add_white_and_black (cube.number[0][0][0][i], cube.number[lst][lst][0][i]);
                        // bows .add_white_and_black (cube.number[0][0][0][i], cube.number[lst][lst][lst][i]);

                        // точки на рёбрах
                        for (st j = 1; j < lst; ++j)
                        {
                            for (st k = 0; k < 4; ++k)
                            {
                                st indx_1 = ((k & 2) >> 1) * lst;
                                st indx_2 = (k & 1) * lst;
                                bows .add_white_and_black (cube.number[0][0][j][i], cube.number[indx_1][indx_2][j][i]);
                                bows .add_white_and_black (cube.number[0][j][0][i], cube.number[indx_1][j][indx_2][i]);
                                bows .add_white_and_black (cube.number[j][0][0][i], cube.number[j][indx_1][indx_2][i]);
                                
                            };
                            // bows .add_white_and_black (cube.number[0][0][j][i], cube.number[0][lst][j][i]);
                            // bows .add_white_and_black (cube.number[0][0][j][i], cube.number[0][lst][j][i]);
                            // bows .add_white_and_black (cube.number[0][0][j][i], cube.number[lst][0][j][i]);
                            // bows .add_white_and_black (cube.number[0][0][j][i], cube.number[lst][lst][j][i]);
                            //
                            // bows .add_white_and_black (cube.number[0][j][0][i], cube.number[0][j][lst][i]);
                            // bows .add_white_and_black (cube.number[0][j][0][i], cube.number[0][j][lst][i]);
                            // bows .add_white_and_black (cube.number[0][j][0][i], cube.number[lst][j][0][i]);
                            // bows .add_white_and_black (cube.number[0][j][0][i], cube.number[lst][j][lst][i]);
                            //
                            // bows .add_white_and_black (cube.number[j][0][0][i], cube.number[j][0][lst][i]);
                            // bows .add_white_and_black (cube.number[j][0][0][i], cube.number[j][0][lst][i]);
                            // bows .add_white_and_black (cube.number[j][0][0][i], cube.number[j][lst][0][i]);
                            // bows .add_white_and_black (cube.number[j][0][0][i], cube.number[j][lst][lst][i]);
                        };

                        // точки на гранях
                        for (st j = 1; j < lst; ++j)
                        {
                            for (st k = 1; k < lst; ++k)
                            {
                                bows .add_white_and_black (cube.number[0][j][k][i], cube.number[lst][j][k][i]);
                                bows .add_white_and_black (cube.number[j][0][k][i], cube.number[j][lst][k][i]);
                                bows .add_white_and_black (cube.number[j][k][0][i], cube.number[j][k][lst][i]);
                            };
                            
                        };
                    };
                    puts("i3");

                    // Добавляем связей в csp между белыми точками и соседями черных

                    auto csp_add = [&csp] (cst white_point, cst black_point_neighbor){
                        csp .add (white_point, black_point_neighbor);
                        csp .add (black_point_neighbor, white_point);
                    };

                    // Пробегаем ячейку в которой ноходится черная точка и добавляем её соседей как соседей белой точки
                    auto tie_w_point_with_b_cell = [&csp, &cube, &bows, &csp_add] (cst wnx, cst wny, cst wnz, cst bnx, cst bny, cst bnz){
                        // (bnx, bny) - левый нижний угол ячейки
                        // пробегаем все точки ячейки
                        for (st i = bnx; i < bnx + 2; ++i)
                            for (st j = bny; j < bny + 2; ++j)
                                for (st k = bnz; k < bnz + 2; ++k)
                                    if (bows .subst(cube.number[i][j][k][0]) != cube.number[wnx][wny][wnz][0]) // если это не та же самая точка
                                        for (st l = 0; l < type_space; ++l)
                                            for (st m = 0; m < type_space; ++m)
                                                csp_add (cube.number[wnx][wny][wnz][l], bows .subst(cube.number[i][j][k][m]));
                    };

                    // точки на углах
                    for (st i = 1; i < 8; ++i)
                    {
                        st indx_1 = ((i & 4) >> 2) * (lst - 1);
                        st indx_2 = ((i & 2) >> 1) * (lst - 1);
                        st indx_3 = (i & 1) * (lst - 1);
                        tie_w_point_with_b_cell(0, 0, 0, indx_1, indx_2, indx_3);
                    };
                    puts("i4");
                    // tie_w_point_with_b_cell(0, 0, 0, 0, 0, lst - 1);
                    // tie_w_point_with_b_cell(0, 0, 0, 0, lst - 1, 0);
                    // tie_w_point_with_b_cell(0, 0, 0, 0, lst - 1, lst - 1);
                    // tie_w_point_with_b_cell(0, 0, 0, lst - 1, 0, 0);
                    // tie_w_point_with_b_cell(0, 0, 0, lst - 1, 0, lst - 1);
                    // tie_w_point_with_b_cell(0, 0, 0, lst - 1, lst - 1, 0);
                    // tie_w_point_with_b_cell(0, 0, 0, lst - 1, lst - 1, lst - 1);

                    // точки на рёбрах
                    // не забываем, что у каждой точки на ребре две соседние ячейки
                    for (st i = 1; i < lst; ++i)
                    {
                        for (st j = 0; j < 4; ++j)
                        {
                            st indx_1 = ((j & 2) >> 1) * (lst - 1);
                            st indx_2 = (j & 1) * (lst - 1);
                            tie_w_point_with_b_cell(0, 0, i, indx_1, indx_2, i - 1);
                            tie_w_point_with_b_cell(0, 0, i, indx_1, indx_2, i);
                            tie_w_point_with_b_cell(0, i, 0, indx_1, i - 1, indx_2);
                            tie_w_point_with_b_cell(0, i, 0, indx_1, i, indx_2);
                            tie_w_point_with_b_cell(i, 0, 0, i - 1, indx_1, indx_2);
                            tie_w_point_with_b_cell(i, 0, 0, i, indx_1, indx_2);
                        };
                    };
                    puts("i5");

                    // точки на гранях
                    for (st i = 1; i < lst; ++i)
                    {
                        for (st j = 1; j < lst; ++j)
                        {
                            // не забываем, что у каждой точки на грани четыре соседние ячейки
                            for (st k = i - 1; k < i + 1; ++k)
                            {
                                for (st l = j - 1; l < j + 1; ++l)
                                {
                                    tie_w_point_with_b_cell(0, i, j, lst - 1, k, l);
                                    tie_w_point_with_b_cell(i, 0, j, k, lst - 1, l);
                                    tie_w_point_with_b_cell(i, j, 0, k, l, lst - 1);
                                };
                            };
                        };
                    };
                    puts("i6");
                };

            private:
                static const u8 dim = 3;

        };
};
#endif /* end of include guard: DOMAIN_LOOPER_TRIVIAL_H_4C8LSHXA */ 
