#ifndef PROBLEM_ON_CELL_PRINT
#define PROBLEM_ON_CELL_PRINT 1

#include <stdio.h>
#include <stdlib.h>
#include <fstream>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/numerics/data_out.h>

// #include "../calculate_meta_coefficients/calculate_meta_coefficients.h"
#include "../stress_calculator/stress_calculator.h"
#include "../deform_calculator/deform_calculator.h"

namespace OnCell
{
    namespace BDStore
    {
        template <st dim>
            void print_coor (const dealii::DoFHandler<dim> &dof_handler, const str path)
            {
                const str file_name = path + "/coor_";
                vec<dealii::Point<dim>> coor(dof_handler.n_dofs() / dim);
                {
                    cu8 dofs_per_cell =  dof_handler.get_fe().dofs_per_cell;

                    std::vector<u32> local_dof_indices (dofs_per_cell);

                    auto cell = dof_handler.begin_active();
                    auto endc = dof_handler.end();
                    for (; cell != endc; ++cell)
                    {
                        cell ->get_dof_indices (local_dof_indices);

                        FOR (i, 0, dofs_per_cell / dim)
                        {
                            coor[local_dof_indices[i*dim] / dim] = cell ->vertex (i);
                        };
                    };
                };

                {
                    arr<str, 3> ort = {"x", "y", "z"};
                    for (st i = 0; i < dim; ++i)
                    {
                        std::ofstream out (file_name + ort[i] + ".bin", std::ios::out | std::ios::binary);
                        for (st j = 0; j < coor.size(); ++j)
                        {
                            out.write ((char *) &(coor[j](i)), sizeof(dbl));
                        };
                        out.close ();
                    };
                };
                {
                    std::ofstream f(path + "/size_solution.bin", std::ios::out | std::ios::binary);
                    cst size = dof_handler.n_dofs() / dim;
                    f.write ((char *) &(size), sizeof(st));
                    f.close ();
                };
            };

        void print_meta_coef(
                const ATools::FourthOrderTensor& meta_coef,
                const str path)
        {
            std::ofstream out (path + "/meta_coef.bin", std::ios::out | std::ios::binary);
            out.write ((char *) &meta_coef, sizeof meta_coef);
            out.close ();
        };

        void print_move (
                const dealii::DoFHandler<3> &dof_handler,
                const ArrayWithAccessToVector<arr<dealii::Vector<dbl>, 3>>& move,
                cst number_of_approx,
                const str path)
        {
            const arr<str, 3> ort = {"x", "y", "z"};
            const arr<str, 6> aprx = {"0", "1", "2", "3", "4", "5"};

            for (st approx_number = 1; approx_number < number_of_approx; ++approx_number)
            {
                const str move_path = path + aprx[approx_number] + "/";
                mkdir(move_path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
                for (i32 i = 0; i < approx_number+1; ++i)
                {
                    for (i32 j = 0; j < approx_number+1; ++j)
                    {
                        for (i32 k = 0; k < approx_number+1; ++k)
                        {
                            if ((i+j+k) == approx_number)
                            {
                                arr<i32, 3> approximation = {i, j, k};
                                const str move_path_approx = move_path + 
                                    aprx[i] + "_" +
                                    aprx[j] + "_" +
                                    aprx[k] + "/";
                                mkdir(move_path_approx.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
                                for (st nu = 0; nu < 3; ++nu)
                                {
                                    const str move_path_approx_nu = move_path_approx + ort[nu] + "/"; 
                                    mkdir(move_path_approx_nu.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
                                    const str move_path_final = move_path_approx_nu + "move/"; 
                                    mkdir(move_path_final.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
                                    for (st alpha = 0; alpha < 3; ++alpha)
                                    {
                                        vec<dbl> move_ort(move[approximation][nu].size() / 3);
                                        for (st l = 0; l < move[approximation][nu].size(); ++l)
                                        {
                                            if ((l % 3) == alpha)
                                            {
                                                move_ort[l / 3] = move[approximation][nu][l];
                                            };
                                        };

                                        std::ofstream out (
                                                move_path_final +
                                                ort[alpha] +".bin",
                                                std::ios::out | std::ios::binary);
                                        for (st l = 0; l < move_ort.size(); ++l)
                                        {
                                            out.write ((char *) &(move_ort[l]), sizeof(dbl));
                                        };
                                        out.close ();
                                    };
                                };
                            };
                        };
                    };
                };
            };
        };

        void print_stress(
                const dealii::DoFHandler<3> &dof_handler,
                const vec<ATools::FourthOrderTensor> &coefficient,
                const dealii::FiniteElement<3> &fe,
                const ArrayWithAccessToVector<arr<dealii::Vector<dbl>, 3>>& move,
                cst number_of_approx,
                const str path)
        {
            StressCalculatorBD stress_calculator (dof_handler, coefficient, fe);
            dealii::Vector<dbl> stress(dof_handler.n_dofs() / 3);

            const arr<str, 3> ort = {"x", "y", "z"};
            const arr<str, 6> aprx = {"0", "1", "2", "3", "4", "5"};

            for (st approx_number = 1; approx_number < number_of_approx; ++approx_number)
            {
                const str stress_path = path + aprx[approx_number] + "/";
                mkdir(stress_path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
                for (i32 i = 0; i < approx_number+1; ++i)
                {
                    for (i32 j = 0; j < approx_number+1; ++j)
                    {
                        for (i32 k = 0; k < approx_number+1; ++k)
                        {
                            if ((i+j+k) == approx_number)
                            {
                                arr<i32, 3> approximation = {i, j, k};
                                const str stress_path_approx = stress_path + 
                                    aprx[i] + "_" +
                                    aprx[j] + "_" +
                                    aprx[k] + "/";
                                mkdir(stress_path_approx.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
                                for (st nu = 0; nu < 3; ++nu)
                                {
                                    const str stress_path_approx_nu = stress_path_approx + ort[nu] + "/"; 
                                    mkdir(stress_path_approx_nu.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
                                    const str stress_path_final = stress_path_approx_nu + "stress/"; 
                                    mkdir(stress_path_final.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
                                    for (st alpha = 0; alpha < 3; ++alpha)
                                    {
                                        for (st beta = 0; beta < 3; ++beta)
                                        {
                                            for (st l = 0; l < stress.size(); ++l)
                                            {
                                                stress[l] = 0.0;
                                            };

                                            stress_calculator .calculate (
                                                    approximation, nu, alpha, beta,
                                                    move, stress);

                                            std::ofstream out (
                                                    stress_path_final +
                                                    ort[alpha] + ort[beta] + ".bin",
                                                    std::ios::out | std::ios::binary);
                                            for (st l = 0; l < stress.size(); ++l)
                                            {
                                                out.write ((char *) &(stress[l]), sizeof(dbl));
                                            };
                                            out.close ();
                                        };
                                    };
                                };
                            };
                        };
                    };
                };
            };
        };

        void print_deform(
                const dealii::DoFHandler<3> &dof_handler,
                const dealii::FiniteElement<3> &fe,
                const ArrayWithAccessToVector<arr<dealii::Vector<dbl>, 3>> &move,
                cst number_of_approx,
                const str path)
        {
            DeformCalculatorBD deform_calculator (dof_handler, fe);
            dealii::Vector<dbl> deform(dof_handler.n_dofs() / 3);

            const arr<str, 3> ort = {"x", "y", "z"};
            const arr<str, 6> aprx = {"0", "1", "2", "3", "4", "5"};

            for (st approx_number = 1; approx_number < number_of_approx; ++approx_number)
            {
                const str deform_path = path + aprx[approx_number] + "/";
                mkdir(deform_path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
                for (i32 i = 0; i < approx_number+1; ++i)
                {
                    for (i32 j = 0; j < approx_number+1; ++j)
                    {
                        for (i32 k = 0; k < approx_number+1; ++k)
                        {
                            if ((i+j+k) == approx_number)
                            {
                                arr<i32, 3> approximation = {i, j, k};
                                const str deform_path_approx = deform_path + 
                                    aprx[i] + "_" +
                                    aprx[j] + "_" +
                                    aprx[k] + "/";
                                mkdir(deform_path_approx.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
                                for (st nu = 0; nu < 3; ++nu)
                                {
                                    const str deform_path_approx_nu = deform_path_approx + ort[nu] + "/"; 
                                    mkdir(deform_path_approx_nu.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
                                    const str deform_path_final = deform_path_approx_nu + "deform/"; 
                                    mkdir(deform_path_final.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
                                    for (st alpha = 0; alpha < 3; ++alpha)
                                    {
                                        for (st beta = 0; beta < 3; ++beta)
                                        {
                                            for (st l = 0; l < deform.size(); ++l)
                                            {
                                                deform[l] = 0.0;
                                            };

                                            deform_calculator .calculate (
                                                    approximation, nu, alpha, beta,
                                                    move, deform);

                                            std::ofstream out (
                                                    deform_path_final +
                                                    ort[alpha] + ort[beta] + ".bin",
                                                    std::ios::out | std::ios::binary);

                                            for (st l = 0; l < deform.size(); ++l)
                                            {
                                                out.write ((char *) &(deform[l]), sizeof(dbl));
                                            };
                                            out.close ();
                                        };
                                    };
                                };
                            };
                        };
                    };
                };
            };
        };
    }
}
#endif
