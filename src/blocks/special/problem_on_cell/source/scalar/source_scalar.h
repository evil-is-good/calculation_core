#ifndef SOURCE_SCALAR_ON_CELL
#define SOURCE_SCALAR_ON_CELL
 
#include "../../../../general/source/scalar/source_scalar.h"
#include "../../problem_on_cell_tools/problem_on_cell_tools.h"
#include <deal.II/grid/tria.h>
#include <deal.II/fe/fe.h>
#include <deal.II/base/quadrature_lib.h>

namespace OnCell
{
//! Элемент вектора праой части уравнения в МКЭ для задачи на ячейке, скалярный случай
/*!
 * Правая часть уравнений на ячейке имеет вид:
 \f[
 \sum_{n=0}^{2}\frac{\partial}{\partial x_n}\lambda_{x_n\nu}; \nu \in \{0, 1, 2\}
 \f]
 где \f$\lambda\f$ - коэфффициент физических свойств материалов входящих 
 в состав композита (например тепропроводность) 
 После интегрированиия по ячейке
 \f[
 \sum_{n=0}^{2}\int_{cell}\frac{\partial}{\partial x_n}\lambda_{x_n\nu} \phi_i; \nu \in \{0, 1, 2\}
 \f]
 Так как производная от прерывной функции бесконечна, то интегрируем по частям чтобы избавится 
 от сингулярности.
 \f[
 \sum_{n=0}^{2}\left(\int_{cell}\frac{\partial}{\partial x_n}(\lambda_{x_n\nu} \phi_i) - \int_{cell}\lambda_{x_n\nu} \frac{\partial}{\partial x_n}\phi_i\right); 
 \nu \in \{0, 1, 2\}
 \f]
 Из-за цикличности границ ячейки первый интеграл равен нулю.
 \f[
 - \sum_{n=0}^{2}\int_{cell}\lambda_{x_n\nu}(q) \frac{\partial}{\partial x_n}\phi_i(q); 
 \nu \in \{0, 1, 2\}
 \f]
  Интеграл расчитывается в квадратурах.
  \f[
  - \sum_{q=0}^{numq} \sum_{n=0}^{2}\int_{cell}\lambda_{x_n\nu}(q) \frac{\partial}{\partial x_n}\phi_i(q) J(q)W(q); \nu \in \{0, 1, 2\}
  \f]
*/
    template <u8 dim>
        class SourceScalar : public ::SourceInterface<dim>
    {
        public:

            SourceScalar (const dealii::FiniteElement<dim> &fe);
            SourceScalar (const vec<arr<dbl, dim>> &coefficient, const dealii::FiniteElement<dim> &fe);

            virtual void update_on_cell (
                    typename dealii::DoFHandler<dim>::active_cell_iterator &cell) override;

            virtual dbl operator() (cst i) override;

            virtual u8 get_dofs_per_cell () override;


            vec<arr<dbl, dim>> coef;
            dealii::QGauss<dim>        quadrature_formula; //!< Формула интегрирования в квадратурах.
            dealii::FEValues<dim, dim> fe_values; //!< Тип функций формы.
            cu8                        dofs_per_cell; //!< Количество узлов в ячейке (зависит от типа функций формы).
            cu8                        num_quad_points; //!< Количество точек по которым считается квадратура.
            u8                         material_id = 0; //!< Идентефикатор материала ячейки.
    };

    template <u8 dim>
        SourceScalar<dim>::SourceScalar (const dealii::FiniteElement<dim> &fe) :
            quadrature_formula (2),
            fe_values (fe, quadrature_formula,
                    dealii::update_gradients | dealii::update_quadrature_points | 
                    dealii::update_JxW_values),
            dofs_per_cell (fe.dofs_per_cell),
            num_quad_points (quadrature_formula.size())
    {};

    template <u8 dim>
        SourceScalar<dim>::SourceScalar (
                const vec<arr<dbl, dim>> &coefficient,
                const dealii::FiniteElement<dim> &fe) :
            SourceScalar<dim>(fe)
    {
        coef .resize (coefficient.size());
        for (st i = 0; i < coefficient.size(); ++i)
        {
            for (st j = 0; j < dim; ++j)
            {
                coef[i][j] = coefficient[i][j];
            };
        };
    };

    template <u8 dim>
        void SourceScalar<dim>::update_on_cell (
                typename dealii::DoFHandler<dim>::active_cell_iterator &cell)
        {
            fe_values .reinit (cell);
            material_id = cell->material_id();
        };

    template <u8 dim>
        dbl SourceScalar<dim>::operator () (cst i)
        {
            const uint8_t num_quad_points = this->quadrature_formula.size();

            dbl res = 0.0;

            for (st q_point = 0; q_point < num_quad_points; ++q_point)
            {
                for(st ort = 0; ort < dim; ++ort)
                {
        // res +=  
        //     fe_values.shape_value (i, q_point) *
        //     fe_values.JxW(q_point);
                        // fe_values.shape_grad (0, 0);// *
                    res += 
                        -(this->fe_values.shape_grad (i, q_point)[ort]) *
                        coef[this->material_id][ort] *
                        this->fe_values.JxW(q_point);
                // res +=  
                //     fe_values.shape_grad (i, q_point)[0] *
                //     fe_values.shape_grad (i, q_point)[0] *
                //     fe_values.JxW(q_point);
                };
            };

            return res;
        };

    template <u8 dim>
        u8 SourceScalar<dim>::get_dofs_per_cell ()
        {
            return dofs_per_cell;
        };




    template <u8 dim>
        class SourceScalarApprox : public ::SourceInterface<dim>
    {
        public:

            SourceScalarApprox (const dealii::FiniteElement<dim> &fe);
            SourceScalarApprox (const arr<i32, 3> approximation,
                    const vec<arr<arr<dbl, 3>, 3>> &coefficient, 
                    const ArrayWithAccessToVector<dbl> &meta_coefficient,
                    // ArrayWithAccessToVector<dealii::Vector<dbl>> *psi_func,
                    ArrayWithAccessToVector<dealii::Vector<dbl>> &psi_func,
                    const dealii::FiniteElement<dim> &fe);

            virtual void update_on_cell (
                    typename dealii::DoFHandler<dim>::active_cell_iterator &cell) override;

            virtual dbl operator() (cst i) override;

            virtual u8 get_dofs_per_cell () override;

            static cst x = 0;
            static cst y = 1;
            static cst z = 2;

            vec<arr<arr<dbl, 3>, 3>> coef; //!< Коэффициенты материла, например теплопроводность.
            ArrayWithAccessToVector<dbl> meta_coef; //!< Массив коэффициентов метаматериала.
            // ArrayWithAccessToVector<dealii::Vector<dbl>> *psi; //!< Массив ячейкивых функций.
            ArrayWithAccessToVector<dealii::Vector<dbl>> psi; //!< Массив ячейкивых функций.
            arr<i32, 3> k; //!< Номер приближения в векторной форме записи.
            dealii::QGauss<dim>        quadrature_formula; //!< Формула интегрирования в квадратурах.
            dealii::FEValues<dim, dim> fe_values; //!< Тип функций формы.
            vec<u32>                   global_dof_indices; //!< Список глобальных значений индексов с ячейки.
            cu8                        dofs_per_cell; //!< Количество узлов в ячейке (зависит от типа функций формы).
            cu8                        num_quad_points; //!< Количество точек по которым считается квадратура.
            u8                         material_id = 0; //!< Идентефикатор материала ячейки.

            dbl tmp = 0.0;
            arr<st, 4> indx;

    };

    template <u8 dim>
        SourceScalarApprox<dim>::SourceScalarApprox (const dealii::FiniteElement<dim> &fe) :
            quadrature_formula (2),
            fe_values (fe, quadrature_formula,
                    dealii::update_values | dealii::update_gradients | dealii::update_quadrature_points | 
                    dealii::update_JxW_values),
            dofs_per_cell (fe.dofs_per_cell),
            global_dof_indices (fe.dofs_per_cell),
            num_quad_points (quadrature_formula.size())
    {};

    template <u8 dim>
        SourceScalarApprox<dim>::SourceScalarApprox (
                    const arr<i32, 3> approximation,
                    const vec<arr<arr<dbl, 3>, 3>> &coefficient, 
                    const ArrayWithAccessToVector<dbl> &meta_coefficient,
                    // ArrayWithAccessToVector<dealii::Vector<dbl>> *psi_func,
                    ArrayWithAccessToVector<dealii::Vector<dbl>> &psi_func,
                    const dealii::FiniteElement<dim> &fe) : 
            SourceScalarApprox<dim>(fe)
               {
        k[x] = approximation[x];
        k[y] = approximation[y];
        k[z] = approximation[z];

        coef .resize (coefficient.size());
        for (st i = 0; i < coefficient.size(); ++i)
        {
            for (st j = 0; j < dim; ++j)
            {
                for (st k = 0; k < dim; ++k)
                {
                    coef[i][j][k] = coefficient[i][j][k];
                };
            };
        };

        meta_coef .reinit (meta_coefficient.size);
        meta_coef = meta_coefficient;

        // cst size = meta_coefficient.size();
        // meta_coef .resize (size);
        // for (st i = 0; i < size; ++i)
        // {
        //     meta_coef[i] .resize (size);
        //     for (st j = 0; j < size; ++j)
        //     {
        //         meta_coef[i][j] .resize (size);
        //         for (st k = 0; k < size; ++k)
        //         {
        //             meta_coef[i][j][k]  = meta_coefficient[i][j][k];
        //
        //         };
        //     };
        // };

        // psi = psi_func;
        psi.reinit(psi_func.content.size());
        // psi.content.resize(psi_func.content.size());
        for (st i = 0; i < psi_func.content.size(); ++i)
        {
            // psi.content[i].resize(psi_func.content[i].size());
            for (st j = 0; j < psi_func.content[i].size(); ++j)
            {
                // psi.content[i][j].resize(psi_func.content[i][j].size());
                for (st k = 0; k < psi_func.content[i][j].size(); ++k)
                {
                    psi.content[i][j][k].reinit(psi_func.content[i][j][k].size());
                    for (st l = 0; l < psi_func.content[i][j][k].size(); ++l)
                    {
                        psi.content[i][j][k][l] = psi_func.content[i][j][k][l];
                    };
                };
            };
        };

    };

    template <u8 dim>
        void SourceScalarApprox<dim>::update_on_cell (
                typename dealii::DoFHandler<dim>::active_cell_iterator &cell)
        {
            fe_values .reinit (cell);
            material_id = cell->material_id();
            cell->get_dof_indices (global_dof_indices);
            indx[0] = cell->vertex_dof_index(0, 0);
            indx[1] = cell->vertex_dof_index(1, 0);
            indx[2] = cell->vertex_dof_index(2, 0);
            indx[3] = cell->vertex_dof_index(3, 0);
        };

    template <u8 dim>
        dbl SourceScalarApprox<dim>::operator () (cst i)
        {
            const uint8_t num_quad_points = this->quadrature_formula.size();

            dbl res = 0.0;
            dbl summ = 0.0;

            dbl f1 = 0.0;
            dbl f2 = 0.0;
            dbl f3 = 0.0;

            for (st q_point = 0; q_point < num_quad_points; ++q_point)
            {
                for (st l = 0; l < dim; ++l)
                {
                    for (st m = 0; m < dim; ++m)
                    {
                        // arr<i32, 3> km = {k[x], k[y], k[z]}; // k - э_m
                        // if (k[m] > 0)
                        // {
                        //     km[m]--;
                        //     for (st n = 0; n < dofs_per_cell; ++n)
                        //     {
                        //         cdbl psi_km = psi.content.at(km[0]).at(km[1]).at(km[2])[global_dof_indices[n]];
                        //         f1 += 
                        //             this->coef[this->material_id][l][m] * 
                        //             psi_km *
                        //             // this->psi[0][km][global_dof_indices[n]] * 
                        //             this->fe_values.shape_value (n, q_point) *
                        //             this->fe_values.shape_grad (i, q_point)[l] *
                        //             this->fe_values.JxW(q_point);
                        //     };
                        //     arr<i32, 3> klm = {km[x], km[y], km[z]}; // k - э_l - э_m
                        //     if (klm[l] > 0)
                        //     {
                        //         klm[l]--;
                        //         for (st n = 0; n < dofs_per_cell; ++n)
                        //         {
                        //             cdbl psi_km = psi.content.at(km[0]).at(km[1]).at(km[2])[global_dof_indices[n]];
                        //             cdbl psi_klm = psi.content.at(klm[0]).at(klm[1]).at(klm[2])[global_dof_indices[n]];
                        //             f3 += 
                        //                 (this->coef[this->material_id][l][m] *
                        //                  (psi_km *
                        //                   this->fe_values.shape_grad (n, q_point)[m] +
                        //                   psi_klm *
                        //                   this->fe_values.shape_value (n, q_point)
                        //                  )
                        //                 ) *
                        //                 this->fe_values.shape_value (i, q_point) *
                        //                 this->fe_values.JxW(q_point);
                        //         };
                        //     };
                        // };
                        for (st n = 0; n < dofs_per_cell; ++n)
                        {
                            arr<i32, 3> km = {k[x], k[y], k[z]}; // k - э_m
                            bool psi_km_exist = true;
                            if (k[m] == 0)
                                psi_km_exist = false;
                            else
                                km[m]--;

                            arr<i32, 3> kl = {k[x], k[y], k[z]}; // k - э_l
                            bool psi_kl_exist = true;
                            if (k[l] == 0)
                                psi_kl_exist = false;
                            else
                                kl[l]--;

                            bool psi_klm_exist = true;
                            arr<i32, 3> klm = {km[x], km[y], km[z]}; // k - э_l - э_m
                            if (km[l] == 0)
                                psi_klm_exist = false;
                            else
                                klm[l]--;
                            // printf("k[m]=%d km=%d k[l]=%d klm=%d %d %d\n",
                            //         k[m], km[m], k[l], klm[l], psi_km_exist, psi_klm_exist);
                            // printf("size %d %d %f\n", psi[0][km].size(), global_dof_indices[n], psi[0][km] != dealii::Vector<dbl>(0) ? psi[0][km][global_dof_indices[n]] : 0.0);
                            dbl psi_km = 0.0;
                            dbl psi_kl = 0.0;
                            dbl psi_klm = 0.0;
                            // printf("%d %d %d %f\n", km[x], km[y], km[z], psi[0][km][0]);
                            // if (psi[0][km] == dealii::Vector<dbl>(1))
                            // if ((km[x] > -1) and (km[y] > -1) and (km[z] > -1))
                            if (psi_km_exist)
                                psi_km = psi[km][global_dof_indices[n]];
                            else
                                psi_km = 0.0;

                            if (psi_kl_exist)
                                psi_kl = psi[kl][global_dof_indices[n]];
                            else
                                psi_kl = 0.0;
                                // psi_km = psi[0][km][indx[n]];
                            // if (psi[0][klm] == dealii::Vector<dbl>(1))
                            // if ((klm[x] > -1) and (klm[y] > -1) and (klm[z] > -1))
                            if (psi_klm_exist)
                                psi_klm = psi[klm][global_dof_indices[n]];
                            else
                                psi_klm = 0.0;
                                // psi_klm = psi[0][klm][indx[n]];

                            // cdbl psi_km = psi[0][km] != dealii::Vector<dbl>(0) ? psi[0][km][global_dof_indices[n]] : 0.0;
                            // cdbl psi_klm = psi[0][klm] != dealii::Vector<dbl>(0) ? psi[0][klm][global_dof_indices[n]] : 0.0;
                            f1 += 
                                this->coef[this->material_id][l][m] * 
                                psi_km *
                                // this->psi[0][km][global_dof_indices[n]] * 
                                this->fe_values.shape_value (n, q_point) *
                                this->fe_values.shape_grad (i, q_point)[l] *
                                this->fe_values.JxW(q_point);
                            // cdbl f2 = 
                            //     (this->meta_coef[k] -
                            //      this->coef[this->material_id][l][m] *
                            //      // (this->psi[0][km][global_dof_indices[n]] *
                            //      (psi_km *
                            //       this->fe_values.shape_grad (n, q_point)[m] +
                            //       // this->psi[0][klm][global_dof_indices[n]] *
                            //       psi_klm *
                            //       this->fe_values.shape_value (n, q_point)
                            //      )
                            //     ) *
                            //     this->fe_values.shape_value (i, q_point) *
                            //     this->fe_values.JxW(q_point);
                            f3 += 
                                 (this->coef[this->material_id][l][m] *
                                 (psi_kl *
                                  this->fe_values.shape_grad (n, q_point)[m] +
                                  psi_klm *
                                  this->fe_values.shape_value (n, q_point)
                                 )
                                ) *
                                this->fe_values.shape_value (i, q_point) *
                                this->fe_values.JxW(q_point);

                                // if (k[x] == 2)
                                // res += 
                                // - f1;
                                // // - (f1 + f2);
                                // else
                               //  res += 
                               //  // - f2;
                               //  // - (f1 + f2);
                               //  // f2 - f1;
                               // -f1-f2-f3;
                               //  tmp += f2 - f3;//-f3;
                               //  summ += f3;
                        };
                    };
                };
                f2 += 1.0 *(
                    this->meta_coef[k] *
                    this->fe_values.shape_value (i, q_point) *
                    this->fe_values.JxW(q_point));
                // tmp += res;
            };
            tmp += f3-f2;
            // tmp += f2-f3;
            // printf("%f %f %f %f %ld\n", f2, f3, tmp, this->meta_coef[k], k[0]);
        //         {
        // // res +=  
        // //     fe_values.shape_value (i, q_point) *
        // //     fe_values.JxW(q_point);
    //                 // fe_values.shape_grad (0, 0);// *
        //             // res += 
        //             //     -(this->fe_values.shape_grad (i, q_point)[ort]) *
        //             //     coef[this->material_id][ort] *
        //             //     this->fe_values.JxW(q_point);
        //         // res +=  
        //         //     fe_values.shape_grad (i, q_point)[0] *
        //         //     fe_values.shape_grad (i, q_point)[0] *
        //         //     fe_values.JxW(q_point);
        //         };
        //     };

            return f3-f2-f1;//-(f1 + f3 - f2);
        };

    template <u8 dim>
        u8 SourceScalarApprox<dim>::get_dofs_per_cell ()
        {
            return dofs_per_cell;
        };
};

#endif
