#ifndef SOURCE_SCALAR_NIKOLA_PROBLEM
#define SOURCE_SCALAR_NIKOLA_PROBLEM
 
#include "../../../../general/source/scalar/source_scalar.h"
#include <deal.II/grid/tria.h>
#include <deal.II/fe/fe.h>
#include <deal.II/base/quadrature_lib.h>

namespace Nikola
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
            typedef std::function<dbl (const dealii::Point<dim>&)> Func;

            SourceScalar (const dealii::FiniteElement<dim> &fe);
            SourceScalar (const vec<arr<Func, dim>> &U, const vec<Func> &tau, const dealii::FiniteElement<dim> &fe);

            virtual void update_on_cell (
                    typename dealii::DoFHandler<dim>::active_cell_iterator &cell) override;

            virtual dbl operator() (cst i) override;

            virtual u8 get_dofs_per_cell () override;


            vec<arr<Func, dim>> U;
            vec<Func> tau;
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
                    dealii::update_values | dealii::update_gradients | dealii::update_quadrature_points | 
                    dealii::update_JxW_values),
            dofs_per_cell (fe.dofs_per_cell),
            num_quad_points (quadrature_formula.size())
    {};

    template <u8 dim>
        SourceScalar<dim>::SourceScalar (
                const vec<arr<Func, dim>> &U_p,
                const vec<Func> &tau_p,
                const dealii::FiniteElement<dim> &fe) :
            SourceScalar<dim>(fe)
    {
        U .resize (U_p.size());
        for (st i = 0; i < U.size(); ++i)
        {
            for (st j = 0; j < dim; ++j)
            {
                U[i][j] = U_p[i][j];
            };
        };

        tau .resize (tau_p.size());
        for (st i = 0; i < tau.size(); ++i)
        {
            tau[i] = tau_p[i];
        };
    };

    template <u8 dim>
        void SourceScalar<dim>::update_on_cell (
                typename dealii::DoFHandler<dim>::active_cell_iterator &cell)
        {
            fe_values .reinit (cell);
            material_id = cell->material_id();
            printf("111 %ld\n", material_id);
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
                    res += 
                        -(this->fe_values.shape_grad (i, q_point)[ort]) *
                        this->U[this->material_id][ort](fe_values.quadrature_point(q_point)) *
                        this->fe_values.JxW(q_point);
                };
                res +=
                    -(this->fe_values.shape_value (i, q_point)) *
                    this->tau[this->material_id](fe_values.quadrature_point(q_point)) *
                    this->fe_values.JxW(q_point);
            };

            return res;
        };

    template <u8 dim>
        u8 SourceScalar<dim>::get_dofs_per_cell ()
        {
            return dofs_per_cell;
        };
};

#endif
