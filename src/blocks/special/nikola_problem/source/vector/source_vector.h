#ifndef SOURCE_Vector_NIKOLA_PROBLEM
#define SOURCE_Vector_NIKOLA_PROBLEM
 
#include "../scalar/source_scalar.h"

namespace Nikola
{
//! Элемент вектора праой части уравнения в МКЭ для задачи Николы, векторный случай
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
        class SourceVector : public ::SourceInterface<dim>
    {
        public:
            // typedef std::function<dbl (const dealii::Point<dim>&)> Func;
            typedef typename ::Nikola::SourceScalar<dim>::Func Func;

            SourceVector (const dealii::FiniteElement<dim> &fe);
            SourceVector (const vec<arr<Func, dim>> &U, const vec<arr<Func, dim>> &tau, const dealii::FiniteElement<dim> &fe);

            virtual void update_on_cell (
                    typename dealii::DoFHandler<dim>::active_cell_iterator &cell) override;

            virtual dbl operator() (cst i) override;

            virtual u8 get_dofs_per_cell () override;


            vec<arr<Func, dim>> U;
            vec<arr<Func, dim>> tau;
            ::Nikola::SourceScalar<dim> source; 
    };

    template <u8 dim>
        SourceVector<dim>::SourceVector (const dealii::FiniteElement<dim> &fe) :
            source (fe)
    {};

    template <u8 dim>
        SourceVector<dim>::SourceVector (
                const vec<arr<Func, dim>> &U_p,
                const vec<arr<Func, dim>> &tau_p,
                const dealii::FiniteElement<dim> &fe) :
            SourceVector<dim>(fe)
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
        for (st i = 0; i < U.size(); ++i)
        {
            for (st j = 0; j < dim; ++j)
            {
                tau[i][j] = tau_p[i][j];
            };
        };
            source.U .resize (U.size());
            source.tau .resize (tau.size());
    };

    template <u8 dim>
        void SourceVector<dim>::update_on_cell (
                typename dealii::DoFHandler<dim>::active_cell_iterator &cell)
        {
            source .update_on_cell (cell);
        };

    template <u8 dim>
        dbl SourceVector<dim>::operator () (cst i)
        {
            for (st m = 0; m < U.size(); ++m)
            {
                for (st n = 0; n < dim; ++n)
                {
                    source.U[m][n] = (i % dim) == n ? U[m][n] : [](const dealii::Point<2> &p){return 0.0;};
                };

                source.tau[m] = tau[m][i % dim];
            };

            return source (i);
        };

    template <u8 dim>
        u8 SourceVector<dim>::get_dofs_per_cell ()
        {
            return source.dofs_per_cell;
        };
};

#endif
