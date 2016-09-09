#ifndef FEATURE_SOURCE_Vector_PROBLEM
#define FEATURE_SOURCE_Vector_PROBLEM
 
#include "../scalar/source_scalar.h"

// namespace Nikola
// {
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
        class SourceVectorFeature : public ::SourceInterface<dim>
    {
        public:
            // typedef std::function<dbl (const dealii::Point<dim>&)> Func;
            typedef typename SourceScalarFeature<dim>::Func Func;

            SourceVectorFeature (const dealii::FiniteElement<dim> &fe);
            SourceVectorFeature (const vec<arr<arr<Func, dim>, dim>> &U_p, const dealii::FiniteElement<dim> &fe);

            virtual void update_on_cell (
                    typename dealii::DoFHandler<dim>::active_cell_iterator &cell) override;

            virtual dbl operator() (cst i) override;

            virtual u8 get_dofs_per_cell () override;


            vec<arr<arr<Func, dim>, dim>> U;
            SourceScalarFeature<dim> source; 
    };

    template <u8 dim>
        SourceVectorFeature<dim>::SourceVectorFeature (const dealii::FiniteElement<dim> &fe) :
            source (fe)
    {};

    template <u8 dim>
        SourceVectorFeature<dim>::SourceVectorFeature (
                const vec<arr<arr<Func, dim>, dim>> &U_p,
                const dealii::FiniteElement<dim> &fe) :
            SourceVectorFeature<dim>(fe)
    {
        U .resize (U_p.size());
        for (st i = 0; i < U.size(); ++i)
        {
            for (st j = 0; j < dim; ++j)
            {
                for (st k = 0; k < dim; ++k)
                {
                    U[i][j][k] = U_p[i][j][k];
                };
            };
        };

        source.U .resize (U.size());
    };

    template <u8 dim>
        void SourceVectorFeature<dim>::update_on_cell (
                typename dealii::DoFHandler<dim>::active_cell_iterator &cell)
        {
            source .update_on_cell (cell);
        };

    template <u8 dim>
        dbl SourceVectorFeature<dim>::operator () (cst i)
        {
            for (st m = 0; m < U.size(); ++m)
            {
                for (st n = 0; n < dim; ++n)
                {
                    source.U[m][n] = U[m][i % dim][n];
                };
            };

            return source (i);
        };

    template <u8 dim>
        u8 SourceVectorFeature<dim>::get_dofs_per_cell ()
        {
            return source.dofs_per_cell;
        };
// };

#endif
