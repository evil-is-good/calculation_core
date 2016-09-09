#ifndef SOURCE_VECTOR
#define SOURCE_VECTOR
 
#include "../scalar/source_scalar.h"

//! Элемент вектора праой части уравнения в МКЭ, векторный случай
/*!
 * Представляет собой интеграл по ячейке:
  \f[
  \int_{cell} f_{(i \mod dim)} \phi_i 
  \f]
  Так как функции векторные, то в интеграле используются компоненты функций. В диле компоненты функций формы
  записываются по такому порядку (двумерный случай): \phi_0 - компонет x нулевой функции формы,
  \phi_1 - компонет y нулевой функции формы, \phi_2 - компонет x первой функции формы и тд.
  Интеграл расчитывается в квадратурах.
  \f[
  \sum_{q=0}^{numq} f_{(i \mod dim)}(q)\phi_i(q)  J(q)W(q)
  \f]
*/
template <u8 dim>
class SourceVector : public SourceInterface<dim>
{
    public:
    typedef std::function<dbl (const dealii::Point<dim>&)> Func;

    SourceVector (
            const arr<Func, dim> func, const dealii::FiniteElement<dim> &fe);

    virtual dbl operator() (cst i) override;

    virtual void update_on_cell (
        typename dealii::DoFHandler<dim>::active_cell_iterator &cell) override;

    virtual u8 get_dofs_per_cell () override;


    arr<Func, dim>    f; //!< Функция из правой части уравнения.
    SourceScalar<dim> source; 
};

template <u8 dim>
SourceVector<dim>::SourceVector (
        const arr<Func, dim> func, const dealii::FiniteElement<dim> &fe) :
    f(func), source (fe)
{

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
    source.f = f[i % dim];

    return source (i);
};

template <u8 dim>
u8 SourceVector<dim>::get_dofs_per_cell ()
{
    return source.get_dofs_per_cell();
};

#endif
