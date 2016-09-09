#ifndef POLY_MATERIALS_SOURCE_SCALAR
#define POLY_MATERIALS_SOURCE_SCALAR
 
#include "../../../general/source/scalar/source_scalar.h"
#include <deal.II/grid/tria.h>
#include <deal.II/fe/fe.h>
#include <deal.II/base/quadrature_lib.h>

//! Элемент вектора праой части уравнения в МКЭ, скалярный случай
/*!
 * Представляет собой интеграл по ячейке:
  \f[
  \int_{cell} f \phi_i 
  \f]
  Интеграл расчитывается в квадратурах.
  \f[
  \sum_{q=0}^{numq} f(q)\phi_i(q)  J(q)W(q)
  \f]
*/
template <u8 dim>
class SourceScalarPolyMaterials : public SourceInterface<dim>
{
    public:
    typedef std::function<dbl (const dealii::Point<dim>&)> Func;

    SourceScalarPolyMaterials (const dealii::FiniteElement<dim> &fe);
    SourceScalarPolyMaterials (const vec<Func> &func, const dealii::FiniteElement<dim> &fe);

    virtual void update_on_cell (
        typename dealii::DoFHandler<dim>::active_cell_iterator &cell) override;

    virtual dbl operator() (cst i) override;

    virtual u8 get_dofs_per_cell () override;


    vec<Func>                  f; //!< Функция из правой части уравнения.
    dealii::QGauss<dim>        quadrature_formula; //!< Формула интегрирования в квадратурах.
    dealii::FEValues<dim, dim> fe_values; //!< Тип функций формы.
    cu8                        dofs_per_cell; //!< Количество узлов в ячейке (зависит от типа функций формы).
    cu8                        num_quad_points; //!< Количество точек по которым считается квадратура.
    u8                         material_id = 0; //!< Идентефикатор материала ячейки.

};

template <u8 dim>
SourceScalarPolyMaterials<dim>::SourceScalarPolyMaterials (const dealii::FiniteElement<dim> &fe) :
    quadrature_formula (2),
    fe_values (fe, quadrature_formula,
            dealii::update_values | dealii::update_quadrature_points | 
            dealii::update_JxW_values),
    dofs_per_cell (fe.dofs_per_cell),
    num_quad_points (quadrature_formula.size())
{};

template <u8 dim>
SourceScalarPolyMaterials<dim>::SourceScalarPolyMaterials (const vec<Func> &func, const dealii::FiniteElement<dim> &fe) :
    SourceScalarPolyMaterials(fe)
{
    f .resize (func.size());
    for (st i = 0; i < f.size(); ++i)
    {
        f[i] = func[i];
    };
};

template <u8 dim>
void SourceScalarPolyMaterials<dim>::update_on_cell (
        typename dealii::DoFHandler<dim>::active_cell_iterator &cell)
{
    fe_values .reinit (cell);
    material_id = cell->material_id();
    // printf("%в\n", material_id);
};

template <u8 dim>
dbl SourceScalarPolyMaterials<dim>::operator () (cst i)
{
    dbl res = 0.0;

    FOR (q_point, 0, num_quad_points)
    {
        res -=  
            this->fe_values.shape_value (i, q_point) *
            this->f[material_id](fe_values.quadrature_point(q_point)) *
            this->fe_values.JxW(q_point);
    };

    return res;
};

template <u8 dim>
u8 SourceScalarPolyMaterials<dim>::get_dofs_per_cell ()
{
    return dofs_per_cell;
};

#endif
