#ifndef LAPLACIAN_SCALAR
#define LAPLACIAN_SCALAR
 
#include "../interface/laplacian_interface.h"
#include "../../additional_tools/types/types.h"
#include <deal.II/grid/tria.h>
#include <deal.II/fe/fe.h>
#include <deal.II/base/quadrature_lib.h>

//! Элемент матрици жескости получиный из лапласиана в МКЭ, скалярный случай
/*!
 * Представляет собой интеграл по ячейке:
  \f[
  \int_{cell}
  \sum_{a,b=0}^2 C_{ab}\frac{\partial\phi_i}{\partial x_a}
  \frac{\partial\phi_j}{\partial x_b} 
  \f]
  Интеграл расчитывается в квадратурах.
  \f[
  \sum_{q=0}^{numq}
  \sum_{a,b=0}^2 C_{ab}\frac{\partial\phi_i(q)}{\partial x_a}
  \frac{\partial\phi_j(q)}{\partial x_b} * J(q)*W(q)
  \f]
  ## Выведение
  ### Лапласиан в вариационной записи
  \f[
    \int_{domain}
    \sum_{a,b=0}^2 C_{ab}\frac{\partial U}{\partial x_a}\frac{\partial U}{\partial x_b}
  \f]
  После разбивки расчетной области на конечные элементы заменяем искомую функцию \f$U\f$ на
  сумму функций формы с параметрами \f$\sum_{m=0}^3\alpha_m\phi_m\f$. Далее уравнение 
  записывается для одного контрольного объёма.
  \f[
    \int_{cell}  
    \sum_{a,b=0}^2 C_{ab}\frac{(\sum_{m=0}^3\alpha_m\partial\phi_m)}
    {\partial x_a}\frac{(\sum_{m=0}^3\alpha_m\partial\phi_m)}{\partial x_b}
  \f]
  Проварьируем это уравнение по \f$\alpha_i\f$, где \f$i\f$ есть номер узла в 
  триангуляции, а также номер строки в системе уравнений полученой после вариации. Далее
  одна строка из матрици для ячейки, а не для всей метрици жескости.
  \f[
    \begin{array}{c}
        \int_{cell}\left( \right.
        C_{00}
        (\frac{\partial\phi_{i}}{\partial x_0}\frac{(\sum_{m=0}^3\alpha_m\partial\phi_m)}{\partial x_0}+
        \frac{(\sum_{m=0}^3\alpha_m\partial\phi_m)}{\partial x_0}\frac{\partial\phi_{i}}{\partial x_0})+ \\ \\

        C_{01}
        (\frac{\partial\phi_{i}}{\partial x_0}\frac{(\sum_{m=0}^3\alpha_m\partial\phi_m)}{\partial x_1}+
        \frac{(\sum_{m=0}^3\alpha_m\partial\phi_m)}{\partial x_0}\frac{\partial\phi_{i}}{\partial x_1})+ \\ \\

        C_{10}
        (\frac{\partial\phi_{i}}{\partial x_1}\frac{(\sum_{m=0}^3\alpha_m\partial\phi_m)}{\partial x_0}+
        \frac{(\sum_{m=0}^3\alpha_m\partial\phi_m)}{\partial x_1}\frac{\partial\phi_{i}}{\partial x_0})+ \\ \\

        C_{11}
        (\frac{\partial\phi_{i}}{\partial x_1}\frac{(\sum_{m=0}^3\alpha_m\partial\phi_m)}{\partial x_1}+
       \left.  \frac{(\sum_{m=0}^3\alpha_m\partial\phi_m)}{\partial x_1}\frac{\partial\phi_{i}}{\partial x_1})\right)
\end{array}
  \f]
  Разбиваем интеграл суммы на сумму интегралов по множителям \f$\alpha_j\f$, которые являются искомыми
  величинами и составляют вектор. Далее запишем \f$j\f$-тый интеграл, который является одним элементом матрици,
  помноженный на \f$\alpha_j\f$.
  \f[
    \begin{array}{c}
        \alpha_{j}* \int_{cell}\left[

        C_{00}
        (\frac{\partial\phi_{i}}{\partial x_0}\frac{\partial\phi_{j}}{\partial x_0}+
        \frac{\partial\phi_{j}}{\partial x_0}\frac{\partial\phi_{i}}{\partial x_0})+ \\ \\

        C_{01}
        (\frac{\partial\phi_{i}}{\partial x_0}\frac{\partial\phi_{j}}{\partial x_1}+
        \frac{\partial\phi_{j}}{\partial x_0}\frac{\partial\phi_{i}}{\partial x_1})+ \\ \\

        C_{10}
        (\frac{\partial\phi_{i}}{\partial x_1}\frac{\partial\phi_{j}}{\partial x_0}+
        \frac{\partial\phi_{j}}{\partial x_1}\frac{\partial\phi_{i}}{\partial x_0})+ \\ \\

        C_{11}
        (\frac{\partial\phi_{i}}{\partial x_1}\frac{\partial\phi_{j}}{\partial x_1}+
        \frac{\partial\phi_{j}}{\partial x_1}\frac{\partial\phi_{i}}{\partial x_1}) \right] = \\ \\

        \int_{cell}
        \alpha_{j}*2*\sum_{a,b=0}^2C_{ab}\frac{\partial\phi_{i}}{\partial x_a}\frac{\partial\phi_{j}}{\partial x_b} \\ \\
\end{array}
  \f]
  Двойка сокращается, так как она есть в правой части функционала энергии \f$F(u)=(Lu,u)-2(u,f)\f$. А
  \f$\alpha_j\f$ идёт в вектор. Получаем нужное выражение элемента ячейковой матрици. 
  \f[
  \int_{cell}
  \sum_{a,b=0}^2 C_{ab}\frac{\partial\phi_i}{\partial x_a}
  \frac{\partial\phi_j}{\partial x_b} 
  \f]
*/
template <u8 dim>
class LaplacianScalar : public LaplacianInterface<dim>
{
    public:
    LaplacianScalar (const dealii::FiniteElement<dim> &fe);

    virtual void update_on_cell (
        typename dealii::DoFHandler<dim>::active_cell_iterator &cell) override;

    virtual dbl operator() (cst i, cst j) override;

    virtual u8 get_dofs_per_cell () override;


    vec<ATools::SecondOrderTensor> C; //!< Матрица коэффициентов, например теплопроводности.
    dealii::QGauss<dim>        quadrature_formula; //!< Формула интегрирования в квадратурах.
    dealii::FEValues<dim, dim> fe_values; //!< Тип функций формы.
    cu8                        dofs_per_cell; //!< Количество узлов в ячейке (зависит от типа функций формы).
    cu8                        num_quad_points; //!< Количество точек по которым считается квадратура.
    u8                         material_id = 0; //!< Идентефикатор материала ячейки.
};

template <u8 dim>
LaplacianScalar<dim>::LaplacianScalar (const dealii::FiniteElement<dim> &fe) :
    quadrature_formula (2),
    fe_values (fe, quadrature_formula,
            dealii::update_gradients | dealii::update_quadrature_points | 
            dealii::update_JxW_values),
    dofs_per_cell (fe.dofs_per_cell),
    num_quad_points (quadrature_formula.size())
{

};

template <u8 dim>
void LaplacianScalar<dim>::update_on_cell (
        typename dealii::DoFHandler<dim>::active_cell_iterator &cell)
{
    fe_values .reinit (cell);
    material_id = cell->material_id();
    // printf("%ld %f\n", material_id, C[material_id][0][0]);
};

template <u8 dim>
dbl LaplacianScalar<dim>::operator () (cst i, cst j)
{
    dbl res = 0.0;

    FOR (q_point, 0, num_quad_points)
        FOR (a, 0, dim)
            FOR (b, 0, dim)
            {
                // printf("%ld %ld %d\n", i, j, material_id);
                res += C[material_id][a][b] * 
                    fe_values.shape_grad (i, q_point)[a] *
                    fe_values.shape_grad (j, q_point)[b] *
                    fe_values.JxW(q_point);
            };

    return res;
};

template <u8 dim>
u8 LaplacianScalar<dim>::get_dofs_per_cell ()
{
    return dofs_per_cell;
};

#endif
