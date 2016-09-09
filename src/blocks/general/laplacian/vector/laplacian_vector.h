#ifndef LAPLACIAN_VECTOR
#define LAPLACIAN_VECTOR
 
#include "../scalar/laplacian_scalar.h"

//! Элемент матрици жескости получиный из лапласиана в МКЭ, векторный случай
/*!
 * Представляет собой интеграл по ячейке:
  \f[
  \int_{cell}
  \sum_{a,b=0}^2 C_{labr}\frac{\partial\phi_{vl}}{\partial x_a}
  \frac{\partial \phi_{wr}}{\partial x_b}
  \f]
  Интеграл расчитывается в квадратурах.
  \f[
  \sum_{q=0}^{numq}
  \sum_{a,b=0}^2 C_{labr}\frac{\partial\phi_{vl}(q)}{\partial x_a}
  \frac{\partial\phi_{wr}(q)}{\partial x_b} J(q)W(q)
  \f]
  Этот случай практически аналогичен LaplacianScalar, только функций
  формы \f$\phi_{n}\f$ больше в размерность раз, тоесть для двумерного
  случая их в два раза больше. Эти функции относятся к разным компанентам 
  вектора решения, то есть для разных направлений. Они чередуются по
  очереди, то есть, например в трёхмерном случае, \f$\phi_0\f$ для
  компоненты \f$x\f$, \f$\phi_1\f$ для \f$y\f$, \f$\phi_2\f$ для \f$z\f$,
  \f$\phi_3\f$ для \f$x\f$ и так далее. Поэтому индекс \f$n\f$ можно
  пердставить в виде:
  \f[
  n = \xi \zeta, \\ 0 \le \xi \le 3, \\ 0 \le \zeta \le dim, \\ \zeta = n \bmod{dim} 
  \f]
  То есть \f$\xi\f$ - это номер вершины, а \f$\zeta\f$ - номер компоненты.
  Таким образом отличие векторного случая от скалярного только в индексах
  \f$l\f$ и \f$r\f$ и в том, что тензор \f$C\f$ четвёртого порядка. Однако
  для конкретных двух функций \f$\phi_{n}\f$, для конкретной точки, этот
  случай совпадает со скалярным, то есть тензор становится второго порядка.
  Поэтому класс LaplacianVector содержит в себе экземпляр класса LaplacianScalar
  и выполняет вспомогательную роль подготавливая тензор.

  ## Выведение 
  Оно, в основном, совпадает с выведением для LaplacianScalar и потому 
  здесь подробно не описывается.
  
  Уравнение упругости.
  \f[
    \sum_{a,b,d=0}^2 \frac{\partial}{\partial x_a}C_{cabd}\frac{\partial U_d}{\partial x_b}
    = f_c
  \f]
  Функционал упругой энергии.
  \f[
    \int_{domain}
    \sum_{a,b,c,d=0}^2 \frac{\partial}{\partial x_a}C_{cabd}\frac{\partial U_d}{\partial x_b}U_c
    = \int_{domain}\sum_c2f_cU_c
  \f]
  \f[
    \int_{domain}
    \sum_{a,b,c,d=0}^2 C_{cabd}\frac{\partial U_c}{\partial x_a}\frac{\partial U_d}{\partial x_b}
  \f]
  \f[
    \begin{array}{c}
    \int_{cell}
    \sum_{c,a,b,d=0}^2 C_{iabl}\frac{(\sum_{m=0}^3\alpha_m\partial\phi_m)_c}
    {\partial x_a}\frac{(\sum_{m=0}^3\alpha_m\partial\phi_m)_d}{\partial x_b}= \\ \\

    \int_{cell}
    \sum_{a,b=0}^2 C_{0ab0}\frac{(\sum_{m=0}^3\alpha_m\partial\phi_m)_0}
    {\partial x_a}\frac{(\sum_{m=0}^3\alpha_m\partial\phi_m)_0}{\partial x_b}+ \\ \\

    \int_{cell}
    \sum_{a,b=0}^2 C_{0ab1}\frac{(\sum_{m=0}^3\alpha_m\partial\phi_m)_0}
    {\partial x_a}\frac{(\sum_{m=0}^3\alpha_m\partial\phi_m)_1}{\partial x_b}+ \\ \\

    \int_{cell}
    \sum_{a,b=0}^2 C_{0ab2}\frac{(\sum_{m=0}^3\alpha_m\partial\phi_m)_0}
    {\partial x_a}\frac{(\sum_{m=0}^3\alpha_m\partial\phi_m)_2}{\partial x_b}+ \\ \\


    \int_{cell}
    \sum_{a,b=0}^2 C_{1ab0}\frac{(\sum_{m=0}^3\alpha_m\partial\phi_m)_1}
    {\partial x_a}\frac{(\sum_{m=0}^3\alpha_m\partial\phi_m)_0}{\partial x_b}+ \\ \\

    \int_{cell}
    \sum_{a,b=0}^2 C_{1ab1}\frac{(\sum_{m=0}^3\alpha_m\partial\phi_m)_1}
    {\partial x_a}\frac{(\sum_{m=0}^3\alpha_m\partial\phi_m)_1}{\partial x_b}+ \\ \\

    \int_{cell}
    \sum_{a,b=0}^2 C_{1ab2}\frac{(\sum_{m=0}^3\alpha_m\partial\phi_m)_1}
    {\partial x_a}\frac{(\sum_{m=0}^3\alpha_m\partial\phi_m)_2}{\partial x_b}+ \\ \\


    \int_{cell}
    \sum_{a,b=0}^2 C_{2ab0}\frac{(\sum_{m=0}^3\alpha_m\partial\phi_m)_2}
    {\partial x_a}\frac{(\sum_{m=0}^3\alpha_m\partial\phi_m)_0}{\partial x_b}+ \\ \\

    \int_{cell}
    \sum_{a,b=0}^2 C_{2ab1}\frac{(\sum_{m=0}^3\alpha_m\partial\phi_m)_2}
    {\partial x_a}\frac{(\sum_{m=0}^3\alpha_m\partial\phi_m)_1}{\partial x_b}+ \\ \\

    \int_{cell}
    \sum_{a,b=0}^2 C_{2ab2}\frac{(\sum_{m=0}^3\alpha_m\partial\phi_m)_2}
    {\partial x_a}\frac{(\sum_{m=0}^3\alpha_m\partial\phi_m)_2}{\partial x_b} \\ \\
\end{array}
  \f]
  Проварьируем по \f$\alpha_{vl}\f$.
  \f[
    \begin{array}{c}
    \int_{cell}
        \sum_{a,b=0}^2 C_{labl}\frac{\partial\phi_{vl}}{\partial x_a}
    \frac{(\sum_{m=0}^3\alpha_m\partial\phi_m)_l}{\partial x_b}+ \\ \\

    \int_{cell}
    \sum_{a,b=0}^2 C_{labl}\frac{(\sum_{m=0}^3\alpha_m\partial\phi_m)_l}{\partial x_a}
    \frac{\partial\phi_{vl}}{\partial x_b}+ \\ \\

    \int_{cell}
    \sum_{a,b=0}^2 C_{lab\overline{l}}\frac{\partial\phi_{vl}}{\partial x_a}
    \frac{(\sum_{m=0}^3\alpha_m\partial\phi_m)_{\overline{l}}}{\partial x_b}+ \\ \\

    \int_{cell}
    \sum_{a,b=0}^2 C_{lab\overline{\overline{l}}}\frac{\partial\phi_{vl}}{\partial x_a}
    \frac{(\sum_{m=0}^3\alpha_m\partial\phi_m)_{\overline{\overline{l}}}}{\partial x_b}+ \\ \\

    \int_{cell}
    \sum_{a,b=0}^2 C_{\overline{l}abl}\frac{(\sum_{m=0}^3\alpha_m\partial\phi_m)_{\overline{l}}}{\partial j}
    \frac{\partial\phi_{vl}}{\partial x_b}+ \\ \\

    \int_{cell}
    \sum_{a,b=0}^2 C_{\overline{\overline{l}}abl}\frac{(\sum_{m=0}^3\alpha_m\partial\phi_m)_{\overline{\overline{l}}}}{\partial x_a}
    \frac{\partial\phi_{vl}}{\partial x_b}= \\ \\
\end{array}
    \f]
    Здесь \f$ l \in \{0, dim\}; \overline{l} \in \{0, dim\}, \overline{l} \neq l;
    \overline{\overline{l}} \in \{0, dim\}, \overline{\overline{l}} \neq l, 
    \overline{\overline{l}} \neq \overline{l} \f$.

    \f[
    \begin{array}{c}
    \int_{cell}
    \sum_{a,b,d=0}^2 C_{dabl}\frac{(\sum_{m=0}^3\alpha_m\partial\phi_m)_d}{\partial x_a}
    \frac{\partial\phi_{vl}}{\partial x_b}+ \\ \\

    \int_{cell}
    \sum_{a,b,d=0}^2 C_{labd}\frac{\partial\phi_{vl}}{\partial x_a}
    \frac{(\sum_{m=0}^3\alpha_m\partial\phi_m)_d}{\partial x_b}= \\ \\

    \int_{cell}
    2*\sum_{a,b,d=0}^2 C_{labd}\frac{\partial\phi_{vl}}{\partial x_a}
    \frac{(\sum_{m=0}^3\alpha_m\partial\phi_m)_d}{\partial x_b} \\ \\
\end{array}
  \f]
  \f[
    \int_{cell}
    2*\alpha_{wr}*\sum_{a,b=0}^2 C_{labr}\frac{\partial\phi_{vl}}{\partial x_a}
    \frac{\partial\phi_{wr}}{\partial x_b} \\ \\
  \f]
  \f[
    \int_{cell}
    \sum_{a,b=0}^2 C_{labr}\frac{\partial\phi_{vl}}{\partial x_a}
    \frac{\partial\phi_{wr}}{\partial x_b} \\ \\
  \f]
  */
template <u8 dim>
class LaplacianVector : public LaplacianInterface<dim>
{
    public:

    LaplacianVector (const dealii::FiniteElement<dim> &fe);

    virtual dbl operator() (cst i, cst j) override;

    virtual void update_on_cell (
        typename dealii::DoFHandler<dim>::active_cell_iterator &cell) override;

    virtual u8 get_dofs_per_cell () override;


    vec<ATools::FourthOrderTensor> C; //!< Тензор коэффициентов, четвёртого порядка, например упругости.
    LaplacianScalar<dim>                         laplacian; 
};

template <u8 dim>
LaplacianVector<dim>::LaplacianVector (const dealii::FiniteElement<dim> &fe) :
    laplacian (fe)
{

};

template <u8 dim>
void LaplacianVector<dim>::update_on_cell (
        typename dealii::DoFHandler<dim>::active_cell_iterator &cell)
{
    laplacian .update_on_cell (cell);
};

template <u8 dim>
dbl LaplacianVector<dim>::operator () (cst i, cst j)
{
    // phi[indx_n] = 
    // A
    // (phi_x[indx_n % dim], phi_y[indx_n % dim], phi_x[indx_n % dim]);
    // indx_n in [dim*dim]

    if (C.size() != laplacian.C.size())
    {
        laplacian.C.resize(C.size());
    };

    cst l = i % dim; 
    cst r = j % dim;

    FOR (n, 0, laplacian.C.size())
        FOR (a, 0, dim)
            FOR (b, 0, dim)
                laplacian.C[n][a][b] = C[n][l][a][b][r];

    return laplacian(i, j);
};

template <u8 dim>
u8 LaplacianVector<dim>::get_dofs_per_cell ()
{
    return laplacian.get_dofs_per_cell();
};

#endif
