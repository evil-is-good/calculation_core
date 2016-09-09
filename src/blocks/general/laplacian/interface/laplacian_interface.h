#ifndef LAPLACIAN_INTERFACE
#define LAPLACIAN_INTERFACE
 
#include "../../../../../../prmt_sintactic_addition/prmt_sintactic_addition.h"
#include <deal.II/dofs/dof_handler.h>

//! Интерфейс к классу являющимся элементом матрици жескости
/*!
 * Класс "элемент матрици жескости" используется в функциях
 * асуществляющих сборку матрици жескости из значений возвращаемых 
 * этим классом.  \n
 * Этот класс принимает на вход через свои методы 
 * ячейку из dealii::DoFHandler, два индекса элемента ячейковой 
 * матрици. \n 
 * Возвращает через оператор () элемент матрици. \n
 * Элементом матрици является интеграл по ячейке от произведения 
 * двух производных функций формы.
  \f[
  \int_{cell}
  \sum_{a,b=0}^2 C_{ab}\frac{\partial\phi_i}{\partial a}
  \frac{\partial\phi_j}{\partial b}
  \f]
*/
template <u8 dim>
class LaplacianInterface
{
    public:

    /*!
     * Обновление данных в классе, единых для одной ячейки.
    */
    virtual void update_on_cell (
        typename dealii::DoFHandler<dim>::active_cell_iterator &cell) = 0;

    /*!
     * Принимает индексы ячейковой матрици, которые являются номерами
     * функций формы, определённых на этой ячейке. \n
     * Возвращает элемент этой матрици.
    */
    virtual dbl operator() (cst i, cst j) = 0;

    virtual u8 get_dofs_per_cell () = 0;
};

#endif
