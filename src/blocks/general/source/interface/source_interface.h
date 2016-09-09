#ifndef SOURCE_INTERFACE
#define SOURCE_INTERFACE
 
#include "../../../../../../prmt_sintactic_addition/prmt_sintactic_addition.h"
#include <deal.II/dofs/dof_handler.h>

//! Интерфейс к классу являющимся элементом вектора правой части
/*!
 * Класс "элементом вектора правой части" используется в функциях
 * асуществляющих сборку вектора правой части из значений возвращаемых 
 * этим классом.  \n
 * Этот класс принимает на вход через свои методы 
 * ячейку из dealii::DoFHandler, индекс элемента ячейкового
 * вектора. \n 
 * Возвращает через оператор () элемент вектора. \n
 * Элементом вектора является интеграл по ячейке от произведения 
 * функции источника и функции формы.
  \f[
  \int_{cell} f \phi
  \f]
*/
template <u8 dim>
class SourceInterface
{
    public:

    /*!
     * Обновление данных в классе, единых для одной ячейки.
    */
    virtual void update_on_cell (
        typename dealii::DoFHandler<dim>::active_cell_iterator &cell) = 0;

    /*!
     * Принимает индекс ячейкового  вектора, который являются номером
     * функций формы, определённой на этой ячейке. \n
     * Возвращает элемент этого мектора.
    */
    virtual dbl operator() (cst i) = 0;

    virtual u8 get_dofs_per_cell () = 0;
};

#endif
