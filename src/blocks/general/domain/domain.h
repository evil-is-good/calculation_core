/*
 * =====================================================================================
 *
 *       Filename:  domain.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  13.09.2012 11:44:17
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */

#ifndef CALCULATION_CORE_DOMAIN
#define CALCULATION_CORE_DOMAIN

#include "../../../../../prmt_sintactic_addition/prmt_sintactic_addition.h"
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe.h>

//! Область решения
template <u8 dim>
struct Domain
{
    //! Инициализация степеней свободы
    void dof_init (const dealii::FiniteElement<dim> &fe) 
    {
        dof_handler .initialize (grid, fe);
    };

    //! Сетка
    dealii::Triangulation <dim> grid;
    //! Степени свободы
    dealii::DoFHandler    <dim> dof_handler;
};

#endif
