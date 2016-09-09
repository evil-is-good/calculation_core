#ifndef SOURCE_VECTOR_ON_CELL
#define SOURCE_VECTOR_ON_CELL
 
#include "../scalar/source_scalar.h"

namespace OnCell
{
    //! Элемент вектора правой части уравнения в МКЭ для задачи на ячейке, векторный случай
    template <u8 dim>
        class SourceVector : public SourceInterface<dim>
    {
        public:
            typedef std::function<arr<dbl, dim> (const dealii::Point<dim>&)> Func;

            SourceVector (
                    const vec<arr<arr<dbl, dim>, dim>> &coefficient, const dealii::FiniteElement<dim> &fe);

            virtual dbl operator() (cst i) override;

            virtual void update_on_cell (
                    typename dealii::DoFHandler<dim>::active_cell_iterator &cell) override;

            virtual u8 get_dofs_per_cell () override;


            vec<arr<arr<dbl, dim>, dim>> coef;
            ::OnCell::SourceScalar<dim> source; 
    };

    template <u8 dim>
        SourceVector<dim>::SourceVector (
                const vec<arr<arr<dbl, dim>, dim>> &coefficient, const dealii::FiniteElement<dim> &fe) :
            source (fe)
    {
        coef .resize (coefficient.size());
        for (st i = 0; i < coefficient.size(); ++i)
        {
            for (st j = 0; j < dim; ++j)
            {
                for (st k = 0; k < dim; ++k)
                {
                    coef[i][j][k] = coefficient[i][j][k];
                };
            };
        };
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
            source.coef .resize (coef.size());
            for (st m = 0; m < coef.size(); ++m)
            {
                for (st n = 0; n < dim; ++n)
                {
                    source.coef[m][n] = coef[m][i % dim][n];
                };
            };

            return source (i);
        };

    template <u8 dim>
        u8 SourceVector<dim>::get_dofs_per_cell ()
        {
            return source.get_dofs_per_cell();
        };

    //! Элемент вектора правой части уравнения в МКЭ для задачи на ячейке для произвольного
    //! номера приближения, векторный случай
    template <u8 dim>
        class SourceVectorApprox : public SourceInterface<dim>
    {
        public:
            SourceVectorApprox (const dealii::FiniteElement<dim> &fe);
            SourceVectorApprox (const arr<i32, 3> approximation,
                    cst nu_indx,
                    const vec<ATools::FourthOrderTensor> &coefficient,
                    const ArrayWithAccessToVector<arr<arr<dbl, 3>, 3>> &meta_coefficient,
                    // ArrayWithAccessToVector<dealii::Vector<dbl>> *psi_func,
                    ArrayWithAccessToVector<arr<dealii::Vector<dbl>, 3>> &cell_function,
                    const dealii::FiniteElement<dim> &fe);

            virtual dbl operator() (cst i) override;

            virtual void update_on_cell (
                    typename dealii::DoFHandler<dim>::active_cell_iterator &cell_p) override;

            virtual u8 get_dofs_per_cell () override;

            static cst x = 0;
            static cst y = 1;
            static cst z = 2;

            vec<ATools::FourthOrderTensor> coef; //!< Коэффициенты материла, например теплопроводность.
            ArrayWithAccessToVector<arr<arr<dbl, 3>, 3>> meta_coef; //!< Массив коэффициентов метаматериала.
            // ArrayWithAccessToVector<dealii::Vector<dbl>> *psi; //!< Массив ячейкивых функций.
            ArrayWithAccessToVector<arr<dealii::Vector<dbl>, 3>> cell_func; //!< Массив ячейкивых функций.
            arr<i32, 3> k; //!< Номер приближения в векторной форме записи.
            dealii::QGauss<dim>        quadrature_formula; //!< Формула интегрирования в квадратурах.
            dealii::FEValues<dim, dim> fe_values; //!< Тип функций формы.
            vec<u32>                   global_dof_indices; //!< Список глобальных значений индексов с ячейки.
            cu8                        dofs_per_cell; //!< Количество узлов в ячейке (зависит от типа функций формы).
            cu8                        num_quad_points; //!< Количество точек по которым считается квадратура.
            u8                         material_id = 0; //!< Идентефикатор материала ячейки.
            st                         nu;
            typename dealii::DoFHandler<dim>::active_cell_iterator cell;
            dbl tmp = 0.0;
    };

    template <u8 dim>
        SourceVectorApprox<dim>::SourceVectorApprox (const dealii::FiniteElement<dim> &fe) :
            quadrature_formula (2),
            fe_values (fe, quadrature_formula,
                    dealii::update_values | dealii::update_gradients | dealii::update_quadrature_points | 
                    dealii::update_JxW_values),
            dofs_per_cell (fe.dofs_per_cell),
            global_dof_indices (fe.dofs_per_cell),
            num_quad_points (quadrature_formula.size())
    {};

    template <u8 dim>
        SourceVectorApprox<dim>::SourceVectorApprox (const arr<i32, 3> approximation,
                    cst nu_indx,
                    const vec<ATools::FourthOrderTensor> &coefficient,
                    const ArrayWithAccessToVector<arr<arr<dbl, 3>, 3>> &meta_coefficient,
                    // ArrayWithAccessToVector<dealii::Vector<dbl>> *psi_func,
                    ArrayWithAccessToVector<arr<dealii::Vector<dbl>, 3>> &cell_function,
                    const dealii::FiniteElement<dim> &fe) :
            SourceVectorApprox<dim>(fe)
    {
        k[x] = approximation[x];
        k[y] = approximation[y];
        k[z] = approximation[z];

        nu = nu_indx;

        coef .resize (coefficient.size());
        for (st i = 0; i < coefficient.size(); ++i)
        {
            for (st j = 0; j < 3; ++j)
            {
                for (st k = 0; k < 3; ++k)
                {
                    for (st l = 0; l < 3; ++l)
                    {
                        for (st m = 0; m < 3; ++m)
                        {
                            coef[i][j][k][l][m] = coefficient[i][j][k][l][m];
                        };
                    };
                };
            };
        };

        meta_coef .reinit (meta_coefficient.size);
        for (st i = 0; i < cell_function.size; ++i)
        {
            for (st j = 0; j < cell_function.size; ++j)
            {
                for (st k = 0; k < cell_function.size; ++k)
                {
                    for (st l = 0; l < 3; ++l)
                    {
                        for (st m = 0; m < 3; ++m)
                        {
                            meta_coef.content[i][j][k][l][m] = meta_coefficient.content[i][j][k][l][m];
                        };
                    };
                };
            };
        };


        cell_func.reinit(cell_function.size);
        for (st i = 0; i < cell_function.size; ++i)
        {
            for (st j = 0; j <cell_function.size; ++j)
            {
                for (st k = 0; k < cell_function.size; ++k)
                {
                    for (st l = 0; l < 3; ++l)
                    {
                        cell_func.content[i][j][k][l].reinit(cell_function.content[i][j][k][l].size());
                        // printf("%dgcc\n",cell_func.content[i][j][k][l].size() );
                        for (st m = 0; m < cell_function.content[i][j][k][l].size(); ++m)
                        {
                            cell_func.content[i][j][k][l][m] = cell_function.content[i][j][k][l][m];
                        };
                    };
                };
            };
        };
        // printf("dfdfvdfv %d\n", cell_func.content[0][0][0][0].size());
        // printf("1111111111111111111\n");
    };

    template <u8 dim>
        void SourceVectorApprox<dim>::update_on_cell (
                typename dealii::DoFHandler<dim>::active_cell_iterator &cell_p)
        {
            fe_values .reinit (cell_p);
            material_id = cell_p->material_id();
            cell_p->get_dof_indices (global_dof_indices);
            cell = cell_p;
        };

    template <u8 dim>
        dbl SourceVectorApprox<dim>::operator () (cst i)
        {
            // const uint8_t num_quad_points = this->quadrature_formula.size();

            dbl res = 0.0;
            dbl summ = 0.0;

            dbl f1 = 0.0;
            dbl f2 = 0.0;
            dbl f3 = 0.0;

            cst alpha = i % dim;

            for (st q_point = 0; q_point < num_quad_points; ++q_point)
            {
                for (st beta = 0; beta < dim; ++beta)
                    // for (st phi = 0; phi < 3; ++phi)
                {
                    for (st m = 0; m < dofs_per_cell; ++m)
                    {
                        cst phi = m % dim;
                            // cst beta = m % 3;
                        for (st psi = 0; psi < dim; ++psi)
                        {
                            arr<i32, 3> k_psi = {k[x], k[y], k[z]}; // k - э_psi
                            bool cell_k_psi_exist = true;
                            if (k[psi] == 0)
                                cell_k_psi_exist = false;
                            else
                                k_psi[psi]--;

                            arr<i32, 3> k_beta = {k[x], k[y], k[z]}; // k - э_beta
                            bool cell_k_beta_exist = true;
                            if (k[beta] == 0)
                                cell_k_beta_exist = false;
                            else
                                k_beta[beta]--;

                            arr<i32, 3> k_psi_beta = {k_psi[x], k_psi[y], k_psi[z]}; // k - э_psi - э_beta
                            bool cell_k_psi_beta_exist = true;
                            if (cell_k_psi_exist)
                            {
                                if (k_psi[beta] == 0)
                                    cell_k_psi_beta_exist = false;
                                else
                                    k_psi_beta[beta]--;
                            };

                            dbl cell_k_psi      = 0.0;
                            dbl cell_k_beta     = 0.0;
                            dbl cell_k_psi_beta = 0.0;

                            if (cell_k_psi_exist)
                                cell_k_psi = cell_func[k_psi][nu][global_dof_indices[m]];
                            else
                                cell_k_psi = 0.0;

                            if (cell_k_beta_exist)
                                cell_k_beta = cell_func[k_beta][nu][global_dof_indices[m]];
                            else
                                cell_k_beta = 0.0;

                            if (cell_k_psi_beta_exist)
                                cell_k_psi_beta = cell_func[k_psi_beta][nu][global_dof_indices[m]];
                            else
                                cell_k_psi_beta = 0.0;
                            // printf("%f\n", cell_k_psi_beta);

                            f1 += 
                                this->coef[this->material_id][alpha][beta][phi][psi] * 
                                cell_k_psi *
                                this->fe_values.shape_value (m, q_point) *
                                this->fe_values.shape_grad (i, q_point)[beta] *
                                this->fe_values.JxW(q_point);

                            f3 += 
                                 (this->coef[this->material_id][alpha][beta][phi][psi] *
                                 (
                                  // 0.0 *
                                  cell_k_beta *
                                  this->fe_values.shape_grad (m, q_point)[psi] +
                                  // 0.0 *
                                  cell_k_psi_beta *
                                  this->fe_values.shape_value (m, q_point)
                                 )
                                ) *
                                this->fe_values.shape_value (i, q_point) *
                                this->fe_values.JxW(q_point);
                        };
                    };
                };
                f2 += 1.0 *(
                        // (std::sin(2*M_PI*this->fe_values.quadrature_point(q_point)(0)) + 1.0) *
                        // (2.0*(cell->vertex(i)(0))) *
                        // 2.0*this->fe_values.quadrature_point(q_point)(0)*
                        // 3.0*this->fe_values.quadrature_point(q_point)(0)*this->fe_values.quadrature_point(q_point)(0)*
                        // 00.0*
                    this->meta_coef[k][nu][alpha] *
                    this->fe_values.shape_value (i, q_point) *
                    this->fe_values.JxW(q_point));
            };
            // tmp += f3-f2;
            // tmp += f2-f3;
            // printf("%f %f %f %f %ld\n", f2, f3, tmp, this->meta_coef[k], k[0]);
            tmp += f3-f2-f1;

            return f3-f2-f1;
            // return -f3;
        };

    template <u8 dim>
        u8 SourceVectorApprox<dim>::get_dofs_per_cell ()
        {
            return dofs_per_cell;
        };
};
#endif
