#ifndef GRID_GENERATOR_e4rt545t4etg
#define GRID_GENERATOR_e4rt545t4etg 1

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_reordering.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_out.h>

namespace GridGenerator
{
    void extrude_triangulation(
            const dealii::Triangulation<2> &tria2d,
            cst n_slices,
            cdbl hight,
            dealii::Triangulation<3> &tria3d
            )
    {
        std::vector<dealii::Point<3> > points(n_slices*tria2d.n_vertices());
        std::vector<dealii::CellData<3> > cells;
        cells.reserve((n_slices-1)*tria2d.n_active_cells());

        // copy the array of points as many times as there will be slices,
        // one slice at a time
        for (st slice = 0; slice < n_slices; ++slice)
        {
            for (st i = 0; i < tria2d.n_vertices(); ++i)
            {
                const dealii::Point<2> &v = tria2d.get_vertices()[i];
                points[slice*tria2d.n_vertices()+i](0) = v(0);
                points[slice*tria2d.n_vertices()+i](1) = v(1);
                points[slice*tria2d.n_vertices()+i](2) = hight * slice / (n_slices-1);
            }
        }

        // then create the cells of each of the slices, one stack at a
        // time
        st count = 0;
        auto cell = tria2d.begin_active();
        auto endc = tria2d.end();
        for (; cell != endc; ++cell)
        {
            ++count;
            for (st slice = 0; slice < n_slices-1; ++slice)
            {
                dealii::CellData<3> this_cell;
                for (st v = 0; v < dealii::GeometryInfo<2>::vertices_per_cell; ++v)
                {
                    this_cell.vertices[v]
                        = cell->vertex_index(v)+slice*tria2d.n_vertices();
                    this_cell.vertices[v+dealii::GeometryInfo<2>::vertices_per_cell]
                        = cell->vertex_index(v)+(slice+1)*tria2d.n_vertices();
                }

                this_cell.material_id = cell->material_id();
                cells.push_back(this_cell);
            }
        }
        tria3d .create_triangulation (points, cells, dealii::SubCellData());
    };

    void gen_cylinder_true_ordered_speciment(
            dealii::Triangulation<3> &triangulation,
            cdbl R,
            cst n_ref_cells,
            cst level_density,
            cst n_slices,
            cdbl H)
    {
        enum {x, y, z};
        cdbl eps = 1.0e-8;

        cst n_cells = 1 << n_ref_cells;
        cst all_n_cells = n_cells * n_cells;
        cdbl start_center = 1.0 / (2 << n_ref_cells);
        cdbl size_cell = 1.0 / n_cells;
        cdbl Rc = R / n_cells;
        vec<dealii::Point<2>> center(all_n_cells);
        {
            st count = 0;
            dealii::Point<2> c(start_center, start_center);
            for (st i = 0; i < n_cells; ++i)
            {
                c(y) = start_center;
                for (st j = 0; j < n_cells; ++j)
                {
                    center[count] = c;
                    c(y) += size_cell;
                    ++count;
                };
                c(x) += size_cell;
            };
        };

        dealii::Triangulation<2> tria2d;
        dealii::GridGenerator ::hyper_cube (tria2d, 0.0, 1.0);
        tria2d .refine_global (level_density);

        cdbl max_dist = std::sqrt(std::pow(1.0 / std::pow(2.0, level_density), 2.0) * 2.0) / 2.0;

        auto face = tria2d.begin_active_face();
        auto endf = tria2d.end_face();
        for (; face != endf; ++face)
        {
            for (st i = 0; i < 2; ++i)
            {
                auto p = face->vertex(i);
                for (st j = 0; j < all_n_cells; ++j)
                {
                    cdbl r = center[j].distance(p); 
                    cdbl dist = std::abs(r - Rc);
                    if (dist < max_dist)
                    {
                        face->vertex(i) -= center[j];
                        face->vertex(i) *= Rc/r;
                        face->vertex(i) += center[j];
                    };
                };
            };
        };

        {
            auto cell = tria2d .begin_active();
            auto end_cell = tria2d .end();
            for (; cell != end_cell; ++cell)
            {
                cell->set_material_id(0);
                for (st i = 0; i < all_n_cells; ++i)
                {
                    if (center[i].distance(cell->center()) < Rc)
                        // if (
                        //         (std::abs(center[i](x) - cell->center()(x)) < Rc) and
                        //         (std::abs(center[i](y) - cell->center()(y)) < Rc)
                        //         )
                    {
                        cell->set_material_id(1);
                    };
                };
            };
        };

        {
            std::ofstream out ("grid-speciment.eps");
            dealii::GridOut grid_out;
            grid_out.write_eps (tria2d, out);
        };

        extrude_triangulation (tria2d, n_slices, H, triangulation);
    };
}

#endif
