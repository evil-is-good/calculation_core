#ifndef problem_on_cell_tools_def
#define problem_on_cell_tools_def 1

#include "../../../../../../prmt_sintactic_addition/prmt_sintactic_addition.h"
#include <deal.II/lac/vector.h>

namespace OnCell
{
//! Массив для хранения метакоэффициентов, ячейковых функций и их производных, с доступом по векторному 
//  индексу в формате записи Горынина.
    template <typename T>
    class ArrayWithAccessToVector
    {
        public:
            ArrayWithAccessToVector () : size(0) {};

            ArrayWithAccessToVector (cst a_size) : size (a_size)//, empty(1)
            {
                content .resize (a_size);
                for (auto &&a : content)
                {
                    a .resize (a_size);
                    for (auto &&b : a)
                    {
                        b .resize (a_size);
                        // for (auto &&c : b)
                        // {
                        //     c = 0.0;
                        // };
                    };
                };
                // content[0][0][0] = 1.0; // Нулевое приближение ячейковой функции равно 1.0
                // for (st i = 0; i < size; ++i)
                // {
                //     for (st j = 0; j < size; ++j)
                //     {
                //         for (st k = 0; k < size; ++k)
                //         {
                //             content[i][j][k] = 0.0;
                //         };
                //     };
                // };
            };
            ArrayWithAccessToVector (ArrayWithAccessToVector &a) : size (a.size)//, empty(1)
            {
                content .resize (a.size);
                for (auto &&a : content)
                {
                    a .resize (a.size);
                    for (auto &&b : a)
                    {
                        b .resize (a.size);
                        // for (auto &&c : b)
                        // {
                        //     c = 0.0;
                        // };
                    };
                };
                for (st i = 0; i < a.size; ++i)
                {
                    for (st j = 0; j < a.size; ++j)
                    {
                        for (st k = 0; k < a.size; ++k)
                        {
                            content[i][j][k] = a.content[i][j][k];
                        };
                    };
                };
            };

            void reinit (cst a_size)
            {
                size = a_size;

                content .resize (a_size);
                for (auto &&a : content)
                {
                    a .resize (a_size);
                    for (auto &&b : a)
                    {
                        b .resize (a_size);
                    };
                };
            };

            ArrayWithAccessToVector<T>& operator= (const ArrayWithAccessToVector<T> &rhs)
            {
                if (size == rhs.size)
                {
                    for (st i = 0; i < size; ++i)
                    {
                        for (st j = 0; j < size; ++j)
                        {
                            for (st k = 0; k < size; ++k)
                            {
                                content[i][j][k] = rhs.content[i][j][k];
                            };
                        };
                    };
                }
                else
                {
                    printf("\x1B[31m ERROR in ArrayWithAccessToVector operator=, size != rhs.size \x1B[0m\n");
                };

                return *this;
            };


            static cst x = 0;
            static cst y = 1;
            static cst z = 2;

            T& operator [] (const arr<i32, 3> k)
            {
                // printf("%d %d %d, %ld\n", k[x], k[y], k[z], size);
                if ((k[x] >=   0) and (k[y] >=   0) and (k[z] >=   0) and
                    (k[x] < size) and (k[y] < size) and (k[z] < size))
                {
                    // puts("da");
                    return content[k[x]][k[y]][k[z]];
                } 
                else
                    return empty;
                // return content.at(k[x]).at(k[y]).at(k[z]);
            };

            const T& operator [] (const arr<i32, 3> k) const
            {
                // printf("%d %d %d, %ld\n", k[x], k[y], k[z], size);
                if ((k[x] >=   0) and (k[y] >=   0) and (k[z] >=   0) and
                    (k[x] < size) and (k[y] < size) and (k[z] < size))
                {
                    // puts("da");
                    return content[k[x]][k[y]][k[z]];
                } 
                else
                    return empty;
                // return content.at(k[x]).at(k[y]).at(k[z]);
            };

            vec<vec<vec<T>>> content;
            st size;
            T empty;
    };

    template <typename T>
    class ArrayWithAccessToVector2
    {
        public:
            ArrayWithAccessToVector2 () : size(0) {};

            ArrayWithAccessToVector2 (cst a_size) : size (a_size)//, empty(1)
            {
                switch (a_size)
                {
                    case 1:
                        content .resize (1);
                        break;
                    case 2:
                        content .resize (4);
                        break;
                    case 3:
                        content .resize (10);
                        break;
                    default:
                    printf("\x1B[31m ERROR in ArrayWithAccessToVector a_size > 3 \x1B[0m\n");
                        break;
                };
            };
            ArrayWithAccessToVector2 (ArrayWithAccessToVector2 &a) : size (a.size)//, empty(1)
            {
                switch (a.size)
                {
                    case 1:
                        content .resize (1);
                        break;
                    case 2:
                        content .resize (4);
                        break;
                    case 3:
                        content .resize (10);
                        break;
                };
                for (st i = 0; i < a.size; ++i)
                {
                    content[i] = a.content[i];
                };
            };

            void reinit (cst a_size)
            {
                switch (a_size)
                {
                    case 1:
                        content .resize (1);
                        break;
                    case 2:
                        content .resize (4);
                        break;
                    case 3:
                        content .resize (10);
                        break;
                };
            };

            ArrayWithAccessToVector2<T>& operator= (const ArrayWithAccessToVector2<T> &rhs)
            {
                if (size == rhs.size)
                {
                    for (st i = 0; i < rhs.size; ++i)
                    {
                        content[i] = rhs.content[i];
                    };
                }
                else
                {
                    printf("\x1B[31m ERROR in ArrayWithAccessToVector operator=, size != rhs.size \x1B[0m\n");
                };

                return *this;
            };


            static cst x = 0;
            static cst y = 1;
            static cst z = 2;

            T& operator [] (const arr<i32, 3> k)
            {
                if ((k[x] >=   0) and (k[y] >=   0) and (k[z] >=   0) and
                    (k[x] < size) and (k[y] < size) and (k[z] < size))
                {
                    for (st i = 0; i < size; ++i)
                    {
                        if ((k[x] == index[i][x]) and (k[y] == index[i][y]) and (k[z] == index[i][z]))
                            return content[i];
                    };
                    return empty;
                } 
                else
                    return empty;
            };

            const T& operator [] (const arr<i32, 3> k) const
            {
                if ((k[x] >=   0) and (k[y] >=   0) and (k[z] >=   0) and
                    (k[x] < size) and (k[y] < size) and (k[z] < size))
                {
                    for (st i = 0; i < size; ++i)
                    {
                        if ((k[x] == index[i][x]) and (k[y] == index[i][y]) and (k[z] == index[i][z]))
                            return content[i];
                    };
                    return empty;
                } 
                else
                    return empty;
            };

            vec<T> content;
            // vec<vec<vec<T>>> content;
            st size;
            T empty;
            ci32 index[10][3] = 
            {
                {0, 0, 0},
                {1, 0, 0},
                {0, 1, 0},
                {0, 0, 1},
                {2, 0, 0},
                {1, 1, 0},
                {1, 0, 1},
                {0, 2, 0},
                {0, 1, 1},
                {0, 0, 2}
            };
    };

};

#endif
