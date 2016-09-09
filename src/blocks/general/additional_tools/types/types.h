#ifndef types_def
#define types_def 1

namespace ATools
{
    using SecondOrderTensor = arr<arr<dbl, 3>, 3>;

    using FourthOrderTensor = arr<arr<arr<arr<dbl, 3>, 3>, 3>, 3>;

    void copy_tensor(const SecondOrderTensor& from, SecondOrderTensor& to)
    {
        for (st i = 0; i < 3; ++i)
        {
            for (st j = 0; j < 3; ++j)
            {
                to[i][j] = from[i][j];
            };
        };
    };

    void copy_tensor(const FourthOrderTensor& from, FourthOrderTensor& to)
    {
        for (st i = 0; i < 3; ++i)
        {
            for (st j = 0; j < 3; ++j)
            {
                for (st k = 0; k < 3; ++k)
                {
                    for (st l = 0; l < 3; ++l)
                    {
                        to[i][j][k][l] = from[i][j][k][l];
                    };
                };
            };
        };
    };
}

#endif
