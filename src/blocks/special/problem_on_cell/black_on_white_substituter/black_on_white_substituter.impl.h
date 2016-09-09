#ifndef BLACK_ON_WHITE_SUBSTITUTER_IMPL_2

#define BLACK_ON_WHITE_SUBSTITUTER_IMPL_2

namespace OnCell
{
    BlackOnWhiteSubstituter::BlackOnWhiteSubstituter()
        :
            size (0)
    {
    };

    BlackOnWhiteSubstituter::BlackOnWhiteSubstituter(cst n)
        :
            white (n), black (n), size (n)
    {
    };

    size_t BlackOnWhiteSubstituter::subst (cst index) const
    {
        if (size > 0)
        {
            for (size_t i = 0; i < size; ++i)
            {
                if (black[i] == index)
                    return white[i];
            };
        };

        return index;
    };

    bool BlackOnWhiteSubstituter::is_black (cst index) const
    {
        if (size > 0)
        {
            for (size_t i = 0; i < size; ++i)
            {
                if (black[i] == index)
                    return true;
            };
        };

        return false;
    };

    void BlackOnWhiteSubstituter::add_white_and_black (cst w, cst b)
    {
        white .push_back (w);
        black .push_back (b);
        ++size;
    };

    BlackOnWhiteSubstituter::~BlackOnWhiteSubstituter()
    {
        white .clear ();
        black .clear ();
    };
};

#endif

