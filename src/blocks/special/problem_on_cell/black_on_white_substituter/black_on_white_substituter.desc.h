#ifndef BLACK_ON_WHITE_SUBSTITUTER_DESC_2

#define BLACK_ON_WHITE_SUBSTITUTER_DESC_2

#include "../../../../../../prmt_sintactic_addition/prmt_sintactic_addition.h"

namespace OnCell
{
    //! Замена черного индекса точки на соответствующий белый
    /*!
     * Черный индекс - это индекс "виртуальной" точки, эта точка лежит на
     * границе на которую нахлёстывается другая граница, поэтому она фактически
     * заменяется другой точкой. \n
     * Белый индекс - индекс точки которая заменяет "виртуальную" точку.
    */
    class BlackOnWhiteSubstituter
    {
        public:
            BlackOnWhiteSubstituter  ();
            BlackOnWhiteSubstituter  (cst n);
            ~BlackOnWhiteSubstituter ();

            //METHODS
        public:
            st   subst    (cst index) const; //!< Если index черный, то заменяет его на белый, иначе оставляет тот же
            bool is_black (cst index) const; //!< Проверить черный index или нет

        // private:
            void add_white_and_black (cst w, cst b); //!< Добавить пару связанных точек

            //FIELDS
        public:
            vec<st> black;
            vec<st> white;
            st  size;

        public:
            template<u8 dim, bool type_space>
                friend class DomainLooper;
            template<u8 dim, u8 type_space>
                friend class DomainLooperTrivial;
    };
};

#endif
