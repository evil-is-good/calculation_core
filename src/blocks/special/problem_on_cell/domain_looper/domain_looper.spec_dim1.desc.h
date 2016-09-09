/*
 * =====================================================================================
 *
 *       Filename:  domain_looper.spec_dim1.desc.h
 *
 *    Description:  for dim = 1
 *
 *        Version:  1.0
 *        Created:  31.08.2012 11:53:41
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */
#ifndef DOMAIN_LOOPER_SPEC_DIM1_DESC

#define DOMAIN_LOOPER_SPEC_DIM1_DESC

//#include "./domain_looper.desc.h"

template<bool type_space>
class DomainLooper<1, type_space>
{
    public:
        DomainLooper ();
        ~DomainLooper ();

        const static uint8_t dim = 1;

    //METHODS
    public:
        void loop_domain (const dealii::DoFHandler<dim> &dof_h,
                            OnCell::BlackOnWhiteSubstituter &bows,
                            dealii::CompressedSparsityPattern &csp);

    OPERATOR_REPORT;
};

#endif
