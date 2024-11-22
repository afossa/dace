/******************************************************************************
*                                                                             *
* DIFFERENTIAL ALGEBRA CORE ENGINE                                            *
*                                                                             *
*******************************************************************************
*                                                                             *
* Copyright 2016 Politecnico di Milano (2014 Dinamica Srl)                    *
* Licensed under the Apache License, Version 2.0 (the "License");             *
* you may not use this file except in compliance with the License.            *
* You may obtain a copy of the License at                                     *
*                                                                             *
*    http://www.apache.org/licenses/LICENSE-2.0                               *
*                                                                             *
* Unless required by applicable law or agreed to in writing, software         *
* distributed under the License is distributed on an "AS IS" BASIS,           *
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.    *
* See the License for the specific language governing permissions and         *
* limitations under the License.                                              *
*                                                                             *
*******************************************************************************/

/*
 * StatExtension.cpp
 *
 *  Created on: Sep. 12, 2024
 *      Author: Alberto Fossa'
 */

#include <cmath>
#include <algorithm>
#include <numeric>
#include <vector>
#include <utility>

#include "dace/DA.h"
#include "dace/AlgebraicVector.h"
#include "dace/AlgebraicVector_t.h"
#ifdef WITH_ALGEBRAICMATRIX
    #include "dace/AlgebraicMatrix.h"
    #include "dace/AlgebraicMatrix_t.h"
#endif /* WITH_ALGEBRAICMATRIX */
#include "dace/StatExtension.h"

namespace DACE {

namespace {

long getProdFactorials(const vectorui& v) {
/*! Compute the product between the factorials of the elements of the input vector.
    \param[in] v vector of non-negative integers.
    \return The product between factorials.
 */

    long double s = 0.0;
    for (unsigned int i : v) {
        s += std::lgammal(i + 1);
    }

    return std::lround(std::exp(s));
}

void getMultiIndices(const unsigned int no, const unsigned int nv, vectorui& curr, const unsigned int sum, const unsigned int idx, matrixui& res) {
/*! Recursive function to generate all multi-indices of order no in nv variables.
    \param[in] no order of the multi-indices.
    \param[in] nv number of variables.
    \param[in] curr current multi-index.
    \param[in] sum current sum of the multi-index.
    \param[in] idx current index.
    \param[out] res vector of multi-indices.
 */

    if (idx == nv) {
        if (sum <= no) {
            res.push_back(curr);
        }
        return;
    }

    for (unsigned int i = 0; i <= no; ++i) {
        if (sum + i <= no) {
            curr[idx] = i;
            getMultiIndices(no, nv, curr, sum + i, idx + 1, res);
        }
    }
}

bool compare(const vectorui& a, const vectorui& b) {
/*! Compare two multi-indices.
    \param[in] a first multi-index.
    \param[in] b second multi-index.
    \return True if the order of a is smaller than that of b, or the orders
     are the same and a is lexicographically greater than b. False otherwise.
 */

    int sa = std::accumulate(a.begin(), a.end(), 0);
    int sb = std::accumulate(b.begin(), b.end(), 0);
    if (sa != sb) {
        return sa < sb;
    }
    return a > b;
}

}

matrixui getMultiIndices(const unsigned int no, const unsigned int nv) {
/*! Get all multi-indices of order no in nv variables.
    \param[in] no order of the multi-indices.
    \param[in] nv number of variables.
    \return A vector of multi-indices.
 */

    matrixui res;
    vectorui curr(nv, 0);
    getMultiIndices(no, nv, curr, 0, 0, res);
    std::sort(res.begin(), res.end(), compare);

    return res;
}

std::pair<matrixui, vectordb> getRawMoments(const DA& mgf, const unsigned int no) {
/*! Compute all raw moments up to order no.
    \param[in] mgf Taylor expansion of the moment generating function around zero.
    \param[in] no maximum order of the moments.
    \return A pair containing the multi-indices and the raw moments.
 */

    matrixui idx = getMultiIndices(no, DA::getMaxVariables());
    vectordb m0(idx.size());
    for (size_t i = 0; i < idx.size(); i++) {
        m0[i] = getProdFactorials(idx[i])*mgf.getCoefficient(idx[i]);
    }

    return std::make_pair(idx, m0);
}

std::pair<matrixui, vectordb> getCentralMoments(const DA& mgf, const unsigned int no) {
/*! Compute all central moments up to order no.
    \param[in] mgf Taylor expansion of the moment generating function around zero.
    \param[in] no maximum order of the moments.
    \return A pair containing the multi-indices and the central moments.
 */

    vectordb mu = getRawMoments(mgf, 1).second;
    vectorDA id = vectorDA::identity();
    DA ip = std::inner_product(id.begin(), id.end(), mu.begin() + 1, DA(0.0));
    DA mt = (-ip).exp()*mgf;

    return getRawMoments(mt, no);
}

#ifdef WITH_ALGEBRAICMATRIX
DA getMGFGaussian(const vectordb& mu, const matrixdb& cov) {
/*! Return the moment generating function of a Gaussian distribution.
    \param[in] mu mean vector.
    \param[in] cov covariance matrix.
    \return The moment generating function of the Gaussian distribution.
 */

    vectorDA t = vectorDA::identity();
    DA mut = std::inner_product(t.begin(), t.end(), mu.begin(), DA(0.0));
    DA tPt = std::inner_product(t.begin(), t.end(), (cov*t).begin(), DA(0.0));

    return (mut + 0.5*tPt).exp();
}
#endif /* WITH_ALGEBRAICMATRIX */

}
