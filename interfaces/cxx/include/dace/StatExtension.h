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
 * StatExtension.h
 *
 *  Created on: Sep. 12, 2024
 *      Author: Alberto Fossa'
 */

#ifndef DINAMICA_STATEXTENSION_H_
#define DINAMICA_STATEXTENSION_H_

#include <vector>
#include <utility>

#include "dace/DA.h"
#include "dace/AlgebraicVector.h"
#ifdef WITH_ALGEBRAICMATRIX
    #include "dace/AlgebraicMatrix.h"
#endif /* WITH_ALGEBRAICMATRIX */

namespace DACE {

typedef std::vector<unsigned int> vectorui; //!< Shorthand notation for std::vector<unsigned int>.
typedef std::vector<std::vector<unsigned int>> matrixui; //!< Shorthand notation for std::vector<std::vector<unsigned int>>.

DACE_API matrixui getMultiIndices(const unsigned int no, const unsigned int nv); //!< Get all multi-indices of order no in nv variables
DACE_API std::pair<matrixui, vectordb> getRawMoments(const DA& mgf, const unsigned int no = DA::getMaxOrder()); //!< Get raw moments up to order no
DACE_API std::pair<matrixui, vectordb> getCentralMoments(const DA& mgf, const unsigned int no = DA::getMaxOrder()); //!< Get central moments up to order no
#ifdef WITH_ALGEBRAICMATRIX
    DACE_API DA getMGFGaussian(const vectordb& mu, const matrixdb& cov); //!< Get the Taylor expansion of the moment generating function of a Gaussian distribution
#endif /* WITH_ALGEBRAICMATRIX */

}

#endif /* DINAMICA_STATEXTENSION_H_ */
