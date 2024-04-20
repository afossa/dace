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
 * Interval.cpp
 *
 *  Created on: Apr 17, 2024
 *      Author: Alberto Fossa'
 */

// C++ stdlib classes used only internally in the implementation
#include <sstream>
#include <iomanip>

// DACE classes
#include "dace/config.h"
#include "dace/Interval.h"
#include "dace/DA.h"

namespace DACE{

std::string Interval::toString() const {
/*! Convert interval to string.
   \return A string representing the interval in human readable form.
 */
    std::ostringstream oss;
    oss << std::setprecision(16) << std::scientific;
    oss << "[" << m_lb << ", " << m_ub << "]" << std::endl;
    return oss.str();
}

std::ostream& operator<< (std::ostream &out, const Interval &m){
/*! Overload of std::operator<< in iostream.
   \param[in] out standard output stream.
   \param[in] m Interval to be printed in the stream
   \return The output stream out.
 */
    out << m.toString();
    return out;
}

}