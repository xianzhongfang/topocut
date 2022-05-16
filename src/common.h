/*
  This file is part of TopoCut (TopoCut: Fast and Robust Planar Cutting of Arbitrary Domains)
  
  Copyright (C) 2022,  Xianzhong Fang
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:

  * Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

  * Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

  * Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
  
  Author: Xianzhong Fang, Email: xzfangcs@163.com
*/


#ifndef COMMON_FXZ_H
#define COMMON_FXZ_H
#include <cstddef>
#include <chrono>
#include <bits/stdc++.h>

namespace fxz {
  
  const size_t INVALID_NUM = (size_t)(-1);
  const size_t UNSIGN_NEG_ONE = (size_t)(-1);
const double kPI = 3.14159265;

#define CALL_FUNC(pn) {                         \
    if(pn){ return __LINE__; }}


#define CALL_FUNC_WITH_CLOCK(pn, func_name) {                           \
    auto start = std::chrono::high_resolution_clock::now();             \
    if (pn) { std::cerr << "# [ ERROR ] " << #pn << std::endl; return __LINE__; } \
    auto end = std::chrono::high_resolution_clock::now();               \
    double time_taken=std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count(); \
    time_taken *= 1e-9;                                                 \
    std::cout << "--- Time taken by program " << func_name << " : " << std::fixed \
              << time_taken << std::setprecision(9) << "  sec" << std::endl; \
  }


#define CALL_SUBPROG(sp, func) {if (prog==sp) {CALL_FUNC(func(argc,argv));}}

#define FASSERT(func) {if(!(func)){                             \
      std::cerr<<"# [ ERROR ] " << " L"                         \
               << __LINE__ << ": assert "<<#func<<std::endl;    \
      exit(1);}}


#define FASSERT2(func, str) {if(!(func)){std::cerr<<"# [ ERROR ] assert "<<#func << " | " << str <<std::endl; return 1;}}

#define ERROR_RETURN { exit(1);return 5; }

} // namespace fxz

#include <Eigen/Dense>

namespace fxz {
  namespace eigen {
    typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> ematrixd;
    typedef Eigen::Matrix<double,Eigen::Dynamic,1> evecd;
    typedef Eigen::Matrix<double,3,1> evec3d;
    typedef Eigen::Map<const ematrixd> cmap_ematrixd;
    typedef Eigen::Map<evec3d> map_evec3d;
    typedef Eigen::Map<const evecd> cmap_evecd;
    typedef Eigen::Map<const evec3d> cmap_evec3d;
  } // namespace eigen
} // namespace fxz


#endif
