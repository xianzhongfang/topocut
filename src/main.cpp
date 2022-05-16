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


#include <iostream>
#include <string>
#include <fstream>
#include "grid.h"
#include "common.h"

namespace fxz {
  namespace grid {
  // transform all triangles in polyhedral mesh into ply
  int polyvtk_to_ply(int argc, char* argv[]);
  }
} // fxz

int main(int argc, char* argv[])
{
  if (argc == 1) {
    std::cout << " There is no any parameter." << std::endl;
    std::cout << "# sub-program: plane, cut" << std::endl;
    return 0;
    
  } else if (argc > 1) {
    
    std::string prog = argv[1];
    std::cout << "# main prog: \n" << argv[0] << std::endl;
    std::cout << "# sub  prog: " << argv[1] << std::endl;
    
    CALL_SUBPROG("plane", fxz::grid::make_cut_plane);
    CALL_SUBPROG("cut", fxz::grid::general_exact_cut);
    CALL_SUBPROG("check", fxz::grid::check_input);
    CALL_SUBPROG("poly2ply", fxz::grid::polyvtk_to_ply);
  }

  return 0;
}
