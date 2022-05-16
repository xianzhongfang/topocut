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



#include "common.h"
#include "grid.h"


namespace fxz {
namespace grid {

int write_ply(
    const std::string& out_file,    
    const std::vector<double>& V,
    const std::vector<size_t>& M);


int read_poly_vtk(
    const std::string& in_poly,
    std::vector<double>& V,
    std::vector<size_t>& C, size_t& cell_num);


int polyvtk_to_ply(
      const std::string& in_poly, const std::string& out_ply)
{

  std::vector<double> V;
  std::vector<size_t> C;
  size_t cell_num(0);
  if (read_poly_vtk(in_poly, V, C, cell_num)) {return 1;}

  std::set<std::tuple<size_t,size_t,size_t> > tris;

  size_t id=0;
  for (size_t ci = 0; ci < cell_num; ++ci) {
    //size_t cell_size = C[id];
    ++id;
    size_t fn = C[id];
    ++id;
    for (size_t fi = 0; fi < fn; ++fi) {
      size_t f_vn = C[id];
      std::vector<size_t> f_vs;
      f_vs.reserve(3);
      ++id;
      for (size_t vi = 0; vi < f_vn; ++vi) {
        f_vs.push_back(C[id+vi]);
      }
      id+=f_vn;
      if (f_vn != 3) {
        std::cerr << "# [ ERROR ] the polyhedral mesh has non triangle face." << std::endl;
        return 1;
      }
      std::sort(f_vs.begin(), f_vs.end());
      tris.insert(std::make_tuple(f_vs[0], f_vs[1], f_vs[2]));
    }
  }
  
  std::vector<size_t> all_tris;
  all_tris.reserve(tris.size()*3);
  for (auto& t : tris) {
    all_tris.push_back(std::get<0>(t));
    all_tris.push_back(std::get<1>(t));
    all_tris.push_back(std::get<2>(t));
  }

  std::cout << "-- start to write ply ..." << std::endl;
  if (write_ply(out_ply, V, all_tris)) {return 1;}
  std::cout << "-- write into file: " << out_ply << std::endl;

  return 0;
}


int polyvtk_to_ply(const pm_map& para)
{
  const std::string& in_poly = para["input"].as<std::string>();
  const std::string& out_ply = para["output"].as<std::string>();

  return polyvtk_to_ply(in_poly, out_ply);
}

int polyvtk_to_ply(int argc, char* argv[])
{
  using namespace boost::program_options;

  try {
    options_description console_desc("TopoCut Console: transform polyhedral mesh (vtk) into triangles as ply format");
    console_desc.add_options()
        ("help,h", "help screen")
        ("input,i", value<std::string>()->notifier(on_input), "the input polyhedral vtk")
        ("output,o", value<std::string>()->notifier(on_output), "the output triangles ply");
    
    fxz::grid::pm_map vm;
    store(parse_command_line(argc, argv, console_desc), vm);
    notify(vm);
    if (vm.count("help")) { std::cout << console_desc << std::endl; }
    else {
      CALL_FUNC(fxz::grid::polyvtk_to_ply(vm));
    }
  } catch (const error &ex) {
    std::cerr << ex.what() << std::endl;
  }
  
  return 0;
}



} // grid
} // fxz

