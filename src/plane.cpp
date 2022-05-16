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
#include <random>

#include "common.h"
#include "grid.h"


using namespace boost::program_options;

namespace fxz {
namespace grid {

int get_plane_by_random(
    const std::vector<double>& V,
    const std::vector<size_t>& M,
    const size_t pn,
    std::vector<double>& P)
{
  P.reserve(4*pn);

  size_t vn = V.size()/3;

  double L[3], R[3];
  for (size_t j = 0; j < 3; ++j) {
    L[j] = R[j] = V[j]; // initialize
  }
  
  for (size_t vi = 0; vi < vn; ++vi) {
    for (size_t j = 0; j < 3; ++j) {
      L[j] = std::min(L[j], V[3*vi+j]);
      R[j] = std::max(R[j], V[3*vi+j]);
    }
  }
  
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> uni_float(-1.0, 1.0);

  double eps = 1e-10;
  
  for (size_t pi = 0; pi < pn; ++pi) {
    double N[3];
    do {
      for (size_t j = 0; j < 3; ++j) {
        N[j] = uni_float(mt);
      }
    } while (N[0]*N[0]+N[1]*N[1]+N[2]*N[2] < eps);
    
    double p[3];
    double k = (uni_float(mt)+1.0)*0.5;
    for (size_t j = 0; j < 3; ++j) {
      p[j] = L[j] + (R[j]-L[j])*k;
    }

    double d = -(N[0]*p[0]+N[1]*p[1]+N[2]*p[2]);

    P.push_back(N[0]);
    P.push_back(N[1]);
    P.push_back(N[2]);
    P.push_back(d);
  }

  std::cout << "-- generate random plane number: " << pn << std::endl;
  
  return 0;
}

int get_plane_by_bound_box(
    const std::vector<double>& V,
    const std::vector<size_t>& M,
    size_t ele_n,
    std::vector<double>& P)
{
  if (ele_n == 0) {
    std::cerr << "# [ ERROR ] get plane by bounding box: please input right ele_n: "
              << ele_n << std::endl;
    return 1;
  }
  
  using namespace fxz::eigen;
  cmap_ematrixd mat_V(&V[0], 3, V.size()/3);
  FASSERT(mat_V(0,1)==V[3]);

  auto p_max = mat_V.rowwise().maxCoeff();
  auto p_min = mat_V.rowwise().minCoeff();

  auto d = p_max-p_min;

  double len = pow(d[0]*d[1]*d[2]/(double)ele_n, 1.0/3.0);

  size_t n[3];
  for (size_t i = 0; i < 3; ++i) {
    n[i] = (size_t)(d[i]/len);
    if (n[i]==0) n[i] = 1;
  }

  for (size_t i = 0; i < 3; ++i) {
    evec3d N;
    N[0] = N[1] = N[2] = 0.0;
    N[i] = 1.0;

    auto dx = d/(double)(n[i]);
    
    for (size_t j = 0; j <= n[i]; ++j) {
      double d_val = -N.dot(p_min+dx*(double)j);
      if (j==n[i]) {
        d_val = -N.dot(p_max);
      }
      P.push_back(N[0]);
      P.push_back(N[1]);
      P.push_back(N[2]);
      P.push_back(d_val);
    }
  }
  
  return 0;
}


int get_plane_by_bound_box_xyz(
    const std::vector<double>& V,
    const std::vector<size_t>& M,
    size_t nx, size_t ny, size_t nz,
    std::vector<double>& P)
{
  if (nx < 1 || ny < 1 || nz < 1) {
    std::cerr << "# [ ERROR ] get plane by bounding box: please input right nx, ny, nz (>1): "
              << nx << ", " << ny << ", " << nz << std::endl;
    return 1;
  }
  
  using namespace fxz::eigen;
  cmap_ematrixd mat_V(&V[0], 3, V.size()/3);
  FASSERT(mat_V(0,1)==V[3]);

  auto p_max = mat_V.rowwise().maxCoeff();
  auto p_min = mat_V.rowwise().minCoeff();
  auto d = p_max-p_min;

  size_t n[3] = {nx, ny, nz};

  for (size_t i = 0; i < 3; ++i) {
    evec3d N;
    N[0] = N[1] = N[2] = 0.0;
    N[i] = 1.0;

    auto dx = d/(double)(n[i]-1);
    
    for (size_t j = 0; j < n[i]; ++j) {
      double d_val = -N.dot(p_min+dx*(double)j);
      if (j==n[i]) {
        d_val = -N.dot(p_max);
      }
      P.push_back(N[0]);
      P.push_back(N[1]);
      P.push_back(N[2]);
      P.push_back(d_val);
    }
  }
  
  return 0;
}


int make_cut_plane(const pm_map& para)
{
  const std::string& in_mesh = para["input"].as<std::string>();
  const std::string& strategy = para["strategy"].as<std::string>();

  const std::string& out_plane = para["output"].as<std::string>();
  
  std::vector<double> V, P; // V: 3*vn, P: 4*pn.
  std::vector<size_t> M; // M: 3*fn.
  CALL_FUNC(read_mesh(in_mesh, V, M));

  if (strategy == "bound_box") {
    const size_t ele_n = para["num"].as<size_t>();
    std::cout << "-- ele num: " << ele_n << std::endl;
    get_plane_by_bound_box(V, M, ele_n, P);
  } else if (strategy == "random") {
    const size_t pn = para["planenum"].as<size_t>();
    std::cout << "-- plane num: " << pn << std::endl;
    get_plane_by_random(V, M, pn, P);
  } else if (strategy == "bound_box_xyz") {
    const size_t nx = para["nx"].as<size_t>();
    const size_t ny = para["ny"].as<size_t>();
    const size_t nz = para["nz"].as<size_t>();
    std::cout << "-- nx, ny, nz: " << nx << ", " << ny << ", " << nz << std::endl;
    get_plane_by_bound_box_xyz(V, M, nx, ny, nz, P);
  } else {
    std::cout << "# There is no such strategy for plane making: "
              << strategy << std::endl;
    return 0;
  }

  CALL_FUNC(write_vector(P, 4, out_plane));
  
  return 0;
}

void on_strategy(const std::string& str)
{
  std::cout << "Making planes, strategy: " << str << std::endl;
}

void on_element_num(size_t ele_n)
{
  std::cout << "Making planes, element number : " << ele_n << std::endl;
}

int make_cut_plane(int argc, char* argv[])
{
  try {
    options_description console_desc("TopoCut Console Setting");
    console_desc.add_options()
        ("help,h", "help screen")
        ("input,i", value<std::string>()->notifier(on_input), "the input mesh")
        ("strategy,s", value<std::string>()->default_value("bound_box"), "strategy for making: bound_box, bound_box_xyz, random, ...")
        ("num,n", value<size_t>()->default_value(100), "expected element number for bounding box strtegy")
        ("nx", value<size_t>()->default_value(10), "expected x element number for bounding box xyz strategy")
        ("ny", value<size_t>()->default_value(10), "expected y element number for bounding box xyz strategy")
        ("nz", value<size_t>()->default_value(10), "expected z element number for bounding box xyz strategy")
        ("planenum,P", value<size_t>()->default_value(20), "expected plane number for random strategy")        
        ("output,o", value<std::string>()->notifier(on_output), "the output plane");
    
    fxz::grid::pm_map vm;
    store(parse_command_line(argc, argv, console_desc), vm);
    
    notify(vm);

    if (vm.count("help")) { std::cout << console_desc << std::endl; }

    CALL_FUNC(fxz::grid::make_cut_plane(vm));
    
  } catch (const error &ex) {
    std::cerr << ex.what() << std::endl;
  }

  return 0;
}

} // grid
} // fxz
