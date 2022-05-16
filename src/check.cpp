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
#include <queue>
#include <list>
#include <iostream>
#include <sstream>
#include <fstream>

#include "common.h"
#include "grid.h"
#include "my_gmp.h"

using namespace boost::program_options;

namespace fxz {
namespace grid {


int check_topology(
    const std::vector<size_t>& M,
    const std::map<std::pair<size_t,size_t>, std::vector<size_t> >& e2f)
{

  for (auto& e : e2f) {
    if ((e.second.size()&1) == 1) {
      std::cerr << "-- [ ERROR ] The mesh is not watertight." << std::endl;
      return 1;
    } else {
      size_t ln(0), rn(0);
      for (auto f : e.second) {
        for (size_t j = 0; j < 3; ++j) {
          size_t va = M[3*f+j];
          size_t vb = M[3*f+(j+1)%3];

          if (va == e.first.first && vb == e.first.second) {
            ++ln;
          } else if (vb == e.first.first && va == e.first.second) {
            ++rn;
          }
        }
      } // for
      if (ln != rn) {
        std::cerr << "-- [ ERROR ] The mesh is not watertight, non-manifold inner edge." << std::endl;
        return 1;
      }
    } // if
  } // for
  return 0;
}

void regularize_plane(double* P, const double norm_eps, const double len_eps)
{
  //assert(eps > 1e-15);
  
  for (size_t i = 0; i < 4; ++i) {
    if ((i < 3 && fabs(P[i]) > norm_eps) || (i==3 && fabs(P[i]) > len_eps)) {
      if (P[i]<0) {
        for (size_t j = i; j < 4; ++j) {
          P[j] *= -1.0;
        }
      }
      break;
    } else {
      P[i] = 0.0;
    }
  }
}

struct plane_type {
 public:
  plane_type(const double* p) {
    for (size_t i=0; i<4; ++i) {
      P_[i] = p[i];
    }
  }
  double operator[] (size_t id) const { return P_[id];}
  double P_[4];
};

double NORM_CMP_EPS;
double LEN_CMP_EPS;

class myc {
 public:
  bool operator()(const plane_type& A, const plane_type& B) const
  {
    if (fabs(A[0]-B[0]) > NORM_CMP_EPS) {
      return (A[0] < B[0]);
    } else if (fabs(A[1]-B[1]) > NORM_CMP_EPS) {
      return (A[1] < B[1]);
    } else if (fabs(A[2]-B[2]) > NORM_CMP_EPS) {
      return (A[2] < B[2]);
    } else if (fabs(A[3]-B[3]) > LEN_CMP_EPS) {
      return (A[3] < B[3]);
    } else {
      return false; // need checking
    }
  }
};

int check_plane(
    const std::vector<double>& P, double norm_eps, double len_eps,
    std::vector<size_t>& invalid_p)
{
  // regularize planes
  std::vector<double> new_P(P);

  size_t pn = P.size()/4;

  {
    if (len_eps < 1e-15) {
      len_eps = 1e-15;
      std::cout << "-- reset len eps for plane checking: len_eps=1e-15" << std::endl;
    }
    if (norm_eps < 1e-15) {
      norm_eps = 1e-15;
      std::cout << "-- reset norm eps for plane checking: norm_eps=1e-15" << std::endl;
    }
  }

  std::vector<char> flag_p(pn, 1);
  for (size_t i = 0; i < pn; ++i) {
    double len = sqrt(P[4*i]*P[4*i]+P[4*i+1]*P[4*i+1]+P[4*i+2]*P[4*i+2]);
    if (len < norm_eps) {
      std::cout << "# [ WARNING ] there exists one plane with degenerated norm: "
                << i << std::endl;
      flag_p[i] = 0; // invalid
    } else {
      for (size_t j = 0; j < 4; ++j) {
        new_P[4*i+j] = P[4*i+j]/len;
      }
    }
  }

  for (size_t pi = 0; pi < pn; ++pi) {
    regularize_plane(&new_P[4*pi], norm_eps, len_eps);
  }

  NORM_CMP_EPS=norm_eps;
  LEN_CMP_EPS=len_eps;
  
  std::set<plane_type, myc> plane_set;

  size_t cnt = 0;

  for (size_t pi = 0; pi < pn; ++pi) {

    if (flag_p[pi] == 0) {
      invalid_p.push_back(pi);
      continue;
    }
    
    auto pt(plane_type(&new_P[4*pi]));
    auto it = plane_set.find(pt);

    if (it == plane_set.end()) {
      plane_set.insert(pt);
    } else {
      std::cout << "-- [ WARNING ] find repeated plane, id: "
                << pi << std::endl;
      invalid_p.push_back(pi);
      ++cnt;
    }
  }

  if (cnt > 0) {
    std::cout << "-- [ WARNING ] find " << cnt
              << " repeated planes." << std::endl;
    return 1;
  }
  
  return 0;
}


int check_topology(const std::vector<size_t>& M)
{
  std::map<std::pair<size_t,size_t>, std::vector<size_t> > e2f;
  build_edge_to_triangle_map(M, e2f);
  return check_topology(M, e2f);
}

int check_topology(const pm_map& para)
{
  const std::string in_mesh = para["input"].as<std::string>();

  std::vector<double> V;
  std::vector<size_t> M;
  CALL_FUNC(read_mesh(in_mesh, V, M));
  CALL_FUNC(check_topology(M));
  
  return 0;
}


double check_volume(
    const std::vector<double>& V,
    const std::vector<size_t>& M)
{
  mpf_class sum_vol(0.0);
  size_t tn = M.size()/3;

  for (size_t ti = 0; ti < tn; ++ti) {
    size_t va = M[3*ti];
    size_t vb = M[3*ti+1];
    size_t vc = M[3*ti+2];
    Eigen::Matrix<mpf_class,3,1> A, B, C;
    A << V[3*va], V[3*va+1], V[3*va+2];
    B << V[3*vb], V[3*vb+1], V[3*vb+2];
    C << V[3*vc], V[3*vc+1], V[3*vc+2];
    sum_vol += A.cross(B).dot(C);
  }
  return sum_vol.get_d();
}

int check_volume(
    const std::vector<double>& V,
    const std::vector<size_t>& M, const double eps, double& sum_vol)
{
  sum_vol = check_volume(V, M);
  if (sum_vol < eps) {
    std::cout << "-- [ ERROR ] The volume of the mesh is smaller than " << eps << std::endl;
    return 1;
  }
  std::cout << "-- The mesh has positive volume." << std::endl;

  return 0;
}

double ave_edge_length(
    const std::vector<double>& V,
    const std::map<std::pair<size_t,size_t>, std::vector<size_t> >& e2f)
{
  typedef Eigen::Matrix<double,3,1> evec3d;
  typedef Eigen::Map<const evec3d> cmap_evec3d;

  double sum_len(0.0);
  
  for (auto& e2 : e2f) {
    size_t va = e2.first.first;
    size_t vb = e2.first.second;
    cmap_evec3d A(&V[3*va], 3, 1);
    cmap_evec3d B(&V[3*vb], 3, 1);
    double len = (A-B).norm();
    sum_len += len;
  }

  return sum_len/(size_t)e2f.size();
}

int check_degenerated_edge(
    const std::vector<double>& V,
    const std::map<std::pair<size_t,size_t>, std::vector<size_t> >& e2f,
    const double eps)
{
  typedef Eigen::Matrix<double,3,1> evec3d;
  typedef Eigen::Map<const evec3d> cmap_evec3d;

  size_t cnt(0);
  
  for (auto& e2 : e2f) {
    size_t va = e2.first.first;
    size_t vb = e2.first.second;

    cmap_evec3d A(&V[3*va], 3, 1);
    cmap_evec3d B(&V[3*vb], 3, 1);

    double len = (A-B).norm();

    if (len < eps) {
      std::cout << "-- [ ERROR ] degenerated edge: (" << va << ", " << vb
                << "), length is smaller than " << eps << std::endl;
      ++cnt;
    }
  }

  std::cout << "-- degenerated edge number: " << cnt << std::endl;

  if (cnt > 0) return 1;
  return 0;
}

int check_degenerated_triangle(
    const std::vector<double>& V,
    const std::vector<size_t>& M, const double eps)
{
  typedef Eigen::Matrix<double,3,1> evec3d;
  typedef Eigen::Map<const evec3d> cmap_evec3d;
  
  size_t tn = M.size()/3;

  size_t cnt(0);

  for (size_t ti = 0; ti < tn; ++ti) {
    size_t va = M[3*ti];
    size_t vb = M[3*ti+1];
    size_t vc = M[3*ti+2];

    cmap_evec3d A(&V[3*va], 3, 1);
    cmap_evec3d B(&V[3*vb], 3, 1);
    cmap_evec3d C(&V[3*vc], 3, 1);
    
    double vol = (B-A).cross(C-A).norm();

    if (vol < eps) {
      std::cout << "-- [ ERROR ] degenerated triangle: "
                << "(" << va << ", " << vb << ", " << vc << "), volume is smaller than "
                << eps << std::endl;
      ++cnt;
    }
  }

  std::cout << "-- degenerated triangle number: " << cnt << std::endl;

  if (cnt > 0) return 1;
  return 0;
}

int check_input(const pm_map& para)
{
  const std::string in_mesh = para["input"].as<std::string>();
  double      eps = para["eps"].as<double>();

  std::cout << "-- mesh: " << in_mesh << std::endl;
  std::cout << "--  eps: " << eps << std::endl;

  std::vector<double> V;
  std::vector<size_t> M;
  CALL_FUNC(read_mesh(in_mesh, V, M));

  if (M.size()%3 != 0) {
    std::cout << "-- [ ERROR ] The mesh is not triangular." << std::endl;
    return 1;
  }

  // check topology
  std::map<std::pair<size_t,size_t>, std::vector<size_t> > e2f;
  CALL_FUNC(build_edge_to_triangle_map(M, e2f));
  CALL_FUNC(check_topology(M, e2f));

  if (e2f.size() == 0) {
    std::cout << "-- [ ERROR ] There is no triangle in the mesh." << std::endl;
    return 1;
  }

  double ave_len = ave_edge_length(V, e2f);

  eps *= ave_len;

  std::cout << "-- average length of   edges: " << ave_len << std::endl;
  std::cout << "-- new eps considering edges: " << eps << std::endl;

  double all_vol(0.0);
  CALL_FUNC(check_volume(V, M, eps, all_vol));

  CALL_FUNC(check_degenerated_triangle(V, M, eps));

  return 0;
}

int check_input(int argc, char* argv[])
{
  try {
    options_description console_desc("TopoCcut Console Setting");
    console_desc.add_options()
        ("help,h", "help screen")
        ("input,i", value<std::string>()->notifier(on_input), "the input mesh")
        ("eps,e", value<double>()->default_value(1e-15), "eps for degeneration checking");
    
    fxz::grid::pm_map vm;
    store(parse_command_line(argc, argv, console_desc), vm);
    notify(vm);
    if (vm.count("help")) { std::cout << console_desc << std::endl; }
    if (fxz::grid::check_input(vm)) {
      std::cout << "The mesh is not valid for cutting." << std::endl;
    } else {
      std::cout << "# [ MESSAGE ] The mesh is valid." << std::endl;
      std::cout << "OK" << std::endl;
    }
  } catch (const error &ex) {
    std::cerr << ex.what() << std::endl;
  }

  return 0;
}


} // grid
} // fxz
