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


#ifndef GRID_FXZ_H
#define GRID_FXZ_H

#include <iostream>
#include <boost/program_options.hpp>
#include "seq.h"

namespace fxz {
namespace grid {

typedef boost::program_options::variables_map pm_map;

int make_cut_plane(const pm_map& para);
int make_cut_plane(int argc, char* argv[]);

double ave_edge_length(
    const std::vector<double>& V,
    const std::map<std::pair<size_t,size_t>, std::vector<size_t> >& e2f);

int check_plane(
    const std::vector<double>& P, const double norm_eps, const double len_eps,
    std::vector<size_t>& repeated_p);

int check_delete_repeated_plane(
    std::vector<double>& P, const double eps);

int general_exact_cut(const pm_map& para);
int general_exact_cut(int argc, char* argv[]);

int general_exact_cut(
    const pm_map& para,
    const std::string& out_file,
    const std::vector<double>& V,
    const std::vector<size_t>& M,
    std::vector<double>& P,
    std::vector<double>& CV,
    std::vector<size_t>& CC, size_t& cn,
    std::vector<double>& new_CV, double& sum_vol);

int triangulate_cut_mesh(
    const std::string& out_ply,
    size_t ori_vn, double all_vol,
    const std::vector<double>& V,
    const std::vector<size_t>& C,
    size_t cn);


int read_mesh(
    const std::string& in_mesh,
    std::vector<double>& V,
    std::vector<size_t>& M);

int read_vector(
    const std::string& file,
    std::vector<double>& P);

int write_obj(
    const std::vector<double>& V,
    const std::vector<std::vector<size_t> >& M,
    const std::string& out_file);


int write_ply(
    const std::vector<double>& V,
    const std::vector<std::vector<size_t> >& M,
    const std::string& out_file);


int write_vector(
    const std::vector<double>& P, size_t nc,
    const std::string& file);

int write_cutcell_to_vtk(
    const std::string& file,
    const std::vector<double>& CV,
    const std::vector<size_t>& CC, size_t cell_num);


int build_edge_to_triangle_map(
    const std::vector<size_t>& M,
    std::map<std::pair<size_t,size_t>, std::vector<size_t> >& e2f);

int check_input(const pm_map& para);
int check_input(int argc, char* argv[]);

int check_degenerated_triangle(
    const std::vector<double>& V,
    const std::vector<size_t>& M, const double eps);

int check_volume(
    const std::vector<double>& V,
    const std::vector<size_t>& M, const double eps, double& sum_vol);

double check_volume(
    const std::vector<double>& V,
    const std::vector<size_t>& M);


int check_topology(
    const std::vector<size_t>& M,
    const std::map<std::pair<size_t,size_t>, std::vector<size_t> >& e2f);

int check_topology(const std::vector<size_t>& M);



inline void on_input(const std::string& str)
{
  std::cout << "Input mesh: " << str << std::endl;
}

inline void on_output(const std::string& str)
{
  std::cout << "Ouput mesh: " << str << std::endl;
}

inline void on_plane(const std::string& str)
{
  std::cout << "Input plane: " << str << std::endl;
}

int generate_cut_cells(
    const pm_map& para,
    const sequence_type& seq, const std::string& out_file,
    std::vector<size_t>& CC, size_t& cn,
    std::vector<double>& new_CV,
    std::vector<std::vector<size_t> >& face_triangles);

} //
} // fxz


#endif

