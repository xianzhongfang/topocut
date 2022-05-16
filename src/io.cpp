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


#include <string>
#include <fstream>
#include <iostream>
#include <unordered_set>
#include <iomanip>
#include "common.h"
#include "grid.h"

namespace fxz {
namespace grid {


/// a simple obj parser. only read vertices coordinates and faces
int read_obj(
    const std::string& file_name,
    std::vector<double>& verts,
    std::vector<size_t>& mesh,
    std::vector<size_t>& mesh_ptr)
{
  std::cout << "# [ Start to read OBJ file ]" << std::endl;
  std::ifstream inf(file_name.c_str());

  if (inf.fail()) {
    std::cerr << "# [ Error: read obj ] Can't open file: "
        << file_name << std::endl;
    return __LINE__;
  }

  std::vector<double>& V  = verts;
  std::vector<size_t>& M  = mesh;
  std::vector<size_t>& MP = mesh_ptr;

  std::string str_line;
  size_t cnt = 0;
  while (std::getline(inf, str_line)) {
    std::string ty;
    std::istringstream ss(str_line);
    ss >> ty;
    if (ty[0] == '#') { continue; }
    else if (ty == "v") {
      double x, y, z;
      ss >> x >> y >> z;
      V.push_back(x); V.push_back(y); V.push_back(z);
    } else if ( ty == "f") {
      std::string fid_str;
      //mesh_ptr.push_back(mesh.size());
      MP.push_back(cnt);
      while (ss >> fid_str) {
        std::istringstream fid_ss(fid_str);
        size_t fid;
        fid_ss >> fid;
        ++cnt;
        M.push_back(fid-1);
      }
    } else if ( ty == "l") {
      std::string fid_str;
      //mesh_ptr.push_back(mesh.size());
      MP.push_back(cnt);
      while (ss >> fid_str) {
        std::istringstream fid_ss(fid_str);
        size_t fid;
        fid_ss >> fid;
        ++cnt;
        M.push_back(fid-1);
      }
    }
  }

  MP.push_back(cnt);

  if (cnt != M.size()) {
    std::cerr << "-- [ ERROR: read obj ] mesh size is not compatible with faces type sum."
              << std::endl;
  }

  std::cout << "-- Verts Num: " << V.size()/3 << std::endl;
  std::cout << "-- Faces Num: " << MP.size()-1 << std::endl;
  if (MP.size() > 1) 
    std::cout << "-- First face type: "<< MP[1]-MP[0] << std::endl;
  std::cout << "# [ End: read obj.]" << std::endl;
  
  return 0;
}


int read_mesh(
    const std::string& in_mesh,
    std::vector<double>& V,
    std::vector<size_t>& M)
{
  std::vector<size_t> MP;
  return read_obj(in_mesh, V, M, MP);
}

int read_vector(
    const std::string& file,
    std::vector<double>& P)
{
  std::ifstream inf(file);

  if (inf) {
    size_t nr(0), nc(0);
    inf >> nr >> nc;
    if (nr*nc == 0) {
      std::cerr << "# [ ERROR ] the first line is not right." << std::endl;
      return 1;
    }
    P.reserve(nc*nr);
    for (size_t i = 0; i < nr; ++i) {
      for (size_t j = 0; j < nc; ++j) {
        double v;
        inf >> v;
        P.push_back(v);
      }
    }
  } else {
    std::cerr << "# [ ERROR ] can not open file: " << file << std::endl;
    return 1;
  }

  return 0;
}

int write_vector(
    const std::vector<double>& P, size_t nc,
    const std::string& file)
{
  if (nc == 0) {
    std::cout << "# [ write vector ] nc = 0." << std::endl;
    return 0;
  }
  size_t nr = P.size()/nc;

  if (nr*nc != P.size()) {
    std::cerr << "# [ write vetor ] nc is not compatible with vector." << std::endl;
    return 1;
  }
  
  std::ofstream outf(file);

  outf << nr << " " << nc << std::endl;

  for (size_t i = 0; i < nr; ++i) {
    outf << P[nc*i];
    for (size_t j = 1; j < nc; ++j) {
      outf << " " << P[nc*i+j];
    }
    outf << std::endl;
  }

  return 0;
}


int write_cutcell_to_vtk(
    const std::string& file,
    const std::vector<double>& CV,
    const std::vector<size_t>& CC, size_t cell_num)
{
  std::ofstream outf(file);
  if (outf) {
    outf << "# vtk DataFile Version 4.1" << std::endl;
    outf << "General exact cutted cells" << std::endl;
    outf << "ASCII\nDATASET UNSTRUCTURED_GRID" << std::endl;
    size_t vn = CV.size()/3;
    outf << "POINTS " << vn << " float" << std::endl;

    for (size_t i = 0; i < CV.size(); i+=3) {
      outf << std::setprecision(17) << CV[i] << " " << CV[i+1] << " " << CV[i+2] << std::endl;
    }

    outf << "CELLS " << cell_num << " " << CC.size() << std::endl;
    size_t cnt = 0;
    size_t id = 0;
    do {
      size_t n = CC[id];
      outf << n; //all size and fn
      for (size_t i = 0; i < n; ++i) {
        outf << " " << CC[id+1+i];
      }
      outf << std::endl;
      id += n+1;
      ++cnt;
    } while (id < CC.size());

    FASSERT(cnt == cell_num);

    outf << "CELL_TYPES " << cell_num << std::endl;
    for (size_t ci = 0; ci < cell_num; ++ci) {
      outf << 42 << std::endl;
    }
    
  } else {
    std::cerr << "# [ ERROR ] can not open file: " << file << std::endl;
    return 1;
  }
  

  return 0;
}

int write_obj(
    const std::vector<double>& V,
    const std::vector<std::vector<size_t> >& M,
    const std::string& out_file)
{
  std::ofstream outf(out_file);
  
  if (outf.fail()) {
    std::cerr << "-- [ ERROR ] can not open file " << out_file << std::endl;
    return 1;
  }

  for (size_t i = 0; i < V.size(); i+=3) {
    outf << std::setprecision(17) << "v " << V[i] << " " << V[i+1] << " " << V[i+2] << std::endl;
  }

  for (auto& f : M) {
    for (size_t i=0; i < f.size(); i+=3) {
      outf << "f " << f[i]+1 << " " << f[i+1]+1 << " " << f[i+2]+1 << std::endl;
    }
  }

  outf.close();

  return 0;
}


int write_ply(
    const std::vector<double>& V,
    const std::vector<std::vector<size_t> >& M,
    const std::string& out_file)
{
  std::ofstream outf(out_file);
  
  if (outf.fail()) {
    std::cerr << "-- [ ERROR ] can not open file " << out_file << std::endl;
    return 1;
  }

  size_t fn(0);
  for (auto& f : M) {
    fn += f.size()/3;
  }

  outf << "ply\n"
       << "format ascii 1.0\n"
       << "comment write ply\n"
       << "element vertex " << V.size()/3 << "\n"
       << "property double x\n"
       << "property double y\n"
       << "property double z\n"
       << "element face " << fn << "\n"
       << "property list uchar int vertex_indices\n"
       << "end_header" << std::endl;
  

  for (size_t i = 0; i < V.size(); i+=3) {
    outf << std::setprecision(17) << V[i] << " " << V[i+1] << " " << V[i+2] << std::endl;
  }

  for (auto& f : M) {
    for (size_t i=0; i < f.size(); i+=3) {
      outf << "3 " << f[i] << " " << f[i+1] << " " << f[i+2] << std::endl;
    }
  }

  outf.close();

  return 0;
}

int write_ply(
    const std::string& out_file,    
    const std::vector<double>& V,
    const std::vector<size_t>& M)
{
  std::ofstream outf(out_file);
  
  if (outf.fail()) {
    std::cerr << "-- [ ERROR ] can not open file " << out_file << std::endl;
    return 1;
  }

  size_t fn = M.size()/3;

  outf << "ply\n"
       << "format ascii 1.0\n"
       << "comment write ply\n"
       << "element vertex " << V.size()/3 << "\n"
       << "property double x\n"
       << "property double y\n"
       << "property double z\n"
       << "element face " << fn << "\n"
       << "property list uchar int vertex_indices\n"
       << "end_header" << std::endl;
  

  for (size_t i = 0; i < V.size(); i+=3) {
    outf << std::setprecision(17) << V[i]
         << " " << V[i+1] << " " << V[i+2] << std::endl;
  }

  for (size_t i=0; i < M.size(); i+=3) {
    outf << "3 " << M[i] << " " << M[i+1]
         << " " << M[i+2] << std::endl;
  }

  outf.close();

  return 0;
}

int read_poly_vtk(
    const std::string& in_poly,
    std::vector<double>& V,
    std::vector<size_t>& C, size_t& cell_num)
{
  std::ifstream inf(in_poly);

  if (inf.fail()) {
    std::cerr << "# [ ERROR ] cannot open file : "
              << in_poly << std::endl;
    return 1;
  }

  std::string str_line;

  size_t vn(0);
  size_t& cn = cell_num;

  while (std::getline(inf, str_line)) {
    std::string ty;
    std::istringstream ss(str_line);
    ss >> ty;
    if (ty == "POINTS") {
      ss >> vn;
      V.resize(3*vn);
      for (size_t vi = 0; vi < vn; ++vi) {
        std::string v_str;
        std::getline(inf, v_str);
        std::istringstream v_ss(v_str);
        v_ss >> V[3*vi] >> V[3*vi+1] >> V[3*vi+2];
      }
    } else if (ty == "CELLS") {
      size_t nn;
      ss >> cn >> nn;
      C.reserve(nn);
      for (size_t ci = 0; ci < cn; ++ci) {
        std::string c_str;
        std::getline(inf, c_str);
        std::istringstream c_ss(c_str);
        size_t cell_s;
        c_ss >> cell_s;
        C.push_back(cell_s);
        for (size_t i = 0; i < cell_s; ++i) {
          size_t id;
          c_ss >> id;
          C.push_back(id);
        }
      }
    } else if (ty == "CELL_TYPES") {
      size_t n;
      ss >> n;
      if (n != cn) {
        std::cerr << "# [ ERROR ] cell-type size is not equal to cell size." << std::endl;
        return 1;
      }
      for (size_t i = 0; i < cn; ++i) {
        size_t id;
        inf >> id;
        if (id != 42) {
          std::cerr << "# [ ERROR ] cell-type is not 42, i.e., for polhedral mesh." << std::endl;
          return 1;
        }
      }
    }
  }

  inf.close();

  return 0;
}


} // grid
} // fxz
