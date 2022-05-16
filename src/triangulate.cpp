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


#include <map>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Polygon_2.h>

#include "common.h"
#include "grid.h"
#include "my_gmp.h"

#include <Eigen/Geometry>

namespace fxz {
namespace grid {

short choose_axis(const Eigen::Matrix<double,3,1>& D);


struct FaceInfo2
{
  FaceInfo2(){}
  int nesting_level;
  bool in_domain(){  return nesting_level%2 == 1; }
};
typedef CGAL::Exact_predicates_inexact_constructions_kernel       K;
typedef CGAL::Triangulation_vertex_base_2<K>                      Vb;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2,K>    Fbb;
typedef CGAL::Constrained_triangulation_face_base_2<K,Fbb>        Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>               TDS;
typedef CGAL::Exact_predicates_tag                                Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag>  CDT;
typedef CDT::Point                                                Point;
typedef CGAL::Polygon_2<K>                                        Polygon_2;


void mark_domains(
    CDT& ct, CDT::Face_handle start, int index, std::list<CDT::Edge>& border)
{
  if(start->info().nesting_level != -1){
    return;
  }
  std::list<CDT::Face_handle> queue;
  queue.push_back(start);
  while(! queue.empty()){
    CDT::Face_handle fh = queue.front();
    queue.pop_front();
    if(fh->info().nesting_level == -1){
      fh->info().nesting_level = index;
      for(int i = 0; i < 3; i++){
        CDT::Edge e(fh,i);
        CDT::Face_handle n = fh->neighbor(i);
        if(n->info().nesting_level == -1){
          if(ct.is_constrained(e)) border.push_back(e);
          else queue.push_back(n);
        }
      }
    }
  }
}
void mark_domains(CDT& cdt)
{
  for(CDT::All_faces_iterator it = cdt.all_faces_begin(); it != cdt.all_faces_end(); ++it){
    it->info().nesting_level = -1;
  }
  std::list<CDT::Edge> border;
  mark_domains(cdt, cdt.infinite_face(), 0, border);
  while(! border.empty()){
    CDT::Edge e = border.front();
    border.pop_front();
    CDT::Face_handle n = e.first->neighbor(e.second);
    if(n->info().nesting_level == -1){
      mark_domains(cdt, n, e.first->info().nesting_level+1, border);
    }
  }
}


int triangulate_with_CGAL(
    const std::vector<double>& UV,
    std::vector<size_t>& tris)
{
  using namespace CGAL;

  CDT cdt;  

  std::map<Point, size_t> c2id;  
  {
    Polygon_2 cgal_polygon;
    for (size_t i = 0; i < UV.size(); i+=2) {
      Point P(UV[i], UV[i+1]);
      cgal_polygon.push_back(P);
      c2id.insert(std::make_pair(P, i/2));
    }
    cdt.insert_constraint(cgal_polygon.vertices_begin(), cgal_polygon.vertices_end(), true);
  }

  mark_domains(cdt);
  
  for (CDT::Finite_faces_iterator fit=cdt.finite_faces_begin();
       fit!=cdt.finite_faces_end();++fit) {
    if ( fit->info().in_domain() ) {
      for (size_t i = 0; i < 3; ++i) {
        tris.push_back(c2id[fit->vertex(i)->point()]);
      }
    }
  }

  
  return 0;
}

int check_triangles_area(
    const std::vector<double>& UV,
    const std::vector<size_t>& tris)
{
  for (size_t i = 0; i < tris.size(); i+=3) {
    size_t va(tris[i]), vb(tris[i+1]), vc(tris[i+2]);
    Eigen::Matrix<mpf_class, 2, 1> pa, pb, pc;
    pa << UV[2*va], UV[2*va+1];
    pb << UV[2*vb], UV[2*vb+1];
    pc << UV[2*vc], UV[2*vc+1];
    auto ab = pb-pa;
    auto ac = pc-pa;
    mpf_class area = ab[0]*ac[1]-ab[1]*ac[0];
    if (sgn(area)<=0) {
      std::cerr << "-- [ WARNING ] negative triangle." << std::endl;
    }
  }

  return 0;
}

int triangulate_face(
    const std::vector<double>& V,
    const std::vector<size_t>& f,
    std::vector<size_t>& tris)
{
  tris.reserve(3*(f.size()-2));
  
  if (f.size() == 3) {
    tris = f;
  }     
  else {
    Eigen::Map<const Eigen::Matrix<double, 3, 1> > O(&V[3*f[0]], 3, 1);
    Eigen::Matrix<double, 3, 1> sum_N;
    sum_N.setZero();
    for (size_t i = 1; i+1 < f.size(); ++i) {
      Eigen::Map<const Eigen::Matrix<double, 3, 1> > A(&V[3*f[i]], 3, 1);
      Eigen::Map<const Eigen::Matrix<double, 3, 1> > B(&V[3*f[i+1]], 3, 1);
      sum_N += (A-O).cross(B-O);
    }
    short direc = choose_axis(sum_N);


    size_t x = abs(direc)%3; // +1
    size_t y = (abs(direc)+1)%3;//+2
    if (direc < 0) std::swap(x, y);
    
    std::vector<double> UV(2*f.size());
    for (size_t i = 0; i < f.size(); ++i) {
      UV[2*i] = V[3*f[i]+x];
      UV[2*i+1] = V[3*f[i]+y];
    }
    
    triangulate_with_CGAL(UV, tris);
    check_triangles_area(UV, tris); //for test
    
    for (auto& v : tris) {
      v = f[v];
    }
  }

  return 0;
}


int triangulate_cut_mesh(
    size_t ori_vn, double all_vol,
    const std::vector<double>& V,
    const std::vector<size_t>& CC,
    size_t cn,
    std::vector<std::vector<size_t> >& faces_tris,
    std::vector<size_t>& tris_CC,
    size_t& new_cn)
{
  std::vector<std::vector<size_t> > all_faces; // faces <verts id list>
  std::vector<std::vector<size_t> > all_cells; // cells <faces id list>

  // <verts num, sum(verts id)> -> faces id list
  std::map<std::pair<size_t,size_t>, std::vector<size_t> > faces_key_map;
  std::vector<size_t> faces_pair; // faces pair flag;

  {
    size_t cnt(0);
    for (size_t i = 0; i < CC.size(); i+=CC[i]+1) {
      size_t fn = CC[i+1];
      cnt += fn;
    }
    all_faces.reserve(cnt);
    faces_pair.resize(cnt, INVALID_NUM);
  }

  all_cells.reserve(cn);  

  size_t face_id(0);

  for (size_t i = 0; i < CC.size(); i+=CC[i]+1) {
    size_t fn = CC[i+1];
    std::vector<size_t> fs;
    fs.reserve(fn);
    for (size_t j = i+2; j < i+CC[i]+1; j+=CC[j]+1) {
      size_t vn = CC[j];
      std::vector<size_t> f;
      f.reserve(vn);
      for (size_t k = 0; k < CC[j]; ++k) {
        f.push_back(CC[j+k+1]);
      }
      all_faces.push_back(f);

      {
        size_t vs_sum = std::accumulate(f.begin(), f.end(), 0);
        auto it = faces_key_map.find(std::make_pair(vn, vs_sum));
        if (it == faces_key_map.end()) {
          std::vector<size_t> list_fs;
          list_fs.reserve(2);
          list_fs.push_back(face_id);
          faces_key_map.insert(std::make_pair(std::make_pair(vn, vs_sum), list_fs));
        } else {
          it->second.push_back(face_id);
        }
      }
      fs.push_back(face_id);
      ++face_id;
    }

    all_cells.push_back(fs);
  } // for

  FASSERT(face_id == all_faces.size());
  FASSERT(all_cells.size() == cn);

  for (auto& s : faces_key_map) {
    const auto& fs_ids = s.second;
    if (fs_ids.size() == 1) {continue;}

    std::vector<std::string> fs_str;
    fs_str.reserve(fs_ids.size());
    for (auto& fid : fs_ids) {
      std::vector<size_t> f = all_faces[fid];
      std::sort(f.begin(), f.end());
      std::stringstream ss;
      std::copy(f.begin(), f.end(), std::ostream_iterator<size_t>(ss, ","));
      fs_str.push_back(ss.str());
    }
    for (size_t i = 0; i < fs_str.size(); ++i) {
      for (size_t j = i+1; j < fs_str.size(); ++j) {
        if (fs_str[i] == fs_str[j]) {
          faces_pair[fs_ids[i]] = fs_ids[j];
          faces_pair[fs_ids[j]] = fs_ids[i];
          break;
        }
      }
    }
  }

  faces_tris.resize(faces_pair.size());
  // triangulate each faces
#pragma omp parallel for
  for (size_t i = 0; i < faces_pair.size(); ++i) {
    std::vector<size_t>& tris = faces_tris[i];
    
    if (faces_pair[i] != INVALID_NUM) {
      FASSERT(faces_pair[faces_pair[i]] == i);
      if (faces_pair[i] > i) { // only one side is triangulated
        triangulate_face(V, all_faces[i], tris);
      }
    } else { // boundary faces, convex, directly triangulate it
      const auto& f = all_faces[i];
      tris.reserve(3*(f.size()-2));
      for (size_t j = 1; j+1 < f.size(); ++j) {
        tris.push_back(f[0]);
        tris.push_back(f[j]);
        tris.push_back(f[j+1]);
      }
    }

    FASSERT(tris.size()%3==0);
  }

  FASSERT(faces_tris.size() == faces_pair.size());

  
  new_cn = cn;
  for (size_t ci = 0; ci < cn; ++ci) {
    size_t tris_fn(0);
    for (auto fid : all_cells[ci]) {
      const auto& f = faces_tris[fid].size()>0 ? faces_tris[fid] : faces_tris[faces_pair[fid]];
      tris_fn += f.size()/3;
      FASSERT(f.size()%3 == 0);
    }
    
    tris_CC.push_back(1+tris_fn+3*tris_fn);
    tris_CC.push_back(tris_fn);

    std::vector<size_t> tris_M;

    bool is_bd = false;
    bool is_connect_inner=false;
    
    for (auto fid : all_cells[ci]) {
      FASSERT(fid < faces_pair.size());
      if (faces_pair[fid] == INVALID_NUM) {
        is_bd = true;
      } else {
        is_connect_inner = true;
      }

      const auto& f = faces_tris[fid].size()>0 ? faces_tris[fid] : faces_tris[faces_pair[fid]];

      FASSERT(f.size() >= 3 && f.size()%3 == 0);
      short direc = faces_tris[fid].size()>0 ? 1 : -1;
      for (size_t i = 0; i < f.size(); i+=3) {
        tris_CC.push_back(3);
        tris_CC.push_back(f[i]);
        
        tris_M.push_back(f[i]);
        
        if (direc > 0) {
          tris_CC.push_back(f[i+1]);
          tris_CC.push_back(f[i+2]);
          
          tris_M.push_back(f[i+1]);
          tris_M.push_back(f[i+2]);
          
        } else {
          tris_CC.push_back(f[i+2]);
          tris_CC.push_back(f[i+1]);
          
          tris_M.push_back(f[i+2]);
          tris_M.push_back(f[i+1]);
        }
      }
    }// for all_cells

    auto r_topo = check_topology(tris_M);
    double cell_vol = check_volume(V, tris_M);
    if (is_bd && is_connect_inner && r_topo==0 && (cell_vol>=0 && (cell_vol) <= 1e-6*all_vol)) {
      size_t all_n = 2+4*tris_fn;
      for (size_t i = 0; i < all_n; ++i) {
        tris_CC.pop_back();
      }
      --new_cn;
    }
  }//for ci

  return 0;
}

int triangulate_cut_mesh(
    const std::string& out_mesh,
    size_t ori_vn, double all_vol,
    const std::vector<double>& V,
    const std::vector<size_t>& CC,
    size_t cn)
{
  std::vector<std::vector<size_t> > faces_tris;
  std::vector<size_t> tris_CC;
  size_t new_cn;

  CALL_FUNC_WITH_CLOCK(
      triangulate_cut_mesh(
          ori_vn, all_vol, V, CC, cn,
          faces_tris, tris_CC, new_cn),
      "triangulate-cut-mesh");

  // CALL_FUNC(write_ply(V, faces_tris, out_mesh+"-tri.ply"));
  // std::cout << "-- write triangles ply into file: " << out_mesh+"-tri.ply" << std::endl;

  CALL_FUNC(write_cutcell_to_vtk(
      out_mesh+"-tri.vtk", V, tris_CC, new_cn));
  std::cout << "-- write triangulated cut-mesh into file: " << out_mesh+"-tri.vtk" << std::endl;

  
  return 0;
}

} // grid
} // fxz
