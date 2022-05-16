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


#include "grid.h"
#include "my_gmp.h"


using namespace boost::program_options;

namespace fxz {
namespace grid {


int general_exact_cut(int argc, char* argv[])
{
  try {
    options_description console_desc("TopoCut Console Setting");
    console_desc.add_options()
        ("help,h", "help screen")
        ("input,i", value<std::string>()->notifier(on_input), "the input mesh")
        ("plane,p", value<std::string>()->notifier(on_plane), "the input planes")
        ("is-triangulate", value<size_t>()->default_value(0), "whether to triangulate cut-faces (0 or 1)")
        ("is-round-and-embed", value<size_t>()->default_value(1), "whether to do rounding and embedding (0 or 1)")
        ("is-connect-hole", value<size_t>()->default_value(1), "whether to connect holes in cut-planes (0 or 1)")
        ("eps,e", value<double>()->default_value(1e-20), "used for detecting degeneration")
        ("line_sort_eps_ratio", value<double>()->default_value(1e-3), "used for coarse line-sorting")
        ("planeeps,E", value<double>()->default_value(1e-2), "used for detecting repeated plane")
        ("detecplane", value<size_t>()->default_value(1), "detec input planes (0 or 1)")
        ("output,o", value<std::string>()->notifier(on_output), "the output mesh");
    
    fxz::grid::pm_map vm;
    store(parse_command_line(argc, argv, console_desc), vm);
    
    notify(vm);

    if (vm.count("help")) { std::cout << console_desc << std::endl; }
    else {
      CALL_FUNC(fxz::grid::general_exact_cut(vm));
    }
    
  } catch (const error &ex) {
    std::cerr << ex.what() << std::endl;
  }

  return 0;
}


void reduce_precision(double& x, size_t p)
{
  std::stringstream ss;
  ss << std::setprecision(p) << x;
  ss >> x;
}

int reduce_precision(std::vector<double>& V, size_t p)
{
  for (auto& x : V) {
    reduce_precision(x, p);
  }

  return 0;
}

int general_exact_cut(const pm_map& para)
{
  // read input triangular mesh
  const std::string& in_mesh = para["input"].as<std::string>();
  const std::string& in_plane = para["plane"].as<std::string>();
  const std::string& out_mesh = para["output"].as<std::string>();
  
  std::vector<double> V, P; // V: 3*vn, P: 4*pn.
  std::vector<size_t> M; // M: 3*fn.
  CALL_FUNC(read_mesh(in_mesh, V, M));
  CALL_FUNC(read_vector(in_plane, P));
  
  std::cout << "-- plane number: " << P.size()/4 << std::endl;

  if (P.size() == 0) {
    std::cout << "-- [ WARNING ] there is no input plane." << std::endl;
    return 2;
  }

  if (M.size() == 0) {
    std::cout << "-- [ WARNING ] there is no input triangle."
              << std::endl;
    return 3;
  }

  
  std::vector<double> CV, new_CV;
  std::vector<size_t> CC;
  size_t cn(0);
  double all_vol(0.0);
  CALL_FUNC_WITH_CLOCK(
      general_exact_cut(para, out_mesh, V, M, P, CV, CC, cn, new_CV, all_vol),
      "generate-cut-mesh");

  if (cn > 0) {
    CALL_FUNC(write_cutcell_to_vtk(out_mesh, new_CV, CC, cn));
    std::cout << "-- write out mesh into file: " << out_mesh << std::endl;

    const size_t is_triangulate = para["is-triangulate"].as<size_t>();
    if (is_triangulate == 1) {
      CALL_FUNC(triangulate_cut_mesh(
          out_mesh, V.size()/3, all_vol, new_CV, CC, cn));
    }
    
  } else {
    std::cout << "-- [ WARNING ] there is no output cut-cell." << std::endl;
    return 1;
  }
  
  return 0;
}

///////////////////////////////////////////////////////////////////
//// general exact cut
///////////////////////////////////////////////////////////////////

int generate_PPP(
    const table_relation& tr,
    std::list<cut_vert>& cut_verts)
{
  size_t tn = tr.tn();
  size_t pn = tr.pn();
  auto& P = tr.plane();

  std::vector<std::list<cut_vert> > plane_cv(pn);
  
  /// PPP
#pragma omp parallel for
  for (size_t pi = 0; pi < pn; ++pi) {//
    auto& pcv = plane_cv[pi];
    
    for (size_t pj = pi+1; pj < pn; ++pj) {//
      for (size_t pk = pj+1; pk < pn; ++pk) {//
        if (tr.pp(pi,pj)==2 && tr.pp(pi,pk)==2
            && tr.pp(pj,pk)==2) {
          int r = vol_sign_with_norms(
              &P[4*pi], &P[4*pj], &P[4*pk]);

          if (r != 0) {
            pcv.push_back(cut_vert(tn+pi,tn+pj,tn+pk));
            pcv.back().set_on_grid();
          } // if
        } // if
      } // for pk
    } // for pj
  } // for pi


  for (auto& cv : plane_cv) {
    cut_verts.insert(cut_verts.end(), cv.begin(), cv.end());
  }
  
  return 0;
}

int generate_TTP(
    const table_relation& tr,
    const std::map<std::pair<size_t,size_t>, std::vector<size_t> >& e2f,
    std::list<cut_vert>& cut_verts)    
{
  size_t pn = tr.pn();
  size_t tn = tr.tn();

  { // simple method
    size_t eid = 0;
    for (const auto& e : e2f) {
      size_t va = e.first.first;
      size_t vb = e.first.second;
      for (size_t pi = 0; pi < pn; ++pi) {
        if (tr.pe(pi, va, vb) == 1) {
          for (size_t k = 0; k+1 < e.second.size(); k+=2) {
            size_t fa = e.second[k];
            size_t fb = e.second[k+1];
            if (fa > fb) std::swap(fa, fb);
            cut_verts.push_back(cut_vert(fa,fb,pi+tn));//TTP
            cut_verts.back().set_on_edge(eid);
          }
        }
      }
      ++eid;
    }
  }
  return 0;
}

int generate_PPT(
    const table_relation& tr,
    std::list<cut_vert>& cut_verts)    
{
  size_t pn = tr.pn();
  size_t tn = tr.tn();
  auto& P = tr.plane();
  auto& V = tr.mesh().verts();
  auto& M = tr.mesh().triangles();

  std::vector<std::list<cut_vert> > plane_cv(pn);
  
  /// PPT
#pragma omp parallel for
  for (size_t pi = 0; pi < pn; ++pi) {
    auto& pcv = plane_cv[pi];
    
    for (size_t pj = pi+1; pj < pn; ++pj) {
      if (tr.pp(pi, pj)==2) {
        std::list<size_t> com_ts;
        {
          // precondition, LA and LB is ordered from small to large.
          const auto& LA = tr.tuple_pt()[pi];
          const auto& LB = tr.tuple_pt()[pj];
          for (size_t id_pi=0,id_pj=0; id_pi<LA.size() && id_pj<LB.size();) {
            if (LA[id_pi] < LB[id_pj]) {
              ++id_pi;
            } else if (LA[id_pi] > LB[id_pj]) {
              ++id_pj;
            } else {
              com_ts.push_back(LA[id_pi]);
              ++id_pi; ++id_pj;
            }
          }
          
        }
        //for (size_t ti = 0; ti < tn; ++ti) { // simple method
        for (auto ti : com_ts) {
          auto si = tr.pt(pi, ti);
          auto sj = tr.pt(pj, ti);
          const size_t* const vp = &M[3*ti];
          
          bool is_add = false;
          if (si==2 && sj==2) {
            is_add = cal_ppt_intersection(
                &P[4*pi], &P[4*pj],
                &V[3*vp[0]], &V[3*vp[1]], &V[3*vp[2]]);
          } else if (si==1 &&
              ( (tr.pv(pj, vp[0])==0&&tr.pv(pi,vp[0])==0)
                || (tr.pv(pj, vp[1])==0&&tr.pv(pi,vp[1])==0)
                || (tr.pv(pj, vp[2])==0&&tr.pv(pi,vp[2])==0) ) )
          { is_add = true; }
          else if (sj==1 &&
              ( (tr.pv(pi, vp[0])==0&&tr.pv(pj,vp[0])==0)
                || (tr.pv(pi, vp[1])==0&&tr.pv(pj,vp[1])==0)
                || (tr.pv(pi, vp[2])==0&&tr.pv(pj,vp[2])==0) ) )
          { is_add = true; }
          
          if (is_add) {
            pcv.push_back(cut_vert(pi+tn, pj+tn, ti));
            pcv.back().set_on_triangle(ti);
          }
        }// for ti
      }// if
    }// for pj
  }// for pi

  for (auto& cv : plane_cv) {
    cut_verts.insert(cut_verts.end(), cv.begin(), cv.end());
  }
  
  return 0;
}

int generate_cut_verts(
    const table_relation& tr,
    const std::map<std::pair<size_t,size_t>, std::vector<size_t> >& e2f,
    std::vector<cut_vert>& cut_verts)
{
  std::vector<std::list<cut_vert> > triples(3);
  std::vector<std::string> triple_name={"PPP","TTP","PPT"};
  
  CALL_FUNC(generate_PPP(tr, triples[0]));  
  CALL_FUNC(generate_TTP(tr, e2f, triples[1]));
  CALL_FUNC(generate_PPT(tr, triples[2]));
      
  cut_verts.reserve(triples[0].size()+triples[1].size()+triples[2].size());
  for (size_t i = 0; i < 3; ++i) {
    cut_verts.insert(cut_verts.end(), triples[i].begin(), triples[i].end());
    std::cout << "-- " << triple_name[i] << " num : " << triples[i].size() << std::endl;
  }
  std::cout << "-- all triples of (PPP,TTP,PPT) : " << cut_verts.size() << std::endl;
  
  return 0;
}

typedef std::pair<size_t, size_t> vert_pair;

int sort_adjacent_triangles(
    const std::vector<double>& V,
    const std::vector<size_t>& M,
    std::map<vert_pair, std::vector<size_t> >& e2f)
{
  for (auto& ts : e2f) {
    if (ts.second.size() == 0) { return 1; }

    auto& faces = ts.second;
    size_t va = ts.first.first;
    size_t vb = ts.first.second;

    FASSERT(va < vb);

    {
      size_t k = INVALID_NUM;
      for (size_t i = 0; i < faces.size(); ++i) {
        size_t fid = faces[i];
        for (size_t j = 0; j < 3; ++j) {
          if (M[3*fid+j]==va && M[3*fid+(j+1)%3]==vb) {//norm outer forward
            k = i;
            break;
          }
        }
        if (k != INVALID_NUM) {
          break;
        }
      }
      FASSERT(k != INVALID_NUM);
      if (k != 0) {
        std::swap(faces[0], faces[k]);
      }
    }

    if (ts.second.size() == 2) { continue; }
    
    std::vector<size_t> third_vs;
    third_vs.reserve(faces.size());

    for (auto& t : faces) {
      for (size_t k = 0; k < 3; ++k) {
        if (M[3*t+k]!=va && M[3*t+k]!=vb) {
          third_vs.push_back(M[3*t+k]);
          break;
        }
      }
    }

    FASSERT(third_vs.size() == faces.size());

    std::vector<short> si(faces.size(), 0);
    for (size_t i = 1; i < si.size(); ++i) {
      si[i] = vol_sgn(&V[3*va], &V[3*third_vs[0]], &V[3*third_vs[i]], &V[3*vb]);
    }

    // Here, there is an assumption that there is no self-intersection.
    std::vector<size_t> new_third_vs, new_fs;
    new_third_vs.reserve(faces.size());
    new_fs.reserve(faces.size());

    // push the first
    new_fs.push_back(faces[0]);
    new_third_vs.push_back(third_vs[0]);
    
    for (size_t i = 1; i < si.size(); ++i) {
      if (si[i] < 0) {
        new_fs.push_back(faces[i]);
        new_third_vs.push_back(third_vs[i]);
      }
    }
    
    for (size_t i = 1; i < new_third_vs.size(); ++i) {
      for (size_t j = i+1; j < new_third_vs.size(); ++j) {
        int r = vol_sgn(&V[3*va], &V[3*new_third_vs[i]], &V[3*new_third_vs[j]], &V[3*vb]);
        if (r > 0) {
          std::swap(new_third_vs[i], new_third_vs[j]);
          std::swap(new_fs[i], new_fs[j]);
        }
        FASSERT(r != 0);

      }
    }
    
    for (size_t i = 1; i < si.size(); ++i) {
      if (si[i] == 0) {
        new_fs.push_back(faces[i]);
        new_third_vs.push_back(third_vs[i]);
      }
    }
    
    size_t last_s = new_third_vs.size();
    for (size_t i = 1; i < si.size(); ++i) {
      if (si[i] > 0) {
        new_fs.push_back(faces[i]);
        new_third_vs.push_back(third_vs[i]);
      }
    }
    
    FASSERT(third_vs.size() == new_third_vs.size());
    
    for (size_t i = last_s; i < new_third_vs.size(); ++i) {
      for (size_t j = i+1; j < new_third_vs.size(); ++j) {
        int r = vol_sgn(&V[3*va], &V[3*new_third_vs[i]], &V[3*new_third_vs[j]], &V[3*vb]);
        if (r > 0) {
          std::swap(new_third_vs[i], new_third_vs[j]);
          std::swap(new_fs[i], new_fs[j]);
        }
        FASSERT(r != 0);
      }
    }
    
    std::swap(ts.second, new_fs);
    
  } // for e2f
  
  return 0;
}

int build_edge_to_triangle_map(
    const std::vector<size_t>& M,
    std::map<vert_pair, std::vector<size_t> >& e2f)
{
  for (size_t i = 0; i+2 < M.size(); i+=3) {
    for (size_t j = 0; j < 3; ++j) {
      vert_pair e(M[i+j], M[i+(j+1)%3]);

      if (e.first > e.second) {
        std::swap(e.first, e.second);
      }

      auto it = e2f.find(e);
      if (it != e2f.end()) {
        it->second.push_back(i/3);
      } else {
        std::vector<size_t> fs;
        fs.reserve(2);
        fs.push_back(i/3);
        e2f.insert(std::make_pair(e, fs));
      }
    }
  }

  return 0;
}

int build_edge_to_triangle_map(
    const std::vector<double>& V,
    const std::vector<size_t>& M,
    std::map<vert_pair, std::vector<size_t> >& e2f)
{
  CALL_FUNC(build_edge_to_triangle_map(M, e2f));
  
  // sort triangles adjacent to each edge with more than 2 adjacent faces.
  CALL_FUNC(sort_adjacent_triangles(V, M, e2f));
  
  return 0;
}

int check_delete_repeated_plane(
    std::vector<double>& P, const double norm_eps, const double len_eps)
{
  std::vector<size_t> invalid_P;
  
  if (check_plane(P, norm_eps, len_eps, invalid_P)) {
    std::vector<double> new_P;
    new_P.reserve(P.size());

    size_t pn = P.size()/4;
    size_t cnt = 0;
    for (size_t pi=0; pi<pn; ++pi) {
      if (cnt < invalid_P.size() && pi==invalid_P[cnt]) {
        ++cnt;
      } else {
        for (size_t j = 0; j < 4; ++j) {
          new_P.push_back(P[4*pi+j]);
        }
      }
    }
    
    std::swap(new_P, P);
  }

  return 0;
}


int general_exact_cut(
    const pm_map& para,
    const std::string& output_file,
    const std::vector<double>& V,
    const std::vector<size_t>& M,
    std::vector<double>& P,
    std::vector<double>& CV,
    std::vector<size_t>& CC, size_t& cn, std::vector<double>& new_CV, double& all_vol)
{

  std::cout << "-- default prec: " << mpf_get_default_prec() << std::endl;
  size_t mpf_prec = 4096; // enough for 32 float-point input
  mpf_set_default_prec(mpf_prec);
  std::cout << "-- set     prec: " << mpf_prec << std::endl;

  auto tmp_bm_start = std::chrono::high_resolution_clock::now(); 
  
  base_grid_mesh bound_mesh(V, M);
  
  std::map<vert_pair, std::vector<size_t> > e2f;
  CALL_FUNC(build_edge_to_triangle_map(V, M, e2f));

  const double ave_len = ave_edge_length(V, e2f);

  std::cout << "-- edges num : " << e2f.size() << std::endl;
  std::cout << "-- average length: " << ave_len << std::endl;

  double eps = para["eps"].as<double>();
  std::cout << "-- eps for detecting degeneration: " << eps
            << std::endl;
  eps *= ave_len;
  std::cout << "-- eps (considering average length): " << eps
            << std::endl;
  {
    CALL_FUNC(check_topology(M, e2f));
    CALL_FUNC(check_volume(V, M, eps, all_vol));
    CALL_FUNC(check_degenerated_triangle(V, M, eps));
  }

  {
    double plane_eps = para["planeeps"].as<double>();
    std::cout << "-- eps for detecting plane degeneration (for norm): " << plane_eps
              << std::endl;
    double len_eps = plane_eps*ave_len;
    std::cout << "-- eps (considering average length): "
              << len_eps << std::endl;
    
    const size_t is_detec_plane = para["detecplane"].as<size_t>();
    if (is_detec_plane == 1) {
      CALL_FUNC(check_delete_repeated_plane(P, plane_eps, len_eps));
    }

    if (P.size() == 0) {
      std::cout << "-- [ ERROR ] there is no input plane after plane checking." << std::endl;
      return 1;
    }
  }


  auto tmp_bm_end = std::chrono::high_resolution_clock::now();
  
  {
    double time_taken=std::chrono::duration_cast<std::chrono::nanoseconds>(tmp_bm_end - tmp_bm_start).count();
    time_taken *= 1e-9;
    std::cout << "--- Time taken by program build-and-check-mesh is "
              << std::fixed << time_taken << std::setprecision(9)
              << "  sec" << std::endl; 
  }

  table_relation tr(bound_mesh, P);
  
  std::vector<cut_vert> cut_verts;
  
  CALL_FUNC(generate_cut_verts(tr, e2f, cut_verts));
  
  std::vector<e2f_map_type::iterator> eid_map;
  eid_map.reserve(e2f.size());
  for (e2f_map_type::iterator it = e2f.begin(); it != e2f.end(); ++it) {
    eid_map.push_back(it);
  }

  // set seq eps
  const double seq_eps_ratio = para["line_sort_eps_ratio"].as<double>();
  std::cout << "-- seq eps ratio: " << seq_eps_ratio << std::endl;
  const double seq_eps = ave_len*seq_eps_ratio;
  std::cout << "-- seq eps: " << seq_eps << std::endl;

  line_sort_type lst(line_sort_type::LINE_COORD);
  sequence_type seq(tr, e2f, eid_map, seq_eps, cut_verts, CV, lst);
  
  CALL_FUNC(seq.run());

  seq.get_cut_verts_type();
  seq.get_cut_verts_id(); // calculate CV id, will be used in cell building

  std::vector<std::vector<size_t> > face_triangles;
  CALL_FUNC(generate_cut_cells(para, seq, output_file+"-cf", CC, cn,
                               new_CV, face_triangles));

  return 0;
}

} // namespace grid
} // fxz
