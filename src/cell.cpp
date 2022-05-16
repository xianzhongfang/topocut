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
#include <queue>
#include <list>
#include <unordered_set>

#include "grid.h"
#include "seq.h"
#include "my_gmp.h"

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/AABB_segment_primitive.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_2.h>

const double EPS_A= pow(2, -30);
const double EPS_B= pow(2, -45);
const double EPS_C= pow(2, -50);
const double EPS_MERGE = pow(2, -52);
const double EPS_MV= pow(2, -53);
const double EPS_LAMBDA = 1e-3;


namespace fxz {
namespace grid {

typedef Eigen::Matrix<mpf_class,3,1> MPF_VEC;
typedef std::tuple<mpf_class, mpf_class, mpf_class> mpf_rational_point2d;



double my_abs(const Eigen::Matrix<double,3,1>& v)
{
  return (fabs(v[0])+fabs(v[1])+fabs(v[2]));
}

double my_abs2(const Eigen::Matrix<double,2,1>& v)
{
  return (fabs(v[0])+fabs(v[1]));
}


int vec2d_cross_sgn(
    const mpf_class& mA_0,
    const mpf_class& mA_1,
    const mpf_class& mA_d,
    const mpf_class& mB_0,
    const mpf_class& mB_1,
    const mpf_class& mB_d,
    const mpf_class& mC_0,
    const mpf_class& mC_1,
    const mpf_class& mC_d)
{
  mpf_class D = (mC_0*mB_d-mC_d*mB_0)*(mA_1*mB_d-mA_d*mB_1)
      -(mC_1*mB_d-mC_d*mB_1)*(mA_0*mB_d-mA_d*mB_0);
  return sgn(D)*sgn(mA_d)*sgn(mC_d);
}


void vector_minus(
    const std::vector<size_t>& A,
    const std::vector<size_t>& B,
    std::vector<size_t>& C)
{
  std::unordered_set<size_t> vs(B.begin(), B.end());

  C.clear();
  C.reserve(A.size());
  for (auto& x : A) {
    if (vs.find(x) == vs.end()) {
      C.push_back(x);
    }
  }
}


int triangulate_with_CGAL(
    int direc,
    const std::vector<double>& V,
    const std::vector<size_t>& poly,
    std::vector<size_t>& tris);


///////////////////////////////////////////////////////////////
class segment_space {
 public:
  segment_space(size_t n) : n_(n)
  {
    seg_.resize(n);
  }

  void set_min_max(double min_val, double max_val) {
    min_val_=min_val; max_val_=max_val;
    if (min_val_ == max_val_) {
      max_val_ += 1e-10;
    }
    FASSERT(max_val_>min_val_);

    dv_ = (max_val_ - min_val_)/n_;
    seg_val_.resize(n_+1);
    for (size_t i = 0; i < n_; ++i) {
      seg_val_[i] = dv_*(double)i+min_val_;
    }
    seg_val_[n_] = max_val_;
  }

  void add_segment(double L, double R, size_t va, size_t vb) {
    size_t L_id = cal_id(L);
    size_t R_id = cal_id(R);
    FASSERT(L_id!=INVALID_NUM && R_id!=INVALID_NUM && L_id <= R_id);
    for (size_t i = L_id; i <= R_id; ++i) {
      seg_[i].push_back(std::make_pair(va, vb));      
    }
  }

  const std::list<std::pair<size_t,size_t> >& get_segments(size_t id) {
    return seg_[id];
  }

  size_t cal_id(double val)
  {
    if (val < min_val_) val = min_val_;
    size_t r = std::floor((val-min_val_)/dv_);

    if (r>0 && r<=n_ && val>=seg_val_[r-1] && val<seg_val_[r]) {
      return r-1;
    } else if (r+1<=n_ && val>=seg_val_[r] && val < seg_val_[r+1]) {
      return r;
    } else if (r+2<=n_ && val>=seg_val_[r+1] && val < seg_val_[r+2]) {
      return r+1;
    } else if (r==n_ && val>=seg_val_[n_]) {
      return n_-1;
    }
    return INVALID_NUM;
  }

  void get_segs(
      double L, double R, std::set<std::pair<size_t,size_t> >& es) {
    if (L > R) {
      std::swap(L, R);
    }
    L -= EPS_C*fabs(L);
    R += EPS_C*fabs(R);
    size_t L_id = cal_id(L);
    size_t R_id = cal_id(R);
    FASSERT(L_id<=R_id);
    if (L_id == INVALID_NUM || R_id == INVALID_NUM) {return;}

    for (size_t i=L_id; i<=R_id; ++i) {
      es.insert(seg_[i].begin(), seg_[i].end());
    }
  }

 private:
  double min_val_, max_val_, dv_;
  size_t n_;
  std::vector<double> seg_val_;
  std::vector<std::list<std::pair<size_t,size_t> > > seg_;
};


typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_2 Point;
typedef K::Segment_2 Segment;
typedef CGAL::Search_traits_2<K> TreeTraits;
typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
typedef Neighbor_search::Tree Point_Tree;

int build_segment_space(
    const std::vector<std::vector<size_t> >& polys, const std::vector<size_t>& ps,
    const std::unordered_map<size_t, std::pair<double, double> >& p2uv_d,    
    segment_space& seg_x, segment_space& seg_y)
{
  std::vector<double> X, Y;
  std::vector<size_t> start_id;
  size_t cnt(0);
  for (auto& pid : ps) {
    const auto& po = polys[pid];
    start_id.push_back(cnt);
    for (size_t i = 0; i < po.size(); ++i) {
      auto it = p2uv_d.find(po[i]);
      const auto& uv = it->second;

      X.push_back(uv.first);
      Y.push_back(uv.second);
    }
    cnt += po.size();
  }

  double min_X = *std::min_element(X.begin(), X.end());
  double max_X = *std::max_element(X.begin(), X.end());
  double min_Y = *std::min_element(Y.begin(), Y.end());
  double max_Y = *std::max_element(Y.begin(), Y.end());
      
  seg_x.set_min_max(min_X, max_X);
  seg_y.set_min_max(min_Y, max_Y);

  cnt = 0;
  for (auto& pid : ps) {
    auto& po = polys[pid];
    for (size_t i = 0; i < po.size(); ++i) {
      size_t j = (i+1)%po.size();
      double a(X[cnt+i]),  b(X[cnt+j]);
      if (a > b) std::swap(a, b);
      seg_x.add_segment(a, b, po[i], po[j]);
          
      a = Y[cnt+i]; b = Y[cnt+j];
      if (a > b) std::swap(a, b);          
      seg_y.add_segment(a, b, po[i], po[j]);
    }
        
    cnt += po.size();
  }

  return 0;
}

template<int T, typename VEC>
int build_segment_space(
    const std::vector<std::vector<size_t> >& polys,
    const VEC& ps,

    const std::unordered_map<size_t, std::pair<double, double> >& p2uv_d,        
    segment_space& seg_y)
{
  std::vector<double> Y;
  FASSERT(ps.size() > 0);

  for (auto pid : ps) {
    if (polys[pid].size() == 1) {continue;}
    for (auto v : polys[pid]) {
      auto it = p2uv_d.find(v);
      const auto& uv = it->second;

      if (T == 1) {
        Y.push_back(uv.second);
      } else {
        Y.push_back(uv.first);
      }
    }
  }

  if (Y.size() == 0) {
    seg_y.set_min_max(-1e100,1e100);
    return 0;
  }

  double min_Y = *std::min_element(Y.begin(), Y.end());
  double max_Y = *std::max_element(Y.begin(), Y.end());


  seg_y.set_min_max(min_Y, max_Y);


  size_t cnt = 0;
  for (auto pid : ps) {
    auto& po = polys[pid];
    if (po.size()==1) {continue;}
    for (size_t i = 0; i < po.size(); ++i) {
      size_t j = (i+1)%po.size();
      FASSERT(cnt+i<Y.size() && cnt+j<Y.size());
      double a = Y[cnt+i], b = Y[cnt+j];
      if (a > b) std::swap(a, b);          
      seg_y.add_segment(a, b, po[i], po[j]);
    }
    cnt += po.size();
  }

  return 0;
}

int build_trees(
    const std::vector<std::vector<size_t> >& polys, const std::vector<size_t>& ps,
    const std::unordered_map<size_t, std::pair<double,double> >& p2uv_d,    
    std::vector<std::shared_ptr<Point_Tree> >& po_trees,
    std::vector<std::map<std::pair<double,double>,size_t > >& p2id)
{
  std::list<Segment> segs;
  p2id.resize(ps.size());
  po_trees.reserve(ps.size());
  size_t cnt=0;
  for (auto id : ps) {
    std::list<Point> points;

    for (auto v : polys[id]) {
      auto it = p2uv_d.find(v);
      Point pc(it->second.first, it->second.second);
      p2id[cnt].insert(std::make_pair(it->second, v));
      points.push_back(pc);

    }

    po_trees.push_back(std::shared_ptr<Point_Tree>(
        new Point_Tree(points.begin(), points.end())));
    ++cnt;
  }

  return 0;
}
///////////////////////////////////////////////////////////////


class half_edge {
 public:
  half_edge(size_t v, size_t f, size_t inter_f)
      : start_v_(v), ori_f_(f), inter_f_(inter_f),
        next_(invalid_ptr()), prev_(invalid_ptr()),
        oppo_(invalid_ptr()), foppo_(invalid_ptr()), ce_id_(INVALID_NUM) { }

  typedef half_edge* HE_PTR;
  typedef const half_edge* const HE_CPTR;
  
  HE_CPTR next() const { return next_; }
  HE_CPTR prev() const { return prev_; }
  HE_CPTR oppo() const { return oppo_; }
  HE_CPTR foppo()const { return foppo_;}

  HE_PTR& next() { return next_; }
  HE_PTR& prev() { return prev_; }
  HE_PTR& oppo() { return oppo_; }
  HE_PTR& foppo(){ return foppo_;}

  size_t ce_id() const { return ce_id_; }
  size_t& ce_id() { return ce_id_; }

  static std::nullptr_t invalid_ptr() { return nullptr; }

  static bool is_valid(half_edge::HE_PTR ptr) {
    return (ptr!=invalid_ptr());
  }

  static bool is_valid_vert(size_t id) {
    return (id != INVALID_NUM);
  }

  /////////
  size_t start_v() const { return start_v_; }
  size_t end_v() const {
    if (oppo() != invalid_ptr()) {
      return oppo()->start_v();
    } else if (next() != invalid_ptr()) {
      return next()->start_v();
    }
    return INVALID_NUM;
  }
  size_t ori_f() const { return ori_f_; }
  size_t& start_v(){ return start_v_; }
  size_t ori_f() { return ori_f_; }

  size_t inter_f() const { return inter_f_; }
  size_t& inter_f() { return inter_f_; }

  char flag() const { return flag_; }
  char& flag() { return flag_; }

 private:
  size_t start_v_, ori_f_, inter_f_;
  HE_PTR next_, prev_, oppo_, foppo_;

  // temp flag
  char flag_;
  size_t ce_id_;
};

class hf_mesh {
 public:
  hf_mesh() {}

  // used for parallel running
  void resize(size_t num) { he_.clear(); he_.resize(num); }

  std::list<half_edge>& ori_face_he(size_t id) {
    FASSERT(id < he_.size());
    return he_[id];
  }
  
  const std::list<half_edge>& ori_face_he(size_t id) const {
    FASSERT(id < he_.size());
    return he_[id];
  }

  half_edge::HE_PTR add_he(size_t va, size_t cid, size_t inter_f)
  {
    assert(cid < he_.size());
    he_[cid].push_back(half_edge(va, cid, inter_f));
    return &(he_[cid].back());
  }

  void set_all_flag(char s)
  {
    for (auto& fh : he_) {
      for (auto& he : fh) {
        he.flag() = s;
      }
    }
  }

  int set_flag(size_t ori_fid, char s)
  {
    FASSERT(ori_fid < he_.size());    
    if (ori_fid >= he_.size()) return 1;

    for (auto& he : he_[ori_fid]) {
      he.flag() = s;
    }
    return 0;
  }

  half_edge& he(half_edge::HE_PTR e)
  {
    return *e;
  }

  const half_edge& he(half_edge::HE_CPTR e) const
  {
    return *e;
  }

  void set_next_prev(half_edge::HE_PTR ea, half_edge::HE_PTR eb)
  {
    FASSERT(ea != half_edge::invalid_ptr() && eb != half_edge::invalid_ptr());
    he(ea).next() = eb;
    he(eb).prev() = ea;
  }

  void set_oppo(half_edge::HE_PTR ea, half_edge::HE_PTR eb)
  {
    FASSERT(ea != half_edge::invalid_ptr() && eb != half_edge::invalid_ptr());
    he(ea).oppo() = eb;
    he(eb).oppo() = ea;
  }

  void set_foppo(half_edge::HE_PTR ea, half_edge::HE_PTR eb)
  {
    FASSERT(ea != half_edge::invalid_ptr() && eb != half_edge::invalid_ptr());
    he(ea).foppo() = eb;
    he(eb).foppo() = ea;
  }

 private:
  std::vector<std::list<half_edge> > he_;
};

class cut_cell {
 public:
  cut_cell(const sequence_type& seq) : seq_(seq) {
    sorted_cv_.resize(seq_.tr().pn());
  }

  int build_faces()
  {
    std::cout << "# start to generate cut-faces ..." << std::endl;
    hfm_.resize(2*seq_.tr().pn()+seq_.tr().tn());
    CALL_FUNC(build_bound_faces());
    CALL_FUNC(build_plane_faces());
    std::cout << "# end of generating cut-faces." << std::endl;

    {
      std::unordered_set<size_t> all_sorted_cv;
      for (auto& C : sorted_cv_) {
        all_sorted_cv.insert(C.begin(), C.end());
      }
      std::cout << "-- cut-verts by using cycle-sort : " << all_sorted_cv.size() << std::endl;
      std::cout << "-- all-cut-verts-number          : " << seq_.pure_cv_num() << std::endl;
    }
    
    return 0;
  }

  int build_cells()
  {
    // build cut edges map
    std::map<std::pair<size_t,size_t>,size_t> ce_map;
    // get map and set he->eid;
    size_t ce_num = get_cutedge_map(ce_map);

    // eid -> he (only for boundary, others can be visited by foppo operator)
    // if the boundary input mesh is not manifold, need more info for each ce
    // if the boundary is manifold, only one he is enough
    std::vector<std::list<half_edge::HE_PTR> > ce2he(ce_num);
    // some of them may be invalid ptr at first, which means that
    // they are not in the current half-face structure.
    for (size_t ci = 0; ci < seq_.tr().tn(); ++ci) { // bound he
      auto& hes = hfm_.ori_face_he(ci);
      for (auto& he : hes) {
        add_he_into_ce2he(ce2he, &he);
      }
    }

    // add cutting plane faces one by one.
    for (size_t pi = 0; pi < seq_.tr().pn(); ++pi) {
      if (add_plane_into_hf_mesh(ce2he, pi)) { return 1; }
    }

    std::cout << "-- cycle-sorted cut-edge number : " << sorted_ce_.size() << std::endl;
    
    return 0;
  }

  int add_he_into_ce2he(
      std::vector<std::list<half_edge::HE_PTR> >& ce2he, half_edge::HE_PTR add_he)
  {
    size_t ceid = add_he->ce_id();

    FASSERT(ceid != INVALID_NUM);
    
    if (ce2he[ceid].empty()) {
      ce2he[ceid].push_back(add_he);
      return 1;
    } else {
      bool is_got = false;
      for (auto& hp : ce2he[ceid]) {//check whether the half-edge can be got from other he.
        auto st = hp;
        do {
          st = st->oppo();
          if (st == add_he) {break;}
          if (!half_edge::is_valid(st)) {break;}
          st = st->foppo();
        } while (half_edge::is_valid(st) && st!=add_he);
        if (st == add_he) {is_got=true; break;}
        
        st = hp->foppo();
        while (half_edge::is_valid(st) && st != add_he) {
          st = st->oppo();
          if (st == add_he) {break;}
          if (!half_edge::is_valid(st)) {break;}
          st = st->foppo();
        }
        if (st == add_he) {is_got=true; break;}
      } // for
      
      if (!is_got) {
        ce2he[ceid].push_back(add_he);
      }
      if (is_got) {return 2;}
    }

    return 0;
  }

  int add_plane_into_hf_mesh(
      std::vector<std::list<half_edge::HE_PTR> >& ce2he,
      size_t pid)
  {
    size_t tn = seq_.tr().tn();
    auto& hes = hfm_.ori_face_he(pid+tn);

    // used for searching boundary bd
    hfm_.set_flag(pid+tn, 0);
    
    for (auto& he : hes) {
      size_t eid = he.ce_id();
      if (!ce2he[eid].empty()) {
        he.flag() = 1;
        if (half_edge::is_valid(he.oppo())) {
          he.oppo()->flag() = 1;
        }
      }
      { // for test
        if (!half_edge::is_valid(he.oppo())) {
          FASSERT(he.flag() == 1);
        }
      }
    }
    
    // get all sub domain of cutted faces that are needed to be inserted.
    std::vector<std::vector<half_edge::HE_PTR> > faces_bd;
    for (auto& he : hes) {
      if (he.flag() == 1) {
        faces_bd.push_back(std::vector<half_edge::HE_PTR>());
        search_bound_he(&he, faces_bd.back());
        FASSERT(faces_bd.back().size() > 0);
      }
    }
    

    std::unordered_set<half_edge::HE_PTR> no_pasted_he;
    
    // insert faces into cell
    for (auto& bd : faces_bd) {
      int r = add_subplane_into_hf_mesh(ce2he, bd);
      if (r == 2) { // no pasted
        for (auto& h : bd) {
          no_pasted_he.insert(h);
        }
      }
    }
    
    // add new ce into ce2he
    for (auto& he : hes) {
      if (no_pasted_he.empty() || (no_pasted_he.find(&he)==no_pasted_he.end())) {
        int r = add_he_into_ce2he(ce2he, &he);
        FASSERT(r==1 || r==2);
      }
      
      { // for test: the oppo half-edge
        FASSERT(half_edge::is_valid(he.oppo()));
      }
    }

    // for test: the oppo half-edge
    for (auto& he : hfm_.ori_face_he(tn+seq_.tr().pn()+pid)) {
      FASSERT(half_edge::is_valid(he.oppo()));
    }
    
    return 0;
  }

  // for search insersion places
  int global_search_loop(half_edge::HE_PTR start_he, std::list<half_edge::HE_PTR>& got_bd)
  {
    std::queue<half_edge::HE_PTR> que;
    std::unordered_set<half_edge::HE_PTR> he_set;
    que.push(start_he);

    while (!que.empty()) {
      auto s = que.front();
      que.pop();

      if (he_set.count(s) != 0) {continue;}
      
      auto h = s;
      do {
        if (h->flag() == 'c') {
          got_bd.push_back(h);
        }
        he_set.insert(h);
        auto oh = h->oppo();
        if (he_set.count(oh) == 0) {
          que.push(oh);
        }
        h = h->next();
      } while (h != s);
    } // while


    return 0;
  }
  
  int add_subplane_into_hf_mesh(
      const std::vector<std::list<half_edge::HE_PTR> >& ce2he,
      const std::vector<half_edge::HE_PTR>& face_bd)
  {
    std::vector<size_t> loop_vs;
    loop_vs.reserve(face_bd.size());
    
    for (size_t i = 0; i < face_bd.size(); ++i) {
      loop_vs.push_back(face_bd[i]->start_v());
    }
    
    size_t pid = face_bd.front()->ori_f()-seq_.tr().tn();
    
    FASSERT(pid < seq_.tr().pn());
    
    std::vector<char> bd_st(face_bd.size(), 's'); // s: split; u: up; d: down
    
    std::vector<std::list<half_edge::HE_PTR> > candidate_paste_bd(face_bd.size());
    std::vector<half_edge::HE_PTR> paste_bd(face_bd.size(),half_edge::invalid_ptr());
    
    
    // for each bound half-edge, find the oppo halfedge in half-face structure
    size_t start_id = 0;
    size_t min_adj_cnt(2);
    
    for (size_t hi = 0; hi < face_bd.size(); ++hi) {
      auto& he = face_bd[hi];
      FASSERT(half_edge::is_valid(he));
      size_t vb = he->end_v();
      assert(he->ce_id() < ce2he.size());
      auto& edge_hes = ce2he[he->ce_id()];
      FASSERT(!edge_hes.empty());
      
      auto& adj_hes = candidate_paste_bd[hi];
      for (auto sub_domain_he : edge_hes) {
        // search all the adjacent he about he by oppo and foppo
        search_adj_he_by_oppo_and_foppo(sub_domain_he, vb, adj_hes);
      }
      
      judge_adj_hes_by_cycle_sort(he, adj_hes);
      
      if (adj_hes.size() == 1) {
        bd_st[hi] = plane_status(adj_hes.front(), pid);
        paste_bd[hi] = adj_hes.front();
      } else if (adj_hes.size() == 0) {
        std::cout << "-- [ MESSAGE ] There is no appropriate half-edge to be pasted for subplane."
                  << std::endl;
        return 2;
      }
      if (adj_hes.size() < min_adj_cnt) {
        start_id = hi;
        min_adj_cnt = adj_hes.size();
      }

    }
    
    std::unordered_map<half_edge::HE_PTR,size_t> candidate_map;
    for (size_t i = 0; i < candidate_paste_bd.size(); ++i) {
      auto& bd = candidate_paste_bd[i];
      for (auto& he : bd) {
        FASSERT(he->start_v() == face_bd[i]->end_v());
        FASSERT(he->end_v() == face_bd[i]->start_v());
        candidate_map.insert(std::make_pair(he, i));
        he->flag() = 'c'; // candidate
      }
    }
    
    size_t valid_loop(0);

    
    for (auto& start_he : candidate_paste_bd[start_id]) {
      auto st = start_he;

      std::list<half_edge::HE_PTR> got_bd;
      
      size_t v_cnt = start_id;
      
      FASSERT(loop_vs[v_cnt] == st->end_v());
      
      do {
        got_bd.push_back(st);
        // find the approrpiate he
        v_cnt = (v_cnt+loop_vs.size()-1)%loop_vs.size();
        size_t cur_vid = loop_vs[v_cnt];

        
        std::list<half_edge::HE_PTR> nex_he;
        cycle_search_he(st, nex_he, cur_vid);
        
        if (nex_he.size() == 1) {
          st = nex_he.front();
        } else if (nex_he.size() > 1) {
          std::cerr << "... [ ERROR ] L" << __LINE__ 
                    << " edge: " << st->start_v() << ", " << st->end_v() << std::endl;
          for (auto& h : nex_he) {
            std::cerr << "   h edge: " << h->start_v() << ", " << h->end_v()
                      << ", h fid: " << h->ori_f() << ", h inter f: " << h->inter_f() << std::endl;
            st = h;
            return 1; // error
          }
        } else {
          //std::cout << "... there is no next he." <<  std::endl; // for test
          break; // there is no next he
        }
      } while (st != start_he);
      
      if (st != start_he || got_bd.size() != face_bd.size()) {
        continue;
      }
      
      FASSERT(v_cnt == start_id);
      
      // add into paste_bd
      for (auto& he : got_bd) {
        auto it = candidate_map.find(he);
        if (it != candidate_map.end()) {
          paste_bd[it->second] = he;
        } else {
          std::cerr << "# [ ERROR ] strange error, paste the map has no such half-edge."
                    << std::endl;
          return 1;
        }
      }
      
      { // for checking
        bool is_ok = true;
        for (size_t i = 0; i < paste_bd.size() && is_ok; ++i) {
          if (!half_edge::is_valid(paste_bd[i])) {
            is_ok = false;
          }
        }
        if (is_ok) {++valid_loop;}
      }
    } // for


    if (valid_loop == 0) {
      std::cout << "# [ MESSAGE ] can not find the related paste loop by simple loop searching, now use global searching."
                << std::endl;
      FASSERT(candidate_paste_bd[start_id].size() == 1);
      std::list<half_edge::HE_PTR> got_bd;
      for (auto h : candidate_paste_bd[start_id]) {
        global_search_loop(h, got_bd);
        if (got_bd.size() == paste_bd.size()) {
          break;
        } else {
          got_bd.clear();
        }
      }

      if (got_bd.size() == paste_bd.size()) {
        for (auto he : got_bd) {
          auto it = candidate_map.find(he);
          if (it != candidate_map.end()) {
            paste_bd[it->second] = he;
          } else {
            std::cerr << "# [ ERROR ] global search ... strange error, paste the map has no such half-edge."
                      << std::endl;
            return 1;
          }
        }//for
      } else {
        std::cerr << "# [ ERROR ] can not find the related paste loop with global searching." << std::endl;
        ERROR_RETURN;
      }//if
    }//if
    
    // reset the flag, in order to use the flag in the next time
    for (size_t i = 0; i < candidate_paste_bd.size(); ++i) {
      auto& bd = candidate_paste_bd[i];
      for (auto& he : bd) {
        he->flag() = 0;
      }
    }

    if (valid_loop > 1) {
        std::cerr << "# [ ERROR ] more than one valid loops for pasting subplane." << std::endl;
        return 1;
    }
    
    // paste half-edges
    for (size_t i = 0; i < face_bd.size(); ++i) {
      auto he = face_bd[i];
      auto op_he = face_bd[i]->foppo();
      auto op_paste_he = paste_bd[i]->oppo();

      if (bd_st[i] == 's' || bd_st[i] == 'u') {

        FASSERT(half_edge::is_valid(paste_bd[i]));
        FASSERT(half_edge::is_valid(op_he));
        FASSERT(he->start_v() == paste_bd[i]->end_v());
        FASSERT(he->end_v() == paste_bd[i]->start_v());
        FASSERT(op_he->start_v() == op_paste_he->end_v());
        
        hfm_.set_oppo(he, paste_bd[i]);
        hfm_.set_oppo(op_he, op_paste_he);
        
      } else {
        if (op_paste_he->ori_f() < seq_.tr().tn()) {
          hfm_.set_oppo(he, paste_bd[i]);
          hfm_.set_oppo(op_he, op_paste_he);
        } else {
          auto new_paste_he = op_paste_he->foppo();
          auto op_new_paste_he = new_paste_he->oppo();
          FASSERT(half_edge::is_valid(new_paste_he));
          FASSERT(half_edge::is_valid(new_paste_he->oppo()));
          FASSERT(new_paste_he->start_v() == he->end_v());
          hfm_.set_oppo(he, new_paste_he);
          hfm_.set_oppo(op_he, op_new_paste_he);
        }
      }
      
    } // for
    
    return 0;
  }
  
  int judge_triangle_and_plane(size_t pid, size_t f)
  {
    FASSERT(f < seq_.tr().tn());
    size_t va = seq_.tr().mesh().tv(f, 0);
    size_t vb = seq_.tr().mesh().tv(f, 1);
    size_t vc = seq_.tr().mesh().tv(f, 2);

    int a = seq_.tr().pv(pid, va);
    int b = seq_.tr().pv(pid, vb);
    int c = seq_.tr().pv(pid, vc);

    int s = a+b+c;

    if (a>=0 && b>=0 && c>=0 && s!=0) { return 1; }
    if (a<=0 && b<=0 && c<=0 && s!=0) { return -1;}
    return 0;
  }

  // here, only remove ambiguity half-edges in edges of original mesh
  void judge_adj_hes_by_cycle_sort(
      half_edge::HE_PTR insert_he,
      std::list<half_edge::HE_PTR>& adj_hes)
  {
    std::vector<std::list<half_edge::HE_PTR>::iterator> del_hes;

    // get the edge id of insert_he
    size_t va = insert_he->start_v();
    size_t vb = insert_he->end_v();
    FASSERT(half_edge::is_valid_vert(va) && half_edge::is_valid_vert(vb));

    char ty_a = seq_.cut_verts_type()[va];
    char ty_b = seq_.cut_verts_type()[vb];
    
    if ((ty_a=='e' || ty_a=='v') && (ty_b=='e' || ty_b=='v')) {
      size_t Lv(INVALID_NUM), Rv(INVALID_NUM);

      if (ty_a=='v' && ty_b=='v') {
        Lv = va; Rv = vb;
        
      } else if (ty_a=='e' && ty_b=='e') {
        size_t Leid = seq_.cv(seq_.cut_verts_id()[va]).ele_id(); //edge id
        size_t Reid = seq_.cv(seq_.cut_verts_id()[vb]).ele_id(); //edge id

        if (Leid == Reid) {
          auto ev = seq_.ev(Leid);
          const auto& Q_tt = seq_.seq_tt(Leid);
          for (auto& v : Q_tt) {
            size_t vid = seq_.cv(v.cvid).id();
            if (vid == va || vid == vb) {
              if (vid == va) {
                Lv = ev.first;
                Rv = ev.second;
              } else {
                Lv = ev.second;
                Rv = ev.first;
              }
              break;
            }//if
          }//for

          FASSERT(Lv != INVALID_NUM && Rv != INVALID_NUM);
        }//if
        
      } else if (ty_a=='e' && ty_b=='v') {
        size_t Leid = seq_.cv(seq_.cut_verts_id()[va]).ele_id(); //edge id        
        auto ev = seq_.ev(Leid);
        if (ev.first == vb) {
          Lv = ev.second;
          Rv = ev.first;
        } else if (ev.second == vb) {
          Lv = ev.first;
          Rv = ev.second;
        }
        
      } else if (ty_a=='v' && ty_b=='e') {
        size_t Reid = seq_.cv(seq_.cut_verts_id()[vb]).ele_id(); //edge id        
        auto ev = seq_.ev(Reid);
        if (ev.first == va) {
          Lv = ev.first;
          Rv = ev.second;
        } else if (ev.second == va) {
          Lv = ev.second;
          Rv = ev.first;
        }
      }
      
      if (Lv!=INVALID_NUM && Rv!=INVALID_NUM) {
        for (auto it = adj_hes.begin(); it != adj_hes.end(); ++it) {
          if (!is_appropriate_he(insert_he, Lv, Rv, *it)) {
            del_hes.push_back(it);
          }
        }

        for (auto& it : del_hes) {
          adj_hes.erase(it);
        }
      }
      
    } // if
    
    return;
  }

  template<typename T>
  void get_face_norm(
      const double* va, const double* vb, const double* vc,
      T& Nt)
  {
    T VA, VB, VC;
    VA << va[0], va[1], va[2];
    VB << vb[0], vb[1], vb[2];
    VC << vc[0], vc[1], vc[2];
    Nt = (VB-VA).cross(VC-VA);
  }

  void get_point_lambda_in_triangle(
      const double* p,
      const double* va, const double* vb, const double* vc,
      Eigen::Matrix<double, 2, 1>& L)
  {
    Eigen::Map<const Eigen::Matrix<double,3,1> > P(p, 3, 1);
    Eigen::Map<const Eigen::Matrix<double,3,1> > A(va, 3, 1);
    Eigen::Map<const Eigen::Matrix<double,3,1> > B(vb, 3, 1);
    Eigen::Map<const Eigen::Matrix<double,3,1> > C(vc, 3, 1);

    Eigen::Matrix<double, 3, 2> M;
    M.col(0) = A-C;
    M.col(1) = B-C;
    Eigen::Matrix<double, 2, 2> MTM = M.transpose()*M;
    L = MTM.inverse()*M.transpose()*(P-C);
  }

  void get_face_norm(size_t fid, MPF_VEC& N)
  {
    if (fid < seq_.tr().tn()) {
      auto& bm = seq_.tr().mesh();
      auto va_ptr = bm.tv(fid);
      get_face_norm(
          bm.vert_coord(va_ptr[0]),
          bm.vert_coord(va_ptr[1]),
          bm.vert_coord(va_ptr[2]), N);
      N *= -1.0;
    } else {
      size_t pid = fid - seq_.tr().tn();
      bool is_oppo = false;
      if (pid >= seq_.tr().pn()) {
        pid -= seq_.tr().pn();
        is_oppo = true;
      }
      auto& P = seq_.tr().plane();
      N << P[4*pid], P[4*pid+1], P[4*pid+2];
      if (is_oppo) { N *= -1.0; }
    }
  }

  void get_face_norm(
      const half_edge::HE_PTR he, MPF_VEC& N)
  {
    size_t fid = he->ori_f();
    get_face_norm(fid, N);
  }
  
  bool is_appropriate_he(
      half_edge::HE_PTR insert_he,
      size_t va, size_t vb,
      half_edge::HE_PTR adj_he)
  {
    FASSERT(insert_he->start_v() == adj_he->end_v());
    auto& bm = seq_.tr().mesh();
    auto& V = bm.verts();
    MPF_VEC VA, VB;
    VA << V[3*va], V[3*va+1], V[3*va+2];
    VB << V[3*vb], V[3*vb+1], V[3*vb+2];
    MPF_VEC E = VB-VA;

    sorted_ce_.insert(insert_he->ce_id());

    MPF_VEC Nf1, Nf2, Np;

    FASSERT(half_edge::is_valid(adj_he->oppo()));    
    get_face_norm(adj_he->oppo(), Nf1);

    get_face_norm(adj_he, Nf2);
    Nf2 *= -1.0;

    FASSERT(insert_he->ori_f() >= seq_.tr().tn());
    FASSERT(insert_he->ori_f() < seq_.tr().pn()+seq_.tr().tn());
    get_face_norm(insert_he, Np);

    
    int s12 = sgn(Nf1.cross(Nf2).dot(E));
    int s1 = sgn(Nf1.cross(Np).dot(E));
    int s2 = sgn(Np.cross(Nf2).dot(E));
    
    if ((s12>=0 && s1>0 && s2>0)
        || (s12<0 && (s1>=0 || s2>=0))) {
      return true;
    }
    return false;
  }

  char plane_status(
      half_edge::HE_PTR he, size_t pid)
  {
    auto op_he = he->oppo();

    size_t fa = he->ori_f();
    size_t fb = op_he->ori_f();

    if (fa == fb) {
      return 's';
    }

    if (fa < seq_.tr().tn() && fb < seq_.tr().tn()) {
      int sa = judge_triangle_and_plane(pid, fa);
      int sb = judge_triangle_and_plane(pid, fb);

      if (sa*sb <= 0) {
        return 's';
      } else {
        if (sa > 0 && sb > 0) { return 'u'; }
        else { return 'd'; }
      }
    }
    
    return 's';
  }

  void cycle_search_he(
      half_edge::HE_PTR start_he,
      std::list<half_edge::HE_PTR>& nex_he,
      size_t cur_vid)
  {
    auto st = start_he;
    
    do {
      st = st->next();
      FASSERT(half_edge::is_valid(st));
      if (st->flag() == 'c' && st->end_v()==cur_vid) {
        nex_he.push_back(st);
      }
      FASSERT(half_edge::is_valid(st->oppo()));
      if (!half_edge::is_valid(st->oppo())) {
        std::cerr << "--- [ ERROR ] st->ori_f(): " << st->ori_f()
                  << ", tn: " << seq_.tr().tn()
                  << ", pn: " << seq_.tr().pn() << std::endl;
      }
      st = st->oppo();
    } while (st != start_he);
  }
  
  void search_adj_he_by_oppo_and_foppo(
      half_edge::HE_PTR start_he, size_t vb,
      std::list<half_edge::HE_PTR>& adj_hes)
  {
    auto st = start_he;
    do {
      if (st->start_v() == vb) {
        adj_hes.push_back(st);
      }
      st = st->oppo();
      FASSERT(half_edge::is_valid(st));
      if (st->start_v() == vb) {
        adj_hes.push_back(st);
      }
      st = st->foppo();
    } while (half_edge::is_valid(st) && st!=start_he);
    
    if (st != start_he) {
      st = start_he->foppo();
      if (half_edge::is_valid(st)) {
        do {
          if (st->start_v() == vb) {
            adj_hes.push_back(st);
          }
          st = st->oppo();
          FASSERT(half_edge::is_valid(st));
          if (st->start_v() == vb) {
            adj_hes.push_back(st);
          }
          st = st->foppo();
        } while (half_edge::is_valid(st) && st!=start_he->foppo());
      }
    }

  }

  void search_bound_he(
      half_edge::HE_PTR start_he,
      std::vector<half_edge::HE_PTR>& bound)
  {
    auto st = start_he;

    do {
      bound.push_back(st);
      st->flag() = 2;
      st = st->next();
      while (st->flag() == 0) {
        st = st->oppo()->next();
      }
    } while (st != start_he);
  }

  size_t get_cutedge_map(std::map<std::pair<size_t,size_t>,size_t>& ce_map)
  {
    size_t pn = seq_.tr().pn();
    size_t tn = seq_.tr().tn();

    std::cout << "-- tn :" << tn  << std::endl;
    std::cout << "-- pn :" << pn  << std::endl;    
    
    size_t cnt = 0;
    for (size_t ci = 0; ci < 2*pn+tn; ++ci) {
      auto& hes = hfm_.ori_face_he(ci);
      for (auto& he : hes) {
        std::pair<size_t,size_t> e(he.start_v(), he.end_v());
        if (e.first > e.second) { std::swap(e.first, e.second); }

        auto it = ce_map.find(e);
        size_t id(INVALID_NUM);
        if (it != ce_map.end()) {
          id = it->second;
        } else {
          id = cnt++;
          ce_map.insert(std::make_pair(e, id));
        }
        FASSERT(id != INVALID_NUM);
        he.ce_id() = id;
      } // for hes
    } // for ci
    
    std::cout << "-- all-cut-edges-number : " << cnt << std::endl;

    return cnt;
  }

  int build_bound_faces()
  {
    size_t tn = seq_.tr().mesh().tn();
    std::vector<half_edge::HE_PTR> tri_he3(3*tn);
    // add boundary of each triangle
#pragma omp parallel for    
    for (size_t ti = 0; ti < tn; ++ti) {
      build_triangle_he(ti, &tri_he3[3*ti]);
    }

    // add oppo
    CALL_FUNC(build_bound_oppo(tri_he3));
    
    // add inner cut-edges into
    std::vector<std::list<size_t> > tri_pt(tn);
    auto& pt_map = seq_.pt_map();
    for (auto& pt : pt_map) {
      FASSERT(pt.first.first >= tn);
      size_t t = pt.first.second;
      FASSERT(t < tn);
      tri_pt[t].push_back(pt.second);
    }

     // clean the same inner half-edges
#pragma omp parallel for          
    for (size_t ti = 0; ti < tn; ++ti) {
      std::set<std::pair<size_t,size_t> > inner_lines_ending;
      std::vector<std::list<size_t>::iterator > del_pt;
      std::list<size_t> new_pt;
      for (auto& ptid : tri_pt[ti]) {
        auto bd = seq_.pt_bd(ptid);
        bd.first = seq_.cv(bd.first).id();
        bd.second = seq_.cv(bd.second).id();
        if (bd.first > bd.second) { std::swap(bd.first, bd.second); }
        if (inner_lines_ending.find(bd) == inner_lines_ending.end()) {
          inner_lines_ending.insert(bd);
          new_pt.push_back(ptid);
        }
      }
      std::swap(new_pt, tri_pt[ti]);
    }

    // can parallel running
#pragma omp parallel for
    for (size_t ti = 0; ti < tn; ++ti) {
      build_triangle_inner_he(ti, tri_pt[ti], tri_he3[3*ti]);
    }
    
    return 0;
  }

  int build_triangle_inner_he(
      size_t ti, const std::list<size_t>& pt, half_edge::HE_PTR start_he)
  {
    std::unordered_map<size_t, half_edge::HE_PTR> v2he;

    {
      auto st = start_he;
      FASSERT(start_he->ori_f() == ti);
      do {
        size_t v = st->start_v();
        v2he[v] = st;

        st = st->next();
      } while (st != start_he);

    }

    // here should be dealed with serial running
    for (auto& sid : pt) { // for all pt line belonging to triangle ti.
      auto& Q = seq_.seq_pt(sid);
      auto& s_bd = seq_.pt_bd(sid);
      
      size_t st_v = seq_.cv(s_bd.first).id();
      size_t en_v = seq_.cv(s_bd.second).id();
      
      if (st_v < seq_.tr().vn() and en_v < seq_.tr().vn()) {
        continue; // belongs to edge
      }
      
      FASSERT(st_v != en_v);
      
      std::vector<half_edge::HE_PTR> line_he;
      
      size_t pre_id = st_v;
      
      for (auto& v : Q) {
        size_t vid = seq_.cv(v.cvid).id();
        if (vid != pre_id) {
          line_he.push_back(hfm_.add_he(pre_id, ti, INVALID_NUM));
          auto oppo_he = hfm_.add_he(vid, ti, INVALID_NUM);
          hfm_.set_oppo(line_he.back(), oppo_he);
          pre_id = vid;

        }
      }
      if (pre_id != en_v) {

        line_he.push_back(hfm_.add_he(pre_id, ti, INVALID_NUM));
        auto oppo_he = hfm_.add_he(en_v, ti, INVALID_NUM);
        hfm_.set_oppo(line_he.back(), oppo_he);
      }

      
      for (size_t i = 0; i+1 < line_he.size(); ++i) {
        hfm_.set_next_prev(line_he[i], line_he[i+1]);
        hfm_.set_next_prev(line_he[i+1]->oppo(), line_he[i]->oppo());
      }
      // used for inserting
      hfm_.set_next_prev(line_he.back(), line_he.back()->oppo());
      hfm_.set_next_prev(line_he.front()->oppo(), line_he.front());
      
      // add into triangle he
      add_line_into_he_loops(ti, v2he, line_he);
    }
    
    return 0;
  }

  int find_adj_he(
      size_t fid,
      half_edge::HE_PTR he,
      std::vector<half_edge::HE_PTR>& adj_vh)
  {
    //auto fid = he->ori_f();
    FASSERT(he->ori_f() == fid);
    auto st = he;
    do {
      if (st->ori_f() == fid) {
        adj_vh.push_back(st);
      }
      st = st->oppo();
      if (!half_edge::is_valid(st)) break;
      if (st->ori_f() != fid) { break; }
      st = st->next();
    } while (half_edge::is_valid(st) && st != he);

    if (st != he) {
      st = he->prev()->oppo();
      while (half_edge::is_valid(st)) {
        if (st->ori_f() == fid) {
          adj_vh.push_back(st);
        } else {
          break;
        }
        FASSERT(half_edge::is_valid(st->prev()));
        st = st->prev()->oppo();
      }
    }

    return 0;
  }

  int choose_loop(
      const std::vector<half_edge::HE_PTR>& adj_vh, size_t end_v,
      half_edge::HE_PTR& start_he, half_edge::HE_PTR& end_he)
  {
    for (auto he : adj_vh) {
      auto st = he;
      //FASSERT(he->start_v() == va);
      size_t cnt = 0;
      do {
        st = st->next();
        FASSERT(half_edge::is_valid(st));
        if (st->start_v() == end_v) {
          start_he = he;
          end_he = st;
          return 0;
        }
        ++cnt;
      } while (st != he);
    }
    
    return 1;
  }

  int add_line_into_he_loops(
      size_t tid,
      std::unordered_map<size_t,half_edge::HE_PTR>& v2he,
      const std::vector<half_edge::HE_PTR>& line_he)
  {
    std::vector<half_edge::HE_PTR> line_ip;
    std::vector<half_edge::HE_PTR> loop_ip;

    {
      auto end_he = line_he.back()->oppo()->next();
      auto he = line_he.front();
      
      do {
        auto it = v2he.find(he->start_v());
        if (it != v2he.end()) {
          line_ip.push_back(he);
          loop_ip.push_back(it->second);
        } else {
          v2he[he->start_v()] = he;
        }
        he = he->next();
      } while (he != end_he);
    }
    
    FASSERT(line_ip.size() > 1);

    for (size_t i = 0; i+1 < line_ip.size(); ++i) {
      size_t vb = line_ip[i+1]->start_v();
      
      std::vector<half_edge::HE_PTR> adj_vh;
      // find adj inserted half edge
      find_adj_he(tid, loop_ip[i], adj_vh);
      half_edge::HE_PTR sp(half_edge::invalid_ptr()), ep(half_edge::invalid_ptr());
      
      if (choose_loop(adj_vh, vb, sp, ep)) {
        std::cerr << "# [ ERROR ] can not find the right loop." << std::endl;
        return 1;
      }

      // merge
      half_edge::HE_PTR sl = line_ip[i];
      half_edge::HE_PTR el = line_ip[i+1]->prev(); // must valid

      FASSERT(half_edge::is_valid(el));
      FASSERT(half_edge::is_valid(sl));
      FASSERT(half_edge::is_valid(sp));
      FASSERT(half_edge::is_valid(ep));
      FASSERT(half_edge::is_valid(sp->prev()));
      FASSERT(half_edge::is_valid(sl->oppo()));
      FASSERT(half_edge::is_valid(el->oppo()));
      FASSERT(half_edge::is_valid(ep->prev()));
      
      hfm_.set_next_prev(sp->prev(), sl);
      hfm_.set_next_prev(sl->oppo(), sp);
      hfm_.set_next_prev(ep->prev(), el->oppo());      
      hfm_.set_next_prev(el, ep);
    }
    
    return 0;
  }

  size_t find_tri_edge(size_t f, size_t va, size_t vb)
  {
    const size_t* v3 = seq_.tr().mesh().tv(f);
    for (size_t k = 0; k+1 < 3; ++k) {
      if ( (v3[k]==va && v3[k+1]==vb)
           || (v3[k]==vb && v3[k+1]==va) ) {
        return k;
      }
    }
    if ( (v3[2]==va && v3[0]==vb)
         || (v3[2]==vb && v3[0]==va) ) {
      return 2;
    }

    return INVALID_NUM;
  }

  // input halfedge of each direc-triangle-edge
  int build_bound_oppo(const std::vector<half_edge::HE_PTR>& tri_he3)
  {
    auto& e2f = seq_.e2f();

    for (auto& ef : e2f) {
      size_t va = ef.first.first;
      size_t vb = ef.first.second;

      // std::cout << "   (va,vb): " << va << ", " << vb << std::endl;
      auto& fs = ef.second;
      
      for (size_t i = 0; i+1 < fs.size(); i+=2) {
        size_t fa = fs[i];
        size_t fb = fs[i+1];
        size_t ea = find_tri_edge(fa, va, vb);
        size_t eb = find_tri_edge(fb, va, vb);
        
        FASSERT(ea < 3 && eb < 3);

        std::vector<half_edge::HE_PTR> ea_he, eb_he;
        ea_he.push_back(tri_he3[3*fa+ea]);
        eb_he.push_back(tri_he3[3*fb+eb]);

        FASSERT(ea_he.front()->start_v() == va || ea_he.front()->start_v() == vb);
        FASSERT(eb_he.front()->start_v() == va || eb_he.front()->start_v() == vb);
        FASSERT(ea_he.front()->start_v() != eb_he.front()->start_v());

        size_t ea_end, eb_end;
        if (ea_he.front()->start_v() == va) {
          ea_end = vb;
        } else if (ea_he.front()->start_v() == vb) {
          ea_end = va;
        } else {
          std::cerr << "# [ ERROR ] invalid half-edge." << std::endl;
          return 1;
        }
        
        if (eb_he.front()->start_v() == va) {
          eb_end = vb;
        } else if (eb_he.front()->start_v() == vb) {
          eb_end = va;
        } else {
          std::cerr << "# [ ERROR ] invalid half-edge." << std::endl;
          return 1;
        }

        while (ea_he.back()->end_v() != ea_end) {
          auto pre = ea_he.back();
          FASSERT(half_edge::is_valid(pre->next()));
          ea_he.push_back(pre->next());
          FASSERT(ea_he.back()->end_v() != INVALID_NUM);
        }

        
        while (eb_he.back()->end_v() != eb_end) {
          auto pre = eb_he.back();
          FASSERT(half_edge::is_valid(pre->next()));
          eb_he.push_back(pre->next());
          FASSERT(eb_he.back()->end_v() != INVALID_NUM);
        }
        
        FASSERT(ea_he.size() == eb_he.size());
        
        
        size_t n = ea_he.size();
        for (size_t i = 0; i < ea_he.size(); ++i) {
          hfm_.set_oppo(ea_he[i], eb_he[n-i-1]);
        }
      }
    }
    std::cout << "-- end of building bound oppo." << std::endl;
    return 0;
  }

  int build_triangle_he(
      size_t tid, half_edge::HE_PTR* tri_he)
  {
    auto& bm = seq_.tr().mesh();
    
    std::vector<size_t> eids(3);
    std::vector<char> is_swap(3, 0);
    
    for (size_t k = 0; k < 3; ++k) {
      std::pair<size_t,size_t> e(bm.tv(tid,k), bm.tv(tid,(k+1)%3));
      if (e.first > e.second) {
        std::swap(e.first, e.second);
        is_swap[k] = 1;
      }
      eids[k] = seq_.edge_id(e);
    }

    std::vector<half_edge::HE_PTR> start_end_he(6);

    // add he of tt edge
    for (size_t k = 0; k < 3; ++k) {
      FASSERT(eids[k] != INVALID_NUM);
      auto& s = seq_.seq_tt(eids[k]);
      size_t pre_id = seq_.ev(eids[k]).first; //also is the cut vert id


      std::vector<half_edge::HE_PTR> temp_he;

      // when adding, revert the order
      for (auto& v : s) {
        size_t vid = seq_.cv(v.cvid).id();
        if (vid != pre_id) {
          if (is_swap[k] == 1) {
            temp_he.push_back(hfm_.add_he(pre_id, tid, INVALID_NUM));
          } else {
            temp_he.push_back(hfm_.add_he(vid, tid, INVALID_NUM));
          }
          pre_id = vid;
        }
      }
      
      { // add the last cut vert
        size_t end_id = seq_.ev(eids[k]).second;
        if (pre_id != end_id) {
          if (is_swap[k] == 1) {
            temp_he.push_back(hfm_.add_he(pre_id, tid, INVALID_NUM));
          } else {
            temp_he.push_back(hfm_.add_he(end_id, tid, INVALID_NUM));
          }
        }
      }

      if (is_swap[k] == 1) {
        start_end_he[2*k] = temp_he.front();
        start_end_he[2*k+1] = temp_he.back();
      } else {
        start_end_he[2*k] = temp_he.back();
        start_end_he[2*k+1] = temp_he.front();
      }
      FASSERT(start_end_he[2*k]->start_v() == bm.tv(tid,(k+1)%3));

      for (size_t i = 0; i+1 < temp_he.size(); ++i) {
        if (is_swap[k] == 1) {
          hfm_.set_next_prev(temp_he[i], temp_he[i+1]);
        } else {
          hfm_.set_next_prev(temp_he[i+1], temp_he[i]);
        }
      }
    }

    // add prev and next at each vertex
    for (size_t k = 0; k < 3; ++k) {
      hfm_.set_next_prev(start_end_he[2*((k+1)%3)+1], start_end_he[2*k]);
      tri_he[k] = start_end_he[2*k];
    }

    return 0;
  }

  size_t get_faces_num(size_t sf, size_t ef)
  {
    size_t cnt(0);
    for (size_t fi = sf; fi < ef; ++fi) {
      hfm_.set_flag(fi, 0);
      auto& face_he = hfm_.ori_face_he(fi);
      for (auto& he : face_he) {
        if (he.flag() == 0) {
          auto st = &he;
          do {
            st->flag() = 1;
            st = st->next();
          } while (st != &he);
          ++cnt;
        }
      }
    }

    return cnt;
  }

  int show_faces_num()
  {
    size_t tn = seq_.tr().tn();
    size_t pn = seq_.tr().pn();
    size_t bd_n = get_faces_num(0, tn);
    size_t p_n = get_faces_num(tn, pn+tn);
    std::cout << "--- bound-cut-faces : " << bd_n << std::endl;
    std::cout << "--- plane-cut-faces : " << p_n << std::endl;
    std::cout << "--- all-cut-faces   : " << bd_n+p_n << std::endl;
    return 0;
  }

  int get_faces(size_t sf, size_t ef, std::vector<size_t>& M, std::vector<size_t>& MP)
  {
    MP.push_back(0);

    for (size_t fi = sf; fi < ef; ++fi) {
      hfm_.set_flag(fi, 0);
      auto& face_he = hfm_.ori_face_he(fi);
      for (auto& he : face_he) {
        if (he.flag() == 0) {
          auto st = &he;
          do {
            st->flag() = 1;
            M.push_back(st->start_v());
            st = st->next();
          } while (st != &he);
          MP.push_back(M.size());
        }
      }
    }

    return 0;
  }

  int get_plane_faces(
      std::vector<size_t>& M, std::vector<size_t>& MP)
  {
    size_t tn = seq_.tr().tn();
    size_t pn = seq_.tr().pn();
    return get_faces(tn, pn+tn, M, MP);
  }
  

  int get_bound_faces(std::vector<size_t>& M, std::vector<size_t>& MP)
  {
    size_t tn = seq_.tr().tn();
    return get_faces(0, tn, M, MP);
  }

  struct pt_edge {
    pt_edge() : pt_id(INVALID_NUM), t(INVALID_NUM) { }
    pt_edge(size_t ptid, size_t tid): pt_id(ptid), t(tid) {}
    pt_edge(const pt_edge& pe) {pt_id = pe.pt_id; t = pe.t;}
    pt_edge& operator=(const pt_edge& pe) {pt_id=pe.pt_id; t = pe.t; return *this;}
    size_t pt_id, t;
  };

  struct pp_edge {
    pp_edge(size_t ppid, size_t fid) : pp_id(ppid), f(fid) { }
    size_t pp_id, f;
  };

  int build_plane_faces()
  {
    size_t pn = seq_.tr().pn();
    size_t tn = seq_.tr().tn();

    std::vector<std::list<pt_edge> > pt_id(pn);
    std::vector<std::list<pp_edge> > pp_id(pn);

    // get the boundary pt    
    for (auto& pt : seq_.pt_map()) {
      FASSERT(pt.first.first >= tn);
      FASSERT(pt.first.second < tn);
      size_t pid = pt.first.first-tn;
      pt_id[pid].push_back(pt_edge(pt.second, pt.first.second));
      
    }
    
    // get the inner pp
    for (auto& pp : seq_.pp_map()) {
      size_t pa = pp.first.first-tn;
      size_t pb = pp.first.second-tn;
      FASSERT(pa<pn && pb<pn);
      pp_id[pa].push_back(pp_edge(pp.second,pb));//push pp line id
      pp_id[pb].push_back(pp_edge(pp.second,pa));
    }

    // clean inner edges
#pragma omp parallel for    
    for (size_t pi = 0; pi < pn; ++pi) {
      std::list<pp_edge> new_pp;
      std::set<std::pair<size_t,size_t> > line_set;
      for (const auto& ppe : pp_id[pi]) {
        size_t ppid = ppe.pp_id;
        size_t a = seq_.cv(seq_.seq_pp(ppid).front().cvid).id();
        size_t b = seq_.cv(seq_.seq_pp(ppid).back().cvid).id();
        if (a > b) {std::swap(a, b);}
        std::pair<size_t,size_t> bd(a,b);
        if (line_set.find(bd) == line_set.end()) {
          line_set.insert(bd);
          new_pp.push_back(ppe);
        }
      }
      std::swap(new_pp, pp_id[pi]);
    }

#pragma omp parallel for
    for (size_t pi = 0; pi < pn; ++pi) {
      if (build_plane_bound(pi, pt_id[pi])) {
        std::cerr << "... build plane bound error, pid " << pi << std::endl;
      }

      hfm_.set_flag(pi+tn, 0); // set not visited
      
      if (build_plane_inner(pi, pp_id[pi])) {
        std::cerr << "... build plane inner error, pid " << pi << std::endl;
      }
    }

    return build_plane_foppo();
  }

  int build_plane_foppo()
  {
    size_t pn = seq_.tr().pn();
    size_t tn = seq_.tr().tn();
    
#pragma omp parallel for
    for (size_t pi = 0; pi < pn; ++pi) {
      auto& hes = hfm_.ori_face_he(pi+tn);
      for (auto& he : hes) {
        auto fop_he_ptr = hfm_.add_he(he.end_v(), tn+pn+pi, he.inter_f());
        hfm_.set_foppo(&he, fop_he_ptr);
      }
      for (auto& he : hes) {
        FASSERT(half_edge::is_valid(he.prev()));
        if (!half_edge::is_valid(he.prev())) { // for test
            std::cout << "... [ ERROR ] he :  " << he.start_v() << ", " << he.end_v()
                      << " | pid " << pi << std::endl;
        }
        
        hfm_.set_next_prev(he.foppo(), he.prev()->foppo());
        if (half_edge::is_valid(he.oppo())) {
          FASSERT(half_edge::is_valid(he.foppo()));
          FASSERT(half_edge::is_valid(he.oppo()->foppo()));
          hfm_.set_oppo(he.foppo(), he.oppo()->foppo());
        }
      }
    }
    return 0;
  }


  int remove_plane_isolated_half_edge(size_t pid)
  {
    auto& hes = hfm_.ori_face_he(pid);

    hfm_.set_flag(pid, 'o'); // set not visited

    for (auto& he : hes) {
      if (he.flag() == 'o') {
        auto st = &he;
        do {
          st->flag() = 'u';
          st = st->next();
          if (!half_edge::is_valid(st)) break;
        } while (st != &he);
        if (st != &he) { // delete
          auto st1 = &he;
          do {
            st1->flag() = 'd';
            st1 = st1->next();
            if (!half_edge::is_valid(st1)) break;
          } while (st1 != &he);
          auto st2 = &he;
          do {
            st2->flag() = 'd';
            st2 = st2->prev();
            if (!half_edge::is_valid(st2)) break;
          } while (st2 != &he);
        }
      }//if
    }//for

    std::vector<std::list<half_edge>::iterator> del_he;
    for (auto it = hes.begin(); it != hes.end(); ++it) {
      if (it->flag() == 'd') {
        del_he.push_back(it);
      }
    }

    if (del_he.size() > 0) {
      std::cout << "-- delete isolated half-edges : " << del_he.size() << std::endl;
      for (auto& it : del_he) {
        hes.erase(it);
      }
    }

    return 0;
  }

  int build_plane_bound(
      size_t pid, const std::list<pt_edge>& pt_ids)
  {
    size_t tn = seq_.tr().tn();
    auto& P = seq_.tr().plane();
    auto& bm = seq_.tr().mesh();
    auto& V = bm.verts();

    size_t vn = seq_.tr().vn();

    // delete half-edge of the common faces between the plane and triangles
    std::map<std::pair<size_t,size_t>, std::list<pt_edge> > bd_loop_list;

    std::vector<std::pair<size_t,size_t> > pt_bd(pt_ids.size());

    std::unordered_map<size_t,size_t> v2e; // edge-type cv to eid

    // build v2he
    for (auto& pte : pt_ids) {
      size_t id = pte.pt_id;
      char d = seq_.pt_dir(id);
      auto& bd = seq_.pt_bd(id);
      size_t va = seq_.cv(bd.first).id();
      size_t vb = seq_.cv(bd.second).id();

      if (va >= vn) {
        FASSERT(seq_.cv(bd.first).is_edge());
        v2e.insert(std::make_pair(va, seq_.cv(bd.first).ele_id()));
      }
      if (vb >= vn) {
        FASSERT(seq_.cv(bd.second).is_edge());
        v2e.insert(std::make_pair(vb, seq_.cv(bd.second).ele_id()));
      }
      
      if (d=='-') { std::swap(va, vb); }

      auto he = std::make_pair(va, vb);
      auto it = bd_loop_list.find(he);
      if (it == bd_loop_list.end()) {
        std::list<pt_edge> pt_l;
        pt_l.push_back(pte);
        bd_loop_list.insert(std::make_pair(he, pt_l));
      } else {
        it->second.push_back(pte);
      }

    }

    // remove common half-edges of plane and triangles, <va,vb,fid,eid>
    std::vector<std::pair<std::pair<size_t,size_t>, std::pair<size_t,size_t> > > common_he;
    std::unordered_map<size_t, std::pair<size_t,size_t> > tri_edges;
    for (auto& bd : bd_loop_list) {
      auto he = bd.first;
      if (he.first < vn && he.second < vn) { // edge
        auto e = he;
        if (e.first > e.second) { std::swap(e.first, e.second); }
        size_t eid = seq_.edge_id(e);
        FASSERT(eid != INVALID_NUM);
        auto& fs = seq_.edge_fs(eid);

        tri_edges[eid] = e;

        for (auto& f : fs) {
          if (seq_.tr().pt(pid, f) == 3) {//coplanar triangle
            bool is_same = is_triangle_and_plane_has_same_norm(
                &V[3*bm.tv(f,0)],&V[3*bm.tv(f,1)],&V[3*bm.tv(f,2)],&P[4*pid]);
            if (is_same) {
              for (size_t i = 0; i<3; ++i) {
                size_t va = bm.tv(f,i);
                size_t vb = bm.tv(f,(i+1)%3);
                if (va==he.first&&vb==he.second) {
                  common_he.push_back(std::make_pair(std::make_pair(va, vb), std::make_pair(f,eid)));
                }
              }
            } else {
              for (size_t i = 0; i<3; ++i) {
                size_t va = bm.tv(f,i);
                size_t vb = bm.tv(f,(i+1)%3);
                if (va==he.second&&vb==he.first) {
                  common_he.push_back(std::make_pair(std::make_pair(vb, va), std::make_pair(f,eid)));
                }
              }

            }
          }
        }
      }
    }// for
    // remove common he
    for (auto& he : common_he) {
      size_t eid = he.second.second;
      size_t fid = he.second.first;
      auto& e = he.first;

      // when the adjacent face in the same domain has the same direction
      std::unordered_map<size_t, std::pair<size_t,size_t> > f2de;

      {
        auto bd_it = bd_loop_list.find(e);
        if (bd_it != bd_loop_list.end()) {
          for (auto it = bd_it->second.begin(); it != bd_it->second.end(); ++it) {
            f2de.insert(std::make_pair(it->t, std::make_pair(e.first, e.second)));
          }//for
        }//if
      }
      {
        auto op_bd_it = bd_loop_list.find(std::make_pair(e.second, e.first));
        if (op_bd_it != bd_loop_list.end()) {
          for (auto it = op_bd_it->second.begin(); it != op_bd_it->second.end(); ++it) {
            f2de.insert(std::make_pair(it->t, std::make_pair(e.second, e.first)));
          }
        }
      }

      auto& fs = seq_.edge_fs(eid);

      std::unordered_set<size_t> del_f;
      for (size_t i = 0; i < fs.size(); ++i) {
        if (fid == fs[i]) {
          size_t adj_id = (i%2==0) ? ((i+1)%fs.size()) : ((i+fs.size()-1)%fs.size());
          size_t adj_f = fs[adj_id];
          auto it = f2de.find(adj_f);
          if (it!=f2de.end() && it->second.first == e.first) { // the same direction, delete the adjacent he
            del_f.insert(adj_f);
          } else {// del other he
            for (size_t j = 0; j < fs.size(); ++j) {
              if (j != i && j != adj_id) {
                del_f.insert(fs[j]);
              }
            }
          }
          break;
        }
      }
      
      {
        auto bd_it = bd_loop_list.find(e);
        if (bd_it != bd_loop_list.end()) {
          std::vector<std::list<pt_edge>::iterator > del_it;
          for (auto it = bd_it->second.begin(); it != bd_it->second.end(); ++it) {
            if (del_f.count(it->t) > 0) del_it.push_back(it);
          }//for
          for (auto& it : del_it) { bd_it->second.erase(it); }
        }//if
      }
      {
        auto op_bd_it = bd_loop_list.find(std::make_pair(e.second, e.first));
        if (op_bd_it != bd_loop_list.end()) {
          std::vector<std::list<pt_edge>::iterator > del_it;
          for (auto it = op_bd_it->second.begin(); it != op_bd_it->second.end(); ++it) {
            if (del_f.count(it->t) > 0) del_it.push_back(it);
          }
          for (auto& it : del_it) { op_bd_it->second.erase(it); }
        }
      }
    }

    // deal with boundary half-edge in triangle edges
    for (auto& te : tri_edges) {
      auto& e = te.second;
      auto& fs = seq_.edge_fs(te.first);
      
      // build f2de
      // when the adjacent face in the same domain has the same direction
      std::unordered_map<size_t, std::pair<size_t,size_t> > f2de;

      {
        auto bd_it = bd_loop_list.find(e);
        if (bd_it != bd_loop_list.end()) {
          for (auto it = bd_it->second.begin(); it != bd_it->second.end(); ++it) {
            f2de.insert(std::make_pair(it->t, std::make_pair(e.first, e.second)));
          }//for
        }//if
      }
      {
        auto op_bd_it = bd_loop_list.find(std::make_pair(e.second, e.first));
        if (op_bd_it != bd_loop_list.end()) {
          for (auto it = op_bd_it->second.begin(); it != op_bd_it->second.end(); ++it) {
            f2de.insert(std::make_pair(it->t, std::make_pair(e.second, e.first)));
          }
        }
      }

      std::unordered_set<size_t> del_f;

      for (size_t i = 0; i+1 < fs.size(); i+=2) {
        auto ita = f2de.find(fs[i]);
        auto itb = f2de.find(fs[i+1]);


        if (ita != f2de.end() && itb != f2de.end() && ita->second.first != itb->second.first) {
          size_t va = seq_.tr().mesh().tv(fs[i], 0);
          size_t vb = seq_.tr().mesh().tv(fs[i], 1);
          size_t vc = seq_.tr().mesh().tv(fs[i], 2);
          std::unordered_set<size_t> vs;
          vs.insert(va);
          vs.insert(vb);
          vs.insert(vc);
          size_t vd(INVALID_NUM);
          for (size_t vj = 0; vj < 3; ++vj) {
            size_t fb_vid = seq_.tr().mesh().tv(fs[i+1], vj);
            if (vs.count(fb_vid)==0) { vd=fb_vid; break; }
          }
          FASSERT(vd != INVALID_NUM);
          int f_sgn = vol_sgn(&V[3*vd],&V[3*va],&V[3*vb],&V[3*vc]);
          if (f_sgn > 0) {
            del_f.insert(fs[i]);
            del_f.insert(fs[i+1]);

          }
        }
      }

      {
        auto bd_it = bd_loop_list.find(e);
        if (bd_it != bd_loop_list.end()) {
          std::vector<std::list<pt_edge>::iterator > del_it;
          for (auto it = bd_it->second.begin(); it != bd_it->second.end(); ++it) {
            if (del_f.count(it->t) > 0) del_it.push_back(it);
          }//for
          for (auto& it : del_it) { bd_it->second.erase(it); }
        }//if
      }
      {
        auto op_bd_it = bd_loop_list.find(std::make_pair(e.second, e.first));
        if (op_bd_it != bd_loop_list.end()) {
          std::vector<std::list<pt_edge>::iterator > del_it;
          for (auto it = op_bd_it->second.begin(); it != op_bd_it->second.end(); ++it) {
            if (del_f.count(it->t) > 0) del_it.push_back(it);
          }
          for (auto& it : del_it) { op_bd_it->second.erase(it); }
        }
      }
    }


    std::map<std::pair<size_t,size_t>, pt_edge > bd_loop;

    for (auto& e : bd_loop_list) {
      if (e.second.size() > 0) {
        bd_loop[e.first] = e.second.front();
      }
    }

    // build graph and add prev, next in each he
    std::unordered_map<size_t, std::vector<half_edge::HE_PTR> > v2he;
    std::map<std::pair<size_t,size_t>, std::pair<half_edge::HE_PTR,half_edge::HE_PTR> > bd_he;
    
    for (const auto& he : bd_loop) {
      auto& bd = he.first;
      auto& Q = seq_.seq_pt(he.second.pt_id);
      char d = seq_.pt_dir(he.second.pt_id);

      size_t pre_id(bd.first);
      if (d == '-') {
        pre_id = bd.second;
      }
      std::vector<half_edge::HE_PTR> hes;
      for (auto& v : Q) {
        size_t vid = seq_.cv(v.cvid).id();
        if (pre_id != vid) {
          if (d == '+') {
            FASSERT(he.second.t < tn);
            hes.push_back(hfm_.add_he(pre_id, tn+pid, he.second.t));
            hes.back()->flag() = 0;
          } else {
            FASSERT(he.second.t < tn);
            hes.push_back(hfm_.add_he(vid, tn+pid, he.second.t));
            hes.back()->flag() = 0;
          }
          pre_id = vid;
        }
      }
      size_t en_v(bd.second);
      if (d == '-') {
        en_v = bd.first;
      }
      
      if (en_v != pre_id) {
        if (d == '+') {
          hes.push_back(hfm_.add_he(pre_id, tn+pid, he.second.t));
          hes.back()->flag() = 0;
        } else {
          hes.push_back(hfm_.add_he(en_v, tn+pid, he.second.t));
          hes.back()->flag() = 0;
        }
      }

      if (d == '+') {      
        for (size_t i = 0; i+1 < hes.size(); ++i) {
          hfm_.set_next_prev(hes[i], hes[i+1]);
        }
      } else {
        for (size_t i = 0; i+1 < hes.size(); ++i) {
          hfm_.set_next_prev(hes[i+1], hes[i]);
        }
      } // if

      std::pair<half_edge::HE_PTR, half_edge::HE_PTR> start_end_he;
      if (d=='+') {
        start_end_he = std::make_pair(hes.front(), hes.back());
      } else {
        start_end_he = std::make_pair(hes.back(), hes.front());
      }
      bd_he[bd] = start_end_he;

      { // for building graph
        size_t va = bd.first;
        size_t vb = bd.second;
        // va
        auto it = v2he.find(va);
        if (it != v2he.end()) {
          it->second.push_back(start_end_he.first);
        } else {
          std::vector<half_edge::HE_PTR> hes;
          hes.push_back(start_end_he.first);
          v2he.insert(std::make_pair(va, hes));
        }
        // vb
        auto it2 = v2he.find(vb);
        if (it2 != v2he.end()) {
          it2->second.push_back(start_end_he.second);
        } else {
          std::vector<half_edge::HE_PTR> hes;
          hes.reserve(2);
          hes.push_back(start_end_he.second);
          v2he.insert(std::make_pair(vb, hes));
        }
      }
    } // for

    {// set oppo if need
      std::vector<std::pair<size_t,size_t> > double_he;
      for (const auto& bd : bd_loop) {
        auto e = bd.first;
        if (e.first < e.second) {
          std::swap(e.first, e.second);
          auto it = bd_loop.find(e);
          if (it != bd_loop.end()) {
            double_he.push_back(e);
          }
        }
      }
      for (const auto& he : double_he) {
        const auto& ptr_a = bd_he[he];
        const auto& ptr_b = bd_he[std::make_pair(he.second, he.first)];
        auto sa = ptr_a.first;
        auto sb = ptr_b.second;
        while (sa != ptr_a.second) {
          hfm_.set_oppo(sa, sb);
          sa = sa->next();
          sb = sb->prev();
        }
        hfm_.set_oppo(sa, sb);
      }
    }

    // add prev, next between he, need sort if needed
    for (const auto& c : v2he) {
      size_t v = c.first;
      auto& adj_he = c.second;

      if (adj_he.size() == 2) {
        if (adj_he.front()->start_v() == v) {
          hfm_.set_next_prev(adj_he.back(), adj_he.front());
        } else if (adj_he.back()->start_v() == v) {
          hfm_.set_next_prev(adj_he.front(), adj_he.back());
        } else {
          std::cerr << "# [ ERROR ] not right." << std::endl;
          return 1;
        }
      } else if (adj_he.size() > 2) { // need sorting
        FASSERT(((adj_he.size()&1)==0));
        
        
        if (v >= vn) { // edge cv, do not need to sort
          std::unordered_map<size_t,half_edge::HE_PTR> f2he;
          for (size_t i = 0; i < adj_he.size(); ++i) {
            auto& he = adj_he[i]; // HE
            // get the triangle id
            size_t tid = he->inter_f();
            FASSERT(tid < tn);
            f2he[tid] = he;
          }

          // get the edge id

          size_t eid = INVALID_NUM;
          {
            auto it = v2e.find(v);
            if (it != v2e.end()) {
              eid = it->second;
            } else {
              std::cerr << "# [ ERROR ] invalid edge id." << std::endl;
              return 1;
            }
          }
          
          auto& fs = seq_.edge_fs(eid);
          FASSERT((fs.size()&1)==0);
          
          for (size_t i = 0; i+1 < fs.size(); i+=2) {
            FASSERT(f2he.find(fs[i])!=f2he.end());
            FASSERT(f2he.find(fs[i+1])!=f2he.end());
            auto a = f2he[fs[i]];
            auto b = f2he[fs[i+1]];
            
            if (a->start_v() == v) {
              hfm_.set_next_prev(b, a);
            } else if (b->start_v() == v) {
              hfm_.set_next_prev(a, b);
            } else {
              std::cerr << "# [ ERROR ] sort cut edge of edge-type cv." << std::endl;
              return 1;
            }
          }

        } else { // vert-type cv, sort the adjacent cut edges
          
          std::vector<size_t> cut_fs;
          for (auto& he : adj_he) {
            FASSERT(he->inter_f() < tn);
            cut_fs.push_back(he->inter_f()); //tid
            FASSERT(seq_.tr().pt(pid, he->inter_f())==2);
          }


          std::unordered_map<size_t,size_t> fs2id;
          for (size_t i = 0; i < cut_fs.size(); ++i) {
            fs2id.insert(std::make_pair(cut_fs[i], i));
          }
          
          
          std::vector<MPF_VEC> dt;
          seq_.cal_pt_outer_direc(v, pid, cut_fs, dt);

          std::vector<char> cp;
          std::vector<short> dir;
          seq_.cycle_sort_triangles(v, pid, dt, fs2id, cut_fs, cp, dir, cycle_sort_type::SIMPLE);

          sorted_cv_[pid].insert(v);
          size_t np_cnt(0);
          size_t n = cut_fs.size();
          for (size_t i = 0; i < n+1; ++i) {
            size_t j = (i+1)%n;
            if (dir[i]==-1 && dir[j]==1) {
              auto ha = adj_he[fs2id[cut_fs[i]]];
              auto hb = adj_he[fs2id[cut_fs[j]]];
              FASSERT(hb->start_v() == v);
              hfm_.set_next_prev(ha, hb);
              ++np_cnt;
            }
          } // for
          FASSERT(np_cnt*2 == n);
        } // if, sorting
      } else if (adj_he.size() == 1) {
        std::cout << "# [ ERROR ] there exists one dangle half-edge." << std::endl;
      } // if, classify adj he num
    }

    return 0;
  }

  int build_plane_inner(
      size_t pid, const std::list<pp_edge>& pp_ids)
  {
    size_t tn = seq_.tr().tn();
    // build v2he
    auto& plane_HE = hfm_.ori_face_he(pid+tn);
    std::unordered_map<size_t, std::vector<half_edge::HE_PTR> > v2he;

    // now all the HE are on the bound
    for (auto& he : plane_HE) {
      auto it = v2he.find(he.start_v());
      if (it != v2he.end()) {
        it->second.push_back(&he);
      } else {
        std::vector<half_edge::HE_PTR> vec;
        vec.reserve(2);
        vec.push_back(&he);
        v2he.insert(std::make_pair(he.start_v(), vec));
      }
    }

    // insert pp lines, each time add one line
    for (auto it = pp_ids.begin(); it != pp_ids.end(); ++it) {
      build_plane_inner(v2he, pid, it->pp_id, it->f);
    }// for
    
    return 0;
  }

  int build_plane_inner(
      std::unordered_map<size_t, std::vector<half_edge::HE_PTR> >& v2he,
      size_t pid, size_t pp_id, size_t line_pid)
  {
    uint8_t pp_direc = (pid < line_pid)? 0 : 1;

    size_t tn = seq_.tr().tn();
    auto& Q = seq_.seq_pp(pp_id);
    auto& S = seq_.seq_pp_s(pp_id);

    std::vector<half_edge::HE_PTR> hes;

    FASSERT(Q.size() == S.size());


    for (size_t i = 0; i+1 < Q.size(); ++i) {//each cv in pp line
      if (S[i] == 0) continue;

      FASSERT(S[i] == 1);
      size_t vid = seq_.cv(Q[i].cvid).id();
      size_t nex_vid = seq_.cv(Q[i+1].cvid).id();


      FASSERT(vid != nex_vid);
      hes.push_back(hfm_.add_he(vid, pid+tn, line_pid+tn));
      hes.back()->flag() = (pp_direc==0)?0:1;
      auto op_he = hfm_.add_he(nex_vid, pid+tn, line_pid+tn);
      op_he->flag() = (pp_direc==0)?1:0;
      hfm_.set_oppo(hes.back(), op_he);
    }

    
    for (size_t i = 0; i+1 < hes.size(); ++i) {
      if(hes[i]->end_v() == hes[i+1]->start_v()) {
        hfm_.set_next_prev(hes[i], hes[i+1]);
        hfm_.set_next_prev(hes[i+1]->oppo(), hes[i]->oppo());
      }
    }

    size_t pre_i = 0;
    for (size_t i = 0; i < hes.size(); ++i) {
      if (!half_edge::is_valid(hes[i]->next())) {
        hfm_.set_next_prev(hes[i], hes[i]->oppo());
        hfm_.set_next_prev(hes[pre_i]->oppo(), hes[pre_i]);
        add_subline_into_plane_loop(v2he, pid, hes[pre_i], hes[i], line_pid);
        pre_i = i+1;
      }
    }

    return 0;
  }


  int add_subline_into_plane_loop(
      std::unordered_map<size_t, std::vector<half_edge::HE_PTR> >& v2he,
      size_t pid, const half_edge::HE_PTR& line_start_he, const half_edge::HE_PTR& line_end_he,
      size_t line_pid)
  {
    typedef std::unordered_map<size_t, std::vector<half_edge::HE_PTR> >::iterator map_iter_type;
    std::vector<half_edge::HE_PTR> line_ip;
    std::vector<map_iter_type> loop_ip;

    FASSERT(line_end_he->oppo() == line_end_he->next());

    {
      auto he = line_start_he;
      do {
        size_t v = he->start_v();
        auto it = v2he.find(v);
        if (it != v2he.end()) {
          line_ip.push_back(he);
          loop_ip.push_back(it);
        } else {
          // cv in the middle of line, only need consider one he
          std::vector<half_edge::HE_PTR> vec(1);
          vec[0] = he;
          v2he.insert(std::make_pair(v, vec));
        }
        he = he->next();
      } while (he != line_end_he->oppo()->next());
    }

    FASSERT(line_ip.size() > 1);

    
    // find all half-edge fan
    for (size_t i = 0; i+1 < line_ip.size(); ++i) {
      std::unordered_set<half_edge::HE_PTR> s_fan, e_fan;

      //each time, need to search it with the current graph
      find_plane_adj_he(loop_ip[i]->second, s_fan);
      find_plane_adj_he(loop_ip[i+1]->second, e_fan);

      half_edge::HE_PTR sl = line_ip[i];
      half_edge::HE_PTR el = line_ip[i+1]->prev();


      // must know the direction of each fan
      half_edge::HE_PTR sp(half_edge::invalid_ptr()), ep(half_edge::invalid_ptr());
      if (choose_plane_fan(s_fan, e_fan, pid, line_pid, sp, ep)) {
        std::cerr << "# [ ERROR ] can not find appropriate fans." << std::endl;
        exit(1);
      }

      hfm_.set_next_prev(sp->prev(), sl);
      hfm_.set_next_prev(sl->oppo(), sp);
      hfm_.set_next_prev(ep->prev(), el->oppo());
      hfm_.set_next_prev(el, ep);
    }

    return 0;
  }

  int choose_plane_fan(
      const std::unordered_set<half_edge::HE_PTR>& s_fan,
      const std::unordered_set<half_edge::HE_PTR>& e_fan,
      size_t pid, size_t line_pid,
      half_edge::HE_PTR& start_he, half_edge::HE_PTR& end_he)
  {
    if (s_fan.size() == 0 || e_fan.size() == 0) {return 1;}

    
    if (s_fan.size() == 1 && e_fan.size() == 1) {
      start_he = *s_fan.begin();
      end_he = *e_fan.begin();
      return 0;
    }

    std::set<std::pair<half_edge::HE_PTR,half_edge::HE_PTR> > valid_fan_pair;
    
    {
      size_t end_v = (*e_fan.begin())->start_v();    
      for (auto& sf : s_fan) {
        std::vector<half_edge::HE_PTR> end_hes;
        FASSERT(half_edge::is_valid(sf->next()));
        for (auto st=sf->next(); st!=sf; st=st->next()) {
          FASSERT(half_edge::is_valid(st->next()));
          
          if (!half_edge::is_valid(st->next())) {//for test
            std::cout << "... error : st verts  ("
                      << st->start_v() << ", "
                      << st->end_v() << ")" << std::endl;
          }
          
          if (st->start_v() == end_v) {
            end_hes.push_back(st);
          }
        }
        for (auto he : end_hes) {
          valid_fan_pair.insert(std::make_pair(sf, he));
        }
      }
    }


    uint8_t pp_direc = (pid<line_pid) ? 0 : 1;

    if (valid_fan_pair.size() == 1) {
      start_he = valid_fan_pair.begin()->first;
      end_he = valid_fan_pair.begin()->second;
      return 0;
    } else {
      for (auto& pf : valid_fan_pair) {
        if (judge_fan(pf.first, pid, line_pid, pp_direc)
            && judge_fan(pf.second, pid, line_pid, (pp_direc^1))) {
          start_he = pf.first;
          end_he = pf.second;
          return 0;
        }
      }
      
    }


    {
      // judge the start he
      start_he = half_edge::invalid_ptr();
      end_he = half_edge::invalid_ptr();

      FASSERT(!s_fan.empty() && !e_fan.empty());

      if (s_fan.size() == 1) {
        start_he = *s_fan.begin();
      } else {
        for (auto& he : s_fan) {
          if (judge_fan(he, pid, line_pid, pp_direc)) {
            start_he = he;
            break;
          }
        }
      }

      if (e_fan.size() == 1) {
        end_he = *e_fan.begin();
      } else {
        for (auto& he : e_fan) {
          if (judge_fan(he, pid, line_pid, (pp_direc^1))) {
            end_he = he;
            break;
          }
        }
      }

      
      if (half_edge::is_valid(start_he)
          && half_edge::is_valid(end_he)) {
        return 0;
      } else {
        return 1;
      }
    }
    
    return 1;
  }

  void cal_line_direc(const MPF_VEC& N, size_t f, MPF_VEC& d)
  {
    if (f < seq_.tr().tn()) {
      size_t va = seq_.tr().mesh().tv(f,0);
      size_t vb = seq_.tr().mesh().tv(f,1);
      size_t vc = seq_.tr().mesh().tv(f,2);
      auto& V = seq_.tr().mesh().verts();
      MPF_VEC A, B, C;
      A << V[3*va], V[3*va+1], V[3*va+2];
      B << V[3*vb], V[3*vb+1], V[3*vb+2];
      C << V[3*vc], V[3*vc+1], V[3*vc+2];
      d = N.cross((B-A).cross(C-A));
    } else {
      MPF_VEC Np;
      auto& P = seq_.tr().plane();
      size_t pid = f-seq_.tr().tn();
      Np << P[4*pid], P[4*pid+1], P[4*pid+2];
      d = N.cross(Np);
    }
  }
  
  bool judge_fan(half_edge::HE_PTR he, size_t pa, size_t pb, uint8_t L_d)
  {
    // for cound the sorted vert by using cycle-sorting
    sorted_cv_[pa].insert(he->start_v());
    
    auto pre_he = he->prev();
    
    int p_d = (pre_he->flag()==0)?1:(-1);
    int d = (he->flag()==0)?1:(-1);
    
    MPF_VEC N, di, dj, dk;
    auto& P = seq_.tr().plane();
    N << P[4*pa], P[4*pa+1], P[4*pa+2];
    
    FASSERT(pre_he->ori_f() == pa+seq_.tr().tn());
    FASSERT(he->ori_f() == pa+seq_.tr().tn());
    
    cal_line_direc(N, pre_he->inter_f(), di);
    cal_line_direc(N, he->inter_f(), dj);
    cal_line_direc(N, pb+seq_.tr().tn(), dk);
    
    int ks = (L_d == 0) ? 1 : -1;
    
    di *= p_d*(-1);
    dj *= d;
    
    int sij = sgn(di.cross(dj).dot(N));
    int dij = sgn(di.dot(dj));

    
    if (sij != 0 || (sij==0 && dij<0)) {
      if (seq_.judge_fan(N, di, dj, dk, ks)==1) {
        return true;
      }
    } else {
      if (seq_.judge_full_fan(N, di, dk, ks)==1) {
        return true;
      }
    }
    
    return false;
  }

  int find_plane_adj_he(
      const std::vector<half_edge::HE_PTR>& hes,
      std::unordered_set<half_edge::HE_PTR>& fan)
  {
    for (auto he : hes) {
      auto st = he;
      do {
        fan.insert(st);
        st = st->oppo();
        if (!half_edge::is_valid(st)) {break;}
        st = st->next();
      } while (half_edge::is_valid(st) && st != he);

      if (st != he) {
        st = he->prev()->oppo();
        while (half_edge::is_valid(st)) {
          fan.insert(st);
          st = st->prev()->oppo();
        }
      }
    }
    
    return 0;
  }

  // format for each cell: cn, fn, v1n, (v1,v2,...,v1n), v2n, (....), v3n, (...), ..., vfn.
  int get_cells(std::vector<size_t>& CC, size_t& cells_num)
  {
    hfm_.set_all_flag(0);

    size_t all_cn = seq_.tr().pn()*2+seq_.tr().tn();

    cells_num = 0;

    for (size_t ci = 0; ci < all_cn; ++ci) {
      auto& C = hfm_.ori_face_he(ci);
      for (auto& he : C) {
        if (he.flag() == 1) continue;
        
        std::queue<half_edge::HE_PTR> que;
        que.push(&he);

        std::vector<size_t> vs;
        std::vector<size_t> f_vn;

        bool is_valid = true;
        
        while (!que.empty()) {
          auto start_he = que.front();
          que.pop();          

          if (!half_edge::is_valid(start_he->oppo())) {
            is_valid = false;
          }
          if (start_he->flag() == 1) {continue;}
          
          auto st = start_he;

          do {
            st->flag() = 1;
            vs.push_back(st->start_v());

            auto op_st = st->oppo();

            if (!half_edge::is_valid(op_st)) {
              is_valid = false;
            }
      
            if (op_st->flag() == 0 && half_edge::is_valid(op_st)) {
              que.push(op_st);
            }
            
            st = st->next();
          } while (st != start_he);


          f_vn.push_back(vs.size());
        } // while

        if (is_valid) {
          // add cells
          CC.push_back(f_vn.size()+vs.size()+1);
          CC.push_back(f_vn.size());

          size_t pre_id = 0;
          for (size_t i = 0; i < f_vn.size(); ++i) {
            CC.push_back(f_vn[i]-pre_id);
            for (size_t j = pre_id; j<f_vn[i]; ++j) {
              CC.push_back(vs[j]);
            }
            pre_id = f_vn[i];
          }
          ++cells_num;
        } else {
          std::cout << "# [ WARNING ] topology degenerated cell." << std::endl;
        }
      }
    }
    
    std::cout << "-- cells num: " << cells_num << std::endl;
    
    return 0;
  }


  int get_all_passing_faces(
      std::vector<std::vector<size_t> >& cv2passF,
      std::vector<size_t>& cv2ele)
  {
    size_t cvn = seq_.ori_cv_num();
    size_t pure_cvn = seq_.pure_cv_num();
    
    cv2passF.resize(pure_cvn);
    cv2ele.resize(pure_cvn);
    
    for (size_t vi = 0; vi < cvn; ++vi) {
      auto& ori_CV = seq_.cv(vi);
      for (size_t j = 0; j < 3; ++j) {
        cv2passF[ori_CV.id()].push_back(ori_CV.f(j));
      }
      cv2ele[ori_CV.id()] = ori_CV.ele_id();
    }

    size_t vn = seq_.tr().vn();
    for (size_t vi = 0; vi < vn; ++vi) {
      auto& fs = seq_.v2f(vi);
      cv2passF[vi].insert(cv2passF[vi].end(), fs.begin(), fs.end());
    }
    
    return 0;
  }

  int merge_all_info_with_ufs(
      const std::vector<std::vector<size_t> >& cv2nearF,
      const std::vector<std::vector<size_t> >& cv2passF,
      union_find_set& ufs,
      std::vector<size_t>& cv2newid,
      std::vector<size_t>& new2ori,
      std::vector<std::vector<size_t> >& new_cv2checkF)
  {
    size_t cvn = cv2nearF.size();
    size_t new_cvn(0);


    cv2newid.resize(cvn, INVALID_NUM);
    new2ori.reserve(cvn);
    for (size_t vi = 0; vi < cvn; ++vi) {

      size_t pa = ufs.find(vi);

      if (cv2newid[pa]==INVALID_NUM) {
        cv2newid[pa] = new_cvn;
        new2ori.push_back(pa);
        ++new_cvn;
      }
      cv2newid[vi] = cv2newid[pa];
    }//for
    std::cout << "-- cut vertices after merging : " << new_cvn << std::endl;

    // get merged cv info
    {
      std::vector<std::vector<size_t> > new_cv2nearF(new_cvn);
      std::vector<std::vector<size_t> > new_cv2passF(new_cvn);
      
      for (size_t vi = 0; vi < cvn; ++vi) {
        size_t new_id = cv2newid[vi];

        
        new_cv2nearF[new_id].insert(
            new_cv2nearF[new_id].end(),
            cv2nearF[vi].begin(), cv2nearF[vi].end());
        
        new_cv2passF[new_id].insert(
            new_cv2passF[new_id].end(),
            cv2passF[vi].begin(), cv2passF[vi].end());
      }//for
      
      new_cv2checkF.resize(new_cvn);
      
      for (size_t vi = 0; vi < new_cvn; ++vi) {
        std::sort(new_cv2nearF[vi].begin(), new_cv2nearF[vi].end());
        new_cv2nearF[vi].erase(std::unique(new_cv2nearF[vi].begin(), new_cv2nearF[vi].end()),
                           new_cv2nearF[vi].end());
        
        std::sort(new_cv2passF[vi].begin(), new_cv2passF[vi].end());
        new_cv2passF[vi].erase(std::unique(new_cv2passF[vi].begin(), new_cv2passF[vi].end()),
                           new_cv2passF[vi].end());
        
        vector_minus(new_cv2nearF[vi], new_cv2passF[vi], new_cv2checkF[vi]);
      }
    }
    
    return 0;
  }

  int move_cv(
      const std::vector<size_t>& checkF,
      const std::pair<Eigen::Matrix<mpf_class,3,1>,mpf_class>& c,
      double* dc)
  {
    size_t tn = seq_.tr().tn();
    auto& bm = seq_.tr().mesh();
    auto& P = seq_.tr().plane();

    Eigen::Map<Eigen::Matrix<double,3,1> > X(dc, 3, 1);

    for (auto& f : checkF) {
      Eigen::Matrix<double, 3, 1> N;
      double df(0.0);
      Eigen::Matrix<mpf_class, 3, 1> mpf_N;
      mpf_class mpf_df(0.0);
      
      if (f >= tn) {
        size_t pi = f-tn;
        N << P[4*pi], P[4*pi+1], P[4*pi+2];
        df = P[4*pi+3];
        mpf_N << P[4*pi], P[4*pi+1], P[4*pi+2];
        mpf_df = P[4*pi+3];
      } else {
        auto va_ptr = bm.tv(f);
        get_face_norm(
            bm.vert_coord(va_ptr[0]),
            bm.vert_coord(va_ptr[1]),
            bm.vert_coord(va_ptr[2]), N);
        Eigen::Map<const Eigen::Matrix<double,3,1> > va(bm.vert_coord(va_ptr[0]), 3, 1);
        df = -N.dot(va);

        get_face_norm(
            bm.vert_coord(va_ptr[0]),
            bm.vert_coord(va_ptr[1]),
            bm.vert_coord(va_ptr[2]), mpf_N);
        
        Eigen::Matrix<mpf_class,3,1> mpf_va;
        mpf_va << va[0],va[1],va[2];
        mpf_df = -mpf_N.dot(mpf_va);
      } // else-if
      
      int s = dsgn(X.dot(N)+df);
      int mpf_s = sgn(c.first.dot(mpf_N)+c.second*mpf_df)*sgn(c.second);
      if (mpf_s == 0) return 0;
      size_t move_cnt(0);
      while (mpf_s != s && move_cnt < 10) {
        X = X + EPS_MV*my_abs(X)*N/N.norm()*static_cast<double>(mpf_s);
        std::cout << "... move one cut vertex " << std::endl;
        s = dsgn(X.dot(N)+df);
        ++move_cnt;
      }
      {// for check
        s = dsgn(X.dot(N)+df);
        if (mpf_s != s) {
          std::cerr << "-- [ WARNING ] move cv fail, mpf_s != s" << std::endl;
        }
      }
    }
    
    return 0;
  }


  int move_cv(
      const std::vector<std::vector<size_t> >& cv2checkF,
      std::vector<double>& new_CV)
  {
    auto& CV_cvid = seq_.cut_verts_id();

    size_t vn = seq_.tr().vn();
    auto& V = seq_.tr().mesh().verts();
    size_t cvn = seq_.pure_cv_num();

#pragma omp parallel for
    for (size_t vid = 0; vid < cvn; ++vid) {
      if (cv2checkF.size() == 0) {continue;}

      std::pair<Eigen::Matrix<mpf_class,3,1>,mpf_class> c;

      if (vid >= vn) {
        size_t ori_cvid = CV_cvid[vid];
        seq_.get_cut_vert_coord(seq_.cv(ori_cvid), c);
      } else {
        c.first << V[3*vid], V[3*vid+1], V[3*vid+2];
        c.second = 1.0;
      }

      move_cv(cv2checkF[vid], c, &new_CV[3*vid]);
    }//for
    
    return 0;
  }

  // try to move the coordinates of cut-vertices
  // under double floating-point to get embedding results
  int round_and_embed(
      std::vector<double>& new_CV)
  {
    auto& P = seq_.tr().plane();

    // for each cut-vertex, find all adjacent triangles/planes under predefined radius
    std::vector<std::vector<size_t> > cv2nearF;
    std::vector<std::vector<size_t> > cv2passF; // all faces passing each cut-vertex.
    std::vector<size_t> cv2ele;
    get_all_passing_faces(cv2passF, cv2ele);
    
    size_t pure_cvn = new_CV.size()/3;
    size_t pn = seq_.tr().pn();
    size_t tn = seq_.tr().tn();

    cv2nearF.resize(pure_cvn);

#pragma omp parallel for
    for (size_t vi = 0; vi < pure_cvn; ++vi) {
      const Eigen::Map<const Eigen::Matrix<double,3,1> > v(&new_CV[3*vi], 3, 1);
      std::vector<size_t> fs;

      size_t start_f = tn;

      for (size_t fi = start_f; fi < pn+tn; ++fi) {
        Eigen::Matrix<double, 3, 1> N={1.0,0.0,0.0};
        double df(0.0);
        
        if (fi >= tn) {
          size_t pi = fi-tn;
          N << P[4*pi], P[4*pi+1], P[4*pi+2];
          df = P[4*pi+3];
        }
        
        if ( fabs(v.dot(N)+df) <= ((N.norm()*my_abs(v)*EPS_A)) ) {
          fs.push_back(fi);
        }
      }

      std::swap(cv2nearF[vi], fs);
      
    }//for


    {
      std::vector<std::vector<size_t> > cv2checkF(pure_cvn);
      {
#pragma omp parallel for        
        for (size_t vi = 0; vi < pure_cvn; ++vi) {
          std::sort(cv2nearF[vi].begin(), cv2nearF[vi].end());
          cv2nearF[vi].erase(std::unique(cv2nearF[vi].begin(), cv2nearF[vi].end()),
                             cv2nearF[vi].end());
        
          std::sort(cv2passF[vi].begin(), cv2passF[vi].end());
          cv2passF[vi].erase(std::unique(cv2passF[vi].begin(), cv2passF[vi].end()),
                             cv2passF[vi].end());
        
          vector_minus(cv2nearF[vi], cv2passF[vi], cv2checkF[vi]);
        }
      }

      
      CALL_FUNC(move_cv(cv2checkF, new_CV));
    }
    
    return 0;
  }

  int move_by_point_line_order(
      int direc,
      const std::vector<mpf_class>& mpf_cut_verts_nm,
      const std::vector<mpf_class>& mpf_cut_verts_dm,
      std::vector<double>& new_CV,
      size_t va, size_t vb, size_t vc)
  {
    size_t x = abs(direc)%3;
    size_t y = (abs(direc)+1)%3;
    if (direc < 0) std::swap(x, y);

    Eigen::Matrix<double, 2, 1> A, B, C;
    A << new_CV[3*va+x], new_CV[3*va+y];
    B << new_CV[3*vb+x], new_CV[3*vb+y];
    C << new_CV[3*vc+x], new_CV[3*vc+y];

    auto BA = A-B;
    auto BC = C-B;

    int s = vec2d_cross_sgn(
        mpf_cut_verts_nm[3*va+x],
        mpf_cut_verts_nm[3*va+y],
        mpf_cut_verts_dm[va],
        mpf_cut_verts_nm[3*vb+x],
        mpf_cut_verts_nm[3*vb+y],
        mpf_cut_verts_dm[vb],
        mpf_cut_verts_nm[3*vc+x],
        mpf_cut_verts_nm[3*vc+y],
        mpf_cut_verts_dm[vc]);

    double c = BC[0]*BA[1]-BC[1]*BA[0];
    
    // find the smallest path between (A, BC)
    // dx*|BC| = |BC x BA|
    Eigen::Matrix<double, 2, 1> BC_p;
    BC_p[0] = -BC[1]; BC_p[1] = BC[0];

    while (s!=0 && dsgn(c)!=s && fabs(c)<EPS_MV*(my_abs2(A)+my_abs2(B)+my_abs2(C))*BC.norm()) {
      A += (double)s*EPS_MV*(my_abs2(A)+my_abs2(B)+my_abs2(C))*BC_p/BC_p.norm();
      // update
      new_CV[3*va+x] = A[0];
      new_CV[3*va+y] = A[1];
      c = BC[0]*(A[1]-B[1])-BC[1]*(A[0]-B[0]);
    }
    
    return 0;
  }

  int move_bound(
      const std::vector<mpf_class>& mpf_cut_verts_nm,
      const std::vector<mpf_class>& mpf_cut_verts_dm,
      std::vector<double>& new_CV)
  {
    size_t pn = seq_.tr().pn();
    size_t tn = seq_.tr().tn();

    // for each cut-plane
    for (size_t fi = 0; fi < pn; ++fi) {
      // find the boundary edge
      // move the vertex according to the two adjacent vertices
      auto& hes = hfm_.ori_face_he(fi+tn); // for cut-plane fi

      MPF_VEC norm;
      get_face_norm(fi+tn, norm);
      auto direc = choose_axis(norm);

      for (auto& he : hes) {
        if (he.oppo()->ori_f() < tn || he.next()->oppo()->ori_f() < tn) {// only choose boundary
          size_t vc = he.start_v();
          size_t va = he.end_v();
          size_t vb = he.next()->end_v();

          move_by_point_line_order(direc, mpf_cut_verts_nm,
                                   mpf_cut_verts_dm, new_CV, va, vb, vc);
        }
      }
    } //for
    
    return 0;
  }

  int triangulate(
      //const std::vector<size_t>& cv2newid,
      const std::vector<double>& new_CV,
      std::vector<std::vector<size_t> >& face_triangles)
  {
    std::cout << "# Triangulate cut-faces ..." << std::endl;
    size_t pn = seq_.tr().pn();
    size_t tn = seq_.tr().tn();

    //// only need to triangulate pn + tn,
    //// the opposite of cut-faces in cut-plane are not needed to be done.

    face_triangles.resize(pn+tn);
    
#pragma omp parallel for
    for (size_t fi = 0; fi < pn+tn; ++fi) {
      auto& hes = hfm_.ori_face_he(fi);
      hfm_.set_flag(fi, 0);

      std::vector<std::vector<size_t> > polys;
      // get polys
      for (auto& he : hes) {
        if (he.flag() == 0) {
          std::vector<size_t> poly;
          auto st = &he;
          do {
            st->flag() = 1;
            poly.push_back(st->start_v());
            st = st->next();
          } while (st != &he);

          
          if (poly.size() > 2) {
            polys.push_back(poly);
          }
        }
      }

      std::list<size_t> new_tris;

      if (fi < tn) { // for bounary triangles, all cut-faces are convex
        for (auto& poly : polys) {
          std::vector<size_t> tris;
          if (poly.size() > 2) { // cut-faces on triangles, all polys are convex
            tris.reserve(3*(poly.size()-2));
            for (size_t i = 1; i+1 < poly.size(); ++i) {
              tris.push_back(poly[0]);
              tris.push_back(poly[i]);
              tris.push_back(poly[i+1]);
            }
          }
          new_tris.insert(new_tris.end(), tris.begin(), tris.end());
        }//for
        
      } else {
        // get norm
        MPF_VEC norm;
        get_face_norm(fi, norm);
        auto direc = choose_axis(norm);
        FASSERT(direc != 0);

        for (auto& poly : polys) {
          std::vector<size_t> tris;
          if (poly.size() == 3) {
            tris = poly;
          } else {
            triangulate_with_CGAL(direc, new_CV, poly, tris);
          }

          new_tris.insert(new_tris.end(), tris.begin(), tris.end());
        }//polys
        
      }//if
      
      face_triangles[fi].resize(new_tris.size());
      std::copy(new_tris.begin(), new_tris.end(), face_triangles[fi].begin());
    }
    
    std::cout << "# Triangulated." << std::endl;
    return 0;
  }

  int print_poly(
      const std::vector<size_t>& poly,
      const std::unordered_map<size_t, std::pair<double,double> >& p2uv_d)
  {

    std::cout << "-- poly n: " << poly.size() << std::endl;
    for (auto& v : poly) {
      auto it = p2uv_d.find(v);
      std::cout << "v " << it->second.first << "  " << it->second.second << "  1.0" << std::endl;
    }

    return 0;
  }

  int get_nearest_right_intersection(
      size_t vid,
      const std::set<std::pair<size_t,size_t> >& segs,
      const std::unordered_map<size_t, mpf_rational_point2d>& p2uv,
      std::pair<size_t,size_t>& near_p)
  {
    using namespace std;
    auto& p = p2uv.find(vid)->second;

    Eigen::Matrix<mpf_class,2,1> near_x;
    near_x << std::numeric_limits<double>::max(), 1.0;

    int r = 0;

    for (auto& e : segs) {
      auto& pa = p2uv.find(e.first)->second;
      auto& pb = p2uv.find(e.second)->second;
      
      // sgn(p_v - pa_v)
      int s_a=sgn(get<2>(p))*sgn(get<2>(pa))*sgn(get<1>(p)*get<2>(pa)-get<2>(p)*get<1>(pa));
      // sgn(p_v - pb_v)
      int s_b=sgn(get<2>(p))*sgn(get<2>(pb))*sgn(get<1>(p)*get<2>(pb)-get<2>(p)*get<1>(pb));
      
      if (s_a*s_b <=0) { // intersect
        mpf_class temp_ab=get<2>(pa)*get<1>(pb) - get<1>(pa)*get<2>(pb);
        int sab = sgn(temp_ab);
        int s = sab*sgn(get<2>(pa))*sgn(get<2>(pb));

        if (s==0) {continue;}

        // get intersection point
        Eigen::Matrix<mpf_class,2,1> x;
        if (sab != 0) {
          x << (get<0>(pb)*get<2>(pa)-get<2>(pb)*get<0>(pa))*
              (get<1>(p)*get<2>(pa)-get<2>(p)*get<1>(pa)) +
              get<0>(pa)*get<2>(p)*(temp_ab),
              get<2>(p)*get<2>(pa)*(temp_ab);

        } else {
          x << get<0>(pa), get<2>(pa);
        }

        // check intersection
        if (sgn(x[0]*get<2>(p)-get<0>(p)*x[1])*sgn(get<2>(p))*sgn(x[1]) > 0) {
          int s_xx = sgn(x[1])*sgn(near_x[1])*sgn(x[0]*near_x[1]-near_x[0]*x[1]);
          if (s_xx < 0 || (s_xx==0 && s>0)) {
            std::swap(near_x, x);
            r = s;
            near_p = e;
          }
        }//if
      }
    }
    
    return r;
  }
  

  int fast_find_parent_positive_loop(
      const std::vector<std::vector<size_t> >& polys,
      const std::vector<size_t>& neg_loops,
      const std::list<size_t>& pos_loops,
      const std::unordered_map<size_t, mpf_rational_point2d>& p2uv,
      const std::unordered_map<size_t, std::pair<double,double> >& p2uv_d,
      std::vector<size_t>& nl_ct)
  {
    segment_space pos_seg_space_y(40);
    segment_space neg_seg_space_y(40);

    build_segment_space<1,std::list<size_t> >(polys, pos_loops, p2uv_d, pos_seg_space_y);
    build_segment_space<1,std::vector<size_t> >(polys, neg_loops, p2uv_d, neg_seg_space_y);
    nl_ct.resize(neg_loops.size(), INVALID_NUM);
    
    std::map<std::pair<size_t,size_t>, size_t> e2p;
    for (size_t pi = 0; pi < polys.size(); ++pi) {
      auto& po = polys[pi];
      if (po.size() == 1) { continue; }
      size_t pre_id = po.back();
      for (auto v : po) {
        e2p.insert(std::make_pair(std::make_pair(pre_id, v), pi));
        std::swap(pre_id, v);
      }
    }//for

#pragma omp parallel for
    for (size_t pi = 0; pi < neg_loops.size(); ++pi) {
          
      size_t nl_id = neg_loops[pi];
      size_t vid = polys[nl_id].front();
      FASSERT(p2uv_d.find(vid)!=p2uv_d.end());
      auto& start_uv_d = p2uv_d.find(vid)->second;

      {
        std::set<std::pair<size_t,size_t> > pos_segs;
        pos_seg_space_y.get_segs(start_uv_d.second, start_uv_d.second, pos_segs);
        std::pair<size_t,size_t> near_pos_R;
        // >0, find the near right point, direction is also right
        // =0, no intersection in the right
        // <0, the right near intersection has wrong direction.
        if (get_nearest_right_intersection(vid, pos_segs, p2uv, near_pos_R)>0) {
          auto eit = e2p.find(near_pos_R);
          nl_ct[pi] = eit->second;
        }
      }
      if (polys[nl_id].size() == 1 && nl_ct[pi]!=INVALID_NUM) {
        std::set<std::pair<size_t,size_t> > neg_segs;
        neg_seg_space_y.get_segs(start_uv_d.second, start_uv_d.second, neg_segs);
        
        std::pair<size_t,size_t> near_neg_R;        
        if (get_nearest_right_intersection(vid, neg_segs, p2uv, near_neg_R)<0) {
          nl_ct[pi] = INVALID_NUM;
        }
      }
      
    }//for
    
    return 0;
  }


  int find_parent_positive_loop(
      const std::vector<std::vector<size_t> >& polys,
      const std::vector<size_t>& neg_loops,
      const std::list<size_t>& pos_loops,
      const std::unordered_map<size_t, mpf_rational_point2d>& p2uv,
      std::vector<size_t>& nl_ct)
  {
    std::vector<poly_bound_box> poly_bb(polys.size());
    for (size_t i = 0; i < polys.size(); ++i) {
      auto& bb = poly_bb[i];
      bb.init();
      for (auto& v : polys[i]) {
        auto it = p2uv.find(v);
        bb.set_max(it->second); bb.set_min(it->second);
      }

    }

    std::unordered_set<size_t> neg_set(neg_loops.begin(), neg_loops.end());
    // judge whether there exists
    // only one smallest bounding box to contain each negative loop
#pragma omp parallel for
    for (size_t i = 0; i < neg_loops.size(); ++i) {
      size_t nl_id = neg_loops[i];
      std::vector<size_t> ps;

      //// a robust one
      if (polys[nl_id].size() > 1) {
        for (auto pj : pos_loops) {
          if (poly_bb[pj].is_contain(poly_bb[nl_id])) {
            ps.push_back(pj);
          }
        }
      } else {
        for (size_t pj = 0; pj < polys.size(); ++pj)
          if (poly_bb[pj].is_contain(poly_bb[nl_id])) {
            ps.push_back(pj);
          }
      }
      
      if (ps.size() > 1) {
        { // find the loop containing one vertex of the negative loop
          std::vector<size_t> temp_ct;
          for (auto& pid : ps) {
            // consider negative loop for single point polygon
            if (fabs(calculate_theta_sum(polys[pid], polys[nl_id].front(), p2uv)) > kPI) {
              temp_ct.push_back(pid);
            }
          }

          std::swap(ps, temp_ct);
        }
        if (ps.size() > 1) { // only keep the smallest one
          for (size_t j = 0; j < ps.size(); ++j) {
            size_t pa = ps[j];
            size_t pa_cnt = 0;
            for (size_t k = 0; k < ps.size(); ++k) {
              if (k == j) continue;
              size_t pb = ps[k];
              if (poly_bb[pb].is_contain(poly_bb[pa])) { ++pa_cnt; }
            }
            if (pa_cnt == ps.size()-1) {

              if (neg_set.count(pa) == 0) {
                nl_ct[i] = pa;
              }
              break;
            }
          }
        } else if (ps.size()==1) {
          nl_ct[i] = ps.front();
        }
      } else if (ps.size() == 1) {
        if (polys[nl_id].size() > 1) {
          nl_ct[i] = ps.front();
        } else {
          if (calculate_theta_sum(polys[ps.front()], polys[nl_id].front(), p2uv) > kPI) {
            nl_ct[i] = ps.front();
          }
        }
      }
      if (polys[nl_id].size()>1) {
        FASSERT(nl_ct[i] != INVALID_NUM);
      }

    }
    return 0;
  }
  

  // connect isolated loops
  int connect_cut_faces(
      const std::vector<mpf_class>& mpf_cut_verts_nm,
      const std::vector<mpf_class>& mpf_cut_verts_dm)
  {
    using namespace std;
    size_t pn = seq_.tr().pn();
    size_t tn = seq_.tr().tn();


    std::vector<std::vector<size_t> > p2v_vec;
    p2v_vec.resize(pn);

    {
      auto& M = seq_.tr().mesh().triangles();
      std::unordered_set<size_t> used_vs(M.begin(), M.end());
      
#pragma omp parallel for
      for (size_t pi = 0; pi < pn; ++pi) {
        std::list<size_t> v_list;
        for (auto vi : used_vs) {
          if (seq_.tr().pv(pi,vi)==0) {
            //v_list.push_back(vi);
            auto& fs = seq_.v2f(vi);
            bool is_add=true;
            for (auto& f : fs)  {
              if (seq_.tr().pt(pi, f)>1) {
                is_add = false; break;
              }
            }
            if (is_add) { v_list.push_back(vi); }
          }
        }
        p2v_vec[pi].resize(v_list.size());
        std::copy(v_list.begin(), v_list.end(), p2v_vec[pi].begin());
      }
    }

    const auto& V = seq_.tr().mesh().verts();

    // only need to deal with cut-faces in cut planes
#pragma omp parallel for
    for (size_t fi = 0; fi < pn; ++fi) {
      auto& hes = hfm_.ori_face_he(fi+tn);
      hfm_.set_flag(fi+tn, 0);


      std::vector<std::vector<size_t> > polys;
      std::vector<half_edge::HE_PTR> polys_he;
      // get polys
      for (auto& he : hes) {
        if (he.flag() == 0) {
          std::vector<size_t> poly;
          auto st = &he;
          do {
            st->flag() = 1;
            poly.push_back(st->start_v());
            st = st->next();
          } while (st != &he);
          polys.push_back(poly);
          polys_he.push_back(&he);
        }
      }

      MPF_VEC norm;
      get_face_norm(fi+tn, norm);
      auto direc = choose_axis(norm);
      FASSERT(direc != 0);


      // <nominator(u,v), denominator>
      std::unordered_map<size_t, mpf_rational_point2d> p2uv;
      std::unordered_map<size_t, std::pair<double,double> > p2uv_d;
      {
        size_t id = abs(direc)-1;

        for (size_t pi = 0; pi < polys.size(); ++pi) {
          auto& po = polys[pi];
          //for (auto& po : polys) {
          for (auto& vid : po) {
            if (p2uv.count(vid) == 0) {
              FASSERT(3*vid+2< mpf_cut_verts_nm.size() && vid < mpf_cut_verts_dm.size());
              mpf_rational_point2d
                               uv(mpf_cut_verts_nm[3*vid+(id+1)%3],
                                  mpf_cut_verts_nm[3*vid+(id+2)%3], mpf_cut_verts_dm[vid]);


              std::pair<double, double> uv_d;
              uv_d.first = std::get<0>(uv).get_d()/std::get<2>(uv).get_d();
              uv_d.second = std::get<1>(uv).get_d()/std::get<2>(uv).get_d();
              
              if (direc < 0) {
                std::swap(get<0>(uv), get<1>(uv));
                std::swap(uv_d.first, uv_d.second);
              }
              p2uv.insert(std::make_pair(vid, uv));
              p2uv_d.insert(std::make_pair(vid, uv_d));
            }
          }
        }

        for (auto& vv : p2v_vec[fi]) {
          if (p2uv.count(vv)==0) {
            std::vector<size_t> vs(1, vv);
            polys.push_back(vs);
            polys_he.push_back(nullptr);
            mpf_rational_point2d uv(
                mpf_cut_verts_nm[3*vv+(id+1)%3],
                mpf_cut_verts_nm[3*vv+(id+2)%3], mpf_cut_verts_dm[vv]);
            p2uv.insert(std::make_pair(vv, uv));

            std::pair<double,double> uv_d(V[3*vv+(id+1)%3], V[3*vv+(id+2)%3]);
            p2uv_d.insert(std::make_pair(vv, uv_d));
          }
        }
      }


      std::vector<size_t> conn_edges, conn_polys;

      connect_cut_faces(polys, p2uv, p2uv_d, conn_edges, conn_polys);


      // update half-face structure
      for (size_t i=0; i < conn_edges.size(); i+=2) {
        auto pa_he = polys_he[conn_polys[i]];
        auto pb_he = polys_he[conn_polys[i+1]];
        size_t va = conn_edges[i];
        size_t vb = conn_edges[i+1];
        std::vector<half_edge::HE_PTR> h1, h2;

        // a
        if (pa_he != nullptr) {
          auto st = pa_he;          
          do {
            if (st->end_v() == va) {
              h1.push_back(st);
            }
            st = st->next();
          } while (st != pa_he);
        }

        // b
        if (pb_he != nullptr) {
          auto st = pb_he; 
          do {
            if (st->end_v() == vb) {
              h2.push_back(st);
            }
            st = st->next();
          } while (st != pb_he);
        }

        
        if (pa_he!=nullptr) {FASSERT(h1.size()>0);}
        if (pb_he!=nullptr) {FASSERT(h2.size()>0);}

        half_edge::HE_PTR ha(nullptr), hb(nullptr), hc(nullptr), hd(nullptr);
        if (h1.size()>0) {
          ha = h1.front();
          hb = ha->next();
        }
        if (h2.size()>0) {
          hc = h2.front();
          hd = hc->next();
        }

        if (h1.size()>1 || h2.size()>1) {
          auto ita = p2uv.find(va);
          auto itb = p2uv.find(vb);
          auto& pA = ita->second;
          auto& pB = itb->second;

          if (h1.size() > 1) {
            bool is_find=false;
            for (auto h : h1) {
              auto itc = p2uv.find(h->start_v());
              auto itd = p2uv.find(h->next()->end_v());
              FASSERT(itc!=p2uv.end() && itd!=p2uv.end());
              auto& pC = itc->second;
              auto& pD = itd->second;
              short bad = tri_area_sign(pB, pA, pD);
              short cab = tri_area_sign(pC, pA, pB);
              short cad = tri_area_sign(pC, pA, pD);


              if ((cad>0 && bad>0 && cab>0) || (cad<0 && (bad>0 || cab>0))
                  || (cad==0 && itc->first!=itd->first && bad>0 && cab>0)
                  || (cad==0 && itc->first==itd->first)) {
                is_find = true;
                ha = h; hb = h->next(); break;
              }
            }

            FASSERT(is_find);
          }// if h1

          if (h2.size() > 1) {
            bool is_find = false;
            for (auto h : h2) {
              auto itc = p2uv.find(h->start_v());
              auto itd = p2uv.find(h->next()->end_v());
              FASSERT(itc!=p2uv.end() && itd!=p2uv.end());
              auto& pC = itc->second;
              auto& pD = itd->second;
              short abd = tri_area_sign(pA, pB, pD);
              short cba = tri_area_sign(pC, pB, pA);
              short cbd = tri_area_sign(pC, pB, pD);

              
              if ((cbd>0 && abd>0 && cba>0) || (cbd<0 && (abd>0 || cba>0))
                  || (cbd==0 && itc->first!=itd->first && abd>0 && cba>0)
                  || (cbd==0 && itc->first==itd->first)) {
                is_find = true;
                hc = h; hd = h->next(); break;
              }
            }//for h2

            FASSERT(is_find);
          }//if h2

        }// if h1 or h2

        // update half-edge
        auto h1_new = hfm_.add_he(va, fi+tn, INVALID_NUM);
        auto h2_new = hfm_.add_he(vb, fi+tn, INVALID_NUM);

        if (ha!=nullptr && hb!=nullptr) {
          hfm_.set_next_prev(ha, h1_new);
          hfm_.set_next_prev(h2_new, hb);
        } else {
          hfm_.set_next_prev(h2_new,h1_new);
        }

        if (hc!=nullptr && hd!=nullptr) {
          hfm_.set_next_prev(h1_new, hd);
          hfm_.set_next_prev(hc, h2_new);
        } else {
          hfm_.set_next_prev(h1_new, h2_new);
        }
        
        hfm_.set_oppo(h1_new, h2_new);        

        auto h1_fop_new = hfm_.add_he(vb, fi+tn+pn, INVALID_NUM);
        auto h2_fop_new = hfm_.add_he(va, fi+tn+pn, INVALID_NUM);

        hfm_.set_foppo(h1_new, h1_fop_new);
        hfm_.set_foppo(h2_new, h2_fop_new);

        if (ha!=nullptr && hb!=nullptr) {
          hfm_.set_next_prev(h1_fop_new, ha->foppo());
          hfm_.set_next_prev(hb->foppo(), h2_fop_new);
        } else {
          hfm_.set_next_prev(h1_fop_new, h2_fop_new);
        }

        if (hc!=nullptr && hd!=nullptr) {
          hfm_.set_next_prev(hd->foppo(), h1_fop_new);
          hfm_.set_next_prev(h2_fop_new, hc->foppo());
        } else {
          hfm_.set_next_prev(h2_fop_new, h1_fop_new);
        }
        
        hfm_.set_oppo(h1_fop_new, h2_fop_new);
      }
    }
    
    return 0;
  }

  double calculate_theta_sum(
      const std::vector<size_t>& poly,
      const std::unordered_map<size_t, std::tuple<mpf_class,mpf_class,mpf_class> >& p2uv)
  {
    using namespace std;

    if (poly.size() > 2) {
      double sum(0.0);
      for (size_t i = 0; i < poly.size(); ++i) {
        size_t pre_id = poly[(i+poly.size()-1)%poly.size()];
        size_t cur_id = poly[i];
        size_t nex_id = poly[(i+1)%poly.size()];

        auto it1 = p2uv.find(pre_id);
        const auto& a = it1->second;
        auto it2 = p2uv.find(cur_id);
        const auto& b = it2->second;
        auto it3 = p2uv.find(nex_id);
        const auto& c = it3->second;

        if (a!=c) {
          mpf_class y = (get<0>(b)*get<2>(a)-get<0>(a)*get<2>(b))*
              (get<1>(c)*get<2>(b)-get<1>(b)*get<2>(c))-
              (get<1>(b)*get<2>(a)-get<1>(a)*get<2>(b))*
              (get<0>(c)*get<2>(b)-get<0>(b)*get<2>(c));
          mpf_class x = (get<0>(b)*get<2>(a)-get<0>(a)*get<2>(b))*
              (get<0>(c)*get<2>(b)-get<0>(b)*get<2>(c))+          
              (get<1>(c)*get<2>(b)-get<1>(b)*get<2>(c))*
              (get<1>(b)*get<2>(a)-get<1>(a)*get<2>(b));

          short s = sgn(get<2>(a))*sgn(get<2>(c));
          y *= s;
          x *= s;
      
          mpf_class xy = sgn(y)*y+sgn(x)*x;
          double dy = y.get_d()/xy.get_d();
          double dx = x.get_d()/xy.get_d();

          double ang = atan2(dy, dx);
          sum += ang;
        } else {
          sum += -kPI;
        }
      }
      return sum;    
    }

    return -2.0*kPI;
  }

  class poly_bound_box {
   public:
    mpf_class u1[2], v1[2], u2[2], v2[2];
    void init()
    {
      u1[0] = std::numeric_limits<double>::max();
      u1[1] = 1.0;
      v1[0] = std::numeric_limits<double>::max();
      v1[1] = 1.0;
      
      u2[0] = -std::numeric_limits<double>::max();
      u2[1] = 1.0;
      v2[0] = -std::numeric_limits<double>::max();
      v2[1] = 1.0;
    }

    void print()
    {
      std::cout << "-- bbox: min p  " <<  u1[0].get_d()/u1[1].get_d() << ", "
                << v1[0].get_d()/v1[1].get_d() << std::endl;
      std::cout << "   bbox: max p  " <<  u2[0].get_d()/u2[1].get_d() << ", "
                << v2[0].get_d()/v2[1].get_d() << std::endl;
    }
    
    void set_max(const std::tuple<mpf_class, mpf_class, mpf_class>& p)
    {
      using namespace std;
      {
        int s = sgn(get<0>(p)*u2[1]-get<2>(p)*u2[0])*sgn(u2[1])*sgn(get<2>(p));
        if (s > 0) {
          u2[0] = get<0>(p);
          u2[1] = get<2>(p);
        }
      }

      {
        int s = sgn(get<1>(p)*v2[1]-get<2>(p)*v2[0])*sgn(v2[1])*sgn(get<2>(p));
        if (s > 0) {
          v2[0] = get<1>(p);
          v2[1] = get<2>(p);
        }
      }
    }

    void set_min(const std::tuple<mpf_class, mpf_class, mpf_class>& p)
    {
      {
        int s = sgn(std::get<0>(p)*u1[1]-std::get<2>(p)*u1[0])*sgn(u1[1])*sgn(std::get<2>(p));
        if (s < 0) {
          u1[0] = std::get<0>(p);
          u1[1] = std::get<2>(p);
        }
      }

      {
        int s = sgn(std::get<1>(p)*v1[1]-std::get<2>(p)*v1[0])*sgn(v1[1])*sgn(std::get<2>(p));
        if (s < 0) {
          v1[0] = std::get<1>(p);
          v1[1] = std::get<2>(p);
        }
      }
    }

    // p>q
    bool is_big(const mpf_class* p, const mpf_class* q)
    {
      return (sgn(p[0]*q[1]-q[0]*p[1])*sgn(p[1])*sgn(q[1])>0);
    }

    bool is_contain(const poly_bound_box& b)
    {
      auto& a = *this;
      // check a > b
      if (is_big(b.u1, a.u1) && is_big(a.u2, b.u2)
          && is_big(b.v1, a.v1) && is_big(a.v2, b.v2)) {
        return true;
      }
      
      return false;
    }

    bool is_intersect(const poly_bound_box& b)
    {
      if ( is_big(b.u1, u2) || is_big(u1,b.u2)
          || is_big(v1,b.v2) || is_big(b.v1,v2) ) return false;
      return true;
    }

  };

  double calculate_theta_sum(
      const std::vector<size_t>& poly, size_t vid,
      const std::unordered_map<size_t, std::tuple<mpf_class,mpf_class,mpf_class> >& p2uv)
  {
    double sum(0.0);

    auto mid_it = p2uv.find(vid);
    const auto& p = mid_it->second;

    using namespace std;
    
    for (size_t i = 0; i < poly.size(); ++i) {
      size_t a_id = poly[i];
      size_t b_id = poly[(i+1)%poly.size()];

      auto it1 = p2uv.find(a_id);
      const auto& a = it1->second;
      auto it2 = p2uv.find(b_id);
      const auto& b = it2->second;

      mpf_class sv = (get<0>(a)*get<2>(p)-get<2>(a)*get<0>(p))*
          (get<1>(b)*get<2>(p)-get<2>(b)*get<1>(p))
          - (get<1>(a)*get<2>(p)-get<2>(a)*get<1>(p))*
          (get<0>(b)*get<2>(p)-get<2>(b)*get<0>(p)); // sin
      
      mpf_class cv = (get<0>(a)*get<2>(p)-get<2>(a)*get<0>(p))*
          (get<0>(b)*get<2>(p)-get<2>(b)*get<0>(p))
          + (get<1>(a)*get<2>(p)-get<2>(a)*get<1>(p))*
          (get<1>(b)*get<2>(p)-get<2>(b)*get<1>(p)); // cos

      short s = sgn(get<2>(a))*sgn(get<2>(b));

      sv *= s;
      cv *= s;

      mpf_class xy = sgn(sv)*sv+sgn(cv)*cv;

      double dy = sv.get_d() / xy.get_d();
      double dx = cv.get_d() / xy.get_d();
      
      sum += atan2(dy, dx); // use angle to avoid small area
    }
    
    return sum;
  }

  short tri_area_sign(
      const mpf_rational_point2d& a,
      const mpf_rational_point2d& b,
      const mpf_rational_point2d& c) {
    using namespace std;
    mpf_class area = (get<2>(a)*get<0>(b)-get<2>(b)*get<0>(a))*(get<2>(a)*get<1>(c)-get<2>(c)*get<1>(a))
        -(get<2>(a)*get<1>(b)-get<2>(b)*get<1>(a))*(get<2>(a)*get<0>(c)-get<2>(c)*get<0>(a));
    return sgn(area)*sgn(get<2>(b))*sgn(get<2>(c));
  }

  short tri_area_sign(
      const std::pair<double, double>& A,
      const std::pair<double, double>& B,
      const std::pair<double, double>& C)
  {
    // (B-A)x(C-A)
    return dsgn((B.first-A.first)*(C.second-A.second)
                -(B.second-A.second)*(C.first-A.first));
  }

  // a, edge(b,c)
  short line_partial_order(
      const mpf_rational_point2d& a,
      const mpf_rational_point2d& b,
      const mpf_rational_point2d& c)
  {
    using namespace std;
    mpf_class po =
        (get<2>(a)*get<0>(b)-get<2>(b)*get<0>(a))*(get<2>(a)*get<0>(c)-get<2>(c)*get<0>(a))
        +(get<2>(a)*get<1>(b)-get<2>(b)*get<1>(a))*(get<2>(a)*get<1>(c)-get<2>(c)*get<1>(a));
    return sgn(po)*sgn(get<2>(b))*sgn(get<2>(c));
  }

  bool is_pass_vert(size_t va, size_t vb, size_t vc,
      const std::unordered_map<size_t, mpf_rational_point2d>& p2uv)
  {
    auto ita = p2uv.find(va);
    auto itb = p2uv.find(vb);
    auto itc = p2uv.find(vc);
    FASSERT(va!=vb);
    FASSERT(ita != p2uv.end() && itb!=p2uv.end() && itc!=p2uv.end());
    
    auto& ta = ita->second;
    auto& tb = itb->second;
    auto& tc = itc->second;

    short cab = line_partial_order(tc, ta, tb);
    if (cab < 0) {return true;}
    return false;
  }

  bool is_edge_intersect(
      size_t va, size_t vb, size_t vc, size_t vd,
      const std::unordered_map<size_t, mpf_rational_point2d>& p2uv)
  {
    auto ita = p2uv.find(va);
    auto itb = p2uv.find(vb);
    auto itc = p2uv.find(vc);
    auto itd = p2uv.find(vd);
    //FASSERT(va!=vb && vc!=vd);
    FASSERT(ita != p2uv.end() && itb!=p2uv.end() && itc!=p2uv.end() && itd!=p2uv.end());
    
    auto& ta = ita->second;
    auto& tb = itb->second;
    auto& tc = itc->second;
    auto& td = itd->second;

    short abc = tri_area_sign(ta, tb, tc);
    short abd = tri_area_sign(ta, tb, td);
    short cda = tri_area_sign(tc, td, ta);
    short cdb = tri_area_sign(tc, td, tb);

    if (abc*abd<0 && cda*cdb<0) {return true;}
    
    if (abc*abd>0 || cda*cdb>0) {return false;}
    
    if (abc==0 && abd!=0) {
      short cab = line_partial_order(tc, ta, tb);
      if (cab < 0) {return true;}
      return false;
    }
    
    if (abd==0 && abc!=0) {
      short dab = line_partial_order(td, ta, tb);
      if (dab < 0) {return true;}
      return false;
    }
    
    if (cda==0 && cdb!=0) {
      short acd = line_partial_order(ta, tc, td);
      if (acd < 0) {return true;}
      return false;
    }
    
    if (cdb==0 && cda!=0) {
      short bcd = line_partial_order(tb, tc, td);
      if (bcd < 0) {return true;}
      return false;
    }
    
    if (abc==0 && abd==0) {
      short acd = line_partial_order(ta, tc, td);
      short bcd = line_partial_order(tb, tc, td);
      short cab = line_partial_order(tc, ta, tb);
      short dab = line_partial_order(td, ta, tb);
      if (acd<0 || bcd<0 || cab<0 || dab<0) {return true;}
    }
    return false;
  }


  int connect_cut_faces(
      const std::vector<std::vector<size_t> >& polys,
      const std::unordered_map<size_t, std::tuple<mpf_class,mpf_class,mpf_class> >& p2uv,
      const std::unordered_map<size_t, std::pair<double,double> >& p2uv_d,  
      std::vector<size_t>& conn_edges, std::vector<size_t>& conn_polys)
  {
    const auto& CV_type = seq_.cut_verts_type();
    std::vector<size_t> neg_loops;
    std::list<size_t> pos_loops;
    for (size_t i = 0; i < polys.size(); ++i) {
      if (polys[i].size()>3) {
        bool is_neg = true;
        for (auto& v : polys[i]) {
          if (CV_type[v]=='g') {
            is_neg = false;
            break;
          }
        }
        if (!is_neg) {
          pos_loops.push_back(i);
          continue;
        }
      }
      double theta = calculate_theta_sum(polys[i], p2uv);
      if (theta < -kPI) {
        neg_loops.push_back(i);
      } else {
        pos_loops.push_back(i);
      }

    }

    if (neg_loops.size() == 0 || pos_loops.size()==0) return 0;

    std::vector<size_t> nl_ct(neg_loops.size(), INVALID_NUM);

    // find parent positive loop
    find_parent_positive_loop(polys, neg_loops, pos_loops, p2uv, nl_ct);

    std::unordered_map<size_t, std::vector<size_t> > p_loop_to_neg_loop;

    // group negative loops under positive loops
    // for positive loops containing negative loops
    for (size_t i = 0; i < nl_ct.size(); ++i) {
      if (nl_ct[i] == INVALID_NUM) {continue;}
      
      auto it = p_loop_to_neg_loop.find(nl_ct[i]);
      if (it == p_loop_to_neg_loop.end()) {
        std::vector<size_t> nls;
        nls.reserve(2);
        nls.push_back(nl_ct[i]); // positive loop
        nls.push_back(neg_loops[i]);
        p_loop_to_neg_loop.insert(std::make_pair(nl_ct[i], nls));
      } else {
        it->second.push_back(neg_loops[i]);
      }
    }

    union_find_set ufs(polys.size());

    // connect negative to other loops
    for (auto& pl : p_loop_to_neg_loop) {
      const auto& ps = pl.second; // one positive loop + all contained negative loops

      // <P, (fid, vid)>
      std::vector<std::map<std::pair<double,double>, size_t> > p2id(ps.size());
      //std::shared_ptr<Seg_Tree> tree;
      std::vector<std::shared_ptr<Point_Tree> > po_trees;
      build_trees(polys, ps, p2uv_d, po_trees, p2id);

      segment_space seg_space_x(20), seg_space_y(20);
      build_segment_space(polys, ps, p2uv_d, seg_space_x, seg_space_y);
      
      size_t conn_cnt(0);

      // build near point-pair for each loops
      for (size_t loop_i = 0; loop_i < ps.size() && conn_cnt+1<ps.size(); ++loop_i) {
        for (size_t loop_j = loop_i+1; loop_j < ps.size() && conn_cnt+1<ps.size(); ++loop_j) {
          size_t pi(ps[loop_i]), pj(ps[loop_j]);

          if (ufs.find(pi) == ufs.find(pj)) {continue;}
          
          size_t tree_i(loop_i), tree_j(loop_j);
          if (polys[pi].size() > polys[pj].size()) {
            std::swap(pi, pj);
            std::swap(tree_i, tree_j);
          }
            
          for (auto& pit : p2id[tree_i]) {
            Point pA(pit.first.first, pit.first.second);
            Point pB;
            {
              Neighbor_search search_B(*po_trees[tree_j], pA, 1);
              pB = search_B.begin()->first;
            }

            auto it_pB = p2id[tree_j].find(std::make_pair(pB.x(),pB.y()));
            FASSERT(it_pB != p2id[tree_j].end());
            size_t va_id = pit.second;
            size_t vb_id = it_pB->second;

            std::set<std::pair<size_t,size_t> > segs_X, segs_Y;
            seg_space_x.get_segs(pA.x(), pB.x(), segs_X);
            seg_space_y.get_segs(pA.y(), pB.y(), segs_Y);

            bool is_ok = true;
            for (auto& e : segs_X) {
              if (segs_Y.find(e) != segs_Y.end()) { // common edge
                if (is_edge_intersect(va_id, vb_id, e.first, e.second, p2uv)) {
                  is_ok = false;
                  break;
                }
              }
            }

            if (is_ok) {//check intersection between previous connected edges.
              for (size_t k = 0; k < conn_edges.size(); k+=2) {
                size_t e_va = conn_edges[k];
                size_t e_vb = conn_edges[k+1];
                if (is_edge_intersect(va_id, vb_id, e_va, e_vb, p2uv)) {
                  is_ok = false;
                  break;
                }// if
              } // for k
            }
            
            if (is_ok) {
              conn_edges.push_back(va_id);
              conn_edges.push_back(vb_id);
              conn_polys.push_back(pi);
              conn_polys.push_back(pj);
              ufs.set_union(pi, pj);
              ++conn_cnt;
              break;
            }
          }
        }//for
      }//for

      if (conn_cnt+1 == ps.size()) {continue;}

      
      for (size_t i = 0; i < ps.size() && conn_cnt+1 < ps.size(); ++i) {
        size_t pa = ps[i];
        for (size_t j = i+1; j < ps.size() && conn_cnt+1<ps.size(); ++j) {
          size_t pb = ps[j];
          if (ufs.find(pa) == ufs.find(pb)) {continue;}

          // given pa, pb
          bool is_conn = false;
          
          for (auto va : polys[pa]) {
            if (is_conn) break;
            for (auto vb : polys[pb]) {

              double pA[2], pB[2];
              {
                auto ita = p2uv.find(va);
                auto itb = p2uv.find(vb);
                auto& pA_uv = ita->second;
                auto& pB_uv = itb->second;
                pA[0] = std::get<0>(pA_uv).get_d()/std::get<2>(pA_uv).get_d();
                pA[1] = std::get<1>(pA_uv).get_d()/std::get<2>(pA_uv).get_d();

                pB[0] = std::get<0>(pB_uv).get_d()/std::get<2>(pB_uv).get_d();
                pB[1] = std::get<1>(pB_uv).get_d()/std::get<2>(pB_uv).get_d();
              }

              std::set<std::pair<size_t,size_t> > segs_X, segs_Y;
              seg_space_x.get_segs(pA[0], pB[0], segs_X);
              seg_space_y.get_segs(pA[1], pB[1], segs_Y);

              bool is_ok = true;
              for (auto& e : segs_X) {
                if (segs_Y.find(e) != segs_Y.end()) { // common edge
                  if (is_edge_intersect(va, vb, e.first, e.second, p2uv)) {
                    is_ok = false;
                    break;
                  }
                }
              }

              if (is_ok) {//check intersection between previous connected edges.
                for (size_t k = 0; k < conn_edges.size(); k+=2) {
                  size_t e_va = conn_edges[k];
                  size_t e_vb = conn_edges[k+1];
                  if (is_edge_intersect(va, vb, e_va, e_vb, p2uv)) {
                    is_ok = false;
                    break;
                  }// if
                } // for k
              }
              
              if (is_ok) {
                conn_edges.push_back(va);
                conn_edges.push_back(vb);
                conn_polys.push_back(pa);
                conn_polys.push_back(pb);
                is_conn = true;
                ufs.set_union(pa, pb);

                break;
              }
              
            }// for vb
          }// for va

          if (is_conn) { ++conn_cnt;}
        }// for j
      }// for i

      FASSERT(conn_cnt == ps.size()-1);
    }
    
    return 0;
  }

  const sequence_type& seq() const { return seq_; }
  

 protected:
  const sequence_type& seq_;

  // half face mesh
  hf_mesh hfm_;

  std::vector<std::unordered_set<size_t> > sorted_cv_;
  std::unordered_set<size_t> sorted_ce_;//used to count the cut-edges of cycle-sorting
};



bool is_contain(
    const std::unordered_set<size_t>& fa_set,
    const std::vector<size_t>& fb)
{
  for (auto& v : fb) {
    if (fa_set.find(v) == fa_set.end()) return false;
  }
  return true;
}


int generate_cut_cells(
    const pm_map& para,
    const sequence_type& seq,
    const std::string& out_file,
    std::vector<size_t>& CC, size_t& cell_num,
    std::vector<double>& new_CV,
    std::vector<std::vector<size_t> >& face_triangles)
{
  cut_cell cc(seq);
  CALL_FUNC(cc.build_faces());

  CALL_FUNC(cc.show_faces_num());

  CALL_FUNC(cc.build_cells());

  {
    //calculate CV before connecting cut-faces    
    size_t vn = seq.tr().vn();
    size_t cvn = seq.pure_cv_num();
    std::vector<mpf_class> mpf_cut_verts_nm;
    std::vector<mpf_class> mpf_cut_verts_dm;
    mpf_cut_verts_nm.resize(3*cvn);
    mpf_cut_verts_dm.resize(cvn, 1.0);

    new_CV.resize(3*cvn);
    
    const auto& V = seq.tr().mesh().verts();
    std::copy(V.begin(), V.end(), new_CV.begin());

    // initialize cut-verts
#pragma omp parallel for
    for (size_t i = 0; i < V.size(); ++i) {
      mpf_cut_verts_nm[i] = V[i];
    }

    {
      const auto& CV_id = seq.cut_verts_id();
      FASSERT(3*vn == V.size());

#pragma omp parallel for
      for (size_t vi=vn; vi < cvn; ++vi) {
        size_t cv_id = CV_id[vi];
        std::pair<Eigen::Matrix<mpf_class, 3, 1>, mpf_class> mpf_c;
        seq.get_cut_vert_coord(seq.cv(cv_id), mpf_c);

        for (size_t j = 0; j < 3; ++j) {
          mpf_cut_verts_nm[3*vi+j] = mpf_c.first[j];
          new_CV[3*vi+j] = mpf_c.first[j].get_d()/mpf_c.second.get_d();
        }
        mpf_cut_verts_dm[vi] = mpf_c.second;
      }
    }
    
    const size_t is_connect_hole = para["is-connect-hole"].as<size_t>();
    if (is_connect_hole == 1) {
      CALL_FUNC(cc.connect_cut_faces(
          mpf_cut_verts_nm, mpf_cut_verts_dm));
    }

    // try to keep vert-vert order
    CALL_FUNC(seq.move_cut_verts(new_CV));
    CALL_FUNC(cc.move_bound(
        mpf_cut_verts_nm, mpf_cut_verts_dm, new_CV));

    const size_t is_round_and_embed = para["is-round-and-embed"].as<size_t>();
    if (is_round_and_embed == 1) {
      CALL_FUNC(cc.round_and_embed(new_CV));
    }

  }

  CALL_FUNC(cc.get_cells(CC, cell_num));

  return 0;
}


} // grid
} // fxz
