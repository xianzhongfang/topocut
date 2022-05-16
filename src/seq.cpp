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


#include <queue>
#include <set>

#include "seq.h"
#include "union_find_set.h"

namespace fxz {
namespace grid {

typedef Eigen::Matrix<mpf_class,3,1> MPF_VEC;

// key: (k0,k1), k0/k1=key.
typedef std::pair<double, double> key_type;
typedef std::pair<key_type,size_t> key_map;

typedef std::pair<mpf_class, mpf_class> key_mpf;
typedef std::pair<key_mpf,size_t> key_mpf_map;

double my_abs(const Eigen::Matrix<double,3,1>&);

// for post-processing: rounding
const double MACHINE_EPS=1e-14;
const double EPS_MERGE = pow(2, -50);
const double EPS_MV = pow(2, -50);

inline void keep_key_valid(key_type& k)
{
  if (k.second == 0) {
    std::cout << "# [ MESSAGE ] k.second = 0, set it as eps." << std::endl;
    k.second = MACHINE_EPS;
  }
}
inline void keep_value_valid(double& x)
{
  if (x == 0) {
    std::cout << "# [ MESSAGE ] k.second = 0, set it as eps." << std::endl;
    x = MACHINE_EPS;
  }
}

class myc {
 public:
  myc(const double eps) : eps_(eps) {}
  myc() : eps_(1e-10) {}
  bool operator() (const std::pair<key_type,size_t>& ka,
                   const std::pair<key_type,size_t>& kb) const {
    // key_a > key_b
    auto& a = ka.first;
    auto& b = kb.first;
    return (dsgn(a.first*b.second-a.second*b.first)*dsgn(a.second)*dsgn(b.second)>0);
  }
  bool small_distance(const key_type& a, const key_type& b) const {
    return (fabs(a.first*b.second-a.second*b.first)<eps_);
  }

  
 private:

  double eps_;
};


class pc_myc {
 public:
  pc_myc(const double eps) : eps_(eps) {}
  pc_myc() : eps_(1e-10) {}
  bool operator() (const std::pair<key_type,size_t>& ka,
                   const std::pair<key_type,size_t>& kb) const {
    // key_a > key_b
    auto& a = ka.first;
    auto& b = kb.first;
    return (dsgn(a.first*b.second-a.second*b.first)*dsgn(a.second)*dsgn(b.second)>0);
  }
  
  bool small_distance(const double a1, const double a2, const double b1, const double b2) const {
    return (fabs(a1*b2-a2*b1) < (eps_*fabs(a2)*fabs(b2)+1e-10));
  }
  
  bool small_distance(const key_type& a, const key_type& b) const {
    return small_distance(a.first, a.second, b.first, b.second);
  }

  
 private:

  double eps_;
};


typedef std::priority_queue<key_map,std::vector<key_map>,myc> pri_queue;


///////////////////////////////////////////////
short sgn(double x) {
  if (x<0) return -1;
  if (x>0) return 1;
  return 0;
}

template<typename T>
short key_cmp(const std::pair<T,T>& a, const std::pair<T,T>& b)
{
  return (sgn(a.first*b.second-a.second*b.first)*sgn(a.second)*sgn(b.second));
}

template<typename T1, typename T2>
short key_cmp(const std::pair<T1,T1>& a, const std::pair<T2,T2>& b)
{
  return (sgn(a.first*b.second-a.second*b.first)*sgn(a.second)*sgn(b.second));  
}

/// for sort vector
template<typename T>
class line_sort_myc {
 public:
  line_sort_myc() {}
  typedef std::pair<T, T> KT;
  bool operator() (const std::pair<KT,size_t>& ka,
                   const std::pair<KT,size_t>& kb) const {
    // key_a > key_b
    auto& a = ka.first;
    auto& b = kb.first;
    return (key_cmp(a, b)<0);
  }
};


///////////////

int sequence_type::build_v2f()
{
  size_t tn = tr_.tn();
  size_t vn = tr_.vn();
  auto& M = tr_.mesh().triangles();

  v2f_.resize(vn);
  for (size_t vi = 0; vi < vn; ++vi) {
    v2f_[vi].reserve(6);
  }

  for (size_t ti = 0; ti < tn; ++ti) {
    for (size_t j = 0; j < 3; ++j) {
      v2f_[M[3*ti+j]].push_back(ti);
    }
  }

  return 0;
}

int sequence_type::run()
{
  tn_ = tr_.mesh().faces_num();
  for (size_t i = 0; i < eid_map_.size(); ++i) {
    e2id_.insert(std::make_pair(eid_map_[i]->first, i));
  }
    
  ufs_.reset(cv_.size()+tr_.vn());
  CALL_FUNC(classify());
  CALL_FUNC(sort_pp());
  CALL_FUNC(sort_pt());
  CALL_FUNC(sort_tt());
  CALL_FUNC(reset_cutvert_type());
  CALL_FUNC(keep_pp_inner());

  return 0;
}

int sequence_type::reset_cutvert_type()
{
  size_t vn = tr_.mesh().verts_num();
  size_t all_cvn = vn+cv_.size();
  
  std::vector<char> cv_type(all_cvn, 'o');
  // get element id, and reset element id into cv_
  std::vector<size_t> cv_ele_id(all_cvn, INVALID_NUM);
  
  // set vert-type cv
  for (size_t i = cv_.size(); i<all_cvn; ++i) {
    cv_type[i] = 'v';
    cv_ele_id[i] = i-cv_.size();
    size_t pa = ufs_.find(i);
    if (pa != i) {
      FASSERT(pa < cv_type.size());
      cv_type[pa] = 'v';
      cv_ele_id[pa] = i-cv_.size();
    }
  }
  
  for (size_t i = 0; i < cv_.size(); ++i) {
    size_t pa = ufs_.find(i);
    if (cv_type[i]=='o' && cv_type[pa]=='v') {
      cv_type[i] = 'v';
      cv_ele_id[i] = cv_ele_id[pa];
    }
  }
  
  // set edge-type cv
  for (size_t i = 0; i < cv_.size(); ++i) {
    if (cv_type[i] == 'o' && cv_[i].is_edge()) {
      cv_type[i] = 'e';
      cv_ele_id[i] = cv_[i].ele_id();
      size_t pa = ufs_.find(i);
      FASSERT(cv_type[pa] == 'o' || cv_type[pa] == 'e')
      if (pa != i && cv_type[pa] == 'o') {
        FASSERT(pa < cv_type.size());

        cv_type[pa] = 'e';
        cv_ele_id[pa] = cv_ele_id[i];
      }
    }
  }
  
  for (size_t i = 0; i < cv_.size(); ++i) {// run the second time
    size_t pa = ufs_.find(i);
    if (cv_type[i] == 'o' && cv_type[pa]=='e') {
      cv_type[i] = 'e';
      cv_ele_id[i] = cv_ele_id[pa];
    }
  }
    
  // set triangle-type cv
  for (size_t i = 0; i < cv_.size(); ++i) {
    if (cv_type[i] == 'o' && cv_[i].is_triangle()) {
      cv_type[i] = 't';
      cv_ele_id[i] = cv_[i].ele_id();
      size_t pa = ufs_.find(i);
      if (pa != i) {
        FASSERT(pa < cv_type.size());
        cv_type[pa] = 't';
        cv_ele_id[pa] = cv_ele_id[i];
      }
    }
  }
  for (size_t i = 0; i < cv_.size(); ++i) {
    size_t pa = ufs_.find(i);
    if (cv_type[i] == 'o' && cv_type[pa] == 't') {
      cv_type[i] = 't';
      cv_ele_id[i] = cv_ele_id[pa];
    }
  }

  // reset cv type about bound cv
  for (size_t i = 0; i < cv_.size(); ++i) {
    if (cv_type[i] == 'v') {
      cv_[i].set_on_vert(cv_ele_id[i]);
    } else if (cv_type[i] == 'e') {
      cv_[i].set_on_edge(cv_ele_id[i]);
    } else if (cv_type[i] == 't') {
      cv_[i].set_on_triangle(cv_ele_id[i]);
    }
  }

  // reset grid cv
  std::vector<char> g_cv_type(all_cvn, 'o');
  for (size_t i = 0; i < cv_.size(); ++i) {
    if (cv_[i].is_grid()) {
      g_cv_type[i] = 'g';
      size_t pa = ufs_.find(i);
      if (pa != i) { g_cv_type[pa] = 'g'; }
    }
  }
  for (size_t i = 0; i < cv_.size(); ++i) {
    size_t pa = ufs_.find(i);
    if (g_cv_type[i] == 'o' && g_cv_type[pa] == 'g') {
      g_cv_type[i] = 'g';
    }
  }

  // reset grid point
  for (size_t i = 0; i < cv_.size(); ++i) {
    if (g_cv_type[i] == 'g') {
      cv_[i].set_on_grid();
    }
  }
    
  {
    std::unordered_map<size_t,size_t> id_map;
    for (size_t i = 0; i < vn; ++i) {//first, add original V
      size_t pa = ufs_.find(i+cv_.size());
      id_map.insert(std::make_pair(pa, i));
    }
      
    size_t cnt = vn;
    for (size_t i = 0; i < cv_.size(); ++i) {//then, add cut V
      size_t pa = ufs_.find(i);
      auto it = id_map.find(pa);
      if (it == id_map.end()) {
        id_map.insert(std::make_pair(pa, cnt));
        ++cnt;
      }
    }

    for (size_t i = 0; i < cv_.size(); ++i) {
      auto it = id_map.find(ufs_.find(i));
      FASSERT(it != id_map.end());
      cv_[i].id() = it->second;
    }
    std::cout << "-- no repeated cut verts number: " << cnt << std::endl;
    pure_cv_num_ = cnt;
  }

  { // print out the number of each type.
    size_t cvn_v(0), cvn_e(0), cvn_t(0);
    for (size_t i = 0; i < all_cvn; ++i) {
      if (cv_type[i] == 'v') ++cvn_v;
      else if (cv_type[i] == 'e') ++cvn_e;
      else if (cv_type[i] == 't') ++cvn_t;
      else if (cv_type[i] != 'o') {
        std::cerr << "# [ ERROR ] not v, e, t, o." << std::endl;
      }
    }
    std::cout << "-- cut verts   " << "\n"
              << "     on vert : " << cvn_v << "\n"
              << "     on edge : " << cvn_e << "\n"
              << " on triangle : " << cvn_t << std::endl;
  }
    
  return 0;
}


// unique id -> cv id
const std::vector<size_t>& sequence_type::get_cut_verts_id()
{
  if (CV_id_.size() == pure_cv_num_) { return CV_id_; }

  CV_id_.resize(pure_cv_num_);

  size_t vn = tr_.mesh().verts_num();  

  for (size_t vi = 0; vi < vn; ++vi) {
    CV_id_[vi] = vi+cv_.size();
  }

  std::vector<char> flag(pure_cv_num_, 0);
  std::fill(&flag[0], &flag[0]+vn, 1);

  size_t cnt(0);
  for (auto& v : cv_) {
    FASSERT(v.id() != INVALID_NUM);
    if (flag[v.id()] == 0) {
      flag[v.id()] = 1;
      CV_id_[v.id()] = cnt;
    }
    ++cnt;
  }
  
  return CV_id_;
}

const std::vector<char>& sequence_type::get_cut_verts_type()
{
  if (CV_type_.size() == pure_cv_num_) { return CV_type_; }

  //CV_.resize(pure_cv_num_*3);
  CV_type_.resize(pure_cv_num_, 'g');

  std::vector<char> flag(pure_cv_num_, 0);

  size_t vn = tr_.mesh().verts_num();
  std::fill(&flag[0], &flag[0]+vn, 1);
  std::fill(&CV_type_[0], &CV_type_[0]+vn, 'v');

  for (auto& v : cv_) {
    FASSERT(v.id() != INVALID_NUM);
    if (flag[v.id()] == 0) {
      flag[v.id()] = 1;
      FASSERT(!v.is_vert());
      if (v.is_triangle()) {
        CV_type_[v.id()] = 't';
        FASSERT(!v.is_edge());
      } else if (v.is_edge()) {
        CV_type_[v.id()] = 'e';
        FASSERT(!v.is_triangle());
      }
    }
  }

  return CV_type_;
}

const std::vector<double>& sequence_type::get_cut_verts()
{
  if (CV_.size() == pure_cv_num_*3) { return CV_; }
  
  CV_.resize(pure_cv_num_*3);

  std::vector<char> flag(pure_cv_num_, 0);

  size_t vn = tr_.mesh().verts_num();

  auto& V = tr_.mesh().verts();
  
  std::copy(V.begin(), V.end(), CV_.begin());
  std::fill(&flag[0], &flag[0]+vn, 1);

  for (auto& v : cv_) {
    FASSERT(v.id() != INVALID_NUM);
    if (flag[v.id()] == 0) {
      cal_cv_coord(v, &CV_[3*v.id()]);
      flag[v.id()] = 1;
      FASSERT(!v.is_vert());
    }
  }

  return CV_;
}

int sequence_type::cal_cv_coord(const cut_vert& v, double* coord) const
{
  auto& P = tr_.plane();
  auto& V = tr_.mesh().verts();
  auto& M = tr_.mesh().triangles();
    
  if (v.f(0)>=tn_ && v.f(1)>=tn_ && v.f(2)>=tn_) {//PPP
    size_t p0(v.f(0)-tn_), p1(v.f(1)-tn_), p2(v.f(2)-tn_);
    cal_ppp_point(&P[4*p0], &P[4*p1], &P[4*p2], coord);
  } else if (v.f(0)>=tn_ && v.f(1)>=tn_ && v.f(2)<tn_) {//PPT
    FASSERT(!v.is_vert());//for guarantee,the planes of ppt should be intersected at one point
    size_t p0(v.f(0)-tn_), p1(v.f(1)-tn_), t(v.f(2));
    size_t vi(M[3*t]), vj(M[3*t+1]), vk(M[3*t+2]);

    cal_ppt_point(
        &P[4*p0], &P[4*p1], &V[3*vi], &V[3*vj], &V[3*vk], coord);

  } else if (v.f(0)<tn_ && v.f(1)<tn_ && v.f(2)>=tn_) {//TTP
    FASSERT(!v.is_vert());// must add
    FASSERT(v.is_edge());
    size_t p = v.f(2)-tn_;
    auto& ev = eid_map_[v.ele_id()]->first;
    FASSERT(ev.first < ev.second);
    cal_ttp_point(&V[3*ev.first], &V[3*ev.second], &P[4*p], coord);
  }
    
  return 0;
}

int sequence_type::get_cut_edges(
    std::vector<char>& cv_type,
    std::vector<size_t>& edge_CE,
    std::vector<size_t>& tri_CE,
    std::vector<size_t>& grid_in_CE,
    std::vector<size_t>& grid_out_CE) const
{
    
  {// add tt (edge chains)
    for (auto& e : e2id_) {
      const auto& seq = seq_tt_[e.second];
      size_t pre_id = e.first.first;//it is equal to the final id
      for (auto& v : seq) {
        size_t vid = cv_[v.cvid].id();
        if (vid == pre_id) continue;
        edge_CE.push_back(pre_id);
        edge_CE.push_back(vid);
        pre_id = vid;
      }
      if (pre_id != e.first.second) {
        edge_CE.push_back(pre_id);
        edge_CE.push_back(e.first.second);
      }
    }
  }
    
  {// add pt (triangle lines)
    for (size_t pti = 0; pti < seq_pt_.size(); ++pti) {
      auto& seq = seq_pt_[pti];
      size_t pre_id = cv_[seq_pt_bound_[pti].first].id();
      FASSERT(pre_id != INVALID_NUM);
      for (auto& v : seq) {
        size_t vid = cv_[v.cvid].id();
        if (vid == pre_id) continue;
        tri_CE.push_back(pre_id);
        tri_CE.push_back(vid);
        pre_id = vid;
      }
      size_t right_id = cv_[seq_pt_bound_[pti].second].id();
      FASSERT(right_id != INVALID_NUM);
      if (pre_id != right_id) {
        tri_CE.push_back(pre_id);
        tri_CE.push_back(right_id);
      }
    }
  }

  {// add pp (grid lines)
    for (size_t si = 0; si < seq_pp_.size(); ++si) {
      auto& seq = seq_pp_[si];
      auto& seq_s = pp_s_[si];
      FASSERT(seq.size() == seq_s.size());

      if (seq.size() == 0){continue;}
      if (seq.size() == 1){std::cout<<"-- seq pp length: 1 " << std::endl;continue;}
      for (size_t i = 0; i+1 < seq.size(); ++i) {
        auto& a = cv_[seq[i].cvid];
        auto& b = cv_[seq[i+1].cvid];
        size_t vi = a.id();
        size_t vj = b.id();
        FASSERT(vi != vj);
        if (seq_s[i]==1) {
          grid_in_CE.push_back(vi);
          grid_in_CE.push_back(vj);

          if (!a.is_triangle() && !a.is_edge() && !a.is_vert()) {
            cv_type[vi] = 'i';
          }
          if (!b.is_triangle() && !b.is_edge() && !b.is_vert()) {
            cv_type[vj] = 'i';
          }
            
        } else {
          grid_out_CE.push_back(vi);
          grid_out_CE.push_back(vj);

          if (!a.is_triangle() && !a.is_edge() && !a.is_vert()) {
            cv_type[vi] = 'o';
          }
          if (!b.is_triangle() && !b.is_edge() && !b.is_vert()) {
            cv_type[vj] = 'o';
          }
            
        }
      }
    }
  }
    
  return 0;
}
  

  

int sequence_type::classify()
{
  size_t tn = tr_.tn();

  seq_tt_.resize(e2f_.size());

  // classify the cutted vertices into several subset.
  size_t cid = 0;
  for (const auto& c : cv_) {
    // without merging, the type of cut-vert is determined when inserting.
    if (c.is_grid()) { //PPP, add into (pi,pj)
      FASSERT(c.f(0)>=tn && c.f(1)>=tn && c.f(2)>=tn);
        
      add_pp(cid, c.f(0), c.f(1), c.f(2)); // not minus tn, i.e. c.f(2)-tn
      add_pp(cid, c.f(0), c.f(2), c.f(1));
      add_pp(cid, c.f(1), c.f(2), c.f(0)); // not minus tn, i.e. c.f(0)-tn
        
    } else if (c.is_triangle()) { //PPT,  add into (p,t)
      FASSERT(c.f(0)>=tn && c.f(1)>=tn && c.f(2)<tn);

      //here, only add PPT type points, not TTP (which can be got by topology info)
      add_pt(cid, c.f(0), c.f(2), c.f(1)-tn);
      add_pt(cid, c.f(1), c.f(2), c.f(0)-tn);
      // add into pp
      add_pp(cid, c.f(0), c.f(1), c.f(2)); // triangles
        
    } else if (c.is_edge()) { //TTP, add into (ti,tj)
      FASSERT(c.f(0)<tn && c.f(1)<tn && c.f(2)>=tn);

      add_tt(cid, c.ele_id(), c.f(2)-tn);
      add_pt(cid, c.f(2), c.f(0), INVALID_NUM);// no plane
      add_pt(cid, c.f(2), c.f(1), INVALID_NUM);// no plane
    }
    ++cid;
  }

  // std::cout << "-- seq tt: " << seq_tt_.size() << std::endl;
  // std::cout << "-- seq pt: " << seq_pt_.size() << std::endl;    
  // std::cout << "-- seq pp: " << seq_pp_.size() << std::endl;    

  return 0;
}

void sequence_type::add_pp(size_t cvid, size_t a, size_t b, size_t c) //(a,b),c
{
  std::pair<size_t, size_t> pp(a, b);
  if (pp.first > pp.second) std::swap(pp.first, pp.second);

  auto s = tr_.pp(a-tn_, b-tn_);
  if (s != 2) return; // not intersected at one line.

  auto it = pp_map_.find(pp);
  if (it != pp_map_.end()) {
    seq_pp_[it->second].push_back(pp_point(cvid, c));
  } else {
    pp_map_.insert(std::make_pair(pp, seq_pp_.size()));
    std::vector<pp_point> vec;
    vec.push_back(pp_point(cvid, c));
    seq_pp_.push_back(vec);
  }
}

void sequence_type::add_pt(size_t cvid, size_t a, size_t b, size_t c)
{
  std::pair<size_t, size_t> pt(a, b);
  size_t s = tr_.pt(a-tn_, b);

  if (s == 1) {
    auto& M = tr_.mesh().triangles();
    size_t va(M[3*b]), vb(M[3*b+1]), vc(M[3*b+2]);
    if (tr_.pv(a-tn_, va)==0) {
      ufs_.set_union(cvid, cv_.size()+va);
    } else if (tr_.pv(a-tn_, vb)==0) {
      ufs_.set_union(cvid, cv_.size()+vb);
    } else if (tr_.pv(a-tn_, vc)==0) {
      ufs_.set_union(cvid, cv_.size()+vc);
    } else {
      std::cerr << "# [ ERROR ] pt == 1, but no vert satisfies the condition.\n"
                << "            (va, vb, vc): " << va << ", " << vb << ", " << vc
                << "            pid: " << a-tn_
                << std::endl;
    }

    return;
  } else if (s == 0 || s == 3) {
    return; // not intersected at one line.
  }
    
  //if (s != 2) return; 
    
  auto it = pt_map_.find(pt);
  if (it != pt_map_.end()) {
    if (c != INVALID_NUM) {
      seq_pt_[it->second].push_back(pt_point(cvid, c));
    }
  } else {
    pt_map_.insert(std::make_pair(pt, seq_pt_.size()));//new id
    std::vector<pt_point> vec;
    if (c != INVALID_NUM) {
      vec.push_back(pt_point(cvid, c));
    }
    seq_pt_.push_back(vec);

  }
}

void sequence_type::add_tt(size_t cvid, size_t eid, size_t fid) {
  seq_tt_[eid].push_back(tt_point(cvid, fid));
}
  
int sequence_type::sort_tt()
{
  for (auto& e : e2id_) {
    auto& seq = seq_tt_[e.second];
    if (seq.size()==0) continue;
    FASSERT(e.second < eid_map_.size());
    if(sort_tt_api(e.first, seq)) {
      std::cerr << "# [ ERROR ] sort tt" << std::endl;
      return 1;
    }
  }
  
  return 0;
}


int sequence_type::sort_tt_with_lc(
    const std::pair<size_t, size_t>& e,
    std::vector<tt_point>& seq)
{
  if (seq.size() == 0) {return 0;}
  const auto& V = tr_.mesh().verts();
  const auto& P = tr_.plane();
    
  // key value and index of seq element
  pri_queue seq_map;

  size_t cnt = 0;
  for (auto& cv : seq) {
    key_type kt;
    FASSERT(cv.plane < P.size()/4);
    cal_ttp_point(&V[3*e.first], &V[3*e.second], &P[4*cv.plane], kt);
    keep_key_valid(kt);
    seq_map.push(std::make_pair(kt, cnt));
    ++cnt;
  }
  FASSERT(cnt == seq.size());
    
  std::vector<tt_point> new_seq;
  std::vector<key_type> seq_key;
  new_seq.reserve(seq.size());
  seq_key.reserve(seq.size());
    
  while (!seq_map.empty()) {
    auto s = seq_map.top();
    seq_map.pop();

    new_seq.push_back(seq[s.second]);//info
    seq_key.push_back(s.first);//key
  }

  std::swap(new_seq, seq);

  
  // fine sorting
  bool is_ok = false;
  while (!is_ok) {
    is_ok = true;
    for (size_t si = 0; si+1 < seq.size(); ++si) {
      if (myc(seq_eps_).small_distance(seq_key[si], seq_key[si+1])) {
        int r = cmp_tt_point(
            &V[3*e.first], &V[3*e.second],
            &P[4*seq[si].plane], &P[4*seq[si+1].plane]);
        // r>0, si>si+1; r<0, si<si+1.
        if (r > 0) {
          std::swap(seq[si], seq[si+1]);
          std::swap(seq_key[si], seq_key[si+1]);
          is_ok = false;
        } else if (r == 0) {
          ufs_.set_union(seq[si].cvid, seq[si+1].cvid);//cut verts index
        }
      }
    }//for
  }//while

  // check whether the starting and ending cutted vertices are equal to vl or vr.
  key_type bd_v;
  bd_v.first = 0.0; bd_v.second = 1.0;
  if (myc(seq_eps_).small_distance(seq_key.front(),bd_v)) {
    int r = cmp_tt_point_to_const(
        &V[3*e.first], &V[3*e.second], &P[4*seq.front().plane], 0.0);
    if (r==0) {
      ufs_.set_union(seq.front().cvid, e.first+cv_.size());
    } else if (r < 0) {
      std::cerr << "# [ ERROR ] this error cannot be happened theoretically (tt left)." << std::endl;
      return 1;
    }
  }
    
  bd_v.first = bd_v.second = 1.0;
  if (myc(seq_eps_).small_distance(seq_key.back(),bd_v)) {
    int r = cmp_tt_point_to_const(
        &V[3*e.first], &V[3*e.second], &P[4*seq.back().plane], 1.0);
    if (r==0) {
      ufs_.set_union(seq.back().cvid, e.second+cv_.size());
    } else if (r > 0) {
      std::cerr << "# [ ERROR ] this error cannot be happened theoretically (tt right)." << std::endl;
      return 1;
    }
  }
    
  return 0;
}

int sequence_type::sort_pp()
{
  size_t tn = tr_.mesh().faces_num();
  for (auto& ppm : pp_map_) {
    auto& seq = seq_pp_[ppm.second];
    size_t ori_size = seq.size();
    if (seq.size() == 0) continue;

    auto& pp = ppm.first;
    FASSERT(pp.first>=tn && pp.second>=tn);
    sort_pp_api(pp.first-tn, pp.second-tn, seq);
    FASSERT(ori_size == seq.size());
  }
    
  return 0;
}



///////////////
#define CAL_PP_GMP_KEY_TYPE(F, K)    {                                  \
    if (F>=tn) {                                                        \
      cal_ppp_point(gmp_NSd,gmp_dir,gmp_Nd,&P[4*(F-tn)],K);             \
    } else {                                                            \
      auto it = f2v.find(F);                                            \
      if (it == f2v.end()) {                                            \
        size_t vi(M[3*F]), vj(M[3*F+1]), vk(M[3*F+2]);                  \
        cal_ppt_point(gmp_NSd,gmp_dir,gmp_Nd,&V[3*vi],&V[3*vj],&V[3*vk], K); \
      } else {                                                          \
        cal_ppv_point(gmp_NSd,gmp_dir,gmp_Nd,&V[3*it->second], K);      \
      }                                                                 \
    }}
///////////////
  
int sequence_type::sort_pp_with_lc(size_t pa, size_t pb, std::vector<pp_point>& seq)
{
  size_t tn = tr_.mesh().faces_num();
  auto& M = tr_.mesh().triangles();
  auto& V = tr_.mesh().verts();
  auto& P = tr_.plane();

  // define a starting point by choosing a axis-plane
  double sP[4];
  sP[0] = sP[1] = sP[2] = sP[3] = 0.0;
  {
    double max = 0.0;
    size_t k = 0;
    for (size_t i = 0; i < 3; ++i) {
      sP[i] = 1.0;
      double vol = fabs(cal_vol(&P[4*pa],&P[4*pb],sP));
      if (vol > max) { max = vol; k = i;}
      sP[i] = 0.0;
    }
    sP[k] = 1.0;
  }

  Eigen::Matrix<double,3,1> NSd, dir;
  double Nd;
  pp_precompute(&P[4*pa], &P[4*pb], sP, NSd, dir, Nd);

  std::unordered_map<size_t,size_t> f2v;
    
  //std::map<key_type, size_t, myc> seq_map;
  pri_queue seq_map;
    
  size_t cnt(0);
  for (auto& v : seq) {
    key_type kt;
    if (v.face >= tn) { // PPP
      cal_ppp_point(NSd,dir,Nd,&P[4*(v.face-tn)], kt);
    } else { // PPT, precondition, there exists only one intersection point between them.
      size_t fid = v.face;

      size_t ppv = INVALID_NUM;
      for (size_t k = 0; k < 3; ++k) {
        if (tr_.pv(pa, M[3*fid+k])==0 && tr_.pv(pb, M[3*fid+k])==0) {
          ppv = M[3*fid+k];
          f2v.insert(std::make_pair(fid, ppv));
          break;
        }
      }
        
      if (ppv == INVALID_NUM) {
        size_t vi(M[3*fid]), vj(M[3*fid+1]), vk(M[3*fid+2]);
        cal_ppt_point(NSd,dir,Nd,&V[3*vi],&V[3*vj],&V[3*vk], kt);
      } else {
        cal_ppv_point(NSd,dir,Nd,&V[3*ppv], kt);
      }
    }
    keep_key_valid(kt);
    seq_map.push(std::make_pair(kt, cnt));
    ++cnt;
  }
    
  std::vector<pp_point> new_seq;
  std::vector<key_type> seq_key;
  new_seq.reserve(seq.size());
  seq_key.reserve(seq.size());
    
  while (!seq_map.empty()) {
    auto s = seq_map.top();
    seq_map.pop();
    new_seq.push_back(seq[s.second]);//info
    seq_key.push_back(s.first);//key
  }

  FASSERT(seq.size() == new_seq.size());
    
  std::swap(new_seq, seq);


  Eigen::Matrix<mpf_class,3,1> gmp_NSd, gmp_dir;
  mpf_class gmp_Nd;
  pp_precompute(&P[4*pa], &P[4*pb], sP, gmp_NSd, gmp_dir, gmp_Nd);

  // fine sorting
  bool is_ok = false;
  std::unordered_map<size_t, std::pair<mpf_class,mpf_class> > face_key;
  while (!is_ok) {
    is_ok = true;
    for (size_t si = 0; si+1 < seq.size(); ++si) {
      if (myc(seq_eps_).small_distance(seq_key[si], seq_key[si+1])) {
        size_t fi = seq[si].face;
        size_t fj = seq[si+1].face;
        auto k1_it = face_key.find(fi);
        auto k2_it = face_key.find(fj);
        if (k1_it == face_key.end()) {
          std::pair<mpf_class,mpf_class> ke;
          CAL_PP_GMP_KEY_TYPE(fi,ke);
          FASSERT(ke.second != 0);
          face_key.insert(std::make_pair(fi,ke));
          k1_it = face_key.find(fi);
        }
        if (k2_it == face_key.end()) {
          std::pair<mpf_class,mpf_class> ke;          
          CAL_PP_GMP_KEY_TYPE(fj,ke);
          FASSERT(ke.second != 0);          
          face_key.insert(std::make_pair(fj,ke));
          k2_it = face_key.find(fj);
        }
        auto& k1 = k1_it->second;
        auto& k2 = k2_it->second;


        int k12s =sgn(k1.second)*sgn(k2.second);
        FASSERT(k12s != 0);
        int r = sgn(k1.first*k2.second - k1.second*k2.first)*k12s;
        
        // r>0, si>si+1; r<0, si<si+1.
        if (r > 0) {
          std::swap(seq[si], seq[si+1]);
          std::swap(seq_key[si], seq_key[si+1]);
          is_ok = false;
        } else if (r == 0) {
          ufs_.set_union(seq[si].cvid, seq[si+1].cvid);//cut verts index
        }
      }
    }//for
  }//while


  return 0;
}

void get_edge_direc(const double* a,  const double* b, MPF_VEC& d)
{
  MPF_VEC A, B;
  A << a[0], a[1], a[2];
  B << b[0], b[1], b[2];
  d = B-A;
}

void get_plane2_direc(const double* pa, const double* pb, MPF_VEC& d)
{
  MPF_VEC Na, Nb;
  Na << pa[0], pa[1], pa[2];
  Nb << pb[0], pb[1], pb[2];
  d = Na.cross(Nb);
}

void get_pt_direc(
    const double* pa, const double* va, const double* vb, const double* vc,
    MPF_VEC& d)
{
  MPF_VEC N, A, B, C;
  N << pa[0], pa[1], pa[2];
  A << va[0], va[1], va[2];
  B << vb[0], vb[1], vb[2];
  C << vc[0], vc[1], vc[2];

  MPF_VEC Nt = (B-A).cross(C-A);

  d = N.cross(Nt);
}


int sequence_type::keep_pp_inner()
{
  //std::cout << "# [ WARNING ] now the adjacent triangles of vertices are got by global search."
  //          << std::endl;
    
  pp_s_.resize(seq_pp_.size());

  for (auto& ppm : pp_map_) {
    FASSERT(ppm.first.first>=tn_ && ppm.first.second>=tn_);
    if (keep_pp_inner(ppm.first.first-tn_, ppm.first.second-tn_, ppm.second)) {
      return 1;
    }
  }
  return 0;
}

size_t sequence_type::find_third_vert(size_t va, size_t vb, const size_t* M)
{
  for (size_t k = 0; k < 3; ++k) {
    if (M[k]==va && M[(k+1)%3]==vb) {
      return M[(k+2)%3];
    }
  }
  return INVALID_NUM;
}

uint8_t sequence_type::passing_triangle(size_t pa, size_t pb, size_t fid)
{
  auto& P = tr_.plane();
  auto& V = tr_.mesh().verts();
  auto& M = tr_.mesh().triangles();
  size_t vi(M[3*fid]), vj(M[3*fid+1]), vk(M[3*fid+2]);
  return judge_triangle_pass(&P[4*pa], &P[4*pb], &V[3*vi], &V[3*vj], &V[3*vk]);
}

uint8_t sequence_type::passing_edge(size_t eid, size_t pa, size_t pb, const std::vector<size_t>& fs)
{
  auto& P = tr_.plane();
  auto& V = tr_.mesh().verts();
  auto& M = tr_.mesh().triangles();

  size_t va = eid_map_[eid]->first.first;
  size_t vb = eid_map_[eid]->first.second;

  if (tr_.pv(pa,va)==0 && tr_.pv(pa,vb)==0
      && tr_.pv(pb,va)==0 && tr_.pv(pb,vb)==0) {
    return 0;
  }

    
  size_t s = 0;
  for (size_t i = 0; i+1 < fs.size(); i+=2) {
    size_t fa(fs[i]), fb(fs[i+1]);
    FASSERT(va < vb);
    FASSERT(fa < tn_ && fb < tn_);

    size_t vc = find_third_vert(va, vb, &M[3*fa]);
    size_t vd = find_third_vert(vb, va, &M[3*fb]);
    FASSERT(vc != INVALID_NUM);
    FASSERT(vd != INVALID_NUM);
    s+=judge_edge_pass(&P[4*pa], &P[4*pb], &V[3*va], &V[3*vb], &V[3*vc], &V[3*vd]);
  }
    
  return (s&1);
}


  
uint8_t sequence_type::passing_vert(size_t vid, size_t pa, size_t pb)
{
  const auto& fs = v2f_[vid];

  std::vector<size_t> cut_fs;
  std::vector<size_t> bd_fs;
  cut_fs.reserve(fs.size());
  for (auto& f : fs) {
    size_t r = tr_.pt(pa, f);
    if (r==2) {
      cut_fs.push_back(f);
    } else if (r==3) {
      bd_fs.push_back(f);
    }
  }

    
  if (cut_fs.size() == 0) { return 0;}

  // sort bound fs
  std::vector<std::vector<size_t> > bd_fans;
  cycle_sort_verts(vid, pa, bd_fs, bd_fans);
  if (bd_fans.size() == 1 && bd_fans[0].front()==bd_fans[0].back()) {

    return 0;
  }


  FASSERT(cut_fs.size()!=1);
  
  std::unordered_map<size_t,size_t> fs2id;
  for (size_t i = 0; i < cut_fs.size(); ++i) {
    fs2id.insert(std::make_pair(cut_fs[i], i));
  }
    
  std::vector<MPF_VEC> dt;
  cal_pt_outer_direc(vid, pa, cut_fs, dt);

  FASSERT(dt.size() == cut_fs.size());

  std::vector<char> cp;
  std::vector<short> dir;

  // if not success, return 5
  if (cycle_sort_triangles(vid, pa, dt, fs2id, cut_fs, cp, dir)) {
    std::cerr << "# [ ERROR ] cycle sort triangles." << std::endl;
    ERROR_RETURN;
  }

  return passing_cycle_fan(vid, bd_fans, pa, pb, dt, fs2id, cut_fs, cp, dir);
}

void sequence_type::cal_pt_outer_direc(
    size_t vid, size_t pid,
    const std::vector<size_t>& fs,
    std::vector<MPF_VEC>& dt) const
{
  auto& V = tr_.mesh().verts();
  auto& M = tr_.mesh().triangles();
  auto& P = tr_.plane();
    
  dt.reserve(fs.size());


  MPF_VEC VO, NP;
  VO << V[3*vid],V[3*vid+1],V[3*vid+2];
  NP << P[4*pid], P[4*pid+1], P[4*pid+2];

  for (auto f : fs) {
    size_t va(INVALID_NUM), vb(INVALID_NUM);
    for (size_t i = 0; i < 3; ++i) {
      if (M[3*f+i] == vid) {
        va = M[3*f+(i+1)%3];
        vb = M[3*f+(i+2)%3];
        break;
      }
    }
    FASSERT(va!=INVALID_NUM && vb!=INVALID_NUM);
    if (tr_.pv(pid, va) < tr_.pv(pid, vb)) {
      std::swap(va, vb);
    }

    MPF_VEC VA,VB;
    VA << V[3*va], V[3*va+1], V[3*va+2];
    VB << V[3*vb], V[3*vb+1], V[3*vb+2];
    dt.push_back(((VA-VO).cross(VB-VO)).cross(NP));
  }
}


// > 2 adjacent triangles
int sequence_type::cycle_sort_triangles(
    size_t vid, size_t pid,
    const std::vector<MPF_VEC>& dt,
    const std::unordered_map<size_t,size_t>& fs2id,
    std::vector<size_t>& fs,
    std::vector<char>& cp, std::vector<short>& dir,
    cycle_sort_type cst) const
{
  FASSERT(fs.size() != 1);
  std::vector<size_t> neg_id, zero_begin_id, zero_end_id, one_id;
  zero_begin_id.push_back(0);

  auto& V = tr_.mesh().verts();
  auto& P = tr_.plane();
  auto& M = tr_.mesh().triangles();
    
  MPF_VEC NP;
  NP << P[4*pid], P[4*pid+1], P[4*pid+2];
    
  for (size_t i = 1; i < fs.size(); ++i) {
    int s = sgn(dt[0].cross(dt[i]).dot(NP));
    if (s < 0) {neg_id.push_back(i);}
    else if (s == 0) {
      int ss = sgn(dt[i].dot(dt[0]));
      if (ss > 0) {
        zero_begin_id.push_back(i);
      } else {
        zero_end_id.push_back(i);
      }
    } else if (s > 0) {one_id.push_back(i);}
  }

  if (neg_id.size() > 1) {
    for (size_t i = 0; i < neg_id.size(); ++i) {
      for (size_t j = i+1; j < neg_id.size(); ++j) {
        int s = sgn(dt[neg_id[i]].cross(dt[neg_id[j]]).dot(NP));
        if (s > 0) {
          std::swap(neg_id[i], neg_id[j]);
        }
      }
    }
  }

  if (one_id.size() > 1) {
    for (size_t i = 0; i < one_id.size(); ++i) {
      for (size_t j = i+1; j < one_id.size(); ++j) {
        auto s = sgn(dt[one_id[i]].cross(dt[one_id[j]]).dot(NP));
        if (s > 0) {
          std::swap(one_id[i], one_id[j]);
        }
      }
    }
  }


  // (vid, va, vb) -> Nt
  std::vector<size_t> new_fs;
  new_fs.reserve(fs.size());
  cp.clear(); cp.reserve(fs.size());

  if (zero_begin_id.size() > 0) {
    for (size_t i = 0; i < zero_begin_id.size(); ++i) {
      cp.push_back('=');
      new_fs.push_back(fs[zero_begin_id[i]]);
    }
    cp.back() = '<';
  }


  if (neg_id.size() > 0) { // sort the negative part
    for (size_t i = 0; i+1 < neg_id.size(); ++i) {
      auto s = sgn(dt[neg_id[i]].cross(dt[neg_id[i+1]]).dot(NP));
      new_fs.push_back(fs[neg_id[i]]);
      if (s < 0) { cp.push_back('<'); }
      else if (s == 0) { cp.push_back('='); }
      else {
        std::cerr << "# [ ERROR ] s > 0." << std::endl;
        return 1;
      }
    }
    new_fs.push_back(fs[neg_id.back()]);
    cp.push_back('<');
  }


  if (zero_end_id.size() > 0) {
    // push the zero part in the middle
    for (size_t i = 0; i < zero_end_id.size(); ++i) {
      cp.push_back('=');
      new_fs.push_back(fs[zero_end_id[i]]);
    }
    cp.back() = '<';
  }



  if (one_id.size() > 0) { // sort the positive part
    for (size_t i = 0; i+1 < one_id.size(); ++i) {
      new_fs.push_back(fs[one_id[i]]);
      auto s = sgn(dt[one_id[i]].cross(dt[one_id[i+1]]).dot(NP));
      if (s < 0) { cp.push_back('<'); }
      else if (s == 0) { cp.push_back('='); }
      else {
        std::cerr << "# [ ERROR ] s > 0." << std::endl;
        return 1;
      }
    }
    new_fs.push_back(fs[one_id.back()]);
    cp.push_back('<');
  }



  if (cst == cycle_sort_type::FULL) {
    FASSERT(new_fs.size()==cp.size());
    // sort equal triangles, and remove unuseful triangles    
    sort_equal_triangles(vid, pid, new_fs, cp);
    FASSERT(new_fs.size()==cp.size());
    std::swap(new_fs, fs);
    FASSERT(fs.size() != 1);

    dir.reserve(fs.size());  

    // set the half-edge direction
    for (auto f : fs) {

      int r = plane_triangle_edge_direc(fs2id, dt, V, M, NP, f);
      if (r > 0) {
        dir.push_back(1);
      } else {
        dir.push_back(-1);
      }
    }
    
  } else { // simple sort, because the boundary redundant half-edges have been cleared before.
    std::swap(new_fs, fs);
    dir.reserve(fs.size());  
    for (auto f : fs) {
      int r = plane_triangle_edge_direc(fs2id, dt, V, M, NP, f);
      dir.push_back(r>0?1:-1);
    }

    FASSERT(fs.size() > 2);
    
    if (cp[0] == '=' && dir[0]==-1) {
      std::swap(fs[0], fs[1]);
      std::swap(dir[0], dir[1]);
    }

    // forward adjust
    for (size_t i = 1; i+1 < fs.size();) {
      if (cp[i] == '=') {
        size_t& fa = fs[i];
        size_t& fb = fs[i+1];
        FASSERT(dir[i] != dir[i+1]);
        if (dir[i]==dir[i-1] && dir[i]==-1) {
          std::swap(fa, fb);
          std::swap(dir[i],dir[i+1]);
        }
        i+=2;
      } else { ++i; }
    }

    // backford adjust
    for (int i = static_cast<int>(fs.size())-3; i>=0;) {
      FASSERT(i>=0);
      if (cp[i] == '=') {
        size_t& fa = fs[i];
        size_t& fb = fs[i+1];
        FASSERT(dir[i] != dir[i+1]);
        if (dir[i+1]==dir[i+2] && dir[i]==-1) {
          std::swap(fa, fb);
          std::swap(dir[i],dir[i+1]);
        }
        i-=2;
      } else { --i; }
    }
    
  }


  return 0;
}

int sequence_type::plane_triangle_edge_direc(
    const std::unordered_map<size_t,size_t>& fs2id,
    const std::vector<MPF_VEC>& dt,
    const std::vector<double>& V,
    const std::vector<size_t>& M,
    const MPF_VEC& NP, size_t f) const
{
  auto it = fs2id.find(f);
  FASSERT(it != fs2id.end() && it->second < dt.size());
  auto& d = dt[it->second];

  size_t va = M[3*f];
  size_t vb = M[3*f+1];
  size_t vc = M[3*f+2];

  MPF_VEC VA, VB, VC;
  VA << V[3*va], V[3*va+1], V[3*va+2];
  VB << V[3*vb], V[3*vb+1], V[3*vb+2];
  VC << V[3*vc], V[3*vc+1], V[3*vc+2];

  return sgn(NP.cross((VB-VA).cross(VC-VA)).dot(d));
}


int sequence_type::sort_equal_triangles(
    size_t vid, size_t pid,
    std::vector<size_t>& fs,
    std::vector<char>& cp) const
{
  FASSERT(fs.size() == cp.size());
  std::vector<size_t> new_fs;
  std::vector<char> new_cp;
  new_fs.reserve(fs.size());
  new_cp.reserve(cp.size());
    
  for (size_t i = 0; i < fs.size();) {
    std::vector<size_t> eq_fs;
    eq_fs.push_back(fs[i]);
    size_t j = i;
    for (; j+1 < cp.size(); ++j) {
      if (cp[j]=='=') {
        eq_fs.push_back(fs[j+1]);
      } else {
        break;
      }
    }

    if (eq_fs.size() > 1) {
      sort_equal_triangles(vid, pid, eq_fs);
    }
      
    FASSERT(eq_fs.size() > 0);
      
    for (auto f : eq_fs) {
      new_fs.push_back(f);
      new_cp.push_back('=');
    }
    new_cp.back() = '<';
      
    i = j+1;
  }
    
  std::swap(new_fs, fs);
  std::swap(new_cp, cp);
    
  return 0;
}

// in the plane and vert
int sequence_type::sort_equal_triangles(
    size_t vid, size_t pid,
    std::vector<size_t>& fs) const
{
  std::vector<size_t> up_fs, down_fs;
  std::vector<size_t> up_v3, down_v3;
  up_fs.reserve(fs.size());
  down_fs.reserve(fs.size());
  up_v3.reserve(fs.size());
  down_v3.reserve(fs.size());

  auto& V = tr_.mesh().verts();
  auto& M = tr_.mesh().triangles();

  for (auto f : fs) {
    int s[3];
    s[0] = tr_.pv(pid, M[3*f]);
    s[1] = tr_.pv(pid, M[3*f+1]);
    s[2] = tr_.pv(pid, M[3*f+2]);
      
    if (s[0]>=0 && s[1]>=0 && s[2]>=0) {
      up_fs.push_back(f);
      for (size_t k = 0; k < 3; ++k) {
        if (s[k] > 0) {up_v3.push_back(M[3*f+k]);break;}
      }
    } else if (s[0]<=0 && s[1]<=0 && s[2]<=0) {
      down_fs.push_back(f);
      for (size_t k = 0; k < 3; ++k) {
        if (s[k] < 0) {down_v3.push_back(M[3*f+k]);break;}
      }
    } else {
      std::cerr << "# [ ERRRO ] sort equal triangles." << std::endl;
      return 1;
    }
  }


  FASSERT(up_v3.size()==up_fs.size() && down_v3.size()==down_fs.size());

  size_t ve(INVALID_NUM);
  for (size_t k = 0; k < 3; ++k) {
    size_t id = M[3*fs[0]+k];
    if (id!=vid && tr_.pv(pid,id)==0) {
      ve = id; break;
    }
  }

  int global_s = 1;

  if ( up_fs.size()==0 || (down_fs.size()<up_fs.size()&&down_fs.size()>0) ) {
    std::swap(up_fs, down_fs);
    std::swap(up_v3, down_v3);
    global_s = -1;
  }

  FASSERT(up_fs.size()>0);

  for (size_t i = 0; i < up_fs.size(); ++i) {
    for (size_t j=i+1; j < up_fs.size(); ++j) {
      int s = global_s*vol_sgn(&V[3*vid],&V[3*ve],&V[3*up_v3[i]],&V[3*up_v3[j]]);
      if (s < 0) {
        std::swap(up_fs[i], up_fs[j]);
        std::swap(up_v3[i], up_v3[j]);
      }
    }
  }

  std::swap(up_fs, fs);

  return 0;
}

int sequence_type::cycle_sort_verts(
    size_t vid, size_t pid, const std::vector<size_t>& bd_fs,
    std::vector<std::vector<size_t> >& bd_fans)
{
  if (bd_fs.size() == 0) return 0;    


  auto& M = tr_.mesh().triangles();
  auto& V = tr_.mesh().verts();
  auto& P = tr_.plane();

  std::unordered_set<size_t> third_vs;
  std::set<std::pair<size_t,size_t> > edges;
  for (auto& f : bd_fs) {
    std::vector<size_t> other_vs;
    other_vs.reserve(2);
    for (size_t i = 0; i < 3; ++i) {
      size_t id = M[3*f+i];
      if (id != vid) {
        other_vs.push_back(id);
        third_vs.insert(id);
      }
    }
    FASSERT(other_vs.size() == 2);
    edges.insert(std::make_pair(other_vs[0], other_vs[1]));
    edges.insert(std::make_pair(other_vs[1], other_vs[0]));
  }
  // cycle-sort verts
  std::vector<size_t> vs(third_vs.begin(), third_vs.end());

  FASSERT(vs.size() != 0);


  std::vector<int> s(vs.size());

  for (size_t i = 1; i < vs.size(); ++i) {
    s[i] = two_edges_order_on_plane(&V[3*vid], &V[3*vs[0]], &V[3*vs[i]], &P[4*pid]);
  }
  std::vector<size_t> new_vs;
  new_vs.reserve(vs.size());
    
  new_vs.push_back(vs[0]);
    
  for (size_t i = 1; i < vs.size(); ++i) {
    if (s[i]<0) {
      new_vs.push_back(vs[i]);
    }
  }
  for (size_t i = 1; i < new_vs.size(); ++i) {
    for (size_t j = i+1; j < new_vs.size(); ++j) {
      int sij = two_edges_order_on_plane(
          &V[3*vid], &V[3*new_vs[i]], &V[3*new_vs[j]], &P[4*pid]);
      if (sij > 0) {
        std::swap(new_vs[i], new_vs[j]);
      }
      FASSERT(sij != 0);
    }
  }
  for (size_t i = 1; i < vs.size(); ++i) {
    if (s[i] == 0) {
      new_vs.push_back(vs[i]);
    }
  }
  size_t last_start = new_vs.size();
  for (size_t i = 1; i < vs.size(); ++i) {
    if (s[i] > 0) {
      new_vs.push_back(vs[i]);
    }
  }
  for (size_t i = last_start; i < new_vs.size(); ++i) {
    for (size_t j = i+1; j < new_vs.size(); ++j) {
      int sij = two_edges_order_on_plane(
          &V[3*vid], &V[3*new_vs[i]], &V[3*new_vs[j]], &P[4*pid]);
      if (sij > 0) {
        std::swap(new_vs[i], new_vs[j]);
      }
      FASSERT(sij != 0);
    }
  }


  std::vector<char> st(new_vs.size(), 'm');
  size_t start_id = INVALID_NUM;
  for (size_t i = 0; i < new_vs.size(); ++i) {
    size_t pre_id = new_vs[(i+new_vs.size()-1)%new_vs.size()];
    size_t id = new_vs[i];
    size_t nex_id = new_vs[(i+1)%new_vs.size()];

    int s_L = two_edges_order_on_plane(
        &V[3*vid], &V[3*pre_id], &V[3*id], &P[4*pid]);
    int s_R = two_edges_order_on_plane(
        &V[3*vid], &V[3*id], &V[3*nex_id], &P[4*pid]);

    
    if ((edges.find(std::make_pair(pre_id,id))==edges.end() || s_L>0)
        && (edges.find(std::make_pair(id,nex_id))!=edges.end() && s_R<0)) {
      st[i] = 's';
      if (start_id == INVALID_NUM) {
        start_id = i;
      }
    } else if ((edges.find(std::make_pair(pre_id,id))!=edges.end()&&s_L<0)
               && (edges.find(std::make_pair(id,nex_id))==edges.end()||s_R>0)) {
      st[i] = 'e';
    } else if (edges.find(std::make_pair(pre_id,id))==edges.end()
               && edges.find(std::make_pair(id,nex_id))==edges.end()) {
      std::cerr << "# [ ERROR ] can not happen!" << std::endl;
    }
    
  }

  if (start_id != INVALID_NUM) { // full bound fan
    std::vector<size_t> temp_vs;
    size_t n = new_vs.size();

    size_t i = start_id;
    do {
      FASSERT(st[i] == 's');
      if (st[i]=='s') {
        temp_vs.clear();
        temp_vs.push_back(new_vs[i]);
        size_t j = (i+1)%n;
        for (; j != start_id; j=(j+1)%n) {
          if (st[j]!='s') {
            temp_vs.push_back(new_vs[j]);
          } else {
            break;
          }
        }
        i = j;
        bd_fans.push_back(temp_vs);
      }
    } while (i != start_id);
  } else {
    bd_fans.push_back(new_vs);
    bd_fans.front().push_back(new_vs.front());
  }

  //std::cout << "--- end cycle sort verts." << std::endl;

  return 0;
}

bool sequence_type::is_bound_fan(
    size_t vid, const std::vector<std::vector<size_t> >& bd_fans,
    size_t pid,
    size_t fa, size_t fb)
{
    
  if (bd_fans.size() == 0) {return false;}
  auto& M = tr_.mesh().triangles();
    
  size_t va(INVALID_NUM);
  for (size_t i = 0; i<3; ++i) {
    size_t id = M[3*fa+i];
    if (id!=vid && tr_.pv(pid, id)==0) {
      va = id;
    }
  }
  size_t vb(INVALID_NUM);
  for (size_t i = 0; i<3; ++i) {
    size_t id = M[3*fb+i];
    if (id!=vid && tr_.pv(pid, id)==0) {
      vb = id;
    }
  }

  for (auto& fan : bd_fans) {
    for (size_t i=0; i < fan.size(); ++i) {
      for (size_t j = i+1; j < fan.size(); ++j) {
        if (fan[i]==va && fan[j]==vb) {
          return true;
        }
      }
    }
  }
    
  return false;
}

uint8_t sequence_type::passing_cycle_fan(
    size_t vid,
    const std::vector<std::vector<size_t> >& bd_fans,
    size_t pa, size_t pb,
    const std::vector<MPF_VEC>& dt,
    std::unordered_map<size_t,size_t>& fs2id,
    const std::vector<size_t>& cut_fs,
    const std::vector<char>& cp,
    const std::vector<short>& dir)
{
  FASSERT(cut_fs.size()!=1);
  FASSERT(cp.size()==dir.size() && cut_fs.size() == dir.size());

  auto& P = tr_.plane();

  int sa(0), sb(0);

  FASSERT(4*pa < P.size() && 4*pb < P.size());

  MPF_VEC NA, NB;
  NA << P[4*pa], P[4*pa+1], P[4*pa+2];
  NB << P[4*pb], P[4*pb+1], P[4*pb+2];
  MPF_VEC L = NA.cross(NB);

    
  for (size_t i = 0; i+1 < cut_fs.size(); ++i) {
    if (dir[i] == -1 && dir[i+1] == 1) {
      if (is_bound_fan(vid, bd_fans, pa, cut_fs[i], cut_fs[i+1])) continue; // boundary fan

      auto& di = dt[fs2id[cut_fs[i]]];
      auto& dj = dt[fs2id[cut_fs[i+1]]];
      if (cp[i]!='=') {
        if (sa == 0) {
          sa |= judge_fan(NA, di, dj, L, 1); //di, dj outerside
        }
        if (sb == 0) {
          sb |= judge_fan(NA, di, dj, L, -1);
        }
      } 
    }
  }

  // deal with the last one
  if (dir.back()==-1 && dir.front()==1) {
    auto& di = dt[fs2id[cut_fs.back()]];
    auto& dj = dt[fs2id[cut_fs.front()]];
    bool is_equal = true;
    for (size_t i = 0; i+1 < cp.size(); ++i) {
      if (cp[i] != '=') { is_equal = false; break; }
    }
    if (!is_equal) {
      // check whether the fan is on the boundary
      if (!is_bound_fan(vid, bd_fans, pa, cut_fs.back(), cut_fs.front())) {
        if (sa == 0) {
          sa |= judge_fan(NA, di, dj, L, 1);
        }
        if (sb == 0) {
          sb |= judge_fan(NA, di, dj, L, -1);
        }
      }

    } else { // the same edge, special case
      if (sa == 0) {
        sa |= judge_full_fan(NA, di, L, 1);
      }
      if (sb == 0) {
        sb |= judge_full_fan(NA, di, L, -1);
      }

    } 
  }
    
  if ((sa>0&&sb==0) || (sa==0&&sb>0)) return 1;
    
  return 0;
}

// 1: inner, 0: bound/outer
// fan: <di, dj>, pass line: dk*ks
int sequence_type::judge_fan(
    const MPF_VEC& N, const MPF_VEC& di, const MPF_VEC& dj,
    const MPF_VEC& dk, int ks) const
{
  int sij = -sgn(di.cross(dj).dot(N));
  int sik = -ks*sgn(di.cross(dk).dot(N));
  int sjk = ks*sgn(dj.cross(dk).dot(N));
  int dij = -sgn(di.dot(dj));

  if ((sij>0&&sik>0&&sjk>0) || (sij<0&&(sik>0||sjk>0))
      || (sij==0&&dij>0&&sik>0)) {
    return 1;
  }
  return 0;
}

int sequence_type::judge_full_fan(
    const MPF_VEC& N, const MPF_VEC& di,
    const MPF_VEC& dk, int ks) const
{
  int sik = sgn(di.cross(dk).dot(N));
  if (sik!=0) return 1;
  int dik = -ks*sgn(di.dot(dk));
  if (dik>0) return 1;
  return 0;
}
  
int sequence_type::keep_pp_inner(size_t pa, size_t pb, size_t seq_id)
{
  auto& seq = seq_pp_[seq_id];
  if (seq.size() == 0) {return 0;}

  std::vector<pp_point> shrink_seq;
  shrink_seq.reserve(seq.size());

  auto& seq_s = pp_s_[seq_id];

    
  seq_s.reserve(seq.size());

  size_t pre_id = INVALID_NUM;

  for (size_t i = 0; i < seq.size(); ++i) {
    auto& v = cv_[seq[i].cvid];
    FASSERT(v.id() != INVALID_NUM);
    if (v.id() == pre_id) {// the same id, not need to recalculate
      continue;
    }

    shrink_seq.push_back(seq[i]);
    
    uint8_t s=0;
    if (v.is_triangle()) {
      size_t fid = v.ele_id();
      s = passing_triangle(pa, pb, fid);

    } else if (v.is_edge()) {
      size_t eid = v.ele_id();
      auto& fs = eid_map_[eid]->second;
      s = passing_edge(eid, pa, pb, fs);

    } else if (v.is_vert()) {
      size_t vid = v.ele_id();
      s = passing_vert(vid, pa, pb);
      
    }

    seq_s.push_back(s);
    pre_id = v.id();
  }

  FASSERT(seq_s.size() == shrink_seq.size());


  size_t sum_s(0);
  for (size_t i = 0; i < seq_s.size(); ++i) {
    sum_s += seq_s[i];
    seq_s[i] = (sum_s&1);
  }

  
  FASSERT(seq_s.back() == 0);

  if (shrink_seq.size() == 1) {
    std::cerr << "-- shrink seq size 1, original size " << seq.size() << std::endl;
  }
    
  std::swap(shrink_seq, seq);

  {  // for test, after shrink
    if (seq_s.back() != 0) {// print detailed information for error
      std::cerr << "... [ ERROR ] seq s back is not zero : " << pa << ", " << pb << "  | "<< std::endl;
      size_t pp_id = seq_id;
      const auto& Q = seq_pp_[pp_id];
      for (size_t i = 0; i < Q.size(); ++i) {
        size_t cvid  = Q[i].cvid;
        std::cout << cv_[cvid].id() << " , ";
      }
      std::cout << std::endl;
      for (size_t i = 0; i < Q.size(); ++i) {
        size_t cvid  = Q[i].cvid;
        char ty('o');
        ty = cv_[cvid].is_triangle() ? 't': ty;
        ty = cv_[cvid].is_edge() ? 'e': ty;        
        ty = cv_[cvid].is_vert() ? 'v': ty;    
        std::cout << cv_[cvid].id() << " (cvid " << cvid << ", " << ty << ") , ";
      }
      std::cout << std::endl;

      const auto& S = pp_s_[pp_id];
      for (size_t i = 0; i < S.size(); ++i) {
        std::cout << (size_t)S[i] << ", ";
      }
      std::cout << std::endl;
      
      for (size_t i = 0; i < Q.size(); ++i) {
        auto& v = cv_[Q[i].cvid];
        std::cout << " (" << v.f0() << ", " << v.f1() << ", " << v.f2() << ") ";
      }
      std::cout << std::endl;
    }
  }
    
  return 0;
}

size_t sequence_type::find_edge_cv(size_t va, size_t vb, size_t pid)
{
  std::pair<size_t,size_t> e(va,vb);
  if (e.first > e.second) std::swap(e.first, e.second);
  size_t id = e2id_[e];
  for (auto& tt : seq_tt_[id]) {
    if (tt.plane == pid) {
      return tt.cvid;
    }
  }
  return INVALID_NUM;
}

int sequence_type::sort_pt()
{
  seq_pt_bound_.resize(seq_pt_.size(), std::make_pair(INVALID_NUM,INVALID_NUM));
  seq_pt_dir_.resize(seq_pt_.size(), 'n');
    
  for (auto& ptm : pt_map_) {
    auto& pt = ptm.first;
    FASSERT(pt.first>=tn_ && pt.second<tn_);
    size_t pt_id = ptm.second;

    sort_pt_api(pt_id, pt.first-tn_, pt.second, seq_pt_[pt_id]);
  }
    
  return 0;
}

int sequence_type::sort_pt_with_lc(
    size_t pt_id, size_t pid, size_t tid, std::vector<pt_point>& seq)
{
  const auto& V = tr_.mesh().verts();
  const auto& M = tr_.mesh().triangles();
  const auto& P = tr_.plane();

  std::vector<short> pv(3);
  for (size_t i = 0; i < 3; ++i) {
    pv[i] = tr_.pv(pid, M[3*tid+i]);
  }
  size_t vi(INVALID_NUM), vj(INVALID_NUM), vk(INVALID_NUM);
  for (size_t i = 0; i < 3; ++i) {
    if ((pv[i]>pv[(i+1)%3] && pv[i]>pv[(i+2)%3])
        || (pv[i]<pv[(i+1)%3] && pv[i]<pv[(i+2)%3])) {
      vi = M[3*tid+i];
      vj = M[3*tid+(i+1)%3];
      vk = M[3*tid+(i+2)%3];

      if (pv[i]>pv[(i+1)%3]) {
        seq_pt_dir_[pt_id] = '+';
      } else {
        seq_pt_dir_[pt_id] = '-';
      }
      
      break;
    }
  }
  FASSERT(vi!=INVALID_NUM && vj!=INVALID_NUM && vk!=INVALID_NUM);

  // find left_edge_point, right_edge_point.
  size_t left_edge_point = find_edge_cv(vi, vj, pid);
  size_t right_edge_point = find_edge_cv(vi, vk, pid);
  FASSERT(left_edge_point != INVALID_NUM);
  FASSERT(right_edge_point != INVALID_NUM);
  seq_pt_bound_[pt_id] = std::make_pair(left_edge_point, right_edge_point);

  if (seq.size() == 0) return 0;

  pri_queue seq_map;

  size_t cnt = 0;
  for (auto& v : seq) {
    key_type k;
    cal_pt_point(&V[3*vi],&V[3*vj],&V[3*vk], &P[4*pid], &P[4*v.plane], k);
    keep_key_valid(k);
    seq_map.push(std::make_pair(k, cnt));
    ++cnt;
  }
  std::vector<pt_point> new_seq;
  std::vector<key_type> seq_key;
  new_seq.reserve(seq.size());
  seq_key.reserve(seq.size());
    

  while (!seq_map.empty()) {
    auto s = seq_map.top();
    seq_map.pop();

    new_seq.push_back(seq[s.second]);//info
    seq_key.push_back(s.first);//key
  }

    
  std::swap(new_seq, seq);


  // fine sorting
  bool is_ok = false;
  while (!is_ok) {
    is_ok = true;
    for (size_t si = 0; si+1 < seq.size(); ++si) {
      if (myc(seq_eps_).small_distance(seq_key[si], seq_key[si+1])) {
        int r = cmp_pt_point(
            &V[3*vi], &V[3*vj], &V[3*vk], &P[4*pid],
            &P[4*seq[si].plane], &P[4*seq[si+1].plane]);
        // r>0, si>si+1; r<0, si<si+1.
        if (r > 0) {
          std::swap(seq[si], seq[si+1]);
          std::swap(seq_key[si], seq_key[si+1]);
          is_ok = false;
        } else if (r == 0) {
          ufs_.set_union(seq[si].cvid, seq[si+1].cvid);//cut verts index
        }
      }
    }//for
  }//while

  key_type bd_v;
  bd_v.first = 0.0; bd_v.second = 1.0;
  if (myc(seq_eps_).small_distance(seq_key.front(),bd_v)) {
    int r = cmp_pt_point_to_const(
        &V[3*vi], &V[3*vj], &V[3*vk], &P[4*pid], &P[4*seq.front().plane], 0.0);
    if (r==0) {
      ufs_.set_union(seq.front().cvid, left_edge_point);
    } else if (r < 0) {
      std::cerr << "# [ ERROR ] this error cannot be happened theoretically (pt left).\n"
                << "            key: " << seq_key.front().first << " / "
                << seq_key.front().second << std::endl;
      return 1;
    }
  }
    
  bd_v.first = bd_v.second = 1.0;
  if (myc(seq_eps_).small_distance(seq_key.back(),bd_v)) {
    int r = cmp_pt_point_to_const(
        &V[3*vi], &V[3*vj], &V[3*vk], &P[4*pid], &P[4*seq.back().plane], 1.0);
    if (r==0) {
      ufs_.set_union(seq.back().cvid, right_edge_point);
    } else if (r > 0) {
      std::cerr << "# [ ERROR ] this error cannot be happened theoretically (pt right)."
                << std::endl;
      return 1;
    }
  }
    
  return 0;
}


template<typename T>
short key_cmp(
    const T& a1, const T& a2,
    const T& b1, const T& b2)
{
  return (sgn(a1*b2-a2*b1)*sgn(a2)*sgn(b2));
}


template<typename T>
class coord_cmp {
 public:
  typedef std::pair<Eigen::Matrix<T,3,1>, T> KT;
  typedef std::pair<KT, size_t> key_type;

  int cmp(const KT& a, const KT& b)
  {
    short r0 = key_cmp(a.first[0], a.second, b.first[0], b.second);
    if (r0 == 0) {
      short r1 = key_cmp(a.first[1], a.second, b.first[1], b.second);
      if (r1 == 0) {
        short r2 = key_cmp(a.first[2], a.second, b.first[2], b.second);
        return r2;
      } else {
        return r1;
      }
    } else {
      return r0;
    }
  }
  
  int cmp(const key_type& ka, const key_type& kb)
  {
    return cmp(ka.first, kb.first);
  }
  
  bool operator() (const key_type& ka, const key_type& kb)
  {
    return (cmp(ka,kb)<0);
  }
};

//////////////////////////////
template<typename T>
struct point_coord_type {
  Eigen::Matrix<T,3,1> nm;
  T dm;
};


int sequence_type::sort_tt_with_coord(
    const std::pair<size_t,size_t>& e, std::vector<tt_point>& seq)
{
  if (seq.size() == 0) {return 0;}
  const auto& V = tr_.mesh().verts();
  const auto& P = tr_.plane();

  MPF_VEC direc;
  get_edge_direc(&V[3*e.first], &V[3*e.second], direc);
  auto direc_type = choose_axis(direc);


  FASSERT(direc_type != 0);

  // key value and index of seq element
  // pri_queue seq_map;
  std::vector<key_map> seq_vec;
  seq_vec.reserve(seq.size());

  {
    size_t cnt = 0;
    for (auto& cv : seq) {
      point_coord_type<double> pc;
      FASSERT(cv.plane < P.size()/4);
      // double
      cal_ttp_point_coord(&V[3*e.first], &V[3*e.second], &P[4*cv.plane], pc.nm, pc.dm);
      seq_vec.push_back(std::make_pair(std::make_pair(pc.nm[abs(direc_type)-1],pc.dm), cnt));
      keep_key_valid(seq_vec.back().first);
      ++cnt;
    }
    FASSERT(cnt == seq.size());
    std::sort(seq_vec.begin(), seq_vec.end(), line_sort_myc<double>());    
  }


  if (direc_type < 0) {
    std::reverse(seq_vec.begin(), seq_vec.end());
  }


  // fine sorting
  bool is_ok = false;
  std::unordered_map<size_t, std::pair<mpf_class,mpf_class> > LVC; //line vert coordinate
  while (!is_ok) {
    is_ok = true;
    for (size_t si = 0; si+1 < seq.size(); ++si) {
      if (pc_myc(seq_eps_).small_distance(seq_vec[si].first,seq_vec[si+1].first)) {
        key_mpf ka;
        size_t id_a = seq_vec[si].second;
        {
          auto ita = LVC.find(id_a);
          if (ita == LVC.end()) {
            point_coord_type<mpf_class> co;
            cal_ttp_point_coord(&V[3*e.first], &V[3*e.second], &P[4*seq[id_a].plane], co.nm, co.dm);
            ka.first = co.nm[abs(direc_type)-1];
            ka.second = co.dm;
            LVC.insert(std::make_pair(id_a, ka));
          } else {
            ka = ita->second;
          }
        }

        key_mpf kb;
        size_t id_b = seq_vec[si+1].second;
        {
          auto itb = LVC.find(id_b);
          if (itb == LVC.end()) {
            point_coord_type<mpf_class> co;
            cal_ttp_point_coord(&V[3*e.first], &V[3*e.second], &P[4*seq[id_b].plane], co.nm, co.dm);
            kb.first = co.nm[abs(direc_type)-1];
            kb.second = co.dm;
            LVC.insert(std::make_pair(id_b, kb));
          } else {
            kb = itb->second;
          }
        }

        // compare, consider direction
        short r = key_cmp(ka, kb)*direc_type;
        if (r == 0) {
          ufs_.set_union(seq[id_a].cvid, seq[id_b].cvid);
          //std::cout << "-- merge " << id_a << ", " << id_b << std::endl;
        } else if (r > 0) {
          std::swap(seq_vec[si], seq_vec[si+1]);
          is_ok = false;
        }
      }//if
    }//for
  }//while


  // double
  key_type left, right;
  left.first = V[3*e.first+abs(direc_type)-1];
  left.second = 1.0;
  right.first = V[3*e.second+abs(direc_type)-1];
  right.second = 1.0;


  if (pc_myc(seq_eps_).small_distance(left, seq_vec.front().first)) {
    key_mpf k;
    size_t id = seq_vec.front().second;
    auto it = LVC.find(id);
    if (it == LVC.end()) {
      point_coord_type<mpf_class> co;
      cal_ttp_point_coord(&V[3*e.first], &V[3*e.second], &P[4*seq[id].plane], co.nm, co.dm);
      k.first = co.nm[abs(direc_type)-1];
      k.second = co.dm;
    } else {
      k = it->second;
    }
    key_mpf L;
    L.first = V[3*e.first+abs(direc_type)-1];
    L.second = 1.0;
    
    short r = key_cmp(L, k)*direc_type;
    if (r==0) {
      ufs_.set_union(seq[id].cvid, e.first+cv_.size());

    }
    FASSERT(r<=0);
  }

  if (pc_myc(seq_eps_).small_distance(right, seq_vec.back().first)) {
    key_mpf k;
    size_t id = seq_vec.back().second;
    auto it = LVC.find(id);
    if (it == LVC.end()) {
      point_coord_type<mpf_class> co;
      cal_ttp_point_coord(&V[3*e.first], &V[3*e.second], &P[4*seq[id].plane], co.nm, co.dm);
      k.first = co.nm[abs(direc_type)-1];
      k.second = co.dm;
    } else {
      k = it->second;
    }
    key_mpf R;
    R.first = V[3*e.second+abs(direc_type)-1];
    R.second = 1.0;
    short r = key_cmp(k, R)*direc_type;
    if (r == 0) {
      ufs_.set_union(seq[id].cvid, e.second+cv_.size());

    }
    FASSERT(r<=0);
  }

  {
    std::vector<tt_point> new_seq;
    new_seq.reserve(seq.size());
    for (auto& s : seq_vec) {
      new_seq.push_back(seq[s.second]);
    }
    std::swap(new_seq, seq);
  }
  
  return 0;
}

    
int sequence_type::sort_pp_with_coord(
    size_t pa, size_t pb, std::vector<pp_point>& seq)
{
  size_t tn = tr_.mesh().faces_num();
  auto& M = tr_.mesh().triangles();
  auto& V = tr_.mesh().verts();
  auto& P = tr_.plane();

  // define a starting point by choosing a axis-plane
  double sP[4];
  sP[0] = sP[1] = sP[2] = sP[3] = 0.0;
  {
    double max = 0.0;
    size_t k = 0;
    for (size_t i = 0; i < 3; ++i) {
      sP[i] = 1.0;
      double vol = fabs(cal_vol(&P[4*pa],&P[4*pb],sP));
      if (vol > max) { max = vol; k = i;}
      sP[i] = 0.0;
    }
    sP[k] = 1.0;
  }
  
  Eigen::Matrix<double,3,1> NSd, dir;
  double Nd;
  pp_precompute(&P[4*pa], &P[4*pb], sP, NSd, dir, Nd);

  MPF_VEC direc;
  get_plane2_direc(&P[4*pa], &P[4*pb], direc);
  auto direc_type = choose_axis(direc);
  FASSERT(direc_type != 0);  

  std::unordered_map<size_t,size_t> f2v;
    
  std::vector<key_map> seq_vec;

  {
    size_t cnt(0);
    for (auto& v : seq) {
      // MPF_VEC nm;
      // mpf_class dm;
      point_coord_type<double> pc;
      if (v.face >= tn) { // PPP
        cal_ppp_point_coord(NSd,dir,Nd,&P[4*(v.face-tn)], pc.nm, pc.dm);
      } else { // PPT, precondition, there exists only one intersection point between them.
        size_t fid = v.face;
        size_t ppv = INVALID_NUM;
        for (size_t k = 0; k < 3; ++k) {
          if (tr_.pv(pa, M[3*fid+k])==0 && tr_.pv(pb, M[3*fid+k])==0) {
            ppv = M[3*fid+k];
            f2v.insert(std::make_pair(fid, ppv));
            break;
          }
        }
        if (ppv == INVALID_NUM) {
          size_t vi(M[3*fid]), vj(M[3*fid+1]), vk(M[3*fid+2]);
          cal_ppt_point_coord(NSd,dir,Nd,&V[3*vi],&V[3*vj],&V[3*vk], pc.nm, pc.dm);
        } else {
          pc.nm << V[3*ppv], V[3*ppv+1], V[3*ppv+2];
          pc.dm = 1.0;
        }
      }
      seq_vec.push_back(std::make_pair(std::make_pair(pc.nm[abs(direc_type)-1],pc.dm), cnt));
      keep_key_valid(seq_vec.back().first);
      ++cnt;
    }
    std::sort(seq_vec.begin(), seq_vec.end(), line_sort_myc<double>());
  }


  if (direc_type < 0) {
    std::reverse(seq_vec.begin(), seq_vec.end());
  }

  Eigen::Matrix<mpf_class,3,1> gmp_NSd, gmp_dir;
  mpf_class gmp_Nd;
  pp_precompute(&P[4*pa], &P[4*pb], sP, gmp_NSd, gmp_dir, gmp_Nd);
  

  std::unordered_map<size_t, std::pair<mpf_class, mpf_class> > LVC; // line vert coordinates
  bool is_ok = false;
  
  while (!is_ok) {
    is_ok = true;

    for (size_t si = 0; si+1 < seq_vec.size(); ++si) {
      if (pc_myc(seq_eps_).small_distance(seq_vec[si].first, seq_vec[si+1].first)) {
        size_t id_a = seq_vec[si].second;
        key_mpf ka;
        
        {// ita
          auto ita = LVC.find(id_a);
          if (ita == LVC.end()) {
            point_coord_type<mpf_class> co;
            auto& lv = seq[id_a];
            if (lv.face >= tn) {
              cal_ppp_point_coord(gmp_NSd,gmp_dir,gmp_Nd,&P[4*(lv.face-tn)], co.nm, co.dm);
            } else {
              auto itv = f2v.find(lv.face);
              if (itv == f2v.end()) {
                size_t vi(M[3*lv.face]), vj(M[3*lv.face+1]), vk(M[3*lv.face+2]);
                cal_ppt_point_coord(gmp_NSd,gmp_dir,gmp_Nd,&V[3*vi],&V[3*vj],&V[3*vk], co.nm, co.dm);
              } else {
                co.nm << V[3*itv->second], V[3*itv->second+1], V[3*itv->second+2];
                co.dm = 1.0;
              }
            } // if
            ka.first = co.nm[abs(direc_type)-1];
            ka.second = co.dm;
            LVC.insert(std::make_pair(id_a, ka));
          } else {
            ka = ita->second;
          }
        }// ita

        key_mpf kb;        
        size_t id_b = seq_vec[si+1].second;
        
        {// itb
          auto itb = LVC.find(id_b);
          if (itb == LVC.end()) {
            point_coord_type<mpf_class> co;
            auto& lv = seq[id_b];
            if (lv.face >= tn) {
              cal_ppp_point_coord(&P[4*pa], &P[4*pb], &P[4*(lv.face-tn)], co.nm, co.dm);
            } else {
              auto itv = f2v.find(lv.face);
              if (itv == f2v.end()) {
                size_t vi(M[3*lv.face]), vj(M[3*lv.face+1]), vk(M[3*lv.face+2]);
                cal_ppt_point_coord(&P[4*pa], &P[4*pb], &V[3*vi], &V[3*vj], &V[3*vk], co.nm, co.dm);
              } else {
                co.nm << V[3*itv->second], V[3*itv->second+1], V[3*itv->second+2];
                co.dm = 1.0;
              }
            } // if
            kb.first = co.nm[abs(direc_type)-1];
            kb.second = co.dm;
            LVC.insert(std::make_pair(id_b, kb));
          } else {
            kb = itb->second;
          }
        }//itb

        // compare, consider direction
        short r = key_cmp(ka, kb)*direc_type;
        if (r == 0) {
          ufs_.set_union(seq[id_a].cvid, seq[id_b].cvid);
        } else if (r > 0) {
          std::swap(seq_vec[si], seq_vec[si+1]);
          is_ok = false;
        }
      }
    }
    
  }

  {
    std::vector<pp_point> new_seq;
    new_seq.reserve(seq.size());
    for (auto& s : seq_vec) {
      new_seq.push_back(seq[s.second]);
    }
    std::swap(new_seq, seq);
  }


  return 0;
}
    
int sequence_type::sort_pt_with_coord(
    size_t pt_id, size_t pid, size_t tid, std::vector<pt_point>& seq)
{
  const auto& V = tr_.mesh().verts();
  const auto& M = tr_.mesh().triangles();
  const auto& P = tr_.plane();

  MPF_VEC direc;
  get_pt_direc(&P[4*pid], &V[3*M[3*tid]], &V[3*M[3*tid+1]], &V[3*M[3*tid+2]], direc);
  auto direc_type = choose_axis(direc);

  FASSERT(direc_type != 0);

  std::vector<short> pv(3);
  for (size_t i = 0; i < 3; ++i) {
    pv[i] = tr_.pv(pid, M[3*tid+i]);
  }
  size_t vi(INVALID_NUM), vj(INVALID_NUM), vk(INVALID_NUM);
  for (size_t i = 0; i < 3; ++i) {
    if ((pv[i]>pv[(i+1)%3] && pv[i]>pv[(i+2)%3])
        || (pv[i]<pv[(i+1)%3] && pv[i]<pv[(i+2)%3])) {
      vi = M[3*tid+i];
      vj = M[3*tid+(i+1)%3];
      vk = M[3*tid+(i+2)%3];

      if (pv[i]>pv[(i+1)%3]) {
        seq_pt_dir_[pt_id] = '+';
      } else {
        seq_pt_dir_[pt_id] = '-';
        direc_type *= -1;
      }
      
      break;
    }
  }
  FASSERT(vi!=INVALID_NUM && vj!=INVALID_NUM && vk!=INVALID_NUM);

  // check whether the starting and ending cutted vertices are equal to vl or vr.
  // find left_edge_point, right_edge_point.
  size_t left_edge_point = find_edge_cv(vi, vj, pid);
  size_t right_edge_point = find_edge_cv(vi, vk, pid);
  FASSERT(left_edge_point != INVALID_NUM);
  FASSERT(right_edge_point != INVALID_NUM);
  seq_pt_bound_[pt_id] = std::make_pair(left_edge_point, right_edge_point);

  
  if (seq.size() == 0) return 0;

  std::vector<key_map> seq_vec;
  seq_vec.reserve(seq.size());

  {
    size_t cnt = 0;
    for (auto& v : seq) {
      point_coord_type<double> pc;
      cal_ppt_point_coord(&P[4*pid], &P[4*v.plane], &V[3*vi],&V[3*vj],&V[3*vk], pc.nm, pc.dm);
      seq_vec.push_back(std::make_pair(std::make_pair(pc.nm[abs(direc_type)-1],pc.dm), cnt));
      keep_key_valid(seq_vec.back().first);
      ++cnt;
    }
    std::sort(seq_vec.begin(), seq_vec.end(), line_sort_myc<double>());
  }


  if (direc_type < 0) {
    std::reverse(seq_vec.begin(), seq_vec.end());
  }

  std::unordered_map<size_t, std::pair<mpf_class,mpf_class> > LVC; //line vert coordinate  

  bool is_ok = false;
  while (!is_ok) {
    is_ok = true;

    for (size_t si = 0; si+1 < seq_vec.size(); ++si) {
      if (pc_myc(seq_eps_).small_distance(seq_vec[si].first, seq_vec[si+1].first)) {
        key_mpf ka, kb;

        size_t id_a = seq_vec[si].second;
        {
          auto ita = LVC.find(id_a);
          if (ita == LVC.end()) {
            point_coord_type<mpf_class> pc;
            // ppt is not the vertex of triangle
            cal_ppt_point_coord(
                &P[4*pid], &P[4*seq[id_a].plane],
                &V[3*vi], &V[3*vj], &V[3*vk], pc.nm, pc.dm);
            ka.first = pc.nm[abs(direc_type)-1];
            ka.second = pc.dm;
            LVC.insert(std::make_pair(id_a, ka));
          } else {
            ka = ita->second;
          }
        }

        size_t id_b = seq_vec[si+1].second;
        {
          auto itb = LVC.find(id_b);
          if (itb == LVC.end()) {
            point_coord_type<mpf_class> pc;
            // ppt is not the vertex of triangle
            cal_ppt_point_coord(
                &P[4*pid], &P[4*seq[id_b].plane],
                &V[3*vi], &V[3*vj], &V[3*vk], pc.nm, pc.dm);
            kb.first = pc.nm[abs(direc_type)-1];
            kb.second = pc.dm;
            LVC.insert(std::make_pair(id_b, kb));
          } else {
            kb = itb->second;
          }
        }

        // compare, consider direction
        short r = key_cmp(ka, kb)*direc_type;
        if (r == 0) {
          ufs_.set_union(seq[id_a].cvid, seq[id_b].cvid);
        } else if (r > 0) {
          std::swap(seq_vec[si], seq_vec[si+1]);
          is_ok = false;
        }
      }
    }
  }

  key_type L, R;
  {
    point_coord_type<double> pc;
    cal_ttp_point_coord(&V[3*vi], &V[3*vj], &P[4*pid], pc.nm, pc.dm);
    L.first = pc.nm[abs(direc_type)-1];
    L.second = pc.dm;
    keep_key_valid(L);
    ////
    cal_ttp_point_coord(&V[3*vi], &V[3*vk], &P[4*pid], pc.nm, pc.dm);
    R.first = pc.nm[abs(direc_type)-1];
    R.second = pc.dm;
    keep_key_valid(R);
  }

  if (pc_myc(seq_eps_).small_distance(L, seq_vec.front().first)) {

    key_mpf k, kL;
    {
      point_coord_type<mpf_class> pc;
      cal_ttp_point_coord(&V[3*vi], &V[3*vj], &P[4*pid], pc.nm, pc.dm);
      kL.first = pc.nm[abs(direc_type)-1];
      kL.second = pc.dm;
    }
    
    size_t id = seq_vec.front().second;
    {
      auto it = LVC.find(id);
      if (it == LVC.end()) {
        point_coord_type<mpf_class> pc;
        cal_ppt_point_coord(&P[4*pid], &P[4*seq[id].plane],
                            &V[3*vi],&V[3*vj],&V[3*vk], pc.nm, pc.dm);
        k.first = pc.nm[abs(direc_type)-1];
        k.second = pc.dm;
        LVC.insert(std::make_pair(id, k));
      } else {
        k = it->second;
      }
    }

    short r = key_cmp(kL, k)*direc_type;
    FASSERT(r<=0);    
    if (r == 0) {
      ufs_.set_union(seq[id].cvid, left_edge_point);
    }
  }

  if (pc_myc(seq_eps_).small_distance(R, seq_vec.back().first)) {
    key_mpf k, kR;
    {
      point_coord_type<mpf_class> pc;
      cal_ttp_point_coord(&V[3*vi], &V[3*vk], &P[4*pid], pc.nm, pc.dm);
      kR.first = pc.nm[abs(direc_type)-1];
      kR.second = pc.dm;
    }
    
    size_t id = seq_vec.back().second;
    {
      auto it = LVC.find(id);
      if (it == LVC.end()) {
        point_coord_type<mpf_class> pc;
        cal_ppt_point_coord(&P[4*pid], &P[4*seq[id].plane],
                            &V[3*vi],&V[3*vj],&V[3*vk], pc.nm, pc.dm);
        k.first = pc.nm[abs(direc_type)-1];
        k.second = pc.dm;
        LVC.insert(std::make_pair(id, k));
      } else {
        k = it->second;
      }
    }

    short r = key_cmp(k, kR)*direc_type;
    FASSERT(r<=0);
    if (r == 0) {
      ufs_.set_union(seq[id].cvid, right_edge_point);
    }
    
  }

  {  
    std::vector<pt_point> new_seq;
    new_seq.reserve(seq.size());  
    for (auto& s : seq_vec) {
      new_seq.push_back(seq[s.second]);
    }
    std::swap(new_seq, seq);
  }

  return 0;
}

int sequence_type::move_cv(
    std::vector<double>& V,
    const std::vector<size_t>& vs,
    short d, bool is_fix_bd,
    const Eigen::Matrix<double, 3, 1>& D,
    char& is_move) const
{
  FASSERT(fabs(D.norm()-1.0) < 1e-4);
  FASSERT(d!=0);


  for (size_t i = 1; i+1 < vs.size(); ++i) {
    Eigen::Map<Eigen::Matrix<double, 3, 1> > p(&V[3*vs[i]], 3, 1);
    Eigen::Map<Eigen::Matrix<double, 3, 1> > p_p(&V[3*vs[i-1]], 3, 1);

    double eps = EPS_MERGE*(my_abs(p)+my_abs(p_p));    

    while ((p-p_p).dot(D) < eps) {
      p += EPS_MV*(my_abs(p)+my_abs(p_p))*D;
      is_move = 1;
    }

  }

  FASSERT(vs.size()>=2);
  for (size_t i = vs.size()-2; i > 0; --i) {
    Eigen::Map<Eigen::Matrix<double, 3, 1> > p(&V[3*vs[i]], 3, 1);
    Eigen::Map<Eigen::Matrix<double, 3, 1> > n_p(&V[3*vs[i+1]], 3, 1);

    double eps = EPS_MERGE*(my_abs(p)+my_abs(n_p));
    
    while ( (n_p-p).dot(D) < eps ) {
      p -= EPS_MV*(my_abs(p)+my_abs(n_p))*D;
      is_move = 1;
    }

  }


  if (!is_fix_bd) {
    Eigen::Map<Eigen::Matrix<double, 3, 1> > pA(&V[3*vs[0]], 3, 1);
    Eigen::Map<Eigen::Matrix<double, 3, 1> > pB(&V[3*vs[1]], 3, 1);

    while ((pB-pA).dot(D) <= EPS_MERGE*(my_abs(pA)+my_abs(pB))) {
      pA -= EPS_MV*(my_abs(pA)+my_abs(pB))*D;
      is_move = 1;
    }

    Eigen::Map<Eigen::Matrix<double, 3, 1> > pC(&V[3*vs[vs.size()-2]], 3, 1);
    Eigen::Map<Eigen::Matrix<double, 3, 1> > pD(&V[3*vs.back()], 3, 1);

    while ((pD-pC).dot(D) <= EPS_MERGE*(my_abs(pC)+my_abs(pD))) {
      pD += EPS_MV*(my_abs(pC)+my_abs(pD))*D;
      is_move = 1;
    }
  }
  
  return 0;
}

int sequence_type::move_cut_verts(std::vector<double>& new_CV) const
{
  for (size_t i = 0; i < 5; ++i) {
    char is_move(0);
    move_cv_on_line(new_CV, is_move);
    if (is_move==0) {break;}
  }
  return 0;
}

int sequence_type::move_cv_on_line(std::vector<double>& new_CV, char& is_move) const
{
  FASSERT(new_CV.size() == 3*pure_cv_num());

  const auto& V = tr_.mesh().verts();
  const auto& M = tr_.mesh().triangles();
  const auto& P = tr_.plane();

  size_t en = e2f().size();
  // get tt seq, each edge is not related
#pragma omp parallel for
  for (size_t ei = 0; ei < en; ++ei) {
    const auto& stt = seq_tt(ei);
    if (stt.size() == 0) {continue;}
    auto& e = ev(ei);
    
    std::vector<size_t> s;
    s.reserve(stt.size()+2);
    
    s.push_back(e.first);
    size_t pre_id = s.front();
    for (auto& v : stt) {
      size_t vid = cv_[v.cvid].id();
      if (pre_id == vid) {continue;}
      s.push_back(vid);
      pre_id = vid;
    }
    if (pre_id != e.second) {
      s.push_back(e.second);
    }
    
    MPF_VEC direc;
    get_edge_direc(&V[3*s.front()], &V[3*s.back()], direc);
    Eigen::Matrix<double, 3, 1> D;
    D[0]=direc[0].get_d();
    D[1]=direc[1].get_d();
    D[2]=direc[2].get_d();
    D/=D.norm();
    move_cv(new_CV, s, choose_axis(direc), true, D, is_move);
  }

  std::vector<std::pair<size_t,size_t> > pt_ids(pt_map_.size());
  for (auto& pt : pt_map_) {
    pt_ids[pt.second] = pt.first;
  }
  // get pt seq
#pragma omp parallel for
  //  for (auto& pt : pt_map_) {
  for (size_t pti = 0; pti < pt_ids.size(); ++pti) {
    //size_t pti = pt.second;
    const auto& seq = seq_pt_[pti];
    size_t pre_id = cv_[seq_pt_bound_[pti].first].id();
    FASSERT(pre_id != INVALID_NUM);
    std::vector<size_t> s;
    s.reserve(seq.size()+2);
    s.push_back(pre_id);
    for (auto& v : seq) {
      size_t vid = cv_[v.cvid].id();
      if (vid == pre_id) continue;
      s.push_back(vid);
      pre_id = vid;
    }
    size_t right_id = cv_[seq_pt_bound_[pti].second].id();
    FASSERT(right_id != INVALID_NUM);
    if (pre_id != right_id) {
      s.push_back(right_id);
    }
    //size_t pid = pt.first.first;
    //size_t tid = pt.first.second;
    size_t pid = pt_ids[pti].first;
    size_t tid = pt_ids[pti].second;
    FASSERT(pid >= tr_.tn());
    FASSERT(tid < tr_.tn());
    pid -= tr_.tn();
    MPF_VEC direc;
    get_pt_direc(&P[4*pid], &V[3*M[3*tid]], &V[3*M[3*tid+1]], &V[3*M[3*tid+2]], direc);
    Eigen::Matrix<double, 3, 1> D;
    D[0] = direc[0].get_d();
    D[1] = direc[1].get_d();
    D[2] = direc[2].get_d();
    D/=D.norm();
    short direc_type = choose_axis(direc);
    if (seq_pt_dir_[pti] == '-') {
      direc_type *= -1;
      D *= -1.0;
    }
    move_cv(new_CV, s, direc_type, true, D, is_move);
  }


  // get pp seq
  for (auto& pp : pp_map_) {
    size_t si = pp.second;
    const auto& seq = seq_pp_[si];

    if (seq.size() == 0){continue;}
    std::vector<size_t> s;
    s.reserve(seq.size());
    
    size_t pre_id = INVALID_NUM;
    for (size_t i = 0; i < seq.size(); ++i) {
      size_t vid = cv_[seq[i].cvid].id();
      if (vid != pre_id) {
        s.push_back(vid);
      }
      pre_id = vid;
    }
    size_t pa = pp.first.first;
    size_t pb = pp.first.second;
    FASSERT(pa >= tr_.tn());
    FASSERT(pb >= tr_.tn());
    pa -= tr_.tn();
    pb -= tr_.tn();

    MPF_VEC direc;
    get_plane2_direc(&P[4*pa], &P[4*pb], direc);
    Eigen::Matrix<double, 3, 1> D;
    D[0] = direc[0].get_d();
    D[1] = direc[1].get_d();
    D[2] = direc[2].get_d();
    D /= D.norm();
    move_cv(new_CV, s, choose_axis(direc), false, D, is_move);
  }

  return 0;
}


} // grid
}// fxz
