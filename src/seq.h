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


#ifndef SEQ_FXZ_H
#define SEQ_FXZ_H


#include <unordered_map>
#include <unordered_set>

#include "common.h"
#include "tr.h"
#include "union_find_set.h"

namespace fxz {
  namespace grid{

  enum class cycle_sort_type : short {
    SIMPLE=0, FULL=1
  };

  enum class line_sort_type : short {
    NONE=0,LINE_COORD=1, POINT_COORD=2,
  };

  class sequence_type {
   public:
   sequence_type(
       const table_relation& tr,
       const e2f_map_type& e2f,
       const std::vector<e2f_map_type::iterator>& eid_map,
       const double seq_eps,
       std::vector<cut_vert>& cut_verts,
       std::vector<double>& CV,
       line_sort_type lst=line_sort_type::POINT_COORD)
       : tr_(tr), e2f_(e2f), eid_map_(eid_map),
         seq_eps_(seq_eps), cv_(cut_verts), CV_(CV), lst_(lst)
      {
        switch (lst_) {
          case line_sort_type::LINE_COORD:
            //std::cout << "-- line sort type: line coordinates." << std::endl;
            break;
          case line_sort_type::POINT_COORD:
            //std::cout << "-- line sort type: point coordinates." << std::endl;
            break;
          default:
            std::cerr << "# [ ERROR ] there is no such line sort type." << std::endl;
            exit(1);
        }
        build_v2f();
      }

      int build_v2f();

      int run();
      int reset_cutvert_type();

      const std::vector<double>& get_cut_verts() ;
      const std::vector<double>& cut_verts_coord() const { return CV_; }
      std::vector<double>& cut_verts_coord() { return CV_; }    

      const std::vector<char>& get_cut_verts_type();
      const std::vector<char>& cut_verts_type() const { return CV_type_; }

      // unique id -> cut vert id.
      // here unique id is the vert id of input mesh M,
      // will return an id out of scope of cut verts.
      const std::vector<size_t>& get_cut_verts_id();
      const std::vector<size_t>& cut_verts_id() const { return CV_id_; }
      
      int cal_cv_coord(const cut_vert& v, double* coord) const;
      int get_cut_edges(
          std::vector<char>& cv_type,
          std::vector<size_t>& edge_CE,
          std::vector<size_t>& tri_CE,
          std::vector<size_t>& grid_in_CE,
          std::vector<size_t>& grid_out_CE) const;

    template<typename T>
    int get_cut_vert_coord(
        const cut_vert& x,
        std::pair<Eigen::Matrix<T,3,1>,T>& c) const
    {
      std::vector<short> t_id, p_id;
      t_id.reserve(3);
      p_id.reserve(3);

      for (size_t i = 0; i < 3; ++i) {
        if (x.f(i) < tr_.tn()) {
          t_id.push_back(i);
        } else {
          p_id.push_back(i);
        }
      }

      FASSERT(t_id.size()<3);
      if (t_id.size()==2) { // ttp
        std::unordered_set<size_t> fa_v;
        for (size_t i = 0; i < 3; ++i) {
          size_t vid = tr_.mesh().tv(x.f(t_id[0]), i); 
          fa_v.insert(vid);
        }
        std::vector<size_t> cm_v;
        cm_v.reserve(2);
        for (size_t i = 0; i < 3; ++i) {
          size_t vid = tr_.mesh().tv(x.f(t_id[1]), i);
          if (fa_v.count(vid) > 0) {
            cm_v.push_back(vid);
          }
        }
        size_t p = x.f(p_id[0]) - tr_.tn();
        FASSERT(cm_v.size() == 2);
        cal_ttp_point_coord(
            &tr_.mesh().verts()[3*cm_v[0]],
            &tr_.mesh().verts()[3*cm_v[1]],
            &tr_.plane()[4*p],
            c.first, c.second);
    
      } else if (t_id.size()==1) { // ppt

        size_t pa = x.f(p_id[0]) - tr_.tn();
        size_t pb = x.f(p_id[1]) - tr_.tn();
        bool is_vert = false;
        std::vector<size_t> v3(3);
        for (size_t i = 0; i < 3; ++i) {
          v3[i] = tr_.mesh().tv(x.f(t_id[0]), i);

          if (tr_.pv(pa, v3[i])==0 && tr_.pv(pb, v3[i])==0) {
            is_vert = true;
            auto vert = tr_.mesh().vert_coord(v3[i]);
            c.first << vert[0], vert[1], vert[2];
            c.second = 1.0;
            break;
          }
        }

        if (!is_vert) {
          cal_ppt_point_coord(
              &tr_.plane()[4*pa],
              &tr_.plane()[4*pb],
              &tr_.mesh().verts()[3*v3[0]],
              &tr_.mesh().verts()[3*v3[1]],
              &tr_.mesh().verts()[3*v3[2]],
              c.first, c.second);
        }
    
      } else { // ppp
    
        size_t pa = x.f(p_id[0]) - tr_.tn();
        size_t pb = x.f(p_id[1]) - tr_.tn();
        size_t pc = x.f(p_id[2]) - tr_.tn();
        cal_ppp_point_coord(
            &tr_.plane()[4*pa],
            &tr_.plane()[4*pb],
            &tr_.plane()[4*pc],
            c.first, c.second);
      }
  
      return 0;
    }


      struct tt_point {
       public:
       tt_point(size_t a, size_t b) : cvid(a), plane(b) { }
        size_t cvid, plane;
      };

      struct pt_point {
       public:
       pt_point(size_t a, size_t b) : cvid(a), plane(b) { }
        size_t cvid, plane;
      };

      struct pp_point {
       public:
       pp_point(size_t a, size_t b) : cvid(a), face(b) { }
        size_t cvid, face; // triangle or plane
      };

      const table_relation& tr() const { return tr_; }

      int classify();
      void add_pp(size_t cvid, size_t a, size_t b, size_t c);
      void add_pt(size_t cvid, size_t a, size_t b, size_t c);
      void add_tt(size_t cvid, size_t eid, size_t fid);
  
      int sort_tt();


      int sort_pp();

      int keep_pp_inner();
      size_t find_third_vert(size_t va, size_t vb, const size_t* M);
      uint8_t passing_triangle(size_t pa, size_t pb, size_t fid);
      uint8_t passing_edge(size_t eid, size_t pa, size_t pb, const std::vector<size_t>& fs);

      

      typedef Eigen::Matrix<mpf_class,3,1> MPF_VEC;
  
      uint8_t passing_vert(size_t vid, size_t pa, size_t pb);
      void cal_pt_outer_direc(
          size_t vid, size_t pid,
          const std::vector<size_t>& fs,
          std::vector<MPF_VEC>& dt) const;

      // > 2 adjacent triangles
      int cycle_sort_triangles(
          size_t vid, size_t pid,
          const std::vector<MPF_VEC>& dt,
          const std::unordered_map<size_t,size_t>& fs2id,
          std::vector<size_t>& fs,
          std::vector<char>& cp, std::vector<short>& dir,
          cycle_sort_type cst=cycle_sort_type::FULL) const;
      int sort_equal_triangles(
          size_t vid, size_t pid,
          std::vector<size_t>& fs,
          std::vector<char>& cp) const;
      
      // in the plane and vert
      int sort_equal_triangles(
          size_t vid, size_t pid,
          std::vector<size_t>& fs) const;

      int plane_triangle_edge_direc(
          const std::unordered_map<size_t,size_t>& fs2id,
          const std::vector<MPF_VEC>& dt,
          const std::vector<double>& V,
          const std::vector<size_t>& M,          
          const MPF_VEC& NP, size_t f) const;
      
      
      int cycle_sort_verts(
          size_t vid, size_t pid, const std::vector<size_t>& bd_fs,
          std::vector<std::vector<size_t> >& bd_fans);
      bool is_bound_fan(
          size_t vid, const std::vector<std::vector<size_t> >& bd_fans,
          size_t pid,
          size_t fa, size_t fb);
      uint8_t passing_cycle_fan(
          size_t vid,
          const std::vector<std::vector<size_t> >& bd_fans,
          size_t pa, size_t pb,
          const std::vector<MPF_VEC>& dt,
          std::unordered_map<size_t,size_t>& fs2id,
          const std::vector<size_t>& cut_fs,
          const std::vector<char>& cp,
          const std::vector<short>& dir);


      // 1: inner, 0: bound/outer
      // fan: <di, dj>, pass line: dk*ks
      int judge_fan(
          const MPF_VEC& N, const MPF_VEC& di, const MPF_VEC& dj,
          const MPF_VEC& dk, int ks) const;


      int judge_full_fan(
          const MPF_VEC& N, const MPF_VEC& di,
          const MPF_VEC& dk, int ks) const;
  
      int keep_pp_inner(size_t pa, size_t pb, size_t seq_id);
      size_t find_edge_cv(size_t va, size_t vb, size_t pid);
      int sort_pt();


      //// interface
      int sort_tt_api(
          const std::pair<size_t, size_t>& e,
          std::vector<tt_point>& seq)
      {
        switch(lst_) {
          case line_sort_type::LINE_COORD:
            return sort_tt_with_lc(e, seq);
          case line_sort_type::POINT_COORD:
            return sort_tt_with_coord(e, seq);
          default:
            std::cerr << "# [ ERROR ] the sort type is not right!" << std::endl;
            exit(1);
        }
        return 1;
      }
      
      int sort_pp_api(size_t pa, size_t pb, std::vector<pp_point>& seq)
      {
        switch(lst_) {
          case line_sort_type::LINE_COORD:
            return sort_pp_with_lc(pa, pb, seq);
          case line_sort_type::POINT_COORD:
            return sort_pp_with_coord(pa, pb, seq);
          default:
            std::cerr << "# [ ERROR ] the sort type is not right!" << std::endl;
            exit(1);
        }
        return 1;
      }
      
      int sort_pt_api(size_t pt_id, size_t pid, size_t tid, std::vector<pt_point>& seq)
      {
        switch(lst_) {
          case line_sort_type::LINE_COORD:
            return sort_pt_with_lc(pt_id, pid, tid, seq);
          case line_sort_type::POINT_COORD:
            return sort_pt_with_coord(pt_id, pid, tid, seq);
          default:
            std::cerr << "# [ ERROR ] the sort type is not right!" << std::endl;
            exit(1);
        }
        return 1;
      }


      //////////
      // sort with parametric coordinates
      int sort_tt_with_lc(
          const std::pair<size_t, size_t>& e,
          std::vector<tt_point>& seq);
      int sort_pp_with_lc(size_t pa, size_t pb, std::vector<pp_point>& seq);
      int sort_pt_with_lc(size_t pt_id, size_t pid, size_t tid, std::vector<pt_point>& seq);
      //////////

      //////////
      // sort with one component of coordinates
      int sort_tt_with_coord(const std::pair<size_t,size_t>& e, std::vector<tt_point>& seq);
      int sort_pp_with_coord(size_t pa, size_t pb, std::vector<pp_point>& seq);
      int sort_pt_with_coord(size_t pt_id, size_t pid, size_t tid, std::vector<pt_point>& seq);
      //////////

    int move_cv_on_line(std::vector<double>& new_CV, char& is_move) const;
    int move_cut_verts(std::vector<double>& new_CV) const;

    int move_cv(
        std::vector<double>& V,
        const std::vector<size_t>& vs,
        short d, bool is_fix_bd, const Eigen::Matrix<double, 3, 1>& D,
        char& is_move) const;
    
      const std::map<std::pair<size_t,size_t>, size_t>& edge_map_to_id() const { return e2id_; }
      size_t edge_id(const std::pair<size_t,size_t>& ev) const
      {
        auto it = e2id_.find(ev);
        if (it != e2id_.end()) {
          return it->second;
        }
        return INVALID_NUM;
      }
      
      const std::pair<size_t,size_t>& edge_v(size_t eid) const { return eid_map_[eid]->first; }

      const std::vector<tt_point>& seq_tt(size_t eid) const { return seq_tt_[eid]; }
      const std::vector<pt_point>& seq_pt(size_t id) const { return seq_pt_[id]; }
      const std::vector<pp_point>& seq_pp(size_t id) const { return seq_pp_[id]; }
      const std::vector<uint8_t>& seq_pp_s(size_t id) const { return pp_s_[id]; }

      const std::map<std::pair<size_t,size_t>,size_t>& pt_map() const { return pt_map_; }
      const std::pair<size_t,size_t>& pt_bd(size_t id) const { return seq_pt_bound_[id]; }

      const std::map<std::pair<size_t,size_t>,size_t>& pp_map() const { return pp_map_; }

      const cut_vert& cv(size_t cvid) const { FASSERT(cvid<cv_.size()); return cv_[cvid]; }
      const e2f_map_type& e2f() const { return e2f_; }
      const std::pair<size_t,size_t>& ev(size_t eid) const { return eid_map_[eid]->first; }
      const std::vector<size_t>& edge_fs(size_t eid) const { return eid_map_[eid]->second; }
      const std::vector<double>& CV() const { return CV_; }
      const char pt_dir(size_t id) const { return seq_pt_dir_[id]; }

      const size_t ori_cv_num() const { return cv_.size(); }
      const size_t pure_cv_num() const { return pure_cv_num_; }

    const std::vector<size_t>& v2f(size_t vid) const { FASSERT(vid<tr_.vn());return v2f_[vid]; }

   private:
      const table_relation& tr_;
      const e2f_map_type& e2f_;
      const std::vector<e2f_map_type::iterator>& eid_map_;

      std::vector<std::vector<size_t> > v2f_;

      double seq_eps_;
  
      // will be modified after sorting
      std::vector<cut_vert>& cv_;

      size_t tn_;
      std::map<std::pair<size_t,size_t>, size_t> e2id_;
      std::vector<std::vector<tt_point> > seq_tt_;

      // pp lines
      std::map<std::pair<size_t,size_t>,size_t> pp_map_;
      std::vector<std::vector<pp_point> > seq_pp_;
      std::vector<std::vector<uint8_t> > pp_s_;

      // pt lines, <pid+tn, tid>
      std::map<std::pair<size_t,size_t>,size_t> pt_map_;
      std::vector<std::vector<pt_point> > seq_pt_;
      std::vector<std::pair<size_t,size_t> > seq_pt_bound_;
      std::vector<char> seq_pt_dir_;//'+':the same;'-':opposite of the standard dir (np x nt)

      union_find_set ufs_;
      size_t pure_cv_num_;


      std::vector<double>& CV_;
      std::vector<char> CV_type_;
      std::vector<size_t> CV_id_;

      line_sort_type lst_;

    };

  } // grid
} // fxz

#endif
