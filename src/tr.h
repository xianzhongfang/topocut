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


#ifndef TR_FXZ_H
#define TR_FXZ_H

#include <map>
#include <iostream>
#include "my_gmp.h"

namespace fxz {
  namespace grid {

    typedef std::map<std::pair<size_t,size_t>, std::vector<size_t> > e2f_map_type;

    class base_grid_mesh
    {
   public:
   base_grid_mesh(const std::vector<double>& V, const std::vector<size_t>& M)
       : V_(V), M_(M) {  vn_=V_.size()/3; tn_=M_.size()/3; }

      const std::vector<size_t>& triangles() const { return M_; }
      const std::vector<double>& verts() const { return V_; }

      size_t tv(size_t fid, size_t id) const { assert(fid<tn_ && id<3); return M_[3*fid+id]; }
      const double* vert_coord(size_t vid) const { assert(vid<vn_); return &V_[3*vid]; }
      const size_t* tv(size_t fid) const { return &M_[3*fid]; }

      size_t tn() const { return tn_; }
      size_t vn() const { return vn_; }
      size_t faces_num() const { return tn_; }
      size_t verts_num() const { return vn_; }
      
   private:
      const std::vector<double>& V_;
      const std::vector<size_t>& M_;
      size_t vn_;
      size_t tn_;
    };


    enum cut_vert_type: char {
      VERT=1, EDGE=2, TRIANGLE=3, GRID=4, INNER=5
          };

    struct base_cut_vert {
     public:
      base_cut_vert() { flag_=0; }

      void set_on_triangle(size_t tid) {
        flag_ = ((tid<<3)|(flag_&1)|2);
      }
  
      void set_on_edge(size_t eid) {
        flag_ = ((eid<<3)|(flag_&1)|4);
      }
  
      void set_on_vert(size_t vid) {
        flag_ = ((vid<<3)|(flag_&1)|6);
      }

      void set_on_grid() { flag_ |= 1; }
  
      bool is_set() const { return (flag_==0)?false:true; }
  
      bool is_grid() const { return ((flag_&1)==1)?true:false; }
      bool is_triangle() const { return ((flag_&6)==2)?true:false; }//01*
      bool is_edge() const { return ((flag_&6)==4)?true:false; } //10*
      bool is_vert() const { return ((flag_&6)==6)?true:false; } //11*
  
      cut_vert_type cv_type() const {
        if (is_vert()) return cut_vert_type::VERT;
        else if (is_edge()) return cut_vert_type::EDGE;
        else if (is_triangle()) return cut_vert_type::TRIANGLE;
        return cut_vert_type::GRID;
      }
  
      size_t ele_id() const { return (flag_>>3); }

     private:
      size_t flag_; // ele_id + bound_flag (2bit, such as vert,edge,triangle) + grid_flag (1bit)
    };

    class cut_vert : public base_cut_vert {
   public:
      cut_vert(size_t fa, size_t fb, size_t fc) // T=>(fa<tn), P=>(fa>=tn).
      { f_[0]=fa; f_[1]=fb; f_[2]=fc; id_ = INVALID_NUM; }

      size_t f0() const { return f_[0]; }
      size_t f1() const { return f_[1]; }
      size_t f2() const { return f_[2]; }
      size_t f(size_t i) const {return f_[i]; }
      size_t& id() { return id_; }
      size_t id() const { return id_; }
  
   private:
      size_t f_[3]; // faces id
      size_t id_;   // cut vert id, there are repeated ones.
    };

    struct cut_face {
      std::vector<size_t> v;
    };


    class table_relation {
   public:
   table_relation(const base_grid_mesh& bm, const std::vector<double>& P)
       : bm_(bm), P_(P)
      {
        pn_ = P.size()/4;
        vn_ = bm.verts_num();
        tn_ = bm.faces_num();
        

        T_pp_.resize(pn_*(pn_-1)/2, 0);
        T_pv_.resize(pn_*vn_, 0);


#pragma omp parallel for
        for (size_t i = 0; i < pn_; ++i) {
          for (size_t j = i+1; j < pn_; ++j) {
            set_pp(i, j);
          }
        }

#pragma omp parallel for        
        for (size_t i = 0; i < pn_; ++i) {
          for (size_t j = 0; j < vn_; ++j) {
            set_pv(i, j);
          }
        }

        tuple_pt_.resize(pn_);

#pragma omp parallel for        
        for (size_t i = 0; i < pn_; ++i) {
          std::list<size_t> ts;
          for (size_t j = 0; j < tn_; ++j) {
            auto s = pt(i, j);
            if (s > 0) ts.push_back(j);
          }
          tuple_pt_[i].resize(ts.size());
          std::copy(ts.begin(), ts.end(), tuple_pt_[i].begin());
        }

      }

      size_t tn() const { return tn_; }
      size_t vn() const { return vn_; }
      size_t pn() const { return pn_; }
      const base_grid_mesh& mesh() const { return bm_; }
      const std::vector<double>& plane() const { return P_; }

      // triangle
      uint8_t pt(size_t p, size_t tid) const {
        if (tid > bm_.triangles().size()/3) { return 5; } //invalid
        size_t va = bm_.triangles()[3*tid];
        size_t vb = bm_.triangles()[3*tid+1];
        size_t vc = bm_.triangles()[3*tid+2];
        return pt(p, va, vb, vc);
      }

      uint8_t pt(size_t p, size_t va, size_t vb, size_t vc) const {
        short sa = pv(p, va);
        short sb = pv(p, vb);
        short sc = pv(p, vc);
        if (abs(sa+sb+sc)==2) return 1;
        else if (sa==sb && sb==sc && sc!=0) return 0;
        else if (sa==sb && sb==sc && sc==0) return 3;
        return 2;
      }



      uint8_t pe(size_t p, size_t va, size_t vb) const {
        short a = pv(p, va);
        short b = pv(p, vb);
        if (a==b && b!=0) { return 0; }
        else if (a==b && b==0) { return 2; }
        return 1;
      }

      void set_value(uint8_t val, size_t id, std::vector<uint16_t>& T) {
        auto& s = T[id>>3];
        uint16_t mask = 65535^(3<<(2*(id&7)));
        s = (s&mask)|(val<<(2*(id&7)));
      }

      uint8_t get_value(size_t id, const std::vector<uint16_t>& T) const {
        return (((T[id>>3]>>((id&7)*2)))&3);
      }

      size_t pv_id(size_t p, size_t v) const
      { return p*vn_+v; }


      void set_pv(size_t p, size_t v) {
        T_pv_[pv_id(p,v)] = cal_pv(&P_[4*p], &(bm_.verts()[3*v]));
      }
          

      size_t pp_id(size_t pa, size_t pb) const
      {
        assert(pa != pb);
        if (pa > pb) { std::swap(pa, pb); }
        return (pn_*pa-(1+pa)*pa/2+(pb-pa)-1);
      }

      void set_pp(size_t pa, size_t pb) {
        if (pa != pb) {
          T_pp_[pp_id(pa,pb)] = cal_pp(&P_[4*pa], &P_[4*pb]);
        } else {
          std::cerr << "# [ ERROR ] pa == pb." << std::endl;
        }
      }

      
      short pv(size_t p, size_t v) const {
        return T_pv_[pv_id(p,v)];
      }

      
      uint8_t pp(size_t pa, size_t pb) const {
        return T_pp_[pp_id(pa, pb)];
      }

      const auto& tuple_pt() const { return tuple_pt_; }

   private:
      const base_grid_mesh& bm_;
      const std::vector<double>& P_;
      size_t pn_, vn_, tn_;

      std::vector<short> T_pv_;
      std::vector<uint8_t> T_pp_;

      std::vector<std::vector<size_t> > tuple_pt_;
    };

  } // grid
} // fxz

#endif
