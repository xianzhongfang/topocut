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


#ifndef GMP_FXZ_H
#define GMP_FXZ_H

#include <stdio.h>
#include <gmp.h>
#include <gmpxx.h>

#include "common.h"

namespace fxz {
namespace grid {
    
short cal_pv(const double* p, const double* v);
uint8_t cal_pp(const double* pa, const double* pb);

int vol_sign_with_norms(const double* na, const double* nb, const double* nc);

bool cal_ppt_intersection(
    const double* pa, const double* pb,
    const double* va, const double* vb, const double* vc);

void cal_ttp_point( // ttp
    const double* va, const double* vb, const double* p,
    std::pair<double, double>& k);

void cal_pt_point( // ptp
    const double* va, const double* vb, const double* vc,
    const double* p, const double* px, std::pair<double,double>& k);

// r>0, p1>p2; r<0, p1<p2.
int cmp_tt_point(
    const double* va, const double* vb, const double* p1, const double* p2);

int cmp_tt_point_to_const(
    const double* va, const double* vb, const double* p, const double kc);

int cmp_pt_point(
    const double* va, const double* vb, const double* vc,
    const double* p, const double* pa, const double* pb);

int cmp_pt_point_to_const(
    const double* va, const double* vb, const double* vc,
    const double* p, const double* px, const double kc);

template<typename T>
void N_star(const T& N, T& NS)
{
  NS(0,0) = N(1,1)*N(2,2)-N(1,2)*N(2,1);
      
  NS(0,1) = -(N(0,1)*N(2,2)-N(0,2)*N(2,1));
      
  NS(0,2) = (N(0,1)*N(1,2)-N(0,2)*N(1,1));
      
  NS(1,0) = -(N(1,0)*N(2,2)-N(1,2)*N(2,0));
      
  NS(1,1) = (N(0,0)*N(2,2)-N(0,2)*N(2,0));
      
  NS(1,2) = -(N(0,0)*N(1,2)-N(0,2)*N(1,0));
      
  NS(2,0) = (N(1,0)*N(2,1)-N(1,1)*N(2,0));
      
  NS(2,1) = -(N(0,0)*N(2,1)-N(0,1)*N(2,0));
      
  NS(2,2) = (N(0,0)*N(1,1)-N(0,1)*N(1,0));
}

template<typename T>
void pp_precompute(
    const double* pa, const double* pb, const double* pc,
    Eigen::Matrix<T,3,1>& NSd, Eigen::Matrix<T,3,1>& dir, T& Nd)
{
  Eigen::Matrix<T, 3, 3> N, NS;
  N << pa[0], pa[1], pa[2],
      pb[0], pb[1], pb[2],
      pc[0], pc[1], pc[2];
      
  N_star(N, NS);
  
  Eigen::Matrix<T, 3, 1> D;
  D << pa[3], pb[3], pc[3];
  
  NSd = NS*D;
  dir = N.row(0).cross(N.row(1));
  Nd = dir.dot(N.row(2));
}

template<typename T>
void cal_ppp_point(
    const Eigen::Matrix<T,3,1>& NSd,
    const Eigen::Matrix<T,3,1>& dir,
    const T& Nd,
    const double* px, std::pair<T,T>& k)
{
  Eigen::Matrix<T,3,1> nx;
  nx << px[0],px[1],px[2];
  k.first = nx.dot(NSd)-Nd*px[3];
  k.second = Nd*(dir.dot(nx));
  FASSERT(Nd != 0);
}

template<typename T>
void cal_ppp_point_coord(
    const Eigen::Matrix<T,3,1>& NSd,
    const Eigen::Matrix<T,3,1>& dir,
    const T& Nd, const double* px,
    Eigen::Matrix<T,3,1>& nm, T& dm)
{
  Eigen::Matrix<T,3,1> nx;
  nx << px[0],px[1],px[2];
  T ka = nx.dot(NSd)-Nd*px[3];
  T kb = Nd*(dir.dot(nx));

  nm = -kb*NSd+Nd*ka*dir;
  dm = Nd*kb;
  
  FASSERT(Nd != 0);
}


template<typename T>
void cal_ppt_point(
    const Eigen::Matrix<T,3,1>& NSd,
    const Eigen::Matrix<T,3,1>& dir,
    const T& Nd,
    const double* vi, const double* vj, const double* vk,
    std::pair<T, T>& k)
{
  using namespace fxz::eigen;
  Eigen::Matrix<T,3,1> VI,VJ,VK;
  VI << vi[0],vi[1],vi[2];
  VJ << vj[0],vj[1],vj[2];
  VK << vk[0],vk[1],vk[2];

  auto nt = (VJ-VI).cross(VK-VI);
  T dt = -VI.dot(nt);

  k.first = nt.dot(NSd)-Nd*dt;
  k.second = Nd*dir.dot(nt);
}

template<typename T>
void cal_ppt_point_coord(
    const Eigen::Matrix<T,3,1>& NSd,
    const Eigen::Matrix<T,3,1>& dir,
    const T& Nd,
    const double* vi, const double* vj, const double* vk,
    Eigen::Matrix<T,3,1>& nm, T& dm)
{
  using namespace fxz::eigen;
  Eigen::Matrix<T,3,1> VI,VJ,VK;
  VI << vi[0],vi[1],vi[2];
  VJ << vj[0],vj[1],vj[2];
  VK << vk[0],vk[1],vk[2];

  auto nt = (VJ-VI).cross(VK-VI);
  T dt = -VI.dot(nt);

  T ka = nt.dot(NSd)-Nd*dt;
  T kb = Nd*dir.dot(nt);

  nm = -kb*NSd+Nd*ka*dir;
  dm = Nd*kb;
}


template<typename T>
void cal_ppv_point(
    const Eigen::Matrix<T,3,1>& NSd,
    const Eigen::Matrix<T,3,1>& dir,
    const T& Nd,
    const double* vi, std::pair<T,T>& k)
{
  using namespace fxz::eigen;

  Eigen::Matrix<T,3,1> VI;
  VI << vi[0], vi[1], vi[2];

  k.first = (Nd*VI+NSd).dot(dir);
  k.second = Nd*dir.dot(dir);
}

size_t judge_edge_pass(
    const double* pa, const double* pb,
    const double* va, const double* vb,
    const double* vc, const double* vd);

double cal_vol(const double* na, const double* nb, const double* nc);

int vol_sgn(
    const double* va, const double* vb, const double* vc, const double* vd);

int cal_ppp_point(const double* p0, const double* p1, const double* p2, double* rv);
int cal_ppt_point(
    const double* p0, const double* p1,
    const double* vi, const double* vj, const double* vk,
    double* rv);
int cal_ttp_point(
    const double* va, const double* vb, const double* p, double* rv);

double dsgn(const double x);

uint8_t judge_triangle_pass(
    const double* pa, const double* pb,
    const double* va, const double* vb, const double* vc);

int two_edges_order_on_plane(const double* v, const double* va, const double* vb, const double* p);

bool is_colinear_line(const double* va, const double* vb, const double* pa, const double* pb);

bool is_triangle_and_plane_has_same_norm(
    const double* va, const double* vb, const double* vc, const double* p);


template<typename T>
int cal_ttp_point_coord(
    const double* va, const double* vb, const double* p,
    Eigen::Matrix<T,3,1>& nm, T& dm)
{
  Eigen::Matrix<T,3,1> N, A, B;
  A << va[0], va[1], va[2];
  B << vb[0], vb[1], vb[2];
  N << p[0], p[1], p[2];
  dm = N.dot(B-A);
  //FASSERT(dm != 0);  
  nm = dm*A-(N.dot(A)+p[3])*(B-A);
  return 0;
}

template<typename T>
int cal_ppp_point_coord(
    const double* pa, const double* pb, const double* pc,
    Eigen::Matrix<T,3,1>& nm, T& dm)    
{
  Eigen::Matrix<T,3,1> D;
  Eigen::Matrix<T, 3, 3> N, NS;
  N << pa[0], pa[1], pa[2],
      pb[0], pb[1], pb[2],
      pc[0], pc[1], pc[2];
  D << pa[3], pb[3], pc[3];
  dm = N.determinant();
  N_star(N, NS);
  nm = -NS*D;
  return 1;
}

template<typename T>
int cal_ppt_point_coord(
    const double* pa, const double* pb,
    const double* va, const double* vb, const double* vc,
    Eigen::Matrix<T,3,1>& nm, T& dm)
{
  Eigen::Matrix<T,3,1> A, B, C, Na, Nb, D, Nt;
  Eigen::Matrix<T, 3, 3> N, NS;
  A << va[0], va[1], va[2];
  B << vb[0], vb[1], vb[2];
  C << vc[0], vc[1], vc[2];
  Nt = (B-A).cross(C-A);
  D << pa[3], pb[3], -Nt.dot(A);
  N << pa[0], pa[1], pa[2],
      pb[0], pb[1], pb[2],
      Nt[0], Nt[1], Nt[2];
  N_star(N, NS);
  dm = N.determinant();
  nm = -NS*D;
  return 1;
}

short choose_axis(const Eigen::Matrix<mpf_class,3,1>& D);
short choose_small_axis(const Eigen::Matrix<mpf_class,3,1>& D);
short choose_axis(const Eigen::Matrix<double,3,1>& D);


} // grid
} // fxz


#endif
