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
#include <stdio.h>
#include <gmp.h>
#include <gmpxx.h>

#include "common.h"
#include "my_gmp.h"

namespace fxz {
namespace grid {

typedef Eigen::Matrix<mpf_class,3,1> MPF_VEC;

// np^T v+dp
short cal_pv(const double* p, const double* v)
{
  mpf_class n1(p[0]),n2(p[1]),n3(p[2]),v1(v[0]),v2(v[1]),v3(v[2]),d(p[3]);
  mpf_class val = n1*v1+n2*v2+n3*v3+d;
  return sgn(val);
}


double cal_vol(const double* na, const double* nb, const double* nc)
{
  using namespace fxz::eigen;
  cmap_evec3d N1(na, 3, 1);
  cmap_evec3d N2(nb, 3, 1);
  cmap_evec3d N3(nc, 3, 1);
  return N1.cross(N2).dot(N3);
}

int vol_sign_with_norms(
    const double* na, const double* nb, const double* nc)
{
  MPF_VEC NA, NB, NC;
  NA << na[0], na[1], na[2];
  NB << nb[0], nb[1], nb[2];
  NC << nc[0], nc[1], nc[2];
  
  return sgn(NA.cross(NB).dot(NC));
}

uint8_t cal_pp(const double* pa, const double* pb)
{
  //typedef Eigen::Matrix<mpf_class,3,1> MPF_VEC;
  MPF_VEC NA, NB;
  NA << pa[0], pa[1], pa[2];
  NB << pb[0], pb[1], pb[2];
  MPF_VEC N = NA.cross(NB);
  int nij = abs(sgn(N[0]))+abs(sgn(N[1]))+abs(sgn(N[2]));
  int d = sgn((pb[3]*NA-pa[3]*NB).dot(NA));
  if (nij==0 && d!=0) return 0;
  else if (nij!=0) return 2;
  else if (nij==0 && d==0) return 3;
  std::cerr << "# [ ERROR ] invalid relation between two planes." << std::endl;
  return 10; // invalid
}

// check whether the intersection is only one point.
bool cal_ppt_intersection(
    const double* pa, const double* pb,
    const double* va, const double* vb, const double* vc)
{
  //typedef Eigen::Matrix<mpf_class,3,1> MPF_VEC;
  MPF_VEC NA, NB, VA, VB, VC;
  NA << pa[0], pa[1], pa[2];
  NB << pb[0], pb[1], pb[2];
  VA << va[0], va[1], va[2];
  VB << vb[0], vb[1], vb[2];
  VC << vc[0], vc[1], vc[2];

  mpf_class b0 = -(pa[3]+NA.dot(VA));
  mpf_class b1 = -(pb[3]+NB.dot(VA));

  mpf_class M00 = NA.dot(VB-VA);
  mpf_class M01 = NA.dot(VC-VA);
  mpf_class M10 = NB.dot(VB-VA);
  mpf_class M11 = NB.dot(VC-VA);

  mpf_class det = M00*M11-M01*M10;

  mpf_class c0 = M11*b0-M01*b1;
  mpf_class c1 = -M10*b0+M00*b1;

  int ds = sgn(det);
  int s0 = sgn(c0);
  int s1 = sgn(c1);
  int s01 = sgn(det-c0-c1);
  
  if (ds == 0) {return false;}
  
  return (ds*s0>=0 && ds*s1>=0 && ds*s01>=0);
}


void cal_ttp_point(
    const double* va, const double* vb, const double* p,
    std::pair<double, double>& k)
{
  using namespace fxz::eigen;
  cmap_evecd VA(va, 3, 1);
  cmap_evecd VB(vb, 3, 1);
  cmap_evecd N(p, 3, 1);
  k.first = -(N.dot(VA)+p[3]);
  k.second = N.dot(VB-VA);
}

int cmp_tt_point(
    const double* va, const double* vb,
    const double* p1, const double* p2)
{
  Eigen::Matrix<mpf_class,3,1> VA,VB,N1,N2;
  VA << va[0], va[1], va[2];
  VB << vb[0], vb[1], vb[2];
  N1 << p1[0], p1[1], p1[2];
  N2 << p2[0], p2[1], p2[2];

  mpf_class ka1 = -(N1.dot(VA)+p1[3]);
  mpf_class ka2 = N1.dot(VB-VA);
  mpf_class kb1 = -(N2.dot(VA)+p2[3]);
  mpf_class kb2 = N2.dot(VB-VA);
  // ka1/ka2 : kb1/kb2
  return sgn(ka1*kb2-ka2*kb1)*sgn(ka2*kb2);
}

int cmp_tt_point_to_const(
    const double* va, const double* vb, const double* p, const double kc)
{
  Eigen::Matrix<mpf_class,3,1> VA, VB, N;
  VA << va[0], va[1], va[2];
  VB << vb[0], vb[1], vb[2];
  N  << p[0], p[1], p[2];
  mpf_class k1 = -(N.dot(VA)+p[3]);
  mpf_class k2 = N.dot(VB-VA);

  return sgn(k1-k2*kc)*sgn(k2);
}

template<typename T,typename R>
inline R pt_point_first(const T& VI,const T& VJ, const T& VK, const T& NP, const T& NX, const double d, const double dx)
{ return -(VK-VI).dot(NP)*(VJ-VI).dot((dx+VI.dot(NX))*NP-(d+VI.dot(NP))*NX); }

template<typename T,typename R>
inline R pt_point_second(const T& VI,const T& VJ, const T& VK, const T& NP, const T& NX, const double d)
{ return (d+VI.dot(NP))*(NP.dot(VJ-VI)*(VI-VK).dot(NX)+NP.dot(VK-VI)*(VJ-VI).dot(NX)); }


void cal_pt_point(
    const double* va, const double* vb, const double* vc,
    const double* p, const double* px, std::pair<double,double>& k)
{
  using namespace fxz::eigen;
  cmap_evecd VI(va, 3, 1);
  cmap_evecd VJ(vb, 3, 1);
  cmap_evecd VK(vc, 3, 1);
  cmap_evecd NP(p,  3, 1);
  cmap_evecd NX(px, 3, 1);

  k.first = pt_point_first<cmap_evecd,double>(VI,VJ,VK,NP,NX,p[3],px[3]);
  k.second = pt_point_second<cmap_evecd,double>(VI,VJ,VK,NP,NX,p[3]);
}

int cmp_pt_point(
    const double* va, const double* vb, const double* vc,
    const double* p, const double* pa, const double* pb)
{
  MPF_VEC VI,VJ,VK,N,NA,NB;
  VI << va[0], va[1], va[2];
  VJ << vb[0], vb[1], vb[2];
  VK << vc[0], vc[1], vc[2];
  N  << p[0], p[1], p[2];
  NA << pa[0], pa[1], pa[2];
  NB << pb[0], pb[1], pb[2];

  mpf_class ka1 = pt_point_first<MPF_VEC,mpf_class>(VI,VJ,VK,N,NA,p[3],pa[3]);
  mpf_class ka2 = pt_point_second<MPF_VEC,mpf_class>(VI,VJ,VK,N,NA,p[3]);
  mpf_class kb1 = pt_point_first<MPF_VEC,mpf_class>(VI,VJ,VK,N,NB,p[3],pb[3]);
  mpf_class kb2 = pt_point_second<MPF_VEC,mpf_class>(VI,VJ,VK,N,NB,p[3]);
  return sgn(ka1*kb2-ka2*kb1)*sgn(ka2*kb2);
}

int cmp_pt_point_to_const(
    const double* va, const double* vb, const double* vc,
    const double* p, const double* px, const double kc)
{
  MPF_VEC VI,VJ,VK,NP,NX;
  VI << va[0], va[1], va[2];
  VJ << vb[0], vb[1], vb[2];
  VK << vc[0], vc[1], vc[2];
  NP  << p[0], p[1], p[2];
  NX << px[0], px[1], px[2];

  mpf_class k1 = pt_point_first<MPF_VEC,mpf_class>(VI,VJ,VK,NP,NX,p[3],px[3]);
  mpf_class k2 = pt_point_second<MPF_VEC,mpf_class>(VI,VJ,VK,NP,NX,p[3]);
  return sgn(k1-k2*kc)*sgn(k2);
}

size_t judge_edge_pass(
    const double* pa, const double* pb,
    const double* va, const double* vb,
    const double* vc, const double* vd)
{
  MPF_VEC NA, NB, VA, VB, VC, VD;
  NA << pa[0], pa[1], pa[2];
  NB << pb[0], pb[1], pb[2];
  VA << va[0], va[1], va[2];
  VB << vb[0], vb[1], vb[2];
  VC << vc[0], vc[1], vc[2];
  VD << vd[0], vd[1], vd[2];

  MPF_VEC ni = (VB-VA).cross(VC-VA);
  MPF_VEC nj = (VD-VA).cross(VB-VA);
  MPF_VEC L = NA.cross(NB);

  int di = sgn(ni.dot(L));
  int dj = sgn(nj.dot(L));
  int r  = sgn((VB-VA).dot(ni.cross(nj)));

  if (di*dj>0 || (di*dj==0 && r<0)) {return 1;}
  
  return 0;
}

int vol_sgn(
    const double* a, const double* b, const double* c, const double* d)
{
  MPF_VEC A, B, C, D;
  A << a[0], a[1], a[2];
  B << b[0], b[1], b[2];
  C << c[0], c[1], c[2];
  D << d[0], d[1], d[2];

  return sgn((B-A).cross(C-A).dot(D-A));
}

int cal_ppp_point(
    const double* p0, const double* p1, const double* p2,
    double* rv)
{
  Eigen::Matrix<double,3,1> NSd, dir;
  double Nd;
  pp_precompute(p0, p1, p2, NSd, dir, Nd);
  FASSERT2(fabs(Nd)>0, "ppp");

  eigen::map_evec3d R(rv, 3, 1);
  R = -NSd/Nd;

  return 0;
}

int cal_ppt_point(
    const double* p0, const double* p1,
    const double* vi, const double* vj, const double* vk,
    double* rv)
{
  using namespace fxz::eigen;
  cmap_evec3d VI(vi, 3, 1);
  cmap_evec3d VJ(vj, 3, 1);
  cmap_evec3d VK(vk, 3, 1);

  evec3d Nt = (VJ-VI).cross(VK-VI);
  double pt[4];
  pt[0] = Nt[0];
  pt[1] = Nt[1];
  pt[2] = Nt[2];
  pt[3] = -VI.dot(Nt);

  FASSERT(pt == &pt[0]);

  Eigen::Matrix<double,3,1> NSd, dir;
  double Nd(0);
  pp_precompute(p0, p1, pt, NSd, dir, Nd);
  
  FASSERT2(fabs(Nd)>0, "ppt");

  map_evec3d R(rv, 3, 1);
  R = -NSd/Nd;

  return 0;
}

int cal_ttp_point(
    const double* va, const double* vb, const double* p,
    double* rv)
{
  std::pair<double, double> k;
  cal_ttp_point(va, vb, p, k);

  FASSERT(fabs(k.second)>0);

  using namespace fxz::eigen;
  cmap_evec3d VA(va, 3, 1);
  cmap_evec3d VB(vb, 3, 1);
  map_evec3d R(rv, 3, 1);
  R = VA+(VB-VA)*k.first/k.second;

  return 0;
}

double dsgn(const double x)
{
  if (x > 0) return 1.0;
  else if (x < 0) return -1.0;
  return 0.0;
}

uint8_t judge_triangle_pass(
    const double* pa, const double* pb,
    const double* va, const double* vb, const double* vc)
{
  //typedef Eigen::Matrix<mpf_class,3,1> MPF_VEC;
  MPF_VEC NA, NB, A, B, C;
  NA << pa[0], pa[1], pa[2];
  NB << pb[0], pb[1], pb[2];
  A << va[0], va[1], va[2];
  B << vb[0], vb[1], vb[2];
  C << vc[0], vc[1], vc[2];

  int r = sgn((B-A).cross(C-A).dot(NA.cross(NB)));
  
  if (r != 0) return 1;
  return 0;
}


int two_edges_order_on_plane(
    const double* v, const double* va, const double* vb, const double* p)
{
  //typedef Eigen::Matrix<mpf_class,3,1> MPF_VEC;
  MPF_VEC N, V, A, B;
  N << p[0] , p[1], p[2];
  V << v[0], v[1], v[2];
  A << va[0], va[1], va[2];
  B << vb[0], vb[1], vb[2];

  return sgn((A-V).cross(B-V).dot(N));
}

// assume the line is passing one point of the edge.
bool is_colinear_line(const double* va, const double* vb, const double* pa, const double* pb)
{
  //typedef Eigen::Matrix<mpf_class,3,1> MPF_VEC;
  MPF_VEC NA, NB, A, B;
  A << va[0], va[1], va[2];
  B << vb[0], vb[1], vb[2];
  NA << pa[0], pa[1], pa[2];
  NB << pb[0], pb[1], pb[2];


  MPF_VEC c = (A-B).cross(NA.cross(NB));

  if (sgn(c[0])==0 && sgn(c[1])==0 && sgn(c[2])==0) {
    return true;
  }
  return false;
}

bool is_triangle_and_plane_has_same_norm(
    const double* va, const double* vb, const double* vc, const double* p)
{

  MPF_VEC A, B, C, N;
  A << va[0], va[1], va[2];
  B << vb[0], vb[1], vb[2];
  C << vc[0], vc[1], vc[2];
  N << p[0], p[1], p[2];

  int s = sgn((B-A).cross(C-A).dot(N));
  
  return (s>0);
}

// |x|<=eps, 0; x > eps, 1; x < -eps, -1;
int my_sgn(const mpf_class& x, const double eps)
{
  // cmp: a>b, 1; a<b, -1; a=b, 0;
  int r = cmp(abs(x), eps);
  if (r <= 0) return 0;
  return sgn(x);
}


short choose_axis(const Eigen::Matrix<mpf_class,3,1>& D)
{
  auto s01 = cmp(abs(D[0]),abs(D[1]));
  auto s02 = cmp(abs(D[0]),abs(D[2]));
  auto s12 = cmp(abs(D[1]),abs(D[2]));

  if (s01>=0 && s02>=0) { return sgn(D[0]); }
  if (s01<=0 && s12>=0) { return sgn(D[1])*2; }

  return sgn(D[2])*3;
}

short choose_small_axis(const Eigen::Matrix<mpf_class,3,1>& D)
{
  auto s01 = cmp(abs(D[0]),abs(D[1]));
  auto s02 = cmp(abs(D[0]),abs(D[2]));
  auto s12 = cmp(abs(D[1]),abs(D[2]));

  if (s01<=0 && s02<=0) { return sgn(D[0]); }
  if ((s01>=0 && s12<=0)) { return sgn(D[1])*2; }

  return sgn(D[2])*3;
}

short choose_axis(const Eigen::Matrix<double,3,1>& D)
{
  if (fabs(D[0])>=fabs(D[1]) && fabs(D[0])>=fabs(D[2])) { return dsgn(D[0]); }
  if (fabs(D[0])<=fabs(D[1]) && fabs(D[1])>=fabs(D[2])) { return dsgn(D[1])*2; }
  return dsgn(D[2])*3;
}


} //grid
} //fxz
