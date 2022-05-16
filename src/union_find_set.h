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


#ifndef UNION_FIND_SET_H
#define UNION_FIND_SET_H

#include <vector>


namespace fxz {

  class union_find_set
  {
 public:
    union_find_set() { n_=0;}
 union_find_set(size_t n):n_(n) { parent_.resize(n_, n_+5);}
    ~union_find_set() { }

    size_t find(size_t a)
    {
      assert(a<n_);
      return (parent_[a]>n_)?a:(parent_[a]=find(parent_[a]));
    }

    void set_union(size_t a, size_t b)
    {
      assert(a<n_ && b<n_);
      a = find(a); b = find(b);
      if (a != b) parent_[a] = b;
    }

    void set_union_by_order(size_t a, size_t b)
    {
      assert(a<n_ && b<n_);
      a = find(a); b = find(b);
      if (a != b) parent_[a] = b;
    }

    bool is_connected(size_t a, size_t b)
    {
      assert(a<n_ && b<n_);
      return (find(a)==find(b));
    }

    void reset(size_t n)
    {
      n_ = n; parent_.clear();
      parent_.resize(n_, n_+5);
    }

    size_t num() { return n_; }

 private:
    size_t n_;
    std::vector<size_t> parent_;
  };

} // namespace fxz


#endif
