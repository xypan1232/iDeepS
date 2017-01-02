// -*- C++ -*-
// Little library of vectors and sparse vectors
// Copyright (C) 2007- Leon Bottou

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111, USA



// $Id: vectors.h,v 1.13 2007/10/02 20:40:06 cvs Exp $

#ifndef VECTORS_H
#define VECTORS_H 1

#include <cstring>
#include <iostream>

#include <sstream>
#include <fstream>

#include <cassert>

#include <vector>
#include <cmath>

#include "wrapper.h"

class FVector;
class SVector;

typedef float VFloat;

class FVector
{
private:
  struct Rep
  {
    int refcount;
    int size;
    VFloat *data;
    Rep() : size(0), data(0) {}
    ~Rep() { delete [] data; }
    void resize(int n);
    Rep *copy();
  };

  Wrapper<Rep> w;
  Rep *rep() { return w.rep(); }
  const Rep *rep() const { return w.rep(); }
  void qset(int i, double v);

public:
  FVector();
  FVector(int n);
  FVector(const SVector &v);
  int size() const { return rep()->size; }

  // these accessors are range-checked.
  // get() returns 0 when i is out-of-range.
  // set() expands the vector.
  double get(int i) const;
  double set(int i, double v);

  // warning: you can write vector[i] but
  // the subscripts are not range-checked!
  // on the other hand, that's fast.
  operator const VFloat* () const { return rep()->data; }
  operator VFloat* () { w.detach(); return rep()->data; }

  void clear();
  void resize(int n);
  void touch(int i);
  FVector slice(int fi, int ti) const;

  void add(double c1);
  void add(const FVector &v2);
  void add(const SVector &v2);
  void add(const FVector &v2, double c2);
  void add(const SVector &v2, double c2);
  void add(const FVector &v2, double c2, const FVector &q2);
  void add(const SVector &v2, double c2, const FVector &q2);
  void scale(double c1);
  void combine(double c1, const FVector &v2, double c2);
  void combine(double c1, const SVector &v2, double c2);


  friend std::ostream& operator<<(std::ostream &f, const FVector &v);
  friend std::istream& operator>>(std::istream &f, FVector &v);
  bool save(std::ostream &f) const;
  bool load(std::istream &f);
};



class SVector
{
public:
  struct Pair
  {
    int i;
    VFloat v;
  };
private:
  struct Rep
  {
    int refcount;
    int npairs;
    int mpairs;
    int size;
    struct Pair *pairs;

    Rep() : npairs(0), mpairs(-1), size(0), pairs(0) {}
    ~Rep() { delete [] pairs; }
    void resize(int n);
    double qset(int i, double v);
    Rep *copy();
  };

  Wrapper<Rep> w;
  Rep *rep() { return w.rep(); }
  const Rep *rep() const { return w.rep(); }

public:
  SVector();
  SVector(const FVector &v);
  int size() const { return rep()->size; }

  // these accessors are range-checked.
  // get() returns 0 when i is out-of-range.
  // set() expands the vector.
  double get(int i) const;
  double set(int i, double v);

  // to quickly iterate over the non-zero coefficients,
  // do for(SVector::Pair *p = x; p->i>=0; p++) { ... }
  int npairs() const { return rep()->npairs; }
  operator const Pair* () const { return rep()->pairs; }

  void clear();
  void trim();
  SVector slice(int fi, int ti) const;

  std::pair<unsigned, double> extract_component(unsigned aFeatureOrderIndex);
  int sparse_size();

  void add(const SVector &v2);
  void add(const SVector &v2, double c2);
  void scale(double c1);
  void combine(double c1, const SVector &v2, double c2);
  void normalize();
  std::vector< std::pair <int,double> >  unpack();
  void binarize();


  friend std::ostream& operator<<(std::ostream &f, const SVector &v);
  friend std::istream& operator>>(std::istream &f, SVector &v);
  bool save(std::ostream &f) const;
  bool load(std::istream &f);

  friend SVector combine(const SVector &v1, double a1,
                         const SVector &v2, double a2);
  friend SVector elementwise_product(const SVector &v1, const SVector &v2, const double hardlimit);

};

double dot(const FVector &v1, const FVector &v2);
double dot(const FVector &v1, const SVector &v2);
double dot(const SVector &v1, const FVector &v2);
double dot(const SVector &v1, const SVector &v2);

SVector combine(const SVector &v1, double a1, const SVector &v2, double a2);
FVector combine(const FVector &v1, double a1, const SVector &v2, double a2);
FVector combine(const SVector &v1, double a1, const FVector &v2, double a2);
FVector combine(const FVector &v1, double a1, const FVector &v2, double a2);

SVector elementwise_product(const SVector &v1, const SVector &v2, const double hardlimit=1e9);

#endif

/* -------------------------------------------------------------
   Local Variables:
   c++-font-lock-extra-types: ("\\sw+_t" "[A-Z]\\sw*[a-z]\\sw*" "std::\\sw+")
   End:
   ------------------------------------------------------------- */
