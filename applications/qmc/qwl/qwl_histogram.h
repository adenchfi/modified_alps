/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2004 by Stefan Wessel <wessel@comp-phys.org>
*
* This software is part of the ALPS Applications, published under the ALPS
* Application License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Application License along with
* the ALPS Applications; see the file LICENSE.txt. If not, the license is also
* available from http://alps.comp-phys.org/.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/

#ifndef ALPS_QWL_HISTOGRAM_H
#define ALPS_QWL_HISTOGRAM_H

#include "alps/scheduler/montecarlo.h"
#include <vector>
#include <algorithm>
#include <valarray>

template<typename T=double> 
class histogram {
 public:
  histogram() : data(), left_(0), right_(0) {}
  void resize(unsigned int newsize,unsigned int newleft=0) {
    data.resize(newsize);
    left_=newleft;
    right_=left_+size()-1;
    fill(0);
  }
  void fill(T x) {
    std::fill(data.begin(), data.end(), x);
  }
  void subtract() {
    T x=data[0];
    for (unsigned int i=0;i<data.size();++i)
      data[i]-=x;
  }
  std::valarray<T> getvalarray(unsigned int min,unsigned int max) const {
    std::valarray<T> dval;
    dval.resize(max-min+1);
    int count=0;
    for (int i=min-left();i<=max-left();++i) {
      dval[count]=data[i];
      ++count;
    }
    return dval;
  }
  T& operator[](unsigned int i) {
    return data[i-left()];
  }
  unsigned int left() const {
    return left_;
  }  
  unsigned int right() const {
    return right_;
  }  
  unsigned int size() const {
    return data.size();
  }
  double flatness() const {  
    double av=data[0];
    for (int i=1;i<size();++i)
      av+=data[i];
    av/=size();
    double diff=fabs(data[0]-av);
    for (int i=1;i<size();++i)
      if (diff<fabs(data[i]-av))  
        diff=fabs(data[i]-av);
    return ( av>0 ? diff/av : 0);
  }  
  T min() const {  
   T min=data[0];
    for (int i=1;i<size();++i)
      if (data[i]<min)
         min=data[i];
    return min;
  }  
  void save(alps::ODump& dump) const {
    dump << data << left_;
  }
  void load(alps::IDump& dump) {
    dump >> data >> left_;
    right_=left_+size()-1;
  }
 private:
  std::vector<T> data;
  unsigned int left_;
  unsigned int right_;
};

#endif


