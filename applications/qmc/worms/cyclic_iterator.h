/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2001-2004 by Matthias Troyer <troyer@comp-phys.org>,
*                            Simon Trebst <trebst@comp-phys.org>
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

/* $Id: cyclic_iterator.h 1345 2004-10-07 04:06:39Z wistaria $ */

#ifndef CYCLIC_ITERATOR_H___
#define CYCLIC_ITERATOR_H___

#include <boost/throw_exception.hpp>
#include <stdexcept>
#include <cassert>
#include <iostream>

#ifdef ICC_CYCLIC_ITERATOR

template <class C>
class cyclic_iterator
{
public:
  typedef C container_type;
  typedef typename container_type::value_type base_value_type;
  typedef typename container_type::iterator   base_iterator;
  
  cyclic_iterator() : c_(0) {}
  cyclic_iterator(container_type& c) : c_iterator(c.end()), c_(&c) {}
  cyclic_iterator(container_type& c, base_iterator it) : c_iterator(it), c_(&c) {}

  C& container() { return *c_;}
  const C& container() const { return *c_;}
  
  const cyclic_iterator& operator=(const base_iterator& it)
  {
    if(c_==0)
      boost::throw_exception(std::logic_error("(1) container not set in cyclic_iterator"));
    c_iterator = it;
    return *this;
  }

  // added function
  const cyclic_iterator& operator=(const cyclic_iterator& cit)
  {
    c_ = cit.c_;
    if(c_==0)
      boost::throw_exception(std::logic_error("(2) container not set in cyclic_iterator"));
    c_iterator = cit.c_iterator;
    return *this;
  }

  void check_container(const container_type& c)
  {
    if (&c!=c_) {
      std::cerr << "Incorrect container: " << &c << " " << c_ << "\n";
      boost::throw_exception(std::runtime_error("Incorrect container"));
    }
  }
  
  bool valid () const    { return c_iterator != c_->end(); }
  operator bool () const { return valid(); }

  void invalidate () { c_iterator = c_->end(); }

  // --- new dereference operators ---
  const base_value_type& operator*() const {
    return *c_iterator;
  }

  const base_value_type* operator->() const {
    return &(*c_iterator);
  }

  operator base_iterator() const {
    return c_iterator;
  } 

  bool operator==(const  base_iterator& it) const {
    return c_iterator==it;
  }

  bool operator!=(const  base_iterator& it) const {
    return c_iterator!=it;
  }

  cyclic_iterator& operator++() {
    c_iterator++; 
    if( c_iterator==c_->end() ) {
      c_iterator = c_->begin();
    }
    return *this;
  }

  cyclic_iterator operator++(int) { 
    cyclic_iterator tmp(*this);
    operator++();
    return tmp;
  }

  cyclic_iterator& operator--() { 
    if( c_iterator==c_->begin() ) {
      c_iterator=c_->end();
    }
    c_iterator--;
    return *this;
  }

  cyclic_iterator operator--(int) { 
    cyclic_iterator old(*this);
    operator--();
    return old;
  }
  
  cyclic_iterator operator+(int n) {
    assert (n==1);
    cyclic_iterator result(*this);
    return ++result;
  }

  cyclic_iterator operator-(int n) {
    assert (n==1);
    cyclic_iterator result(*this);
    return --result;
  }

 public:
  base_iterator c_iterator;

 private:
  container_type* c_;
};

#else

template <class C>
class cyclic_iterator : public C::iterator
{
public:
  typedef C container_type;
  typedef typename container_type::iterator base_iterator;
  cyclic_iterator() : c_(0) {}
  cyclic_iterator(container_type& c) : base_iterator(c.end()), c_(&c) {}
  cyclic_iterator(container_type& c, base_iterator it) : base_iterator(it), c_(&c) {}

  C& container() { return *c_;}
  const C& container() const { return *c_;}
  
  const cyclic_iterator& operator=(const base_iterator& it)
  {
    if (c_==0)
      boost::throw_exception(std::logic_error("container not set in cyclic_iterator"));
    base_iterator::operator=(it);
    return *this;
  }

  void check_container(const container_type& c)
  {
    if (&c!=c_) {
      std::cerr << "Incorrect container: " << &c << " " << c_ << "\n";
      boost::throw_exception(std::runtime_error("Incorrect container"));
      }
  }
  bool valid () const {return *this != c_->end();}
  operator bool () const {return valid();}

  void invalidate () { (*this) = c_->end(); }

  cyclic_iterator& operator++() { 
    base_iterator::operator++();
    if ((*this)==c_->end()) {
      (*this) = c_->begin();
    }
    return *this;
  }

  cyclic_iterator operator++(int) { 
    cyclic_iterator tmp(*this);
    operator++();
    return tmp;
  }

  cyclic_iterator& operator--() { 
    if((*this)==c_->begin()) {
      (*this)=c_->end();
    }
    base_iterator::operator--();
    return *this;
  }

  cyclic_iterator operator--(int) { 
    cyclic_iterator old(*this);
    operator--();
    return old;
  }
  
  cyclic_iterator operator+(int n) {
    assert (n==1);
    cyclic_iterator result(*this);
    return ++result;
  }

  cyclic_iterator operator-(int n) {
    assert (n==1);
    cyclic_iterator result(*this);
    return --result;
  }

private:
  container_type* c_;
};

#endif


#endif

