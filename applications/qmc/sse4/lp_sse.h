/*****************************************************************************
 *
 * ALPS Project Applications
 *
 * Copyright (C) 2010 by Lode Pollet <pollet@itp.phys.ethz.ch>
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

#ifndef ALPS_APPLICATIONS_LP_SSE
#define ALPS_APPLICATIONS_LP_SSE

#include <iostream>
#include <vector>

using namespace std;

class LP_Solve_SSE {
public:
  void init(const vector<double>& );
  void calc_solution(vector<double>& );
  double calc_loc_opt_bounce(); // NOT NEEDED FOR Linear Programming !!!!
private:
  void do_pivot(const int, const int);
  bool check_optimal();
  void print_tableau() const;
  void copy_solution(vector<double>& );
  int find_pivot_column() const;
  int find_pivot_row(const int) const;
  int Nvar;
  int Nslack;
  vector<int> mVar;
  vector<double> mWeight;
  vector<double> C; // optimization expression
  vector<vector<double> > A; // Matrix of variables
  vector<double> B; // right hand side of inequalies;
  vector<double> mSolution;  
  double mBounce;
};

#endif
