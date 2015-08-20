/*
 * Copyright (c) 2015, Georgia Tech Research Corporation
 * All rights reserved.
 *
 * Author(s): Michael X. Grey <mxgrey@gatech.edu>
 *
 * This file is provided under the following "BSD-style" License:
 *   Redistribution and use in source and binary forms, with or
 *   without modification, are permitted provided that the following
 *   conditions are met:
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
 *   CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
 *   INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *   MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *   DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
 *   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 *   USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *   AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *   LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *   ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *   POSSIBILITY OF SUCH DAMAGE.
 */

#include <dart/dart.h>
#include <osgDart/osgDart.h>

using namespace dart::dynamics;

size_t factorial(size_t m)
{
  size_t result = 1u;
  for(size_t i=2; i <= m; ++i)
    result *= i;
  return result;
}

class BezierTau : dart::optimizer::Function
{
  BezierTau(const Eigen::VectorXd& alpha, size_t M = 4)
    : mAlpha(alpha), mM(M) { }

  double computeTau(const Eigen::VectorXd& _x)
  {
    // TODO: implement
  }

  double eval(const Eigen::VectorXd& _x) override final
  {
    double result = 0;
    for(size_t k=0; k < mM; ++k)
    {
      double tau = computeTau(_x);
      result += mAlpha[k] * factorial(mM)/(factorial(k)*factorial(mM-k))
          * pow(tau,k) * pow(1-tau, mM-k);
    }

    return result;
  }


protected:

  Eigen::VectorXd mAlpha;
  size_t mM;
};

class LinearDofCombo : dart::optimizer::Function
{
public:

  struct Term
  {
    Term(size_t _index = INVALID_INDEX, double _coeff = 0.0)
      : index(_index), coeff(_coeff) { }

    size_t index;
    double coeff;
  };

  LinearDofCombo(const std::vector<Term>& terms = std::vector<Term>())
    : mTerms(terms) { }

  double eval(const Eigen::VectorXd& _x) override final
  {
    double result = 0;
    for(const Term& term : mTerms)
      result +=
  }

  void evalGradient(const Eigen::VectorXd& _x,
                    Eigen::Map<Eigen::VectorXd> _grad) override final
  {

  }

protected:

  std::vector<Term> mTerms;

  BezierTau mTau;

};

int main()
{

}
