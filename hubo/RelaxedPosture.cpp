/*
 * Copyright (c) 2016, Georgia Tech Research Corporation
 * All rights reserved.
 *
 * Author(s): Michael X. Grey <mxgrey@gatech.edu>
 *
 * Georgia Tech Graphics Lab and AMBER Lab
 *
 * Directed by Prof. C. Karen Liu and Prof. Aaron Ames
 * <karenliu@cc.gatech.edu> <ames@gatech.edu>
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

#include <dart/common/Console.h>

#include "hubo/RelaxedPosture.hpp"

namespace hubo {

//==============================================================================
RelaxedPosture::RelaxedPosture(
    const Eigen::VectorXd& idealPosture,
    const Eigen::VectorXd& lowerIdeal, const Eigen::VectorXd& upperIdeal,
    const Eigen::VectorXd& lowerOkay, const Eigen::VectorXd& upperOkay,
    const Eigen::VectorXd& weights, bool enforceIdeal)
  : enforceIdealPosture(enforceIdeal),
    mIdeal(idealPosture),
    mIdealLower(lowerIdeal),
    mIdealUpper(upperIdeal),
    mLower(lowerOkay),
    mUpper(upperOkay),
    mWeights(weights)
{
  int dofs = mIdeal.size();
  if(mLower.size() != dofs || mWeights.size() != dofs || mUpper.size() != dofs)
  {
    dterr << "[RelaxedPose::RelaxedPose] Dimension mismatch:\n"
          << "  ideal:   " << mIdeal.size()   << "\n"
          << "  lower:   " << mLower.size()   << "\n"
          << "  upper:   " << mUpper.size()   << "\n"
          << "  weights: " << mWeights.size() << "\n";
  }
  mResultVector.setZero(dofs);
}

//==============================================================================
double RelaxedPosture::eval(const Eigen::VectorXd& _x)
{
  computeResultVector(_x);
  return 0.5 * mResultVector.dot(mResultVector);
}

//==============================================================================
void RelaxedPosture::evalGradient(
    const Eigen::VectorXd& _x,
    Eigen::Map<Eigen::VectorXd> _grad)
{
  computeResultVector(_x);

  _grad.setZero();
  int smaller = std::min(mResultVector.size(), _grad.size());
  for(int i=0; i < smaller; ++i)
    _grad[i] = mResultVector[i];
}

//==============================================================================
void RelaxedPosture::computeResultVector(const Eigen::VectorXd& _x)
{
  mResultVector.setZero();

  if(enforceIdealPosture)
  {
    for(int i=0; i < _x.size(); ++i)
    {
      if(mIdeal.size() <= i)
        break;

      if(mIdealLower[i] < _x[i] && _x[i] < mIdealUpper[i])
        continue;

      mResultVector[i] = mWeights[i]*(_x[i] - mIdeal[i]);
    }
  }
  else
  {
    for(int i=0; i < _x.size(); ++i)
    {
      if(mIdeal.size() <= i)
        break;

      if(_x[i] < mLower[i])
        mResultVector[i] = mWeights[i]*(_x[i] - mLower[i]);
      else if(mUpper[i] < _x[i])
        mResultVector[i] = mWeights[i]*(_x[i] - mUpper[i]);
    }
  }
}


} // namespace hubo
