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

#ifndef HUBO_RELAXEDPOSTURE_HPP
#define HUBO_RELAXEDPOSTURE_HPP

#include <dart/optimizer/Function.h>

namespace hubo {

class RelaxedPosture : public dart::optimizer::Function
{
public:

  RelaxedPosture(const Eigen::VectorXd& idealPosture,
                 const Eigen::VectorXd& lowerIdeal, const Eigen::VectorXd& upperIdeal,
                 const Eigen::VectorXd& lowerOkay, const Eigen::VectorXd& upperOkay,
                 const Eigen::VectorXd& weights, bool enforceIdeal = false);

  double eval(const Eigen::VectorXd& _x) override;

  void evalGradient(const Eigen::VectorXd& _x,
                    Eigen::Map<Eigen::VectorXd> _grad) override;

  void computeResultVector(const Eigen::VectorXd& _x);

  bool enforceIdealPosture;

  Eigen::VectorXd mResultVector;
  Eigen::VectorXd mIdeal;
  Eigen::VectorXd mIdealLower;
  Eigen::VectorXd mIdealUpper;
  Eigen::VectorXd mLower;
  Eigen::VectorXd mUpper;
  Eigen::VectorXd mWeights;

};


} // namespace hubo

#endif // HUBO_RELAXEDPOSTURE_HPP
