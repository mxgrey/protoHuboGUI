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

#ifndef HUBO_UTILS_HPP
#define HUBO_UTILS_HPP

#include <Eigen/Geometry>

namespace hubo {

const double DefaultIsometryTolerance = 1e-6;
const double DefaultVectorTolerance = 1e-6;

//==============================================================================
static inline bool checkDist(Eigen::Vector3d& p, double a, double b)
{
  double d = p.norm();
  double dmax = a+b;
  double dmin = fabs(a-b);

  if (d > dmax)
  {
    p *= dmax/d;
    return false;
  }
  else if (d < dmin)
  {
    p *= dmin/d;
    return false;
  }
  else
  {
    return true;
  }
}

//==============================================================================
static inline void clamp_sincos(double& sincos, bool& valid)
{
  if (sincos < -1)
  {
    valid = false;
    sincos = -1;
  }
  else if (sincos > 1)
  {
    valid = false;
    sincos = 1;
  }
}

//==============================================================================
static inline Eigen::Vector3d flipEuler3Axis(const Eigen::Vector3d& u)
{
  Eigen::Vector3d v;
  v[0] = u[0] - M_PI;
  v[1] = M_PI - u[1];
  v[2] = u[2] - M_PI;
  return v;
}

//==============================================================================
inline double diff(const Eigen::Isometry3d& tf0,
                   const Eigen::Isometry3d& tf1)
{
  double cost = (tf0.translation() - tf1.translation()).norm();
  Eigen::AngleAxisd aa(tf0.linear().transpose()*tf1.linear());
  cost += aa.angle();

  return cost;
}

//==============================================================================
inline bool equal(const Eigen::Isometry3d& tf0,
                  const Eigen::Isometry3d& tf1,
                  const double tolerance = DefaultIsometryTolerance)
{
  return (diff(tf0, tf1) < tolerance);
}

//==============================================================================
inline bool equal(const Eigen::VectorXd& v0,
                  const Eigen::VectorXd& v1,
                  const double tolerance = DefaultVectorTolerance)
{
  return ((v1-v0).norm() < tolerance);
}

//==============================================================================
inline Eigen::AngleAxisd findRotation(
    const Eigen::Vector3d& x1, const Eigen::Vector3d& x0)
{
  Eigen::Vector3d N = x1.cross(x0);
  double s = N.norm();
  double theta = asin(s);
  // Deal with numerical imprecision issues
  if(s < -1.1 || 1.1 < s)
  {
    std::cout << "[findRotation] Strange value for s: " << s
              << " | Results in theta: " << theta << std::endl;
    return Eigen::AngleAxisd(0, Eigen::Vector3d::UnitZ());
  }
  else if(1.0 <= s)
  {
    theta = M_PI/2.0;
  }
  else if(s <= -1.0)
  {
    theta = -M_PI/2.0;
  }

  N.normalize();

  const double wrap = x1.dot(x0);
  if(wrap < 0)
    theta = M_PI - theta;

  if(std::abs(theta-M_PI) < 1e-4
     || std::abs(theta+M_PI) < 1e-4)
    return Eigen::AngleAxisd(M_PI, Eigen::Vector3d::UnitZ());

  return Eigen::AngleAxisd(theta, N);
}



} // namespace hubo

#endif // HUBO_UTILS_HPP
