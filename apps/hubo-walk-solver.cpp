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
using namespace dart::optimizer;

size_t factorial(size_t m)
{
  size_t result = 1u;
  for(size_t i=2; i <= m; ++i)
    result *= i;
  return result;
}

double computeBezier(double tau, const Eigen::VectorXd& alpha, double M = 4)
{
  double result = 0.0;
  for(size_t k = 0; k < M; ++k)
  {
    result += alpha[k] * factorial(M) / (factorial(k) * factorial(M-k))
        * pow(tau, k) * pow(1-tau, M-k);
  }

  return result;
}

class BezierDependentFunc
{
public:

  void setTau(double tau)
  {
    mTau = tau;
  }

protected:

  double mTau;
};

typedef std::shared_ptr<BezierDependentFunc> BezierFuncPtr;
typedef std::vector<BezierFuncPtr> BezierFuncArray;

class LinearCombo : public Function,
                       public BezierDependentFunc
{
public:

  struct Term
  {
    Term(size_t _index = INVALID_INDEX, double _coeff = 0.0)
      : index(_index), coeff(_coeff) { }

    size_t index;
    double coeff;
  };

  LinearCombo(const std::vector<Term>& terms,
              const Eigen::VectorXd& alphas)
    : mTerms(terms), mAlphas(alphas) { }

  double computeCost(const Eigen::VectorXd& _x)
  {
    double cost = 0;
    for(const Term& term : mTerms)
      cost += term.coeff * _x[term.index];

    cost -= computeBezier(mTau, mAlphas);
    return cost;
  }

  double eval(const Eigen::VectorXd& _x) override final
  {
    const double cost = computeCost(_x);
    return cost*cost;
  }

  void evalGradient(const Eigen::VectorXd& _x,
                    Eigen::Map<Eigen::VectorXd> _grad) override final
  {
    _grad.setZero();
    const double cost = computeCost(_x);
    for(const Term& term : mTerms)
      _grad[term.index] = 2 * term.coeff * cost;
  }

protected:

  std::vector<Term> mTerms;

  Eigen::VectorXd mAlphas;
};

class EndEffectorConstraint : public Function,
                              public BezierDependentFunc
{
public:

  /// alphas should have six entries, one corresponding to each component of the
  /// EndEffector pose. Unused alphas can be set to all zeros.
  EndEffectorConstraint(JacobianNode* node,
                        const std::vector<Eigen::VectorXd>& alphas)
    : mAlphaArray(alphas)
  {
    mIK = InverseKinematics::create(node);

    const SkeletonPtr& skel = node->getSkeleton();

    std::vector<size_t> indices;
    indices.reserve(skel->getNumDofs());
    for(size_t i=0; i < skel->getNumDofs(); ++i)
      indices.push_back(i);
  }

  double eval(const Eigen::VectorXd& _x) override final
  {
    Eigen::Vector6d screw;
    for(size_t i=0; i < 6; ++i)
      screw[i] = computeBezier(mTau, mAlphaArray[i]);

    Eigen::Isometry3d tf(Eigen::Isometry3d::Identity());
    // TODO: Find out how to convert the screw into an Isometry3d

    mIK->getTarget()->setTransform(tf);

    const Eigen::Vector6d& error = mIK->getErrorMethod().evalError(_x);
    return error.norm();
  }

  void evalGradient(const Eigen::VectorXd& _x,
                    Eigen::Map<Eigen::VectorXd> _grad) override final
  {
    mIK->getGradientMethod().evalGradient(_x, _grad);
  }

  InverseKinematicsPtr mIK;

protected:

  std::vector<Eigen::VectorXd> mAlphaArray;

};

void build_alphas(std::vector<double>& /*alphas*/)
{
  // Terminate template recursion
}

template <typename ... Args>
void build_alphas(std::vector<double>& alphas, double nextAlpha, Args... args)
{
  alphas.push_back(nextAlpha);
  build_alphas(alphas, args...);
}

template <typename ... Args>
Eigen::VectorXd make_alphas(Args... args)
{
  std::vector<double> alphas;
  build_alphas(alphas, args...);
  Eigen::VectorXd vec(alphas.size());
  for(size_t i=0; i < alphas.size(); ++i)
    vec[i] = alphas[i];

  return vec;
}

void build_terms(std::vector<LinearCombo::Term>& /*terms*/)
{
  // Terminate template recursion
}

template <typename ... Args>
void build_terms(std::vector<LinearCombo::Term>& terms,
                 size_t index, double coeff, Args... args)
{
  terms.push_back(LinearCombo::Term(index, coeff));
  build_terms(terms, args...);
}

template <typename ... Args>
std::vector<LinearCombo::Term> make_terms(Args... args)
{
  std::vector<LinearCombo::Term> terms;
  build_terms(terms, args...);
  return terms;
}

Eigen::VectorXd solve(const std::shared_ptr<Solver>& solver,
                      const BezierFuncArray& bezierFuncs,
                      double tau)
{
  for(const BezierFuncPtr& func : bezierFuncs)
    func->setTau(tau);

  if(!solver->solve())
    std::cerr << "Could not solve at tau = " << tau << std::endl;

  return solver->getProblem()->getOptimalSolution();
}

class TrajectoryDisplayWorld : public osgDart::WorldNode
{
public:

  TrajectoryDisplayWorld(dart::simulation::WorldPtr world,
                         const std::vector<Eigen::VectorXd>& traj)
    : osgDart::WorldNode(world), mTrajectory(traj), count(0)
  {
    hubo = world->getSkeleton(0);
  }

  void customPreRefresh() override
  {
    if(mTrajectory.size() == 0 && count == 0)
    {
      std::cerr << "No trajectory was generated!" << std::endl;
      ++count;
      return;
    }

    hubo->setPositions(mTrajectory[count]);
    ++count;
    if(count >= mTrajectory.size())
      count = 0;
  }

protected:

  SkeletonPtr hubo;
  std::vector<Eigen::VectorXd> mTrajectory;
  size_t count;
};

int main()
{
  dart::utils::DartLoader loader;
  loader.addPackageDirectory("drchubo", DART_DATA_PATH"/urdf/drchubo");
  SkeletonPtr hubo =
      loader.parseSkeleton(DART_DATA_PATH"/urdf/drchubo/drchubo.urdf");

  std::shared_ptr<Solver> solver = std::make_shared<GradientDescentSolver>();
  std::shared_ptr<Problem> problem = std::make_shared<Problem>();
  BezierFuncArray bezierFuncs;

  Eigen::VectorXd alphas = make_alphas();
  std::vector<LinearCombo::Term> terms = make_terms();
  std::shared_ptr<LinearCombo> f = std::make_shared<LinearCombo>(terms, alphas);
  bezierFuncs.push_back(f);
  problem->addEqConstraint(f);

  JacobianNode* swing_foot = hubo->getBodyNode("Body_RAR");
  std::vector<Eigen::VectorXd> alphaArray;
  alphaArray.push_back(Eigen::VectorXd::Constant(4, 0.0));
  alphaArray.push_back(Eigen::VectorXd::Constant(4, 0.0));
  alphaArray.push_back(Eigen::VectorXd::Constant(4, 0.0));
  alphaArray.push_back(make_alphas());
  alphaArray.push_back(make_alphas());
  alphaArray.push_back(make_alphas());
  std::shared_ptr<EndEffectorConstraint> g =
      std::make_shared<EndEffectorConstraint>(swing_foot, alphaArray);
  g->mIK->getErrorMethod().setLinearBounds(
        Eigen::Vector3d::Constant(-std::numeric_limits<double>::infinity()),
        Eigen::Vector3d::Constant( std::numeric_limits<double>::infinity()));
  bezierFuncs.push_back(g);
  problem->addEqConstraint(g);

  JacobianNode* support_foot = hubo->getBodyNode("Body_LAR");
  InverseKinematicsPtr support = InverseKinematics::create(support_foot);
  problem->addEqConstraint(support->getProblem()->getEqConstraint(0));

  std::vector<Eigen::VectorXd> trajectory;
  double tau = 0.0;
  while(tau <= 1.0)
  {
    trajectory.push_back(solve(solver, bezierFuncs, tau));
    tau += hubo->getTimeStep();
  }

  dart::simulation::WorldPtr world = std::make_shared<dart::simulation::World>();
  world->addSkeleton(hubo);
  osg::ref_ptr<TrajectoryDisplayWorld> display = new TrajectoryDisplayWorld(world, trajectory);

  osgDart::Viewer viewer;
  viewer.addWorldNode(display);
  viewer.allowSimulation(false);

  viewer.run();
}
