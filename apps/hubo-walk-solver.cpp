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

#include <iostream>
#include <fstream>

#include <dart/dart.h>
#include <osgDart/osgDart.h>
#include <yaml-cpp/yaml.h>

#include <HuboPath/Operator.hpp>

#include <Eigen/QR>
#include <array>

#include <osg/Timer>

using namespace dart::dynamics;
using namespace dart::optimizer;

const double frequency = 1000.0;

//static double sign(double value)
//{
//  return value > 0? 1.0 : (value < 0? -1.0 : 0.0);
//}

class WalkSolver : public Solver
{
public:

  WalkSolver(const std::shared_ptr<Problem>& problem = nullptr)
    : Solver(problem)
  {
    // Do nothing
  }

  bool solve() override final
  {
    const std::shared_ptr<Problem>& problem = mProperties.mProblem;

    Eigen::VectorXd x = problem->getInitialGuess();
    Eigen::VectorXd dx = Eigen::VectorXd::Zero(x.size());
    Eigen::Map<Eigen::VectorXd> dxmap(dx.data(), dx.size());

    bool solved = false;
    mLastNumIterations = 0;

    while(!solved && mLastNumIterations < mProperties.mNumMaxIterations)
    {
      for(size_t i=0; i < problem->getNumEqConstraints(); ++i)
      {
        const FunctionPtr& constraint = problem->getEqConstraint(i);
//        double cost = constraint->eval(x);
        dxmap.setZero();
        constraint->evalGradient(x, dxmap);
        x -= dxmap/* * sign(cost)*/;
      }

      solved = true;
      for(size_t i=0; i < problem->getNumEqConstraints(); ++i)
      {
        const FunctionPtr& constraint = problem->getEqConstraint(i);
        double cost = constraint->eval(x);
        if(std::abs(cost) > mProperties.mTolerance)
        {
//          std::cout << "[" << i << "] " << cost << ": ";
          constraint->evalGradient(x, dxmap);

//          for(size_t j=0; j < dxmap.size(); ++j)
//          {
//            if(std::abs(dxmap[j]) > mProperties.mTolerance)
//              std::cout << "(" << j << ") " << dxmap[j] << ", ";
//          }
//          std::cout << " \t|\t ";
          solved = false;
        }
      }

//      if(!solved)
//        std::cout << std::endl;

      ++mLastNumIterations;
    }

    problem->setOptimalSolution(x);
    problem->setOptimumValue(0.0);
    return solved;
  }

  std::string getType() const override final
  {
    return "WalkSolver";
  }

  std::shared_ptr<Solver> clone() const override
  {
    return std::make_shared<WalkSolver>(mProperties.mProblem);
  }

  size_t getLastNumIterations() const
  {
    return mLastNumIterations;
  }

protected:

  size_t mLastNumIterations;

};

class HuboLegIK : public InverseKinematics::Analytical
{
public:

  /// baseLink should be Body_LHY or Body_RHY
  HuboLegIK(InverseKinematics* _ik, const std::string& baseLinkName)
    : Analytical(_ik, "HuboLegIK_"+baseLinkName, Analytical::Properties()),
      configured(false),
      mBaseLinkName(baseLinkName)
  {
    // Do nothing
  }

  std::unique_ptr<GradientMethod> clone(InverseKinematics* _newIK) const override
  {
    return std::unique_ptr<GradientMethod>(new HuboLegIK(_newIK, mBaseLinkName));
  }

  const std::vector<Solution>& computeSolutions(
      const Eigen::Isometry3d& _desiredBodyTf) override
  {
    mSolutions.clear();
    mSolutions.reserve(8);

    if(!configured)
    {
      configure();

      if(!configured)
      {
        dtwarn << "[HuboLegIK::computeSolutions] This analytical IK was not able "
              << "to configure properly, so it will not be able to compute "
              << "solutions\n";
        return mSolutions;
      }
    }

    const BodyNodePtr& base = mBaseLink.lock();
    if(nullptr == base)
    {
      dterr << "[HuboLegIK::computeSolutions] Attempting to perform IK on a "
            << "limb that no longer exists!\n";
      assert(false);
      return mSolutions;
    }

    double nx, ny, sx, sy, ax, ay, az, px, py, pz;
    double q1, q2, q3, q4, q5, q6;
    double S2, S4, S6;
    double C2, C4, C5, C6;
    double C45, psi, q345;
    std::complex<double> radical;
    std::complex<double> sqrt_radical;
    Eigen::Isometry3d B, Binv;

    Eigen::Vector6d testQ;

    B = (base->getParentBodyNode()->getWorldTransform() * waist).inverse()
        * _desiredBodyTf * footTfInv;
    Binv = B.inverse();

    nx = Binv(0,0); sx = Binv(0,1); ax = Binv(0,2); px = Binv(0,3);
    ny = Binv(1,0); sy = Binv(1,1); ay = Binv(1,2); py = Binv(1,3);
                                    az = Binv(2,2); pz = Binv(2,3);

    for(size_t i=0; i < 8; ++i)
    {
      bool isValid = true;

      C4 = ((px+L6)*(px+L6) - L4*L4 - L5*L5 + py*py + pz*pz)/(2*L4*L5);
      radical = 1-C4*C4;
      sqrt_radical = std::sqrt(radical);
      if(sqrt_radical.imag() != 0)
          isValid = false;
      q4 = atan2(alternatives(i,0)*sqrt_radical.real(), C4);

      S4 = sin(q4);
      psi = atan2(S4*L4, C4*L4+L5);
      radical = (px+L6)*(px+L6) + py*py;
      sqrt_radical = std::sqrt(radical);
      if(sqrt_radical.imag() != 0)
          isValid = false;

      q5 = dart::math::wrapToPi(atan2(-pz, alternatives(i,1)*sqrt_radical.real())-psi);

      q6 = atan2(py, -(px+L6));
      C45 = cos(q4+q5);
      C5 = cos(q5);
      if( C45*L4 + C5*L5 < 0 )
          q6 = dart::math::wrapToPi(q6+M_PI);

      S6 = sin(q6);
      C6 = cos(q6);

      S2 = C6*ay + S6*ax;
      radical = 1-S2*S2;
      sqrt_radical = std::sqrt(radical);
      if(sqrt_radical.imag() != 0)
          isValid = false;
      q2 = atan2(S2, alternatives(i,2)*sqrt_radical.real());

      q1 = atan2(C6*sy + S6*sx, C6*ny + S6*nx);
      C2 = cos(q2);
      if( C2 < 0 )
          q1 = dart::math::wrapToPi(q1+M_PI);

      q345 = atan2(-az/C2, -(C6*ax - S6*ay)/C2);
      q3 = dart::math::wrapToPi(q345 - q4 - q5);

      testQ[0]=q1; testQ[1]=q2; testQ[2]=q3; testQ[3]=q4; testQ[4]=q5; testQ[5]=q6;

      for(int k=0; k<testQ.size(); ++k)
          if( fabs(testQ[k]) < zeroSize )
              testQ[k] = 0;

      int validity = isValid? VALID : OUT_OF_REACH;
      mSolutions.push_back(Solution(testQ, validity));
    }

    checkSolutionJointLimits();

    return mSolutions;
  }

  const std::vector<size_t>& getDofs() const override
  {
    if(!configured)
      configure();

    return mDofs;
  }

  const double zeroSize = 1e-8;

protected:

  void configure() const
  {
    configured = false;

    mBaseLink = mIK->getNode()->getSkeleton()->getBodyNode(mBaseLinkName);

    BodyNode* base = mBaseLink.lock();
    if(nullptr == base)
    {
      dterr << "[HuboLegIK::configure] base link is a nullptr\n";
      assert(false);
      return;
    }

    const SkeletonPtr& skel = mIK->getNode()->getSkeleton();
    BodyNode* pelvis = skel->getBodyNode("Body_TSY");
    if(nullptr == pelvis)
    {
      dterr << "[HuboLegIK::configure] Could not find Hubo's pelvis "
            << "(Body_TSY)\n";
      assert(false);
      return;
    }

    Eigen::Vector6d saved_q;

    DegreeOfFreedom* dofs[6];
    BodyNode* bn = base;
    for(size_t i=0; i < 6; ++i)
    {
      Joint* joint = bn->getParentJoint();
      if(joint->getNumDofs() != 1)
      {
        dterr << "[HuboLegIK::configure] Invalid number of DOFs ("
              << joint->getNumDofs() << ") in the Joint [" << joint->getName()
              << "]\n";
        assert(false);
        return;
      }

      dofs[i] = joint->getDof(0);
      saved_q[i] = dofs[i]->getPosition();
      dofs[i]->setPosition(0.0);

      if(bn->getNumChildBodyNodes() > 0)
        bn = bn->getChildBodyNode(0);
    }

    // Thigh length
    L4 = std::abs(dofs[3]->getChildBodyNode()->
        getRelativeTransform().translation()[2]);

    // Calf length
    L5 = std::abs(dofs[4]->getChildBodyNode()->
        getRelativeTransform().translation()[2]);

    // This offset will be taken care of with footTfInv
    L6 = 0.0;

    hipRotation = Eigen::Isometry3d::Identity();
    hipRotation.rotate(Eigen::AngleAxisd(90*M_PI/180.0,
                                         Eigen::Vector3d::UnitZ()));

    waist = dofs[2]->getChildBodyNode()->getTransform(
              dofs[0]->getParentBodyNode()) * hipRotation;

    footTfInv = Eigen::Isometry3d::Identity();
    footTfInv.rotate(Eigen::AngleAxisd(-90*M_PI/180.0,
                                       Eigen::Vector3d::UnitY()));
    footTfInv = footTfInv * mIK->getNode()->getTransform(dofs[5]->getChildBodyNode());
    footTfInv = footTfInv.inverse();

    alternatives <<
         1,  1,  1,
         1,  1, -1,
         1, -1,  1,
         1, -1, -1,
        -1,  1,  1,
        -1,  1, -1,
        -1, -1,  1,
        -1, -1, -1;

    for(size_t i=0; i < 6; ++i)
    {
      dofs[i]->setPosition(saved_q[i]);
      mDofs.push_back(dofs[i]->getIndexInSkeleton());
    }

    configured = true;
  }

  mutable double L4, L5, L6;
  mutable Eigen::Isometry3d waist;
  mutable Eigen::Isometry3d hipRotation;
  mutable Eigen::Isometry3d footTfInv;
  mutable Eigen::Matrix<int, 8, 3> alternatives;

  mutable std::vector<size_t> mDofs;

  mutable bool configured;

  std::string mBaseLinkName;

  mutable WeakBodyNodePtr mBaseLink;

};

size_t factorial(size_t m)
{
  size_t result = 1u;
  for(size_t i=2; i <= m; ++i)
    result *= i;
  return result;
}

double computeBezier(double tau, const Eigen::VectorXd& alpha, double M = 4)
{
  if(alpha.size() != M+1)
    std::cerr << "Invalid number of alphas: " << alpha.size() << std::endl;
  assert(alpha.size() == M+1);
  double result = 0.0;
  for(size_t k = 0; k <= M; ++k)
  {
    result += alpha[k] * factorial(M) / (factorial(k) * factorial(M-k))
        * pow(tau, k) * pow(1-tau, M-k);
  }

  return result;
}

class BezierDependentFunc
{
public:

  virtual void setTau(double tau)
  {
    mTau = tau;
  }

protected:

  double mTau;
};

typedef std::shared_ptr<BezierDependentFunc> BezierFuncPtr;
typedef std::vector<BezierFuncPtr> BezierFuncArray;

struct Term
{
  Term(size_t _index = INVALID_INDEX, double _coeff = 0.0)
    : index(_index), coeff(_coeff) { }

  size_t index;
  double coeff;
};

class LinearCombo : public Function,
                    public BezierDependentFunc
{
public:

  LinearCombo(const std::vector<Term>& terms,
              const Eigen::VectorXd& alphas)
    : mTerms(terms), mAlphas(alphas)
  {
    assert(mAlphas.size() == 5);
  }

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

template <int eqns, int dofs>
class LinearComboSystem : public Function,
                          public BezierDependentFunc
{
public:

  LinearComboSystem()
  {
    A.setZero();
  }

  void recomputeNullspace()
  {
    dart::math::computeNullSpace(A, NS);
    if(NS.rows() > 0 && NS.cols() > 0)
      Anull = NS*NS.transpose();
    else
      Anull.resize(0,0);
  }

  void setEqn(size_t row,
              const std::vector<Term>& terms,
              const Eigen::VectorXd& alphas)
  {
    assert(row < eqns);

    for(size_t i=0; i < terms.size(); ++i)
    {
      A(row, terms[i].index) = terms[i].coeff;
    }

    mAlphas[row] = alphas;

    recomputeNullspace();
  }

  void setTauSatisfaction(const std::vector<Term>& terms,
                          double p_plus_, double p_minus_)
  {
    p_plus = p_plus_;
    p_minus = p_minus_;

    for(size_t i=0; i < terms.size(); ++i)
    {
      A(0, terms[i].index) = terms[i].coeff;
    }

    recomputeNullspace();
  }

  void setTau(double tau) override
  {
    mTau = tau;

    computeX();
  }

  void computeX()
  {
//    std::cout << "computing x" << std::endl;
    b[0] = p_plus + mTau*(p_minus - p_plus);

    for(int i=1; i < eqns; ++i)
    {
      b[i] = computeBezier(mTau, mAlphas[i]);
    }

    x = A.colPivHouseholderQr().solve(b);
  }

  double eval(const Eigen::VectorXd& _x) override final
  {
    Eigen::Map<Eigen::VectorXd> dxmap(dx.data(), dofs);
    evalGradient(_x, dxmap);
    return dxmap.norm();
  }

  void evalGradient(const Eigen::VectorXd& _x,
                    Eigen::Map<Eigen::VectorXd> _grad) override final
  {
    _grad = _x - x;
    if(Anull.cols() > 0 && Anull.rows() > 0)
      _grad -= Anull*_grad;
  }

protected:

  Eigen::Matrix<double, eqns, dofs> A;
  Eigen::Matrix<double, eqns, 1> b;
  Eigen::Matrix<double, dofs, 1> x;
  Eigen::Matrix<double, dofs, 1> dx;

  double p_plus;
  double p_minus;

  Eigen::Matrix<double, dofs, Eigen::Dynamic> NS;
  Eigen::MatrixXd Anull;

  std::array<Eigen::VectorXd, eqns> mAlphas;
};

class Holonomic : public Function
{
public:

  Holonomic(BodyNode* node)
    : bn(node)
  {
    skel = bn->getSkeleton();
    mBnTf = Eigen::Isometry3d::Identity();
  }

  void setTargetTransform(const Eigen::Isometry3d& tf)
  {
    mBnTf = tf;
  }

  double eval(const Eigen::VectorXd& _x) override final
  {
    skel->setPositions(_x);
    const Eigen::Isometry3d& currentBnTf = bn->getWorldTransform();
    const Eigen::Isometry3d& currentRootTf = skel->getBodyNode(0)->getWorldTransform();

    const Eigen::Isometry3d& deltaTf = mBnTf * currentBnTf.inverse();
    const Eigen::Isometry3d& newRootTf = deltaTf * currentRootTf;

    const Eigen::Vector6d& newRootPos = FreeJoint::convertToPositions(newRootTf);
    rootGrad = _x.head<6>() - newRootPos;

    return rootGrad.norm();
  }

  void evalGradient(const Eigen::VectorXd& /*_x*/,
                    Eigen::Map<Eigen::VectorXd> _grad) override final
  {
    _grad.setZero();

    // Note: This function relies on the fact that GradientDescentSolver calls
    // eval() before it calls evalGradient()
    _grad.head<6>() = rootGrad;
  }

protected:

  SkeletonPtr skel;

  BodyNode* bn;

  Eigen::Isometry3d mBnTf;

  Eigen::Vector6d rootGrad;
};

class EndEffectorConstraint : public Function,
                              public BezierDependentFunc
{
public:

  /// alphas should have six entries, one corresponding to each component of the
  /// EndEffector pose. Unused alphas can be set to all zeros.
  EndEffectorConstraint(JacobianNode* node,
                        const std::vector<Eigen::VectorXd>& alphas,
                        const std::vector<size_t>& ignoreDofs = std::vector<size_t>())
    : mAlphaArray(alphas), mIgnoreDofs(ignoreDofs)
  {
    for(size_t i=0; i < mAlphaArray.size(); ++i)
    {
      if(mAlphaArray[i].size() != 5)
        std::cerr << "EndEffectorConstraint: Invalid alpha size for index "
                  << i << ": " << mAlphaArray[i].size() << std::endl;
      assert(mAlphaArray[i].size() == 5);
    }

    mIK = InverseKinematics::create(node);

    if(node->getName() == "Body_LAR")
    {
      mIK->setGradientMethod<HuboLegIK>("Body_LHY");
    }
    else if(node->getName() == "Body_RAR")
    {
      mIK->setGradientMethod<HuboLegIK>("Body_RHY");
    }
    else
    {
      dterr << "[EndEffectorConstraint::constructor] Unsupported node: "
            << node->getName() << "\n";
      assert(false);
    }

    const SkeletonPtr& skel = node->getSkeleton();

    std::vector<size_t> indices;
    indices.reserve(skel->getNumDofs());
    for(size_t i=0; i < skel->getNumDofs(); ++i)
      indices.push_back(i);

    mIK->setDofs(indices);
  }

  void setTau(double tau) override final
  {
    mTau = tau;

    mIK->clearCaches();
  }

  double eval(const Eigen::VectorXd& _x) override final
  {
    Eigen::Vector6d screw;
    for(size_t i=0; i < 6; ++i)
      screw[i] = computeBezier(mTau, mAlphaArray[i]);

    Eigen::Isometry3d tf(Eigen::Isometry3d::Identity());
    // TODO: Find out how to convert the screw into an Isometry3d
    tf.translate(screw.head<3>());
    for(size_t i=0; i < 3; ++i)
    {
      double angle = screw[i+3];
      Eigen::Vector3d axis(Eigen::Vector3d::Zero());
      axis[i] = 1.0;
      tf.rotate(Eigen::AngleAxisd(angle, axis));
    }

    mIK->getTarget()->setRelativeTransform(tf);

    const Eigen::Vector6d& error = mIK->getErrorMethod().evalError(_x);
    return error.norm();
  }

  void evalGradient(const Eigen::VectorXd& _x,
                    Eigen::Map<Eigen::VectorXd> _grad) override final
  {
    mIK->getGradientMethod().evalGradient(_x, _grad);

    for(size_t dof : mIgnoreDofs)
    {
      _grad[dof] = 0.0;
    }
  }

  InverseKinematicsPtr mIK;

protected:

  std::vector<Eigen::VectorXd> mAlphaArray;
  std::vector<size_t> mIgnoreDofs;

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

void build_terms(const SkeletonPtr& /*hubo*/,
                 std::vector<Term>& /*terms*/)
{
  // Terminate template recursion
}

template <typename ... Args>
void build_terms(const SkeletonPtr& hubo, std::vector<Term>& terms,
                 double coeff, const std::string& name, Args... args)
{
  terms.push_back(Term(hubo->getDof(name)->getIndexInSkeleton(), coeff));
  build_terms(hubo, terms, args...);
}

template <typename ... Args>
void build_terms(const SkeletonPtr& hubo, std::vector<Term>& terms,
                 double coeff, size_t index, Args... args)
{
  terms.push_back(Term(index, coeff));
  build_terms(hubo, terms, args...);
}

template <typename ... Args>
std::vector<Term> make_terms(const SkeletonPtr& hubo, Args... args)
{
  std::vector<Term> terms;
  build_terms(hubo, terms, args...);
  return terms;
}

Eigen::VectorXd solve(const std::shared_ptr<Solver>& solver,
                      const BezierFuncArray& bezierFuncs,
                      double tau)
{
  for(const BezierFuncPtr& func : bezierFuncs)
    func->setTau(tau);

  if(!solver->solve())
  {
    std::cerr << "Could not solve at tau = " << tau << " in "
              << solver->getNumMaxIterations() << " steps" << std::endl;

    const std::shared_ptr<Problem> problem = solver->getProblem();
    for(size_t i=0; i < problem->getNumEqConstraints(); ++i)
    {
      const FunctionPtr& f = problem->getEqConstraint(i);
      double cost = f->eval(problem->getOptimalSolution());
      if(std::abs(cost) > solver->getTolerance())
        std::cerr << "Cost of (" << i+1 << "): " << cost << std::endl;
    }
  }

  return solver->getProblem()->getOptimalSolution();
}

class TrajectoryDisplayWorld : public osgDart::WorldNode
{
public:

  TrajectoryDisplayWorld(dart::simulation::WorldPtr world,
                         const std::vector<Eigen::VectorXd>& traj)
    : osgDart::WorldNode(world), mTrajectory(traj)
  {
    hubo = world->getSkeleton(0);
    LSR = hubo->getDof("LSR")->getIndexInSkeleton();
    RSR = hubo->getDof("RSR")->getIndexInSkeleton();

    firstLoop = true;
  }

  void customPreRefresh() override
  {
    if(mTrajectory.size() == 0)
    {
      if(firstLoop)
      {
        std::cerr << "No trajectory was generated!" << std::endl;
        firstLoop = false;
      }
      return;
    }

    if(firstLoop)
    {
      mTimer.setStartTick();
      firstLoop = false;
    }

    double time = mTimer.time_s();
    size_t count = (size_t)(floor(frequency*time)) % mTrajectory.size();

    double dt = 1.0/frequency;

    size_t last_1 = count <= 0? 0 : count - 1;
    size_t last_2 = count <= 1? 0 : count - 2;
    size_t last_3 = count <= 2? 0 : count - 3;
    size_t last_4 = count <= 3? 0 : count - 4;

    size_t next_1 = count < mTrajectory.size() - 1? count + 1 : count;
    size_t next_2 = count < mTrajectory.size() - 2? count + 2 : count;
    size_t next_3 = count < mTrajectory.size() - 3? count + 3 : count;
    size_t next_4 = count < mTrajectory.size() - 4? count + 4 : count;

    Eigen::VectorXd velocities =
        (0.5*mTrajectory[next_1] - 0.5*mTrajectory[last_1])/dt;

    Eigen::Matrix<double, 5, 1> c;
    c << -205.0/72.0, 8.0/5.0, -1.0/5.0, 8.0/315.0, -1.0/560.0;

    Eigen::VectorXd accelerations =
        (  c[1]*mTrajectory[next_1] + c[2]*mTrajectory[next_2] + c[3]*mTrajectory[next_3] + c[4]*mTrajectory[next_4]
         + c[0]*mTrajectory[count]
         + c[1]*mTrajectory[last_1] + c[2]*mTrajectory[last_2] + c[3]*mTrajectory[last_3] + c[4]*mTrajectory[last_4])/pow(dt,2);


//    accelerations.head<6>().setZero();

    hubo->setPositions(mTrajectory[count]);
    hubo->setVelocities(velocities);
    hubo->setAccelerations(accelerations);

//    hubo->setPositions(mMapping, positions);
//    clone->setPositions(mMapping, mRaw[count]);
  }

protected:

  SkeletonPtr hubo;
  std::vector<Eigen::VectorXd> mTrajectory;
  size_t LSR;
  size_t RSR;

  bool firstLoop;
  osg::Timer mTimer;
};

Eigen::VectorXd getAlphas(const YAML::Node& a, size_t index)
{
  const YAML::Node entry = a[index];
  Eigen::VectorXd alphas(entry.size());
  for(size_t i=0; i < entry.size(); ++i)
    alphas[i] = entry[i].as<double>();

  return alphas;
}

#define ADD_LINEAR_OUTPUT( index, T ) \
{\
  alphas = getAlphas(a, index-1);\
  terms = T;\
  system->setEqn(index, terms, alphas);\
}

#define ADD_DOUBLE_SUPPORT_OUTPUT( index1, index2, index3, index4, index5, index6, body, relativeTo )\
{\
  JacobianNode* ee = body;\
  std::vector<Eigen::VectorXd> alphaArray;\
  alphaArray.push_back(getAlphas(a, index1-1));\
  alphaArray.push_back(getAlphas(a, index2-1));\
  alphaArray.push_back(getAlphas(a, index3-1));\
  alphaArray.push_back(getAlphas(a, index4-1));\
  alphaArray.push_back(getAlphas(a, index5-1));\
  std::shared_ptr<EndEffectorConstraint> f =\
      std::make_shared<EndEffectorConstraint>(ee, alphaArray);\
  f->mIK->getTarget()->setParentFrame(relativeTo);\
  bezierFuncs.push_back(f);\
  problem->addEqConstraint(f);\
}


#define ADD_SWING_FOOT_ORIENTATION_OUTPUT( index1, index2, index3, body )\
{\
  JacobianNode* ee = body;\
  std::vector<Eigen::VectorXd> alphaArray;\
  alphaArray.push_back(Eigen::VectorXd::Constant(5, 0.0));\
  alphaArray.push_back(Eigen::VectorXd::Constant(5, 0.0));\
  alphaArray.push_back(Eigen::VectorXd::Constant(5, 0.0));\
  alphaArray.push_back(getAlphas(a, index1-1));\
  alphaArray.push_back(getAlphas(a, index2-1));\
  alphaArray.push_back(getAlphas(a, index3-1));\
  std::shared_ptr<EndEffectorConstraint> f =\
      std::make_shared<EndEffectorConstraint>(ee, alphaArray);\
  f->mIK->getErrorMethod().setLinearBounds(\
      Eigen::Vector3d::Constant(-std::numeric_limits<double>::infinity()),\
      Eigen::Vector3d::Constant( std::numeric_limits<double>::infinity()));\
  bezierFuncs.push_back(f);\
  problem->addEqConstraint(f);\
}

#define ZERO_OUT_ALPHAS(index) \
  for(size_t i=0; i < a[index-1].size(); ++i)\
    a[index-1][i] = 0.0;

double computeP(const SkeletonPtr& hubo, const std::string& st,
                double thighLength, double calfLength)
{
//  std::cout << "coeffs for actual: "
//            << "(" << hubo->getDof(st+"AP")->getIndexInSkeleton() << ") " << -calfLength << " \t|\t "
//            << "(" << hubo->getDof(st+"")
  return calfLength*(-hubo->getDof(st+"AP")->getPosition())
      + thighLength*(-hubo->getDof(st+"AP")->getPosition()
                     -hubo->getDof(st+"KP")->getPosition());
}

double computeP(double t, double p0, double pdot0, double vd)
{
  double eps = 1.0;

  return vd*t + p0 + ((1.0-exp(-eps*t))*pdot0 + (exp(-eps*t)-1.0)*vd)/eps;
}

std::vector<Eigen::VectorXd> setupAndSolveProblem(
    const SkeletonPtr& hubo, YAML::Node params,
    const std::string& st/*ance*/, const std::string& sw/*ing*/,
    bool double_support = false)
{
  Eigen::VectorXd lastPositions = hubo->getPositions();
  hubo->resetPositions();

  double p_plus = params["p"][1].as<double>();
  double p_minus = params["p"][0].as<double>();

  Eigen::VectorXd x_plus(2*hubo->getNumDofs()), x_minus(2*hubo->getNumDofs());
  Eigen::VectorXd q_plus(hubo->getNumDofs()), q_minus(hubo->getNumDofs());
  YAML::Node node_x_plus = params["x_plus"];
  for(size_t i=0; i < node_x_plus.size(); ++i)
    x_plus[i] = node_x_plus[i].as<double>();
  YAML::Node node_x_minus = params["x_minus"];
  for(size_t i=0; i < node_x_minus.size(); ++i)
    x_minus[i] = node_x_minus[i].as<double>();

  q_plus = x_plus.block(0,0,hubo->getNumDofs(),1);
  q_minus = x_minus.block(0,0,hubo->getNumDofs(),1);

//  YAML::Node p = params["p"];
//  std::cout << "p: ";
//  for(size_t i=0; i < p.size(); ++i)
//    std::cout << p[i].as<double>() << "\t";
//  std::cout << "\n\n";

  YAML::Node a = params["a"];
//  std::cout << "a:\n";
//  for(size_t i=0; i < a.size(); ++i)
//  {
//    std::cout << i << ") ";
//    YAML::Node entry = a[i];
//    for(size_t j=0; j < entry.size(); ++j)
//      std::cout << entry[j].as<double>() << "\t";
//    std::cout << "\n";
//  }
//  std::cout << "\n";

//  std::shared_ptr<GradientDescentSolver> solver = std::make_shared<GradientDescentSolver>();
//  solver->setStepSize(1.0);
  std::shared_ptr<WalkSolver> solver = std::make_shared<WalkSolver>();
  std::shared_ptr<Problem> problem = std::make_shared<Problem>();
  problem->setDimension(hubo->getNumDofs());
  solver->setProblem(problem);
  BezierFuncArray bezierFuncs;
  Eigen::VectorXd alphas;
  std::vector<Term> terms;

  Eigen::VectorXd lower(hubo->getNumDofs());
  Eigen::VectorXd upper(hubo->getNumDofs());

  for(size_t i=0; i < hubo->getNumDofs(); ++i)
  {
    lower[i] = hubo->getDof(i)->getPositionLowerLimit();
    upper[i] = hubo->getDof(i)->getPositionUpperLimit();
  }

  problem->setLowerBounds(lower);
  problem->setUpperBounds(upper);

  std::shared_ptr<LinearComboSystem<24,33>> system =
      std::make_shared<LinearComboSystem<24,33>>();
  problem->addEqConstraint(system);
  bezierFuncs.push_back(system);

  ADD_LINEAR_OUTPUT(1,  make_terms(hubo, 1.0, st+"KP"));
  ADD_LINEAR_OUTPUT(2,  make_terms(hubo, -1.0, st+"AP", -1.0, st+"KP", -1.0, st+"HP"));
  ADD_LINEAR_OUTPUT(3,  make_terms(hubo, 1.0, st+"AR"));
  ADD_LINEAR_OUTPUT(4,  make_terms(hubo, -1.0, st+"AR", -1.0, st+"HR"));
  ADD_LINEAR_OUTPUT(5,  make_terms(hubo, 1.0, st+"HY"));
  ADD_LINEAR_OUTPUT(6,  make_terms(hubo, 1.0, st+"SP"));
  ADD_LINEAR_OUTPUT(7,  make_terms(hubo, 1.0, st+"SR"));
  ADD_LINEAR_OUTPUT(8,  make_terms(hubo, 1.0, st+"SY"));
  ADD_LINEAR_OUTPUT(9,  make_terms(hubo, 1.0, st+"EP"));
  ADD_LINEAR_OUTPUT(10, make_terms(hubo, 1.0, st+"WY"));
  ADD_LINEAR_OUTPUT(11, make_terms(hubo, 1.0, st+"WP"));
  ADD_LINEAR_OUTPUT(12, make_terms(hubo, 1.0, st+"WR"));
  ADD_LINEAR_OUTPUT(13, make_terms(hubo, 1.0, "TSY"));
  ADD_LINEAR_OUTPUT(14, make_terms(hubo, 1.0, sw+"SP"));
  ADD_LINEAR_OUTPUT(15, make_terms(hubo, 1.0, sw+"SR"));
  ADD_LINEAR_OUTPUT(16, make_terms(hubo, 1.0, sw+"SY"));
  ADD_LINEAR_OUTPUT(17, make_terms(hubo, 1.0, sw+"EP"));
  ADD_LINEAR_OUTPUT(18, make_terms(hubo, 1.0, sw+"WY"));
  ADD_LINEAR_OUTPUT(19, make_terms(hubo, 1.0, sw+"WP"));
  ADD_LINEAR_OUTPUT(20, make_terms(hubo, 1.0, sw+"WR"));

//  BodyNode* knee = hubo->getBodyNode("Body_"+st+"KP");
//  double thighLength = std::abs(knee->getRelativeTransform().translation()[2]);
//  double thighLength = std::abs(knee->getTransform(knee->getParentBodyNode()->getParentBodyNode()->getParentBodyNode()).translation()[2]);
  double hLength = 0.1401;
  double thighLength = 0.3299;
//  std::cout << "thigh: " << thighLength << std::endl;
//  BodyNode* ankle = hubo->getBodyNode("Body_"+st+"AP");
//  double calfLength = std::abs(ankle->getRelativeTransform().translation()[2]);
  double calfLength = 0.33;

//  std::cout << "calf: " << calfLength << std::endl;
//  double legRatio = calfLength/(thighLength+calfLength);
  double legRatio = calfLength/(thighLength+calfLength+hLength);

  std::shared_ptr<Holonomic> h = std::make_shared<Holonomic>(
        hubo->getBodyNode("Body_"+st+"AR"));
  problem->addEqConstraint(h);

  if(double_support)
  {
    ADD_DOUBLE_SUPPORT_OUTPUT(21, 22, 23, 24, 25, 26,
                              hubo->getBodyNode("Body_"+sw+"AR"),
                              hubo->getBodyNode("Body_"+st+"AR"));
  }
  else
  {
    ADD_LINEAR_OUTPUT(21, make_terms(hubo, 1.0, sw+"KP"));
    ADD_LINEAR_OUTPUT(22, make_terms(hubo, -1.0, st+"AP", -1.0, st+"KP", -1.0, st+"HP",
                                     1.0, sw+"HP", legRatio, sw+"KP"));

    ADD_LINEAR_OUTPUT(23, make_terms(hubo, 1.0, st+"HR", -1.0, sw+"HR"));

    ADD_SWING_FOOT_ORIENTATION_OUTPUT( 24, 25, 26, hubo->getBodyNode("Body_"+sw+"AR") );
  }

  system->setTauSatisfaction(
        make_terms(hubo, -calfLength-thighLength, st+"AP",
                                    -thighLength, st+"KP"),
        p_plus, p_minus);

  std::cout << "thigh: " << thighLength << "\n" << "calf: " << calfLength << std::endl;


  solver->setNumMaxIterations(1000);
  solver->setTolerance(1e-6);

//  hubo->setPositions(q_plus);
  hubo->setPositions(lastPositions);
  h->setTargetTransform(hubo->getBodyNode("Body_"+st+"AR")->getWorldTransform());

  std::vector<Eigen::VectorXd> trajectory;

  hubo->setPositions(lastPositions);

  double p0 = params["p"][1].as<double>();
//  double p0 = params["p0"].as<double>();

  YAML::Node pdot0_node = params["pdot0"];
  double pdot0 = pdot0_node? pdot0_node.as<double>() : 0.0;
  double vd = params["v"].as<double>();
  std::cout << "\n\n" << st << " Foot Trajectory (p0 " << p0 << ") :\n";
  double time = 0.0;
  double tau = 0.0;
  do
  {
    double p = computeP(time, p0, pdot0, vd);
    tau = (p - p_plus)/(p_minus - p_plus);
    if(tau > 1.0)
    {
      std::cout << "stopping" << std::endl;
      break;
    }

    trajectory.push_back(solve(solver, bezierFuncs, tau));

//    double actualP = computeP(hubo, st, thighLength, calfLength);
//    double actualTau = (actualP - p_plus)/(p_minus - p_plus);
//    std::cout << "tau: " << tau << " \t actualTau: " << actualTau
//              << " \t p: " << p << " \t actualP: " << actualP << std::endl;

    time += hubo->getTimeStep();

//    std::cout << "solved in " << solver->getLastNumIterations() << " steps" << std::endl;
  } while(tau <= 1.0);

  trajectory.push_back(solve(solver, bezierFuncs, tau));

  return trajectory;
}

SkeletonPtr createHubo()
{
  dart::utils::DartLoader loader;
  loader.addPackageDirectory("drchubo", DART_DATA_PATH"/urdf/drchubo");

  SkeletonPtr hubo =
      loader.parseSkeleton(DART_DATA_PATH"/urdf/drchubo/drchubo.urdf");

  for(size_t i = 0; i < hubo->getNumBodyNodes(); ++i)
  {
    BodyNode* bn = hubo->getBodyNode(i);
    if(bn->getName().substr(0, 7) == "Body_LF"
       || bn->getName().substr(0, 7) == "Body_RF"
       || bn->getName().substr(0, 7) == "Body_NK")
    {
      bn->remove();
      --i;
    }
  }

  hubo->getDof("REP")->setPositionUpperLimit(0.0);
  hubo->getDof("LEP")->setPositionUpperLimit(0.0);

//  std::cout << "REP: " << hubo->getDof("REP")->getPositionLowerLimit()
//            << " : " << hubo->getDof("REP")->getPositionUpperLimit() << std::endl;

//  std::cout << "LEP: " << hubo->getDof("LEP")->getPositionLowerLimit()
//            << " : " << hubo->getDof("LEP")->getPositionUpperLimit() << std::endl;

  return hubo;
}

void dumpTrajectory(const std::vector<Eigen::VectorXd>& traj,
                    const std::string& file)
{
  std::cout << "dumping trajectory of size " << traj.size() << " to " << file << std::endl;
  std::ofstream dump;
  dump.open(file, std::ofstream::out);
  for(size_t i=0; i < traj.size(); ++i)
  {
    const Eigen::VectorXd& pos = traj[i];
    for(int j=0; j < pos.size(); ++j)
      dump << pos[j] << "\t";
    dump << "\n";
  }
  dump.close();

}

int main()
{
  dart::simulation::WorldPtr world = std::make_shared<dart::simulation::World>();
  SkeletonPtr hubo = createHubo();
  world->addSkeleton(hubo);
  world->setTimeStep(1.0/frequency);

//  std::string file = "/home/grey/projects/protoHuboGUI/params_2015-08-27T07-01-0400.yaml";
//  std::string yaml = "/home/grey/projects/protoHuboGUI/params_2015-08-29T16-11-0400.yaml";
  std::string yaml = "/home/grey/projects/protoHuboGUI/params_2015-09-01T15-35-0400.yaml";

  bool loadfile = false;
//  loadfile = true;

  std::string dump_name = "/home/grey/projects/protoHuboGUI/trajectory.dat";

  std::vector<Eigen::VectorXd> raw_trajectory;
  if(loadfile)
  {
    std::ifstream file;
    file.open(dump_name);
    if(file.is_open())
    {
      while(!file.eof())
      {
        raw_trajectory.push_back(Eigen::VectorXd(hubo->getNumDofs()));
        Eigen::VectorXd& q = raw_trajectory.back();
        for(size_t i=0; i < hubo->getNumDofs(); ++i)
          file >> q[i];
      }
    }
    else
    {
      std::cerr << "Could not open file: " << dump_name << std::endl;
    }
  }
  else
  {
    // Start by testing left stance leg
    YAML::Node config = YAML::LoadFile(yaml);
    YAML::Node domain = config["domain"];

    YAML::Node leftStartParams, leftWalkParams, rightStartParams, rightWalkParams;
    for(size_t i=0; i < domain.size(); ++i)
    {
      const std::string& paramName = domain[i]["name"].as<std::string>();

      if("Left3DFlatStarting" == paramName)
      {
        std::cout << "Found left starting params" << std::endl;
        leftStartParams  = domain[i];
      }
      else if("Left3DFlatWalking" == paramName)
      {
        std::cout << "Found left walking params" << std::endl;
        leftWalkParams   = domain[i];
      }
      else if("Right3DFlatStarting" == paramName)
      {
        std::cout << "Found right starting params" << std::endl;
        rightStartParams = domain[i];
      }
      else if("Right3DFlatWalking" == paramName)
      {
        std::cout << "Found right walking params" << std::endl;
        rightWalkParams  = domain[i];
      }
      else
        std::cerr << "Unknown parameters: " << domain[i]["name"] << std::endl;
    }

    osg::Timer timer;
    timer.setStartTick();

    bool startWithLeft = true;
    startWithLeft = false;

    std::vector<Eigen::VectorXd> leftStart;
    if(startWithLeft)
      leftStart = setupAndSolveProblem(hubo, leftStartParams, "L", "R", true);

    std::vector<Eigen::VectorXd> rightWalk =
        setupAndSolveProblem(hubo, rightWalkParams, "R", "L");
    std::vector<Eigen::VectorXd> leftWalk =
        setupAndSolveProblem(hubo, leftWalkParams, "L", "R");

    std::cout << "Computation Time: " << timer.time_s() << std::endl;


    if(startWithLeft)
    {
      for(const Eigen::VectorXd& pos : leftStart)
        raw_trajectory.push_back(pos);
    }

    for(const Eigen::VectorXd& pos : rightWalk)
      raw_trajectory.push_back(pos);
    for(const Eigen::VectorXd& pos : leftWalk)
      raw_trajectory.push_back(pos);

//    Eigen::VectorXd diff = leftStart[0] - leftWalk[0]

    std::cout << "Trajectory Time:  " << (double)(raw_trajectory.size())*hubo->getTimeStep() << std::endl;

    dumpTrajectory(raw_trajectory, dump_name);
  }


  bool operate = false;
//  operate = true;

  if(operate)
  {
    HuboPath::Operator mOperator;
    std::vector<std::string> indexNames;
    std::vector<size_t> mOperatorIndices;
    for(size_t i=6; i < hubo->getNumDofs(); ++i)
    {
      DegreeOfFreedom* dof = hubo->getDof(i);
      mOperatorIndices.push_back(i);
      indexNames.push_back(dof->getName());

      dof->setPosition(raw_trajectory[0][i]);
    }
    mOperator.setJointIndices(indexNames);

    for(size_t i=0; i < raw_trajectory.size(); ++i)
    {
      hubo->setPositions(raw_trajectory[i]);
      mOperator.addWaypoint(hubo->getPositions(mOperatorIndices));
    }

    mOperator.setInterpolationMode(HUBO_PATH_RAW);

    HuboPath::Trajectory traj = mOperator.getCurrentTrajectory();
    traj.elements.erase(traj.elements.begin());
    if(!traj.check_limits())
      std::cout << "Limits violated" << std::endl;
    else
      std::cout << "Limits are okay" << std::endl;


    std::cout << "trajectory size: " << mOperator.getWaypoints().size() << std::endl;
    mOperator.sendNewTrajectory();
  }
  else
  {
    osg::ref_ptr<TrajectoryDisplayWorld> display =
        new TrajectoryDisplayWorld(world, raw_trajectory);

    osgDart::Viewer viewer;
    viewer.addWorldNode(display);
    viewer.allowSimulation(false);
//    viewer.record("/home/grey/dump/");

//    osg::ref_ptr<osgDart::SupportPolygonVisual> support =
//        new osgDart::SupportPolygonVisual(hubo);
    viewer.addAttachment(new osgDart::SupportPolygonVisual(hubo, -0.97+0.02));

    viewer.setUpViewInWindow(0, 0, 1280, 960);

    // Set up the default viewing position
    viewer.getCameraManipulator()->setHomePosition(osg::Vec3( 5.34,  3.00, 1.00),
                                                   osg::Vec3( 0.00,  0.00, 0.00),
                                                   osg::Vec3(-0.20, -0.08, 0.98));

    // Reset the camera manipulator so that it starts in the new viewing position
    viewer.setCameraManipulator(viewer.getCameraManipulator());

    viewer.run();
  }
}
