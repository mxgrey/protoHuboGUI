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

using namespace dart::dynamics;
using namespace dart::optimizer;

class HuboLegIK : public InverseKinematics::Analytical
{
public:

  /// baseLink should be Body_LHY or Body_RHY
  HuboLegIK(InverseKinematics* _ik, const std::string& baseLinkName)
    : Analytical(_ik, "HuboLegIK_"+baseLinkName),
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

  void setTau(double tau)
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

    dart::math::computeNullSpace(A, NS);
    if(NS.rows() > 0 && NS.cols() > 0)
      Anull = NS*NS.transpose();
    else
      Anull.resize(0,0);
  }

  void computeX()
  {
    for(int i=0; i < eqns; ++i)
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
    computeX();

    _grad = _x - x;
    if(Anull.cols() > 0 && Anull.rows() > 0)
      _grad -= Anull*_grad;
//    std::cout << "system: (10) " << _x[10] << "," << x[10] << " | (11) " << _x[11] << "," << x[11] << std::endl;
  }

protected:

  Eigen::Matrix<double, eqns, dofs> A;
  Eigen::Matrix<double, eqns, 1> b;
  Eigen::Matrix<double, dofs, 1> x;
  Eigen::Matrix<double, dofs, 1> dx;

  Eigen::Matrix<double, dofs, Eigen::Dynamic> NS;
  Eigen::MatrixXd Anull;

  std::array<Eigen::VectorXd, eqns> mAlphas;
};

class SatisfyTau : public Function,
                   public BezierDependentFunc
{
public:

  SatisfyTau(const std::vector<Term>& p_terms_,
             double p_plus_, double p_minus_)

    : p_terms(p_terms_),
      p_plus(p_plus_),
      p_minus(p_minus_)
  {
    // Do nothing
  }

  double computeP(const Eigen::VectorXd& _x)
  {
    double p = 0;
    for(const Term& term : p_terms)
      p += term.coeff * _x[term.index];

    return p;
  }

  double computeCost(const Eigen::VectorXd& _x)
  {
    double p = computeP(_x);
    return (p - p_plus) - mTau*(p_minus - p_plus);
  }

  double eval(const Eigen::VectorXd& _x) override final
  {
    double cost = computeCost(_x);
    return cost*cost;
  }

  void evalGradient(const Eigen::VectorXd& _x,
                    Eigen::Map<Eigen::VectorXd> _grad) override final
  {
    _grad.setZero();
    const double cost = computeCost(_x);
    for(const Term& term : p_terms)
      _grad[term.index] = 2 * term.coeff * cost;
  }

protected:

  std::vector<Term> p_terms;

  double p_plus;
  double p_minus;
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
//    std::cout << "hip grad:   " << _grad[mIK->getNode()->getSkeleton()->getDof("LHP")->getIndexInSkeleton()] << "\n"
//              << "knee grad:  " << _grad[mIK->getNode()->getSkeleton()->getDof("LKP")->getIndexInSkeleton()] << "\n"
//              << "ankle grad: " << _grad[mIK->getNode()->getSkeleton()->getDof("LAP")->getIndexInSkeleton()] << "\n"
//              << std::endl;
//    std::cout << "grad: " << _grad.transpose() << std::endl;

//    for(size_t i=0; i < _grad.size(); ++i)
//    {
//      if(_grad[i] != 0.0)
//        std::cout << "(" << i << ")" << mIK->getNode()->getSkeleton()->getDof(i)->getName() << " " << _grad[i] << " | ";
//    }
//    std::cout << std::endl;
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
  else
  {
    const std::shared_ptr<Problem> problem = solver->getProblem();
    for(size_t i=0; i < problem->getNumEqConstraints(); ++i)
    {
      const FunctionPtr& f = problem->getEqConstraint(i);
    }
  }

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
    LSR = hubo->getDof("LSR")->getIndexInSkeleton();
    RSR = hubo->getDof("RSR")->getIndexInSkeleton();
  }

  void customPreRefresh() override
  {
    if(mTrajectory.size() == 0 && count == 0)
    {
      std::cerr << "No trajectory was generated!" << std::endl;
      ++count;
      return;
    }

//    Eigen::VectorXd positions = mTrajectory[count];
//    positions[LSR] += 90.0*M_PI/180.0;
//    positions[RSR] -= 90.0*M_PI/180.0;


    double dt = mWorld->getTimeStep();

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

//    std::cout << hubo->getZMP().transpose() << " \t:\t ";
    std::cout << accelerations.block<8,1>(0,0).transpose() << std::endl;

//    hubo->setPositions(mMapping, positions);
//    clone->setPositions(mMapping, mRaw[count]);

    ++count;
    if(count >= mTrajectory.size())
      count = 0;
  }

protected:

  SkeletonPtr hubo;
  std::vector<Eigen::VectorXd> mTrajectory;
  size_t count;
  size_t LSR;
  size_t RSR;
};

Eigen::VectorXd getAlphas(const YAML::Node& a, size_t index)
{
  const YAML::Node entry = a[index];
  Eigen::VectorXd alphas(entry.size());
  for(size_t i=0; i < entry.size(); ++i)
    alphas[i] = entry[i].as<double>();

  return alphas;
}

//#define ADD_LINEAR_OUTPUT( index, T ) \
//{\
//  alphas = getAlphas(a, index-1);\
//  terms = T;\
//  std::shared_ptr<LinearCombo> f = std::make_shared<LinearCombo>(terms, alphas);\
//  bezierFuncs.push_back(f);\
//  problem->addEqConstraint(f);\
//}

#define ADD_LINEAR_OUTPUT( index, T ) \
{\
  alphas = getAlphas(a, index-1);\
  terms = T;\
  system->setEqn(index-1, terms, alphas);\
}

#define ADD_EE_ORIENTATION_OUTPUT( index1, index2, index3, body )\
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
  return calfLength*(-hubo->getDof(st+"AP")->getPosition())
      + thighLength*(-hubo->getDof(st+"AP")->getPosition()
                     -hubo->getDof(st+"KP")->getPosition());
}

double computeP(const YAML::Node& params, double t, double pdot0, double vd)
{
  double p0 = params["p0"].as<double>();

  double eps = 1.0;

  return vd*t + p0 + ((1.0-exp(-eps*t))*pdot0 + (exp(-eps*t)-1.0)*vd)/eps;
}

std::vector<Eigen::VectorXd> setupAndSolveProblem(
    const SkeletonPtr& hubo, const std::string& file,
    const std::string& st/*ance*/, const std::string& sw/*ing*/, double* pdot0, double* vd)
{
  Eigen::VectorXd lastPositions = hubo->getPositions();
  hubo->resetPositions();

  // Start by testing left stance leg
  YAML::Node config = YAML::LoadFile(file);

  YAML::Node domain = config["domain"];
  YAML::Node params = st == "L"? domain[1] : domain[0];
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

  std::shared_ptr<GradientDescentSolver> solver = std::make_shared<GradientDescentSolver>();
  std::shared_ptr<Problem> problem = std::make_shared<Problem>();
  problem->setDimension(hubo->getNumDofs());
  solver->setProblem(problem);
//  solver->setStepSize(0.1);
//  solver->setStepSize(1.0);
  solver->setStepSize(1.0);
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

  std::shared_ptr<LinearComboSystem<23,33>> system =
      std::make_shared<LinearComboSystem<23,33>>();
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
  ZERO_OUT_ALPHAS(9);
  ADD_LINEAR_OUTPUT(9,  make_terms(hubo, 1.0, st+"EP"));
//  ADD_LINEAR_OUTPUT(24,  make_terms(hubo, 1.0, st+"EP"));
  ADD_LINEAR_OUTPUT(10, make_terms(hubo, 1.0, st+"WY"));
  ADD_LINEAR_OUTPUT(11, make_terms(hubo, 1.0, st+"WP"));
  ADD_LINEAR_OUTPUT(12, make_terms(hubo, 1.0, st+"WR"));
  ADD_LINEAR_OUTPUT(13, make_terms(hubo, 1.0, "TSY"));
  ADD_LINEAR_OUTPUT(14, make_terms(hubo, 1.0, sw+"SP"));
  ADD_LINEAR_OUTPUT(15, make_terms(hubo, 1.0, sw+"SR"));
  ADD_LINEAR_OUTPUT(16, make_terms(hubo, 1.0, sw+"SY"));
  ZERO_OUT_ALPHAS(17);
  ADD_LINEAR_OUTPUT(17, make_terms(hubo, 1.0, sw+"EP"));
//  ADD_LINEAR_OUTPUT(24, make_terms(hubo, 1.0, sw+"EP"));
  ADD_LINEAR_OUTPUT(18, make_terms(hubo, 1.0, sw+"WY"));
  ADD_LINEAR_OUTPUT(19, make_terms(hubo, 1.0, sw+"WP"));
  ADD_LINEAR_OUTPUT(20, make_terms(hubo, 1.0, sw+"WR"));
  ADD_LINEAR_OUTPUT(21, make_terms(hubo, 1.0, sw+"KP"));

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
  ADD_LINEAR_OUTPUT(22, make_terms(hubo, -1.0, st+"AP", -1.0, st+"KP", -1.0, st+"HP",
                                   1.0, sw+"HP", legRatio, sw+"KP"));

  ADD_LINEAR_OUTPUT(23, make_terms(hubo, 1.0, st+"HR", -1.0, sw+"HR"));

  ADD_EE_ORIENTATION_OUTPUT( 24, 25, 26, hubo->getBodyNode("Body_"+sw+"AR") );

  std::cout << "thigh: " << thighLength << "\n" << "calf: " << calfLength << std::endl;

  {
    // Constraint #25
    terms = make_terms(hubo, -calfLength, st+"AP", -thighLength, st+"AP", -thighLength, st+"KP");
    std::shared_ptr<SatisfyTau> f = std::make_shared<SatisfyTau>(terms, p_plus, p_minus);
    problem->addEqConstraint(f);
    bezierFuncs.push_back(f);
  }

//  // Constraint #26
//  JacobianNode* support_foot = hubo->getBodyNode("Body_"+st+"AR");
//  InverseKinematicsPtr support = InverseKinematics::create(support_foot);
//  std::vector<size_t> allDofs;
//  for(size_t i=0; i < hubo->getNumDofs(); ++i)
//    allDofs.push_back(i);
//  support->setDofs(allDofs);
//  Eigen::VectorXd weights = Eigen::VectorXd::Zero(allDofs.size());
////  weights.head<6>() = Eigen::Vector6d::Constant(5.0);
//  weights.head<6>() = Eigen::Vector6d::Constant(10.0);
//  support->getGradientMethod().setComponentWeights(weights);
//  support->getGradientMethod().setComponentWiseClamp(std::numeric_limits<double>::infinity());
//  problem->addEqConstraint(support->getProblem()->getEqConstraint(0));

  solver->setNumMaxIterations(1000);
  solver->setTolerance(1e-6);

//  hubo->setPositions(q_plus);
  hubo->setPositions(lastPositions);
//  support->getTarget()->setTransform(support_foot->getTransform());
  std::vector<Eigen::VectorXd> trajectory;

//  hubo->setPositions(solve(solver, bezierFuncs, 0.0));
//  double p_plus = computeP(hubo, st, thighLength, calfLength);
//  hubo->setPositions(solve(solver, bezierFuncs, 1.0));
//  double p_minus = computeP(hubo, st, thighLength, calfLength);

  hubo->setPositions(lastPositions);
//  double p = computeP(hubo, st, thighLength, calfLength);

//  double tau = 0.0;
//  if(tau < 0)
//  {
//    std::cout << "clamping tau from " << tau << std::endl;
//    tau = 0;
//  }

  double p0 = params["p0"].as<double>();
  std::cout << "\n\n" << st << " Foot Trajectory (p0 " << p0 << ") :\n";
  double time = 0.0;
  double tau = 0.0;
  double lastTau = 0.0;
  do
  {
    double p = computeP(params, time,
                        pdot0==nullptr? params["pdot0"].as<double>() : *pdot0,
                        vd==nullptr? params["v"].as<double>() : *vd);
    tau = (p - p_plus)/(p_minus - p_plus);
    if(tau > 1.0)
    {
      std::cout << "stopping" << std::endl;
      break;
    }

//    std::cout << "p: " << p << " \t p_plus: " << p_plus << " \t p_minus: " << p_minus << std::endl;
    trajectory.push_back(solve(solver, bezierFuncs, tau));

    double actualP = computeP(hubo, st, thighLength, calfLength);
    double actualTau = (actualP - p_plus)/(p_minus - p_plus);

    std::cout << "tau: " << tau << " \t actualTau: " << actualTau
              << " \t p: " << p << " \t actualP: " << actualP << std::endl;

    if(std::abs(tau-lastTau) < 1e-6)
      break;

    lastTau = tau;

    time += hubo->getTimeStep();

    std::cout << "solved in " << solver->getLastNumIterations() << " steps" << std::endl;
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
  world->setTimeStep(1.0/200.0);

  std::cout << "dofs: " << hubo->getNumDofs() << std::endl;

//  std::string file = "/home/grey/projects/protoHuboGUI/params_2015-08-27T07-01-0400.yaml";
  std::string yaml = "/home/grey/projects/protoHuboGUI/params_2015-08-29T16-11-0400.yaml";

  bool loadfile = false;
//  loadfile = true;

  std::string filename = "/home/grey/projects/protoHuboGUI/trajectory.dat";

  std::vector<Eigen::VectorXd> raw_trajectory;
  if(loadfile)
  {
    std::ifstream file;
    file.open(filename);
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
      std::cerr << "Could not open file: " << filename << std::endl;
    }
  }
  else
  {
    double pdot0 = 0.0;
    double vd = 0.0;
//    std::vector<Eigen::VectorXd> leftRamp =
//        setupAndSolveProblem(hubo, yaml, "L", "R", &pdot0, nullptr);
    std::vector<Eigen::VectorXd> rightStance =
        setupAndSolveProblem(hubo, yaml, "R", "L", nullptr, nullptr);
    std::vector<Eigen::VectorXd> leftStance =
        setupAndSolveProblem(hubo, yaml, "L", "R", nullptr, nullptr);
//    std::vector<Eigen::VectorXd> rightRamp =
//        setupAndSolveProblem(hubo, yaml, "R", "L", nullptr, &vd);


//    for(const Eigen::VectorXd& pos : leftRamp)
//      raw_trajectory.push_back(pos);
    for(const Eigen::VectorXd& pos : rightStance)
      raw_trajectory.push_back(pos);
    for(const Eigen::VectorXd& pos : leftStance)
      raw_trajectory.push_back(pos);
//    for(const Eigen::VectorXd& pos : rightRamp)
//      raw_trajectory.push_back(pos);

    dumpTrajectory(raw_trajectory, filename);
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

    mOperator.addWaypoint(hubo->getPositions(mOperatorIndices));
    mOperator.setInterpolationMode(HUBO_PATH_SPLINE);

    mOperator.sendNewTrajectory();

    mOperator.clearWaypoints();

    sleep(10);

    mOperator.update();
    hubo_player_state_t player = mOperator.getPlayerState();
    while(player.trajectory_size == 0 || player.current_index < player.trajectory_size - 1)
    {
      mOperator.update();
      player = mOperator.getPlayerState();
    }


    for(size_t i=0; i < raw_trajectory.size(); ++i)
    {
      hubo->setPositions(raw_trajectory[i]);
      mOperator.addWaypoint(hubo->getPositions(mOperatorIndices));
    }

    mOperator.setInterpolationMode(HUBO_PATH_RAW);

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
