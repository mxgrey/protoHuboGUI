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

#include <memory>

#include <dart/dart.h>
#include <osgDart/osgDart.h>

#include <HuboCan/HuboDescription.hpp>
#include <HuboCan/VirtualPump.hpp>
#include <HuboCmd/Aggregator.hpp>
#include <HuboState/State.hpp>

#include <osg/Timer>

using namespace dart::dynamics;
using namespace dart::simulation;
using namespace osgDart;
using namespace HuboCan;

class SimulationWorld : public osgDart::WorldNode
{
public:

  SimulationWorld(WorldPtr world)
    : osgDart::WorldNode(world)
  {
    mHubo = world->getSkeleton("drchubo");

    mDesc.parseFile("/opt/hubo/devices/DrcHubo.dd");
    mDesc.broadcastInfo();

    for(size_t i=0; i < mDesc.getJointCount(); ++i)
    {
      const hubo_joint_info_t& info = mDesc.joints[i]->info;
      DegreeOfFreedom* dof = mHubo->getDof(info.name);
      if(dof)
        mIndexMapping.push_back(dof->getIndexInSkeleton());
      else
        mIndexMapping.push_back(dart::dynamics::INVALID_INDEX);
    }

    mState = std::unique_ptr<HuboState::State>(new HuboState::State(mDesc));
    mAgg = std::unique_ptr<HuboCmd::Aggregator>(new HuboCmd::Aggregator(mDesc));

    if(!mState->initialized())
    {
      std::cout << "State was not initialized correctly, so we are quitting.\n"
                << " -- Either your ach channels are not open"
                << " or your HuboDescription was not valid!\n" << std::endl;
      return;
    }

    mAgg->run();

    mNumTicks = 0;
    mTimer.setStartTick();
  }

  void customPreStep() override
  {
    const HuboCmd::JointCmdArray& commands = mAgg->update();
    assert(commands.size() == mDesc.getJointCount());
    for(size_t i=0; i < commands.size(); ++i)
    {
      size_t index = mIndexMapping[i];
      if(dart::dynamics::INVALID_INDEX == index)
        continue;

      const hubo_joint_cmd_t& cmd = commands[i];
      DegreeOfFreedom* dof = mHubo->getDof(index);
      double velocity = cmd.position - dof->getPosition();
      velocity /= mHubo->getTimeStep();

      dof->setCommand(velocity);

      hubo_joint_state_t& state = mState->joints[i];
      state.duty = dof->getCommand();
      state.current = velocity;
    }
  }

  void customPostStep() override
  {
    double time = mTimer.time_s();
    double nextTickTime = (double)(mNumTicks)/mDesc.params.frequency;

    if(time < nextTickTime)
      return;

    ++mNumTicks;

    const HuboCmd::JointCmdArray& commands = mAgg->last_commands();
    for(size_t i=0; i < mState->joints.size(); ++i)
    {
      size_t index = mIndexMapping[i];
      if(dart::dynamics::INVALID_INDEX == index)
        continue;

      const hubo_joint_cmd_t& cmd = commands[i];
      hubo_joint_state_t& state = mState->joints[i];
      const DegreeOfFreedom* dof = mHubo->getDof(index);
      state.reference = cmd.position;
      state.position = dof->getPosition();
    }

    mState->publish();
  }

protected:

  SkeletonPtr mHubo;

  HuboDescription mDesc;
  std::vector<size_t> mIndexMapping;
  std::unique_ptr<HuboState::State> mState;
  std::unique_ptr<HuboCmd::Aggregator> mAgg;

  osg::Timer mTimer;

  size_t mNumTicks;
};

SkeletonPtr createHubo()
{
  dart::utils::DartLoader loader;
  loader.addPackageDirectory("drchubo", DART_DATA_PATH"/urdf/drchubo");
  SkeletonPtr hubo =
      loader.parseSkeleton(DART_DATA_PATH"/urdf/drchubo/drchubo.urdf");

  hubo->setPosition(5, 0.97);

  for(size_t i=1; i < hubo->getNumJoints(); ++i)
  {
    hubo->getJoint(i)->setActuatorType(Joint::VELOCITY);
  }

  for(size_t i=0; i < hubo->getNumBodyNodes(); ++i)
  {
    BodyNode* bn = hubo->getBodyNode(i);
    for(size_t j=0; j < bn->getNumVisualizationShapes(); ++j)
    {
      const ShapePtr& shape = bn->getVisualizationShape(j);
      shape->setColor(Eigen::Vector3d(0.2, 0.2, 0.2));

      if(MeshShapePtr mesh = std::dynamic_pointer_cast<MeshShape>(shape))
        mesh->setColorMode(MeshShape::SHAPE_COLOR);
    }
  }

  hubo->setName("drchubo");
  return hubo;
}

SkeletonPtr createGround()
{
  // Create a Skeleton to represent the ground
  SkeletonPtr ground = Skeleton::create("ground");
  Eigen::Isometry3d tf(Eigen::Isometry3d::Identity());
  double thickness = 0.01;
  tf.translation() = Eigen::Vector3d(0,0,-thickness/2.0);
  WeldJoint::Properties joint;
  joint.mT_ParentBodyToJoint = tf;
  ground->createJointAndBodyNodePair<WeldJoint>(nullptr, joint);
  ShapePtr groundShape =
      std::make_shared<BoxShape>(Eigen::Vector3d(10,10,thickness));
  groundShape->setColor(dart::Color::Blue(0.2));

  ground->getBodyNode(0)->addVisualizationShape(groundShape);
  ground->getBodyNode(0)->addCollisionShape(groundShape);

  return ground;
}

int main()
{
  WorldPtr world(new dart::simulation::World);

  world->addSkeleton(createHubo());
  world->addSkeleton(createGround());

  osg::ref_ptr<SimulationWorld> node = new SimulationWorld(world);
  node->setNumStepsPerCycle(30);

  osgDart::Viewer viewer;
  viewer.addWorldNode(node);
  viewer.simulate(true);

  viewer.setUpViewInWindow(0, 0, 640, 480);

  // Set up the default viewing position
  viewer.getCameraManipulator()->setHomePosition(osg::Vec3( 5.34,  3.00, 1.91),
                                                 osg::Vec3( 0.00,  0.00, 0.50),
                                                 osg::Vec3(-0.20, -0.08, 0.98));

  viewer.setCameraManipulator(viewer.getCameraManipulator());

  viewer.run();
}
