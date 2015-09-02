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
#include <iostream>
#include <fstream>

#include <dart/dart.h>
#include <osgDart/osgDart.h>

#include <osg/Timer>

using namespace dart::dynamics;
using namespace dart::simulation;
using namespace osgDart;

class SimulationWorld : public osgDart::WorldNode
{
public:

  SimulationWorld(WorldPtr world, const std::vector<Eigen::VectorXd>& trajectory)
    : osgDart::WorldNode(world), mTrajectory(trajectory)
  {
    mHubo = world->getSkeleton("drchubo");

    double height = mHubo->getPosition(5);
    mHubo->setPositions(mTrajectory[0]);
    mHubo->setPosition(5, height);
//    index = 0;
  }

  void customPreStep() override
  {
    if(mWorld->getTime() < t0)
      return;

    double time = mWorld->getTime() - t0;
    size_t index = floor(200*(time - t0));

    if(index >= mTrajectory.size())
      return;


    std::cout << "Index: " << index << std::endl;

    const Eigen::VectorXd& qd = mTrajectory[index];
    for(size_t j=6; j < mHubo->getNumDofs(); ++j)
    {
      DegreeOfFreedom* dof = mHubo->getDof(j);
      double velocity = qd[j] - dof->getPosition();
      velocity /= mHubo->getTimeStep();

      dof->setCommand(velocity);
    }

  }

protected:

  const double frequency = 200;

  SkeletonPtr mHubo;

  std::vector<Eigen::VectorXd> mTrajectory;

  const double t0 = 1.0;

//  size_t index;
};

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
//  world->setTimeStep(1.0/200.0);

  SkeletonPtr hubo = createHubo();
  world->addSkeleton(hubo);
  world->addSkeleton(createGround());

  std::string name = "/home/ayonga/protoHuboGUI/trajectory.dat";

  std::ifstream file;
  file.open(name);
  std::vector<Eigen::VectorXd> trajectory;
  if(file.is_open())
  {
    while(!file.eof())
    {
      trajectory.push_back(Eigen::VectorXd(hubo->getNumDofs()));
      Eigen::VectorXd& q = trajectory.back();
      for(size_t i=0; i < hubo->getNumDofs(); ++i)
        file >> q[i];
    }
  }
  else
  {
    std::cerr << "Could not open file: " << name << std::endl;
  }

  osg::ref_ptr<SimulationWorld> node = new SimulationWorld(world, trajectory);
  node->setNumStepsPerCycle(30);

  osgDart::Viewer viewer;
  viewer.addWorldNode(node);
//  viewer.simulate(true);

  viewer.setUpViewInWindow(0, 0, 1280, 960);

  viewer.addAttachment(new osgDart::SupportPolygonVisual(hubo));

  // Set up the default viewing position
  viewer.getCameraManipulator()->setHomePosition(osg::Vec3( 5.34,  3.00, 1.91),
                                                 osg::Vec3( 0.00,  0.00, 0.50),
                                                 osg::Vec3(-0.20, -0.08, 0.98));

  viewer.setCameraManipulator(viewer.getCameraManipulator());
//  viewer.record("/home/ayonga/dump");

  viewer.run();
}
