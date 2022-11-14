/*
 * SPDX-FileCopyrightText: 2022-2025 Smart Robotics Lab / Imperial College London
 * SPDX-FileCopyrightText: 2022-2025 Smart Robotics Lab / Technical University of Munich
 * SPDX-FileCopyrightText: 2022-2025 Sotiris Papatheodorou
 * SPDX-License-Identifier: BSD-3-Clause
 */

#include <gtest/gtest.h>

#include "srl/projection/EquirectangularCamera.hpp"

// All tests follow the same structure, just using different (back-)projection methods:
// 1. Back-project a random image point.
// 2. Randomize the 3D point distance.
// 3. Project the 3D point.
// 4. Ensure the two image points are the same.
class EquirectangularCamera : public ::testing::Test {
  protected:
    EquirectangularCamera() :
      camera_(srl::projection::EquirectangularCamera::createTestInstance()),
      imagePoints_(2, 1000)
    {
    }

    void SetUp()
    {
      // The constructor can't have asserts.
      ASSERT_TRUE(camera_);
      ASSERT_FALSE(camera_->hasMask());
      for (int i = 0; i < imagePoints_.cols(); ++i) {
         imagePoints_.col(i) = camera_->createRandomImagePoint();
      }
    }

    srl::Vector3f randomizeDistance(const srl::Vector3f& ray)
    {
      srl::Vector3f point = ray;
      point.normalize();
      point *= 0.2f + 8 * (srl::Vector2f::Random()[0] + 1);
      return point;
    }

    srl::Vector4f randomizeDistance(const srl::Vector4f& ray)
    {
      return randomizeDistance(ray.template head<3>().eval()).homogeneous();
    }

    static constexpr double imageThres_ = 0.01;
    std::shared_ptr<srl::projection::ProjectionBase> camera_;
    srl::Matrix2Xf imagePoints_;

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

constexpr double EquirectangularCamera::imageThres_;



TEST_F(EquirectangularCamera, backProjectAndProject)
{
  for (int i = 0; i < imagePoints_.cols(); ++i) {
    const srl::Vector2f imagePoint = imagePoints_.col(i);

    srl::Vector3f ray;
    ASSERT_TRUE(camera_->backProject(imagePoint, &ray))
      << "imagePoint: " << imagePoint.transpose();

    const srl::Vector3f point = randomizeDistance(ray);
    srl::Vector2f imagePoint2;
    ASSERT_EQ(camera_->project(point, &imagePoint2), srl::projection::ProjectionStatus::Successful)
      << "point: " << point.transpose();

    ASSERT_LT((imagePoint - imagePoint2).norm(), imageThres_)
      << "imagePoint:  " << imagePoint.transpose()
      << "\nimagePoint2: " << imagePoint2.transpose();
  }
}



TEST_F(EquirectangularCamera, backProjectAndProjectHomogeneous)
{
  for (int i = 0; i < imagePoints_.cols(); ++i) {
    const srl::Vector2f imagePoint = imagePoints_.col(i);

    srl::Vector4f ray;
    ASSERT_TRUE(camera_->backProjectHomogeneous(imagePoint, &ray))
      << "imagePoint: " << imagePoint.transpose();

    const srl::Vector4f point = randomizeDistance(ray);
    srl::Vector2f imagePoint2;
    ASSERT_EQ(camera_->projectHomogeneous(point, &imagePoint2), srl::projection::ProjectionStatus::Successful)
        << "point: " << point.transpose();

    ASSERT_LT((imagePoint - imagePoint2).norm(), imageThres_)
      << "imagePoint:  " << imagePoint.transpose()
      << "\nimagePoint2: " << imagePoint2.transpose();
  }
}



TEST_F(EquirectangularCamera, backProjectAndProjectBatch)
{
  srl::Matrix3Xf rays(3, imagePoints_.cols());
  std::vector<bool> success;
  ASSERT_TRUE(camera_->backProjectBatch(imagePoints_, &rays, &success));
  for (size_t i = 0; i < success.size(); ++i) {
    ASSERT_TRUE(success[i])
      << "at success[" << i << "]";
  }

  srl::Matrix3Xf points(3, rays.cols());
  for (int i = 0; i < points.cols(); ++i) {
    points.col(i) = randomizeDistance(rays.col(i).eval());
  }
  srl::Matrix2Xf imagePoints2(2, points.cols());
  std::vector<srl::projection::ProjectionStatus> stati(points.cols());
  camera_->projectBatch(points, &imagePoints2, &stati);
  for (size_t i = 0; i < stati.size(); ++i) {
    ASSERT_EQ(stati[i], srl::projection::ProjectionStatus::Successful)
      << "at stati[" << i << "]";
  }

  for (int i = 0; i < imagePoints_.cols(); ++i) {
    ASSERT_LT((imagePoints_.col(i) - imagePoints2.col(i)).norm(), imageThres_)
      << "imagePoint:  " << imagePoints_.col(i).transpose()
      << "\nimagePoint2: " << imagePoints2.col(i).transpose();
  }
}



TEST_F(EquirectangularCamera, backProjectAndProjectHomogeneousBatch)
{
  srl::Matrix4Xf rays(4, imagePoints_.cols());
  std::vector<bool> success;
  ASSERT_TRUE(camera_->backProjectHomogeneousBatch(imagePoints_, &rays, &success));
  for (size_t i = 0; i < success.size(); ++i) {
    ASSERT_TRUE(success[i])
      << "at success[" << i << "]";
  }

  srl::Matrix4Xf points(4, rays.cols());
  for (int i = 0; i < points.cols(); ++i) {
    points.col(i) = randomizeDistance(rays.col(i).eval());
  }
  srl::Matrix2Xf imagePoints2(2, points.cols());
  std::vector<srl::projection::ProjectionStatus> stati(points.cols());
  camera_->projectHomogeneousBatch(points, &imagePoints2, &stati);
  for (size_t i = 0; i < stati.size(); ++i) {
    ASSERT_EQ(stati[i], srl::projection::ProjectionStatus::Successful)
      << "at stati[" << i << "]";
  }

  for (int i = 0; i < imagePoints_.cols(); ++i) {
    ASSERT_LT((imagePoints_.col(i) - imagePoints2.col(i)).norm(), imageThres_)
      << "imagePoint:  " << imagePoints_.col(i).transpose()
      << "\nimagePoint2: " << imagePoints2.col(i).transpose();
  }
}
