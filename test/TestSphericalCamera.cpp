/*
 * SPDX-FileCopyrightText: 2022-2025 Smart Robotics Lab / Imperial College London
 * SPDX-FileCopyrightText: 2022-2025 Smart Robotics Lab / Technical University of Munich
 * SPDX-FileCopyrightText: 2022-2025 Sotiris Papatheodorou
 * SPDX-License-Identifier: BSD-3-Clause
 */

#include <gtest/gtest.h>

#include "srl/projection/SphericalCamera.hpp"

TEST(SphericalCamera, functions)
{
  const auto camera = srl::projection::SphericalCamera::createTestInstance();

  ASSERT_TRUE(camera);
  ASSERT_FALSE(camera->hasMask());

  for (size_t i = 0; i < 100; ++i) {
    // Back-project a random point.
    const srl::Vector2f imagePoint = camera->createRandomImagePoint();
    srl::Vector3f ray;
    ASSERT_TRUE(camera->backProject(imagePoint, &ray)) << "imagePoint: " << imagePoint.transpose();

    // Randomise the point distance and project it.
    ray.normalize();
    ray *= (0.2f + 8.0f * (srl::Vector2f::Random()[0] + 1.0f));
    srl::Vector2f imagePoint2;
    ASSERT_EQ(camera->project(ray, &imagePoint2), srl::projection::ProjectionStatus::Successful)
        << "ray: " << ray.transpose();

    // Ensure the two image points are the same.
    ASSERT_LT((imagePoint - imagePoint2).norm(), 0.01) << "imagePoint: " << imagePoint.transpose()
      << "\nimagePoint2: " << imagePoint2.transpose();
  }
}
