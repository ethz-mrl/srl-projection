/*
 * SPDX-FileCopyrightText: 2023 Smart Robotics Lab, Imperial College London, Technical University of Munich
 * SPDX-FileCopyrightText: 2023-2026 Sotiris Papatheodorou
 * SPDX-License-Identifier: BSD-3-Clause
 */

#include <fstream>
#include <gtest/gtest.h>
#include <srl/projection/NoDistortion.hpp>
#include <srl/projection/PinholeCamera.hpp>
#include <srl/projection/RaycastingCamera.hpp>
#include <sstream>

TEST(RaycastingCamera, backProjectAndProject)
{
  typedef srl::projection::PinholeCamera<srl::projection::NoDistortion> Sensor;

  static constexpr float_t pi = 3.14159265358979323846264338;
  static constexpr int width = 640;
  static constexpr int height = 480;
  static constexpr float_t f = 554.25; // ~30° horizontal FoV
  static constexpr float_t cx = width / 2 - 0.5;
  static constexpr float_t cy = height / 2 - 0.5;
  static constexpr float_t near_plane = 0.2;
  static constexpr float_t far_plane = 1.0;
  const srl::Isometry3f T_BS = srl::Translation3f(0.1, 0, 0)
    * srl::AngleAxisf(pi / 2, srl::Vector3f::UnitY())
    * srl::AngleAxisf(-pi / 2, srl::Vector3f::UnitZ())
    * srl::AngleAxisf(-pi / 18, srl::Vector3f::UnitX());

  // We need a camera without distortion so we can't use
  // srl::projection::PinholeCamera::createTestObject().
  const srl::projection::NoDistortion distortion;
  const Sensor sensor (width, height, f, f, cx, cy, distortion);
  const srl::projection::RaycastingCamera<Sensor> raycasting_sensor(
    Eigen::Vector2i(36, 10), sensor, near_plane, far_plane, T_BS);

  for (int y = 0; y < raycasting_sensor.raycastingHeight(); y++) {
    for (int x = 0; x < raycasting_sensor.raycastingWidth(); x++) {
      const srl::Vector2f p(x, y);
      std::stringstream error_msg;
      error_msg << "for pixel [" << p.transpose() << "]";
      // Back-project
      const cv::Point3f ray_Br_des = raycasting_sensor.raysBr().at<cv::Point3f>(y, x);
      srl::Vector3f ray_Br = srl::Vector3f::Zero();
      ASSERT_TRUE(raycasting_sensor.backProject(p, &ray_Br)) << error_msg.str();
      EXPECT_FLOAT_EQ(ray_Br.x(), ray_Br_des.x) << error_msg.str();
      EXPECT_FLOAT_EQ(ray_Br.y(), ray_Br_des.y) << error_msg.str();
      EXPECT_FLOAT_EQ(ray_Br.z(), ray_Br_des.z) << error_msg.str();
      error_msg << " back-projected to ray_Br [" << ray_Br.transpose() << "]";
      // Project
      srl::Vector2f p2 = srl::Vector2f::Zero();
      const auto status = raycasting_sensor.project(ray_Br, &p2);
      ASSERT_EQ(status, srl::projection::ProjectionStatus::Successful) << error_msg.str();
      // The same image coordinates must be obtained. Due to numerical issues with floating
      // point numbers the resulting values can have a rather large error.
      EXPECT_NEAR(p.x(), p2.x(), 1e-2) << error_msg.str();
      EXPECT_NEAR(p.y(), p2.y(), 1e-2) << error_msg.str();
    }
  }
}
