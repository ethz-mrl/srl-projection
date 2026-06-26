/*
 * SPDX-FileCopyrightText: 2023 Smart Robotics Lab, Imperial College London
 * SPDX-FileCopyrightText: 2023-2026 Sotiris Papatheodorou
 * SPDX-License-Identifier: BSD-3-Clause
 */

#ifndef SRL_PROJECTION_RAYCASTINGCAMERA_HPP
#define SRL_PROJECTION_RAYCASTINGCAMERA_HPP

#include <Eigen/Geometry>
#include <cmath>
#include <opencv2/core/core.hpp>

namespace srl {
namespace projection {

/** A virtual camera used for 360° raycasting using some projection model \p ProjectionT derived
 * from srl::projection::ProjectionBase. The srl::projection::RaycastingCamera instance observes the
 * same space as that obtained by rotating \p ProjectionT by 360° around the z-axis of the body
 * frame B. Unlike other cameras it operates on the raycasting-body frame Br which has the same
 * orientation as the world frame W and the same origin as the body frame B.
 *
 * \bug Only the plain, 2-argument project()/backProject() methods are implemented. The others will
 * give garbage results.
 */
template<typename ProjectionT>
struct RaycastingCamera : public ProjectionT {
  /** Construct an instance producing 360° images with dimensions \p raycastingResolution from a \p
   * sensor and a transformation \p T_BS from the z-forward, x-right sensor frame S of \p sensor to
   * the x-forward, z-up body frame B.
   *
   * \warning Assumes \p T_BS doesn't have a roll component or a translation in the y-axis.
   */
  RaycastingCamera(const Eigen::Vector2i& raycastingResolution,
                   const ProjectionT& sensor,
                   const float_t nearPlane,
                   const float_t farPlane,
                   const Isometry3f& T_BS);

  ProjectionStatus project(const Vector3f& pointB, Vector2f* imagePoint) const override;

  bool backProject(const Vector2f& imagePoint, Vector3f* directionB) const override;

  int raycastingWidth() const;

  int raycastingHeight() const;

  const cv::Mat& raysBr() const;

  const float_t nearPlane;

  const float_t farPlane;

  const Isometry3f T_BS;

  const Isometry3f T_SB;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  protected:
  const int raycastingWidth_;
  const int raycastingHeight_;
  const float_t pixelHFov_;
  const float_t invPixelHFov_;
  cv::Mat raysBr_;
};

} // namespace projection
} // namespace srl

#include "implementation/RaycastingCamera.hpp"

#endif // SRL_PROJECTION_RAYCASTINGCAMERA_HPP
