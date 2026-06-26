/*
 * SPDX-FileCopyrightText: 2023 Smart Robotics Lab, Imperial College London
 * SPDX-FileCopyrightText: 2023-2026 Sotiris Papatheodorou
 * SPDX-License-Identifier: BSD-3-Clause
 */

#ifndef SRL_PROJECTION_RAYCASTINGCAMERA_IMPL_HPP
#define SRL_PROJECTION_RAYCASTINGCAMERA_IMPL_HPP

namespace srl {
namespace projection {

template<typename ProjectionT>
RaycastingCamera<ProjectionT>::RaycastingCamera(const Eigen::Vector2i& raycastingResolution,
                                                const ProjectionT& sensor,
                                                const float_t nearPlane,
                                                const float_t farPlane,
                                                const Isometry3f& T_BS) :
    ProjectionT(sensor),
    nearPlane(nearPlane),
    farPlane(farPlane),
    T_BS(T_BS),
    T_SB(T_BS.inverse()),
    raycastingWidth_(raycastingResolution.x()),
    raycastingHeight_(raycastingResolution.y()),
    pixelHFov_(6.28318530717958647692528677 / raycastingWidth_), // tau / raycastingWidth_
    invPixelHFov_(1 / pixelHFov_),
    raysBr_(cv::Size(raycastingWidth_, raycastingHeight_), CV_32FC3)
{
  // The assumption that a 360° image can be produced by rotation of the camera around the z-axis
  // of the body frame B doesn't hold when T_BS contains a translation in the y-axis.
  assert(std::abs(T_BS.translation().y()) < 1e-5f);

  for (int y = 0; y < raysBr_.rows; y++) {
    for (int x = 0; x < raysBr_.cols; x++) {
      // Everything needed to call backProject has been initialized by now.
      Vector3f ray;
      backProject(Vector2f(x, y), &ray);
      raysBr_.at<cv::Point3f>(y, x) = cv::Point3f(ray.x(), ray.y(), ray.z());
    }
  }
}



template<typename ProjectionT>
ProjectionStatus RaycastingCamera<ProjectionT>::project(const Vector3f& pointB,
                                                        Vector2f* imagePoint) const
{
  // The body ray frame Br has the same origin and z-axis as the body frame B but its x-axis
  // points in the direction of pointB.
  const float_t theta = std::atan2(pointB.y(), pointB.x());
  const Isometry3f T_BrB(AngleAxisf(-theta, Vector3f::UnitZ()));
  // Transform pointB to the sensor ray frame Sr, i.e. the frame of an
  // srl::projection::PinholeCamera mounted on Br, and project it.
  const Isometry3f& T_SrBr = T_SB;
  const Vector3f pointSr = T_SrBr * T_BrB * pointB;
  Vector2f originalImagePoint;
  const auto status = ProjectionT::project(pointSr, &originalImagePoint);
  // Apply the inverse of T_BrB on the image coordinates.
  imagePoint->x() = raycastingWidth() / 2.0f - theta * invPixelHFov_;
  // Ensure coordinates ∈[dim-0.5, dim) wrap-around to [-0.5, 0).
  if (imagePoint->x() >= raycastingWidth() - 0.5f) {
    imagePoint->x() -= raycastingWidth();
  }
  // Convert the y coordinate of originalImagePoint from the original image to the raycasting
  // image by going through UV coordinates.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
  const float_t originalImagePointV = (originalImagePoint.y() + 0.5f) / this->imageHeight();
#pragma GCC diagnostic pop
  imagePoint->y() = originalImagePointV * raycastingHeight() - 0.5f;
  return status;
}



template<typename ProjectionT>
bool RaycastingCamera<ProjectionT>::backProject(const Vector2f& imagePoint,
                                                Vector3f* directionB) const
{
  // The projection model considers pixel coordinates ∈[-0.5, dim-0.5) to be inside the image thus
  // the midpoint is at (dim-1)/2.
  const float_t originalImageMidpointX = (this->imageWidth() - 1) / 2.0f;
  // Convert the y coordinate of imagePoint from the raycasting image to the original image by
  // going through UV coordinates.
  const float_t imagePointV = (imagePoint.y() + 0.5f) / raycastingHeight();
  const float_t originalImagePointY = imagePointV * this->imageHeight() - 0.5f;
  // Apply T_BrB on the image coordinates, i.e. center imagePoint horizontally on an original
  // camera image.
  const Vector2f originalImagePoint(originalImageMidpointX, originalImagePointY);
  Vector3f directionSr;
  const bool success = ProjectionT::backProject(originalImagePoint, &directionSr);
  // Transform the backprojected point to the body frame B.
  const float_t theta = pixelHFov_ * (raycastingWidth() / 2.0f - imagePoint.x());
  const Isometry3f T_BBr(AngleAxisf(theta, Vector3f::UnitZ()));
  const Isometry3f& T_BrSr = T_BS;
  // Normalize the ray instead of setting the z component to 1 since this is no longer a pinhole
  // camera.
  *directionB = T_BBr * T_BrSr * directionSr.normalized();
  return success;
}



template<typename ProjectionT>
int RaycastingCamera<ProjectionT>::RaycastingCamera::raycastingWidth() const
{
  return raycastingWidth_;
}



template<typename ProjectionT>
int RaycastingCamera<ProjectionT>::raycastingHeight() const
{
  return raycastingHeight_;
}



template<typename ProjectionT>
const cv::Mat& RaycastingCamera<ProjectionT>::raysBr() const
{
  return raysBr_;
}

} // namespace projection
} // namespace srl

#endif // SRL_PROJECTION_RAYCASTINGCAMERA_IMPL_HPP
