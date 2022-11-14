/*
 * SPDX-FileCopyrightText: 2022-2025 Smart Robotics Lab / Imperial College London
 * SPDX-FileCopyrightText: 2022-2025 Smart Robotics Lab / Technical University of Munich
 * SPDX-FileCopyrightText: 2022-2025 Sotiris Papatheodorou
 * SPDX-License-Identifier: BSD-3-Clause
 */

#ifndef INCLUDE_SRL_PROJECTION_EQUIRECTANGULARCAMERA_HPP_
#define INCLUDE_SRL_PROJECTION_EQUIRECTANGULARCAMERA_HPP_

#include <Eigen/Geometry>
#include <cassert>
#include <memory>
#include <vector>

#include "srl/projection/ProjectionBase.hpp"

namespace srl {
namespace projection {

/// \brief A standard equirectangular projection model using spherical coordinates.
///
/// The camera's frame is x-forward, y-left, z-up. Spherical coordinates are mapped to an image,
/// with azimuth (angle from the x axis) increasing towards the left and inclination (angle from the
/// z-axis) increasing downwards.
/// \code{.txt}
///              m ┌───────────────────┐
///                │                   │
/// Inclination    │       Image       │
///                │                   │
///       m + vFoV └───────────────────┘
///              hFoV/2      0°     -hFoV/2
///                       Azimuth
/// \endcode
/// where `m = (π - vFoV) / 2`.
///
/// - https://en.wikipedia.org/wiki/Equirectangular_projection
/// - https://en.wikipedia.org/wiki/Spherical_coordinate_system
class EquirectangularCamera : public ProjectionBase
{
 public:
  /// \param[in] imageWidth The image width in pixels.
  /// \param[in] imageHeight The image height in pixels.
  /// \param[in] horizontalFov The horizontal field of view in radians in the interval (0, 2π].
  /// \param[in] verticalFov The vertical field of view in radians in the interval (0, π].
  inline EquirectangularCamera(int imageWidth, int imageHeight, float_t horizontalFov, float_t verticalFov);

  /// \brief Write the intrinsics into a vector.
  /// \param[out] intrinsics The intrinsics as a vector where `intrinsics[0]` is the horizontal and
  /// `intrinsics[1]` is the vertical FoV in radians.
  inline void getIntrinsics(VectorXf& intrinsics) const override;

  /// \brief Overwrite all intrinsics from a vector.
  /// \param[in] intrinsics The intrinsics as a vector where `intrinsics[0]` is the horizontal and
  /// `intrinsics[1]` is the vertical FoV in radians.
  /// \return Whether the intrinsics were updated successfully.
  inline bool setIntrinsics(const VectorXf& intrinsics) override;

  /// \brief Return the number of intrinsic parameters.
  inline int numIntrinsicsParameters() const override;

  /// \brief Return the horizontal field of view in radians.
  inline float_t horizontalFov() const;

  /// \brief Return the vertical field of view in radians.
  inline float_t verticalFov() const;

  //////////////////////////////////////////////////////////////
  /// \name Methods to project points
  /// @{

  /// \brief Project a 3D point to 2D image coordinates.
  /// \param[in]  point      The 3D point in Euclidean coordinates.
  /// \param[out] imagePoint The 2D image coordinates.
  /// \return The projection status.
  inline ProjectionStatus project(const Vector3f& point, Vector2f* imagePoint) const override;

  /// \brief Project a 3D point to 2D image coordinates.
  /// \param[in]  point              The 3D point in Euclidean coordinates.
  /// \param[out] imagePoint         The 2D image coordinates.
  /// \param[out] pointJacobian      The Jacobian of the projection function w.r.t. the point.
  /// \param[out] intrinsicsJacobian The Jacobian of the projection function w.r.t. the intrinsics.
  /// \return The projection status.
  inline ProjectionStatus project(const Vector3f& point,
                                  Vector2f* imagePoint,
                                  Matrixf<2, 3>* pointJacobian,
                                  Matrix2Xf* intrinsicsJacobian = nullptr) const override;

  /// \brief Project a 3D point to 2D image coordinates.
  /// \param[in]  point              The 3D point in Euclidean coordinates.
  /// \param[in]  parameters         The intrinsic parameters.
  /// \param[out] imagePoint         The 2D image coordinates.
  /// \param[out] pointJacobian      The Jacobian of the projection function w.r.t. the point.
  /// \param[out] intrinsicsJacobian The Jacobian of the projection function w.r.t. the intrinsics.
  /// \return The projection status.
  inline ProjectionStatus projectWithExternalParameters(
      const Vector3f& point,
      const VectorXf& parameters,
      Vector2f* imagePoint,
      Matrixf<2, 3>* pointJacobian,
      Matrix2Xf* intrinsicsJacobian = nullptr) const override;

  /// \brief Project a batch of 3D points to 2D image coordinates.
  /// \param[in]  points      The 3D points in Euclidean coordinates, one per column.
  /// \param[out] imagePoints The 2D image coordinates, one point per column.
  /// \param[out] stati       The projection status for each point.
  inline void projectBatch(const Matrix3Xf& points,
                           Matrix2Xf* imagePoints,
                           std::vector<ProjectionStatus>* stati) const override;

  /// \brief Project a 3D point in homogeneous coordinates to 2D image coordinates.
  /// \param[in]  point      The 3D point in homogeneous Euclidean coordinates.
  /// \param[out] imagePoint The 2D image coordinates.
  /// \return The projection status.
  inline ProjectionStatus projectHomogeneous(const Vector4f& point, Vector2f* imagePoint) const override;

  /// \brief Project a 3D point in homogeneous coordinates to 2D image coordinates.
  /// \param[in]  point              The 3D point in homogeneous Euclidean coordinates.
  /// \param[out] imagePoint         The 2D image coordinates.
  /// \param[out] pointJacobian      The Jacobian of the projection function w.r.t. the point.
  /// \param[out] intrinsicsJacobian The Jacobian of the projection function w.r.t. the intrinsics.
  /// \return The projection status.
  inline ProjectionStatus projectHomogeneous(const Vector4f& point,
                                             Vector2f* imagePoint,
                                             Matrixf<2, 4>* pointJacobian,
                                             Matrix2Xf* intrinsicsJacobian = nullptr) const override;

  /// \brief Project a 3D point in homogeneous coordinates to 2D image coordinates.
  /// \param[in]  point              The 3D point in homogeneous Euclidean coordinates.
  /// \param[in]  parameters         The intrinsic parameters.
  /// \param[out] imagePoint         The 2D image coordinates.
  /// \param[out] pointJacobian      The Jacobian of the projection function w.r.t. the point.
  /// \param[out] intrinsicsJacobian The Jacobian of the projection function w.r.t. the intrinsics.
  /// \return The projection status.
  inline ProjectionStatus projectHomogeneousWithExternalParameters(
      const Vector4f& point,
      const VectorXf& parameters,
      Vector2f* imagePoint,
      Matrixf<2, 4>* pointJacobian = nullptr,
      Matrix2Xf* intrinsicsJacobian = nullptr) const override;

  /// \brief Project a batch of 3D points in homogeneous coordinates to 2D image coordinates.
  /// \param[in]  points      The 3D points in homogeneous Euclidean coordinates, one per column.
  /// \param[out] imagePoints The 2D image coordinates, one point per column.
  /// \param[out] stati       The projection status for each point.
  inline void projectHomogeneousBatch(const Matrix4Xf& points,
                                      Matrix2Xf* imagePoints,
                                      std::vector<ProjectionStatus>* stati) const override;
  /// @}

  //////////////////////////////////////////////////////////////
  /// \name Methods to backproject points
  /// @{

  /// \brief Back-project 2D image coordinates into a 3D direction vector.
  /// \param[in]  imagePoint The 2D image coordinates.
  /// \param[out] direction  The 3D direction vector in Euclidean coordinates.
  /// \return                Whether the back-projection succeeded.
  inline bool backProject(const Vector2f& imagePoint, Vector3f* direction) const override;

  /// \brief Back-project 2D image coordinates into a 3D direction vector.
  /// \param[in]  imagePoint    The 2D image coordinates.
  /// \param[out] direction     The 3D direction vector in Euclidean coordinates.
  /// \param[out] pointJacobian The Jacobian of the back-projection function  w.r.t. the point.
  /// \return                   Whether the back-projection succeeded.
  inline bool backProject(const Vector2f& imagePoint,
                          Vector3f* direction,
                          Matrixf<3, 2>* pointJacobian) const override;

  /// \brief Back-project a batch of 2D image coordinates into 3D direction vectors.
  /// \param[in]  imagePoints The 2D image coordinates, one per column.
  /// \param[out] directions  The 3D direction vectors in Euclidean coordinates, one per column.
  /// \param[out] success     Whether each back-projection succeeded.
  /// \return                 Always true.
  inline bool backProjectBatch(const Matrix2Xf& imagePoints,
                               Matrix3Xf* directions,
                               std::vector<bool>* success) const override;

  /// \brief Back-project 2D image coordinates into a 3D direction vector in homogeneous coordinates.
  /// \param[in]  imagePoint The 2D image coordinates.
  /// \param[out] direction  The 3D direction vector in homogeneous Euclidean coordinates.
  /// \return                Whether the back-projection succeeded.
  inline bool backProjectHomogeneous(const Vector2f& imagePoint, Vector4f* direction) const override;

  /// \brief Back-project 2D image coordinates into a 3D direction vector in homogeneous coordinates.
  /// \param[in]  imagePoint    The 2D image coordinates.
  /// \param[out] direction     The 3D direction vector in homogeneous Euclidean coordinates.
  /// \param[out] pointJacobian The Jacobian of the back-projection function  w.r.t. the point.
  /// \return                   Whether the back-projection succeeded.
  inline bool backProjectHomogeneous(const Vector2f& imagePoint,
                                     Vector4f* direction,
                                     Matrixf<4, 2>* pointJacobian) const override;

  /// \brief Back-project a batch of 2D image coordinates into 3D direction vectors in homogeneous coordinates.
  /// \param[in]  imagePoints The 2D image coordinates, one per column.
  /// \param[out] directions  The 3D direction vectors in homogeneous Euclidean coordinates, one per column.
  /// \param[out] success     Whether each back-projection succeeded.
  /// \return                 Always true.
  inline bool backProjectHomogeneousBatch(const Matrix2Xf& imagePoints,
                                          Matrix4Xf* directions,
                                          std::vector<bool>* success) const override;
  /// @}

  /// \brief Return a test instance.
  inline static std::shared_ptr<ProjectionBase> createTestInstance();

  /// \brief Return the string `"EquirectangularCamera"`.
  inline std::string type() const;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

 protected:

  Vector2f intrinsics_;
};

}  // namespace projection
}  // namespace srl

#include "implementation/EquirectangularCamera.hpp"

#endif /* INCLUDE_SRL_PROJECTION_EQUIRECTANGULARCAMERA_HPP_ */
