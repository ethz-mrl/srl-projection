/*
 * SPDX-FileCopyrightText: 2022-2025 Smart Robotics Lab / Imperial College London
 * SPDX-FileCopyrightText: 2022-2025 Smart Robotics Lab / Technical University of Munich
 * SPDX-FileCopyrightText: 2022-2025 Sotiris Papatheodorou
 * SPDX-License-Identifier: BSD-3-Clause
 */

namespace srl {
namespace projection {

SphericalCamera::SphericalCamera(int imageWidth, int imageHeight, float_t horizontalFov, float_t verticalFov)
    : ProjectionBase(imageWidth, imageHeight)
{
  intrinsics_[0] = horizontalFov;
  intrinsics_[1] = verticalFov;
}

void SphericalCamera::getIntrinsics(VectorXf& intrinsics) const {
  intrinsics = intrinsics_;
}

bool SphericalCamera::setIntrinsics(const VectorXf& intrinsics) {
  if (intrinsics.size() != intrinsics_.size()) {
    return false;
  }
  intrinsics_ = intrinsics;
  return true;
}

int SphericalCamera::numIntrinsicsParameters() const
{
  return intrinsics_.size();
}

float_t SphericalCamera::horizontalFov() const
{
  return intrinsics_[0];
}

float_t SphericalCamera::verticalFov() const
{
  return intrinsics_[1];
}

//////////////////////////////////////////
// Methods to project points

ProjectionStatus SphericalCamera::project(const Vector3f& point, Vector2f* imagePoint) const
{
  const float_t norm = point.norm();
  if (norm < 1e-5f || (std::fabs(point.x()) < 1e-5f && std::fabs(point.y()) < 1e-5f)) {
    return ProjectionStatus::Invalid;
  }
  const float_t azimuth = std::atan2(point.y(), point.x());
  if (azimuth > horizontalFov() / 2 || azimuth < -horizontalFov() / 2) {
    return ProjectionStatus::OutsideImage;
  }
  const float_t inclination = std::acos(point.z() / norm);
  const float_t min_inclination = (M_PI - verticalFov()) / 2;
  if (inclination < min_inclination || inclination > min_inclination + verticalFov()) {
    return ProjectionStatus::OutsideImage;
  }
  imagePoint->x() = imageWidth_ * (-azimuth / horizontalFov() + 0.5f);
  imagePoint->y() = (imageHeight_ - 1) * (inclination - min_inclination) / verticalFov();
  if (ProjectionBase::isMasked(*imagePoint)) {
    return ProjectionStatus::Masked;
  }
  return ProjectionStatus::Successful;
}

ProjectionStatus SphericalCamera::project(const Vector3f& point,
                                          Vector2f* imagePoint,
                                          Matrixf<2, 3>* /*pointJacobian*/,
                                          Matrix2Xf* /*intrinsicsJacobian*/) const
{
  // TODO: Write into pointJacobian, intrinsicsJacobian.
  return project(point, imagePoint);
}

ProjectionStatus SphericalCamera::projectWithExternalParameters(const Vector3f& point,
                                                                const VectorXf& parameters,
                                                                Vector2f* imagePoint,
                                                                Matrixf<2, 3>* pointJacobian,
                                                                Matrix2Xf* intrinsicsJacobian) const
{
  if (parameters.size() != numIntrinsicsParameters()) {
    return ProjectionStatus::Invalid;
  }
  return SphericalCamera(imageWidth(), imageHeight(), parameters[0], parameters[1]).project(point,
      imagePoint, pointJacobian, intrinsicsJacobian);
}

void SphericalCamera::projectBatch(const Matrix3Xf& points,
                                   Matrix2Xf* imagePoints,
                                   std::vector<ProjectionStatus>* stati) const
{
  for (int i = 0; i < points.cols(); ++i) {
    Vector2f imagePoint;
    const ProjectionStatus status = project(points.col(i), &imagePoint);
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
    imagePoints->col(i) = imagePoint;
#pragma GCC diagnostic pop
    if (stati) {
      stati->push_back(status);
    }
  }
}

ProjectionStatus SphericalCamera::projectHomogeneous(const Vector4f& point,
                                                     Vector2f* imagePoint) const
{
  if (point.w() < 0) {
    return project(-point.head<3>(), imagePoint);
  } else {
    return project(point.head<3>(), imagePoint);
  }
}

ProjectionStatus SphericalCamera::projectHomogeneous(const Vector4f& point,
                                                     Vector2f* imagePoint,
                                                     Matrixf<2, 4>* pointJacobian,
                                                     Matrix2Xf* intrinsicsJacobian) const
{
  Matrixf<2, 3> pointJacobian3;
  ProjectionStatus status;
  if (point.w() < 0) {
    status = project(-point.head<3>(), imagePoint, &pointJacobian3, intrinsicsJacobian);
  } else {
    status = project(point.head<3>(), imagePoint, &pointJacobian3, intrinsicsJacobian);
  }
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
  pointJacobian->template leftCols<3>() = pointJacobian3;
#pragma GCC diagnostic pop
  pointJacobian->template rightCols<1>() = Vector2f::Zero();
  return status;
}

ProjectionStatus SphericalCamera::projectHomogeneousWithExternalParameters(
    const Vector4f& point,
    const VectorXf& parameters,
    Vector2f* imagePoint,
    Matrixf<2, 4>* pointJacobian,
    Matrix2Xf* intrinsicsJacobian) const
{
  if (parameters.size() != numIntrinsicsParameters()) {
    return ProjectionStatus::Invalid;
  }
  return SphericalCamera(imageWidth(), imageHeight(), parameters[0],
      parameters[1]).projectHomogeneous(point, imagePoint, pointJacobian, intrinsicsJacobian);
}

void SphericalCamera::projectHomogeneousBatch(const Matrix4Xf& points,
                                              Matrix2Xf* imagePoints,
                                              std::vector<ProjectionStatus>* stati) const
{
  for (int i = 0; i < points.cols(); ++i) {
    Vector2f imagePoint;
    const ProjectionStatus status = projectHomogeneous(points.col(i), &imagePoint);
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
    imagePoints->col(i) = imagePoint;
#pragma GCC diagnostic pop
    if (stati) {
      stati->push_back(status);
    }
  }
}

//////////////////////////////////////////
// Methods to backproject points

bool SphericalCamera::backProject(const Vector2f& imagePoint, Vector3f* direction) const
{
  // Azimuth is horizontalFov at the left edge of the image and -horizontalFov at the right edge.
  // Normalize with imageWidth_ so that the azimuth is in the interval
  // (-horizontalFov, horizontalFov].
  const float_t azimuth = -horizontalFov() * (imagePoint.x() / imageWidth_ - 0.5f);
  // Inclination is min_inclination at the top of the image and min_inclination + verticalFov at
  // the bottom. Normalize with imageHeight_ - 1 so that the inclination is in the interval
  // [min_inclination, min_inclination + verticalFov].
  const float_t min_inclination = (M_PI - verticalFov()) / 2;
  const float_t inclination = verticalFov() * imagePoint.y() / (imageHeight_ - 1) + min_inclination;
  direction->x() = std::sin(inclination) * std::cos(azimuth);
  direction->y() = std::sin(inclination) * std::sin(azimuth);
  direction->z() = std::cos(inclination);
  return true;
}

inline bool SphericalCamera::backProject(const Vector2f& imagePoint,
                                         Vector3f* direction,
                                         Matrixf<3, 2>* /*pointJacobian*/) const
{
  // TODO: Write into pointJacobian
  return backProject(imagePoint, direction);
}

bool SphericalCamera::backProjectBatch(const Matrix2Xf& imagePoints,
                                       Matrix3Xf* directions,
                                       std::vector<bool>* success) const
{
  for (int i = 0; i < imagePoints.cols(); ++i) {
    Vector3f point;
    const bool suc = backProject(imagePoints.col(i), &point);
    directions->col(i) = point;
    if (success) {
      success->push_back(suc);
    }
  }
  return true;
}

bool SphericalCamera::backProjectHomogeneous(const Vector2f& imagePoint, Vector4f* direction) const
{
  Vector3f ray;
  const bool success = backProject(imagePoint, &ray);
  *direction = ray.homogeneous();
  return success;
}

bool SphericalCamera::backProjectHomogeneous(const Vector2f& imagePoint,
                                             Vector4f* direction,
                                             Matrixf<4, 2>* pointJacobian) const
{
  Vector3f ray;
  Matrixf<3, 2> pointJacobian3;
  const bool success = backProject(imagePoint, &ray, &pointJacobian3);
  *direction = ray.homogeneous();
  pointJacobian->template topRows<3>() = pointJacobian3;
  pointJacobian->template bottomRows<1>() = Vector2f::Zero();
  return success;
}

bool SphericalCamera::backProjectHomogeneousBatch(const Matrix2Xf& imagePoints,
                                                  Matrix4Xf* directions,
                                                  std::vector<bool>* success) const
{
  for (int i = 0; i < imagePoints.cols(); ++i) {
    Vector3f point;
    const bool suc = backProject(imagePoints.col(i), &point);
    directions->col(i) = point.homogeneous();
    if (success) {
      success->push_back(suc);
    }
  }
  return true;
}

std::shared_ptr<ProjectionBase> SphericalCamera::createTestInstance()
{
  return std::shared_ptr<ProjectionBase>(new SphericalCamera(1024, 512, 2 * M_PI, M_PI));
}

std::string SphericalCamera::type() const
{
  return "SphericalCamera";
}

}  // namespace projection
}  // namespace srl
