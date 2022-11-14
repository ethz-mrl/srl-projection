/*
 * SPDX-FileCopyrightText: 2022-2025 Smart Robotics Lab / Imperial College London
 * SPDX-FileCopyrightText: 2022-2025 Smart Robotics Lab / Technical University of Munich
 * SPDX-FileCopyrightText: 2022-2025 Sotiris Papatheodorou
 * SPDX-License-Identifier: BSD-3-Clause
 */

namespace srl {
namespace projection {

EquirectangularCamera::EquirectangularCamera(int imageWidth, int imageHeight, float_t horizontalFov, float_t verticalFov)
    : ProjectionBase(imageWidth, imageHeight)
{
  assert(imageWidth > 0);
  assert(imageHeight > 0);
  assert(horizontalFov > 0);
  assert(horizontalFov <= 2 * M_PI);
  assert(verticalFov > 0);
  assert(verticalFov <= M_PI);
  intrinsics_[0] = horizontalFov;
  intrinsics_[1] = verticalFov;
}

void EquirectangularCamera::getIntrinsics(VectorXf& intrinsics) const {
  intrinsics = intrinsics_;
}

bool EquirectangularCamera::setIntrinsics(const VectorXf& intrinsics) {
  if (intrinsics.size() != intrinsics_.size()) {
    return false;
  }
  intrinsics_ = intrinsics;
  return true;
}

int EquirectangularCamera::numIntrinsicsParameters() const
{
  return intrinsics_.size();
}

float_t EquirectangularCamera::horizontalFov() const
{
  return intrinsics_[0];
}

float_t EquirectangularCamera::verticalFov() const
{
  return intrinsics_[1];
}

//////////////////////////////////////////
// Methods to project points

ProjectionStatus EquirectangularCamera::project(const Vector3f& point, Vector2f* imagePoint) const
{
  return project(point, imagePoint, nullptr, nullptr);
}

ProjectionStatus EquirectangularCamera::project(const Vector3f& point,
                                          Vector2f* imagePoint,
                                          Matrixf<2, 3>* pointJacobian,
                                          Matrix2Xf* intrinsicsJacobian) const
{
  assert(imagePoint);
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

  // The equations the derivatives were computed from:
  //   imagePoint.x = imageWidth * (-atan2(point.y, point.x) / horizontalFov + 0.5)
  //   imagePoint.y = (imageHeight - 1) * acos(point.z / point.norm) / verticalFov
  //                  - (imageHeight - 1) * (π - verticalFov) / (2 * verticalFov)
  if (pointJacobian) {
    const float_t norm2_sq = point.head<2>().squaredNorm();
    // Derivatives of imagePoint.x w.r.t. point:
    // J[0,0] =  point.y * imageWidth / (horizontalFov * (point.x^2 + point.y^2))
    // J[0,1] = -point.x * imageWidth / (horizontalFov * (point.x^2 + point.y^2))
    const float_t temp0 = imageWidth_ / (horizontalFov() * norm2_sq);
    (*pointJacobian)(0, 0) = point.y() * temp0;
    (*pointJacobian)(0, 1) = -point.x() * temp0;
    (*pointJacobian)(0, 2) = 0;
    // Derivatives of imagePoint.y w.r.t. point:
    // J[1,0] = point.x * point.z * (imageHeight - 1) / (verticalFov * point.norm^2 * sqrt(point.x^2 + point.y^2))
    // J[1,1] = point.y * point.z * (imageHeight - 1) / (verticalFov * point.norm^2 * sqrt(point.x^2 + point.y^2))
    // J[1,2] = (imageHeight_ - 1) * sqrt(point.x^2 + point.y^2) / (verticalFov * point.norm^2)
    const int h = imageHeight_ - 1;
    const float_t norm2 = std::sqrt(norm2_sq);
    const float_t temp1 = verticalFov() * point.squaredNorm();
    const float_t temp2 = point.z() * h / (temp1 * norm2);
    (*pointJacobian)(1, 0) = point.x() * temp2;
    (*pointJacobian)(1, 1) = point.y() * temp2;
    (*pointJacobian)(1, 2) = h * norm2 / temp1;
  }
  if (intrinsicsJacobian) {
    // Derivatives of imagePoint.x w.r.t. horizontalFov:
    // J[0,0] = imageWidth * atan2(point.y, point.x) / horizontalFov^2
    (*intrinsicsJacobian)(0, 0) = imageWidth_ * azimuth / (horizontalFov() * horizontalFov());
    (*intrinsicsJacobian)(0, 1) = 0;
    // Derivatives of imagePoint.x w.r.t. verticalFov:
    // J[1,1] = (imageHeight - 1) * (π - 2 * acos(point.z / point.norm)) / (2 * verticalFov^2)
    (*intrinsicsJacobian)(1, 0) = 0;
    (*intrinsicsJacobian)(1, 1) = (imageHeight_ - 1) * (M_PI - inclination)
      / (2 * verticalFov() * verticalFov());
    intrinsicsJacobian->resize(2, numIntrinsicsParameters());
  }

  if (ProjectionBase::isMasked(*imagePoint)) {
    return ProjectionStatus::Masked;
  }
  return ProjectionStatus::Successful;
}

ProjectionStatus EquirectangularCamera::projectWithExternalParameters(const Vector3f& point,
                                                                const VectorXf& parameters,
                                                                Vector2f* imagePoint,
                                                                Matrixf<2, 3>* pointJacobian,
                                                                Matrix2Xf* intrinsicsJacobian) const
{
  if (parameters.size() != numIntrinsicsParameters()) {
    return ProjectionStatus::Invalid;
  }
  return EquirectangularCamera(imageWidth(), imageHeight(), parameters[0], parameters[1]).project(point,
      imagePoint, pointJacobian, intrinsicsJacobian);
}

void EquirectangularCamera::projectBatch(const Matrix3Xf& points,
                                   Matrix2Xf* imagePoints,
                                   std::vector<ProjectionStatus>* stati) const
{
  assert(imagePoints);
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

ProjectionStatus EquirectangularCamera::projectHomogeneous(const Vector4f& point,
                                                     Vector2f* imagePoint) const
{
  if (point.w() < 0) {
    return project(-point.head<3>(), imagePoint);
  } else {
    return project(point.head<3>(), imagePoint);
  }
}

ProjectionStatus EquirectangularCamera::projectHomogeneous(const Vector4f& point,
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
  if (pointJacobian) {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
    pointJacobian->template leftCols<3>() = pointJacobian3;
#pragma GCC diagnostic pop
    pointJacobian->template rightCols<1>() = Vector2f::Zero();
  }
  return status;
}

ProjectionStatus EquirectangularCamera::projectHomogeneousWithExternalParameters(
    const Vector4f& point,
    const VectorXf& parameters,
    Vector2f* imagePoint,
    Matrixf<2, 4>* pointJacobian,
    Matrix2Xf* intrinsicsJacobian) const
{
  if (parameters.size() != numIntrinsicsParameters()) {
    return ProjectionStatus::Invalid;
  }
  return EquirectangularCamera(imageWidth(), imageHeight(), parameters[0],
      parameters[1]).projectHomogeneous(point, imagePoint, pointJacobian, intrinsicsJacobian);
}

void EquirectangularCamera::projectHomogeneousBatch(const Matrix4Xf& points,
                                              Matrix2Xf* imagePoints,
                                              std::vector<ProjectionStatus>* stati) const
{
  assert(imagePoints);
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

bool EquirectangularCamera::backProject(const Vector2f& imagePoint, Vector3f* direction) const
{
  return backProject(imagePoint, direction, nullptr);
}

inline bool EquirectangularCamera::backProject(const Vector2f& imagePoint,
                                         Vector3f* direction,
                                         Matrixf<3, 2>* pointJacobian) const
{
  assert(direction);
  // Azimuth is horizontalFov/2 at the left edge of the image and -horizontalFov/2 at the right
  // edge. Normalize with imageWidth_ so that the azimuth is in the interval
  // (-horizontalFov/2, horizontalFov/2].
  const float_t azimuth = -horizontalFov() * (imagePoint.x() / imageWidth_ - 0.5f);
  // Inclination is min_inclination at the top of the image and min_inclination + verticalFov at
  // the bottom. Normalize with imageHeight_ - 1 so that the inclination is in the interval
  // [min_inclination, min_inclination + verticalFov].
  const float_t min_inclination = (M_PI - verticalFov()) / 2;
  const float_t inclination = verticalFov() * imagePoint.y() / (imageHeight_ - 1) + min_inclination;
  const float_t sin_azimuth = std::sin(azimuth);
  const float_t cos_azimuth = std::cos(azimuth);
  const float_t sin_inclination = std::sin(inclination);
  const float_t cos_inclination = std::cos(inclination);
  direction->x() = sin_inclination * cos_azimuth;
  direction->y() = sin_inclination * sin_azimuth;
  direction->z() = cos_inclination;

  // The equations the derivatives were computed from:
  //   point.x = sin(inclination) * cos(azimuth)
  //   point.y = sin(inclination) * sin(azimuth)
  //   point.z = cos(inclination)
  // with
  //   inclination = verticalFov * imagePoint.y / (imageHeight_ - 1) + (π - verticalFov) / 2
  //   azimuth = -horizontalFov * (imagePoint.x / imageWidth - 0.5)
  if (pointJacobian) {
    const float_t dinclination_dy = verticalFov() / (imageHeight_ - 1);
    const float_t dazimuth_dx = -horizontalFov() / imageWidth_;
    // Derivatives of point.x w.r.t. imagePoint:
    // J[0,0] = -sin(azimuth) * ∂azimuth/∂x
    // J[0,1] = cos(inclination) * ∂inclination/∂y
    (*pointJacobian)(0, 0) = -sin_azimuth * dazimuth_dx;
    (*pointJacobian)(0, 1) = cos_inclination * dinclination_dy;
    // Derivatives of point.y w.r.t. imagePoint:
    // J[1,0] = cos(inclination) * ∂inclination/∂y
    // J[1,1] = cos(inclination) * ∂inclination/∂y
    (*pointJacobian)(1, 0) = cos_azimuth * dazimuth_dx;
    (*pointJacobian)(1, 1) = (*pointJacobian)(0, 1);
    // Derivatives of point.z w.r.t. imagePoint:
    // J[2,1] = -sin(inclination) * ∂inclination/∂y
    (*pointJacobian)(2, 0) = 0;
    (*pointJacobian)(2, 1) = -sin_inclination * dinclination_dy;
  }
  return true;
}

bool EquirectangularCamera::backProjectBatch(const Matrix2Xf& imagePoints,
                                       Matrix3Xf* directions,
                                       std::vector<bool>* success) const
{
  assert(directions);
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

bool EquirectangularCamera::backProjectHomogeneous(const Vector2f& imagePoint, Vector4f* direction) const
{
  assert(direction);
  Vector3f ray;
  const bool success = backProject(imagePoint, &ray);
  *direction = ray.homogeneous();
  return success;
}

bool EquirectangularCamera::backProjectHomogeneous(const Vector2f& imagePoint,
                                             Vector4f* direction,
                                             Matrixf<4, 2>* pointJacobian) const
{
  assert(direction);
  Vector3f ray;
  Matrixf<3, 2> pointJacobian3;
  const bool success = backProject(imagePoint, &ray, &pointJacobian3);
  *direction = ray.homogeneous();
  if (pointJacobian) {
    pointJacobian->template topRows<3>() = pointJacobian3;
    pointJacobian->template bottomRows<1>() = Vector2f::Zero();
  }
  return success;
}

bool EquirectangularCamera::backProjectHomogeneousBatch(const Matrix2Xf& imagePoints,
                                                  Matrix4Xf* directions,
                                                  std::vector<bool>* success) const
{
  assert(directions);
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

std::shared_ptr<ProjectionBase> EquirectangularCamera::createTestInstance()
{
  return std::shared_ptr<ProjectionBase>(new EquirectangularCamera(1024, 512, 2 * M_PI, M_PI));
}

std::string EquirectangularCamera::type() const
{
  return "EquirectangularCamera";
}

}  // namespace projection
}  // namespace srl
