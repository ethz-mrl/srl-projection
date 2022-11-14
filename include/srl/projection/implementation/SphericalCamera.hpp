/*********************************************************************************
 *  Copyright (c) 2022, Smart Robotics Lab / Imperial College London
 *  Copyright (c) 2022, Smart Robotics Lab / Technical University of Munich
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *   * Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *   * Neither the name of Autonomous Systems Lab / ETH Zurich nor the names of
 *     its contributors may be used to endorse or promote products derived from
 *     this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *  Created on: Nov 14, 2022
 *      Author: Sotiris Papatheodorou (s.leutenegger@imperial.ac.uk)
 *********************************************************************************/

/**
 * @file implementation/SphericalCamera.hpp
 * @brief Header implementation file for the SphericalCamera class.
 * @author Sotiris Papatheodorou
 */

#include <Eigen/Geometry>

namespace srl {
namespace projection {

SphericalCamera::SphericalCamera(int imageWidth, int imageHeight, float horizontalFov, float verticalFov)
    : ProjectionBase(imageWidth, imageHeight),
      horizontalFov_(horizontalFov),
      verticalFov_(verticalFov)
{
}

void SphericalCamera::getIntrinsics(VectorXf&) const {
}

bool SphericalCamera::setIntrinsics(const VectorXf&) {
  return true;
}

int SphericalCamera::numIntrinsicsParameters() const
{
  return 0;
}

float SphericalCamera::horizontalFov() const
{
  return horizontalFov_;
}

float SphericalCamera::verticalFov() const
{
  return verticalFov_;
}

//////////////////////////////////////////
// Methods to project points

ProjectionStatus SphericalCamera::project(const Vector3f& point, Vector2f* imagePoint) const
{
  const float_t norm = point.norm();
  if (norm < 1e-5 || (std::fabs(point.x()) < 1e-5 && std::fabs(point.y()) < 1e-5)) {
    return ProjectionStatus::Invalid;
  }
  const float_t azimuth = std::atan2(point.y(), point.x());
  if (azimuth > horizontalFov_ / 2 || azimuth < -horizontalFov_ / 2) {
    return ProjectionStatus::OutsideImage;
  }
  const float_t inclination = std::acos(point.z() / norm);
  const float_t min_inclination = (M_PI - verticalFov_) / 2;
  if (inclination < min_inclination || inclination > min_inclination + verticalFov_) {
    return ProjectionStatus::OutsideImage;
  }
  imagePoint->x() = imageWidth_ * (-azimuth / horizontalFov_ + 0.5);
  imagePoint->y() = (imageHeight_ - 1) * (inclination - min_inclination) / verticalFov_;
  if (ProjectionBase::isMasked(*imagePoint)) {
    return ProjectionStatus::Masked;
  }
  return ProjectionStatus::Successful;
}

ProjectionStatus SphericalCamera::project(const Vector3f& /*point*/,
                                          Vector2f* /*imagePoint*/,
                                          Matrixf<2, 3>* /*pointJacobian*/,
                                          Matrix2Xf* /*intrinsicsJacobian*/) const
{
  throw std::runtime_error("not implemented");
}

ProjectionStatus SphericalCamera::projectWithExternalParameters(const Vector3f& /*point*/,
                                                                const VectorXf& /*parameters*/,
                                                                Vector2f* /*imagePoint*/,
                                                                Matrixf<2, 3>* /*pointJacobian*/,
                                                                Matrix2Xf* /*intrinsicsJacobian*/) const
{
  throw std::runtime_error("not implemented");
}

void SphericalCamera::projectBatch(const Matrix3Xf& points,
                                   Matrix2Xf* imagePoints,
                                   std::vector<ProjectionStatus>* stati) const
{
  const int numPoints = points.cols();
  for (int i = 0; i < numPoints; ++i) {
    Vector3f point = points.col(i);
    Vector2f imagePoint;
    ProjectionStatus status = project(point, &imagePoint);
    imagePoints->col(i) = imagePoint;
    if (stati) {
      stati->push_back(status);
    }
  }
}

ProjectionStatus SphericalCamera::projectHomogeneous(const Vector4f& point,
                                                     Vector2f* imagePoint) const
{
  Vector3f head = point.head<3>();
  if (point[3] < 0) {
    return project(-head, imagePoint);
  } else {
    return project(head, imagePoint);
  }
}

ProjectionStatus SphericalCamera::projectHomogeneous(const Vector4f& point,
                                                     Vector2f* imagePoint,
                                                     Matrixf<2, 4>* pointJacobian,
                                                     Matrix2Xf* intrinsicsJacobian) const
{
  Vector3f head = point.head<3>();
  Matrixf<2, 3> pointJacobian3;
  ProjectionStatus status;
  if (point[3] < 0) {
    status = project(-head, imagePoint, &pointJacobian3, intrinsicsJacobian);
  } else {
    status = project(head, imagePoint, &pointJacobian3, intrinsicsJacobian);
  }
  pointJacobian->template bottomRightCorner<2, 1>() = Vector2f::Zero();
  pointJacobian->template topLeftCorner<2, 3>() = pointJacobian3;
  return status;
}

ProjectionStatus SphericalCamera::projectHomogeneousWithExternalParameters(
    const Vector4f& point,
    const VectorXf& parameters,
    Vector2f* imagePoint,
    Matrixf<2, 4>* pointJacobian,
    Matrix2Xf* intrinsicsJacobian) const
{
  Vector3f head = point.head<3>();
  Matrixf<2, 3> pointJacobian3;
  ProjectionStatus status;
  if (point[3] < 0) {
    status = projectWithExternalParameters(-head, parameters, imagePoint, &pointJacobian3,
        intrinsicsJacobian);
  } else {
    status = projectWithExternalParameters(head, parameters, imagePoint, &pointJacobian3,
        intrinsicsJacobian);
  }
  pointJacobian->template bottomRightCorner<2, 1>() = Vector2f::Zero();
  pointJacobian->template topLeftCorner<2, 3>() = pointJacobian3;
  return status;
}

void SphericalCamera::projectHomogeneousBatch(const Matrix4Xf& points,
                                              Matrix2Xf* imagePoints,
                                              std::vector<ProjectionStatus>* stati) const
{
  const int numPoints = points.cols();
  for (int i = 0; i < numPoints; ++i) {
    Vector4f point = points.col(i);
    Vector2f imagePoint;
    ProjectionStatus status = projectHomogeneous(point, &imagePoint);
    imagePoints->col(i) = imagePoint;
    if (stati) {
      stati->push_back(status);
    }
  }
}

//////////////////////////////////////////
// Methods to backproject points

bool SphericalCamera::backProject(const Vector2f& imagePoint, Vector3f* direction) const
{
  // Azimuth is horizontalFov_ at the left edge of the image and -horizontalFov_ at the right edge.
  // Normalize with imageWidth_ so that the azimuth is in the interval
  // (-horizontalFov_, horizontalFov_].
  const float_t azimuth = -horizontalFov_ * (imagePoint.x() / imageWidth_ - 0.5f);
  // Inclination is min_inclination at the top of the image and min_inclination + verticalFov_ at
  // the bottom. Normalize with imageHeight_ - 1 so that the inclination is in the interval
  // [min_inclination, min_inclination + verticalFov_].
  const float_t min_inclination = (M_PI - verticalFov_) / 2;
  const float_t inclination = verticalFov_ * imagePoint.y() / (imageHeight_ - 1) + min_inclination;
  direction->x() = std::sin(inclination) * std::cos(azimuth);
  direction->y() = std::sin(inclination) * std::sin(azimuth);
  direction->z() = std::cos(inclination);
  return true;
}

inline bool SphericalCamera::backProject(const Vector2f& /*imagePoint*/,
                                         Vector3f* /*direction*/,
                                         Matrixf<3, 2>* /*pointJacobian*/) const
{
  throw std::runtime_error("not implemented");
}

bool SphericalCamera::backProjectBatch(const Matrix2Xf& imagePoints,
                                       Matrix3Xf* directions,
                                       std::vector<bool>* success) const
{
  const int numPoints = imagePoints.cols();
  directions->row(3) = VectorXf::Ones(numPoints);
  for (int i = 0; i < numPoints; ++i) {
    Vector2f imagePoint = imagePoints.col(i);
    Vector3f point;
    bool suc = backProject(imagePoint, &point);
    if (success) {
      success->push_back(suc);
    }
    directions->col(i) = point;
  }
  return true;
}

bool SphericalCamera::backProjectHomogeneous(const Vector2f& imagePoint, Vector4f* direction) const
{
  Vector3f ray;
  bool success = backProject(imagePoint, &ray);
  direction->template head<3>() = ray;
  (*direction)[3] = 1.0;  // arbitrary
  return success;
}

bool SphericalCamera::backProjectHomogeneous(const Vector2f& imagePoint,
                                             Vector4f* direction,
                                             Matrixf<4, 2>* pointJacobian) const
{
  Vector3f ray;
  Matrixf<3, 2> pointJacobian3;
  bool success = backProject(imagePoint, &ray, &pointJacobian3);
  direction->template head<3>() = ray;
  (*direction)[3] = 1.0;  // arbitrary
  pointJacobian->template bottomRightCorner<1,2>() = Vector2f::Zero();
  pointJacobian->template topLeftCorner<3, 2>() = pointJacobian3;
  return success;
}

bool SphericalCamera::backProjectHomogeneousBatch(const Matrix2Xf& imagePoints,
                                                  Matrix4Xf* directions,
                                                  std::vector<bool>* success) const
{
  const int numPoints = imagePoints.cols();
  directions->row(3) = VectorXf::Ones(numPoints);
  for (int i = 0; i < numPoints; ++i) {
    Vector2f imagePoint = imagePoints.col(i);
    Vector3f point;
    bool suc = backProject(imagePoint, &point);
    if (success) {
      success->push_back(suc);
    }
    directions->template block<3, 1>(0, i) = point;
  }
  return true;
}

std::shared_ptr<ProjectionBase> SphericalCamera::createTestObject()
{
  return std::shared_ptr<ProjectionBase>(new SphericalCamera(1024, 512, 2 * M_PI, M_PI));
}

SphericalCamera SphericalCamera::testObject()
{
  return SphericalCamera(1024, 512, 2 * M_PI, M_PI);
}

std::string SphericalCamera::type() const
{
  return "SphericalCamera";
}

}  // namespace projection
}  // namespace srl
