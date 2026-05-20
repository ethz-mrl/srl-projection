//
// Created by boche on 5/5/22.
//
/**
 * @file projection/Lidar.hpp
 * @brief Header file for the Lidar class.
 * @author Simon Boche
 */

#ifndef INCLUDE_SRL_PROJECTION_LIDAR_HPP_
#define INCLUDE_SRL_PROJECTION_LIDAR_HPP_

#include <vector>
#include <memory>
#include <stdint.h>
#include <Eigen/Core>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnon-virtual-dtor"
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/core/core.hpp> // Code that causes warning goes here
#pragma GCC diagnostic pop
#include "srl/projection/ProjectionBase.hpp"
#include "srl/projection/DistortionBase.hpp"
#include "srl/projection/NoDistortion.hpp"

/// \brief Main namespace of this package.
namespace srl {
/// \brief Namespace for camera-related functionality.
namespace projection {

/// \class Lidar
/// \brief This implements the Lidar projection model.
struct Lidar : public ProjectionBase
{
    /// \brief Constructor that will figure out the type of distortion
    /// @param[in] imageWidth The width in pixels.
    /// @param[in] imageHeight The height in pixels.
    /// @param[in] beamAzimuthAngles The azimuth start angles per scan row.
    /// @param[in] beamElevationAngles The elevation angle per scan row.
    inline Lidar(const int imageWidth, const int imageHeight);

    /// \brief Destructor.
    virtual ~Lidar() = default;

    /// \brief Get the intrinsics as a concatenated vector.
    /// \param[out] intrinsics The intrinsics as a concatenated vector.
    inline void getIntrinsics(VectorXf & intrinsics) const override;

    /// \brief overwrite all intrinsics - use with caution !
    /// \param[in] intrinsics The intrinsics as a concatenated vector.
    inline bool setIntrinsics(const VectorXf & intrinsics) override;

    /// \brief Get the total number of intrinsics.
    /// \return Number of intrinsics parameters.
    int numIntrinsicsParameters() const override
    {
        return -1;
    }

    //////////////////////////////////////////////////////////////
    /// \name Methods to project points
    /// @{

    /// \brief Projects a Euclidean point to a 2d image point (projection).
    ///        Uses projection including distortion models.
    /// @param[in]  point      The point in Euclidean coordinates.
    /// @param[out] imagePoint The image point.
    /// @return     Get information about the success of the projection. See
    ///             \ref ProjectionStatus for more information.
    inline ProjectionStatus project(
        const Vector3f & point, Vector2f * imagePoint) const override;

    /// \brief Projects a Euclidean point to a 2d image point (projection).
    ///        Uses projection including distortion models.
    /// @param[in]  point              The point in Euclidean coordinates.
    /// @param[out] imagePoint         The image point.
    /// @param[out] pointJacobian      The Jacobian of the projection function w.r.t. the point..
    /// @param[out] intrinsicsJacobian The Jacobian of the projection function w.r.t. the intinsics.
    /// @return     Get information about the success of the projection. See
    ///             \ref ProjectionStatus for more information.
    inline ProjectionStatus project(
        const Vector3f & point, Vector2f * imagePoint,
        Matrixf<2, 3> * pointJacobian,
        Matrix2Xf * intrinsicsJacobian = nullptr) const override;

    inline ProjectionStatus projectSphere(
        const Vector3f & center, float radius, Vector2f * imageCenter, float& imageRadius) const;

    /// \brief Projects a Euclidean point to a 2d image point (projection).
    ///        Uses projection including distortion models.
    /// @param[in]  point              The point in Euclidean coordinates.
    /// @param[in]  parameters         The intrinsics.
    /// @param[out] imagePoint         The image point.
    /// @param[out] pointJacobian      The Jacobian of the projection function w.r.t. the point..
    /// @param[out] intrinsicsJacobian The Jacobian of the projection function w.r.t. the intinsics.
    /// @return     Get information about the success of the projection. See
    ///             \ref ProjectionStatus for more information.
    inline ProjectionStatus projectWithExternalParameters(
        const Vector3f & point, const VectorXf & parameters,
        Vector2f * imagePoint, Matrixf<2, 3> * pointJacobian,
        Matrix2Xf * intrinsicsJacobian = nullptr) const override;

    /// \brief Projects a point in homogenous coordinates to a 2d image point (projection).
    ///        Uses projection including distortion models.
    /// @param[in]  point      The point in Homogeneous coordinates.
    /// @param[out] imagePoint The image point.
    /// @return     Get information about the success of the projection. See
    ///             \ref ProjectionStatus for more information.
    inline ProjectionStatus projectHomogeneous(
        const Vector4f & point, Vector2f * imagePoint) const override;

    /// \brief Projects a point in homogenous coordinates to a 2d image point (projection).
    ///        Uses projection including distortion models.
    /// @param[in]  point              The point in Homogeneous coordinates.
    /// @param[out] imagePoint         The image point.
    /// @param[out] pointJacobian      The Jacobian of the projection function w.r.t. the point.
    /// @param[out] intrinsicsJacobian The Jacobian of the projection function w.r.t. the intrinsics.
    /// @return     Get information about the success of the projection. See
    ///             \ref ProjectionStatus for more information.
    inline ProjectionStatus projectHomogeneous(
        const Vector4f & point, Vector2f * imagePoint,
        Matrixf<2, 4> * pointJacobian,
        Matrix2Xf * intrinsicsJacobian = nullptr) const override;

    /// \brief Projects a point in homogenous coordinates to a 2d image point (projection).
    ///        Uses projection including distortion models.
    /// @param[in]  point              The point in Homogeneous coordinates.
    /// @param[in]  parameters         The intrinsics.
    /// @param[out] imagePoint         The image point.
    /// @param[out] pointJacobian      The Jacobian of the projection function w.r.t. the point.
    /// @param[out] intrinsicsJacobian The Jacobian of the projection function w.r.t. the intrinsics.
    /// @return     Get information about the success of the projection. See
    ///             \ref ProjectionStatus for more information.
    inline ProjectionStatus projectHomogeneousWithExternalParameters(
        const Vector4f & point, const VectorXf & parameters,
        Vector2f * imagePoint,
        Matrixf<2, 4> * pointJacobian = nullptr,
        Matrix2Xf * intrinsicsJacobian = nullptr) const override;
    /// @}

    //////////////////////////////////////////////////////////////
    /// \name Methods to backproject points
    /// @{

    /// \brief Back-project a 2d image point into Euclidean space (direction vector).
    /// @param[in]  imagePoint The image point.
    /// @param[out] direction  The Euclidean direction vector.
    /// @return     true on success.
    inline bool backProject(const Vector2f & imagePoint,
                            Vector3f * direction) const override;

    /// \brief Back-project a 2d image point into Euclidean space (direction vector).
    /// @param[in]  imagePoint         The image point.
    /// @param[out] direction          The Euclidean direction vector.
    /// @param[out] pointJacobian      Jacobian of the back-projection function  w.r.t. the point.
    /// @return     true on success.
    inline bool backProject(const Vector2f & imagePoint,
                            Vector3f * direction,
                            Matrixf<3, 2> * pointJacobian) const override;

    /// \brief Back-project a 2d image point into homogeneous point (direction vector).
    /// @param[in]  imagePoint The image point.
    /// @param[out] direction  The homogeneous point as direction vector.
    /// @return     true on success.
    inline bool backProjectHomogeneous(const Vector2f & imagePoint,
                                       Vector4f * direction) const override;

    /// \brief Back-project a 2d image point into homogeneous point (direction vector).
    /// @param[in]  imagePoint         The image point.
    /// @param[out] direction          The homogeneous point as direction vector.
    /// @param[out] pointJacobian      Jacobian of the back-projection function.
    /// @return     true on success.
    inline bool backProjectHomogeneous(
        const Vector2f & imagePoint, Vector4f * direction,
        Matrixf<4, 2> * pointJacobian) const override;
    /// @}

    /// \brief get a test instance
    static inline std::shared_ptr<ProjectionBase> createTestObject();
    /// \brief get a test instance
    static inline Lidar testObject();

    /// \brief Obtain the projection type
    std::string type() const
    {
        return "Lidar";
    }

    inline float azimuthResolutionAngle() const;
    inline void setAzimuthResolutionAngle(const float azimuthResolutionAngle);

    inline float elevationResolutionAngle() const;
    inline void setElevationResolutionAngle(const float elevationResolutionAngle);

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    protected:

    /// \brief No default constructor.
    Lidar() = delete;

    float azimuthResolution_;
    float elevationResolution_;
};

}  // namespace projection
}  // namespace srl

#include "implementation/Lidar.hpp"

#endif /* INCLUDE_SRL_PROJECTION_LIDAR_HPP_ */
