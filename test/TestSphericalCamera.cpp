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

#include <gtest/gtest.h>

#include "srl/projection/SphericalCamera.hpp"

TEST(SphericalCamera, functions)
{
  const std::shared_ptr<srl::projection::ProjectionBase> camera
      = srl::projection::SphericalCamera::createTestObject();

  ASSERT_FALSE(camera->hasMask());

  constexpr size_t num_points = 100;
  for (size_t i = 0; i < num_points; ++i) {
    // generate a random point and back project
    const srl::Vector2f imagePoint = camera->createRandomImagePoint();
    srl::Vector3f ray;
    ASSERT_TRUE(camera->backProject(imagePoint, &ray)) << "unsuccessful back projection";

    // randomise distance and project
    ray.normalize();
    ray *= (0.2f + 8 * (srl::Vector2f::Random()[0] + 1));
    srl::Vector2f imagePoint2;
    ASSERT_EQ(camera->project(ray, &imagePoint2), srl::projection::ProjectionStatus::Successful)
        << "unsuccessful projection";

    // check they are the same
    ASSERT_LT((imagePoint2 - imagePoint).norm(), 0.01) << "project/unproject failure";
  }
}
