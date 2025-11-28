# SPDX-FileCopyrightText: 2025 Smart Robotics Lab, Imperial College London, Technical University of Munich
# SPDX-License-Identifier: CC0-1.0
cmake_minimum_required(VERSION 3.1...3.16)

# Disable warnings for third-party code
add_compile_options(-w)

if(SRL_PROJECTION_BUILD_TESTS AND NOT TARGET gtest_main)
    set(CMAKE_POLICY_DEFAULT_CMP0069 NEW)
    option(BUILD_GMOCK OFF)
    option(INSTALL_GTEST OFF)
    add_subdirectory(googletest)
endif()
