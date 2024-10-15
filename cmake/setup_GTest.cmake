# This file includes logic for setting up GoogleTest
find_package(GTest)

if (NOT GTest_FOUND)
    message(STATUS
    "Preparing to download GTest from GitHub & configure its build. Don't "
    "worry about any mentions of Python in the following text -- that is "
    "only used internally by the Google test build"
    )

    # NOTE: it is idiomatic to use FetchContent with a git hash rather than the
    # name of a tag or branch (since the latter 2 can freely change). If we use
    # the name of a tag/branch, then CMake will query GitHub on every
    # subsequent build to check whether the tag/branch is still the same

    include(FetchContent)
    FetchContent_Declare(
      googletest # the url contains the hash for v1.13.0
      URL https://github.com/google/googletest/archive/b796f7d44681514f58a683a3a71ff17c94edb0c1.zip
    )

    # For Windows: Prevent overriding the parent project's compiler/linker settings
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)

elseif("${CMAKE_VERSION}" VERSION_LESS "3.20")
    # CMake's built-in `FindGTest` module imported targets have different names
    # in earlier CMake versions
    add_library(GTest::gtest ALIAS GTest::GTest)
    add_library(GTest::gtest_main ALIAS GTest::Main)
endif()

if (NOT TARGET GTest::gtest_main)
    message("Target GTest:: stuff MISSING")
endif()
