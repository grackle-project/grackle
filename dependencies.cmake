# Handle external dependencies
# -> locate all existing prebuilt dependencies
# -> prepare other dependencies that must be built from source

if (GRACKLE_BUILD_TESTS)  # deal with testing dependencies

  # the only testing dependency is googletest. If we can't find a pre-installed
  # version of the library, we will fetch it and build from source
  find_package(GTest)

  if (NOT GTest_FOUND)
    message(STATUS
      "GTest not found. Fetching via FetchContent and configuring its build"
    )

    # NOTE: it is idiomatic to use FetchContent with a git hash rather than the
    # name of a tag or branch (since the latter 2 can freely change). If we use
    # the name of a tag/branch, then CMake will query GitHub on every
    # subsequent build to check whether the tag/branch is still the same

    include(FetchContent)
    FetchContent_Declare(
      googletest # the url contains the hash for v1.15.2
      URL https://github.com/google/googletest/archive/b514bdc898e2951020cbdca1304b75f5950d1f59.zip
    )

    # Tell GoogleTest's build-system not to define installation rules (since we
    # only use it to run tests from the build-directory)
    set(INSTALL_GTEST OFF)
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
endif() # GRACKLE_BUILD_TESTS
