# declare system tests (aka end2end tests)
# ========================================

# NOTE: the testdeps target is defined within ../grtestutils/CMakeLists.txt

add_executable(runSystemTests test_ghost_zone.cpp)
target_link_libraries(runSystemTests testdeps)
gtest_discover_tests(runSystemTests)
