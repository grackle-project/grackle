#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <memory>  // unique_ptr
#include <string>

#include "grackle.h"  // GR_FAIL
#include "status_reporting.h"

#include "grtestutils/os.hpp"

TEST(StatusReportingTest, PrintErrMsg) {
  // setup capture of stderr
  std::unique_ptr<grtest::CaptureSink> capstderr =
      grtest::CaptureSink::create(stderr);
  if (capstderr.get() == nullptr) {
    GTEST_SKIP() << "Stream redirection doesn't seem to be working. It may "
                 << "not be implemented for the current platform";
  }

  GrPrintErrMsg("message with string: `%s` & int: %d", "str", 4);
  // end the capture
  std::string msg = capstderr->end_capture();

  EXPECT_THAT(msg, testing::EndsWith("\n"));
  EXPECT_THAT(msg, testing::HasSubstr("message with string: `str` & int: 4"));
}

TEST(StatusReportingTest, PrintAndReturnErr) {
  // setup capture of stderr
  std::unique_ptr<grtest::CaptureSink> capstderr =
      grtest::CaptureSink::create(stderr);
  if (capstderr.get() == nullptr) {
    GTEST_SKIP() << "Stream redirection doesn't seem to be working. It may "
                 << "not be implemented for the current platform";
  }

  int ec = GrPrintAndReturnErr("message with string: `%s` & int: %d", "str", 4);

  // end the capture
  std::string msg = capstderr->end_capture();

  EXPECT_EQ(ec, GR_FAIL) << "GrPrintAndReturnErr's exit code, " << ec
                         << ", doesn't match GR_FAIL (" << GR_FAIL << ')';
  EXPECT_THAT(msg, testing::EndsWith("\n"));
  EXPECT_THAT(msg, testing::HasSubstr("message with string: `str` & int: 4"));
}

TEST(StatusReportingDeathTest, InternalError) {
  // right now, we are mostly just confirm that the program aborts
  // -> we are checking that the specified error is included
  // -> in the future, we have the option to test the formatting of the message
  EXPECT_DEATH(
      { GR_INTERNAL_ERROR("This is a dummy error message %s", "dummy-str"); },
      ".*This is a dummy error message dummy-str.*");
}

TEST(StatusReportingDeathTest, InternalRequireFalse) {
  // right now, we are mostly just confirm that the program aborts
  // -> we are checking that the specified error is included
  // -> in the future, we have the option to test the formatting of the message
  EXPECT_DEATH(
      { GR_INTERNAL_REQUIRE(1 == 0, "This is a dummy error message %d", 10); },
      ".*This is a dummy error message 10.*");
}

TEST(StatusReportingDeathTest, InternalRequireTrue) {
  GR_INTERNAL_REQUIRE(1 == 1, "This program shouldn't abort");
}
