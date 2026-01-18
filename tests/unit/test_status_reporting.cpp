#include <gtest/gtest.h>

// the following 2 headers are the c versions of the headers (rather than the
// C++ versions) since it seems more likely that posix-specific functions are
// provided in this version
//#include <stdio.h> // fflush, tmpfile, fread, fileno
//#include <stdlib.h> // mkstemp, getenv

#include <memory>
#include <string>

#include "grackle.h" // GR_FAIL
#include "status_reporting.h"

#include "grtestutils/os.hpp"


testing::AssertionResult ContainsFormattedMessage_(int n) {
  if ((n % 2) == 0)
    return testing::AssertionSuccess();
  else
    return testing::AssertionFailure() << n << " is odd";
}


TEST(StatusReportingTest, PrintErrMsg) {

  // setup capture of stderr
  std::unique_ptr<grtest::CaptureSink> capstderr =
    grtest::CaptureSink::create(stderr);
  if (capstderr.get() == nullptr) {
    GTEST_SKIP() << "Stream redirection doesn't seem to be working. It may "
                 << "not be implemented for the current platform";
  }

  GrPrintErrMsg("message with string: `%s` & int: %d", "str", 4);
  const char* expected_msg = "message with string: `str` & int: 4";

  // end the capture
  std::string msg = capstderr->end_capture();

  // now do some checking (it would be cleaner to use gmock's matchers)
  ASSERT_GT(msg.size(), 0) << "nothing seems to have printed";
  EXPECT_EQ(msg[msg.size()-1], '\n') << "we expect a trailing newline";
  ASSERT_NE(msg.find(expected_msg), std::string::npos)
    << "The error message doesn't contain the substring \"" << expected_msg
    << "\". The contents of the message is: \"" << msg << '"';
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
  const char* expected_msg = "message with string: `str` & int: 4";

  // end the capture
  std::string msg = capstderr->end_capture();

  EXPECT_EQ(ec, GR_FAIL) << "GrPrintAndReturnErr's exit code, " << ec
                         << ", doesn't match GR_FAIL (" << GR_FAIL << ')';

  // now do some checking (it would be cleaner to use gmock's matchers)
  ASSERT_GT(msg.size(), 0) << "nothing seems to have printed";
  EXPECT_EQ(msg[msg.size()-1], '\n') << "we expect a trailing newline";
  ASSERT_NE(msg.find(expected_msg), std::string::npos)
    << "The error message doesn't contain the substring \"" << expected_msg
    << "\". The contents of the message is: \"" << msg << '"';
}

TEST(StatusReportingDeathTest, InternalError) {
  // right now, we are mostly just confirm that the program aborts
  // -> we are checking that the specified error is included
  // -> in the future, we have the option to test the formatting of the message
  EXPECT_DEATH({
    GR_INTERNAL_ERROR("This is a dummy error message %s", "dummy-str");
  }, ".*This is a dummy error message dummy-str.*");
}

TEST(StatusReportingDeathTest, InternalRequireFalse) {
  // right now, we are mostly just confirm that the program aborts
  // -> we are checking that the specified error is included
  // -> in the future, we have the option to test the formatting of the message
  EXPECT_DEATH({
    GR_INTERNAL_REQUIRE(1 == 0, "This is a dummy error message %d", 10);
  }, ".*This is a dummy error message 10.*");
}

TEST(StatusReportingDeathTest, InternalRequireTrue) {
  GR_INTERNAL_REQUIRE(1 == 1, "This program shouldn't abort");
}
