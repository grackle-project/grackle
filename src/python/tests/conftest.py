# this defines some basic utilities shared among all of the tests
#
# DO NOT import pygrackle into global scope in this file (or in any file
# imported by this file). Some tests need to be runable without installing
# pygrackle

import os
from typing import NamedTuple

import pytest


# this hook is used to add more command line flags to the pytest launcher
def pytest_addoption(parser):
    parser.addoption(
        "--answer-skip",
        action="store_true",
        help="Indicates that we should skip all answer tests."
    )
    parser.addoption(
        "--answer-store",
        action="store_true",
        help="Indicates that we should generate test results.",
    )


    # in the future, I think we should revisit the default behavior when this
    # is omitted. Instead of providing a default location, I think we should
    # just skip all answer tests (and maybe raise an error if the
    # --answer-store flag is provided without this flag)
    parser.addoption(
        "--local-dir",
        action="store",
        default=os.path.join(os.path.abspath(__file__), "test_answers"),
        help="Path to directory where answers are/will be stored.",
    )


class AnswerTestSpec(NamedTuple):
    generate_answers: bool
    answer_dir: str


@pytest.fixture(scope="session")
def answertestspec(request):
    """
    Return an object specifying all user-specified directives regarding
    answer tests (whether to run them and where to store/find results)
    """
    if request.config.getoption("--answer-skip"):
        if request.config.getoption("--answer-store"):
            raise RuntimeError(
                "It is an error to specify both --answer-skip and --answer-store"
            )
        pytest.skip("--answer-skip was specified")

    generate_answers = (
        request.config.getoption("--answer-store")
        or int(os.environ.get("GENERATE_PYGRACKLE_TEST_RESULTS", 0)) == 1
    )
    answer_dir = request.config.getoption("--local-dir")


    if (not os.path.isdir(answer_dir)) and (not generate_answers):
        pytest.skip(f"the directory of test answers can't be found, {answer_dir}")
    os.makedirs(answer_dir, exist_ok=True)

    return AnswerTestSpec(generate_answers=generate_answers, answer_dir=answer_dir)
