import os
from typing import NamedTuple

import pytest

from testing_common import generate_test_results, test_answers_dir


# this hook is used to add more command line flags to the pytest launcher
def pytest_addoption(parser):
    parser.addoption(
        "--answer-store",
        action="store_true",
        help="Indicates whether to generate test results.",
    )

    # in the future, I think we should revisit the default behavior when this
    # is omitted. Instead of providing a default location, I think we should
    # just skip all answer tests (and maybe raise an error if the
    # --answer-store flag is provided without this flag)
    parser.addoption(
        "--local-dir",
        action="store",
        default=test_answers_dir,
        help="Path to directory where answers are/will be stored.",
    )


class AnswerTestSpec(NamedTuple):
    generate_answers: bool
    answer_dir: str


@pytest.fixture(scope="session")
def answertestspec(pytestconfig):
    """
    Return an object specifying all user-specified directives regarding
    answer tests (whether to run them and where to store/find results)
    """
    generate_answers = pytestconfig.getoption("--answer-store") or generate_test_results
    answer_dir = pytestconfig.getoption("--local-dir", default=test_answers_dir)

    if (not os.path.isdir(answer_dir)) and (not generate_answers):
        pytest.skip(f"the directory of test answers can't be found, {answer_dir}")
    os.makedirs(answer_dir, exist_ok=True)

    return AnswerTestSpec(generate_answers=generate_answers, answer_dir=answer_dir)
