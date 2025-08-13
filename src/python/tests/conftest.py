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
        "--answer-dir",
        action="store",
        default=None,
        help=(
            "Path to directory used for storing answer-tests answers. All "
            "logic related to answer-tests is disabled unless this is "
            "specified. When specified alonside --answer-store, logic to "
            "generate the answers will be run (the answers are stored to this "
            "directory). When specified WITHOUT --answer-store, the results "
            "of each answer-test is compared against the answers that were "
            "previously stored in this directory."
        )
    )

    parser.addoption(
        "--answer-store",
        action="store_true",
        help=(
            "Indicates that we should store answer-test results. It is an "
            "error to specifiy this arg without --answer-dir."
        ),
    )

    parser.addoption(
        "--model-comparison-dir",
        action="store",
        default=None,
        help=(
            "Path to directory used to store additional comparison of model "
            "tests between two versions of the code. This is only used when "
            "running in comparison mode, i.e., when --answer-store has not "
            "been given."
        ),
    )


class AnswerTestSpec(NamedTuple):
    generate_answers: bool
    answer_dir: str
    model_comparison_dir: str


@pytest.fixture(scope="session")
def answertestspec(request):
    """
    Return an object specifying all user-specified directives regarding
    answer tests (whether to run them and where to store/find results)
    """

    generate_answers = request.config.getoption("--answer-store")

    answer_dir = request.config.getoption("--answer-dir")

    model_comparison_dir = request.config.getoption("--model-comparison-dir")

    if (answer_dir is None) and generate_answers:
        raise ValueError(
            "--answer-store option can't be specified without --answer-dir"
        )
    elif answer_dir is None:
        pytest.skip("no --answer-dir option found")
    elif (not os.path.isdir(answer_dir)) and (not generate_answers):
        pytest.skip(f"directory of test answers can't be found, {answer_dir}")
    else:
        os.makedirs(answer_dir, exist_ok=True)

    if model_comparison_dir is not None:
        if generate_answers:
            raise ValueError(
                "--answer-store option can't be specified with "
                "--model-comparison-dir"
            )
        os.makedirs(model_comparison_dir, exist_ok=True)

    return AnswerTestSpec(generate_answers=generate_answers, answer_dir=answer_dir,
                          model_comparison_dir=model_comparison_dir)
