# define some hook-functions that will customize pytest's behavior

from pygrackle.utilities.data_path import _download_all_datafiles


def pytest_sessionstart(session):
    # this is a hook that is called just before collecting tests and entering
    # the test loop.

    # All we want to do is make sure that we have all of the data files that we
    # need downloaded (This might not be the right place to put this logic)
    _download_all_datafiles()
