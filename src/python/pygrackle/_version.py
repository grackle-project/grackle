# pyproject.toml currently reads in this information
# -> in a simple project, the setuptools_scm plugin could be used to infer the
#    version number from git and then generate a _version.py file
#    https://scikit-build-core.readthedocs.io/en/latest/configuration.html#dynamic-metadata
# -> currently, this isn't viable for pygrackle since it lives in a "monorepo"
#    and it has a different version number from everything else in the repository.
#    There was some recent discussion discussion about addressing this scenario
#    in setuptools_scm, but it probably won't happen very soon
#    https://github.com/pypa/setuptools_scm/issues/1056

__all__ = ["__version__"]

__version__ = "1.1.dev"
