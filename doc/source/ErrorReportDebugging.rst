.. _err_report_debug:

Error Reporting and Debugging
=============================

Like with most computational software, one of the first steps of dealing with most (non-obvious) Grackle-related issue is to identify a minimal example that reproduces the issue (whether you want to simply get help OR you try to debug the underlying issue).

Unfortunately, a few factors have historically complicated this process:

  * Grackle has a **LOT** of input parameters

  * Grackle is usually integrated inside an application that is being executed across many processors (e.g. with MPI)


To help with identifying a minimal example that reproduces your issue, we have introduced the :c:func:`grunstable_h5dump_state`.
This function is most useful when you can reliably and repeatedly reproduce the underlying issues.

This function essentially dumps the current state of Grackle (and the input fields) to an hdf5 file.
Therefore, you should make sure to dump the state right before the error gets triggered.
