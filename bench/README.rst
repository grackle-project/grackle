Benchmarking
============

This page documents benchmarks for Grackle.

**WARNING:** These benchmarks are all fairly experimental (we are mostly in the proof-of-concept phase).

These benchmarks are constructed using `Google Benchmark <https://github.com/google/benchmark/tree/main>`__ and make use of parts of the test-harness machinery for the C++ tests (i.e. namely the C++ wrappers and machinery to aide initialization).

At the time of writing, the benchmarks **ONLY** test Grackle's public API.
See the "Unanswered Questions" section at the end of this page for an extended discussion of this choice (and the ramifications for changing this decision).

Compiling the benchmarks
------------------------

**WARNING:** Since this is all experimental, the details may change!

When configuring the grackle-build, enable the ``GRACKLE_BUILD_BENCHMARK`` cmake-option.

If your build directory is called `<my-build-dir>`, then you would invoke:

.. code:: shell-session

   $ <my-build-dir>/bench/grackle-benchmarks

The result will look something look something like the following:

::

   Unable to determine clock rate from sysctl: hw.cpufrequency: No such file or directory
   This does not affect benchmark measurements, only the metadata output.
   ***WARNING*** Failed to set thread affinity. Estimated CPU frequency may be incorrect.
   2026-02-11T12:17:32-05:00
   Running ./build/bench/grackle-benchmarks
   Run on (12 X 24 MHz CPU s)
   CPU Caches:
     L1 Data 64 KiB
     L1 Instruction 128 KiB
     L2 Unified 4096 KiB (x12)
   Load Average: 2.42, 2.36, 2.25
   ---------------------------------------------------------------
   Benchmark                     Time             CPU   Iterations
   ---------------------------------------------------------------
   solve_chemistry/1         12744 ns        12747 ns        54470
   solve_chemistry/8         29175 ns        29173 ns        23949
   solve_chemistry/64       154456 ns       154450 ns         4518
   solve_chemistry/256      573804 ns       573801 ns         1219
   solve_chemistry/512     1167628 ns      1167652 ns          600
   solve_chemistry/1024    2322132 ns      2322095 ns          304
   solve_chemistry/2048    5147197 ns      5147243 ns          136
   solve_chemistry/4096   10334208 ns     10334281 ns           64

Obviously, the details in the heading and timings will change from machine to machine.
The listed benchmarks will probably evolve over time


Writing a good benchmark (and Grackle-specific challenges)
----------------------------------------------------------

Writing a good benchmark is a famously challenging problem.
Essentially, we want a nice repeatable operation that is representative of real-world behavior.

Lucky for us we mostly just care about profiling a small handful of Grackle's functions:

- most important: ``local_solve_chemsitry``
- secondary importance: ``local_calculate_temperature`` and ``local_calculate_cooling_time``
- tertiary importance: other ``local_calculate_<quan>`` functions

The challenge is that performance depends a lot upon inputs.
While there are lots of configuration parameters, that's managable.
The real issue concerns the actual field values.

TODO
----

- (obviously) address the unanswered questions
- Work on coming up with a nice minimal set of representative benchmarks
- On Linux, make it possible to gather hardware counters to measure:
  - instruction counts: apparently, the number of instruction counts invoked by a program is often more robust than wallclock time for detecting performance regression
  - CPU cycles
  - we may also want to track things like number of cache misses, number of branches, number of branch misses
- Introduce some kind of CI solution to regularly run benchmarks for every CI test
  - if nothing else, we should run it to make sure it doesn't break
  - ideally, we would make a GitHub action that is executed on PRs.
    - On PR creation it would run benchmarks on the current branch and the destination branch and post a comment with a table summarizing benchmark differences
    - Every time a new commit is pushed, the action would ideally rerun the benchmark and update the comment
    - CI should probably focus on comparing instruction count (since we can't control what else is running on the machine during the benchmark)
  - an example of a project that does something like this is `AirspeedVelocity.jl <https://github.com/MilesCranmer/AirspeedVelocity.jl>`__

Unanswered Questions
--------------------

What is our big-picture testing strategy?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We essentially have 2 options:
1. Only make benchmarks compatible with the current Grackle version
2. Generalize the benchmark executable so that it's compatible with older Grackle versions (in practice, we would probably select just a small subset of versions...)

The tradeoffs:
- Choice number 1 is obviously the easier thing to do and maintain
  - in practice, choice number 2 probably isn't *that* bad...
  - the test-harness machinery was explicitly designed to make it straight-forward to support option 2 (in practice, we would just need to replace occurrences of ``GR_INTERNAL_ERRROR`` & ``GR_INTERNAL_REQUIRE``). Any ifdef statements would be very concentrated
  - the build-system logic would get slightly more complex, but it wouldn't be terrible (it would probably be a little simpler the than the logic that makes it possible to compile ``gracklepy`` with an embedded Grackle-build or link against an external Grackle library)
- Choice number 2 is a lot more "forgiving" when it comes to benchmark-design
  - designing robust benchmarks for a library like Grackle is a little tricky ()
