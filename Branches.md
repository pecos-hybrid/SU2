Branches
=======

Active branches
---------------

+ _pecos-dev_ The "develop" branch for the pecos-hybrid group.  This is
  meant to be a stable branch, with relatively no unfinished work and as
  few bugs as possible.
+ _pecos-improved-RANS-averaging_ The new implementation of averaging, using
  a whole CVariable vector for averages instead of separate "average
  solution" fields.
+ _refactor-viscous-numerics_ Changes the viscous numerics to allow more
  general stress tensors and heat flux Jacobians, with the focus being
  on the model-split hybridization.
+ _pecos-rk-fixes_ Various improvements in the Runge-Kutta methods.
   This branch was working, but other problems were noticed with the RK
   time-stepping.  These problems were believed to be due to the
   approximations made in calculating Jacobians.  Focus has shifted away
   from RK methods due to these problems, and this branch was abandoned.
+ _pecos-model-split_ The main development branch for the fully integrated
  model-split hybridization.  While other branches may contain pieces of
  the model-split hybridization, this branch is used to combine all the
  pieces together.
+ _pecos-new-forcing_ The branch containing the stream-function version
  of the model-split forcing.  This branch has not been tested, and may
  have some issues.

Branches based off su2code/SU2
------------------------------

+ _fix-cpp-inlet-with-periodic_ Fixes some issues when using the c++
  version of inlet profile specification, especially with periodic
  boundaries.
+ _fix_unify-FD-NTS-blending_ This branch contains some improvements and
  bug fixes for the FD/NTS central/upwind blending.

Deprecated branches
-------------------

+ _model-split_ The "original" model-split branch.  This branch has an early
  implementation of the forcing and a 4/3 precomputed power of the resolution
  tensor. Last updated May 16, 2018.
+ _pecos-hybrid-forcing_ Contains an early version of the forcing. Last
  updated Feb 5, 2018
+ _pecos-old-forcing_ Contains an early version of the forcing. Last
  updated on May 16, 2018.
+ _pecos-zetaf_
+ _pecos-zetaf_freestream_debug_
