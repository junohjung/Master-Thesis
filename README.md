Master-Thesis
=============
Development of a fully coupled zonal RANS/LES method for the simulation of turbulent backward-facing step flow


ZFS (Zonal Flow Solver)
=======================

ZFS is a multi-physics PDE solver framework with a focus on fluid dynamics
equations. It is developed at the Institute of Aerodynamics of the RWTH Aachen
University.

Zonal Flow Solver Tutorial for developers
=========================================

Quickstart
----------

1. Get ZFS from repository: `svn co http://svn.aia.rwth-aachen.de/zfs/trunk`
2. Run `./configure.py ? ?` in the source directory
3. Run `make` to build ZFS
4. Create property and geometry configuration files, and translate them using
   `ncmpigen -o geometry geometry.cdl` and `ncmpigen -o properties
   properties.cdl`
5. Execute `./zfs` for serial or `mpiexec -n NP ./zfs` for parallel execution


Getting ZFS
-----------

The best way to get ZFS is to check out the current trunk from the Subversion
repository. The trunk is considered to be *release*-grade, thus the latest
version should always be the greatest.

In order to access the repository, you need to be within the RWTH Aachen
University network (either physically or using VPN). Then, you can check out the
trunk by running

    svn co http://svn.aia.rwth-aachen.de/zfs/trunk

Note for Windows users: The Visual Studio project files have to be obtained
separately. After checking out a version of ZFS from the Subversion repository,
enter the working copy and get the project files by executing

    svn export http://svn.aia.rwth-aachen.de/zfs/windows/VS


Install
-------

### Prerequisites

There are a number of hosts which are natively supported by ZFS, so you can skip
this section if you want to use ZFS on any of them:

* AIA (Institute of Aerodynamics)
* Hazel Hen (HLRS Stuttgart)
* JUQUEEN (Forschungszentrum Juelich)
* JURECA (Forschungszentrum Juelich)
* Nehalem (HLRS Stuttgart)
* RWTH Cluster (RWTH Aachen University)

If you would like to use ZFS on another system, these are the minimum
requirements for building ZFS:

* FFTW >= 3.3.2
* Parallel NetCDF >= 1.6.0
* OpenMPI >= 1.8.7

The following components are optional and only need to be installed when using a
part of ZFS that requires them:

* zlib >= 1.2.8 (only needed when running `configure.py` with the `--with-zlib`
  flag)
* doxygen >= 1.6 (only needed when running `make doc` to auto-generate HTML
  documentation)

Once you have satisfied all dependencies, you also need to create a custom host
file. Copy the template file `aux/hosts/Host.cmake.in` to a place of your choice
and uncomment/modify all necessary lines to accommodate your host-specific
settings. Only the fields `HOST_SUPPORTED_COMPILERS`, `HOST_DEFAULT_COMPILER`,
and `HOST_COMPILER_EXE_compilername` are mandatory, but you usually also want
to set the `INCLUDE_DIRS`, `LIBRARY_DIRS`, and `LIBRARY_NAMES` properties as
well. Once you created your custom host file, make sure that CMake can find it
by setting the environment variable `ZFS_HOST_FILE` to the full path of the
host file, e.g. by putting

    export ZFS_HOST_FILE=/home/walle/path/to/Host.cmake

in your `.bashrc` file. If you are not using bash, you need to use your
shell-specific method of setting environment variables.

### Build

Building ZFS requires Python 2.x and CMake >= 2.8, so make sure that both are
on your `PATH` before proceeding. Create and enter an empty build directory (or
go to the ZFS source directory) and run

    path/to/zfs/configure.py [compiler] [build_type]
    make

where `compiler` is the selected compiler and `build_type` the set of
optimization flags to use. If unsure which compilers/build types are supported
on the current host, run `configure.py` without any arguments. If you would like
to just use the host-specific defaults, use `?` instead of the compiler and/or
the build type. To get a full list of options, run `configure.py -h`.

There is no installation step required. You can freely move around the `zfs`
executable on the system, as long as you do not change the paths of the shared
libraries that were used to build it.

The following compilers are currently supported:

* GCC >= 5.1
* Clang >= 3.8
* Intel >= 15.0.1
* PGI >= 15.10

Older versions might work, but there is no guarantee.


Getting started
---------------

The easiest way to get started is by downloading and running a test case from
the ZFS regression tests. You can get the full set of testcases (recommended) by
executing

    svn co http://svn.aia.rwth-aachen.de/zfs/testcases

For this to work, you must be on the RWTH Aaachen University network. If you are
on an AIA host, just compile ZFS manually and then run the full test suite by
executing

    cry test -z path/to/zfs/binary

in the `testcases` root directory. If you are not on an AIA host or want to run
a test case manually, go to one of the test case subdirectories (to be found
under `<blockname>/<test_case_name>`) and execute the following commands:

    ncmpigen -o geometry geometry.cdl
    ncmpigen -o properties properties_grid.cdl
    mpiexec -n NP path/to/zfs
    ncmpigen -o properties properties_run.cdl
    mpiexec -n NP path/to/zfs

The `ncmpigen` commands (part of the Parallel netCDF package) generate the
necessary binary files with runtime properties that are read by ZFS. The first
exexcution of ZFS generates the grid (since we used `properties_grid.cdl` as the
property file) and the second starts the solver. The number of MPI tasks to use
`NP` should be set to something reasonable - if in doubt, check the
corresponding `noDomains` property in the corresponding `.cdl` file (note that
you usually can deviate from the `noDomains` setting by a few numbers without
causing errors).

-----------------------------------------------------------------------------

Contact
-------

If you have questions about zonal RANS/LES connection and zonal exchange procedure, 

please get in contact with:

Junoh Jung, M.Sc. Aerospace Engineering     
Mail: junoh.jung@rwth-aachen.de    
Phone: +49 177 7428216    

Pontstr.135
52062 Aachen  
Germany  
