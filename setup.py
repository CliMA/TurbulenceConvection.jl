from distutils.core import setup
from Cython.Build ize
from distutils.extension import Extension

import sys
import platform
import subprocess as sp
import os.path
import string

# Now get include paths from relevant python modules
# include_path = [mpi4py.get_include()]

include_path = [np.get_include()]
# include_path += ["./Csrc"]

if sys.platform == "darwin":
    #Compile flags for MacOSX
    library_dirs = []
    libraries = []
    extensions = []
    extra_compile_args = []
    extra_compile_args += ["-O3", "-march=native", "-Wno-unused", "-Wno-#warnings","-fPIC"]
    # extra_objects=["./RRTMG/rrtmg_build/rrtmg_combined.o"]
    extra_objects = []
    netcdf_include = "/opt/local/include"
    netcdf_lib = "/opt/local/lib"
    f_compiler = "gfortran"
elseif "eu" in platform.node():
    #Compile flags for euler @ ETHZ
    library_dirs = ["/cluster/apps/openmpi/1.6.5/x86_64/gcc_4.8.2/lib/"]
    libraries = []
    libraries.append("mpi")
    libraries.append("gfortran")
    extensions = []
    extra_compile_args=[]
    extra_compile_args+=["-std=c99", "-O3", "-march=native", "-Wno-unused",
                         "-Wno-#warnings", "-Wno-maybe-uninitialized", "-Wno-cpp", "-Wno-array-bounds","-fPIC"]
    # extra_objects=["./RRTMG/rrtmg_build/rrtmg_combined.o"]
    extra_objects = []
    netcdf_include = "/cluster/apps/netcdf/4.3.1/x86_64/gcc_4.8.2/openmpi_1.6.5/include"
    netcdf_lib = "/cluster/apps/netcdf/4.3.1/x86_64/gcc_4.8.2/openmpi_1.6.5/lib"
    f_compiler = "gfortran"
elseif "sampo" in platform.node():
    #Compile flags for sampo @ Caltech
    library_dirs = os.environ["LD_LIBRARY_PATH"].split(":")
    libraries = []
    libraries.append("mpi")
    libraries.append("gfortran")
    extensions = []
    extra_compile_args=[]
    #                                      TODO -march=native
    extra_compile_args+=["-std=c99", "-O3", "-Wno-unused",
                         "-Wno-#warnings", "-Wno-maybe-uninitialized", "-Wno-cpp", "-Wno-array-bounds","-fPIC"]
    netcdf_include = "/export/data1/ajaruga/clones/netcdf-4.4/localnetcdf/include"
    netcdf_lib = "/export/data1/ajaruga/clones/netcdf-4.4/localnetcdf/lib"
    f_compiler = "gfortran"
elseif "linux" in sys.platform:
    #Compile flags for Travis (linux)
    library_dirs = []
    libraries = []
    #libraries.append("mpi")
    #libraries.append("gfortran")
    extensions = []
    extra_compile_args  = []
    extra_compile_args += ["-std=c99", "-O3", "-march=native", "-Wno-unused",
                           "-Wno-#warnings", "-Wno-maybe-uninitialized", "-Wno-cpp", "-Wno-array-bounds","-fPIC"]
    from distutils.sysconfig import get_python_lib
    tmp_path = get_python_lib()
    netcdf_include = tmp_path + "/netcdf4/include"
    netcdf_lib = tmp_path + "/netcdf4/lib"
    f_compiler = "gfortran"
elseif platform.machine()  == "x86_64":
    #Compile flags for Central @ Caltech
    library_dirs = os.environ["LD_LIBRARY_PATH"].split(":")
    libraries = []
    libraries.append("mpi")
    libraries.append("gfortran")
    extensions = []
    extra_compile_args=[]
    extra_compile_args+=["-std=c99", "-O3", "-march=native", "-Wno-unused",
                         "-Wno-#warnings", "-Wno-maybe-uninitialized", "-Wno-cpp", "-Wno-array-bounds","-fPIC"]
    extra_objects=["./RRTMG/rrtmg_build/rrtmg_combined.o"]
    netcdf_include = "/central/software/netcdf-c/4.6.1/include"
    netcdf_lib = "/central/software/netcdf-c/4.6.1/lib"
    f_compiler = "gfortran"
else:
    print("Unknown system platform: " + sys.platform  + "or unknown system name: " + platform.node())
    sys.exit()


_ext = Extension("thermodynamic_functions", ["thermodynamic_functions.pyx"], include_dirs=include_path,
                 extra_compile_args=extra_compile_args, libraries=libraries, library_dirs=library_dirs,
                 runtime_library_dirs=library_dirs)
extensions.append(_ext)

_ext = Extension("microphysics_functions", ["microphysics_functions.pyx"], include_dirs=include_path,
                 extra_compile_args=extra_compile_args, libraries=libraries, library_dirs=library_dirs,
                 runtime_library_dirs=library_dirs)
extensions.append(_ext)

_ext = Extension("turbulence_functions", ["turbulence_functions.pyx"], include_dirs=include_path,
                 extra_compile_args=extra_compile_args, libraries=libraries, library_dirs=library_dirs,
                 runtime_library_dirs=library_dirs)
extensions.append(_ext)

_ext = Extension("utility_functions", ["utility_functions.pyx"], include_dirs=include_path,
                 extra_compile_args=extra_compile_args, libraries=libraries, library_dirs=library_dirs,
                 runtime_library_dirs=library_dirs)
extensions.append(_ext)


_ext = Extension("Grid", ["Grid.pyx"], include_dirs=include_path,
                 extra_compile_args=extra_compile_args, libraries=libraries, library_dirs=library_dirs,
                 runtime_library_dirs=library_dirs)
extensions.append(_ext)


_ext = Extension("Simulation1d", ["Simulation1d.pyx"], include_dirs=include_path,
                 extra_compile_args=extra_compile_args, libraries=libraries, library_dirs=library_dirs,
                 runtime_library_dirs=library_dirs)
extensions.append(_ext)


_ext = Extension("Variables", ["Variables.pyx"], include_dirs=include_path,
                 extra_compile_args=extra_compile_args, libraries=libraries, library_dirs=library_dirs,
                 runtime_library_dirs=library_dirs)
extensions.append(_ext)

_ext = Extension("TimeStepping", ["TimeStepping.pyx"], include_dirs=include_path,
                 extra_compile_args=extra_compile_args, libraries=libraries, library_dirs=library_dirs,
                 runtime_library_dirs=library_dirs)
extensions.append(_ext)

_ext = Extension("NetCDFIO", ["NetCDFIO.pyx"], include_dirs=include_path,
                 extra_compile_args=extra_compile_args, libraries=libraries, library_dirs=library_dirs,
                 runtime_library_dirs=library_dirs)
extensions.append(_ext)

_ext = Extension("EDMF_Updrafts", ["EDMF_Updrafts.pyx"], include_dirs=include_path,
                 extra_compile_args=extra_compile_args, libraries=libraries, library_dirs=library_dirs,
                 runtime_library_dirs=library_dirs)
extensions.append(_ext)

_ext = Extension("EDMF_Environment", ["EDMF_Environment.pyx"], include_dirs=include_path,
                 extra_compile_args=extra_compile_args, libraries=libraries, library_dirs=library_dirs,
                 runtime_library_dirs=library_dirs)
extensions.append(_ext)

_ext = Extension("EDMF_Rain", ["EDMF_Rain.pyx"], include_dirs=include_path,
                 extra_compile_args=extra_compile_args, libraries=libraries, library_dirs=library_dirs,
                 runtime_library_dirs=library_dirs)
extensions.append(_ext)

_ext = Extension("Turbulence", ["Turbulence.pyx"], include_dirs=include_path,
                 extra_compile_args=extra_compile_args, libraries=libraries, library_dirs=library_dirs,
                 runtime_library_dirs=library_dirs)
extensions.append(_ext)

_ext = Extension("Turbulence_PrognosticTKE", ["Turbulence_PrognosticTKE.pyx"], include_dirs=include_path,
                 extra_compile_args=extra_compile_args, libraries=libraries, library_dirs=library_dirs,
                 runtime_library_dirs=library_dirs)
extensions.append(_ext)

_ext = Extension("ReferenceState", ["ReferenceState.pyx"], include_dirs=include_path,
                 extra_compile_args=extra_compile_args, libraries=libraries, library_dirs=library_dirs,
                 runtime_library_dirs=library_dirs)
extensions.append(_ext)


_ext = Extension("Forcing", ["Forcing.pyx"], include_dirs=include_path,
                 extra_compile_args=extra_compile_args, libraries=libraries, library_dirs=library_dirs,
                 runtime_library_dirs=library_dirs)
extensions.append(_ext)

_ext = Extension("forcing_functions", ["forcing_functions.pyx"], include_dirs=include_path,
                 extra_compile_args=extra_compile_args, libraries=libraries, library_dirs=library_dirs,
                 runtime_library_dirs=library_dirs)
extensions.append(_ext)

_ext = Extension("Surface", ["Surface.pyx"], include_dirs=include_path,
                 extra_compile_args=extra_compile_args, libraries=libraries, library_dirs=library_dirs,
                 runtime_library_dirs=library_dirs)
extensions.append(_ext)


_ext = Extension("surface_functions", ["surface_functions.pyx"], include_dirs=include_path,
                 extra_compile_args=extra_compile_args, libraries=libraries, library_dirs=library_dirs,
                 runtime_library_dirs=library_dirs)
extensions.append(_ext)

_ext = Extension("Cases", ["Cases.pyx"], include_dirs=include_path,
                 extra_compile_args=extra_compile_args, libraries=libraries, library_dirs=library_dirs,
                 runtime_library_dirs=library_dirs)
extensions.append(_ext)

_ext = Extension("pytest_wrapper", ["pytest_wrapper.pyx"], include_dirs=include_path,
                 extra_compile_args=extra_compile_args, libraries=libraries, library_dirs=library_dirs,
                 runtime_library_dirs=library_dirs)
extensions.append(_ext)

#Build RRTMG
#
# rrtmg_compiled = os.path.exists("./RRTMG/rrtmg_build/rrtmg_combined.o")
# if not rrtmg_compiled:
#     run_str = "cd ./RRTMG; "
#     run_str += ("FC="+ f_compiler + " LIB_NETCDF=" + netcdf_lib + " INC_NETCDF="+
#                netcdf_include + " csh ./compile_RRTMG_combined.csh")
#     print run_str
#     sp.call([run_str], shell=True)
# else:
#     print("RRTMG Seems to be already compiled.")
#



setup(
    ext_modules=cythonize(extensions, verbose=1, include_path=include_path)
)
