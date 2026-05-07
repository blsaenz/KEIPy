
# set compile environment
# ----------------------------------------------------------------------
FC = "gfortran"

# Interpreter that owns NumPy (and thus ``python -m numpy.f2py``). Example with env ``py313``:
#   mamba activate py313
#   make kei PYTHON=$(which python)
# NumPy 2 f2py on Python >=3.12 uses Meson; use conda-forge only (override defaults).
# Fortran + LAPACK (SGTSV, etc.) + extension link: conda-forge ``openblas`` is preferred because
# ``-framework Accelerate`` commonly misses symbols f2py’s link step expects—see repo ``environment.yml``.
# One-liner (after ``mamba activate <env>``):
#   mamba install --override-channels -c conda-forge compilers liblapack meson ninja charset-normalizer
# Pip cannot install the conda ``compilers`` or ``openblas`` packages; declare them via ``environment.yml``.
PYTHON ?= python3

# Meson/Ninja live next to PYTHON in conda envs; f2py spawns ``meson`` by name. Adding this avoids
# FileNotFoundError when ``make`` runs without shell activation (still pass full path via PYTHON).
_PYTHON_EXEC := $(shell $(PYTHON) -c "import sys; print(sys.executable)")
PYTHON_BINDIR := $(dir $(_PYTHON_EXEC))
export PATH := $(PYTHON_BINDIR):$(PATH)

# Conda env root: prefer CONDA_PREFIX; else infer from PYTHON if …/<env>/conda-meta exists (no ``activate``).
CONDA_ENV_PREFIX :=
ifneq ($(strip $(CONDA_PREFIX)),)
CONDA_ENV_PREFIX := $(strip $(CONDA_PREFIX))
else
CONDA_ENV_PREFIX := $(shell $(PYTHON) -c 'import pathlib,sys; p=pathlib.Path(sys.executable).resolve(); m=p.parent.parent/"conda-meta"; print(str(m.parent) if m.is_dir() else "")')
endif
HAVE_CONDA_OPENBLAS :=
ifneq ($(strip $(CONDA_ENV_PREFIX)),)
ifneq ($(wildcard $(CONDA_ENV_PREFIX)/lib/libopenblas*),)
HAVE_CONDA_OPENBLAS := 1
endif
endif

# Fortran ``-c`` builds do not need BLAS libraries; LAPACK_MK is legacy (kept for non-conda only).
LAPACK_MK :=
LAPACK_F90_PY :=
F2PY_EXT_LIBS :=
F2PY_LINK_LAPACK := --link-lapack_opt

# determine compile options from environment set above
# ----------------------------------------------------------------------

ifeq ($(FC),ifort)
	COMP_SWITCH := -Difort
	#BUG := -g -O0 -traceback -fpe:0 -debug all -C -fp-model source -assume minus0
	BUG := -g -O2 -traceback -fpe:0 -debug all -xSSE4.2 -ip -fp-model source -align -assume minus0

	MKL_INC = /opt/intel/Compiler/11.1/073/mkl/include
	MKL_LIB = /opt/intel/Compiler/11.1/073/mkl/lib/em64t
  LAPACK = -L$(MKL_LIB) -I$(MKL_INC) -lmkl_lapack95_lp64 $(MKL_LIB)/libmkl_intel_lp64.a $(MKL_LIB)/libmkl_sequential.a $(MKL_LIB)/libmkl_core.a $(MKL_LIB)/libmkl_sequential.a $(MKL_LIB)/libmkl_core.a $(MKL_LIB)/libmkl_sequential.a $(MKL_LIB)/libmkl_core.a -lpthread
	LAPACK_F90_PY := $(LAPACK)
	LAPACK_MK := $(LAPACK)

endif
ifeq ($(FC),"gfortran")

  # mac debug:
  # lldb <executable>
  # >> process launch -i <file>
  # where <file> is std/carrot pipe in

	COMP_SWITCH :=
	BUG := -gdwarf-2 -O0 -fbounds-check -ffpe-trap=invalid,overflow
	#BUG := -gdwarf-2 -O2 -funsafe-math-optimizations
	#BUG := -gdwarf-2 -O2 -ffpe-trap=invalid,overflow

	# for gfortran on 10.9 macs
	#LAPACK = -L/System/Library/Frameworks/Accelerate.framework/Versions/Current/Frameworks/vecLib.framework/Versions/Current/ -lBLAS -lLAPACK
	LAPACK = -framework Accelerate

	LAPACK_F90_PY := $(LAPACK)
	LAPACK_MK := $(LAPACK)

	#
	MKL_INC = /usr/local/opt/lapack/include
	MKL_LIB = /usr/local/opt/lapack/lib
  #LAPACK = -L$(MKL_LIB) -I$(MKL_INC) -lblas -llapack

  # Ubuntu linux, installed apt packages: libblas-dev, liblapack-dev
  #BLAS_LIB = /usr/lib/aarch64-linux-gnu/blas
  #LAPACK_LIB = /usr/lib/aarch64-linux-gnu/lapack
  #LAPACK = -L$(BLAS_LIB) -L$(LAPACK_LIB) -lblas -llapack

endif

ifeq ($(HAVE_CONDA_OPENBLAS),1)
ifeq ($(FC),"gfortran")
# Pass ``-L``/``-l`` as standalone f2py args so Meson records them as link libraries (NumPy parses them;
# embedding in ``--f90flags`` is compile-only). On macOS conda-forge often has no unversioned ``libopenblas.dylib``,
# but ``liblapack.{dylib,so}`` resolves to OpenBLAS; use ``-llapack`` so ``ld`` succeeds.
ifneq ($(wildcard $(CONDA_ENV_PREFIX)/lib/liblapack.dylib $(CONDA_ENV_PREFIX)/lib/liblapack.so),)
F2PY_EXT_LIBS := -L$(CONDA_ENV_PREFIX)/lib -llapack
else
F2PY_EXT_LIBS := -L$(CONDA_ENV_PREFIX)/lib -lopenblas
endif
F2PY_LINK_LAPACK :=
LAPACK_MK :=
LAPACK_F90_PY :=
endif
endif


BD = build

SIESTA_OBJ := kei_kinds.o siesta_constants.o kei_hacks.o siesta_parameters.o siesta_state.o \
    siesta_types.o siesta_grid.o siesta_flux_heat.o siesta_desalination.o

OBJ1D  := kei_kinds.o siesta_constants.o kei_hacks.o siesta_parameters.o siesta_state.o \
  siesta_types.o siesta_grid.o siesta_flux_heat.o siesta_desalination.o\
  kei_ecocommon.o \
  kei_ice.o kei_eco.o kei_init.o kei_fluxes.o \
  kei_subs1D.o kei_atm.o kei_kpp.o kei_ocn.o \
  kei_sw.o

# f2py must compile these .f90 files to expose ``kei.<module>`` in Python (see ``kei.py`` YAML).
# Omit ``kei_ecocommon.o`` from the link list so f2py compiles ``kei_ecocommon.f90`` once (no duplicate symbols).
OBJ1D_F2PY := $(filter-out kei_ecocommon.o,$(OBJ1D))

F901D  := kei_parameters.f90 kei_common.f90 kei_icecommon.f90 kei_ecocommon.f90 kei_link.f90

MACMODS_OBJ := macmods_mod.o macmods_param_mod.o macmods_kinds_mod.o \
  macmods_util_mod.o macmods_calc.o

# Meson f2py compiles listed .f90 in a temp dir; prebuilt $(OBJ1D) do not expose .mod there.
# GNU Fortran searches -I paths for module files built by ``make kei`` prerequisites in this dir.
F2PY_F90_MOD_INC := -I$(CURDIR)

FFLAGS :=  $(BUG)


# compile directives
# ----------------------------------------------------------------------

# keep ice_common (siesta_routines) and macmods out of f2py linking, since it cannot
# handle variables with the same names across different modules
kei: $(OBJ1D)
# Non-conda: Accelerate inside ``--f90flags`` and numpy ``--link-lapack_opt``. Conda: ``F2PY_EXT_LIBS``.
	$(PYTHON) -m numpy.f2py -c --f90exec="gfortran" --f90flags="$(FFLAGS) $(LAPACK_F90_PY) $(F2PY_F90_MOD_INC)" \
		$(F2PY_EXT_LIBS) $(F2PY_LINK_LAPACK) $(OBJ1D_F2PY) $(F901D) $(MACMODS_OBJ) -m kei

# compile routines
kei_fluxes.o: kei_hacks.o
	$(FC) -c $(FFLAGS) kei_fluxes.f90
kei_kpp.o: kei_subs1D.o kei_common.o kei_parameters.o kei_kinds.o kei_icecommon.o
	$(FC) -c $(FFLAGS) kei_kpp.f90
kei_ocn.o: kei_kpp.o kei_common.o kei_parameters.o kei_kinds.o
	$(FC) -c $(FFLAGS) kei_ocn.f90
kei_subs1D.o:
	$(FC) -c $(FFLAGS) kei_subs1D.f90
kei_atm.o	: kei_atm.f90
	$(FC) -c $(FFLAGS) kei_atm.f90
kei_init.o	: kei_ecocommon.o kei_parameters.o kei_common.o kei_icecommon.o
	$(FC) -c $(FFLAGS) kei_init.f90

# compile modules
kei_link.o	: $(OBJ1D) kei_parameters.o kei_common.o kei_icecommon.o
	$(FC) -c $(FFLAGS) kei_link.f90 $(LAPACK_MK)
kei_ice.o	: kei_ice.f90 kei_hacks.o kei_fluxes.o
	$(FC) -c $(FFLAGS) kei_ice.f90
kei_eco.o	: kei_ecocommon.o kei_parameters.o kei_common.o kei_hacks.o
	$(FC) -c $(FFLAGS) kei_eco.f90
kei_ecocommon.o	:
	$(FC) -c $(FFLAGS) kei_ecocommon.f90
kei_icecommon.o	: kei_icecommon.f90 kei_kinds.o
	$(FC) -c $(FFLAGS) kei_icecommon.f90
kei_common.o	: kei_icecommon.o
	$(FC) -c $(FFLAGS) kei_common.f90
kei_parameters.o	: kei_kinds.o
	$(FC) -c $(FFLAGS) kei_parameters.f90

siesta_flux_heat.o	: siesta_constants.o siesta_parameters.o siesta_types.o
	$(FC) -c $(FFLAGS) siesta_flux_heat.f90
siesta_grid.o	: siesta_constants.o siesta_parameters.o siesta_types.o siesta_desalination.o
	$(FC) -c $(FFLAGS) siesta_grid.f90
siesta_desalination.o	: siesta_constants.o siesta_parameters.o siesta_state.o siesta_types.o
	$(FC) -c $(FFLAGS) siesta_desalination.f90
siesta_state.o	: siesta_parameters.o siesta_types.o
	$(FC) -c $(FFLAGS) siesta_state.f90
siesta_types.o	: siesta_constants.o siesta_parameters.o
	$(FC) -c $(FFLAGS) siesta_types.f90
siesta_parameters.o	: siesta_constants.o kei_icecommon.o
	$(FC) -c $(FFLAGS) siesta_parameters.f90
kei_hacks.o	: kei_parameters.o kei_common.o kei_subs1D.o kei_ecocommon.o
	$(FC) -c $(FFLAGS) kei_hacks.f90
kei_kinds.o	: kei_kinds.f90
	$(FC) -c $(FFLAGS) kei_kinds.f90
siesta_constants.o	: kei_kinds.o
	$(FC) -c $(FFLAGS) siesta_constants.f90

kei_sw.o	: $(MACMODS_OBJ)
	$(FC) -c $(FFLAGS) kei_sw.f90
macmods_mod.o    : macmods_param_mod.o macmods_kinds_mod.o macmods_util_mod.o macmods_calc.o
	$(FC) -c $(FFLAGS) macmods_mod.f90
macmods_calc.o    : macmods_param_mod.o macmods_kinds_mod.o macmods_util_mod.o
	$(FC) -c $(FFLAGS) macmods_calc.f90
macmods_util_mod.o    : macmods_param_mod.o macmods_kinds_mod.o
	$(FC) -c $(FFLAGS) macmods_util_mod.f90
macmods_param_mod.o    : macmods_kinds_mod.o
	$(FC) -c $(FFLAGS) macmods_param_mod.f90
macmods_kinds_mod.o    : kei_kinds.o
	$(FC) -c $(FFLAGS) macmods_kinds_mod.f90

#	-------
#	CLEANUP
#	-------
clean :
	rm -f *.o
	rm -f *.mod
	rm -f *.so
	rm -f -r KEI.run.dSYM

# end
