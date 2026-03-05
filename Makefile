#======================================================================
# FDPACK - Finite Difference Package
#======================================================================

# compiler selection
FC = ifort
#FC = ifx
#FC = gfortran

# build mode
MODE = debug
#MODE = release

# LAPACK support (set to 1 to
# enable diag_num_hess)
USE_LAPACK = 1

# LAPACK/BLAS libraries
LAPACK_PATH = ~/.local/lib/lapack-3.12.0/lib/$(TC)
LAPACK = $(LAPACK_PATH)/liblapack.a
BLAS = $(LAPACK_PATH)/librefblas.a

# directory structure
SRC   = src
BIN   = bin
BUILD = build
DIST  = dist

#======================================================================
# compiler flags
#======================================================================

ifeq ($(FC),gfortran)
  TC     = gnu
  OFLAGS = -O2
  DFLAGS = -O0 -g -fcheck=all
  GFLAGS = -cpp -fopenmp -J$(BUILD) -Wunused-variable \
           -fbacktrace -ffpe-summary=none -fno-stack-arrays -I.
else ifneq ($(filter ifort ifx,$(FC)),)
  TC = intel
  ifeq ($(FC),ifort)
    CHFLAG  = -check all
    DIAFLAG = -diag-disable=10448
  else ifeq ($(FC),ifx)
    CHFLAG  = -check all,nouninit
    DIAFLAG =
  endif
  OFLAGS = -O3 $(DIAFLAG) -fp-model strict
  DFLAGS = -O0 -traceback $(DIAFLAG) $(CHFLAG) -debug -warn all
  GFLAGS = -cpp -module $(BUILD)
else
  $(error Unsupported Fortran compiler: FC=$(FC))
endif

# linker flags
LDFLAGS = -Wl,-z,execstack,--no-warn-execstack

# combined flags for selected mode
ifeq ($(MODE),release)
  FFLAGS = $(OFLAGS) $(GFLAGS)
else ifeq ($(MODE),debug)
  FFLAGS = $(DFLAGS) $(GFLAGS)
else
  $(error Unsupported MODE: MODE=$(MODE))
endif

# enable diag_num_hess if USE_LAPACK is set
ifeq ($(USE_LAPACK),1)
  FFLAGS += -DUSE_LAPACK
  ifeq ($(wildcard $(LAPACK)),)
    $(error USE_LAPACK=1 but LAPACK not found: $(LAPACK))
  endif
  ifeq ($(wildcard $(BLAS)),)
    $(error USE_LAPACK=1 but BLAS not found: $(BLAS))
  endif
endif

#======================================================================
# targets
#======================================================================

# find test_*.f90 sources; derive executable names
TSRC = $(shell find $(SRC) -name 'test_*.f90')
TARG = $(addprefix $(BIN)/, $(notdir $(TSRC:.f90=.x)))

# 'all' must be the first target
all: $(TARG) dist

#======================================================================
# object generation system
#======================================================================

include mk/obgen.mk

#----------------------------------------------------------------------
# utility modules
#----------------------------------------------------------------------

# error handling
SFILE = mod_catch.f90
DEP   = NONE
CMD   = $(FC) $(FFLAGS) -c $$< -o $$@
$(eval $(call obgen, $(SFILE), $(DEP), $(CMD)))

# swap utilities for bracket management
SFILE = mod_swap.f90
DEP   = NONE
CMD   = $(FC) $(FFLAGS) -c $$< -o $$@
$(eval $(call obgen, $(SFILE), $(DEP), $(CMD)))

#----------------------------------------------------------------------
# core modules
#----------------------------------------------------------------------

# finite difference routines (grad, hess)
SFILE = mod_fdbase.f90
DEP   = mod_catch.o mod_swap.o
CMD   = $(FC) $(FFLAGS) -c $$< -o $$@
$(eval $(call obgen, $(SFILE), $(DEP), $(CMD)))

# offset control and noise estimation
SFILE = mod_fdtune.f90
DEP   = mod_catch.o mod_swap.o mod_fdbase.o
CMD   = $(FC) $(FFLAGS) -c $$< -o $$@
$(eval $(call obgen, $(SFILE), $(DEP), $(CMD)))

#----------------------------------------------------------------------
# test programs
#----------------------------------------------------------------------

SFILE = test_demo.f90
DEP   = mod_catch.o mod_swap.o mod_fdbase.o mod_fdtune.o
CMD   = $(FC) $(FFLAGS) -c $$< -o $$@
$(eval $(call obgen, $(SFILE), $(DEP), $(CMD)))

SFILE = test_func_noise.f90
DEP   = mod_catch.o mod_swap.o mod_fdbase.o mod_fdtune.o
CMD   = $(FC) $(FFLAGS) -c $$< -o $$@
$(eval $(call obgen, $(SFILE), $(DEP), $(CMD)))

SFILE = test_grad_offsets.f90
DEP   = mod_catch.o mod_swap.o mod_fdbase.o mod_fdtune.o
CMD   = $(FC) $(FFLAGS) -c $$< -o $$@
$(eval $(call obgen, $(SFILE), $(DEP), $(CMD)))

SFILE = test_hess_offsets.f90
DEP   = mod_catch.o mod_swap.o mod_fdbase.o mod_fdtune.o
CMD   = $(FC) $(FFLAGS) -c $$< -o $$@
$(eval $(call obgen, $(SFILE), $(DEP), $(CMD)))

#======================================================================
# build rules
#======================================================================

# module objects only (no test mains)
OBJS  = $(filter-out $(BUILD)/test_%.o,\
        $(addprefix $(BUILD)/,\
        $(notdir $(SFILES:.f90=.o))))

# test objects
TOBJS = $(addprefix $(BUILD)/, $(notdir $(TSRC:.f90=.o)))

# stamp: all objects compiled before linking begins
$(BUILD)/.cmp_stamp: $(OBJS) $(TOBJS)
	@touch $@

$(BIN)/%.x: $(BUILD)/.cmp_stamp $(BUILD)/%.o
	@echo "--> building $@ ..."
	@$(FC) $(FFLAGS) $(LDFLAGS) $(OBJS) \
	$(BUILD)/$(notdir $(@:.x=.o)) \
	$(if $(filter 1,$(USE_LAPACK)),$(LAPACK) $(BLAS),) -o $@

#======================================================================
# distribution
#======================================================================

# static library containing only the public modules
$(BUILD)/libfdpack.a: $(BUILD)/mod_fdbase.o $(BUILD)/mod_fdtune.o
	@echo "--> building $@ ..."
	@ar rcs $@ $^

$(BUILD)/.dist_stamp: $(BUILD)/libfdpack.a
	@echo "--> populating $(DIST)/$(TC) ..."
	@cp $(BUILD)/libfdpack.a $(DIST)/$(TC)
	@cp $(BUILD)/mod_fdbase.mod $(DIST)/$(TC)
	@cp $(BUILD)/mod_fdtune.mod $(DIST)/$(TC)
	@touch $@

dist: $(BUILD)/.dist_stamp

#======================================================================
# run tests
#======================================================================

test: $(TARG)
	@echo "Available tests:"; \
	select x in \
		$(notdir $(TARG)) \
		"test_hess_offsets approx" \
		all quit; do \
		case $$x in \
			"test_hess_offsets approx") \
				$(BIN)/test_hess_offsets.x approx; \
				break;; \
			all) \
				for t in $(TARG); do \
					echo "--- $$t ---"; \
					$$t; \
				done; \
				echo "--- test_hess_offsets" \
				"approx ---"; \
				$(BIN)/test_hess_offsets.x approx; \
				break;; \
			quit) \
				break;; \
			'') \
				echo "Invalid choice";; \
			*) \
				$(BIN)/$$x; \
				break;; \
		esac; \
	done

#======================================================================
# cleanup
#======================================================================

.PHONY: all test dist clean

clean:
	rm -f $(BUILD)/* $(BIN)/* $(DIST)/$(TC)/*
