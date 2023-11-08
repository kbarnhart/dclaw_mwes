"""
Module to set up run time parameters for Clawpack.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.

"""
import os

import numpy as np
import yaml
from pyclaw import data

with open("params.yaml") as file:
    params = yaml.safe_load(file)

tfinal = 30
nout = 30

xlower = params["xlower"]
ylower = params["ylower"]
yupper = params["yupper"]
xupper = params["xupper"]

mx = int(xupper-xlower)
my = int(yupper-ylower)

m_crit = 0.64
m0=0.62
kappita = 1e-10
g=9.81
mannings=0.025
alpha_seg=1
chi_init_val=0.5
kappita_diff = 0.5
topo_file = params["topo_file"]
q1_file = params["q1_file"]
aux5_file = params["aux5_file"]


# ------------------------------
def setrun(claw_pkg="digclaw"):
    # ------------------------------

    """
    Define the parameters used for running Clawpack.

    INPUT:
        claw_pkg expected to be "geoclaw" for this setrun.

    OUTPUT:
        rundata - object of class ClawRunData

    """

    # assert claw_pkg.lower() == 'digclaw',  "Expected claw_pkg = 'digclaw'"
    ndim = 2
    rundata = data.ClawRunData(claw_pkg, ndim)

    # ------------------------------------------------------------------
    # GeoClaw specific parameters:
    # ------------------------------------------------------------------

    rundata = setgeo(rundata)  # Defined below

    # ------------------------------------------------------------------
    # DigClaw specific parameters:
    # ------------------------------------------------------------------

    rundata = setdig(rundata)  # Defined below

    # ------------------------------------------------------------------
    # Standard Clawpack parameters to be written to claw.data:
    #   (or to amr2ez.data for AMR)
    # ------------------------------------------------------------------

    clawdata = rundata.clawdata  # initialized when rundata instantiated

    # Set single grid parameters first.
    # See below for AMR parameters.

    # ---------------
    # Spatial domain:
    # ---------------

    # Number of space dimensions:
    clawdata.ndim = ndim

    # Lower and upper edge of computational domain:

    clawdata.xlower = xlower
    clawdata.xupper = xupper

    clawdata.ylower = ylower
    clawdata.yupper = yupper

    # Number of grid cells:
    clawdata.mx = mx
    clawdata.my = my

    # ---------------
    # Size of system:
    # ---------------

    # Number of equations in the system:
    clawdata.meqn = 7

    # Number of auxiliary variables in the aux array (initialized in setaux)
    clawdata.maux = 10

    # Index of aux array corresponding to capacity function, if there is one:
    clawdata.mcapa = 0

    # -------------
    # Initial time:
    # -------------

    clawdata.t0 = 0.0

    # -------------
    # Output times:
    # --------------

    # Specify at what times the results should be written to fort.q files.
    # Note that the time integration stops after the final output time.
    # The solution at initial time t0 is always written in addition.

    clawdata.outstyle = 1

    if clawdata.outstyle == 1:
        # Output nout frames at equally spaced times up to tfinal:
        clawdata.nout = nout  # HERE
        clawdata.tfinal = tfinal  # HERE

    # ---------------------------------------------------
    # Verbosity of messages to screen during integration:
    # ---------------------------------------------------

    # The current t, dt, and cfl will be printed every time step
    # at AMR levels <= verbosity.  Set verbosity = 0 for no printing.
    #   (E.g. verbosity == 2 means print only on levels 1 and 2.)
    clawdata.verbosity = 0

    # --------------
    # Time stepping:
    # --------------

    # if dt_variable==1: variable time steps used based on cfl_desired,
    # if dt_variable==0: fixed time steps dt = dt_initial will always be used.
    clawdata.dt_variable = 1

    # Initial time step for variable dt.
    # If dt_variable==0 then dt=dt_initial for all steps:
    clawdata.dt_initial = 0.01  # HERE

    # Max time step to be allowed if variable dt used:
    clawdata.dt_max = 1.0  # HERE

    # Desired Courant number if variable dt used, and max to allow without
    # retaking step with a smaller dt:
    clawdata.cfl_desired = 0.25
    clawdata.cfl_max = 0.8

    # Maximum number of time steps to allow between output times:
    clawdata.max_steps = 2147483647  # likely largest integer

    # ------------------
    # Method to be used:
    # ------------------

    # Order of accuracy:  1 => Godunov,  2 => Lax-Wendroff plus limiters
    clawdata.order = 1

    # Transverse order for 2d or 3d (not used in 1d):
    clawdata.order_trans = 0

    # Number of waves in the Riemann solution:
    clawdata.mwaves = 5

    # List of limiters to use for each wave family:
    # Required:  len(mthlim) == mwaves
    clawdata.mthlim = [4, 4, 4, 4, 4]

    # Source terms splitting:
    #   src_split == 0  => no source term (src routine never called)
    #   src_split == 1  => Godunov (1st order) splitting used,
    #   src_split == 2  => Strang (2nd order) splitting used,  not recommended.
    clawdata.src_split = 1

    # --------------------
    # Boundary conditions:
    # --------------------

    # Number of ghost cells (usually 2)
    clawdata.mbc = 2

    # Choice of BCs at xlower and xupper:
    #   0 => user specified (must modify bcN.f to use this option)
    #   1 => extrapolation (non-reflecting outflow)
    #   2 => periodic (must specify this at both boundaries)
    #   3 => solid wall for systems where q(2) is normal velocity

    clawdata.mthbc_xlower = 1
    clawdata.mthbc_xupper = 1

    clawdata.mthbc_ylower = 1
    clawdata.mthbc_yupper = 1

    # ---------------
    # AMR parameters:
    # ---------------

    # max number of refinement levels:

    clawdata.mxnest = -1  # negative ==> anisotropic refinement in x,y,t

    # List of refinement ratios at each level (length at least mxnest-1)
    clawdata.inratx = [5]
    clawdata.inraty = [5]
    clawdata.inratt = [5]

    # Specify type of each aux variable in clawdata.auxtype.
    # This must be a list of length maux, each element of which is one of:
    #   'center',  'capacity', 'xleft', or 'yleft'  (see documentation).

    clawdata.auxtype = [
        "center",
        "center",
        "yleft",
        "center",
        "center",
        "center",
        "center",
        "center",
        "center",
        "center",
    ]

    clawdata.tol = -1.0  # negative ==> don't use Richardson estimator
    clawdata.tolsp = 0.5  # used in default flag2refine subroutine
    # (Not used in geoclaw!)

    clawdata.kcheck = 2  # how often to regrid (every kcheck steps)
    clawdata.ibuff = 4  # width of buffer zone around flagged points
    clawdata.cutoff = 0.7  # efficiency cutoff for grid generator
    clawdata.checkpt_iousr = 2147483647
    clawdata.restart = False
    # More AMR parameters can be set -- see the defaults in pyclaw/data.py

    return rundata
    # end of function setrun
    # ----------------------


def setgeo(rundata):
    """
    Set GeoClaw specific runtime parameters.
    For documentation see ....
    """

    try:
        geodata = rundata.geodata
    except:
        print("*** Error, this rundata has no geodata attribute")
        raise AttributeError("Missing geodata attribute")

    geodata.variable_dt_refinement_ratios = True

    # == setgeo.data values ==
    R1 = 6357.0e3  # polar radius
    R2 = 6378.0e3  # equatorial radius
    Rearth = 0.5 * (R1 + R2)
    geodata.igravity = 1
    geodata.gravity = g
    geodata.icoordsys = 1
    geodata.icoriolis = 0
    geodata.Rearth = Rearth

    # == settsunami.data values ==
    geodata.sealevel = -10000000.0
    geodata.drytolerance = (
        1e-5  # drytol < h in flowgrades (only dryed if derefined naturally)
    )
    geodata.wavetolerance = 5.0e-2
    geodata.depthdeep = 1.0e2
    geodata.maxleveldeep = 1
    geodata.ifriction = 1
    geodata.coeffmanning = mannings
    geodata.frictiondepth = 10000.0

    # == settopo.data values ==
    # set a path variable for the base topo directory for portability
    geodata.topofiles = []
    #   [topotype, minlevel,maxlevel,startTime,endTime,fname]
    geodata.topofiles.append([2, 1, 1, 0.0, 1.0e10, topo_file])  # HERE

    # == setdtopo.data values ==
    # == setdtopo.data values ==
    geodata.dtopofiles = []
    # for moving topography, append lines of the form:
    #   [topotype, minlevel,maxlevel,fname]

    # == setqinit.data values ==
    geodata.qinitfiles = []
    geodata.qinitfiles.append([2, 1, 1, 1, q1_file])
    # for qinit perturbations append lines of the form
    #   [qinitftype,iqinit, minlev, maxlev, fname]

    # qinitftype: file-type, same as topo files, ie: 1, 2 or 3
    # The following values are allowed for iqinit:
    # n=1,meqn perturbation of q(i,j,n)
    # n=meqn+1: surface elevation eta is defined by the file and results in h=max(eta-b,0)

    # == setauxinit.data values ==
    geodata.auxinitfiles = []
    geodata.auxinitfiles.append([2, 5, 1, 1, aux5_file])
    # for auxinit perturbations append lines of the form
    #   [auxinitftype,iauxinit, minlev, maxlev, fname]

    # auxinitftype: file-type, same as topo files, ie: 1, 2 or 3
    # The following values are allowed for iauxinit:
    # n=1,maux perturbation of aux(i,j,n)

    # == setregions.data values ==
    geodata.regions = []

    # == setgauges.data values ==
    geodata.gauges = []
    # for gauges append lines of the form  [gaugeno, x, y, t0, tf]

    # == setfixedgrids.data values ==
    geodata.fixedgrids = []
    # for fixed grids append lines of the form
    # [t1,t2,noutput,x1,x2,y1,y2,xpoints,ypoints,\
    #  ioutarrivaltimes,ioutsurfacemax]
    # geodata.fixedgrids.append([54.e3,55.e3,100,-101.,-96.,14.,19.,1000,1000,0,0])

    # == setflowgrades.data values ==
    geodata.flowgrades = []
    # for using flowgrades for refinement append lines of the form
    # [flowgradevalue, flowgradevariable, flowgradetype, flowgrademinlevel]
    # where:
    # flowgradevalue: floating point relevant flowgrade value for following measure:
    # flowgradevariable: 1=depth, 2= momentum, 3 = sign(depth)*(depth+topo) (0 at sealevel or dry land).
    # flowgradetype: 1 = norm(flowgradevariable), 2 = norm(grad(flowgradevariable))
    # flowgrademinlevel: refine to at least this level if flowgradevalue is exceeded.
    geodata.keep_fine = True
    geodata.flowgrades.append([1.0e-6, 2, 1, 1])
    geodata.flowgrades.append([1.0e-6, 1, 1, 1])

    return rundata


def setdig(rundata):
    """
    Set DigClaw specific runtime parameters.
    For documentation see ....
    """

    try:
        digdata = rundata.digdata
    except:
        print("*** Error, this rundata has no digdata attribute")
        raise AttributeError("Missing digdata attribute")

    # set non-default values if needed
    digdata.c1 = 1.0
    digdata.rho_f = 1100.0
    digdata.rho_s = 2700.0
    digdata.phi_bed = 36.0
    digdata.phi_int = 36.0
    digdata.theta_input = 0.0
    digdata.mu = 0.005
    digdata.m0 = m0
    digdata.m_crit = m_crit
    digdata.kappita = kappita
    digdata.kappita_diff = kappita_diff
    digdata.chi_init_val=chi_init_val
    digdata.alpha_c = 0.05
    digdata.alpha_seg = alpha_seg
    digdata.phi_seg_coeff = 0.0
    digdata.delta = 0.001
    digdata.bed_normal = 0
    digdata.entrainment = 1
    digdata.entrainment_rate = 0.2
    digdata.sigma_0 = 1.0e3
    digdata.mom_autostop = True
    digdata.momlevel = 1
    digdata.mom_perc = 0.0
    digdata.phys_tol = rundata.geodata.drytolerance

    digdata.init_ptype = 0
    digdata.init_pmax_ratio = 0.00e0
    digdata.init_ptf = 0.0
    digdata.init_ptf2 = 0.0

    # -1 =0, 0 = hydro, 1,2 = failure or average failure, 3,4= p(t) to failure or average failure
    # to reduce to shallow water equations, uncomment the following
    # digdata.c1= 0.0
    # digdata.phi_int = 0.0
    # digdata.phi_bed = 0.0
    # digdata.kappita = 0.0
    # digdata.mu = 0.0

    return rundata


if __name__ == "__main__":
    # Set up run-time parameters and write all data files.
    import sys

    if len(sys.argv) == 2:
        rundata = setrun(sys.argv[1])
    else:
        rundata = setrun()

    rundata.write()
