import numpy as np
from pathlib import Path
import dclaw.topotools as gt
import yaml 

# create topo and q1
V = 1000
alpha1 = np.deg2rad(1) # fan angle
alpha2 = np.deg2rad(7) # channel angle
alpha3 = np.deg2rad(30) # hillslope angle

# Calculation of the wedge
tanalpha2 = np.tan(alpha2)
tanalpha3 = np.tan(alpha3)
tantheta = tanalpha2/tanalpha3

prefactor = 3/(tanalpha2 * tantheta)

L = float((prefactor*V)**(1/3))
H = float(L * tanalpha2)
W = float(L * tantheta)

aspect = 3
half_width = np.ceil(6*W)
half_length = aspect * half_width

def topo_src(X,Y):
    """
    make topo
    """
    Z = np.zeros(X.shape)

    # where Z>0 slope down at alpha1
    # where Z<0, slope up at alpha2
    # Where Z<0 make two side slopes at alpha3

    fan = X>=0
    upslope = X<0

    Z[fan]= -np.tan(alpha1) * X[fan]
    Z[upslope]= np.tan(alpha2) * np.abs(X[upslope]) + np.tan(alpha3) * np.abs(Y[upslope])

    return Z

def a_src(X,Y):
    """
    make entrainable
    """
    a = np.zeros(X.shape)

    upslope = X<0

    a[upslope]= 5
    return a

def q1_src(X,Y):
    """
    make q1
    """
    Z=topo_src(X,Y)

    # downstream point of dambreak
    xmax = -20*H
    z_at_xmax = -xmax*tanalpha2

    # top of dam break
    z_top_of_flow = z_at_xmax + H
    dx = X[0,1]-X[0,0]

    volume_adjusting = True
    while volume_adjusting:
        q1 = np.zeros(Z.shape)
        
        # set inundated cells to difference between 
        # flat plane and topography
        sel = (X<xmax)&(Z<z_top_of_flow)
        assert np.sum(sel)>0

        q1[sel]= z_top_of_flow - Z[sel]
        assert np.all(q1>=0)

        # to get exact value for volume, we will need
        # to adjust a little (L was rounded down to an integer)
        missing_V = V-np.sum(q1)*dx*dx
        if missing_V>0:
            ncells = np.sum(sel)
            missing_depth = missing_V/(dx*dx*ncells)
            q1[sel]+= missing_depth
            volume_adjusting = False
        else:
            # lower z top of flow by the equivalent of 
            # one horizontal grid cell at the slope of 
            # alpha2 or alpha3, whichever is smaller
            # dz = dx tan alpha2
            # 
            adjustment = dx*min(tanalpha3, tanalpha2)
            z_top_of_flow -= adjustment    
    print(V-np.sum(q1)*dx*dx)
    return q1

dx = 1
xlower = -half_length
xupper = half_length
ylower = -half_width
yupper = half_width

print(xlower, xupper, ylower, yupper)
nxpoints = int((xupper-xlower)/dx) + 1
nypoints = int((yupper-ylower)/dx) + 1


outfile= 'topo.tt2'
gt.topo2writer(outfile,topo_src,xlower,xupper,ylower,yupper,nxpoints,nypoints)

outfile= 'q1.tt2'
gt.topo2writer(outfile,q1_src,xlower,xupper,ylower,yupper,nxpoints,nypoints)

outfile= 'aux5.tt2'
gt.topo2writer(outfile,a_src,xlower,xupper,ylower,yupper,nxpoints,nypoints)

# write out parameters
setup_params = dict(
topo_file='topo.tt2',
q1_file='q1.tt2',
aux5_file = 'aux5.tt2',
xlower=float(xlower),
ylower=float(ylower),
yupper=float(yupper),
xupper=float(xupper),
L=L, 
W=W, 
H=H)

with open("params.yaml", "w") as file:
    yaml.dump(setup_params, file, sort_keys=True, default_flow_style=False)
