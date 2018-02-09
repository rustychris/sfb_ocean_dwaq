import xarray as xr
import matplotlib.pyplot as plt
from stompy.plot import plot_utils
from stompy.grid import unstructured_grid
##

nc=xr.open_dataset('dwaqA01/dwaq_map.nc')

g=unstructured_grid.UnstructuredGrid.from_ugrid(nc)

##
potw_surface=nc.potw.isel(time=-1,nFlowMesh_layers=0)
cont_surface=nc.continuity.isel(time=-1,nFlowMesh_layers=0)
storm_surface=nc.stormwater.isel(time=-1,nFlowMesh_layers=0)

##

plt.figure(1).clf()
fig,ax=plt.subplots(num=1)

ccoll=g.plot_cells(values=potw_surface.values,ax=ax,lw=0.5)
ccoll.set_clim([0,1e-2])

#ccoll=g.plot_cells(values=cont_surface.values,ax=ax,lw=0.5)
#ccoll.set_clim([0.9,1.1])

#ccoll=g.plot_cells(values=storm_surface.values,ax=ax,lw=0.5)
#ccoll.set_clim([0,0.1])

ccoll.set_edgecolor('face')

##

# Why are the POTW exchanges not being forced?
# from the input file, EBMUD gets listed as the name for
# boundary 57, 297, 485, 661, and 823
# Boundary numbers start with 1, up to 2925
# That matches pointers, which has negative from entries down to -2925.
# 

poi1=scen.hydro.pointers # happen to have this lying around in the namespace


bc_seg=-823 # as reported in the inp file, so starting with -1
ebmud_bed_exchs=np.nonzero( poi1[:,0]==-823 )[0]
assert len(ebmud_bed_exchs)==1
ebmud_bed_exch=ebmud_bed_exchs[0]

##

hydro=scen.hydro
flo=hydro.flows(hydro.t_secs[10])

# This should have some flow, but it's coming up 0.0
flo[ebmud_bed_exch]

##

# Look at the segments associated with EBMUD,
# then crosswalk that to all exchanges in that column
#  do any of them look like our exchange?
# from the inp file:
ebmud_boundaries=np.array([57, 297, 485, 661, 823])
ebmud_bc_segs=-ebmud_boundaries # as they should appear in 1-based pointers

ebmud_bc_exchs=[]

hydro.infer_2d_elements()

for ebmud_bc_seg in ebmud_bc_segs:
    print("BC seg %d"%(ebmud_bc_seg))
    
    exchs=np.nonzero( poi1[:,0]==ebmud_bc_seg )[0]
    # Each of these bc segs maps to exactly one exchange, no surprise
    ebmud_bc_exchs.append( exchs )

    for exch in exchs:
        inside_seg1=poi1[exch,1]
        elt=hydro.seg_to_2d_element[inside_seg1-1]
        print("    exch0=%6d  =>  inside_seg1=%5d  elt=%5d  flo=%12.5f"%(exch,inside_seg1,elt,flo[exch]))

        hydro.check_volume_conservation_incr(seg_select=[inside_seg1-1],
                                             tidx_select=slice(9,11))

    print()
    print()
        
## For the bed segment:

# All of those map to the same 2d element: 22145
# which are nicely spaced by hydro.n_2d_elements.  no surprise with dense output.
# But none have any flow
hydro.check_volume_conservation_incr(seg_select=[inside_seg1-1],
                                     tidx_select=slice(9,11))

# Error is 2.43 m3/s.  Sounds a lot like ebmud flow.
# 
hydro.check_volume_conservation_incr(tidx_select=slice(9,11))

# INFO:HydroFiles:Bad segments: [0]
# WARNING:HydroFiles:  Worst segment is index 248593
# WARNING:HydroFiles:  Vlast=373514.187500  Vpred=369137.316772  Vnow=373514.187500
# WARNING:HydroFiles:  z error=0.069397 m
# WARNING:HydroFiles:  Condition of dV: Qmag=732294.125000 Qnet=-4376.875000 mag/net=-167.309814
# WARNING:HydroFiles:****************BAD Volume Conservation*************
# WARNING:HydroFiles:  t=    959400   RMSE: 1.769324e-03    Max rel. err: 1.769324e-03
# WARNING:HydroFiles:  1 segments above rel err tolerance

##

bad_segs=[  3265,  15623,  15624,  15625,  15626,  15627,  15634,  15635,  15643,  15646,
            16067,  16068,  16069,  16070,  16071,  16172,  16178,  16179,  16180,  16181,
            16182,  17508,  17509,  25200,  25314,  25966,  26189,  26190,  26191,  46213,
            46527,  46528,  46529,  46531,  46532,  47728,  47780,  47785,  47786,  47807,
            48347,  77008, 117962, 126829, 126831, 134906, 138515, 138589, 138807, 140053,
            156413, 200640, 210975, 248593, 249129, 253190, 256203, 269129, 308527, 365752]

# Are all conservation errors linked to a discharge BC?

for bad_seg in bad_segs:
    exchs=np.nonzero( (poi1[:,1]-1==bad_seg)&(poi1[:,0]<0) )[0]
    print("Bad seg %d ==>   BC exchs %s"%(bad_seg,exchs))

# maybe 29 of those line up with a BC
# Other might be evaporation?

##

six.moves.reload_module(waq)
hydro=waq.HydroFiles(hyd_path='../../delft/sfb_ocean/sfb_ocean/runs/short_test_13/global/com-merged.hyd')

##

hydro.adjust_boundaries_for_conservation()

##

