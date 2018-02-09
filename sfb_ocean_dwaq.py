"""
Try some DWAQ runs with the splice SFB/ocean grid
"""

import os
from collections import defaultdict
import datetime

import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()

import numpy as np
import xarray as xr
import six

from stompy import utils
import stompy.model.delft.waq_scenario as waq

here="."
six.moves.reload_module(waq)

##

PC=waq.ParameterConstant
Sub=waq.Substance
IC=waq.Initial

class Scen(waq.Scenario):
    name="sfb_ocean_A00"
    desc=('sfb_ocean_A00',
          '201708',
          'full')
    integration_option="""16.62 ;
    LOWER-ORDER-AT-BOUND NODISP-AT-BOUND
    BALANCES-OLD-STYLE BALANCES-GPP-STYLE 
    BAL_NOLUMPPROCESSES BAL_NOLUMPLOADS BAL_NOLUMPTRANSPORT
    BAL_NOSUPPRESSSPACE BAL_NOSUPPRESSTIME
    """
    base_path='auto'

    time_step=3000 # matches the hydro

    # multigrid_block now set by default in waq_scenario.py
    sfbay_potw_fn='sfbay_potw/outputs/sfbay_potw.nc'

    def init_substances(self):
        subs=super(Scen,self).init_substances()

        # TODO: unify the handling of src_tags so that it is
        # accessible across the class hierarchy, more standardized.
        self.src_tags=[] # list of dicts
        
        subs['continuity']=Sub(initial=IC(default=1.0))

        subs['potw']=Sub(initial=IC(default=0.0))
        subs['stormwater']=Sub(initial=IC(default=0.0))
        
        link_groups=self.hydro.group_boundary_links()

        ocean=[]
        river=[]
        potw=[]
        delta=[]
        
        for link_group in link_groups:
            if link_group['id']<0:
                continue
            
            name=link_group['name']
            
            if name.startswith('oce'):
                ocean.append(link_group['name'])
            elif name.endswith('_flow'):
                if ('RioVista' in name) or ('Jersey' in name):
                    delta.append(name)
                else:
                    river.append(name)
            else:
                potw.append(name)

        self.src_tags.append(dict(tracer='potw',items=potw))
        self.src_tags.append(dict(tracer='stormwater',items=river))

        self.src_tags.append(dict(tracer='continuity',
                                  items=river+ocean+potw+delta,
                                  value=1.0))

        return subs

    def init_parameters(self):
        # choose which processes are enabled.  Includes some
        # parameters which are not currently used.
        params=super(Scen,self).init_parameters()
        
        params['NOTHREADS']=PC(28) # better to keep it to real cores?

        # z-layer specific items
        params['ZThreshold']=PC(0.005)
        params['ACTIVE_Emersion']=PC(1)

        params['ACTIVE_DYNDEPTH']=PC(1)
        params['ACTIVE_TOTDEPTH']=PC(1)

        # looser than this and the errors are quite visible.  This already feels lenient,
        # but in 2016-06 tighter tolerances led to non-convergence.
        # had been 1.0e-5
        params['Tolerance']=PC(1.0e-6)
        
        # if convergence becomes an issue, there *might* be more info here:
        # params['Iteration Report']=PC(1)

        return params
        
    def cmd_default(self):
        """ Write hydro, DWAQ input files, and run the simulation"""
        self.cmd_write_hydro()
        self.cmd_write_inp()
        self.cmd_delwaq1()
        self.cmd_delwaq2()
        self.cmd_write_nc()

    def __init__(self,*a,**k):
        super(Scen,self).__init__(*a,**k)

        extra_fields=('salinity',
                      'temp',
                      'TotalDepth',
                      'VertDisper',
                      'volume',
                      'depth')
        self.map_format=['binary']
        self.map_output+=extra_fields
        self.hist_output+=extra_fields
        
        self.map_time_step=10000 # hourly
        self.mon_time_step=3000 # half-hourly (had been daily)
        self.hist_time_step=3000 # half-hourly (had been daily)
##

hydro=waq.HydroFiles(hyd_path="../../delft/sfb_ocean/sfb_ocean/runs/short_test_13/global/com-merged.hyd",
                     enable_write_symlink=True)


sec=datetime.timedelta(seconds=1)
start_time=hydro.time0+hydro.t_secs[ 0]*sec
stop_time=hydro.time0+hydro.t_secs[-4]*sec 
    
scen=Scen(hydro=hydro,
          start_time=start_time,
          stop_time=stop_time,
          base_path='dwaqA01')

##

scen.cmd_default()

##

if __name__=='__main__':
   scen.main()


        
