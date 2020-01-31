#!/usr/bin/env python

## \file gradients.py
#  \brief python package for gradients
#  \author T. Lukaczyk, F. Palacios
#  \version 6.2.0 "Falcon"
#
# The current SU2 release has been coordinated by the
# SU2 International Developers Society <www.su2devsociety.org>
# with selected contributions from the open-source community.
#
# The main research teams contributing to the current release are:
#  - Prof. Juan J. Alonso's group at Stanford University.
#  - Prof. Piero Colonna's group at Delft University of Technology.
#  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
#  - Prof. Rafael Palacios' group at Imperial College London.
#  - Prof. Vincent Terrapon's group at the University of Liege.
#  - Prof. Edwin van der Weide's group at the University of Twente.
#  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
#
# Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
#                      Tim Albring, and the SU2 contributors.
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import os, sys, shutil, copy
from .. import run  as su2run
from .. import io   as su2io
from .. import util as su2util
from .functions import function, update_mesh
from ..io import redirect_folder, redirect_output
from SU2.eval import functions

# ----------------------------------------------------------------------
#  Main Gradient Interface
# ----------------------------------------------------------------------

def gradient( func_name, method, config, state=None ):    
    """ val = SU2.eval.grad(func_name,method,config,state=None)

        Evaluates the aerodynamic gradients.

        Wraps:
            SU2.eval.adjoint()
            SU2.eval.findiff()

        Assumptions:
            Config is already setup for deformation.
            Mesh need not be deformed.
            Updates config and state by reference.
            Redundancy if state.GRADIENTS has the key func_name.

        Executes in:
            ./ADJOINT_* or ./FINDIFF

        Inputs:
            func_name - SU2 objective function name
            method    - 'CONTINUOUS_ADJOINT' or 'FINDIFF' or 'DISCRETE_ADJOINT'
            config    - an SU2 config
            state     - optional, an SU2 state

        Outputs:
            A list of floats of gradient values
    """

    # Initialize
    grads = {}
    state = su2io.State(state)
    if func_name == 'ALL':
        raise Exception("func_name = 'ALL' not yet supported")
    func_output = func_name
    if (type(func_name)==list):
        if (config.OPT_COMBINE_OBJECTIVE=="YES"):
            func_output = 'COMBO'
        else:
            func_name = func_name[0]
    else:
        config.OPT_COMBINE_OBJECTIVE="NO"
        config.OBJECTIVE_WEIGHT = "1.0"

    # redundancy check
    if not func_output in state['GRADIENTS']:

        # Adjoint Gradients
        if any([method == 'CONTINUOUS_ADJOINT', method == 'DISCRETE_ADJOINT']):
              
            # Aerodynamics
            if func_output in su2io.optnames_aero + su2io.optnames_turbo:
                grads = adjoint( func_name, config, state )

            elif func_name in su2io.optnames_aero + su2io.optnames_turbo:
                grads = adjoint( func_name, config, state )
                
            # Stability
            elif func_output in su2io.optnames_stab:
                grads = stability( func_name, config, state )
            
            # Multipoint
            elif func_output in su2io.optnames_multi:
              grads = multipoint( func_name, config, state )

            # Geometry (actually a finite difference)
            elif func_output in su2io.optnames_geo:
                grads = geometry( func_name, config, state )

            else:
                raise Exception('unknown function name: %s' % func_name)

        # Finite Difference Gradients
        elif method == 'FINDIFF':
            grads = findiff( config, state )

        elif method == 'DIRECTDIFF':
            grad = directdiff (config , state )

        else:
            raise Exception('unrecognized gradient method')
        
        # store
        state['GRADIENTS'].update(grads)

    # if not redundant

    # prepare output
    grads_out = state['GRADIENTS'][func_output]

    return copy.deepcopy(grads_out)

#: def gradient()


# ----------------------------------------------------------------------
#  Adjoint Gradients
# ----------------------------------------------------------------------

def adjoint( func_name, config, state=None ):
    """ vals = SU2.eval.adjoint(func_name,config,state=None)

        Evaluates the aerodynamics gradients using the 
        adjoint methodology with:
            SU2.eval.func()
            SU2.run.deform()
            SU2.run.direct()
            SU2.run.adjoint()

        Assumptions:
            Config is already setup for deformation.
            Mesh may or may not be deformed.
            Updates config and state by reference.
            Adjoint Redundancy if state.GRADIENTS has key func_name.
            Direct Redundancy if state.FUNCTIONS has key func_name.

        Executes in:
            ./ADJOINT_<func_name>

        Inputs:
            func_name - SU2 objective function name
            config    - an SU2 config
            state     - optional, an SU2 state

        Outputs:
            A Bunch() with keys of objective function names
            and values of list of floats of gradient values
    """

    # ----------------------------------------------------
    #  Initialize    
    # ----------------------------------------------------

    # initialize
    state = su2io.State(state)
    special_cases = su2io.get_specialCases(config)
    
    # When a list of objectives is used, they are combined 
    # and the output name is 'COMBO'
    multi_objective = (type(func_name)==list)
    func_output = func_name
    if multi_objective:   func_output = 'COMBO'

    ADJ_NAME = 'ADJOINT_'+func_output

    # console output
    if config.get('CONSOLE','VERBOSE') in ['QUIET','CONCISE']:
        log_adjoint = 'log_Adjoint.out'
    else:
        log_adjoint = None   

    # ----------------------------------------------------
    #  Redundancy Check
    # ----------------------------------------------------    

    # master redundancy check
    if func_output in state['GRADIENTS']:
        grads = state['GRADIENTS']
        return copy.deepcopy(grads)

    # ----------------------------------------------------
    #  Direct Solution    
    # ----------------------------------------------------        

    # run (includes redundancy checks)
    function( func_name, config, state )   

    # ----------------------------------------------------    
    #  Adaptation (not implemented)
    # ----------------------------------------------------

    #if not state.['ADAPTED_ADJOINT']:
    #    config = su2run.adaptation(config)
    #    state['ADAPTED_FUNC'] = True

    # ----------------------------------------------------    
    #  Adjoint Solution
    # ----------------------------------------------------        

    # files to pull
    files = state['FILES']
    pull = []; link = []    

    # files: mesh
    name = files['MESH']
    name = su2io.expand_part(name,config)
    link.extend(name)

    # files: direct solution
    name = files['DIRECT']
    name = su2io.expand_zones(name,config)
    name = su2io.expand_time(name,config)
    link.extend(name)

    # files: adjoint solution
    if ADJ_NAME in files:
        name = files[ADJ_NAME]
        name = su2io.expand_zones(name,config)
        name = su2io.expand_time(name,config)
        link.extend(name)       
    else:
        config['RESTART_SOL'] = 'NO'

    # files: target equivarea adjoint weights
    if 'EQUIV_AREA' in special_cases:
        pull.append(files['WEIGHT_NF'])   
        pull.append(files['TARGET_EA'])

    # files: target pressure coefficient
    if 'INV_DESIGN_CP' in special_cases:
        pull.append(files['TARGET_CP'])

    # files: target heat flux coefficient
    if 'INV_DESIGN_HEATFLUX' in special_cases:
        pull.append(files['TARGET_HEATFLUX'])

    # output redirection
    with redirect_folder( ADJ_NAME, pull, link ) as push:
        with redirect_output(log_adjoint):        

            # Format objective list in config
            if multi_objective:
                config['OBJECTIVE_FUNCTION'] = ", ".join(func_name)
            else:
                config['OBJECTIVE_FUNCTION'] = func_name

            # # RUN ADJOINT SOLUTION # #
            info = su2run.adjoint(config)
            su2io.restart2solution(config,info)
            state.update(info)

            # Gradient Projection
            info = su2run.projection(config,state)
            state.update(info)

            # solution files to push
            name = state.FILES[ADJ_NAME]
            name = su2io.expand_zones(name,config)
            name = su2io.expand_time(name,config)
            push.extend(name)

    #: with output redirection

    # return output 
    grads = su2util.ordered_bunch()
    grads[func_output] = state['GRADIENTS'][func_output]
    return grads

#: def adjoint()


# ----------------------------------------------------------------------
#  Stability Functions
# ----------------------------------------------------------------------

def stability( func_name, config, state=None, step=1e-2 ):


    folder = 'STABILITY' # os.path.join('STABILITY',func_name) #STABILITY/D_MOMENT_Y_D_ALPHA/

    # ----------------------------------------------------
    #  Initialize    
    # ----------------------------------------------------

    # initialize
    state = su2io.State(state)
    if not 'MESH' in state.FILES:
        state.FILES.MESH = config['MESH_FILENAME']
    special_cases = su2io.get_specialCases(config)

    # find base func name
    matches = [ k for k in su2io.optnames_aero if k in func_name ]
    if not len(matches) == 1:
        raise Exception('could not find stability function name')
    base_name = matches[0]    

    ADJ_NAME = 'ADJOINT_'+base_name

    # console output
    if config.get('CONSOLE','VERBOSE') in ['QUIET','CONCISE']:
        log_direct = 'log_Direct.out'
    else:
        log_direct = None

    # ----------------------------------------------------    
    #  Update Mesh
    # ----------------------------------------------------

    # does decomposition and deformation
    info = update_mesh(config,state) 

    # ----------------------------------------------------    
    #  CENTRAL POINT
    # ----------------------------------------------------    

    # will run in ADJOINT/
    grads_0 = gradient(base_name,'CONTINUOUS_ADJOINT',config,state)


    # ----------------------------------------------------    
    #  Run Forward Point
    # ----------------------------------------------------   

    # files to pull
    files = state.FILES
    pull = []; link = []

    # files: mesh
    name = files['MESH']
    name = su2io.expand_part(name,config)
    link.extend(name)

    # files: direct solution
    ## DO NOT PULL DIRECT SOLUTION, use the one in STABILITY/ 

    # files: adjoint solution
    if ADJ_NAME in files:
        name = files[ADJ_NAME]
        name = su2io.expand_time(name,config)
        link.extend(name)       
    else:
        config['RESTART_SOL'] = 'NO'        

    # files: target equivarea adjoint weights
    ## DO NOT PULL EQUIVAREA WEIGHTS, use the one in STABILITY/


    # pull needed files, start folder
    with redirect_folder( folder, pull, link ) as push:
        with redirect_output(log_direct):     

            konfig = copy.deepcopy(config)
            ztate  = copy.deepcopy(state)

            # TODO: GENERALIZE
            konfig.AOA = konfig.AOA + step

            # let's start somethin somthin
            del ztate.GRADIENTS[base_name]
            #ztate.find_files(konfig)

            # the gradient
            grads_1 = gradient(base_name,'CONTINUOUS_ADJOINT',konfig,ztate)


    # ----------------------------------------------------    
    #  DIFFERENCING
    # ----------------------------------------------------

    grads = [ ( g_1 - g_0 ) / step 
              for g_1,g_0 in zip(grads_1,grads_0) ]

    state.GRADIENTS[func_name] = grads
    grads_out = su2util.ordered_bunch()
    grads_out[func_name] = grads

    return grads_out


# ----------------------------------------------------------------------
#  Multipoint Functions
# ----------------------------------------------------------------------

def multipoint( func_name, config, state=None, step=1e-2 ):
  
    mach_list = config['MULTIPOINT_MACH_NUMBER'].replace("(", "").replace(")", "").split(',')
    reynolds_list = config['MULTIPOINT_REYNOLDS_NUMBER'].replace("(", "").replace(")", "").split(',')
    freestream_temp_list = config['MULTIPOINT_FREESTREAM_TEMPERATURE'].replace("(", "").replace(")", "").split(',')    
    freestream_press_list = config['MULTIPOINT_FREESTREAM_PRESSURE'].replace("(", "").replace(")", "").split(',')
    aoa_list = config['MULTIPOINT_AOA'].replace("(", "").replace(")", "").split(',')
    sideslip_list = config['MULTIPOINT_SIDESLIP_ANGLE'].replace("(", "").replace(")", "").split(',')
    target_cl_list = config['MULTIPOINT_TARGET_CL'].replace("(", "").replace(")", "").split(',')
    weight_list = config['MULTIPOINT_WEIGHT'].replace("(", "").replace(")", "").split(',')

    grads = []
    folder = []
    for i in range(len(weight_list)):
      grads.append(0)
      folder.append(0)

    for i in range(len(weight_list)):
      folder[i] = 'MULTIPOINT_' + str(i)
    
    # ----------------------------------------------------
    #  Initialize
    # ----------------------------------------------------
    
    # initialize
    state = su2io.State(state)
    if not 'MESH' in state.FILES:
      state.FILES.MESH = config['MESH_FILENAME']
    special_cases = su2io.get_specialCases(config)
    
    # find base func name
    matches = [ k for k in su2io.optnames_aero if k in func_name ]
    if not len(matches) == 1:
        raise Exception('could not find multipoint function name')
    base_name = matches[0]
    
    ADJ_NAME = 'ADJOINT_' + base_name
    
    # console output
    if config.get('CONSOLE','VERBOSE') in ['QUIET','CONCISE']:
      log_direct = 'log_Direct.out'
    else:
      log_direct = None
  
#    # ----------------------------------------------------
#    #  Update Mesh
#    # ----------------------------------------------------
#    
#    # does decomposition and deformation
#    info = update_mesh(config,state)
    
    # ----------------------------------------------------
    #  FIRST POINT
    # ----------------------------------------------------
    
    # will run in ADJOINT/

    config.AOA = aoa_list[0]
    config.SIDESLIP_ANGLE = sideslip_list[0]
    config.MACH_NUMBER = mach_list[0]
    config.REYNOLDS_NUMBER = reynolds_list[0]
    config.FREESTREAM_TEMPERATURE = freestream_temp_list[0]
    config.FREESTREAM_PRESSURE = freestream_press_list[0]
    config.TARGET_CL = target_cl_list[0]   

    grads[0] = gradient(base_name,'DISCRETE_ADJOINT',config,state)

    src = os.getcwd()
    src = os.path.abspath(src).rstrip('/') + '/' + ADJ_NAME + '/'


    # ----------------------------------------------------
    #  Run Forward Point
    # ----------------------------------------------------
    
    # files to pull
    files = state.FILES
    pull = []; link = []
    
    # files: mesh
    name = files['MESH']
    name = su2io.expand_part(name,config)
    link.extend(name)
    
    # files: direct solution
    ## DO NOT PULL DIRECT SOLUTION, use the one in MULTIPOINT/
    
    # files: adjoint solution
    if ADJ_NAME in files:
      name = files[ADJ_NAME]
      name = su2io.expand_time(name,config)
      link.extend(name)
    else:
      config['RESTART_SOL'] = 'NO'

    # files: target equivarea adjoint weights
    ## DO NOT PULL EQUIVAREA WEIGHTS, use the one in MULTIPOINT/

    # pull needed files, start folder
    with redirect_folder( folder[0], pull, link ) as push:
        with redirect_output(log_direct):
    
          konfig = copy.deepcopy(config)
          ztate  = copy.deepcopy(state)
        
          dst = os.getcwd()
          dst = os.path.abspath(dst).rstrip('/')+'/'

        # make unix link
          string = "ln -s " + src + " " + dst
          os.system(string)

    for i in range(len(weight_list)-1):

      # pull needed files, start folder
      with redirect_folder( folder[i+1], pull, link ) as push:
        with redirect_output(log_direct):
    
          konfig = copy.deepcopy(config)
          ztate  = copy.deepcopy(state)
        
          konfig.AOA = aoa_list[i+1]
          konfig.SIDESLIP_ANGLE = sideslip_list[i+1]
          konfig.MACH_NUMBER = mach_list[i+1]
          konfig.REYNOLDS_NUMBER = reynolds_list[i+1]
          konfig.FREESTREAM_TEMPERATURE = freestream_temp_list[i+1]
          konfig.FREESTREAM_PRESSURE = freestream_press_list[i+1]
          konfig.TARGET_CL = target_cl_list[i+1]         
 
          # let's start somethin somthin
          del ztate.GRADIENTS[base_name]
          #ztate.find_files(konfig)
            
          # the gradient
          grads[i+1] = gradient(base_name,'DISCRETE_ADJOINT',konfig,ztate)
        
    # ----------------------------------------------------
    #  WEIGHT FUNCTIONS
    # ----------------------------------------------------
    
    grad = []
    for variable in range(len(grads[0])):
      grad.append(0)

    for variable in range(len(grads[0])):
      grad[variable] = 0.0
      for point in range(len(weight_list)):
        grad[variable] = grad[variable] + float(weight_list[point])*grads[point][variable]
      
    state.GRADIENTS[func_name] = grad
    grads_out = su2util.ordered_bunch()
    grads_out[func_name] = grad
                 
    return grads_out


# ----------------------------------------------------------------------
#  Finite Difference Gradients
# ----------------------------------------------------------------------

def findiff( config, state=None ):
    """ vals = SU2.eval.findiff(config,state=None)

        Evaluates the aerodynamics gradients using 
        finite differencing with:
            SU2.eval.func()
            SU2.run.deform()
            SU2.run.direct()

        Assumptions:
            Config is already setup for deformation.
            Mesh may or may not be deformed.
            Updates config and state by reference.
            Gradient Redundancy if state.GRADIENTS has the key func_name.
            Direct Redundancy if state.FUNCTIONS has key func_name.

        Executes in:
            ./FINDIFF

        Inputs:
            config - an SU2 config
            state  - optional, an SU2 state

        Outputs:
            A Bunch() with keys of objective function names
            and values of list of floats of gradient values
    """    

    # ----------------------------------------------------
    #  Initialize    
    # ----------------------------------------------------

    # initialize
    state = su2io.State(state)
    special_cases = su2io.get_specialCases(config)
    Definition_DV = config['DEFINITION_DV']

    # console output
    if config.get('CONSOLE','VERBOSE') in ['QUIET','CONCISE']:
        log_findiff = 'log_FinDiff.out'
    else:
        log_findiff = None

    # evaluate step length or set default value
    if 'FIN_DIFF_STEP' in config:
        step = float(config.FIN_DIFF_STEP)
    else:
        step = 0.001 

    # ----------------------------------------------------
    #  Redundancy Check
    # ----------------------------------------------------    

    # master redundancy check
    opt_names = su2io.optnames_aero + su2io.optnames_geo 
    findiff_todo = all([key in state.GRADIENTS for key in opt_names])
    if findiff_todo:
        grads = state['GRADIENTS']
        return copy.deepcopy(grads)

    # ----------------------------------------------------
    #  Zero Step  
    # ----------------------------------------------------   

    # run
    func_base = function( 'ALL', config, state )      

    # ----------------------------------------------------
    #  Plot Setup
    # ----------------------------------------------------          

    grad_filename  = config['GRAD_OBJFUNC_FILENAME']
    grad_filename  = os.path.splitext( grad_filename )[0]
    output_format  = config['OUTPUT_FORMAT']
    plot_extension = su2io.get_extension(output_format)    
    grad_filename  = grad_filename + '_findiff' + plot_extension

    # ----------------------------------------------------
    #  Finite Difference Steps
    # ----------------------------------------------------  

    # local config
    konfig = copy.deepcopy(config)

    # check deformation setup
    n_dv = sum(Definition_DV['SIZE'])
    deform_set = konfig['DV_KIND'] == Definition_DV['KIND']
    if not deform_set: 
        dvs_base = [0.0] * n_dv
        konfig.unpack_dvs(dvs_base,dvs_base)    
    else:
        dvs_base = konfig['DV_VALUE_NEW']

    # initialize gradients
    func_keys = list(func_base.keys())
    func_keys = ['VARIABLE'] + func_keys + ['FINDIFF_STEP']
    grads = su2util.ordered_bunch.fromkeys(func_keys)
    for key in grads.keys(): grads[key] = []

    # step vector
    if isinstance(step,list):
        assert n_dv == len(step) , 'unexpected step vector length'
    else:
        step = [step] * n_dv

    # files to pull
    files = state['FILES']
    pull = []; link = []    
    # files: mesh
    name = files['MESH']
    name = su2io.expand_part(name,konfig)
    link.extend(name)
    # files: direct solution
    if 'DIRECT' in files:
        name = files['DIRECT']
        name = su2io.expand_time(name,config)
        link.extend(name)

    # files: target equivarea distribution
    if 'EQUIV_AREA' in special_cases and 'TARGET_EA' in files:
        pull.append(files['TARGET_EA'])

    # files: target pressure distribution
    if 'INV_DESIGN_CP' in special_cases and 'TARGET_CP' in files:
        pull.append(files['TARGET_CP'])

    # files: target heat flux distribution
    if 'INV_DESIGN_HEATFLUX' in special_cases and 'TARGET_HEATFLUX' in files:
        pull.append(files['TARGET_HEATFLUX'])

       
    # output redirection
    with redirect_folder('FINDIFF',pull,link) as push:
        with redirect_output(log_findiff):

            # iterate each dv    
            for i_dv in range(n_dv):

                this_step = step[i_dv]
                temp_config_name = 'config_FINDIFF_%i.cfg' % i_dv 

                this_dvs    = copy.deepcopy(dvs_base)
                this_konfig = copy.deepcopy(konfig)
                this_dvs[i_dv] = this_dvs[i_dv] + this_step

                this_state = su2io.State()
                this_state.FILES = copy.deepcopy( state.FILES )
                this_konfig.unpack_dvs(this_dvs,dvs_base)

                this_konfig.dump(temp_config_name)

                # Direct Solution, findiff step
                func_step = function( 'ALL', this_konfig, this_state )

                # remove deform step files
                meshfiles = this_state.FILES.MESH
                meshfiles = su2io.expand_part(meshfiles,this_konfig)
                for name in meshfiles: os.remove(name)

                # calc finite difference and store
                for key in grads.keys():
                    if key == 'VARIABLE': 
                        grads[key].append(i_dv)
                    elif key == 'FINDIFF_STEP': 
                        grads[key].append(this_step)
                    else:
                        this_grad = ( func_step[key] - func_base[key] ) / this_step
                        grads[key].append(this_grad)
                        
                   
                #: for each grad name
                    
                su2util.write_plot(grad_filename,output_format,grads)
                os.remove(temp_config_name)

            #: for each dv

    #: with output redirection

    # remove plot items
    del grads['VARIABLE']
    del grads['FINDIFF_STEP']
    state.GRADIENTS.update(grads)

    # return results
    grads = copy.deepcopy(grads)
    return grads

#: def findiff()


# ----------------------------------------------------------------------
#  Geometric Gradients
# ----------------------------------------------------------------------

def geometry( func_name, config, state=None ):
    """ val = SU2.eval.geometry(config,state=None)

        Evaluates geometry with the following:
            SU2.run.deform()
            SU2.run.geometry()

        Assumptions:
            Config is already setup for deformation.
            Mesh may or may not be deformed.
            Updates config and state by reference.
            Redundancy if state.FUNCTIONS does not have func_name.

        Executes in:
            ./GEOMETRY

        Inputs:
            config    - an SU2 config
            state     - optional, an SU2 state

        Outputs:
            Bunch() of functions with keys of objective function names
            and values of objective function floats.
    """

    # ----------------------------------------------------
    #  Initialize    
    # ----------------------------------------------------

    # initialize
    state = su2io.State(state)
    if not 'MESH' in state.FILES:
        state.FILES.MESH = config['MESH_FILENAME']
    special_cases = su2io.get_specialCases(config)

    # console output
    if config.get('CONSOLE','VERBOSE') in ['QUIET','CONCISE']:
        log_geom = 'log_Geometry.out'
    else:
        log_geom = None

    # ----------------------------------------------------    
    #  Update Mesh (check with Trent)
    # ----------------------------------------------------

    # does decomposition and deformation
    # info = update_mesh(config,state)


    # ----------------------------------------------------
    #  Geometry Solution
    # ----------------------------------------------------    

    # redundancy check
    geometry_done = func_name in state.GRADIENTS
    #geometry_done = all([key in state.FUNCTIONS for key in su2io.optnames_geo])
    if not geometry_done:    

        # files to pull
        files = state.FILES
        pull = []; link = []

        # files: mesh
        name = files['MESH']
        name = su2io.expand_part(name,config)
        link.extend(name)

        # update function name
        ## TODO

        # output redirection
        with redirect_folder( 'GEOMETRY', pull, link ) as push:
            with redirect_output(log_geom):     

                # setup config
                config.GEO_PARAM = func_name
                config.GEO_MODE  = 'GRADIENT'

                # # RUN GEOMETRY SOLUTION # #
                info = su2run.geometry(config)
                state.update(info)

                # no files to push

        #: with output redirection

    #: if not redundant 


    # return output 
    grads = su2util.ordered_bunch()
    for key in su2io.optnames_geo:
        if key in state['GRADIENTS']:
            grads[key] = state['GRADIENTS'][key]
    return grads    

#: def geometry()


# ----------------------------------------------------------------------
#  Direct Differentiation Gradients
# ----------------------------------------------------------------------

def directdiff( config, state=None ):
    """ vals = SU2.eval.directdiff(config,state=None)

        Evaluates the aerodynamics gradients using
        direct differentiation with:
            SU2.eval.func()
            SU2.run.deform()
            SU2.run.direct()

        Assumptions:
            Config is already setup for deformation.
            Mesh may or may not be deformed.
            Updates config and state by reference.
            Gradient Redundancy if state.GRADIENTS has the key func_name.
            Direct Redundancy if state.FUNCTIONS has key func_name.

        Executes in:
            ./DIRECTDIFF

        Inputs:
            config - an SU2 config
            state  - optional, an SU2 state
            step   - finite difference step size, as a float or
                     list of floats of length n_DV

        Outputs:
            A Bunch() with keys of objective function names
            and values of list of floats of gradient values
    """

    # ----------------------------------------------------
    #  Initialize
    # ----------------------------------------------------

    # initialize
    state = su2io.State(state)
    special_cases = su2io.get_specialCases(config)
    Definition_DV = config['DEFINITION_DV']

    # console output
    if config.get('CONSOLE','VERBOSE') in ['QUIET','CONCISE']:
        log_directdiff = 'log_DirectDiff.out'
    else:
        log_directdiff = None

    # ----------------------------------------------------
    #  Redundancy Check
    # ----------------------------------------------------

    # master redundancy check
    opt_names = su2io.optnames_aero + su2io.optnames_geo
    directdiff_todo = all([key in state.GRADIENTS for key in opt_names])
    if directdiff_todo:
        grads = state['GRADIENTS']
        return copy.deepcopy(grads)

    # ----------------------------------------------------
    #  Plot Setup
    # ----------------------------------------------------

    grad_filename  = config['GRAD_OBJFUNC_FILENAME']
    grad_filename  = os.path.splitext( grad_filename )[0]
    output_format  = config['OUTPUT_FORMAT']
    plot_extension = su2io.get_extension(output_format)
    grad_filename  = grad_filename + '_directdiff' + plot_extension

    # ----------------------------------------------------
    # Direct Differentiation Evaluation
    # ----------------------------------------------------

    # local config
    konfig = copy.deepcopy(config)

    n_dv = sum(Definition_DV['SIZE'])

    # initialize gradients
    func_keys = list(su2io.grad_names_map.keys())
    func_keys = ['VARIABLE'] + func_keys
    grads = su2util.ordered_bunch.fromkeys(func_keys)
    for key in grads.keys(): grads[key] = []

    # files to pull
    files = state['FILES']
    pull = []; link = []
    # files: mesh
    name = files['MESH']
    name = su2io.expand_part(name,konfig)
    link.extend(name)
    # files: direct solution
    if 'DIRECT' in files:
        name = files['DIRECT']
        name = su2io.expand_time(name,config)
        link.extend(name)

    # files: target equivarea distribution
    if 'EQUIV_AREA' in special_cases and 'TARGET_EA' in files:
        pull.append(files['TARGET_EA'])

    # files: target pressure distribution
    if 'INV_DESIGN_CP' in special_cases and 'TARGET_CP' in files:
        pull.append(files['TARGET_CP'])

    # files: target heat flux distribution
    if 'INV_DESIGN_HEATFLUX' in special_cases and 'TARGET_HEATFLUX' in files:
        pull.append(files['TARGET_HEATFLUX'])

    # output redirection
    with redirect_folder('DIRECTDIFF',pull,link) as push:
        with redirect_output(log_directdiff):

            # iterate each dv
            for i_dv in range(n_dv):

                temp_config_name = 'config_DIRECTDIFF_%i.cfg' % i_dv

                this_konfig = copy.deepcopy(konfig)

                this_dvs = [0.0]*n_dv
                this_dvs[i_dv] = 1.0
                this_dvs_old = [0.0]*n_dv
                this_dvs_old[i_dv] = 1.0
                this_state = su2io.State()
                this_state.FILES = copy.deepcopy( state.FILES )
                this_konfig.unpack_dvs(this_dvs, this_dvs_old)

                this_konfig.dump(temp_config_name)

                # Direct Solution
                func_step = function( 'ALL', this_konfig, this_state )

                # store
                for key in grads.keys():
                    if key == 'VARIABLE':
                        grads[key].append(i_dv)
                    else:
                        if su2io.grad_names_map[key] in func_step:
                          this_grad = func_step[su2io.grad_names_map[key]]
                        else:
                          this_grad = 0.0
                        grads[key].append(this_grad)
                #: for each grad name

                su2util.write_plot(grad_filename,output_format,grads)
                os.remove(temp_config_name)

            #: for each dv

    #: with output redirection

    # remove plot items
    del grads['VARIABLE']
    state.GRADIENTS.update(grads)
    state.update(this_state)

    # return results
    grads = copy.deepcopy(grads)
    return grads

#: def directdiff()

