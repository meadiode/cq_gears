#! /usr/bin/python3

'''
CQ_Gears - CadQuery based involute profile gear generator

Copyright 2021 meadiode@github

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
'''


class GearBase:
    ka = 1.0  # Addendum coefficient
    kd = 1.25 # Dedendum coefficient

    curve_points = 20 # Number of points to approximate a curve
    surface_splines = 5 # Number of curve splines to approximate a surface
    
    wire_comb_tol = 1e-2 # Wire combining tolerance
    spline_approx_tol = 1e-2 # Surface spline approximation tolerance
    shell_sewing_tol = 1e-2 # Tolerance to assembly a shell out of faces
    isection_tol = 1e-7 # Tolerance to find intersections between two surfaces
    spline_approx_min_deg = 3 # Minimum surface spline degree
    spline_approx_max_deg = 8 # Maximum surface spline degree

    working_plane = 'XY'
    axial_plane = 'XZ'
    rotation_axis = 'Z'

    tooth_trace_curve_segments = (1.0,)

    def __init__(self, *args, **kv_args):
        raise NotImplementedError('Constructor is not defined')

    
    def build(self, **kv_params):
        params = {**self.build_params, **kv_params}
        
        return self._build(**params)
