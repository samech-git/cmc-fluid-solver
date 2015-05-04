# General Overview #

Configuration file should be written in plain text with 1 line per option. <br>
Some of the options are required for successfull run, some have default values so you can bypass them in the config. <br>
Standard config files for tests are available in the same directories as data files in SVN. <br> <br>
Example: \trunk\data\3D\example_tests\box_pipe\box_pipe_2D_config.txt. <br>

<h1>Details</h1>

<table><thead><th> <b>name</b> </th><th> <b>required</b> </th><th> <b>value</b> </th><th> <b>description</b> </th></thead><tbody>
<tr><td> dimension </td><td> yes </td><td> <code>2D,3D</code> </td><td> Problem dimension </td></tr>
<tr><td> in_fmt </td><td> yes </td><td> <code>Shape2D,Shape3D,SeaNetCDF</code> </td><td> Input data format </td></tr>
<tr><td> viscosity </td><td> no </td><td> <code>&lt;float&gt;</code> </td><td> Viscosity of the fluid (default is 0.05) </td></tr>
<tr><td> density </td><td> no </td><td> <code>&lt;float&gt;</code> </td><td> Density of the fluid (default is 1000.0) </td></tr>
<tr><td> Re </td><td> no </td><td> <code>&lt;float&gt;</code> </td><td> Reynolds number. Must specify if you use normalized parameters </td></tr>
<tr><td> Pr </td><td> no </td><td> <code>&lt;float&gt;</code> </td><td> Prandtl number. Must specify if you use normalized parameters </td></tr>
<tr><td> lambda </td><td> no </td><td> <code>&lt;float&gt;</code> </td><td> Heat capacity ratio. Must specify if you use normalized parameters </td></tr>
<tr><td> bc_type </td><td> no </td><td> <code>NoSlip,Slip</code> </td><td> Type of boundary conditions: slip or no-slip (default is <code>NoSlip</code>) </td></tr>
<tr><td> bc_strength </td><td> no </td><td> <code>&lt;float&gt;</code> </td><td> Boundary slip strength from 0..1 (default is 0.5) </td></tr>
<tr><td> bc_initv </td><td> no </td><td> <code>&lt;float3&gt;</code> </td><td> Input stream velocity (for sea NetCDF tests) </td></tr>
<tr><td> grid_dx </td><td> yes </td><td> <code>&lt;float&gt;</code> </td><td> Grid step size for X dimension </td></tr>
<tr><td> grid_dy </td><td> yes </td><td> <code>&lt;float&gt;</code> </td><td> Grid step size for Y dimension </td></tr>
<tr><td> grid_dz </td><td> 3D </td><td> <code>&lt;float&gt;</code> </td><td> Grid step size for Z dimension </td></tr>
<tr><td> cycles </td><td> no </td><td> <code>&lt;integer&gt;</code> </td><td> Number of cycles (related to the heart beating, default is 1) </td></tr>
<tr><td> frame_time </td><td> SeaNetCDF </td><td> <code>&lt;float&gt;</code> </td><td> Frame time used for NetCDF input format </td></tr>
<tr><td> time_steps </td><td> no </td><td> <code>&lt;integer&gt;</code> </td><td> Number of time steps per one frame used for compute (default is 50) </td></tr>
<tr><td> out_time_steps </td><td> no </td><td> <code>&lt;integer&gt;</code> </td><td> Output time step results % out_time_steps (default is 10) </td></tr>
<tr><td> out_gridx </td><td> no </td><td> <code>&lt;integer&gt;</code> </td><td> Output grid size for X dimension (default is 50) </td></tr>
<tr><td> out_gridy </td><td> no </td><td> <code>&lt;integer&gt;</code> </td><td> Output grid size for Y dimension (default is 50) </td></tr>
<tr><td> out_gridz </td><td> no </td><td> <code>&lt;integer&gt;</code> </td><td> Output grid size for Z dimension (default is 50) </td></tr>
<tr><td> out_fmt </td><td> yes </td><td> <code>MultiVox, NetCDF</code> </td><td> Output format: standard NetCDF in plain text (<code>NetCDF</code>) or specific output for MultiVox product (<code>MultiVox</code>) </td></tr>
<tr><td> depth </td><td> Shape2D </td><td> <code>&lt;float&gt;</code> </td><td> Depth of the 3D model in Z dimension (only for Shape2D format in 3D mode). </td></tr>
<tr><td> solver </td><td> yes </td><td> <code>Explicit</code> <br> <code>ADI</code> <br> <code>Stable</code> </td><td> <code>Explicit</code> - simple explicit finite-difference scheme  <br> <code>ADI</code> - alternating direction implicit method with non-linear iterations <br> <code>Stable</code> - standard advection & poisson equation for pressure </td></tr>
<tr><td> num_global </td><td> no </td><td> <code>&lt;integer&gt;</code> </td><td> Global non-linear iterations (default is 2) </td></tr>
<tr><td> num_local </td><td> no </td><td> <code>&lt;integer&gt;</code> </td><td> Local non-linear iterations (default is 1) </td></tr></tbody></table>

Notes:<br>
<ul><li>You must use either viscosity/density (non-normalized) or Re/Pr/lambda (normalized) parameters for the flow