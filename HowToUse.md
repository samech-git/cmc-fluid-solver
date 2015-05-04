# Setup #

Requirements:
  * Windows XP/Vista/7
  * Visual Studio 2008
  * CUDA toolkit 4.0

Main solutions:
  * \trunk\src\FluidSolver2D.sln for 2D solver <br>
<ul><li>\trunk\src\FluidSolver3D.sln for 3D solver <br></li></ul>

After building corresponding solutions you should get all executables in \trunk\bin\Release or \trunk\bin\Debug folders.<br>
<br>
<h1>How to Run</h1>

Use Release executables for large tests.<br>
<br>
<ul><li>2D solver<br>
<pre><code> FluidSolver2D &lt;data file&gt; &lt;output file&gt; &lt;config file&gt; </code></pre>
</li></ul><blockquote><code>&lt;data file&gt;</code> - represents input data file with initial and boundary condifitions <br>
<code>&lt;output file&gt;</code> - will contain results in specified format <br>
<code>&lt;config file&gt;</code> - used for different parameters of the solver/grids/formats/etc <br></blockquote>

<ul><li>3D solver<br>
<pre><code> FluidSolver3D &lt;data file&gt; &lt;project name&gt; &lt;config file&gt; [GPU] [align] [transpose] </code></pre>
</li></ul><blockquote><code>&lt;data file&gt;</code>  - represents input data file with initial and boundary condifitions <br>
<code>&lt;project name&gt;</code>  - results are stored in <project name><code>_</code>res.nc using specified format, the grid is saved in <project name><code>_</code>3d folder <br>
<code>&lt;config file&gt;</code>  - used for different parameters of the solver/grids/formats/etc <br>
<code>[GPU]</code>  - (optional) use GPU hardware instead of CPU <br>
<code>[align]</code>  - (optional) align grid sizes to be multiple of 32 <br>
<code>[transpose]</code>  - (optional) transpose data arrays, GPU-specific optimization<br></blockquote>

<h1>Sample Command Lines</h1>

<pre><code> FluidSolver2D ..\..\data\2D\box_pipe\box_pipe_data.txt output.txt ..\..\data\2D\box_pipe\box_pipe_config.txt </code></pre>
Runs a 2D solver on CPU for a sample 2D box_pipe test with the corresponding configuration file. <br> <pre><code> FluidSolver3D ..\..\data\3D\example_tests\box_pipe\box_pipe_2D_data.txt box_pipe_example ..\..\data\3D\example_tests\box_pipe\box_pipe_2D_config.txt align GPU transpose</code></pre>
Runs a 3D solver on GPU for an example test that uses box pipe 2D model extruded along Z direction with corresponding configuration file, also using aligned grids and GPU-only transpose optimization.<br>
<br>
There are a couple of example bat files you can run once you build the project. These batch files test simple scenarios in data folder and usually take just a few minutes to complete:<br>
<pre><code> trunk\bin\Release\run_examples_CPU.bat<br>
trunk\bin\Release\run_examples_GPU.bat </code></pre>

Also you can run full large white sea test on CPU/GPU with the following bat files in svn. Please note that this will run much longer since it uses high-res grids and more time steps for better results.<br>
<pre><code> trunk\bin\Release\run_white_sea_CPU.bat<br>
trunk\bin\Release\run_white_sea_GPU.bat </code></pre>

<h1>How to Visualize</h1>

Both 2D and 3D solvers can output axis and variables data (velocity and temperature) in standard NetCDF binary format. <br>
This format is wide-known and supported in many scientific data visualization programs.<br>
<br>
