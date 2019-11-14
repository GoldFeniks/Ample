# Acoustic
Wide-anlge mode parabolic equation (WAMPE) solver for Helmholtz equation.

## Requirements
* [ALGLIB](https://www.alglib.net/)
* [nlohmann::json](https://github.com/nlohmann/json)

## Building
You can build ACOUSTIC on Linux with [CMake](https://cmake.org/) as follows
```bash
$ cd build/cmake
$ mkdir ACOUSTIC && cd ACOUSTIC
$ cmake -DCMAKE_BUILD_TYPE=Release ..
```
## Usage
`ACOUSTIC [solution|modes] [options]`
* `solution` Compute WAMPE solution (default)
* `modes` Compute and save acoustic modes
### Generic options
* `-h [ --help ]` Print help message.
* `-v [ --verbosity ] arg` Set verbosity level to `arg`. Currently three levels are supported: 
  * `0` report nothing (default value),
  * `1` show execution time,
  * `2` print progress.
* `-r [ --report ] k` If verbosity level > 0 report every `k` computed rows. If set to zero reports nothing (default).
* `-c [ --config ] arg` Specifies path to config file. Default value: `config.json`.
### Output options
* `-o [ --output ] arg` Specifies path to output file. Default value: `output.txt`.
* `-s [ --step ] k` Output every `k`-th computed row. By default set to `100`.
* `--binary` If specified uses binary intead of text output.
### Computation options
* `-w [ --workers ] arg` Sets number of threads to be used for computation. By default uses one thread.
* `-b [ --buff ] arg` Sets buffer size for multithreaded computation. Default size is `100`.
## Configuration
Configuration is a file in JSON format containing any of the following fields.
### Floating point fields
* `"mode_subset"` Value in range `[-1, 1]`, used to truncate computed modes.
* `"x0", "x1"` Domain borders over `x` coordinate (`x0` is currently unused).
* `"y0", "y1"` Domain borders over `y` coordinate.
* `"f"` Source frequency.
* `"z_s"` Source depth.
* `"z_r"` Receiver depth.
* `"y_s"` `y` coordinate of the source.
* `"bottom_rho"` Bottom layers density.
### Integer fields
* `"max_mode"` Maximum number of modes to use.
* `"nx", "ny"` Number of points over `x` and `y` coordinates.
* `"ppm"` Number of points for modes computation.
* `"mnx", "mny"` Number of points over `x` and `y` coordinates for modes computation.
* `"ordRich"` Order of Richardson extrapolation.
* `"n_layers"` Number of water layers.
* `"past_n"` History length for transparent boundary conditions.
* `"border_width"` Width of smoothed areas over left and rights domain borders.
### Boolean fiedls
* `"complex_modes"` Uses complex-valued modes (accounts for attenuation).
* `"const_modes"` Modes are assumed to be `x`-independent.
* `"additive_depth"` Add bottom layer depth instead of setting it (see [Array fields](https://github.com/GoldFeniks/Acoustic#array-fields)).
### Array fields
* `"betas"` Attenuation coeffitients for all layers.
* `"bottom_layers"` Depths of bottom layers.
* `"bottom_c1s", "bottom_c2s"` Sound speed at the top and bottom of each bottom layer.
* `"k0", "phi_s"` Wavenumbers and modal functions of the source.
### Bathymetry
Can be specified in one of the following forms.
#### Text file
```json
"bathymetry": [
  "text_file",
  "filename"
]
```
Assumes that data is stored in a text file as a table.
```
 0  y0  ... yM
x0  d00 ... d0M
... ... ... ...
xN  dN0 ... dNM
```
#### Binary file
```json
"bathymetry": [
  "binary_file",
  "filename"
]
```
Assumes that data is stored in a binary file.

Type | Length | Description
-----|--------|------------
`uint32` | `1`  | Unsigned number `N`
`uint32` | `1`  | Unsigned number `M`
`double` | `N`  | `x` coordinates
`double` | `M`  | `y` coordinates
`double` | `NM` | Depth data

#### Explicit specification
```json
"bathymetry": [
  "values",
  {
    "x": [],
    "y": [],
    "values": [ [], [], ]
  }
]
```
### Hydrology
The text and explicit representation is analogous to bathymetry, except `"y"` is replaced with `"z"` and missing values can be specifeid as `-1`.
#### Binary file
Assumes that data is stored in a following format.

Type | Length | Description
-----|--------|------------
`uint32` | `1` | Unsigned number `N`
⁠| `N` | Profile data
`uitn32` | `1` | Unsigned number `M`
⁠| `M` | Sound speed data
`double` | `1` | Depth
`double` | `1` | Sound speed
### Modes
In each representation modes data complex values are stored as two consecutive real values. For binary data integer values are stored as `uint32` and real values as `double`.
#### Text and binary files
The data for `x`-independent modes is stored in a following format (`k` is wavenumber and `p` is modal function).
```
M K
y0 .. yM
k00 ... k0M
... ... ...
kK0 ... kKM
p00 ... p0M
... ... ...
pK0 ... pKM
```
And for `x`-dependent.
```
N M K
x0 ... xN
y0 ... yN
k000 ... k00M
k100 ... k10M
.... ... ....
kN0K ... kNMK
p000 ... p00M
.... ... ....
pN0K ... pNMK
```
#### Explicit specification
```json
"modes": [
  "values",
  {
    "x": [],
    "y": [],
    "k": [ [ [], [], ], [], ],
    "phi": [ [ [], [], ], [], ]
  }
]
```
For `x`-independent modes "y" can be omitted, "k" and "phi" are two-dimensional.
## Output file
Similar to [modes](https://github.com/GoldFeniks/Acoustic#modes) output file consists of two integer values (number of points over `x` and `y` coordinates) followed by respective points specifications. The rest of the file contains `NM` pairs of real values — WAMPE solution.
