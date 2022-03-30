# Moving Level-of-Detail Surfaces

# Installation instructions

This code was tested under Linux using Ubuntu 18.04 and 20.04.
It might compile under Windows 10, however, performance might be worse than under Linux.

Code dependencies:
- Qt5
- CUDA (version 10.1 tested)
- cmake (>= 3.12)
- openmp

```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j
cd ../bin
```

# Usage
```
./MLoDS -i NAME_OF_INPUT_FILE [OPTIONS]
Or for manual use: ./MLoDS
Points will be automatically projected on the input pointset file
The output folder determines where the output pointset and the timings will be saved
Formats supported :
- ply
- pn
--------------Options--------------
-h: display the help
-o: path of output folder
-f: path of output file
-c: use CPU only
-d: only compute the dual-contouring
-a knn: use APSS with knn nearest neighbors
-e r: use APSS with ball of radius r, r<0 => automatic radius
-k cutoff sigma1 factor number_of_sigmas: use the multiple gaussian kernel
-r cutoff exponent epsilon: use the rational kernel
-s cutoff sigma: use the single gaussian kernel
-g: use global APSS instead of MLoDS
-m depth: adaptive octree depth
-b: use octree built on input points
-----------------------------------
```
# Commands

Keyboard shortcuts:
- 1: display/hide input point set (white)
- 3: display/hide projected points (green)
- 4: display/hide the dual-contoured surface if computed (red)
- X, Y and Z: align the camera with one axis
- D: compute dual contouring
- R: recompile shaders
- Arrow up/down: increase/decrease point size

# Example data

Some data can be found in the folder data.
They are accompanied by .params file with the same names that contain parameters suited to the model (they can be overwritten by quitting the software using File/Quit)

# License
Source code provided FOR REVIEW ONLY, as part of the submission entitled
"Moving Level-of-Detail Surfaces".

A proper version of this code will be released if the paper is accepted
with the proper licence, documentation and bug fix.
Currently, this material has to be considered confidential and shall not
be used outside the review process.

All right reserved. The Authors
