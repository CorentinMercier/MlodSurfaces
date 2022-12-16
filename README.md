# Moving Level-of-Detail Surfaces


>**Moving Level-Of-Detail Surfaces** <br/>
*Corentin Mercier, Thibault Lescoat, Pierre Roussillon, Tamy Boubekeur and Jean-Marc Thiery.*<br/>
ACM Transaction On Graphics 2022<br/>
DOI: 10.1145/3528223.3530151<br/>

You can find the paper and the video presentation here: https://perso.telecom-paristech.fr/boubek/papers/MLoDSurfaces/

This is not the exact same version of the source code that is used to measure performance for the corresponding paper. Performance might have been slightly affected during refactoring.

Copyright(C) 2022 Corentin Mercier, Thibault Lescoat, Pierre Roussillon, Tamy Boubekeur and Jean-Marc Thiery

All right reserved

# Installation instructions

This code was tested under Linux using Ubuntu 18.04 and 20.04.
It compiles and runs under Windows 10 using Visual Studio 2019, however, performance might be worse than under Linux.

Code dependencies:
- Qt5
- CUDA (version 10.1 tested)
- cmake (>= 3.12)
- openmp

Linux compilation instructions:

```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j
cd ../bin
```

A stand alone version of the code, compiled for x64 and tested on windows 10 can be directly executed in the subfolder 
.../StandaloneBinary/MLoDS.exe

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

Source code for the submission:
   Moving Level-of-Detail Surfaces.
   C. Mercier, T. Lescoat, P. Roussillon, T. Boubekeur and J-M. Thiery
   ACM Transaction On Graphics 2022
   DOI: 10.1145/3528223.3530151

All rights reserved. Use of this source code is governed by a
MIT license that can be found in the LICENSE file.