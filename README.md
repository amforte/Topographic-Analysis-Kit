## Topographic Analyis Kit for TopoToolbox ##
These are a series of MATLAB functions written by Adam M. Forte [aforte8 'at' lsu.edu] that build upon the functionality of TopoToolbox [https://github.com/wschwanghart/topotoolbox]. Each function contains a header with basic functionality info along with expected inputs and possible outputs. There is now a detailed user guide included in the repository as both a PDF ('TAKmanual_reduce.pdf') and a LaTeX file (folder 'TAKmanual_texFiles' which inludes the tex file and associated images).

# Release
v.1.2.0 has been released! The release includes all files in the 'master' repository along with a tarball bundling the example data sets referenced in the manual. New releases are created periodically when major changes in the codes occur.

# GUI
Starting with v.1.2.0, TAK comes with a GUI which incorporates almost all of the standard TAK functions. This can be run either in MATLAB or through the Matlab Runtime Environment.

# Tutorial Files
Starting with v.1.2.0, the release includes a detailed tutorial for standard usage of TAK. The tutorial comes in two versions, one for the MATLAB function version and one for the GUI version. The tutorial is only available through downloading a release (because of file size limitations for files in the repository). 

# Old Compiled Versions
The old compiled functions are deprecated and no longer being supported. They are still included in the repository for reference, but issues raised with these will be ignored. These have been replaced with a GUI version of TAK that runs in both MATLAB and through the use of the Matlab Runtime Environment

# Branches
There are two branches on this repository, 'master' and 'Development'. Functions in the 'master' branch are stable (at least I think they are!), fully documented in the manual, and have corresponding functions in the complied versions of TAK. Functions in the 'Development' branch (including versions of functions in 'master' that have been modified) are still in the testing phase and are either only partially documented or completely undocumented in the manual and compiled versions have not yet been implemented. If there is no difference between 'master' and 'Development' branches, this indicates that there are no functions under development that I've released to the public yet.

# Attribution 
If you use or modify these codes and use the results in a publication, please cite the corresponding publication:
A.M. Forte, K.X. Whipple. Short communication: The Topographic Analysis Kit (TAK) for TopoToolbox. Earth Surface Dynamics, 2019, v. 7, p. 87-95, doi: 10.5194/esurf-7-87-2019. 
[https://www.earth-surf-dynam.net/7/87/2019/]

# Error Reporting and Feature Request
If you encounter a bug or have a suggestion for a new feature / improvement the preferred method of communication is to use the 'Issues' function built into GitHub. You can also email Adam [aforte8 'at' lsu.edu]. If you encounter an issue that you know how to fix and are comfortable with how git works, please feel free to fork the code and submit a pull request with fixes and improvements, collaboration is welcome!
