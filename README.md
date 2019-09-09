## Topographic Analyis Kit for TopoToolbox ##
These are a series of MATLAB functions written by Adam M. Forte [aforte8 'at' lsu.edu] that build upon the functionality of TopoToolbox [https://github.com/wschwanghart/topotoolbox]. Each function contains a header with basic functionality info along with expected inputs and possible outputs. There is now a detailed user guide included in the repository as both a PDF ('TAKmanual_reduce.pdf') and a LaTeX file (folder 'TAKmanual_texFiles' which inludes the tex file and associated images).

# Release
v.1.1.0 has been released! The release includes all files in the 'master' repository along with a tarball bundling the example data sets referenced in the manual. New releases are created periodically when major changes in the codes occur.

# Compiled Versions
Functions within the 'master' branch and releases have corresponding 'compiled' versions of functions that can be run using the free MATLAB Runtime Environment for users who do not have access to MATLAB or all the required toolboxes. The '.m' files used to compile along with the binaries for Mac and Windows are found in the 'Compiled_Versions' folder. See manual for more details on the use of the compiled versions and important differences in inputs and outputs. It is NOT recommended that you use the '.m' files in the 'Compiled_Versions' folder in MATLAB as these have been modified from the original functions to make them function properly when compiled, they are included for reference and in case someone else wanted to recompile the functions.

# READ THIS IF YOU ARE USING THE COMPILED VERSIONS
You must run the compiled versions from the Terminal (Mac OS X) or Command Prompt (Windows). If you simply double click on the TAK executable file created at the end of the installation process this will result in an error, specifically an error on line 9 of TAK.m. This is because the program is expecting inputs that you have to provide in the command line. The user manual contains detailed instructions of how to successfully run the compiled versions, including examples for all compiled functions in the appendix.

# Branches
There are two branches on this repository, 'master' and 'Development'. Functions in the 'master' branch are stable (at least I think they are!), fully documented in the manual, and have corresponding functions in the complied versions of TAK. Functions in the 'Development' branch (including versions of functions in 'master' that have been modified) are still in the testing phase and are either only partially documented or completely undocumented in the manual and compiled versions have not yet been implemented. If there is no difference between 'master' and 'Development' branches, this indicates that there are no functions under development that I've released to the public yet.

# Attribution 
If you use or modify these codes and use the results in a publication, please cite the corresponding publication:
A.M. Forte, K.X. Whipple. Short communication: The Topographic Analysis Kit (TAK) for TopoToolbox. Earth Surface Dynamics, 2019, v. 7, p. 87-95, doi: 10.5194/esurf-7-87-2019. 
[https://www.earth-surf-dynam.net/7/87/2019/]

# Error Reporting and Feature Request
If you encounter a bug or have a suggestion for a new feature / improvement the preferred method of communication is to use the 'Issues' function built into GitHub. You can also email Adam [aforte8 'at' lsu.edu]. If you encounter an issue that you know how to fix and are comfortable with how git works, please feel free to fork the code and submit a pull request with fixes and improvements, collaboration is welcome!
