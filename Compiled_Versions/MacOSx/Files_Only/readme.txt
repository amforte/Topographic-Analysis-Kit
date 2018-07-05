TAK Executable

1. Prerequisites for Deployment 

Verify that version 9.4 (R2018a) of the MATLAB Runtime is installed.   
If not, you can run the MATLAB Runtime installer.
To find its location, enter
  
    >>mcrinstaller
      
at the MATLAB prompt.
NOTE: You will need administrator rights to run the MATLAB Runtime installer. 

Alternatively, download and install the Macintosh version of the MATLAB Runtime for R2018a 
from the following link on the MathWorks website:

    http://www.mathworks.com/products/compiler/mcr/index.html
   
For more information about the MATLAB Runtime and the MATLAB Runtime installer, see 
Package and Distribute in the MATLAB Compiler documentation  
in the MathWorks Documentation Center.

2. Files to Deploy and Package

Files to Package for Standalone 
================================
-run_TAK.sh (shell script for temporarily setting environment variables and executing the 
 application)
   -to run the shell script, type
   
       ./run_TAK.sh <mcr_directory> <argument_list>
       
    at Linux or Mac command prompt. <mcr_directory> is the directory 
    where version 9.4 of the MATLAB Runtime is installed or the directory where 
    MATLAB is installed on the machine. <argument_list> is all the 
    arguments you want to pass to your application. For example, 

    If you have version 9.4 of the MATLAB Runtime installed in 
    /mathworks/home/application/v94, run the shell script as:
    
       ./run_TAK.sh /mathworks/home/application/v94
       
    If you have MATLAB installed in /mathworks/devel/application/matlab, 
    run the shell script as:
    
       ./run_TAK.sh /mathworks/devel/application/matlab
-MCRInstaller.zip 
    Note: if end users are unable to download the MATLAB Runtime using the
    instructions in the previous section, include it when building your 
    component by clicking the "Runtime included in package" link in the
    Deployment Tool.
-The Macintosh bundle directory structure TAK.app 
    Note: this can be stored in an archive file with the zip command 
    zip -r TAK.zip TAK.app
    or the tar command 
    tar -cvf TAK.tar TAK.app
-This readme file 



3. Definitions

For information on deployment terminology, go to
http://www.mathworks.com/help and select MATLAB Compiler >
Getting Started > About Application Deployment >
Deployment Product Terms in the MathWorks Documentation
Center.

4. Appendix 

A. Mac systems:
In the following directions, replace MR/v94 by the directory on the target machine where 
   MATLAB is installed, or MR by the directory where the MATLAB Runtime is installed.

If the environment variable DYLD_LIBRARY_PATH is undefined, set it to the following 
   string:

MR/v94/runtime/maci64:MR/v94/sys/os/maci64:MR/v94/bin/maci64

If it is defined, set it to the following:

${DYLD_LIBRARY_PATH}:MR/v94/runtime/maci64:MR/v94/sys/os/maci64:MR/v94/bin/maci64

    For more detailed information about setting the MATLAB Runtime paths, see Package and 
   Distribute in the MATLAB Compiler documentation in the MathWorks Documentation Center.


     
        NOTE: To make these changes persistent after logout on Linux 
              or Mac machines, modify the .cshrc file to include this  
              setenv command.
        NOTE: The environment variable syntax utilizes forward 
              slashes (/), delimited by colons (:).  
        NOTE: When deploying standalone applications, you can
              run the shell script file run_TAK.sh 
              instead of setting environment variables. See 
              section 2 "Files to Deploy and Package".    



5. Launching application using Macintosh finder

If the application is purely graphical, that is, it doesn't read from standard in or 
write to standard out or standard error, it may be launched in the finder just like any 
other Macintosh application.



