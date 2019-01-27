% Topographic-Analysis-Kit (TAK)
% Version 1.0.2  9-Jan-2019
%
% These are a series of MATLAB functions that build upon the functionality 
% of TopoToolbox [https://github.com/wschwanghart/topotoolbox]. Each function
% contains a header with basic functionality info along with expected inputs 
% and possible outputs. See the included manual for additional information
% on tool functionality and discussion of underlying methodology.
%
% If you encounter errors or have suggestions please contact:
%
% Adam M. Forte
% aforte8 'at' lsu.edu
%
% Paper documenting the tools is now published in ESurf:
% Forte, A.M., Whipple, K.X., 2019, Short communication:
% 	The Topographic Analysis Kit (TAK) for TopoToolbox, 
%	Earth Surface Dynamics, v. 7, p. 87-95,
% 	DOI: 10.5194/esurf-7-87-2019.
%
% URL: https://www.earth-surf-dynam.net/7/87/2019/esurf-7-87-2019.html
%
% If you use or modify these tools for use in a publication, please
% cite the above paper.
%
% Required MATLAB Toolbox list:
%
% Image Processing Toolbox
% Mapping Toolbox
% Statistics and Machine Learning Toolbox
% Curve Fitting Toolbox
% Optimization Toolbox
%
% If you do not have all the required toolboxes (can check with 
% CheckTAKDependencies), consider using the compiled versions.
%
% To fully utilize all of TAK, it is recommended that you have
% at least MATLAB 2017b installed.
%
%
% Function List: 
%
% Basin2Raster
% Basin2Shape
% BasinPicker
% BasinStatsPlots
% CatPoly2GRIDobj
% CheckTAKDepedencies
% ClassifyKnicks
% CompileBasinStats
% ConditionDEM
% DippingBedFinder
% FindBasinKnicks
% FindCentroid
% FindThreshold
% KsnChiBatch
% ksncolor
% KsnProfiler
% MakeCombinedSwath
% MakeStreams
% MakeTopoSwath
% Mat2Arc
% PlotIndividualBasins
% ProcessRiverBasins
% ProjectOntoSwath
% RemoveFlats
% SegmentPicker
% SegmentPlotter
% SegmentProjector
% SubDivideBigBasins
%
