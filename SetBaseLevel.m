function [SC]=SetBaseLevel(DEM,FD,A,S,method,varargin)
	% Function to adjust base level of streams within a network via either an elevation or drainage area condition.
	% 
	% Reqiured Inputs:
	% 	DEM - DEM Grid Object (assumes unconditioned DEM)
	% 	FD - FLOW object
	% 	A - GRID object of flow accumulations
	%	S - STREAM object
	%	method - controls how stream network base level is adjusted. Options for control of base level are:
	%		'elevation' - extract streams only above a given elevation (provided by the user using the 'min_elevation' parameter) to ensure that base level
	%			elevation for all streams is uniform. If the provided elevation is too low (i.e. some outlets of the unaltered stream network are above this
	%			elevation) then a warning will be displayed, but the code will still run.
	%		'drain_area' - extract streams only below a given maximum drainage area (provided by the user using the 'max_drainage_area' parameter) to ensure
	%			that the outlets of all extracted streams have the same drainage areas. If the provided maximum drainage area is too large (i.e. some outlets
	%			have drainage areas smaller than this maximum) then a warning will be displayed, but the code will still run.
	%		'max_out_elevation' - uses the maximum elevation of all stream outlets to extract streams only above this elevation
	%		'min_out_drain_area' - uses the minimum drainage area of all stream outlets to extract streams only below this drainage area
	%
	% Optional Inputs:
	%	min_elevation [] - parameter to set minimum elevation for base level, required if 'method' is set to 'elevation'
	%	max_drainage_area [] - parameter to set maximum drainage area for base level, required if 'method' is set to 'drain_area'
	%
	% Example:
	%	[SN]=SetBaseLevel(DEM,FD,A,S,'max_out_elevation');
	%	[SN]=SetBaseLevel(DEM,FD,A,S,'drain_area','max_drainage_area',1e8);
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Last Revised Winter 2017 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	% Parse Inputs
	p = inputParser;         
	p.FunctionName = 'SetBaseLevel';
	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'FD', @(x) isa(x,'FLOWobj'));
	addRequired(p,'A', @(x) isa(x,'GRIDobj'));
	addRequired(p,'S',@(x) isa(x,'STREAMobj'));
	addRequired(p,'method',@(x) ischar(validatestring(x,{'elevation','drain_area','max_out_elevation','min_out_drain_area'})));

	addParamValue(p,'min_elevation',[],@(x) isnumeric(x));
	addParamValue(p,'max_drainage_area',[],@(x) isnumeric(x));

	parse(p,DEM,FD,A,S,method,varargin{:});
	DEM=p.Results.DEM;
	FD=p.Results.FD;
	A=p.Results.A;
	S=p.Results.S;
	method=p.Results.method;

	me=p.Results.min_elevation;
	ma=p.Results.max_drainage_area;

	% Check that cut off values have been provided
	if strcmp(method,'elevation') & isempty(me)
		error('Selected method "elevation" requires that you provide an input for parameter "min_elevation"');
	elseif strcmp(method,'drain_area') & isempty(ma)
		error('Selected method "drain_area" requires that you provide an input for parameter "max_drainage_area"');
	end

	switch method
	case 'elevation'
		coix=streampoi(S,'outlets','ix');
		coel=DEM.Z(coix);
		max_coel=max(coel);
		if sum(coel>me)~=0
			warning(['One or more stream outlets are above the provided elevation, maximum outlet elevation is ' num2str(max_coel)]);
		end
		st_el=getnal(S,DEM);
		idx=st_el>=me;
		IX=S.IXgrid(idx);
		W=GRIDobj(DEM);
		W.Z(IX)=1;
		W.Z=logical(W.Z);
		SC=STREAMobj(FD,W);
	case 'drain_area'
		coix=streampoi(S,'outlets','ix');
		DA=A.*(A.cellsize^2);
		coda=DA.Z(coix);
		min_coda=min(coda);
		if sum(coda<ma)~=0
			warning(['One or more steam outlets have drainage areas less than provided maximum drainage area, minimum outlet drainage area is ' num2str(min_coda)]);
		end
		st_da=getnal(S,DA);
		idx=st_da<=ma;
		IX=S.IXgrid(idx);
		W=GRIDobj(DEM);
		W.Z(IX)=1;
		W.Z=logical(W.Z);
		SC=STREAMobj(FD,W);	
	case 'max_out_elevation'
		coix=streampoi(S,'outlets','ix');
		coel=DEM.Z(coix);
		max_coel=max(coel);
		st_el=getnal(S,DEM);
		idx=st_el>=max_coel;
		IX=S.IXgrid(idx);
		W=GRIDobj(DEM);
		W.Z(IX)=1;
		W.Z=logical(W.Z);
		SC=STREAMobj(FD,W);
	case 'min_out_drain_area'
		coix=streampoi(S,'outlets','ix');
		DA=A.*(A.cellsize^2);
		coda=DA.Z(coix);
		min_coda=min(coda);
		st_da=getnal(S,DA);
		idx=st_da<=min_coda;
		IX=S.IXgrid(idx);
		W=GRIDobj(DEM);
		W.Z(IX)=1;
		W.Z=logical(W.Z);
		SC=STREAMobj(FD,W);	
	end
end