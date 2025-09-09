function [theta_out] = MinimizeDisorder(DEM,FD,A,S,varargin)
	%
	% Usage:
	%	[theta_out] = MinimizeDisorder(DEM,FD,A,S);
	%	[theta_out] = MinimizeDisorder(DEM,FD,A,S,'name',value,...);
	%
	% Description:
	%	Implementation of the methods to estimate a best fit concavity index for a single
	%	watershed descibed in Mudd et al., (2018) which in turn relies on methods from
	%	Hergarten et al., (2016) in the calculation of a "disorder" metric for a watershed.
	%	The minimization of this disorder metric is a method for estimating the concavity.
	%	As described in Mudd et al., (2018), this implements either a simple minimization of 
	%	the disorder metric, whereby a single estimate of the concavity index is produced or 
	% 	can use a bootstrap approach where the trunk stream plus every possible combination of
	% 	three tributaries to that trunk stream are used to calculate a suite of concavity indices.
	%	The latter is more computationally time intensive, but provides a range of estimates that
	%	gives a sense of uncertainty. Mudd et al., (2018) recommends using the median and interquartile
	%	range of the resultant values.
	%
	% Required Inputs:
	%	DEM - Digital Elevation Model as GRIDobj
	% 	FD - Flow direction as FLOWobj
	%	A - GRIDobj that is the result of flowacc(FD)
	%	S - STREAMobj, assumed to be that of a single network, i.e., with a single outlet. If a STREAMobj
	%		with multiple outlets is provided, the code will extract the largest single network for the
	%		calculationa and warn the user.
	%
	% Optional Inputs:
	%	bounded [false] - flag to indicate whether the minimization should be bounded or not in terms of 
	%		whether restrictions are placed on the value of the concavity index (theta) that are permitted.
	%	lower_bound [0] - lower bound on the concavity index allowed if bounded is true.
	%	upper_bound [1.5] - upper bound on the concavity index allowed if bounded is true.
	%	start_val [0.5] - starting value for theta if the minimization is unbounded.
	%	bootstrap [false] - flag to use the bootstrap procedure described in Mudd et al., (2018). If set to false
	%		the output 'theta_out' will be a 1 x 1 array with the single value of theta that is estimated from minizming
	%		disorder across the entire watershed provided. If bootstrap is set to true, then 'theta_out' will be n x 1,
	%		where n is the number of unique set of 3 tributaries that drain to the trunk. Note, that if you are supplying 
	%		a large network or a network where channel definition used a low threshold area (i.e., in a scenario where 
	%		there are a lot of tributaries), turning on bootstrapping may take a long time to complete.
	%	tributary_numbers [3] - the number of tributaries (along with the trunk) to include, along with the trunk
	%		in the bootstrap routine (if bootstrap is set to true). The default, 3, is what is used in Mudd et al., 
	%		(2018). As this number increases, the number of repititions within the bootstrap routine will also increase
	%		as will the computation time.
	%
	% Examples:
	%	[theta_out] = MinimizeDisorder(DEM,FD,A,S,'bounded',true,'lower_bound',0.1,'upper_bound',0.9);
	%	[theta_out] = MinimizeDisorder(DEM,FD,A,S,'bootstrap',true);
	%	best_theta = median(theta_out);
	%
	% References:
	%	Mudd, S.M., Clubb, F.J, Gailleton, B., Hurst, M.D., 2018, How concave are river channels?, Earth
	%		Surface Dynamics, v.6, pg. 505-523, doi: 10.5194/esurf-6-505-2018
	%	Hergarten, S., Robl, J., Stuwe, K., 2016, Tectonic geomorphology at small catchment sizes - extensions
	%		of the stream-power approach and the chi method, Earth Surface Dynamics, v. 4, pg. 1-9,
	%		doi: 10.5194/esurf-4-1-2016
	%	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Created : 09/08/25 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	p = inputParser;
	p.FunctionName = 'MinimizeDisorder';

	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'FD',@(x) isa(x,'FLOWobj'));
	addRequired(p,'A',@(x) isa(x,'GRIDobj'));
	addRequired(p,'S',@(x) isa(x,'STREAMobj'));

	addParameter(p,'start_val',0.5,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'bootstrap',false,@(x) isscalar(x) && islogical(x));
	addParameter(p,'bounded',false,@(x) isscalar(x) && islogical(x));
	addParameter(p,'lower_bound',0.0,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'upper_bound',1.5,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'tributary_numbers',3,@(x) isscalar(x) && isnumeric(x));

	parse(p,DEM,FD,A,S,varargin{:});

	% Check for single outlet
	oix = streampoi(S,'outlets','ix');
	if numel(oix)>1
		disp('Warning: This function can only operate on a stream network with a single outlet, proceeding with the largest single network.')
		S = klargestconncomps(S,1);
	end

	if p.Results.bootstrap
		% Find outlets of tributaries draining to the trunk
		ST = trunk(S);
		Strib = modify(S,'tributaryto',ST);
		trib_oix = streampoi(Strib,'outlets','ix');

		% The bootstrap routine will only produce a non-zero output if the
		% number of tributaries is equal to or greater than the number of
		% tributaries to use in the bootstrap routine. If there is a small
		% number of tributaries, this reverts to a single determination of 
		% theta based on minimization of the disorder.
		if numel(trib_oix) >= p.Results.tributary_numbers
			% Create 3-n list of tributaries to iterate through
			n = nchoosek(trib_oix,p.Results.tributary_numbers);
			num_trib_pairs = size(n,1);

			% Generate a container for the theta results
			theta_out = zeros(num_trib_pairs,1);

			% Begin iterating through tributary pairs
			for ii=1:num_trib_pairs
				% Generate a stream network of just the three selected tributaries
				Sup = modify(S,'upstreamto',n(ii,:));
				% Create logical raster
				L = GRIDobj(DEM);
				L.Z = logical(L.Z);
				% Populate with the trunk, selected tribs, and their outlets
				L.Z(ST.IXgrid) = 1;
				L.Z(Sup.IXgrid) = 1;
				L.Z(n(ii,:)) = 1;
				% Recreate this candidate stream network
				Stest = STREAMobj(FD,L);
				% Do the optimization
				if p.Results.bounded
					theta_out(ii) = fminbnd(@(theta) disorder(Stest,DEM,A,theta),p.Results.lower_bound,p.Results.upper_bound);
				else
					theta_out(ii) = fminsearch(@(theta) disorder(Stest,DEM,A,theta),p.Results.start_val);
				end
			end
		else
			if p.Results.bounded
				theta_out = fminbnd(@(theta) disorder(S,DEM,A,theta),p.Results.lower_bound,p.Results.upper_bound);
			else
				theta_out = fminsearch(@(theta) disorder(S,DEM,A,theta),p.Results.start_val);
			end
		end


	else
		if p.Results.bounded
			theta_out = fminbnd(@(theta) disorder(S,DEM,A,theta),p.Results.lower_bound,p.Results.upper_bound);
		else
			theta_out = fminsearch(@(theta) disorder(S,DEM,A,theta),p.Results.start_val);
		end
	end
end

function [D] = disorder(S,DEM,A,theta_ref)
	% Calculates scaled disorder, sensu Hergarten et al., 2016, ESurf.

	% Calculate chi
	c = chitransform(S,A,'a0',1,'mn',theta_ref);
	% Extract z
	z = getnal(S,DEM);
	% Sort by ascending elevation
	[zs,idx] = sort(z);
	% Apply sorting to chi
	cs = c(idx);
	% Calculate non-normalized disorder
	S = sum(abs(diff(cs)));
	% Calculate normalized disorder
	D = (S - max(c))/max(c);
end