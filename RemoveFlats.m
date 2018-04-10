function [DEMn,MASK]=RemoveFlats(DEM,strength)
	% Function takes DEM and attempts a semi-automated routine to remove flat areas with some input from
	% the user to select areas considred to be flat. This function sometimes works reliably, but will never
	% produce as clean a result as manually clipping out flat areas.
	%
	% Required Inputs:
	% 		DEM - GRIDobj of the digital elevation model of your area loaded into the workspace
	%		strength - integer value between 1 and 4 that controls how aggressively the function defines
	%			flat areas, specifically related to the size of the neighborhood the function uses to
	%			connect ares of similar elevation. A strength of 1 = a 3x3 neighborhood, 2=5x5, 3=7x7, and 
	%			4=9x9. If the results of the function do not capture enough of the flat areas in the MASK,
	%			increase the strength and rerun. Similarly, if the function erroneously includes areas that
	%			are not part of what you consider the flats, try decreasing the strength.
	%
	% Outputs:
	%		DEMn - Version of the DEM with idenitifed flat areas masked out (values set to nan)
	%		MASK - Logical GRIDobj, true where area was identified as a flat.
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Last Revised Spring 2018 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	% Parse Inputs
	p = inputParser;
	p.FunctionName = 'RemoveFlats';
	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'strength',@(x) isnumeric(x) && x<=4 && mod(x,1)==0);

	parse(p,DEM,strength);
	DEM=p.Results.DEM;
	strength=p.Results.strength;	

	disp('Preparing Grids...')

	% Define kernels
	nhood1=ones(3,3);
	if strength==1
		nhood2=ones(3,3);
	elseif strength==2
		nhood2=ones(5,5);
	elseif strength==3
		nhood2=ones(7,7);
	elseif strength==4
		nhood2=ones(9,9);
	else
		error('Input to "strength" is not recognized, must be an integer between 1 and 4')
	end

	% Fill sinks
	[DEMs]=fillsinks(DEM);

	% Identify flats
	FL=GRIDobj(DEMs);
	FL=erode(DEMs,nhood1)==DEMs;
	CFL=dilate(FL,nhood2);

	% Label flats
	L=bwlabel(CFL.Z);
	[xdim,ydim]=getcoordinates(DEM);
	L=GRIDobj(xdim,ydim,L);

	disp('Select areas that you consider sinks and then press enter')
	% Prompt user to choose sinks
	f1=figure(1);
	hold on 
	imagesc(DEM);
	hold off
	[x,y]=ginput;
	close(f1)

	% Find label for flats
	[ix]=coord2ind(DEM,x,y);
	l=L.Z(ix);
	l=unique(l);

	% Generate mask
	MASK=GRIDobj(DEM);
	for ii=1:numel(l)
		MASK.Z(L.Z==l(ii))=1;
	end
	MASK.Z=logical(MASK.Z);

	% Produce new DEM
	DEMn=DEM;
	DEMn.Z(MASK.Z)=nan;
end