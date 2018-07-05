function cmpRemoveFlats(wdir,dem,strength,file_name)
	% 	Function takes DEM and attempts a semi-automated routine to remove flat areas with some input from 
	% 	the user to select areas considred to be flat. This function sometimes works reliably, but will 
	% 	never produce as clean a result as manually clipping out flat areas in gis software (but it's
	%	a lot faster!)
	%
	% Required Inputs:
	%		wdir - full path of working directory
	% 		dem - either full path of dem file as either an ascii text file or geotiff 
	%		strength - integer value between 1 and 4 that controls how aggressively the function defines
	%			flat areas, specifically related to the size of the neighborhood the function uses to
	%			connect ares of similar elevation. A strength of 1 = a 3x3 neighborhood, 2=5x5, 3=7x7, and 
	%			4=9x9. If the results of the function do not capture enough of the flat areas in the MASK,
	%			increase the strength and rerun. Similarly, if the function erroneously includes areas that
	%			are not part of what you consider the flats, try decreasing the strength.
	%		file_name - name of output ascii file (without a file type suffix)
	%
	% Outputs:
	%		georeferenced ascii file of the processed dem and mask (mask will have '_mask' appended to the
	%		file name you provide).
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Updated : 07/02/18 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if isdeployed
		if ~isempty(varargin)
			varargin=varargin{1};
		end
	end

	% Parse Inputs
	p = inputParser;
	p.FunctionName = 'cmpRemoveFlats';
	addRequired(p,'wdir',@(x) ischar(x));
	addRequired(p,'dem',@(x) ~isempty(regexp(x,regexptranslate('wildcard','*.txt'))) || ~isempty(regexp(x,regexptranslate('wildcard','*.tif'))));
	addRequired(p,'strength',@(x) isnumeric(x) && x<=4 && mod(x,1)==0);
	addRequired(p,'file_name',@(x) ischar(x));

	parse(p,wdir,dem,strength,file_name);
	wdir=p.Results.wdir;
	dem=p.Results.dem;
	strength=p.Results.strength;
	file_name=p.Results.file_name;	

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

	% Generate GRIDobj
	disp('Loading DEM')
	demLoc=fullfile(wdir,dem);
	DEM=GRIDobj(demLoc);

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

	% Generate gradient for visualizing
	G=gradient8(DEM);

	% Prompt user to choose sinks
	f1=figure(1);
	hold on 
	title('Gradient of DEM : Select areas that you consider sinks and then press enter')
	imageschs(DEM,G,'colormap','jet','caxis',[0 1]);
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

	% Export the processed DEM
	fn=fullfile(wdir,[file_name '.txt']);
	GRIDobj2ascii(DEMn,fn);
	maskfn=fullfile(wdir,[file_name '_mask.txt']);
	GRIDobj2ascii(MASK,maskfn);

	% Generate output to test results
	f1=figure(1);
	hold on
	imageschs(DEMn,MASK);
	hold off
end