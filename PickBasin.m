function [DEMc,FDc,Ac,Sc,drainage_area]=PickBasin(DEM,FD,S,thresh_area)
	% Function to select a single drainage basin from an existing stream network and output
	%	a clipped DEM and the drainage area of the drainage basin in km^2 along with a FLOWobj,
	%	flow accumulation grid, and STREAMobj for the clipped basin.
	%
	% Inputs:
	%	DEM - Digital Elevation as a GRIDobj
	%	FD - Flow direction grid as FLOWobj
	%	S - Stream network as STREAMobj
	%	thresh_area - Threshold drainage area for defining streams in m^2
	%
	% Outputs:
	%	DEMc - Clipped DEM as GRIDobj
	%	FDc - FLOWobj calculated from clipped GRIDobj
	%	Ac - Flow accumulation from clipped GRIDobj
	%	Sc - STREAMobj from clipped GRIDobj
	%	drainage_area - Drainage area in km^2 of the clipped basin
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Last Revised Spring 2016 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	str='N';

	while strcmp(str,'N')==1;	
		close all
		f1=figure(1);
		set(f1,'Units','normalized','Position',[0.05 0.1 0.45 0.8],'renderer','painters');

		hold on
		imagesc(DEM);
		plot(S,'-k');
		hold off

		disp('Choose mouth of desired drainage basin')
		[x,y]=ginput(1);

		% Build dependenc map and clip out drainage basins
		[xn,yn]=snap2stream(S,x,y);
		I=dependencemap(FD,xn,yn);
		DEMc=crop(DEM,I,nan);

		% Calculate drainage area
		dep_map=GRIDobj2mat(I);
		num_pix=sum(sum(dep_map));
		drainage_area=(num_pix*DEMc.cellsize*DEMc.cellsize)/(1e6);

		clf
		hold on
		imageschs(DEM,I);
		plot(S,'-k');
		scatter(xn,yn,20,'r','filled');
		hold off

		prompt='Is this the basin you wanted? Y/N [Y]: ';
		str=input(prompt,'s');
		if isempty(str)
			str = 'Y';
		end
	end

	FDc=FLOWobj(DEMc,'preprocess','carve');
	Ac=flowacc(FDc);
	DEM_res=DEMc.cellsize;
	min_area=floor(thresh_area/(DEM_res*DEM_res));
	isstream=Ac>min_area;
	Sc=STREAMobj(FDc,isstream);
% Function End
end


