function [ExpandedGRID]=GridExpand(GRID_small,GRID_large,fill_val);
	% Function takes a small grid and expands it to the dimensions of a larger grid. Smaller grid must be completely contained within bounds
	% of larger grids. Two grids must be in same projection and have the same resolution (see resample if this is not the case)
	%
	%Inputs:
	% GRID_small - smaller GRIDobj
	% GRID_larger - larger of the two GRIDobjs
	% fill_val - value to fill in the extra space
	%Outputs:
	% ExpandedGRID - smaller grid now the same dimensions as larger grid
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Last Revised Fall 2015 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	[xs,ys]=getoutline(GRID_small);
	[xl,yl]=getoutline(GRID_large);

	[xx,yy]=getcoordinates(GRID_large);
	grid_size=GRID_small.cellsize;

	extra_xl=single((xs(1)-xl(1))/grid_size);
	extra_xr=single((xl(4)-xs(4))/grid_size);

	extra_yt=single((yl(3)-ys(3))/grid_size);
	extra_yb=single((ys(4)-yl(4))/grid_size);

	[Gmat,~,~]=GRIDobj2mat(GRID_small);
	num_rows=size(Gmat,1);
	
	% Add Width
	left_add=ones(num_rows,extra_xl)*fill_val;
	right_add=ones(num_rows,extra_xr)*fill_val;

	Gmat_AddW=horzcat(left_add,Gmat,right_add);


	% Add Height
	num_cols=size(Gmat_AddW,2);

	top_add=ones(extra_yt,num_cols)*fill_val;
	bot_add=ones(extra_yb,num_cols)*fill_val;

	Gmat_AddW_AddH=vertcat(top_add,Gmat_AddW,bot_add);

	% Build New Grid
	ExpandedGRID=GRIDobj(xx,yy,Gmat_AddW_AddH);
end