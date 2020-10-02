function JunctionSegments(S,IX,junctions)
	% Usage:
	%	JunctionSegments(S,IX,junctions);
	%
	% Description:
	%	Function to create a shapefile(s) showing the segments up and downstream
	%	of the confluence used to measure the junction angle. Will produce
	%	one shapefile for each fit distance you provided when running JunctionAngle 
	%
	% Required Inputs:
	%	S - STREAMobj
	%	IX - IX cell array output from Junction Angle
	%	junctions - junctions table output from JunctionAngle, if you provided multiple
	%		fit_distances to JunctionAngle, the output 'links' will be a cell
	%		array of similar dimensions as the input 'junctions'.
	%
	% Examples:
	%	JunctionSegments(S,IX,junctions);
	%
	% Related Functions:
	% JunctionAngle, JunctionLinks, InspectJunction
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Created : 10/02/20 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	% Parse Inputs
	p = inputParser;
	p.FunctionName = 'JunctionSegments';
	addRequired(p,'S',@(x) isa(x,'STREAMobj'));
	addRequired(p,'IX',@(x) iscell(x));
	addRequired(p,'junctions',@(x) istable(x) | iscell(x));

	parse(p,S,IX,junctions);
	S=p.Results.S;
	IX=p.Results.IX;
	JIN=p.Results.junctions;

	if istable(JIN)
		J=JIN;
		ms=struct;

		% Thin the junction table and IX table removing any non-viable junctions
		ja=J.junction_angle;
		idx=~isnan(ja);
		J=J(idx,:);
		IX=IX(idx,:);

		w1=waitbar(0,'Building shapefile');
		% Start loop
		num_juncs=size(J,1);
		ind=1;
		for ii=1:num_juncs
			ds=IX{ii,2};
			us1=IX{ii,3};
			us2=IX{ii,4};

			ms(ind,1).Geometry='PolyLine';
			ms(ind,1).X=S.x(us1);
			ms(ind,1).Y=S.y(us1);
			ms(ind,1).Confl_ID=J.junction_number(ii);
			ms(ind,1).Part='upstream_1';			
			ind=ind+1;
			ms(ind,1).Geometry='PolyLine';
			ms(ind,1).X=S.x(us2);
			ms(ind,1).Y=S.y(us2);
			ms(ind,1).Confl_ID=J.junction_number(ii);
			ms(ind,1).Part='upstream_2';	
			ind=ind+1;
			ms(ind,1).Geometry='PolyLine';
			ms(ind,1).X=S.x(ds);
			ms(ind,1).Y=S.y(ds);
			ms(ind,1).Confl_ID=J.junction_number(ii);
			ms(ind,1).Part='downstream';
			ind=ind+1;
			waitbar(ii/num_juncs);
		end
		close(w1);
		shapewrite(ms,'junction_segments.shp');
	else
		
		num_lengths=size(JIN,2);

		for jj=1:num_lengths
			% Grab junction table
			J=JIN{2,jj};
			fit_distance=JIN{1,jj};
			ms=struct;

			% Thin the junction table and IX table removing any non-viable junctions
			ja=J.junction_angle;
			idx=~isnan(ja);
			J=J(idx,:);
			IXN=IX(idx,:);

			w1=waitbar(0,['Building shapefile ' num2str(jj) ' of ' num2str(num_lengths)]);
			% Start loop
			num_juncs=size(J,1);
			ind=1;
			for ii=1:num_juncs
				ds=IXN{ii,2};
				us1=IXN{ii,3};
				us2=IXN{ii,4};

				% Thin out based on distances
				n=floor(fit_distance/hypot(S.cellsize,S.cellsize));
				if n<1
					n=1;
				end

				if n<numel(ds)
					ds=ds(1:n);
				end

				if n<numel(us1)
					us1=us1(1:n);
				end

				if n<numel(us2)
					us2=us2(1:n);
				end

				ms(ind,1).Geometry='PolyLine';
				ms(ind,1).X=S.x(us1);
				ms(ind,1).Y=S.y(us1);
				ms(ind,1).Confl_ID=J.junction_number(ii);
				ms(ind,1).Part='upstream_1';			
				ind=ind+1;
				ms(ind,1).Geometry='PolyLine';
				ms(ind,1).X=S.x(us2);
				ms(ind,1).Y=S.y(us2);
				ms(ind,1).Confl_ID=J.junction_number(ii);
				ms(ind,1).Part='upstream_2';	
				ind=ind+1;
				ms(ind,1).Geometry='PolyLine';
				ms(ind,1).X=S.x(ds);
				ms(ind,1).Y=S.y(ds);
				ms(ind,1).Confl_ID=J.junction_number(ii);
				ms(ind,1).Part='downstream';
				ind=ind+1;
				waitbar(ii/num_juncs);
			end
			close(w1);
			shapewrite(ms,['junction_segments_' num2str(jj) '.shp']);
		end
	end