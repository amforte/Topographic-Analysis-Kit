function cmpInspectJunction(wdir,MatFile,JunctionMatFile,num,varargin)
	% Description:
	%	Function for visualizing the fit of up- and downstream links to evaluate
	%	the values for the junction angles. It will display both a zoomed view of
	%	the junction, displayed in polar coordinates, and the location of the junction
	%	on a map of the stream network. 
	%
	% Required Inputs:
	%	wdir - full path of working directory
	%	MatFile - Full path of matfile output from either 'cmpMakeStreams' or the name of a single basin mat 
	%		file from 'cmpProcessRiverBasins', in this case it must be the same matfile you provided to JunctionAngle
	%	JunctionMatFile - Full path of matfile output from JunctionAngle
	%	num - junction number, referenced from JunctionAngle table. If
	%		an empty array is provided, you will be prompted with a map 
	%		of the stream network to select the junction of interest.
	%	
	% Optional Inputs:
	%	fit_distance - distance in map units along stream to
	%		do fit of selected junction. Distance must be
	%		less than the maximum distance used when JunctionAngle
	%		was run. This input (and the related 'num_nodes') is provided
	%		if you want to visually explore the sensitivity of the
	%		fit distance to junction angles
	%	num_nodes - number of stream nodes to do fit of selected
	%		junction. Number of nodes must be less than the number
	%		of nodes extracted when JunctionAngle was run.
	%	save_fig [false] - logical flag to toggle saving the figure output of 
	%		this function as a PDF.
	%
	% Note: You cannot provide entries to both fit distance and num_nodes
	%
	%
    % Examples if running for the command line, minus OS specific way of calling main TAK function:
	%	InspectJunction /path/to/wdir Topo.mat Topo_junctions.mat 50
	%	InspectJunction /path/to/wdir Topo.mat Topo_junctions.mat 50 fit_distance 500
	%	InspectJunction /path/to/wdir Topo.mat Topo_junctions.mat 50 num_nodes 10	
	%
	% Related Functions:
	% JunctionAngle, JunctionLinks
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Updated : 06/26/19 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if isdeployed
		if ~isempty(varargin)
			varargin=varargin{1};
		end
	end

	% Parse Inputs
	p = inputParser;
	p.FunctionName = 'cmpInspectJunction';
	addRequired(p,'wdir',@(x) ischar(x));
	addRequired(p,'MatFile',@(x) ~isempty(regexp(x,regexptranslate('wildcard','*.mat'))));
	addRequired(p,'JunctionMatFile',@(x) ~isempty(regexp(x,regexptranslate('wildcard','*.mat'))));	
	addRequired(p,'num',@(x) isnumeric(x) && isscalar(x) || isempty(x));

	addParameter(p,'fit_distance',[],@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'num_nodes',[],@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'save_fig',false,@(x) isscalar(x) && islogical(x));

	parse(p,wdir,MatFile,JunctionMatFile,num,varargin{:});
	wdir=p.Results.wdir;
	MatFile=p.Results.MatFile;
	JMat=p.Results.JunctionMatFile;
	num=p.Results.num;

	fit_distance=p.Results.fit_distance;
	num_nodes=p.Results.num_nodes;
	save_fig=p.Results.save_fig;

	if ~isempty(fit_distance) & ~isempty(num_nodes)
		error('Provide a value to either fit_distance of num_nodes, not both')
	end

	% Determine the type of input
	MatFile=fullfile(wdir,MatFile);
	D=whos('-file',MatFile);
	VL=cell(numel(D),1);
	for ii=1:numel(D);
		VL{ii}=D(ii,1).name;
	end

	if any(strcmp(VL,'S'))
		load(MatFile,'S');
	elseif any(strcmp(VL,'Sc'))
		load(MatFile,'Sc');
		S=Sc;
	end

	JMat=fullfile(wdir,JMat);
	load(JMat,'IX');	

	if isempty(num)
		cix=cell2mat(IX(:,1));
		axc=S.x(cix);
		ayc=S.y(cix);

		f1=figure(1);
		hold on
		plot(S,'-k');
		axis equal
		title('Click near junction of interest');
		hold off

		[xpick,ypick]=ginput(1);
		close(f1);
		[~,num]=min(hypot(axc-xpick,ayc-ypick));
	end

	c=IX{num,1};
	ds=IX{num,2};
	us1=IX{num,3};
	us2=IX{num,4};

	% Accepts variable argument of a fit distance
	if ~isempty(fit_distance)
		n=floor(fit_distance/hypot(S.cellsize,S.cellsize));
		if n<1
			warning('Input fit_distance is too short, will not average across more than one node, setting to minimum distance');
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
	elseif ~isempty(num_nodes)
		n=floor(num_nodes);
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
	end

	if isempty(S.georef)
		xl=S.x; yl=S.y;
		projcoord=false;
	else
		try
			[yl,xl]=projinv(S.georef.mstruct,S.x,S.y);
			projcoord=true;
		catch
			xl=S.x; yl=S.y;
			projcoord=false;
		end
	end

	% Confluence point
	xc=xl(c); yc=yl(c);
	% Upstream links
	xu1=xl(us1); yu1=yl(us1);
	xu2=xl(us2); yu2=yl(us2);
	% Downstream link
	xd=xl(ds); yd=yl(ds);

	xu1=xu1-xc; xu2=xu2-xc;
	yu1=yu1-yc; yu2=yu2-yc;
	xd=xd-xc; yd=yd-yc;

	theta1ap=median(atan2(yu1,xu1)); 
	theta2ap=median(atan2(yu2,xu2));
	theta3ap=median(atan2(yd,xd));

	% Rotate to horizontal
	[xu1r,yu1r]=RotCoord(xu1,yu1,-theta1ap,0,0); 
	[xu2r,yu2r]=RotCoord(xu2,yu2,-theta2ap,0,0);	
	[xdr,ydr]=RotCoord(xd,yd,-theta3ap,0,0);

	% Simple approximation of linear segment on
	% rotated positions
	a1=xu1r\yu1r;
	a2=xu2r\yu2r; 
	a3=xdr\ydr;

	% Find projected y positions
	yp1r=xu1r*a1;
	yp2r=xu2r*a2;
	yp3r=xdr*a3;

	% Rotate back
	[xp1,yp1]=RotCoord(xu1r,yp1r,theta1ap,0,0);
	[xp2,yp2]=RotCoord(xu2r,yp2r,theta2ap,0,0);
	[xp3,yp3]=RotCoord(xdr,yp3r,theta3ap,0,0);	

	% Convert to polar
	[theta1,rho1]=cart2pol(xp1,yp1);
	[theta2,rho2]=cart2pol(xp2,yp2);
	[theta3,rho3]=cart2pol(xp3,yp3);

	[theta1o,rho1o]=cart2pol(xu1,yu1);
	[theta2o,rho2o]=cart2pol(xu2,yu2);
	[theta3o,rho3o]=cart2pol(xd,yd);

	%Find maximum radii and thetas for those radii
	[mrho1,ix1]=max(rho1);
	[mrho2,ix2]=max(rho2);
	[mrho3,ix3]=max(rho3);
	theta1=theta1(ix1);
	theta2=theta2(ix2);
	theta3=theta3(ix3);


	f1=figure(1);
	set(f1,'unit','normalized','position',[0.1 0.1 0.8 0.4],'renderer','painters');
	clf 
	subplot(1,2,1);
	polarplot([theta1 theta1],[0 mrho1],'-r','LineWidth',2);
	hold on 
	polarplot([theta2 theta2],[0 mrho2],'-b','LineWidth',2);
	polarplot([theta3 theta3],[0 mrho3],'-k','LineWidth',2);
	if theta3>=0
		theta3proj=mod(theta3,-pi);
	else
		theta3proj=mod(theta3,pi);
	end

	polarplot([theta3proj theta3proj],[0 mrho3],':k','LineWidth',2);
	polarscatter(theta1o,rho1o,20,'r');
	polarscatter(theta2o,rho2o,20,'b');
	polarscatter(theta3o,rho3o,20,'k');
	ax=gca;
	ax.ThetaTickLabel={'90','60','30','0','330','300','270','240','210','180','150','120'};
	ax.RTickLabel={''};
	title(['Junction ' num2str(num)]);
	legend('Trib 1 Projected','Trib 2 Projected','Downstream Projected','Upstream Plane','Trib 1','Trib 2','Downstream');
	hold off

	subplot(1,2,2)
	hold on 
	plot(S,'-k');
	scatter(S.x(us1),S.y(us1),10,'r');
	scatter(S.x(us2),S.y(us2),10,'b');
	scatter(S.x(ds),S.y(ds),10,'k');
	scatter(S.x(c),S.y(c),20,'g','filled');
	axis equal
	hold off

	if save_fig
		figFile=fullfile(wdir,['Junction_' num2str(num) '.pdf']);
		orient(f1,'Landscape');
		print(f1,'-dpdf',figFile,'-fillpage');
	end	

end

function [n_x,n_y]=RotCoord(x,y,theta,x0,y0)
	n_x=(x-x0).*cos(theta)-(y-y0).*sin(theta);
	n_y=(x-x0).*sin(theta)+(y-y0).*cos(theta);
end
