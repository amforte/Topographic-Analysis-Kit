function StreamProjector(DEM,FD,A,S,varargin);
	% Function to interactively select segments of a channel profile you wish to project (e.g. projecting a portion of the profile with a different ksn).
	%	Expects that the STREAMobj you provide will only have one one channel head (i.e. a single profile). You can use the 'SegmentPicker' function to 
	%	interactively choose channels to provide to the StreamProjector function.  
	%	
	% Required Inputs:
	%	DEM - Digital Elevation as a GRIDobj, assumes unconditioned DEM (e.g. DEMoc from ProcessRiverBasins)
	%	FD - Flow direction as FLOWobj
	%	S - Stream network as STREAMobj, expects a channel with a single channel head, will produce an error otherwise.
	%	A - Flow accumulation GRIDobj
	%
	% Optional Inputs:
	%	theta_method ['ref']- options for concavity
	%		'ref' - uses a reference concavity, the user can specify this value with the reference concavity option (see below)
	%		'auto' - function finds a best fit concavity for the provided stream
	%	pick_method ['chi'] - choice of how you want to pick the stream segment to be projected:
	%		'chi' - select segments on a chi - z plot
	%		'stream' - select segments on a longitudinal profile
	%	ref_concavity [0.45] - refrence concavity used if 'theta_method' is set to 'auto'
	%	refit_streams [false] - option to recalculate chi based on the concavity of the picked segment (true), useful if you want to try to precisely 
	%		match the shape of the picked segment of the profile. Only used if 'theta_method' is set to 'auto'
	%		
	% Examples:
	%	StreamProjector(DEM,FD,A,S)
	%	StreamProjector(DEM,FD,A,S,'ref_concavity',0.55);
	%	StreamProjector(DEM,FD,A,S,'theta_method','auto','pick_method','stream','refit_streams',true);
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Last Revised Winter 2017 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Parse Inputs
	p = inputParser;
	p.FunctionName = 'StreamProjector';
	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'FD',@(x) isa(x,'FLOWobj'));
	addRequired(p,'A',@(x) isa(x,'GRIDobj'));
	addRequired(p,'S',@(x) isa(x,'STREAMobj'));

	addParamValue(p,'theta_method','ref',@(x) ischar(validatestring(x,{'ref','auto'})));
	addParamValue(p,'pick_method','chi',@(x) ischar(validatestring(x,{'chi','stream'})));
	addParamValue(p,'ref_concavity',0.45,@(x) isscalar(x) && isnumeric(x));
	addParamValue(p,'refit_streams',false,@(x) isscalar(x) && islogical(x));

	parse(p,DEM,FD,A,S,varargin{:});
	DEM=p.Results.DEM;
	FD=p.Results.FD;
	S=p.Results.S;
	A=p.Results.A;

	theta_method=p.Results.theta_method;
	ref_theta=p.Results.ref_concavity;
	refit_streams=p.Results.refit_streams;
	pick_method=p.Results.pick_method;

	% Make sure channel has only one channel head
	ix=streampoi(S,'channelheads','ix');
	if numel(ix)>1
		error('The provided STREAMobj has more than one channel head, you can use "SegmentPicker" to extract a single channel');
	end

	% Create hydrologically conditioned DEM;
	DEMc=imposemin(FD,DEM,0.001);

	% Parse Switches
	if strcmp(theta_method,'ref') & strcmp(pick_method,'chi');
		method=1;
	elseif strcmp(theta_method,'auto') & strcmp(pick_method,'chi');
		method=2;
	elseif strcmp(theta_method,'ref') & strcmp(pick_method,'stream');
		method=3;
	elseif strcmp(theta_method,'auto') & strcmp(pick_method,'stream');
		method=4;
	end

	close all

	switch method
	case 1
		C=ChiCalc(S,DEMc,A,1,ref_theta);

		f1=figure(1);
		set(f1,'Units','normalized','Position',[0.05 0.1 0.45 0.8],'renderer','painters');

		clf
		subplot(2,1,2);
		hold on
		plotdz(S,DEM,'dunit','km','Color',[0.5 0.5 0.5]);
		plotdz(S,DEMc,'dunit','km','Color','k');
		xlabel('Distance from Mouth (km)')
		ylabel('Elevation (m)')
		legend('Unconditioned DEM','Conditioned DEM','location','best');
		title('Long Profile')
		hold off

		subplot(2,1,1);
		hold on
		plot(C.chi,C.elev,'-k');
		scatter(C.chi,C.elev,10,'k');
		xlabel('Chi','Color','r')
		ylabel('Elevation (m)','Color','r')
		title(['Chi - Z : \theta = ' num2str(C.mn) ' : Pick Segment to Project'],'Color','r')
		hold off

		disp('    Select bounds of stream segment you want to project')
		[cv,e]=ginput(2);

		% Sort knickpoint list and construct bounds list
		cvs=sortrows(cv);

		rc=C.chi;
		rx=C.x;
		ry=C.y;
		re=C.elev;

		lb=cvs(1);
		rb=cvs(2);

		% Clip out stream segment
		lb_chidist=sqrt(sum(bsxfun(@minus, rc, lb).^2,2));
		rb_chidist=sqrt(sum(bsxfun(@minus, rc, rb).^2,2));

		lbx=rx(lb_chidist==min(lb_chidist));
		lby=ry(lb_chidist==min(lb_chidist));

		rbx=rx(rb_chidist==min(rb_chidist));
		rby=ry(rb_chidist==min(rb_chidist));	

		lix=coord2ind(DEM,lbx,lby);
		LIX=GRIDobj(DEM);
		LIX.Z(lix)=1;
		[lixmat,X,Y]=GRIDobj2mat(LIX);
		lixmat=logical(lixmat);
		LIX=GRIDobj(X,Y,lixmat);	

		rix=coord2ind(DEM,rbx,rby);
		RIX=GRIDobj(DEM);
		RIX.Z(rix)=1;
		[rixmat,X,Y]=GRIDobj2mat(RIX);
		rixmat=logical(rixmat);
		RIX=GRIDobj(X,Y,rixmat);	

		Seg=modify(S,'downstreamto',RIX);
		Seg=modify(Seg,'upstreamto',LIX);

		Cseg=ChiCalc(Seg,DEMc,A,1,C.mn);

		lbc=rc(lb_chidist==min(lb_chidist));
		lbe=re(lb_chidist==min(lb_chidist));

		y_int=lbe-(Cseg.ks)*lbc;
		pred_el=(Cseg.ks.*C.chi)+y_int;

		subplot(2,1,1)
		hold on
		plot(C.chi,pred_el,'-r','LineWidth',2);
		hold off

		subplot(2,1,2)
		hold on
		plot(C.distance./1000,pred_el,'-r','LineWidth',2);
		hold off

		f2=figure(2);
		set(f2,'Units','normalized','Position',[0.5 0.1 0.45 0.8],'renderer','painters');
		subplot(2,1,1)
		hold on
		plot([0,max(C.chi)],[0,0],'--k');
		scatter(C.chi,pred_el-C.elev,10,'k')
		xlabel('Chi');
		ylabel('Difference between projection and true profile (m)')
		hold off

		subplot(2,1,2)
		hold on
		plot([0,max(C.distance)/1000],[0,0],'--k');
		scatter(C.distance./1000,pred_el-C.elev,10,'k')
		xlabel('Distance (km)');
		ylabel('Difference between projection and true profile (m)')
		hold off
	case 2
		C=ChiCalc(S,DEMc,A,1);

		f1=figure(1);
		set(f1,'Units','normalized','Position',[0.05 0.1 0.45 0.8],'renderer','painters');

		clf
		subplot(2,1,2);
		hold on
		plotdz(S,DEM,'dunit','km','Color',[0.5 0.5 0.5]);
		plotdz(S,DEMc,'dunit','km','Color','k');
		xlabel('Distance from Mouth (km)')
		ylabel('Elevation (m)')
		legend('Unconditioned DEM','Conditioned DEM','location','best');
		title('Long Profile')
		hold off

		subplot(2,1,1);
		hold on
		plot(C.chi,C.elev,'-k');
		scatter(C.chi,C.elev,10,'k');
		xlabel('Chi','Color','r')
		ylabel('Elevation (m)','Color','r')
		title(['Chi - Z : \theta = ' num2str(C.mn) ' : Pick Segment to Project'],'Color','r')
		hold off

		disp('    Select bounds of stream segment you want to project')
		[cv,e]=ginput(2);

		% Sort knickpoint list and construct bounds list
		cvs=sortrows(cv);

		rc=C.chi;
		rx=C.x;
		ry=C.y;
		re=C.elev;

		lb=cvs(1);
		rb=cvs(2);

		% Clip out stream segment
		lb_chidist=sqrt(sum(bsxfun(@minus, rc, lb).^2,2));
		rb_chidist=sqrt(sum(bsxfun(@minus, rc, rb).^2,2));

		lbx=rx(lb_chidist==min(lb_chidist));
		lby=ry(lb_chidist==min(lb_chidist));

		rbx=rx(rb_chidist==min(rb_chidist));
		rby=ry(rb_chidist==min(rb_chidist));	

		lix=coord2ind(DEM,lbx,lby);
		LIX=GRIDobj(DEM);
		LIX.Z(lix)=1;
		[lixmat,X,Y]=GRIDobj2mat(LIX);
		lixmat=logical(lixmat);
		LIX=GRIDobj(X,Y,lixmat);	

		rix=coord2ind(DEM,rbx,rby);
		RIX=GRIDobj(DEM);
		RIX.Z(rix)=1;
		[rixmat,X,Y]=GRIDobj2mat(RIX);
		rixmat=logical(rixmat);
		RIX=GRIDobj(X,Y,rixmat);	

		Seg=modify(S,'downstreamto',RIX);
		Seg=modify(Seg,'upstreamto',LIX);

		Cseg=ChiCalc(Seg,DEMc,A,1,C.mn);
		Csegrf=ChiCalc(Seg,DEMc,A,1);

		switch refit_streams
		case false
			lbc=rc(lb_chidist==min(lb_chidist));
			lbe=re(lb_chidist==min(lb_chidist));

			y_int=lbe-(Cseg.ks)*lbc;
			pred_el=(Cseg.ks.*C.chi)+y_int;

			subplot(2,1,1)
			hold on
			plot(C.chi,pred_el,'-r','LineWidth',2);
			hold off

			subplot(2,1,2)
			hold on
			plot(C.distance./1000,pred_el,'-r','LineWidth',2);
			hold off

			f2=figure(2);
			set(f2,'Units','normalized','Position',[0.5 0.1 0.45 0.8],'renderer','painters');
			subplot(2,1,1)
			hold on
			plot([0,max(C.chi)],[0,0],'--k');
			scatter(C.chi,pred_el-C.elev,10,'k')
			xlabel('Chi');
			ylabel('Difference between projection and true profile (m)')
			hold off

			subplot(2,1,2)
			hold on
			plot([0,max(C.distance)/1000],[0,0],'--k');
			scatter(C.distance./1000,pred_el-C.elev,10,'k')
			xlabel('Distance (km)');
			ylabel('Difference between projection and true profile (m)')
			hold off
		case true
			CN=ChiCalc(S,DEMc,A,1,Csegrf.mn);

			lbnc=CN.chi(lb_chidist==min(lb_chidist));
			lbe=re(lb_chidist==min(lb_chidist));

			y_int=lbe-(Csegrf.ks)*lbnc;
			pred_el=(Csegrf.ks.*CN.chi)+y_int;

			f1=figure(1);
			set(f1,'Units','normalized','Position',[0.05 0.1 0.45 0.8],'renderer','painters');

			clf
			subplot(2,1,2);
			hold on
			plotdz(S,DEM,'dunit','km','Color',[0.5 0.5 0.5]);
			plotdz(S,DEMc,'dunit','km','Color','k');
			plot(CN.distance./1000,pred_el,'-r','LineWidth',2);
			xlabel('Distance from Mouth (km)')
			ylabel('Elevation (m)')
			legend('Unconditioned DEM','Conditioned DEM','location','best');
			title('Long Profile')
			hold off

			subplot(2,1,1);
			hold on
			plot(CN.chi,CN.elev,'-k');
			scatter(CN.chi,CN.elev,10,'k');
			plot(CN.chi,pred_el,'-r','LineWidth',2);
			xlabel('Chi')
			ylabel('Elevation (m)')
			title(['Chi - Z : \theta = ' num2str(CN.mn)],'Color','r')
			hold off

			f2=figure(2);
			set(f2,'Units','normalized','Position',[0.5 0.1 0.45 0.8],'renderer','painters');
			subplot(2,1,1)
			hold on
			plot([0,max(CN.chi)],[0,0],'--k');
			scatter(CN.chi,pred_el-CN.elev,10,'k')
			xlabel('Chi');
			ylabel('Difference between projection and true profile (m)')
			hold off

			subplot(2,1,2)
			hold on
			plot([0,max(CN.distance)/1000],[0,0],'--k');
			scatter(CN.distance./1000,pred_el-CN.elev,10,'k')
			xlabel('Distance (km)');
			ylabel('Difference between projection and true profile (m)')
			hold off
		end
	case 3
		C=ChiCalc(S,DEMc,A,1,ref_theta);

		f1=figure(1);
		set(f1,'Units','normalized','Position',[0.05 0.1 0.45 0.8],'renderer','painters');

		clf

		subplot(2,1,1);
		hold on
		plot(C.chi,C.elev,'-k');
		scatter(C.chi,C.elev,10,'k');
		xlabel('Chi')
		ylabel('Elevation (m)')
		title('Chi - Z')
		hold off

		subplot(2,1,2);
		hold on
		plotdz(S,DEM,'dunit','km','Color',[0.5 0.5 0.5]);
		plotdz(S,DEMc,'dunit','km','Color','k');
		xlabel('Distance from Mouth (km)','Color','r')
		ylabel('Elevation (m)','Color','r')
		legend('Unconditioned DEM','Conditioned DEM','location','best');
		title(['Long Profile : \theta = ' num2str(C.mn) ' : Pick Segment to Project'],'Color','r')
		hold off

		disp('    Select bounds of stream segment you want to project')
		[d,e]=ginput(2);
		d=d*1000;

		% Sort knickpoint list and construct bounds list
		ds=sortrows(d);

		rd=S.distance;
		rx=S.x;
		ry=S.y;
		% rc=flipud(C.chi);
		rc=chitransform(S,A.*(A.cellsize^2),'mn',C.mn,'a0',1);
		re=getnal(S,DEMc);

		lb=ds(1);
		rb=ds(2);

		lb_dist=sqrt(sum(bsxfun(@minus, rd, lb).^2,2));
		rb_dist=sqrt(sum(bsxfun(@minus, rd, rb).^2,2));

		lbx=rx(lb_dist==min(lb_dist));
		lby=ry(lb_dist==min(lb_dist));

		rbx=rx(rb_dist==min(rb_dist));
		rby=ry(rb_dist==min(rb_dist));	

		lix=coord2ind(DEM,lbx,lby);
		LIX=GRIDobj(DEM);
		LIX.Z(lix)=1;
		[lixmat,X,Y]=GRIDobj2mat(LIX);
		lixmat=logical(lixmat);
		LIX=GRIDobj(X,Y,lixmat);	

		rix=coord2ind(DEM,rbx,rby);
		RIX=GRIDobj(DEM);
		RIX.Z(rix)=1;
		[rixmat,X,Y]=GRIDobj2mat(RIX);
		rixmat=logical(rixmat);
		RIX=GRIDobj(X,Y,rixmat);	

		Seg=modify(S,'downstreamto',RIX);
		Seg=modify(Seg,'upstreamto',LIX);

		Cseg=ChiCalc(Seg,DEMc,A,1,C.mn);

		lbc=rc(lb_dist==min(lb_dist));
		lbe=re(lb_dist==min(lb_dist));

		y_int=lbe-(Cseg.ks)*lbc;
		pred_el=(Cseg.ks.*C.chi)+y_int;

		subplot(2,1,1)
		hold on
		plot(C.chi,pred_el,'-r','LineWidth',2);
		hold off

		subplot(2,1,2)
		hold on
		plot(C.distance./1000,pred_el,'-r','LineWidth',2);
		hold off

		f2=figure(2);
		set(f2,'Units','normalized','Position',[0.5 0.1 0.45 0.8],'renderer','painters');
		subplot(2,1,1)
		hold on
		plot([0,max(C.chi)],[0,0],'--k');
		scatter(C.chi,pred_el-C.elev,10,'k')
		xlabel('Chi');
		ylabel('Difference between projection and true profile (m)')
		hold off

		subplot(2,1,2)
		hold on
		plot([0,max(C.distance)/1000],[0,0],'--k');
		scatter(C.distance./1000,pred_el-C.elev,10,'k')
		xlabel('Distance (km)');
		ylabel('Difference between projection and true profile (m)')
		hold off
	case 4	
		C=ChiCalc(S,DEMc,A,1);

		f1=figure(1);
		set(f1,'Units','normalized','Position',[0.05 0.1 0.45 0.8],'renderer','painters');

		clf

		subplot(2,1,1);
		hold on
		plot(C.chi,C.elev,'-k');
		scatter(C.chi,C.elev,10,'k');
		xlabel('Chi')
		ylabel('Elevation (m)')
		title('Chi - Z')
		hold off

		subplot(2,1,2);
		hold on
		plotdz(S,DEM,'dunit','km','Color',[0.5 0.5 0.5]);
		plotdz(S,DEMc,'dunit','km','Color','k');
		xlabel('Distance from Mouth (km)','Color','r')
		ylabel('Elevation (m)','Color','r')
		legend('Unconditioned DEM','Conditioned DEM','location','best');
		title(['Long Profile : \theta = ' num2str(C.mn) ' : Pick Segment to Project'],'Color','r')
		hold off

		disp('    Select bounds of stream segment you want to project')
		[d,e]=ginput(2);
		d=d*1000;

		% Sort knickpoint list and construct bounds list
		ds=sortrows(d);

		rd=S.distance;
		rx=S.x;
		ry=S.y;
		% rc=flipud(C.chi);
		rc=chitransform(S,A.*(A.cellsize^2),'mn',C.mn,'a0',1);
		re=getnal(S,DEMc);

		lb=ds(1);
		rb=ds(2);

		lb_dist=sqrt(sum(bsxfun(@minus, rd, lb).^2,2));
		rb_dist=sqrt(sum(bsxfun(@minus, rd, rb).^2,2));

		lbx=rx(lb_dist==min(lb_dist));
		lby=ry(lb_dist==min(lb_dist));

		rbx=rx(rb_dist==min(rb_dist));
		rby=ry(rb_dist==min(rb_dist));	

		lix=coord2ind(DEM,lbx,lby);
		LIX=GRIDobj(DEM);
		LIX.Z(lix)=1;
		[lixmat,X,Y]=GRIDobj2mat(LIX);
		lixmat=logical(lixmat);
		LIX=GRIDobj(X,Y,lixmat);	

		rix=coord2ind(DEM,rbx,rby);
		RIX=GRIDobj(DEM);
		RIX.Z(rix)=1;
		[rixmat,X,Y]=GRIDobj2mat(RIX);
		rixmat=logical(rixmat);
		RIX=GRIDobj(X,Y,rixmat);	

		Seg=modify(S,'downstreamto',RIX);
		Seg=modify(Seg,'upstreamto',LIX);

		Cseg=ChiCalc(Seg,DEMc,A,1,C.mn);
		Csegrf=ChiCalc(Seg,DEMc,A,1);
		switch refit_streams
		case false
			lbc=rc(lb_dist==min(lb_dist));
			lbe=re(lb_dist==min(lb_dist));

			y_int=lbe-(Cseg.ks)*lbc;
			pred_el=(Cseg.ks.*C.chi)+y_int;

			subplot(2,1,1)
			hold on
			plot(C.chi,pred_el,'-r','LineWidth',2);
			hold off

			subplot(2,1,2)
			hold on
			plot(C.distance./1000,pred_el,'-r','LineWidth',2);
			hold off

			f2=figure(2);
			set(f2,'Units','normalized','Position',[0.5 0.1 0.45 0.8],'renderer','painters');
			subplot(2,1,1)
			hold on
			plot([0,max(C.chi)],[0,0],'--k');
			scatter(C.chi,pred_el-C.elev,10,'k')
			xlabel('Chi');
			ylabel('Difference between projection and true profile (m)')
			hold off

			subplot(2,1,2)
			hold on
			plot([0,max(C.distance)/1000],[0,0],'--k');
			scatter(C.distance./1000,pred_el-C.elev,10,'k')
			xlabel('Distance (km)');
			ylabel('Difference between projection and true profile (m)')
			hold off

		case true
			CN=ChiCalc(S,DEMc,A,1,Csegrf.mn);
			rcn=chitransform(S,A.*(A.cellsize^2),'mn',Csegrf.mn,'a0',1);

			lbnc=rcn(lb_dist==min(lb_dist));
			lbe=re(lb_dist==min(lb_dist));

			y_int=lbe-(Csegrf.ks)*lbnc;
			pred_el=(Csegrf.ks.*CN.chi)+y_int;

			f1=figure(1);
			set(f1,'Units','normalized','Position',[0.05 0.1 0.45 0.8],'renderer','painters');

			clf
			subplot(2,1,2);
			hold on
			plotdz(S,DEM,'dunit','km','Color',[0.5 0.5 0.5]);
			plotdz(S,DEMc,'dunit','km','Color','k');
			plot(CN.distance./1000,pred_el,'-r','LineWidth',2);
			xlabel('Distance from Mouth (km)')
			ylabel('Elevation (m)')
			legend('Unconditioned DEM','Conditioned DEM','location','best');
			title(['Long Profile : \theta = ' num2str(CN.mn)],'Color','r')
			hold off

			subplot(2,1,1);
			hold on
			plot(CN.chi,CN.elev,'-k');
			scatter(CN.chi,CN.elev,10,'k');
			plot(CN.chi,pred_el,'-r','LineWidth',2);
			xlabel('Chi')
			ylabel('Elevation (m)')
			title('Chi - Z ')
			hold off

			f2=figure(2);
			set(f2,'Units','normalized','Position',[0.5 0.1 0.45 0.8],'renderer','painters');
			subplot(2,1,1)
			hold on
			plot([0,max(CN.chi)],[0,0],'--k');
			scatter(CN.chi,pred_el-CN.elev,10,'k')
			xlabel('Chi');
			ylabel('Difference between projection and true profile (m)')
			hold off

			subplot(2,1,2)
			hold on
			plot([0,max(CN.distance)/1000],[0,0],'--k');
			scatter(CN.distance./1000,pred_el-CN.elev,10,'k')
			xlabel('Distance (km)');
			ylabel('Difference between projection and true profile (m)')
			hold off
		end
	end
end

function [OUT]=ChiCalc(S,DEM,A,a0,varargin)
% Modified version of chiplot function by Wolfgang Schwanghart to remove unused options/outputs and 
% to interpolate chi-z relationship so that chi values are equally spaced to avoid biasing of ksn fit
% by clustering at high drainage areas 

betamethod='ls'; % Default is to find best fit slope via least squares, can
				 % change to 'lad' to calculate via least absolute deviations
mnmethod='ls'; % Similar switch for calculating best fit concavity (if reference not provided)

% nr of nodes in the entire stream network
nrc = numel(S.x);
M   = sparse(double(S.ix),double(S.ixc),true,nrc,nrc);
% find outlet
outlet = sum(M,2) == 0 & sum(M,1)'~=0;
if nnz(outlet)>1
    % there must not be more than one outlet (constraint could be removed
    % in the future).
    error('The stream network must not have more than one outlet');
end


% elevation values at nodes
zx   = double(DEM.Z(S.IXgrid));
% elevation at outlet
zb   = double(DEM.Z(S.IXgrid(outlet)));
% a is the term inside the brackets of equation 6b 
a    = double(a0./(A.Z(S.IXgrid)*(A.cellsize.^2)));
% x is the cumulative horizontal distance in upstream direction
x    = S.distance;
Lib = true(size(x));

% Find bestfit concavity if no concavity provided or use provided concavity
if isempty(varargin)
    mn0  = 0.5; % initial value
    mn   = fminsearch(@mnfit,mn0);
else
	mn=varargin{1};
end

% calculate chi
chi = netcumtrapz(x,a.^mn,S.ix,S.ixc);

% Resample chi-elevation relationship using cubic spline interpolation
chiF=chi(Lib);
zabsF=zx(Lib)-zb;

chiS=linspace(0,max(chiF),numel(chiF)).';
zS=spline(chiF,zabsF,chiS);

% Calculate beta
switch betamethod
    case 'ls'
        % least squares
        BETA = chiS\(zS);
    case 'lad'
        % least absolute deviations
        BETA = fminsearch(@(b) sum(abs(b*chiS - (zS))),0.0334);
end

OUT=struct;
OUT.ks   = BETA*a0^mn;
OUT.mn   = mn;
[OUT.x,...
 OUT.y,...
 OUT.chi,...
 OUT.elev,...
 OUT.elevbl,...
 OUT.distance,...
 OUT.pred,...
 OUT.area] = STREAMobj2XY(S,chi,DEM,zx-zb,S.distance,BETA*chi,A.*(A.cellsize^2));
 OUT.res = OUT.elevbl - OUT.pred;

	% Nested function for calculating mn ratio
	function sqres = mnfit(mn)
		% calculate chi with a given mn ratio
		% and integrate in upstream direction
		CHI = netcumtrapz(x(Lib),a(Lib).^mn,S.ix,S.ixc);%*ab.^mn
		% normalize both variables
		CHI = CHI ./ max(CHI);
		z   = zx(Lib)-zb;
		z   = z./max(z);
        switch mnmethod
            case 'ls' % least squares
                sqres = sum((CHI - z).^2);
            case 'lad' % least absolute devation
                sqres = sum(sqrt(abs(CHI-z)));
        end
	end

end

function z = netcumtrapz(x,y,ix,ixc)
% cumtrapz along upward direction in a directed tree network

z = zeros(size(x));
for lp = numel(ix):-1:1;
    z(ix(lp)) = z(ixc(lp)) + (y(ixc(lp))+(y(ix(lp))-y(ixc(lp)))/2) *(abs(x(ixc(lp))-x(ix(lp))));
end
end
