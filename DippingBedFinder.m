function [BED]=DippingBedFinder(dem_location,x_coord,y_coord,hght_abv_base,thickness,strike,dip);
% Function to determine the expected location of a planar dipping bed within a landscape based on an input coordinate
%
% Required Inputs:
% 	dem_location - file path of the DEM file, must be either a geotiff or ascii text file and be in a projected coordinate 
%					system with horizontal and vertical units in meters
% 	x_coord - x (longitude/easting) coordinate of the location from which you wish to project bed locations, must be in the 
%				same coordinate system as the dem
% 	y_coord - y (latitude/northing) coordinate of the location from which you with to project bed locations
% 	hght_abv_base - height of the outcrop of interest above the base of the bed of interest (i.e. positin of the outcrop in the section)
%   thickness - thickness of the bed (hght_abv_base must be smaller than total thickness)
% 	strike - strike of bed, report with right hand rule
% 	dip - dip of bed
%
% Output:
% 	Code will produce a figure showing expected location of bed and will also output a binary ascii text file. This ascii
% 	is importable into Arc or another GIS program. It will be in the same projection and the exact same dimensions as your input DEM
% 	and will have a value of 0 wherever the bed does not outcrop and a value of 1 where the bed does outcrop.
%
% Example Input:
% 	DippingBedFinder('C:\Users\Adam\Project\location_dem.tif',45325.23,1024567.2,10,50,270,40);
%	For a bed striking east-west and dipping to the north at 40 degrees, that is 50 meters thick and the outcrop is near the bottom of the bed
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function Written by Adam M. Forte - Last Revised Fall 2016 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert dem to GRIDobj
DEM=GRIDobj(dem_location);

% Get info regarding DEM dimensions
gs=DEM.cellsize;
[xvec,yvec]=getcoordinates(DEM);
dim=DEM.size;

% Find row and column position and elevation of input coordinate
x_dif=abs(xvec-x_coord);
y_dif=abs(yvec-y_coord);

[~,col]=min(x_dif);
[~,row]=min(y_dif);

sam_el=DEM.Z(row,col);	

% Generate Grid
xmax=dim(2)*gs;
ymax=dim(1)*gs;

x_spacing=(xmax)/(gs);
y_spacing=(ymax)/(gs);

x=linspace(0,xmax,x_spacing);
y=linspace(0,ymax,y_spacing);

[xi,yi]=meshgrid(x,y);

% Calculate plane geometry
de=PlaneOrient(strike,dip,xmax,ymax,xi,yi);

% Recalculate plane elevations based on coordinate
orig_plane_el=de(row,col);
el_dif=sam_el-orig_plane_el;
bed_el=de+el_dif;

% Determine thickness in bed
theta=90-dip;

if hght_abv_base>thickness
    error('Height above base cannot be greater than the total bed thickness')
end

tu=thickness-hght_abv_base;
tb=thickness-tu;
atu=tu/sin(deg2rad(theta));
atb=tb/sin(deg2rad(theta));

% Determine locations of intersections
up_int=bed_el+atu;
lb_int=bed_el-atb;

z=DEM.Z;
int=z<=up_int & z>=lb_int;
int=double(int);

BED=GRIDobj(xvec,yvec,int);

% Plot
f1=figure(1);
set(f1,'Units','normalized','Position',[0.1 0.1 0.8 0.8],'renderer','painters');
hold on
imageschs(DEM,BED);
scatter(x_coord,y_coord,20,'w','filled');
title('Projected Intersection of Bed');
hold off

% Output
GRIDobj2ascii(BED,'Bed_Location.txt');


end

function [de]=PlaneOrient(strike,dip,xmax,ymax,xi,yi)
   strike=-strike; 

  % Check for vertical beds
  if dip==90
    error('Beds cannot be vertical');
  end

  if strike>=0 & strike<90;
    BetaX=90-strike;
    BetaY=strike;

    AlphaX=atan(sin(deg2rad(BetaX))*tan(deg2rad(dip)));
    AlphaY=atan(sin(deg2rad(BetaY))*tan(deg2rad(dip)));

    % Calc normal vector
    XVec=[xmax,0,-xmax*tan(AlphaX)];
    YVec=[0,-ymax,-ymax*tan(AlphaY)];
    normal=cross(XVec,YVec);

    a=normal(:,1); b=normal(:,2); c=normal(:,3);
    d=-1*(a*0+b*ymax+c*0);

    de=(a.*xi + b.*yi+d)./-c;

  elseif strike>=90 & strike<180;
    BetaX=strike-90;
    BetaY=180-strike;

    AlphaX=atan(sin(deg2rad(BetaX))*tan(deg2rad(dip)));
    AlphaY=atan(sin(deg2rad(BetaY))*tan(deg2rad(dip)));

    % Calc normal vector
    XVec=[-xmax,0,-xmax*tan(AlphaX)];
    YVec=[0,-ymax,-ymax*tan(AlphaY)];
    normal=cross(XVec,YVec);

    a=normal(:,1); b=normal(:,2); c=normal(:,3);
    d=-1*(a*xmax+b*ymax+c*0);

    de=(a.*xi + b.*yi+d)./-c;
  elseif strike>=180 & strike<270;
    BetaX=270-strike;
    BetaY=strike-180;

    AlphaX=atan(sin(deg2rad(BetaX))*tan(deg2rad(dip)));
    AlphaY=atan(sin(deg2rad(BetaY))*tan(deg2rad(dip)));

    % Calc normal vector
    XVec=[-xmax,0,-xmax*tan(AlphaX)];
    YVec=[0,ymax,-ymax*tan(AlphaY)];
    normal=cross(XVec,YVec);

    a=normal(:,1); b=normal(:,2); c=normal(:,3);
    d=-1*(a*xmax+b*0+c*0);

    de=(a.*xi + b.*yi+d)./-c;
  else 
    BetaX=strike-270;
    BetaY=360-strike;

    AlphaX=atan(sin(deg2rad(BetaX))*tan(deg2rad(dip)));
    AlphaY=atan(sin(deg2rad(BetaY))*tan(deg2rad(dip)));

    % Calc normal vector
    XVec=[xmax,0,-xmax*tan(AlphaX)];
    YVec=[0,ymax,-ymax*tan(AlphaY)];
    normal=cross(XVec,YVec);

    a=normal(:,1); b=normal(:,2); c=normal(:,3);
    d=-1*(a*0+b*0+c*0);

    de=(a.*xi + b.*yi+d)./-c;
  end


end 