function [BED]=DippingBedFinder(DEM,xy,hght_abv_base,thickness,strike,dip);
    %
    % Usage:
    %   [GRIDobj]=DippingBedFinder(DEM,xy,hght_abv_base,thickness,strike,dip); 
    %
    % Description:
    %   Function to determine the expected location of a planar dipping bed within a landscape based on an input coordinate
    %
    % Required Inputs:
    % 	DEM - DEM GRIDobj
    %   xy - 1 x 2 vector with the x and y coordinate (i.e. easting and northing) of the location of interest, if you provide an empty vector
    %        you will be given the opportunity to pick a location on the DEM
    % 	hght_abv_base - height of the outcrop of interest above the base of the bed of interest (i.e. positin of the outcrop in the section)
    %   thickness - thickness of the bed (hght_abv_base must be smaller than total thickness)
    % 	strike - strike of bed, report with right hand rule
    % 	dip - dip of bed
    %
    % Output:
    % 	Code will produce a figure showing expected location of bed and will produce a binary GRIDobj with expected location of the bed 
    %   (1 where the bed should appear, 0 where it should not)
    %
    % Example Input:
    %   To make a DEM GRIDobj:
    %   [DEM]=GRIDobj('name_of_tif.tif'); or [DEM]=GRIDobj('name_of_ascii.txt');
    %   then
    % 	[BED]=DippingBedFinder(DEM,[45325.23 1024567.2],10,50,270,40); If you want to provide a coordinate
    %   [BED]=DippingBedFinder(DEM,[],10,50,270,40); If you want to pick a spot on the DEM
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Function Written by Adam M. Forte - Updated : 06/18/18 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if isempty(xy);
        f1=figure(1);
        set(f1,'Units','normalized','Position',[0.1 0.1 0.8 0.8],'renderer','painters');
        imageschs(DEM,DEM);
        [x_coord,y_coord]=ginput(1);
        close(f1)
    else
        x_coord=xy(1);
        y_coord=xy(2);
    end


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

    % Calculate general plane geometry
    de=PlaneEq(strike,dip,xi,yi);

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
    atu=tu/sind(theta);
    atb=tb/sind(theta);

    % Determine locations of intersections
    up_int=bed_el+atu;
    lb_int=bed_el-atb;

    z=DEM.Z;
    inter=z<=up_int & z>=lb_int;
    inter=double(inter);

    BED=GRIDobj(xvec,yvec,inter);
    BED.georef=DEM.georef;

    % Plot
    f1=figure(1);
    set(f1,'Units','normalized','Position',[0.1 0.1 0.8 0.8],'renderer','painters');
    hold on
    imageschs(DEM,BED);
    scatter(x_coord,y_coord,20,'w','filled');
    title('Projected Intersection of Bed');
    hold off

end


function [de]=PlaneEq(strike,dip,xi,yi);
    % Equation for normal vector comes from Charles Ammon:
    % http://eqseis.geosc.psu.edu/~cammon/HTML/UsingMATLAB/PDF/ML1%20FaultNormals.pdf
    %
    % Calculate normal vector
    a=sind(dip)*cosd(strike);
    b=-sind(dip)*sind(strike);
    c=-cosd(dip);

    % Convert to equation for plane centered on origin
    de=(a.*xi + b.*yi)./-c;

    % Correct for coordinate system flip
    de=fliplr(de);
end

