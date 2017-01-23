function [Outlets]=DetritalSamplePicker(DEM,FD,A,S,varargin)
    % Function takes results of makes streams and allows for interactive picking of basins for detrital
    % analyses (e.g. Be-10 cosmo). Displays two panel figure with topography colored by elevation and local relief
    % on which to pick individual basins. After the figure displays, it will wait until you press enter to begin the sample
    % picking process. This is to allow you to zoom, pan, etc to find a stream you are interested in. When you click enter, 
    % cross hairs will appear in the elevation map so you can select a pour point. Once you select a pour point, a new figure 
    % will display this basin and stream to confirm that's the watershed you wanted (it will also display the drainage area). If 
    % You can either accept this basin or reject it if it was misclick. If you accept it will then display a new figure with the
    % chi-z and longitudinal profiles for that basin. It will then give you a choice to either save the choice or discard it.
    % Finally it will ask if you want to keep picking streams, if you choose yes (the default) it will start the process over.
    % Note that any selected (and saved) pour point will be displayed as a white circle on the main figure. As you pick basins
    % the funciton saves a file called 'Outlets.mat' that contains the outlets you've picked so far. If you exit out of the function
    % and restart it later, it looks for this Oulets file so you can pick up where you left off.
    %
    %
    % Required Inputs:
    %       DEM - GRIDobj of the DEM
    %       FD - FLOWobj from the supplied DEM
    %       A - Flow accumulation grid (GRIDobj)
    %       S - STREAMobj derived from the DEM
    %
    % Optional Inputs:
    %       theta_ref [0.45]- reference concavity for chi-Z plots
    %       rlf_radius [2500] - radius in map units for calculating local relief
    %       plot_type ['native'] - expects either 'native' or 'downsample', default is 'native'. Controls whether all streams are drawn as individual lines ('native') or if
    %           the stream network is plotted as a grid and downsampled ('downsample'). The 'downsample' option is much faster for large datasets, 
    %           but can result in inaccurate site selections. The 'native' option is easier to see, but can be very slow to load and interact with.
    %       threshold_area [1e6] - used to redraw downsampled stream network if 'plot_type' is set to 'downsample'
    % 
    % Outputs:
    %       Outlets - n x 3 matrix of sample locations with columns basin number, x coordinate, and y coordinate
    %
    % Examples:
    %       [Outs]=DetritalSamplePicker(DEM,FD,A,S);
    %       [Outs]=DetritalSamplePicker(DEM,FD,A,S,'rlf_radius',5000);
    %  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Function Written by Adam M. Forte - Last Revised Winter 2017 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Parse Inputs
    p = inputParser;
    p.FunctionName = 'DetritalSamplePicker';
    addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
    addRequired(p,'FD',@(x) isa(x,'FLOWobj'));
    addRequired(p,'A',@(x) isa(x,'GRIDobj'));
    addRequired(p,'S',@(x) isa(x,'STREAMobj'));

    addParamValue(p,'theta_ref',0.45,@(x) isscalar(x) && isnumeric(x));
    addParamValue(p,'rlf_radius',2500,@(x) isscalar(x) && isnumeric(x));
    addParamValue(p,'plot_type','native',@(x) ischar(validatestring(x,{'native','downsample'})));
    addParamValue(p,'threshold_area',1e6,@(x) isnumeric(x));   

    parse(p,DEM,FD,A,S,varargin{:});
    DEM=p.Results.DEM;
    FD=p.Results.FD;
    S=p.Results.S;
    A=p.Results.A;

    theta_ref=p.Results.theta_ref;
    rlf_radius=p.Results.rlf_radius;
    plot_type=p.Results.plot_type;
    threshold_area=p.Results.threshold_area;

    % Check for outlets file from previous run
    if exist('Outlets.mat','file')==2;
        load('Outlets.mat','Outlets');
        ii=max(Outlets(:,3))+1;
    else 
       ii=1; 
       Outlets=zeros(0,3);
    end
             
    % Hydrologically condition dem
    zc=mincosthydrocon(S,DEM,'interp',0.1);
    DEMf=GRIDobj(DEM);
    DEMf.Z(DEMf.Z==0)=NaN;
    DEMf.Z(S.IXgrid)=zc;

    switch plot_type
    case 'native'
        disp('Calculating local relief')
        RLF=localtopography(DEM,rlf_radius);
          
        % Plot main figure
        f1=figure(1);
        clf
        set(f1,'Units','normalized','Position',[0.05 0.1 0.45 0.8]);
        
        ax(2)=subplot(2,1,2);
        hold on
        imageschs(DEM,RLF);
        plot(S,'-k','LineWidth',1.5);
        scatter(Outlets(:,1),Outlets(:,2),20,'w','filled');
        title('Local Relief')
        hold off
            
        ax(1)=subplot(2,1,1);
        hold on
        imageschs(DEM,DEM);
        plot(S,'-k','LineWidth',1.5);
        scatter(Outlets(:,1),Outlets(:,2),20,'w','filled');
        title('Elevation')
        hold off
        
        linkaxes(ax,'xy');
    case 'downsample'
        disp('Downsampling datasets for display purposes')  
        % Redo flow direction   
        DEMr=resample(DEM,DEM.cellsize*4);
        FDr=FLOWobj(DEMr,'preprocess','carve');
        % True outlets
        out_T_xy=streampoi(S,'outlets','xy');
        % Downsampled total stream network
        Sr_temp=STREAMobj(FDr,'minarea',threshold_area,'unit','mapunits');
        out_D_xy=streampoi(Sr_temp,'outlets','xy');
        out_D_ix=streampoi(Sr_temp,'outlets','ix');
        % Find if outlet list is different
        dists=pdist2(out_T_xy,out_D_xy);
        [~,s_out_ix]=min(dists,[],2);
        out_D_ix=out_D_ix(s_out_ix);
        % Rebuild downsampled network
        Sr=STREAMobj(FDr,'minarea',threshold_area,'unit','mapunits','outlets',out_D_ix);
        % Turn it into a grid
        SG=STREAMobj2GRIDobj(Sr);
        % Redo local relief
        RLFr=localtopography(DEMr,rlf_radius);

        % Plot main figure
        f1=figure(1);
        clf
        set(f1,'Units','normalized','Position',[0.05 0.1 0.45 0.8]);
        
        ax(2)=subplot(2,1,2);
        hold on
        imageschs(DEMr,RLFr);
        plot(Sr,'-k','LineWidth',1.5);
        scatter(Outlets(:,1),Outlets(:,2),20,'w','filled');
        title('Local Relief')
        hold off
            
        ax(1)=subplot(2,1,1);
        hold on
        imageschs(DEMr,DEMr);
        plot(Sr,'-k','LineWidth',1.5);
        scatter(Outlets(:,1),Outlets(:,2),20,'w','filled');
        title('Elevation')
        hold off
        
        linkaxes(ax,'xy');
    end
    
    % Start sample selection process
    str1='N';
    str2='Y';

    while strcmp(str2,'Y')==1;
        while strcmp(str1,'N')==1;	
            
            disp('Zoom and pan to area of interest and press "return/enter" when ready to pick');
            pause()

    
            disp('    Choose sample point on elevation map')
            [x,y]=ginput(1);

            % Build logical raster
            [xn,yn]=snap2stream(S,x,y);
            ix=coord2ind(DEM,xn,yn);
            IX=GRIDobj(DEM);
            IX.Z(ix)=1;
            [ixmat,X,Y]=GRIDobj2mat(IX);
            ixmat=logical(ixmat);
            IX=GRIDobj(X,Y,ixmat);

            % Clip out stream network
            Sn=modify(S,'upstreamto',IX);
            
            % Clip out DEM
            I=dependencemap(FD,xn,yn);
            DEMc=crop(DEM,I,nan);
                        
            % Calculate drainage area
            dep_map=GRIDobj2mat(I);
            num_pix=sum(sum(dep_map));
            drainage_area=(num_pix*DEMc.cellsize*DEMc.cellsize)/(1e6);
        
            f2=figure(2);
            clf
            set(f2,'Units','normalized','Position',[0.5 0.1 0.45 0.45])
            hold on
            title(['Drainage Area = ' num2str(round(drainage_area)) 'km^2']);
            imagesc(DEMc)
            plot(Sn,'-r','LineWidth',2);
            scatter(xn,yn,20,'w','filled');
            hold off

            prompt='    Is this the stream segment you wanted? Y/N [Y]: ';
            str5=input(prompt,'s');
            if isempty(str5) | strcmp(str5,'Y')==1;
                str1 = 'Y';
            elseif strcmp(str5,'N')
                str1='N';
            end
            close figure 2
        end
        
        Sn=klargestconncomps(Sn,1);
        C=chiplot(Sn,DEM,A,'a0',1,'mn',theta_ref,'plot',false);

        f2=figure(2);
        clf
        set(f2,'Units','normalized','Position',[0.5 0.1 0.45 0.8],'renderer','painters');
        subplot(2,1,1);
        hold on
        plot(C.chi,C.elev);
        xlabel('Chi')
        ylabel('Elevation (m)')
        title('Chi - Z')
        hold off

        subplot(2,1,2);
        hold on
        plotdz(Sn,DEMf,'dunit','km','Color','k');
        xlabel('Distance from Mouth (km)')
        ylabel('Elevation (m)')
        title('Long Profile')
        hold off
       
        prompt='    Keep this stream? Y/N [Y]: ';
        str2=input(prompt,'s');
        if isempty(str2) | strcmp(str2,'Y')==1;
            Outlets(ii,3)=ii;
            Outlets(ii,1)=xn;
            Outlets(ii,2)=yn;
            ii=ii+1;
            
            % Plot selected point on figures
            figure(1);
            subplot(2,1,2);
            hold on
            scatter(xn,yn,20,'w','filled');
            hold off
                      
            subplot(2,1,1);
            hold on
            scatter(xn,yn,20,'w','filled');
            hold off
            
            save('Outlets.mat','Outlets');
        end
        
        prompt='    Keep picking streams? Y/N [Y]: ';
        str3=input(prompt,'s');
        if isempty(str3) | strcmp(str3,'Y')==1;
           str2 = 'Y';
           str1 = 'N';
        end
        
        close figure 2;

    end
    
    
    
