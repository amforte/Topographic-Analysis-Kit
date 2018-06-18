function [Outlets]=BasinPicker(DEM,FD,A,S,varargin)
    %
    % Usage:
    %   [Outlets]=BasinPicker(DEM,FD,A,S);
    %   [Outlets]=BasinPicker(DEM,FD,A,S,'name',value,...);
    %
    % Description:
    %   Function takes results of makes streams and allows for interactive picking of basins (watersheds). Function was
    %   designed intially for choosing basins suitable for detrital analyses (e.g. Be-10 cosmo). Displays two panel figure 
    %   with topography colored by elevation and local relief on which to pick individual basins. After the figure displays,
    %   it will wait until you press enter to begin the watershed picking process. This is to allow you to zoom, pan, etc to 
    %   find a stream you are interested in. When you click enter, cross hairs will appear in the elevation map so you can 
    %   select a pour point. Once you select a pour point, a new figure will display this basin and stream to confirm that's 
    %   the watershed you wanted (it will also display the drainage area). You can either accept this basin or reject it if 
    %   it was misclick. If you accept it will then display a new figure with the chi-z and longitudinal profiles for that basin. 
    %   It will then give you a choice to either save the choice or discard it. Finally it will ask if you want to keep picking 
    %   streams, if you choose yes (the default) it will start the process over. Note that any selected (and saved) pour point 
    %   will be displayed as a white circle on the main figure. As you pick basins the funciton saves a file called 'Outlets.mat'
    %   that contains the outlets you've picked so far. If you exit out of the function and restart it later, it looks for this 
    %   Outlets file in the current working directory so you can pick up where you left off.
    %
    %
    % Required Inputs:
    %       DEM - GRIDobj of the DEM
    %       FD - FLOWobj from the supplied DEM
    %       A - Flow accumulation grid (GRIDobj)
    %       S - STREAMobj derived from the DEM
    %
    % Optional Inputs:
    %       ref_concavity [0.50]- reference concavity for chi-Z plots
    %       rlf_radius [2500] - radius in map units for calculating local relief OR
    %       rlf_grid [] - if you already have a local relief grid that you've calculated, you can provide it with 'rlf_grid', it must be a GRIDobj and must be the same
    %           coordinates and dimensions as the provided DEM.
    %       extra_grid [] - sometimes it can be useful to also view an additional grid (e.g. georeferenced road map, precipitation grid, etc) along with the DEM and relief.
    %           This grid can be a different size or have a different cellsize than the underlying dem (but still must be the same projection and coordinates system!), it will be
    %           resampled to match the provided DEM. 
    %       cmap ['jet'] - colormap to use for the displayed maps. Input can be the name of a standard colormap or a nx3 array of rgb values
    %           to use as a colormap.  
    %       conditioned_DEM [] - option to provide a hydrologically conditioned DEM for use in this function (do not provide a conditoned DEM
    %           for the main required DEM input!) which will be used for extracting elevations. See 'ConditionDEM' function for options for making a 
    %           hydrological conditioned DEM. If no input is provided the code defaults to using the mincosthydrocon function.
    %       interp_value [0.1] - value (between 0 and 1) used for interpolation parameter in mincosthydrocon (not used if user provides a conditioned DEM)
    %       plot_type ['vector'] - expects either 'vector' or 'grid', default is 'vector'. Controls whether all streams are drawn as individual lines ('vector') or if
    %           the stream network is plotted as a grid and downsampled ('grid'). The 'grid' option is much faster for large datasets, 
    %           but can result in inaccurate site selections. The 'vector' option is easier to see, but can be very slow to load and interact with.
    %       threshold_area [1e6] - used to redraw downsampled stream network if 'plot_type' is set to 'grid'
    %
    % 
    % Outputs:
    %       Outlets - n x 3 matrix of sample locations with columns basin number, x coordinate, and y coordinate (valid input to 'ProcessRiverBasins'
    %           as 'river_mouths' parameter)
    %
    % Examples:
    %       [Outs]=DetritalSamplePicker(DEM,FD,A,S);
    %       [Outs]=DetritalSamplePicker(DEM,FD,A,S,'rlf_radius',5000);
    %       [Outs]=DetritalSamplePicker(DEM,FD,A,S,,'rlf_grid',RLF);
    %  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Function Written by Adam M. Forte - Updated : 06/18/18 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Parse Inputs
    p = inputParser;
    p.FunctionName = 'DetritalSamplePicker';
    addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
    addRequired(p,'FD',@(x) isa(x,'FLOWobj'));
    addRequired(p,'A',@(x) isa(x,'GRIDobj'));
    addRequired(p,'S',@(x) isa(x,'STREAMobj'));

    addParameter(p,'ref_concavity',0.5,@(x) isscalar(x) && isnumeric(x));
    addParameter(p,'rlf_radius',2500,@(x) isscalar(x) && isnumeric(x));
    addParameter(p,'plot_type','vector',@(x) ischar(validatestring(x,{'vector','grid'})));
    addParameter(p,'cmap','jet',@(x) ischar(x) || isnumeric(x) & size(x,2)==3);
    addParameter(p,'threshold_area',1e6,@(x) isnumeric(x));   
    addParameter(p,'rlf_grid',[],@(x) isa(x,'GRIDobj'));
    addParameter(p,'extra_grid',[],@(x) isa(x,'GRIDobj'));
    addParameter(p,'conditioned_DEM',[],@(x) isa(x,'GRIDobj'));
    addParameter(p,'interp_value',0.1,@(x) isnumeric(x) && x>=0 && x<=1);

    parse(p,DEM,FD,A,S,varargin{:});
    DEM=p.Results.DEM;
    FD=p.Results.FD;
    S=p.Results.S;
    A=p.Results.A;

    theta_ref=p.Results.ref_concavity;
    rlf_radius=p.Results.rlf_radius;
    plot_type=p.Results.plot_type;
    threshold_area=p.Results.threshold_area;
    cmap=p.Results.cmap;
    RLF=p.Results.rlf_grid;
    EG=p.Results.extra_grid;
    DEMf=p.Results.conditioned_DEM;
    iv=p.Results.interp_value;


    % Check for outlets file from previous run
    if exist('Outlets.mat','file')==2;
        load('Outlets.mat','Outlets');
        ii=max(Outlets(:,3))+1;
    else 
       ii=1; 
       Outlets=zeros(0,3);
    end
             
    % Hydrologically condition dem
    if isempty(DEMf)
        zc=mincosthydrocon(S,DEM,'interp',iv);
        DEMf=GRIDobj(DEM);
        DEMf.Z(DEMf.Z==0)=NaN;
        DEMf.Z(S.IXgrid)=zc;
    end

    switch plot_type
    case 'vector'


        if isempty(RLF)
            disp('Calculating local relief')
            RLF=localtopography(DEM,rlf_radius);
        end

        % Plot main figure
        if isempty(EG)
            f1=figure(1);
            clf
            set(f1,'Units','normalized','Position',[0.05 0.1 0.45 0.8]);        

            ax(2)=subplot(2,1,2);
            hold on
            imageschs(DEM,RLF,'colormap',cmap);
            plot(S,'-k','LineWidth',1.5);
            scatter(Outlets(:,1),Outlets(:,2),20,'w','filled');
            title('Local Relief')
            hold off
                
            ax(1)=subplot(2,1,1);
            hold on
            imageschs(DEM,DEM,'colormap',cmap);
            plot(S,'-k','LineWidth',1.5);
            scatter(Outlets(:,1),Outlets(:,2),20,'w','filled');
            title('Elevation')
            hold off
            
            linkaxes(ax,'xy');
        else

            if ~validatealignment(EG,DEM);
                disp('Resampling extra grid');
                EG=resample(EG,DEM,'bicubic');
            end

            % Plot main figure
            f1=figure(1);
            clf
            set(f1,'Units','normalized','Position',[0.05 0.1 0.45 0.8]);
            
            ax(3)=subplot(3,1,3);
            hold on
            imageschs(DEM,EG,'colormap',cmap)
            plot(S,'-k','LineWidth',1.5);
            scatter(Outlets(:,1),Outlets(:,2),20,'w','filled');
            title('User Provided Extra Grid')
            hold off

            ax(2)=subplot(3,1,2);
            hold on
            imageschs(DEM,RLF,'colormap',cmap);
            plot(S,'-k','LineWidth',1.5);
            scatter(Outlets(:,1),Outlets(:,2),20,'w','filled');
            title('Local Relief')
            hold off
                
            ax(1)=subplot(3,1,1);
            hold on
            imageschs(DEM,DEM,'colormap',cmap);
            plot(S,'-k','LineWidth',1.5);
            scatter(Outlets(:,1),Outlets(:,2),20,'w','filled');
            title('Elevation')
            hold off
            
            linkaxes(ax,'xy');
        end

    case 'grid'
        disp('Downsampling datasets for display purposes')  
        % Redo flow direction   
        DEMr=resample(DEM,DEM.cellsize*4,'bicubic');
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
        if isempty(RLF)
            RLFr=localtopography(DEMr,rlf_radius);
        else
            RLFr=resample(RLF,DEMr,'bicubic');
        end

        if isempty(EG)
            % Plot main figure
            f1=figure(1);
            clf
            set(f1,'Units','normalized','Position',[0.05 0.1 0.45 0.8]);
            
            ax(2)=subplot(2,1,2);
            hold on
            imageschs(DEMr,RLFr,'colormap',cmap);
            plot(Sr,'-k','LineWidth',1.5);
            scatter(Outlets(:,1),Outlets(:,2),20,'w','filled');
            title('Local Relief')
            hold off
                
            ax(1)=subplot(2,1,1);
            hold on
            imageschs(DEMr,DEMr,'colormap',cmap);
            plot(Sr,'-k','LineWidth',1.5);
            scatter(Outlets(:,1),Outlets(:,2),20,'w','filled');
            title('Elevation')
            hold off
            
            linkaxes(ax,'xy');
        else

            EGr=resample(EG,DEMr,'bicubic');

            % Plot main figure
            f1=figure(1);
            clf
            set(f1,'Units','normalized','Position',[0.05 0.1 0.45 0.8]);
            
            ax(3)=subplot(3,1,3);
            hold on
            imageschs(DEMr,EGr,'colormap',cmap)
            plot(S,'-k','LineWidth',1.5);
            scatter(Outlets(:,1),Outlets(:,2),20,'w','filled');
            title('User Provided Extra Grid')
            hold off

            ax(2)=subplot(3,1,2);
            hold on
            imageschs(DEMr,RLFr,'colormap',cmap);
            plot(Sr,'-k','LineWidth',1.5);
            scatter(Outlets(:,1),Outlets(:,2),20,'w','filled');
            title('Local Relief')
            hold off
                
            ax(1)=subplot(3,1,1);
            hold on
            imageschs(DEMr,DEMr,'colormap',cmap);
            plot(Sr,'-k','LineWidth',1.5);
            scatter(Outlets(:,1),Outlets(:,2),20,'w','filled');
            title('Elevation')
            hold off
            
            linkaxes(ax,'xy');
        end

    end
    
    % Start sample selection process
    str1='N'; % Watershed choice
    str2='Y'; % Continue choosing

    while strcmpi(str2,'Y');
        while strcmpi(str1,'N');	
            
            if isempty(EG)
                subplot(2,1,1);
                hold on
                title('Zoom and pan to area of interest and press "return/enter" when ready to pick')
                hold off
                pause()

                subplot(2,1,1);
                hold on
                title('Choose sample point on elevation map')
                hold off
            else
                subplot(3,1,1);
                hold on
                title('Zoom and pan to area of interest and press "return/enter" when ready to pick')
                hold off
                pause()

                subplot(3,1,1);
                hold on
                title('Choose sample point on elevation map')
                hold off
            end

            [x,y]=ginput(1);

            % Build logical raster
            [xn,yn]=snap2stream(S,x,y);
            ix=coord2ind(DEM,xn,yn);
            IX=GRIDobj(DEM,'logical');
            IX.Z(ix)=true;

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
            colormap(cmap);
            imagesc(DEMc)
            plot(Sn,'-r','LineWidth',2);
            scatter(xn,yn,20,'w','filled');
            hold off

            if isempty(EG)
                figure(1);

                subplot(2,1,2);
                hold on
                pl(1)=plot(Sn,'-r','LineWidth',2);
                sc(1)=scatter(xn,yn,20,'w','filled');
                hold off
                          
                subplot(2,1,1);
                hold on
                pl(2)=plot(Sn,'-r','LineWidth',2);        
                sc(2)=scatter(xn,yn,20,'w','filled');      
                hold off
            else
                figure(1);

                subplot(3,1,3);
                hold on
                pl(1)=plot(Sn,'-r','LineWidth',2);
                sc(1)=scatter(xn,yn,20,'w','filled');
                hold off

                subplot(3,1,2);
                hold on
                pl(2)=plot(Sn,'-r','LineWidth',2);
                sc(2)=scatter(xn,yn,20,'w','filled');
                hold off
                          
                subplot(3,1,1);
                hold on
                pl(3)=plot(Sn,'-r','LineWidth',2);
                sc(3)=scatter(xn,yn,20,'w','filled');
                hold off
            end


            qa=questdlg('Is this the watershed you wanted?','Basin Selection','No','Yes','Yes');

            switch qa
            case 'Yes'
                str1 = 'Y';
            case 'No'
                str1 = 'N';
            end

            delete(pl);
            delete(sc);
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
        xlabel('\chi')
        ylabel('Elevation (m)')
        title('\chi - Z')
        hold off

        subplot(2,1,2);
        hold on
        plotdz(Sn,DEMf,'dunit','km','Color','k');
        xlabel('Distance from Mouth (km)')
        ylabel('Elevation (m)')
        title('Long Profile')
        hold off

        qa2=questdlg('Keep this basin?','Basin Selection','No','Yes','Yes');       
        % prompt='    Keep this stream? Y/N [Y]: ';
        % str2=input(prompt,'s');
        % if isempty(str2) | strcmpi(str2,'Y');
        switch qa2
        case 'Yes'
            Outlets(ii,1)=xn;
            Outlets(ii,2)=yn;
            Outlets(ii,3)=ii;
            ii=ii+1;
            
            % Plot selected point on figures

            if isempty(EG)
                figure(1);
                subplot(2,1,2);
                hold on
                scatter(xn,yn,20,'w','filled');
                hold off
                          
                subplot(2,1,1);
                hold on
                scatter(xn,yn,20,'w','filled');
                hold off
            else
                figure(1);

                subplot(3,1,3);
                hold on
                scatter(xn,yn,20,'w','filled');
                hold off

                subplot(3,1,2);
                hold on
                scatter(xn,yn,20,'w','filled');
                hold off
                          
                subplot(3,1,1);
                hold on
                scatter(xn,yn,20,'w','filled');
                hold off
            end
            
            save('Outlets.mat','Outlets');
        end
        

        qa3=questdlg('Keep choosing basins?','Basin Selection','No','Yes','Yes'); 
        switch qa3
        case 'Yes'
            str2='Y';
            str1='N';
        case 'No'
            str2='N';
        end  

        
        close figure 2;
    end
end   
    
    
