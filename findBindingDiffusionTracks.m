function findBindingDiffusionTracks(tracks, params)
% Stephan Uphoff. 14.09.11

sigmaNoise = params.sigmaNoise;
pixel = params.pixel;
dT = params.dT;
DThresh = params.DThresh;
minStep = params.plotTracksMinSteps;
holdFigure = params.holdFigureCheckbox;
if holdFigure
    figureHandle = params.figureHandle;    
end

nMolecules = max(tracks(:,4)); % number of tracks

D1 = zeros(nMolecules,1);
D2 = zeros(nMolecules,1);

if holdFigure
    figure(figureHandle);
else
    figure;
end
hold all

kk = 1;

for jj = 1:nMolecules % loop over tracks
    
    % indexes of this track in the tracks data
    xx = find(tracks(:,4)==jj);
    
    %plot(tracks(xx,1),tracks(xx,2),...
    %        '-','Color',[0.5 0.5 0.5])
    
    if numel(xx)>minStep % include only tracks with at least minSteps
        
        % compare D between two sections of a track
        MSD1 = 0; % reset MSD1
        MSD2 = 0; % reset MSD2
        
        u = [];
        v = [];
        
        % sum all squared displacement
        for ii = 1:floor(numel(xx)/2)-1
            
            MSD1 = MSD1 + ((tracks(xx(ii+1),1) - tracks(xx(ii),1))^2 +...
                (tracks(xx(ii+1),2) - tracks(xx(ii),2))^2);
            
            u(ii) = tracks(xx(ii+1),1) - tracks(xx(ii),1);
            v(ii) = tracks(xx(ii+1),2) - tracks(xx(ii),2);
            
        end
        
        for ll = floor(numel(xx)/2):numel(xx)-1
            
            MSD2 = MSD2 + ((tracks(xx(ll+1),1) - tracks(xx(ll),1))^2 +...
                (tracks(xx(ll+1),2) - tracks(xx(ll),2))^2);
            
            u(ll) = tracks(xx(ll+1),1) - tracks(xx(ll),1);
            v(ll) = tracks(xx(ll+1),2) - tracks(xx(ll),2);
            
        end        
        
        if MSD1==0
            pause()
        end
        MSD1 = (MSD1/ii) * pixel^2; % mean square displacement
        MSD2 = (MSD2/ll) * pixel^2; % mean square displacement        
 
        % calculate D from MSD
        D1(kk) = MSD1/(4*dT) - (sigmaNoise*pixel)^2/dT;
        D2(kk) = MSD2/(4*dT) - (sigmaNoise*pixel)^2/dT;

        if (D1(kk)<DThresh(1) && D2(kk)>DThresh(2)) ||...
                (D2(kk)<DThresh(1) && D1(kk)>DThresh(2))   

        quiver(tracks(xx(1:end-1),1),tracks(xx(1:end-1),2),u',v',0,'LineWidth',2)
        
        end
        
        kk = kk+1;  
        
    end
    
end

D1(kk:end) = [];
D2(kk:end) = [];

axis image; % same scale on x and y axis.
hold off;

end