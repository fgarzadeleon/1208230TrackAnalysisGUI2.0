function plotTracksColourCoded(tracks,params)
% function plotTracksColourCoded(tracks,minStep,plotTracksColour) plots all tracks that have
% at least minStep steps. Tracks is the output of the tracking software:
% tracks(:,1) = x coordinates
% tracks(:,2) = y coordinates
% tracks(:,3) = frame
% tracks(:,4) = track id
%
% plotTracksColour defines the one letter colour of the tracks.
% plotTracksColour = 'random' plots tracks in random colours
% plotTracksColour = 'time' plot tracks coloured coded by order of appearance in the
% movie. Tracks early in the movie are plotted blue, later tracks are
% plotted red. 
%
% Stephan Uphoff. 21.09.11


minStep = params.plotTracksMinSteps;
plotTracksColour = params.plotTracksColour;
holdFigure = params.holdFigureCheckbox;
if holdFigure
    figureHandle = params.figureHandle;    
end

nMolecules = max(tracks(:,4)); % number of tracks

cmap = colormap(jet); % load colour map

if holdFigure
    figure(figureHandle);
else
    figure;
end
hold all


for jj = 1:nMolecules % loop over tracks
    
    % indexes of this track in the tracks data
    xx = find(tracks(:,4)==jj);
    
    if numel(xx)>minStep % include only tracks with at least minSteps
        
        if strcmp(plotTracksColour, 'time')
        
        % choose the colour by the track id
        colorindex = ceil( length(cmap) * tracks(xx(1),4) / nMolecules );
        
        % plot the track in the chosen colour
        plot(tracks(xx,1),tracks(xx,2),...
            '-','Color',cmap(colorindex,:),'lineWidth', 1.5)
        
        elseif strcmp(plotTracksColour,'random') 
        plot(tracks(xx,1),tracks(xx,2),...
            '-','Color',rand(1,3),'lineWidth', 1.5)            
        
        else
         plot(tracks(xx,1),tracks(xx,2),...
            '-','Color',plotTracksColour,'lineWidth', 1.5)
        
        end
       % plot(mean(tracks(xx,1)), mean(tracks(xx,2)),'-')
    end
    
end

axis image; % same scale on x and y axis.
hold off;

end