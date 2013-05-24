function [diffusionFraction, meanD1, meanD2] = twoSpeciesMSD2Threshold_vStephan(tracks, params)
% MSD diffusion analysis

cmap = colormap(jet);
pixel = params.pixel;
dT = params.dT;
nMolecules = max(tracks(:,4));
sigmaNoise = params.sigmaNoise;
threshold1 = params.DThresh(1);
threshold2 = params.DThresh(2);
minSteps = params.DhistMinSteps;
rangeD = params.rangeD;
holdFigure = params.holdFigureCheckbox;
if holdFigure
    figureHandle = params.figureHandle;
end

tracks(:,1:2) = tracks(:,1:2) * pixel;

D1 = zeros(nMolecules,1);
D2 = zeros(nMolecules,1);
meanPos = zeros(nMolecules,2);
kk = 1;
ll = 1;

if holdFigure
    figure(figureHandle);
    axis image;
else
    figure;
end
hold all

immobileTracks = [];

for ii = 1:nMolecules
    
    xx = find(tracks(:,4)==ii);
    
    if numel(xx)>1
%          plot(tracks(xx,1)/pixel,tracks(xx,2)/pixel,...
%              'b-','lineWidth', 1)
                        colorindex = ceil( length(cmap) * tracks(xx(1),4) / nMolecules );
            
            %plot the track in the chosen colour
                        plot(tracks(xx,1)/pixel,tracks(xx,2)/pixel,...
                            '-','Color',cmap(colorindex,:),'lineWidth', 2)
    end
    
    if numel(xx)>minSteps
        
        MSD1 = 0; % reset MSD
        MSD2 = 0;
        
        % sum all squared displacement
        for jj = 1:numel(xx)-1
            
            MSD1 = MSD1 + ((tracks(xx(jj+1),1) - tracks(xx(jj),1))^2 +...
                (tracks(xx(jj+1),2) - tracks(xx(jj),2))^2);
            
        end
        MSD1 = MSD1/jj; % mean square displacement
        
        for jj = 1:numel(xx)-2
            
            MSD2 = MSD2 + ((tracks(xx(jj+2),1) - tracks(xx(jj),1))^2 +...
                (tracks(xx(jj+2),2) - tracks(xx(jj),2))^2);
            
        end       
        MSD2 = MSD2/jj; % mean square displacement
        
        % calculate D from MSD and correct for localization noise
        D1(kk) = MSD1/(4*dT) - sigmaNoise^2*pixel^2/dT;
        D2(kk) = MSD2/(8*dT) - sigmaNoise^2*pixel^2/dT;
        
        if D1(kk) < threshold1 && D2(kk) < threshold2
            
%                         colorindex = ceil( length(cmap) * tracks(xx(1),4) / nMolecules );
%             
%             %plot the track in the chosen colour
%                         plot(tracks(xx,1)/pixel,tracks(xx,2)/pixel,...
%                             '-','Color',cmap(colorindex,:))
            
            immobileTracks = [immobileTracks; tracks(xx,:)];
            %             plot(tracks(xx,1)/pixel,tracks(xx,2)/pixel,...
            %                 'r-')
            %
            meanPos(ll,:) = [mean(tracks(xx,1)/pixel),mean(tracks(xx,2)/pixel)];
            ll = ll + 1;
            
        else
%                    plot(tracks(xx,1)/pixel,tracks(xx,2)/pixel,...
%            'b-','lineWidth', 2)
        end
        
        % reset MSD and update index
        kk = kk + 1;
        
    end
    
end


D1(kk:end) = []; % delete unused rows
D2(kk:end) = []; % delete unused rows
meanPos(ll:end,:) = [];

for jj = 1:nMolecules
    xx = find(immobileTracks(:,4)==jj);
    
    if ~isempty(xx)
%        colorindex = ceil( length(cmap) * immobileTracks(xx(1),4) / nMolecules );
        
        % plot the track in the chosen colour
%                 plot(immobileTracks(xx,1)/pixel,immobileTracks(xx,2)/pixel,...
%                     '-','Color',cmap(colorindex,:),'lineWidth', 1.5)
        
%                 plot(immobileTracks(xx,1)/pixel,immobileTracks(xx,2)/pixel,...
%                     '-','Color',rand(1,3),'lineWidth', 1.5)
        
%       plot(immobileTracks(xx,1)/pixel,immobileTracks(xx,2)/pixel,...
%            'r-','lineWidth', 2)
        
    end
    
end

%plot mean positions of immobile molecules
% plot(meanPos(:,1),meanPos(:,2),'ko');

hold off

% save('immobileMoleculesMeanPos','meanPos');

D1c = histc(D1,rangeD'); %normalized histogram count
D2c = histc(D2,rangeD'); %normalized histogram count

figure;
hist(D1,rangeD');
xlim([min(rangeD) max(rangeD)]);
hold all
stem(threshold1,max(D1c),'marker','none');
xlabel('diffusion coefficient [um^2/s]');
ylabel('histogram count');
hold off

figure;
hist(D2,rangeD');
xlim([min(rangeD) max(rangeD)]);
hold all
stem(threshold2,max(D2c),'marker','none');
xlabel('diffusion coefficient [um^2/s]');
ylabel('histogram count');
hold off

diffusionFraction = (ll-1)/(kk-1)

observedFraction = numel(D1)/nMolecules

meanD1 = NaN;
meanD2 = NaN;

%keyboard

end