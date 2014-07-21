function D = histD2(tracks,params)
% function D = histD(tracks, params) calculates the diffusion coefficient
% for each track in "tracks" and plots a histogram of D.
%
% Stephan Uphoff. 09.09.11
%
% params contains the following parameters:
%
%modified from histD to ignore blinked frames when calculating the MSD

pixel = params.pixel; % length per pixel
dT = params.dT; % time per frame
sigmaNoise = params.sigmaNoise; % localization noise
DhistMinSteps = params.DhistMinSteps; % minimum number of steps for a track to be analyzed
rangeD = params.rangeD; % D range for the histogram

nMolecules = max(tracks(:,4));

MSD = zeros(nMolecules,1);
kk = 1;
track_lengths = zeros(nMolecules,1);


for ii = 1:nMolecules
    
    xx = find(tracks(:,4)==ii);
    num_frames = 0;
    
    if numel(xx) >= DhistMinSteps
        
        % sum all squared displacement in the track
        for jj = 1:numel(xx)-1
            
            if tracks(xx(jj+1),3) == tracks(xx(jj),3) + 1
            MSD(kk) = MSD(kk) + ((tracks(xx(jj+1),1) - tracks(xx(jj),1))^2 +...
                (tracks(xx(jj+1),2) - tracks(xx(jj),2))^2);
            num_frames = num_frames +1;
            end
            
        end
        
%         for jj = 1:numel(xx)-2
%             
%             if tracks(xx(jj+2),3) == tracks(xx(jj),3) + 2
%             MSD(kk) = MSD(kk) + ((tracks(xx(jj+2),1) - tracks(xx(jj),1))^2 +...
%                 (tracks(xx(jj+2),2) - tracks(xx(jj),2))^2);
%             num_frames = num_frames +1;
%             end
%         end       
        
        
        if num_frames >=DhistMinSteps
        MSD(kk) = MSD(kk)/num_frames; % mean square displacement
        track_lengths(kk) = num_frames;
        kk = kk + 1;
        else
        MSD(kk) = [];    
        end
    end

    
end
track_lengths(kk:end) = [];
MSD(kk:end) = []; % delete unused rows
MSD = MSD * pixel^2; % convert from pixel to length units

% calculate D from MSD and correct for localization noise
D = MSD/(4*dT) - sigmaNoise^2*pixel^2/dT;
%D = MSD/(8*dT) - sigmaNoise^2*pixel^2/dT;
noise_correction = sigmaNoise^2*pixel^2/dT;

mean(D)
% Dhist_plot1 = figure;
% axes1 = axes('Parent',Dhist_plot1,'LineWidth',3,'FontSize',16);
% box(axes1,'off');
% hold(axes1,'all');
% 
% % plot histogram of single-molecule diffusion coefficients over rangeD
% 
% hist(D,rangeD);
% xlim([min(rangeD), max(rangeD)]);
% xlabel('Apparent diffusion coefficient in \mum^{2}s^{-1}','FontSize',16);
% ylabel('Fraction of molecules','FontSize',16);


binSpacing = rangeD(2)-rangeD(1);
histDCount = histc(D,rangeD + binSpacing/2);
histDCount = histDCount./numel(D); %normalize area


Dhist_plot2 = figure;
axes2 = axes('Parent',Dhist_plot2,'LineWidth',3,'FontSize',16);
box(axes2,'off');
hold(axes2,'all');
bar1 = bar(rangeD+binSpacing,histDCount,'BarWidth',1);
baseline1 = get(bar1,'BaseLine');
set(baseline1,'LineWidth',3);

%
xlim([min(rangeD), max(rangeD)]);
xlabel('Apparent diffusion coefficient in \mum^{2}s^{-1}','FontSize',16);
ylabel('Fraction of molecules','FontSize',16);

x = rangeD+binSpacing; 
x = x';
n =DhistMinSteps;
mean_track_length = mean(track_lengths)


disp([' number of molecules = ' num2str(kk)]);

% if params.fitting
% c = Dhist_fitting(histDCount,x,n, noise_correction,'three_constrained',0.2459);
% keyboard;
% end   

 % c
%hold on;
%kde(D);
%sum(histDCount)


end