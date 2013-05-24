function D = histD(tracks,params)
% function D = histD(tracks, params) calculates the diffusion coefficient
% for each track in "tracks" and plots a histogram of D.
%
% Stephan Uphoff. 09.09.11
%
% params contains the following parameters:

pixel = params.pixel; % length per pixel
dT = params.dT; % time per frame
sigmaNoise = params.sigmaNoise; % localization noise
DhistMinSteps = params.DhistMinSteps; % minimum number of steps for a track to be analyzed
rangeD = params.rangeD; % D range for the histogram

nMolecules = max(tracks(:,4));

MSD = zeros(nMolecules,1);
kk = 1;

for ii = 1:nMolecules
    
    xx = find(tracks(:,4)==ii);
    
    if numel(xx) > DhistMinSteps
        
        % sum all squared displacement in the track
        for jj = 1:numel(xx)-1
            
            MSD(kk) = MSD(kk) + ((tracks(xx(jj+1),1) - tracks(xx(jj),1))^2 +...
                (tracks(xx(jj+1),2) - tracks(xx(jj),2))^2);
            
        end
        
        MSD(kk) = MSD(kk)/jj; % mean square displacement
        kk = kk + 1;
        
    end
    
    ii
    
end

MSD(kk:end) = []; % delete unused rows
MSD = MSD * pixel^2; % convert from pixel to length units

% calculate D from MSD and correct for localization noise
D = MSD/(4*dT) - sigmaNoise^2*pixel^2/dT;

% plot histogram of single-molecule diffusion coefficients over rangeD
figure;
hist(D,rangeD);
xlim([min(rangeD), max(rangeD)]);
xlabel('diffusion coefficient [um^2/s]');
ylabel('histogram count');

end