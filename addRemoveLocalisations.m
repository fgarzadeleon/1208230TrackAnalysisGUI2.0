function [pos] = addRemoveLocalisations(pos,fig,loadname,localizationThresh,testDetectionFrame)
flag = 1;
while flag
    if ~any(isnan(pos))
        
        hold all
        h_plot = plot(pos(:,1),pos(:,2),'ro','MarkerSize',10);
        xlabel('x [pixels]');
        ylabel('y [pixels]');
        axis image;
    end
    figure(fig)
    dcm_obj = datacursormode(fig);
    set(dcm_obj,'DisplayStyle','datatip',...
        'SnapToDataVertex','on','Enable','on')
    
    disp('Click on the image to display a data tip, then press Return.')
    disp('To Select multiple press Alt+LeftClick.')
    pause                            % Wait while the user does this.
    
    c_info = getCursorInfo(dcm_obj);
    %set(c_info.Target,'LineWidth',2)  % Make selected line wider
    if isempty(c_info)
        flag=0;
    end
    
    
    for indexStruct = 1:length(c_info)
       
        posDif = [];
        posDif(:,1) = abs(pos(:,1)-c_info(indexStruct).Position(1));
        posDif(:,2) = abs(pos(:,2)-c_info(indexStruct).Position(2));
        posDifLow = posDif==0;
        indexDelete = sum(posDifLow,2)==2;
        if ~isempty(indexDelete)
            if any(indexDelete)
                pos(indexDelete,:) = [];
            else
                [posTemp(:,1),posTemp(:,2)] = gaussStormSingleParticle(loadname,...,
                    localizationThresh, 3,...,
                    testDetectionFrame,c_info(indexStruct).Position,0);
                pos = [pos; posTemp];
            end
        end
        
    end
    delete(findall(fig,'Type','hggroup','HandleVisibility','off'));
    x = [];
    y = [];
    set(h_plot, 'xdata', x, 'ydata', y);
end



