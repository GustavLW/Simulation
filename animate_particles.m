clear all 
close all
clc
dd      = dir('Datasets');  
dd      = dd(3:end);
nFiles  = length(dd);
save    = 0;
for f = 9:nFiles
    close all
    load(['Datasets\' dd(f).name])
    sample_cell = observed_cells{1};
    K           = length(sample_cell.location);
    clear sample_cell
    N           = length(observed_cells)-2;
    ECEAT       = NaN*ones(N,2,K);

    for i = 1:N
        for k = 1:K
            if k >= observed_cells{i}.b_time
                ECEAT(i,:,k) = observed_cells{i}.location(:,k);
            end
        end
    end
    
    h = figure('units','normalized','outerposition',[0 0 1 1]);
    axis tight manual; % this ensures that getframe() returns a consistent size
    filename = strcat(dd(f).name(1:end-4),'.gif');
    %
    skip = 2;
    for k = skip:skip:K
        clf
        picture = squeeze(ECEAT(:,:,k));
        viscircles(picture,0.5*ones(length(picture),1))
        hold on
        plot([0 50 50 0 0],[0 0 50 50 0],'k')
        axis equal
        axis([-5 55 -5 55])
        drawnow;
        % Capture the plot as an image
        if save == 1
            frame = getframe(h);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            % Write to the GIF File
            if k == skip
                imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
            else
                imwrite(imind,cm,filename,'gif','WriteMode','append');
            end
        end
        pause(0.05)
    end
    pause(1)
end