N = length(observed_cells) - 1;
K = length(observed_cells{1}.location);

for k = 1:K
    current_cells = [];
    for i = 1:N
         if k >= observed_cells{i}.b_time && k < observed_cells{i}.d_time
             current_cells = [current_cells observed_cells{i}.location(:,k)];
         end
    end
    scatter(current_cells(1,:),current_cells(2,:),'filled')
    axis equal
    axis([0 50 0 50])
    drawnow;
end
