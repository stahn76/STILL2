function [ ] = setColorbar( rng1, rng2 )



colorbar('off');

subplots = get(gcf,'Children'); % Get each subplot in the figure
for i=1:length(subplots) % for each subplot
    caxis(subplots(i),[rng1,rng2]); % set the clim
end

h=colorbar('location','Manual', 'position', [0.93 0.43 0.03 0.2]);
% colormap(colormap(brewermap([],'*Spectral')));
% colormap(brewermap([],'*RdYlBu'))
% colormap('hot');

% set(get(h,'title'),'string','WPLI');


end
