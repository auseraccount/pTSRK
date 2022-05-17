% fig_init initializes the axis and figure to my preferable standard

function fig_init(fig_pos)
if nargin < 1
    fig_pos = [0 0 8 6];
end
set(gcf,'units','inches')
set(gcf,'position',fig_pos,'paperpositionmode','auto')
set(gca,'FontSize',20,'LineWidth',1.5)
end

