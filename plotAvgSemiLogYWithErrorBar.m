function [linehandle, e_handle] = plotAvgSemiLogYWithErrorBar(avg, std, freq, color, linewidth)
% plotAvgWithSD: Plot two lines between Y1 (lower) and Y2 (upper) in
% semilog scale (return p for handling the legend selectively)
X = 1 : size(avg, 2);
X_bar =  0 : freq : size(avg, 2); 
X_bar = X_bar( 2 : (end - 1));
Y_bar = avg(X_bar);
Error_bar = std(X_bar);
% default values
if nargin < 5; linewidth = 1.5; end
if nargin < 4; color = 'b'; end % default color is blue
%% Plot and fill
linehandle = plot(X, avg, color, 'LineWidth', linewidth); hold on;
e_handle = errorbar(X_bar, Y_bar, Error_bar);  
set(gca,'YScale','log');

e_handle.Marker = 'x';
e_handle.MarkerSize = 5;
e_handle.Color = color;
e_handle.CapSize = 5;
e_handle.LineStyle = 'None';

end
