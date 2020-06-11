function ColorSet=varycolor(NumberOfPlots)
% VARYCOLOR Produces colors with maximum variation on plots with multiple
% lines.
%
%     VARYCOLOR(X) returns a matrix of dimension X by 3.  The matrix may be
%     used in conjunction with the plot command option 'color' to vary the
%     color of lines.  
%
%     Uses colors based on the MATLAB 'lines' colour set for plots.
% 
%     Example Usage:
%         NumberOfPlots=50;
%
%         ColorSet=varycolor(NumberOfPlots);
% 
%         figure
%         hold on;
% 
%         for m=1:NumberOfPlots
%             plot(ones(20,1)*m,'Color',ColorSet(m,:))
%         end

% Based on the concept of VARYCOLOR by Daniel Helmick 8/12/2008
% Written Steven Johnson 13/03/2015


narginchk(1,1)  %Check number of input arguments.
nargoutchk(0,1) %Check number of output arguments.

% Based on Standard MATLAB 'lines' colors
%      0    0.4470    0.7410
% 0.8500    0.3250    0.0980
% 0.9290    0.6940    0.1250
% 0.4940    0.1840    0.5560
% 0.4660    0.6740    0.1880
% 0.3010    0.7450    0.9330
% 0.6350    0.0780    0.1840

R = [     0    0.8500    0.9290    0.4940    0.4660    0.3010    0.6350];
G = [0.4470    0.3250    0.6940    0.1840    0.6740    0.7450    0.0780];
B = [0.7410    0.0980    0.1250    0.5560    0.1880    0.9330    0.1840];

spaces = length(R) - 1;
lines = 0:spaces;
nPlotValues = 0:spaces/(NumberOfPlots-1):spaces;

valueR = interp1(lines, R, nPlotValues);
valueG = interp1(lines, G, nPlotValues);
valueB = interp1(lines, B, nPlotValues);

ColorSet = [valueR;valueG;valueB]';