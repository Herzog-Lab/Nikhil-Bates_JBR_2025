% Part2 of rayleighplotcolor to overlay points

% This was how I export the plot from the original rayleighplotcolor
angles = PeakBottom; %febLpeak2 was a dataset I input
color = 'black';
lightcycle = 1;
output=rayleighplotcolor(angles,color,lightcycle);

%I am essentially overwriting the plot now with a new color
color = 'red';
angles = PeakMiddle; %new dataset is febLtrough2

anglesnew = zeros(length(angles),1);
for i = 1:length(angles)
    anglesnew(i,1) = angles(i);
end
angles = anglesnew;

radians = circ_ang2rad(angles);
[pval, z] = circ_rtest(radians);
mu = circ_mean(radians);
r = circ_r(radians);


scatr = .95 + (1-.95) .* rand(size(angles,1),1); 
hold on

polarscatter(radians.',scatr,400,color,'filled','MarkerFaceAlpha',0.75,'MarkerEdgeColor','none');
%you can play with these two commands to change the way the marker and line of confidence looks if
%you dont want to change color
polarplot([mu mu],[0 r],'Color',color,'LineWidth',4)

%adding title
title('Peak and Trough FebR Quadrants E16.5-E17.5')




