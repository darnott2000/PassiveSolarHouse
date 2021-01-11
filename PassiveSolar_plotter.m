firstSweepRange = range1;   
secondSweepRange = range2;
%Matrix of Heats from Optimize Function
Heats = allheat;
%number of steps between the ranges (gets bumped to 201 from 
%200 because of matlab quirks)
numSteps = 201; 

%This will give us the positions of the data from the array
Xwidth = linspace(firstSweepRange(1), firstSweepRange(2),numSteps);
Ywidth = linspace(secondSweepRange(1), secondSweepRange(2),numSteps);

%this line is not necessary since our data is in this format
% F = scatteredInterpolant(Xwidth,Ywidth,Heats);

[xq,yq] = meshgrid(Xwidth,Ywidth);
grid on;

%create a surface with the heat data
surf(xq,yq,Heats,"EdgeColor","none");
hold on;
%show meaning from colors
colorbar
colormapeditor
title("Average Temperature over 30 years period for Tile & Insulation Thicknesses")
ylabel("Thickness of Insulation (m)")
xlabel("Thickness of Tile (m)")
zlabel("Temperature (degrees Celcius)")
scatter3(.0373,.4,23.7119,'rx')
text(.0373,.4,23.7119,"Optimal","HorizontalAlignment","right")
hold off;

