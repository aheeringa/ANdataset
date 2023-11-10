function [differenceTwo,oor] = bandwidthBF(x_axis, y_axis, peak, offset)
%%
% input:
%   x_axis: vector of the x axis values
%   y_axis: vector of the y axis values
%   peak: scalar or vector of the peak(s) of the best frequency (max peak will be figured out)
%   offset: scalar which gets subtracted from the amplitude, e.g. mean spontaneous rate
%output:
%   differenceTwo: scalar of the difference in Hz between the 50% values
%   oor: out of range, oor = 0 when the two scalars fall inside the x-axis
%      values, oor = 1 when the outer bandwidth edge could not be faithfully
%      determined.

% By: Amarins Heeringa
%%

%value to find for the 50% value for the bandwidth
value_half = round(max(peak)-(max(peak)-offset)/2); %subtracts the offset of the peak to get the absolute amplitude; 
                                                    %divides that by 2 to get the 50% value which is subtracted from the peak to get the exact 50% excluding the offset
% find the index of the peak
[~,xPeak]=max(y_axis-max(peak)); 
ind1=xPeak;
while y_axis(ind1)>value_half && ind1>1
    ind1=ind1-1;
end
ind2=xPeak;
while y_axis(ind2)>value_half && ind2<length(y_axis)
    ind2=ind2+1;
end
indx=[ind1 ind2];

% determine out of range
if ind1==1 && y_axis(ind1)>value_half
    oor=1;
elseif ind2==length(y_axis) && y_axis(ind2)>value_half
    oor=1;
else 
    oor=0;
end

%calculating the bandwidth
differenceTwo = diff(x_axis(indx));
    
%calculating the points of the 50% points to draw a line
y = [value_half,value_half];
x = x_axis(indx);
lw = 1.5;

%plotting the line for bandwidth
plot(x,y,'b-', 'LineWidth',lw)
strmin = ['Bandwidth ', num2str(differenceTwo), ' HZ'];  %creating a string box
text(x_axis(indx(2))+70,y_axis(indx(2)),strmin,'HorizontalAlignment','left'); %plotting the text box on the points

end