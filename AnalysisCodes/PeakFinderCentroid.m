function [xy] = PeakFinderCentroid(I,threshold,edge_pad,bg_subtract, radius)

% Find maxima in image using built-in MATLAB functions for centroid finding


xdim = size(I,2);
ydim = size(I,1);

if strcmp(bg_subtract,'y') == 1
    I = I - imopen(I,strel('disk',10));
end

I = im2double(I);
I = I/max(max(I));

Ilinear = reshape(I,1,size(I,1)*size(I,2));
Istd = std(Ilinear);
% figure (2);
% imshow(I);

pause(0.01);

hold on

for i = 1:size(I,1)
    for j = 1:size(I,2)
        if I(i,j)<0
            I(i,j)=0;
        end
    end
end

ISupp = imhmax(I,Istd*threshold, 8);
IrMax = imregionalmax(ISupp);
regions = bwlabel(IrMax);
centroids = regionprops(regions, 'centroid');
RegArr = struct2cell(centroids);
RegMat = cell2mat(RegArr.');

N=size(RegMat,1);
molecules = zeros(0,2);

for n = 1:N
    if RegMat(n,1) > edge_pad && RegMat(n,1) < xdim-edge_pad && RegMat(n,2) > edge_pad && RegMat(n,2) < ydim-edge_pad
        molecules=cat(1,molecules,RegMat(n,:));
    end
end


%% Sujay Mod.
stats = regionprops(regions, 'MajorAxisLength','MinorAxisLength');
a1=struct2cell(stats);
a2=cell2mat(a1.');

b1=a2(:,1)-a2(:,2);
b2=b1<radius;
a3=mean(a2,2);
a4=a3<2*radius;
a5=b2.*a4;
NewList=RegMat((a5==1),:);

N=size(NewList,1);
NewMolecules = zeros(0,2);

for n = 1:N
    if NewList(n,1) > edge_pad && NewList(n,1) < xdim-edge_pad && NewList(n,2) > edge_pad && NewList(n,2) < ydim-edge_pad
        NewMolecules=cat(1,NewMolecules,NewList(n,:));
    end
end

%%


% plot(molecules(:,1),molecules(:,2), 'ro');
% hold on;
% plot(NewMolecules(:,1),NewMolecules(:,2), 'go');
% hold off

xy = cat(1,NewMolecules(:,2)',NewMolecules(:,1)');

end

