function chains = PointProcess(input_image,net,varargin)
% 
% %Load in image to label
% Im = rgb2gray(imread(input_image));%Test image
% Im = Im(19:530,38:549);
Im = input_image;
if ~isempty(varargin) && strcmp(varargin,'filter')
    if mean(Im(:)) < 90
        Im = Im + (100-mean(Im(:)));
        se = strel('square',5);
        Im = (imtophat(Im,se)+100)*.8;
    elseif length(find(Im > 250)) > 1000 || std(double(Im(:))) > 30
        se = strel('square',9);
        Im = (imtophat(Im,se)+100)*.7;
    elseif mean(Im(:)) < 120 && std(double(Im(:))) > 18
        se = strel('square',5);
        Im = (imtophat(Im,se)+100)*.8;
    elseif length(find(Im > 250)) < 1000 && kurtosis(Im(:)) < 0.1
        se = strel('square',12);
        Im = imtophat(Im,se)*2;
    end
end
Im([1 512],:) = 0; Im(:,[1 512]) = 0;
%Load network and label image
% netv5 = load('NetworkGen2.mat');%Better labeling and changed learning rate
I = semanticseg(Im,net);% Classify pixels with network
B = labeloverlay(Im,I);% overlay classification with image

%Turn the categorical label matrix into a binary image
%   then remove small and large objects
bw = logical(ismember(I,'Chains')+ismember(I,'Overlap'));
%  for n = 1:3
if strcmp(varargin,'actin')
    bw = logical(ismember(I,'Chains')+ismember(I,'Overlap')+ismember(I,'Noise'));
end
%     for i=2:length(bw)-1
%         for j=2:length(bw)-1
%             ul = bw(i,j-1) + bw(i+1,j) + bw(i+1,j-1);
%             ur = bw(i,j+1) + bw(i+1,j) + bw(i+1,j+1);
%             bl = bw(i,j-1) + bw(i-1,j) + bw(i-1,j-1);
%             br = bw(i,j+1) + bw(i-1,j) + bw(i-1,j+1);
%             if ul < 3 && ur < 3 && bl < 3 &&  br < 3
%                 bw(i,j) = 0;
%             end
%         end
%     end
%  end
if isempty(varargin) || strcmp(varargin,'filter')
    bw = bwareaopen(bw,40);
elseif strcmp(varargin{1},'EM')
    bw = bwareaopen(bw,20);
    se = strel('disk',4);
    bw = imclose(bw,se);
    bw = bwareaopen(bw,200);

end
chains = regionprops(bw,'Area');
cc = bwconncomp(bw,8);
L = labelmatrix(cc);
Lbw = ismember(L,find([chains.Area] > 400));
bw = ismember(L,find([chains.Area] <= 400));
sbw = ismember(L,find([chains.Area] <= 150));
se = strel('disk',4);
sc = regionprops(sbw,'Circularity','Image','BoundingBox');

for i = 1:length(sc)
    box = floor(sc(i).BoundingBox);
    if sc(i).Circularity > 0.4 
        sbw((box(2)+1:box(2)+box(4)),(box(1)+1:box(1)+box(3))) = ...
            sbw((box(2)+1:box(2)+box(4)),(box(1)+1:box(1)+box(3))) - sc(i).Image;
    end
end
sbw = imclose(sbw,se); sbw = bwareaopen(sbw,200);
bw = bw+sbw; bw(bw==2) = 1; bw = logical(bw);
bw = bwareaopen(bw,100);
chains = regionprops(bw,'Area','BoundingBox');
%Lbw = bwmorph(Lbw,'open');
if isempty(varargin) || strcmp(varargin,'filter')
    over = ismember(I,'Overlap');
    over = bwareaopen(over,10);
    oc = regionprops(over,'Circularity','Eccentricity','Image','BoundingBox');

    for i = 1:length(oc)
        box = floor(oc(i).BoundingBox);
        if oc(i).Circularity < 0.7 || oc(i).Eccentricity > 0.75
            over((box(2)+1:box(2)+box(4)),(box(1)+1:box(1)+box(3))) = ...
                over((box(2)+1:box(2)+box(4)),(box(1)+1:box(1)+box(3))) - oc(i).Image;
        end
    end

    oc = regionprops(over,'BoundingBox');
    temp = zeros(size(over)+14);
    for i = 1:length(oc)
        bnd = ceil(oc(i).BoundingBox);

        temp(bnd(2)+4:bnd(2)+bnd(4)+10,bnd(1)+4:bnd(1)+bnd(3)+10) = 1;
    end
    mask = temp(8:end-7,8:end-7);
    % bw = (bw-mask);
    % bw(bw<0)=0;
    % bw = logical(bw);
    Lbw = (Lbw-mask); Lbw(Lbw<0)=0; Lbw = logical(Lbw);
end
%Remove pixels with less than 3 adjacent white pixels
 for n = 1:2
    for i=2:length(bw)-1
        for j=2:length(bw)-1
            ul = bw(i,j-1) + bw(i+1,j) + bw(i+1,j-1);
            ur = bw(i,j+1) + bw(i+1,j) + bw(i+1,j+1);
            bl = bw(i,j-1) + bw(i-1,j) + bw(i-1,j-1);
            br = bw(i,j+1) + bw(i-1,j) + bw(i-1,j+1);
            if ul < 3 && ur < 3 && bl < 3 &&  br < 3
                bw(i,j) = 0;
            end
        end
    end
 end
 
  for n = 1:3
    for i=2:length(Lbw)-1
        for j=2:length(Lbw)-1
            ul = Lbw(i,j-1) + Lbw(i+1,j) + Lbw(i+1,j-1);
            ur = Lbw(i,j+1) + Lbw(i+1,j) + Lbw(i+1,j+1);
            bl = Lbw(i,j-1) + Lbw(i-1,j) + Lbw(i-1,j-1);
            br = Lbw(i,j+1) + Lbw(i-1,j) + Lbw(i-1,j+1);
            if ul < 3 && ur < 3 && bl < 3 &&  br < 3
                Lbw(i,j) = 0;
            end
        end
    end
 end

%Remove small objects again and get object properties
% Circ = regionprops(bw,'Circularity','Image','BoundingBox');
% 
% for i = 1:length(Circ)
%     box = floor(Circ(i).BoundingBox);
%     if Circ(i).Circularity > 0.4
%         bw((box(2)+1:box(2)+box(4)),(box(1)+1:box(1)+box(3))) = ...
%             bw((box(2)+1:box(2)+box(4)),(box(1)+1:box(1)+box(3))) - Circ(i).Image;
%     end
% end
bw = bwareaopen(bw,80,4);
chains = regionprops(bw,'BoundingBox','Image');

for i = 1:length(chains)
    skel = bwskel(chains(i).Image,'MinBranchLength',12);
    [y x] = find(skel);
    ep = bwmorph(skel,'endpoints');
    [p k] = find(ep);
    bp = bwmorph(bwmorph(skel,'thin'),'branchpoints');
    [r c] = find(bp);
    if ~isempty(r)
        for j = 1:length(r)
            temp = zeros(size(chains(i).Image)+14);
            temp(8:end-7,8:end-7) = chains(i).Image;
            temp(r(j)+2:r(j)+12,c(j)+2:c(j)+12) = 0;
            chains(i).Image = temp(8:end-7,8:end-7);
        end
        box = floor(chains(i).BoundingBox);
        bw((box(2)+1:box(2)+box(4)),(box(1)+1:box(1)+box(3))) = chains(i).Image;
    end
    if isempty(p) && ~isempty(y)
        temp = zeros(size(chains(i).Image)+14);
        temp(8:end-7,8:end-7) = chains(i).Image;
        temp(y(1)+2:y(1)+12,x(1)+2:x(1)+12) = 0;
        chains(i).Image = temp(8:end-7,8:end-7);
    end
end

ml = 15; maskL = 11;
if ~isempty(varargin) && strcmp(varargin{1},'EM')
    ml = 30; maskL = 20;
    Lbw = imclose(Lbw,strel('disk',5));
end
    % if isempty(varargin)

    for n = 1:2
        Lchains = regionprops(Lbw,'Area','BoundingBox','Image');
        for i = 1:length(Lchains)
            skel = bwskel(Lchains(i).Image,'MinBranchLength',ml);
            bp = bwmorph(bwmorph(skel,'thin'),'branchpoints');
            [r c] = find(bp);            
            box = floor(Lchains(i).BoundingBox);
    %         temp = zeros(size(Lbw)+14);
            if ~isempty(r)
                for j = 1:length(r)
                    temp = zeros(size(Lbw)+14);
    %                 temp(8:end-7,8:end-7) = Lchains(i).Image;
                    temp(box(2)+r(j)+3:box(2)+r(j)+maskL,box(1)+c(j)+3:box(1)+c(j)+maskL) = 1;
                    temp = temp(8:end-7,8:end-7);
                    temp((box(2)+1:box(2)+box(4)),(box(1)+1:box(1)+box(3))) = ...
                        temp((box(2)+1:box(2)+box(4)),(box(1)+1:box(1)+box(3))) + Lchains(i).Image;
                    temp(temp == 1) = 0; temp(temp == 2) = 1;
                    Lbw = (Lbw - temp); Lbw(Lbw < 0) = 0; Lbw = logical(Lbw);
    %                 Lchains(i).Image = temp(8:end-7,8:end-7);
                end

    %             Lbw((box(2)+1:box(2)+box(4)),(box(1)+1:box(1)+box(3))) = Lchains(i).Image;
            end
        end
        se = strel('square',3);Lbw = imclose(Lbw,se);ml = 3;
    end
% end
bw = bwareaopen(bw,80,8);
Lbw = bwareaopen(Lbw,80,8);

bw = logical(bw+Lbw);
if ~isempty(varargin) && strcmp(varargin{1},'EM')
    bw = Lbw;
    bw = imclose(bw,strel('disk',5));
end
chains = regionprops(bw,'BoundingBox','Image');
MBL = 15;
if ~isempty(varargin) && strcmp(varargin{1},'EM')
    MBL = 30;
end
EdgeFlag = zeros(1,length(chains));
for i = 1:length(chains)
    [r1 c1] = size(chains(i).Image);
    pic = chains(i).Image;
    box = chains(i).BoundingBox;

    %Get the skeleton of the chain
    skel = bwskel(pic,'MinBranchLength',MBL);    
    chains(i).skeleton = skel;
    [r c] = find(skel);
    N = length(r);
    bp = bwmorph(bwmorph(skel,'thin'),'branchpoints');
    endpnts = bwmorph(skel,'endpoints');
    [y x] = find(endpnts); 
    if isempty(x) && ~isempty(r)
        skel(r(1),c(1)) = 0;
        endpnts = bwmorph(skel,'endpoints');
        [y x] = find(endpnts);
    elseif length(x)==1 && length(bp) == 1
        skel = skel - bp;
    end

    if ~isempty(x) && all(bp(:)==0) && length(r) > 30 && length(r) < 300
        for j = 1:N/2
                [r c] = find(skel);
                endpnts = bwmorph(skel,'endpoints');
                if j == 1
                    if std(r) <= 2 || std(c) <= 2
                         EdgeFlag(i) = 1;
                    end
                end
                if j == 1
                    [y x] = find(endpnts);    
                    chains(i).points(1,j) = x(1);
                    chains(i).points(2,j) = y(1);
                    chains(i).points(1,N-j) = x(2);
                    chains(i).points(2,N-j) = y(2);
                    skel = skel - endpnts;
                elseif j > 1
                    [y x] = find(endpnts);

                    if length(x) == 1
                        
                        
                    elseif abs(x(1) - chains(i).points(1,j-1)) > 2 || ...
                            abs(y(1) - chains(i).points(2,j-1)) > 2 ...
                            abs(x(2) - chains(i).points(1,N-j+1)) > 2 ...
                            || abs(y(2) - chains(i).points(2,N-j+1)) > 2;

                        chains(i).points(1,j) = x(2);
                        chains(i).points(2,j) = y(2);
                        chains(i).points(1,N-j) = x(1);
                        chains(i).points(2,N-j) = y(1);
                        skel = skel - endpnts;
                    else
                        chains(i).points(1,j) = x(1);
                        chains(i).points(2,j) = y(1);
                        chains(i).points(1,N-j) = x(2);
                        chains(i).points(2,N-j) = y(2);
                        skel = skel - endpnts;
                    end
                end
        end
        box = floor(box);
        chains(i).points(1,:) = chains(i).points(1,:) + box(1);
        chains(i).points(2,:) = chains(i).points(2,:) + box(2);
        if EdgeFlag(i) == 1
            chains(i).points = [];
        end
    else
        chains(i).points = [];
    end
end



end

