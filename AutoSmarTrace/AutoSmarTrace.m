function AutoSmarTrace(input_path,out_filename,nmperpix,varargin)
% --- Auto tracing changes to be made:
% Replace all user input with one input at the beginning
%   Ideas: Change to a function that takes a string and datastore as input.
%           The string is the filename for saving, datastore contains all
%           images to be traced.
%         The algorithm should remain the same, but loop through each chain
%           in each image using the initial points given by the skeleton
%           from 'PointProcess.m'.
% Remove most most or all of the GUI.  No user input means no graphics
% required. Maybe a message box or something to show user it is running.


%---------------------------------------------------------------
% copyright 2010-2016 Guillaume Lamour, Naghmeh Rezaei
% *******Do NOT distribute*******
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Edits by Mathew Schneider:
%Lots
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%---------------------------------------------------------------

% Edit the above text to modify the response to help SmarTrace


%{
flag: Decides when to stop or adjust tracing parameters (incremented in steps of 0.5 or 1)
    Loop breaks at flag = 2
pt_sep: Controls how many points skipped in each chain's skeleton
varargin: optional parameter, if EM passes skip a larger step (40) cuz EM images might be huge or need sparser sampling
%}

flag = 0;
pt_sep = 5;
if ~isempty(varargin) && strcmp(varargin{1},'EM')
    pt_sep = 40;
end

%{
Loop the script until flag <= 2. When traced chain has sections with curvature below certain threshold, tracing too many details/overfit
-> change the seperation distance (pt_sep 1 to 3) and increment flag by 0.5.
When parameters not found, more comprehensive analyis run (SmTr_ChainSamplingTracing and SmTr_AnalysisTracing) and flag incremented 
by 1.
When chain traces successfully and no more optimization to be done flag set to 2 and while look ended.
%}
while flag < 2

clc
handles = struct; %Handles store any data throughout processing (current image, resolution, output file name, results etc.)
handles.nmperpx = nmperpix; %Appends the size of one pixel into the handles struct
imds = imageDatastore(input_path); %Automatically manages collections of images in the input folder. [Find, Track, provide access to iterate etc.]
FilePaths = imds.Files; %Pulls the file paths of images from imds struct into the variable
N_images = length(imds.Files); %Number of images found
N = 0; %Counter to enumerate extracted chains across the images
ifield = 'chain';
ivalue = {[];[];[]};
directData = struct(ifield,ivalue); %Create a struct that holds the name 'chain' and 3x1 empty array
ImDir = [out_filename, 'Traces']; %Folder name to save traced images
mkdir(ImDir); %Create the save image folder if not existing
netv5 = load('NetworkGen2.mat'); %Load the pretrained neural net to analyze and extract the chain skeletons from the images.


%Run the for loop for all images found in the directory 
for im_num = 1:N_images
    disp(['Tracing Image #' num2str(im_num) '...']) %Track the current image number
    %Pull data from imds (imageDatastore), output directory, handles (image info), current image, image filepaths and load it to handles
    handles = load_data(imds, out_filename, handles,im_num,FilePaths);
    chains = struct; %Create structure to process chains on the image

    %{
    Call custom pointProcess function to process chains from the images.
    Different behaviour depending on the optional parameter. If parameter EM, PointProcess takes in Em as a parameter If it does not 
    exist, PointProcess runs without it and if it is named actin, runs actin in PointProcess
    %}
    if ~isempty(varargin) && strcmp(varargin{1},'EM') 
        chains1 = PointProcess(handles.A,netv5.net,'EM');
    elseif isempty(varargin)
        chains1 = PointProcess(handles.A,netv5.net);
    elseif ~isempty(varargin) && strcmp(varargin{1},'actin') 
        chains1 = PointProcess(handles.A,netv5.net,'actin');
    end
    
    %Filters out empty chains and assigns the nonempty ones to the chains struct
    N_chains = length(chains1);
    j=1;
    for i2=1:N_chains
        if ~length(chains1(i2).points) == 0
            chains(j).points = chains1(i2).points;
            j=j+1;
        end
    end
    N_chains = length(chains);

%Start plotting the images 
    %Framework for the image
    figure(1)
    imshow(handles.A);
    colormap(bone(255));
    if im_num == 1
        set(gcf,'Position',[300, 150, 700, 700]);
    end
    hold on
    
    for chain_num = 1:N_chains
        figure(1)
        hold on
        if length(chains(chain_num).points) == 0 || ~isfield(chains,'points')
            N=N;
        else
        handles = getpoints(handles,directData,chains,[chain_num N],pt_sep);
        handles = combo_calc(handles,chain_num + N);
        if any(abs(handles.curv(chain_num+N).curvature) > 1) && pt_sep == 1
            pt_sep = 3;flag = flag + .5;
            break
        end

    %gif code
        % --- BEGIN GIF ANIMATION BLOCK ---

        % Force MATLAB to finish drawing the plot updates
        drawnow; 

        % Capture the current figure window as a frame
        frame = getframe(figure(1));
        im = frame2im(frame);

        % Convert the captured RGB image to an indexed image, which is necessary for the GIF format
        [imind,cm] = rgb2ind(im,256); 

        % Define a unique filename for the GIF within the 'Traces' directory
        gif_filename = fullfile(ImDir, ['trace_animation_' handles.filename '.gif']);

        % Write the frame to the GIF file
        if chain_num == 1
            % For the first chain, create a new GIF file. 
            % 'Loopcount',inf makes the GIF loop forever.
            % 'DelayTime' sets the time in seconds between frames.
            imwrite(imind,cm,gif_filename,'gif', 'Loopcount',inf, 'DelayTime',0.2);
        else
            % For all subsequent chains, append the new frame to the existing file
            imwrite(imind,cm,gif_filename,'gif','WriteMode','append', 'DelayTime',0.2);
        end

            % --- END GIF ANIMATION BLOCK ---
            % ================================================================





        end
    end


    if flag == 1.5 
        flag = 1;
        break
    end
    N = N + N_chains;
    handles.figure = gcf;
    if im_num == 3 && flag == 0 && isempty(varargin)
        sampledChains = SmTr_ChainSamplingTracing(handles);
        P = SmTr_AnalysisTracing(sampledChains,30,0);
        pt_sep = ceil(.0015*exp((P-18)/20))+ceil(4*log((P-13)/15));
        if pt_sep < 1
            pt_sep = 1;
        end
        flag = flag + 1;
%         if P <= 66
%             flag = flag + 1;
%             pt_sep = 1;
%         elseif P >= 90 && P <= 110
%             flag = flag + 1;
%             pt_sep = 12;
%         elseif P > 110 && P <= 160
%             flag = flag + 1;
%             pt_sep = 18;
%         elseif P > 160
%             flag = flag + 1;
%             pt_sep = 25;
%         else
%             flag = flag + 2;
%         end
        if flag == 1
            break
        end
    end
    if mod(im_num,10) == 0 || im_num == N_images
        savedata(directData,handles);
    end
    cd(ImDir)
%     if mod(im_num,20) == 0
        saveimage(handles,im_num);
%     end
    cd '..'
%     pause(2)
    hold off
    if im_num == N_images
        flag = flag + 1;
    end
    if ~isempty(varargin) && strcmp(varargin{1},'EM')
        flag = 3;
    end
end
close all
% savedata(directData,handles);

%msgbox('Tracing Complete','Success')
end

function handles = load_data(input_datastore, out_pathname, handles,im_num,paths)
    
    handles.ans_filepath = paths{im_num};
%     bslash = find(handles.ans_filepath == '\');
%     fslash = find(handles.ans_filepath == '/');
%     if isempty(fslash) && bslash(end) ~= length(handles.ans_filepath)
%         handles.filename = handles.ans_filepath(bslash(end)+1:end);
%     elseif isempty(fslash) && bslash(end) == length(handles.ans_filepath)
%         handles.filename = handles.ans_filepath(bslash(end-1)+1:end);
%     elseif isempty(bslash) && fslash(end) ~= length(handles.ans_filepath)
%         handles.filename = handles.ans_filepath(bslash(end)+1:end);
%     elseif isempty(bslash) && fslash(end) == length(handles.ans_filepath)
%         handles.filename = handles.ans_filepath(bslash(end-1)+1:end);
%     end
    pathparts = strsplit(handles.ans_filepath,filesep);
    handles.filename = pathparts{end};
    A = readimage(input_datastore,im_num); %Reduntant? just assign handles.A from start, annoying to write maybe
    if size(A, 3) > 1
        A = rgb2gray(A(:,:,1:3));
    end
    if size(A,1) > 512
%        A = A(19:530,38:549);
         A = A(24:535,47:558);
    end
    % get resolution and store it
    image_resol1 = size(A,1);
    handles.image_resol = image_resol1;
    handles.txt_resol1 = num2str(image_resol1);
    image_resol2 = size(A,2);
    handles.txt_resol2 = num2str(image_resol2);

    % normalize heights to make the image readable
    C = round((A-min(min(A)))./(max(max(A))-min(min(A)))*255-min(min(A)));

    handles.output_name = out_pathname;
    handles.A = A;
    handles.C = C;
end


function d = angdiff(th1, th2)
d = th1 - th2;
d = mod(d+pi, 2*pi) - pi;
end

% --------------------------------------------------------------------
% Main part of the algorithm, should be mostly untouched save for point
% input
% Gets points on the fibril from user input, traces the
% fibril that corresponds to these points, and fits a spline
% to the resulting points.
function handles = getpoints(handles,directData,chains,count,pt_sep)

%Adding debugging code to check the zoom on each LL--------------------
DEBUG_PLOTS = true;
PLOT_EVERY = 1;
%Added above-------------------------------------------------------------

cn = count(1);
Ntxt = count(2);

X = chains(cn).points(1,2:pt_sep:end-1);
Y = chains(cn).points(2,2:pt_sep:end-1);

if mod(length(chains(cn).points) , pt_sep) ~= 1 && pt_sep ~= 1
    X(end+1) = chains(cn).points(1,end);
    Y(end+1) = chains(cn).points(2,end);
end

if length(X) == 1
    X(end+1) = chains(cn).points(1,end);
    Y(end+1) = chains(cn).points(2,end);
end

nmperpx = handles.nmperpx;

%b-spline fit to the user defined points 
%returns spline values (sp_val) and its derivatives (dir_der)
%modified interparc fuction fit returns equially distanced points
%note dir_der is not unit length
[sp_val, dir_der] = mod_interparc(1e-6:600,X,Y,'spline'); %Doesnt look like a proper b-spline but more of a interploated cubic spline
Xsp = sp_val(:,1);
Ysp = sp_val(:,2);


%initial b_spline to the user defiend points in nm:
L2 = length(Xsp)*nmperpx*1.2;
[sp_val1, dir_der1,xxx, kappa1] = mod_interparc2(1e-6:L2,X*nmperpx,Y*nmperpx,'spline');

Xsp1 = sp_val1(:,1);
Ysp1 = sp_val1(:,2);

handles.Xsp1 = Xsp1;
handles.Ysp1 = Ysp1;
handles.dir_der1 = dir_der1;
handles.kappa1 = kappa1;

% Remove background
A = handles.A;
se = strel('disk', 12); %default 12
im2 = imtophat(A, se);

% median filter to remove speckles, use odd sizes to prevent shifting
im1 = medfilt2(im2, [3 3]);

[hh ww] = size(im1);
imcenter = im1(round(hh*0.1):round(hh*0.9),round(ww*0.1):round(ww*0.9));
THRESH = double(prctile(reshape(imcenter,1,[]), 80));
MAX = double(prctile(reshape(imcenter,1,[]), 99.9));

XT = Xsp;
YT = Ysp;

% search range
steps = 0.25; %default 0.2
trange = -6:steps:6; %default 3
prange = -2:steps:2; %default -2

tr=repmat(trange,length(prange),1);
pr=repmat(prange',1,length(trange));


lwcv = zeros(0,0,0); % value (score) at a specific length, width & centre compare to where the initial spline is.

for LL=1:length(Xsp)
  % tangent direction
    dx = dir_der(LL,1); % dir-der from mod_interparc
    dy = dir_der(LL,2);
    dl = sqrt(dx^2+dy^2);
    dx = dx / dl; % normal tangent vector  
    dy = dy / dl;



    
%Adding debugging code=====================================================
    if DEBUG_PLOTS && mod(LL,PLOT_EVERY)==0
        % radius of the visual crop ~ search extents
        Rn = max(abs(trange)); 
        Rp = max(abs(prange));
        R  = ceil(1.5*(Rn+Rp));
    
        cx = Xsp(LL);  cy = Ysp(LL);
        x1 = max(1, round(cx - R));  x2 = min(size(im1,2), round(cx + R));
        y1 = max(1, round(cy - R));  y2 = min(size(im1,1), round(cy + R));
    
        patchImg = im1(y1:y2, x1:x2);
    
        figure(101); clf
        imagesc([x1 x2], [y1 y2], patchImg); 
        axis image ij; colormap gray; hold on
        plot(cx, cy, 'r+', 'MarkerSize', 8, 'LineWidth', 1.5);
        % show local frame: tangent (green) and normal (cyan)
        quiver(cx, cy, dx*R,  dy*R,  0, 'g', 'LineWidth', 1.2);
        quiver(cx, cy, dy*R, -dx*R,  0, 'c', 'LineWidth', 1.2);
        title(sprintf('Pre-grid zoom @ LL=%d  (center=(%.1f, %.1f))', LL, cx, cy));
    end

%End first debug block=====================================================




    
    % can use meshgrid instead
    Px = Xsp(LL) + tr * dy + pr * dx; % x corrdinates of search area (grid) 
    Py = Ysp(LL) + tr * -dx + pr * dy;
    
    % remove out-of-bound points
    Px(Px<1)=1;
    Py(Py<1)=1;
    Px(Px>size(im1,1)-1)=size(im1,1)-1; % repeats border px for out of bound coords
    Py(Py>size(im1,2)-1)=size(im1,2)-1;
    
    pts = [];
    %         pts = get_subpixel(im1, [Px ; Py]'-1, 'cubic')';
    
    % % TODO: remove this part of the code, sub pixel calculations
    Pxflat = reshape(Px,1,[]); %grid x points into a list
    Pyflat = reshape(Py,1,[]);
    % NRtest 
    PxPyflat = num2cell([Pxflat ; Pyflat]',1); PxPyflat = PxPyflat([2,1,3:end]);
    pts = interpn(double(im1),PxPyflat{:},'cubic');
%     pts = pkget_subpixel(im1, [Pxflat ; Pyflat]', 'cubic')'; 
%     pts = pkget_subpixel(im1, [Pxflat ; Pyflat]', 'nearest')'; 
    pts = reshape(pts,length(prange),length(trange)); % intensities at the grid points




% 2nd debug block=========================================================
    if DEBUG_PLOTS && mod(LL,PLOT_EVERY)==0
        % overlay the actual queried subpixel grid on the same crop
        figure(101); hold on
        plot(Px(:), Py(:), '.', 'MarkerSize', 4);  % sampling lattice
        drawnow;
    
        % also visualize the sampled intensity patch in (prange x trange) coords
        figure(102); clf
        imagesc(trange, prange, pts);
        axis image; axis xy; colormap gray
        xlabel('normal offset t (px)'); ylabel('tangent offset p (px)');
        title(sprintf('Sampled patch (pts) @ LL=%d', LL));
        drawnow;
    end

% End 2nd block===========================================================




        
    if LL == 20
        wcv = centroid06(pts, 0, THRESH, MAX);
    else
        wcv = centroid06(pts, 0, THRESH, MAX);
    end
        
    
    if (max(size(lwcv))==0)
        lwcv=zeros(length(Xsp),size(wcv,1),size(wcv,2));
    end
    lwcv(LL,:,:,:) = wcv;
        
end

lwcv(lwcv<0.01)=0.01;

% chain points
chpX = zeros(length(Xsp),2);
chpY = zeros(length(Xsp),2);

% find widths with max values (scores) along chain
wl = zeros(1,length(Xsp));
for LL=1:length(Xsp)
    wcv = squeeze(lwcv(LL,:,:));
    
    %findpeaksn finds local maxima and ignores nearby peaks
    [cw,cc]=ind2sub(size(wcv), find(findpeaksn(wcv,true(1, ndims(wcv)), 0.1, 5)));
    if isempty(cw) 
        warning('no peaks at length=%d', LL);
        wl(LL) = 3/steps; % NRtest wl(LL-1);
        continue;
    end
    wl(LL) = mean(cw);
    
end
% smooth widths, moving average
norm_filter = fspecial('gauss', [1 10], 3);
wl = imfilter(wl, norm_filter, 'replicate');
for LL=1:length(Xsp)
    % TODO:: refactor this for faster speed
    dx = dir_der(LL,1);
    dy = dir_der(LL,2);
    dl = sqrt(dx^2+dy^2);
    dx = dx / dl;
    dy = dy / dl;
    
    % give more weight to those near the smooth width at that point of the
    % chain
    wcv = squeeze(lwcv(LL,:,:));
    W = repmat(normpdf(1:size(wcv,1), wl(LL), size(wcv,1)*0.1)', 1, size(wcv,2));
    W = W/max(max(W));
    wcv = W.*wcv;
    
    % more weight for points in the same direction as the last point
    if LL>2
        tdx = XT(LL-1) - XT(LL-2);
        tdy = YT(LL-1) - YT(LL-2);
        tdl = sqrt(tdx.^2+tdy.^2);
        tdx = tdx / tdl;
        tdy = tdy / tdl;
        
        xn = XT(LL-1) + tdx;
        yn = YT(LL-1) + tdy;
        
        Wc=[];
        for cc=1:size(wcv,2)
            c = trange(cc);
            
            x2 = Xsp(LL) + dy * c;
            y2 = Ysp(LL) + -dx * c;
            L2 = (x2-xn).^2+(y2-yn).^2;
            Wc(cc) = 5/(5+L2);
        end
        W=repmat(Wc, size(wcv,1), 1);
        wcv = W.*wcv;
    end
    
    [cw, cc] = ind2sub(size(wcv), find(findpeaksn(wcv,true(1, ndims(wcv)))));
    if length(cw) == 0
        error('No peaks found')
    end
    cv=[];
    for ii=1:length(cw)
        cv(end+1)=wcv(cw(ii),cc(ii));
    end
    %[cw';cc';cv]
    [a,k]=max(cv);
    ii=cc(k);
    ww=cw(k);
    vv=cv(k);
    
    t = trange(ii);
    
    XT(LL) = Xsp(LL) + dy * t;
    YT(LL) = Ysp(LL) + -dx * t;
    
    XJ = [XT(LL)-dy*steps*ww/2 XT(LL)+dy*steps*ww/2];
    YJ = [YT(LL)+dx*steps*ww/2 YT(LL)-dx*steps*ww/2];
    chpX(LL,:) = XJ;
    chpY(LL,:) = YJ;

end

% NRtest
cm = lines;
rgb = squeeze(ind2rgb(round(vv * 255), cm));

Xsp = XT;
Ysp = YT;

%pause(0.5);


figure(1)
hold on
plt = plot(chpX(:,1), chpY(:,1), '--b', 'LineWidth', 1.5);
plt = plot(chpX(:,2), chpY(:,2), '--b', 'LineWidth', 1.5);

XT=XT'; %the final traced points on the chain
YT=YT';

handles.X = X;
handles.Y = Y;
handles.XT = XT;
handles.YT = YT;

% ---------------------------------------------------------------
% --- Compile spline fitting and image data (--> normalize lengths vs. pixels)
%---- Executes on button press in btn_compile_and_check.

%nmperpx = image_width*1000/image_resol % lengths are converted from microns to nanometers (*1000)

DOTS = [XT ; YT];

% DOTS = DOTS(:,find(sum(diff(DOTS').^2,2)~=0));
L2 = length(Xsp)*nmperpx*1.2;
% [sp_val, dir_der] = mod_interparc(1e-6:L2,DOTS(1,:)*nmperpx,DOTS(2,:)*nmperpx,'spline'); %note dir_der is not unit length

%return curvature, derivatives and spline
[sp_val, dir_der, xxx, kappa] = mod_interparc2(1e-6:L2,DOTS(1,:)*nmperpx,DOTS(2,:)*nmperpx,'spline');

x_norm = sp_val(:,1);
y_norm = sp_val(:,2);

Xsp1 = handles.Xsp1;
Ysp1 = handles.Ysp1;

figure(1)
hold on
% plt = plot(X, Y, 'ko','MarkerSize',1,'linewidth',.1);
plt = plot(Xsp1/nmperpx, Ysp1/nmperpx, 'k--','linewidth',1);
plt = plot(x_norm/nmperpx, y_norm/nmperpx, 'r-','linewidth',1.5);
txt = num2str(cn+Ntxt);
text(Xsp1(1)/nmperpx,Ysp1(1)/nmperpx,txt,'HorizontalAlignment','right','Color','g');

handles.dir_der = dir_der;
handles.kappa = kappa;

% plt = text(max(x_norm)+5, max(y_norm)+5, num2str(chain_n), 'Color', 'black', 'BackgroundColor', 'white', 'Margin', 0.1, 'fontsize', 14);
% handles.chain_plots{chain_n,plot_nn} = plt;
% plot_nn = plot_nn+1;


bspline_norm_coord(1,:) = x_norm;
bspline_norm_coord(2,:) = y_norm;
handles.bspline_norm_coord = bspline_norm_coord;


handles.x_norm = x_norm; % used in interval/check increment function
handles.y_norm = y_norm;

handles.nmperpx = nmperpx;

arclen = arclength(x_norm,y_norm);

field = 'chain';
tpoints = [0:length(XT)-1];
value = {XT;YT;tpoints};
tempStruct = struct(field,value);
directData = [directData,tempStruct];



% ---------------------------------------------------------------
% Check knots spline increments over the fibril
%INTERVAL FUNCTION
x_norm = handles.x_norm;
y_norm = handles.y_norm;

for i=1:length(x_norm)-1
    dist = ( x_norm(i) - x_norm(i+1) ).^2 + ( y_norm(i) - y_norm(i+1) ).^2;
    intervals(i) = sqrt(dist);
end

% average_interval between each knot over all the spline
mean_itv = mean (intervals(1,:));
mean_itv = round (mean_itv);
% standard deviation
sd_itv = std(intervals(1,:));
sd_itv = round (sd_itv);


handles.mean_itv = mean_itv;
handles.sd_itv = sd_itv;
handles.intervals = intervals;


end


function handles = combo_calc(handles,n)

%TANTAN-COREL
%--------------------------------------------------------------------
%-------Calculate the decay of tangent-tangent correlations over the fibril
x_norm = handles.x_norm;
y_norm = handles.y_norm;
nmperpx = handles.nmperpx;

arclen = arclength(x_norm,y_norm);

seps_mean = [];
cosines_mean = [];
cosines_std = [];

R2s_mean = [];
R2s_std = [];

cosines = [];
seps = [];
R2s = [];

angles2pi = [];
seps_angles = [];

thetas = [];


separations = 1:1:max(arclen);
for sep = separations
    
    %LL = 0:sep:max(arclen);
    %[x_points, y_points, xder,yder, idx_points] = point_at_length(handles,LL);
    
    x_points = handles.x_norm(1:sep:end);
    y_points = handles.y_norm(1:sep:end);
    der = normr(handles.dir_der(1:sep:end,:));
    
    %correlation (cos(tetha)) calculations
    %der = [xder, yder];
    
    raw_theta = (atan2(der(:,2),der(:,1)));
    
    theta = abs(diff(raw_theta));
    thetas = [thetas;theta];
    
    angle = abs(mod(diff(raw_theta)+pi,2*pi)-pi);
    angles2pi = [angles2pi;angle];
    seps_angles = [seps_angles;ones(length(angle),1)*sep];
    
    cosine = dot(der(1:end-1,:),der(2:end,:),2);
    cosines = [cosines;cosine];
    seps = [seps;ones(length(cosine),1)*sep];
    cosines_mean = [cosines_mean;mean(cosine)];
    cosines_std = [cosines_std;std(cosine)];
    seps_mean = [seps_mean,sep];
    
    %end-to-end (R^2) calculations
    rr = [x_points,y_points];
    R2 = sum(diff(rr).^2,2);
    R2s = [R2s;R2];
    
    R2s_mean = [R2s_mean;mean(R2)];
    R2s_std = [R2s_std;std(R2)];
end

handles.angles2pi = angles2pi;
handles.cosines = cosines;
handles.contours = seps;
handles.thetas = thetas;

corel_single = [cosines,seps];
handles.corel_single = corel_single;

angles2pi_single = [angles2pi,seps_angles];
handles.angles2pi_single = angles2pi_single;

thetas_single = [thetas,seps_angles];
handles.thetas_single = thetas_single;

handles.end2end = R2s;

end2cont_single = [R2s,seps];
handles.end2cont_single = end2cont_single;

goodcontour = arclen;
handles.goodcontour = goodcontour ;

handles.total_end2end = sqrt (( y_norm(end) - y_norm(1) ).^2 +...
    (( x_norm(end) - x_norm(1) ).^2 ));


%--------------------------------------------------------------------
%BELOW is the code that fill in the structure with all the elements, pertaining
% to one fibril, calculated above in this very same function

% IMPORTANT to keep this

handles.mat_length_struct(n).contourL = handles.goodcontour;
handles.mat_intervals_struct(n).meanitv_fib = handles.mean_itv;
handles.mat_sd_struct(n).meansd_fib = handles.sd_itv;

handles.bsplines_struct(n).bspline = handles.bspline_norm_coord;

handles.correlations_struct(n).corelfib = handles.corel_single;
handles.wormlike_struct(n).wormfib = handles.end2cont_single;

handles.length_conformation(n, 1) = handles.goodcontour;
handles.length_conformation(n, 2) = handles.total_end2end;

handles.points(n).X = handles.X;
handles.points(n).Y = handles.Y;
handles.points(n).XT = handles.XT;
handles.points(n).YT = handles.YT;
handles.points(n).filename = handles.filename;
handles.points(n).filepath = handles.ans_filepath;

handles.der(n).Xder = handles.dir_der(:,1);
handles.der(n).Yder = handles.dir_der(:,2);

handles.curv(n).curvature = handles.kappa;

handles.angles2pi_struc(n).angles2pi = handles.angles2pi_single;
handles.thetas_struc(n).thetas = handles.thetas_single;

handles.info_struc(n).image_width = nmperpx;
handles.info_struc(n).image_resol = handles.image_resol;

handles.userspline(n).Xsp1 = handles.Xsp1;
handles.userspline(n).Ysp1 = handles.Ysp1;
handles.userspline(n).dir_der1 = handles.dir_der1;
handles.userspline(n).kappa1 = handles.kappa1;


% Set the string and color of text to confirm chain has been added
%disp('Chain added!');
end


% --- This just needs to run once automatically at the end
% --- Executes on button press in btn_savedata.
function savedata(directData, handles)

% Important Stuff
contour_lengths = handles.mat_length_struct;
intervals_means = handles.mat_intervals_struct;
intervals_sd = handles.mat_sd_struct;

bsplines_norm = handles.bsplines_struct;

mat_correlations = handles.correlations_struct;
mat_wormlike = handles.wormlike_struct;

length_conformation = handles.length_conformation;
points = handles.points;
der = handles.der;
angles2pi = handles.angles2pi_struc;
thetas = handles.thetas_struc;

curv = handles.curv;

info = handles.info_struc;

userspline = handles.userspline;
save('dataForVectorAnalysis.mat','directData');
output_base = handles.output_name;
output_sup = '_SmTr';
updated_filename=[output_base, output_sup];

%newfilename = fullfile(pathname, updated_filename);
save(updated_filename,'contour_lengths','intervals_means',...
    'intervals_sd','bsplines_norm',...
    'mat_correlations','mat_wormlike','length_conformation','points','angles2pi','thetas', 'info', 'der','curv','userspline');

end


% --- Have this auto run once per image
% --- Executes on button press in btn_saveimage.
function saveimage(handles,n)
% hObject    handle to btn_saveimage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% fname = [ 'SmarTrace_', handles.output_name, '_image', num2str(n)];
 fname = ['Traced_', handles.filename];
saveas(handles.figure, fname)
saveas(handles.figure, fname,'png')
end
% --- Should not need to touch this
function [x,y,xder,yder,idx] = point_at_length(handles, LL)

x_norm = handles.x_norm;
y_norm = handles.y_norm;
dir_der = handles.dir_der;

    if ~isfield(handles, 'length_on_chain')

        x_diff = diff(x_norm);
        y_diff = diff(y_norm);

        r_diff = sqrt(x_diff.^2 + y_diff.^2);

        handles.length_on_chain = cumsum(r_diff);
    end
length_on_chain = handles.length_on_chain;

idx = [];
    for ll2=LL
        idx(end+1) = find(length_on_chain>ll2,1);
    end
x = x_norm(idx);
y = y_norm(idx);
xder = dir_der(idx,1);
yder = dir_der(idx,2);

end

end
