% prototype for retroconversion of picked particles (subtomo) to segmentations and model volumes
% hacky mess for I/O right now, put everything in a working folder
% have a /tomos folder for inputs(just for size/resolution really), /densities folder for particle maps
% and a folder for picked .star files. auto-fetches identical X.mrc densities from X.star selections.
% outputs atlas files into a /segmentations folder, along with a record of which labels each output includes

%% iter 2 for chlamy picked particles to segmentation script

% needed variables:
% bin of source (currently 4)
% shrink val (to replace water density opposing creep of particles outward)
% controls for separate inputs of tomograms?

% select parent folder, assume a specific structure
workingfolder = uigetdir();

%% select target star files to process
[list, path] = uigetfile({'*.star'},'Select coordinate file(s)',workingfolder,'MultiSelect','on');
clear starfiles
if ~iscell(list), starfiles{1}=list; else starfiles=list; end
% pre-extract tomo IDs to make it faster? not sure how much time saved, probably not a lot
%might just make code cleaner later on, fewer checks needed?
clear rec
for i=1:numel(starfiles)
    [~,base,ext] = fileparts(starfiles{i});
    star = readstar(fullfile(path,starfiles{i}));
    
    densmrc = fullfile(workingfolder,append('densities\',base,'.mrc'));
    [~,dat] = densityprep(densmrc);
    % density might need work, the trim process might be shifting away from 0-centering
    rec(i).name = base; rec(i).dat = star; rec(i).atoms = dat;
end
clear star starfiles dat list

%% fetch tomograms from fixed folder
list2 = dir(append(workingfolder,'\tomos\*.mrc'));
tomolist = cell(numel(list2),1);
for i=1:numel(list2)
    tomolist{i} = append(workingfolder,'\tomos\',list2(i).name);
end
clear list2

%% loop over tomos and generate segmentations
clear vol split picks pickrecord idx atlas tmp t
%tmp = randn(3,1000); tmp = 25*(tmp./vecnorm(tmp))'; tmp(:,4) = 10;
%rec(3).atoms = tmp;
outfolder = append(workingfolder,'\segmentations\');
mkdir(outfolder)
pickrecord = zeros(numel(tomolist),numel(rec),'single');
for i=1:numel(tomolist)
    [~,tomobase,~] = fileparts(tomolist{i});
    tomobase = pad(tomobase,4,'left','0');
    [~,srchead] = ReadMRC(tomolist{i}); srcsize = [srchead.nx,srchead.ny,srchead.nz];
    
    for j=1:numel(rec)
        split.(rec(j).name) = zeros(0,4,'single');
        idx = contains(rec(j).dat.tomo,tomobase(1:end));
        pickrecord(i,j) = single(any(idx)); %array to record which tomos have particles from which classes
        if any(idx>0)
            picks = rec(j).dat(idx,:); % particles picked of class in the tomogram
            n = size(picks,1); m = size(rec(j).atoms,1);
            split.(rec(j).name) = zeros(n*m,4,'single');
            
            for k=1:n
                co = (picks.coords(k,:)+0)*srchead.pixA/4/1; %convert bin4 pixels to angstroms (need variable)
                rot = deg2rad(picks.rots(k,:)); tmp = rec(j).atoms;
                
                [mat,tform] = eul2rot(rot);
                tmp(:,1:3) = transformPointsForward(affine3d(tform),tmp(:,1:3))+co; %rotate and translate
                split.(rec(j).name)(1+(k-1)*m:m*k,:) = tmp;
            end
        end
    end
    [vol,atlas] = trunc_atoms2vol(srchead.pixA,split,srcsize*srchead.pixA/1);
    outf = append(outfolder,tomobase,'_seg.mrc');
    WriteMRC(atlas,srchead.pixA,outf,0)
    fprintf('%i,',i); if rem(i,25)==0; fprintf('\n'); end
end
fprintf('\n')
t = table(tomolist,pickrecord,'VariableNames',{'source tomo','particle'});
t = splitvars(t,'particle','NewVariableNames',{rec.name});
writetable(t,append(outfolder,'atlas_records.csv'))


%% internal functions
function [dtrim,dat] = densityprep(mrc)
[dvol,dhead] = ReadMRC(mrc);

dvol(dvol<mean(dvol,'all')) = 0;
trimr = any(dvol,[2 3]); trimc = any(dvol,[1 3]); triml = any(dvol,[1 2]); 
dtrim = single(dvol(trimr,trimc,triml)); % zeros trimmed to shrink volume
%dtrim = dvol;
dat = zeros(numel(dtrim),4,'single');
[xx,yy,zz] = ind2sub(size(dtrim),1:size(dat,1));
dat = [xx',yy',zz',dtrim(:)];

ix = dat(:,4)>(mean(dat(:,4))/2); dat = dat(ix,:); %coordinate array w/ associated density
n = size(dat,1); ix = randperm(n);
ix = ix(1:round(n/10*dhead.pixA));dat = dat(ix,:); % randomly keep pix/10 of pts for faster run
cen = size(dtrim)/2-0.5; dat(:,1:3) = (dat(:,1:3)-cen)*dhead.pixA;
end

function star = readstar(starfile)
text = textscan(fileread(starfile),'%s','delimiter','\n'); text = text{1};

headstart = find(strncmp(text,'_rlnCoordinateX #1',15)); %header id start
headend = find(strncmp(text,'_rlnParticleName #9',15)); %header id end
head = erase(text(headstart:headend),{'_rln',' #'}); % extract header lines, prune prefixes

text = text(headend+1:end); % extract data, after header lines
text(cellfun(@isempty,text)) = []; %remove empty cells, leaving only data

q = cellfun(@(x) strsplit(x, '\t'), text, 'UniformOutput', false); %split cells by tab
q = vertcat(q{:}); % concatenate into single array of cells
dat = cellfun(@str2double,q(:,1:6));
srcix = contains(head,'TomoName'); % which col has tomo name, because it's not static
%tm = cellfun(@convertCharsToStrings,q(:,srcix)); % super slow
tm = cellfun(@(x) x(1:end),q(:,srcix), 'UniformOutput', false);
star = table(dat(:,1:3),dat(:,4:6),tm,'VariableNames',{'coords','rots','tomo'});
%f = contains(q(:,srcix),tomo); q = q(f,:); % filter down to records matching input tomogram
%if size(q,1)<1, error('looks like this tomogram has no picked particles'), end
end

function [mat,tform] = eul2rot(e) %currently hardcoded zyz convention - relion standard, EMAN2 is zxz
r1 = axrot2mat([0,0,1],e(1));
r2 = axrot2mat([0,1,0],e(2));
r3 = axrot2mat([0,0,1],e(3));
mat = r1*r2*r3;
tform = eye(4); tform(1:3,1:3) = mat;
end

function [mat,tform] = axrot2mat(ax, rad) %axis angle to 3x3 rotation matrix
  bw = [0, -ax(3), ax(2); ax(3), 0, -ax(1); -ax(2), ax(1), 0];
  mat = eye(3) + sin(rad)*bw + (1-cos(rad))*bw*bw;
  tform = eye(4); tform(1:3,1:3) = mat;
  %tform = makehgtform("axisrotate",ax,rad); mat = tform(1:3,1:3);
end

function [vol,atlas,split] = trunc_atoms2vol(pix,pts,sz,offset)
if isstruct(pts)
    names = fieldnames(pts); pts = struct2cell(pts); %convert to cell for easy looping
else
    names = 0;
end
if ~iscell(pts) && numel(size(pts))==2
    pts = {pts}; %convert single array to cell for ease of use
end

s = numel(pts); %temp patch for old loop definition code and cell/array switcher

if nargin<3 %no box inputs, output tight bounds
    dd = vertcat(pts{:});
    offset = min(dd(:,1:3),[],1)-pix;
    sz = max(dd(:,1:3),[],1)+pix-offset;
elseif all(sz==[0,0,0]) %if box 0, keep the object centered in the volume
    dd = vertcat(pts{:});
    offset = -max(abs(dd(:,1:3)))-pix/2;
    sz = offset*-2;
elseif nargin<4 %box limit only, output corner starting at 0
    offset = [0,0,0];
end
% rough constants - need improved values, per-atom vol especially
%avol = 4/3*pi*(1.9^3); %eyeballed volume of the average organic atom (radii approx 1.8A)- get per-atom measure?
%h20 = 3.041/2; %computed scatter factor for H2O - /2 for similarity to vol and simulate defaults
%wd = 6.022e23/18/(1e8)^3; %molecules of water per a^3 - ~1/30 for liquid water
%wvol = 32; %eyeballed volume of amorphous ice molecules in angstroms

emsz = floor(sz/pix); 

%solv = (rand(emsz)-0.6)*1.5*pix^2+(pix^3); %set initial solvent density
%solv = imgaussfilt3(solv,0.5); % smoother solvation test
%solv = 0;
%acount = zeros(emsz,'single');
sptmp = cell(1,s);
for j=1:s
    p = single(pts{j});
    
    if size(p,2)<4, p(:,4)=1; end %intensity==1 if not provided in 4th column
    mag = p(:,4); p = p(:,1:3); p = round( (p-offset)/pix+0.5 );
    
    [sptmp{j}] = internal_accumarray(p,mag,emsz);
end

vol = zeros(emsz,'single');
atlas = vol;
for i=1:numel(sptmp)%:-1:1
    hit = sptmp{i}>vol; 
    atlas(hit)=i;
    vol = vol+sptmp{i};
end

%{
tmp = cat(4,zeros(emsz,'single'),sptmp{:}); %memory limitations
[~,atlas] = max(tmp,[],4); atlas = atlas-1; % memory limitations
atlas = single(atlas);
vol = sum(tmp,4);
%}

if iscell(names)
    for i=1:s
        split.(names{i}) = sptmp{i};%tmp(:,:,:,i+1);
    end
elseif s>0
    split = cell(1,s);
    for i=1:s
        split{i} = sptmp{i};%tmp(:,:,:,i+1);
    end
else
    split = 0;
end
end

function [tmpvol] = internal_accumarray(p,mag,emsz)
if nargin<6, acount = zeros(emsz,'single'); end
ixf = ones(size(mag),'single');
for i=1:3
    ix = p(:,i) <= emsz(i) & p(:,i) >= 1; %index points inside the box
    %p = p(ix,:); mag=mag(ix); %drop points outside the box
    ixf = ixf & ix; %ixf = ixf.*ix;
end
ixf = ixf>0;
p = p(ixf,:); 
mag = single(mag);
mag=mag(ixf); %out of memory error - problematic reuse of same array? is mag too high precision?
tmpvol = accumarray(p,mag,emsz);
end