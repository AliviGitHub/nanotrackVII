% Last updated by: Michelle Crook on 07/2019, at University of
% California,Berkeley
     %Adapted from Matt Hauwiller
     %Cropping code adapted from Amy McKeown-Green

% VARIABLES TO CHECK FOR EACH ANALYSIS %
fps=1;
x=1; %number of frames to average over


%%%CROPPING CODE%%%
%gets cropping dimensions from DM
[filename, pathname] = uigetfile('*.*', '*SELECT TEXT FILE WITH CROPPING DIMENSIONS*', 'Select text file with cropping dimensions'); %Obtain desired file and pathname
fullpathname = strcat(pathname, filename); %Create full file pathname

%%Open file and scan
fileID = fopen(fullpathname, 'rt'); %Opening user selected file
fileData = textscan(fileID, '%s', 'delimiter', '\n', 'whitespace', ' '); % Reading information inside a file

%%Isolate number data and convert from cell to double
TotHeadInfoStrcell = cellstr(fileData{1,1}(1:1));
cellstrdimensions = (strsplit(TotHeadInfoStrcell{1,1},' ')); % Cuts string based on presence of " "
ROIdimensions = str2double(cellstrdimensions); % converts cell array to an array of doubles

%%Assign dimensions to variables (may seem redundant but reduces overhead for parfor loop)
left = ROIdimensions(1,1);
right=ROIdimensions(1,2);
down=ROIdimensions(1,3);
up=ROIdimensions(1,4);

%selects where data will be saved 
wheretosavename = uigetdir('','*SELECT LOCATION TO SAVE DATA*'); 


%%%CREATE STRUCTURE TO STORE ALL FRAMES FROM DATA TREE%%%
start_min_path=fullfile(matlabroot, '\Users\chellecrook\Documents\Grad School\TEM\MatLab Code\Hour_00');
hour00=uigetdir(start_min_path, 'SELECT HOUR_00 FOLDER OF DATASET');
minFolders=dir(hour00);
minFolders=minFolders(~ismember({minFolders.name},{'.','..'}));
numberOfMinFolders = length(minFolders);

allSubFolders = genpath(hour00);

remain = allSubFolders;
listOfFolderNames = {};
while true
	[singleSubFolder, remain] = strtok(remain, ';');
	if isempty(singleSubFolder)
		break;
	end
	listOfFolderNames = [listOfFolderNames singleSubFolder];
end

for a=1:numberOfMinFolders-1
vid(1).min=listOfFolderNames(:,[2 a+61]); 
end

listOfSecFolders=listOfFolderNames;
lngthSecFolderName=max(length(listOfFolderNames{1,3}));
listOfSecFolders(cellfun('length',listOfSecFolders)<lngthSecFolderName)=[];
numberOfSecFolders=length(listOfSecFolders);
%%
vid(1).sec=listOfSecFolders;
chr='\';

for bb= 1:numberOfSecFolders
    frames=listOfSecFolders{:,bb};
    frameb=dir(frames);
    frameb=frameb(~ismember({frameb.name},{'Thumbs.db','.','..'}));%to remove thumbs.db, '.', '..' files
    frameFolder=frameb.folder;
    for l=1:length(frameb)
        frameName=frameb(l,:).name;
        frameFull=strcat(frameFolder,chr,frameName);
        vid(1).frame(l,:)=frameFull;
    end
    vid(1).compframe((fps*(bb-1))+1:fps*bb,:)=vid(1).frame;
end

corefilenamecell = (strsplit(frameName,'_Hour')); %cuts it so there is just the core file name (no hour, minute, or second data)
corefilenamestr = corefilenamecell{1,1};

%%%DETERMINE TIME ZERO FOR ANALYSIS(when dose rate is increased)%%% 
    counts=zeros(fps*7,1);
    for       i = 1:(fps*7)
    fr=vid(1).compframe(i,:);
    frs=dmread(fr);
    frvec=frs.ImageList.Unnamed0.ImageData.Data.Value;
    xdim =frs.ImageList.Unnamed0.ImageData.Dimensions.Unnamed0.Value;
    ydim =frs.ImageList.Unnamed0.ImageData.Dimensions.Unnamed1.Value;
    pic=reshape(frvec,[xdim, ydim]);
    
    avCount=mean(mean(pic));
    counts(i,:)=avCount;
    end
    
    TF=ischange(counts,'MaxNumChanges',1);
    if TF==0
        nstart=1;
    else
        nstart=find(TF)+1;
    end
%%    
%define the location struct
locations=struct('ftr',{});
f=0;

%Conversion of Pixel to nm (changes for each mag)
ptonm=284.25/2048; %97KX 

%%%START MAIN LOOP OF ANALYSIS%%%
 
 for       i = nstart:x:length(vid(1).compframe)       
 
     %reads in frames from vid structure and reshapes pixels into appropiate
%sized image
    
    for ii=nstart:nstart+x-1   
        fr=vid(1).compframe(ii,:);
        frs=dmread(fr);
        frvec=frs.ImageList.Unnamed0.ImageData.Data.Value;
        xdim =frs.ImageList.Unnamed0.ImageData.Dimensions.Unnamed0.Value;
        ydim =frs.ImageList.Unnamed0.ImageData.Dimensions.Unnamed1.Value;
        pic=reshape(frvec,[xdim, ydim]);
        picture=transpose(pic);
        
        picCell=mat2cell(picture,xdim,ydim);
        picArray(1,(ii-nstart+1))=picCell;
    end
    
    centrd(1,1)=xdim/2;
    centrd(1,2)=ydim/2;
    %change nstart so that it will loop over all frames
    nstart=nstart+x;
    %sums x number of frames for better contrast
    sumPic=sum(cat(3,picArray{:}),3); 
    
    %crops image
    imagecropped = sumPic(left:right,down:up); % crops image 
    imagegrysc=uint16(65536*mat2gray(imagecropped)); %normalizes --> max pixel value = 1 and min =0 --> converts to integer necessary for tiff
    imagefiletiff= strcat(wheretosavename,'\', corefilenamestr, '_edit','\', corefilenamestr,'_edit_', num2str(i),'.tif'); %writes new file name
    cropFoldName=strcat(corefilenamestr, '_edit');
    mkdir(fullfile(pathname, cropFoldName));
    imwrite(imagegrysc, imagefiletiff); %writes a tiff file   
    
    %bins image to give better contrast
    binim=imresize(imagegrysc,[1024 1024]);
    
    % Read in current image/frame; reverse the image since image process
    % toolbox identifies bright areas as features
    maxint= max(max(max(imagegrysc)));
    frame = maxint-imagegrysc;
    frame8=im2uint8(frame);
    
    % show up the figure frame to check the accuracy of tracking
    figure
    imagesc(frame);colormap(gray);
    drawnow;
    
    % Use an averaging filter to remove noise 
    % NOTE: averaging and disk filters work the best
    % for this applicatiodn
    %H = fspecial('average',5);
    
    % Temporary variable for image processing
    %ftemp = imfilter(frame,H,'replicate');
    
    %Convert image to black and white
    %takes 8bit image input
    Kmean = imsegkmeans( frame8,2); %K-means
    BW = false(size(Kmean));
    if mean(frame8(Kmean==1))>mean(frame8(Kmean==2)) %assign 0s and 1s
       BW(Kmean==1) = true;
    else
       BW(Kmean==2) = true;
    end
    
    % Fills in brightest pixel clusters
    bw = imfill(BW,'holes');
    
    
    se=strel('disk',5);
    BW=imopen(bw,se);
    
    % Keep only pixel clusters that fall between mincluster and maxcluster
    %bw=xor(bwareaopen(fbw,mincluster), bwareaopen(fbw,maxcluster));
  
    % Identify the boundaries of candidate particles
    [B,L] = bwboundaries(BW);
    
    %delete any outlined background noise
    %for ii=1:length(B)
        %[M,N]=size(B{ii,:});
        %if M<100
            %B{ii,:}=[0,0];
        %end
    %end
    
    % extract all the relevant information
    point = regionprops(L,'Centroid','Area','Orientation','MajorAxisLength','MinorAxisLength','Perimeter');
    
    hold on;
        for k = 1:length(B);
               boundary = B{k};
               plot(boundary(:,2), boundary(:,1), 'b', 'LineWidth', 3);
               drawnow;
        end
        f=f+1;
        hold off;
         
      if isempty(B)==1
         continue  
     end
%%    
    % create the struct for the calculation of ends
   % feature = struct('ftr',{});
    % number of candidate particles
    npmax = size(B,1);
   
    for index=1:npmax
        locations(i).ftr(index,1:7)=[point(index,1).Centroid(1,1) point(index,1).Centroid(1,2) point(index,1).Area(1,1) point(index,1).Orientation(1,1) point(index,1).MajorAxisLength(1,1) point(index,1).MinorAxisLength(1,1) point(index,1).Perimeter(1,1)];
        Areas(1,index)=point(index,1).Area(1,1);
        Data(index,1)=point(index,1).Centroid(1,1);
        Data(index,2)=point(index,1).Centroid(1,2);
        Data(index,3)=point(index,1).Area(1,1);
    end
    
    Data=sortrows(Data,-3);
    
%     assignin('base','Data',Data)
    %if i==1;

%         centrd(1,1)=Data(1,1);
%         centrd(1,2)=Data(1,2);
   % end
    assignin('base','centrd',centrd)
     %This is the "screening" section
    
     b=0;%Use b as the checkpoint to determine whether another iteration is necessary
    while b<0.5
        s=size(Data(:,3)); %Get the size of the remaining matrix

        %The first possibility is that the particle has fully
        %disappeared in the frame, and there are no elements in the
        %matrix fitting the correct centroid (screener has deleted all
        %elements)  In this case, just put in all zeros and end the
        %screening process
        if s(1)<1.5
            s(1);
            for g=1:7
                c(1,g)=0;
            end
            b=1;
             %If the centroid of the line is the same as the centroid of
        %particle we are looking for, the screener can exit
        
        elseif Data(1,1)<centrd(1,1)+40 && Data(1,1)>centrd(1,1)-40 && Data(1,2)<centrd(1,2)+40 && Data(1,2)>centrd(1,2)-40
            b=1;

        %For all other cases, delete the row in question and repeat
        %screener loop
        else
             Data(1,:)=[];
            b=0;
        end
    end
   
        centrd(1,1)=Data(1,1);
        centrd(1,2)=Data(1,2);    
        
%     assignin('base','Areas',Areas)
    m=Data(1,3);
    [m1,n1]=find(Areas==m);
    boundarynp=B{n1};
    Outlines{i}=boundarynp;
    index=n1(1);
    endpt(1,2)= point(index,1).Centroid(1,1)+cos(point(index,1).Orientation(1)*pi/180)*point(index,1).MajorAxisLength/2;
    endpt(1,1)= point(index,1).Centroid(1,2)-sin(point(index,1).Orientation(1)*pi/180)*point(index,1).MajorAxisLength/2;
    endpt(2,2)= point(index,1).Centroid(1,1)-cos(point(index,1).Orientation(1)*pi/180)*point(index,1).MajorAxisLength/2;
	endpt(2,1)= point(index,1).Centroid(1,2)+sin(point(index,1).Orientation(1)*pi/180)*point(index,1).MajorAxisLength/2;
    endpt(3,2)= point(index,1).Centroid(1,1);
    endpt(3,1)= point(index,1).Centroid(1,2);
    endpt(4,2)= point(index,1).Centroid(1,1)+sin(point(index,1).Orientation(1)*pi/180)*point(index,1).MinorAxisLength/2;
    endpt(4,1)= point(index,1).Centroid(1,2)+cos(point(index,1).Orientation(1)*pi/180)*point(index,1).MinorAxisLength/2;
    endpt(5,2)= point(index,1).Centroid(1,1)-sin(point(index,1).Orientation(1)*pi/180)*point(index,1).MinorAxisLength/2;
	endpt(5,1)= point(index,1).Centroid(1,2)-cos(point(index,1).Orientation(1)*pi/180)*point(index,1).MinorAxisLength/2;
    
    
    plot(boundarynp(:,2), boundarynp(:,1), 'b', 'LineWidth', 3);
    hold on
    plot(endpt(:,2),endpt(:,1),'g')
    hold off
    drawnow;
   % figure
    assignin('base','boundarynp',boundarynp)
    minFolders=sqrt((endpt(2,1)-endpt(1,1))^2+(endpt(2,2)-endpt(1,2))^2);
    lngth=length(boundarynp);
    bndryxy=zeros(lngth,2);
    for j=1:lngth
        r=sqrt((boundarynp(j,1)-endpt(1,1))^2+(boundarynp(j,2)-endpt(1,2))^2);
        c=sqrt((boundarynp(j,1)-endpt(2,1))^2+(boundarynp(j,2)-endpt(2,2))^2);
        angl=acos((r^2+minFolders^2-c^2)/(2*r*minFolders));
        
        bndryxy(j,1)=r*cos(angl);
        bndryxy(j,2)=r*sin(angl);
    end
    assignin('base','bndryxy',bndryxy)

    mn1=min(bndryxy(1:round(lngth/2),2));
    [min1x,min1y]=find(bndryxy(1:round(lngth/2),2)==mn1);
    if min1x==1
        mn1=min(bndryxy(1:round(lngth*3/5),2));
        [min1x,min1y]=find(bndryxy(1:round(lngth*3/5),2)==mn1);
    end;
    lngth2=lngth;
    if min1x<0.15*lngth
        lngth2=0.9*lngth;
    end
    mn2=min(bndryxy(min1x+round(lngth/4):round(lngth2),2));
    [min2x,min2y]=find(bndryxy(min1x+round(lngth/4):round(lngth2),2)==mn2);
    min2x=min2x+min1x+round(lngth/4)-1;
    for j=1:min1x
        hf1(min1x-j+1,1)=bndryxy(j,1)*ptonm;
        hf1(min1x-j+1,2)=bndryxy(j,2)*ptonm;
    end
    for j=1:min2x-min1x
        hf2(j,1)=bndryxy(min1x+j,1)*ptonm;
        hf2(j,2)=bndryxy(min1x+j,2)*ptonm;
    end
    assignin('base','hf2',hf2)
    for j=1:lngth-min2x
        hf1(min1x+j,2)=bndryxy(lngth+1-j,2)*ptonm;
    end
    assignin('base','hf1',hf1)

    
   % plot(hf1(:,1), hf1(:,2), 'b', 'LineWidth', 3)
   % hold on
   % plot(hf2(:,1), hf2(:,2), 'g', 'LineWidth', 3)
   % hold off
   % drawnow;
    Vol=0;
    lngthhf1=length(hf1);
    for j=2:lngthhf1
        slc=(hf1(j,1)-hf1(j-1,1))*pi*hf1(j,2)^2;
        Vol=Vol+slc;
    end
    
    %Put First Volume in
    Total(i,1)=abs(Vol);
    %Put First Major Axis in
    mjraxis(1,1)=0;
    mjraxis(1,2)=0;
    for j=1:1+round(lngthhf1/100)
        mjraxis(1,1)=mjraxis(1,1)+hf1(j,1);
        mjraxis(1,2)=mjraxis(1,2)+hf1(lngthhf1-j+1);
    end
    Total(i,3)=abs(mjraxis(1,1)/j-mjraxis(1,2)/j);
    
    rnd=round(lngthhf1/2);
    mnraxis(1,1)=hf1(rnd,2);
    
    for j=1:1+round(lngthhf1/100)
        mnraxis(1,1)=mnraxis(1,1)+hf1(j+rnd,2)+hf1(-j+rnd,2);
    end
    
    Total(i,5)=2*mnraxis(1,1)/(2*j+1);
    
    lngthhf2=length(hf2);
    Vol=0;
    for j=2:lngthhf2
        slc=(hf2(j,1)-hf2(j-1,1))*pi*hf2(j,2)^2;
        Vol=Vol+slc;
    end
    
    %Put Seccond Volume in
    Total(i,2)=abs(Vol);
    %Put Second Major Axis in
    mjraxis(1,1)=0;
    mjraxis(1,2)=0;
    for j=1:1+round(lngthhf2/100)
        mjraxis(1,1)=mjraxis(1,1)+hf2(j,1);
        mjraxis(1,2)=mjraxis(1,2)+hf2(lngthhf2-j+1);
    end
    Total(i,4)=abs(mjraxis(1,1)/j-mjraxis(1,2)/j);
    
    rnd=round(lngthhf2/2);
    mnraxis(1,1)=hf2(rnd,2);
    
    for j=1:1+round(lngthhf2/100)
        mnraxis(1,1)=mnraxis(1,1)+hf2(j+rnd,2)+hf2(-j+rnd,2);
    end
    
    Total(i,6)=2*mnraxis(1,1)/(2*j+1);
    
    Surface=1;
    
    if Surface >0.5
        lnseg=0;
        rd=0;
        SArea=0;
        for j=2:lngthhf1
            lnseg=sqrt((hf1(j,1)-hf1(j-1,1))^2+(hf1(j,2)-hf1(j-1,2))^2);
            rd=(hf1(j,2)+hf1(j-1,2))/2;
            SArea=SArea+lnseg*rd*2*pi;
        end
        Total(i,7)=SArea;
        
        lnseg=0;
        rd=0;
        SArea=0;
        for j=2:lngthhf2
            lnseg=sqrt((hf2(j,1)-hf2(j-1,1))^2+(hf2(j,2)-hf2(j-1,2))^2);
            rd=(hf2(j,2)+hf2(j-1,2))/2;
            SArea=SArea+lnseg*rd*2*pi;
        end
        Total(i,8)=SArea;
    end
 
   i
   
end
   assignin('base','Outlines',Outlines)
   assignin('base','Total',Total) 
    
    % Clear up workspace

clear allSubFolders acCount bb binim cellsrdimensions centrd chr counts
clear cropFoldName filename fps fr frame frame8 frameb frameFolder frameFull
clear frameName frames frs frvec fullpathname g hour00 i ii imagecropped
clear imagefiletiff imagegrysc index listOfFolderNames listOfSecFolders
clear lngthSecfolderNames maxint numberOfMinFolders numberOfSecFolders pathname
clear pic picArray picCell picture remain ROIdimensions se singleSubFolder sumPic TF 
clear TotHeadInfoStrcell up down left right xdim ydim wheretosavename