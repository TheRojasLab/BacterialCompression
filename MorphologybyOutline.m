% This is a script that takes the normalisation from Ecoli count and adds
% it to the colony counts script


clear, close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUT
recrunch=1;
basename='4_9_19_100X';
%dirname=['E:\Google Drive\Data Share\Guy\Computer Sync\' basename '\' basename '_2_a'];
dirname=['C:\Users\GuyMa\Documents\Rojas Lab\Experiments\'  basename '\' basename '_2_a'];
cd(dirname);

if recrunch == 1
 %   load (dirname '\' basename '_CNT');
 load ('4_9_19_100X_CNT_PersisterMorphs.mat');
else 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%onecolony=1.2185e+03; %this was found by dividing the total #of counted pixels
% by the # of counted bacteria in the image - this is for 100X
%onecolony= 50;

%determine number of frames

directory=dir('*tif'); %directory of all the .tif files
T=length(directory); %make a string of length # of images



    
%%%%%%%%%%%%%%%%%%%%%%%
tscale=180; %frame rate in seconds
phagein=80;
persistersleft=115;
%Thresholding
%%%%%%%%%%%%%%%%%%%%%%%
ppix=0.2; %This is what % the of the top and bottom boundaries of the pixel values will be
cell_thresh=0.9; %This applies the threshold for upper limit of pixels to count
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
time=[0:T-1]*tscale/60; %create time vector from tscale for # of mins

imagename=directory(1).name; %imagename becomes first image
im=imread(imagename); %im becomes vector of pixels
[m,n] = size(im); %dimensions of im
%z=m*n*T; %z is total size of vector to be filled (1 image dims x # of images)
z=m*n;
imsize=m*n; % imsize is one image dims
pixel_int=zeros(z,1); %vector of values to be filled in, in one column
     
[n,m]=size(im);
morphs=[]; % this is for the morphology tracking


%Normalising
% normally 1:T
for t=1:10; % This for loop is used to find all the pixel values to generate thresholds
    t % image
    imagename=directory(t).name; %As above, taking image
im=imread(imagename); % reading the image into pixels
    imvec=im(:); %Creates a vector with the above values in 1 column
    pixel_int((t-1)*imsize+1:t*imsize,1)=imvec; %This decides where to place the values in
    % the pixel intensity vector
end;




  %Saturate some defined % of the image
    [im,minint,maxint]=norm16bit(pixel_int,ppix);
    im=norm16bit(pixel_int,ppix);
    %Saturate some % of the pixel_int vector and define the min/max int
    %values
    lowin=minint/2^16; %Convert minint value to a value between 0-1
    highin=maxint/2^16; %Convert maxint value to a value between 0-1

    

   st=strel('disk',1);    
    %This for loop applies the threshold to each image
    for i=1:T; %this for loop is to apply the new thresholds to each image
    i
     imagename=directory(i).name; %start at the first image
   im=imread(imagename) ;  % read the pixels into im


   imalt=imadjust(im,[lowin;highin],[0;1]); 

   %Scales the image to between 0-1 using the defined thresholds
   % by saturating the bottom/top values as 0/1
   

   %figure
   %imhist(imalt)
   %pause
     bw=im2bw(imalt,cell_thresh); %create matrix of 0 and 1 values, by thresholding images
     
    
 CC = bwconncomp(bw);
 numberofpersisters(i)=CC.NumObjects; %gives number of connected things at each time point
 %plot this and see when best time point for calculating no of persisters
     
 stats=regionprops(CC,'Centroid'); % measures properties of image regions
 xycoords=[stats.Centroid]; %gives x and y coords
 xcoords=xycoords(1:2:end); %x coordinates of all connected components (skips y coords)
 ycoords=xycoords(2:2:end);

 
f=0;
count=0;


    if i==phagein
        
        imshow(imalt); %shows the image I will track the morphology of
hold on;

load ('4_9_19_100X_CNT_NPs.mat','morphs');
                plot(morphs(:,1),morphs(:,2),'r');
                hold on
while f~=1



%imshow(imalt); %shows the image I will track the morphology of
%hold on;

load ('NPmorphs.mat','morphs');
                plot(morphs(1,:),morphs(2,:),'r');
                hold on

G=0; % this is a loop to continue the cell counts
while G==0  ;
    
    disp('Track cells?');
    disp ('1 for yes');
    disp ('2 for no');
    disp ('3 for zoom');
    P = input(' ');
    switch P;
        case 1;
            xmorphs_vec=[];
            ymorphs_vec=[];
            k=0;
            while k~=1;
               % zoom on               
                k=waitforbuttonpress;
                [xmorphs,ymorphs]=ginput(1); %Pick point/s of cell outline
                xmorphs_vec=[xmorphs_vec xmorphs];
                ymorphs_vec=[ymorphs_vec ymorphs];
                plot(xmorphs_vec,ymorphs_vec,'r');
                G=0  ;               
            end
             
             
             morphs{count+1,1}=xmorphs_vec;
             morphs{count+1,2}=ymorphs_vec;
             count=count+1;
             
        case 2;
            G=G+1;
            f=1
        case 3
            zoom on;
            G=0;
            pause;
            zoom off;
    end
    
end



end
 
end


if i==phagein;
    firstvalue=xcoords;
end;

if i==persistersleft;
   secondvalue=xcoords; 
end;
 

 % want histogram of x positions pre and post phage, put them on same graph
 % to compare them. Bins need to be same. Define centre and width of
 % histograms. Zero to max histo, not all on one side
 
  im_dil=imdilate(bw,st); %dilate the image using 3x3 structure to create more 1s
  dil_ends=im_dil-bw; %Matrix - matrix to leave only dilated 1s
  
  overlaid=imoverlay(imalt,dil_ends,'red'); %Overlay red lines using dilution
  
%if vis==1;
    imshow(overlaid);
    
i

%location=xycoords(i)
%figure
%plot(i)(xycoords
    end
    
bins=8;

[k,l]=size(im);
numbins=round(l/bins);
    xbins=(1:numbins:l);
    
    xbins=(1:bins+1)*l/(bins+1);
xbins=(xbins(2:end)+xbins(1:end-1))/2;

end % this is the end of recrunch

% compare ccs at phage in vs persisters
 %figure
 %histogram(firstvalue,xbins, 'FaceColor','b')
 %hold on
 %histogram(secondvalue,xbins, 'Facecolor','k')
 %hold off
    
 figure;
 imshow(imalt);
 hold on;
 %axis on;
 
 idx=[];
 idy=[];
 
 for ii = 1:length(morphs)
     
idx= morphs{ii,1};
meanx=mean(idx);
normalx=(idx-meanx);

idy= morphs{ii,2};
meany=mean(idy);
normaly=(idy-meany)+(ii*50);

%px and py would be the x and y locations of the cell boundaries. 
%take a look at the documentation for csaps.


%dS=sqrt(diff(idx).^2+diff(idy).^2);
%S=[0 cumsum(dS)'];
%px=csaps(S,idx,0.05,S);
%py=csaps(S,idy,0.05,S);


dS=sqrt(diff(normalx).^2+diff(normaly).^2);
S=[0 cumsum(dS)'];
px=csaps(S,normalx,0.05,S);
py=csaps(S,normaly,0.05,S);

 %plot(px,py);
 end
% xlabel('uM')
axis equal
 
hold off



