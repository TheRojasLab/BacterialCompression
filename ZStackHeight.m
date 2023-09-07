%  This is a script to get the average height profiles across Z stacks by
%  using a fluorescent dye in a CellASIC chip


%INSTRUCTIONS FOR USE:
%Use BacTrack2ROI %and BTfluoROI (Region of Interest) to generate
%tiled data across the whole chip
%
%INPUT:
%basename: name of the .mat directory.
%dirname:the full pathname of the directory where you saved the ROI .mat files
%lscale:microscope calibration in microns per pixels.
%trappos:location where cells begin to get trapped in uM. Tallest ->
%Shortest e.g. 1500 -> 0uM 

%OUTPUT:
%z: x dimension of chamber. Chamber position 
%averages: height of chamber


%Calls on the following m-files:
%norm16bit.m
%fig2pretty.m
%movingaverage


clear, close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%User Input
basename='14_6_21';
dirname={['C:\Users\GuyMa\Documents\Rojas Lab\Exported Movies\' basename '\']};
regname=['C:\Users\GuyMa\Documents\Rojas Lab\Exported Movies\' basename '\'];

% look at where cells begin to get trapped, this is the xposition
trappos=400;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(regname)


recrunch=1;
 if recrunch==1
      load([basename '_z_stack.mat'])
 else
 
directory=dir('*tif'); %directory of all the .tif files
imagename=directory(1).name; %imagename becomes first image
im=imread(imagename); %im becomes vector of pixels

level=graythresh(im);
    bw=im2bw(im,level);
 stats_lab=cell(1,1);

  %De-speckle
        im=medfilt2(im);
        
        %Normalize images
        ppix=0.5;
        im2=norm16bit(im,ppix);

windows=18;
t=length(im2);
rollwindow=t/windows;
rollwindow=round(rollwindow);
mwindows=[];
avg=[];
B=[];
%%%%% #of images across
totalimages=18;
%%%%%
 end


recrunch2=1;
 if recrunch2==1
      load([basename '_z_stack.mat'], 'smooth_im')
 else
 smooth_im=movingaverage(im2,1200);
 end

 
 
 
imshow(smooth_im,[]);

%%lets get the average column values and plot them
averages=[];
averages=mean(smooth_im);
%averages=averages/1000;
x=length(averages);
position=(1:x);
position=position/10;

y=1500/x;
z=(0:y:1500);
z(end)=[];

a=1.5-0.75;
y=a/x;

   onez=length(z)/1500;
   ztrap=trappos*onez; %place in the chip where cells are trapped
   ztrap=round(ztrap); %make it an integer
   % equation for fluorescence is h= 1um/fluor(1um)-f(0) * (f-f(0)) 
   f0= 127;
   f1um=averages(ztrap);
   height=[];
   
   for i=1:length(averages) %Find the height of each spot from the fluorescence
   height(i)=1/(f1um-f0)*(averages(i)-f0);
   end
averages=height;
   
figure
plot(z,averages)
xlabel('Chamber Position (\mum)')
%set(gca,'XTick',[])

hold on
grid on
fig2pretty
xlabel('Chamber Position (\mum)','FontSize', 20)
ylabel('Height (\mum)','FontSize',20)


%save([basename '_z_stack'])

