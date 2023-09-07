%BTfluo.m
%Guy Mason, updated 7/9/23
%Calculates the average cytoplasmic fluorescence intensity from cell
%tracked with BacTrack2ROI.m.  

clear, close all

%INSTRUCTIONS FOR USE:
%Remove frames with poor contrast and save fluorescent image stacks
%directories by themselves. 

%INPUT
%basename: name to save the results to.
%channels: list of directories containing fluorescent image stacks to quantify.
%analysisdist: window of CellASIC chip to analyse

%OUTPUT:
%intensities_temp:  Matrix where each entry is an average of the fluorescent intensities of each
        %cell, where rows are the cell and columns are time points.
%intensities_all:  Cell array where each entry is a vector containing the
%values of the single-cell fluorescent intensities.

% Remember to move BTfluo and lab to _2_a

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%USER INPUT
basename='14_6_21';
channels={['C:\Users\GuyMa\Documents\Rojas Lab\Exported Movies\' basename '\' basename '_2_a']};
recrunch=0;
fluor=0; %0=Skip the filtering step if there isn't much fluorescence in the movie
analysisdist=200;
viz=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(channels{1});

if recrunch==0;

curdir=cd;

for i=1:length(channels)
    cd(channels{i}); 
    fluo_directory{i}=dir('*.tif');
end

analysisdist=num2str(analysisdist);
load([basename '_' analysisdist '_BT'])
load([basename '_' analysisdist '_BTlab'])

intensities=cell(length(channels),1);

for i=1:length(channels)
    cd(channels{i}); 
    intensities_temp=zeros(size(lcell));
    for t=1:T
       % t
        imagename=fluo_directory{i}(t).name;
        im=imread(imagename);
          x=size(im,2);
    windowsize=round(x/15); %1370
    analysisdistance=analysisregion/100; %2
    viewpoint=analysisdistance*windowsize; %274
    viewpoints=[viewpoint-windowsize viewpoint+windowsize];
    
   if analysisdistance == 1
                    viewpoint=windowsize; %274
    viewpoints=[1 windowsize];
       im=im(:,1:windowsize);
       
       
    elseif analysisdistance == 14
     %   disp('hi')
            viewpoint=14*windowsize; %274
    viewpoints=[viewpoint-windowsize size(im,2)];
            im=im(:,viewpoints(1):end);
            
            
    else
    im=im(:,viewpoints(1):viewpoints(2));

   end
    
        for j=1:ncells
            intensities_all{j,t}=im(pixels{j,t});
            intensities_temp(j,t)=mean(im(pixels{j,t})); %This point gets the pix values
        end
    end
    intensities_temp(intensities_temp==0)=NaN;
    icell{i}=intensities_temp;
    
    Dims=size(intensities_temp);
    multi=Dims(1)*Dims(2);
    cellids=(1:multi);
    cellids=reshape(cellids,Dims(1),Dims(2));
    idcell{i}=cellids;
    
end

holdcells=icell{1};
[q r]=size(icell{1});

for i = 1:Dims(1)
    for d = 1:Dims(2) 
        verify = holdcells(i,d);
        query = isnan(verify);
        if query ==1;
   cellids(i,d)=NaN;
   pixels{i,d}=NaN;
        end
    
    end
end
  idcell{1}=cellids;

  
  
  mask=[];
  mask2=[];
  %Now lets check that pixels and hold cells are the same
  for i = 1:Dims(1)
    for d = 1:Dims(2) 
        
        verify = holdcells(i,d);
        query = isnan(verify);
        if query ==1;
   mask(i,d)=1;
        else
            mask(i,d)=0;
        end
    verify2 = pixels{i,d};
    query2 = isnan(verify2);
    if query2 ==1;
    mask2(i,d)=1;
    else
        mask(i,d)=0;
    end
end
  end
  
  isequal(mask,mask2)
  

%Lets remove background
background=min(im,[],'all');

for aa=1:length(holdcells)
    for bb=1:size(holdcells,2)
        validate=isnan(holdcells(aa,bb));
 if validate ==0;
     holdcells(aa,bb)=holdcells(aa,bb)-background;
 end
    end
end

elseif recrunch==1
    load ([basename '_BTfluo'])
end

icell_av=[];
for i=1:size(holdcells,2)
    icell_av(i)=nanmean(holdcells(:,i));
end
%end


forplot=B;
 for k = 1:size(B,1)
     for t= 1:size(B,2)
         
 % Get rid of the cells in intensities_all that aren't in forplot
 reference=forplot{k,t};
 verify=isempty(reference);
 if verify==1;
 intensities_all{k,t}=[];
 
 else
  %now I'm going to filter by aspect ratio
        dim1=range(B{k,t}(:,1));
         dim2=range(B{k,t}(:,2));
 absvals{k,t} = abs(dim1-dim2);
 if  absvals{k,t} < 5;
     % This is only at one place, should really be the whole row
    B{k,t}=[]; 
   %  intensities_all{k,t}=[];  % Now make pixels empty where B is
 end
 end

 end
 end

 
 for k = 1:size(B,1)
     for t= 1:size(B,2)
         
 % Get rid of the cells in pixels that aren't in forplot
 reference=intensities_all{k,t};
 verify=isempty(reference);
 if verify==1;
 B{k,t}=[];     
 forplot{k,t}=[];
 end
     end
 end
        
        
 %recalc intensities_temp
  for i=1:length(intensities_all)
     for t=1:size(intensities_all,2)
 intensities_temp(i,t)=nanmean(intensities_all{i,t});

     end
  end
 
 
%This is a loop to generate the diff of the fluorescence
dfluor=[];

for i = 1:length(intensities_temp)
    i;
fluorchange(i,:)=diff(intensities_temp(i,:));
%************** this might have to be changed
dfluor(i,:)=(fluorchange(i,:)/time(2)); %dfluor per minute for all cells at all time points

%now the averages need it only of the same time points
check=isnan(nanmean(dfluor(i,:)));
if check==0
finalpoint= size(dfluor,2);
tbmeaned{i}=dfluor(i,1:finalpoint-1); % means of all cells at all time points
meandfluor(i)=nanmean(tbmeaned{i}); %These are the average rate increases
end
end


%Figure for plotting rate vs time
figure
hold on
plot(time(1:end-1),dfluor)%(:,i)(1:21))
xlabel('Time (Mins)')
ylabel('Rate of Induction (Mins^-^1)')
hold off
 

%Now get mean rates and filter by extremes
S=nanstd(meandfluor);
M=nanmean(meandfluor);
SDP=M+S+S;
SDN=M-S-S;
meandfluor1=meandfluor;
for i=1:length(meandfluor)
   if meandfluor1(i) < SDN ;
    meandfluor1(i)=NaN;
   elseif meandfluor1(i) > SDP;
       meandfluor1(i) = NaN;

       
end
end

                  figure
            imshow(im,[])
            hold on
            for t=1:size(B,2)
            for k=1:length(B)
                verify = isnan(B{k,t});
                if verify == 0
                plot(B{k,t}(:,1),B{k,t}(:,2),'-r')
                end 
            end
            end
 


% histogram of rates of all cells over every time point
figure
hist(meandfluor1)

     
% Now to see what's being tracked

  for k = 1:size(B,1)
     for t= 1:size(B,2)
         
 % Get rid of the cells in intensities_all that aren't in forplot
 reference=intensities_all{k,t};
 verify=isempty(reference);
 if verify==1;
 forplot{k,t}=[];
 end
     end
  end

          figure
            imshow(im,[])
            hold on
            for k=1:length(forplot)
                verify = isnan(forplot{k,t});
                if verify == 0
                plot(forplot{k,t}(:,1),forplot{k,t}(:,2),'-r')
                end 
            end
    

basename=num2str(basename);
save([basename '_' analysisdist '_BTfluo' ])