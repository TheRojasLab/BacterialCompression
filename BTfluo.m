%BTfluo.m
%Rico Rojas, updated 1/21/19
%Calculates the average cytoplasmic fluorescence intensity from cell
%tracked with BacTrack.m.  

clear, close all

%INSTRUCTIONS FOR USE:
%Remove frames with poor contrast and save fluorescent image stacks
%directories by themselves. 

%INPUT
%basename: name to save the results to.
%channels: list of directories containing fluorescent image stacks to quantify.

%OUTPUT:
%icell: Cell array with length equal to the number of fluorescent
        %channels.  Each entry is a matrix (ncellxT) with the fluorescent intensities of each
        %cell, where rows are the cells and columns are time points.
%icell_av:  Cell array with length equal to the number of fluorescent
        %channels.  Each entry is a vector containing the population-
        %average of the single-cell fluorescent intensities.

% Remember to move BTfluo and lab to _2_a

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%USER INPUT
basename='18_4_23_1';
channels={['C:\Users\GuyMa\Documents\Rojas Lab\Exported Movies\' basename '\' basename '_2_a']};
recrunch=0;
fluor=0; %0=Skip the filtering step if there isn't much fluorescence in the movie
%nobrightcell=0 %if there are super high bright cells messing with error bars
viz=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(channels{1});

if recrunch==0;

curdir=cd;

for i=1:length(channels)
    cd(channels{i}); 
    fluo_directory{i}=dir('*.tif');
end

load([basename '_BT'])
load([basename '_BTlab'])

intensities=cell(length(channels),1);

T = length(fluo_directory);
for i=1:length(channels)
    cd(channels{i}); 
    intensities_temp=zeros(size(lcell));
    for t=1:T
       % t
        imagename=fluo_directory{i}(t).name;
        im=imread(imagename);   
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
 verify = abs(dim1-dim2);
 if verify < 1;
     % This is only at one place, should really be the whole row
     B{k,t}=[]; 
     intensities_all{k,t}=[];  % Now make pixels empty where B is
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
 forplot{k,t}=[];
 end
     end
 end

 % I want to try and filter by total change again
 for i=1:size(intensities_all,1) % This filters out the posts by value of intensities_all
    tbcleaned=intensities_all(i,:);
     idk=tbcleaned(~cellfun('isempty',tbcleaned));
     verify=isempty(idk);
     if verify == 0;
     idk=idk';
         allvals=cell2mat(idk);
        uintmin = min(allvals) ;
        uintmin=double(uintmin);
           uintmax = max(allvals) ;
        uintmax=double(uintmax);
        fluordiffs(i) = abs(uintmin-uintmax);
     end
end

         thresholdscore = (mean(fluordiffs)*0.5);%;
        for i=1:size(intensities_all,1)
     if fluordiffs(i) < thresholdscore;
        forplot(i,:)={[]}; 
        intensities_all(i,:)={[]};
      end
        end
 
 %recalc intensities_temp
  for i=1:length(intensities_all)
     for t=1:size(intensities_all,2)
 intensities_temp(i,t)=nanmean(intensities_all{i,t});

     end
 end
 
 if T ~= 1
%This is a loop to generate the diff of the fluorescence
dfluor=[];
for i = 1:length(intensities_temp)
    i;
fluorchange=diff(intensities_temp(i,:));
%************** this might have to be changed
dfluor(i,:)=(fluorchange/time(2)); %dfluor per minute for all cells at all time points

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
       
  % check=isempty(intensities_all{i,1})
  % elseif check==1
  %      meandfluor1(i) = NaN;
       
end
end

% histogram of rates of all cells over every time point
hist(meandfluor1)

% Now clear all rows in intensities_all that aren't in meandfluor1
%because meandfluor is the mean of the whole thing, if it's empty
%I should clear whole row 
 for t= 1:size(meandfluor1,2)
     
     reference=meandfluor1(t);
 verify=isnan(reference);
 
 if verify==1;
     [X,Y]=size(intensities_all(t,:));
     for i=1:Y 
         intensities_all{t,i}=[];
         forplot{t,i}=[];
 end
 end
 end
 
 
 % also get intensities_mean
 for i=1:length(intensities_all)
     for t=1:Y
 intensities_mean(i,t)=nanmean(intensities_all{i,t});
 
     end
 end
 figure
 plot(time,intensities_mean)
 hold on fig2pretty
  xlabel('Time (s)')
ylabel('Intensity (A.U.)')
     
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
 end
  
 elseif recrunch==1
    load ([basename '_BTfluo'])
end
 
         extracted=forplot(:,1);
 bounplot=extracted(~cellfun('isempty',extracted));
 figure
     imshow(im,[])
    hold on
 for i=1:size(bounplot,1)
    plot(bounplot{i,1}(:,1),bounplot{i,1}(:,2),'-r')   
 end
    hold off  
    


save([basename '_BTfluo' ])
