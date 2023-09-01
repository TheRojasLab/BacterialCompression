clear, close all

basename='14_6_21';
channels=['C:\Users\GuyMa\Documents\Rojas Lab\MatLAB Workspaces\' basename '\'];
recrunch=0;

if recrunch ==1
    cd(['C:\Users\GuyMa\Documents\Rojas Lab\MatLAB WorkSpaces\' basename]);
    load([basename '_intensity.mat']);
   
else
cd(channels)
      directory=dir('*.mat');
T=length(directory);

filenames=cell(T,1);
for i=1:length(directory) % This extracts the file names
  
  [pathstr, name, ext] = fileparts(directory(i).name);
  filenames{i}=name;
end

files = string(filenames) + '.mat'; % adds the extenstion
  %%%%

l=length(files);
intensities_mean=[];
for i = 1:l-1 %loading mean intensities and 
Ls(i)=load(files(i), 'analysisregion');
locations(i)=Ls(i).analysisregion;
    
cellholder=load(files(i), 'intensities_mean');
intensities_mean{i}=struct2cell(cellholder);
intensities_mean{i}=intensities_mean{i}{1};

rateholder=load(files(i), 'meandfluor1');
mratefluor{i}=struct2cell(rateholder);
mratefluor{i}=mratefluor{i}{1}; % average rate increase for every cell
end
    

z=load([files(end)], 'z');
z=z.z;
averages=load([files(end)], 'averages');
averages=averages.averages;
[q r]=size(intensities_mean);
r=1:r;

end


for i = 1:length(intensities_mean)
    i
dfluo{i}=diff(intensities_mean{i}');
dfluo{i}=dfluo{i}';
meanfluo=nanmean(dfluo{i}');
totalSD(i)=nanstd(meanfluo);
totalmean(i)=nanmean(meanfluo); %meanfluo is rate of every cell, total mean is rates of all

Hlimit=totalmean(i)+(2*totalSD);
Llimit=totalmean(i)-(2*totalSD);
Above_h=meanfluo>Hlimit(i);
Below_l=meanfluo<Llimit(i);
Invalid=Above_h + Below_l;
meanfluo=meanfluo(~Invalid);
cellmeans{i}=meanfluo;


% This is now filtering the traces by the rates excluded

Valid=(~Invalid)
    traces{i} = intensities_mean{i}(Valid,:);
    
    
    % Now I'm going to filter by time averaged rate to get rid of extreme
    % cells
    timemean=nanmean(dfluo{i}); % This is the mean rate at each time point
    timeSD=nanstd(dfluo{i});
    Hlimit=timemean(i)+(2*timeSD);
Llimit=timemean(i)-(2*timeSD);
Ncells = length(dfluo{i});
Hlims=ones(Ncells,1)*Hlimit;
Llims=ones(Ncells,1)*Llimit;
InvalidH=dfluo{i}>Hlims;
InvalidL=dfluo{i}<Llims;
Invalid=InvalidH + InvalidL;
Invalid_ind=find(Invalid);
filtrate{i}=dfluo{i};
filtrate{i}(Invalid_ind)=NaN;

Valid2=(~Invalid);
Valid2=[ones(1,Ncells) Valid];
    traces{i} = intensities_mean{i}(Valid,:);



% now filter again
meanfluo=nanmean(filtrate{i}');
totalSD(i)=nanstd(meanfluo);
totalmean(i)=nanmean(meanfluo);

Hlimit=totalmean(i)+(2*totalSD);
Llimit=totalmean(i)-(2*totalSD);
Above_h=meanfluo>Hlimit(i);
Below_l=meanfluo<Llimit(i);
Invalid=Above_h + Below_l;
meanfluo=meanfluo(~Invalid);
cellmeans{i}=meanfluo;


end


[Dimy Dimx] =  size(intensities_mean{1});
figure
hold on
for i=1:length(dfluo)
plot(dfluo{i})
end
hold off

figure
hold on
for i=1:length(dfluo)
plot(filtrate{i})
end
hold off


for i = 1:length(traces)
   tracemeans(i,:)=nanmean(traces{i});
end
time=0:length(tracemeans)-1;


figure
hold on
for i=1:l-1
   plot(time,tracemeans(i,:))
end
legendStrings = string(locations);
legend(legendStrings,'location','NorthWest')
xlabel('Time (mins)')
ylabel('Average Fluorescence')
hold off


figure
hold on
for i=1:l-1
   plot(time,intensities_mean{i})
end
xlabel('Time (mins)')
ylabel('Average Fluorescence')
hold off



figure
hold on
for i =1:l-1
    hist(cellmeans{i})
end
xlabel('Rate of increase')
ylabel('Count')
hold off
hold off


h_locs=interp1(z,averages,locations); %averages is fluor values of dye, z is x dimension
% rates v height
figure
plot(h_locs,totalmean) %totalmean is rate of every cell averaged - graph looks weird
xlabel('Interpolated height (\mum)')
ylabel('Mean Rate at Location')


alldiffs=diff(tracemeans'); %tracemeans is the mean change at every time point per experiment
meanymeans=nanmean(alldiffs); %mean change for all time points per experiment

figure
plot(h_locs,meanymeans)  
hold on
xlabel('Interpolated height (\mum)')
ylabel('Mean Rate at Location')
hold off

for i=1:l-1
allSDs(i)=nanstd(cellmeans{i});
N=length(isnan(cellmeans{i}));
stderror(i) = (allSDs(i)/sqrt(N));
end


% now I want to add the mean rate of change from 0 to 70 mins
basename='9_2_22';
channels=['C:\Users\GuyMa\Documents\Rojas Lab\MatLAB Workspaces\FLUO\' basename '\'];

cd(channels);
      directory=dir('*.mat');
T=length(directory);

filenames=cell(T,1);
for i=1:length(directory) % This extracts the file names
  
  [pathstr, name, ext] = fileparts(directory(i).name);
  filenames{i}=name;
end

files = string(filenames) + '.mat'; % adds the extenstion
  %%%%
l=length(files);
for i = 1 %loading mean intensities and 
Rs(i)=load(files(i), 'rates');
Er(i)=load(files(i), 'rateerr');
Ms(i)=load(files(i), 'means');
aSDs(i)=load(files(i), 'allctrlSD');
MOE(i)=load(files(i), 'moe');
moe=MOE.moe;
rates=Rs(i).rates;
Ctrlmeans=Ms(i).means;
allSD=aSDs(i).allctrlSD;
rateerr=Er(i).rateerr;
end

dirname= ['C:\Users\Guyma\Documents\Rojas Lab\MatLAB WorkSpaces\FLUO\WT and Low'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 cd(dirname);
directory=dir('*.mat');
T=length(directory);
filenames=cell(T,1);
for i=1:length(directory) % This extracts the file names
  [pathstr, name, ext] = fileparts(directory(i).name);
  filenames{i}=name;
end
files = string(filenames) + '.mat'; % adds the extenstion
  for tt=2;
NR=load(files(tt), 'meanintensities'); %loads the  data
NRall=load(files(tt), 'holdcells');
noreportfluos{1}=NRall.holdcells;
noreport=NR.meanintensities;
  end

%for the dashed line from final data point to uncompressed

fillinex=[h_locs(end) 1.1];
filliney=[meanymeans(end) rates];


figure
hold on
   errorbar(h_locs(2:end),meanymeans(2:end),stderror(2:end),'Color','b')
xlabel('Device Height (\mum)')
ylabel('Rate of Induction (M^-^1)')
%yline(rates,'-.b','Rate when cells are not trapped')
%errorbar(1.1,rates,)
errorbar(1.1,rates,rateerr,'Color','b')
    plot(fillinex,filliney,'--','Color','b')
xlim([0.5 1.2])
xticks([0.5 0.6 0.7 0.8 0.9 1 1.1])
xticklabels({'0.5','0.6','0.7','0.8','0.9','1','<1'})

plot(h_locs(2),meanymeans(2),'r*','MarkerSize',10)
plot(h_locs(3),meanymeans(3),'b*','MarkerSize',10)
plot(h_locs(4),meanymeans(4),'g*','MarkerSize',10)
plot(h_locs(5),meanymeans(5),'k*','MarkerSize',10)
plot(h_locs(6),meanymeans(6),'m*','MarkerSize',10)
fig2pretty
%set(gca,'fontsize', 20)
text(0.85-0.15,rates+10,'Rate of Untrapped Cells','FontSize',15)
hold off

%saveloc=['C:\Users\Guyma\Documents\Rojas Lab\MatLAB WorkSpaces\FLUO'];
%cd(saveloc);
%save([basename '_intensity'])

SEM=[];
means=[];
std_dev=[];
CI95=[];
altdev=[];
%nanmean(fluo)
%nanstd(fluo)
%nanstd(fluo/sqrt(N)) N is number of non NaNs 
%. before divide means element by element


for i=1:length(traces) %getting all fluor values of all cells
means{i}=nanmean(traces{i});
std_dev{i}=nanstd(traces{i});

N=sum(~isnan(traces{i}(:,1)));
std_err{i}=nanstd(traces{i}./sqrt(N));

for a=1:size(traces{i},2) %For each time point generate CIs
SEM{i} = std_err{i}; %/sqrt(nnz(~isnan(c_counts1{i}(:,a))));
nonzeros=nnz(~isnan(traces{i}(:,end))); %Find out how many cells per time point for DoF
CI95= tinv([0.025 0.975], nonzeros-1);    % Calculate 95% Probability Intervals Of t-Distribution
yCI95 = bsxfun(@times, SEM{i}(end), CI95); 
CI95n{i,a}=yCI95(1);
CI95p{i,a}=yCI95(2);

end
end

%need mean per time point
for i=1:length(traces)
   meanpertp{i}=nanmean(traces{i});
end


%Add CIs to data
for i=1:length(meanpertp) %not cellmeans
B=meanpertp{i};
C=[CI95n{i,:}];
CI1{i}=B+C;


B=meanpertp{i};
C=[CI95p{i,:}];
CI2{i}=B+C;
line_color = ['r' 'b' 'g' 'k' 'm' 'y' 'c'];
end

for i=1:length(meanpertp)
x2{i} = 0:length(meanpertp{i})-1;
end

   figure % Fig 9
%ciplot(CI2{1},CI1{1},x2{1},'r')
%hold on
ciplot(CI1{2},CI2{2},x2{2},'r')
hold on
ciplot(CI1{3},CI2{3},x2{3},'b')
hold on
ciplot(CI1{4},CI2{4},x2{4},'g')
hold on
ciplot(CI1{5},CI2{5},x2{5},'k')
hold on   
ciplot(CI1{6},CI2{6},x2{6},'m')
hold on 

 for i=2:length(meanpertp)
plot(x2{i},meanpertp{i},'Color',line_color(i-1))
 end
 
 
 
Locstrings=round(h_locs*100)/100;
Locstrings=Locstrings(2:end);
Locstrings=num2str(Locstrings);
nopressmeans=zeros(1,length(time));
nopressmeans(1)=Ctrlmeans(1);
for i =1:length(time)-1
    nopressmeans(i+1)=nopressmeans(i)+rates;
end
 
%now make the CIs for nopressmeans
NPCI1=nopressmeans+moe;
NPCI2=nopressmeans-moe;

% Now I want the CIs for the controls with no reporter

 %getting all fluor values of all time points not cells individually
    noreportmeans=nanmean(noreportfluos{1});
%means{i}=nanmean(traces{i});
noreportstd_dev=nanstd(noreportfluos{1});

N=sum(~isnan(noreportfluos{1}(:,1)));
noreportstd_err=nanstd(noreportfluos{1}./sqrt(N));

for a=1:size(noreportfluos{1},2) %For each time point generate CIs
noreportSEM = noreportstd_err; %/sqrt(nnz(~isnan(c_counts1{i}(:,a))));
noreportnonzeros=nnz(~isnan(noreportfluos{1}(:,end))); %Find out how many cells per time point for DoF
CI95= tinv([0.025 0.975], noreportnonzeros-1);    % Calculate 95% Probability Intervals Of t-Distribution
yCI95 = bsxfun(@times, noreportSEM(end), CI95); 
ctrlCI95n=yCI95(1);
ctrlCI95p=yCI95(2);
end
%end

B=noreportmeans;
C=ctrlCI95n;
ctrlCI1=B+C;
B=noreportmeans;
C=ctrlCI95p;
ctrlCI2=B+C;
noreporttime=0:length(noreport)-1;

% I want to stretch out the no reporter data
interpednoreport=interp1(noreporttime,noreport,[31:60],'linear','extrap');
noreport = [noreport interpednoreport];
noreporttime2=0:length(noreport)-1;

%for a=1:size(noreport) %For each time point generate CIs

cd('C:\Users\GuyMa\Documents\Rojas Lab\MatLAB WorkSpaces\FLUO\RcsF Del');
datastruc=load('RcsAdel');
RcsFdel=datastruc.c_counts{6,1};
RcsFtime=1:length(RcsFdel);
% extrapolate RcsF data
interpedRcsF=interp1(RcsFtime,RcsFdel,[31:60],'linear','extrap');
RcsFdel=[RcsFdel interpedRcsF];

plot(time,nopressmeans,'--','Color','r')
ciplot(NPCI1,NPCI2,time,'r')
ciplot(ctrlCI1,ctrlCI2,noreporttime)
plot(noreporttime2,noreport,'--','Color','b')
plot(time,RcsFdel)
fig2pretty
set(gca,'fontsize', 15)
text(40,nopressmeans(end)+650,'No Compression','FontSize',12)
text(60,noreport(end)+100,'No Reporter','FontSize',12)
text(60,RcsFdel(end)+300,'RcsF deletion','FontSize',12)
ylabel('Fluorescence (A.U.)')
xlabel('Time (Mins)')
legend('Location','northwest');
lgd=legend('0.63\mum','0.71\mum','0.73\mum','0.75\mum','0.98\mum');
title(lgd,'Device Height','FontSize',12)
hold off