clear, close all
% 10% w/v PDMS script
%Tracks bacterial fluoresence in the torture chamber.
%Customized for E. coli.


%INSTRUCTIONS FOR USE:
% Have directory of BTfluo files individually partitioned 
%
%INPUT:
%channels: name of the file containing BTfluo .mat-s
%
%
%OUTPUT:
%times: cell array containg number of time points.
%normexp: vector of normalized mean cell fluorescence change.
%avxchange: vector of mean cell fluorescence change.
%avx: vector of mean cell fluorescence at each time point

%Calls on the following m-files:
%shadedErrorBar.m
%fig2pretty.m



channels=['C:\Users\GuyMa\Documents\Rojas Lab\MatLAB Workspaces\FLUO/Smush\10%\0'];
recrunch=0;
cd(channels)
      directory=dir('*.mat');
T=length(directory);
filenames=cell(T,1);
for i=1:length(directory) % This extracts the file names
  [pathstr, name, ext] = fileparts(directory(i).name);
  filenames{i}=name;
end
files = string(filenames) + '.mat'; % adds the extenstion
l=length(files);
for i = 1:l %loading data
intenses{i}=load(files(i), 'icell_av');
holder=intenses{i}.icell_av;
avmaxintensities{i}=holder{1};
 
icells{i}=load(files(i), 'icell');
allcellis=icells{i}.icell;
allintensities{i}=allcellis{1};
SDs{i}=nanstd(allintensities{i});

Ts(i)=load(files(i), 'time');
times{i}=Ts(i).time;
end

channels=['C:\Users\GuyMa\Documents\Rojas Lab\MatLAB Workspaces\FLUO/Smush\10%\1.25'];
recrunch=0;
cd(channels)
      directory=dir('*.mat');
T=length(directory);
filenames=cell(T,1);

for i=1:length(directory) % This extracts the file names
  [pathstr, name, ext] = fileparts(directory(i).name);
  filenames{i}=name;
end
files = string(filenames) + '.mat'; % adds the extenstion
l=length(files);
for i = 1:l %loading data
intenses{i}=load(files(i), 'icell_av');
holder=intenses{i}.icell_av;
avmaxintensities=[avmaxintensities holder{1}];
 
icells{i}=load(files(i), 'icell');
allcellis=icells{i}.icell;
allintensities=[allintensities allcellis{1}];
SDs=[SDs nanstd(allintensities{i})];

Ts(i)=load(files(i), 'time');
times=[times Ts(i).time];
end

channels=['C:\Users\GuyMa\Documents\Rojas Lab\MatLAB Workspaces\FLUO/Smush\10%\1.5'];
recrunch=0;
cd(channels)
      directory=dir('*.mat');
T=length(directory);
filenames=cell(T,1);

for i=1:length(directory) % This extracts the file names
  [pathstr, name, ext] = fileparts(directory(i).name);
  filenames{i}=name;
end
files = string(filenames) + '.mat'; % adds the extenstion
l=length(files);
for i = 1:l %loading data
intenses{i}=load(files(i), 'icell_av');
holder=intenses{i}.icell_av;
avmaxintensities=[avmaxintensities holder{1}];
 
icells{i}=load(files(i), 'icell');
allcellis=icells{i}.icell;
allintensities=[allintensities allcellis{1}];
SDs=[SDs nanstd(allintensities{i})];

Ts(i)=load(files(i), 'time');
times=[times Ts(i).time];
end

channels=['C:\Users\GuyMa\Documents\Rojas Lab\MatLAB Workspaces\FLUO/Smush\10%\2'];
recrunch=0;
cd(channels)
      directory=dir('*.mat');
T=length(directory);
filenames=cell(T,1);

for i=1:length(directory) % This extracts the file names
  [pathstr, name, ext] = fileparts(directory(i).name);
  filenames{i}=name;
end
files = string(filenames) + '.mat'; % adds the extenstion
l=length(files);
for i = 1:l %loading data
intenses{i}=load(files(i), 'icell_av');
holder=intenses{i}.icell_av;
avmaxintensities=[avmaxintensities holder{1}];
 
icells{i}=load(files(i), 'icell');
allcellis=icells{i}.icell;
allintensities=[allintensities allcellis{1}];
SDs=[SDs nanstd(allintensities{i})];

Ts(i)=load(files(i), 'time');
times=[times Ts(i).time];
end

channels=['C:\Users\GuyMa\Documents\Rojas Lab\MatLAB Workspaces\FLUO/Smush\10%\2.25'];
recrunch=0;
cd(channels)
      directory=dir('*.mat');
T=length(directory);
filenames=cell(T,1);

for i=1:length(directory) % This extracts the file names
  [pathstr, name, ext] = fileparts(directory(i).name);
  filenames{i}=name;
end
files = string(filenames) + '.mat'; % adds the extenstion
l=length(files);
for i = 1:l %loading data
intenses{i}=load(files(i), 'icell_av');
holder=intenses{i}.icell_av;
avmaxintensities=[avmaxintensities holder{1}];
 
icells{i}=load(files(i), 'icell');
allcellis=icells{i}.icell;
allintensities=[allintensities allcellis{1}];
SDs=[SDs nanstd(allintensities{i})];

Ts(i)=load(files(i), 'time');
times=[times Ts(i).time];
end

% going to try filtering out cells found at only one or two time points
%maybe filter those that have diffs above a value
for i=1:length(allintensities)
   for t=1:length(allintensities{i}) 
    celltrace=allintensities{i}(t,:);
    celltracenan = celltrace(~isnan(celltrace));
    tps=length(celltracenan);
    verify= tps >= 1;
    if verify == 0
        allintensities{i}(t,:) = NaN;
    end  
    
    changecheck=diff(allintensities{i}(t,:));
    maxdiff=max(changecheck);
    mindiff=min(changecheck);
    verify = maxdiff > 700;
    if verify == 1;
         allintensities{i}(t,:) = NaN;
    end
       verify = mindiff > -300;
    if verify == 0;
         allintensities{i}(t,:) = NaN;
    end
   
    end
end

%Partition data from same experiments
zero=[allintensities{1}];
one25=[allintensities{2} ; allintensities{3}];
one5=[allintensities{4} ; allintensities{5} ; allintensities{6}];
two={allintensities{7} allintensities{8}  allintensities{9} allintensities{10}};
two25={allintensities{11} ; allintensities{12} ; allintensities{13}}; %allintensities{16} allintensities{17}};

% time too long, add nans to data that's too short

for i =1:length(two)
    framelen=size(two{4},2);
    celllen=size(two{i},2);
    verify = celllen == framelen;
    if verify == 0 
   newmat=NaN(size(two{i},1),size(two{4},2));
        newmat(1:size(two{i},1),1:size(two{i},2))=two{i};
        two{i}=newmat;
    end
end

%Partition data from same experiments with new lengths
two=[two{1} ; two{2} ; two{3} ; two{4}];
two25=[two25{1} ; two25{2};  two25{3}];

%tidy 1.25 by remove slope outliers and negatives
one25changes=diff(one25');
one25changes=one25changes'; 
one25change=nanmean(one25changes'); %these are the, slope change per cell
avone25change=nanmean(one25change);
for i =1:length(one25)
    z=one25change(i)/avone25change;
    verify = z > 2.5;
        if verify == 1;
         one25(i,:) = NaN;
        end
            verify = one25change(i) < 0;
        if verify == 1;
         one25(i,:) = NaN;
        end
end

 % I want to tidy 1.5
one5changes=diff(one5');
one5changes=one5changes'; 
one5change=nanmean(one5changes'); %these are the slope change per cell
avone5change=nanmean(one5change);
for i =1:length(one5)
    z=one5change(i)/avone5change;
    verify = z > 2.5;
        if verify == 1;
         one5(i,:) = NaN;
        end
                    verify = one5change(i) < 0;
        if verify == 1;
         one5(i,:) = NaN;
        end
end

% tidy 2 too
twochanges=diff(two');
twochanges=twochanges'; 
twochange=nanmean(twochanges'); %these are the slope change per cell
avtwochange=nanmean(twochange);
for i =1:length(two)
    z=twochange(i)/avtwochange;
    verify = z > 2.5;
        if verify == 1;
         two(i,:) = NaN;
        end
                    verify = twochange(i) < 0;
        if verify == 1;
         two(i,:) = NaN;
        end
end

%2.25
two25changes=diff(two25');
two25changes=two25changes'; 
two25change=nanmean(two25changes'); %these are the slope change per cell
avtwo25change=nanmean(two25change);
for i =1:length(two25)
    z=two25change(i)/avtwo25change;
    verify = z > 2.5;
        if verify == 1;
         two25(i,:) = NaN;
        end
                    verify = two25change(i) < 0;
        if verify == 1;
         two25(i,:) = NaN;
        end
end



zerosd=nanstd(zero);
one25sd=nanstd(one25);
one5sd=nanstd(one5);
twosd=nanstd(two);
two25sd=nanstd(two25);

avzero=nanmean(zero);
avone25=nanmean(one25);
avone5=nanmean(one5);
avtwo=nanmean(two);
avtwo25=nanmean(two25);

figure
shadedErrorBar(times{1},avzero,zerosd,'lineprops','r')
hold on
shadedErrorBar(times{6},avone25,one25sd,'lineprops','g')
shadedErrorBar(times{6},avone5,one5sd,'lineprops','b')
shadedErrorBar(times{6},avtwo,twosd,'lineprops','y')
shadedErrorBar(times{6},avtwo25,two25sd,'lineprops','h')
fig2pretty
hold off
legend('0','1.25','1.5','2','2.25')

normexp1=normalize(avzero,"scale","first");
normexp2=normalize(avone25,"scale","first");
normexp3=normalize(avone5,"scale","first");
normexp4=normalize(avtwo,"scale","first");
normexp5=normalize(avtwo25,"scale","first");

%normalize the error - divide by each time point data?
normsd1=nanstd(zero)./avzero;
normsd2=nanstd(one25)./avone25;
normsd3=nanstd(one5)./avone5;
normsd4=nanstd(two)./avtwo;
normsd5=nanstd(two25)./avtwo25;

figure
shadedErrorBar(times{1},normexp1,normsd1,'lineprops','r')
hold on
shadedErrorBar(times{6},normexp2,normsd2,'lineprops','g')
shadedErrorBar(times{6},normexp3,normsd3,'lineprops','b')
shadedErrorBar(times{6},normexp4,normsd4,'lineprops','m')
shadedErrorBar(times{6},normexp5,normsd5,'lineprops','k')
fig2pretty
hold off
ylabel('Normalized Fluorescence (A.U.)')
xlabel('Time (mins)')
legend('0','1.25','1.5','2','2.25','location','northwest')

zerochanges=diff(zero');
zerochanges=zerochanges'; 
zerochange=nanmean(zerochanges'); %these are the slope change per cell
zerochange = zerochange(~isnan(zerochange));
avzchange=nanmean(zerochange)/16.333; %this is change per minute
normzchange=avzchange/avzero(1);
avzsd=nanstd(zerochange)/16.333;

% lets do sem - need accurate numbers of cells to square by
stderrszero=nanstd(zero)/sqrt(length(zero));
normseszero=stderrszero/avzero(1);

for i = 1:size(one25,2)
    col=one25(:,i);
    col = col(~isnan(col));
    lengthcol=length(col);
stderrsone25(i)=nanstd(col)/sqrt(lengthcol);
end
normsesone25=stderrsone25/avone25(1);

for i = 1:size(one5,2)
    col=one25(:,i);
    col = col(~isnan(col));
    lengthcol=length(col);
stderrsone5(i)=nanstd(col)/sqrt(lengthcol);
end
normsesone5=stderrsone5/avone5(1);

for i = 1:size(two,2)
    col=one25(:,i);
    col = col(~isnan(col));
    lengthcol=length(col);
stderrstwo(i)=nanstd(col)/sqrt(lengthcol);
end
normsestwo=stderrstwo/avtwo(1);

for i = 1:size(two25,2)
    col=one25(:,i);
    col = col(~isnan(col));
    lengthcol=length(col);
stderrstwo25(i)=nanstd(col)/sqrt(lengthcol);
end
normsestwo25=stderrstwo25/avtwo25(1);

scaler=150/30;
scale1=0:50:150;
scale2=0:scaler:150;
%normexp1=normexp1.*scale1;
%normseszero=normseszero.*scale1;
%normexp2=normexp2.*scale2;
%normsesone25=normsesone25.*scale2;
%normexp3=normexp3.*scale2;
%normsesone5=normsesone5.*scale2;
%normexp4=normexp4.*scale2;
%normsestwo=normsestwo.*scale2;
%normexp5=normexp5.*scale2;
%normsestwo25=normsestwo25.*scale2;


figure
shadedErrorBar(times{1},normexp1,normseszero,'lineprops','r')
hold on
shadedErrorBar(times{6},normexp2,normsesone25,'lineprops','g')
shadedErrorBar(times{6},normexp3,normsesone5,'lineprops','b')
shadedErrorBar(times{6},normexp4,normsestwo,'lineprops','m')
shadedErrorBar(times{6},normexp5,normsestwo25,'lineprops','k')
fig2pretty
hold off
ylabel('Single cell fluorescence (A.U.)')
xlabel('Time (mins)')
lgd = legend('0','0.63','0.75','1','1.25','location','northwest');
title(lgd,'Distance piston lowered (cm)')



one25changes=diff(one25');
one25changes=one25changes'; 
one25change=nanmean(one25changes'); %these are the, slope change per cell
one25change = one25change(~isnan(one25change));
avone25change=nanmean(one25change)/2;
normone25change=avone25change/avone25(1);
avone25sd=nanstd(one25change)/2;

one5changes=diff(one5');
one5changes=one5changes'; 
one5change=nanmean(one5changes'); %these are the slope change per cell
one5change = one5change(~isnan(one5change));
avone5change=nanmean(one5change)/2;
normone5change=avone5change/avone5(1);
avone5sd=nanstd(one5change)/2;

twochanges=diff(two');
twochanges=twochanges'; 
twochange=nanmean(twochanges'); %these are the slope change per cell
twochange = twochange(~isnan(twochange));
avtwochange=nanmean(twochange)/2;
normtwochange=avtwochange/avtwo(1);
avtwosd=nanstd(twochange)/2;

two25changes=diff(two25');
two25changes=two25changes'; 
two25change=nanmean(two25changes'); %these are the slope change per cell
two25change = two25change(~isnan(two25change));
avtwo25change=nanmean(two25change)/2;
normtwo25change=avtwo25change/avtwo25(1);
avtwo25sd=nanstd(two25change)/2;

figure
errorbar(1,avzchange,avzsd)
hold on
errorbar(2,avone25change,avone25sd)
errorbar(3,avone5change,avone5sd)
errorbar(4,avtwochange,avtwosd)
errorbar(5,avtwo25change,avtwo25sd)
xlim([0 6])
ylabel(['change per min'])
hold off

% now I want normalized rate of increase with error
%std of all changes - this is affected by time change
normzsd=nanstd(zerochange)/avzero(1)/16.333;
normone25sd=nanstd(one25change)/avone25(1)/2;
normone5sd=nanstd(one5change)/avone5(1)/2;
normtwosd=nanstd(twochange)/avtwo(1)/2;
normtwo25sd=nanstd(two25change)/avtwo25(1)/2;


figure
errorbar(1,normzchange,normzsd,'-o')
hold on
errorbar(2,normone25change,normone25sd,'-o')
errorbar(3,normone5change,normone5sd,'-o')
errorbar(4,normtwochange,normtwosd,'-o')
errorbar(5,normtwo25change,normtwo25sd,'-o')
xlim([0 6])
%ylim([-0.06 0.07])
hold off

% standard error
stderrzero=nanstd(zerochange)/sqrt(length(zerochange));
normsezero=stderrzero/avzero(1);
stderrone25=nanstd(one25change)/sqrt(length(one25change));
normseone25=stderrone25/avone25(1);
stderrone5=nanstd(one5change)/sqrt(length(one5change));
normseone5=stderrone5/avone5(1);
stderrtwo=nanstd(twochange)/sqrt(length(twochange));
normsetwo=stderrtwo/avtwo(1);
stderrtwo25=nanstd(two25change)/sqrt(length(two25change));
normsetwo25=stderrtwo25/avtwo25(1);

%normzchange=normzchange*10000;
%normsezero=normsezero*10000;
%normone25change=normone25change*10000;
%normseone25=normseone25*10000;
%normone5change=normone5change*10000;
%normseone5=normseone5*10000;
%normtwochange=normtwochange*10000;
%normsetwo=normsetwo*10000;
%normtwo25change=normtwo25change*10000;
%normsetwo25=normsetwo25*10000;

set(groot,'defaultLineMarkerSize',8)
figure
errorbar(1,normzchange,normsezero,'-ok')
hold on
aa=scatter(1,normzchange, 'o', "MarkerEdgeColor","k", ...
'MarkerFaceColor', 'r')
errorbar(2,normone25change,normseone25,'-ok')
bb=scatter(2,normone25change,'o', "MarkerEdgeColor","k", ...
    'MarkerFaceColor', 'g')
errorbar(3,normone5change,normseone5,'-ok')
cc=scatter(3,normone5change,'o', "MarkerEdgeColor","k", ...
    'MarkerFaceColor', 'b')
errorbar(4,normtwochange,normsetwo,'-ok')
dd=scatter(4,normtwochange,'o', "MarkerEdgeColor","k", ...
    'MarkerFaceColor', 'm')
errorbar(5,normtwo25change,normsetwo25,'-ok')
ee=scatter(5,normtwo25change,'o', "MarkerEdgeColor","k", ...
    'MarkerFaceColor', 'k')
xlim([0 6])
xticks([1 2 3 4 5])
xticklabels({'0','0.63','0.75','1','1.25'})
fig2pretty
plot(1:2,[normzchange normone25change],'--k')
plot(2:5,[normone25change normone5change normtwochange normtwo25change],'k')
ylabel('Reporter induction rate (A.U.)')
xlabel('Distance piston lowered (cm)')
uistack(aa, 'top')
uistack(bb, 'top')
uistack(cc, 'top')
uistack(dd, 'top')
uistack(ee, 'top')
hold off

%correct for photobleaching here
photobleach1=avone25(1)-avone25(2);
avone252=[avone25(1) avone25(2:end)+photobleach1];
normexp2=normalize(avone252,"scale","first");
photobleach2=avone5(1)-avone5(2);
avone52=[avone5(1) avone5(2:end)+photobleach2];
normexp32=normalize(avone52,"scale","first");

figure
shadedErrorBar(times{1},normexp1,normseszero,'lineprops','r')
hold on
shadedErrorBar(times{6},normexp2,normsesone25,'lineprops','g')
shadedErrorBar(times{6},normexp32,normsesone5,'lineprops','b')
shadedErrorBar(times{6},normexp4,normsestwo,'lineprops','m')
shadedErrorBar(times{6},normexp5,normsestwo25,'lineprops','k')
fig2pretty
hold off
ylabel('Single cell fluorescence (A.U.)')
xlabel('Time (mins)')
lgd = legend('0','0.63','0.75','1','1.25','location','northwest');
title(lgd,'Distance piston lowered (cm)')