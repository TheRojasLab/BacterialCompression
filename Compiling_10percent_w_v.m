% Script for compiling all the torture chamber data
clear, close all

% This is 
channels=['C:\Users\GuyMa\Documents\Rojas Lab\MatLAB Workspaces\FLUO/Smush\25_4_23'];
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
maxintenses={};
for i = 1:l %loading data
maxintenses{i}=load(files(i), 'icell_av');
holder=maxintenses{i}.icell_av;
avmaxintensities{i}=holder{1};
 
icells{i}=load(files(i), 'icell');
allcellis=icells{i}.icell;
allintensities{i}=allcellis{1};
SDs{i}=nanstd(allintensities{i});

Ts(i)=load(files(i), 'time');
times{i}=Ts(i).time;
end

% going to try filtering out cells found at only one or two time points
%maybe filter those that have diffs above a value
for i=1:length(allintensities)
   for t=1:length(allintensities{i}) 
    celltrace=allintensities{i}(t,:);
    celltracenan = celltrace(~isnan(celltrace));
    tps=length(celltracenan);
    verify= tps > 3;
    if verify == 0
        allintensities{i}(t,:) = NaN;
    end  
    %custom for 2.25 turns
    if i==11
        verify= tps > 6;
    if verify == 0
        allintensities{i}(t,:) = NaN;
    end 
    end
     if i==10
        verify= tps > 6;
    if verify == 0
        allintensities{i}(t,:) = NaN;
    end 
     end
         if i==8
        verify= tps > 6;
    if verify == 0
        allintensities{i}(t,:) = NaN;
    end 
    end
    
    changecheck=diff(allintensities{i}(t,:));
    maxdiff=max(changecheck);
    mindiff=min(changecheck);
    verify = maxdiff < 400;
    if verify == 0;
         allintensities{i}(t,:) = NaN;
    end
       verify = mindiff > -300;
    if verify == 0;
         allintensities{i}(t,:) = NaN;
    end
   end
end


%Partition data from same experiments
exp1=[allintensities{1} ; allintensities{2} ; allintensities{3}];
exp2=[allintensities{4} ; allintensities{5}];
exp3=[allintensities{6} ; allintensities{7}];
exp4=[allintensities{8}];
exp5=[allintensities{9}];
exp6=[allintensities{10}; allintensities{11}];


exp1sd=nanstd(exp1);
exp2sd=nanstd(exp2);
exp3sd=nanstd(exp3);
exp4sd=nanstd(exp4);
exp6sd=nanstd(exp6);

avexp1=nanmean(exp1);
avexp2=nanmean(exp2);
avexp3=nanmean(exp3);
avexp4=nanmean(exp4);
avexp5=nanmean(exp5);
avexp6=nanmean(exp6);

figure
hold on
for i = 1:length(times)
errorbar(times{i},avmaxintensities{i},SDs{i})
hold on
legend('1.1','1.2','1.3','2.1','2.2','3.1','3.2','4','5','6')
end
hold off

figure
hold on
errorbar(times{1},avexp1,exp1sd)
errorbar(times{4},avexp2,exp2sd)
errorbar(times{6},avexp3,exp3sd)
errorbar(times{8},avexp4,exp4sd)
errorbar(times{9},avmaxintensities{9},SDs{9})
errorbar(times{10},avexp6,exp6sd)
legend('1.5','2','1.25','2.25','0','2.25')
hold off

% now get absolute change
normexp1=normalize(avexp1,"scale","first");
normexp2=normalize(avexp2,"scale","first");
normexp3=normalize(avexp3,"scale","first");
normexp4=normalize(avexp4,"scale","first");
normexp5=normalize(avmaxintensities{9},"scale","first");
normexp6=normalize(avexp6,"scale","first");


figure
plot(times{4},normexp2)
hold on
plot(times{10},normexp6)
plot(times{6},normexp3)
plot(times{8},normexp4)
plot(times{1},normexp1)
plot(times{9},normexp5)
legend('2','2.25','1.25','2.25','1.5','0')
hold off

% induction rate - get slopes from individual traces for error bars
% get slope from normalized data

for i=1:length(allintensities)
changes{i}=diff(allintensities{i}');
changes{i}=changes{i}'; 
avchange{i}=nanmean(changes{i}'); %these are the sds, slope change per cell
%sdavchange(i)=nanstd(avchange{i});
end

%combine same traces - these are individual slopes
exp1traces=[avchange{1} avchange{2} avchange{3}];
exp2traces=[avchange{4} avchange{5}];
exp3traces=[avchange{6} avchange{7}];
exp4traces=[avchange{8}];
exp5traces=[avchange{9}];
exp6traces=[avchange{10} avchange{11}];

alltraces{1}=exp1traces;
alltraces{2}=exp2traces;
alltraces{3}=exp3traces;
alltraces{4}=exp4traces;
alltraces{5}=exp5traces;
alltraces{6}=exp6traces;

% some of these are too big, so lets filter the extremes
for i=1:length(alltraces)
   y=nanmean(alltraces{i}); 
    for t=1:length(alltraces{i})
      verify =  abs(alltraces{i}(t)/y); %y > alltraces{i}(t);
        check= verify > 3;
        if check == 1
        %   if i==6
          % exp6(t,:)= NaN;
           alltraces{i}(t) = NaN; 
        %   end
        end
           %Filter negatives - there are no negatives in the data of 6
             verify = alltraces{i}(t);   
            check = verify < 0;
            if check == 1
                if i==6
           exp6(t,:)= NaN;
           alltraces{i}(t) = NaN; 
                end
            end
       end
end

    for i = 1:length(alltraces)
          alltraces{i} = alltraces{i}(~isnan(alltraces{i})) ;
    end
    %think I need to redefine data from cleaned - just remove nans from
%alltraces

%these are the SDs
slopesd1=nanstd(alltraces{1});
slopesd2=nanstd(alltraces{2});
slopesd3=nanstd(alltraces{3});
slopesd4=nanstd(alltraces{4});
slopesd5=nanstd(alltraces{5});
slopesd6=nanstd(alltraces{6});

%normalize by inital value - need to combine traces from same experiment
%sdavchange(1)=sdavchange(1)/
ratesofchange1=diff(normexp1)./diff(times{1});
rate1=nanmean(ratesofchange1);
ratesofchange2=diff(normexp2)./diff(times{4});
rate2=nanmean(ratesofchange2);
ratesofchange3=diff(normexp3)./diff(times{6});
rate3=nanmean(ratesofchange3);
ratesofchange4=diff(normexp4)./diff(times{8});
rate4=nanmean(ratesofchange4);
ratesofchange5=diff(normexp5)./diff(times{9});
rate5=nanmean(ratesofchange5);
ratesofchange6=diff(normexp6)./diff(times{10});
rate6=nanmean(ratesofchange6);

normsd1=slopesd1/avexp1(1);
normsd2=slopesd2/avexp2(1);
normsd3=slopesd3/avexp3(1);
normsd4=slopesd4/avexp4(1);
normsd5=slopesd5/avexp5(1);
normsd6=slopesd6/avexp6(1);
 

x=1:6;
figure
errorbar(x(1),rate5,normsd5) % this is 0
hold on
errorbar(x(2),rate3,normsd3) %1.25
errorbar(x(3),rate1,normsd1) %1.5
errorbar(x(4),rate2,normsd2) % this is 2
errorbar(x(5),rate4,normsd4) %2.25
errorbar(x(6),rate6,normsd6) %2.25
xlim([0 7])
ylim([-0.05 0.08])
%ylim([-0.1 0.1])
xticks([1 2 3 4 5 6])
xticklabels({'0','1.25','1.5','2','2.25','2.25'})
hold off

%Combine both two turns
exp4=[allintensities{8}; allintensities{10}; allintensities{11}];
avexp4=nanmean(exp4);
alltraces{4}=[alltraces{4}  alltraces{6}];
slopesd4=nanstd(alltraces{4});
normexp4=normalize(avexp4,"scale","first");
ratesofchange4=diff(normexp4)./diff(times{8});
rate4=nanmean(ratesofchange4);
normsd4=slopesd4/avexp4(1);

figure
errorbar(x(1),rate5,normsd5,'*') % this is 0
hold on
errorbar(x(2),rate3,normsd3,'*') %1.25
errorbar(x(3),rate1,normsd1,'*') %1.5
errorbar(x(4),rate2,normsd2,'*') % this is 2
errorbar(x(5),rate4,normsd4,'*') %2.25
xlim([0 6])
ylim([-0.05 0.08])
%ylim([-0.1 0.1])
xticks([1 2 3 4 5])
xticklabels({'0','0.63','0.75','1','1.25'})
ylabel('Rate of Induction (min^-^1)')
xlabel('Distance Screw Lowered (cm)')
hold off
fig2pretty
 
%this is the sd 
exp1sd=nanstd(exp1);
exp2sd=nanstd(exp2);
exp3sd=nanstd(exp3);
exp4sd=nanstd(exp4);
exp5sd=nanstd(exp5);

%normalize it
norm1sd=exp1sd/avexp1(1);
norm2sd=exp2sd/avexp1(2);
norm3sd=exp3sd/avexp1(3);
norm4sd=exp4sd/avexp1(4);
norm5sd=exp5sd/avexp1(5);

figure
errorbar(times{9},normexp5,norm5sd,'-o')
hold on
errorbar(times{6},normexp3,norm3sd,'-o')
errorbar(times{1},normexp1,norm1sd,'-o')
errorbar(times{4},normexp2,norm2sd,'-o')
errorbar(times{8},normexp4,norm4sd,'-o')
xlabel('Time (Mins)')
ylabel('Normalized Fluorescence')
hold off
fig2pretty
legend('0cm','0.625cm','0.75cm','1cm','1.25cm','Location','northwest')

%Lets try the standard error - std( data ) / sqrt( length( data ))
stderr1=nanstd(alltraces{1})/sqrt(length(alltraces{1}));
normse1=stderr1/avexp1(1);
stderr2=nanstd(alltraces{2})/sqrt(length(alltraces{2}));
normse2=stderr2/avexp2(1);
stderr3=nanstd(alltraces{3})/sqrt(length(alltraces{3}));
normse3=stderr3/avexp3(1);
stderr4=nanstd(alltraces{4})/sqrt(length(alltraces{4}));
normse4=stderr4/avexp4(1);
stderr5=nanstd(alltraces{5})/sqrt(length(alltraces{5}));
normse5=stderr5/avexp5(1);

figure
errorbar(x(1),rate5,normse5,'*') % this is 0
hold on
errorbar(x(2),rate3,normse3,'*') %1.25
errorbar(x(3),rate1,normse1,'*') %1.5
errorbar(x(4),rate2,normse2,'*') % this is 2
errorbar(x(5),rate4,normse4,'*') %2.25
xlim([0 6])
%ylim([-0.05 0.08])
%ylim([-0.1 0.1])
xticks([1 2 3 4 5])
xticklabels({'0','0.63','0.75','1','1.25'})
ylabel('Rate of Induction (min^-^1)')
xlabel('Distance Screw Lowered (cm)')
hold off
fig2pretty



