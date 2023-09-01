clear, close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%USER INPUT
dirname= ['C:\Users\Guyma\Documents\Rojas Lab\MatLAB WorkSpaces\FLUO\13_8_21'];
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
  %%%%
  pixvalues=[];
  c_counts1=[];
  
  for t=1:length(files);

  load(files(t), 'meanintensities'); %loads the  data
  c_counts{t,1} = meanintensities; % Allocates it into a cell
  load(files(t), 'holdcells');
  c_counts1{t,1} = holdcells;
y=c_counts{t}; % for easier manipulation, moving it into a matrix

  end

  dirname= ['C:\Users\Guyma\Documents\Rojas Lab\MatLAB WorkSpaces\FLUO\14_6_21_z_stack'];
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
  %%%%

l=length(files);
load(files, 'z');
load(files, 'averages');
stack=z;
averages=averages/1000;

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

  for tt=1:length(files);

  load(files(tt), 'meanintensities'); %loads the  data
  c_counts{t+tt,1} = meanintensities; % Allocates it into a cell
  load(files(tt), 'holdcells');
  c_counts1{t+tt,1} = holdcells;
y=c_counts{t+tt}; % for easier manipulation, moving it into a matrix


  end

    dirname= ['C:\Users\Guyma\Documents\Rojas Lab\MatLAB WorkSpaces\FLUO\16_11_21'];
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

  for ttt=1:length(files);

  load(files(ttt), 'meanintensities'); %loads the  data
  c_counts{t+tt+ttt,1} = meanintensities; % Allocates it into a cell
  load(files(ttt), 'holdcells');
  c_counts1{t+tt+ttt,1} = holdcells;
y=c_counts{t+tt+ttt}; % for easier manipulation, moving it into a matrix
  end

  
  
  
figure % Fig 1   
yyaxis left      
plot(z,averages)

xlabel('X Position (\mum)')
%set(gca,'XTick',[])
ylabel('Height (\mum)')
hold on
grid on  
set(gca,'FontSize',16)  
  %Get the errors from c_counts1
u=[];
error=[];
for i=1:length(c_counts1)
u = c_counts1{i};
    ui=u(:,end);
        stderror(i)=nanstd(ui);
end

%xspots=[13487 13313 13028 12591 12284 12190 11987]; This is the data from
%the metadata
xspots = [13313 13028 12591 12284 12190];  %This is the data I want to use right now
xspots=flip(xspots);
xspots = xspots-11987;

xspots=[xspots 100 100];
%Get the mean values for the final time point of the shortest array

for i=1:length(c_counts) %Ccounts has the mean values already
   Dims=c_counts{i};
   Lengths(i)= length(Dims);
    q=min(Lengths);
    aufluors(i)=c_counts{i}(q);
end


idx=xspots;
yyaxis right
for i=1:length(c_counts{5})
errorbar(idx(1:5),aufluors(1:5),stderror(1:5),'-o');
hold on
end
ylabel('Final Fluorescence (A.U.)')
%title('Plot of Mean Fluorescence After 1 Hour of Growth')
hold off

figure
plot(z,averages)
xlabel('X Position (\mum)')
%set(gca,'XTick',[])
ylabel('Height (\mum)')
hold on
grid on
x0=10;
y0=10;
width=1850;
height=400
set(gcf,'position',[x0,y0,width,height])
hold off


names={'200\mum';'300\mum';'600\mum';'1050\mum';'1350\mum';'RcsF Deletion';'WT'};
figure %Fig 2
errorbar(aufluors,stderror)
set(gca,'xtick',[1:7],'xticklabel',names)
hold on
ylabel('Fluorescence (Arbitrary Units)')
grid on 
title('Plot of Mean Fluorescence After 1 Hour of Growth')
hold off


time=(1:q);
tscale=2;
time=(time*tscale);

% c_counts1 is values of each cell at all time points, get CI from those
%z=[];
SEM=[];
means=[];
std_dev=[];
CI95=[];
altdev=[];
%nanmean(fluo)
%nanstd(fluo)
%nanstd(fluo/sqrt(N)) N is number of non NaNs 
%. before divide means element by element


for i=1:length(c_counts1) %getting all fluor values of all cells
means{i}=nanmean(c_counts1{i});
std_dev{i}=nanstd(c_counts1{i});

N=sum(~isnan(c_counts1{i}(:,1)));
std_err{i}=nanstd(c_counts1{i}./sqrt(N));

for a=1:size(c_counts{i},2) %For each time point generate CIs
SEM{i} = std_err{i}%/sqrt(nnz(~isnan(c_counts1{i}(:,a))));
nonzeros=nnz(~isnan(c_counts1{i}(:,a))); %Find out how many cells per time point for DoF
CI95= tinv([0.025 0.975], nonzeros-1);    % Calculate 95% Probability Intervals Of t-Distribution
yCI95 = bsxfun(@times, SEM{i}(a), CI95); 
CI95n{i,a}=yCI95(1);
CI95p{i,a}=yCI95(2);

end
end


%Add CIs to data
for i=1:length(c_counts)
B=c_counts{i}
C=[CI95n{i,:}]
CI1{i}=B+C;


B=c_counts{i};
C=[CI95p{i,:}];
CI2{i}=B+C;
line_color = ['r' 'b' 'g' 'k' 'm' 'y' 'c'];
end

for i=1:length(c_counts)
x2{i} = 1:length(c_counts{i});
end

   figure % Fig 3
ciplot(CI2{1},CI1{1},x2{1},'r')
hold on
ciplot(CI1{2},CI2{2},x2{2},'b')
hold on
ciplot(CI1{3},CI2{3},x2{3},'g')
hold on
ciplot(CI1{4},CI2{4},x2{4},'k')
hold on
ciplot(CI1{5},CI2{5},x2{5},'m')
hold on
%ciplot(CI1{6},CI2{6},x2{6},'m')
%hold on
%ciplot(CI1{7},CI2{7},x2{7},'m')

%Fig 3

   figure % Fig 3
ciplot(CI2{1},CI1{1},x2{1},'r')
hold on
ciplot(CI1{2},CI2{2},x2{2},'b')
hold on
ciplot(CI1{3},CI2{3},x2{3},'g')
hold on
ciplot(CI1{4},CI2{4},x2{4},'k')
hold on
ciplot(CI1{5},CI2{5},x2{5},'m')
hold on
  for i=1:length(line_color)
plot(x2{i},c_counts{i},'Color',line_color(i))
end
fig2pretty
ylabel('Fluorescence (A.U.)')
xlabel('Time (Mins)')
text(30,c_counts{6}(end)+500,'RcsF Knockout')
text(30,c_counts{7}(end)+300,'Untrapped Cells')
ylabel('Fluorescence (A.U.)')
xlabel('Time (Mins)')
legend('Location','northwest');
lgd=legend('0.63\mum','0.71\mum','0.73\mum','0.75\mum','0.98\mum');
title(lgd,'Device Height','FontSize',12)
hold off

%normalized fluorescence over time
% divide by normalized data
for i =1:length(CI2)
normcounts{i}=normalize(c_counts{i},"scale","first");
normCI2{i}=CI2{i}./c_counts{i}(1); 
normCI1{i}=CI1{i}./c_counts{i}(1);
end

figure
ciplot(normCI2{1},normCI1{1},x2{1},'r')
hold on
ciplot(normCI1{2},normCI2{2},x2{2},'b')
hold on
ciplot(normCI1{3},normCI2{3},x2{3},'g')
hold on
ciplot(normCI1{4},normCI2{4},x2{4},'k')
hold on
ciplot(normCI1{5},normCI2{5},x2{5},'m')
hold on
  for i=1:length(line_color)
plot(x2{i},normcounts{i},'Color',line_color(i))
end
fig2pretty
ylabel('Fluorescence (A.U.)')
xlabel('Time (Mins)')
text(30,c_counts{6}(end)+500,'RcsF Knockout')
text(30,c_counts{7}(end)+300,'Untrapped Cells')
ylabel('Fluorescence (A.U.)')
xlabel('Time (Mins)')
legend('Location','northwest');
lgd=legend('0.63\mum','0.71\mum','0.73\mum','0.75\mum','0.98\mum');
title(lgd,'Device Height','FontSize',12)
hold off


%Now I want to work out some line slopes
fits=[];
ms=[];
lines=(1:37);
generated=[];

%dif function - F is 1xN vector


%df(rate)=dif(f)/dt  - df = 1xN-1 vecot
%plot rate vs time
% [rate 1 rate2...] / rate average = mean

%Thisis a loop to generate the diff of the fluorescence
dfluor=[];

for i = 1:5%(c_counts1)
    i
dfluo{i}=diff(c_counts1{i}');
dfluo{i}=dfluo{i}';
meanfluo=nanmean(dfluo{i}');
totalSD(i)=nanstd(meanfluo);
totalmean(i)=nanmean(meanfluo);

Hlimit=totalmean(i)+(2*totalSD);
Llimit=totalmean(i)-(2*totalSD);
Above_h=meanfluo>Hlimit(i);
Below_l=meanfluo<Llimit(i);
Invalid=Above_h + Below_l;
meanfluo=meanfluo(~Invalid);
cellmeans{i}=meanfluo;


% This is now filtering the traces by the rates excluded

Valid=(~Invalid)
    traces{i} = c_counts1{i}(Valid,:);
    
    
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
    traces{i} = c_counts1{i}(Valid,:);



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


for i = 1:length(traces)
   tracemeans(i,:)=nanmean(traces{i});
end
time=1:37;
figure
hold on
for i=1:5
   plot(time,tracemeans(i,:))
end
hold off


figure
hold on
for i=1:5
   plot(time,c_counts{i})
end




figure
for i =1:5
    hist(cellmeans{i})
end

h_locs=interp1(z,averages,idx(1:5));
% rates v height
plot(h_locs,totalmean)

figure
hold on


alldiffs=diff(tracemeans');
meanymeans=nanmean(alldiffs);

for i=1:5
   plot(h_locs,meanymeans)
end
hold off

for i=1:5
allSDs(i)=nanstd(cellmeans{i})
N=length(isnan(cellmeans{i}));
stderror(i) = (allSDs(i)/sqrt(N));
end

figure
hold on

   errorbar(h_locs,meanymeans,stderror(1:5))
ylim([0 250])
hold off




for i=1:length(traces) %getting all fluor values of all cells
means{i}=nanmean(traces{i});
std_dev{i}=nanstd(traces{i});

N=sum(~isnan(traces{i}(:,1)));
std_err{i}=nanstd(traces{i}./sqrt(N));

for a=1:size(traces{i},2) %For each time point generate CIs
SEM{i} = std_err{i}%/sqrt(nnz(~isnan(c_counts1{i}(:,a))));
nonzeros=nnz(~isnan(c_counts1{i}(:,a))); %Find out how many cells per time point for DoF
CI95= tinv([0.025 0.975], nonzeros-1);    % Calculate 95% Probability Intervals Of t-Distribution
yCI95 = bsxfun(@times, SEM{i}(a), CI95); 
CI95n{i,a}=yCI95(1);
CI95p{i,a}=yCI95(2);

end
end



%Add CIs to data




%Figure for plotting rate vs time
figure
hold on
for i=1:5
 hist(cellmeans{i})
    %      plot(time(1:36),cellmeans{i}(1:36))
%hold on 
end
%xlabel('Time (Mins)')
%ylabel('Rate of Induction (Mins^-^1)')
%legend('200','300','600','1000','1200','Location','northwest')
hold off




for i=1:length(c_counts) %This gives me the slope of the lines
ms(i) = (c_counts{i}(1,20)-c_counts{i}(1,10)/(time(20)-time(10)));
generated(i,:)=time*ms(i);
end

names = {'200'; '300'; '600'; '1050'; '1350'} %; 'WT'; 'RcsF'};
%Plot induction vs height
height=[1 2];
induction=ms(1:5);
%induction=flip(induction);
figure % Fig 4
plot(induction)
set(gca,'XTick',[], 'YTick', [])
xlabel('Height (\mum)')
ylabel('Rate of Induction')
set(gca,'xtick',[1:5],'xticklabel',names)
%fig2pretty





%replacing the final time point of the 70 mins WT with the first point
aufluors(8)=c_counts{8}(1);
%Getting the std error for the correct time point
for i=8
u = c_counts1{i};
    ui=u(:,1);
        stderror(i)=nanstd(ui);
end



