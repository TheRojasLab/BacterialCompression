% This is a script to get the average height profiles across Z stacks

clear, close all

dirname= ['C:\Users\Guyma\Documents\Rojas Lab\MatLAB WorkSpaces\Z stacks'];
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
lfl=length(files); 
av_counts=cell(lfl,1); %creates cell array, is more lenient than a matrix

for i = 1:lfl;
av_counts{i}=load([files{i}],'averages');
zz_counts{i}=load([files{i}],'z');
   avg_counts{i}=av_counts{i}.averages;
   z_counts{i}=zz_counts{i}.z;
end
%editing first z stack to get uM
avg_counts{1}=avg_counts{1}/1000;

% Need to interpolate to fit the axis #1 is the shortest so will use that as
% the x axis, going to fit them to the tallest ends?
xlength=length(z_counts{1});
onecutoff=length(z_counts{2});
z_counts{2}=z_counts{2}((onecutoff-xlength)+1:end);
avg_counts{2}=avg_counts{2}((onecutoff-xlength)+1:end);
twocutoff=length(z_counts{3});
z_counts{3}=z_counts{3}((twocutoff-xlength)+1:end);
avg_counts{3}=avg_counts{3}((twocutoff-xlength)+1:end);
threecutoff=length(z_counts{4});
z_counts{4}=z_counts{4}((threecutoff-xlength)+1:end);
avg_counts{4}=avg_counts{4}((threecutoff-xlength)+1:end);

heightavs=[];
for i =1:lfl
    height(i,:)=avg_counts{i};
end

heightavs=mean(height);
heightsd=std(height);

for i =1:lfl
   plot(z_counts{1},avg_counts{i}) 
   hold on
   plot(z_counts{1},heightavs,'o-')
end
hold off