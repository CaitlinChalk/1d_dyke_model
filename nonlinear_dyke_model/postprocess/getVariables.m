function [z,h,Pe,zf,hmax,t] = getVariables(res_dir, filename,dt)

%function to load variables from output .txt files
%input - res_dir - string of directory containing results
%      - filename - filename of results, prior to the saving index

files = dir([res_dir filename '*.txt']);
n = length(files);

%intialise arrays
t = zeros(n,1);
hmax = zeros(n,1);
zf = zeros(n,1);
h = cell(1,n);
Pe = cell(1,n);
z = cell(1,n);
sizeA = [3 Inf];

for i = 1:n
    thisfile = files(i).name;
    %extract saving index at the end of the file name
    endIndex = regexp(thisfile,'.txt')-1;
    startIndex = regexp(thisfile,'_');
    startIndex = startIndex(2)+1;
    %startIndex = regexp(thisfile,'_i'); change to this later
    isave = str2double(thisfile(startIndex:endIndex));
    t(i) = isave*dt;
    fileID = fopen([res_dir thisfile],'r');
    A = fscanf(fileID,'%f',sizeA);
    A = A'; %z,h,Pe
    %save dyke height, width and pressure
    ztemp = A(:,1);
    z{i} = ztemp;
    zf(i) = ztemp(end);
    h{i} = A(:,2);
    hmax(i) = max(A(:,2));
    Pe{i} = A(:,3);
end

%reorder arrays in ascending order of time
[t, ind] = sort(t);
z = z(ind);
zf = zf(ind);
h = h(ind);
hmax = hmax(ind);
Pe = Pe(ind);

end
