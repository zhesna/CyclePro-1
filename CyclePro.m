function CyclePro(input_signal,sampling_frequency)


data = load(input_signal);
x = data(:,1);
y = data(:,2);
z = data(:,3);


%%sensor fusion
svm = sqrt(x.^2+y.^2+z.^2);
L = length(svm);
%figure;
%plot(svm);
index = [1:L];

%%salient points
salient = zeros(1,L);
for i = 1:L-1
    mi = 0;
    for j = i+1:L
        if svm(j) > svm(i)
            mi = mi+1;
        else
            break
        end
    end
    salient(i) = mi;
end

n = 0;
for i = 1:L
    if salient(i) > floor((60.*sampling_frequency)/100)
        n = n + 1;
    end
end

minipeak = zeros(1,n);
indpeak = zeros(1,n);
m = 1;
for i = 1:L
    if salient(i) > floor((60.*sampling_frequency)/100)
        minipeak(m) = salient(i);
        indpeak(m) = index(i);
        m = m+1;
    end
end

eliminateP = 0;
for i = 2:n
    if indpeak(i) - indpeak(i-1) < floor((60.*sampling_frequency)/100)
        indpeak(i) = indpeak(i-1);
        minipeak(i) = 0;
        eliminateP = eliminateP+1;
    end
end
k = n - eliminateP - 1;
finalpeak = zeros(1,k);
finalIndex = zeros(1,k);
q = 1;
for i = 2:n
    if minipeak(i)>0
        finalpeak(q) = minipeak(i);
        finalIndex(q) = indpeak(i);
        q = q+1;
    end
end

x = length(finalIndex);
%%template generation
answer = TemplateGeneration(svm,finalIndex,x);
aryPeak = zeros(1,1);
cadMean = zeros(1,1);
for i = 1:length(answer)
    temp = svm(answer(i,1):answer(i,2));
%%pattern identification
    corResult = correlation(temp, svm);
%%optimal stride detection
    [numOfPeak,Peak,Index] = newPeak(corResult,sampling_frequency);
    aryPeak(i) = numOfPeak;
    cadMean(i) = cadenceCal(Peak,Index,sampling_frequency);
end
Number_of_Strides = round(mean(aryPeak));
Cadence = mean(cadMean);


disp('Number of Strides=')
disp(Number_of_Strides)
disp('cadence=')
disp(Cadence)




%%%function to convert csv file to dat file

% function loadFile(filename)
% data = csvread(filename);
%  x = data(:,1);%6 or 9 or 1
%  y = data(:,2);%7 or 10 or 1
%  z = data(:,3);%9 or 11 or 1
%  dlmwrite('example.dat',[x y z]);
% end




function answer = TemplateGeneration(signal,cutAry,x)

numList = 1:1:x;
if x > 4
    numList = numList(2:end-2);
    m = cutAry(1)+1;
else
    numList = numList(1:end-1);
    m = 1;
end
diff = zeros(length(numList),3);

for i=1:length(numList)
    k = numList(i);
    if cutAry(k)<0 || cutAry(k)>length(signal)
        break;
    end
    a = signal(m:cutAry(k));
    
    minD = abs(min(a(1:3))) - abs(min(a(end-3:end)));
    diff(i,1) = (std(signal) - std(a))^2 + minD^2;
    diff(i,2) = m;
    diff(i,3) = cutAry(k);
    m = cutAry(k)+1;
end

[value,index] = sort(diff(:,1),'ascend');
answer = zeros(1,2);
for opt =1:length(index)
    answer(opt,1) = floor(diff(index(opt),2));
    answer(opt,2) = floor(diff(index(opt),3));
    sel = signal(answer(opt,1):answer(opt,2));
    if opt == 3
        break;
    end
    %figure;
    %plot(sel);
end




function cor = correlation(window, ary)

lw = length(window);
la = length(ary);
raw_cor = zeros(1,la-lw+1);
lastP = la-lw+1;

for i = 1:lastP
    aryInWindow = ary(i:(i+lw-1));
    score = sum(times(window,aryInWindow));
    raw_cor(i) = score;
end
cor = raw_cor/(max(abs(raw_cor)));










function [minNum,finalP,finalI] = newPeak(signal,s)

L = length(signal);
index = [1:L];

tempPeak = [];
tempIndex = [];
n = 1;
%%local maximum points
for i = 2:L-1
    if signal(i) > signal(i-1) && signal(i) > signal(i+1)
        tempPeak(n) = signal(i);
        tempIndex(n) = index(i);
        n = n+1;
    end
end

eliminateP = 0;
i = 1;
n = n-1;

%%anchor window
while i < n
    j = i+1;
    maxT = tempPeak(i);
    maxI = i;
    while ( j <= n && tempIndex(j) - tempIndex(i) < floor((60.*s)/100))
        if maxT <= tempPeak(j)
            maxT = tempPeak(j);
            tempPeak(maxI) = 0;
            eliminateP = eliminateP+1;
            maxI = j;
        else
            tempPeak(j) = 0;
            eliminateP = eliminateP + 1;
        end 
        j = j+1;
    end
    i = j;
end

tPeak = [];
tIndex = [];
k = 1;
for i = 1:n
    if tempPeak(i) > 0
        tPeak(k) = tempPeak(i);
        tIndex(k) = tempIndex(i);
        k = k+1;
    end
end
k = k-1;

%%distance filter

value = tPeak(1);
in = tIndex(1);
for i = 2:k
    if tIndex(i) - in < floor((60.*s)/70)
        if tPeak(i) > value
            value = tPeak(i);
            in = tIndex(i);
            tPeak(i-1) = 0;
        else
            value = tPeak(i-1);
            in = tIndex(i-1);
            tPeak(i) = 0;
        end
    else
        in = tIndex(i);
        value = tPeak(i);
    end
end
          
Peak = [];
Index = [];
p = 1;
for i = 1:k
    if tPeak(i) > 0
        Peak(p) = tPeak(i);
        Index(p) = tIndex(i);
        p = p+1;
    end
end
p = p-1;

%%optimal peak separator
[valInd,valVal] = findValley(signal,Index);
valueBar = min(valVal);

[value,index] = sort(Peak,'descend');
valueTemp = [value,valueBar];
indexTemp = [index,p+1];
minV = var(valueTemp);
minNum = 1;

for i=1:p
    if var(valueTemp(1:i))+var(valueTemp(i+1:end)) <= minV
        minV = var(valueTemp(1:i))+var(valueTemp(i+1:end));
        minNum = i;
    end
end

finalP = zeros(minNum,1);
finalI = zeros(minNum,1);
for j = 1:minNum
    finalP(j) = Peak(indexTemp(j));
    finalI(j) = Index(indexTemp(j));
end

if minNum < p  
    deNum = p-minNum;
    deleteP = zeros(deNum,1);
    deleteI = zeros(deNum,1);
    for j=minNum+1:p
        deleteP(j-minNum) = Peak(indexTemp(j));
        deleteI(j-minNum) = Index(indexTemp(j));
    end
else
    deleteP = [];
    deleteI = [];
end


figure;
plot(signal);
hold on;
plot(finalI,finalP,'ro','MarkerSize',6);
hold on;
plot(deleteI,deleteP,'ko','MarkerSize',6);
hold on;
plot(valInd,valVal,'m+','MarkerSize',8);



function [valleyIndex, valleyValue] = findValley(signal,index)
L = length(index);

valleyIndex = zeros(1,L-1);
valleyValue = zeros(1,L-1);
for i=1:L-1
    valleyValue(i) = min(signal(index(i):index(i+1)));
    valleyIndex(i) = find(signal(index(i):index(i+1)) == valleyValue(i),1) + index(i)-1;
end




function cad = cadenceCal(peak,index,s)
indexO = sort(index);
l = length(peak);
array = [];
for i = 2:l
    step = times(60,s)/(indexO(i)-indexO(i-1));
    tempC = step.*2;
    array(i-1) = roundNum(tempC);
end
cad = mean(array);%%%cadence bar asase 2 peake avalie select shode(mean)????



function rNum = roundNum(num)
temp = num2str(num,4);
rNum = str2num(temp);


