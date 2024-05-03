function [C,D] = BryBoschan(data_bench,data_candi)
% Calculate the turning point matching rate
% Input parameters
% data_candi: Original sequence of candidate indicators
% data_bench: Original sequence of benchmark indicators
% Return parameters
% C: Each column data represents: matching turning points, missing turning points, excess turning points, no data turning points, turning point matching rate, excess rate, average lead order, standard deviation of lead
% D: Each column data represents: original sequence of candidate indicators, original sequence of benchmark indicators, smoothed sequence of candidate indicators, smoothed sequence of benchmark indicators ...
% Turning points on the original sequence of candidate indicators, turning points on the original sequence of benchmark indicators, turning points on the smoothed sequence of candidate indicators, turning points on the smoothed sequence of benchmark indicators
% Note: Please manually adjust the legend position for the fourth subplot
% Example: [C,D] = BryBoschan(y,x)
%% Parameter definition
pp = 16; % Distance between peak (valley) and peak (valley)
pt = 6; % Distance between peak and valley
N = 6; % Step 1 defines the maximum (minimum) value within the range [-N,N] of the sequence as a local maximum (minimum)
N_ad = 4; % Step 3 adjusts the turning points of the smoothed sequence to the original sequence (or secondary smoothed sequence) within [-N_ad,N_ad]
lead = 15; % Step 4 searches for matching turning points on the candidate indicator within lead months before the benchmark indicator
lag = 8; % Step 4 searches for matching turning points on the candidate indicator within lag months after the benchmark indicator
name_candi = 'Candidate Indicator';
name_bench = 'Benchmark Indicator';

%% Step 1: Initial identification of turning points
data_candi1 = smooth(data_candi,12); % Using a 12-item moving average to smooth the sequence
[IndPeaks_candi1,IndTroughs_candi1] = Check(data_candi1,N);

data_bench1 = smooth(data_bench,12); % Using a 12-item moving average to smooth the sequence
[IndPeaks_bench1,IndTroughs_bench1] = Check(data_bench1,N);

% Step 2: Turning point filtering
[IndPeaks_candi2,IndTroughs_candi2] = PointsClean(IndPeaks_candi1,IndTroughs_candi1,pt,pp,data_candi1);

[IndPeaks_bench2,IndTroughs_bench2] = PointsClean(IndPeaks_bench1,IndTroughs_bench1,pt,pp,data_bench1);

% Step 3: Turning point adjustment
% The adjustment method is to find the maximum (minimum) value in the original sequence within [-N_ad,N_ad] months of the smoothed sequence turning points
[IndPeaks_candi3_1, IndTroughs_candi3_1] = PointsAdjust(IndPeaks_candi2, IndTroughs_candi2, data_candi, N_ad);
[IndPeaks_candi3_2,IndTroughs_candi3_2] = PointsClean(IndPeaks_candi3_1,IndTroughs_candi3_1,pt,pp,data_candi);

[IndPeaks_bench3_1, IndTroughs_bench3_1] = PointsAdjust(IndPeaks_bench2, IndTroughs_bench2, data_bench, N_ad);
[IndPeaks_bench3_2,IndTroughs_bench3_2] = PointsClean(IndPeaks_bench3_1,IndTroughs_bench3_1,pt,pp,data_bench);

% Plotting
% Candidate Indicator
figure;
subplot(2,2,1);
MyPlot(IndPeaks_candi1, IndTroughs_candi1, data_candi1,[255 0 0]/255,[255 128 128]/255,[255 178 178]/255);
title([name_candi ' Smoothed Sequence Potential Turning Points Identification']);
legend({[name_candi ' Smoothed Sequence'],'Maximum Points','Minimum Points'},'Orientation','horizontal','Location','North');
subplot(2,2,2);
MyPlot(IndPeaks_candi2, IndTroughs_candi2, data_candi1,[255 0 0]/255,[255 128 128]/255,[255 178 178]/255);
title([name_candi ' Smoothed Sequence Turning Points Filtering']);
legend({[name_candi ' Smoothed Sequence'],'Maximum Points','Minimum Points'},'Orientation','horizontal','Location','North');
subplot(2,2,3);
MyPlot(IndPeaks_candi3_2, IndTroughs_candi3_2, data_candi,[255 0 0]/255,[255 128 128]/255,[255 178 178]/255);
title([name_candi ' Turning Point Adjustment']);
legend({[name_candi ' Original Sequence'],'Adjusted Maximum Points','Adjusted Minimum Points'},'Orientation','horizontal','Location','North');
subplot(2,2,4);
[h1, h2, h3] = MyPlot(IndPeaks_candi2, IndTroughs_candi2, data_candi1,[255 0 0]/255,[255 128 128]/255,[255 178 178]/255); % Smoothed sequence is red, dark red, deep green
[h4, h5, h6] = MyPlot(IndPeaks_candi3_2, IndTroughs_candi3_2, data_candi,[4 78 126]/255,[129 166 190]/255,[180 202 216]/255); % Adjusted to original sequence is blue, dark blue, light green
title([name_candi ' Smoothed Sequence & Original Sequence']);
legend([h1,h2,h3],[name_candi ' Smoothed Sequence'],'Maximum Points','Minimum Points','Orientation','horizontal','Location','North');
ah = axes('position',get(gca,'position'),'visible','off');
legend(ah,[h4,h5,h6],[name_candi ' Original Sequence'],'Adjusted Maximum Points','Adjusted Minimum Points','Orientation','horizontal','Location','North');

% Benchmark Indicator
figure;
subplot(2,2,1);
MyPlot(IndPeaks_bench1, IndTroughs_bench1,data_bench1,[255 0 0]/255,[255 128 128]/255,[255 178 178]/255);
title([name_bench ' Smoothed Sequence Potential Turning Points Identification']);
legend({[name_bench ' Smoothed Sequence'],'Maximum Points','Minimum Points'},'Orientation','horizontal','Location','North');
subplot(2,2,2);
MyPlot(IndPeaks_bench2, IndTroughs_bench2,data_bench1,[255 0 0]/255,[255 128 128]/255,[255 178 178]/255);
title([name_bench ' Smoothed Sequence Turning Points Filtering']);
legend({[name_bench ' Smoothed Sequence'],'Maximum Points','Minimum Points'},'Orientation','horizontal','Location','North');
subplot(2,2,3);
MyPlot(IndPeaks_bench3_2, IndTroughs_bench3_2, data_bench,[255 0 0]/255,[255 128 128]/255,[255 178 178]/255);
title([name_bench ' Benchmark Indicator Turning Point Adjustment']);
legend({[name_bench ' Original Sequence'],'Adjusted Maximum Points','Adjusted Minimum Points'},'Orientation','horizontal','Location','North');
subplot(2,2,4);
[h1,h2,h3] = MyPlot(IndPeaks_bench2, IndTroughs_bench2,data_bench1,[255 0 0]/255,[255 128 128]/255,[255 178 178]/255); % Smoothed sequence is red, dark red, deep green
[h4,h5,h6] = MyPlot(IndPeaks_bench3_2, IndTroughs_bench3_2, data_bench,[4 78 126]/255,[129 166 190]/255,[180 202 216]/255); % Adjusted to original sequence is blue, dark blue, light green
title([name_bench ' Smoothed Sequence & Original Sequence']);
legend([h1,h2,h3],[name_bench ' Smoothed Sequence'],'Maximum Points','Minimum Points','Orientation','horizontal','Location','North');
ah = axes('position',get(gca,'position'),'visible','off');
legend(ah,[h4,h5,h6],[name_bench ' Original Sequence'],'Adjusted Maximum Points','Adjusted Minimum Points','Orientation','horizontal','Location','North');

% Candidate Indicator & Benchmark Indicator
figure;
[h1,h2,h3] = MyPlot(IndPeaks_bench3_2, IndTroughs_bench3_2, data_bench,[255 0 0]/255,[255 128 128]/255,[129 166 190]/255); % Benchmark indicator is blue, dark blue, light green
[h4, ~, ~] = MyPlot(IndPeaks_candi3_2, IndTroughs_candi3_2, data_candi,[4 78 126]/255,[255 128 128]/255,[129 166 190]/255); % Candidate indicator is red, dark red, deep green
title([name_candi ' and ' name_bench ' Turning Point Identification and Matching']);
legend([h1,h4,h2,h3],name_bench,name_candi,'Maximum Points','Minimum Points','Orientation','horizontal','Location','North');

% Step 4: Turning point matching
[PiPei, QueShi, DuoYu, WuShu, LeadLag] = PointsMatch(IndPeaks_candi3_2, IndTroughs_candi3_2, data_candi, IndPeaks_bench3_2,IndTroughs_bench3_2 , lead, lag);
ratio_pipei = length(PiPei)/(length(PiPei)+length(QueShi)+length(DuoYu));
ratio_duoyu = length(DuoYu)/(length(PiPei)+length(DuoYu));
inflectionpoint=length(IndPeaks_bench3_2)+length(IndTroughs_bench3_2);
fprintf('Matched turning points = %d\n',length(PiPei));
%fprintf('Missing turning points = %d\n',length(QueShi));
fprintf('Missing turning points = %d\n',(inflectionpoint-length(PiPei)));
fprintf('Excess turning points = %d\n',length(DuoYu));
fprintf('No data turning points = %d\n',length(WuShu));
fprintf('Turning point matching rate = %.2f\n',ratio_pipei);
fprintf('Excess rate = %.2f\n',ratio_duoyu);
fprintf('Average lead order = %.2f\n',mean(LeadLag));
fprintf('Standard deviation of lead = %.2f\n',std(LeadLag));

% Return parameter organization
% Turning points on the original sequence of candidate indicators
mark_candi = zeros(size(data_candi,1),1);
mark_candi(IndPeaks_candi3_2) = 1;
mark_candi(IndTroughs_candi3_2) = -1;
% Turning points on the original sequence of benchmark indicators
mark_bench = zeros(size(data_bench,1),1);
mark_bench(IndPeaks_bench3_2) = 1;
mark_bench(IndTroughs_bench3_2) = -1;
% Turning points on the smoothed sequence of candidate indicators
mark_candi1 = zeros(size(data_candi1,1),1);
mark_candi1(IndPeaks_candi2) = 1;
mark_candi1(IndTroughs_candi2) = -1;
% Turning points on the smoothed sequence of benchmark indicators
mark_bench1 = zeros(size(data_bench1,1),1);
mark_bench1(IndPeaks_bench2) = 1;
mark_bench1(IndTroughs_bench2) = -1;
%C = [length(PiPei),length(QueShi),length(DuoYu),length(WuShu),ratio_pipei,ratio_duoyu,mean(LeadLag),std(LeadLag)];
C = [length(PiPei),(inflectionpoint-length(PiPei)),length(DuoYu),length(WuShu),ratio_pipei,ratio_duoyu,mean(LeadLag),std(LeadLag)];
D = [data_candi,data_bench,data_candi1,data_bench1,mark_candi,mark_bench,mark_candi1,mark_bench1];
end

%% Identify potential turning points
% Find peaks and troughs
% Input data is a column vector (M*1), and N (local turning points in the previous and following N months)
% Returns the processed peak array IndPeaks and trough array IndTroughs
function [IndPeaks,IndTroughs] = Check(data, N)
dataSize = size(data,1);
% Versions above r2017b can use
% IndPeaks = islocalmax(data, 'MinSeparation', 6);
% IndTroughs = islocalmin(data, 'MinSeparation', 6);

% Since my version does not support it, I wrote the function myself
% Find local turning points one month before and after
if N == 1
    [~, IndPeaks] = findpeaks(data);
    [~, IndTroughs] = findpeaks(-data);
    % Find local turning points N months before and after (N>1)
elseif N > 1
    % First select local turning points one month before and after (to speed up)
    [~, IndPeaks0] = findpeaks(data);
    [~, IndTroughs0] = findpeaks(-data);
    m = size(IndPeaks0,1);
    n = size(IndTroughs0,1);
    % Then filter based on N months before and after
    for i = 1:m
        for j = -N:N
            if IndPeaks0(i)+j > 0 && IndPeaks0(i)+j < dataSize+1 && data(IndPeaks0(i)+j) > data(IndPeaks0(i))
                IndPeaks0(i) = 0;
                break;
            end
        end
    end
    for i = 1:n
        for j = -N:N
            if IndTroughs0(i)+j > 0 &&  IndTroughs0(i)+j < dataSize+1 && data( IndTroughs0(i)+j) < data( IndTroughs0(i))
                IndTroughs0(i) = 0;
                break;
            end
        end
    end
    IndPeaks = Delete(IndPeaks0);
    IndTroughs = Delete(IndTroughs0);
else
    error('An error occurred, returning!');
end
end

%% Turning point standardization and filtering
% Input IndPeaks and IndTroughs are row vectors (N*1), data is a column vector (N*1)
% Returns the processed peak array IndPeaks2 and trough array IndTroughs2, same format
function [IndPeaks2,IndTroughs2] = PointsClean(IndPeaks1,IndTroughs1,pt,pp,data)
IndPeaks2_1 = CheckCirclePeaks(IndPeaks1, pp, data); % Cycle test
IndTroughs2_1 = CheckCircleTroughs(IndTroughs1, pp, data); % Cycle test
[IndPeaks2_2, IndTroughs2_2] = CheckCircleHalf(IndPeaks2_1, IndTroughs2_1, pt); % Half-cycle test
[IndPeaks2, IndTroughs2] = CheckAlternate(IndPeaks2_2, IndTroughs2_2, data); % Alternation test
end

% Exclude turning points within 6 months of the beginning and end of the data, and restrict the interval between adjacent peaks (cycle) to at least pp months
% Input IndPeaks1, pp, data (pp generally takes 16 months)
% Returns the processed peak array IndPeaks2
% The method of deleting turning points is to delete the smaller peak between adjacent peaks
function IndPeaks2 = CheckCirclePeaks(IndPeaks1, pp, data)
m = size(IndPeaks1,1);
prePeak = 0;
IndPeaks2 = nan(m,1);
for i = 1:m
    % Exclude turning points within 6 months of the beginning and end of the data
    if prePeak == 0 && IndPeaks1(i)  < 6 || size(data,1) - IndPeaks1(i) < 6
        continue;
    end
    % Restrict the interval between adjacent peaks (cycle) to at least pp months
    if prePeak ~= 0 && IndPeaks1(i) - prePeak < pp
        % If the peak value of the latter is less than or equal to the former, exclude the latter
        if data(IndPeaks1(i)) <= data(prePeak)
            continue;
            % If the peak value of the latter is greater than the former, exclude the former
        else
            IndPeaks2(length(IndPeaks2(~isnan(IndPeaks2)))) = nan; % Remove the last turning point in the peak array
        end
    end
    prePeak = IndPeaks1(i);
    IndPeaks2(length(IndPeaks2(~isnan(IndPeaks2)))+1) = prePeak; % Add the peak turning point to the peak array
end
IndPeaks2(isnan(IndPeaks2)) = []; % Remove nan values
end

% Exclude turning points within 6 months of the beginning and end of the data, and restrict the interval between adjacent valleys (cycle) to at least pp months
% Input IndTroughs1, pp, data (pp generally takes 16 months)
% Returns the processed trough array IndTroughs2
% The method of deleting turning points is to delete the larger trough between adjacent troughs
function IndTroughs2 = CheckCircleTroughs(IndTroughs1, pp, data)
m = size(IndTroughs1,1);
preTrough = 0;
IndTroughs2 = nan(m,1);
for i = 1:m
    % Exclude turning points within 6 months of the beginning and end of the data
    if preTrough == 0 && IndTroughs1(i) < 6 || size(data,1) - IndTroughs1(i) < 6
        continue;
    end
    % Restrict the interval between adjacent valleys (cycle) to at least pp months
    if preTrough ~= 0 && IndTroughs1(i) - preTrough < pp
        % If the trough value of the latter is less than or equal to the former, exclude the latter
        if data(IndTroughs1(i)) >= data(preTrough)
            continue;
            % If the trough value of the latter is greater than the former, exclude the former
        else
            IndTroughs2(length(IndTroughs2(~isnan(IndTroughs2)))) = nan; % Remove the last turning point in the trough array
        end
    end
    preTrough = IndTroughs1(i);
    IndTroughs2(length(IndTroughs2(~isnan(IndTroughs2)))+1) = preTrough; % Add the trough turning point to the trough array
end
IndTroughs2(isnan(IndTroughs2)) = []; % Remove nan values
end

% Restrict the interval between adjacent peaks and valleys (half-cycle) to at least pt months
% Input IndPeaks1, IndTroughs1, pt (pt generally takes 6 months)
% Returns the processed peak array IndPeaks2, trough array IndTroughs2
% The method of deleting turning points is to delete the earlier turning point between peaks and troughs
function [IndPeaks2, IndTroughs2] = CheckCircleHalf(IndPeaks1, IndTroughs1, pt)
m = size(IndPeaks1,1);
n = size(IndTroughs1,1);
for i = 1:m
    for j = 1:n
        if IndTroughs1(j) < IndPeaks1(i) && abs(IndTroughs1(j)-IndPeaks1(i)) < pt
            IndTroughs1(j) = 0;
        end
        if IndTroughs1(j) > IndPeaks1(i) && abs(IndTroughs1(j)-IndPeaks1(i)) < pt
            IndPeaks1(i) = 0;
        end
    end
end
% Remove non-conforming turning points
IndPeaks2 = Delete(IndPeaks1);
IndTroughs2 = Delete(IndTroughs1);
end

% Process peaks and troughs alternation
% Input IndPeaks1, IndTroughs1, data
% Returns the processed peak array IndPeaks2, trough array IndTroughs2
% The method of deleting turning points is to delete the smaller (larger) turning point between adjacent peaks (troughs)
function [IndPeaks2, IndTroughs2] = CheckAlternate(IndPeaks1, IndTroughs1, data)
dataSize = size(data,1);
mark = zeros(dataSize,1);
mark(IndPeaks1) = 1;
mark(IndTroughs1) = -1;
pre = 0; % Mark the position of the previous turning point
for i = 1:dataSize
    if mark(i) == 0
        continue;
    else
        if pre == 0  % If there was no turning point before
            pre = i;
        else
            if mark(i)+mark(pre) == 0 % If the current turning point alternates with the previous one
                pre = i;
                % If the current turning point does not alternate with the previous one, and they are consecutive peaks, the smaller peak must be deleted
            elseif mark(i) == 1
                if data(i) >= data(pre)
                    mark(pre) = 0;
                    pre = i;
                else
                    mark(i) = 0;
                end
                % If the current turning point does not alternate with the previous one, and they are consecutive troughs, the larger trough must be deleted
            else
                if data(i) <= data(pre)
                    mark(pre) = 0;
                    pre = i;
                else
                    mark(i) = 0;
                end
            end
        end
    end
end
IndPeaks2 = find(mark==1);
IndTroughs2 = find(mark==-1);
end

% Remove non-conforming turning points
% Input IndPoints1, non-conforming turning points are marked as 0
% Returns the processed turning point array IndPoints2
function IndPoints2 = Delete(IndPoints1)
IndPoints2_0 = sort(IndPoints1,1);
i = 1;
while i< length(IndPoints2_0)+1 && IndPoints2_0(i) == 0
    i = i+1;
end
if i < length(IndPoints2_0)+1
    IndPoints2 = IndPoints2_0(i:end);
else
    IndPoints2 = [];
end
end

%% Turning point adjustment
% Input IndPeaks1, IndTroughs1, data, N_ad (data is the sequence to be adjusted)
% Returns the processed peak array IndPeaks2, trough array IndTroughs2
function [IndPeaks2, IndTroughs2] = PointsAdjust(IndPeaks1, IndTroughs1, data, N_ad)
dataSize = size(data,1);
m = size(IndPeaks1,1);
n = size(IndTroughs1,1);
IndPeaks2 = IndPeaks1;
IndTroughs2 = IndTroughs1;
for i = 1:m
    max = intmax;
    for j = -N_ad:N_ad
        if IndPeaks1(i)+j > 0 && IndPeaks1(i)+j < dataSize+1
            if max == intmax
                max = data(IndPeaks1(i)+j);
                IndPeaks2(i) = IndPeaks1(i)+j;
            else
                if data(IndPeaks1(i)+j) > max
                    max = data(IndPeaks1(i)+j);
                    IndPeaks2(i) = IndPeaks1(i)+j;
                end
            end
        end
    end
end
for i = 1:n
    min = intmin;
    for j = -N_ad:N_ad
        if IndTroughs1(i)+j > 0 && IndTroughs1(i)+j < dataSize+1
            if min == intmin
                min = data(IndTroughs1(i)+j);
                IndTroughs2(i) = IndTroughs1(i)+j;
            else
                if data(IndTroughs1(i)+j) < min
                    min = data(IndTroughs1(i)+j);
                    IndTroughs2(i) = IndTroughs1(i)+j;
                end
            end
        end
    end
end
end

%% Turning point matching
% Input IndPeaks_candi, IndTroughs_candi, data_candi, IndPeaks_bench, IndTroughs_bench, lead, lag
% Output PiPei, QueShi, DuoYu, WuShu, records the positions of each category in the respective sequence, all in N*1 format
function [PiPei, QueShi, DuoYu, WuShu,LeadLag] = PointsMatch(IndPeaks_candi, IndTroughs_candi, data_candi, IndPeaks_bench, IndTroughs_bench, lead, lag)
m_candi = size(IndPeaks_candi,1);
n_candi = size(IndTroughs_candi,1);
m_bench = size(IndPeaks_bench,1);
n_bench = size(IndTroughs_bench,1);
PiPei0 = nan(m_bench+n_bench,1);
PiPei_candi = nan(m_candi+n_candi,1);
QueShi0 = nan(m_bench+n_bench,1);
DuoYu0 = nan(m_candi+n_candi,1);
WuShu0 = nan(m_bench+n_bench,1);
LeadLag0 = nan(m_bench+n_bench,1); % If an element is 1, it means that the point on the candidate indicator leads the point on the benchmark indicator by 1 period
% Marking of benchmark indicator turning points: PiPei0 (matching), QueShi0 (missing), WuShu0 (no data)
% Marking of candidate indicator turning points: PiPei_candi (matching), DuoYu (excess)
for i = 1:m_bench
    isNotNull = 0;
    for j = -lead:lag
        if IndPeaks_bench(i)+j > 0 && IndPeaks_bench(i)+j < size(data_candi,1)+1 && isnan(data_candi(IndPeaks_bench(i)+j))
        else
            isNotNull = 1;
            IdElem = Search(IndPeaks_candi, IndPeaks_bench(i)+j);
            if IdElem ~= 0
                PiPei0(length(PiPei0(~isnan(PiPei0)))+1) = IndPeaks_bench(i);
                PiPei_candi(length(PiPei_candi(~isnan(PiPei_candi)))+1) = IndPeaks_bench(i)+j;
                LeadLag0(length(LeadLag0(~isnan(LeadLag0)))+1) = -j;
                break;
            end
        end
    end
    if isNotNull == 0
        WuShu0(length(WuShu0(~isnan(WuShu0)))+1) = IndPeaks_bench(i);
    else
        if j == lag
            QueShi0(length(QueShi0(~isnan(QueShi0)))+1) = IndPeaks_bench(i);
        end
    end
end
for i = 1:n_bench
    isNotNull = 0;
    for j = -lead:lag
        if IndTroughs_bench(i)+j > 0 && IndTroughs_bench(i)+j < size(data_candi,1)+1 && isnan(data_candi(IndTroughs_bench(i)+j))
        else
            isNotNull = 1;
            IdElem = Search(IndTroughs_candi, IndTroughs_bench(i)+j);
            if IdElem ~= 0
                PiPei0(length(PiPei0(~isnan(PiPei0)))+1) = IndTroughs_bench(i);
                PiPei_candi(length(PiPei_candi(~isnan(PiPei_candi)))+1) = IndTroughs_bench(i)+j;
                LeadLag0(length(LeadLag0(~isnan(LeadLag0)))+1) = -j;
                break;
            end
        end
    end
    if isNotNull == 0
        WuShu0(length(WuShu0(~isnan(WuShu0)))+1) = IndTroughs_bench(i);
    else
        if j == lag
            QueShi0(length(QueShi0(~isnan(QueShi0)))+1) = IndTroughs_bench(i);
        end
    end
end

% Marking of candidate indicator: DuoYu0 (excess)
for i = 1:m_candi
    IdElem = Search(PiPei_candi, IndPeaks_candi(i));
    if IdElem == 0
        DuoYu0(length(DuoYu0(~isnan(DuoYu0)))+1) = IndPeaks_candi(i);
    end
end
for i = 1:n_candi
    IdElem = Search(PiPei_candi, IndTroughs_candi(i));
    if IdElem == 0
        DuoYu0(length(DuoYu0(~isnan(DuoYu0)))+1) = IndTroughs_candi(i);
    end
end
% Remove nan values
PiPei0(isnan(PiPei0)) = [];
% PiPei_candi(isnan(PiPei_candi)) = [];
QueShi0(isnan(QueShi0)) = [];
DuoYu0(isnan(DuoYu0)) = [];
WuShu0(isnan(WuShu0)) = [];
LeadLag0(isnan(LeadLag0)) = [];

% Sort the marked turning point sequences
[PiPei,b] = sort(PiPei0,1);
QueShi = sort(QueShi0,1);
DuoYu = sort(DuoYu0,1);
WuShu = sort(WuShu0,1);
LeadLag = LeadLag0(b,:);
end

% Search for an element elem in the sequence arr
% Input arr,elem
% Output IdElem, the position of elem in arr; if not found, then IdElem = 0
function IdElem  = Search(arr,elem)
IdElem = 0;
for i = 1:length(arr)
    if arr(i) == elem
        IdElem = i;
    end
end
end

% Plotting
function [h1, h2, h3] = MyPlot(IndPeaks, IndTroughs, data , lineColor, peaksColor, troughsColor)
x=1:size(data,1);
hold on
h1 = plot(x,data,'-','Color',lineColor);
h2 = plot(x(IndPeaks),data(IndPeaks),'.','Color',peaksColor,'MarkerSize',15);
h3 = plot(x(IndTroughs),data(IndTroughs),'.','Color',troughsColor,'MarkerSize',15);
hold off
end
