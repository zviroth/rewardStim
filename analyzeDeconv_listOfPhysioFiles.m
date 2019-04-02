dataFolder = '~/data/rwdFmri/';
% subFolders = {'001920190214/'};
% roiNames = {'lh_V1_exvivo','rh_V1_exvivo'};
samplerate=500;
subFolders = {'007420190221/'};
subFolders = {'007520190318/'};
subFolders = {'008320190320/'};
subFolders = {'009020190325/'};
% subFolders = {'007420190221_OC/'};
% physioNames = {[
% roiNames = {'left','right'};
% roiNames = {'l_v1','r_v1'};
% roiNames = {'l_v1_ecc4','r_v1_ecc4'};
% roiNames = {'l_v1','r_v1','l_v1_ecc4','r_v1_ecc4'};
roiNames = {'lV1_eccen4','rV1_eccen4'};
% roiNames = {'lV1_eccen0to1','rV1_eccen0to1','lV1_eccen1to4','rV1_eccen1to4','lV1_eccen4to10','rV1_eccen4to10','lV1_eccen10','rV1_eccen10'};
frames=170;
junkedFrames = 10;
TR=1.5;

physioSampleRate = 50; %Hz
numSubs = length(subFolders);
numRois = length(roiNames);
curFolder = pwd;
onlyCorrect=0;
cd(dataFolder);
mrQuit;
if exist('v')
    deleteView(v);
end
rows = numRois;
cols=numSubs;
trialsPerRun=16;
trialLength = 10;
nVolumes = trialsPerRun*trialLength;
numContrasts = 2;
numFreqs = 5;

plotColors = {[0 0 1], [1 0 0], [0 1 0], [0.5 1 0.2]};
plotStyles = {'-','--',':','-.','-','--',':','-.'};
clear d concatInfo sumResponse subResponse roiMeanTseries meanResponse roiTC subRunResponse trialCorrectness trialResponse trialRT propCorrect
clear contrast spatFreq nullTrial contrastConcat spatFreqConcat nullTrialConcat designMatContrast designMatFreq nullDeconv deconvMatNull meanDeconvNull
clear taskResp taskPrediction meanPupil
clear ecg ecgPeaks ecgPeaksDiff ecgPeaksAmp resp respPeaks respPeaksDiff respPeaksAmp
clear designMatContrastFreq designMatNull taskBaseline
% figure(1)
% clf;

%set canonical HRF
params = struct;
sampleDuration = 1.5;
sampleDelay=sampleDuration/2;
defaultParams=1;
[params hrf] = hrfDoubleGamma(params,sampleDuration,sampleDelay,defaultParams);
deconvLength = 10;



for iSub = 1:numSubs
    cd(subFolders{iSub});
    v=newView;
    % switch to the concatenation group
    v = viewSet(v, 'curGroup', 'Concatenation');
    nScans = viewGet(v, 'nscans');
    
    for iScan = 1:nScans%2 concatenations, 1 for each reward type
        d{iScan} = viewGet(v, 'd', iScan);
        %             concatInfo{iSub,iScan} = viewGet(v, 'concatInfo', iScan);
        s = viewGet(v, 'stimfile', iScan);
        
        rwdType = s{1}.myscreen.stimulus.rewardType;
        if strcmp(rwdType, 'H')
            rwdTypeNum = 1;
        elseif strcmp(rwdType, 'L')
            rwdTypeNum = 2;
        else
            disp('wtf');
            keyboard
        end
        concatRwdTypeNum(iScan) = rwdTypeNum;
        sacRwd{rwdTypeNum}=[];
        for r=1:length(s)
            trialCorrectness{iSub,rwdTypeNum}(r,:) = s{r}.task{1}{1}.correctness;
            trialResponse{iSub,rwdTypeNum,r} = s{r}.task{1}{1}.response;
            propCorrect{iSub,rwdTypeNum}(r) = sum(trialCorrectness{iSub,rwdTypeNum}(r,:)) / length(trialCorrectness{iSub,rwdTypeNum}(r,:));
            stimfile = s{r}.filename;
            e{r} = myGetTaskEyeTraces(stimfile, 'removeBlink=3');
            runEyeSize{rwdTypeNum}(r,:) = size(e{r}.eye.pupil);
            sacRun{r} = [];
            for itrial=1:e{r}.nTrials
                eyePos = [e{r}.eye.xPos(itrial,:)'  e{r}.eye.yPos(itrial,:)'];
               eyevel{r,itrial} = vecvel(eyePos,samplerate, 2);
               sac{r,itrial} = microsacc(eyePos, eyevel{r,itrial}, 1, 5);
               sacRun{r} = [sacRun{r}; sac{r,itrial}];
            end
            sacRwd{rwdTypeNum} = [sacRwd{rwdTypeNum}; sacRun{r}];

       
        end
        [sacRwd{rwdTypeNum}(:,1), ind] = sort(sacRwd{rwdTypeNum}(:,1));%saccade onset
        sacRwd{rwdTypeNum}(:,2:7) = sacRwd{rwdTypeNum}(ind,2:7);%saccade end, velocity, horiz/vert component, horiz/vert amplitude
        [binnedSac{rwdTypeNum},edges] = histcounts(sacRwd{rwdTypeNum}(:,1),max(sacRwd{rwdTypeNum}(:,1)));

        trialLengthEye(rwdTypeNum) = max(runEyeSize{rwdTypeNum}(:,2));
        numTrials(rwdTypeNum) = sum(runEyeSize{rwdTypeNum}(:,1));
        rwdPupil{rwdTypeNum} = NaN(numTrials(rwdTypeNum), trialLengthEye(rwdTypeNum));
        trialCounter=0;
        for r=1:length(s)
            rwdPupil{rwdTypeNum}(trialCounter+1:trialCounter+runEyeSize{rwdTypeNum}(r,1), 1:runEyeSize{rwdTypeNum}(r,2)) = e{r}.eye.pupil;
            trialCounter = trialCounter + runEyeSize{rwdTypeNum}(r,1);
        end
        %plot mean pupil size
        meanPupil{rwdTypeNum} = nanmean(rwdPupil{rwdTypeNum})';
        
        
        
        %extract conditions for all scans within concatenation
        % doing this manually because the volNum is messed up, missed
        % triggers, instead of using getTaskParameters
        % names_: {'contrast'  'spatFreq'  'nullTrial'}
        contrastConcat{iScan} = [];
        spatFreqConcat{iScan} = [];
        nullTrialConcat{iScan} = [];
        for r=1:length(s)%per run
            
            contrastConcat{iScan} = [contrastConcat{iScan} s{r}.task{2}{1}.randVars.contrast(2:trialsPerRun+1)];
            spatFreqConcat{iScan} = [spatFreqConcat{iScan} s{r}.task{2}{1}.randVars.spatFreq(2:trialsPerRun+1)];
            nullTrialConcat{iScan} = [nullTrialConcat{iScan} s{r}.task{2}{1}.randVars.nullTrial(2:trialsPerRun+1)];
        end
        
        contrastConcat{iScan} = contrastConcat{iScan}.*(1-nullTrialConcat{iScan});
        spatFreqConcat{iScan} = spatFreqConcat{iScan}.*(1-nullTrialConcat{iScan});
        
        totalTrials = length(contrastConcat{iScan});
        %contrast
        smallContrastMat = zeros(numContrasts,totalTrials);%number per trial
        for c=1:numContrasts
            smallContrastMat(c,contrastConcat{iScan}==c) = ones;
        end
        contrastMat = upsample(smallContrastMat,trialLength);%number per volume
        for c=1:numContrasts
             temp = conv(contrastMat(c,:),hrf);
             designMatContrast{iScan}(c,:) =  temp(1:totalTrials*trialLength);
        end
        designMatContrast{iScan}(numContrasts+1,:) = ones;
        %spatial frequency
        smallFreqMat = zeros(numFreqs,totalTrials);%number per trial
        for c=1:numFreqs
            smallFreqMat(c,spatFreqConcat{iScan}==c) = ones;
        end
        sfMat = upsample(smallFreqMat,trialLength);%number per volume
        for c=1:numFreqs
             temp = conv(sfMat(c,:),hrf);
             designMatFreq{iScan}(c,:) =  temp(1:totalTrials*trialLength);
        end        
        designMatFreq{iScan}(numFreqs+1,:) = ones;
        
        %spatial frequency & contrast!
        smallContrastSfMat = zeros(numFreqs*numContrasts,totalTrials);%number per trial
        for c=1:numFreqs*numContrasts
            curCond = (contrastConcat{iScan}-1)*2+spatFreqConcat{iScan};
            smallContrastSfMat(c,curCond==c) = ones;
        end
        contrastSfMat = upsample(smallContrastSfMat,trialLength);%number per volume
        for c=1:numFreqs*numContrasts
             temp = conv(contrastSfMat(c,:),hrf);
             designMatContrastFreq{iScan}(c,:) =  temp(1:totalTrials*trialLength);
        end        
        designMatContrastFreq{iScan}(numFreqs*numContrasts+1,:) = ones;
        
        %null vs stimulated trials
        smallNullMat = zeros(2,totalTrials);%number per trial
        for c=0:1
            smallNullMat(c+1,nullTrialConcat{iScan}==c) = ones; %1=stim, 2=null
        end
        nullTrialMat = upsample(smallNullMat,trialLength);%number per volume
        for c=1:2
             temp = conv(nullTrialMat(c,:),hrf);
             designMatNull{iScan}(c,:) =  temp(1:totalTrials*trialLength);
        end        
        designMatNull{iScan}(2+1,:) = ones;
        
        
        
        %nullTrial deconvolution matrix
        deconvMatNull{iScan} = zeros((size(nullTrialMat,1)-1)*deconvLength,size(nullTrialMat,2));
        for c=1:2
            for j=1:deconvLength
                deconvMatNull{iScan}((c-1)*deconvLength+j,:) = [zeros(1,j-1) nullTrialMat(c,1:end-j+1)];%shift and pad with zeros
            end
        end
        %do we need the constant predictor? Or does it mess up the
        %subtraction?
        deconvMatNull{iScan}(end+1,:) = ones;
        
        
        
        
        %create  vector with a 1 for each trial, to be convolved with the
        %task response (null trials)
        taskSmallMat = ones(1,totalTrials);
        taskMat = upsample(taskSmallMat,trialLength);%add zeros in between

        for iRoi=1:numRois
            roiTC{iRoi,rwdTypeNum} = loadROITSeries(v, roiNames{iRoi}, iScan, [], 'keepNAN',true);
            concatInfo{iSub,rwdTypeNum} = viewGet(v, 'concatInfo', iScan);
            %average trial
            trialTC{iRoi}(rwdTypeNum,:) = nanmean(reshape(nanmean(roiTC{iRoi,rwdTypeNum}.tSeries), trialLength, totalTrials),2);
            %contrast betas
            contrastBetasWithTask{iRoi,rwdTypeNum} = designMatContrast{iScan}'\roiTC{iRoi,rwdTypeNum}.tSeries';
%             residuals = roiTC{iRoi,rwdTypeNum} - designMatContrast{iScan}' * contrastBetas{iRoi,rwdTypeNum};
            %sf betas
            sfBetasWithTask{iRoi,rwdTypeNum} = designMatFreq{iScan}'\roiTC{iRoi,rwdTypeNum}.tSeries';
            sfContrastBetasWithTask{iRoi,rwdTypeNum}= designMatContrastFreq{iScan}'\roiTC{iRoi,rwdTypeNum}.tSeries';
            
            %deconvolution with stim/null predictors
            nullDeconv{iRoi,rwdTypeNum} = deconvMatNull{iScan}'\roiTC{iRoi,rwdTypeNum}.tSeries';
            for c=1:2
                for t=1:deconvLength
                    meanDeconvNull(iRoi,rwdTypeNum,c,t) = nanmean(nullDeconv{iRoi,rwdTypeNum}((c-1)*deconvLength+t,:));%mean over voxels
                end
            end
            
%             %get the null response
            taskResp{iRoi}(rwdTypeNum,:,:) = nullDeconv{iRoi,rwdTypeNum}(deconvLength+1:2*deconvLength,:); %for each voxel
            taskBaseline{iRoi}(rwdTypeNum,:) = nullDeconv{iRoi,rwdTypeNum}(2*deconvLength+1,:); %for each voxel
            %get the stim response
%             taskResp{iRoi}(rwdTypeNum,:,:) = nullDeconv{iRoi,rwdTypeNum}(1:deconvLength,:); %for each voxel
            
%             meanTaskResp{iRoi}(rwdTypeNum,:) = meanDeconvNull(iRoi,rwdTypeNum,2,:); %mean over voxels
            %convole with task events
            nvox = roiTC{iRoi,rwdTypeNum}.n;
            for i=1:nvox
                temp = conv(taskMat, taskResp{iRoi}(rwdTypeNum,:,i));
                taskPrediction{iRoi,rwdTypeNum}(i,:) = temp(1:totalTrials*trialLength);
            end
            
            
            %subtract the task response
            tcMinusTask{iRoi,rwdTypeNum} = roiTC{iRoi,rwdTypeNum}.tSeries - taskPrediction{iRoi,rwdTypeNum};
%             tcMinusTask{iRoi,rwdTypeNum} = zscore(tcMinusTask{iRoi,rwdTypeNum},0,2);
            
            %contrast betas
            contrastBetas{iRoi,rwdTypeNum} = designMatContrast{iScan}'\tcMinusTask{iRoi,rwdTypeNum}';
            %sf betas
            sfBetas{iRoi,rwdTypeNum} = designMatFreq{iScan}'\tcMinusTask{iRoi,rwdTypeNum}';
            %betas for combination of contrast and sf
            sfContrastBetas{iRoi,rwdTypeNum}= designMatContrastFreq{iScan}'\tcMinusTask{iRoi,rwdTypeNum}';
        end
        
    end
    deleteView(v);
    
    %% load list of original files - for physio
    load('origFilenames.mat','listOfOriginalFilenames');
%     listOfOriginalFilenames(1) = [];%this is just the blip up file
    numRuns(iSub) = concatInfo{iSub,1}.n+concatInfo{iSub,1}.n;
    if numRuns(iSub)~=length(listOfOriginalFilenames)
       sprintf([num2str(numRuns(iSub)) ' runs, but ' num2str(length(listOfOriginalFilenames)) ' physio files!']); 
    end
%     listOfOriginalFilenames{1}(6:7)
%     keyboard
    %% reach 'realtime' folder for physio data
    origFolder = pwd;
    for j=1:2
        dirList = dir;
        nameLength=0;
        for i=1:length(dirList)
            if dirList(i).isdir & length(dirList(i).name) > nameLength
                nameLength = length(dirList(i).name);
                mrnFolder = dirList(i).name;
            end
        end
        cd(mrnFolder)
    end
    %% ECG
    selectECG=0.03;
    selectResp=0.1;
    display=0;
    clear d
    %the first file will be for the blip?
    for i=1:length(listOfOriginalFilenames)
        d(i)=dir(['realtime/ECG_epiRTnih_scan_00' listOfOriginalFilenames{i}(6:7) '*.*']);
        ecg{i} = load(['realtime/' d(i).name]); 
        [ecgPeaks{i},criterion] = pickpeaks(ecg{i},selectECG,display);
        ecgPeaksDiff{i} = diff(ecgPeaks{i});
        ecgPeaksAmp{i} = ecg{i}(ecgPeaks{i});  
        runMeanEcgPeakDiff(i) = mean(ecgPeaksDiff{i});
        runMeanEcgPeakAmp(i) = mean(ecgPeaksAmp{i});
        
        d(i)=dir(['realtime/Resp_epiRTnih_scan_00' listOfOriginalFilenames{i}(6:7) '*.*']);
        resp{i} = load(['realtime/' d(i).name]); 
        [respPeaks{i},criterion] = pickpeaks(resp{i},selectResp,display);
        respPeaksDiff{i} = diff(respPeaks{i});
        respPeaksAmp{i} = resp{i}(respPeaks{i});
        runMeanRespPeakDiff(i) = mean(respPeaksDiff{i});
        runMeanRespPeakAmp(i) = mean(respPeaksAmp{i});
    end
    %assign physio files to fMRI runs
    %first, get scan numbers for each reward type
    for iScan=1:nScans
        scanNums = [];
        for i=1:length(concatInfo{iScan}.filename)
            scanNums = [scanNums str2num(concatInfo{iScan}.filename{i}(9:10))];
        end
        concatRuns{concatRwdTypeNum(iScan)} = scanNums;
%         concatRwdTypeNum(iScan)
    end
    %then, grab physio data for each reward type

    physioSamples = physioSampleRate*(frames-junkedFrames)*TR;
    for rwd=1:2
        allEcg{rwd}=[];
        allResp{rwd}=[];
        for i=1:length(concatRuns{rwd})
            if length(ecg(concatRuns{rwd}(i)))~=physioSamples
               sprintf(['ecg for run ' num2str(concatRuns{rwd}(i)) ' has only ' num2str(length(ecg(concatRuns{rwd}(i)))) ' samples, instead of ' num2str(physioSamples)]);
            end
            if length(resp(concatRuns{rwd}(i)))~=physioSamples
               sprintf(['respiration for run ' num2str(concatRuns{rwd}(i)) ' has only ' num2str(length(resp(concatRuns{rwd}(i)))) ' samples, instead of ' num2str(physioSamples)]);
            end
            allEcg{rwd} = [allEcg{rwd}; ecg{concatRuns{rwd}(i)}];
            allResp{rwd} = [allResp{rwd}; resp{concatRuns{rwd}(i)}];
        end
        [concatEcgPeaks{rwd},criterion] = pickpeaks(allEcg{rwd},selectECG,display);
        [concatRespPeaks{rwd},criterion] = pickpeaks(allResp{rwd},selectResp,display);
        
        concatEcgPeaksDiff{rwd} = diff(concatEcgPeaks{rwd});
        concatEcgPeaksAmp{rwd} = allResp{rwd}(concatEcgPeaks{rwd});
        meanEcgPeakDiff(rwd) = mean(concatEcgPeaksDiff{rwd});
        meanEcgPeakAmp(rwd) = mean(concatEcgPeaksAmp{rwd});
        
        concatRespPeaksDiff{rwd} = diff(concatRespPeaks{rwd});
        concatRespPeaksAmp{rwd} = allResp{rwd}(concatRespPeaks{rwd});
        meanRespPeakDiff(rwd) = mean(concatRespPeaksDiff{rwd});
        meanRespPeakAmp(rwd) = mean(concatRespPeaksAmp{rwd});
        
    end
%     ecgMinBytes = 40000;
%     respMinBytes = 50000;
%     d=dir('realtime/ECG*.*');
%     select=0.03;
%     display=0;
%     goodEcgRuns = [];
%     for i=1:length(d)
%         sz = d(i).bytes;
%         if sz>ecgMinBytes
%             goodEcgRuns = [goodEcgRuns i];
%         end
%     end
%     for i=1:length(goodEcgRuns)
%         filename = d(goodEcgRuns(i)).name;
%         ecg{i}=load(['realtime/' filename]);
%         [ecgPeaks{i},criterion] = pickpeaks(ecg{i},select,display);
%         ecgPeaksDiff{i} = diff(ecgPeaks{i});
%         ecgPeaksAmp{i} = ecg{i}(ecgPeaks{i});        
%     end
    
%     %% RESP
%     d=dir('realtime/Resp*.*');
%     select=0.1;
%     display=0; 
%     goodRespRuns = [];
%     for i=1:length(d)
%         sz = d(i).bytes;
%         if sz>respMinBytes
%             goodRespRuns = [goodRespRuns i];
%         end
%     end
%     for i=1:length(goodRespRuns)
%         filename = d(goodRespRuns(i)).name;
%         resp{i}=load(['realtime/' filename]);
%         [respPeaks{i},criterion] = pickpeaks(resp{i},select,display);
%         respPeaksDiff{i} = diff(respPeaks{i});
%         respPeaksAmp{i} = resp{i}(respPeaks{i});
%     end
%     
%     cd('..'); %cd(origFolder);
    %% return to home directory
    cd('..');
end
%% FIGURES
for i=1:11
    figure(i)
    clf
end

figure(1)
for rwdTypeNum=1:2
    subplot(4,1,rwdTypeNum)
%     meanPupil{rwdTypeNum} = nanmean(rwdPupil{rwdTypeNum})';
    plot(meanPupil{rwdTypeNum}, 'linewidth', 1)
    hold all
    stdRwd = std(rwdPupil{rwdTypeNum},0,1,'omitnan');
    stdRwd(isnan(stdRwd)) = ones;
%     meanPupil{rwdTypeNum}(isnan(meanPupil{rwdTypeNum})) = zeros;
    dsErrorsurface(1:length(stdRwd), meanPupil{rwdTypeNum}, stdRwd, [0.8 0.8 0.8], 0.2);
    if rwdTypeNum==1
        title('high')
    else
        title('low')
    end
%     title(rwdType);
    cueTime = s{r}.fixStimulus.cueTime * samplerate;
    interTime = s{r}.fixStimulus.interTime * samplerate;
    responseTime = s{r}.fixStimulus.responseTime * samplerate;
    allTimes = [cueTime interTime cueTime responseTime];
    plotRwdStimSegments(meanPupil{rwdTypeNum}, allTimes)
end

figure(1)%pupil
subplot(4,1,3)
for rwdTypeNum=1:2
    plot(meanPupil{rwdTypeNum}, 'Color', plotColors{rwdTypeNum}, 'linewidth', 1)
    hold on
end
legend('high', 'low')
subplot(4,1,4)
l = min(length(meanPupil{1}), length(meanPupil{2}));
l=l-50;
diffPupil = meanPupil{1}(1:l) - meanPupil{2}(1:l);
plot(zeros(l,1),'--');
hold all
plot(diffPupil, 'linewidth', 1)
plotRwdStimSegments(diffPupil, allTimes)
     


figure(2)%microsaccades
subplot(2,1,1)
for rwdTypeNum=1:2
   plot(binnedSac{rwdTypeNum})
   hold on
end
subplot(2,1,2)
L = 60;
for rwdTypeNum=1:2
    smoothSac{rwdTypeNum} = conv(binnedSac{rwdTypeNum}, ones(L,1)/L, 'same');
   plot(smoothSac{rwdTypeNum})
   hold on
end
legend('high', 'low')
        
rows=2;
cols=numRois;
rwdTypes = {'high','low'};

for iRoi=1:numRois
    i=3;
    figure(i)
    i=i+1;
    subplot(1,numRois,iRoi)
    plot(trialTC{iRoi}', 'linewidth', 1);
    xlabel('time (TR)');
    title([roiNames{iRoi}])
    if iRoi==1; legend('high','low'); end
    for rwdTypeNum=1:2
        j=i;
        
        figure(j)
        j=j+1;
        subplot(rows,cols,iRoi + (rwdTypeNum-1)*cols)
        plot(squeeze(meanDeconvNull(iRoi, rwdTypeNum,:,:))', 'linewidth', 1);
        title([roiNames{iRoi} ': ' rwdTypes{rwdTypeNum}]);
        if iRoi==1 && rwdTypeNum==1;  legend('stim', 'null'); end
        
        figure(j)
        j=j+1;
        subplot(rows,cols,iRoi)
        plot(nanmean(contrastBetas{iRoi,rwdTypeNum}(1:end-1,:),2), 'linewidth', 1);
        hold all
        title([roiNames{iRoi}  ': task removed']);
        
        subplot(rows,cols,cols+iRoi)
        plot(nanmean(contrastBetasWithTask{iRoi,rwdTypeNum}(1:end-1,:),2), 'linewidth', 1);
        hold all
        title([roiNames{iRoi}]);
        if iRoi==1 && rwdTypeNum==2;legend('high', 'low'); end
        
        figure(j)
        j=j+1;
        subplot(rows,cols,iRoi)
        plot(nanmean(sfBetas{iRoi,rwdTypeNum}(1:end-1,:),2), 'linewidth', 1);
        hold all
        title([roiNames{iRoi} ': task removed']);
        subplot(rows,cols,cols+iRoi)
        plot(nanmean(sfBetasWithTask{iRoi,rwdTypeNum}(1:end-1,:),2), 'linewidth', 1);
        hold all
        title([roiNames{iRoi}]);
        if iRoi==1 && rwdTypeNum==2;legend('high', 'low'); end
        
    end
end

%% ECG figure
rows = ceil(length(ecg)/3);
cols = ceil(length(ecg)/rows);
figure(j)
j=j+1;
for i=1:length(ecg)
    subplot(rows,cols,i)
    plot(ecg{i})
    hold on
    plot(ecgPeaks{i},ecgPeaksAmp{i},'ro','MarkerSize',3,'LineWidth',2)
end
figure(j)
j=j+1;
for i=1:length(ecgPeaksDiff)
    subplot(rows,cols,i)
    plot(ecgPeaksDiff{i})
    hold all
    title(['ecg ' num2str(mean(ecgPeaksDiff{i}))])
end

%% Resp figure
figure(j)
j=j+1;
for i=1:length(ecg)
    subplot(rows,cols,i)
    plot(resp{i})
    hold on
    plot(respPeaks{i},respPeaksAmp{i},'ro','MarkerSize',3,'LineWidth',2)
end
figure(j)
j=j+1;
for i=1:length(respPeaksDiff)
    subplot(rows,cols,i)
    plot(respPeaksDiff{i})
    hold all
    title(['resp ' num2str(mean(respPeaksDiff{i}))])
end

%% timecourse before and after removing null trial response
figure(11)
clf
% subplot(1,2,1); 
rows=2;
cols = numRois;
for iRoi=1:numRois
    for rwdTypeNum=1:2
        subplot(rows,cols,iRoi+(rwdTypeNum-1)*cols)
        plot(nanmean(roiTC{iRoi,rwdTypeNum}.tSeries),'linewidth', 1)
        hold on
        plot(nanmean(tcMinusTask{iRoi,rwdTypeNum}),'linewidth', 1)
        % subplot(1,2,2);
%         plot(1+nanmean(taskPrediction{iRoi,rwdTypeNum}), 'linewidth', 1)
        plot(1+nanmean(taskPrediction{iRoi,rwdTypeNum}), 'linewidth', 1)
        title([roiNames{iRoi} ': ' rwdTypes{rwdTypeNum}]);
    end
end
legend('original', 'null removed', 'null model');
% subplot(1,3,3); plot(mean(tcMinusTask{iRoi,rwdTypeNum}))




%%
%%%%%%%%%%%%%%%%%%%%%%%%%
% increases the sample rate of x by inserting n - 1 zeros between samples. If x is a matrix, the function treats each column as a separate sequence.
% inserts additional n-1 zeros after the last sample.
%%%%%%%%%%%%%%%%%%%%%%%%%%
function bigMatrix = upsample(smallMatrix,volumesPerTrial)
numTrials = size(smallMatrix,2);
numConditions = size(smallMatrix,1);
bigMatrix = zeros(numConditions, volumesPerTrial*numTrials);
bigMatrix(:,1:volumesPerTrial:end) = smallMatrix;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = plotRwdStimSegments(meanPupil, allTimes)

itime = 1;
startT=1;
endT = startT+allTimes(itime);
plot(startT:endT-1,meanPupil(startT:endT-1),'lineWidth',4,'color',[0 1 1]);
itime = itime+1;
startT=endT;
endT = startT+allTimes(itime);
plot(startT:endT-1,meanPupil(startT:endT-1),'lineWidth',4,'color',[0.1 0.7 0.6]);
itime = itime+1;
startT=endT;
endT = startT+allTimes(itime);
plot(startT:endT-1,meanPupil(startT:endT-1),'lineWidth',4,'color',[0 1 1]);
itime = itime+1;
startT=endT;
endT = startT+allTimes(itime);
plot(startT:endT-1,meanPupil(startT:endT-1),'lineWidth',4,'color',[1 1 0.0]);

end