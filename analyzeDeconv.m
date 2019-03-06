dataFolder = '~/data/rwdFmri/';
% subFolders = {'001920190214/'};
% roiNames = {'lh_V1_exvivo','rh_V1_exvivo'};
samplerate=500;
subFolders = {'007420190221/'};
roiNames = {'left','right'};

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
%             [sacRun{r}(:,1), ind] = sort(sacRun{r}(:,1));%saccade onset
%             sacRun{r}(:,2:7) = sacRun{r}(ind,2:7);%saccade end, velocity, horiz/vert component, horiz/vert amplitude
%             [binnedSac{r},edges] = histcounts(sacRun{r}(:,1),max(sacRun{r}(:,1)));
       
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
            nullDeconv{iRoi,rwdTypeNum} = deconvMatNull{iScan}'\roiTC{iRoi,rwdTypeNum}.tSeries';
            for c=1:2
                for t=1:deconvLength
                    meanDeconvNull(iRoi,rwdTypeNum,c,t) = nanmean(nullDeconv{iRoi,rwdTypeNum}((c-1)*deconvLength+t,:));%mean over voxels
                end
            end
            %get the null response
            taskResp{iRoi}(rwdTypeNum,:,:) = nullDeconv{iRoi,rwdTypeNum}(deconvLength+1:2*deconvLength,:); %for each voxel
%             meanTaskResp{iRoi}(rwdTypeNum,:) = meanDeconvNull(iRoi,rwdTypeNum,2,:); %mean over voxels
            %convole with task events
            nvox = roiTC{iRoi,rwdTypeNum}.n;
            for i=1:nvox
                temp = conv(taskMat, taskResp{iRoi}(rwdTypeNum,:,i));
                taskPrediction{iRoi,rwdTypeNum}(i,:) = temp(1:totalTrials*trialLength);
            end
            %subtract the task response
            
            
            
            % This has to be normalized in some way - the amplitude is too
            % small!!
            tcMinusTask{iRoi,rwdTypeNum} = roiTC{iRoi,rwdTypeNum}.tSeries - taskPrediction{iRoi,rwdTypeNum};
            
            
            
            
            
            %contrast betas
            contrastBetas{iRoi,rwdTypeNum} = designMatContrast{iScan}'\tcMinusTask{iRoi,rwdTypeNum}';
            %sf betas
            sfBetas{iRoi,rwdTypeNum} = designMatFreq{iScan}'\tcMinusTask{iRoi,rwdTypeNum}';
            %betas for combination of contrast and sf
            sfContrastBetas{iRoi,rwdTypeNum}= designMatContrastFreq{iScan}'\tcMinusTask{iRoi,rwdTypeNum}';
        end
        
    end
    deleteView(v);
    
    
    %% ECG
    %% RESP
    
    
end
%% FIGURES
for i=1:6
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
    title(rwdType);
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
    subplot(rows,cols,iRoi)
    plot(trialTC{iRoi}');
    legend('high','low');
    for rwdTypeNum=1:2
        j=i;
        
        figure(j)
        j=j+1;
        subplot(rows,cols,iRoi + (rwdTypeNum-1)*2)
        plot(squeeze(meanDeconvNull(iRoi, rwdTypeNum,:,:))', 'linewidth', 1);
        title([roiNames{iRoi} ': ' rwdTypes{rwdTypeNum}]);
        legend('stim', 'null')
        
        figure(j)
        j=j+1;
        subplot(rows,cols,iRoi)
        plot(nanmean(contrastBetas{iRoi,rwdTypeNum}(1:end-1,:),2), 'linewidth', 1);
        hold all
        title([roiNames{iRoi}]);
        
        subplot(rows,cols,cols+iRoi)
        plot(nanmean(contrastBetasWithTask{iRoi,rwdTypeNum}(1:end-1,:),2), 'linewidth', 1);
        hold all
        title([roiNames{iRoi}]);
        %         legend('high', 'low')
        
        figure(j)
        j=j+1;
        subplot(rows,cols,iRoi)
        plot(nanmean(sfBetas{iRoi,rwdTypeNum}(1:end-1,:),2), 'linewidth', 1);
        hold all
        title([roiNames{iRoi}]);
        subplot(rows,cols,cols+iRoi)
        plot(nanmean(sfBetasWithTask{iRoi,rwdTypeNum}(1:end-1,:),2), 'linewidth', 1);
        hold all
        title([roiNames{iRoi}]);
        
    end
end

%%



%%
%%%%%%%%%%%%%%%%%%%%%%%%%
% increases the sample rate of x by inserting n ? 1 zeros between samples. If x is a matrix, the function treats each column as a separate sequence.
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