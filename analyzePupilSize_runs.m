dataFolder = '~/data/reward/';
% datafiles = {'181012_stim06','181012_stim07','181012_stim08'};
datafiles = {'181018_stim03','181018_stim04','181018_stim05'};%me, very tired
samplerate=500;

runs = length(datafiles);
for i=1:runs
    stimfile = [dataFolder datafiles{i} ];%'181003_stim07.mat'];
    % stimfile = [dataFolder '181010_stim01.mat'];
    e = getTaskEyeTraces(stimfile, 'removeBlink=7');
    s = load(stimfile);
    % plot(e.eye.pupil')
    %%
    meanPupil{i} = nanmean(e.eye.pupil)';
%     meanPupil = nanmean(e.eye.pupil(1:16,:))';
    figure(1)
    subplot(1,runs,i)
    plot(meanPupil{i})
    hold all
    stimTime = s.fixStimulus.stimTime * samplerate;
    interTime = s.fixStimulus.interTime * samplerate;
    responseTime = s.fixStimulus.responseTime * samplerate;
    allTimes = [stimTime interTime stimTime interTime responseTime];
    itime = 1;
    startT=1;
    endT = startT+allTimes(itime);
    plot(startT:endT-1,meanPupil{i}(startT:endT-1),'lineWidth',4);
    itime = itime+1;
    startT=endT;
    endT = startT+allTimes(itime);
    plot(startT:endT-1,meanPupil{i}(startT:endT-1),'lineWidth',4);
    itime = itime+1;
    startT=endT;
    endT = startT+allTimes(itime);
    plot(startT:endT-1,meanPupil{i}(startT:endT-1),'lineWidth',4);
    itime = itime+1;
    startT=endT;
    endT = startT+allTimes(itime);
    plot(startT:endT-1,meanPupil{i}(startT:endT-1),'lineWidth',4);
    itime = itime+1;
    startT=endT;
    endT = startT+allTimes(itime);
    plot(startT:endT-1,meanPupil{i}(startT:endT-1),'lineWidth',4);
    %%
    figure(2)
    subplot(1,runs,i)
    plot(e.eye.pupil')
    hold on
    plot(meanPupil{i},'color','k','linewidth',3);
    %%
    mean1 = nanmean(e.eye.pupil(s.task{1}{1}.correctResponse==1,:));
    std1 = std(e.eye.pupil(s.task{1}{1}.correctResponse==1,:),0,1,'omitnan');
    sem1 = std1./sqrt(size(e.eye.pupil(s.task{1}{1}.correctResponse==1,:),1));
    mean2 = nanmean(e.eye.pupil(s.task{1}{1}.correctResponse==2,:));
    std2 = std(e.eye.pupil(s.task{1}{1}.correctResponse==2,:),0,1,'omitnan');
    ste2 = std2./sqrt(size(e.eye.pupil(s.task{1}{1}.correctResponse==2,:),1));
    figure(3)
    subplot(1,runs,i)
    plot(mean1);
    hold all
    plot(mean2)
    % s.fixStimulus.stimTime;
    % s.fixStimulus.interTime;
    % s.fixStimulus.responseTime;
    
    %%
    % figure
    % mean1(isnan(mean1)) = nanmean(mean1);
    % std1(isnan(std1)) = 0;
    % sem1(isnan(sem1)) = 0;
    % % nanmean1 = mean1(~isnan(mean1));
    % % % nanstd1 = std1(~isnan(std1));
    % % pp_ = dsErrorsurface(1:length(nanmean1),nanmean1, nanstd1, [0 0 1], [0.5])
    % % nanstd1 = std1(~isnan(std1));
    % pp_ = dsErrorsurface(1:length(mean1),mean1, sem1, [0 0 1], [0.5])
    
end

%%
minTrialLength = 7000;
meanAll = meanPupil{1}(1:minTrialLength)+meanPupil{2}(1:minTrialLength)+meanPupil{3}(1:minTrialLength);

figure(4)
plot(meanAll)
hold all
itime = 1;
startT=1;
endT = startT+allTimes(itime);
plot(startT:endT-1,meanAll(startT:endT-1),'lineWidth',4);
itime = itime+1;
startT=endT;
endT = startT+allTimes(itime);
plot(startT:endT-1,meanAll(startT:endT-1),'lineWidth',4);
itime = itime+1;
startT=endT;
endT = startT+allTimes(itime);
plot(startT:endT-1,meanAll(startT:endT-1),'lineWidth',4);
itime = itime+1;
startT=endT;
endT = startT+allTimes(itime);
plot(startT:endT-1,meanAll(startT:endT-1),'lineWidth',4);
itime = itime+1;
startT=endT;
endT = startT+allTimes(itime);
plot(startT:endT-1,meanAll(startT:endT-1),'lineWidth',4);
