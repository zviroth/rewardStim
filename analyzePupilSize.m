dataFolder = '/Users/rothzn/data/reward/';
% datafiles = {'181012_stim07','181012_stim07','181012_stim08'};
samplerate=500;

% runs = length(datafiles);

% stimfile = [dataFolder datafiles{i} ];%'181003_stim07.mat'];
stimfile = [dataFolder '181003_stim07.mat'];
% stimfile = [dataFolder '181010_stim01.mat'];

e = getTaskEyeTraces(stimfile, 'removeBlink=3');
s = load(stimfile);
%edf = mglEyelinkEDFRead(sprintf('%s%s.edf', dataFolder, s.myscreen.eyetracker.datafilename));

% plot(e.eye.pupil')
%%
meanPupil = nanmean(e.eye.pupil)';

figure(1)
plot(meanPupil)
hold all
stimTime = s.fixStimulus.stimTime * samplerate;
interTime = s.fixStimulus.interTime * samplerate;
responseTime = s.fixStimulus.responseTime * samplerate;
allTimes = [stimTime interTime stimTime interTime responseTime];
itime = 1;
startT=1;
endT = startT+allTimes(itime);
plot(startT:endT-1,meanPupil(startT:endT-1),'lineWidth',4);
itime = itime+1;
startT=endT;
endT = startT+allTimes(itime);
plot(startT:endT-1,meanPupil(startT:endT-1),'lineWidth',4);
itime = itime+1;
startT=endT;
endT = startT+allTimes(itime);
plot(startT:endT-1,meanPupil(startT:endT-1),'lineWidth',4);
itime = itime+1;
startT=endT;
endT = startT+allTimes(itime);
plot(startT:endT-1,meanPupil(startT:endT-1),'lineWidth',4);
itime = itime+1;
startT=endT;
endT = startT+allTimes(itime);
plot(startT:endT-1,meanPupil(startT:endT-1),'lineWidth',4);
%%
figure
plot(e.eye.pupil')
hold on
plot(meanPupil,'color','k','linewidth',3);
%%
mean1 = nanmean(e.eye.pupil(s.task{1}{1}.correctResponse==1,:));
std1 = std(e.eye.pupil(s.task{1}{1}.correctResponse==1,:),0,1,'omitnan');
sem1 = std1./sqrt(size(e.eye.pupil(s.task{1}{1}.correctResponse==1,:),1));
mean2 = nanmean(e.eye.pupil(s.task{1}{1}.correctResponse==2,:));
std2 = std(e.eye.pupil(s.task{1}{1}.correctResponse==2,:),0,1,'omitnan');
ste2 = std2./sqrt(size(e.eye.pupil(s.task{1}{1}.correctResponse==2,:),1));
figure(2)
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

