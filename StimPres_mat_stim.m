 function StimPres_mat_stim
 % modified from Stimulus_presentation and face_1000_para_Mar.m in Manabu's Lab and to add dio1(parallel
 % port) output to TDT ; this version with small size and rotation postition variation 110411 

close all; clear all;
fprintf('\n\n================================\n%s\n%s\n\n',...
    mfilename,datestr(now));

%% ---------------------------user definition 


condition = 'Takayuki104';

switch condition
    case 'Takayuki104'
        prm.stim = [1:104];
        prm.Stimmatfile  = 'D:\Takayuki\Takayuki104.mat';
        prm.FrameStim  = 42; % 42/85 ~=  500ms 
        prm.FrameBlank = 90; % 90/85 ~= 1050ms
        
    case 'Fragment300'
        prm.stim = [1:300];
        prm.Stimmatfile  = 'D:\Takayuki\2016101401\20161014_FragmentStim250.mat';
        prm.FrameStim  = 17; % 17/85~= 200ms
        prm.FrameBlank = 34; % 34/85~= 400ms
end


num.blocks = 12;
prm.DestStimSeqPath = 'D:\Takayuki\StimSeq';

% prm.Stimmatfile  = 'D:\Takayuki\Go890.mat';
% prm.Stimmatfile  = 'D:\Takayuki\2016100603\201610060302_stim14_10.mat';

prm.BGColor = 128; % 0=black  128=gray background 170=background of Cheston's stimuli

% prm.FrameStim  = 17;% 17/85~= 200ms 
% prm.FrameBlank = 34;% 34/85~= 400ms
% prm.FrameBlank = 17;% 17/85~= 200ms




prm.RefreshRate = 85; % [Hz]
prm.ScreenSize = [800 600];  % width height [pixels]
prm.PresenPos = prm.ScreenSize/2; % stimulus center
blankscreentime =5;  % [sec.]
prm.TextureOrient = 0;%-90;

%% main start ==========================================================
num_stim = length(prm.stim);
totalstimulinum = num.blocks*num_stim;

%% Genarate stimulus sequence
flag = 1;
while flag
    stimSeq = rand(num_stim,num.blocks);
    [~, stimSeq] = sort(stimSeq);
    flag = ~min(abs(diff(stimSeq(:))));
end
startTime = datestr(now,'yymmddHHMM');
stimSeqFile = fullfile(prm.DestStimSeqPath,['StimSeq_',startTime]);
save(stimSeqFile,'stimSeq','prm','num','startTime');
fprintf('Save %s\n\n',stimSeqFile);

%% Parallel port setting
% DIO1= digitalio('parallel','LPT1');
% addline(DIO1,0:7,0,'out');
% addline(DIO1,[0,1,3],2,'out'); % for 11 bits
% addline(DIO1,2,2,'out');    % for trigger
% sf_Put_value_Stim(DIO1, 0, 0);


% Parallel port setting
DIO1= digitalio('parallel','LPT1');
addline(DIO1,0:7,0,'out');
addline(DIO1,[0,1,3],2,'out'); % added by Takayuki for 11 bits
sf_Put_value_Stim(DIO1, 0, 0);

%% open main screen
AssertOpenGL;
screenNumber=max(Screen('Screens')); % max number of screen is usually the main control monitor
oldResolution = Screen('Resolution',screenNumber,prm.ScreenSize(1),prm.ScreenSize(2),prm.RefreshRate);
[win1, screenRect]=Screen('OpenWindow',screenNumber, prm.BGColor);
Screen('Flip', win1); % Flip buffer to show initial gray background:
HideCursor;

% confirm display parameters
disp_set = Screen('Resolution',screenNumber);
fprintf('intended screen setting: %4d x %4d  %3dHz\n',prm.ScreenSize(1),prm.ScreenSize(2),prm.RefreshRate);
fprintf('    real screen setting: %4d x %4d  %3dHz\n\n',disp_set.width,disp_set.height,disp_set.hz);
if ~isequal([prm.ScreenSize prm.RefreshRate],[disp_set.width disp_set.height disp_set.hz])
    error('Screen setting is wrong.');
end

%calling sync patch                 
% Screen('FillOval', win1, 0, photo_rectangle); % load sync patch onto bitmap                   
Screen('Flip', win1);  % show bitmap in window

%% road stimulus

fprintf('Read image file  %s\n', prm.Stimmatfile);
load(prm.Stimmatfile,'images','stimulus_name');

StimTex = zeros(1,num_stim);

for m=1:num_stim
    StimTex(m) = Screen('MakeTexture', win1, images{m});
    
    stim_size(:,m) = size(images{m});
    stim_shift(:,m) = [...
        -stim_size(2,m)/2,  -stim_size(1,m)/2,...
         stim_size(2,m)/2,  stim_size(1,m)/2];
     fprintf('%4d [%3d %3d %d] %s\n',m,stim_size(:,m),stimulus_name{m});
end
clear images
disp('Finish reading images');


stim_pos = repmat([...
    prm.PresenPos(1), prm.PresenPos(2),...
    prm.PresenPos(1), prm.PresenPos(2)]',1,num_stim)...
    + stim_shift;

% total time
one_trial = (prm.FrameStim+prm.FrameBlank)/prm.RefreshRate;
total_sec = one_trial*totalstimulinum;
total_min = fix(total_sec/60);
total_secm = ceil(mod(total_sec,60));
fprintf('\n\nEstimated presentation time: %dmin %dsec\n\n',...
    total_min,total_secm);

% wait
input('\nPress Enter for start\n');

%%
disp('Presentation start');
a = GetSecs;

counter = 0;
for repeat_num = 1:num.blocks
    for Stim = 1:num_stim
        
        counter = counter + 1;
%         StimID = S_seq(repeat_num,Stim);
        StimID = stimSeq(Stim,repeat_num);
        stim_posX = stim_pos(:,StimID); 
        remaining = total_sec-one_trial*counter;
        fprintf('%3d:%02d   %5d (%2d%%)  %2d %3d   StimID=%4d (%s)\n',...
            fix(remaining/60),fix(mod(remaining,60)),counter,...
            fix(counter/totalstimulinum*100),repeat_num,Stim,...
            StimID,stimulus_name{prm.stim(StimID)});
        
        % Blank presentation ===================
%         Screen('FillOval', win1, 0, photo_rectangle);
        Screen('Flip', win1, 0, 0, 0, 0);
        sf_Put_value_Stim(DIO1, StimID, 0);
        for frames = 2:prm.FrameBlank
%             Screen('FillOval', win1, 0, photo_rectangle);
            Screen('Flip', win1, 0, 0, 0, 0);            
        end
        
        % Stimulus presentation ================
        Screen('DrawTexture', win1, StimTex(StimID),[],stim_posX,prm.TextureOrient);
%         Screen('FillOval', win1,0, photo_rectangle);
%         Screen('FillOval', win1,255, photo_circle);
        Screen('Flip', win1, 0, 0, 0, 0);
        sf_Put_value_Stim(DIO1, StimID, 1);
        for frames = 2:prm.FrameStim
            pplace=round(rand(1)*2.5)+2; %these two commands are to randomly gitter the stimulus.
            stim_posY=stim_posX+pplace;  %these two commands are to randomly gitter the stimulus.
%             Screen('DrawTexture', win1, StimTex(StimID));
            Screen('DrawTexture', win1, StimTex(StimID),[],stim_posY,prm.TextureOrient);
%             Screen('FillOval', win1,  0, photo_rectangle);
            Screen('Flip', win1, 0, 0, 0, 0);
        end
        
    end
end

% END_OF_BLOCK blank screen

% Screen('FillOval', win1, 0, photo_rectangle);
Screen('Flip', win1, 0, 0, 0, 0);
sf_Put_value_Stim(DIO1, 0, 0);

elapsed_time = GetSecs-a;
fprintf('Presentation finish\n\n');
fprintf('%f\n',elapsed_time);

fprintf('\nwait for %d sec.\n\n',blankscreentime);
pause(blankscreentime);

Screen('CloseAll')
delete(DIO1)

fprintf('\n--------------------------------Fine.\n\n');


function sf_Put_value_Stim(dio_line,data,stim)


bvdata = dec2binvec(data,11);

bvdata(11) = logical(stim);

mask = logical([0 0 0 0 0 0 0 0 1 1 1]);
putvalue(dio_line,xor(mask,bvdata));


% bvdata = dec2binvec(data,11);
% 
% bvdata(12) = logical(stim);
% 
% mask = logical([0 0 0 0 0 0 0 0 1 1 1 0]);
% putvalue(dio_line,xor(mask,bvdata));

% function sf_Put_value(dio_line,data)
% mask = logical([0 0 0 0 0 0 0 0 1 1 1]);
% bvdata = dec2binvec(data,11);
% putvalue(dio_line,xor(mask,bvdata));

