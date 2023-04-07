%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%Krait Array Trials
%Nuno Pessanha Santos - santos.naamp@academiamilitar.pt
%Victor Lobo - vlobo@novaims.unl.pt
%André Dias 
%Last update - 05/04/2023
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialization - Pre-Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NUMPOINTS = 30000/16; %Number of 
%NUMUDPPACKETS = 29141;
NUMUDPPACKETS = 6000;
NUMPOINTS = NUMUDPPACKETS*15;
NCHANNELS = 16; %Number of channels - [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
%NCHANNELS = 32; %Number of channels - [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
SAVE_WAV = 1; %SAVE == 1 // NO_SAVE == 0
SAVE_RESULTS = 1; %SAVE == 1 // NO_SAVE == 0
SAMPLING_RATE_CALCULATION = 20; %Packets used to calculate the sampling rate per channel
DIVIDE_CHANNELS = 20; %Divide each channel 24

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Aux Variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COUNT_POINTS = 1000; %Count cycle points
temp_data=zeros(SAMPLING_RATE_CALCULATION,1); %Initialize frequency calculation variable
temp = 1; %Concatenate array - Temporary variable

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initilization - Variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hydro_1 = zeros(NUMPOINTS,NCHANNELS); %Initialize Hydrophones variable - Each column is one hydrophone
Hydro_2 = zeros(NUMPOINTS*NCHANNELS,1); %Initialize Hydrophones variable - All saved sequentially
%".wav" files
FILENAME_CH_1 = 'Channel_1.wav'; %Channel 1
FILENAME_CH_2 = 'Channel_2.wav'; %Channel 2
FILENAME_CH_3 = 'Channel_3.wav'; %Channel 3
FILENAME_CH_4 = 'Channel_4.wav'; %Channel 4
FILENAME_CH_5 = 'Channel_5.wav'; %Channel 5
FILENAME_CH_6 = 'Channel_6.wav'; %Channel 6
FILENAME_CH_7 = 'Channel_7.wav'; %Channel 7
FILENAME_CH_8 = 'Channel_8.wav'; %Channel 8
FILENAME_CH_9 = 'Channel_9.wav'; %Channel 9
FILENAME_CH_10 = 'Channel_10.wav'; %Channel 10
FILENAME_CH_11 = 'Channel_11.wav'; %Channel 11
FILENAME_CH_12 = 'Channel_12.wav'; %Channel 12
FILENAME_CH_13 = 'Channel_13.wav'; %Channel 13
FILENAME_CH_14 = 'Channel_14.wav'; %Channel 14
FILENAME_CH_15 = 'Channel_15.wav'; %Channel 15
FILENAME_CH_16 = 'Channel_16.wav'; %Channel 16
FILENAME_CH_ALL = 'Channel_ALL.wav'; %All Channels in a single ".wav" file
FILENAME_CH_ALL_NORM = 'Channel_ALL_NORM.wav'; %All Channels in a single ".wav" file
FILENAME_CH_HYDR = 'Hydro_Values.mat'; %Save Hydrophone values
FILENAME_CH_HYDR_CONC = 'Hydro_Values_concatenated.mat'; %Save Hydrophone values

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Communications
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u = udpport("LocalPort",51001); %UDP connection
u.InputBufferSize = 15040000000; % UDP Buffer
fopen(u) %Open connection
u.flush();  %Clean Buffer

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Obtain automatically the Sampling Rate - Difference beween 2 packets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u.flush();  %Clean Buffer
%Adquire some packets - "Time of day in ns"
for i=1:SAMPLING_RATE_CALCULATION
    count_read=0;
    disp('estou a escuta');
    while count_read==0
        [data, count_read] = fread(u, 1464, 'uint8');
    end
    disp('comecei a ler');
    temp_data(i,1)=data(21)*(2^24)+data(22)*(2^16)+data(23)*(2^8)+data(24);
end

%Calculate mean sampling rate
SAMPLE_RATE_calculated = round(15/(mean(diff(temp_data))*10^-9));
SAMPLE_RATE = 10000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Main - Hydrophone Values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u.flush();  %Clean Buffer

NUMPOINTS = 0;
index_jump = 96;

%for k=1:NUMPOINTS %Number of samples
for k=1:NUMUDPPACKETS %Number of UDP Packets
   data = fread(u, 1464, 'uint8');

    index = 0;

   %NUMPOINTS = NUMPOINTS + 15

   %Show point counter
   %if(mod(k,COUNT_POINTS) == 0)
   %    k
   %end

   for sets=1:15

    NUMPOINTS = NUMPOINTS + 1;

   for j=1:NCHANNELS

       %Save the packet content to a variable - Each column is one hydro
       m=(j-1)*3;
       
 

       %Hydro_1(k,j) = data(25+m+index)*65536+data(26+m+index)*256+data(27+m+index);
       Hydro_1(NUMPOINTS,j) = data(25+m+index)*65536+data(26+m+index)*256+data(27+m+index);
       
       %Concatenate values - To save a single ".wav" file
       %Hydro_2(temp,1) = Hydro_1(k,j);
       Hydro_2(temp,1) = Hydro_1(NUMPOINTS,j);
    
       % Convert from signed int to number
       %if(Hydro_1(k,j) >= 2^(23))
       %  %Hydrophone values
       %  Hydro_1(k,j) = - (2^(24) - Hydro_1(k,j));
       %  %Concatenate
       %  Hydro_2(temp,1) =  Hydro_1(k,j);
       %end
       if(Hydro_1(NUMPOINTS,j) >= 2^(23))
         %Hydrophone values
         Hydro_1(NUMPOINTS,j) = - (2^(24) - Hydro_1(NUMPOINTS,j));
         %Concatenate
         Hydro_2(temp,1) =  Hydro_1(NUMPOINTS,j);
       end
       %Aux variable - concatenate values
       temp=temp+1;
   end
    index = index + index_jump

   end % NO SETS
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
plot(Hydro_1)

figure()
plot(Hydro_2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Save data into "wav" file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Save all data into ".wav" file
if(SAVE_WAV == 1)
    %Delete older files - Check if it is important later
    if(exist('Channel_1.wav')==2)
    delete Channel_1.wav
    end
    if(exist('Channel_2.wav')==2)
    delete Channel_2.wav
    end
    if(exist('Channel_3.wav')==2)
    delete Channel_3.wav
    end
    if(exist('Channel_4.wav')==2)
    delete Channel_4.wav
    end
    if(exist('Channel_5.wav')==2)
    delete Channel_5.wav
    end
    if(exist('Channel_6.wav')==2)
    delete Channel_6.wav
    end
    if(exist('Channel_7.wav')==2)
    delete Channel_7.wav
    end
    if(exist('Channel_8.wav')==2)
    delete Channel_8.wav
    end
    if(exist('Channel_9.wav')==2)
    delete Channel_9.wav
    end
    if(exist('Channel_10.wav')==2)
    delete Channel_10.wav
    end
    if(exist('Channel_11.wav')==2)
    delete Channel_11.wav
    end
    if(exist('Channel_12.wav')==2)
    delete Channel_12.wav
    end
    if(exist('Channel_13.wav')==2)
    delete Channel_13.wav
    end
    if(exist('Channel_14.wav')==2)
    delete Channel_14.wav
    end
    if(exist('Channel_15.wav')==2)
    delete Channel_15.wav
    end
    if(exist('Channel_16.wav')==2)
    delete Channel_16.wav
    end
    if(exist('Channel_ALL.wav')==2)
    delete Channel_ALL.wav
    end

%Save new files
switch NCHANNELS
    case 1
    audiowrite(FILENAME_CH_1,int16(Hydro_1(:,1)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_ALL,int16(Hydro_2(:,1)),SAMPLE_RATE*NCHANNELS);
    case 2
    audiowrite(FILENAME_CH_1,int16(Hydro_1(:,1)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_2,int16(Hydro_1(:,2)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_ALL,int16(Hydro_2(:,1)),SAMPLE_RATE*NCHANNELS);
    case 3
    audiowrite(FILENAME_CH_1,int16(Hydro_1(:,1)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_2,int16(Hydro_1(:,2)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_3,int16(Hydro_1(:,3)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_ALL,int16(Hydro_2(:,1)),SAMPLE_RATE*NCHANNELS);
    case 4
    audiowrite(FILENAME_CH_1,int16(Hydro_1(:,1)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_2,int16(Hydro_1(:,2)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_3,int16(Hydro_1(:,3)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_4,int16(Hydro_1(:,4)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_ALL,int16(Hydro_2(:,1)),SAMPLE_RATE*NCHANNELS);
    case 5
    audiowrite(FILENAME_CH_1,int16(Hydro_1(:,1)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_2,int16(Hydro_1(:,2)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_3,int16(Hydro_1(:,3)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_4,int16(Hydro_1(:,4)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_5,int16(Hydro_1(:,5)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_ALL,int16(Hydro_2(:,1)),SAMPLE_RATE*NCHANNELS);
    case 6
    audiowrite(FILENAME_CH_1,int16(Hydro_1(:,1)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_2,int16(Hydro_1(:,2)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_3,int16(Hydro_1(:,3)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_4,int16(Hydro_1(:,4)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_5,int16(Hydro_1(:,5)),SAMPLE_RATE);  
    audiowrite(FILENAME_CH_6,int16(Hydro_1(:,6)),SAMPLE_RATE);  
    audiowrite(FILENAME_CH_ALL,int16(Hydro_2(:,1)),SAMPLE_RATE*NCHANNELS);
    case 7
    audiowrite(FILENAME_CH_1,int16(Hydro_1(:,1)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_2,int16(Hydro_1(:,2)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_3,int16(Hydro_1(:,3)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_4,int16(Hydro_1(:,4)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_5,int16(Hydro_1(:,5)),SAMPLE_RATE);  
    audiowrite(FILENAME_CH_6,int16(Hydro_1(:,6)),SAMPLE_RATE);  
    audiowrite(FILENAME_CH_7,int16(Hydro_1(:,7)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_ALL,int16(Hydro_2(:,1)),SAMPLE_RATE*NCHANNELS);
    case 8
    audiowrite(FILENAME_CH_1,int16(Hydro_1(:,1)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_2,int16(Hydro_1(:,2)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_3,int16(Hydro_1(:,3)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_4,int16(Hydro_1(:,4)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_5,int16(Hydro_1(:,5)),SAMPLE_RATE);  
    audiowrite(FILENAME_CH_6,int16(Hydro_1(:,6)),SAMPLE_RATE);  
    audiowrite(FILENAME_CH_7,int16(Hydro_1(:,7)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_8,int16(Hydro_1(:,8)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_ALL,int16(Hydro_2(:,1)),SAMPLE_RATE*NCHANNELS);
    case 9
    audiowrite(FILENAME_CH_1,int16(Hydro_1(:,1)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_2,int16(Hydro_1(:,2)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_3,int16(Hydro_1(:,3)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_4,int16(Hydro_1(:,4)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_5,int16(Hydro_1(:,5)),SAMPLE_RATE);  
    audiowrite(FILENAME_CH_6,int16(Hydro_1(:,6)),SAMPLE_RATE);  
    audiowrite(FILENAME_CH_7,int16(Hydro_1(:,7)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_8,int16(Hydro_1(:,8)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_9,int16(Hydro_1(:,9)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_ALL,int16(Hydro_2(:,1)),SAMPLE_RATE*NCHANNELS);
    case 10
    audiowrite(FILENAME_CH_1,int16(Hydro_1(:,1)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_2,int16(Hydro_1(:,2)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_3,int16(Hydro_1(:,3)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_4,int16(Hydro_1(:,4)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_5,int16(Hydro_1(:,5)),SAMPLE_RATE);  
    audiowrite(FILENAME_CH_6,int16(Hydro_1(:,6)),SAMPLE_RATE);  
    audiowrite(FILENAME_CH_7,int16(Hydro_1(:,7)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_8,int16(Hydro_1(:,8)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_9,int16(Hydro_1(:,9)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_10,int16(Hydro_1(:,10)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_ALL,int16(Hydro_2(:,1)),SAMPLE_RATE*NCHANNELS);
    case 11
    audiowrite(FILENAME_CH_1,int16(Hydro_1(:,1)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_2,int16(Hydro_1(:,2)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_3,int16(Hydro_1(:,3)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_4,int16(Hydro_1(:,4)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_5,int16(Hydro_1(:,5)),SAMPLE_RATE);  
    audiowrite(FILENAME_CH_6,int16(Hydro_1(:,6)),SAMPLE_RATE);  
    audiowrite(FILENAME_CH_7,int16(Hydro_1(:,7)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_8,int16(Hydro_1(:,8)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_9,int16(Hydro_1(:,9)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_10,int16(Hydro_1(:,10)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_11,int16(Hydro_1(:,11)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_ALL,int16(Hydro_2(:,1)),SAMPLE_RATE*NCHANNELS);
    case 12
    audiowrite(FILENAME_CH_1,int16(Hydro_1(:,1)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_2,int16(Hydro_1(:,2)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_3,int16(Hydro_1(:,3)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_4,int16(Hydro_1(:,4)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_5,int16(Hydro_1(:,5)),SAMPLE_RATE);  
    audiowrite(FILENAME_CH_6,int16(Hydro_1(:,6)),SAMPLE_RATE);  
    audiowrite(FILENAME_CH_7,int16(Hydro_1(:,7)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_8,int16(Hydro_1(:,8)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_9,int16(Hydro_1(:,9)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_10,int16(Hydro_1(:,10)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_11,int16(Hydro_1(:,11)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_12,int16(Hydro_1(:,12)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_ALL,int16(Hydro_2(:,1)),SAMPLE_RATE*NCHANNELS);
    case 13
    audiowrite(FILENAME_CH_1,int16(Hydro_1(:,1)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_2,int16(Hydro_1(:,2)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_3,int16(Hydro_1(:,3)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_4,int16(Hydro_1(:,4)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_5,int16(Hydro_1(:,5)),SAMPLE_RATE);  
    audiowrite(FILENAME_CH_6,int16(Hydro_1(:,6)),SAMPLE_RATE);  
    audiowrite(FILENAME_CH_7,int16(Hydro_1(:,7)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_8,int16(Hydro_1(:,8)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_9,int16(Hydro_1(:,9)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_10,int16(Hydro_1(:,10)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_11,int16(Hydro_1(:,11)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_12,int16(Hydro_1(:,12)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_13,int16(Hydro_1(:,13)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_ALL,int16(Hydro_2(:,1)),SAMPLE_RATE*NCHANNELS);
    case 14
    audiowrite(FILENAME_CH_1,int16(Hydro_1(:,1)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_2,int16(Hydro_1(:,2)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_3,int16(Hydro_1(:,3)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_4,int16(Hydro_1(:,4)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_5,int16(Hydro_1(:,5)),SAMPLE_RATE);  
    audiowrite(FILENAME_CH_6,int16(Hydro_1(:,6)),SAMPLE_RATE);  
    audiowrite(FILENAME_CH_7,int16(Hydro_1(:,7)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_8,int16(Hydro_1(:,8)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_9,int16(Hydro_1(:,9)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_10,int16(Hydro_1(:,10)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_11,int16(Hydro_1(:,11)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_12,int16(Hydro_1(:,12)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_13,int16(Hydro_1(:,13)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_14,int16(Hydro_1(:,14)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_ALL,int16(Hydro_2(:,1)),SAMPLE_RATE*NCHANNELS);
    case 15
    audiowrite(FILENAME_CH_1,int16(Hydro_1(:,1)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_2,int16(Hydro_1(:,2)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_3,int16(Hydro_1(:,3)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_4,int16(Hydro_1(:,4)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_5,int16(Hydro_1(:,5)),SAMPLE_RATE);  
    audiowrite(FILENAME_CH_6,int16(Hydro_1(:,6)),SAMPLE_RATE);  
    audiowrite(FILENAME_CH_7,int16(Hydro_1(:,7)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_8,int16(Hydro_1(:,8)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_9,int16(Hydro_1(:,9)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_10,int16(Hydro_1(:,10)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_11,int16(Hydro_1(:,11)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_12,int16(Hydro_1(:,12)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_13,int16(Hydro_1(:,13)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_14,int16(Hydro_1(:,14)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_15,int16(Hydro_1(:,15)),SAMPLE_RATE);
    audiowrite(FILENAME_CH_ALL,int16(Hydro_2(:,1)),SAMPLE_RATE*NCHANNELS);
    case 16
    audiowrite(FILENAME_CH_1,int16(Hydro_1(:,1)/DIVIDE_CHANNELS),SAMPLE_RATE);
    audiowrite(FILENAME_CH_2,int16(Hydro_1(:,2)/DIVIDE_CHANNELS),SAMPLE_RATE);
    audiowrite(FILENAME_CH_3,int16(Hydro_1(:,3)/DIVIDE_CHANNELS),SAMPLE_RATE);
    audiowrite(FILENAME_CH_4,int16(Hydro_1(:,4)/DIVIDE_CHANNELS),SAMPLE_RATE);
    audiowrite(FILENAME_CH_5,int16(Hydro_1(:,5)/DIVIDE_CHANNELS),SAMPLE_RATE);  
    audiowrite(FILENAME_CH_6,int16(Hydro_1(:,6)/DIVIDE_CHANNELS),SAMPLE_RATE);  
    audiowrite(FILENAME_CH_7,int16(Hydro_1(:,7)/DIVIDE_CHANNELS),SAMPLE_RATE);
    audiowrite(FILENAME_CH_8,int16(Hydro_1(:,8)/DIVIDE_CHANNELS),SAMPLE_RATE);
    audiowrite(FILENAME_CH_9,int16(Hydro_1(:,9)/DIVIDE_CHANNELS),SAMPLE_RATE);
    audiowrite(FILENAME_CH_10,int16(Hydro_1(:,10)/DIVIDE_CHANNELS),SAMPLE_RATE);
    audiowrite(FILENAME_CH_11,int16(Hydro_1(:,11)/DIVIDE_CHANNELS),SAMPLE_RATE);
    audiowrite(FILENAME_CH_12,int16(Hydro_1(:,12)/DIVIDE_CHANNELS),SAMPLE_RATE);
    audiowrite(FILENAME_CH_13,int16(Hydro_1(:,13)/DIVIDE_CHANNELS),SAMPLE_RATE);
    audiowrite(FILENAME_CH_14,int16(Hydro_1(:,14)/DIVIDE_CHANNELS),SAMPLE_RATE);
    audiowrite(FILENAME_CH_15,int16(Hydro_1(:,15)/DIVIDE_CHANNELS),SAMPLE_RATE);
    audiowrite(FILENAME_CH_16,int16(Hydro_1(:,16)/DIVIDE_CHANNELS),SAMPLE_RATE);
    audiowrite(FILENAME_CH_ALL,int16(Hydro_2(:,1)/DIVIDE_CHANNELS),SAMPLE_RATE*NCHANNELS);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Save results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(SAVE_RESULTS == 1)
save(FILENAME_CH_HYDR, 'Hydro_1');
save(FILENAME_CH_HYDR_CONC, 'Hydro_2');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Communications
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delete(u);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Resume - Debug
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LogicalStr = {'false', 'true'};
fprintf('Number of adquired Channels: %d \n',NCHANNELS);
fprintf('Number of adquired Points: %d \n',NUMPOINTS);
fprintf('Sample Rate (calculated): %d Hz\n',SAMPLE_RATE);
fprintf('Sample Rate x16 (calculated): %d Hz\n',SAMPLE_RATE*16);
fprintf('Save the results (.mat): %s \n',LogicalStr{SAVE_RESULTS + 1});
fprintf('Save the results (.wav): %s \n',LogicalStr{SAVE_WAV + 1});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Normalize each Channel - Concatenate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEAN_CHANNELS = mean(abs(Hydro_1));
% HYDRO_NORM = Hydro_1./MEAN_CHANNELS;
% %DEBUG_MEAN_CHANNEL = mean(abs(HYDRO_NORM));
% 
% %Concatenate on the same matrix
% tmp = HYDRO_NORM(1,:);
% for(i=2:NUMPOINTS)
% tmp=[tmp HYDRO_NORM(i,:)];
% end

% tmp = tmp/max(abs(tmp));
% audiowrite(FILENAME_CH_ALL_NORM,double(tmp),SAMPLE_RATE*NCHANNELS);
% 
% 
% figure()
% plot(HYDRO_NORM)


