%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%Nuno Pessanha Santos - santos.naamp@academiamilitar.pt
%Victor Lobo - vlobo@novaims.unl.pt
%André Dias 
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

function [angle]=Calc_Angle(FS, C, timedifference, S1, S2, d)

%FS - Sampling frequency (Hz)
%C - Sound water velocity (m/s)
%timedifference - Time difference between the arrival of the sound (n times FS)
%S1 - Hydrohpne number 1
%S1 - Hydrohpne number 2
%d - Physical distance between hydrophones (meters)

%Returns the angle in degrees
angle=acosd((timedifference*1/(FS*15)*C)/((S2-S1)*d));

end