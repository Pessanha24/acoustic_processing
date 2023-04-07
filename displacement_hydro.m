%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%Nuno Pessanha Santos - santos.naamp@academiamilitar.pt
%Victor Lobo - vlobo@novaims.unl.pt
%André Dias 
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

function [timedifference]=timedifference(A,B)

%Size of the vector
%Vector size A = B
N=length(A);

%Obtain correlation between vectors
x=xcorr(A,B);

%Obtain the maximum
[m,i]=max(x);

%Displacement of the maximum
%Time interval between the hydrohpne sound arrival
timedifference = i-N;

end
