function kertemp = kertemp(k,A,ker)

% Computes integration of a potential matrix A by succesive convolution against
% the heat kernel in all dimensions except the k-th

% More specifically, let 
% - A be a matrix of potentials of dimension (N+1)*S 
% - ker a S*S kernel. 
% 
% Define P the S^*(N+1) path measure :
%   P(i_0,...,i_N) = ker(i_1-i_0)...ker(i_N-i_{N-1})
%
% Denote as XA the tensor product of A over paths :
%   XA(i_0,...,i_N) = A(0,i_0)*...*A(N,i_N)
%
% => When ker is the Heat Kernel, P is the discretized Wiener measure.

% This function computes :
% Sum_{ i_0, ... i_{k-1}, i_{k+1}, ... , N } (XA*P)(i_0,...,i_k,i,i_{k+1},...i_N)
% As a vector in i.

[Nplus,~] = size(A) ;
N = Nplus - 1 ;

if k == 1
    % When k=1, integration can be computed by backwards convolution
    kertemp = sum( ker .* A(N+1,:),2)' ;
    for i = 1:N-1
        kertemp = sum(ker .* A(N+1-i,:) .* kertemp,2)' ; 
    end 
end

if (k >= 2) && (k <= N)
    
    ker_backward = sum( ker .* A(N+1,:),2)' ;
    for i = 1:N-k
        ker_backward= sum(ker .* A(N+1-i,:) .* ker_backward,2)' ; 
    end 
    
    ker_forward = sum( ker' .* A(1,:),2)' ;
    for i = 2:k-1
        ker_forward = sum(ker' .* A(i,:) .* ker_forward,2)' ; 
    end 
    
    kertemp = ker_forward.*ker_backward ;
end

if k == N+1
    kertemp = sum( ker' .* A(1,:),2)' ;
    for i = 2:N
        kertemp = sum(ker' .* A(i,:) .* kertemp,2)' ; 
    end 
end

end