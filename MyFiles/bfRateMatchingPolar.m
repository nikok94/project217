function out = bfRateMatchingPolar(in, K, E, varargin)
%%
 narginchk(3,4);
 if nargin==3
    iBIL = false;
 else
    iBIL = varargin{1};
 end
%% 
 N = length(in);

 if (floor(log2(N))~=log2(N))
     error('bfRateMatchingPolar_InvalidInputRMLength');
 end
 
 if iBIL % Uplink
   if ((K < 18) || (K > N) || (K>25 && K<31))
     error('bfRateMatchingPolar_UnsupportedKforUL');
   end 

 else % Downlink
   if ((K < 36) || (K > N))
     error('bfRateMatchingPolar_UnsupportedKforDL');
   end 
 end
 
 if (K>=18 && K<=25) % PC-Polar
    nPC = 3;
 else
    nPC = 0;
 end

 if (E > 8192 || E < K+nPC)
    rror('bfRateMatchingPolar_InvalidELength');
 end
%% Sub-block interleaving, Section 5.4.1.1
SbIP = [0 1 2 4 ...
        3 5 6 7 ...
        8 16 9 17 ...
        10 18 11 19 ...
        12 20 13 21 ...
        14 22 15 23 ...
        24 25 26 28 ...
        27 29 30 31]';
    
jn = zeros(N,1);
    for n = 0:N-1
        i = floor(32*n/N);
        jn(n+1) = SbIP(i+1)*(N/32)+ mod(n,N/32);
    end
    y = in(jn+1);

 %% Bit selection, Section 5.4.1.2
    N = length(in);
    outE = zeros(E,1,class(in));
    if E >= N
        % Bit repetition
        for k = 0:E-1
            outE(k+1) = y(mod(k,N)+1);
        end
    else
        if K/E <= 7/16
            % puncturing (take from the end)
            outE = y(end-E+1:end);
        else
            % shortening (take from the start)
            outE = y(1:E);
        end
    end
%% Interleaving, Section 5.4.1.3
    if iBIL
        % Specified for uplink only
      T =  ceil((-1+sqrt(1+8*length(outE)))/2);
      v = -1*ones(T,T);   % <NULL> bits
      k = 0;
      for i = 0:T-1
        for j = 0:T-1-i
            if k < length(outE)
                v(i+1,j+1) = outE(k+1);
            end
            k = k+1;
        end
      end
      out = zeros(size(outE));
      k = 0;
      for j = 0:T-1
          for i = 0:T-1-j
              if v(i+1,j+1) ~= -1
                  out(k+1) = v(i+1,j+1);
                  k = k+1;
              end
          end
      end     
    else
        % No interleaving
        out = outE;
    end
end

  