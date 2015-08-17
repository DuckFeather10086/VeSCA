%{
Copyright (c) 2015, Christian Wuerslin, University of Tuebingen and University of Stuttgart, Germany
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in
      the documentation and/or other materials provided with the distribution
    * Neither the name of the Stanford University nor the names
      of its contributors may be used to endorse or promote products derived
      from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
%}

function dF = fLiveWireGetCostFcn(dImg, dWz, dWg, dWd)

% -------------------------------------------------------------------------
% Calculate the cost function

% THE GRADIENT DIRECTION COST FD
dFd = zeros(2047);
[dX, dY] = gradient(dImg);

dX = interleave_zeros(dX);
dY = interleave_zeros(dY);

for x = 1:2:size(dX, 2)-2
    for y = 1:2:size(dX, 1)-2
        
        % right
        p = [x; y];
        q = [x+2; y];
        
        Dp = [dY(x, y); -dX(x ,y)];
        Dq = [dY(x+2, y); -dX(x+2 ,y)];
        
        if Dp' * (q - p) >= 0
            Lpq = q - p;
        else
            Lpq = p - q;
        end
        
        dppq = Dp' * Lpq;
        dqpq = Lpq' * Dp;
        
        fDpq = 1/pi * (1/cos(dppq) + 1/cos(dqpq));
        
        dFd(x+1, y) = fDpq;
        
        % down
        p = [x; y];
        q = [x; y+2];
        
        Dp = [dY(x, y); -dX(x ,y)];
        Dq = [dY(x, y+2); -dX(x, y+2)];
        
        if Dp' * (q - p) >= 0
            Lpq = q - p;
        else
            Lpq = p - q;
        end
        
        dppq = Dp' * Lpq;
        dqpq = Lpq' * Dp;
        
        fDpq = 1/pi * (1/cos(dppq) + 1/cos(dqpq));
        
        dFd(x, y+1) = fDpq;
        
        % diag (down, right)
        p = [x; y];
        q = [x+2; y+2];
        
        Dp = [dY(x, y); -dX(x ,y)];
        Dq = [dY(x+2, y+2); -dX(x+2, y+2)];
        
        if Dp' * (q - p) >= 0
            Lpq = q - p;
        else
            Lpq = p - q;
        end
        
        dppq = Dp' * Lpq;
        dqpq = Lpq' * Dp;
        
        fDpq = 1/pi * (1/cos(dppq) + 1/cos(dqpq));
        
        dFd(x+1, y+1) = fDpq;
    end
end

% normalize
norm = max(max(dFd)) - min(min(dFd));
dFd = dFd ./ norm;
dFd = dFd - min(min(dFd));

% THE GRADIENT STRENGTH COST FG
dImg = double(dImg);
[Gmag, Gdir] = imgradient(dImg);
dFg = 1 - Gmag./max(Gmag(:));

% THE ZERO-CROSSING COST FZ
lFz = ~edge(dImg, 'zerocross');

% THE SUM:

% first add Fg and Fz
dF = dWz.*double(lFz)+ dWg.*dFg;

% then interleave zeros and fill in the gaps with Fd
dF = interleave_zeros(dF);
dF = dF + dWd .* dFd;

disp('hi')
% -------------------------------------------------------------------------