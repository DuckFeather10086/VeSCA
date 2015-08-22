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

% the indices of dX and dY which are not newly interleaved zeros
ind = logical(interleave_zeros(ones(size(dImg))));

% gradient direction cost right
Dp = [dY(ind), -dX(ind)];
Dp(end-1023:end, :) = [];
Dq = [dY(ind), -dX(ind)];
Dq(1:1024, :) = [];

Lpq = [ones(length(Dq), 1) zeros(length(Dq), 1)];
dppq = sum(Lpq .* Dp, 2);
dppq(dppq<0) = dppq(dppq<0) * -1;
dqpq = sum(Lpq .* Dq, 2);
dqpq(dqpq<0) = dqpq(dqpq<0) * -1;

fDpq = 1/pi .* (ones(size(dppq)) ./ cos(dppq) + ones(size(dqpq)) ./ cos(dqpq));
fDpq = reshape(fDpq, 1024, 1023);
fDpq = interleave_zeros(fDpq);
fDpq = [zeros(2047, 1) fDpq zeros(2047, 1)];

dFd = dFd + fDpq;

% gradient direction cost down

Dp = [dY(ind), -dX(ind)];
Dp(mod(1:1024^2, 1024)==0, :) = [];
Dq = [dY(ind), -dX(ind)];
Dq(mod(1:1024^2, 1024)==1, :) = [];

Lpq = [zeros(length(Dq), 1) ones(length(Dq), 1)];
dppq = sum(Lpq .* Dp, 2);
dppq(dppq<0) = dppq(dppq<0) * -1;
dqpq = sum(Lpq .* Dq, 2);
dqpq(dqpq<0) = dqpq(dqpq<0) * -1;

fDpq = 1/pi .* (ones(size(dppq)) ./ cos(dppq) + ones(size(dqpq)) ./ cos(dqpq));
fDpq = reshape(fDpq, 1023, 1024);
fDpq = interleave_zeros(fDpq);
fDpq = [zeros(1, 2047); fDpq; zeros(1, 2047)];

dFd = dFd + fDpq;

% gradient direction cost diagonally down+right
Dp = [dY(ind), -dX(ind)];
Dp(mod(1:1024^2, 1024)==0, :) = [];
Dp(end-1022:end, :) = [];

Dq = [dY(ind), -dX(ind)];
Dq(mod(1:1024^2, 1024)==1, :) = [];
Dq(1:1023, :) = [];

Lpq = ones(length(Dq), 2);
dppq = sum(Lpq .* Dp, 2);
dppq(dppq<0) = dppq(dppq<0) * -1;
dqpq = sum(Lpq .* Dq, 2);
dqpq(dqpq<0) = dqpq(dqpq<0) * -1;

fDpq = 1/pi .* (ones(size(dppq)) ./ cos(dppq) + ones(size(dqpq)) ./ cos(dqpq));
fDpq = reshape(fDpq, 1023, 1023);
fDpq = interleave_zeros(fDpq);
fDpq = [zeros(1, 2045); fDpq; zeros(1, 2045)];
fDpq = [zeros(2047, 1) fDpq zeros(2047, 1)];

dFd = dFd + fDpq;

% normalize
dFd(dFd==0) = min(min(dFd));
dFd = dFd - min(min(dFd));
dFd = dFd ./ max(max(dFd));

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
% -------------------------------------------------------------------------