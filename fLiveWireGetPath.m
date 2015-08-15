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

function [iX, iY] = fLiveWireGetPath(iPX, iPY, iXS, iYS)
%FLIVEWIREGETPATH Traces the cheapest path from (IXS, IYS)^T through the
%pathmaps IPX, IPY back to the seed (where both, IPX and IPY are 0).
%
%   See also LIVEWIRE, FLIVEWIRECALCP, FLIVEWIREGETCOSTFCN.
%
%
%   Copyright 2013 Christian Würslin, University of Tübingen and University
%   of Stuttgart, Germany. Contact: christian.wuerslin@med.uni-tuebingen.de

iMAXPATH = 1000;

% -------------------------------------------------------------------------
% Initialize the variables
iPX  = int16(iPX);
iPY  = int16(iPY);
iXS = int16(iXS);
iYS = int16(iYS);

iX = zeros(iMAXPATH, 1, 'int16');
iY = zeros(iMAXPATH, 1, 'int16');

iLength = 1;
iX(iLength) = iXS;
iY(iLength) = iYS;
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% While not at the seed point: march back in the direction indicated by the
% path maps iPX (x-direction) and iPY (y-direction).
while (iPX(iYS, iXS) ~= 0) || (iPY(iYS, iXS) ~= 0) % We're not at the seed
    iXS = iXS + iPX(iYS, iXS);
    iYS = iYS + iPY(iYS, iXS);
    iLength = iLength + 1;
    iX(iLength) = iXS;
    iY(iLength) = iYS;
end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% revert vectors (to make it a forward path) and don't return the seed point.
iX = iX(iLength - 1:-1:1);
iY = iY(iLength - 1:-1:1);
% -------------------------------------------------------------------------