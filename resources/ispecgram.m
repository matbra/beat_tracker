function x = ispecgram(d, ftsize, sr, win, nov)
% X = ispecgram(D, F, SR, WIN, NOV)           Inverse specgram
%    Overlap-add the inverse of the output of specgram
%    ftsize is implied by sizeof d, sr is ignored, nov defaults to ftsize/2
% dpwe 2005may16.  after istft

[nspec,ncol] = size(d);

if nargin < 2
  ftsize = 2*(nspec-1);
end
if nargin < 3
  % who cares?
end
if nargin < 4
  win = ftsize;  % doesn't matter either - assume it added up OK
end
if nargin < 5
  nov = ftsize/2;
end

hop = win - nov;

if nspec ~= (ftsize/2)+1
  error('number of rows should be fftsize/2+1')
end

window = sqrt(hanning(ftsize,'periodic'));

xlen = ftsize + (ncol-1) * hop;
x = zeros(xlen,1);

halff = ftsize/2;   % midpoint of win

% No reconstruction win (for now...)

for c = 1:ncol
  ft = d(:,c);
  ft = [ft(1:(ftsize/2+1)); conj(ft([(ftsize/2):-1:2]))];

  if max(imag(ifft(ft))) > 1e-5
    disp('imag oflow');
  end
  
  px = real(ifft(ft));  % no shift in specgram
  
%  figure;plot(px);pause;
  
  b = (c-1)*hop;
  x(b+[1:ftsize]) = x(b+[1:ftsize]) + window.*px;
end;

x = x * win/ftsize;  % scale amplitude

end

% Copyright (c) 2002-2010 Daniel P.W. Ellis & Columbia University
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
%
% Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.