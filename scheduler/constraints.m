%   Copyright 2023 Zuse Institute Berlin
%
%   Licensed under the Apache License, Version 2.0 (the "License");
%   you may not use this file except in compliance with the License.
%   You may obtain a copy of the License at
%
%       http://www.apache.org/licenses/LICENSE-2.0
%
%   Unless required by applicable law or agreed to in writing, software
%   distributed under the License is distributed on an "AS IS" BASIS,
%   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%   See the License for the specific language governing permissions and
%   limitations under the License.


function [c, ceq] = constraints(x)

global conflict;
global optx;
global window;
% global msg;

% local_msg = '';

K = 999999999;

% constraints
A = @(xi, si, xj, sj) xi + si - xj - sj < 0;
cij = @(xi, si, wi, xj, sj) wi + si + xi - xj - sj - K * (1 - A(xi, si, xj, sj));

%c = zeros(size(conflict, 2) * (size(conflict, 2) - 1) + 2, 1);
n = size(conflict, 2);
c = zeros(n*(n + 1), 1);
%ceq = zeros(n);
zz = 1;
for ii = 1:size(conflict, 2)
    for jj = 1:size(conflict, 2)
        if(ii ~= jj)
            c(zz) = cij(x(ii), conflict{ii}.time, conflict{ii}.duration, x(jj), conflict{jj}.time);
            zz = zz + 1;
%             newmsg = strcat('i: ', num2str(ii), ', j: ', num2str(jj), ', xi: ', num2str(x(ii)), ', xj: ', ...
%                 num2str(x(jj)), ' si: ', num2str(conflict{ii}.time), ', wi: ',...
%                 num2str(conflict{ii}.duration), ', sj: ', num2str(conflict{jj}.time), ', c: ', num2str(c(zz - 1)));
%             arr = {local_msg, newmsg};
%             local_msg = strjoin(arr, '\n');
            if(strcmp(conflict{ii}.owner.id, conflict{jj}.owner.id) == 1 && conflict{ii}.index == conflict{jj}.index - 1)
                % cp s of same application should avoid overlaps
                c(zz) = x(ii) - x(jj) - optx{ii};
                zz = zz + 1;
            end
        end
    end
    % x < opt
    c(zz) = abs(x(ii)) - optx{ii} + 1;
    zz = zz + 1;
    
%     newmsg = strcat('i: ', int2str(ii), ' < opt', ', c: ', num2str(c(zz - 1)));
%     arr = {local_msg, newmsg};
%     local_msg = strjoin(arr, '\n');
 
    % x in window
    c(zz) = window{1} - x(ii) - conflict{ii}.time;
    zz = zz + 1;
    c(zz) = x(ii) + conflict{ii}.time + conflict{ii}.duration - window{2};
    zz = zz + 1;
    
    %l = double(conflict{ii}.owner.natural_intv);
    %r = double(rem(x(ii), l));
    %ceq(ii) = sin(pi * (r/l + mod(floor(x(ii)/l), 2)));
    %zz = zz + 1;
end
% no = 0;
% for ii = 1:(zz - 1)
%     if(c(ii) > 0)
%         no = 1;
%     end
% end
% 
% if(no == 0)
%     msg = strjoin({msg, '=====new test=====', local_msg}, '\n');
% end

ceq = [];
