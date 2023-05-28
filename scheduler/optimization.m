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



d = 1;
fu = @(x, w, t, R, MTTI) abs(MTTI * exp(R / MTTI) * ...
        (exp((t + x + w) / MTTI) - 1) / (t + x));
f = @(x, w, t, R, MTTI) t * fu(x, w, t, R, MTTI);

global conflict; 
global optx; 
global msg;
global window;

msg = '';
conflict = cell(p_conflict);
optx = cell(p_optx);
mx = cell(p_mx);
rx = cell(p_rx);
window = cell(p_window);

window{1} = double(window{1});
window{2} = double(window{2});

diff_to_prev_index = zeros(size(conflict, 2));

for ii = 1:size(conflict, 2)
    conflict{ii}.time = double(conflict{ii}.time);
    conflict{ii}.duration = double(conflict{ii}.duration);
    optx{ii} = double(optx{ii});
    mx{ii} = double(mx{ii});
    rx{ii} = double(rx{ii});
    
    for jj = 1:size(conflict, 2)
        if(strcmp(conflict{ii}.owner.id, conflict{jj}.owner.id) == 1 && conflict{jj}.index == conflict{ii}.index - 1)
            diff_to_prev_index(ii) = jj;
            break;
        end
    end
end

g = @(x) f(x(1) - x(diff_to_prev_index(1) + (diff_to_prev_index(1) == 0)) * (diff_to_prev_index(1) > 0), ...
    double(conflict{1}.duration), double(optx{1}), double(rx{1}), double(mx{1}));
for ii = 2:size(conflict, 2)
    g = @(x) f(x(ii) - x(diff_to_prev_index(ii) + (diff_to_prev_index(ii) == 0)) * (diff_to_prev_index(ii) > 0), ...
        double(conflict{ii}.duration), double(optx{ii}), double(rx{ii}), double(mx{ii})) ^ d + g(x); 
end

x0 = zeros(size(conflict, 2), 1);
options = optimoptions(@fmincon,'Algorithm','sqp');
[x,fval]=fmincon(g,x0,[],[],[],[],[],[],@constraints, []);

