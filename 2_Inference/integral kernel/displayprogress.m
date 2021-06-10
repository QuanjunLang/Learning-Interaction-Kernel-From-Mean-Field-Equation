
%------------------- display progress -------------------------------------
%{
- put this function in your large for-loop folder; then
- 1. before the for-loop 
   fprintf('Progress: '); reverseStr = []; 
- 2. at the beginnng of for i=1:end: 
   reverseStr = displayprogress(100*i/(allsteps-1), reverseStr); 
%}

function reverseStr = displayprogress(perc,reverseStr)
msg = sprintf('%3.1f', perc);
fprintf([reverseStr, msg, '%%']);
reverseStr = repmat(sprintf('\b'), 1, length(msg)+1);
end
