function index = convertindex(i,j,s,nL, nW, nR, h, R, H)
curve = @(r) -(H/R)*r+H;
numj = @(i) round((curve(i*h)/h)-0.5);
if s==-1
    if j==0
        index = (nL-1)*(nW-1)+i;
    else
        index = (i-1)*(nL-1) +j;
    end
    
    
else
    if j==0
        index = (nL-1)*(nW-1)+i;
    else
        index = 0;
        prev = (nL-1)*(nW-1) + nR -1;
        for k = 1: i-1
            index = index + numj(k);
        end
        index = index + j + prev;
    end
end
% if index == 35402
%     disp([i,j,s]);
% end
%if index> 24481
%    disp([i,j,s])
%end