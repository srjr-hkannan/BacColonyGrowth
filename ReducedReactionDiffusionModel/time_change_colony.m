% Data_A = load('/Users/dge/Desktop/colony/code/colony files/Harish data/Data_1+1/Cs_500/Data_C.mat');
% radius=Data_A.Radius_list;
% height=Data_A.Height_list;

time = 15:5:60;
radius = zeros(10,1);
height = zeros(10,1);
path = '/Users/dge/Desktop/colony/Quasi-static codes/';
for i =1:10
    radius(i) = 27*time(i)-5;
    height(i) = 100 +time(i)*1.24;
end
index=1;
R = radius(index);
H = height(index);

W=512*5;
L=512*2;
Ks=20;
Q=10;
h=16;
t=0.1;
cs=500;
nL = round(L/h);
nW = round(W/h);
nR = round(R/h);
nH = round(H/h);
curve = @(r) -(H/R)*r+H;
numj = @(i) round((curve(i*h)/h)-0.5);
numk = (nL-1)*(nW-1) + nR -1;
for i = 1:nR-1
    numk = numk +numj(i);
end
u  = zeros(numk, 1);
% u_A  = zeros(numk, 1);
N = (nW-1)*(nL-1);
for i = 1:N
    u(i) = cs;
end
for i= (nL-1)*(nW-1)+1:N
    u(i) = cs*Dm/(Dp+Dm);
end
% sa = '/Users/dge/Desktop/colony/code/colony files/agent-comp/Deep-L/cs500_c_hr_'+string(5)+ ' 2du_'+ string(2)+' t= '+string(t)+' h='+string(h)+'.mat';
% uu = load(sa);
% u=uu.u;
u_1dp=u;
u_2dp=u;
% uA_1dp = u_A;
% uA_2dp = u_A;
t1=0;
t2=0;

for index= 1:10
    R =radius(index);
    H = height(index);
    nL = round(L/h);
    nW = round(W/h);
    nR = round(R/h);
    curve = @(r) -(H/R)*r+H;
    numj = @(i) round((curve(i*h)/h)-0.5);
    numk = (nL-1)*(nW-1) + nR -1;
    for i = 1:nR-1
        numk = numk +numj(i);
    end

    if index == 1
        vol_1d = (radius(index)*height(index))/2;
        vol_2d = ((radius(index)^2)*height(index))*pi/3;
        [u_1d, t_1d] = time_fd1d_main(W,L,Ks,R,H,Q, u_1dp, vol_1d, time(index), t, h, 1, path);
        [u_2d, t_2d] = time_fd2d_main(W,L,Ks,R,H,Q, u_2dp, vol_2d, time(index), t, h, path);

    else
        vol_1d = (R*H-radius(index-1)*height(index-1))/2;
        vol_2d = ((R^2)*H-(radius(index-1)^2)*height(index-1))*pi/3;
        [u_1d, t_1d] = time_fd1d_main(W,L,Ks,R,H,Q, u_1dp, vol_1d, time(index), t, h, 1, path);
        [u_2d, t_2d] = time_fd2d_main(W,L,Ks,R,H,Q, u_2dp, vol_2d, time(index), t, h, path);
    end
%     t1=t1+t_1d;
    t2=t2+t_2d;
    nRp=nR;
    Rp=R;
    Hp=H;
    curveP=curve;
    numjp = numj;
    R =radius(index+1);
    H = height(index+1);
    nR = round(R/h);
    nH = round(H/h);
    curve = @(r) -(H/R)*r+H;
    numj = @(i) round((curve(i*h)/h)-0.5);
    numk = (nL-1)*(nW-1) + nR -1;
    for i = 1:nR-1
        numk = numk +numj(i);
    end
%     u_1dp = zeros(numk, 1);
    u_2dp = zeros(numk, 1);
%     uA_1dp = zeros(numk, 1);
%     uA_2dp = zeros(numk, 1);
    Np= (nL-1)*(nW-1) + nRp -1;
    for i=1:Np
%         u_1dp(i) = u_1d(i);
        u_2dp(i) = u_2d(i);
%         uA_1dp(i) = uA_1d(i);
%         uA_2dp(i) = uA_2d(i);
    end
    for i = 1:nRp-1
        for j = 1 :numjp(i)
            curr= convertindex(i,j,-1,nL, nW, nRp, h, Rp, Hp);
            prev= convertindex(i,j,-1,nL, nW, nR, h, R, H);
%             u_1dp(prev) = u_1d(curr);
            u_2dp(prev) = u_2d(curr);
%             uA_1dp(prev) = uA_1d(curr);
%             uA_2dp(prev) = uA_2d(curr);
        end
    end
end
