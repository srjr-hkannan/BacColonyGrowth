function [u, t_min] = time_fd2d_main(W,L,Ks,R,H,Q,u,vol,time_hour, t, h, path)
%path = '/Users/dge/Desktop/colony/code/colony files/agent-comp/exp/middle-time/';
% L=10;
% W=20;
% R=5;
% H=2;
% h=1;
% t=1;
% cs=1000;

r=1.8*10^(-7);
lambda_s = 0.6/(3600);
Y=0.5;
rho =150;
Dp=110;
Dm=740;
q=28000;
p=16000;
Q0=Q*1000/3600;
sig = lambda_s*rho*q/Ks;
% sig = lambda_s*rho/(Y*Ks*r);
beta = Q0/(q*lambda_s);
alpha = sig*Q0/(q*lambda_s);
glu=vol*q*rho/(Ks);
nL = round(L/h-0.5);
nW = round(W/h-0.5);
nR = round(R/h-0.5);
nH = round(H/h-0.5);
curve = @(r) -(H/R)*r+H;
numj = @(i) round((curve(i*h)/h)-0.5);
numk = (nL-1)*(nW-1) + nR -1;
for i = 1:nR-1
    numk = numk +numj(i);
end
A = sparse(numk, numk);
% B = sparse(numk, numk);
mark=0;
for i = 1: nW-1
    ri = i*h;
    for j = 1:nL-1
        
        curr = convertindex(i,j,-1,nL, nW, nR, h, R, H);
        if i==1
            A(curr,convertindex(i,j-1,-1,nL, nW, nR, h, R, H)) = -Dm/(2*h^2);
            A(curr,convertindex(i+1,j,-1,nL, nW, nR, h, R, H)) = -Dm/(2*h^2)*(1+h/(2*ri));
%             B(curr,convertindex(i,j-1,-1,nL, nW, nR, h, R, H)) = Dm/(2*h^2);
%             B(curr,convertindex(i+1,j,-1,nL, nW, nR, h, R, H)) = Dm/(2*h^2);
            if j==nL-1 
                A(curr,curr) = (1/t)+(2+h/(2*ri))*Dm/(2*h^2);
%                 B(curr,curr) = (1/t)-2*Dm/(2*h^2);
            else 
                A(curr,convertindex(i,j+1,-1,nL, nW, nR, h, R, H)) = -Dm/(2*h^2);
%                 B(curr,convertindex(i,j+1,-1,nL, nW, nR, h, R, H)) = Dm/(2*h^2);
                A(curr,curr) = (1/t)+(3+h/(2*ri))*Dm/(2*h^2);
%                 B(curr,curr) = (1/t)-3*Dm/(2*h^2);
            end
        elseif i==nW-1
            A(curr,convertindex(i-1,j,-1,nL, nW, nR, h, R, H)) = -Dm/(2*h^2)*(1-h/(2*ri));
%             B(curr,convertindex(i-1,j,-1,nL, nW, nR, h, R, H)) = Dm/(2*h^2);
            if j==1
                A(curr, convertindex(i,j+1,-1,nL, nW, nR, h, R, H)) = -Dm/(2*h^2);
                A(curr, curr) = (1/t)+(2-h/(2*ri))*Dm/(2*h^2);
%                 B(curr, convertindex(i,j+1,-1,nL, nW, nR, h, R, H)) = Dm/(2*h^2);
%                 B(curr, curr) = (1/t)-2*Dm/(2*h^2);
            elseif j==nL-1
                A(curr, convertindex(i,j-1,-1,nL, nW, nR, h, R, H)) = -Dm/(2*h^2);
                A(curr, curr) = (1/t)+(2-h/(2*ri))*Dm/(2*h^2);
%                 B(curr, convertindex(i,j-1,-1,nL, nW, nR, h, R, H)) = Dm/(2*h^2);
%                 B(curr, curr) = (1/t)-2*Dm/(2*h^2);
            else
                A(curr, convertindex(i,j+1,-1,nL, nW, nR, h, R, H)) = -Dm/(2*h^2);
                A(curr, convertindex(i,j-1,-1,nL, nW, nR, h, R, H)) = -Dm/(2*h^2);
                A(curr, curr) = (1/t)+(3-h/(2*ri))*Dm/(2*h^2);
%                 B(curr, convertindex(i,j+1,-1,nL, nW, nR, h, R, H)) = Dm/(2*h^2);
%                 B(curr, convertindex(i,j-1,-1,nL, nW, nR, h, R, H)) = Dm/(2*h^2);
%                 B(curr, curr) = (1/t)-3*Dm/(2*h^2);
            end
        else
            A(curr,convertindex(i+1,j,-1,nL, nW, nR, h, R, H)) = -Dm/(2*h^2)*(1+h/(2*ri));
            A(curr,convertindex(i-1,j,-1,nL, nW, nR, h, R, H)) = -Dm/(2*h^2)*(1-h/(2*ri));
%             B(curr,convertindex(i+1,j,-1,nL, nW, nR, h, R, H)) = Dm/(2*h^2);
%             B(curr,convertindex(i-1,j,-1,nL, nW, nR, h, R, H)) = Dm/(2*h^2);
            if j==nL-1
                A(curr, convertindex(i,j-1,-1,nL, nW, nR, h, R, H)) = -Dm/(2*h^2);
                A(curr, curr) = (1/t)+3*Dm/(2*h^2);
%                 B(curr, convertindex(i,j-1,-1,nL, nW, nR, h, R, H)) = Dm/(2*h^2);
%                 B(curr, curr) = (1/t)-3*Dm/(2*h^2);
            elseif j==1 && i>=nR
                A(curr, convertindex(i,j+1,-1,nL, nW, nR, h, R, H)) = -Dm/(2*h^2);
                A(curr, curr) = (1/t)+3*Dm/(2*h^2);
%                 B(curr, convertindex(i,j+1,-1,nL, nW, nR, h, R, H)) = Dm/(2*h^2);
%                 B(curr, curr) = (1/t)-3*Dm/(2*h^2);
%                 disp([curr, i]);
            else
                A(curr, convertindex(i,j-1,-1,nL, nW, nR, h, R, H)) = -Dm/(2*h^2);
                A(curr, convertindex(i,j+1,-1,nL, nW, nR, h, R, H)) = -Dm/(2*h^2);
                A(curr, curr) = (1/t)+4*Dm/(2*h^2);
%                 B(curr, convertindex(i,j-1,-1,nL, nW, nR, h, R, H)) = Dm/(2*h^2);
%                 B(curr, convertindex(i,j+1,-1,nL, nW, nR, h, R, H)) = Dm/(2*h^2);
%                 B(curr, curr) = (1/t)-4*Dm/(2*h^2);
            end
        end
        
        
    end
end
reached = 0;
for i = 1:nR-1
    if numj(i)>=1
        curr = convertindex(i,0,-1,nL, nW, nR, h, R, H);
        A(curr,convertindex(i,1,1,nL, nW, nR, h, R, H)) = Dp;
        A(curr,convertindex(i,1,-1,nL, nW, nR, h, R, H)) = Dm;
        A(curr, curr) = -(Dp+Dm);
    else
        if reached ==0
            reached =i-1;
        end
         
        curr = convertindex(i,0,-1,nL, nW, nR, h, R, H);
        A(curr,convertindex(i,1,-1,nL, nW, nR, h, R, H)) = Dm;
%         [i1,j1] = reflex(i*h,h,R,H, h);
        reflexpos = convertindex(reached,1,1,nL, nW, nR, h, R, H);
        A(curr, curr) = -(Dp+Dm);
        A(curr, reflexpos) = A(curr, reflexpos) + Dp;
    end
end
for j = 1:numj(1)
    i = 1;
    ri=h;
    curr = convertindex(i,j,1,nL, nW, nR, h, R, H);
    A(curr,convertindex(i,j-1,1,nL, nW, nR, h, R, H)) = -Dp/(h^2);
    if j==numj(1) 
        [i1, j1] = reflex(i*h,(j+1)*h,R,H, h);
        reflexpos = convertindex(i1,j1,1,nL, nW, nR, h, R, H);
        A(curr, reflexpos) = A(curr, reflexpos) -Dp/(h^2);
        A(curr,curr) = (1/t)+((3+h/(2*ri))*Dp)/(h^2);
    else 
        A(curr,convertindex(i,j+1,1,nL, nW, nR, h, R, H)) = -Dp/(h^2);
        A(curr,curr) = (1/t)+((3+h/(2*ri))*Dp)/(h^2);
    end
    if curve((i+1)*h)<=j*h
        [i1, j1] = reflex((i+1)*h,j*h,R,H,h);
        reflexpos = convertindex(i1,j1,1,nL, nW, nR, h, R, H);
        A(curr, reflexpos) = A(curr, reflexpos) -(1 + h/(2*ri))*Dp/(h^2) ;
    else
           
       pos = convertindex(i+1,j,1,nL, nW, nR, h, R, H);
       A(curr,pos) = A(curr,pos) -(1 + h/(2*ri))*Dp/(h^2) ;
    end
end
for i = 2:nR-1
    ri=i*h;
    for j = 1: numj(i)
%         disp(i);
        curr = convertindex(i,j,1,nL, nW, nR, h, R, H);
        A(curr,convertindex(i-1,j,1,nL, nW, nR, h, R, H)) = -Dp/(h^2)*(1 - h/(2*ri));
        A(curr, convertindex(i,j-1,1,nL, nW, nR, h, R, H)) = -Dp/(h^2);
        A(curr, curr) = (1/t)+(4*Dp)/(h^2);
        if curve(i*h)<=(j+1)*h
            [i1, j1] = reflex(i*h,(j+1)*h,R,H, h);
%               disp([i,j+1,i1,j1]);
            reflexpos = convertindex(i1,j1,1,nL, nW, nR, h, R, H);
            A(curr, reflexpos) = A(curr, reflexpos) -Dp/(h^2);
%             if curr == reflexpos
%                 disp('y');
%             end
        else
            pos = convertindex(i,j+1,1,nL, nW, nR, h, R, H);
            A(curr, pos) = A(curr, pos) -Dp/(h^2);
        end
        if curve((i+1)*h)<=j*h
            [i1, j1] = reflex((i+1)*h,j*h,R,H,h);
%              disp([i+1,j,i1,j1]);
            reflexpos = convertindex(i1,j1,1,nL, nW, nR, h, R, H);
            A(curr, reflexpos) = A(curr, reflexpos) -Dp/(h^2)*(1 + h/(2*ri)) ;
        else
           
            pos = convertindex(i+1,j,1,nL, nW, nR, h, R, H);
            A(curr,pos) = A(curr,pos) -Dp/(h^2)*(1 + h/(2*ri)) ;
        end
    end
end
N = (nL-1)*(nW-1) + nR -1;
% u  = zeros(numk, 1);
% for i = 1:N
%     u(i) = cs;
% end
% for i= (nL-1)*(nW-1)+1:N
%     u(i) = cs*Dm/(Dp+Dm);
% end


% Data_A = load('/Users/dge/Desktop/colony/code/colony files/Harish data/Data_1+1/Cs_1000/Data_A.mat');
% M12 = squeeze(Data_A.Glucose(10,:,:));
% u_h = mat_to_u(M12, nW,nL);
% total_glu_h = r*(1000*W*L-totcl_c(u_h,nW,nL,h));
% hr=10;

% cs_r=1000/cs;
% sa = '/Users/dge/Desktop/colony/code/colony files/time/u_hr'+ string(hr)+' Ks= '+string(Ks)+'.mat';
% uu=load(sa);
% u=uu.u;
growth_glu=0;
%T = time_min*60000;
t_min= 0;
[L_M,U,P] = lu(A);
n_reached = true;
while growth_glu< glu || n_reached
    t_min = t_min + 1;
    if growth_glu>= glu
        n_reached = false;
        sa = string(path)+'cs1000_c_hr_'+string(time_hour)+ ' 2du_iter'+ string(t_min)+' t= '+string(t)+' h='+string(h)+'.mat';
        save(sa,'u');
    end

    b = zeros(numk,1);
    b_M = zeros(nW-1,nL-1);
    u_M = zeros(nW-1,nL-1);
    i=1;
    j=1;
    for curr = 1: (nW-1)*(nL-1)
        u_M(i,j) = u(curr);
        if rem(curr,nL-1)==0
           i=i+1;
           j=1;
       else
           j=j+1;
       end
    end
    for i=1:nW-1
        ri = i*h;
        for j=1:nL-1
            if j==1 && i<nR
                if i==1
                    b_M(i,j)= u(convertindex(i,j-1,-1,nL, nW, nR, h, R, H)) *Dm/(2*h^2);
                    b_M(i,j)= b_M(i,j) + u_M(i+1,j) *Dm/(2*h^2)*(1+h/(2*ri));
                    b_M(i,j)= b_M(i,j)+u_M(i,j+1) *Dm/(2*h^2);
                    b_M(i,j)= b_M(i,j)+u_M(i,j)*((1/t)-(3+h/(2*ri))*Dm/(2*h^2));
                else
                    b_M(i,j)= u_M(i+1,j) *(1+h/(2*ri))*Dm/(2*h^2);
                    b_M(i,j)= b_M(i,j) +u_M(i-1,j) *(1-h/(2*ri))*Dm/(2*h^2);
                    b_M(i,j)= b_M(i,j)+u(convertindex(i,j-1,-1,nL, nW, nR, h, R, H)) *Dm/(2*h^2);
                    b_M(i,j)= b_M(i,j)+u_M(i,j+1) *Dm/(2*h^2);
                    b_M(i,j)= b_M(i,j)+u_M(i,j)*((1/t)-4*Dm/(2*h^2));
                end
                continue
            end
            if i==1
                b_M(i,j)= u_M(i,j-1) *Dm/(2*h^2);
                b_M(i,j)= b_M(i,j) + u_M(i+1,j) *Dm/(2*h^2)*(1+h/(2*ri));
                if j==nL-1 
                    b_M(i,j)= b_M(i,j)+u_M(i,j)*((1/t)-(2+h/(2*ri))*Dm/(2*h^2));
                else 
                    b_M(i,j)= b_M(i,j)+u_M(i,j+1) *Dm/(2*h^2);
                    b_M(i,j)= b_M(i,j)+u_M(i,j)*((1/t)-(3+h/(2*ri))*Dm/(2*h^2));
                end
            elseif i==nW-1
                b_M(i,j)= u_M(i-1,j) *Dm/(2*h^2)*(1-h/(2*ri));
                if j==1
                    b_M(i,j)= b_M(i,j)+u_M(i,j+1) *Dm/(2*h^2);
                    b_M(i,j)= b_M(i,j)+u_M(i,j)*((1/t)-(2-h/(2*ri))*Dm/(2*h^2));
                elseif j==nL-1
                    b_M(i,j)= b_M(i,j)+u_M(i,j-1) *Dm/(2*h^2);
                    b_M(i,j)= b_M(i,j)+u_M(i,j)*((1/t)-(2-h/(2*ri))*Dm/(2*h^2));
                else
                    b_M(i,j)= b_M(i,j)+u_M(i,j+1) *Dm/(2*h^2);
                    b_M(i,j)= b_M(i,j)+u_M(i,j-1) *Dm/(2*h^2);
                    b_M(i,j)= b_M(i,j)+u_M(i,j)*((1/t)-(3-h/(2*ri))*Dm/(2*h^2));
                end
            else
                b_M(i,j)= u_M(i+1,j) *(1+h/(2*ri))*Dm/(2*h^2);
                b_M(i,j)= b_M(i,j) +u_M(i-1,j) *(1-h/(2*ri))*Dm/(2*h^2);
                if j==nL-1
                    b_M(i,j)= b_M(i,j)+u_M(i,j-1) *Dm/(2*h^2);
                    b_M(i,j)= b_M(i,j)+u_M(i,j)*((1/t)-3*Dm/(2*h^2));
                elseif j==1 && i>=nR
                    b_M(i,j)= b_M(i,j)+u_M(i,j+1) *Dm/(2*h^2);
                    b_M(i,j)= b_M(i,j)+u_M(i,j)*((1/t)-3*Dm/(2*h^2));
                else
                    b_M(i,j)= b_M(i,j)+u_M(i,j-1) *Dm/(2*h^2);
                    b_M(i,j)= b_M(i,j)+u_M(i,j+1) *Dm/(2*h^2);
                    b_M(i,j)= b_M(i,j)+u_M(i,j)*((1/t)-4*Dm/(2*h^2));
                end
            end
       
       end
    end

    i = 1;
    j = 1;
    for curr = 1: (nW-1)*(nL-1)
        b(curr) = b_M(i,j);
        if rem(curr,nL-1)==0
           i=i+1;
           j=1;
       else
           j=j+1;
       end
    end


    % Acetate calculation
%     for i = 1: nW-1
%         ri=i*h;
%         for j = 1:nL-1
%         
%             curr = convertindex(i,j,-1,nL, nW, nR, h, R, H);
%             if i==1
%                 b_A(curr)= b_A(curr)+u_A(convertindex(i,j-1,-1,nL, nW, nR, h, R, H)) *Dm/(2*h^2);
%                 b_A(curr)= b_A(curr)+u_A(convertindex(i+1,j,-1,nL, nW, nR, h, R, H)) *Dm/(2*h^2)*(1+h/(2*ri));
%                 if j==nL-1 
%                     b_A(curr)= b_A(curr)+u_A(curr)*((1/t)-(2+h/(2*ri))*Dm/(2*h^2));
%                 else 
%                     b_A(curr)= b_A(curr)+u_A(convertindex(i,j+1,-1,nL, nW, nR, h, R, H)) *Dm/(2*h^2);
%                     b_A(curr)= b_A(curr)+u_A(curr)*((1/t)-(3+h/(2*ri))*Dm/(2*h^2));
%                 end
%             elseif i==nW-1
%                 b_A(curr)= b_A(curr)+u_A(convertindex(i-1,j,-1,nL, nW, nR, h, R, H)) *Dm/(2*h^2)*(1-h/(2*ri));
%                 if j==1
%                     b_A(curr)= b_A(curr)+u_A(convertindex(i,j+1,-1,nL, nW, nR, h, R, H)) *Dm/(2*h^2);
%                     b_A(curr)= b_A(curr)+u_A(curr)*((1/t)-(2-h/(2*ri))*Dm/(2*h^2));
%                 elseif j==nL-1
%                     b_A(curr)= b_A(curr)+u_A(convertindex(i,j-1,-1,nL, nW, nR, h, R, H)) *Dm/(2*h^2);
%                     b_A(curr)= b_A(curr)+u_A(curr)*((1/t)-(2-h/(2*ri))*Dm/(2*h^2));
%                 else
%                     b_A(curr)= b_A(curr)+u_A(convertindex(i,j+1,-1,nL, nW, nR, h, R, H)) *Dm/(2*h^2);
%                     b_A(curr)= b_A(curr)+u_A(convertindex(i,j-1,-1,nL, nW, nR, h, R, H)) *Dm/(2*h^2);
%                     b_A(curr)= b_A(curr)+u_A(curr)*((1/t)-(3-h/(2*ri))*Dm/(2*h^2));
%                 end
%             else
%                 b_A(curr)= b_A(curr)+u_A(convertindex(i+1,j,-1,nL, nW, nR, h, R, H)) *(1+h/(2*ri))*Dm/(2*h^2);
%                 b_A(curr)= b_A(curr)+u_A(convertindex(i-1,j,-1,nL, nW, nR, h, R, H)) *(1-h/(2*ri))*Dm/(2*h^2);
%                 if j==nL-1
%                     b_A(curr)= b_A(curr)+u_A(convertindex(i,j-1,-1,nL, nW, nR, h, R, H)) *Dm/(2*h^2);
%                     b_A(curr)= b_A(curr)+u_A(curr)*((1/t)-3*Dm/(2*h^2));
%                 elseif j==1 && i>=nR
%                     b_A(curr)= b_A(curr)+u_A(convertindex(i,j+1,-1,nL, nW, nR, h, R, H)) *Dm/(2*h^2);
%                     b_A(curr)= b_A(curr)+u_A(curr)*((1/t)-3*Dm/(2*h^2));
%                 else
%                     b_A(curr)= b_A(curr)+u_A(convertindex(i,j-1,-1,nL, nW, nR, h, R, H)) *Dm/(2*h^2);
%                     b_A(curr)= b_A(curr)+u_A(convertindex(i,j+1,-1,nL, nW, nR, h, R, H)) *Dm/(2*h^2);
%                     b_A(curr)= b_A(curr)+u_A(curr)*((1/t)-4*Dm/(2*h^2));
%                 end
%             end
%         end
%     end
    maintain = 0;
    for i=N+1:numk
        ui=u(i);
%         uAi=u_A(i);
        m=1;
        if ui/(beta)<1
            m=ui/beta;
        end
        l = lambda_c(lambda_s, Q0, q,  ui);
        b(i) = (1/t)*ui-q*l*((ui)/(ui+1))*rho-alpha*m;
%         b_A(i) = (1/t)*uAi+p*l*((ui)/(ui+1))*rho/Ks;
    end
    pb = P*b;
    y = L_M\pb;
    u = U\y;
    for i=N+1:numk
        if u(i) < 0
%             disp(i)
%             disp(u(i))
             u(i) = 0;
            
        end
    end
%     u_A = A\b_A;

%   maintain
    for i = 1:nR-1
        for j = 1:numj(i)
            curr = convertindex(i,j,1,nL, nW, nR, h, R, H);
            ui=u(curr);
            m=1;
            if ui/(beta)<1
                m=ui/beta;
            end
            if j==0
                maintain = maintain+alpha*m*t*(h^3)*pi*(2*i-1)/2;
            else
                maintain = maintain+alpha*m*t*(h^3)*pi*(2*i-1);
            end
        end
    end


%     if total_glu>total_glu_h
%         sa = '/Users/dge/Desktop/colony/code/colony files/time/u_hr'+ string(hr)+' Ks= '+string(Ks)+'.mat';
%         save(sa,'u');
%         hr=hr+1;
%         M12 = squeeze(Data_A.Glucose(hr,:,:));
%         u_h = mat_to_u(M12, nW,nL);
%         total_glu_h = r*(1000*W*L-totcl_c(u_h,nW,nL,h));
%     end
    grow =0;

    for i = 1:nR-1
        r = t*Dm;
        interface_index = convertindex(i,0,-1,nL, nW, nR, h, R, H);
        index = convertindex(i,1,-1,nL, nW, nR, h, R, H);
        grow = grow + (2*i*h*pi)*r*(u(index) - u(interface_index));
        if i==1
            grow = grow + (1/4)*r*h*pi*(u(index) - u(interface_index));
        end
        if i==nR-1
            grow = grow + r*(R^2 - ((nR-0.5)*h)^2)*pi*(u(index) - u(interface_index))/h;
        end
    end
%     mark = mark + grow - maintain;
    growth_glu = growth_glu + grow-maintain;
    if rem(t_min,125*5)==0
        fprintf(string(t_min)+'\n')
        fprintf(string(u(1))+'  '+string(u(N-nR+1))+'\n')
        fprintf(string(growth_glu/glu)+'\n')
        if rem(t_min,625*5)==0
            sa = string(path)+'cs1000_c_hr_'+string(time_hour)+ ' 2du_iter'+ string(t_min)+' t= '+string(t)+' h='+string(h)+'.mat';
            save(sa,'u');
        end
    end

end
sa = string(path)+'cs1000_c_hr_'+string(time_hour)+ ' 2du_'+ string(2)+' t= '+string(t)+' h='+string(h)+'.mat';
save(sa,'u');
