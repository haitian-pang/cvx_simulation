clear;
%% Solve the optimization problem using cvx
q = 5;
Q = [q q q q q q q q q q q q q q q q q q q q q q q q q q q q q q q q q q q q q q q q q q q q q q q q q q];  % data quota of users
%pred_watch = rem(randn(N,T),100);  % prediction result of users' watching videos 
%% 2. Content Parameters
B = 1100;  % The budget of CP for sponsoring content
%C = 1000;        % The number of total contents
U = 50;
T = 14;
C = 5;
L = 55;
Bb = 0;
%S = 5;
Cm = load(['cmat.csv']);
Lm = load(['lmat.csv']);
%P = zeros(C,L);
%P = load(['pmat5.csv']);
%P = load(['pmat.csv']);
%P = load(['pmat30.csv']);
%P = load(['pmat50.csv']);
%P = load(['pmat54.csv']);
%P = load(['pmat100.csv']);
%P = ones(C,L);

S = ones(C,1);
Sc = 3*S;
G = [5,4,3,6,7] ;

cvx_begin quiet 
    variable x(U, T) integer  % x(t, k): 
    variable y(U, T) integer
	variable z(U, T) integer
    variable P(C, L) integer
	%%
	for k = 1:U
		for i = 1:T
			sponsor(k,i) = z(k,i) + P(Cm(k,i),Lm(k,i));
		end
	end
	%%
	ccache = 0;
	for c = 1:C
		for l = 1:L
			ccache = ccache + P(c,l) * Sc(c,1);
		end
    end
     ccellular = x(1,1)*0;
	for k = 1:U
		for i = 1:T
			ccellular = ccellular + x(k,i) * S(Cm(k,i),1);
		end
	end
	%%
    ubudget = x(1,1)*zeros(50,1);
    for k = 1:U
        for i = 1:T
            ubudget(k,1) = ubudget(k,1) + (z(k,i)-x(k,i)) * S(Cm(k,i),1); %(y(k,i) - x(k,i) - P(Cm(k,i),Lm(k,i)))
        end
    end
	%%
    ubudgetr = x(1,1)*zeros(50,14);
	for k = 1:U
		for i = 1:T
            ubudgetr(k,i) = Q(k);
            if(i > 1)
                for r = 1:(i-1)
                    ubudgetr(k,i) = ubudgetr(k,i) - S(Cm(k,r))*(z(k,r) - x(k,r));
                end
            end
			ubudgetr(k,i) = ubudgetr(k,i)/Q(k);
		end
	end
	%%
    bandcap =  x(1,1)*zeros(14,55);
	for i = 1:T
		for k = 1:U
			bandcap(i,Lm(k,i)) = bandcap(i,Lm(k,i)) + (y(k,i) - z(k,i) - P(Cm(k,i),Lm(k,i)))*S(Cm(k,i));
		end
	end
	%%
	
	%%
    subject to
        1 >= y >= z >= x >= 0;
		1 >= sponsor >= 0;
        1 >= P >= 0;
		0<= ccellular + ccache <= B;
		0 <= ubudget(:,1) <= Q';
		y >= ubudgetr >= 0;
		0 <= bandcap <= Bb;
     
    g1 = x(k,i)*0;
	for k = 1:U
        for i = 1:T
			g1 = g1 + G(Cm(k,i))*y(k,i) - S(Cm(k,i))*x(k,i);
		end
    end
    g2 = x(1,1)*0;
    for c = 1:C
        for l = 1:L
            g2 = g2 - P(c,l) * Sc(c);
        end
    end
	maximize(g1 + g2)   
cvx_end


%g1 + g2

%% plot the result
figure;
plot(sum(y,2),'.:r');hold on;
plot(sum(z,2), 'o-g');
plot(sum(x,2), '<:b');
plot(sum(y,2)-sum(z,2), 'v:k');
legend('serve','cellular' ,'cellular sponsor','cache');

fprintf('\n -- Total welfare = %0.6f, cellular cost = %d, cache cost = %d --\n', g1+g2, ccellular, ccache); 

%figure;
%plot(ubudget,'.:r');
%figure;
%plot([2800, 2920, 2985,3195 ,3285, 2675 ],'.:g');