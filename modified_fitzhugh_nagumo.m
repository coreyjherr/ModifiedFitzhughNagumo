%% Initialize variables
clear all

a1 = .167;
a2 = 16.67;
a3 = 167;
a4 = 1.44;
a5 = 1.47;
a = [a1 a2 a3 a4 a5];

c1 = .1;
c2 = 3.9;
c = [c1 c2];

Du = .1;
Dv = 1;
eps = .52;

dx = .05;
N = 1000;
L = 600;
parmas.L = L;
r = 0;

u0 = zeros(L,1);
v0 = 2*ones(L,1);
y0 = [u0;v0];
t0 = 0;
t_wave = 5;
t_final = 30;


%% Integrate dynamical equations

wrapper = @(t,y)fitzhugh_nagumo(y, r, a, c, Du, Dv, dx, eps, L);
[ts1, y1] = ode45(wrapper, [t0, t_wave], y0);

y00 = y1(end,:);
y00(L/2) = 1;
[ts2, y2] = ode45(wrapper, [t_wave, t_final], y00);

%% Plot kymograph
us = cat(1, y1(:,1:L), y2(:,1:L));
vs = cat(1, y1(:,L+1:end), y2(:,L+1:end));
ys = cat(1, y1, y2);

figure(1)
imagesc(us')

%% Plot nullclines
x = 240;

for t = 3000:100:20000

us = -1:.01:4;
vs = 0:.01:6;

d2u = lap(ys(t,1:L)')/dx^2;
d2v = lap(ys(t,L+1:end)')/dx^2;

inhib = 1/c(1)*(Dv/eps*d2v(x)+c(2)*us);
activ = r + 1/a(2)*(Du*d2u(x)./us-a(1)+a(3)*us.^2./(us.*(a(4)+us.^2))+a(5));

figure(2)
plot(us, activ)
hold on
plot(us, inhib)
plot(ys(t,x), ys(t, L+x), 'b*')
hold off
xlim([-1,4])
ylim([-1,6])
drawnow
pause(.03)
end

%% Functions to evolve equation

function next_state = fitzhugh_nagumo(current_state, r, a, c, Du, Dv, dx, eps, L)

u = current_state(1:L);
v = current_state(L+1:end);
next_state = zeros(size(current_state));

next_state(1:L) = Du*lap(u)/dx^2 - a(1)*u - a(2)*u.*(v-r) + a(3)*u.^2./(a(4)+u.^2) + a(5);
next_state(L+1:end) = Dv*lap(v)/dx^2 + eps*(-c(1)*v+c(2)*u);
end


function out = lap(x)
out = cat(1,x(2:end), x(1)) - 2*x + cat(1,x(end), x(1:end-1));
end
