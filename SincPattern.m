function[] = SincPattern()
close all

[Teta, Phi] = meshgrid(-pi/2:0.02:pi/2,-pi:0.02:pi);
L1 = 2;
L2 = 4;
P1 = sin(pi*L1*Teta)./(pi*L1*Teta);
Pattern = sin(pi*L1*Teta)./(pi*L1*Teta) .*  sin(pi*L2*Phi)./(pi*L2*Phi);
TetaDeg= Teta*180/pi; PhiDeg= Phi*180/pi;


figure(1)% title(' Linear Antenna Pattern ')
surf(PhiDeg,TetaDeg,Pattern,'EdgeColor','none');
ax1 = gca;
axis([-180, 180, -90, 90, -0.5, 1]);
grid(ax1,'on')
ax1.XTick =[-180:30:180]; ax1.YTick =[-90:30:90]; ax1.ZTick = [-0.5:0.25:1]
colormap(jet)
xlabel('\phi','fontsize',12);ylabel('\theta','fontsize',12);zlabel(' Pattern [Linear]','fontsize',12);

PatternLog= 10*log10(Pattern.^2);
PatternLog(find(PatternLog <= -60)) = -60; % cut the minimum to -60 dB

figure(2) % DB
title(' Antenna Pattern ')
surf(PhiDeg,TetaDeg,PatternLog,'EdgeColor','none');
ax2 = gca;
axis([-180, 180, -90, 90, -60, 0]);
ax2.XTick =[-180:30:180]; ax2.YTick =[-90:30:90]; ax2.ZTick = [-60:10:0];
colormap(jet)
xlabel('\phi','fontsize',12);ylabel('\theta','fontsize',12);zlabel('Pattern [dB]','fontsize',12)

%%%%
L1 = 8;
L2 = 4;
P2 = sin(pi*L1*Teta)./(pi*L1*Teta);
Pattern2 = sin(pi*L1*Teta)./(pi*L1*Teta) .*  sin(pi*L2*Phi)./(pi*L2*Phi);
Pattern2Log= 10*log10(Pattern.^2);
PatternLog(find(Pattern2Log <= -60)) = -60; % cut the minimum to -60 dB

figure(3) % DB
title(' Antenna Pattern ')
surf(PhiDeg,TetaDeg,Pattern2Log,'EdgeColor','none');
ax2 = gca;
axis([-180, 180, -90, 90, -60, 0]);
ax2.XTick =[-180:30:180]; ax2.YTick =[-90:30:90]; ax2.ZTick = [-60:10:0];
colormap(jet)
xlabel('\phi','fontsize',12);ylabel('\theta','fontsize',12);zlabel('Pattern [dB]','fontsize',12)

end

