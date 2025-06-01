function[] = SincPattern()

Teta = [-pi/2, 0.01, i/2];
Phi = [-pi, 0.01, pi];
L1 = 2;
L2 = 4;
Pattern = sin(pi*L1*Teta)/(pi*L1*Teta) .*  sin(pi*L2*Phi)/(pi*L2*Phi);
surf(Teta,Phi,Patten);
colormap(jet)
