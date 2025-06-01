function[]  = TestBatAnimate()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
AnimateSpeed = 1;
h1 = animatedline('Marker','o');
h2 = animatedline('Marker','+', 'color','r');
t=1;
x = 1:5000;
y = 1:5000;
z = max(x)*sin(0.01*x);
axis([0, max(x), 0, max(y)]);

for k =1:AnimateSpeed:length(x)-AnimateSpeed
    
    k2= max(1, k-5*AnimateSpeed);
    addpoints(h1,x(k2:k+AnimateSpeed-1),y(k2:k+AnimateSpeed-1))
    addpoints(h2,x(k:k+AnimateSpeed-1),z(k:k+AnimateSpeed-1))
    drawnow %limitrate
    clearpoints(h1)
% % %     DrawFlag = ~mod(k,100);
% % %     if DrawFlag
% % %         drawnow
% % %     end
end
end