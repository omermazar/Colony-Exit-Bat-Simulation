function [] = my_swarm_xyplot(swarm_strcut_full)

figure
hold on

for BatNum = 1:swarm_strcut_full.num_of_bats
   plot(swarm_strcut_full.all_bats.BAT(BatNum).xBat, ...
       swarm_strcut_full.all_bats.BAT(BatNum).yBat,'.') 
   plot(swarm_strcut_full.all_bats.BAT(BatNum).start_pos(1), ...
       swarm_strcut_full.all_bats.BAT(BatNum).start_pos(2), 'sr')
end
xlabel('xPosition [m]')
ylabel('yPosition [m]')
title('Bats Flight')
