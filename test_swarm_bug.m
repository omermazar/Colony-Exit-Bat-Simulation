nbats = curr_strct.num_of_bats; 
all_bats_num = curr_strct.AllParams.SimParams.TotalBatsNumber;
all_x = [curr_strct.all_bats.BAT.xBat];
all_y = [curr_strct.all_bats.BAT.yBat];

all_x = reshape(all_x, [], nbats)';
all_y = reshape(all_y, [], nbats)';

figure
plot(all_x', all_y','.')

all_teta0 = [curr_strct.all_bats.BAT.start_angle];

start_xy = [curr_strct.all_bats.BAT.start_pos];
start_xy = reshape(start_xy,2,[])';
start_x = start_xy(:,1);
start_y = start_xy(:,2);

histogram(start_x,100)
histogram(start_y,100)

x= [1 2 1 2 1 2]
xy = reshape(x,2,[])'


histogram(all_teta0)