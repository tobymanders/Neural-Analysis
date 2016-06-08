
difference = diff(board_dig_in_data);
test = find (difference==1);
time = test/30000;
time = transpose(time);
xlswrite ('laserOn.xlsx',time);