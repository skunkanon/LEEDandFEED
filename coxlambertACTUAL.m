%% 4/11: Cleaning up first attempt. First, making a plot of coverage vs temperature. 

x_0p8_raw = fig5_0p8(:,1);
y_0p8_raw = fig5_0p8(:,2);

scatter(x_0p8_raw, y_0p8_raw);

bg_pre_raw = y_0p8_raw(1:51);
bg_post_raw = y_0p8_raw(85:100);


bg_pre = mean(bg_pre_raw);
bg_post = mean(bg_post_raw);

