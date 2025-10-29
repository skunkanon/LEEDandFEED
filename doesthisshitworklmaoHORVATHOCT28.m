%OCT 28, RUN HORVARTHSCAN FIRST 
Mg_timeSPAN = Mg_scan(:,1);
Mg_rateSPAN = Mg_scan(:,2);
Mg_tempSPAN = Mg_scan(:,3);

figure;
plot(Mg_tempSPAN, Mg_rateSPAN);

%%
figure;
plot(Mg_timeSPAN, Mg_rateSPAN);