scp_so_finitewell;

Pe123 = Pe;
clear Pe;
load fin_well_orig;
plot(E, real(Pe));
hold on;
plot(E, real(Pe123),'r--');
