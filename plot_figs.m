clear all
close all

load('model_generated_data/sim_runs_20240128');
load('model_generated_data/PACE_blood_summary_20240109');
load('model_generated_data/sim_runs_20240128_organ_sum');
load('model_generated_data/PACE_FACS_summary_time_20240109');
load('model_generated_data/PACE_FACS_summary_20240109');

color = flip({'r' 'm' 'g' 'c' 'b'});
dose_rel=flip([0.1 0.5 2 2.5 3.5]);

ind = 1:12000;
times = ind/12000*48;
axis_fs = 15;
tick_fs=12;

figure('units','normalized','outerposition',[0 0 0.4 1])

%s=subplot(4,2,1);
s=figure()
for i=1:length(dose_rel)
blood_time=blood_summary.time(blood_summary.dose==dose_rel(i));
blood_mean=blood_summary.mean(blood_summary.dose==dose_rel(i));
blood_sd=blood_summary.sd(blood_summary.dose==dose_rel(i));

load(['model_generated_data/sim_runs_ ',num2str(dose_rel(i)),' _20240118']);

plot(times(ind),OP.Blood_mn(ind),[char(color(i)), '-'],'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');hold on;
plot(times(ind),OP.Blood_mn(ind)+OP.Blood_sd(ind)/10,[char(color(i)), '--'],'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');
plot(times(ind),OP.Blood_mn(ind)-OP.Blood_sd(ind)/10,[char(color(i)), '--'],'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');
plot(blood_time,blood_mean,[char(color(i)), '.'],'MarkerSize',20);
errorbar(blood_time,blood_mean,blood_sd/sqrt(3),[char(color(i)), '.'],'LineWidth',2);

end
hold off;
set(gca, 'box','off','FontSize',tick_fs,'XTick',0:8:48,'YTick',0:500:2000);
xlim([0 48]);
ylim([0 2000]);
xlabel("Time (hr)",'FontSize',axis_fs);
ylabel("NP Conc. (\mu{g}/ml)",'FontSize',axis_fs);
axis square
title('Blood');

   saveas(gcf,'figures/tot_sims_blood.fig')
   saveas(gcf,'figures/tot_sims_blood','epsc')
   saveas(gcf,'figures/tot_sims_blood','jpeg')

%s=subplot(4,2,3);
s=figure()
for i=1:length(dose_rel)
load(['model_generated_data/sim_runs_ ',num2str(dose_rel(i)),' _20240118']);
plot(times(ind),OP.Liver_mn(ind),[char(color(i)), '-'],'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');hold on;
plot(times(ind),OP.Liver_mn(ind)+OP.Liver_sd(ind)/10,[char(color(i)), '--'],'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');
plot(times(ind),OP.Liver_mn(ind)-OP.Liver_sd(ind)/10,[char(color(i)), '--'],'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');
end
hold off;
set(gca, 'box','off','FontSize',tick_fs,'XTick',0:8:48,'YTick',0:50:200);
xlim([0 48]);
ylim([0 200]);
xlabel("Time (hr)",'FontSize',axis_fs);
ylabel("NP Conc. (\mu{g}/ml)",'FontSize',axis_fs);
axis square
title('Liver'); %s,

   saveas(gcf,'figures/tot_sims_liver.fig')
   saveas(gcf,'figures/tot_sims_liver','epsc')
   saveas(gcf,'figures/tot_sims_liver','jpeg')

%s=subplot(4,2,4);
s=figure()
for i=1:length(dose_rel)
load(['model_generated_data/sim_runs_ ',num2str(dose_rel(i)),' _20240118']);
plot(times(ind),OP.Spleen_mn(ind),[char(color(i)), '-'],'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');hold on;
plot(times(ind),OP.Spleen_mn(ind)+OP.Spleen_sd(ind)/10,[char(color(i)), '--'],'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');
plot(times(ind),OP.Spleen_mn(ind)-OP.Spleen_sd(ind)/10,[char(color(i)), '--'],'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');
end
hold off;
set(gca, 'box','off','FontSize',tick_fs,'XTick',0:8:48,'YTick',0:50:200);
xlim([0 48]);
ylim([0 200]);
xlabel("Time (hr)",'FontSize',axis_fs);
ylabel("NP Conc. (\mu{g}/ml)",'FontSize',axis_fs);
axis square
title('Spleen');

   saveas(gcf,'figures/tot_sims_spleen.fig')
   saveas(gcf,'figures/tot_sims_spleen','epsc')
   saveas(gcf,'figures/tot_sims_spleen','jpeg')

%s=subplot(4,2,2);
s=figure()
for i=1:length(dose_rel)
load(['model_generated_data/sim_runs_ ',num2str(dose_rel(i)),' _20240118']);
plot(times(ind),OP.Heart_mn(ind),[char(color(i)), '-'],'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');hold on;
plot(times(ind),OP.Heart_mn(ind)+OP.Heart_sd(ind)/10,[char(color(i)), '--'],'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');
plot(times(ind),OP.Heart_mn(ind)-OP.Heart_sd(ind)/10,[char(color(i)), '--'],'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');
end
hold off;
set(gca, 'box','off','FontSize',tick_fs,'XTick',0:8:48);
xlim([0 48]);
ylim([0 0.1]);
xlabel("Time (hr)",'FontSize',axis_fs);
ylabel("NP Conc. (\mu{g}/ml)",'FontSize',axis_fs);
axis square
title('Heart');

   saveas(gcf,'figures/tot_sims_heart.fig')
   saveas(gcf,'figures/tot_sims_heart','epsc')
   saveas(gcf,'figures/tot_sims_heart','jpeg')

%s=subplot(4,2,5);
s=figure()
for i=1:length(dose_rel)
load(['model_generated_data/sim_runs_ ',num2str(dose_rel(i)),' _20240118']);
plot(times(ind),OP.Kidneys_mn(ind),[char(color(i)), '-'],'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');hold on;
plot(times(ind),OP.Kidneys_mn(ind)+OP.Kidneys_sd(ind)/10,[char(color(i)), '--'],'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');
plot(times(ind),OP.Kidneys_mn(ind)-OP.Kidneys_sd(ind)/10,[char(color(i)), '--'],'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');
end
hold off;
set(gca, 'box','off','FontSize',tick_fs,'XTick',0:8:48,'YTick',0:20:120);
xlim([0 48]);
ylim([0 120]);
xlabel("Time (hr)",'FontSize',axis_fs);
ylabel("NP Conc. (\mu{g}/ml)",'FontSize',axis_fs);
axis square
title('Kidneys');

   saveas(gcf,'figures/tot_sims_kidney.fig')
   saveas(gcf,'figures/tot_sims_kidney','epsc')
   saveas(gcf,'figures/tot_sims_kidney','jpeg')

%s=subplot(4,2,6);
s=figure()
for i=1:length(dose_rel)
load(['model_generated_data/sim_runs_ ',num2str(dose_rel(i)),' _20240118']);
plot(times(ind),OP.Lungs_mn(ind),[char(color(i)), '-'],'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');hold on;
plot(times(ind),OP.Lungs_mn(ind)+OP.Lungs_sd(ind)/10,[char(color(i)), '--'],'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');
plot(times(ind),OP.Lungs_mn(ind)-OP.Lungs_sd(ind)/10,[char(color(i)), '--'],'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');
end
hold off;
set(gca, 'box','off','FontSize',tick_fs,'XTick',0:8:48,'YTick',0:20:120);
xlim([0 48]);
ylim([0 120]);
xlabel("Time (hr)",'FontSize',axis_fs);
ylabel("NP Conc. (\mu{g}/ml)",'FontSize',axis_fs);
axis square
title('Lungs');

   saveas(gcf,'figures/tot_sims_lung.fig')
   saveas(gcf,'figures/tot_sims_lung','epsc')
   saveas(gcf,'figures/tot_sims_lung','jpeg')

%s=subplot(4,2,7);
s=figure()
for i=1:length(dose_rel)
load(['model_generated_data/sim_runs_ ',num2str(dose_rel(i)),' _20240118']);
plot(times(ind),OP.Bone_mn(ind),[char(color(i)), '-'],'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');hold on;
plot(times(ind),OP.Bone_mn(ind)+OP.Bone_sd(ind)/10,[char(color(i)), '--'],'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');
plot(times(ind),OP.Bone_mn(ind)-OP.Bone_sd(ind)/10,[char(color(i)), '--'],'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');
end
hold off;
set(gca, 'box','off','FontSize',tick_fs,'XTick',0:8:48,'YTick',0:20:120);
xlim([0 48]);
ylim([0 120]);
xlabel("Time (hr)",'FontSize',axis_fs);
ylabel("NP Conc. (\mu{g}/ml)",'FontSize',axis_fs);
axis square
title('Bone');

   saveas(gcf,'figures/tot_sims_bone.fig')
   saveas(gcf,'figures/tot_sims_bone','epsc')
   saveas(gcf,'figures/tot_sims_bone','jpeg')

%s=subplot(4,2,8);
s=figure()
for i=1:length(dose_rel)
load(['model_generated_data/sim_runs_ ',num2str(dose_rel(i)),' _20240118']);
Brain_Rest_mn=OP.Brain_mn(ind)+OP.Rest_mn(ind);
Brain_Rest_sd=sqrt(OP.Brain_sd(ind).^2+OP.Rest_mn(ind).^2);
plot(times(ind),Brain_Rest_mn,[char(color(i)), '-'],'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');hold on;
plot(times(ind),Brain_Rest_mn+Brain_Rest_sd/10,[char(color(i)), '--'],'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');
plot(times(ind),Brain_Rest_mn-Brain_Rest_sd/10,[char(color(i)), '--'],'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');
end
hold off;
set(gca, 'box','off','FontSize',tick_fs,'XTick',0:8:48,'YTick',0:0.05:0.1);
xlim([0 48]);
ylim([0 0.1]);
xlabel("Time (hr)",'FontSize',axis_fs);
ylabel("NP Conc. (\mu{g}/ml)",'FontSize',axis_fs);
axis square
title('Brain + Body');

   saveas(gcf,'figures/tot_sims_body.fig')
   saveas(gcf,'figures/tot_sims_body','epsc')
   saveas(gcf,'figures/tot_sims_body','jpeg')

%figure('units','normalized','outerposition',[0 0 1 0.6])
axis_fs = 15;
tick_fs=12;

color_mat = [0 0 1; 1 0 0;0 0 1; 1 0 0;0 0 1; 1 0 0;0 0 1; 1 0 0];

   y=[0 0 0 0 0; ...
       0 0 0 0 0; ...
       0 0 0 0 0; ...
       0 0 0 0 0; ...
       0 0 0 0 0;...
       0 0 0 0 0; ...
       0 0 0 0 0;...
       0 0 0 0 0];

   FACS_indx = [2 5 3 6 1];

y(1,:)=[OP_organ_sum3.Heart_mn OP_organ_sum3.Lungs_mn ...
    OP_organ_sum3.Kidneys_mn OP_organ_sum3.Spleen_mn OP_organ_sum3.Bone_mn];

FACS_mean=FACS_summary_time.mean(FACS_summary_time.time==3);
y(2,:)=FACS_mean(FACS_indx);

y(3,:)=[OP_organ_sum6.Heart_mn OP_organ_sum6.Lungs_mn ...
    OP_organ_sum6.Kidneys_mn OP_organ_sum6.Spleen_mn OP_organ_sum6.Bone_mn];

FACS_mean=FACS_summary_time.mean(FACS_summary_time.time==6);
y(4,:)=FACS_mean(FACS_indx);

y(5,:)=[OP_organ_sum24.Heart_mn OP_organ_sum24.Lungs_mn ...
    OP_organ_sum24.Kidneys_mn OP_organ_sum24.Spleen_mn OP_organ_sum24.Bone_mn];

FACS_mean=FACS_summary_time.mean(FACS_summary_time.time==24);
y(6,:)=FACS_mean(FACS_indx);

y(7,:)=[OP_organ_sum48.Heart_mn OP_organ_sum48.Lungs_mn ...
    OP_organ_sum48.Kidneys_mn OP_organ_sum48.Spleen_mn OP_organ_sum48.Bone_mn];

FACS_mean=FACS_summary.mean(FACS_summary.Dose==0.5);
y(8,:)=FACS_mean(FACS_indx);

%y(:,1)=y(:,1)*1e3;
y=y';

   err=[0 0 0 0 0; ...
       0 0 0 0 0;...
       0 0 0 0 0; ...
       0 0 0 0 0; ...
       0 0 0 0 0;...
       0 0 0 0 0; ...
       0 0 0 0 0;...
       0 0 0 0 0];
err(1,:)=[OP_organ_sum3.Heart_sd OP_organ_sum3.Lungs_sd ...
    OP_organ_sum3.Kidneys_sd OP_organ_sum3.Spleen_sd OP_organ_sum3.Bone_sd];

FACS_sd=FACS_summary_time.sd(FACS_summary_time.time==3);
err(2,:)=FACS_sd(FACS_indx);

err(3,:)=[OP_organ_sum6.Heart_sd OP_organ_sum6.Lungs_sd ...
    OP_organ_sum6.Kidneys_sd OP_organ_sum6.Spleen_sd OP_organ_sum6.Bone_sd];

FACS_sd=FACS_summary_time.sd(FACS_summary_time.time==6);
err(4,:)=FACS_sd(FACS_indx);

err(5,:)=[OP_organ_sum24.Heart_sd OP_organ_sum24.Lungs_sd ...
    OP_organ_sum24.Kidneys_sd OP_organ_sum24.Spleen_sd OP_organ_sum24.Bone_sd];

FACS_sd=FACS_summary_time.sd(FACS_summary_time.time==24);
err(6,:)=FACS_sd(FACS_indx);

err(7,:)=[OP_organ_sum48.Heart_sd OP_organ_sum48.Lungs_sd ...
    OP_organ_sum48.Kidneys_sd OP_organ_sum48.Spleen_sd OP_organ_sum48.Bone_sd];

FACS_sd=FACS_summary.sd(FACS_summary.Dose==0.5);
err(8,:)=FACS_sd(FACS_indx);

%err(:,1)=err(:,1)*1e3;
err=err';

Z3=(y(:,1)-y(:,2))./sqrt(err(:,1).^2./100 + err(:,2).^2./3);
Z6=(y(:,3)-y(:,4))./sqrt(err(:,3).^2./100 + err(:,4).^2./3);
Z24=(y(:,5)-y(:,6))./sqrt(err(:,5).^2./100 + err(:,6).^2./3);
Z48=(y(:,7)-y(:,8))./sqrt(err(:,7).^2./100 + err(:,8).^2./3);

y=y(:,3:6);
err=err(:,3:6);

%subplot(1,3,1)
figure()
b=bar(y);hold on;

b(1).FaceColor=color_mat(1,:);
b(2).FaceColor=color_mat(2,:);
b(3).FaceColor=color_mat(3,:);
b(4).FaceColor=color_mat(4,:);
%b(5).FaceColor=color_mat(5,:);
%b(6).FaceColor=color_mat(6,:);
%b(7).FaceColor=color_mat(7,:);
%b(8).FaceColor=color_mat(8,:);

b(1).FaceAlpha=1;
b(2).FaceAlpha=1;
b(3).FaceAlpha=0.5;
b(4).FaceAlpha=0.5;
%b(5).FaceAlpha=0.4;
%b(6).FaceAlpha=0.4;
%b(7).FaceAlpha=0.1;
%b(8).FaceAlpha=0.1;

ngroups = size(y, 1);
nbars = size(y, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, y(:,i), 0*err(:,i), err(:,i), 'k.');
end
     %text(0.1,1.7e-1,'(x10^{-3})')
%     text(-0.5,7e-1,'Z (6hr) = ','Color','red')
%     text(-0.5,6e-1,'Z (24hr) = ')
%     text(-0.5,5e-1,'Z (48hr) = ')
% for i = 1:length(Z3)
%     text(i-0.2,8e-1,num2str(round(Z3(i),1)))
%     text(i-0.2,7e-1,num2str(round(Z6(i),1)),'Color','red')
%     text(i-0.2,6e-1,num2str(round(Z24(i),1)),'Color','red')
%     text(i-0.2,5e-1,num2str(round(Z48(i),1)))
% end
hold off
ylim([0 0.5]);
xlim([0 5.6]);
set(gca, 'box','off','XTickLabel',FACS_summary_time.Organ(FACS_indx)','Fontsize',tick_fs);%,'YTick',5:1:11);
ylabel("Organ NP Mass (Relative to Liver)",'FontSize',axis_fs);
legend('Model, 6hr','Data, 6hr',...
    'Model, 24hr','Data, 24hr',...
    'Location','northwest');
axis square

   saveas(gcf,'figures/mod_verify1.fig')
   saveas(gcf,'figures/mod_verify1','epsc')
   saveas(gcf,'figures/mod_verify1','jpeg')

   y=[0 0 0 0 0 0 0; ...
       0 0 0 0 0 0 0];

times=(1:12000)/12000*48;

blood_time=blood_summary.time(blood_summary.dose==0.5);
blood_mean=blood_summary.mean(blood_summary.dose==0.5);
blood_sd=blood_summary.sd(blood_summary.dose==0.5);
load(['model_generated_data/sim_runs_',num2str(0.5),'_20240116']);

y(1,:)=[OP.Blood_mn(abs(times-0.03)==min(abs(times-0.03)))...
    OP.Blood_mn(times==1) OP.Blood_mn(times==2) ...
    OP.Blood_mn(times==8) OP.Blood_mn(times==10)...
    OP.Blood_mn(times==24) OP.Blood_mn(times==48)];

y(2,:)=[blood_mean(blood_time==0.03) ...
        blood_mean(blood_time==1) ...
    blood_mean(blood_time==2) ...
    blood_mean(blood_time==8)...
        blood_mean(blood_time==10)...
     blood_mean(blood_time==24)...
    blood_mean(blood_time==48)];
y=y';

   err=[0 0 0 0 0 0 0; ...
       0 0 0 0 0 0 0];

err(1,:)=[OP.Blood_sd(abs(times-0.03)==min(abs(times-0.03))) OP.Blood_sd(times==1) OP.Blood_sd(times==2) ...
    OP.Blood_sd(times==8) OP.Blood_sd(times==10) OP.Blood_sd(times==24) OP.Blood_sd(times==48)];

err(2,:)=[blood_sd(blood_time==0.03) blood_sd(blood_time==1) ...
    blood_sd(blood_time==2) blood_sd(blood_time==8)...
    blood_sd(blood_time==10)...
    blood_sd(blood_time==24) blood_sd(blood_time==48)];
err=err';

Z=(y(:,1)-y(:,2))./sqrt(err(:,1).^2./100 + err(:,2).^2./3);

% subplot(1,3,1)
% b=bar(y);hold on;
% ngroups = size(y, 1);
% nbars = size(y, 2);
% % Calculating the width for each bar group
% groupwidth = min(0.8, nbars/(nbars + 1.5));
% for i = 1:nbars
%     x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
%     errorbar(x, y(:,i), 0*err(:,i), err(:,i), 'k.');
% end
% %     text(0,450,'Z = ')
% % 
% %     colorZ = [0 0 0;0 0 0;1 0 0;0 0 0;1 0 0;0 0 0;0 0 0]; 
% % for i = 1:length(Z)-1
% %     text(i-0.2,450,num2str(round(Z(i),1)),'Color',colorZ(i,:))
% % end
% b(1).FaceColor=color_mat(1,:);
% b(2).FaceColor=color_mat(2,:);
% b(1).FaceAlpha=0.6;
% b(2).FaceAlpha=0.6;
% ylim([0 420]);
% xlim([-0.2 7.6]);
% hold off
% set(gca, 'box','off','XTickLabel',{'0.03' '1' '2' '8' '10' '24' '48'},'Fontsize',tick_fs);%,'YTick',5:1:11);
% legend('Model','Data','Location','northeast');
% xlabel("Time (hr)")
% ylabel("Blood NP Concentration (\mu{g}/ml)",'FontSize',axis_fs);
% axis square

y=[Z6'; Z24'; 0*Z6'];

%subplot(1,3,3)
figure()
b1=bar(abs(y));hold on;
ngroups = size(y, 1);
nbars = size(y, 2);
% Calculating the width for each bar group
ylim([0 30]);
%xlim([0 14]);
b1(1).FaceAlpha=0.8;
b1(2).FaceAlpha=0.8;
b1(3).FaceAlpha=0.8;
b1(4).FaceAlpha=0.8;
b1(5).FaceAlpha=0.8;
b=bar(abs(Z([3 5 6])));
b.FaceColor='r';
b.FaceAlpha=0.8;
b.XData=[2.85 3 3.15];
set(gca, 'box','off','XTickLabel',{'6' '24' '2,10,24'},'Fontsize',tick_fs);%,'YTick',5:1:11);
legend([FACS_summary_time.Organ(FACS_indx)' 'Blood'],'Location','northeast');
xlabel("Time (hr)")
ylabel("Statistical Model Error (|Z|)",'FontSize',axis_fs);
axis square

   saveas(gcf,'figures/mod_verify2.fig')
   saveas(gcf,'figures/mod_verify2','epsc')
   saveas(gcf,'figures/mod_verify2','jpeg')

   y1=[0 0 0 0 0; ...
       0 0 0 0 0;
       0 0 0 0 0];

   FACS_indx = [2 5 3 6 1];

FACS_mean=FACS_summary_time.mean(FACS_summary_time.time==6);
y1(1,:)=[OP_organ_sum6.Heart_mn OP_organ_sum6.Lungs_mn ...
    OP_organ_sum6.Kidneys_mn OP_organ_sum6.Spleen_mn OP_organ_sum6.Bone_mn]./FACS_mean(FACS_indx)';

FACS_mean=FACS_summary_time.mean(FACS_summary_time.time==24);
y1(2,:)=[OP_organ_sum24.Heart_mn OP_organ_sum24.Lungs_mn ...
    OP_organ_sum24.Kidneys_mn OP_organ_sum24.Spleen_mn OP_organ_sum24.Bone_mn]./FACS_mean(FACS_indx)';

%y1=y1';

y=[OP.Blood_mn(times==2)/blood_mean(blood_time==2)...
    OP.Blood_mn(times==2)/blood_mean(blood_time==10)...
    OP.Blood_mn(times==2)/blood_mean(blood_time==24)];

%subplot(1,3,2)
figure()
b1=bar(abs(y1));hold on;
% Calculating the width for each bar group
ylim([0 3]);
%xlim([0 14]);
b1(1).FaceAlpha=0.8;
b1(2).FaceAlpha=0.8;
b1(3).FaceAlpha=0.8;
b1(4).FaceAlpha=0.8;
b1(5).FaceAlpha=0.8;
b=bar(abs(y));
b.FaceColor='r';
b.FaceAlpha=0.8;
b.XData=[2.85 3 3.15];
set(gca, 'box','off','XTickLabel',{'6' '24' '2,10,24'},'Fontsize',tick_fs);%,'YTick',5:1:11);
legend([FACS_summary_time.Organ(FACS_indx)' 'Blood'],'Location','northeast');
xlabel("Time (hr)")
ylabel("Relative Model Error (Model/Data)",'FontSize',axis_fs);
axis square

   saveas(gcf,'figures/mod_verify3.fig')
   saveas(gcf,'figures/mod_verify3','epsc')
   saveas(gcf,'figures/mod_verify3','jpeg')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()  
y_keep=[0 0];
for i=1:length(dose_rel)
y=[0 0 0 0 0 0 0; ...
       0 0 0 0 0 0 0];

times=(1:12000)/12000*48;

blood_time=blood_summary.time(blood_summary.dose==dose_rel(i));
blood_mean=blood_summary.mean(blood_summary.dose==dose_rel(i));
blood_sd=blood_summary.sd(blood_summary.dose==dose_rel(i));
load(['model_generated_data/sim_runs_',num2str(dose_rel(i)),'_20240116']);

y(1,:)=[OP.Blood_mn(abs(times-0.03)==min(abs(times-0.03)))...
    OP.Blood_mn(times==1) OP.Blood_mn(times==2) ...
    OP.Blood_mn(times==8) OP.Blood_mn(times==10)...
    OP.Blood_mn(times==24) OP.Blood_mn(times==48)];

y(2,:)=[blood_mean(blood_time==0.03) ...
        blood_mean(blood_time==1) ...
    blood_mean(blood_time==2) ...
    blood_mean(blood_time==8)...
        blood_mean(blood_time==10)...
     blood_mean(blood_time==24)...
    blood_mean(blood_time==48)];
y=y';

y_keep=[y_keep; y];
scatter(y(:,2),y(:,1),'k*'); hold on;
end
rl=refline(1,0);
rl.Color='k';
rl.LineStyle='--';

xlim([0 1750]);
ylim([0 1750]);
xlabel("Experimental Blood NP Concentration (\mu{g}/ml)")
ylabel("Predicted Blood NP Concentration (\mu{g}/ml)",'FontSize',axis_fs);
axis square
y_keep=y_keep(2:end,:);
mdl1 = fitlm(y_keep(:,1),y_keep(:,2),'intercept',false);
% %     text(i-0.2,450,num2str(round(Z(i),1)),'Color',colorZ(i,:))
text(1100, 400, ['R^2 = ', num2str(round(mdl1.Rsquared.Ordinary,2))]);
text(1100, 250, ['RMSE = ', num2str(round(rmse(y_keep(:,1),y_keep(:,2)),1))]);

   saveas(gcf,'figures/blood_R2.fig')
   saveas(gcf,'figures/blood_R2','epsc')
   saveas(gcf,'figures/blood_R2','jpeg')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()

   y=[0 0 0 0 0 ...
       0 0 0 0 0 ...
       0 0 0 0 0 ...
       0 0 0 0 0; ...
       0 0 0 0 0 ...
       0 0 0 0 0 ...
       0 0 0 0 0 ...
       0 0 0 0 0];

   FACS_indx = [2 5 3 6 1];

y(1,:)=[OP_organ_sum3.Heart_mn OP_organ_sum3.Lungs_mn ...
    OP_organ_sum3.Kidneys_mn OP_organ_sum3.Spleen_mn OP_organ_sum3.Bone_mn...
    OP_organ_sum6.Heart_mn OP_organ_sum6.Lungs_mn ...
    OP_organ_sum6.Kidneys_mn OP_organ_sum6.Spleen_mn OP_organ_sum6.Bone_mn...
    OP_organ_sum24.Heart_mn OP_organ_sum24.Lungs_mn ...
    OP_organ_sum24.Kidneys_mn OP_organ_sum24.Spleen_mn OP_organ_sum24.Bone_mn...
    OP_organ_sum48.Heart_mn OP_organ_sum48.Lungs_mn ...
    OP_organ_sum48.Kidneys_mn OP_organ_sum48.Spleen_mn OP_organ_sum48.Bone_mn];

FACS_mean1=FACS_summary_time.mean(FACS_summary_time.time==3);
FACS_mean2=FACS_summary_time.mean(FACS_summary_time.time==6);
FACS_mean3=FACS_summary_time.mean(FACS_summary_time.time==24);
FACS_mean4=FACS_summary.mean(FACS_summary.Dose==0.5);

y(2,:)=[FACS_mean1(FACS_indx)' FACS_mean2(FACS_indx)' ...
    FACS_mean3(FACS_indx)' FACS_mean4(FACS_indx)'];

y=y';

y_keep=[y_keep; y];
scatter(y(:,2),y(:,1),'k*');
rl=refline(1,0);
rl.Color='k';
rl.LineStyle='--';

xlim([0 0.3]);
ylim([0 0.3]);
xlabel("Experimental Organ NP Mass (Relative to Liver)")
ylabel("Predicted Organ NP Mass (Relative to Liver)",'FontSize',axis_fs);
axis square

mdl1 = fitlm(y(:,1),y(:,2),'intercept',false);
% %     text(i-0.2,450,num2str(round(Z(i),1)),'Color',colorZ(i,:))
text(0.2, 0.1, ['R^2 = ', num2str(round(mdl1.Rsquared.Ordinary,2))]);
text(0.2, 0.05, ['RMSE = ', num2str(round(rmse(y(:,1),y(:,2)),2))]);

   saveas(gcf,'figures/organs_R2.fig')
   saveas(gcf,'figures/organs_R2','epsc')
   saveas(gcf,'figures/organs_R2','jpeg')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ind = 10:12000;
% 
% dose_rel=flip([0.1 0.5 2 2.5 3.5]);
% color = flip({'r' 'm' 'g' 'c' 'b'});
% 
% figure('units','normalized','outerposition',[0 0 0.8 0.6])
% s=subplot(1,2,1);
% for i=1:length(dose_rel)
% blood_time=blood_summary.time(blood_summary.dose==dose_rel(i));
% blood_mean=blood_summary.mean(blood_summary.dose==dose_rel(i));
% blood_sd=blood_summary.sd(blood_summary.dose==dose_rel(i));
% 
% load(['model_generated_data/sim_runs_',num2str(dose_rel(i)),'_20240116']);
% 
% plot(times(ind),OP.Blood_mn(ind),[char(color(i)), '-'],'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');hold on;
% plot(times(ind),OP.Blood_mn(ind)+OP.Blood_sd(ind),[char(color(i)), '--'],'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');
% plot(times(ind),OP.Blood_mn(ind)-OP.Blood_sd(ind),[char(color(i)), '--'],'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');
% plot(blood_time,blood_mean,[char(color(i)), '.'],'MarkerSize',30);
% errorbar(blood_time,blood_mean,blood_sd/sqrt(3),[char(color(i)), '.'],'LineWidth',2);
% end
% 
% hold off;
% set(gca, 'box','off','FontSize',tick_fs,'XTick',0:6:24);
% xlim([0 48]);
% ylim([0 2000]);
% legend('3.5mg','','','','',...
%     '2.5mg','','','','',...
%     '2mg','','','','',...
%     '0.5mg','','','','',...
%     '0.1mg','','','','',...
%     'Location','northeast');
% xlabel("Time (hr)",'FontSize',axis_fs);
% ylabel("NP Conc. (\mu{g}/ml)",'FontSize',axis_fs);
% axis square
% title(s,'Blood');
% 
%    times=(1:12000)/12000*48;
%     time_common=[1 2 8 10 24];
% rel_tab=0*[time_common; time_common; time_common; time_common; time_common];
% Z_tab=0*[time_common; time_common; time_common; time_common; time_common];
% 
% for i=1:length(dose_rel)
% blood_time=blood_summary.time(blood_summary.dose==dose_rel(i));
% blood_mean=blood_summary.mean(blood_summary.dose==dose_rel(i));
% blood_sd=blood_summary.sd(blood_summary.dose==dose_rel(i));
% 
% load(['model_generated_data/sim_runs_',num2str(dose_rel(i)),'_20240116']);
% 
% rel_tab(i,:)=(OP.Blood_mn(ismember(times,time_common))-blood_mean(ismember(blood_time,time_common)))./...
%     dose_rel(i);
% 
% Z_tab(i,:)=(OP.Blood_mn(ismember(times,time_common))-blood_mean(ismember(blood_time,time_common)))./...
%     sqrt(OP.Blood_sd(ismember(times,time_common)).^2/100 + blood_sd(ismember(blood_time,time_common)).^2/3);
% 
% end
% 
% rel_tab1=0*rel_tab;
% Z_tab1=0*rel_tab;
% 
% for i=1:length(dose_rel)
%     rel_tab1(i,:)=rel_tab(i,:)-rel_tab(dose_rel==0.5,:);
%     Z_tab1(i,:)=Z_tab(i,:)-Z_tab(dose_rel==0.5,:);
% end
% 
% clims=[min(min([Z_tab1])) max(max([Z_tab1]))];
% 
% s=subplot(1,2,2);
% imagesc(Z_tab1,clims)
% colorbar
% axis square
% 
% saveas(gcf,'figures/blood_comp.fig')
% saveas(gcf,'figures/blood_comp','epsc')
% saveas(gcf,'figures/blood_comp','jpeg')


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %clear all
% %close all
% 
% load('model_generated_data/cui_blood_meas_20231210');
% load('model_generated_data/PACE_blood_summary_20231103');
% load('model_generated_data/organ_FACS_cui_20231210');
% 
% ind = 1:12:12000;
% times = ind/12000*48;
% axis_fs = 15;
% tick_fs=12;
% 
% color = {'red' "#EDB120" 'green' 'cyan' 'blue'};
% 
% %figure('units','normalized','outerposition',[0 0 0.6 0.6])
% figure()
% plot(times,OP.Blood_mn(ind)/max(OP.Blood_mn)*100,...
%     'k','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');hold on;
% plot(blood_conc_cui.time,blood_conc_cui.mean,'r.','MarkerSize',15);hold on;
% plot(blood_summary.time,blood_summary.mean/max(blood_summary.mean)*100,'b.','MarkerSize',15);hold on;
% errorbar(blood_conc_cui.time,blood_conc_cui.mean,blood_conc_cui.sd,'r.','LineWidth',2); hold on;
% errorbar(blood_summary.time,blood_summary.mean/max(blood_summary.mean)*100,...
%     blood_summary.sd/max(blood_summary.mean)*100,'b.','LineWidth',2); hold on;
% plot(times,(OP.Blood_mn(ind)+OP.Blood_sd(ind))/max(OP.Blood_mn)*100,...
%     'k--','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');hold on;
% plot(times,(OP.Blood_mn(ind)-OP.Blood_sd(ind))/max(OP.Blood_mn)*100,...
%     'k--','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');
% legend('Model','Cui et al. 2019','Fit Data','','','','');
% hold off;
% set(gca, 'box','off','FontSize',tick_fs);%,'YTick',5:1:11);
% xlim([0 48]);
% ylim([0 120]);
% xlabel("Time (hr)",'FontSize',axis_fs);
% ylabel("% of Initial NP Dose",'FontSize',axis_fs);
% axis square
% 
% 
%    saveas(gcf,'figures/cui_validate.fig')
%    saveas(gcf,'figures/cui_validate','epsc')
%    saveas(gcf,'figures/cui_validate','jpeg')
% 
%    y=[0 0 0 0 0 ; ...
%        0 0 0 0 0 ; ...
%        0 0 0 0 0 ];
% 
% times=(1:12000)/12000*48;
% 
% y(1,:)=[OP.Blood_mn(abs(times-0.0333)==min(abs(times-0.0333))) OP.Blood_mn(times==2) OP.Blood_mn(times==4) ...
%      OP.Blood_mn(times==10) OP.Blood_mn(times==24)];
% y(2,:)=[blood_summary.mean(blood_summary.time==0.03) ...
%         blood_summary.mean(blood_summary.time==2) ...
%     blood_summary.mean(blood_summary.time==4) ...
%         blood_summary.mean(blood_summary.time==10)...
%      blood_summary.mean(blood_summary.time==24)];
% 
% y(3,:)=[blood_conc_cui.mean(blood_conc_cui.time==blood_conc_cui.time(1)) ...
%         blood_conc_cui.mean(blood_conc_cui.time==2) ...
%     blood_conc_cui.mean(blood_conc_cui.time==4) ...
%         blood_conc_cui.mean(blood_conc_cui.time==10)...
%      blood_conc_cui.mean(blood_conc_cui.time==24)];
% y=y';
% 
%    err=[0 0 0 0 0 ; ...
%    0 0 0 0 0 ; ...    
%    0 0 0 0  0];
% 
% err(1,:)=[OP.Blood_sd(abs(times-0.0333)==min(abs(times-0.0333))) OP.Blood_mn(times==2) OP.Blood_mn(times==4) ...
%     OP.Blood_sd(times==10) OP.Blood_sd(times==24)];
% err(2,:)=[blood_summary.sd(blood_summary.time==0.03) blood_summary.sd(blood_summary.time==2) ...
%     blood_summary.sd(blood_summary.time==4)...
%     blood_summary.sd(blood_summary.time==10)...
%     blood_summary.sd(blood_summary.time==24)];
% 
% err(3,:)=[blood_conc_cui.sd(blood_conc_cui.time==blood_conc_cui.time(1)) ...
%     blood_conc_cui.sd(blood_conc_cui.time==2) ...
%     blood_conc_cui.sd(blood_conc_cui.time==4) ...
%     blood_conc_cui.sd(blood_conc_cui.time==10)...
%     blood_conc_cui.sd(blood_conc_cui.time==24)];
% err=err';
% 
% Z=(y(:,1)-y(:,3))./sqrt(err(:,1).^2./100 + err(:,3).^2./3)
% 
% %subplot(1,2,2)
% figure()
% bar(y);hold on;
% ngroups = size(y, 1);
% nbars = size(y, 2);
% % Calculating the width for each bar group
% groupwidth = min(0.8, nbars/(nbars + 1.5));
% for i = 1:nbars
%     x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
%     errorbar(x, y(:,i), 0*err(:,i), err(:,i), 'k.');
% end
%     text(0,250,'Z = ')
% 
% for i = 1:length(Z)
%     text(i-0.2,250,num2str(round(Z(i),1)))
% end
% ylim([0 350]);
% xlim([-0.2 7.6]);
% hold off
% set(gca, 'box','off','XTickLabel',{'0.03' '2' '4' '10' '24'});%,'YTick',5:1:11);
% xlabel("Time (hr)")
% ylabel("Blood NP Concentration (\mu{g}/ml)",'FontSize',axis_fs);
% axis square
% legend('Model','Fit Data','Cui et al. 2019','','Location','northwest');
% 


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %clear all
% %close all
% 
% load('model_generated_data/kdeg_sim_runs_20231212');
% load('model_generated_data/PACE_blood_summary_20231103');
% 
% ind = 1:12:12000;
% times = ind/12000*48;
% axis_fs = 12;
% tick_fs=10;
% 
% color = {'red' "#EDB120" 'green' 'cyan' 'blue'};
% 
% figure('units','normalized','outerposition',[0 0 0.4 1])
% 
% s=subplot(4,2,1);
% plot(times,OP.Blood_mn(ind),'k','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');hold on;
% plot(times,OP.Blood_mn(ind)+OP.Blood_sd(ind),'k--','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');hold on;
% plot(times,OP.Blood_mn(ind)-OP.Blood_sd(ind),'k--','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');hold on;
% plot(blood_summary.time,blood_summary.mean,'b.','MarkerSize',15);hold on;
% errorbar(blood_summary.time,blood_summary.mean,blood_summary.sd,'b.','LineWidth',2); 
% hold off;
% set(gca, 'box','off','FontSize',tick_fs);%,'YTick',5:1:11);
% xlim([0 48]);
% ylim([0 200]);
% xlabel("Time (hr)",'FontSize',axis_fs);
% ylabel("NP Concentration (\mu{g}/ml)",'FontSize',axis_fs);
% axis square
% title(s,'Blood');
% 
% s=subplot(4,2,2);
% plot(times,OP.Liver_mn(ind),'k','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');hold on;
% plot(times,OP.Liver_mn(ind)+OP.Liver_sd(ind),'k--','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');hold on;
% plot(times,OP.Liver_mn(ind)-OP.Liver_sd(ind),'k--','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');hold on;
% hold off;
% set(gca, 'box','off','FontSize',tick_fs);%,'YTick',5:1:11);
% xlim([0 48]);
% ylim([0 500]);
% xlabel("Time (hr)",'FontSize',axis_fs);
% ylabel("NP Concentration (\mu{g}/ml)",'FontSize',axis_fs);
% axis square
% title(s,'Liver');
% 
% s=subplot(4,2,3);
% plot(times,OP.Spleen_mn(ind),'k','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');hold on;
% plot(times,OP.Spleen_mn(ind)+OP.Spleen_sd(ind),'k--','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');hold on;
% plot(times,OP.Spleen_mn(ind)-OP.Spleen_sd(ind),'k--','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');hold on;
% hold off;
% set(gca, 'box','off','FontSize',tick_fs);%,'YTick',5:1:11);
% xlim([0 48]);
% ylim([0 15]);
% xlabel("Time (hr)",'FontSize',axis_fs);
% ylabel("NP Concentration (\mu{g}/ml)",'FontSize',axis_fs);
% axis square
% title(s,'Spleen');
% 
% s=subplot(4,2,4);
% plot(times,OP.Heart_mn(ind),'k','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');hold on;
% plot(times,OP.Heart_mn(ind)+OP.Heart_sd(ind),'k--','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');hold on;
% plot(times,OP.Heart_mn(ind)-OP.Heart_sd(ind),'k--','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');hold on;
% hold off;
% set(gca, 'box','off','FontSize',tick_fs);%,'YTick',5:1:11);
% xlim([0 48]);
% ylim([0 0.2]);
% xlabel("Time (hr)",'FontSize',axis_fs);
% ylabel("NP Concentration (\mu{g}/ml)",'FontSize',axis_fs);
% axis square
% title(s,'Heart');
% 
% s=subplot(4,2,5);
% plot(times,OP.Kidneys_mn(ind),'k','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');hold on;
% plot(times,OP.Kidneys_mn(ind)+OP.Kidneys_sd(ind),'k--','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');hold on;
% plot(times,OP.Kidneys_mn(ind)-OP.Kidneys_sd(ind),'k--','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');hold on;
% hold off;
% set(gca, 'box','off','FontSize',tick_fs);%,'YTick',5:1:11);
% xlim([0 48]);
% ylim([0 0.2]);
% xlabel("Time (hr)",'FontSize',axis_fs);
% ylabel("NP Concentration (\mu{g}/ml)",'FontSize',axis_fs);
% axis square
% title(s,'Kidneys');
% 
% s=subplot(4,2,6);
% plot(times,OP.Lungs_mn(ind),'k','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');hold on;
% plot(times,OP.Lungs_mn(ind)+OP.Lungs_sd(ind),'k--','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');hold on;
% plot(times,OP.Lungs_mn(ind)-OP.Lungs_sd(ind),'k--','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');hold on;
% hold off;
% set(gca, 'box','off','FontSize',tick_fs);%,'YTick',5:1:11);
% xlim([0 48]);
% ylim([0 15]);
% xlabel("Time (hr)",'FontSize',axis_fs);
% ylabel("NP Concentration (\mu{g}/ml)",'FontSize',axis_fs);
% axis square
% title(s,'Lungs');
% 
% s=subplot(4,2,7);
% plot(times,OP.Bone_mn(ind),'k','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');hold on;
% plot(times,OP.Bone_mn(ind)+OP.Lungs_sd(ind),'k--','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');hold on;
% plot(times,OP.Bone_mn(ind)-OP.Lungs_sd(ind),'k--','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');hold on;
% hold off;
% set(gca, 'box','off','FontSize',tick_fs);%,'YTick',5:1:11);
% xlim([0 48]);
% ylim([0 15]);
% xlabel("Time (hr)",'FontSize',axis_fs);
% ylabel("NP Concentration (\mu{g}/ml)",'FontSize',axis_fs);
% axis square
% title(s,'Bone');
% 
%    saveas(gcf,'figures/kdeg_sims.fig')
%    saveas(gcf,'figures/kdeg_sims','epsc')
%    saveas(gcf,'figures/kdeg_sims','jpeg')
