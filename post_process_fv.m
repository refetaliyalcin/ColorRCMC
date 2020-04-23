% close all
cie_inc
lamda_ub=max(lamda)*10^6;
T_sur=300;
T_amb=300;
Solar_l=I_solar(lamda)';
BB_Tamb_l=I_bb(lamda,T_amb)';
BB_Tsur_l=I_bb(lamda,T_sur)';
trans_atm_pure=emis_atm_new(lamda);
emit_atm=1-trans_atm_pure.^(1./cos(polar_angle_rad));

emittance=abs_lamda;

P_rad=zeros(length(f_v_2),length(f_v_1));
P_atm=zeros(length(f_v_2),length(f_v_1));
P_sol=zeros(length(f_v_2),length(f_v_1));
P_net=zeros(length(f_v_2),length(f_v_1));
for k=1:length(f_v_1)
    for j=1:length(f_v_2)
        emittance_rad=trapz(lamda,BB_Tsur_l.*squeeze(emittance(:,j,k,:))',2);
        emittance_atm=trapz(lamda,BB_Tamb_l.*squeeze(emittance(:,j,k,:))'.*emit_atm',2);
        P_rad(j,k)=trapz(polar_angle_rad,2*cos(polar_angle_rad).*sin(polar_angle_rad).*emittance_rad');
        P_atm(j,k)=trapz(polar_angle_rad,2*cos(polar_angle_rad).*sin(polar_angle_rad).*emittance_atm');
        P_sol(j,k)=trapz(lamda,Solar_l.*squeeze(emittance(:,j,k,1))');
        P_net(j,k)=P_rad(j,k)-P_atm(j,k)- P_sol(j,k);
    end
end


figure
RGBtoXYZ=inv([0.4124564  0.3575761  0.1804375
 0.2126729  0.7151522  0.0721750
 0.0193339  0.1191920  0.9503041]);
xyz=zeros(length(f_v_2),length(f_v_1),3);
rgb=zeros(length(f_v_2),length(f_v_1),3);
for k=1:length(f_v_1)
    for j=1:length(f_v_2)
        X=trapz(col_lamda,I_.*x_.*ref_lamda(61:531,j,k,1));
        Y=trapz(col_lamda,I_.*y_.*ref_lamda(61:531,j,k,1));
        Z=trapz(col_lamda,I_.*z_.*ref_lamda(61:531,j,k,1));
        toplam=X+Y+Z;
        x=X/toplam;
        y=Y/toplam;
        z=Z/toplam;
        if ~isnan(x)
%             rgb(k,j,:)=xyz2rgb([x y z]);
            rgb(k,j,:)=RGBtoXYZ*[x;y;z];
            rgb(k,j,1)=gamma_correct(rgb(k,j,1));
            rgb(k,j,2)=gamma_correct(rgb(k,j,2));
            rgb(k,j,3)=gamma_correct(rgb(k,j,3));
            L=116*f_s(100*Y/(Yn*bolu_Y))-16;
            a=500*(f_s(100*X/(Xn*bolu_Y))-f_s(100*Y/(Yn*bolu_Y)));
            b=500*(f_s(100*Y/(Yn*bolu_Y))-f_s(100*Z/(Zn*bolu_Y)));
            Sat=sqrt(a^2+b^2)/L;
            [a1,a2]=max(rgb(k,j,:));
            rgb(k,j,:)=rgb(k,j,:)/rgb(k,j,a2);
        else
            L=0;
            Sat=0;
        end
        
        
        rectangle('Position',[f_v_2(k)-(f_v_2(2)-f_v_2(1))/2,f_v_1(j)-(f_v_1(2)-f_v_1(1))/2,f_v_2(2)-f_v_2(1),f_v_1(2)-f_v_1(1)],'FaceColor',[rgb(k,j,1)  rgb(k,j,2)  rgb(k,j,3)])
        text(f_v_2(k)-(f_v_2(2)-f_v_2(1)+450*10^-9)/4+2.4e-7,f_v_1(j),{strjoin({'P_n_e_t^''^''=',num2str(round(P_net(j,k))),'Wm^-^2'}),strjoin({'[x y]=',num2str(x,'%.2f'),num2str(y,'%.2f')}),strjoin({'S=',num2str(Sat,'%.2f')})},'FontSize',9,'HorizontalAlignment','center')
    end
end
box on
set(gcf,'renderer','painters')
xlabel('Volume Fraction (f_v)','FontSize',12)
ylabel('Volume Fraction (f_v)','FontSize',12)
ylim([min(f_v_2)-0.25*10^-6 max(f_v_2)+0.25*10^-6])
xlim([min(f_v_1)-0.25*10^-6 max(f_v_1)+0.25*10^-6])

