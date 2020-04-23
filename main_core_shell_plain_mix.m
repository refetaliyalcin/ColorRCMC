clc
clear all
close all

photon_number=10^3;
lamda=(300:1:25000)'*10^-9;
lamda_um=lamda*10^6;
lamda_color=(360:830)'*10^-9;
polar_angle=linspace(0,89.99999,30); %30 makul
polar_angle_rad=polar_angle*pi/180;

c=299792458;
w=2*pi*c./lamda; %1/s
wp=1.375*10^16; %silver
V_f=1.39*10^6;%silver m/s
g_free=3.12*10^13; % silver
C=1; %silver
n_pigment=fn_silver_n_2(lamda); %silver
k_pigment=fn_silver_k_2(lamda); %silver


f_v_1=2*[0,0.25,0.5,0.75,1]*1e-6;
pigment1_ro=15*10^-9; %green yapiyor
pigment1_ri=12*10^-9; % green yapiyor
% pigment1_ro=12*10^-9; %red yapiyor
% pigment1_ri=8.5*10^-9; %red yapiyor
pigment1_r=pigment1_ro;
e_pigment_=n_pigment.^2-k_pigment.^2; % eq. 2.30
e_pigment__=n_pigment.*k_pigment*2;
e_pigment=e_pigment_+e_pigment__*1i;

r_eff=(4/3)*(pigment1_ro^3-pigment1_ri^3)/(pigment1_ro^2+pigment1_ri^2);
delta_e_free=0;
for n=1:100
    delta_e_free_new=(-1)^n*((1i*w*C*V_f/r_eff).^n)./((w.^2+1i*w*g_free).^(n+1));
    if isfinite(delta_e_free_new)==1
        delta_e_free=delta_e_free+delta_e_free_new;
    else
        break
    end
end
delta_e_free=-delta_e_free*wp^2;
e_new_pigment=e_pigment+delta_e_free;
e_new_pigment_=real(e_new_pigment);
e_new_pigment__=imag(e_new_pigment);
n_new_pigment_1=sqrt(0.5*(e_new_pigment_+sqrt(e_new_pigment_.^2+e_new_pigment__.^2)));
k_new_pigment_1=sqrt(0.5*(-e_new_pigment_+sqrt(e_new_pigment_.^2+e_new_pigment__.^2)));


f_v_2=2*[0,0.25,0.5,0.75,1]*1e-6;
pigment2_r=25*10^-9;
e_pigment_=n_pigment.^2-k_pigment.^2; % eq. 2.30
e_pigment__=n_pigment.*k_pigment*2;
e_pigment=e_pigment_+e_pigment__*1i;

r_eff_2=pigment2_r;
delta_e_free=0;
for n=1:100
    delta_e_free_new=(-1)^n*((1i*w*C*V_f/r_eff_2).^n)./((w.^2+1i*w*g_free).^(n+1));
    if isfinite(delta_e_free_new)==1
        delta_e_free=delta_e_free+delta_e_free_new;
    else
        break
    end
end
delta_e_free=-delta_e_free*wp^2;
e_new_pigment=e_pigment+delta_e_free;
e_new_pigment_=real(e_new_pigment);
e_new_pigment__=imag(e_new_pigment);
n_new_pigment_2=sqrt(0.5*(e_new_pigment_+sqrt(e_new_pigment_.^2+e_new_pigment__.^2)));
k_new_pigment_2=sqrt(0.5*(-e_new_pigment_+sqrt(e_new_pigment_.^2+e_new_pigment__.^2)));


solar=I_solar(lamda);
thickness=500*10^-6;
thickness_um=thickness*10^6;

n_silver=fn_silver_n_2(lamda);
k_silver=fn_silver_k_2(lamda);
n_sio2=sio2_n(lamda);
k_sio2=sio2_k(lamda);

n_silver_gpu=gpuArray(n_silver);
k_silver_gpu=gpuArray(k_silver);
n_medium_gpu=gpuArray(n_sio2);
k_medium_gpu=gpuArray(k_sio2);

teta_prime=zeros(length(lamda),length(polar_angle));
sur_reflection=zeros(length(lamda),length(polar_angle));
for o=1:length(polar_angle_rad)
    for i=1:length(lamda)
        teta_prime(i,o)=F_fresnel_2(n_sio2(i),k_sio2(i),polar_angle_rad(o))*180/pi;
        cos_teta=cosd(polar_angle(o));
        sin_teta=sqrt(1-cos_teta*cos_teta);
        carpan2=1/(n_sio2(i)-1i*k_sio2(i));
        sin_x2=sin_teta*carpan2;
        cos_x2=sqrt(1-sin_x2*sin_x2);
        carpan1=cos_teta/cos_x2;
        carpan3=cos_x2/cos_teta;
        E_parallel=(carpan1-carpan2)/(carpan1+carpan2);
        R_parallel=E_parallel*conj(E_parallel);
        E_orth=(carpan3-carpan2)/(carpan3+carpan2);
        R_orth=E_orth*conj(E_orth);
        reflectance=real(R_parallel+R_orth)*0.5;
        sur_reflection(i,o)=reflectance;
    end
end

ref_lamda=zeros(length(lamda),length(polar_angle));

mu_tot_arr_1=zeros(length(lamda),1);
scat_prob_arr_1=zeros(length(lamda),1);
g_arr_1=zeros(length(lamda),1);

mu_tot_arr_2=zeros(length(lamda),1);
scat_prob_arr_2=zeros(length(lamda),1);
g_arr_2=zeros(length(lamda),1);

V_1=4*pi*pigment1_r^3/3;
V_2=4*pi*pigment2_r^3/3;
for k=1:length(f_v_1)
    for j=1:length(f_v_2)
        tic
        for i=1:length(lamda)
            n_medium_sc=n_sio2(i);
            x_1=2*pi*n_medium_sc*pigment1_ri/lamda(i);
            y_1=2*pi*n_medium_sc*pigment1_ro/lamda(i);
            m1_1=n_sio2(i)/n_medium_sc;
            m2_1=(n_new_pigment_1(i)+1i*k_new_pigment_1(i))/n_medium_sc;
            fonksiyon_1=Miecoated(m1_1,m2_1,x_1,y_1,1);
            Qsca_1=fonksiyon_1(2);
            Qabs_1=fonksiyon_1(3);
            Csca_1=pi*pigment1_r^2*Qsca_1;
            Cabs_1=pi*pigment1_r^2*Qabs_1;
            alfa_1=f_v_1(k)*Csca_1/V_1;
            beta_1=f_v_1(k)*Cabs_1/V_1+4*pi*k_sio2(i)/lamda(i);
            mu_tot_arr_1(i)=alfa_1+beta_1;
            scat_prob_arr_1(i)=alfa_1/(alfa_1+beta_1);
            g_arr_1(i)=fonksiyon_1(5);

            x_2=2*pi*n_medium_sc*pigment2_r/lamda(i);
            m_2=(n_new_pigment_2(i)+1i*k_new_pigment_2(i))/n_medium_sc;
            fonksiyon_2=Mie(m_2,x_2);
            Qsca_2=fonksiyon_2(2);
            Qabs_2=fonksiyon_2(3);
            Csca_2=pi*pigment2_r^2*Qsca_2;
            Cabs_2=pi*pigment2_r^2*Qabs_2;
            alfa_2=f_v_2(j)*Csca_2/V_2;
            beta_2=f_v_2(j)*Cabs_2/V_2+4*pi*k_sio2(i)/lamda(i);
            mu_tot_arr_2(i)=alfa_2+beta_2;
            scat_prob_arr_2(i)=alfa_2/mu_tot_arr_2(i);
            g_arr_2(i)=fonksiyon_2(5);

            if i==length(lamda)
                g_gpu_1=gpuArray(g_arr_1);
                mu_tot_gpu_1=gpuArray(mu_tot_arr_1);
                scat_prob_gpu_1=gpuArray(scat_prob_arr_1);
                g_gpu_2=gpuArray(g_arr_2);
                mu_tot_gpu_2=gpuArray(mu_tot_arr_2);
                scat_prob_gpu_2=gpuArray(scat_prob_arr_2);
                for o=1:length(polar_angle)
                    teta_prime_arr=teta_prime(:,o);
                    sur_reflection_arr=sur_reflection(:,o);
                    teta_prime_gpu=gpuArray(teta_prime_arr);
                    sur_reflection_gpu=gpuArray(sur_reflection_arr);
                    r_no=arrayfun(@monte_carlo,photon_number,sur_reflection_gpu,cosd(teta_prime_gpu),thickness,scat_prob_gpu_1,mu_tot_gpu_1,g_gpu_1,scat_prob_gpu_2,mu_tot_gpu_2,g_gpu_2,n_medium_gpu,k_medium_gpu,n_silver_gpu,k_silver_gpu);
                    ref_lamda(:,j,k,o)=gather(r_no);
                end
            end
        end
        toc
    end
end
abs_lamda=1-ref_lamda;


post_process_fv