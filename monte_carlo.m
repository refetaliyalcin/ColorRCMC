function [r_tot]=monte_carlo(photon_number,s_ref,cos_gecen,h,scat_prob_1,mu_tot_1,g_1,scat_prob_2,mu_tot_2,g_2,n_medium,k_medium,n_subs,k_subs)
r_tot=0;
for i=1:photon_number
    r_no=0;
%     a_no=0;
    Rastgele=rand();
    if Rastgele>s_ref
        %gecer
        alive=1;
    else
        %yansir
        alive=0;
        r_no=1;
    end
    x=0;
    y=0;
    z=0;
    s_x=0;s_y=sqrt(1-cos_gecen*cos_gecen);s_z=cos_gecen; %direction vector
    mu_tot=mu_tot_1+mu_tot_2;
    l_beta=-log(rand())/mu_tot; %ext length
    while alive   
        if (s_z>0)
            l_w = (h - z)/s_z; %distance to lower boundary
        else
            l_w = -z/s_z; %distance to upper boundary
        end
        if l_w<l_beta
            min_index=1;
            min_l=l_w;
        else
            min_index=2;
            min_l=l_beta;
        end
    %         [min_l,min_index]=min([l1 l2]); %select minimum
        x=x+min_l*s_x;
        y=y+min_l*s_y;
        z=z+min_l*s_z;
        if (min_index==1)
    %            disp('hit boundary');
            alive=snell(s_z,n_medium,k_medium,n_subs,k_subs);
            if (alive==0)
                if s_z>0
                    %absorbed in silver
    %                 a_no=1;
                else
                    r_no=1;
                end
            else
                l_beta=l_beta-l_w;
                s_z=-s_z;
            end
        else
            random_no=rand();
            pigment_1_prob=mu_tot_1/mu_tot;
            if random_no<pigment_1_prob
                g_=g_1;
                scat_prob_=scat_prob_1;
            else
                g_=g_2;
                scat_prob_=scat_prob_2;
            end
            random_no=rand();
            if random_no<scat_prob_
    %               disp('scattering');
                [s_x,s_y,s_z]=scatter_mc(g_,s_x,s_y,s_z);
                l_beta=-log(rand())/mu_tot;
            else
    %               disp('absorption');
                alive=0;
            end
        end
    end
    r_tot=r_tot+r_no;
end
r_tot=r_tot/photon_number;