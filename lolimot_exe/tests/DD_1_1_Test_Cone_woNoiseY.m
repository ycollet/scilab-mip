function res = DD_1_1_Test_Cone_woNoiseY(x)
Ei = zeros(1, 30);
Li = zeros(1,30);
 U = 0.0;
 V = 0.0;
 dU1 = 0.0;
 dV1 = 0.0;
 ddU = 0.0;
 ddV = 0.0;

    %%-------------------
    %% List of partitions
    %%-------------------

    %% Partition n 1 partition name : part_uputx - [-3, -2.28][-3, -2.52](0.011376 , -0.69038 , -0.72087)
    %% Partition n 2 partition name : part_ptttx - [0.6, 3][-3, -2.28](0.0919529 , 0.553387 , -0.817762)
    %% Partition n 3 partition name : part_ilmcu - [-3, -1.56][-1.8, 0.12](0.106908 , -0.921427 , -0.331011)
    %% Partition n 4 partition name : part_oerky - [0.6, 1.08][0.6, 2.04](0.0277233 , 0.546941 , 0.830828)
    %% Partition n 5 partition name : part_kgufb - [-3, -2.28][0.12, 1.272](0.0298963 , -0.964738 , 0.245693)
    %% Partition n 6 partition name : part_myhab - [-3, -2.28][1.272, 3](0.0409174 , -0.774302 , 0.62503)
    %% Partition n 7 partition name : part_dykuk - [0.6, 0.696][-2.28, 0.6](0.715037 , -0.0454748 , -0.624699)
    %% Partition n 8 partition name : part_uniqs - [-1.56, -0.264][-1.8, -1.416](0.0804118 , -0.482947 , -0.843472)
    %% Partition n 9 partition name : part_khash - [-2.28, -0.552][0.12, 1.272](0.0671749 , -0.877798 , 0.440825)
    %% Partition n 10 partition name : part_lffwb - [-2.28, -1.704][1.272, 3](0.0347556 , -0.684149 , 0.723978)
    %% Partition n 11 partition name : part_qfwyv - [1.08, 3][-2.28, -1.704](0.0603422 , 0.710476 , -0.687386)
    %% Partition n 12 partition name : part_iqlij - [-0.552, 0.6][0.12, 1.272](0.153484 , 0.0617728 , 0.887686)
    %% Partition n 13 partition name : part_lksqh - [-0.264, 0.6][-1.8, 0.12](0.12561 , 0.262044 , -0.888176)
    %% Partition n 14 partition name : part_yywwy - [-3, -2.28][-2.52, -1.8](0.0228256 , -0.774231 , -0.625282)
    %% Partition n 15 partition name : part_milyr - [1.08, 3][-1.704, -0.3216](0.0651112 , 0.880671 , -0.450174)
    %% Partition n 16 partition name : part_pndif - [-1.704, -1.2432][1.272, 3](0.0294862 , -0.562695 , 0.820673)
    %% Partition n 17 partition name : part_eykcv - [-2.28, -1.704][-2.52, -1.8](0.0107674 , -0.673465 , -0.736364)
    %% Partition n 18 partition name : part_efepf - [0.6, 3][2.04, 3](0.0669785 , 0.566438 , 0.819294)
    %% Partition n 19 partition name : part_rpbuy - [1.08, 3][-0.3216, 0.6](0.0317827 , 0.98981 , 0.0719228)
    %% Partition n 20 partition name : part_nxfxb - [0.696, 1.08][-2.28, -1.704](0.00389737 , 0.40846 , -0.912237)
    %% Partition n 21 partition name : part_rqgwh - [-2.28, -1.128][-3, -2.52](0.0144478 , -0.518567 , -0.85395)
    %% Partition n 22 partition name : part_hdjjx - [-1.2432, -0.13728][1.272, 3](0.0368685 , -0.306911 , 0.945357)
    %% Partition n 23 partition name : part_kyjpx - [0.696, 0.7728][-1.704, 0.6](-0.299738 , 1.44507 , -0.500821)
    %% Partition n 24 partition name : part_gaumw - [-1.704, -1.2432][-2.52, -1.8](0.00999558 , -0.564968 , -0.822611)
    %% Partition n 25 partition name : part_lfnup - [-1.56, -0.264][-1.416, 0.12](0.115359 , -0.784627 , -0.552661)
    %% Partition n 26 partition name : part_qbpbr - [0.7728, 1.08][-1.704, 0.6](0.245741 , 0.777714 , -0.433307)
    %% Partition n 27 partition name : part_octjo - [1.08, 3][0.6, 2.04](0.0598539 , 0.834163 , 0.536205)
    %% Partition n 28 partition name : part_yuqmp - [-1.128, 0.6][-3, -2.52](0.0395103 , -0.093775 , -0.995116)
    %% Partition n 29 partition name : part_keqbg - [-1.2432, 0.6][-2.52, -1.8](0.0862158 , -0.136163 , -0.97762)
    %% Partition n 30 partition name : part_lhngr - [-0.13728, 0.6][1.272, 3](0.0180442 , 0.113578 , 0.989991)

    %%------------------------------------------
    %% Liste of variables name and of boundaries
    %%------------------------------------------

%% X1 - [-3  3]
%% X2 - [-3  3]
%% Measure - [0.042855 4.24264]

    %%------------------------------------------
    %% Pour chaque partition i, on calcule Ei(x)
    %%------------------------------------------

 Ei(1) = exp(-0.5 * ( + ((x(1) - -2.64)^2.0) / 0.0250906 + ((x(2) - -2.76)^2.0) / 0.0111514));
 Ei(2) = exp(-0.5 * ( + ((x(1) - 1.8)^2.0) / 0.278784 + ((x(2) - -2.64)^2.0) / 0.0250906));
 Ei(3) = exp(-0.5 * ( + ((x(1) - -2.28)^2.0) / 0.100362 + ((x(2) - -0.84)^2.0) / 0.178422));
 Ei(4) = exp(-0.5 * ( + ((x(1) - 0.84)^2.0) / 0.0111514 + ((x(2) - 1.32)^2.0) / 0.100362));
 Ei(5) = exp(-0.5 * ( + ((x(1) - -2.64)^2.0) / 0.0250906 + ((x(2) - 0.696)^2.0) / 0.0642318));
 Ei(6) = exp(-0.5 * ( + ((x(1) - -2.64)^2.0) / 0.0250906 + ((x(2) - 2.136)^2.0) / 0.144522));
 Ei(7) = exp(-0.5 * ( + ((x(1) - 0.648)^2.0) / 0.000446055 + ((x(2) - -0.84)^2.0) / 0.401449));
 Ei(8) = exp(-0.5 * ( + ((x(1) - -0.912)^2.0) / 0.0812934 + ((x(2) - -1.608)^2.0) / 0.00713687));
 Ei(9) = exp(-0.5 * ( + ((x(1) - -1.416)^2.0) / 0.144522 + ((x(2) - 0.696)^2.0) / 0.0642318));
 Ei(10) = exp(-0.5 * ( + ((x(1) - -1.992)^2.0) / 0.016058 + ((x(2) - 2.136)^2.0) / 0.144522));
 Ei(11) = exp(-0.5 * ( + ((x(1) - 2.04)^2.0) / 0.178422 + ((x(2) - -1.992)^2.0) / 0.016058));
 Ei(12) = exp(-0.5 * ( + ((x(1) - 0.0240002)^2.0) / 0.0642318 + ((x(2) - 0.696)^2.0) / 0.0642318));
 Ei(13) = exp(-0.5 * ( + ((x(1) - 0.168)^2.0) / 0.0361304 + ((x(2) - -0.84)^2.0) / 0.178422));
 Ei(14) = exp(-0.5 * ( + ((x(1) - -2.64)^2.0) / 0.0250906 + ((x(2) - -2.16)^2.0) / 0.0250906));
 Ei(15) = exp(-0.5 * ( + ((x(1) - 2.04)^2.0) / 0.178422 + ((x(2) - -1.0128)^2.0) / 0.0924939));
 Ei(16) = exp(-0.5 * ( + ((x(1) - -1.4736)^2.0) / 0.0102771 + ((x(2) - 2.136)^2.0) / 0.144522));
 Ei(17) = exp(-0.5 * ( + ((x(1) - -1.992)^2.0) / 0.016058 + ((x(2) - -2.16)^2.0) / 0.0250906));
 Ei(18) = exp(-0.5 * ( + ((x(1) - 1.8)^2.0) / 0.278784 + ((x(2) - 2.52)^2.0) / 0.0446054));
 Ei(19) = exp(-0.5 * ( + ((x(1) - 2.04)^2.0) / 0.178422 + ((x(2) - 0.1392)^2.0) / 0.0411084));
 Ei(20) = exp(-0.5 * ( + ((x(1) - 0.888)^2.0) / 0.00713687 + ((x(2) - -1.992)^2.0) / 0.016058));
 Ei(21) = exp(-0.5 * ( + ((x(1) - -1.704)^2.0) / 0.0642318 + ((x(2) - -2.76)^2.0) / 0.0111514));
 Ei(22) = exp(-0.5 * ( + ((x(1) - -0.69024)^2.0) / 0.0591961 + ((x(2) - 2.136)^2.0) / 0.144522));
 Ei(23) = exp(-0.5 * ( + ((x(1) - 0.7344)^2.0) / 0.000285475 + ((x(2) - -0.552)^2.0) / 0.256927));
 Ei(24) = exp(-0.5 * ( + ((x(1) - -1.4736)^2.0) / 0.0102771 + ((x(2) - -2.16)^2.0) / 0.0250906));
 Ei(25) = exp(-0.5 * ( + ((x(1) - -0.912)^2.0) / 0.0812934 + ((x(2) - -0.648)^2.0) / 0.11419));
 Ei(26) = exp(-0.5 * ( + ((x(1) - 0.9264)^2.0) / 0.0045676 + ((x(2) - -0.552)^2.0) / 0.256927));
 Ei(27) = exp(-0.5 * ( + ((x(1) - 2.04)^2.0) / 0.178422 + ((x(2) - 1.32)^2.0) / 0.100362));
 Ei(28) = exp(-0.5 * ( + ((x(1) - -0.264)^2.0) / 0.144522 + ((x(2) - -2.76)^2.0) / 0.0111514));
 Ei(29) = exp(-0.5 * ( + ((x(1) - -0.3216)^2.0) / 0.164434 + ((x(2) - -2.16)^2.0) / 0.0250906));
 Ei(30) = exp(-0.5 * ( + ((x(1) - 0.23136)^2.0) / 0.0263094 + ((x(2) - 2.136)^2.0) / 0.144522));

    %%------------------------------------------
    %% Pour chaque partition i, on calcule Li(x)
    %%------------------------------------------

    Li(1) = 0.011376 + -0.69038 * x(1) + -0.72087 * x(2);
    Li(2) = 0.0919529 + 0.553387 * x(1) + -0.817762 * x(2);
    Li(3) = 0.106908 + -0.921427 * x(1) + -0.331011 * x(2);
    Li(4) = 0.0277233 + 0.546941 * x(1) + 0.830828 * x(2);
    Li(5) = 0.0298963 + -0.964738 * x(1) + 0.245693 * x(2);
    Li(6) = 0.0409174 + -0.774302 * x(1) + 0.62503 * x(2);
    Li(7) = 0.715037 + -0.0454748 * x(1) + -0.624699 * x(2);
    Li(8) = 0.0804118 + -0.482947 * x(1) + -0.843472 * x(2);
    Li(9) = 0.0671749 + -0.877798 * x(1) + 0.440825 * x(2);
    Li(10) = 0.0347556 + -0.684149 * x(1) + 0.723978 * x(2);
    Li(11) = 0.0603422 + 0.710476 * x(1) + -0.687386 * x(2);
    Li(12) = 0.153484 + 0.0617728 * x(1) + 0.887686 * x(2);
    Li(13) = 0.12561 + 0.262044 * x(1) + -0.888176 * x(2);
    Li(14) = 0.0228256 + -0.774231 * x(1) + -0.625282 * x(2);
    Li(15) = 0.0651112 + 0.880671 * x(1) + -0.450174 * x(2);
    Li(16) = 0.0294862 + -0.562695 * x(1) + 0.820673 * x(2);
    Li(17) = 0.0107674 + -0.673465 * x(1) + -0.736364 * x(2);
    Li(18) = 0.0669785 + 0.566438 * x(1) + 0.819294 * x(2);
    Li(19) = 0.0317827 + 0.98981 * x(1) + 0.0719228 * x(2);
    Li(20) = 0.00389737 + 0.40846 * x(1) + -0.912237 * x(2);
    Li(21) = 0.0144478 + -0.518567 * x(1) + -0.85395 * x(2);
    Li(22) = 0.0368685 + -0.306911 * x(1) + 0.945357 * x(2);
    Li(23) = -0.299738 + 1.44507 * x(1) + -0.500821 * x(2);
    Li(24) = 0.00999558 + -0.564968 * x(1) + -0.822611 * x(2);
    Li(25) = 0.115359 + -0.784627 * x(1) + -0.552661 * x(2);
    Li(26) = 0.245741 + 0.777714 * x(1) + -0.433307 * x(2);
    Li(27) = 0.0598539 + 0.834163 * x(1) + 0.536205 * x(2);
    Li(28) = 0.0395103 + -0.093775 * x(1) + -0.995116 * x(2);
    Li(29) = 0.0862158 + -0.136163 * x(1) + -0.97762 * x(2);
    Li(30) = 0.0180442 + 0.113578 * x(1) + 0.989991 * x(2);

    %%-----------------------------
    %% Computation of the numerator
    %%-----------------------------

 U = 0.0;
 U = U + Li(1) * Ei(1);
 U = U + Li(2) * Ei(2);
 U = U + Li(3) * Ei(3);
 U = U + Li(4) * Ei(4);
 U = U + Li(5) * Ei(5);
 U = U + Li(6) * Ei(6);
 U = U + Li(7) * Ei(7);
 U = U + Li(8) * Ei(8);
 U = U + Li(9) * Ei(9);
 U = U + Li(10) * Ei(10);
 U = U + Li(11) * Ei(11);
 U = U + Li(12) * Ei(12);
 U = U + Li(13) * Ei(13);
 U = U + Li(14) * Ei(14);
 U = U + Li(15) * Ei(15);
 U = U + Li(16) * Ei(16);
 U = U + Li(17) * Ei(17);
 U = U + Li(18) * Ei(18);
 U = U + Li(19) * Ei(19);
 U = U + Li(20) * Ei(20);
 U = U + Li(21) * Ei(21);
 U = U + Li(22) * Ei(22);
 U = U + Li(23) * Ei(23);
 U = U + Li(24) * Ei(24);
 U = U + Li(25) * Ei(25);
 U = U + Li(26) * Ei(26);
 U = U + Li(27) * Ei(27);
 U = U + Li(28) * Ei(28);
 U = U + Li(29) * Ei(29);
 U = U + Li(30) * Ei(30);

    %%-------------------------------
    %% Computation of the denominator
    %%-------------------------------

 V = 0.0;
 V = V + Ei(1);
 V = V + Ei(2);
 V = V + Ei(3);
 V = V + Ei(4);
 V = V + Ei(5);
 V = V + Ei(6);
 V = V + Ei(7);
 V = V + Ei(8);
 V = V + Ei(9);
 V = V + Ei(10);
 V = V + Ei(11);
 V = V + Ei(12);
 V = V + Ei(13);
 V = V + Ei(14);
 V = V + Ei(15);
 V = V + Ei(16);
 V = V + Ei(17);
 V = V + Ei(18);
 V = V + Ei(19);
 V = V + Ei(20);
 V = V + Ei(21);
 V = V + Ei(22);
 V = V + Ei(23);
 V = V + Ei(24);
 V = V + Ei(25);
 V = V + Ei(26);
 V = V + Ei(27);
 V = V + Ei(28);
 V = V + Ei(29);
 V = V + Ei(30);

    %%-------------------------------------------------
    %% Computation of the derivative of the denominator
    %%-------------------------------------------------

 dV1 = 0.0;
 dV1 = dV1 - ((x(2) - -2.76) / 0.0111514) * Ei(1);
 dV1 = dV1 - ((x(2) - -2.64) / 0.0250906) * Ei(2);
 dV1 = dV1 - ((x(2) - -0.84) / 0.178422) * Ei(3);
 dV1 = dV1 - ((x(2) - 1.32) / 0.100362) * Ei(4);
 dV1 = dV1 - ((x(2) - 0.696) / 0.0642318) * Ei(5);
 dV1 = dV1 - ((x(2) - 2.136) / 0.144522) * Ei(6);
 dV1 = dV1 - ((x(2) - -0.84) / 0.401449) * Ei(7);
 dV1 = dV1 - ((x(2) - -1.608) / 0.00713687) * Ei(8);
 dV1 = dV1 - ((x(2) - 0.696) / 0.0642318) * Ei(9);
 dV1 = dV1 - ((x(2) - 2.136) / 0.144522) * Ei(10);
 dV1 = dV1 - ((x(2) - -1.992) / 0.016058) * Ei(11);
 dV1 = dV1 - ((x(2) - 0.696) / 0.0642318) * Ei(12);
 dV1 = dV1 - ((x(2) - -0.84) / 0.178422) * Ei(13);
 dV1 = dV1 - ((x(2) - -2.16) / 0.0250906) * Ei(14);
 dV1 = dV1 - ((x(2) - -1.0128) / 0.0924939) * Ei(15);
 dV1 = dV1 - ((x(2) - 2.136) / 0.144522) * Ei(16);
 dV1 = dV1 - ((x(2) - -2.16) / 0.0250906) * Ei(17);
 dV1 = dV1 - ((x(2) - 2.52) / 0.0446054) * Ei(18);
 dV1 = dV1 - ((x(2) - 0.1392) / 0.0411084) * Ei(19);
 dV1 = dV1 - ((x(2) - -1.992) / 0.016058) * Ei(20);
 dV1 = dV1 - ((x(2) - -2.76) / 0.0111514) * Ei(21);
 dV1 = dV1 - ((x(2) - 2.136) / 0.144522) * Ei(22);
 dV1 = dV1 - ((x(2) - -0.552) / 0.256927) * Ei(23);
 dV1 = dV1 - ((x(2) - -2.16) / 0.0250906) * Ei(24);
 dV1 = dV1 - ((x(2) - -0.648) / 0.11419) * Ei(25);
 dV1 = dV1 - ((x(2) - -0.552) / 0.256927) * Ei(26);
 dV1 = dV1 - ((x(2) - 1.32) / 0.100362) * Ei(27);
 dV1 = dV1 - ((x(2) - -2.76) / 0.0111514) * Ei(28);
 dV1 = dV1 - ((x(2) - -2.16) / 0.0250906) * Ei(29);
 dV1 = dV1 - ((x(2) - 2.136) / 0.144522) * Ei(30);

    %%-----------------------------------------------
    %% Computation of the derivative of the numerator
    %%-----------------------------------------------

 dU1 = 0.0;
 dU1 = dU1 + (-0.72087 - ((x(2) - -2.76) / 0.0111514) * Li(1)) * Ei(1);
 dU1 = dU1 + (-0.817762 - ((x(2) - -2.64) / 0.0250906) * Li(2)) * Ei(2);
 dU1 = dU1 + (-0.331011 - ((x(2) - -0.84) / 0.178422) * Li(3)) * Ei(3);
 dU1 = dU1 + (0.830828 - ((x(2) - 1.32) / 0.100362) * Li(4)) * Ei(4);
 dU1 = dU1 + (0.245693 - ((x(2) - 0.696) / 0.0642318) * Li(5)) * Ei(5);
 dU1 = dU1 + (0.62503 - ((x(2) - 2.136) / 0.144522) * Li(6)) * Ei(6);
 dU1 = dU1 + (-0.624699 - ((x(2) - -0.84) / 0.401449) * Li(7)) * Ei(7);
 dU1 = dU1 + (-0.843472 - ((x(2) - -1.608) / 0.00713687) * Li(8)) * Ei(8);
 dU1 = dU1 + (0.440825 - ((x(2) - 0.696) / 0.0642318) * Li(9)) * Ei(9);
 dU1 = dU1 + (0.723978 - ((x(2) - 2.136) / 0.144522) * Li(10)) * Ei(10);
 dU1 = dU1 + (-0.687386 - ((x(2) - -1.992) / 0.016058) * Li(11)) * Ei(11);
 dU1 = dU1 + (0.887686 - ((x(2) - 0.696) / 0.0642318) * Li(12)) * Ei(12);
 dU1 = dU1 + (-0.888176 - ((x(2) - -0.84) / 0.178422) * Li(13)) * Ei(13);
 dU1 = dU1 + (-0.625282 - ((x(2) - -2.16) / 0.0250906) * Li(14)) * Ei(14);
 dU1 = dU1 + (-0.450174 - ((x(2) - -1.0128) / 0.0924939) * Li(15)) * Ei(15);
 dU1 = dU1 + (0.820673 - ((x(2) - 2.136) / 0.144522) * Li(16)) * Ei(16);
 dU1 = dU1 + (-0.736364 - ((x(2) - -2.16) / 0.0250906) * Li(17)) * Ei(17);
 dU1 = dU1 + (0.819294 - ((x(2) - 2.52) / 0.0446054) * Li(18)) * Ei(18);
 dU1 = dU1 + (0.0719228 - ((x(2) - 0.1392) / 0.0411084) * Li(19)) * Ei(19);
 dU1 = dU1 + (-0.912237 - ((x(2) - -1.992) / 0.016058) * Li(20)) * Ei(20);
 dU1 = dU1 + (-0.85395 - ((x(2) - -2.76) / 0.0111514) * Li(21)) * Ei(21);
 dU1 = dU1 + (0.945357 - ((x(2) - 2.136) / 0.144522) * Li(22)) * Ei(22);
 dU1 = dU1 + (-0.500821 - ((x(2) - -0.552) / 0.256927) * Li(23)) * Ei(23);
 dU1 = dU1 + (-0.822611 - ((x(2) - -2.16) / 0.0250906) * Li(24)) * Ei(24);
 dU1 = dU1 + (-0.552661 - ((x(2) - -0.648) / 0.11419) * Li(25)) * Ei(25);
 dU1 = dU1 + (-0.433307 - ((x(2) - -0.552) / 0.256927) * Li(26)) * Ei(26);
 dU1 = dU1 + (0.536205 - ((x(2) - 1.32) / 0.100362) * Li(27)) * Ei(27);
 dU1 = dU1 + (-0.995116 - ((x(2) - -2.76) / 0.0111514) * Li(28)) * Ei(28);
 dU1 = dU1 + (-0.97762 - ((x(2) - -2.16) / 0.0250906) * Li(29)) * Ei(29);
 dU1 = dU1 + (0.989991 - ((x(2) - 2.136) / 0.144522) * Li(30)) * Ei(30);

    %%--------------------------------------------------------
    %% Computation of the second derivative of the denominator
    %%--------------------------------------------------------

 ddV = 0.0;
 ddV = ddV - (1.0 / 0.0111514 - ((x(2) - -2.76) / 0.0111514)^2.0) * Ei(1);
 ddV = ddV - (1.0 / 0.0250906 - ((x(2) - -2.64) / 0.0250906)^2.0) * Ei(2);
 ddV = ddV - (1.0 / 0.178422 - ((x(2) - -0.84) / 0.178422)^2.0) * Ei(3);
 ddV = ddV - (1.0 / 0.100362 - ((x(2) - 1.32) / 0.100362)^2.0) * Ei(4);
 ddV = ddV - (1.0 / 0.0642318 - ((x(2) - 0.696) / 0.0642318)^2.0) * Ei(5);
 ddV = ddV - (1.0 / 0.144522 - ((x(2) - 2.136) / 0.144522)^2.0) * Ei(6);
 ddV = ddV - (1.0 / 0.401449 - ((x(2) - -0.84) / 0.401449)^2.0) * Ei(7);
 ddV = ddV - (1.0 / 0.00713687 - ((x(2) - -1.608) / 0.00713687)^2.0) * Ei(8);
 ddV = ddV - (1.0 / 0.0642318 - ((x(2) - 0.696) / 0.0642318)^2.0) * Ei(9);
 ddV = ddV - (1.0 / 0.144522 - ((x(2) - 2.136) / 0.144522)^2.0) * Ei(10);
 ddV = ddV - (1.0 / 0.016058 - ((x(2) - -1.992) / 0.016058)^2.0) * Ei(11);
 ddV = ddV - (1.0 / 0.0642318 - ((x(2) - 0.696) / 0.0642318)^2.0) * Ei(12);
 ddV = ddV - (1.0 / 0.178422 - ((x(2) - -0.84) / 0.178422)^2.0) * Ei(13);
 ddV = ddV - (1.0 / 0.0250906 - ((x(2) - -2.16) / 0.0250906)^2.0) * Ei(14);
 ddV = ddV - (1.0 / 0.0924939 - ((x(2) - -1.0128) / 0.0924939)^2.0) * Ei(15);
 ddV = ddV - (1.0 / 0.144522 - ((x(2) - 2.136) / 0.144522)^2.0) * Ei(16);
 ddV = ddV - (1.0 / 0.0250906 - ((x(2) - -2.16) / 0.0250906)^2.0) * Ei(17);
 ddV = ddV - (1.0 / 0.0446054 - ((x(2) - 2.52) / 0.0446054)^2.0) * Ei(18);
 ddV = ddV - (1.0 / 0.0411084 - ((x(2) - 0.1392) / 0.0411084)^2.0) * Ei(19);
 ddV = ddV - (1.0 / 0.016058 - ((x(2) - -1.992) / 0.016058)^2.0) * Ei(20);
 ddV = ddV - (1.0 / 0.0111514 - ((x(2) - -2.76) / 0.0111514)^2.0) * Ei(21);
 ddV = ddV - (1.0 / 0.144522 - ((x(2) - 2.136) / 0.144522)^2.0) * Ei(22);
 ddV = ddV - (1.0 / 0.256927 - ((x(2) - -0.552) / 0.256927)^2.0) * Ei(23);
 ddV = ddV - (1.0 / 0.0250906 - ((x(2) - -2.16) / 0.0250906)^2.0) * Ei(24);
 ddV = ddV - (1.0 / 0.11419 - ((x(2) - -0.648) / 0.11419)^2.0) * Ei(25);
 ddV = ddV - (1.0 / 0.256927 - ((x(2) - -0.552) / 0.256927)^2.0) * Ei(26);
 ddV = ddV - (1.0 / 0.100362 - ((x(2) - 1.32) / 0.100362)^2.0) * Ei(27);
 ddV = ddV - (1.0 / 0.0111514 - ((x(2) - -2.76) / 0.0111514)^2.0) * Ei(28);
 ddV = ddV - (1.0 / 0.0250906 - ((x(2) - -2.16) / 0.0250906)^2.0) * Ei(29);
 ddV = ddV - (1.0 / 0.144522 - ((x(2) - 2.136) / 0.144522)^2.0) * Ei(30);

    %%------------------------------------------------------
    %% Computation of the second derivative of the numerator
    %%------------------------------------------------------

 ddU = 0.0;
 ddU = ddU - (Li(1) / 0.0111514 + 2 * -0.72087 * ((x(2) - -2.76) / 0.0111514) - Li(1) * ((x(2) - -2.76) / 0.0111514)^2.0) * Ei(1);
 ddU = ddU - (Li(2) / 0.0250906 + 2 * -0.817762 * ((x(2) - -2.64) / 0.0250906) - Li(2) * ((x(2) - -2.64) / 0.0250906)^2.0) * Ei(2);
 ddU = ddU - (Li(3) / 0.178422 + 2 * -0.331011 * ((x(2) - -0.84) / 0.178422) - Li(3) * ((x(2) - -0.84) / 0.178422)^2.0) * Ei(3);
 ddU = ddU - (Li(4) / 0.100362 + 2 * 0.830828 * ((x(2) - 1.32) / 0.100362) - Li(4) * ((x(2) - 1.32) / 0.100362)^2.0) * Ei(4);
 ddU = ddU - (Li(5) / 0.0642318 + 2 * 0.245693 * ((x(2) - 0.696) / 0.0642318) - Li(5) * ((x(2) - 0.696) / 0.0642318)^2.0) * Ei(5);
 ddU = ddU - (Li(6) / 0.144522 + 2 * 0.62503 * ((x(2) - 2.136) / 0.144522) - Li(6) * ((x(2) - 2.136) / 0.144522)^2.0) * Ei(6);
 ddU = ddU - (Li(7) / 0.401449 + 2 * -0.624699 * ((x(2) - -0.84) / 0.401449) - Li(7) * ((x(2) - -0.84) / 0.401449)^2.0) * Ei(7);
 ddU = ddU - (Li(8) / 0.00713687 + 2 * -0.843472 * ((x(2) - -1.608) / 0.00713687) - Li(8) * ((x(2) - -1.608) / 0.00713687)^2.0) * Ei(8);
 ddU = ddU - (Li(9) / 0.0642318 + 2 * 0.440825 * ((x(2) - 0.696) / 0.0642318) - Li(9) * ((x(2) - 0.696) / 0.0642318)^2.0) * Ei(9);
 ddU = ddU - (Li(10) / 0.144522 + 2 * 0.723978 * ((x(2) - 2.136) / 0.144522) - Li(10) * ((x(2) - 2.136) / 0.144522)^2.0) * Ei(10);
 ddU = ddU - (Li(11) / 0.016058 + 2 * -0.687386 * ((x(2) - -1.992) / 0.016058) - Li(11) * ((x(2) - -1.992) / 0.016058)^2.0) * Ei(11);
 ddU = ddU - (Li(12) / 0.0642318 + 2 * 0.887686 * ((x(2) - 0.696) / 0.0642318) - Li(12) * ((x(2) - 0.696) / 0.0642318)^2.0) * Ei(12);
 ddU = ddU - (Li(13) / 0.178422 + 2 * -0.888176 * ((x(2) - -0.84) / 0.178422) - Li(13) * ((x(2) - -0.84) / 0.178422)^2.0) * Ei(13);
 ddU = ddU - (Li(14) / 0.0250906 + 2 * -0.625282 * ((x(2) - -2.16) / 0.0250906) - Li(14) * ((x(2) - -2.16) / 0.0250906)^2.0) * Ei(14);
 ddU = ddU - (Li(15) / 0.0924939 + 2 * -0.450174 * ((x(2) - -1.0128) / 0.0924939) - Li(15) * ((x(2) - -1.0128) / 0.0924939)^2.0) * Ei(15);
 ddU = ddU - (Li(16) / 0.144522 + 2 * 0.820673 * ((x(2) - 2.136) / 0.144522) - Li(16) * ((x(2) - 2.136) / 0.144522)^2.0) * Ei(16);
 ddU = ddU - (Li(17) / 0.0250906 + 2 * -0.736364 * ((x(2) - -2.16) / 0.0250906) - Li(17) * ((x(2) - -2.16) / 0.0250906)^2.0) * Ei(17);
 ddU = ddU - (Li(18) / 0.0446054 + 2 * 0.819294 * ((x(2) - 2.52) / 0.0446054) - Li(18) * ((x(2) - 2.52) / 0.0446054)^2.0) * Ei(18);
 ddU = ddU - (Li(19) / 0.0411084 + 2 * 0.0719228 * ((x(2) - 0.1392) / 0.0411084) - Li(19) * ((x(2) - 0.1392) / 0.0411084)^2.0) * Ei(19);
 ddU = ddU - (Li(20) / 0.016058 + 2 * -0.912237 * ((x(2) - -1.992) / 0.016058) - Li(20) * ((x(2) - -1.992) / 0.016058)^2.0) * Ei(20);
 ddU = ddU - (Li(21) / 0.0111514 + 2 * -0.85395 * ((x(2) - -2.76) / 0.0111514) - Li(21) * ((x(2) - -2.76) / 0.0111514)^2.0) * Ei(21);
 ddU = ddU - (Li(22) / 0.144522 + 2 * 0.945357 * ((x(2) - 2.136) / 0.144522) - Li(22) * ((x(2) - 2.136) / 0.144522)^2.0) * Ei(22);
 ddU = ddU - (Li(23) / 0.256927 + 2 * -0.500821 * ((x(2) - -0.552) / 0.256927) - Li(23) * ((x(2) - -0.552) / 0.256927)^2.0) * Ei(23);
 ddU = ddU - (Li(24) / 0.0250906 + 2 * -0.822611 * ((x(2) - -2.16) / 0.0250906) - Li(24) * ((x(2) - -2.16) / 0.0250906)^2.0) * Ei(24);
 ddU = ddU - (Li(25) / 0.11419 + 2 * -0.552661 * ((x(2) - -0.648) / 0.11419) - Li(25) * ((x(2) - -0.648) / 0.11419)^2.0) * Ei(25);
 ddU = ddU - (Li(26) / 0.256927 + 2 * -0.433307 * ((x(2) - -0.552) / 0.256927) - Li(26) * ((x(2) - -0.552) / 0.256927)^2.0) * Ei(26);
 ddU = ddU - (Li(27) / 0.100362 + 2 * 0.536205 * ((x(2) - 1.32) / 0.100362) - Li(27) * ((x(2) - 1.32) / 0.100362)^2.0) * Ei(27);
 ddU = ddU - (Li(28) / 0.0111514 + 2 * -0.995116 * ((x(2) - -2.76) / 0.0111514) - Li(28) * ((x(2) - -2.76) / 0.0111514)^2.0) * Ei(28);
 ddU = ddU - (Li(29) / 0.0250906 + 2 * -0.97762 * ((x(2) - -2.16) / 0.0250906) - Li(29) * ((x(2) - -2.16) / 0.0250906)^2.0) * Ei(29);
 ddU = ddU - (Li(30) / 0.144522 + 2 * 0.989991 * ((x(2) - 2.136) / 0.144522) - Li(30) * ((x(2) - 2.136) / 0.144522)^2.0) * Ei(30);
  clear Ei;
  clear Li;
    res = (((ddU*V - ddV*U)*V*V - (dU1*V-dV1*U)*2*V*dV1)/(V*V*V*V));
