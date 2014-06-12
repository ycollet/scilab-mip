function res = DD_0_1_SinusY(x)
Ei = zeros(1, 30);
Li = zeros(1,30);
 U = 0.0;
 V = 0.0;
 dU1 = 0.0;
 dV1 = 0.0;
 dU2 = 0.0;
 dV2 = 0.0;
 ddU = 0.0;
 ddV = 0.0;

    %%-------------------
    %% List of partitions
    %%-------------------

    %% Partition n 1 partition name : part_uputx - [0, 111.111][-1, 1](0.0108676 , 0 , 0.991035)
    %% Partition n 2 partition name : part_ptttx - [333.333, 407.407][-1, -0.333333](-0.031786 , 0 , 0.967875)
    %% Partition n 3 partition name : part_ilmcu - [333.333, 555.556][-0.333333, 0.111111](0.012677 , 0 , 1.00248)
    %% Partition n 4 partition name : part_oerky - [333.333, 555.556][0.111111, 1](0.0130032 , 0 , 0.993099)
    %% Partition n 5 partition name : part_kgufb - [111.111, 185.185][-1, 1](-0.0341052 , 0 , 1.03233)
    %% Partition n 6 partition name : part_myhab - [555.556, 703.704][-1, -0.777778](-0.03914 , 0 , 0.953868)
    %% Partition n 7 partition name : part_dykuk - [555.556, 703.704][-0.777778, -0.62963](-0.0237299 , 0 , 0.974689)
    %% Partition n 8 partition name : part_uniqs - [555.556, 703.704][-0.333333, -0.185185](-0.009378 , 0 , 0.972573)
    %% Partition n 9 partition name : part_khash - [555.556, 604.938][0.111111, 1](0.0228878 , 0 , 0.980401)
    %% Partition n 10 partition name : part_lffwb - [555.556, 1000][-0.62963, -0.596708](0.0551995 , 0 , 1.09248)
    %% Partition n 11 partition name : part_qfwyv - [555.556, 703.704][-0.185185, 0.111111](-0.00118627 , 0 , 0.990434)
    %% Partition n 12 partition name : part_iqlij - [555.556, 703.704][-0.530864, -0.333333](0.00028682 , 0 , 1.0109)
    %% Partition n 13 partition name : part_lksqh - [703.704, 802.469][-0.530864, -0.333333](-0.0138932 , 0 , 0.993938)
    %% Partition n 14 partition name : part_yywwy - [703.704, 802.469][-0.185185, 0.111111](-0.0125165 , 0 , 0.999673)
    %% Partition n 15 partition name : part_milyr - [802.469, 1000][-0.530864, -0.333333](0.0156713 , 0 , 1.01348)
    %% Partition n 16 partition name : part_pndif - [703.704, 1000][0.111111, 1](-0.013556 , 0 , 1.00644)
    %% Partition n 17 partition name : part_eykcv - [802.469, 1000][-0.185185, 0.111111](0.0121893 , 0 , 0.99922)
    %% Partition n 18 partition name : part_efepf - [604.938, 703.704][0.111111, 1](-0.0248321 , 0 , 1.02342)
    %% Partition n 19 partition name : part_rpbuy - [407.407, 555.556][-1, -0.333333](0.0169308 , 0 , 1.01202)
    %% Partition n 20 partition name : part_nxfxb - [703.704, 802.469][-0.333333, -0.185185](-0.0129954 , 0 , 0.996558)
    %% Partition n 21 partition name : part_rqgwh - [703.704, 802.469][-1, -0.777778](-0.0239252 , 0 , 0.979475)
    %% Partition n 22 partition name : part_hdjjx - [802.469, 868.313][-1, -0.777778](-0.032689 , 0 , 0.969202)
    %% Partition n 23 partition name : part_kyjpx - [802.469, 1000][-0.333333, -0.185185](0.0133166 , 0 , 1.0064)
    %% Partition n 24 partition name : part_gaumw - [703.704, 1000][-0.777778, -0.728395](-0.0215838 , 0 , 0.972939)
    %% Partition n 25 partition name : part_lfnup - [868.313, 1000][-1, -0.777778](0.0369446 , 0 , 1.03568)
    %% Partition n 26 partition name : part_qbpbr - [703.704, 802.469][-0.728395, -0.62963](-0.017022 , 0 , 0.988424)
    %% Partition n 27 partition name : part_octjo - [185.185, 333.333][-1, 1](-0.00869546 , 0 , 1.00141)
    %% Partition n 28 partition name : part_yuqmp - [802.469, 1000][-0.728395, -0.62963](0.0210843 , 0 , 1.02471)
    %% Partition n 29 partition name : part_keqbg - [555.556, 1000][-0.596708, -0.57476](-0.00266084 , 0 , 0.996492)
    %% Partition n 30 partition name : part_lhngr - [555.556, 1000][-0.57476, -0.530864](-0.0332969 , 0 , 0.942638)

    %%------------------------------------------
    %% Liste of variables name and of boundaries
    %%------------------------------------------

%% T - [0  1000]
%% S - [-1  1]
%% Measure - [-1 1]

    %%------------------------------------------
    %% Pour chaque partition i, on calcule Ei(x)
    %%------------------------------------------

 Ei(1) = exp(-0.5 * ( + ((x(1) - 55.5556)^2.0) / 1344.44 + ((x(2) - 0)^2.0) / 0.4356));
 Ei(2) = exp(-0.5 * ( + ((x(1) - 370.37)^2.0) / 597.531 + ((x(2) - -0.666667)^2.0) / 0.0484));
 Ei(3) = exp(-0.5 * ( + ((x(1) - 444.444)^2.0) / 5377.78 + ((x(2) - -0.111111)^2.0) / 0.0215111));
 Ei(4) = exp(-0.5 * ( + ((x(1) - 444.444)^2.0) / 5377.78 + ((x(2) - 0.555556)^2.0) / 0.0860445));
 Ei(5) = exp(-0.5 * ( + ((x(1) - 148.148)^2.0) / 597.531 + ((x(2) - 0)^2.0) / 0.4356));
 Ei(6) = exp(-0.5 * ( + ((x(1) - 629.63)^2.0) / 2390.12 + ((x(2) - -0.888889)^2.0) / 0.00537778));
 Ei(7) = exp(-0.5 * ( + ((x(1) - 629.63)^2.0) / 2390.12 + ((x(2) - -0.703704)^2.0) / 0.00239012));
 Ei(8) = exp(-0.5 * ( + ((x(1) - 629.63)^2.0) / 2390.12 + ((x(2) - -0.259259)^2.0) / 0.00239012));
 Ei(9) = exp(-0.5 * ( + ((x(1) - 580.247)^2.0) / 265.569 + ((x(2) - 0.555556)^2.0) / 0.0860445));
 Ei(10) = exp(-0.5 * ( + ((x(1) - 777.778)^2.0) / 21511.1 + ((x(2) - -0.613169)^2.0) / 0.000118031));
 Ei(11) = exp(-0.5 * ( + ((x(1) - 629.63)^2.0) / 2390.12 + ((x(2) - -0.037037)^2.0) / 0.00956049));
 Ei(12) = exp(-0.5 * ( + ((x(1) - 629.63)^2.0) / 2390.12 + ((x(2) - -0.432099)^2.0) / 0.00424911));
 Ei(13) = exp(-0.5 * ( + ((x(1) - 753.086)^2.0) / 1062.28 + ((x(2) - -0.432099)^2.0) / 0.00424911));
 Ei(14) = exp(-0.5 * ( + ((x(1) - 753.086)^2.0) / 1062.28 + ((x(2) - -0.037037)^2.0) / 0.00956049));
 Ei(15) = exp(-0.5 * ( + ((x(1) - 901.235)^2.0) / 4249.11 + ((x(2) - -0.432099)^2.0) / 0.00424911));
 Ei(16) = exp(-0.5 * ( + ((x(1) - 851.852)^2.0) / 9560.5 + ((x(2) - 0.555556)^2.0) / 0.0860445));
 Ei(17) = exp(-0.5 * ( + ((x(1) - 901.235)^2.0) / 4249.11 + ((x(2) - -0.037037)^2.0) / 0.00956049));
 Ei(18) = exp(-0.5 * ( + ((x(1) - 654.321)^2.0) / 1062.28 + ((x(2) - 0.555556)^2.0) / 0.0860445));
 Ei(19) = exp(-0.5 * ( + ((x(1) - 481.481)^2.0) / 2390.12 + ((x(2) - -0.666667)^2.0) / 0.0484));
 Ei(20) = exp(-0.5 * ( + ((x(1) - 753.086)^2.0) / 1062.28 + ((x(2) - -0.259259)^2.0) / 0.00239012));
 Ei(21) = exp(-0.5 * ( + ((x(1) - 753.086)^2.0) / 1062.28 + ((x(2) - -0.888889)^2.0) / 0.00537778));
 Ei(22) = exp(-0.5 * ( + ((x(1) - 835.391)^2.0) / 472.123 + ((x(2) - -0.888889)^2.0) / 0.00537778));
 Ei(23) = exp(-0.5 * ( + ((x(1) - 901.235)^2.0) / 4249.11 + ((x(2) - -0.259259)^2.0) / 0.00239012));
 Ei(24) = exp(-0.5 * ( + ((x(1) - 851.852)^2.0) / 9560.5 + ((x(2) - -0.753086)^2.0) / 0.00026557));
 Ei(25) = exp(-0.5 * ( + ((x(1) - 934.156)^2.0) / 1888.49 + ((x(2) - -0.888889)^2.0) / 0.00537778));
 Ei(26) = exp(-0.5 * ( + ((x(1) - 753.086)^2.0) / 1062.28 + ((x(2) - -0.679012)^2.0) / 0.00106228));
 Ei(27) = exp(-0.5 * ( + ((x(1) - 259.259)^2.0) / 2390.12 + ((x(2) - 0)^2.0) / 0.4356));
 Ei(28) = exp(-0.5 * ( + ((x(1) - 901.235)^2.0) / 4249.11 + ((x(2) - -0.679012)^2.0) / 0.00106228));
 Ei(29) = exp(-0.5 * ( + ((x(1) - 777.778)^2.0) / 21511.1 + ((x(2) - -0.585734)^2.0) / 5.24581e-05));
 Ei(30) = exp(-0.5 * ( + ((x(1) - 777.778)^2.0) / 21511.1 + ((x(2) - -0.552812)^2.0) / 0.000209833));

    %%------------------------------------------
    %% Pour chaque partition i, on calcule Li(x)
    %%------------------------------------------

    Li(1) = 0.0108676 + 0 * x(1) + 0.991035 * x(2);
    Li(2) = -0.031786 + 0 * x(1) + 0.967875 * x(2);
    Li(3) = 0.012677 + 0 * x(1) + 1.00248 * x(2);
    Li(4) = 0.0130032 + 0 * x(1) + 0.993099 * x(2);
    Li(5) = -0.0341052 + 0 * x(1) + 1.03233 * x(2);
    Li(6) = -0.03914 + 0 * x(1) + 0.953868 * x(2);
    Li(7) = -0.0237299 + 0 * x(1) + 0.974689 * x(2);
    Li(8) = -0.009378 + 0 * x(1) + 0.972573 * x(2);
    Li(9) = 0.0228878 + 0 * x(1) + 0.980401 * x(2);
    Li(10) = 0.0551995 + 0 * x(1) + 1.09248 * x(2);
    Li(11) = -0.00118627 + 0 * x(1) + 0.990434 * x(2);
    Li(12) = 0.00028682 + 0 * x(1) + 1.0109 * x(2);
    Li(13) = -0.0138932 + 0 * x(1) + 0.993938 * x(2);
    Li(14) = -0.0125165 + 0 * x(1) + 0.999673 * x(2);
    Li(15) = 0.0156713 + 0 * x(1) + 1.01348 * x(2);
    Li(16) = -0.013556 + 0 * x(1) + 1.00644 * x(2);
    Li(17) = 0.0121893 + 0 * x(1) + 0.99922 * x(2);
    Li(18) = -0.0248321 + 0 * x(1) + 1.02342 * x(2);
    Li(19) = 0.0169308 + 0 * x(1) + 1.01202 * x(2);
    Li(20) = -0.0129954 + 0 * x(1) + 0.996558 * x(2);
    Li(21) = -0.0239252 + 0 * x(1) + 0.979475 * x(2);
    Li(22) = -0.032689 + 0 * x(1) + 0.969202 * x(2);
    Li(23) = 0.0133166 + 0 * x(1) + 1.0064 * x(2);
    Li(24) = -0.0215838 + 0 * x(1) + 0.972939 * x(2);
    Li(25) = 0.0369446 + 0 * x(1) + 1.03568 * x(2);
    Li(26) = -0.017022 + 0 * x(1) + 0.988424 * x(2);
    Li(27) = -0.00869546 + 0 * x(1) + 1.00141 * x(2);
    Li(28) = 0.0210843 + 0 * x(1) + 1.02471 * x(2);
    Li(29) = -0.00266084 + 0 * x(1) + 0.996492 * x(2);
    Li(30) = -0.0332969 + 0 * x(1) + 0.942638 * x(2);

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
 dV1 = dV1 - ((x(1) - 55.5556) / 1344.44) * Ei(1);
 dV1 = dV1 - ((x(1) - 370.37) / 597.531) * Ei(2);
 dV1 = dV1 - ((x(1) - 444.444) / 5377.78) * Ei(3);
 dV1 = dV1 - ((x(1) - 444.444) / 5377.78) * Ei(4);
 dV1 = dV1 - ((x(1) - 148.148) / 597.531) * Ei(5);
 dV1 = dV1 - ((x(1) - 629.63) / 2390.12) * Ei(6);
 dV1 = dV1 - ((x(1) - 629.63) / 2390.12) * Ei(7);
 dV1 = dV1 - ((x(1) - 629.63) / 2390.12) * Ei(8);
 dV1 = dV1 - ((x(1) - 580.247) / 265.569) * Ei(9);
 dV1 = dV1 - ((x(1) - 777.778) / 21511.1) * Ei(10);
 dV1 = dV1 - ((x(1) - 629.63) / 2390.12) * Ei(11);
 dV1 = dV1 - ((x(1) - 629.63) / 2390.12) * Ei(12);
 dV1 = dV1 - ((x(1) - 753.086) / 1062.28) * Ei(13);
 dV1 = dV1 - ((x(1) - 753.086) / 1062.28) * Ei(14);
 dV1 = dV1 - ((x(1) - 901.235) / 4249.11) * Ei(15);
 dV1 = dV1 - ((x(1) - 851.852) / 9560.5) * Ei(16);
 dV1 = dV1 - ((x(1) - 901.235) / 4249.11) * Ei(17);
 dV1 = dV1 - ((x(1) - 654.321) / 1062.28) * Ei(18);
 dV1 = dV1 - ((x(1) - 481.481) / 2390.12) * Ei(19);
 dV1 = dV1 - ((x(1) - 753.086) / 1062.28) * Ei(20);
 dV1 = dV1 - ((x(1) - 753.086) / 1062.28) * Ei(21);
 dV1 = dV1 - ((x(1) - 835.391) / 472.123) * Ei(22);
 dV1 = dV1 - ((x(1) - 901.235) / 4249.11) * Ei(23);
 dV1 = dV1 - ((x(1) - 851.852) / 9560.5) * Ei(24);
 dV1 = dV1 - ((x(1) - 934.156) / 1888.49) * Ei(25);
 dV1 = dV1 - ((x(1) - 753.086) / 1062.28) * Ei(26);
 dV1 = dV1 - ((x(1) - 259.259) / 2390.12) * Ei(27);
 dV1 = dV1 - ((x(1) - 901.235) / 4249.11) * Ei(28);
 dV1 = dV1 - ((x(1) - 777.778) / 21511.1) * Ei(29);
 dV1 = dV1 - ((x(1) - 777.778) / 21511.1) * Ei(30);

    %%-----------------------------------------------
    %% Computation of the derivative of the numerator
    %%-----------------------------------------------

 dU1 = 0.0;
 dU1 = dU1 + (0 - ((x(1) - 55.5556) / 1344.44) * Li(1)) * Ei(1);
 dU1 = dU1 + (0 - ((x(1) - 370.37) / 597.531) * Li(2)) * Ei(2);
 dU1 = dU1 + (0 - ((x(1) - 444.444) / 5377.78) * Li(3)) * Ei(3);
 dU1 = dU1 + (0 - ((x(1) - 444.444) / 5377.78) * Li(4)) * Ei(4);
 dU1 = dU1 + (0 - ((x(1) - 148.148) / 597.531) * Li(5)) * Ei(5);
 dU1 = dU1 + (0 - ((x(1) - 629.63) / 2390.12) * Li(6)) * Ei(6);
 dU1 = dU1 + (0 - ((x(1) - 629.63) / 2390.12) * Li(7)) * Ei(7);
 dU1 = dU1 + (0 - ((x(1) - 629.63) / 2390.12) * Li(8)) * Ei(8);
 dU1 = dU1 + (0 - ((x(1) - 580.247) / 265.569) * Li(9)) * Ei(9);
 dU1 = dU1 + (0 - ((x(1) - 777.778) / 21511.1) * Li(10)) * Ei(10);
 dU1 = dU1 + (0 - ((x(1) - 629.63) / 2390.12) * Li(11)) * Ei(11);
 dU1 = dU1 + (0 - ((x(1) - 629.63) / 2390.12) * Li(12)) * Ei(12);
 dU1 = dU1 + (0 - ((x(1) - 753.086) / 1062.28) * Li(13)) * Ei(13);
 dU1 = dU1 + (0 - ((x(1) - 753.086) / 1062.28) * Li(14)) * Ei(14);
 dU1 = dU1 + (0 - ((x(1) - 901.235) / 4249.11) * Li(15)) * Ei(15);
 dU1 = dU1 + (0 - ((x(1) - 851.852) / 9560.5) * Li(16)) * Ei(16);
 dU1 = dU1 + (0 - ((x(1) - 901.235) / 4249.11) * Li(17)) * Ei(17);
 dU1 = dU1 + (0 - ((x(1) - 654.321) / 1062.28) * Li(18)) * Ei(18);
 dU1 = dU1 + (0 - ((x(1) - 481.481) / 2390.12) * Li(19)) * Ei(19);
 dU1 = dU1 + (0 - ((x(1) - 753.086) / 1062.28) * Li(20)) * Ei(20);
 dU1 = dU1 + (0 - ((x(1) - 753.086) / 1062.28) * Li(21)) * Ei(21);
 dU1 = dU1 + (0 - ((x(1) - 835.391) / 472.123) * Li(22)) * Ei(22);
 dU1 = dU1 + (0 - ((x(1) - 901.235) / 4249.11) * Li(23)) * Ei(23);
 dU1 = dU1 + (0 - ((x(1) - 851.852) / 9560.5) * Li(24)) * Ei(24);
 dU1 = dU1 + (0 - ((x(1) - 934.156) / 1888.49) * Li(25)) * Ei(25);
 dU1 = dU1 + (0 - ((x(1) - 753.086) / 1062.28) * Li(26)) * Ei(26);
 dU1 = dU1 + (0 - ((x(1) - 259.259) / 2390.12) * Li(27)) * Ei(27);
 dU1 = dU1 + (0 - ((x(1) - 901.235) / 4249.11) * Li(28)) * Ei(28);
 dU1 = dU1 + (0 - ((x(1) - 777.778) / 21511.1) * Li(29)) * Ei(29);
 dU1 = dU1 + (0 - ((x(1) - 777.778) / 21511.1) * Li(30)) * Ei(30);

    %%-------------------------------------------------
    %% Computation of the derivative of the denominator
    %%-------------------------------------------------

 dV2 = 0.0;
 dV2 = dV2 - ((x(2) - 0) / 0.4356) * Ei(1);
 dV2 = dV2 - ((x(2) - -0.666667) / 0.0484) * Ei(2);
 dV2 = dV2 - ((x(2) - -0.111111) / 0.0215111) * Ei(3);
 dV2 = dV2 - ((x(2) - 0.555556) / 0.0860445) * Ei(4);
 dV2 = dV2 - ((x(2) - 0) / 0.4356) * Ei(5);
 dV2 = dV2 - ((x(2) - -0.888889) / 0.00537778) * Ei(6);
 dV2 = dV2 - ((x(2) - -0.703704) / 0.00239012) * Ei(7);
 dV2 = dV2 - ((x(2) - -0.259259) / 0.00239012) * Ei(8);
 dV2 = dV2 - ((x(2) - 0.555556) / 0.0860445) * Ei(9);
 dV2 = dV2 - ((x(2) - -0.613169) / 0.000118031) * Ei(10);
 dV2 = dV2 - ((x(2) - -0.037037) / 0.00956049) * Ei(11);
 dV2 = dV2 - ((x(2) - -0.432099) / 0.00424911) * Ei(12);
 dV2 = dV2 - ((x(2) - -0.432099) / 0.00424911) * Ei(13);
 dV2 = dV2 - ((x(2) - -0.037037) / 0.00956049) * Ei(14);
 dV2 = dV2 - ((x(2) - -0.432099) / 0.00424911) * Ei(15);
 dV2 = dV2 - ((x(2) - 0.555556) / 0.0860445) * Ei(16);
 dV2 = dV2 - ((x(2) - -0.037037) / 0.00956049) * Ei(17);
 dV2 = dV2 - ((x(2) - 0.555556) / 0.0860445) * Ei(18);
 dV2 = dV2 - ((x(2) - -0.666667) / 0.0484) * Ei(19);
 dV2 = dV2 - ((x(2) - -0.259259) / 0.00239012) * Ei(20);
 dV2 = dV2 - ((x(2) - -0.888889) / 0.00537778) * Ei(21);
 dV2 = dV2 - ((x(2) - -0.888889) / 0.00537778) * Ei(22);
 dV2 = dV2 - ((x(2) - -0.259259) / 0.00239012) * Ei(23);
 dV2 = dV2 - ((x(2) - -0.753086) / 0.00026557) * Ei(24);
 dV2 = dV2 - ((x(2) - -0.888889) / 0.00537778) * Ei(25);
 dV2 = dV2 - ((x(2) - -0.679012) / 0.00106228) * Ei(26);
 dV2 = dV2 - ((x(2) - 0) / 0.4356) * Ei(27);
 dV2 = dV2 - ((x(2) - -0.679012) / 0.00106228) * Ei(28);
 dV2 = dV2 - ((x(2) - -0.585734) / 5.24581e-05) * Ei(29);
 dV2 = dV2 - ((x(2) - -0.552812) / 0.000209833) * Ei(30);

    %%-----------------------------------------------
    %% Computation of the derivative of the numerator
    %%-----------------------------------------------

 dU2 = 0.0;
 dU2 = dU2 + (0.991035 - ((x(2) - 0) / 0.4356) * Li(1)) * Ei(1);
 dU2 = dU2 + (0.967875 - ((x(2) - -0.666667) / 0.0484) * Li(2)) * Ei(2);
 dU2 = dU2 + (1.00248 - ((x(2) - -0.111111) / 0.0215111) * Li(3)) * Ei(3);
 dU2 = dU2 + (0.993099 - ((x(2) - 0.555556) / 0.0860445) * Li(4)) * Ei(4);
 dU2 = dU2 + (1.03233 - ((x(2) - 0) / 0.4356) * Li(5)) * Ei(5);
 dU2 = dU2 + (0.953868 - ((x(2) - -0.888889) / 0.00537778) * Li(6)) * Ei(6);
 dU2 = dU2 + (0.974689 - ((x(2) - -0.703704) / 0.00239012) * Li(7)) * Ei(7);
 dU2 = dU2 + (0.972573 - ((x(2) - -0.259259) / 0.00239012) * Li(8)) * Ei(8);
 dU2 = dU2 + (0.980401 - ((x(2) - 0.555556) / 0.0860445) * Li(9)) * Ei(9);
 dU2 = dU2 + (1.09248 - ((x(2) - -0.613169) / 0.000118031) * Li(10)) * Ei(10);
 dU2 = dU2 + (0.990434 - ((x(2) - -0.037037) / 0.00956049) * Li(11)) * Ei(11);
 dU2 = dU2 + (1.0109 - ((x(2) - -0.432099) / 0.00424911) * Li(12)) * Ei(12);
 dU2 = dU2 + (0.993938 - ((x(2) - -0.432099) / 0.00424911) * Li(13)) * Ei(13);
 dU2 = dU2 + (0.999673 - ((x(2) - -0.037037) / 0.00956049) * Li(14)) * Ei(14);
 dU2 = dU2 + (1.01348 - ((x(2) - -0.432099) / 0.00424911) * Li(15)) * Ei(15);
 dU2 = dU2 + (1.00644 - ((x(2) - 0.555556) / 0.0860445) * Li(16)) * Ei(16);
 dU2 = dU2 + (0.99922 - ((x(2) - -0.037037) / 0.00956049) * Li(17)) * Ei(17);
 dU2 = dU2 + (1.02342 - ((x(2) - 0.555556) / 0.0860445) * Li(18)) * Ei(18);
 dU2 = dU2 + (1.01202 - ((x(2) - -0.666667) / 0.0484) * Li(19)) * Ei(19);
 dU2 = dU2 + (0.996558 - ((x(2) - -0.259259) / 0.00239012) * Li(20)) * Ei(20);
 dU2 = dU2 + (0.979475 - ((x(2) - -0.888889) / 0.00537778) * Li(21)) * Ei(21);
 dU2 = dU2 + (0.969202 - ((x(2) - -0.888889) / 0.00537778) * Li(22)) * Ei(22);
 dU2 = dU2 + (1.0064 - ((x(2) - -0.259259) / 0.00239012) * Li(23)) * Ei(23);
 dU2 = dU2 + (0.972939 - ((x(2) - -0.753086) / 0.00026557) * Li(24)) * Ei(24);
 dU2 = dU2 + (1.03568 - ((x(2) - -0.888889) / 0.00537778) * Li(25)) * Ei(25);
 dU2 = dU2 + (0.988424 - ((x(2) - -0.679012) / 0.00106228) * Li(26)) * Ei(26);
 dU2 = dU2 + (1.00141 - ((x(2) - 0) / 0.4356) * Li(27)) * Ei(27);
 dU2 = dU2 + (1.02471 - ((x(2) - -0.679012) / 0.00106228) * Li(28)) * Ei(28);
 dU2 = dU2 + (0.996492 - ((x(2) - -0.585734) / 5.24581e-05) * Li(29)) * Ei(29);
 dU2 = dU2 + (0.942638 - ((x(2) - -0.552812) / 0.000209833) * Li(30)) * Ei(30);

    %%--------------------------------------------------------
    %% Computation of the second derivative of the denominator
    %%--------------------------------------------------------

 ddV = 0.0;
 ddV = ddV + (((x(2) - 0) / 0.4356) * ((x(1) - 55.5556) / 1344.44)) * Ei(1);
 ddV = ddV + (((x(2) - -0.666667) / 0.0484) * ((x(1) - 370.37) / 597.531)) * Ei(2);
 ddV = ddV + (((x(2) - -0.111111) / 0.0215111) * ((x(1) - 444.444) / 5377.78)) * Ei(3);
 ddV = ddV + (((x(2) - 0.555556) / 0.0860445) * ((x(1) - 444.444) / 5377.78)) * Ei(4);
 ddV = ddV + (((x(2) - 0) / 0.4356) * ((x(1) - 148.148) / 597.531)) * Ei(5);
 ddV = ddV + (((x(2) - -0.888889) / 0.00537778) * ((x(1) - 629.63) / 2390.12)) * Ei(6);
 ddV = ddV + (((x(2) - -0.703704) / 0.00239012) * ((x(1) - 629.63) / 2390.12)) * Ei(7);
 ddV = ddV + (((x(2) - -0.259259) / 0.00239012) * ((x(1) - 629.63) / 2390.12)) * Ei(8);
 ddV = ddV + (((x(2) - 0.555556) / 0.0860445) * ((x(1) - 580.247) / 265.569)) * Ei(9);
 ddV = ddV + (((x(2) - -0.613169) / 0.000118031) * ((x(1) - 777.778) / 21511.1)) * Ei(10);
 ddV = ddV + (((x(2) - -0.037037) / 0.00956049) * ((x(1) - 629.63) / 2390.12)) * Ei(11);
 ddV = ddV + (((x(2) - -0.432099) / 0.00424911) * ((x(1) - 629.63) / 2390.12)) * Ei(12);
 ddV = ddV + (((x(2) - -0.432099) / 0.00424911) * ((x(1) - 753.086) / 1062.28)) * Ei(13);
 ddV = ddV + (((x(2) - -0.037037) / 0.00956049) * ((x(1) - 753.086) / 1062.28)) * Ei(14);
 ddV = ddV + (((x(2) - -0.432099) / 0.00424911) * ((x(1) - 901.235) / 4249.11)) * Ei(15);
 ddV = ddV + (((x(2) - 0.555556) / 0.0860445) * ((x(1) - 851.852) / 9560.5)) * Ei(16);
 ddV = ddV + (((x(2) - -0.037037) / 0.00956049) * ((x(1) - 901.235) / 4249.11)) * Ei(17);
 ddV = ddV + (((x(2) - 0.555556) / 0.0860445) * ((x(1) - 654.321) / 1062.28)) * Ei(18);
 ddV = ddV + (((x(2) - -0.666667) / 0.0484) * ((x(1) - 481.481) / 2390.12)) * Ei(19);
 ddV = ddV + (((x(2) - -0.259259) / 0.00239012) * ((x(1) - 753.086) / 1062.28)) * Ei(20);
 ddV = ddV + (((x(2) - -0.888889) / 0.00537778) * ((x(1) - 753.086) / 1062.28)) * Ei(21);
 ddV = ddV + (((x(2) - -0.888889) / 0.00537778) * ((x(1) - 835.391) / 472.123)) * Ei(22);
 ddV = ddV + (((x(2) - -0.259259) / 0.00239012) * ((x(1) - 901.235) / 4249.11)) * Ei(23);
 ddV = ddV + (((x(2) - -0.753086) / 0.00026557) * ((x(1) - 851.852) / 9560.5)) * Ei(24);
 ddV = ddV + (((x(2) - -0.888889) / 0.00537778) * ((x(1) - 934.156) / 1888.49)) * Ei(25);
 ddV = ddV + (((x(2) - -0.679012) / 0.00106228) * ((x(1) - 753.086) / 1062.28)) * Ei(26);
 ddV = ddV + (((x(2) - 0) / 0.4356) * ((x(1) - 259.259) / 2390.12)) * Ei(27);
 ddV = ddV + (((x(2) - -0.679012) / 0.00106228) * ((x(1) - 901.235) / 4249.11)) * Ei(28);
 ddV = ddV + (((x(2) - -0.585734) / 5.24581e-05) * ((x(1) - 777.778) / 21511.1)) * Ei(29);
 ddV = ddV + (((x(2) - -0.552812) / 0.000209833) * ((x(1) - 777.778) / 21511.1)) * Ei(30);

    %%------------------------------------------------------
    %% Computation of the second derivative of the numerator
    %%------------------------------------------------------

 ddU = 0.0;
 ddU = ddU - (0.991035 * ((x(1) - 55.5556) / 1344.44) + (0 - ((x(1) - 55.5556) / 1344.44) * Li(1)) * ((x(2) - 0) / 0.4356)) * Ei(1);
 ddU = ddU - (0.967875 * ((x(1) - 370.37) / 597.531) + (0 - ((x(1) - 370.37) / 597.531) * Li(2)) * ((x(2) - -0.666667) / 0.0484)) * Ei(2);
 ddU = ddU - (1.00248 * ((x(1) - 444.444) / 5377.78) + (0 - ((x(1) - 444.444) / 5377.78) * Li(3)) * ((x(2) - -0.111111) / 0.0215111)) * Ei(3);
 ddU = ddU - (0.993099 * ((x(1) - 444.444) / 5377.78) + (0 - ((x(1) - 444.444) / 5377.78) * Li(4)) * ((x(2) - 0.555556) / 0.0860445)) * Ei(4);
 ddU = ddU - (1.03233 * ((x(1) - 148.148) / 597.531) + (0 - ((x(1) - 148.148) / 597.531) * Li(5)) * ((x(2) - 0) / 0.4356)) * Ei(5);
 ddU = ddU - (0.953868 * ((x(1) - 629.63) / 2390.12) + (0 - ((x(1) - 629.63) / 2390.12) * Li(6)) * ((x(2) - -0.888889) / 0.00537778)) * Ei(6);
 ddU = ddU - (0.974689 * ((x(1) - 629.63) / 2390.12) + (0 - ((x(1) - 629.63) / 2390.12) * Li(7)) * ((x(2) - -0.703704) / 0.00239012)) * Ei(7);
 ddU = ddU - (0.972573 * ((x(1) - 629.63) / 2390.12) + (0 - ((x(1) - 629.63) / 2390.12) * Li(8)) * ((x(2) - -0.259259) / 0.00239012)) * Ei(8);
 ddU = ddU - (0.980401 * ((x(1) - 580.247) / 265.569) + (0 - ((x(1) - 580.247) / 265.569) * Li(9)) * ((x(2) - 0.555556) / 0.0860445)) * Ei(9);
 ddU = ddU - (1.09248 * ((x(1) - 777.778) / 21511.1) + (0 - ((x(1) - 777.778) / 21511.1) * Li(10)) * ((x(2) - -0.613169) / 0.000118031)) * Ei(10);
 ddU = ddU - (0.990434 * ((x(1) - 629.63) / 2390.12) + (0 - ((x(1) - 629.63) / 2390.12) * Li(11)) * ((x(2) - -0.037037) / 0.00956049)) * Ei(11);
 ddU = ddU - (1.0109 * ((x(1) - 629.63) / 2390.12) + (0 - ((x(1) - 629.63) / 2390.12) * Li(12)) * ((x(2) - -0.432099) / 0.00424911)) * Ei(12);
 ddU = ddU - (0.993938 * ((x(1) - 753.086) / 1062.28) + (0 - ((x(1) - 753.086) / 1062.28) * Li(13)) * ((x(2) - -0.432099) / 0.00424911)) * Ei(13);
 ddU = ddU - (0.999673 * ((x(1) - 753.086) / 1062.28) + (0 - ((x(1) - 753.086) / 1062.28) * Li(14)) * ((x(2) - -0.037037) / 0.00956049)) * Ei(14);
 ddU = ddU - (1.01348 * ((x(1) - 901.235) / 4249.11) + (0 - ((x(1) - 901.235) / 4249.11) * Li(15)) * ((x(2) - -0.432099) / 0.00424911)) * Ei(15);
 ddU = ddU - (1.00644 * ((x(1) - 851.852) / 9560.5) + (0 - ((x(1) - 851.852) / 9560.5) * Li(16)) * ((x(2) - 0.555556) / 0.0860445)) * Ei(16);
 ddU = ddU - (0.99922 * ((x(1) - 901.235) / 4249.11) + (0 - ((x(1) - 901.235) / 4249.11) * Li(17)) * ((x(2) - -0.037037) / 0.00956049)) * Ei(17);
 ddU = ddU - (1.02342 * ((x(1) - 654.321) / 1062.28) + (0 - ((x(1) - 654.321) / 1062.28) * Li(18)) * ((x(2) - 0.555556) / 0.0860445)) * Ei(18);
 ddU = ddU - (1.01202 * ((x(1) - 481.481) / 2390.12) + (0 - ((x(1) - 481.481) / 2390.12) * Li(19)) * ((x(2) - -0.666667) / 0.0484)) * Ei(19);
 ddU = ddU - (0.996558 * ((x(1) - 753.086) / 1062.28) + (0 - ((x(1) - 753.086) / 1062.28) * Li(20)) * ((x(2) - -0.259259) / 0.00239012)) * Ei(20);
 ddU = ddU - (0.979475 * ((x(1) - 753.086) / 1062.28) + (0 - ((x(1) - 753.086) / 1062.28) * Li(21)) * ((x(2) - -0.888889) / 0.00537778)) * Ei(21);
 ddU = ddU - (0.969202 * ((x(1) - 835.391) / 472.123) + (0 - ((x(1) - 835.391) / 472.123) * Li(22)) * ((x(2) - -0.888889) / 0.00537778)) * Ei(22);
 ddU = ddU - (1.0064 * ((x(1) - 901.235) / 4249.11) + (0 - ((x(1) - 901.235) / 4249.11) * Li(23)) * ((x(2) - -0.259259) / 0.00239012)) * Ei(23);
 ddU = ddU - (0.972939 * ((x(1) - 851.852) / 9560.5) + (0 - ((x(1) - 851.852) / 9560.5) * Li(24)) * ((x(2) - -0.753086) / 0.00026557)) * Ei(24);
 ddU = ddU - (1.03568 * ((x(1) - 934.156) / 1888.49) + (0 - ((x(1) - 934.156) / 1888.49) * Li(25)) * ((x(2) - -0.888889) / 0.00537778)) * Ei(25);
 ddU = ddU - (0.988424 * ((x(1) - 753.086) / 1062.28) + (0 - ((x(1) - 753.086) / 1062.28) * Li(26)) * ((x(2) - -0.679012) / 0.00106228)) * Ei(26);
 ddU = ddU - (1.00141 * ((x(1) - 259.259) / 2390.12) + (0 - ((x(1) - 259.259) / 2390.12) * Li(27)) * ((x(2) - 0) / 0.4356)) * Ei(27);
 ddU = ddU - (1.02471 * ((x(1) - 901.235) / 4249.11) + (0 - ((x(1) - 901.235) / 4249.11) * Li(28)) * ((x(2) - -0.679012) / 0.00106228)) * Ei(28);
 ddU = ddU - (0.996492 * ((x(1) - 777.778) / 21511.1) + (0 - ((x(1) - 777.778) / 21511.1) * Li(29)) * ((x(2) - -0.585734) / 5.24581e-05)) * Ei(29);
 ddU = ddU - (0.942638 * ((x(1) - 777.778) / 21511.1) + (0 - ((x(1) - 777.778) / 21511.1) * Li(30)) * ((x(2) - -0.552812) / 0.000209833)) * Ei(30);
  clear Ei;
  clear Li;
    res =  (((ddU*V + dU1*dV2 - ddV*U - dV1*dU2)*V*V - (dU1*V-dV1*U)*2*V*dV2)/(V*V*V*V));
