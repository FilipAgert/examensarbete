clear
close all
Eexc = [13, 14 , 15, 17, 20, 22, 25, 30 ,35 ,40 ,45, 50, 55, 60, 65];
fusion_prob_martin = [0.293,0.485,0.578,0.61,0.538,0.46,0.352,0.238,0.21,0.221,0.249,0.272,0.296,0.322,0.343];
fusion_prob_markovMM33 = [0.6384451564, 0.8358019905, 0.7476546942, 0.6113791779, 0.5644451736, 0.4540475598, 0.378766928, 0.2670794886, 0.203474034, 0.2431601931, 0.2481580806, 0.2567266904, 0.332821196, 0.338954117, 0.3443735955];
fusion_prob_markovMM40 = 1/100*[0.28,86.06,66.24,48.12,39.93,28.84,24.50,21.58,19.35,24.85,25.74,26.66,34.18,34.81,35.38];

plot(Eexc, fusion_prob_martin, 'b-o')
hold on
plot(Eexc, fusion_prob_markovMM33, 'r-o')
plot(Eexc, fusion_prob_markovMM40, 'k-o')
legend('Random walk', 'Markov chain M\_Max = 33', 'Markov chain M\_Max = 40')
title('Result comparison of random walk- and markov chain approach. Z = 102, A = 256')
xlabel('Excitation energy [MeV]')
ylabel('Chance of fusion')

xlim([0, max(Eexc)*1.1])
ylim([0,1])