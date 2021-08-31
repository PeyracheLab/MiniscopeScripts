neuron.options.gSig = 3;
neuron.options.gSiz = 12;
neuron.options.ring_radius = 30;

[cn, pnr] = neuron.correlation_pnr_parallel([1, 5000]);


neuron.options.min_corr = 0.8;
neuron.options.min_pnr = 8;


%find all local maximum as initialization point
tmp_d = max(1,round(gSiz/4));
v_max = ordfilt2(cn.*pnr, tmp_d^2, true(tmp_d));
ind = (v_max==cn.*pnr);

figure('papersize', [12, 3]);
init_fig;
subplot(131);
imagesc(cn, [0, 1]);
title('local corr. image');
axis equal off tight;
subplot(132);
pnr_vmax = max(pnr(:))*0.8;
imagesc(pnr)%, [3, pnr_vmax]);
axis equal off tight;
title('PNR image');

subplot(133);
imagesc(cn.*pnr)%, [3, pnr_vmax]);
hold on;
tmp_ind = ind & (cn>=neuron.options.min_corr) & (pnr>=neuron.options.min_pnr);
[r, c] = find(tmp_ind);
ax_seeds = plot(c, r, '.m', 'markersize', 10);
axis equal off tight;
title('candidate seed pixels');
ylabel('PNR');
xlabel('Cn');

