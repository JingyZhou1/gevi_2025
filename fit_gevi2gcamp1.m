function r2 = fit_gevi2gcamp1(gevi, prm, power, transbasis, gcamp, rng)

    prd = gevi2gcamp_fun1(gevi, prm, power, transbasis, gcamp, rng);

    r2 = sum((prd - gcamp).^2);

end