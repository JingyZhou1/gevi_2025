function r2 = fit_gevi2gcamp(gevi, prm, transbasis, gcamp, rng)

    prd = gevi2gcamp_fun(gevi, prm, transbasis, gcamp, rng);


    r2 = sum((prd - gcamp).^2);

end