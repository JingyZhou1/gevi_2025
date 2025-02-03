function r2 = compute_r2(data, prd)

r2 =  1 - (sum((data - prd).^2))./(sum((prd - mean(prd)).^2));

end
