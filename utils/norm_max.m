function tc = norm_max(tc)

tc = tc ./abs(max(tc(:)));

end