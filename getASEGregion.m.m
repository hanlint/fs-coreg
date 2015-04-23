function region = getASEGregion(ras, aseg)

ras = shiftdim(ras);
v = round(inv(aseg.tkrvox2ras)*ras);

region = aseg.vol(v(2),v(1),v(3));

end