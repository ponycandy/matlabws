function [D_X,D_J] = s_curv_model(X,U,D_U,s,Track)
e = X(1);
erro_PHI = X(2);

D_PHI_ds = U(1);

D2_PHI_ds2 = D_U(1);

curv = interp1(Track.S,Track.curv,s);

D_e = (1 - curv*e)*tan(erro_PHI);
D_erro_PHI = -curv + D_PHI_ds;

D_X = [D_e;D_erro_PHI];

D_J = 1*(D_PHI_ds)^2 + 1*e^2 + 1*(erro_PHI)^2 - 0*(D2_PHI_ds2)^2;

end