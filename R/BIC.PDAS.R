BIC.PDAS <- function(X, y, method, sel="vote", al=0, tau, mu=dim(X)[1]/log(dim(X)[1]),
del=0, weight = 0.5, MaxIt=1, Lmax=1, Lmin=1e-3, N=100) {



    out = pdas_path(X, y, method, sel, al, tau, mu, del, weight, MaxIt, Lmax, Lmin, N);
    ithist_as = out$ithist_as;
    ithist_it = out$ithist_it;
    ithist_x = out$ithist_x;
    ithist_res = out$ithist_res;
    ithist_bic = out$ithist_bic;
    Lam = out$Lam;
    id = out$id;

    if  (sel == "vote"){
        #ii == which(ithist_as == mode(ithist_as));
        freq=table(ithist_as);
        value=as.numeric(names(freq)[which.max(freq)]);
        ii = which(ithist_as == value);
        id = ii[length(ii)];
        x =  ithist_x[,id];
        lam = Lam[id]
    }
    else if (sel == "bic"){
        ii = which(ithist_bic == min(ithist_bic));
        id = ii[length(ii)];
        x =  ithist_x[,id];
        lam = Lam[id]
    }
    obj = list(x = x, lam = lam, id = id, Lam = Lam,
    ithist_as = ithist_as, ithist_it = ithist_it, ithist_x = ithist_x,
    ithist_res = ithist_res, ithist_bic = ithist_bic)
    class(obj) = "BIC.PDAS"
    obj
}