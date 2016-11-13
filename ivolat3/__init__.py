import bs


pc_dict = {0: 0, 'c': 0, 'C': 0,
        1: 1, 'p': 1, 'P': 1,
        2: 2, 'f': 2, 'F': 2, 'm': 2, 'M': 2}


def porc(x):
    if x not in pc_dict: return None
    return pc_dict[x]


def prem_put(s, k, r, q, t, sigma):
    return bs.bs_prem_put(s, k, r, q, t, sigma)


def prem_call(s, k, r, q, t, sigma):
    return bs.bs_prem_call(s, k, r, q, t, sigma)


def delta_put(s, k, r, q, t, sigma):
    return bs.bs_delta_put(s, k, r, q, t, sigma)


def delta_call(s, k, r, q, t, sigma):
    return bs.bs_delta_call(s, k, r, q, t, sigma)


def vega(s, k, r, q, t, sigma, pc=None):
    return bs.bs_vega(s, k, r, q, t, sigma)


def theta_put(s, k, r, q, t, sigma):
    return bs.bs_theta_put(s, k, r, q, t, sigma)


def theta_call(s, k, r, q, t, sigma):
    return bs.bs_theta_call(s, k, r, q, t, sigma)


def rho_put(s, k, r, q, t, sigma):
    return bs.bs_rho_put(s, k, r, q, t, sigma)


def rho_call(s, k, r, q, t, sigma):
    return bs.bs_rho_call(s, k, r, q, t, sigma)


def gamma(s, k, r, q, t, sigma, pc=None):
    return bs.bs_gamma(s, k, r, q, t, sigma)


def vanna(s, k, r, q, t, sigma, pc=None):
    return bs.bs_vanna(s, k, r, q, t, sigma)


def charm_put(s, k, r, q, t, sigma):
    return bs.bs_charm_put(s, k, r, q, t, sigma)


def charm_call(s, k, r, q, t, sigma):
    return bs.bs_charm_call(s, k, r, q, t, sigma)


def speed(s, k, r, q, t, sigma, pc=None):
    return bs.bs_speed(s, k, r, q, t, sigma)


def zomma(s, k, r, q, t, sigma, pc=None):
    return bs.bs_zomma(s, k, r, q, t, sigma)


def color(s, k, r, q, t, sigma, pc=None):
    return bs.bs_color(s, k, r, q, t, sigma)


def DvegaDtime(s, k, r, q, t, sigma, pc=None):
    return bs.bs_DvegaDtime(s, k, r, q, t, sigma)


def vomma(s, k, r, q, t, sigma, pc=None):
    return bs.bs_vomma(s, k, r, q, t, sigma)


def ultima(s, k, r, q, t, sigma, pc=None):
    return bs.bs_ultima(s, k, r, q, t, sigma)


def dualdelta_put(s, k, r, q, t, sigma):
    return bs.bs_dualdelta_put(s, k, r, q, t, sigma)


def dualdelta_call(s, k, r, q, t, sigma):
    return bs.bs_dualdelta_call(s, k, r, q, t, sigma)


def dualgamma(s, k, r, q, t, sigma, pc=None):
    return bs.bs_dualgamma(s, k, r, q, t, sigma)


def ivolat_put(s, k, r, q, t, p):
    return bs.bs_ivolat_put(s, k, r, q, t, p)


def ivolat_call(s, k, r, q, t, p):
    return bs.bs_ivolat_call(s, k, r, q, t, p)


prem_dict = {0: prem_call, 1: prem_put}
delta_dict = {0: delta_call, 1: delta_put}
theta_dict = {0: theta_call, 1: theta_put}
rho_dict = {0: rho_call, 1: rho_put}
charm_dict = {0: charm_call, 1: charm_put}
dualdelta_dict = {0: dualdelta_call, 1: dualdelta_put}
ivolat_dict = {0: ivolat_call, 1: ivolat_put}


def prem(s, k, r, q, t, sigma, pc):
    pc = porc(pc)
    if pc == 2: return
    return prem_dict[pc](s, k, r, q, t, sigma)


def delta(s, k, r, q, t, sigma, pc):
    pc = porc(pc)
    if pc == 2: return 1
    return delta_dict[pc](s, k, r, q, t, sigma)


def theta(s, k, r, q, t, sigma, pc):
    pc = porc(pc)
    if pc == 2: return 0
    return theta_dict[pc](s, k, r, q, t, sigma)


def rho(s, k, r, q, t, sigma, pc):
    pc = porc(pc)
    if pc == 2: return 0
    return rho_dict[pc](s, k, r, q, t, sigma)


def charm(s, k, r, q, t, sigma, pc):
    pc = porc(pc)
    if pc == 2: return 0
    return charm_dict[pc](s, k, r, q, t, sigma)


def dualdelta(s, k, r, q, t, sigma, pc):
    pc = porc(pc)
    if pc == 2: return 0
    return dualdelta_dict[pc](s, k, r, q, t, sigma)


def ivolat(s, k, r, q, t, p, pc):
    pc = porc(pc)
    if pc == 2: return
    #return ivolat_dict[pc](s, k, r, q, t, p)
    return bs.bs_ivolat(s, k, r, q, t, p, pc)


if __name__ == "__main__":
    s = 15420
    k = 15500
    r = 0.001
    q = 0.0
    t = 10.5/365
    p = 100
    sigma = 0.15

    print(prem_put(s, k, r, q, t, sigma))
    print(prem_call(s, k, r, q, t, sigma))
    print(delta_put(s, k, r, q, t, sigma))
    print(delta_call(s, k, r, q, t, sigma))
    print(vega(s, k, r, q, t, sigma))
    print(theta_put(s, k, r, q, t, sigma))
    print(theta_call(s, k, r, q, t, sigma))
    print(rho_put(s, k, r, q, t, sigma))
    print(rho_call(s, k, r, q, t, sigma))
    print(gamma(s, k, r, q, t, sigma))
    print(vanna(s, k, r, q, t, sigma))
    print(charm_put(s, k, r, q, t, sigma))
    print(charm_call(s, k, r, q, t, sigma))
    print(speed(s, k, r, q, t, sigma))
    print(zomma(s, k, r, q, t, sigma))
    print(color(s, k, r, q, t, sigma))
    print(DvegaDtime(s, k, r, q, t, sigma))
    print(vomma(s, k, r, q, t, sigma))
    print(ultima(s, k, r, q, t, sigma))
    print(dualdelta_put(s, k, r, q, t, sigma))
    print(dualdelta_call(s, k, r, q, t, sigma))
    print(dualgamma(s, k, r, q, t, sigma))
    print(ivolat_put(s, k, r, q, t, p))
    print(ivolat_call(s, k, r, q, t, p))
