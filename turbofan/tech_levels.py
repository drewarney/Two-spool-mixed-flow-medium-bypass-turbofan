"""Component technology level 4 (Mattingly, *Elements of Propulsion*,
Table 6.2).  Values transcribed for the component types used in design
problems 7.D1 / 7.D2 - check them against your edition of the table.

    diffuser  pi_d,max   type A (subsonic, nacelle)          0.995
                         type C (supersonic, embedded)       0.96
    fan       e_f                                            0.89
    compressor e_c                                           0.90
    burner    pi_b                                           0.96
              eta_b                                          0.999
    turbine   e_t        cooled                              0.89
    mixer     pi_M,max   (given in problem 7.D2)             0.96
    nozzle    pi_n       type D (fixed convergent)           0.995
                         type F (variable C-D, supersonic)   0.97
"""
LEVEL4 = {
    "pi_d_max": {"A": 0.995, "B": 0.97, "C": 0.96},
    "e_f": 0.89,
    "e_c": 0.90,
    "pi_b": 0.96,
    "eta_b": 0.999,
    "e_t": {"uncooled": 0.91, "cooled": 0.89},
    "pi_n": {"D": 0.995, "E": 0.98, "F": 0.97},
    "Tt4_max_K": 2000.0,
}


def hp1_losses():
    """Problem 7.D1 (HP-1): cooled turbine, type A diffuser, type D nozzle."""
    L = LEVEL4
    return dict(pi_d_max=L["pi_d_max"]["A"], pi_b=L["pi_b"], pi_n=L["pi_n"]["D"], pi_fn=L["pi_n"]["D"],
                e_c=L["e_c"], e_f=L["e_f"], e_t=L["e_t"]["cooled"], eta_b=L["eta_b"], eta_m=0.99)


def hf1_losses():
    """Problem 7.D2 (HF-1): cooled turbine, type C diffuser, type F nozzle,
    M6 = 0.4, pi_M,max = 0.96."""
    L = LEVEL4
    return dict(pi_d_max=L["pi_d_max"]["C"], pi_b=L["pi_b"], pi_n=L["pi_n"]["F"],
                e_c=L["e_c"], e_f=L["e_f"], e_t=L["e_t"]["cooled"], eta_b=L["eta_b"], eta_m=0.99,
                M6=0.4, pi_M_max=0.96)
