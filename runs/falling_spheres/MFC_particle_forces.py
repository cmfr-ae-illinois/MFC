import numpy as np 

def Loth_CD(Re, M, gamma):
    if (Re >= 45.0): # compression dominated regime
        if (M < 1.5):
            CM = 1.65 + 0.65 * np.tanh(4.0 * M - 3.4)
        else:
            CM = 2.18 - 0.13 * np.tanh(0.9 * M - 2.7)

        if (M < 0.8):
            GM = 166.0 * M**3.0 + 3.29 * M**2.0 - 10.9 * M + 20.0
        else:
            GM = 5.0 + 40.0 * M**(-3.0) 
            
        if (M < 1.0):
            HM = 0.0239 * M**3.0 + 0.212 * M**2.0 - 0.074 * M + 1.0
        else:
            HM = 0.93 + 1.0 / (3.5 + M**5.0)

        CD = (24.0 / Re) * (1.0 + 0.15 * Re**0.687) * HM + (0.42 * CM) / (1.0 + (42500.0 / (Re**(1.16 * CM))) + (GM / (Re**0.5)))

    else: # rarefaction dominated regime
        Kn = np.sqrt(np.pi * gamma / 2) * M / Re

        fKn = 1.0 / (1.0 + Kn * (2.514 + 0.8 * np.exp(-0.55 / Kn)))

        CDKnRe = 24.0 / Re * (1.0 + 0.15 * Re**0.687) * fKn

        s = M * np.sqrt(gamma / 2.0)

        CDfm = (1.0 + 2.0 * s**2.0) * np.exp(-s**2.0) / (s**3.0 * np.sqrt(np.pi)) + (4.0 * s**4.0 + 4.0 * s**2.0 - 1.0) * erf(s) / (2.0 * s**4.0) + 2.0 / (3.0 * s) * np.sqrt(np.pi)

        if (M <= 1.0):
            JM = 2.26 - 0.1 / M + 0.14 / (M**3.0)
        else:
            JM = 1.6 + 0.25 / M + 0.11 / (M**2.0) + 0.44 / (M**3.0)

        CDfmRe = CDfm / (1.0 + (CDfm / JM - 1.0) * np.sqrt(Re / 45.0))

        CD = CDKnRe / (1.0 + M**4.0) + M**4.0 * CDfmRe / (1.0 + M**4.0)

    return CD

def Osnes_CD(alpha, Re, M, gamma):
    CDLoth = Loth_CD(Re, M, gamma)

    b1 = 5.81  * alpha / ((1.0 - alpha)**2.0) + 0.48 * alpha**(1.0/3.0) / ((1.0 - alpha)**3)

    b2 = (1.0 - alpha)**2.0 * alpha**3.0 * Re * (0.95 + 0.61 * alpha**3.0 / ((1.0 - alpha)**2.0))

    b3 = np.min((np.sqrt(20.0 * M), 1.0)) * (5.65 * alpha - 22.0 * alpha**2.0 + 23.4*alpha**3.0) * (1.0 + np.tanh((M - (0.65 - 0.24 * alpha)) / 0.35))

    CD = CDLoth / (1.0 - alpha) + b3 + (24.0 * (1.0 - alpha) / Re) * (b1 + b2)

    return CD

def Singh_CD(Re, M, gamma):
    n = 0.74 
    a0 = 0.3555
    delta0 = 9.4 
    C0 = 24.0 / (delta0**2)

    if (M <= 1.0):
        Ms = M
        U_ratio = 1.0
        T_ratio = 1.0 
        a = 1.0 
        C1 = 1.0
    else:
        U_ratio = (2.0 + (gamma - 1.0) * M**2.0) / ((gamma + 1.0) * M**2.0)
        T_ratio = ((gamma - 1.0) * M**2.0 + 2.0) * (2.0 * gamma * M**2.0 - (gamma - 1.0)) / ((gamma + 1.0)**2.0 * (M**2.0))
        Ms = np.sqrt(((gamma - 1.0) * M**2.0 + 2.0) / (2.0 * gamma * M**2.0 - (gamma - 1.0)))
        a = 1.0 / (a0 * (M - 1.0) + 1.0)
        ainf = 1.0 / (a0 * M)
        C1 = (0.9 - C0 * (1.0 + (gamma - 1.0)**2.0 / (4.0 * gamma))**(gamma / (gamma - 1.0))) / (1.0 - ainf * (gamma - 1.0) / (gamma + 1.0))

    Theta = (1.0 + (gamma - 1.0) * Ms**2.0 / 2.0)**(gamma / (gamma - 1.0))
    Rehat = Re * (1.0 / (a**2.0 * T_ratio))**n * Theta**((gamma + 1.0) / (2.0 * gamma) - (gamma - 1.0) * n / gamma)
    CD = (C0 * Theta) * (1.0 + delta0 / np.sqrt(Rehat))**2.0 + C1 * (1.0 - a * U_ratio)

    return CD

def Khalloufi_CD(Re, M, phi, gamma):
    Me = M * (1.0 + 0.33 * (1.0 - np.exp(-34.0 * phi)))

    CDiso = Singh_CD(Re, Me, gamma)

    CD = CDiso * (1.0 + 2.0 * phi) / ((1.0 - phi)**2)

    return CD
