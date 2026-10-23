# -*- coding: utf-8 -*-
"""
Created on Thu Sep 08 16:03:00 2026
@author: Amanda Gabriela

---- MODELE SECHOIR SOLAIRE - MODELE FICK ----

"""
import numpy as np
from scipy.optimize import least_squares
import pandas as pd
import json

# LECTURE DE DONNEES FICHIER .JSON
with open("donnees_diff.json","r") as f:
    donnees_diff = json.load(f)
    
donnees_diff["Aplaque"] = donnees_diff["Lar_plaque"] * donnees_diff["Lon_plaque"]
donnees_diff["Aent"]    = donnees_diff["dist_verre_ardoise"] * donnees_diff["Lar_plaque"]
donnees_diff["Lplaque"] = donnees_diff["Lon_plaque"]

# LECTURE DE DONNEES FICHIER EXCEL EXPERIMENTAL    
donnees_manip = pd.read_excel("manip_exp.xlsx", sheet_name= "Manip 2")
t_exp = donnees_manip["t"].values * 60.0

Tamb_exp = donnees_manip["Tamb"].values
Rayonnement_exp = donnees_manip["S"].values
Tentree_0 = donnees_manip["Tentree"].values[0]
Tplaque_0 = donnees_manip["Tplaque"].values[0]
Tsortie_0 = donnees_manip["Tsortie"].values[0]

def donnees_experimentales(t):
    Tamb_t = np.interp(t, t_exp, Tamb_exp)
    Ray_t = np.interp(t, t_exp, Rayonnement_exp)
    return Tamb_t, Ray_t

# SUBROUTINE D'INTERPOLATION DES PROPRIÉTÉS DE L'AIR
propes_air = pd.read_excel("proprietes_air.xlsx", sheet_name= "Air")
T_tab = propes_air["T"].values
k_tab = propes_air["k"].values
v_tab = propes_air["v"].values
Pr_tab = propes_air["Pr"].values

def proprietes_air(Tmoy):
    kair = np.interp(Tmoy, T_tab, k_tab)
    vair = np.interp(Tmoy, T_tab, v_tab)
    Pr = np.interp(Tmoy, T_tab, Pr_tab)
    return kair, vair, Pr

# SUBROUTINE CALCUL HMIX
def calc_Hmix(Tplaque,Tair_ref,g,angle_plaque,Lplaque,v):
    Tmoy = (Tair_ref + Tplaque)/2.0
    deltaT = Tplaque - Tair_ref
    beta = 1/(Tmoy + 273.15)
    kair, vair, Pr = proprietes_air(Tmoy)
    Gr = (beta * g * np.cos(np.radians(90.0 - angle_plaque)) * deltaT * Lplaque**3) / vair**2
    Ra = Gr * Pr

    if Ra < 1.0e9:
      Nu_nat = 0.678 * (max(0.0, Ra) / (0.952 + Pr))**0.25
    else:
      Nu_nat = 0.021 * max(0.0, Ra)**0.4
    Rey = v * Lplaque / vair
    if Rey < 5.0e5:
       Nu_forc = 0.664 * np.sqrt(max(0.0, Rey)) * Pr**(1.0 / 3.0)
    else:
       Nu_forc = 0.037 * max(0.0, Rey)**0.8 * Pr**(1.0 / 3.0)
    Nu_mix = (Nu_forc**3 + Nu_nat**3)**(1.0 / 3.0)
    h_mix = Nu_mix * kair / Lplaque
    return h_mix

# SUBROUTINE CALCUL Hlv
def calc_Hlv(Tentree, Tsortie, Cpliq_eau, Cpvap_eau, DHvap, Tebu):
    Tmoy = (Tentree + Tsortie) / 2.0
    Hlv = DHvap + (Cpvap_eau - Cpliq_eau) * (Tmoy - Tebu)
    return Hlv

# FUNCTION CALCUL P_SAT
def calc_psat(T_degC):
    T_K = T_degC + 273.15
    Psat = (100.0 * np.exp(75.051 - (7293.3 / T_K) - 8.593 * np.log(T_K) + 0.00617 * T_K))
    return Psat

D0_banane = 9.151e-5
Ea_banane = 34290.0
Rgas_banane = 8.314

def calc_D_banane(T_degC):
    T_k = T_degC + 273.15
    return D0_banane * np.exp(-Ea_banane/(Rgas_banane*(T_k)))

def calc_dMRdt(Theta, D_local, L_demi, n_termos=6):
    MR = 0.0
    somme_deriv = 0.0

    for n in range(n_termos):
        m = 2*n+1
        expo = np.exp(-(m**2)*(np.pi**2)*Theta/(4*L_demi**2))
        MR += (8.0/(np.pi**2)) * (1.0/m**2) * expo
        somme_deriv += expo
    dMRdt = -2.0*D_local/(L_demi**2) * somme_deriv
    return MR, dMRdt

# SUBROUTINE CALCUL Hair
def calcul_Hair(Tamb, HR, Patm, Fm, Aproduit, mdot_air):
    Pvap_amb = HR * calc_psat(Tamb)
    Hair_amb = 0.622 * Pvap_amb / (Patm - Pvap_amb)
    Hair_sortie = Hair_amb + (Fm * Aproduit) / mdot_air
    return Hair_amb, Hair_sortie

# SUBROUTINE RESOLUTION SYSTEME ALGEBRIQUE
def equations(x, Tplaque, Tamb, Fm_interne, Aplaque, Aent, Aproduit, Cpair, g, H_sous, H_sur, Ccheminee, Patm, Rair,
              angle_plaque, Lplaque, Cpliq_eau, Cpvap_eau, DHvap, Tebu, UA_chambre, HR):

    Q, Tentree, Tsortie, debit_mair, v = x

    mdot_evap = Fm_interne * Aproduit

    Tair_ref = (Tamb + Tentree) / 2.0
    h_mix = calc_Hmix(Tplaque, Tair_ref, g, angle_plaque, Lplaque, v)
    rho_air = Patm / (Rair * (Tair_ref + 273.15))
    Hlv = calc_Hlv(Tentree, Tsortie, Cpliq_eau, Cpvap_eau, DHvap, Tebu)

    Eq1 = Q - debit_mair * Cpair * (Tentree - Tamb)
    Eq2 = Q - h_mix * Aplaque * (Tplaque - Tair_ref)
    
    H_coll = Lplaque * np.sin(np.radians(angle_plaque))
    H_effectif = (H_coll * (Tentree - Tamb) / 2.0 + H_sous * (Tentree - Tamb) + H_sur  * (Tsortie - Tamb)) / (Tamb + 273.15)
    
    Eq3 = v - Ccheminee * np.sqrt(max(0.0, 2.0 * g * H_effectif))
    Eq4 = debit_mair - rho_air * v * Aent
    
    T_moy_chambre = (Tentree + Tsortie) / 2.0
    Q_perte_chambre = UA_chambre * (T_moy_chambre - Tamb)
    
    Eq5 = debit_mair * Cpair * (Tentree - Tsortie) - mdot_evap * Hlv - Q_perte_chambre

    return [Eq1, Eq2, Eq3, Eq4, Eq5]

# FUNCTION TPLAQUE AU COURS DU TEMPS
def Tplaque_t(Tplaque, Tamb, Rayonnement, Q, Aplaque, UL, tau_alpha, rho_ardoise, Cp_ardoise, e_plaque, deltat):
    mCp_plaque = rho_ardoise * Aplaque * e_plaque * Cp_ardoise
    dTdt = (Aplaque * Rayonnement * tau_alpha - Aplaque * UL * (Tplaque - Tamb) - Q) / mCp_plaque
    return Tplaque + deltat * dTdt

# SIMULATION SECHOIR SOLAIRE
def simulation_sechoir(p):
    t = 0.0
    Tamb0, _ = donnees_experimentales(0.0)
    Tplaque = Tplaque_0
    xsolution = np.array([0.0, Tentree_0, Tsortie_0, 0.01, 0.3])
    X = p["X0"]                      
    Theta = 0.0                      
    L_demi = p["e_produit"]/2.0

    results = []
    while t <= p["t_fin"]:
        Tamb_t, Ray_t = donnees_experimentales(t)
        limite_inf = [-np.inf, Tamb_t, Tamb_t, 0.0, 0.0]
        limite_sup = [np.inf, 200.0, 200.0, np.inf, np.inf]
        xsolution = np.clip(xsolution, limite_inf, limite_sup)

        T_local = (xsolution[1]+xsolution[2])/2.0
        D_local = calc_D_banane(T_local)
        MR, dMRdt = calc_dMRdt(Theta, D_local, L_demi)
        Fm_interne = -(p["msec"]/p["Aproduit"]) * (p["X0"]-p["Xeq"]) * dMRdt
        
        if t <= 1e-6 :
            Q,Tentree,Tsortie = 0.0, Tentree_0, Tsortie_0 
            H_coll0 = p["Lplaque"]*np.sin(np.radians(p["angle_plaque"]))
            
            dT0 = max(Tplaque_0 - Tamb_t, 0.5)
            
            H_eff0 = (H_coll0*dT0/2.0 + p["H_sous_grille"]*dT0 + p["H_sur_grille"]*dT0)/(Tamb_t+273.15)
            
            v_est = p["Ccheminee"]*np.sqrt(2.0*p["g"]*max(H_eff0,1e-6))
            
            rho_est = p["Patm"]/(p["Rair"]*(Tamb_t+273.15))
            mdot_air = rho_est*v_est*p["Aent"]
            
            v = v_est
            Q = mdot_air*p["Cpair"]*(Tentree_0-Tamb_t)
            sol = np.array([Q, Tentree, Tsortie, mdot_air, v])
            xsolution = sol

        else:
            escala = [10.0, 30.0, 30.0, 1e-3, 0.3]
            result = least_squares(equations, xsolution, bounds=(limite_inf, limite_sup), x_scale=escala,
                max_nfev=150, ftol=1e-5, xtol=1e-5, gtol=1e-5,
                args=(Tplaque, Tamb_t, Fm_interne, p["Aplaque"], p["Aent"], p["Aproduit"],
                      p["Cpair"], p["g"], p["H_sous_grille"], p["H_sur_grille"], p["Ccheminee"], p["Patm"],
                      p["Rair"], p["angle_plaque"], p["Lplaque"], p["Cpliq_eau"],
                      p["Cpvap_eau"], p["DHvap"], p["Tebu"], p["UA_chambre"], p["HR"]))

            sol = result.x
            xsolution = sol
            Q, Tentree, Tsortie, mdot_air, v = sol

        mdot_evap = Fm_interne * p["Aproduit"]

        Hair_amb, Hair_sortie = calcul_Hair(Tamb_t, p["HR"], p["Patm"], Fm_interne, p["Aproduit"], mdot_air)

        results.append([t, Ray_t, Tamb_t, Tentree, Tplaque, Tsortie, Hair_sortie, mdot_air, mdot_evap, X, v])

        Tplaque = Tplaque_t(Tplaque, Tamb_t, Ray_t, Q, p["Aplaque"], p["UL"], p["tau_alpha"],
                             p["rho_ardoise"], p["Cp_ardoise"], p["e_plaque"], p["deltat"])
       
        X = p["Xeq"] + MR*(p["X0"]-p["Xeq"])
        Theta = Theta + D_local*p["deltat"]

        t += p["deltat"]

    col = ["t [s]", "Ray [W/m2]", "Tamb[°C]", "Tentree[°C]", "Tplaque[°C]", "Tsortie[°C]", "Umidite",
       "Debit_air [Kg/s]", "Debit_evap [Kg/s]", "X", "v [m/s]"]
    
    return pd.DataFrame(results, columns=col)

if __name__ == "__main__":
    df = simulation_sechoir(donnees_diff)
    df.to_excel("results_diff.xlsx", index=False)
    print(df.to_string(index=False))