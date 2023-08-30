#Import libraries:
import sympy as sym
import os
import numpy as np
import math  #For importing the value of Pi
import pandas as pd


def calc_ionic_strength_TE(concTE, settings):
    """
    This function calculates the ionic strength of a Tris-EDTA buffer solution.
    Input:
        concTE [Array] - Concentration of TE buffer [x] in the range of 0.001 to 50
        Inside settings dictionary:
            with_BME [Bool] - Includes betamercaptoethanol (BME) in the calculation
            c_BME_percent [Float] - Concentration of BME [% (v/v) or mL BME/mL total]
            save_to_file [Bool] - Saves the data frame to a csv file
            showAllSol [Bool] - Shows all the solutions, even the invalid ones in each loop
            showEq [Bool] - Displays the equations used to solve the system  in each loop
            showConc [Bool] - Displays the concentrations found in each loop
            showCoeffs [Bool] - Displays the activity coefficients found in each loop
    """

    df = get_species_data_frame(with_BME = settings['with_BME'])
    
    df_all, I_s, pH_s = iterate_func(df, concTE, settings)
    
    if settings['save_to_file']:
        save_to_file_fun(df_all, concTE, I_s, pH_s, settings, extra_str='')



def get_species_data_frame(with_BME = False, c_BME_percent = 0):
    """Get the data frame with the species data"""

    #Get dissociation constants
    pK_w, pK_b_Tris, pK_a_EDTA, pK_a_B = get_diss_constants()

    #Get BME concentration
    if with_BME:
        C_B = get_C_BME(c_BME_percent)


    #Add the dissociation constants, the chemical species and valences into lists:
    pK_a = [0, 0, pK_a_EDTA[2], pK_a_EDTA[1],pK_a_EDTA[3], pK_a_B,0, pK_w, 0]
    pK_b = [pK_b_Tris, '', '', '', '', '','', '', '']
    species = ['[TH^{+}]','[T]','[HE^{3-}]','[H_2 E^{2-}]','[E^{4-}]','[B^{-}]','[HB]','[H^{+}]','[OH^{-}]']
    species_f = ['$[TH^{+}]$','$[T]$','$[HE^{3-}]$','$[H_2 E^{2-}]$','$[E^{4-}]$','$[B^{-}]$','$[HB]$','$[H^{+}]$','$[OH^{-}]$']
    keys = ['THP', 'T', 'HE3M', 'H2E2M', 'E4M', 'BM', 'HB', 'HP', 'OHM']
    valences = np.array([1,0,-3,-2,-4,1,0,1,-1])   

    #Add the equilibrium constants
    K_a = np.zeros(len(pK_a))
    for i,x in enumerate(pK_a):
        y = float(x)
        if y > 0:
            K_a[i] = 10**-y
        else:
            K_a[i] = np.nan
    K_b_Tris = 10**-pK_b_Tris
    K_b = [K_b_Tris, '', '', '', '', '','', '', '']
        
    g = np.ones((len(species)))  #initiate acticity coefficients: start with all equal to 1
    c = np.zeros((len(species)))  #initiate list of concentrations
    df = pd.DataFrame({'species_f':species_f,'species':species, 'val':valences,'pK_a':pK_a, 'K_a':K_a, 'pK_b':pK_b, 'K_b':K_b, 'g':g, 'c':c})
    df.index = keys #set the short names as data frame keys
    df.index.name = 'inx'

    if not with_BME: #If BME is not included in the calculation, drop the associated terms from the data frame  
        df = df.drop("HB", axis=0)
        df = df.drop("BM", axis=0)
    print('Initial values:')
    print(df.head(10)) 
    return df

def get_C_BME(c_BME_percent):
    """Get the concentration of BME [M or mol/L]"""

    c_BME_percent = 0.03 # Concentration of BME [% (v/v) or mL BME/mL total]
    M_BME = 78.13# Molar mass for BME [g/mol], Source: Sigma-Aldrich
    rho_BME = 1.114 #Density for BME [g/mL] at 25 °C, Source: Thomas Scientific
    C_B = c_BME_percent * rho_BME * 1e3 / M_BME #Concentration of BME [M or mol/L] unit calc = [mL BME / mL * g/mL * 1000 * mol/g]
    print('Conc. BME'+f' = {c_BME_percent*100}% (v/v) or {np.round(C_B*1e3,3)} mM')
    return C_B

def get_diss_constants():
    """Get the dissociation constants"""

    pK_w = 14.0 #Dissociation constant for water at T=25C
    pK_b_Tris = 5.94 #Tris base dissociation constant, For T=25C from [Iarko 2015]
    pK_a_EDTA = np.array([1.99,2.67,6.16,10.26]) #EDTA acid dissociation constants, For T=25C from [Iarko 2015]
    pK_a_B = 9.6 #BME acid dissociation constant, For T=25C from [Iarko 2015]
    return pK_w, pK_b_Tris, pK_a_EDTA, pK_a_B

def iterate_func(df, concTE, settings):
    """Iterates the calculation until the ionic strength converges"""

    A = 0.51 #Constant for the Davies equation
    c_BME_percent = settings['c_BME_percent'] #BME concentration in percent
    #Get dissociation constants
    pK_w, pK_b_Tris, pK_a_EDTA, pK_a_B = get_diss_constants()


    #Get BME concentration
    if settings['with_BME']:
        C_B = get_C_BME(settings['c_BME_percent'])

    n_s = len(concTE) #number of samples
    # n_s = 1
    I_s = [] #initialize array of ionic strengths
    pH_s = []
    #for j in np.arange(4,5,1):
    for j in range(n_s):
        if settings['with_BME']:
            print('\nsample '+str(j+1)+': '+str(concTE[j])+'x TE + '+str(100*c_BME_percent)+'% BME')
        else:
            print('\nsample '+str(j+1)+': '+str(concTE[j])+'x TE')
        
        #Concentrations
        N = concTE[j] #Numeral prefix to TE or TE buffers
        C_T = N * 0.01 #Concentration of Tris [M]
        C_E = N * 0.001 #Concentration of EDTA [M]
        print('Conc. Tris'+f' = {np.round(C_T*1e3,3)} mM')
        print('Conc. EDTA'+f' = {np.round(C_E*1e3,3)} mM')
        df['N_TE'] = N
        I_list = [] #List of ionic strengths
        n = len(df) #number of equations/variables

        for k in range(7):
            if k == 0: #initiate acticity coefficients: start with all equal to 1
                df['g'] = 1
                g_list = df['g'].to_list()

            sym.init_printing() #Use sympy's "pretty printing"

            # Definition of symbols (unknown concentrations):
            # Water:
            c_HP = sym.Symbol(df.loc['HP','species'])
            c_OHM = sym.Symbol(df.loc['OHM','species'])
            c_THP = sym.Symbol(df.loc['THP','species'])
            c_T = sym.Symbol(df.loc['T','species'])
            c_HE3M = sym.Symbol(df.loc['HE3M','species'])
            c_H2E2M = sym.Symbol(df.loc['H2E2M','species'])
            c_E4M = sym.Symbol(df.loc['E4M','species'])

            print('\n')
            print('<Loop no. '+str(k+1)+'>')
            if settings['showCoeffs']:
                print('Initial values of activity coefficients are:')
                print(df.loc[:,['species', 'g']])

            #Updating activity coefficients (with g as in gamma):
            g_THP = df.loc['THP','g']
            g_T = df.loc['T','g']
            g_HE3M = df.loc['HE3M','g']
            g_H2E2M = df.loc['H2E2M','g']
            g_E4M = df.loc['E4M','g']
            g_HP = df.loc['HP','g']
            g_OHM = df.loc['OHM','g']

            #Set of equations:

            #Tris:
            eq1 = sym.Eq(c_THP*g_THP * c_OHM*g_OHM / (c_T*g_T), 10**(-pK_b_Tris))
            eq6 = sym.Eq(c_T + c_THP, C_T)

            #EDTA:
            eq2 = sym.Eq(c_HE3M*g_HE3M * c_HP*g_HP / (c_H2E2M*g_H2E2M), 10**(-pK_a_EDTA[2]))
            eq3 = sym.Eq(c_E4M*g_E4M*c_HP*g_HP / (c_HE3M*g_HE3M),10**(-pK_a_EDTA[3]))  
            eq7 = sym.Eq(c_H2E2M + c_HE3M + c_E4M, C_E)
            
            #Water:
            eq5 = sym.Eq(c_HP*g_HP * c_OHM*g_OHM, 10**(-pK_w))
            
            if settings['with_BME']: #With BME
                g_BM = df.loc['BM','g']
                g_HB = df.loc['HB','g']
                c_BM = sym.Symbol(df.loc['BM','species'])        
                c_HB = sym.Symbol(df.loc['HB','species'])
                eq4 = sym.Eq(c_BM*g_BM * c_HP*g_HP / (c_HB*g_HB),10**(-pK_a_B))
                eq8 = sym.Eq(c_HB + c_BM, C_B)
                eq9 = sym.Eq(c_THP + c_HP, c_OHM + 2*c_H2E2M + 3*c_HE3M + 4*c_E4M + c_BM)
                system = [eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9]
                
                #Solve the set of equation using Sympy's solver:
                s = list(sym.solve(system,(c_THP, c_T, c_HE3M, c_H2E2M, c_E4M, c_BM, c_HB, c_HP, c_OHM),manual=True, dict=True))
            else: #Without BME:
                eq9 = sym.Eq(c_THP + c_HP, c_OHM + 2*c_H2E2M + 3*c_HE3M + 4*c_E4M)
                system = [eq1, eq2, eq3, eq5, eq6, eq7, eq9]
                
                #Solve the set of equation using Sympy's solver:
                s = list(sym.solve(system,(c_THP, c_T, c_HE3M, c_H2E2M, c_E4M, c_HP, c_OHM),manual=True, dict=True))
                #     if you want include warnings, add warn=True
            
            if settings['showEq']: #For displaying the equations
                show_eq_fun(system)
            if settings['showAllSol']: # For displaying the all possible solutions (incl. negative ones)
                showAllSol_fun(s)
                        
            #Extract the information from the solution:
            n_sol = 0 # Number of solutions with positive values
            for i in range(len(s)):
                a = dict(s[i])
                keys = a.keys()
                vals = list(a.values())
                if all(x > 0 for x in vals):
                    #print(f'Solution {i+1} has no negative variables')
                    inx = i #the index of the solution with positive values
                    n_sol = n_sol+1
                    
            if settings['showConc']: #For displaying the found concentrations:
                print('\nFound '+str(n_sol)+' set of solutions where all concentrations > 0:')
            if n_sol > 1:
                print(f'ERROR: found {n_sol} solutions, should only be 1!!!')
                
            # The solution with positive values
            a = dict(s[inx]) #Make the solution into a dictionary
            keys = list(a.keys())
            vals = list(a.values()) #concentration values

            # Calculate the ionic strength, I from the ionic species concentrations and valances
            c = np.zeros(n) #initiate the array
            I = 0 #initiate the ionic strength value
            for i in range(len(a)):
                key = str(keys[i])
                df_inx = df[df['species'] == key].index
                conc = float(vals[i])
                df.loc[df_inx, 'c'] = conc
                prefix = float((df.loc[df_inx, 'val']**2))
                cont = conc*prefix
                I = I + cont
                print(f'i = {i}, I = {np.round(0.5*I*1000,3)} mM, {key}, contribution = \t0.5*{prefix}*{np.round(1000*conc,4)} mM = \t{np.round(0.5*1000*cont,4)} mM ')
            I = I / 2
            I = float(I) # Convert from sympy float data structure
            I_list.append(I) #Add ionic strength to the list to observe the convergence
            
            if settings['showConc']:
                print('\nA clearer picture of the found concentrations:')
                print(df.loc[:,['species', 'c']]) #display the concentrations
                print('\n')
            pH = -math.log10(df.loc['HP','c'])
            
            if k > 0:
                I_ratio = I/I_list[k-1]
                print(f'I = {np.round(I*1e3,3)} mM {np.round(I_ratio,5)} x I_list[k-1) and the pH is {np.round(pH,2)}')
                if abs(I_ratio - 1) < 0.01:
                    I_s.append(I)                
                    pH_s.append(-math.log10(df.loc['HP','c']))
                    df['I'] = I * 1000
                    df['pH'] = pH
                    if j == 0:
                        df_all = df
                    else:
                        df_all = df_all.append(df, ignore_index=True)
                    break    # break the loop if the difference between last value is too small.
            else:
                print(f'I = {np.round(I*1e3,3)} mM and the pH is {np.round(pH,2)}')        

            df['g'] = 0#set all activity coefficients to zero.
            
            # run the Davies Equation for all ionic species:
            for index, row in df.iterrows():
                df.loc[index,'g'] =  10**(-A*(row['val']**2)*(np.sqrt(I)/(1+np.sqrt(I))-0.3*I))
                #a = row['val']**2
                #print(f'{index}: z^2 = {a}')
                #print(' '+str(row['species'])+', g = \t'+ str(np.round(row['g'],5)))
        #print(df.loc[:,['species', 'g']])    
        if settings['with_BME']:
            print('for a buffer of '+str(N)+'x TE and '+str(c_BME_percent*100)+'% BME.')
        else:
            print('for a buffer of '+str(N)+'x TE.')
    return df_all, I_s, pH_s

def show_eq_fun(system):
    """Displaying of the equations used to solve the system  in each loop"""
    print(f'\nA set of {len(system)} equations:')
    for x in system:
        print(x)

def showAllSol_fun(s):
    """Displaying of the all possible solutions (incl. negative ones)"""

    print(f'\nFound {len(s)} sets of solutions:')
    for i in range(len(s)):
        print(f'Solution {i+1}:')
        dict_s = dict(s[i])
        for key, value in dict_s.items():            
            print('  ',key, ':', value)

# Save the data:
def get_I_contribution(c,val):
    """Calculate the contribution to the ionic strength from a species"""
    return 0.5*c*(val**2)

def save_to_file_fun(df_all, concTE, I_s, pH_s, settings, extra_str=''):
    """Save the data to two csv files.
    One file contains the concentrations of all species (df_all) and the other contains the ionic strength and pH. (df)"""
    c_BME_percent = settings['c_BME_percent'] #BME concentration in percent
    basepath = settings['basepath'] #Path to the folder where the results should be saved

    n_s = len(concTE) #Number of samples
    df_all["c_mM"] = df_all['c']*1000 
    df_all["cont_mM"] = df_all.apply(lambda df_all: get_I_contribution(df_all.c_mM, df_all.val), axis=1)
    df_all["c_BME_percent"] = c_BME_percent*100 
    
    if settings['with_BME']:
        df_BME = pd.DataFrame({'conc_TE':concTE,'conc_BME':np.ones(n_s)*c_BME_percent, 'I':I_s, 'pH':pH_s})
        df_BME['I'] = df_BME['I'] * 1000
        str_BME = str(c_BME_percent)
        str_BME = str_BME.replace('.', '-')

        #Save to file:
        file_name_df = 'df_TE_'+str_BME+'_BME.csv'
        file_name_df_all = 'df_all_TE_'+str_BME+'_BME.csv'
        df_BME.to_csv(os.path.join(basepath, file_name_df))
        df_all.to_csv(os.path.join(basepath, file_name_df_all))
        print(df_BME.head(10))
    else:
        df_0_BME = pd.DataFrame({'conc_TE':concTE,'conc_BME':np.zeros(n_s), 'I':I_s, 'pH':pH_s})
        df_0_BME['I'] = df_0_BME['I'] * 1000
        
        #Save to file:
        file_name_df = extra_str+'df_TE_0_BME.csv'
        file_name_df_all = extra_str+'df_all_TE_0_BME.csv'
        df_0_BME.to_csv(os.path.join(basepath, file_name_df))
        df_all.to_csv(os.path.join(basepath, file_name_df_all))
        print(df_0_BME.head(10))
    print('Summary:')
    for j in range(n_s):
        if settings['with_BME']:
            print('Conc. '+str(j+1)+': '+str(concTE[j])+'x TE + '+str(100*c_BME_percent)+'% BME: I = '+str(np.round(1000*I_s[j],3))+' mM'+f' and the pH is {np.round(pH_s[j],2)}')
        else:
            print('Conc. '+str(j+1)+': '+str(concTE[j])+'x TE'+f' and the pH is {np.round(pH_s[j],2)}')