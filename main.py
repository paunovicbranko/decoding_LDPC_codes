import numpy as np
import random
import matplotlib.pyplot as plt
import time

# --------------------------------------------------------------------------
def calculateSyndrome(codeword, parityCheck):
    # Izračunavanje sindroma parnosti
    syndrome = 0

    for bitPosition in parityCheck:
        syndrome ^= codeword[bitPosition]

    return syndrome



# ---------------------------------------------------------------------------

def calculateSyndromeGDBFM(codeword, parityCheck):
    # Izračunavanje sindroma parnosti
    syndrome = 1

    for bitPosition in parityCheck:
        syndrome *= codeword[bitPosition]

    return syndrome


# --------------------------------------------------------------------------
def findBitFlipPositions(codeword, syndrome, parityCheck):
    # Pronalaženje pozicija bitova koje treba invertovati
    flipPositions = []
    for bitPosition in parityCheck:
        if syndrome != codeword[bitPosition]:
            flipPositions.append(bitPosition)

    return flipPositions

# --------------------------------------------------------------------------

def findParityPositions(LDPCmatrix):
    return [list(np.nonzero(row)[0]) for row in LDPCmatrix]  # vraca se LDPC Code u daljnjem

# --------------------------------------------------------------------------

def recnik(LDPCCode, H):
    dict = {}
    values = []
    for i in range(len(H[0])):
        for parity in LDPCCode:
            if i in parity:
                values.append(LDPCCode.index(parity))
        dict.update({i:values})
        values = []
    return dict

# -------------------------------------------------------------------------

def decisionMaking(codeword, recnik, iterator, LDPCCode):
    syndromeList = []

    parity = recnik[iterator]
    for i in parity:
        syndromeList.append(calculateSyndrome(codeword, LDPCCode[i]))
    
    if syndromeList.count(1) > syndromeList.count(-1):
        return 1
    else:
        return -1


# --------------------------------------------------------------------------

def bitFlippingDecoder(codeword, maxIterations, recnik, LDPCCode):
    # Inicijalizacija
    decodedWord = codeword.copy()
    numIterations = 0

    while numIterations < maxIterations:
        converged = True

        # Iteracija kroz sve parnosti u kodu
        for i in range(len(decodedWord)):
            newSyndrome = decisionMaking(decodedWord, recnik, i, LDPCCode)

            # Provera i ispravka grešaka
            if newSyndrome == -1:
                decodedWord[i] *= -1
            
                converged = False
        numIterations += 1

        if converged: break

    return decodedWord

# -------------------------------------------------------------------------

def gallager_b_algorithm(H, codeword, recnik, LDPCCode, max_iterations=50, alpha=1.0):
    # Inicijalizacija poruka
    rows, cols = len(H), len(H[0])
    lambda_msgs = [[0] * cols for _ in range(rows)]
    rho_msgs = [[1] * cols for _ in range(rows)]
    

    for keys in range(rows):
        for j in LDPCCode[keys]:
            lambda_msgs[keys][j] = codeword[j]

    for iteration in range(max_iterations):
        # Prenos poruka iz čvorova u provere
        n = 1

        for i in range(rows):   # od 0 do 92
            for j in LDPCCode[i]:
                for k in LDPCCode[i]:
                    
                    if j != k:
                        n *= lambda_msgs[i][k]
                rho_msgs[i][j] = n
                n = 1

        pomocni_niz = []

        for i in range(cols):   # od 0 do 155
            for j in recnik[i]:
                for k in recnik[i]:
                    if j != k:
                        pomocni_niz.append(rho_msgs[k][i])
                if pomocni_niz.count(1) >= pomocni_niz.count(-1): lambda_msgs[j][i] = 1
                else: lambda_msgs[j][i] = -1
            pomocni_niz = []

    decoded = []
    temp = []
    for i in range(len(lambda_msgs[0])):
        for j in range(len(lambda_msgs)):
            temp.append(lambda_msgs[j][i])
        if temp.count(1) >= temp.count(-1): decoded.append(1)
        else: decoded.append(-1)
        temp = []

    return decoded


# -------------------------------------------------------------------------

def sign(x):
    return (x > 0) - (x < 0)


# -------------------------------------------------------------------------

def GDBFM(codeword, L, alpha, rho, delta, Itermax, H, LDPCCode, recnik):
    M, N = len(H), len(H[0])
    syndromeList = []

    ln = []
    xn = []

    for i in range(N):
        xn.append(sign(codeword[i]))
        ln.append(L + 1)


    for i in range(Itermax):
        syndromeList = []
        E = []
        
        for parity in LDPCCode:
            syndromeList.append(calculateSyndromeGDBFM(xn, parity))

        # Ukoliko nema losih sindroma, izlazi iz petlje
        if syndromeList.count(-1) == 0: break

        for n in range(N):
            ln[n] = min(ln[n], L) + 1

            # -------------------------------------------------------------------
            # Formula za izracunavanje energije
            temp = 0
            for i in recnik[n]:
                temp += syndromeList[i]

            E.append(alpha * xn[n] * codeword[n] + temp + rho[ln[n]])
            # -------------------------------------------------------------------

        Eth = min(E) + delta

        for n in range(N):
            # Sklonjen uslov za verovatnocu bit-flippinga iz pseudo koda
            if E[n] <= Eth and random.random() < 0.9:
                xn[n] *= -1
                ln[n] = 0

    return xn


# -------------------------------------------------------------------------

def GDBF(codeword, alpha, delta, Itermax, H, LDPCCode, recnik):
    M, N = len(H), len(H[0])
    syndromeList = []

    xn = []

    for i in range(N):
        xn.append(sign(codeword[i]))


    for i in range(Itermax):
        syndromeList = []
        E = []
        
        for parity in LDPCCode:
            syndromeList.append(calculateSyndromeGDBFM(xn, parity))
       
        # Ukoliko nema losih sindroma, izlazi iz petlje
        if syndromeList.count(-1) == 0: break

        for n in range(N):
            
            # Izvlacenje pozicije info bita koji ucestvuju u izracunavanju energije
            
            # -------------------------------------------------------------------
            # Formula za izracunavanje energije
            temp = 0
            for i in recnik[n]:
                temp += syndromeList[i]

            E.append(alpha * xn[n] * codeword[n] + temp )
            # -------------------------------------------------------------------

        Eth = min(E) + delta

        for n in range(N):
            # Sklonjen uslov za verovatnocu bit-flippinga iz pseudo koda
            if E[n] <= Eth and random.random() < 0.9:
                xn[n] *= -1

    return xn


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

matrica = np.loadtxt("matrica.txt")
x = findParityPositions(matrica)        # x JE NIZ U KOJI SE STAVLJAJU POZICIJE JEDINICA U KONTROLNOJ MATRICI
codeword = [0] * len(matrica[0])
codewordCopy = codeword.copy()

for i in range(len(codeword)):      #
    if codeword[i] == 0:            #
        codeword[i] = 1             # PREBACIVANJE NULA U JEDINICE I JEDINICA U MINUS JEDINICE
    else:                           #
        codeword[i] = -1            #

dict = recnik(x, matrica)
verovatnoce_greske = [10,30,50,70,90,110]
FER_bit_flipping = []
FER_GDBF = []
FER_GDBFM =[]
FER_GB =[]

broj_testova = 100000

for j in verovatnoce_greske: # range(10, 130, 20):
    
    lista_gresaka_BF = []                                
    lista_gresaka_GDBF = []             
    lista_gresaka_GDBFM = []
    lista_gresaka_GB = []

    for k in range(broj_testova):
        codewordCopy = codeword.copy() 
        for i in range(len(codewordCopy)): 
            if random.random() < j / 2000:                          
                codewordCopy[i] = codewordCopy[i] * (-1)

# ---------------------------------------------------------------------------

        dekodovana_sekvenca_BF = bitFlippingDecoder(codewordCopy, 4, dict, x)
        dekodovana_sekvenca_GDBF = GDBF(codewordCopy, 1, 0, 4, matrica, x, dict)
        dekodovana_sekvenca_GDBFM = GDBFM(codewordCopy, 3, 1, [3, 3, 2, 2, 1], 0, 7, matrica, x, dict)
        dekodovana_sekvenca_GB = gallager_b_algorithm(matrica, codewordCopy, dict, x, 7)

# ---------------------------------------------------------------------------

        if dekodovana_sekvenca_GB.count(-1) == 0:
            lista_gresaka_GB.append(0)
        else:
            lista_gresaka_GB.append(1)

        if dekodovana_sekvenca_BF.count(-1) == 0:
            lista_gresaka_BF.append(0)
        else:
            lista_gresaka_BF.append(1)
        
        if dekodovana_sekvenca_GDBF.count(-1) == 0:
            lista_gresaka_GDBF.append(0)
        else:
            lista_gresaka_GDBF.append(1)
        
        if dekodovana_sekvenca_GDBFM.count(-1) == 0:
            lista_gresaka_GDBFM.append(0)
        else:
            lista_gresaka_GDBFM.append(1)

    FER_bit_flipping.append(lista_gresaka_BF.count(1) / broj_testova)
    FER_GB.append(lista_gresaka_GB.count(1) / broj_testova)
    FER_GDBF.append(lista_gresaka_GDBF.count(1) / broj_testova)
    FER_GDBFM.append(lista_gresaka_GDBFM.count(1) / broj_testova)



for i in range(len(verovatnoce_greske)):    
    verovatnoce_greske[i] /= 2000
print(len(FER_GB))
plt.plot(verovatnoce_greske, FER_GB, label = 'GB')
plt.plot(verovatnoce_greske, FER_bit_flipping, label = 'BF')
plt.plot(verovatnoce_greske, FER_GDBF, label = 'GDBF')
plt.plot(verovatnoce_greske, FER_GDBFM, label = 'GDBFM')
plt.ylim(10e-6, 10e-2)
#plt.xlim(10e-3, 10e-2)
plt.xscale('log')
plt.yscale('log')
plt.legend()                   
plt.show()

