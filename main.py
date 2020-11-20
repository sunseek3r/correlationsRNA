import tkinter as tk
from tkinter import filedialog
import fastaparser as fp
import numpy as np
import matplotlib.pyplot as plt


def delta (a, alpha):
    if (a == alpha):
        return 1
    else:
        return 0

def p(alpha, sequence):
    res = 0
    for a in sequence:
        res += delta(a, alpha)
    res /= float(len(sequence))
    return res



def correlation(alpha, beta, r, sequence):
    sum = 0
    for i in range (0, len(sequence) - r):
        sum += delta(sequence[i], alpha) * delta(sequence[i+r], beta)
    sum /= float(len(sequence)-r)
    return sum - p(alpha, sequence) * p(beta, sequence)

def corGraph(alpha, beta, r, sequence):
    x = []
    for i in range (1, r):
        x.append(i)
    y = []
    for i in range(1, r):
        y.append(correlation(alpha, beta, i, sequence))

    plt.plot(x, y)
    plt.xlabel("r")
    plt.ylabel("Correlation")
    plt.show()


def funK(i, L, subsum):
    return subsum[i+L] - subsum[i]

def dispersion(L, sequence):
    res = 0
    tempK1 = 0
    tempK2 = 0
    for i in range(0, len(sequence)-L):
        k = funK(i, L, sequence)
        tempK1 += k**2
        tempK2 += k
    res = (tempK1 / float(len(sequence)-L) - (tempK2 / float(len(sequence)-L))**2)
    #print (len(sequence) / float(len(sequence)-L))
    return res 

def sequenceToBinary(zeroCharacter1, zeroCharacter2, zeroCharacter3, sequence): #если надо чтобы два символа переходили в ноль, zeroCharacter3 = ""
    modifiedSequence = ""
    for i in sequence:
        if (i == zeroCharacter1 or i == zeroCharacter2 or i == zeroCharacter3):
            modifiedSequence += "0"
        else:
            modifiedSequence += "1"
    return modifiedSequence


def subsum (sequence):
    res = []
    temp = 0
    for i in sequence:
        temp += int(i)
        res.append(temp)
    return res
    

def dispGraph(L, sequence):
    x = []
    y = []
    i_L = 1
    while True:
        i = int(1.2**i_L)
        if i > L / 5:
            break
        x.append(i)
        y.append(dispersion(i, sequence))
        i_L += 1
    plt.plot(x, y)
    plt.yscale("log")
    plt.xscale("log")
    plt.xlabel("L")
    plt.ylabel("D")
    plt.show()


def corellationPreCalc(sequence, N):
    res = []
    for i in range(0, N+1):
        res.append(correlation("1", "1", i, sequence))
    return res



def memFunction(sequence, N):
    corArray = corellationPreCalc(sequence, N)
    print ("precalc done")
    tempArray = corArray.copy()
    tempArray.pop(0)
    rightPart = np.array(tempArray)
    equationSystem = []
    
    for r in range (1, N+1):
        temp = []
        for r_temp in range (1, N+1):
            temp.append(corArray[abs(r-r_temp)])
        equationSystem.append(temp)
    leftPart = np.array(equationSystem)
    print (len(equationSystem))
    memFunctionArray = np.linalg.solve(leftPart, rightPart)
    return memFunctionArray.tolist()

def memGraph(seqeunce, N):
    memFuncArray = memFunction(seqeunce, N)
    x = []
    for i in range(1, N+1):
        x.append(i)
    plt.plot(x, memFuncArray)
    plt.xlabel("r")
    plt.ylabel("F")
    plt.show()

root = tk.Tk()
root.withdraw()
filepath = filedialog.askopenfilename()
with open(filepath) as fasta_file:
    reader = fp.Reader(fasta_file)
    for seq in reader:
        #print ('ID: ', seq.id)
        #print ('Description:', seq.description)
        #print ('Sequence: ', seq.sequence_as_string())
        sequence = seq.sequence_as_string()
        #print (len(sequence))
        #corGraph("A", "G", 20, seq.sequence_as_string())
        #dispGraph(len(sequence), subsum(sequenceToBinary("A", "G", "", sequence)))
        corGraph("1", "1", 40, sequenceToBinary("A", "G", "", sequence))
        #memGraph(sequenceToBinary("A", "G", "", sequence), 60)
        print ("\n")
    fasta_file.close()
