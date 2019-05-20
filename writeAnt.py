import tellurium as te

def writeBirthDeath(initials, rates, N):
    ant = open('birthdeath.ant','w')
    ant.write("model birthdeath\n")
    rxn_id = 0
    for i in range(0,N):
        ant.write("R"+str(rxn_id+1)+": -> P"+str(i)+"; kb"+str(i)+" \n")
        ant.write("R"+str(rxn_id+2)+": P"+str(i)+" -> ; kd"+str(i)+"*P"+str(i)+" \n")
        ant.write("R"+str(rxn_id+3)+": P"+str(i)+" -> Ss; kout"+str(i)+"*P"+str(i)+" \n")
        ant.write("R"+str(rxn_id+4)+": Ss -> P"+str(i)+"; kin"+str(i)+"*Ss \n\n")
        rxn_id = rxn_id + 4

    for i in range(N):
        ant.write("kb"+str(i)+" = "+str(rates[i][0])+";\n")
        ant.write("kd"+str(i)+" = "+str(rates[i][1])+";\n")
        ant.write("kout"+str(i)+" = "+str(rates[i][2])+";\n")
        ant.write("kin"+str(i)+" = "+str(rates[i][3])+";\n")
        ant.write("P"+str(i)+" = "+str(initials[i][0])+";\n")
        
    ant.write("Ss = "+str(initials[0][1])+";\n")
    ant.write("end")
    ant.close()
    return

def writeAutocatalytic(initials, rates, N):
    ant = open('autocatalytic.ant','w')
    ant.write("model autocatalytic\n")
    rxn_id = 0
    for i in range(0,N):
        ant.write("R"+str(rxn_id+1)+": -> A"+str(i)+"; kb"+str(i)+" \n")
        ant.write("R"+str(rxn_id+2)+": A"+str(i)+" -> A"+str(i)+"+ A"+str(i)+"; ka"+str(i)+"*A"+str(i)+" \n")
        ant.write("R"+str(rxn_id+3)+": A"+str(i)+" -> ; kd"+str(i)+"*A"+str(i)+" \n")
        ant.write("R"+str(rxn_id+4)+": A"+str(i)+" -> Ss; kout"+str(i)+"*A"+str(i)+" \n")
        ant.write("R"+str(rxn_id+5)+": Ss -> A"+str(i)+"; kin"+str(i)+"*Ss \n\n")
        rxn_id = rxn_id + 5

    for i in range(N):
        ant.write("kb"+str(i)+" = "+str(rates[i][0])+";\n")
        ant.write("ka"+str(i)+" = "+str(rates[i][1])+";\n")
        ant.write("kd"+str(i)+" = "+str(rates[i][2])+";\n")
        ant.write("kout"+str(i)+" = "+str(rates[i][3])+";\n")
        ant.write("kin"+str(i)+" = "+str(rates[i][4])+";\n")
        ant.write("A"+str(i)+" = "+str(initials[i][0])+";\n")
        
    ant.write("Ss = "+str(initials[0][1])+";\n")
    ant.write("end")
    ant.close()
    return

def writeDimerRepression(initials, rates, N):
    ant = open('dimerrepression.ant','w')
    ant.write("model dimerrepression\n")
    rxn_id = 0
    for i in range(0,N):
        ant.write("R"+str(rxn_id+1)+": A"+str(i)+" -> A"+str(i)+"+ B"+str(i)+"; kb"+str(i)+"*A"+str(i)+" \n")
        ant.write("R"+str(rxn_id+2)+": B"+str(i)+" + B"+str(i)+" -> C"+str(i)+"; kc"+str(i)+"*B"+str(i)+"*(B"+str(i)+"-1) \n")
        ant.write("R"+str(rxn_id+3)+": A"+str(i)+" + C"+str(i)+" -> D"+str(i)+"; ka"+str(i)+"*A"+str(i)+"*C"+str(i)+" \n")
        ant.write("R"+str(rxn_id+4)+": D"+str(i)+" -> A"+str(i)+" + C"+str(i)+"; kd"+str(i)+"*D"+str(i)+" \n")
        ant.write("R"+str(rxn_id+5)+": B"+str(i)+" -> Ss; kout"+str(i)+"*B"+str(i)+" \n")
        ant.write("R"+str(rxn_id+6)+": Ss -> B"+str(i)+"; kin"+str(i)+"*Ss \n\n")
        rxn_id = rxn_id + 6

    for i in range(N):
        ant.write("kb"+str(i)+" = "+str(rates[i][0])+";\n")
        ant.write("kc"+str(i)+" = "+str(rates[i][1])+";\n")
        ant.write("ka"+str(i)+" = "+str(rates[i][2])+";\n")
        ant.write("kd"+str(i)+" = "+str(rates[i][3])+";\n")
        ant.write("kout"+str(i)+" = "+str(rates[i][4])+";\n")
        ant.write("kin"+str(i)+" = "+str(rates[i][5])+";\n")
        ant.write("A"+str(i)+" = "+str(initials[i][0])+";\n")
        ant.write("B"+str(i)+" = "+str(initials[i][1])+";\n")
        ant.write("C"+str(i)+" = "+str(initials[i][2])+";\n")
        ant.write("D"+str(i)+" = "+str(initials[i][3])+";\n")
        
    ant.write("Ss = "+str(initials[0][4])+";\n")
    ant.write("end")
    ant.close()
    return

def writeDimerization(initials, rates, N):
    ant = open('dimerization.ant','w')
    ant.write("model dimerization\n")
    rxn_id = 0
    for i in range(0,N):
        ant.write("R"+str(rxn_id+1)+": -> P"+str(i)+"; kb"+str(i)+" \n")
        ant.write("R"+str(rxn_id+2)+": P"+str(i)+"+P"+str(i)+" -> D"+str(i)+"; kd"+str(i)+"*P"+str(i)+"*(P"+str(i)+"-1) \n")
        ant.write("R"+str(rxn_id+3)+": P"+str(i)+" -> Ss; kout"+str(i)+"*P"+str(i)+" \n")
        ant.write("R"+str(rxn_id+4)+": Ss -> P"+str(i)+"; kin"+str(i)+"*Ss \n\n")
        rxn_id = rxn_id + 4

    for i in range(N):
        ant.write("kb"+str(i)+" = "+str(rates[i][0])+";\n")
        ant.write("kd"+str(i)+" = "+str(rates[i][1])+";\n")
        ant.write("kout"+str(i)+" = "+str(rates[i][2])+";\n")
        ant.write("kin"+str(i)+" = "+str(rates[i][3])+";\n")
        ant.write("P"+str(i)+" = "+str(initials[i][0])+";\n")
        ant.write("D"+str(i)+" = "+str(initials[i][1])+";\n")
        
    ant.write("Ss = "+str(initials[0][2])+";\n")
    ant.write("end")
    ant.close()
    return

def writeGeneExp(initials, rates, N):
    ant = open('geneexp.ant','w')
    ant.write("model geneexp\n")
    rxn_id = 0
    for i in range(0,N):
        ant.write("R"+str(rxn_id+1)+": D"+str(i)+" -> D"+str(i)+" + M"+str(i)+"; km"+str(i)+"*D"+str(i)+" \n")
        ant.write("R"+str(rxn_id+2)+": M"+str(i)+" -> M"+str(i)+" + P"+str(i)+"; kp"+str(i)+"*M"+str(i)+" \n")
        ant.write("R"+str(rxn_id+3)+": M"+str(i)+" -> ;dm"+str(i)+"*M"+str(i)+" \n")
        ant.write("R"+str(rxn_id+4)+": P"+str(i)+" -> ; dp"+str(i)+"*P"+str(i)+" \n")
        ant.write("R"+str(rxn_id+5)+": P"+str(i)+" -> Ss; kout"+str(i)+"*P"+str(i)+" \n")
        ant.write("R"+str(rxn_id+6)+": Ss -> P"+str(i)+"; kin"+str(i)+"*Ss \n\n")
        rxn_id = rxn_id + 6

    for i in range(N):
        ant.write("km"+str(i)+" = "+str(rates[i][0])+";\n")
        ant.write("kp"+str(i)+" = "+str(rates[i][1])+";\n")
        ant.write("dm"+str(i)+" = "+str(rates[i][2])+";\n")
        ant.write("dp"+str(i)+" = "+str(rates[i][3])+";\n")
        ant.write("kout"+str(i)+" = "+str(rates[i][4])+";\n")
        ant.write("kin"+str(i)+" = "+str(rates[i][5])+";\n")
        ant.write("D"+str(i)+" = "+str(initials[i][0])+";\n")
        ant.write("M"+str(i)+" = "+str(initials[i][1])+";\n")
        ant.write("P"+str(i)+" = "+str(initials[i][2])+";\n")
        
    ant.write("Ss = "+str(initials[0][3])+";\n")
    ant.write("end")
    ant.close()
    return

def writeFeedback(initials, rates, n):
    ant = open('feedback.ant','w')
    ant.write("model feedback\n")
    rxn_id = 0
    for i in range(0,n):
        ant.write("R"+str(rxn_id+1)+": D"+str(i)+" -> D"+str(i)+" + M"+str(i)+"; km"+str(i)+"*D"+str(i)+" \n")
        ant.write("R"+str(rxn_id+2)+": M"+str(i)+" -> M"+str(i)+" + P"+str(i)+"; kp"+str(i)+"*M"+str(i)+" \n")
        ant.write("R"+str(rxn_id+3)+": D"+str(i)+" + P"+str(i)+" -> DP"+str(i)+"; ka"+str(i)+"*D"+str(i)+"*P"+str(i)+" \n")
        ant.write("R"+str(rxn_id+4)+": DP"+str(i)+" -> D"+str(i)+" + P"+str(i)+"; kd"+str(i)+"*DP"+str(i)+" \n")
        ant.write("R"+str(rxn_id+5)+": DP"+str(i)+" -> DP"+str(i)+" + M"+str(i)+"; ki"+str(i)+"*DP"+str(i)+" \n")
        ant.write("R"+str(rxn_id+6)+": M"+str(i)+" -> ;dm"+str(i)+"*M"+str(i)+" \n")
        ant.write("R"+str(rxn_id+7)+": P"+str(i)+" -> ; dp"+str(i)+"*P"+str(i)+" \n")
        ant.write("R"+str(rxn_id+8)+": P"+str(i)+" -> Ss; kout"+str(i)+"*P"+str(i)+" \n")
        ant.write("R"+str(rxn_id+9)+": Ss -> P"+str(i)+"; kin"+str(i)+"*Ss \n\n")
        rxn_id = rxn_id + 9

    for i in range(n):
        ant.write("km"+str(i)+" = "+str(rates[i][0])+";\n")
        ant.write("kp"+str(i)+" = "+str(rates[i][1])+";\n")
        ant.write("ka"+str(i)+" = "+str(rates[i][2])+";\n")
        ant.write("kd"+str(i)+" = "+str(rates[i][3])+";\n")
        ant.write("ki"+str(i)+" = "+str(rates[i][4])+";\n")
        ant.write("dm"+str(i)+" = "+str(rates[i][5])+";\n")
        ant.write("dp"+str(i)+" = "+str(rates[i][6])+";\n")
        ant.write("kout"+str(i)+" = "+str(rates[i][7])+";\n")
        ant.write("kin"+str(i)+" = "+str(rates[i][8])+";\n")
        ant.write("D"+str(i)+" = "+str(initials[i][0])+";\n")
        ant.write("M"+str(i)+" = "+str(initials[i][1])+";\n")
        ant.write("P"+str(i)+" = "+str(initials[i][2])+";\n")
        ant.write("DP"+str(i)+" = "+str(initials[i][3])+";\n")
        
    ant.write("Ss = "+str(initials[0][4])+";\n")
    ant.write("end")
    ant.close()
    return

def writeSimpleFeedback(initials, rates, n):
    ant = open('simplefeedback.ant','w')
    ant.write("model simplefeedback\n")
    rxn_id = 0
    for i in range(0,n):
        ant.write("R"+str(rxn_id+1)+": D"+str(i)+" -> D"+str(i)+" + P"+str(i)+"; kp"+str(i)+"*D"+str(i)+" \n")
        ant.write("R"+str(rxn_id+2)+": D"+str(i)+" + P"+str(i)+" -> DP"+str(i)+"; ka"+str(i)+"*D"+str(i)+"*P"+str(i)+" \n")
        ant.write("R"+str(rxn_id+3)+": DP"+str(i)+" -> D"+str(i)+" + P"+str(i)+"; kd"+str(i)+"*DP"+str(i)+" \n")
        ant.write("R"+str(rxn_id+4)+": DP"+str(i)+" -> DP"+str(i)+" + P"+str(i)+"; ki"+str(i)+"*DP"+str(i)+" \n")
        ant.write("R"+str(rxn_id+5)+": P"+str(i)+" -> ; dp"+str(i)+"*P"+str(i)+" \n")
        ant.write("R"+str(rxn_id+6)+": P"+str(i)+" -> Ss; kout"+str(i)+"*P"+str(i)+" \n")
        ant.write("R"+str(rxn_id+7)+": Ss -> P"+str(i)+"; kin"+str(i)+"*Ss \n\n")
        rxn_id = rxn_id + 7

    for i in range(n):
        ant.write("kp"+str(i)+" = "+str(rates[i][0])+";\n")
        ant.write("ka"+str(i)+" = "+str(rates[i][1])+";\n")
        ant.write("kd"+str(i)+" = "+str(rates[i][2])+";\n")
        ant.write("ki"+str(i)+" = "+str(rates[i][3])+";\n")
        ant.write("dp"+str(i)+" = "+str(rates[i][4])+";\n")
        ant.write("kout"+str(i)+" = "+str(rates[i][5])+";\n")
        ant.write("kin"+str(i)+" = "+str(rates[i][6])+";\n")
        ant.write("D"+str(i)+" = "+str(initials[i][0])+";\n")
        ant.write("P"+str(i)+" = "+str(initials[i][1])+";\n")
        ant.write("DP"+str(i)+" = "+str(initials[i][2])+";\n")
        
    ant.write("Ss = "+str(initials[0][3])+";\n")
    ant.write("end")
    ant.close()
    return