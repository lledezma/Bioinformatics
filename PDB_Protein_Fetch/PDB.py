import math
import requests
import os

class PDB:

    def __init__(self, pdb):
            self.protein = pdb[:4]
            self.file = self.protein+'.txt'
            #try to read a protein from a local file
            if (os.path.isfile(self.file)):
                self.pdb = open(self.file, 'r') #opens pdb file 3UTS.text
                self.pdb.close()
            #if the protein file doesnt exist, we fetch it from online
            else:
                url = 'https://files.rcsb.org/view/' + self.protein + '.pdb'
                myfile = requests.get(str(url).rstrip())        #copy the information of the url webpage
                open(self.file, 'a').write(str(myfile.content).replace("\\n","\n"))  #dump the information into <protein name>.txt
            self.chain = pdb[4:].strip()           #read the last character of the input which is the chain


    def q2(self):
        aminoAcids = []
        aminoAcidscount = "0"
        missingcount = 0
        self.pdb = open(self.file, 'r')
        lines = self.pdb.readlines()
        self.pdb.close()
        for line in lines:
# how many are missing structuralinformation
            if(line[:6] == 'REMARK' and line[7:10].strip() == '465' and line[19:20].strip() == self.chain and line[20:27].strip() != 'SSSEQI' and line[10:15].strip() == ""):
                missingcount+=1
# Show the number of amino acids in the chain
            if(line[:6] == 'SEQRES' and line[10:12].strip() == self.chain):
                aminoAcidscount = (line[12:17].strip())
                break   
        print("total number of amino acids:", aminoAcidscount)
        print("total number of missing amino acids:", missingcount)



    def q3(self,aaSequenceNumber):
        countAtoms = 0
        self.pdb = open(self.file, 'r')
        lines = self.pdb.readlines()
        self.pdb.close()
        sidecountAtoms = 0
        atomCB = 0
# Show the number of atoms in a given aa (aa sequence number should be provided)
        for line in lines:
            if (line[:4] == 'ATOM') and (line[21:22] == self.chain):
                if(line[23:26].strip() == str(aaSequenceNumber)):
                    countAtoms += 1
# Show how many atoms in the side chain
                #ATOM                           #chain                          #sequence number
            if (line[:4] == 'ATOM') and (line[21:22] == self.chain) and (line[23:26].strip() == str(aaSequenceNumber)):
                    #atom name
                if(line[13:16].strip() == 'CB'):
                    atomCB = int(line[7:11].strip())
                    sidecountAtoms +=1
                if(atomCB > 0):
                    if (int(line[7:11].strip()) > atomCB ):
                        sidecountAtoms+=1
        print("Number of atoms:", countAtoms)
        print("number of atoms in the side chain:", sidecountAtoms)



    def q4(self):
        helixcount = 0
        helixType = {'1': 0, '2': 0, '3': 0, '4': 0, '5': 0, '6': 0,
                     '7': 0, '8': 0, '9': 0, '10': 0}
        helixClass = ['Right-handed alpha', 'Right-handed omega', 'Right-handed pi', 'Right-handed gamma', 'Right-handed 3 - 10',
        'Left-handed alpha', 'Left-handed omega', 'Left-handed gamma', '2 - 7 ribbon/helix', 'Polyproline']
        self.pdb = open(self.file, 'r')
        lines = self.pdb.readlines()
        self.pdb.close()
        for line in lines:
            if(line[:6].strip() == "HELIX") and (line[19:20].strip() == self.chain):
# Show the total number of helices
                helixcount +=1
# how many helices for each type
                helixType[line[38:40].strip()] +=1
        print("Total number of Helices: " , helixcount)
        print("\n'TYPE'", " : ", "NUMBER")
        for clss,number in zip(helixClass,helixType):
            print(number, "-",clss, ": ", helixType[number])


    def q5(self,ID):
        sheetcount = 0
        self.pdb = open(self.file, 'r')
        lines = self.pdb.readlines()
        self.pdb.close()
        for line in lines:
            if((line[:6].strip() == "SHEET") and (line[21:22].strip() == self.chain) and line[11:14].strip() == ID):
                sheetcount = int(line[14:16].strip())
                break
        if(sheetcount == 0 ):
            print("Sheet does not exist")
        else:
            print("Number of strands that form sheet",ID, ":", sheetcount)


    def q6(self, atom):
        aanumber, atomname = atom.split('.')
        X = ""
        Y = ""
        Z = ""
        self.pdb = open(self.file, 'r')
        lines = self.pdb.readlines()
        self.pdb.close()
        for line in lines:
            if(  line[:4] == 'ATOM'  and line[21:22] == self.chain and line[12:16].strip() == atomname and  line[23:26].strip() == aanumber):
# Return the coordinates of a particular atom
                X = line[30:38].strip()
                Y = line[38:46].strip()
                Z = line[46:54].strip()
                break
        if(X == "" and Y == "" and Z == ""):
            print("No coordinates for the given atom")
        else:
            print("x: ", X)
            print("y: ", Y)
            print("z: ", Z)


    def q7(self, AA):
        sheetID,number, position = AA.split('.')
        self.pdb = open(self.file, 'r')
        lines = self.pdb.readlines()
        self.pdb.close()
        for line in lines:
            if (line[:6].strip() == "SHEET" and line[11:14].strip() == sheetID and line[21:22].strip() == self.chain and line[7:10].strip() == number  and position == 'begin'):
                print("name of the amino acid at the beginning of a given strand:", line[17:20])
                return
            elif (line[:6].strip() == "SHEET" and line[11:14].strip() == sheetID and line[32:33].strip() == self.chain and line[7:10].strip() == number  and position == 'end'):
                print("name of the amino acid at the end of a given strand:",line[28:31])
                return
        print("No amino acid to show.")



    def q8(self, atom1, atom2):
        aanumber1, atomname1 = atom1.split('.')
        X1 = ""
        Y1 = ""
        Z1 = ""
        aanumber2, atomname2 = atom2.split('.')
        X2 = ""
        Y2 = ""
        Z2 = ""
        self.pdb = open(self.file, 'r')
        lines = self.pdb.readlines()
        self.pdb.close()
        for line in lines:
            if(  line[:4] == 'ATOM'  and line[21:22] == self.chain and line[12:16].strip() == atomname1 and line[23:26].strip() == aanumber1):
                X1 = line[30:38].strip()
                Y1 = line[38:46].strip()
                Z1 = line[46:54].strip()
            if(  line[:4] == 'ATOM'  and line[21:22] == self.chain and line[12:16].strip() == atomname2 and line[23:26].strip() == aanumber2):
                X2 = line[30:38].strip()
                Y2 = line[38:46].strip()
                Z2 = line[46:54].strip()
        try:
            distance = math.sqrt(   ((float(X2)-float(X1))**2) + ((float(Y2)-float(Y1))**2 ) + ((float(Z2)-float(Z1))**2 )  )
        except:
            print("One or more of the given atoms do not exist in the chain.")
            return
        print("Distance:", round(distance,3))



    def q9(self, aa):
        aaNumber, angletype =  aa.split('.')
        Cphi = []
        N = []
        CA = []
        C = []
        Npsi = []
        self.pdb = open(self.file, 'r')
        lines = self.pdb.readlines()
        self.pdb.close()
        for line in lines:
            if(  line[:4] == 'ATOM'  and line[21:22] == self.chain and line[12:16].strip() == 'C' and int(line[23:26].strip()) == int(aaNumber)-1):
                Cphi.append(float(line[30:38].strip()))
                Cphi.append(float(line[38:46].strip()))
                Cphi.append(float(line[46:54].strip()))
            if(  line[:4] == 'ATOM'  and line[21:22] == self.chain and line[12:16].strip() == 'N' and line[23:26].strip() == aaNumber):
                N.append(float(line[30:38].strip())) #append x
                N.append(float(line[38:46].strip())) #append y
                N.append(float(line[46:54].strip())) #append z
            if(  line[:4] == 'ATOM'  and line[21:22] == self.chain and line[12:16].strip() == 'CA' and line[23:26].strip() == aaNumber):
                CA.append(float(line[30:38].strip())) #append x
                CA.append(float(line[38:46].strip())) #append y
                CA.append(float(line[46:54].strip())) #append z
            if(  line[:4] == 'ATOM'  and line[21:22] == self.chain and line[12:16].strip() == 'C' and line[23:26].strip() == aaNumber):
                C.append(float(line[30:38].strip())) #append x
                C.append(float(line[38:46].strip())) #append y
                C.append(float(line[46:54].strip())) #append z
            if(  line[:4] == 'ATOM'  and line[21:22] == self.chain and line[12:16].strip() == 'N' and int(line[23:26].strip()) == int(aaNumber)+1):
                Npsi.append(float(line[30:38].strip()))
                Npsi.append(float(line[38:46].strip()))
                Npsi.append(float(line[46:54].strip()))
        try:
            if angletype == 'phi':
                # A1 = y1 ( z2 - z3 )        +      y2 ( z3 - z1 )            +   y3 ( z1 - z2 )
                A1 =  Cphi[1]*(N[2]-CA[2])   +      N[1]*( CA[2]-Cphi[2] )    +   CA[1]*(Cphi[2]- N[2])
                # B1 = z1 ( x2 - x3 )        +      z2 (   x3 - x1 )          +      z3 ( x1 - x2 )
                B1 = Cphi[2]*(N[0]-CA[0])    +      N[2]*( CA[0]-Cphi[0] )    +   CA[2]*(Cphi[0]- N[0])
                # C1= x1 ( y2 - y3 )         +       x2 ( y3 - y1 )           +  x3 ( y1 - y2 )
                C1 = Cphi[0]*(N[1]-CA[1])    +      N[0]*( CA[1]-Cphi[1] )    +   CA[0]*(Cphi[1]- N[1])
                # ////////////////////////////
                # A1 = y2 ( z3 - z4 )        +      y3 ( z4 - z2 )          +   y4 ( z2 - z3 )
                A2 =  ( N[1]*(CA[2]-C[2]) )  +      (CA[1]*( C[2]-N[2] ))   +  (C[1]*(N[2]- CA[2]))
                # B2 = z2 ( x3 - x4 )        +      z3 (   x4 - x2 )        +   z4 ( x2 - x3 )
                B2 = (N[2]*(CA[0]-C[0]))     +      (CA[2]*( C[0]-N[0] ))   +  (C[2]*(N[0]- CA[0]))
                # C2 = x2 ( y3 - y4 )        +       x3 ( y4 - y2 )         +  x4 ( y2 - y3 )
                C2 = (N[0]*(CA[1]-C[1]))     +      (CA[0]*( C[1]-N[1] ))   +  (C[0]*(N[1]- CA[1]))
                # //////////////////////////
                angle = (A1*A2 + B1*B2 + C1*C2)/ (math.sqrt(A1*A1 + B1*B1 + C1*C1) * math.sqrt(A2*A2 + B2*B2 + C2*C2))
                # //////////////////  vNormal = the cross product of v1 and v2
                V1 = [Cphi[0]-N[0], Cphi[1]-N[1], Cphi[2]-N[2]]
                V2 = [CA[0]-N[0], CA[1]-N[1], CA[2]-N[2]]
                V3 = [C[0]-CA[0], C[1]-CA[1], C[2]-CA[2]]
                vNormal = []
                vNormal.append( (V1[1]*V2[2] - V1[2]*V2[1]) )
                vNormal.append( -(V1[0]*V2[2] - V1[2]*V2[0]) )
                vNormal.append( (V1[0]*V2[1] - V1[1]*V2[0]) )
                dotProduct = (V3[0]*vNormal[0] + V3[1]*vNormal[1] + V3[2]*vNormal[2])
                # //////////////// final angle /// (acos(angle)* 180)/PI
                finalangle = round(((math.acos(angle)*180) / math.pi),4)
                if(dotProduct > 0):
                    print("phi angle:",-1*finalangle)
                else:
                    print("phi angle:",finalangle)
    # ////////////////////////////////////////////////////////////////////////////////////////////////
            if angletype == 'psi':
                # A1 = y1 ( z2 - z3 )     +      y2 ( z3 - z1 )         +    y3 ( z1 - z2 )
                A1 =  N[1]*(CA[2]-C[2])   +      CA[1]*( C[2]-N[2] )    +   C[1]*(N[2]- CA[2])
                # B1 = z1 ( x2 - x3 )     +      z2 (   x3 - x1 )       +    z3 ( x1 - x2 )
                B1 = N[2]*(CA[0]-C[0])    +      CA[2]*( C[0]-N[0] )    +   C[2]*(N[0]- CA[0])
                # C1= x1 ( y2 - y3 )      +       x2 ( y3 - y1 )        +     x3 ( y1 - y2 )
                C1 = N[0]*(CA[1]-C[1])    +      CA[0]*( C[1]-N[1] )    +   C[0]*(N[1]- CA[1])
                # ////////////////////////////
                # A2 = y2 ( z3 - z4 )        +        y3 ( z4 - z2 )          +   y4 ( z2 - z3 )
                A2 =  (CA[1]*(C[2]-Npsi[2])) +      (C[1]*( Npsi[2]-CA[2] ))  +  (Npsi[1]*(CA[2]- C[2]))
                # B2 = z2 ( x3 - x4 )        +        z3 (   x4 - x2 )        +   z4 ( x2 - x3 )
                B2 = (CA[2]*(C[0]-Npsi[0]))  +      (C[2]*( Npsi[0]-CA[0] ))  +  (Npsi[2]*(CA[0]- C[0]))
                # C2 = x2 ( y3 - y4 )        +         x3 ( y4 - y2 )         +   x4 ( y2 - y3 )
                C2 = (CA[0]*(C[1]-Npsi[1]))  +      (C[0]*( Npsi[1]-CA[1] ))  +  (Npsi[0]*(CA[1]- C[1]))
                # //////////////////////////
                angle = (A1*A2 + B1*B2 + C1*C2)/ (math.sqrt(A1*A1 + B1*B1 + C1*C1) * math.sqrt(A2*A2 + B2*B2 + C2*C2))
                # //////////////////  vNormal = the cross product of v1 and v2
                V1 = [CA[0]-N[0], CA[1]-N[1], CA[2]-N[2]]
                V2 = [C[0]-CA[0], C[1]-CA[1], C[2]-CA[2]]
                V3 = [Npsi[0]-C[0], Npsi[1]-C[1], Npsi[2]-C[2]]
                vNormal = []
                vNormal.append( (V1[1]*V2[2] - V1[2]*V2[1]) )
                vNormal.append( -(V1[0]*V2[2] - V1[2]*V2[0]) )
                vNormal.append( (V1[0]*V2[1] - V1[1]*V2[0]) )
                dotProduct = (V3[0]*vNormal[0] + V3[1]*vNormal[1] + V3[2]*vNormal[2])
                # //////////////// final angle /// (acos(angle)* 180)/PI
                finalangle = round(((math.acos(angle)*180) / math.pi),4)
                if(dotProduct < 0):
                    print("psi angle:",-1*finalangle)
                else:
                    print("psi angle:",finalangle)
        except:
            print("Give amino acid or chain doesn't exist")



PDBf = PDB('3UTSA')          # 1. Able to read the protein from local file and fetch it from online.

# print("")
# PDBf.q2()                  # 2. Show the number of amino acids in the chain and how many are missing structural information.

# print("")
# PDBf.q3(2)                 # 3. Show the number of atoms in a given aa (aa sequence number should be provided) and show how many atoms in the side chain.

# print("")
# PDBf.q4()                  # 4. Show the total number of helices and then list how many helices for each type (alpha, 3- 10,...etc).

# print("")
# PDBf.q5('AB')              # 5. Given a sheet ID, show how many strands form that sheet. The user may provide the input as A or AA. Your method will report if the sheet does not exist

# print("")
# PDBf.q6('42.N')            # 5. Return the coordinates of a particular atom (The user should provide aa number and the name of the atom as follow: 56.CA).

# print("")
# PDBf.q7('X.2.begin')       # 6. Show the name (3-letter code) of the amino acid at the beginning or end of a given strand by its number (the user should provide the input as: A.1.end or E.3.begin).

# print("")
# PDBf.q8('9.C', '99.N')     # 7. Find the distance between given two atoms (the user should provide aa numbers and atom names as follow, comma separated: 55.N, 78.CG2).

# print("")
# PDBf.q9('11.psi')          # 8. Calculate phi or psi for a given amino acid (the user provides the input as 45.phi or 23.psi).
