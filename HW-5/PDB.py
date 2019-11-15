import math
class PDB:

    def __init__(self, pdb):
        self.protein = pdb[:4]
        self.file = self.protein+'.txt'
        self.pdb = open(self.file, 'r') #opens pdb file 3UTS.text
        self.pdb.close()
        self.chain = pdb[4:5]           #read the last character of the input which is the chain

    #Show the number of atoms in a given amino acid (by aa sequence number).
    def NumberOfAtoms(self,aaSequenceNumber):
        countAtoms = 0
        for line in lines:
            if (line[:4] == 'ATOM') and (line[21:22] == self.chain):
                #from 23-26
                if(line[23:26].strip() == str(aaSequenceNumber)):
                    countAtoms += 1
                else:
                    pass
        print("Number of atoms:", countAtoms)

    # Show the number of atoms in the side chain for a given aa (aa sequence number should be provided).
    def SideChainAtoms(self, aaSequenceNumber):
        countAtoms = 0
        atomCB = 0
        for line in lines:
                #ATOM                           #chain                          #sequence number
            if (line[:4] == 'ATOM') and (line[21:22] == self.chain) and (line[23:26].strip() == str(aaSequenceNumber)):
                    # atom name
                if(line[13:16].strip() == 'CB'):
                    atomCB = int(line[7:11].strip())
                    countAtoms +=1
                if(atomCB > 0):
                    if (int(line[7:11].strip()) > atomCB ):
                        countAtoms+=1
        print("number of atoms in the side chain:", countAtoms)




    # Show the total number of amino acids.
    def TotalNumberOfAA(self):
        aminoAcids = []
        for line in lines:
            if (line[:4] == 'ATOM') and (line[21:22] == self.chain):
                if(line[17:20].strip() not in aminoAcids):
                    aminoAcids.append(line[17:20].strip())
                else:
                    pass
        print("total number of amino acids:", len(aminoAcids))


    # Show the number of helices and the number of sheets.
    def NumberOfHelicesAndSheets(self):
        helixcount = 0
        sheetcount = 0
        for line in lines:
            if(line[:6].strip() == "HELIX") and (line[19:20].strip() == self.chain):
                helixcount +=1
                # print(line)
            if((line[:6].strip() == "SHEET") and (line[21:22].strip() == self.chain)  and (line[38:40].strip() == '0')):
                sheetcount+=1

        print("Total number of Helices: " , helixcount)
        print("Total number of Sheets: " , sheetcount)


    # Return the coordinates of a particular atom (The user should provide aa number and the name of the atom as
    # follow: 56.CA).
    def getCoordinates(self, atom):
        aanumber, atomname = atom.split('.')
        X = ""
        Y = ""
        Z = ""
        for line in lines:
            if(  line[:4] == 'ATOM'  and line[21:22] == self.chain and line[12:16].strip() == atomname and  line[23:26].strip() == aanumber):
            # print(line[12:16].strip())
                X = line[30:38].strip()
                Y = line[38:46].strip()
                Z = line[46:54].strip()
        print("x: ", X)
        print("y: ", Y)
        print("z: ", Z)


    # Show the name (3-letter code) of the amino acid at the beginning or end of a given helix by its
    # number (the user should provide the input as: 2.end or 4.begin).
    def getCodeName(self, helix):
        number, position = helix.split('.')
        for line in lines:
            if (line[:6].strip() == "HELIX" and line[19:20].strip() == self.chain and line[7:10].strip() == number  and position == 'end'):
                print("name of the amino acid at the end of a given helix:", line[27:30])
                return
            elif (line[:6].strip() == "HELIX" and line[19:20].strip() == self.chain and line[7:10].strip() == number  and position == 'begin'):
                print("name of the amino acid at the beginning of a given helix:",line[15:18])
                return
        print("amino acid doesn't exist in the helix.")


    # Find the distance between given two atoms (the user should provide aa numbers and atom names as follow,
    # comma separated: 55.N, 78.CG2).
    def getDistance(self, atom1, atom2):
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
            print("One of the atoms doesn't exist in the chain.")
            return
        print("Distance:", round(distance,3))

    # Calculate phi or psi for a given amino acid (the user provide the input as 45.phi or 23.psi).
    def get_Phi_or_Psi(self, aa):
        aaNumber, angletype =  aa.split('.')
        Cphi = []
        N = []
        CA = []
        C = []
        Npsi = []
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

        if angletype == 'phi':
            # A1 = y1 ( z2 - z3 )        +      y2 ( z3 - z1 )            +   y3 ( z1 - z2 )
            A1 =  Cphi[1]*(N[2]-CA[2])   +      N[1]*( CA[2]-Cphi[2] )    +   CA[1]*(Cphi[2]- N[2])
            # B1 = z1 ( x2 - x3 )       +      z2 (   x3 - x1 )          +      z3 ( x1 - x2 )
            B1 = Cphi[2]*(N[0]-CA[0])   +      N[2]*( CA[0]-Cphi[0] )    +   CA[2]*(Cphi[0]- N[0])
            # C1= x1 ( y2 - y3 )        +       x2 ( y3 - y1 )           +  x3 ( y1 - y2 )
            C1 = Cphi[0]*(N[1]-CA[1])   +      N[0]*( CA[1]-Cphi[1] )    +   CA[0]*(Cphi[1]- N[1])
    # ////////////////////////////
            # A1 = y2 ( z3 - z4 )        +      y3 ( z4 - z2 )            +   y4 ( z2 - z3 )
            A2 =  ( N[1]*(CA[2]-C[2]) )    +      (CA[1]*( C[2]-N[2] )  )      +   (C[1]*(N[2]- CA[2]))
            # B2 = z2 ( x3 - x4 )       +      z3 (   x4 - x2 )          +      z4 ( x2 - x3 )
            B2 = (N[2]*(CA[0]-C[0]))   +      (CA[2]*( C[0]-N[0] ) )   +   (C[2]*(N[0]- CA[0]))
            # C2 = x2 ( y3 - y4 )        +       x3 ( y4 - y2 )           +  x4 ( y2 - y3 )
            C2 = (N[0]*(CA[1]-C[1]))   +      (CA[0]*( C[1]-N[1] ) )   +  ( C[0]*(N[1]- CA[1]))
        # //////////////////////////
            angle = (A1*A2 + B1*B2 + C1*C2)/ (math.sqrt(A1*A1 + B1*B1 + C1*C1) * math.sqrt(A2*A2 + B2*B2 + C2*C2))
            # //////////////////  vNormal = the cross product of v1 and v2
            vNormal = []
            vNormal.append((Cphi[1]*CA[2] - Cphi[2]*CA[1]))
            vNormal.append(-1*(Cphi[0]*CA[2]-Cphi[2]*CA[0]))
            vNormal.append(-1*(Cphi[0]*CA[1]-Cphi[1]*CA[0]))
            dotProduct = (C[1]*vNormal[2] - C[2]*vNormal[1]) - (C[0]*vNormal[2]-C[2]*vNormal[0]) - (C[0]*vNormal[1]-C[1]*vNormal[0])
            # //////////////// final angle /// (acos(angle)* 180)/PI
            finalangle = round(((math.acos(angle)*180) / math.pi),3)
            if ( dotProduct > 0):
                print("phi angle:",-1*finalangle)
            else:
                print("phi angle:",finalangle)
# ////////////////////////////////////////////////////////////////////////////////////////////////
        if angletype == 'psi':
            # A1 = y1 ( z2 - z3 )        +      y2 ( z3 - z1 )            +   y3 ( z1 - z2 )
            A1 =  N[1]*(CA[2]-C[2])   +      CA[1]*( C[2]-N[2] )    +   C[1]*(N[2]- CA[2])
            # B1 = z1 ( x2 - x3 )       +      z2 (   x3 - x1 )          +      z3 ( x1 - x2 )
            B1 = N[2]*(CA[0]-C[0])   +      CA[2]*( C[0]-N[0] )    +   C[2]*(N[0]- CA[0])
            # C1= x1 ( y2 - y3 )        +       x2 ( y3 - y1 )           +  x3 ( y1 - y2 )
            C1 = N[0]*(CA[1]-C[1])   +      CA[0]*( C[1]-N[1] )    +   C[0]*(N[1]- CA[1])
    # ////////////////////////////
            # A1 = y2 ( z3 - z4 )        +              y3 ( z4 - z2 )            +   y4 ( z2 - z3 )
            A2 =  ( CA[1]*(C[2]-Npsi[2]) )    +      (C[1]*( Npsi[2]-CA[2] )  )      +   (Npsi[1]*(CA[2]- C[2]))
            # B2 = z2 ( x3 - x4 )       +           z3 (   x4 - x2 )          +      z4 ( x2 - x3 )
            B2 = (CA[2]*(C[0]-Npsi[0]))   +      (C[2]*( Npsi[0]-CA[0] ) )   +   (Npsi[2]*(CA[0]- C[0]))
            # C2 = x2 ( y3 - y4 )        +            x3 ( y4 - y2 )           +  x4 ( y2 - y3 )
            C2 = (CA[0]*(C[1]-Npsi[1]))   +      (C[0]*( Npsi[1]-CA[1] ) )   +  ( Npsi[0]*(CA[1]- C[1]))
        # //////////////////////////
            angle = (A1*A2 + B1*B2 + C1*C2)/ (math.sqrt(A1*A1 + B1*B1 + C1*C1) * math.sqrt(A2*A2 + B2*B2 + C2*C2))
            # //////////////////  vNormal = the cross product of v1 and v2
            vNormal = []
            vNormal.append((N[1]*C[2] - N[2]*C[1]))
            vNormal.append(-1*(N[0]*C[2]-N[2]*C[0]))
            vNormal.append(-1*(N[0]*C[1]-N[1]*C[0]))
            dotProduct = (Npsi[1]*vNormal[2] - Npsi[2]*vNormal[1]) - (Npsi[0]*vNormal[2]-Npsi[2]*vNormal[0]) - (Npsi[0]*vNormal[1]-Npsi[1]*vNormal[0])
            # //////////////// final angle /// (acos(angle)* 180)/PI
            finalangle = round(((math.acos(angle)*180) / math.pi),3)
            if ( dotProduct < 0):
                print("psi angle:",-1*finalangle)
            else:
                print("psi angle:",finalangle)

fetch = PDB('3UTSA')
fetch.NumberOfAtoms(3)
print("")
fetch.SideChainAtoms(3)
print("")
fetch.TotalNumberOfAA()
print("")
fetch.NumberOfHelicesAndSheets()
print("")
fetch.getCoordinates('6.N')
print("")
fetch.getCodeName('3.begin')
print("")
fetch.getDistance('1.N', '276.N')
print("")
fetch.get_Phi_or_Psi('45.phi')
