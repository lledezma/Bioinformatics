import math
import requests
import os
import concurrent.futures
import json
from numpy import *

class test:

	def __init__(self):
		self.storeAllAngles = {}
		self.storeProteinAngles = []
		self.errors = []
		self.alltorsions = []
		self.secondary = []


	def get_helices(self,line):
		if line[:4].strip() == 'IDs':
			return
		else:
			cProtein = line[:9].strip()			#name of protein with chain included
			protein = cProtein[:4].strip()				#name of protein without chain
			chain = cProtein[4:9].strip()		#name of chain
		try:
			if(os.path.isfile('proteins/'+protein+'.txt')):
				pdb = open('proteins/'+protein+'.txt', 'r')
				Plines = pdb.readlines()
				pdb.close()
			else:
				pass
		except:
			cProtein = line[:9].strip()			#name of protein with chain included
			protein = cProtein[:4].strip()				#name of protein without chain
			chain = cProtein[4:9].strip()		#name of chain
		if(os.path.isfile('proteins/'+protein+'.txt')):
			pdb = open('proteins/'+protein+'.txt', 'r')
			Plines = pdb.readlines()
			pdb.close()
		else:
			pass
		try:
			lastAA = 0
			firstAA = 0
			for PPline in Plines:
				if(PPline[:6] == 'SEQRES' and chain == PPline[10:12].strip()):
					if(lastAA == 0):
						lastAA = int(PPline[12:17].strip())
					else:
						pass
				if(PPline[:6].strip() == "HELIX") and (PPline[19:20].strip() == chain):
					self.secondary.append(PPline)
				if(PPline[:6].strip() == "SHEET" and chain == PPline[21:22].strip()):
					self.secondary.append(PPline)
				if (PPline[:4] == 'ATOM') and (PPline[21:22] == chain):
					firstAA = int(PPline[22:26].strip())
					for x in range(0,int(lastAA)):
						if (self.get_coord(protein,chain,str(x+int(firstAA))) != None):
							self.alltorsions.append(self.get_coord(protein,chain,str(x+int(firstAA))))
					self.storeAllAngles[cProtein] = self.alltorsions
					self.secondary = []
					return self.storeAllAngles
		except:
			pass
		self.alltorsions = []
		self.secondary = []



	def get_coord(self,protein,chain,aaNumber):
		pdb = open('proteins/'+protein+ '.txt', 'r')
		lines = pdb.readlines()
		pdb.close()
		Cphi = []
		N = []
		CA = []
		C = []
		Npsi = []
		aaName = ""
		for line in lines:
			if(  line[:4] == 'ATOM'  and line[21:22] == chain):
				if(  line[:4] == 'ATOM'  and line[21:22] == chain and line[12:16].strip() == 'C' and int(line[22:26].strip()) == int(aaNumber)-1):
					Cphi.append(float(line[30:38].strip()))
					Cphi.append(float(line[38:46].strip()))
					Cphi.append(float(line[46:54].strip()))
				if(  line[:4] == 'ATOM'  and line[21:22] == chain and line[12:16].strip() == 'N' and line[22:26].strip() == aaNumber):
				    N.append(float(line[30:38].strip())) #append x
				    N.append(float(line[38:46].strip())) #append y
				    N.append(float(line[46:54].strip())) #append z
				if(  line[:4] == 'ATOM'  and line[21:22] == chain and line[12:16].strip() == 'CA' and line[22:26].strip() == aaNumber):
				    CA.append(float(line[30:38].strip())) #append x
				    CA.append(float(line[38:46].strip())) #append y
				    CA.append(float(line[46:54].strip())) #append z
				if(  line[:4] == 'ATOM'  and line[21:22] == chain and line[12:16].strip() == 'C' and line[22:26].strip() == aaNumber):
				    C.append(float(line[30:38].strip())) #append x
				    C.append(float(line[38:46].strip())) #append y
				    C.append(float(line[46:54].strip())) #append z
				    aaName = line[17:20].strip()
				if(  line[:4] == 'ATOM'  and line[21:22] == chain and line[12:16].strip() == 'N' and int(line[22:26].strip()) == int(aaNumber)+1):
				    Npsi.append(float(line[30:38].strip()))
				    Npsi.append(float(line[38:46].strip()))
				    Npsi.append(float(line[46:54].strip()))
				    return self.get_torsion(Cphi,N,CA,C,Npsi, protein, aaNumber,chain,aaName)
			else:
				pass
		return self.get_torsion(Cphi,N,CA,C,Npsi, protein, aaNumber,chain,aaName)



	def get_torsion(self,Cphi,N,CA,C,Npsi, protein, aaNumber,chain,aaName):
			torsions = {}
			torsions['phi'] = 0
			torsions['psi'] = 0
			torsions['protein'] = protein
			torsions['AANumber'] = aaNumber
			torsions['AAName'] = aaName
			torsions['chain'] = chain
			torsions['structure'] = ""
			for i in self.secondary:
				if i[:5] == "HELIX":
					if int(i[21:26].strip()) <= int(aaNumber) <= int(i[33:38].strip()):
						torsions['structure'] = "HELIX"
				else:
					if int(i[23:27].strip()) <= int(aaNumber) <= int(i[33:37].strip()):
						torsions['structure'] = "SHEET"
			for angletype in torsions:
				try:
					if angletype == 'phi':
						A1 =  Cphi[1]*(N[2]-CA[2])   +      N[1]*( CA[2]-Cphi[2] )    +   CA[1]*(Cphi[2]- N[2])
						B1 = Cphi[2]*(N[0]-CA[0])   +      N[2]*( CA[0]-Cphi[0] )    +   CA[2]*(Cphi[0]- N[0])
						C1 = Cphi[0]*(N[1]-CA[1])   +      N[0]*( CA[1]-Cphi[1] )    +   CA[0]*(Cphi[1]- N[1])
						A2 =  ( N[1]*(CA[2]-C[2]) )    +      (CA[1]*( C[2]-N[2] )  )      +   (C[1]*(N[2]- CA[2]))
						B2 = (N[2]*(CA[0]-C[0]))   +      (CA[2]*( C[0]-N[0] ) )   +   (C[2]*(N[0]- CA[0]))
						C2 = (N[0]*(CA[1]-C[1]))   +      (CA[0]*( C[1]-N[1] ) )   +  ( C[0]*(N[1]- CA[1]))
						angle = (A1*A2 + B1*B2 + C1*C2)/ (math.sqrt(A1*A1 + B1*B1 + C1*C1) * math.sqrt(A2*A2 + B2*B2 + C2*C2))
						V1 = [Cphi[0]-N[0], Cphi[1]-N[1], Cphi[2]-N[2]]
						V2 = [CA[0]-N[0], CA[1]-N[1], CA[2]-N[2]]
						V3 = [C[0]-CA[0], C[1]-CA[1], C[2]-CA[2]]
						vNormal = []
						vNormal.append( (V1[1]*V2[2] - V1[2]*V2[1]) )
						vNormal.append( -(V1[0]*V2[2] - V1[2]*V2[0]) )
						vNormal.append( (V1[0]*V2[1] - V1[1]*V2[0]) )
						dotProduct = (V3[0]*vNormal[0] + V3[1]*vNormal[1] + V3[2]*vNormal[2])
						finalangle = round(((math.acos(angle)*180) / math.pi),4)
						if(dotProduct > 0):
							torsions['phi'] = (-1*finalangle)
						else:
							torsions['phi'] = (finalangle)
					if angletype == 'psi':
						A1 =  N[1]*(CA[2]-C[2])   +      CA[1]*( C[2]-N[2] )    +   C[1]*(N[2]- CA[2])
						B1 = N[2]*(CA[0]-C[0])   +      CA[2]*( C[0]-N[0] )    +   C[2]*(N[0]- CA[0])
						C1 = N[0]*(CA[1]-C[1])   +      CA[0]*( C[1]-N[1] )    +   C[0]*(N[1]- CA[1])
						A2 =  ( CA[1]*(C[2]-Npsi[2]) )    +      (C[1]*( Npsi[2]-CA[2] )  )      +   (Npsi[1]*(CA[2]- C[2]))
						B2 = (CA[2]*(C[0]-Npsi[0]))   +      (C[2]*( Npsi[0]-CA[0] ) )   +   (Npsi[2]*(CA[0]- C[0]))
						C2 = (CA[0]*(C[1]-Npsi[1]))   +      (C[0]*( Npsi[1]-CA[1] ) )   +  ( Npsi[0]*(CA[1]- C[1]))
						angle = (A1*A2 + B1*B2 + C1*C2)/ (math.sqrt(A1*A1 + B1*B1 + C1*C1) * math.sqrt(A2*A2 + B2*B2 + C2*C2))
						V1 = [CA[0]-N[0], CA[1]-N[1], CA[2]-N[2]]
						V2 = [C[0]-CA[0], C[1]-CA[1], C[2]-CA[2]]
						V3 = [Npsi[0]-C[0], Npsi[1]-C[1], Npsi[2]-C[2]]
						vNormal = []
						vNormal.append( (V1[1]*V2[2] - V1[2]*V2[1]) )
						vNormal.append( -(V1[0]*V2[2] - V1[2]*V2[0]) )
						vNormal.append( (V1[0]*V2[1] - V1[1]*V2[0]) )
						dotProduct = (V3[0]*vNormal[0] + V3[1]*vNormal[1] + V3[2]*vNormal[2])
						finalangle = round(((math.acos(angle)*180) / math.pi),3)
						if(dotProduct < 0):
							torsions['psi'] = -1*finalangle
						else:
							torsions['psi'] = finalangle
				except:
					pass
			return torsions


		
fetch = test()

proteins = open('proteins.txt', 'r')
lines = proteins.readlines()
proteins.close()

with concurrent.futures.ProcessPoolExecutor() as executor:
	result = executor.map(fetch.get_helices,lines)
		with open('Q2.json') as fh:
			out_file = open("Q2.json", "a+")
		json.dump(i, out_file, indent = 4)
		out_file.close()