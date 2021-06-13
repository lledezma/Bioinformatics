import math
import requests
import os
import concurrent.futures
import json
from numpy import *


class Test:

	def __init__(self):
		self.storeAllAngles = {}
		self.storeProteinAngles = []
		self.storeAlldistances = {}
		self.errors = []

	def fetchall(self,line):
		if line[:4].strip() == 'IDs':
			return
		else:
			cProtein = line[:9].strip()			#name of protein with chain included
			protein = cProtein[:4]				#name of protein without chain
			file = "proteins/"+protein+".txt"
			#try to read a protein from a local file
			if (os.path.isfile(file)):
				pdb = open(file, 'r') #opens pdb file 
				pdb.close()
				return
			#if the protein file doesnt exist, we fetch it from online
			else:
				url = 'https://files.rcsb.org/view/' +protein + '.pdb'
				myfile = requests.get(str(url).rstrip())		#copy the information of the url webpage
				open(file, 'a').write(str(myfile.content).replace("\\n","\n"))  #dump the information into <protein name>.txt
			 # self.chain = pdb[4:9].strip()           #read the last character of the input which is the chain

	def get_q3(self,line):
		if line[:4].strip() == 'IDs':
			pass
		else:
			cProtein = line[:9].strip()			#name of protein with chain included
			protein = cProtein[:4]				#name of protein without chain
			chain = cProtein[4:9].strip()		#name of chain
			if(os.path.isfile('proteins/'+protein+'.txt')):
				pdb = open('proteins/'+protein+'.txt', 'r')
				Plines = pdb.readlines()
				pdb.close()
			else:
				return
			try:
				for PPline in Plines:
					allcords = []		#Store all coordinates of CA values
					allcenters = []		#Store all coordinates of Center values 
					if(PPline[:6].strip() == "HELIX") and (PPline[19:20].strip() == chain):
						try:
							for x in range(int(PPline[21:26].strip()),int(PPline[33:38].strip())+1):
								#get the coordinates of all CA amino acids in the helix
								allcords.append((self.get_coord(protein,chain,str(x))))	
						except:
							self.errors.append(PPline)
							pass
				#get all the centers in the helix
						allcords = [x for x in allcords if x != []]
				#get the center of mass of three consecutive CA amino acids
						accum= []
						for j in range(0,len(allcords)-2):
							accum.append(allcords[j])
							accum.append(allcords[j+1])
							accum.append(allcords[j+2])
							allcenters.append(self.get_center(accum))
							accum = []
				#get torsion angle formed by (CAx, CENTERx, CENTERx+1, CAx+3)
						for z in range(0,len(allcords)-3):
							self.get_torsion(allcords[z],allcenters[z],allcenters[z+1],allcords[z+3])
					elif(PPline[:6].strip() == "SHEET"):
						break
			except:
				pass
			if(os.path.isfile('proteins/'+protein+'.txt') and cProtein not in self.storeAllAngles ):
				self.storeAllAngles[cProtein] = self.storeProteinAngles
				return self.storeAllAngles
			else:
				pass
			self.storeProteinAngles = []

	def get_q4(self,line):
		if line[:4].strip() == 'IDs':
			return
		else:
			cProtein = line[:9].strip()			#name of protein with chain included
			protein = cProtein[:4]				#name of protein without chain
			chain = cProtein[4:9].strip()		#name of chain
			#check to see if the protein exists within our folder of proteins
		if(os.path.isfile('proteins/'+protein+'.txt')):
			pdb = open('proteins/'+protein+'.txt', 'r')
			Plines = pdb.readlines()
			pdb.close()
		else:
			return
		allhelices = []
		allTypes = []
		allinfo = {}
		proteinHelices = []
		tick_label = ['Right-handed alpha', 'Right-handed omega', 'Right-handed pi', 'Right-handed gamma', 'Right-handed 3 - 10', 
			'Left-handed alpha', 'Left-handed omega','Left-handed gamma','2 - 7 ribbon/helix','Polyproline'] 
		for PPline in Plines:
			try:
			#append all helix lines into allhelices array
				if(PPline[:6].strip() == "HELIX" and PPline[18:20].strip() == chain):
					allhelices.append(PPline)
					allTypes.append(PPline[38:40])
			#if "LINK", it means no more helices to append, so we start analysing each helix
				elif(PPline[:4].strip() == "LINK"):
					for i,j in zip(allhelices,allTypes):
						# print(i)
						try:
							lengt=(int(i[33:37].strip())+ 1 - int(i[21:25].strip()))
							proteinHelices.append({"length": lengt, "class": tick_label[int(j)-1] })
						except:
							pass
					self.storeAllAngles[cProtein] = proteinHelices
					return self.storeAllAngles
			except:
				pass

		self.storeAllAngles[cProtein] = proteinHelices
		return self.storeAllAngles

	def get_q5(self,line):
		if line[:4].strip() == 'IDs':
			return
		else:
			cProtein = line[:9].strip()			#name of protein with chain included
			protein = cProtein[:4]				#name of protein without chain
			chain = cProtein[4:9].strip()		#name of chain
		#check to see if the protein exists within our folder of proteins
			if(os.path.isfile('proteins/'+protein+'.txt')):
				pdb = open('proteins/'+protein+'.txt', 'r')
				Plines = pdb.readlines()
				pdb.close()
			else:
				return
			allhelices = []
			allcords = []
			alldistances = [] 
			for PPline in Plines:
				try:
				#append all helix lines into allhelices array
					if(PPline[:6].strip() == "HELIX" and PPline[18:20].strip() == chain):
						allhelices.append(PPline)
				#if "LINK", it means no more helices to append, so we start analysing each helix
					elif(PPline[:4].strip() == "LINK"):
				#get the all of the CA coordinates and distances for each helix (i)
						for i in allhelices:
							try:
								for x in range(int(i[21:25].strip()), int(i[33:37].strip())+1):
									allcords.append((self.get_coord(protein,chain,str(x))))
							#remove any missing amino acids	
								allcords = [x for x in allcords if x != []]
							#number of CA atom coordinates in the helix
								allcordsL = len(allcords)
							#if the number of CA atom coordinates is < 5 we pass, else we get the center line
								if allcordsL < 5:
									alldistances.append("less than five amino acids")
									pass
								else:
									centerLine = self.get_line(allcords[0],allcords[1],allcords[2],allcords[3], allcords[allcordsL-4],
									allcords[allcordsL-3],allcords[allcordsL-2],allcords[allcordsL-1])
							#getting all of the rise distances for the helix
								for j in range(0,allcordsL-1):
									alldistances.append(self.get_distance(allcords[j],allcords[j+1], centerLine))
							#clear the array to fill it with the coordinates of the next helix
								allcords = []
							except:
							#if there was an error trying to read the values from a helix, we append the line to errors array
								self.errors.append(PPline)
								pass
						self.storeAlldistances[cProtein] = alldistances
						return self.storeAlldistances
				except:
					pass
		self.storeAlldistances[cProtein] = alldistances
		return self.storeAlldistances

	def get_q6_3(self,line):
		counter = 0
		if line[:4].strip() == 'IDs':
			pass
		else:
			cProtein = line[:9].strip()			#name of protein with chain included
			protein = cProtein[:4]				#name of protein without chain
			chain = cProtein[4:9].strip()		#name of chain
			if(os.path.isfile('proteins/'+protein+'.txt')):
				pdb = open('proteins/'+protein+'.txt', 'r')
				Plines = pdb.readlines()
				pdb.close()
			else:
				return
			try:
				allsheets = [] 
				allcords = []
				allcenters = []
				for PPline in Plines:
					if(PPline[:6].strip() == "SHEET"):
						if allsheets == [] and chain == PPline[21:22].strip():
							counter+=1
							if counter >1:
								allsheets.append(missedline)
							else:
								allsheets.append(PPline)
							pass
						if(PPline[10:14].strip() == allsheets[0][10:14].strip() and PPline != allsheets[0] and chain == PPline[21:22].strip()):
							allsheets.append(PPline)
							pass
						if(PPline[10:14].strip() != allsheets[0][10:14].strip()):
							try:
								for x in range(int(allsheets[0][23:27].strip()), int(allsheets[len(allsheets)-1][33:37].strip())+1):
					#get the coordinates of all CA amino acids in the helix
									allcords.append((self.get_coord(protein,chain,str(x))))
							except:
								self.errors.append(PPline)
								pass
						#remove missing cordinates
							allcords = [x for x in allcords if x != []]
						#The number of CA coordinates
							allcordsL = len(allcords) 
							accum= []
							for j in range(0,len(allcords)-2):
								accum.append(allcords[j])
								accum.append(allcords[j+1])
								accum.append(allcords[j+2])
								allcenters.append(self.get_center(accum))
								accum = []
					 #get torsion angle formed by (CAx, CENTERx, CENTERx+1, CAx+3)
							for z in range(0,len(allcords)-3):
								self.get_torsion(allcords[z],allcenters[z],allcenters[z+1],allcords[z+3])
							allsheets = []
							allcords = []
							allcenters = []
							missedline = PPline
							if(chain != PPline[21:22].strip()):
								self.storeAllAngles[cProtein] = self.storeProteinAngles
								return self.storeAllAngles
			except:
				pass

	def get_q6_4(self,line):
		counter = 0
		if line[:4].strip() == 'IDs':
			pass
		else:
			cProtein = line[:9].strip()			#name of protein with chain included
			protein = cProtein[:4]				#name of protein without chain
			chain = cProtein[4:9].strip()		#name of chain
			if(os.path.isfile('proteins/'+protein+'.txt')):
				pdb = open('proteins/'+protein+'.txt', 'r')
				Plines = pdb.readlines()
				pdb.close()
			else:
				return
			try:
				allsheets = [] 
				allcords = []
				allcenters = []	
				sheetInfo = {}	
				for PPline in Plines:
					if(PPline[:6].strip() == "SHEET" ):
						if allsheets == [] and chain == PPline[21:22].strip():
							counter+=1
							if counter >1:
								allsheets.append(missedline)
							else:
								allsheets.append(PPline)
							pass
						if(PPline[10:14].strip() == allsheets[0][10:14].strip() and PPline != allsheets[0] and chain == PPline[21:22].strip()):
							allsheets.append(PPline)
							pass
						if(PPline[10:14].strip() != allsheets[0][10:14].strip() ):
							sheetlength= 0
							sheetype = 0
							for h in allsheets:
								sheetype+= int(h[38:40].strip())
								sheetlength += (int(h[33:37].strip())) - (int(h[22:27].strip())) +1
							if sheetype >= 0:
								sheetype = 1
							else:
								sheetype = -1
							sheetInfo['Type'] = sheetype
							sheetInfo['Length'] = sheetlength
							allcords.append(sheetInfo)
							allsheets = []
							sheetInfo = {}
							missedline = PPline
						if(chain != PPline[21:22].strip()):
							self.storeAllAngles[cProtein] = allcords
							return self.storeAllAngles
			except:
				pass

	def get_q6_5(self,line):
		if line[:4].strip() == 'IDs':
			return
		else:
			cProtein = line[:9].strip()			#name of protein with chain included
			protein = cProtein[:4]				#name of protein without chain
			chain = cProtein[4:9].strip()		#name of chain
			if(os.path.isfile('proteins/'+protein+'.txt')):
				pdb = open('proteins/'+protein+'.txt', 'r')
				Plines = pdb.readlines()
				pdb.close()
			else:
				return
			if(cProtein not in self.trackdistances):
				alldistances = []
				allsheets = [] 
				allcords = []		#Store all coordinates of CA values
				for PPline in Plines:
					try:
						if(PPline[:6].strip() == "SHEET" and chain == PPline[21:22].strip()):
							if allsheets == []:
								allsheets.append(PPline)
								pass
							if(PPline[10:14].strip() == allsheets[0][10:14].strip() and PPline != allsheets[0]):
								allsheets.append(PPline)
								pass
							if(PPline[10:14].strip() != allsheets[0][10:14].strip()):
								try:
									for x in range(int(allsheets[0][23:27].strip()), int(allsheets[len(allsheets)-1][33:37].strip())+1):
						#get the coordinates of all CA amino acids in the helix
										allcords.append((self.get_coord(protein,chain,str(x))))
								except:
									self.errors.append(PPline)
									pass
						#remove missing cordinates
								allcords = [x for x in allcords if x != []]
								allcordsL = len(allcords) 
						#if we have more than 5 amino acid coordinates we make a center line
								if(allcordsL >= 5):
									centerLine = self.get_line(allcords[0],allcords[1],allcords[2],allcords[3], allcords[allcordsL-4],
										allcords[allcordsL-3],allcords[allcordsL-2],allcords[allcordsL-1])
								else:
									alldistances.append("less than 5 Amino Acids")
									pass
						#get the distances
								for j in range(0,len(allcords)-1):
									alldistances.append(self.get_distance(allcords[j],allcords[j+1], centerLine))
								allsheets = []
								allcords = []
						if(PPline[:4].strip() == "LINK"):
							break
					except:
						pass
				self.storeAlldistances[cProtein] = alldistances
			else:
				return
			return self.storeAlldistances

	def get_coord(self,protein,chain,aaNumber):
		pdb = open('proteins/'+protein+ '.txt', 'r')
		lines = pdb.readlines()
		pdb.close()
		CA = []
		for line in lines:
			if(  line[:4] == 'ATOM'  and line[21:22] == chain and line[12:16].strip() == 'CA' and line[22:26].strip() == aaNumber):
				CA.append(float(line[30:38].strip())) #append x
				CA.append(float(line[38:46].strip())) #append y
				CA.append(float(line[46:54].strip())) #append z
				break
		return CA

	def get_center(self,atoms):
		center = []
		center.append((atoms[0][0] + atoms[1][0] + atoms[2][0]) / 3)
		center.append((atoms[0][1] + atoms[1][1] + atoms[2][1]) / 3)
		center.append((atoms[0][2] + atoms[1][2] + atoms[2][2]) / 3)
		center[0] = round(center[0],3)
		center[1] = round(center[1],3)
		center[2] = round(center[2],3)
		return center

	def get_distance(self, atom1, atom2,centerLine):
		fp = (self.get_projection(atom1,centerLine[0],centerLine[1]))
		sp = (self.get_projection(atom2,centerLine[0],centerLine[1]))
		try:
			distance = math.sqrt(   ((sp[0]-fp[0])**2) + ((sp[1]-fp[1])**2 ) + ((sp[2]-fp[2])**2 )  )
		except:
		    distance = 0
		return round(distance,3)

	def get_line(self,atom1,atom2, atom3, atom4,atom5,atom6,atom7,atom8):
		maincenter = []
		center1 = []
		center2 = [] 
		#get center 1
		center1.append( (atom1[0] + atom2[0] + atom3[0]  + atom4[0] )  / 4)
		center1.append( (atom1[1] + atom2[1] + atom3[1]  + atom4[1]  )  / 4)
		center1.append( (atom1[2] + atom2[2] + atom3[2]  + atom4[2] )   / 4)
		center1[0] = round(center1[0],3)
		center1[1] = round(center1[1],3)
		center1[2] = round(center1[2],3)
		#get center 2
		center2.append(  (atom5[0] + atom6[0] + atom7[0]  + atom8[0]) / 4)
		center2.append(  (atom5[1] + atom6[1] + atom7[1]  + atom8[1]) / 4)
		center2.append(  (atom5[2] + atom6[2] + atom7[2]  + atom8[2]) / 4)
		center2[0] = round(center2[0],3)
		center2[1] = round(center2[1],3)
		center2[2] = round(center2[2],3)
		maincenter.append(center1)
		maincenter.append(center2)
		return maincenter

	def store_tor_angles(self, TorsionAngel):
		self.storeProteinAngles.append(TorsionAngel)
		return

	def get_torsion(self,Calpha1,center1,center2,Calpha2):
			RadandDeg = {'Radian':0.0, 'Degree':0.0}
			angletype = 'phi'
			if(angletype == 'phi'):
				Cphi = Calpha1
				N = center1
				CA = center2
				C = Calpha2
			try:
				if angletype == 'phi':
					A1 =  Cphi[1]*(N[2]-CA[2])   +      N[1]*( CA[2]-Cphi[2] )    +   CA[1]*(Cphi[2]- N[2])
					B1 = Cphi[2]*(N[0]-CA[0])    +      N[2]*( CA[0]-Cphi[0] )    +   CA[2]*(Cphi[0]- N[0])
					C1 = Cphi[0]*(N[1]-CA[1])    +      N[0]*( CA[1]-Cphi[1] )    +   CA[0]*(Cphi[1]- N[1])
					A2 =  ( N[1]*(CA[2]-C[2]) )  +      (CA[1]*( C[2]-N[2] ))   +  (C[1]*(N[2]- CA[2]))
					B2 = (N[2]*(CA[0]-C[0]))     +      (CA[2]*( C[0]-N[0] ))   +  (C[2]*(N[0]- CA[0]))
					C2 = (N[0]*(CA[1]-C[1]))     +      (CA[0]*( C[1]-N[1] ))   +  (C[0]*(N[1]- CA[1]))
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
						finalangle = -1*finalangle
					radian = math.acos(angle)
					RadandDeg['Radian'] = radian
					RadandDeg['Degree'] = finalangle
			except:
				pass
			self.store_torAngles(RadandDeg)
			return  RadandDeg

	def get_projection(self,atom, center1, center2):
		a = center1 
		b = center2
		p= atom
		ap = []
		ap.append(p[0] - a[0])
		ap.append(p[1] - a[1])
		ap.append(p[2] - a[2])
		ab = []
		ab.append(b[0] - a[0])
		ab.append(b[1] - a[1])
		ab.append(b[2] - a[2])
		last = []
		last.append(dot(ab,ab)*ab[0])
		last.append(dot(ab,ab)*ab[1])
		last.append(dot(ab,ab)*ab[2])
		Nlast = []
		if(last[0] != 0 ):
			Nlast.append(dot(ap,ab)/last[0])
		else:
			Nlast.append(0)
		if(last[1] != 0 ):
			Nlast.append(dot(ap,ab)/last[1])
		else:
			Nlast.append(0)
		if(last[2] != 0 ):
			Nlast.append(dot(ap,ab)/last[2])
		else:
			Nlast.append(0)
		newlast = []
		newlast.append(a[0]+Nlast[0])
		newlast.append(a[1]+Nlast[1])
		newlast.append(a[2]+Nlast[2])
		return newlast


Fetch = Test()

proteins = open('proteins.txt', 'r')
lines = proteins.readlines()
proteins.close()

if __name__ == '__main__':
	with concurrent.futures.ProcessPoolExecutor() as executor:
		result = executor.map(Fetch.get_q4,lines)
		for i in result:
			with open('Q4.json') as fh:
				out_file = open("Q4.json", "a+") 
			json.dump(i, out_file, indent = 4) 
		    out_file.close()

