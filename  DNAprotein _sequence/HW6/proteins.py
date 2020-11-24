import math 

class HW4:

	def __init__(self):
		print("")

	def validate(self,sequence):	#QUESTION 1
		letters = ['A', 'G', 'C', 'T'];
		for i in sequence:
			if i not in letters:
				return False
		return True

	def reverse(self,sequence):		#QUESTION 2
		if self.validate(sequence) == False:
			print("sequence is not valid")
			return False
		else:
			n = len(sequence)
			s = list(sequence)
			for i in range(0,n):
				if (s[i] == 'A'):
					s[i] = 'T'
				elif (s[i] == 'T'):
					s[i] = 'A'
				elif (s[i] == 'C'):
					s[i] = 'G'
				elif (s[i] == 'G'):
					s[i] = 'C'
		sequence = "".join(s)
		sequence = sequence[::-1]
		return sequence

	def transcription(self,sequence):	#QUESTION 3
		if self.validate(sequence) == False:
			print("sequence is not valid")
			return False
		else:
			n = len(sequence)
			s = list(sequence)
			for i in range(0,n):
				if (s[i] == 'T'):
					s[i] = 'U'	
		sequence = "".join(s)
		return sequence

	def translation(self,sequence): #QUESTION 4
		F = ["UUU", "UUC"] 								#phenylalanie
		L = ["UUA", "UUG","CUU", "CUC", "CUA", "CUG"] 	#leucine
		I = ["AUU", "AUC", "AUA"] 						#isoleucine
		M = ["AUG"]										#Methionine, START
		V = ["GUU", "GUC", "GUA", "GUG"] 				#valine
		S = ["UCU", "UCC", "UCA", "UCG", "AGU", "AGC"]	#serina
		P = ["CCU", "CCC", "CCA", "CCG"]				#proline
		T = ["ACU", "ACC", "ACA", "ACG"]				#Threonine
		A = ["GCU", "GCC", "GCA", "GCG"]				#Alanine
		Y = ["UAU", "UAC"]								#Tyrosine
		STOP = ["UAA", "UAG", "UGA"]								#STOP
		H = ["CAU", "CAC"]								#histidine
		Q = ["CAA", "CAG"]								#glutamine
		N = ["AAU", "AAC"]								#asparagine
		K = ["AAA", "AAG"]								#lysine
		D = ["GAU", "GAC"]								#aspartic acid
		E = ["GAA", "GAG"]								#glutamic acid
		C = ["UGU", "UGC"]								#cysteine
		W = ["UGG"]										#tryptophan
		R = ["CGU", "CGC", "CGA", "CGG", "AGA", "AGG"]	#arginine
		G = ["GGU", "GGC", "GGA", "GGG"]				#glycine
		codons = [F,L,I,M,V,S,P,T,A,Y,STOP,H,Q,N,K,D,E,C,W,R,G]
		codonsl = ['F','L','I', 'M','V','S','P','T','A','Y','STOP','H','Q','N','K','D','E','C','W','R','G']
		if(self.transcription(sequence) == False):
			print("No amino acids")
			return False
		else:
			tsequence = self.transcription(sequence)
		protein = []
		n = len(tsequence)
		n = math.floor(n/3)
		if tsequence != "":
			for y in range(n):
				for j in codons:
					if (tsequence[(y*3):(y*3+3)]  == "AUG"):
						if('M' not in protein):
							protein.append(codonsl[codons.index(M)])
							break
						else:
							protein.append(codonsl[codons.index(M)])
							break   
					elif('M' in protein):
						if j == codons[10]:
							if ( tsequence[(y*3):(y*3+3)] in j):
								protein = "".join(protein)
								return protein
						elif ( tsequence[(y*3):(y*3+3)] in j):
							protein.append(codonsl[codons.index(j)])
							break
					else:
						pass
		else:
			print("No amino acids")
			return False
		protein = "".join(protein)
		return protein



	def frames(self,sequence):  #QUESTION 5
		rsequence = sequence

		print("REGULAR SEQUENCE")
		if(self.translation(sequence) == ""):
			print("1st frame:", "No amino acids found")
		else:
			print("1st frame:", self.translation(sequence))
		sequence  = sequence[1:]
		if(self.translation(sequence) == ""):
			print("2nd frame:", "No amino acids found")
		else:
			print("2nd frame:", self.translation(sequence))
		sequence  = sequence[1:]
		if(self.translation(sequence) == ""):
			print("3rd frame:", "No amino acids found")
		else:
			print("3rd frame:", self.translation(sequence))
		print("\n")
		rsequence = self.reverse(rsequence)
		print("REVERSE SEQUENCE")
		if(self.translation(rsequence) == ""):
			print("4th frame:", "No amino acids found")
		else:
			print("4th frame:", self.translation(rsequence))
		rsequence  = rsequence[1:]
		if(self.translation(rsequence) == ""):
			print("5th frame:", "No amino acids found")
		else:
			print("5th frame:", self.translation(rsequence))
		rsequence  = rsequence[1:]
		if(self.translation(rsequence) == ""):
			print("6th frame:", "No amino acids found")
		else:
			print("6th frame:", self.translation(rsequence))




test = HW4()
# print(test.validate("GCAGTCA"))			            				   #QUESTION 1
# print(test.reverse("GCAGTCA"))                       					   #QUESTION 2
# print(test.transcription("GCAGTCA"))                  				   #QUESTION 3
# print(test.translation("AATGGCGCCGATATTATGACGGTCCTTCCTTGATGATAAGGTAA"))                #QUESTION 4
# test.frames("AATGGCGCCGATATTATGACGGTCCTTCCTTGATGATAAGGTAA")  			   #QUESTION 5