import math 

class Test:

	def __init__(self):
		pass

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
			return sequence.replace('T','U')

	def translation(self,sequence): #QUESTION 4
		codons_dict = [ 
			("F", ["UUU", "UUC"]), 								#phenylalanie
			("L" , ["UUA", "UUG","CUU", "CUC", "CUA", "CUG"]), 	#leucine
			("I ", ["AUU", "AUC", "AUA"]),						#isoleucine
			("M" , ["AUG"]),									#Methionine, START
			("V" , ["GUU", "GUC", "GUA", "GUG"]), 				#valine
			("S" , ["UCU", "UCC", "UCA", "UCG", "AGU", "AGC"]),	#serina
			("P" , ["CCU", "CCC", "CCA", "CCG"]),				#proline
			("T" , ["ACU", "ACC", "ACA", "ACG"]),				#Threonine
			("A" , ["GCU", "GCC", "GCA", "GCG"]),				#Alanine
			("Y" , ["UAU", "UAC"]),								#Tyrosine
			("STOP" , ["UAA", "UAG", "UGA"]),					#STOP
			("H" , ["CAU", "CAC"]),								#histidine
			("Q" , ["CAA", "CAG"]),								#glutamine
			("N" , ["AAU", "AAC"]),								#asparagine
			("K" , ["AAA", "AAG"]),								#lysine
			("D" , ["GAU", "GAC"]),								#aspartic acid
			("E" , ["GAA", "GAG"]),								#glutamic acid
			("C" , ["UGU", "UGC"]),								#cysteine
			("W" , ["UGG"]),									#tryptophan
			("R" , ["CGU", "CGC", "CGA", "CGG", "AGA", "AGG"]),	#arginine
			("G" , ["GGU", "GGC", "GGA", "GGG"]),				#glycine
		]
		tsequence = self.transcription(sequence)
		if tsequence != False and tsequence != "":
			protein = ""
			n = len(tsequence)
			n = math.floor(n/3)
			for y in range(n):
				for j in codons_dict:
					if (tsequence[(y*3):(y*3+3)]  == "AUG"):
						protein+='M'
						break
					elif('M' in protein):
						if j[0] == "STOP":
							if ( tsequence[(y*3):(y*3+3)] in j[1]):
								return protein
						elif ( tsequence[(y*3):(y*3+3)] in j[1]):
							protein+=j[0]
							break
					else:
						continue
		else:
			print("No amino acids")
			return False
		return protein

	def frames(self,sequence):  #QUESTION 5
		rsequence = self.reverse(sequence)
		nums = ["1st", "2nd", "3rd", "4th", "5th", "6th"]
		for i in range(6):
			if i < 3:
				if i == 0:
					print("REGULAR SEQUENCE")
				tsequence = self.translation(sequence)
				if(tsequence):
					print( nums[i], " frame:", tsequence)
				else:
					print(nums[i], " frame: No amino acids found")
				sequence  = sequence[1:]

			else:
				if i == 3:
					print("\n")
					print("REVERSE SEQUENCE")
				tsequence = self.translation(rsequence)
				if(tsequence):
					print( nums[i], " frame:", tsequence)
				else:
					print(nums[i], " frame: No amino acids found")
				rsequence  = rsequence[1:]






test = Test()
# print(test.validate("GCAGTCA"))			            				   #QUESTION 1
# print(test.reverse("GCAGTCA"))                       					   #QUESTION 2
# print(test.transcription("GCAGTCA"))                  				   #QUESTION 3
# print(test.translation("AATGGCGCCGATATTATGACGGTCCTTCCTTGATGATAAGGTAA"))  #QUESTION 4
# test.frames("AATGGCGCCGATATTATGACGGTCCTTCCTTGATGATAAGGTAA")  			   #QUESTION 5