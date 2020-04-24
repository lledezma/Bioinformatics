import json
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import concurrent.futures

class parse:

	def __init__(self):
		print("")

	def Q6_5(self):
		with open('Q6_5.json') as f:
			data = json.load(f)

		labels = ['0-0.5', '0.5-1', '1-1.5', '1.5-2', '2-2.5', '2.5-3', '3-3.5',
				  '3.5-4', '4-4.5', '4.5-5', '>5', 'less-than-5-AA']

		lengs = {}
		for i in labels:
			lengs[i] = 0

		avg = 0.0
		count = 0
		for i in data:
			for j in i:
				if j == []:
					pass
				for h in i[j]:
					try:
						if 0 < h < 0.5:
							lengs['0-0.5'] +=1
						elif 0.5 < h < 1:
							lengs['0.5-1'] +=1
						elif 1 < h < 1.5:
							lengs['1-1.5'] +=1
						elif 1.5 < h < 2:
							lengs['1.5-2'] +=1
						elif 2 < h < 2.5:
							lengs['2-2.5'] +=1
						elif 2.5 < h < 3:
							lengs['2.5-3'] +=1
						elif 3 < h < 3.5:
							lengs['3-3.5'] +=1
						elif 3.5 < h < 4:
							lengs['3.5-4'] +=1
						elif 4 < h < 4.5:
							lengs['4-4.5'] +=1
						elif 4.5 < h < 5:
							lengs['4.5-5'] +=1
						elif h >= 5.0:
							lengs['>5'] +=1
						avg+=h
						count+=1
					except:
						lengs['less-than-5-AA'] +=1
		men_means = []

		for i in lengs:
			men_means.append(lengs[i])

		x = np.arange(len(labels))  # the label locations
		width = 0.35  # the width of the bars

		fig, ax = plt.subplots()
		rects1 = ax.bar(x - width/2, men_means, width, label='Average:'+str(round(avg/count, 3)))

		# Add some text for labels, title and custom x-axis tick labels, etc.
		ax.set_ylabel('# of Distances')
		ax.set_xlabel('Rise Distance')
		ax.set_title('Beta Sheets Data')
		ax.set_xticks(x)
		ax.set_xticklabels(labels)
		ax.legend()


		def autolabel(rects):
		    """Attach a text label above each bar in *rects*, displaying its height."""
		    for rect in rects:
		        height = rect.get_height()
		        ax.annotate('{}'.format(height),
		                    xy=(rect.get_x() + rect.get_width() / 2, height),
		                    xytext=(0, 3),  # 3 points vertical offset
		                    textcoords="offset points",
		                    ha='center', va='bottom')
		autolabel(rects1)
		fig.tight_layout()
		plt.show()
		return

	def Q6_4(self):
		with open('/Q6_4.json') as f:
			data = json.load(f)

		negativeType = 0;
		positiveType = 0;

		avgPosLen = 0;
		avgNegLen = 0;

		positiveValues = []
		negativeValues = []


		labels = []

		for i in data:
			for j in i:
				if j == []:
					pass
				if j[4:5] not in labels:
					labels.append(j[4:5])
				for h in i[j]:
					if h['Type'] == 1:
						positiveType +=1
						avgPosLen += h['Length']
						positiveValues.append(h['Length'])
					else:
						negativeType+=1
						avgNegLen += h['Length']
						negativeValues.append(int(h['Length']))

		avgPosLen = round(avgPosLen/positiveType)
		avgNegLen = round(avgNegLen/negativeType)


		Type = ['Anti-Parallel, Average:'+str(avgNegLen), 'Parallel, Average:'+str(avgPosLen)]
		# portion covered by each label
		slices = [negativeType,positiveType]
		# color for each label
		colors = ['r', 'b']
		# plotting the pie chart
		plt.pie(slices, labels = Type, colors=colors,
		        startangle=90, shadow = False, explode = (0, 0),
		        radius = 1.2, autopct = '%1.1f%%')
		# plotting legend
		plt.legend()
		# showing the plot
		plt.show()



		labels = ['Average', 'Largest', 'Smallest']
		men_means = [avgPosLen, max(positiveValues), min(positiveValues) ]
		women_means = [avgNegLen, max(negativeValues), min(negativeValues)]

		x = np.arange(len(labels))  # the label locations
		width = 0.35  # the width of the bars

		fig, ax = plt.subplots()
		rects1 = ax.bar(x - width/2, men_means, width, label='Parallel')
		rects2 = ax.bar(x + width/2, women_means, width, label='Anti-Parallel')

		# Add some text for labels, title and custom x-axis tick labels, etc.
		ax.set_ylabel('Length')
		ax.set_title('Beta Sheets Data')
		ax.set_xticks(x)
		ax.set_xticklabels(labels)
		ax.legend()


		def autolabel(rects):
		    """Attach a text label above each bar in *rects*, displaying its height."""
		    for rect in rects:
		        height = rect.get_height()
		        ax.annotate('{}'.format(height),
		                    xy=(rect.get_x() + rect.get_width() / 2, height),
		                    xytext=(0, 3),  # 3 points vertical offset
		                    textcoords="offset points",
		                    ha='center', va='bottom')


		autolabel(rects1)
		autolabel(rects2)

		fig.tight_layout()

		plt.show()


		return

	def Q6_3(self):
		with open('Q6_3.json') as f:
			data = json.load(f)

		Degrees = []
		Radians = []
		number = 0
		for i in data:
			for j in i:
				if j == []:
					pass
				for h in i[j]:
					Degrees.append(h['Degree'])
					# print(h['Degree'])
					Radians.append(h['Radian'])
		x = np.array(Radians)
		y = np.array(Degrees)
		plt.scatter(x, y)
		plt.show()

	def Q5(self):
		with open('Q5.json') as f:
			data = json.load(f)
		labels = ['0-0.5', '0.5-1', '1-1.5', '1.5-2', '2-2.5', '2.5-3', '3-3.5',
				  '3.5-4', '4-4.5', '4.5-5', '>5', 'less-than-5-AA']

		lengs = {}
		for i in labels:
			lengs[i] = 0

		avg = 0.0
		count = 0
		for i in data:
			for j in i:
				if j == []:
					pass
				for h in i[j]:
					try:
						if 0 < h < 0.5:
							lengs['0-0.5'] +=1
						elif 0.5 < h < 1:
							lengs['0.5-1'] +=1
						elif 1 < h < 1.5:
							lengs['1-1.5'] +=1
						elif 1.5 < h < 2:
							lengs['1.5-2'] +=1
						elif 2 < h < 2.5:
							lengs['2-2.5'] +=1
						elif 2.5 < h < 3:
							lengs['2.5-3'] +=1
						elif 3 < h < 3.5:
							lengs['3-3.5'] +=1
						elif 3.5 < h < 4:
							lengs['3.5-4'] +=1
						elif 4 < h < 4.5:
							lengs['4-4.5'] +=1
						elif 4.5 < h < 5:
							lengs['4.5-5'] +=1
						elif h >= 5.0:
							lengs['>5'] +=1
						avg+=h
						count+=1
					except:
						lengs['less-than-5-AA'] +=1
		men_means = []

		for i in lengs:
			men_means.append(lengs[i])

		x = np.arange(len(labels))  # the label locations
		width = 0.35  # the width of the bars

		fig, ax = plt.subplots()
		rects1 = ax.bar(x - width/2, men_means, width, label='Average:'+str(round(avg/count, 3)))

		# Add some text for labels, title and custom x-axis tick labels, etc.
		ax.set_ylabel('# of Distances')
		ax.set_xlabel('Rise Distance')
		ax.set_title('Helices Data')
		ax.set_xticks(x)
		ax.set_xticklabels(labels)
		ax.legend()


		def autolabel(rects):
		    """Attach a text label above each bar in *rects*, displaying its height."""
		    for rect in rects:
		        height = rect.get_height()
		        ax.annotate('{}'.format(height),
		                    xy=(rect.get_x() + rect.get_width() / 2, height),
		                    xytext=(0, 3),  # 3 points vertical offset
		                    textcoords="offset points",
		                    ha='center', va='bottom')
		autolabel(rects1)
		fig.tight_layout()
		plt.show()

	def Q4(self):
		with open('Q4.json') as f:
			data = json.load(f)

		tick_label = {'Right-handed alpha':0, 'Right-handed omega':0, 'Right-handed pi':0, 'Right-handed gamma':0, 'Right-handed 3 - 10':0, 
			'Left-handed alpha':0, 'Left-handed omega':0,'Left-handed gamma':0,'2 - 7 ribbon/helix':0,'Polyproline':0}
		RHA = 0
		RH3 = 0
		Pol = 0
		for i in data:
			for j in i:
				if j == []:
					pass
				# print (j)
				for h in i[j]:
					tick_label[h['class']] += h['length']

					if h['class'] == 'Right-handed alpha':
						RHA+=1
					elif h['class'] == 'Right-handed 3 - 10':
						RH3+=1
					elif h['class'] == 'Polyproline':
						Pol+=1
		total = RHA + RH3 + Pol
		x = round(tick_label['Right-handed alpha']/RHA,0)
		y = round(tick_label['Right-handed 3 - 10']/RH3,0)
		z = round(tick_label['Polyproline']/Pol,0)

		Type = ['Right-handed alpha, Average:'+str(x), 
		'Right-handed 3 - 10, Average:'+str(y),
		'Polyproline, Average:'+str(z)]
		# portion covered by each label
		slices = [RHA,RH3,Pol]
		# color for each label
		colors = ['r', 'b','g']
		# plotting the pie chart
		plt.pie(slices, labels = Type, colors=colors,
		        startangle=90, shadow = False, explode = (0, 0, 0),
		        radius = 1.2, autopct = '%1.1f%%')
		# plotting legend
		plt.legend()
		# showing the plot
		plt.show()



		labels = ['Right-handed alpha', 'Right-handed omega', 'Right-handed pi', 'Right-handed gamma', 'Right-handed 3 - 10', 
			'Left-handed alpha', 'Left-handed omega','Left-handed gamma','2 - 7 ribbon/helix','Polyproline']

		men_means = [tick_label['Right-handed alpha'], tick_label['Right-handed omega'],
		tick_label['Right-handed pi'], tick_label['Right-handed gamma'], tick_label['Right-handed 3 - 10'],
		tick_label['Left-handed alpha'],tick_label['Left-handed omega'],tick_label['Left-handed gamma'],
		 tick_label['2 - 7 ribbon/helix'], tick_label['Polyproline']]
		
		women_means = [RHA,0, 0,0,RH3,0,0,0,0,Pol]

		x = np.arange(len(labels))  # the label locations
		width = 0.35  # the width of the bars

		fig, ax = plt.subplots()
		rects1 = ax.bar(x - width/2, men_means, width, label='Total Number of Amino Acids')
		rects2 = ax.bar(x + width/2, women_means, width, label='Total number of Helices')

		# Add some text for labels, title and custom x-axis tick labels, etc.
		ax.set_ylabel('class')
		ax.set_title('Alpha Helices Data')
		ax.set_xticks(x)
		ax.set_xticklabels(labels)
		ax.legend()


		def autolabel(rects):
		    """Attach a text label above each bar in *rects*, displaying its height."""
		    for rect in rects:
		        height = rect.get_height()
		        ax.annotate('{}'.format(height),
		                    xy=(rect.get_x() + rect.get_width() / 2, height),
		                    xytext=(0, 3),  # 3 points vertical offset
		                    textcoords="offset points",
		                    ha='center', va='bottom')


		autolabel(rects1)
		autolabel(rects2)

		fig.tight_layout()

		plt.show()

		return

	def Q3(self):
		with open('Q3.json') as f:
			data = json.load(f)

		Degrees = []
		Radians = []
		number = 0
		for i in data:
			for j in i:
				if j == []:
					pass
				for h in i[j]:
					Degrees.append(h['Degree'])
					# print(h['Degree'])
					Radians.append(h['Radian'])
		x = np.array(Radians)
		y = np.array(Degrees)
		plt.scatter(x, y)
		plt.show()

	def Q2(self):
		with open('Q2.json') as f:
			data = json.load(f)

		helicesPhi = []
		helicesPsi = []
		SheetsPhi = []
		SheetsPsi = []
		NonPhi = []
		NonPsi = []

		AminoAcids = ['ASN', 'LEU', 'THR', 'GLU', 'LYS', 'PRO', 'VAL', 'SER', 'ILE', 'GLY', 'ALA', 
		'ARG', 'GLN', 'ASP', 'PHE', 'HIS', 'TYR', 'TRP', 'MET', 'CYS', 'SEC', 'PYL', 'UNK']
		
		for y in AminoAcids:
			for i in data:
				for j in i:
					if j == []:
						pass
					for h in i[j]:
						if h['AAName'] == y:
							if (h['structure'] == "HELIX"):
								helicesPhi.append(h['phi'])
								helicesPsi.append(h['psi'])
							elif (h['structure'] == "SHEET"):
								SheetsPhi.append(h['phi'])
								SheetsPsi.append(h['psi'])
							else:
								NonPhi.append(h['phi'])
								NonPsi.append(h['psi'])
	
			print(y)
		# Draw a single graph
			plt.plot(NonPhi,NonPsi,'k.')
			plt.plot(helicesPhi,helicesPsi,'y.')
			plt.plot(SheetsPhi,SheetsPsi,'g.')
			plt.axis([-180,180,-180,180])
			plt.show()

			# Create a dataset:
			dfhelices=pd.DataFrame({'x': np.array(helicesPhi), 'y': np.array(helicesPsi) })
			dfsheets=pd.DataFrame({'x': np.array(SheetsPhi), 'y': np.array(SheetsPsi) })
			dfNon = pd.DataFrame({'x': np.array(NonPhi), 'y': np.array(NonPsi) })

			# show helices

			# plt.plot( 'x', 'y', color='blue',data=dfhelices, linestyle='none', marker='.')
			# Helices_patch = mpatches.Patch(color='yellow', label='Helices')
			# plt.legend(handles=[Helices_patch])
			# plt.xlabel("phi")
			# plt.ylabel("psi")
			# plt.show()

			# # show sheets
			# plt.plot( 'x', 'y', color='red', data=dfsheets, linestyle='none', marker='.')
			# Sheets_patch = mpatches.Patch(color='green', label='Sheets')
			# plt.legend(handles=[Sheets_patch])
			# plt.xlabel("phi")
			# plt.ylabel("psi")
			# plt.show()

			# # show non-structure
			# plt.plot( 'x', 'y', color='black', data=dfNon, linestyle='none', marker='.')
			# Non_patch = mpatches.Patch(color='blue', label='No Structure')
			# plt.legend(handles=[Non_patch])
			# plt.xlabel("phi")
			# plt.ylabel("psi")
			# plt.show()
			helicesPhi = []
			helicesPsi = []
			SheetsPhi = []
			SheetsPsi = []
			NonPhi = []
			NonPsi = []

test = parse()
