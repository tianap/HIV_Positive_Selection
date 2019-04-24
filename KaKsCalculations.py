#!/usr/bin/env python3
# Name: Tiana Pereira (tmpereir)
# Group Members: n/a

import sys
class FastaReader :
	"""
	Reads in FASTA files and yields header and sequence.

	Input: FASTA file
	Output: yields header, sequence
	"""
	
	def __init__ (self, fname=''):
		'''contructor: saves attribute fname '''
		self.fname = fname
			
	def doOpen (self):
		if self.fname is '':
			return sys.stdin
		else:
			return open(self.fname)
		
	def readFasta (self):
		
		header = ''
		sequence = ''
		
		with self.doOpen() as fileH:
			
			header = ''
			sequence = ''
			
			# skip to first fasta header
			line = fileH.readline()
			while not line.startswith('>') :
				line = fileH.readline()
			header = line[1:].rstrip()

			for line in fileH:
				if line.startswith ('>'):
					yield header,sequence
					header = line[1:].rstrip()
					sequence = ''
				else :
					sequence += ''.join(line.rstrip().split()).upper()

		yield header,sequence

class mutCharStorage:
	"""
	Categorizes and stores all mutations found between two given sequences.

	"""
	
	dnaCodonTable = {
		# DNA codon table
		# T
		'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',  # TxT
		'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',  # TxC
		'TTA': 'L', 'TCA': 'S', 'TAA': '-', 'TGA': '-',  # TxA
		'TTG': 'L', 'TCG': 'S', 'TAG': '-', 'TGG': 'W',  # TxG
		# C
		'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R',  # CxT
		'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
		'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
		'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
		# A
		'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S',  # AxT
		'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
		'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
		'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
		# G
		'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G',  # GxT
		'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
		'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
		'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'   # GxG
		}

	baseStructure = {'C': 'pyr', 'T': 'pyr', 'G': 'pur', 'A': 'pur'}
	ambigCodes = {'Y': ['C', 'T'],
				  'R': ['A', 'G'],
				  'W': ['A', 'T'],
				  'S': ['G', 'C'],
				  'K': ['T', 'G'],
				  'M': ['C', 'A'],
				  'D': ['A', 'G', 'T'],
				  'V': ['A', 'C', 'G'],
				  'H': ['A', 'C', 'T'], 
				  'B': ['C', 'G', 'T']
	}
	bases = ['A', 'C', 'T', 'G']

	def __init__(self, t0Seq, tfSeq):
		
		self.seqt0 = t0Seq
		self.seqtf = tfSeq

		# Key = Codon, value [nonsyn_t, nonsyn_v, syn_t, syn_v, possYt, possYv, possSt, possSv]
		# nonsyn_t: 
		self.mutCharDict = {}

	def findMutations(self):
		"""
		Iterates through each sequence, finds all observed and possible mutations, and categorizes them.

		The method begins by finding a mutation, then categorizing it as synonymous or nonsynonymous. 
		It then determines whether the mutation is a transition or transversion, adding one to both the observed
		and possible mutation counters in the value list for each key.
		Then the method counts all of the possible mutations and stores each count in the value list.

		It finally returns the mutation dictionary for the patient with all of the categorized mutations.
		"""
		codon = 0 # to start off
		for pos in range(0, 1143, 3):
			t0codon = self.seqt0[pos:pos+3]
			tfcodon = self.seqtf[pos:pos+3]
			possYt = 0
			possYv = 0
			possSt = 0
			possSv = 0
			codon += 1 # we want to know which codon it occurs in
			if t0codon != tfcodon:
				if t0codon in self.dnaCodonTable.keys() and tfcodon in self.dnaCodonTable.keys():
					self.mutCharDict.update({codon:[0, 0, 0, 0, 0, 0, 0, 0]}) 
					if self.dnaCodonTable[t0codon] != self.dnaCodonTable[tfcodon]: #nonsyn
						for i in range(3):
							if t0codon[i] != tfcodon[i]: #looks at bases
								if self.baseStructure[t0codon[i]] != self.baseStructure[tfcodon[i]]: # Yv
									self.mutCharDict[codon][1] += 1
									self.mutCharDict[codon][5] += 1

								else: # Yt
									self.mutCharDict[codon][0] += 1
									self.mutCharDict[codon][4] += 1
							
							self.bases.remove(tfcodon[i])
							for base in self.bases:
								possMutCodon = tfcodon[0:i] + base + tfcodon[i+1:len(tfcodon)] # replaces i place with base
								if self.dnaCodonTable[t0codon] != self.dnaCodonTable[possMutCodon]: 
									if self.baseStructure[t0codon[i]] != self.baseStructure[possMutCodon[i]]:
										self.mutCharDict[codon][5] += 1 # possYv
									else:
										self.mutCharDict[codon][4] += 1 # possYt
								else:
									if self.baseStructure[t0codon[i]] != self.baseStructure[possMutCodon[i]]: 
										self.mutCharDict[codon][7] += 1 # possSv
									else:
										self.mutCharDict[codon][6] += 1 # possSt
							self.bases.append(tfcodon[i])
					else: # synonymous
						for i in range(3):
							if t0codon[i] != tfcodon[i]: #looks at bases
								if self.baseStructure[t0codon[i]] != self.baseStructure[tfcodon[i]]: # Sv
									self.mutCharDict[codon][3] += 1
									self.mutCharDict[codon][7] += 1

								else: # St
									self.mutCharDict[codon][2] += 1
									self.mutCharDict[codon][6] += 1
							
							self.bases.remove(tfcodon[i]) # removes the original mutation
							for base in self.bases:
								possMutCodon = tfcodon[0:i] + base + tfcodon[i+1:len(tfcodon)] # replaces i place with base
								if self.dnaCodonTable[t0codon] != self.dnaCodonTable[possMutCodon]: # nonsyn
									if self.baseStructure[t0codon[i]] != self.baseStructure[possMutCodon[i]]:
										self.mutCharDict[codon][5] += 1 # possYv
									else:
										self.mutCharDict[codon][4] += 1 # possYt
								else: # syn
									if self.baseStructure[t0codon[i]] != self.baseStructure[possMutCodon[i]]: 
										self.mutCharDict[codon][7] += 1 # possSv
									else:
										self.mutCharDict[codon][6] += 1 # possSt
							self.bases.append(tfcodon[i]) # resets bases list

		return self.mutCharDict

class PatientProfile :
	"""
	Parses the header of patient FASTA files, stores mutation information.

	This class is used to create a Patient object that holds attributes such as their ID,
	an individual mutation database, and a count of each mutation's characterization.

	Within all patient headers, there is information about whether the patient took drugs, and if 
	they did, what type of drug it was. This can be utilized to find trends among administered
	drugs and associated drug resistance mutations.

	Attributes: patientID, mutations, transverse count, transition count, synonymous count, 
	nonsynonymous count.
	Methods: N/A
	"""
	def __init__ (self, header, mutDatabase):
		thisHeader = header.split('_')
		self.drugs = []
		# records all drugs taken by patient
		for i in range(6, len(thisHeader)-1):
			if len(thisHeader[i]) > 1: 
				self.drugs.append(thisHeader[i])

		
		self.patientID = thisHeader
		self.mutations = mutDatabase
		self.transverseCount = 0
		self.transitionCount = 0
		self.synCount = 0
		self.nonsynCount = 0
		
		# Saves count for t, v to be used in calculating t,v frequencies
		for mutInfo in self.mutations.values():
			if mutInfo[2] == 1 or mutInfo[4]:
				self.transverseCount += 1
			elif mutInfo[1] == 1 or mutInfo[3] == 1:
				self.transitionCount += 1

		for mutInfo in self.mutations.values():
			if mutInfo[1] == 1 or mutInfo[2] == 1:
				self.nonsynCount += 1
			elif mutInfo[3] == 0 or mutInfo[4] == 1:
				self.synCount += 1

class KaKsCalculations:
	"""
	Stores all mutation characterizations across all patients, and 
	performs the calculations for each codon's Ka/Ks ratio.

	Attributes: mutationDict, sampleSize, sequenceLen, Ytmut, Yvmut, Stmut, Svmut, t_Freq, v_Freq.
	"""
	def __init__(self, mutDict, patientSize, obsT, obsV):
		
		self.mutationDict = mutDict
		self.sampleSize = patientSize
		self.sequenceLen = 1143
		
		# all possible comb. of mutations across sample


		# observed transversion and transition mutations across all patients
		self.obsT = obsT
		self.obsV = obsV
		
		# transition/transversion frequencies
		self.t_Freq = self.obsT / (self.sequenceLen * self.sampleSize)
		self.v_Freq = self.obsV / (2*self.sequenceLen * self.sampleSize)

	def KaKsCalc(self, codonInfoList, possYt, possYv, possSt, possSv):
		"""
		Performs the KaKs ratio calculation.

		Input: All mutation characterizations found in a specified codon.
		Output: Ka/Ks ratio 
		"""
		#Key = Codon, value [nonsyn_t, nonsyn_v, syn_t, syn_v, possYt, possYv, possSt, possSv]
		yMutAtCodon = codonInfoList[0] + codonInfoList[1]
		sMutAtCodon = codonInfoList[2] + codonInfoList[3]
		KaKsval = float(0)

		# the following conditional prevents DivisionbyZero errors
		if self.t_Freq == 0 or self.v_Freq == 0 or sMutAtCodon == 0:
			pass
		else:
			numerator = yMutAtCodon/sMutAtCodon
			bottomDenom = (possSt * self.t_Freq) + (possSv * self.v_Freq)
			topDenom = (possYt*self.t_Freq) + (possYv*self.v_Freq)
			fullDenom = topDenom / bottomDenom
			KaKsval = numerator / fullDenom
		return KaKsval
		


import glob, fileinput
from collections import OrderedDict

def main():
	"""
	Goes through each patient file, collects all mutations and their characterizations,
	prints each codon whose Ka/Ks value is greater than 0.
	
	Input: FASTA files for all patients
	Output: Prints all codons that express Ka/Ks values greater than 1.

	"""
	file_list = glob.glob('*.fasta')
	patient_data = []
	myReader = FastaReader()
	patient_list = []
	all_patient_mutations = []
	total_mutation_database = {} #Key = Codon, value [nonsyn_t, nonsyn_v, syn_t, syn_v, possYt, possYv, possSt, possSv]
	obsTMut = 0
	obsVMut = 0

	t0 = 0
	tf = 1

	for file in file_list:
		myReader = FastaReader(file)
		for header, sequence in myReader.readFasta():
			patient_data.append([header, sequence])
		mutationStorage = mutCharStorage(patient_data[0][1], patient_data[1][1])
		totalMutations = mutationStorage.findMutations()
		thisPatient = PatientProfile(patient_data[1][0], totalMutations)
		all_patient_mutations.append(thisPatient)
		patient_data = [] # resets patient_data
		
		for key in thisPatient.mutations.keys(): # mutated codons
			# if the codon position has not yet been recorded in total mutation database
			mutInfo = thisPatient.mutations[key]
			if key not in total_mutation_database.keys():
				total_mutation_database.update({key:mutInfo})
			else:
				for item in range(len(mutInfo)):
					total_mutation_database[key][item] += thisPatient.mutations[key][item]
	# sorts the total mutation database
	sorted_total_mutation_database = {}
	for i in sorted(total_mutation_database.keys()):
		sorted_total_mutation_database.update({i:total_mutation_database[i]})

	#print(sorted_total_mutation_database)
	# calculating observed transition and transversion mutations:
	for key in sorted_total_mutation_database.keys():
		obsVMut += sorted_total_mutation_database[key][1] + sorted_total_mutation_database[key][3]
		obsTMut += sorted_total_mutation_database[key][0] + sorted_total_mutation_database[key][2]

	allMutKaKs = KaKsCalculations(total_mutation_database, len(file_list), obsTMut, obsVMut)

	
	for key in sorted_total_mutation_database.keys():
		codonInfo = sorted_total_mutation_database[key]
		KaKsValue = allMutKaKs.KaKsCalc(codonInfo, codonInfo[4], codonInfo[5], codonInfo[6], codonInfo[7])
		sorted_total_mutation_database[key].append(KaKsValue)
		if KaKsValue > 1:
			print("Codon "+str(key)+" KaKs value: "+str(KaKsValue))





if __name__ == '__main__':
	main()

	

