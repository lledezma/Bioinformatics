class KMP:

    def KMPSearch(self,pat, txt):
        count = 0
        occurences = []
        M = len(pat)
        N = len(txt)
        lps = [0]*M
        j = 0 # index for pat[]
        self.computeLPSArray(pat, M, lps)

        i = 0 # index for txt[]
        while i < N:
            if pat[j] == txt[i]:
                i += 1
                j += 1

            if j == M:
                count+=1
                occurences.append(str(i-j))
                j = lps[j-1]
            elif i < N and pat[j] != txt[i]:
                if j != 0:
                    j = lps[j-1]
                else:
                    i += 1
        print(pat, '-->', count, 'occurences')
        for x in occurences:
            print(x, end=' ')
        print('')

    def computeLPSArray(self, pat, M, lps):
        len = 0 
        lps[0] 
        i = 1
        while i < M:
            if pat[i]== pat[len]:
                len += 1
                lps[i] = len
                i += 1
            else:
                if len != 0:
                    len = lps[len-1]
                else:
                    lps[i] = 0
                    i += 1


# //Read text
file = open('SampleText.txt', 'r')
Text = file.readline()          #Text = SampleText
file.close()

file2 = open('SamplePatterns.txt', 'r')
lines = file2.readlines()
file2.close()
Patterns = []                       #contains the pattern [AAA, AGCG, AGCGTA, GATA]
for line in lines:
    Patterns.append(line.replace('\n', ''))



g = KMP()
g.KMPSearch(Patterns[0],Text)
print('')
g.KMPSearch(Patterns[1],Text)
print('')
g.KMPSearch(Patterns[2],Text)
g.KMPSearch(Patterns[3],Text)

