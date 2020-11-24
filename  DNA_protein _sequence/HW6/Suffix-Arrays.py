class SuffixArray:
    def SuffixArray(self, text, pattern):
        x = 0
        suffixs = []
        orderedSuffix = {}
        occurences = []
        for i in text:
            suffixs.append(text[x:])
            x+=1
        x=0
        for i in sorted(suffixs):
            orderedSuffix[suffixs.index(i)] = i
        for i in orderedSuffix:
            if orderedSuffix[i][0:len(pattern)] == pattern:
                occurences.append(i)
                x+=1
        print(pattern, '-->', x, 'occurences')
        for x in sorted(occurences):
            print(x, end=' ')
        print('')

# //Read text
file = open('SampleText.txt', 'r')
Text = file.readline()                 #Text = SampleText
file.close()

file2 = open('SamplePatterns.txt', 'r')
lines = file2.readlines()
file2.close()
Patterns = []                   #contains the pattern [AAA, AGCG, AGCGTA, GATA]
for line in lines:
    Patterns.append(line.replace('\n', ''))

g = SuffixArray()
g.SuffixArray(Text, Patterns[0])
print('')
g.SuffixArray(Text, Patterns[1])
print('')
g.SuffixArray(Text, Patterns[2])
g.SuffixArray(Text, Patterns[3])
