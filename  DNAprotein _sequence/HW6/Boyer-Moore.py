class BoyerMoore:

    def BoyerMooreHorspool(self,T, P):
        count = 0
        occurences = []
        m = len(P)
        n = len(T)
        if m > n:
            return -1
        jump = []
        for k in range(len(T)):
            jump.append(m)
        for k in range(m - 1):
            jump[ord(P[k])] = m - k - 1
        jump = tuple(jump)
        k = m - 1
        for letters in T:
            while k < n:
                j = m - 1; i = k
                while j >= 0 and T[i] == P[j]:
                    j -= 1; i -= 1
                if j == -1:
                    count+=1
                    occurences.append(i+1)
                k += jump[ord(T[k])]
        print(P, '-->', count, 'occurences')
        for x in occurences:
            print(x, end=' ')
        print('')
        return -1




print('')
# //Read text
file = open('SampleText.txt', 'r')
Text = file.readline()                  #Text == SampleText
file.close()

file2 = open('SamplePatterns.txt', 'r')
lines = file2.readlines()
file2.close()
Patterns = []                   #contains the pattern [AAA, AGCG, AGCGTA, GATA]
for line in lines:
    Patterns.append(line.replace('\n', ''))

g = BoyerMoore()
g.BoyerMooreHorspool(Text, Patterns[0])
print('')
g.BoyerMooreHorspool(Text, Patterns[1])
print('')
g.BoyerMooreHorspool(Text, Patterns[2])
g.BoyerMooreHorspool(Text, Patterns[3])
