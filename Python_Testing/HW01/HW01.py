import math


class Scritp:

    def __init__(self):
        print("")

    def Blinking(self):
        number_of_blinks = input("Enter the number of blinks observed: ")
        blinks = int(number_of_blinks)
        obsarvation = input("Enter the duration of the observation(in seconds): ")
        time = int(obsarvation)
        total = blinks * (86400/time)
        print("")
        print("The extrapolates to ", round(total), " blinks in a day.")


    def Payment(self):
        loan_amount = input("Enter the loan amount(in dollars): ")
        loan = float(loan_amount)
        loan_term = input("Enter the loan term(in years): ")
        months = int(loan_term)
        months = months*12
        year_rate = input("Enter the yearly interest rate(as a percentage): ")
        rate = float(year_rate)
        rate = ((rate/100)/12)
        upper = (1-(math.pow(1+rate,months*-1)))
        monthly_payment = loan/(upper / rate)
        print("")
        print("Your monthly payment will be $",round(monthly_payment, 2), sep='')
        print("After 30 years, your payments will total $",round(monthly_payment,2)*months, sep='')



    def DNA(self):
        nucleotides = ['A', 'C', 'G', 'T']
        file = open('dna.txt', 'r')
        sequence = (file.readline())
        file.close()
        for x in nucleotides:
            print("nucleotide ", x, ": ",  sequence.count(x))


    def Hamming(self, P, distance):
        file = open('dna.txt', 'r')
        sequence = (file.readline())
        file.close()
        arraycounter =  -1;
        counter2=-1
        discounter = 0
        validatecounter = 0
        array = sequence
        arrayP = P
        for x in array:
            arraycounter+=1
            counter2 = arraycounter
            counter= 0
            for y in arrayP:
                try:
                    if(array[counter2] != arrayP[counter]):
                        discounter+=1
                    counter+=1
                    counter2+=1
                except:
                    discounter =0
                    pass
            if(distance == discounter):
                validatecounter+=1
            discounter = 0
        print("Pattern", P, "occurs", validatecounter, "times")


    def Edit(self, DNA1, DNA2):
        counter = 0
        counter2 = 0
        discounter = 0
        array = DNA1
        arrayP = DNA2
        for y in arrayP:
            try:
                if(array[counter2] != arrayP[counter]):
                    discounter += 1
                counter += 1
                counter2 += 1
            except:
                pass
        print("The distance is:",discounter+(abs(len(DNA1)-len(DNA2))))



test = Scritp()
test.Blinking()
test.Payment()
test.DNA()
test.Hamming("ACT", 1)
test.Edit("ACCACTGTC", "ACGTCAG")
