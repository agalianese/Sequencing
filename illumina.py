#for loop to find the shortest sequece length 
for i in SeqIO.parse('illumina_4000.fastq', 'fastq'):
    a = len(i)
    versus = 100
    if a < versus:
        min = len(i)      
print('The shortest sequence is {}'.format(min))

#for loop to find the longest sequence
for b in SeqIO.parse('illumina_4000.fastq', 'fastq'):
    new = len(b)
    max = 0
    if new > max:
        max = len(b)      
print("The longest sequence is {}".format(max))

#how many sequences have a length shorter than the max
count = 0
other = 0
for i in SeqIO.parse('illumina_4000.fastq', 'fastq'):
    if len(i) < max:
        count +=1
    else: other += 1
print(count)
oldtotal = count + other

#calculate the average quality score per sequence and then calculate the average of those scores 
from statistics import mean
list1 = []
list2 = []
list3 = []
counter1 = 0
counter2 = 0
for record in SeqIO.parse("illumina_4000.fastq","fastq"):
    qlist = record.letter_annotations["phred_quality"]
    mu1 = (mean(qlist))
    if len(record) < 100:
        list1.append(mu1)
    else: 
        list2.append(mu1)
    list3.append(mu1)

        
#print(list1)
#print(list2)
print(mean(list1))
print(mean(list2))


#open files and being counters
fast_in = open("failed.fasta", 'w')
fast_out = open('pass.fastq', 'w')
counter1 = 0
counter2 = 0

#create for loop that will open the file and compute the averages
for record in SeqIO.parse('illumina_4000.fastq', 'fastq'):
    qlist = record.letter_annotations["phred_quality"]
    mu1 = (mean(qlist))
    if len(record) < 100:
        list1.append(mu1)
    else: 
        list2.append(mu1)

#write files, print statements and close files        
fast_in.write(str(list1))
fast_out.write(str(list2))
print("There were a total of " + str(oldtotal) + " sequences.")
print(str(other) + " sequences have max length of 100 basepairs and average quality score of " + str(mean(list2)))
print(str(count) + " sequences were shorter than 100 basepairs. Their average quality score is " + str(mean(list1)))
fast_in.close()
fast_out.close()
