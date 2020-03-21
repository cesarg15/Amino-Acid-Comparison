# Amino-Acid-Comparison
This project is based off of my first project, which tried to do something similar, but clearly didn't succeed.
In my first project, I had a lot of repetitive code, obviously because I wasn't comfortable using for loops.
This is a problem that I've tried to remedy here, as for loops constitute the important body of my code.
Everything that I've learned in this course is completely new to me as I had absolutely no experience with coding.
For this project, I spent most of my time trying to figure out the best way to achieve its goals and there were many versions of it that were completely unsuccessful.
I can now say that I understand how to deal with various types of errors that I may get in my future coding projects.
I'm not sure exactly how long these projects were supposed to be, but I do feel like my code is very small.
This in no way reflects the amount of time I spent working on this project as I can comfortably say this is the most time I've spent on any project this quarter.
Although my expertise in python is still quite limited, I feel like this project accomplishes something pretty neat and I am proud of the work I've put in.
This in no way means it's entirely finished, as I think there is a lot of room for improvement and perfection.
Even though it is being turned in, I will continue to work on it so that I may eventually be able to use this myself in an applicable setting.
As you can tell, my skills with python also have a lot of room to improve and there's plenty of stuff out there for me to learn.
I think this is a good step in that direction.

So here is a break down of the code itself and what it does:

The purpose of this code is to provide a comparison between the amino acid content of two (potentially more) organisms.
Genomic data is input as a fasta file.
I define a function in which: 
  Using the SeqIO interface from BioPython, the code will parse through the fasta file and make sure sequences are converted to the RNA alphabet.
  There is some optional code that can print the sequences and their length which I believe can help with any "quality check"     required.
  I created a large dictionary for codons which is used to conduct the counts.
  The RNA sequence for each gene id is then indexed by sets of three, elimintating the risk of overlap in my counts.
  These indexed codons are then appended to an empty list and prepared for counting.
  The next line of code then uses the keys of the codon dictionary to count each amino acid and add the total tally for that one gene id to another dictionary (using the defultdict function from collections) with the amino acids and their corresponding values.
  Then these dictionaries for each gene id are combined into one single defaultdict that represents the entire fasta file.
This function is then applied to both input files which then prints the results to the screen.
I wasn't sure what kind of statistical test to use, and in retrospect I might've tried a monte carlo simulation as well.
I decided on a paired t-test, importing stats from SciPy, even though it is most commonly used in observing changes after a certain condition is applied.
In this case, although perhaps not as crucial, I figured that a test on how similar or different the means in total amino acids between the both organisms could provide some perspective on the data being observed, and provide some control over what organisms are better suited for what is trying to be accomplished.
Finally, the dictionaries that resulted from running each file through my function is used to create a simple bar graph using numpy and matplotlib.
The count for each amino acid from both imput files are displayed side by side for easier comparison.
The x axis identifies which amino acid it is and the y axis displays the count.
The data in the image is then clearly related to its input file by the color legend in the top corner.

Although I only had one function in my project code, I could not figure out how to make it work with a unittest.
I understand this was in the rubric, but I hope that my pseudo test code is at least something.
In my "test code" I took the main code of my project and ran it with two simple and short RNA sequences.
The amino acids in each of these are easily accounted for.
This basically works to confirm that the sequences are being indexed in sets of 3 and that the code is reading them and translating the codons into the proper amino acids.
I make sure to do this for both rna sequences and populate their results into two separate dictionaries.
I then run these same dictionaries through the same code used to create my bar graph to make sure that the data is being charted correctly.

I have also supplied three fasta files which can be used with this code.
These files were downloaded directly from https://ftp.ncbi.nih.gov/genomes/ .
The three fasta files include genomes for; honey bee, aphid, and lettuce.
